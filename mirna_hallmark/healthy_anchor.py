"""Healthy (GTEx) per-arm baseline for healthy-anchored pressure.

Motivation: the canonical pressure z-term ``z(m,s) = (x - mean_tumor)/sd_tumor`` is
centred on each arm's *tumor-cohort* mean, so an arm uniformly elevated in
tumour-vs-healthy (but flat across tumours) carries ~0 z and looks inert. The
biologically meaningful signal is deviation from the **healthy** baseline. This
module supplies ``median_healthy(m)`` so the engine's ``softmax_devhealthy_logrpm``
mode can score ``x(m,s) - median_healthy(m)`` instead.

Healthy reference = **GTEx true-healthy breast** (NAT is field-effect-contaminated:
it already sits tumour-ward, so anchoring on it would cancel the very signal we
want). GTEx is cross-platform (TPM/URS vs TCGA RPM/MIMAT), so its per-arm medians
are bridged onto the TCGA scale by **quantile normalisation**: each GTEx arm is
assigned the TCGA-tumour arm-median value at the rank it occupies in the GTEx
healthy abundance distribution. This keeps the healthy abundance *ordering* (which
arms are high/low in healthy) while removing the cross-platform scale offset, so
``dev`` reflects rank-relative acquired elevation rather than an un-comparable
absolute fold-change.

Arms absent in GTEx -> no baseline -> the engine floors them at 0, i.e. their full
expression is treated as tumour-acquired (correct for tumour-induced arms such as
miR-210 / miR-301a-b). A QC column flags arms detected in NAT but not GTEx (likely
cross-platform mapping gaps that nonetheless get the floor-0 treatment under the
GTEx-only policy).

Run:
  .venv/bin/python3 -m mirna_hallmark.healthy_anchor --build
  .venv/bin/python3 -m mirna_hallmark.healthy_anchor --demo      # z-spine vs healthy-anchored
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _gtex_all_arm_matrix
from mirna_hallmark.family_normal_reference import (
    GTEX_MIRNA_TPM, GTEX_SMALLRNA_ANNOT, _gtex_breast_samples, _load_full_mirna, _split_types,
)

OUT_DIR = C.TISSUE_REFERENCE_DIR / "healthy_anchor"
ABUND_FLOOR = 0.05      # log2(PM+1) median below this = "not detected"
DETECT_THRESH = 1.0     # linear PM > 1 ...
DETECT_MIN_N = 2        # ... in >= this many samples = "detected"

# DIANA-miTED true-healthy breast reference (user-provided export). Expected: a
# tab/CSV table with a miRNA-name column (miRBase mature, e.g. 'hsa-miR-10a-5p') and a
# value column named with 'rpm' (linear) or 'log2' (log2 RPM). Download from
# http://www.microrna.gr/mited filtered to tissue=Breast, status=Healthy.
HEALTHY_REF_DIR = C.MIRNA_EXPRESSION.parent.parent / "mirna_healthy_reference"
MITED_BREAST_HEALTHY = HEALTHY_REF_DIR / "mited_breast_healthy.tsv"


def _mimat_to_arm() -> Dict[str, str]:
    from analysis.expression.mirna_target_integration import load_mimat_to_arm
    return load_mimat_to_arm(C.MIRNA_MATURE_LOCI)


def _gtex_mimat_breast_median() -> pd.Series:
    """GTEx healthy-breast median log2(TPM+1) keyed by **MIMAT accession** (canonical
    cross-platform key; bypasses display-name / locus-suffix mismatches)."""
    annot = pd.read_csv(GTEX_SMALLRNA_ANNOT, sep="\t", usecols=["id", "MIMAT"], low_memory=False)
    annot = annot.dropna(subset=["MIMAT"])
    annot = annot.loc[annot["MIMAT"].astype(str).str.startswith("MIMAT")]
    id_to_mimat = dict(zip(annot["id"].astype(str), annot["MIMAT"].astype(str)))
    breast = set(_gtex_breast_samples("SMLRNA"))
    rows = []
    for chunk in pd.read_csv(GTEX_MIRNA_TPM, sep="\t", index_col=0, chunksize=400, low_memory=False):
        keep = chunk.loc[chunk.index.astype(str).isin(id_to_mimat)]
        if not keep.empty:
            rows.append(keep)
    if not rows:
        return pd.Series(dtype=float)
    g = pd.concat(rows)
    g.index = g.index.astype(str).map(id_to_mimat)
    cols = [c for c in g.columns if c in breast]
    g = np.log2(g[cols].groupby(level=0).mean().clip(lower=0) + 1)
    return g.median(axis=1)


def _gtex_arm_median() -> pd.Series:
    """GTEx breast median in TCGA **arm-name** space (MIMAT -> arm via mimat_to_arm).
    Arm names that never resolved on the TCGA side stay as raw MIMAT and still join."""
    g = _gtex_mimat_breast_median()
    m2a = _mimat_to_arm()
    g.index = [m2a.get(m, m) for m in g.index]
    return g.groupby(level=0).mean()


def _mited_breast_median() -> pd.Series:
    """DIANA-miTED true-healthy breast median log2(RPM+1), keyed by miRBase mature name.

    miTED export is a WIDE samples x miRNA matrix (rows = healthy breast samples, columns
    = metadata + `hsa-*` RPM columns). Returns empty if the export is absent (mode holds).
    """
    if not MITED_BREAST_HEALTHY.is_file():
        return pd.Series(dtype=float)
    df = pd.read_csv(MITED_BREAST_HEALTHY, sep="\t", low_memory=False)
    df.columns = [str(c).lstrip("\ufeff") for c in df.columns]
    if "Health_state" in df.columns:  # defensive: keep only healthy rows
        df = df[df["Health_state"].astype(str).str.lower().str.contains("healthy")]
    mir_cols = [c for c in df.columns if c.startswith("hsa-")]
    if not mir_cols:
        return pd.Series(dtype=float)
    vals = df[mir_cols].apply(pd.to_numeric, errors="coerce")
    med_rpm = vals.median(axis=0, skipna=True)            # per-miRNA median RPM over healthy samples
    logv = np.log2(med_rpm.clip(lower=0) + 1)
    return logv.dropna()


def _reconcile_to_tcga(named: pd.Series, target_arms) -> pd.Series:
    """Map a miRBase-name-keyed healthy series onto TCGA arm names (exact, then
    locus-suffix-stripped match)."""
    if named.empty:
        return pd.Series(dtype=float)
    src = {str(k): float(v) for k, v in named.items()}
    out = {}
    for a in target_arms:
        if a in src:
            out[a] = src[a]
        else:
            base = _strip_locus(a)
            if base in src:
                out[a] = src[base]
    return pd.Series(out, dtype=float)


def gtex_qn_baseline(*, out_dir: Path = OUT_DIR, fallback: str = "gtex_family", force: bool = False) -> pd.Series:
    """Per-arm true-healthy baseline on the TCGA scale, anchored on **GTEx breast only**
    (346 samples, correct tissue). GTEx per-arm medians are quantile-normalised onto the
    TCGA scale (rank-relative acquired elevation, cross-platform-safe).

    GTEx v10's miRNA pipeline (893 uniquely-mappable features) zeros/omits the canonical
    member of several seed families (let-7a, miR-30a, miR-10a, miR-16, miR-17, miR-200b…)
    via multi-mapping collapse, while its seed-mates quantify fine.
    `fallback` (for those gap arms), default "gtex_family" runs two layers then floor0:
      - Layer 1 (PRIMARY) anchor-calibrated DIANA-miTED rank-transfer (`source=
        "mited_anchor"`): co-measured GTEx∩miTED arms = anchors define a monotonic
        miTED->GTEx quantile map; a gap arm's miTED *rank* (not absolute RPM) is mapped to
        a GTEx-equivalent level. Arm-specific and validated (LOO vs real GTEx: r=0.90,
        median err 0.86 log2). Blood markers excluded.
      - Layer 2 (FALLBACK) GTEx-internal seed-family (`source="gtex_family"`): for gap arms
        absent from miTED but with a GTEx seed-mate -> MAX seed-mate (dominant member).
      - Combined coverage 182->213 / 222 (96%) at >=15 RPM.
      - "nat": DEPRECATED (field-effect; audit only). "none": floor0.
    Arms with no GTEx, no miTED, and no seed-mate (or blood) stay floor0 (true orphans,
    e.g. low-abundance housekeeping/blood -> conservative full-acquired credit).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "gtex_qn_baseline.tsv"
    if cache.is_file() and not force:
        df = pd.read_csv(cache, sep="\t", index_col=0)
        return df["healthy_baseline_tcga"].dropna()

    split = _split_types(_load_full_mirna())
    tumor, nat = split["tumor"], split["nat"]
    tcga_med = tumor.median(axis=1)
    nat_med = nat.median(axis=1)
    nat_detect = (nat > DETECT_THRESH).sum(axis=1)
    tumor_detect = (tumor > DETECT_THRESH).sum(axis=1)
    tumor_prev = tumor_detect / max(tumor.shape[1], 1)
    gtex_raw = _gtex_arm_median()
    gtex = gtex_raw[gtex_raw > ABUND_FLOOR].copy()   # reliably-measured GTEx arms
    gtex_reliable = gtex.copy()                       # snapshot (anchor fit uses real arms only)

    # Gap arms = any arm missing from reliable GTEx. No tumor-expression gate: the
    # anchor-miTED transfer is validated (r=0.90) and abundance-weighting makes low arms
    # harmless, so we impute every gap arm with a healthy reference (completeness); arms
    # with neither miTED nor a seed-mate still fall to floor0.
    gap_targets = [a for a in tcga_med.index if a not in gtex.index]

    imp_family, imp_mited = set(), set()
    if fallback == "gtex_family":
        # Layer 1 (PRIMARY): anchor-calibrated miTED rank-transfer. Arm-specific and
        # validated (leave-one-out vs real GTEx: Pearson 0.90, median err 0.86 log2).
        # Preferred over seed-family, which assigns a degenerate family-MAX to every
        # collapsed member (e.g. all let-7 arms -> let-7b's level; miR-103a -> its lone
        # low seed-mate miR-107). Blood markers are excluded inside the transfer.
        mited = _mited_breast_median()
        if not mited.empty:
            mited_arm = _reconcile_to_tcga(mited, tcga_med.index.union(gtex_raw.index))
            mited_arm = mited_arm[mited_arm > ABUND_FLOOR]
            pred = _mited_anchor_transfer(gtex_reliable, mited_arm, gap_targets)
            for a, v in pred.items():
                gtex.loc[a] = float(v)
                imp_mited.add(a)
        # Layer 2 (FALLBACK): GTEx-internal seed-family for gap arms NOT recoverable from
        # miTED (absent or blood-excluded) but with a GTEx-reliable seed-mate. Collapsed
        # arms are the dominant member -> MAX seed-mate.
        fam_map = _family_map()
        fam_levels: Dict[str, list] = {}
        for a, v in gtex_reliable.items():
            fam = _fam_of(a, fam_map)
            if fam:
                fam_levels.setdefault(fam, []).append(float(v))
        for a in gap_targets:
            if a in imp_mited:
                continue
            mates = fam_levels.get(_fam_of(a, fam_map))
            if mates:
                gtex.loc[a] = float(np.max(mates))
                imp_family.add(a)

    ref = np.sort(tcga_med.dropna().values.astype(float))
    gpct = gtex.rank(pct=True)
    base_gtex = pd.Series(np.quantile(ref, gpct.values), index=gtex.index)

    idx = tcga_med.index.union(gtex.index)
    baseline = base_gtex.reindex(idx)
    source = pd.Series(np.where(baseline.notna(), "gtex", None), index=idx)
    if imp_family:
        source.loc[sorted(imp_family)] = "gtex_family"
    if imp_mited:
        source.loc[sorted(imp_mited)] = "mited_anchor"
    print(f"[healthy_anchor] imputed: family={len(imp_family)} (GTEx seed-mate), "
          f"mited_anchor={len(imp_mited)} (rank-transfer).")

    if fallback == "nat":   # DEPRECATED (field-effect, audit only)
        need = baseline.isna() & (nat_detect.reindex(idx).fillna(0) >= DETECT_MIN_N)
        baseline.loc[need] = nat_med.reindex(idx).loc[need]
        source.loc[need] = "nat_fallback"
    source.loc[baseline.isna()] = "floor0"

    table = pd.DataFrame({
        "tcga_tumor_median": tcga_med.reindex(idx),
        "nat_median": nat_med.reindex(idx),
        "nat_detect_n": nat_detect.reindex(idx),
        "tumor_prev": tumor_prev.reindex(idx),
        "gtex_median": gtex.reindex(idx),
        "gtex_pct": gpct.reindex(idx),
        "healthy_baseline_tcga": baseline,
        "source": source,
    })
    # NAT-fallback arms that are tumor-UP vs NAT: NAT (field-effect-contaminated) under-
    # credits their acquired pressure -> flag (interpret dev as a conservative lower bound).
    table["nat_under_credit"] = (
        table["source"].eq("nat_fallback")
        & ((table["tcga_tumor_median"] - table["nat_median"]) > 0.5)
    )
    # floor0 arms with real tumor expression DO get full acquired credit (correct for a
    # tumor-induced arm with no healthy reference); the rest are negligible.
    table["floor0_tumor_expressed"] = table["source"].eq("floor0") & (table["tumor_prev"] >= 0.10)
    # imputed arms carry a proxy baseline (seed-mate or miTED rank-transfer) -> lower
    # precision, flag for audit / optional downstream exclusion
    table["imputed"] = table["source"].isin(["gtex_family", "mited_anchor"])
    table.index.name = "arm"
    table.to_csv(cache, sep="\t")

    # expression-threshold coverage table (how the GTEx-covered / gap split moves with the
    # tumor-expression cutoff used to scope "arms worth a healthy baseline")
    g_ok = set(gtex_reliable.index)
    g_ok_imp = set(gtex.index)  # incl. seed-family + miTED-anchor imputed
    cov_rows = []
    for cut in [0, 1, 2, 3, 4, 5, 6, 8]:
        sel = set(tcga_med[tcga_med > cut].dropna().index)
        cov_rows.append({
            "log2_cutoff": cut, "approx_rpm": round(2 ** cut - 1, 1), "n_arms": len(sel),
            "gtex_reliable": len(sel & g_ok), "gtex_plus_imputed": len(sel & g_ok_imp),
            "gap_after_impute": len(sel - g_ok_imp),
        })
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(out_dir / "expression_threshold_coverage.tsv", sep="\t", index=False)

    vc = source.value_counts()
    bl = baseline.dropna()
    print(f"[healthy_anchor] baseline arms={len(bl)} range [{bl.min():.2f},{bl.max():.2f}]  "
          f"sources: gtex={int(vc.get('gtex',0))}, gtex_family={int(vc.get('gtex_family',0))}, "
          f"mited_anchor={int(vc.get('mited_anchor',0))}, nat_fallback={int(vc.get('nat_fallback',0))}, "
          f"floor0={int(vc.get('floor0',0))} ({int(table['floor0_tumor_expressed'].sum())} tumor-expressed→full acquired credit)")
    return baseline.dropna()


def load_baseline(*, out_dir: Path = OUT_DIR) -> pd.Series:
    return gtex_qn_baseline(out_dir=out_dir)


_LOCUS_SUFFIX = re.compile(r"\.\d+$")


def _strip_locus(arm: str) -> str:
    """TCGA multi-locus arm 'hsa-miR-101-3p.1' -> mirBase mature 'hsa-miR-101-3p'."""
    return _LOCUS_SUFFIX.sub("", str(arm))


# Blood/platelet/RBC-contamination markers: excluded from the miTED anchor fit and from
# miTED-imputation targets (miTED breast n=12 is blood-contaminated — miR-451a ranks #1 —
# and these are not breast-epithelial miRNAs, so a miTED-derived breast baseline is wrong).
BLOOD_MIRNAS = {
    "hsa-miR-451a", "hsa-miR-486-5p", "hsa-miR-486-3p", "hsa-miR-144-3p", "hsa-miR-144-5p",
    "hsa-miR-92a-3p", "hsa-miR-25-3p", "hsa-miR-363-3p", "hsa-miR-223-3p", "hsa-miR-142-3p",
    "hsa-miR-142-5p", "hsa-miR-150-5p", "hsa-miR-20b-5p",
}


def _mited_anchor_transfer(gtex_lvl: pd.Series, mited_lvl: pd.Series,
                           targets) -> pd.Series:
    """Anchor-calibrated rank transfer (user idea): use arms measured in BOTH GTEx and
    miTED as anchors to learn the monotonic miTED->GTEx mapping (quantile-quantile), then
    predict a gap arm's GTEx-equivalent healthy level from its miTED value. Uses miTED only
    for *relative position among reference miRNAs* (platform-robust), not absolute RPM.
    """
    if mited_lvl.empty:
        return pd.Series(dtype=float)
    anchors = [a for a in gtex_lvl.index.intersection(mited_lvl.index) if a not in BLOOD_MIRNAS]
    if len(anchors) < 30:
        return pd.Series(dtype=float)
    xm = np.sort(mited_lvl.loc[anchors].values.astype(float))   # sorted miTED anchor levels
    yg = np.sort(gtex_lvl.loc[anchors].values.astype(float))    # sorted GTEx anchor levels
    out = {}
    for a in targets:
        if a in mited_lvl.index and a not in BLOOD_MIRNAS:
            out[a] = float(np.interp(float(mited_lvl[a]), xm, yg))  # QQ map at matched quantile
    return pd.Series(out, dtype=float)


def _family_map() -> Dict[str, str]:
    """miRBase mature ID -> TargetScan seed family (human)."""
    df = pd.read_csv(C.MIR_FAMILY_INFO, sep="\t")
    h = df.loc[df["Species ID"].eq(9606), ["miR family", "MiRBase ID"]].dropna()
    return dict(zip(h["MiRBase ID"].astype(str), h["miR family"].astype(str)))


def _fam_of(arm: str, fam_map: Dict[str, str]) -> Optional[str]:
    base = _strip_locus(arm)
    return fam_map.get(base, fam_map.get(str(arm)))


def _harmonize_to_gtex(tcga_arms, gtex_arms) -> Dict[str, str]:
    """TCGA arm -> GTEx arm name. exact match, else locus-suffix-stripped match."""
    g = set(gtex_arms)
    out: Dict[str, str] = {}
    for a in tcga_arms:
        if a in g:
            out[a] = a
        else:
            base = _strip_locus(a)
            if base != a and base in g:
                out[a] = base
    return out


def _n_expressed(mat: pd.DataFrame, arm: str, *, thresh: float = 1.0, min_samples: int = 2) -> int:
    """# samples with log2(PM+1) > thresh (i.e. linear PM > 1) for this arm."""
    if mat is None or arm not in mat.index:
        return 0
    v = pd.to_numeric(mat.loc[arm], errors="coerce")
    return int((v > thresh).sum())


def gap_audit(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """Quantify the GTEx-unmapped/NAT-present gap: name recovery + real-expression triage."""
    out_dir.mkdir(parents=True, exist_ok=True)
    split = _split_types(_load_full_mirna())
    tumor, nat = split["tumor"], split["nat"]
    gtex = _gtex_all_arm_matrix()

    tcga_arms = sorted(set(tumor.index) | set(nat.index))
    gtex_arms = set(gtex.index)
    hmap = _harmonize_to_gtex(tcga_arms, gtex_arms)

    exact = {a for a, g in hmap.items() if a == g}
    recovered = {a for a, g in hmap.items() if a != g}
    unmapped = [a for a in tcga_arms if a not in hmap]

    print("=== GTEx coverage of TCGA arms ===")
    print(f"TCGA arms (tumor∪NAT): {len(tcga_arms)}   GTEx arms: {len(gtex_arms)}")
    print(f"  exact-name matched : {len(exact)}")
    print(f"  recovered by stripping locus suffix (.N): {len(recovered)}")
    print(f"  still unmapped     : {len(unmapped)}")
    if recovered:
        ex = sorted(recovered)[:8]
        print(f"  recovered examples : {', '.join(ex)}")

    rows = []
    for a in unmapped:
        nt = _n_expressed(tumor, a)
        nn = _n_expressed(nat, a)
        ng = _n_expressed(gtex, a)  # 0 by construction (unmapped) but kept for symmetry
        real_types = int(nt >= 2) + int(nn >= 2) + int(ng >= 2)
        rows.append({
            "arm": a, "n_expr_tumor": nt, "n_expr_nat": nn, "n_expr_gtex": ng,
            "real_in_tcga": (nt >= 2) or (nn >= 2), "n_real_types": real_types,
            "tumor_median": float(pd.to_numeric(tumor.loc[a], errors="coerce").median()) if a in tumor.index else np.nan,
            "nat_median": float(pd.to_numeric(nat.loc[a], errors="coerce").median()) if a in nat.index else np.nan,
        })
    gap = pd.DataFrame(rows)
    real = gap.loc[gap["real_in_tcga"]] if not gap.empty else gap
    print("\n=== of the still-unmapped arms, how many are genuinely expressed? ===")
    print("(real = linear PM > 1 in >= 2 samples)")
    print(f"  unmapped total              : {len(gap)}")
    print(f"  real in TCGA (tumor OR NAT) : {len(real)}  ({100*len(real)/max(len(gap),1):.0f}%)")
    if not gap.empty:
        print(f"  real in tumor (>=2 smpl)    : {int((gap['n_expr_tumor']>=2).sum())}")
        print(f"  real in NAT   (>=2 smpl)    : {int((gap['n_expr_nat']>=2).sum())}")
        top = real.sort_values("tumor_median", ascending=False).head(12)
        print("\n  most-abundant unmapped-but-real arms (would be floored to 0):")
        print(top[["arm", "n_expr_tumor", "n_expr_nat", "tumor_median", "nat_median"]].to_string(index=False))
    gap.to_csv(out_dir / "gtex_gap_audit.tsv", sep="\t", index=False)
    print(f"\n[healthy_anchor] wrote {out_dir / 'gtex_gap_audit.tsv'}")
    return gap


def _demo(out_dir: Path = OUT_DIR) -> None:
    """Compare tumor pressure-budget arm shares: z-spine vs healthy-anchored."""
    from mirna_hallmark.hallmark_sets import HallmarkSets
    from mirna_hallmark.pressure_build import compute_gene_pressure_contributions, load_mirtar_edges

    base = gtex_qn_baseline(out_dir=out_dir)
    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)

    def _share(mode: str, **kw) -> pd.Series:
        c = compute_gene_pressure_contributions(genes, edges=edges, expr_mode=mode,  # type: ignore[arg-type]
                                                resolve_arms=False, **kw)
        m = c.groupby("miRNA")["mean_abs_contribution"].sum()
        return m / m.sum()

    z_share = _share("softmax_z_logrpm").rename("share_zspine")
    h_share = _share("softmax_devhealthy_logrpm", healthy_baseline=base).rename("share_healthy")
    df = pd.concat([z_share, h_share], axis=1).fillna(0.0)
    df["delta"] = df["share_healthy"] - df["share_zspine"]
    tab = pd.read_csv(out_dir / "gtex_qn_baseline.tsv", sep="\t", index_col=0)
    df["baseline_source"] = df.index.map(tab["source"]).fillna("floor0")
    df = df.sort_values("delta", ascending=False)

    print("\n=== pressure-budget share: healthy-anchored vs z-spine (tumor) ===")
    print("Top GAINERS under healthy anchor (tumor-acquired / healthy-low):")
    print(df.head(12)[["share_zspine", "share_healthy", "delta", "baseline_source"]].to_string())
    print("\nTop LOSERS under healthy anchor (constitutively high in healthy):")
    print(df.tail(12)[["share_zspine", "share_healthy", "delta", "baseline_source"]].to_string())
    df.to_csv(out_dir / "demo_share_zspine_vs_healthy.tsv", sep="\t")
    print(f"\n[healthy_anchor] wrote {out_dir / 'demo_share_zspine_vs_healthy.tsv'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build", action="store_true", help="compute + cache the GTEx QN baseline")
    ap.add_argument("--demo", action="store_true", help="compare z-spine vs healthy-anchored arm shares")
    ap.add_argument("--gap-audit", action="store_true", help="name-recovery + real-expression triage of the GTEx gap")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.gap_audit:
        gap_audit(out_dir=args.out_dir)
        return
    if args.build or not args.demo:
        gtex_qn_baseline(out_dir=args.out_dir, force=args.force)
    if args.demo:
        _demo(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
