"""Per-subtype CN depth: sample distributions, cross-subtype ranges, neighborhood co-alteration.

Focuses on flagged / curated miRNA arms (sections 1–2 of genome-map report). Uses
participant-level CN with MIMAT aggregation via cohort ``w_norm`` weights.

Outputs under ``output/mirna_locus_cnv/maps/``:
  - ``mirna_cnv_focus_subtype_distribution.tsv``
  - ``mirna_cnv_focus_subtype_comparison.tsv``
  - ``mirna_cnv_focus_pairwise_subtype.tsv``
  - ``mirna_cnv_neighborhood_by_subtype.tsv``
  - ``mirna_cnv_neighborhood_summary.tsv``  (one row per anchor × neighbor × subtype)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps import (
    CURATED_ARMS,
    DIPLOID_CN,
    PAM50_ORDER,
    TUMOR_PAM50,
    _load_arm_labels,
    _load_locus_catalog,
    _load_weight_reference,
)

# Arms from genome-map sections 1 (extreme amp) and 2 (subtype biology).
FOCUS_ARMS: Tuple[str, ...] = (
    "hsa-miR-21-5p",
    "hsa-miR-21-3p",
    "hsa-miR-151a-3p",
    "hsa-miR-151a-5p",
    "hsa-miR-4661-5p",
    "hsa-miR-4661-3p",
    "hsa-miR-30b-5p",
    "hsa-miR-30b-3p",
    "hsa-miR-30d-5p",
    "hsa-miR-30d-3p",
    "hsa-miR-135b-5p",
    "hsa-miR-135b-3p",
    "hsa-miR-205-5p",
    "hsa-miR-205-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-10a-5p",
    "hsa-miR-10a-3p",
    "hsa-miR-454-3p",
    "hsa-miR-454-5p",
    "hsa-miR-301a-3p",
    "hsa-miR-301a-5p",
    "hsa-let-7a-5p",
)

NEIGHBORHOOD_WINDOW_BP = 2_000_000
PAIRWISE_SUBTYPES = [
    ("LumB", "LumA"),
    ("Her2", "LumA"),
    ("Her2", "Normal"),
    ("LumB", "Normal"),
    ("Basal", "LumA"),
    ("Basal", "Normal"),
]


def _median_bootstrap_ci(
    values: np.ndarray,
    *,
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 42,
) -> Tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.array(
        [np.median(rng.choice(x, size=len(x), replace=True)) for _ in range(n_boot)]
    )
    return float(np.percentile(boots, 100 * alpha / 2)), float(np.percentile(boots, 100 * (1 - alpha / 2)))


def _distribution_row(
    cn: pd.Series,
    states: Optional[pd.Series],
    *,
    arm_label: str,
    mimat: str,
    pam50: str,
    healthy_ref: float,
) -> dict:
    v = pd.to_numeric(cn, errors="coerce").dropna()
    n = int(len(v))
    if n == 0:
        return {}
    ci_lo, ci_hi = _median_bootstrap_ci(v.values)
    st = states.reindex(v.index).astype(str) if states is not None else pd.Series(dtype=str)
    row = {
        "arm_label": arm_label,
        "mimat": mimat,
        "pam50": pam50,
        "n_participants": n,
        "healthy_ref_cn": round(healthy_ref, 4),
        "mean_cn": round(float(v.mean()), 4),
        "median_cn": round(float(v.median()), 4),
        "std_cn": round(float(v.std(ddof=0)), 4) if n > 1 else 0.0,
        "min_cn": round(float(v.min()), 4),
        "max_cn": round(float(v.max()), 4),
        "cn_range": f"{int(v.min())}-{int(v.max())}",
        "q25_cn": round(float(v.quantile(0.25)), 4),
        "q75_cn": round(float(v.quantile(0.75)), 4),
        "iqr_cn": round(float(v.quantile(0.75) - v.quantile(0.25)), 4),
        "median_ci_lo": round(ci_lo, 4) if np.isfinite(ci_lo) else np.nan,
        "median_ci_hi": round(ci_hi, 4) if np.isfinite(ci_hi) else np.nan,
        "delta_median_vs_healthy": round(float(v.median()) - healthy_ref, 4),
    }
    if len(st):
        row.update(
            {
                "pct_diploid": round(100.0 * (v == 2).mean(), 2),
                "pct_gain": round(100.0 * (st == "gain").mean(), 2),
                "pct_amp": round(100.0 * (st == "amp").mean(), 2),
                "pct_loss": round(100.0 * (st == "loss").mean(), 2),
                "pct_deep_del": round(100.0 * (st == "deep_del").mean(), 2),
                "pct_altered": round(100.0 * st.isin(["gain", "amp", "loss", "deep_del"]).mean(), 2),
            }
        )
    return row


def _participant_mimat_cn(
    locus_long: pd.DataFrame,
    weight_ref: pd.DataFrame,
    mimat: str,
) -> pd.DataFrame:
    """Participant × PAM50 table with expression-weighted MIMAT CN."""
    w = weight_ref.loc[weight_ref["mimat"] == mimat].copy()
    if w.empty:
        return pd.DataFrame()
    lids = set(w["locus_id"].astype(str))
    sub = locus_long.loc[locus_long["entity_id"].astype(str).isin(lids)].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = sub.drop_duplicates(subset=["participant", "entity_id"], keep="first")
    sub = sub.merge(w[["locus_id", "w_norm"]], left_on="entity_id", right_on="locus_id", how="left")
    sub["copy_number"] = pd.to_numeric(sub["copy_number"], errors="coerce")
    sub["w_norm"] = pd.to_numeric(sub["w_norm"], errors="coerce").fillna(0.0)

    rows: List[dict] = []
    for pid, grp in sub.groupby("participant"):
        cn = grp["copy_number"].values
        wt = grp["w_norm"].values
        m = np.isfinite(cn)
        if not m.any():
            continue
        cn_m = cn[m]
        wt_m = wt[m]
        if wt_m.sum() <= 0:
            val = float(np.median(cn_m))
        else:
            val = float(np.average(cn_m, weights=wt_m))
        rows.append(
            {
                "participant": pid,
                "PAM50_final": grp["PAM50_final"].iloc[0],
                "copy_number": val,
                "cn_state": _classify_cn(val),
            }
        )
    return pd.DataFrame(rows)


def _classify_cn(cn: float) -> str:
    from pipeline.CNV.features import _classify_cn_state

    return str(_classify_cn_state(cn))


def build_focus_participant_cn(
    long: pd.DataFrame,
    focus_arms: Sequence[str] = FOCUS_ARMS,
) -> pd.DataFrame:
    """Long participant-level weighted CN for focus arms."""
    weight_ref = _load_weight_reference()
    arm_labels = _load_arm_labels()
    label_to_mimat = arm_labels.set_index("arm_label")["mimat"].to_dict()
    locus_long = long.loc[long["entity_level"] == "locus"].copy()

    parts: List[pd.DataFrame] = []
    for arm in focus_arms:
        mimat = label_to_mimat.get(arm)
        if not mimat:
            continue
        pdf = _participant_mimat_cn(locus_long, weight_ref, str(mimat))
        if pdf.empty:
            continue
        pdf["arm_label"] = arm
        pdf["mimat"] = mimat
        healthy_ref = DIPLOID_CN
        pdf["healthy_ref_cn"] = healthy_ref
        pdf["delta_vs_healthy"] = pdf["copy_number"] - healthy_ref
        parts.append(pdf)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def build_subtype_distribution(participant_cn: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    for arm, grp in participant_cn.groupby("arm_label"):
        mimat = grp["mimat"].iloc[0]
        href = float(grp["healthy_ref_cn"].iloc[0])
        for pam in PAM50_ORDER:
            sub = grp.loc[grp["PAM50_final"].astype(str) == pam]
            if sub.empty:
                continue
            row = _distribution_row(
                sub.set_index("participant")["copy_number"],
                sub.set_index("participant")["cn_state"],
                arm_label=str(arm),
                mimat=str(mimat),
                pam50=pam,
                healthy_ref=href,
            )
            if row:
                rows.append(row)
    return pd.DataFrame(rows)


def build_subtype_comparison(dist: pd.DataFrame, participant_cn: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    for arm, grp in participant_cn.groupby("arm_label"):
        mimat = grp["mimat"].iloc[0]
        href = float(grp["healthy_ref_cn"].iloc[0])
        by_pam = {
            pam: pd.to_numeric(grp.loc[grp["PAM50_final"] == pam, "copy_number"], errors="coerce").dropna().values
            for pam in PAM50_ORDER
        }
        samples = [v for v in by_pam.values() if len(v) >= 5]
        kw_stat, kw_p = (np.nan, np.nan)
        if len(samples) >= 2:
            try:
                kw_stat, kw_p = kruskal(*samples)
            except ValueError:
                pass
        dsub = dist.loc[dist["arm_label"] == arm]
        row = {
            "arm_label": arm,
            "mimat": mimat,
            "healthy_ref_cn": href,
            "kruskal_stat": round(float(kw_stat), 4) if np.isfinite(kw_stat) else np.nan,
            "kruskal_p": kw_p,
        }
        medians = []
        for pam in PAM50_ORDER:
            ps = dsub.loc[dsub["pam50"] == pam]
            if ps.empty:
                row[f"median_cn_{pam}"] = np.nan
                row[f"cn_range_{pam}"] = ""
                row[f"n_{pam}"] = 0
                row[f"median_ci_{pam}"] = ""
            else:
                r = ps.iloc[0]
                row[f"median_cn_{pam}"] = r["median_cn"]
                row[f"cn_range_{pam}"] = r["cn_range"]
                row[f"n_{pam}"] = int(r["n_participants"])
                lo, hi = r.get("median_ci_lo", np.nan), r.get("median_ci_hi", np.nan)
                row[f"median_ci_{pam}"] = (
                    f"[{lo:.2f},{hi:.2f}]" if pd.notna(lo) and pd.notna(hi) else ""
                )
                medians.append(float(r["median_cn"]))
        if medians:
            row["median_spread_subtypes"] = round(max(medians) - min(medians), 4)
            row["max_median_subtype"] = PAM50_ORDER[medians.index(max(medians))]
            row["min_median_subtype"] = PAM50_ORDER[medians.index(min(medians))]
        rows.append(row)
    return pd.DataFrame(rows)


def build_pairwise_subtype(participant_cn: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    for arm, grp in participant_cn.groupby("arm_label"):
        for a, b in PAIRWISE_SUBTYPES:
            va = pd.to_numeric(grp.loc[grp["PAM50_final"] == a, "copy_number"], errors="coerce").dropna()
            vb = pd.to_numeric(grp.loc[grp["PAM50_final"] == b, "copy_number"], errors="coerce").dropna()
            if len(va) < 5 or len(vb) < 5:
                continue
            try:
                stat, p = mannwhitneyu(va, vb, alternative="two-sided")
            except ValueError:
                continue
            rows.append(
                {
                    "arm_label": arm,
                    "mimat": grp["mimat"].iloc[0],
                    "pam50_a": a,
                    "pam50_b": b,
                    "n_a": len(va),
                    "n_b": len(vb),
                    "median_cn_a": round(float(va.median()), 4),
                    "median_cn_b": round(float(vb.median()), 4),
                    "cn_range_a": f"{int(va.min())}-{int(va.max())}",
                    "cn_range_b": f"{int(vb.min())}-{int(vb.max())}",
                    "median_diff_a_minus_b": round(float(va.median() - vb.median()), 4),
                    "mw_p": p,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    from mirna_hallmark import stats as S

    out["mw_q"] = S.bh_fdr(out["mw_p"].values)
    return out.sort_values(["arm_label", "mw_p"])


def _loci_within_window(catalog: pd.DataFrame, anchor_locus: str, window_bp: int) -> pd.DataFrame:
    meta = catalog.loc[catalog["locus_id"].astype(str) == anchor_locus]
    if meta.empty:
        return pd.DataFrame()
    row = meta.iloc[0]
    chrom = row["chrom"]
    mid = float(row["midpoint"])
    sub = catalog.loc[catalog["chrom"] == chrom].copy()
    sub["dist_bp"] = (sub["midpoint"] - mid).abs()
    sub = sub.loc[(sub["dist_bp"] <= window_bp) & (sub["locus_id"].astype(str) != anchor_locus)]
    return sub.sort_values("dist_bp")


def _anchor_locus_for_arm(arm: str, arm_labels: pd.DataFrame, catalog: pd.DataFrame) -> Optional[str]:
    """Dominant locus for arm: highest w_norm, else first locus on same precursor."""
    mimat = arm_labels.loc[arm_labels["arm_label"] == arm, "mimat"]
    if mimat.empty:
        return None
    mimat = str(mimat.iloc[0])
    w = _load_weight_reference()
    w = w.loc[w["mimat"] == mimat].sort_values("w_norm", ascending=False)
    if not w.empty:
        return str(w.iloc[0]["locus_id"])
    if "mature_accessions" in catalog.columns:
        hit = catalog.loc[catalog["mature_accessions"].astype(str).str.contains(mimat, na=False)]
        if not hit.empty:
            return str(hit.iloc[0]["locus_id"])
    return None


def _load_locus_catalog_with_accessions() -> pd.DataFrame:
    catalog = _load_locus_catalog()
    raw = pd.read_csv(C.MIRNA_PRECURSOR_LOCI)
    acc = raw[["gene_id", "mature_accessions"]].rename(columns={"gene_id": "locus_id"})
    return catalog.merge(acc, on="locus_id", how="left")


def build_neighborhood_by_subtype(
    long: pd.DataFrame,
    focus_arms: Sequence[str] = FOCUS_ARMS,
    *,
    window_bp: int = NEIGHBORHOOD_WINDOW_BP,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Neighborhood co-alteration for anchor loci: same ASCAT3 segment + ±window miRNA loci.

    Returns (long per subtype, summary collapsed).
    """
    catalog = _load_locus_catalog_with_accessions()
    if "mature_accessions" not in catalog.columns:
        raw = pd.read_csv(C.MIRNA_PRECURSOR_LOCI)
        acc = raw[["gene_id", "mature_accessions"]].rename(columns={"gene_id": "locus_id"})
        catalog = catalog.merge(acc, on="locus_id", how="left")
    arm_labels = _load_arm_labels()
    loc = long.loc[long["entity_level"] == "locus"].drop_duplicates(
        subset=["participant", "entity_id"], keep="first"
    ).copy()
    loc["copy_number"] = pd.to_numeric(loc["copy_number"], errors="coerce")
    loc["entity_id"] = loc["entity_id"].astype(str)

    long_rows: List[dict] = []
    for arm in focus_arms:
        anchor_locus = _anchor_locus_for_arm(arm, arm_labels, catalog)
        if not anchor_locus:
            continue
        anchor_meta = catalog.loc[catalog["locus_id"] == anchor_locus].iloc[0]
        anchor_cn = loc.loc[loc["entity_id"] == anchor_locus, ["participant", "PAM50_final", "copy_number", "cn_state", "overlap_segment"]]
        anchor_cn = anchor_cn.rename(columns={"copy_number": "anchor_cn", "cn_state": "anchor_state", "overlap_segment": "anchor_segment"})

        # --- same-segment miRNA neighbors (per participant) ---
        seg_neighbors = (
            loc.merge(anchor_cn[["participant", "anchor_segment"]], on="participant")
            .loc[lambda d: (d["overlap_segment"] == d["anchor_segment"]) & (d["entity_id"] != anchor_locus)]
        )
        seg_neighbor_ids = seg_neighbors.groupby("entity_id").size().sort_values(ascending=False).head(12).index.tolist()

        # --- window neighbors (static genomic) ---
        win = _loci_within_window(catalog, anchor_locus, window_bp)
        win_ids = win["locus_id"].astype(str).tolist()[:20]

        neighbor_specs: List[Tuple[str, str, str, float]] = []
        for nid in seg_neighbor_ids:
            nm = catalog.loc[catalog["locus_id"] == nid, "locus_name"]
            neighbor_specs.append(("same_segment", str(nid), str(nm.iloc[0]) if len(nm) else nid, np.nan))
        for _, wr in win.iterrows():
            nid = str(wr["locus_id"])
            if any(nid == x[1] for x in neighbor_specs):
                continue
            neighbor_specs.append(("window_2mb", nid, str(wr["locus_name"]), float(wr["dist_bp"])))

        for ntype, nid, nlabel, dist in neighbor_specs:
            ncn = loc.loc[loc["entity_id"] == nid, ["participant", "PAM50_final", "copy_number", "cn_state"]]
            ncn = ncn.rename(columns={"copy_number": "neighbor_cn", "cn_state": "neighbor_state"})
            merged = anchor_cn.merge(ncn, on=["participant", "PAM50_final"], how="inner")
            if merged.empty:
                continue
            for pam, sub in merged.groupby("PAM50_final"):
                sub = sub.dropna(subset=["anchor_cn", "neighbor_cn"])
                if len(sub) < 5:
                    continue
                a = sub["anchor_cn"].values
                b = sub["neighbor_cn"].values
                rho, rho_p = spearmanr(a, b) if len(sub) >= 10 else (np.nan, np.nan)
                long_rows.append(
                    {
                        "anchor_arm": arm,
                        "anchor_locus_id": anchor_locus,
                        "anchor_locus_name": anchor_meta["locus_name"],
                        "anchor_chrom": anchor_meta["chrom"],
                        "neighbor_type": ntype,
                        "neighbor_locus_id": nid,
                        "neighbor_locus_name": nlabel,
                        "neighbor_dist_bp": round(dist, 0) if np.isfinite(dist) else np.nan,
                        "pam50": str(pam),
                        "n_participants": len(sub),
                        "anchor_median_cn": round(float(np.median(a)), 4),
                        "anchor_cn_range": f"{int(np.min(a))}-{int(np.max(a))}",
                        "neighbor_median_cn": round(float(np.median(b)), 4),
                        "neighbor_cn_range": f"{int(np.min(b))}-{int(np.max(b))}",
                        "pct_both_amp": round(100.0 * ((a >= 4) & (b >= 4)).mean(), 2),
                        "pct_both_gain_or_amp": round(100.0 * ((a >= 3) & (b >= 3)).mean(), 2),
                        "pct_both_diploid": round(100.0 * ((a == 2) & (b == 2)).mean(), 2),
                        "pct_anchor_amp_neighbor_loss": round(
                            100.0 * ((a >= 4) & (b <= 1)).mean(), 2
                        ),
                        "spearman_cn": round(float(rho), 4) if np.isfinite(rho) else np.nan,
                        "spearman_p": rho_p,
                        "median_segment_mb": round(
                            _segment_length_mb(str(sub["anchor_segment"].mode().iloc[0]))
                            if "anchor_segment" in sub.columns and sub["anchor_segment"].notna().any()
                            else np.nan,
                            2,
                        ),
                    }
                )

    detail = pd.DataFrame(long_rows)
    if detail.empty:
        return detail, detail

    # Summary: rank neighbors by subtype-specific co-amp with anchor
    summary = (
        detail.sort_values(["anchor_arm", "neighbor_type", "pam50", "pct_both_amp"], ascending=[True, True, True, False])
        .groupby(["anchor_arm", "neighbor_type", "neighbor_locus_id", "pam50"], as_index=False)
        .first()
    )
    return detail, summary


def _segment_length_mb(seg_label: str) -> float:
    """Parse chr:start-end → length in Mb."""
    try:
        part = str(seg_label).split(":", 1)[1]
        s, e = part.split("-")
        return (int(e) - int(s) + 1) / 1e6
    except (ValueError, IndexError):
        return np.nan


def write_subtype_depth(
    long: pd.DataFrame,
    out_dir: Path,
) -> Dict[str, Path]:
    maps_dir = out_dir / "maps"
    maps_dir.mkdir(parents=True, exist_ok=True)

    participant_cn = build_focus_participant_cn(long)
    dist = build_subtype_distribution(participant_cn)
    comp = build_subtype_comparison(dist, participant_cn)
    pairwise = build_pairwise_subtype(participant_cn)
    nb_detail, nb_summary = build_neighborhood_by_subtype(long)

    paths = {
        "focus_participant_cn": maps_dir / "mirna_cnv_focus_participant_cn.tsv.gz",
        "focus_subtype_distribution": maps_dir / "mirna_cnv_focus_subtype_distribution.tsv",
        "focus_subtype_comparison": maps_dir / "mirna_cnv_focus_subtype_comparison.tsv",
        "focus_pairwise_subtype": maps_dir / "mirna_cnv_focus_pairwise_subtype.tsv",
        "neighborhood_by_subtype": maps_dir / "mirna_cnv_neighborhood_by_subtype.tsv",
        "neighborhood_summary": maps_dir / "mirna_cnv_neighborhood_summary.tsv",
    }
    participant_cn.to_csv(paths["focus_participant_cn"], sep="\t", index=False, compression="gzip")
    dist.to_csv(paths["focus_subtype_distribution"], sep="\t", index=False)
    comp.to_csv(paths["focus_subtype_comparison"], sep="\t", index=False)
    pairwise.to_csv(paths["focus_pairwise_subtype"], sep="\t", index=False)
    nb_detail.to_csv(paths["neighborhood_by_subtype"], sep="\t", index=False)
    nb_summary.to_csv(paths["neighborhood_summary"], sep="\t", index=False)

    print(
        f"[mirna_cnv_subtype_depth] focus arms: {participant_cn['arm_label'].nunique() if not participant_cn.empty else 0}; "
        f"distribution rows: {len(dist)}; neighborhood rows: {len(nb_detail)}"
    )
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--long-cnv", type=Path, default=C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz")
    ap.add_argument("--out-dir", type=Path, default=C.MIRNA_LOCUS_CNV_DIR)
    args = ap.parse_args()
    long = pd.read_csv(args.long_cnv, sep="\t", low_memory=False)
    write_subtype_depth(long, args.out_dir)


if __name__ == "__main__":
    main()
