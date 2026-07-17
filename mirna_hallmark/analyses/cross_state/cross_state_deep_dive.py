"""Deep dive across the three tissue states (tumor / NAT / GTEx-healthy).

Two mirror reports that assemble *everything* we have around one anchor:

GENE deep dive (``--gene PTEN``):
  For one target gene, every miRTarBase regulator arm with, in each of the three
  states, its abundance (median + within-state percentile rank), its realized
  pressure share on the gene (no-z ``softmax_logrpm`` attribution), and its raw
  Spearman coupling to the gene. Then a TargetScan-seed-family roll-up (family
  share + abundance) and a within-family resolution (which arm of a multi-arm
  family actually carries the family's effect on this gene).

miRNA deep dive (``--mirna hsa-miR-200c-3p``):
  The mirror, anchored on one miRNA arm: every Hallmark target it regulates with
  the arm's realized share + coupling per state, the arm's abundance trajectory
  across states, its TargetScan seed-family roster, and a within-family target
  resolution (for targets shared by the family, which arm couples).

Coupling is raw Spearman in every state (the only thing computable in NAT/GTEx),
so it is comparable across states; the tumor *partial* (CN/CPE/HRD-adjusted) rho
lives in the dedicated coupling modules, not here.

Cross-platform: GTEx miRNA (TPM/URS) vs TCGA (RPM/MIMAT) -> read GTEx as
rank/share/sign, not absolute magnitude. mRNA is TPM in both.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_deep_dive --gene PTEN
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_deep_dive --mirna hsa-miR-200c-3p
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _gtex_all_arm_matrix
from mirna_hallmark.family_normal_reference import (
    _gtex_available,
    _gtex_breast_rna,
    _load_full_mirna,
    _load_full_rna,
    _participant,
    _split_types,
)
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import compute_gene_pressure_contributions, load_mirtar_edges

OUT_DIR = C.TISSUE_REFERENCE_DIR / "cross_state_deep_dive"
ATTR_MODE = C.PRESSURE_ATTRIBUTION_EXPR_MODE
MIN_COUPLING_N = 25
_ARM_SUFFIX = re.compile(r"\.\d+$")


# --------------------------------------------------------------------------- #
# Family map (TargetScan seed families, human)
# --------------------------------------------------------------------------- #
def _family_map() -> Dict[str, str]:
    df = pd.read_csv(C.MIR_FAMILY_INFO, sep="\t")
    h = df.loc[df["Species ID"].eq(9606), ["miR family", "MiRBase ID"]].dropna()
    return dict(zip(h["MiRBase ID"].astype(str), h["miR family"].astype(str)))


def _fam_of(arm: str, fam_map: Dict[str, str]) -> str:
    base = _ARM_SUFFIX.sub("", str(arm))
    return fam_map.get(base, fam_map.get(arm, f"<{base}>"))


# --------------------------------------------------------------------------- #
# State bundles: per state, miRNA arm matrix + gene RNA matrix on shared IDs
# --------------------------------------------------------------------------- #
def _state_bundles(genes: Sequence[str]) -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame]]:
    full_mir = _load_full_mirna()
    mir_split = _split_types(full_mir)
    rna_split = _split_types(_load_full_rna(genes))

    def _rename(df: pd.DataFrame) -> pd.DataFrame:
        df = df.rename(columns={c: _participant(c) for c in df.columns})
        return df.loc[:, ~df.columns.duplicated()]

    bundles: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]] = {
        "tumor": (_rename(mir_split["tumor"]), _rename(rna_split["tumor"])),
        "nat": (_rename(mir_split["nat"]), _rename(rna_split["nat"])),
    }
    if _gtex_available():
        gmir = _gtex_all_arm_matrix()
        grna = _gtex_breast_rna(genes)
        if not gmir.empty and not grna.empty:
            bundles["gtex"] = (gmir, grna)
    return bundles


def _arm_pct(mir: pd.DataFrame) -> pd.Series:
    """Within-state abundance percentile rank for every arm (1=highest)."""
    return mir.median(axis=1).rank(pct=True)


def _raw_rho(mir: pd.DataFrame, rna: pd.DataFrame, arm: str, gene: str) -> Tuple[float, float, int]:
    if arm not in mir.index or gene not in rna.index:
        return np.nan, np.nan, 0
    shared = sorted(set(mir.columns) & set(rna.columns))
    x = pd.to_numeric(mir.loc[arm, shared], errors="coerce")
    y = pd.to_numeric(rna.loc[gene, shared], errors="coerce")
    ok = x.notna() & y.notna()
    n = int(ok.sum())
    if n < MIN_COUPLING_N or x[ok].nunique() < 2 or y[ok].nunique() < 2:
        return np.nan, np.nan, n
    rho, p = spearmanr(x[ok], y[ok])
    return float(rho), float(p), n


def _shares_for_genes(
    genes: Sequence[str], bundles: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]], edges: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """Per state -> per-(gene, miRNA) realized share table (global_abs_share within gene)."""
    out: Dict[str, pd.DataFrame] = {}
    for state, (mir, _rna) in bundles.items():
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=mir, expr_mode=ATTR_MODE,  # type: ignore[arg-type]
            resolve_arms=False,
        )
        out[state] = contrib
    return out


def _acquired_tumor_shares(genes: Sequence[str], tumor_mir: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    """Per-(gene, miRNA) HEALTHY-ANCHORED (acquired-vs-healthy) share in the tumor state."""
    from mirna_hallmark.healthy_anchor import load_baseline
    return compute_gene_pressure_contributions(
        genes, edges=edges, mirna=tumor_mir, expr_mode=C.PRESSURE_HEALTHY_ANCHOR_MODE,  # type: ignore[arg-type]
        resolve_arms=False, healthy_baseline=load_baseline(),
    )


def _lookup_share(tbl: pd.DataFrame, gene: str, arm: str) -> float:
    if tbl is None or tbl.empty:
        return np.nan
    r = tbl.loc[tbl["gene"].eq(gene) & tbl["miRNA"].eq(arm)]
    return float(r["global_abs_share"].iloc[0]) if len(r) else np.nan


# --------------------------------------------------------------------------- #
# GENE deep dive
# --------------------------------------------------------------------------- #
def gene_deep_dive(gene: str, *, out_dir: Path = OUT_DIR, healthy_anchor: bool = False) -> Dict[str, pd.DataFrame]:
    fam_map = _family_map()
    edges = load_mirtar_edges([gene], resolve_arms=True)
    if edges.empty:
        raise SystemExit(f"[deep_dive] no miRTarBase regulators for {gene}")
    bundles = _state_bundles([gene])
    pcts = {s: _arm_pct(mir) for s, (mir, _r) in bundles.items()}
    shares = _shares_for_genes([gene], bundles, edges)
    # gated until a reliable true-healthy reference (DIANA-miTED) is ingested
    acquired = _acquired_tumor_shares([gene], bundles["tumor"][0], edges) if healthy_anchor else pd.DataFrame()

    regs = sorted(edges["miRNA"].unique())
    rows: List[dict] = []
    for arm in regs:
        ev = float(edges.loc[edges["miRNA"].eq(arm), "evidence_score"].max())
        row = {"gene": gene, "miRNA": arm, "family": _fam_of(arm, fam_map), "evidence_score": ev}
        for state, (mir, rna) in bundles.items():
            med = float(mir.loc[arm].median()) if arm in mir.index else np.nan
            row[f"abund_{state}"] = med
            row[f"pct_{state}"] = float(pcts[state].get(arm, np.nan))
            sh = shares[state]
            srow = sh.loc[sh["miRNA"].eq(arm)] if not sh.empty else sh
            row[f"share_{state}"] = float(srow["global_abs_share"].iloc[0]) if len(srow) else np.nan
            rho, p, n = _raw_rho(mir, rna, arm, gene)
            row[f"rho_{state}"] = rho
            row[f"rho_p_{state}"] = p
            row[f"n_{state}"] = n
        # healthy-anchored: how much of the ACQUIRED (vs healthy) pressure this arm carries
        row["acquired_share_tumor"] = _lookup_share(acquired, gene, arm)
        rows.append(row)
    reg = pd.DataFrame(rows)
    sort_state = "share_tumor" if "share_tumor" in reg.columns else "evidence_score"
    reg = reg.sort_values(sort_state, ascending=False).reset_index(drop=True)

    # family roll-up
    fam_rows: List[dict] = []
    for fam, grp in reg.groupby("family"):
        frow = {"gene": gene, "family": fam, "n_arms": len(grp),
                "arms": ";".join(grp.sort_values(sort_state, ascending=False)["miRNA"])}
        for state in bundles:
            frow[f"family_share_{state}"] = float(grp[f"share_{state}"].sum(skipna=True))
            # arm carrying the family's pressure on this gene + its coupling
            sub = grp.dropna(subset=[f"share_{state}"])
            if len(sub):
                top = sub.loc[sub[f"share_{state}"].idxmax()]
                frow[f"carrier_{state}"] = top["miRNA"]
                frow[f"carrier_rho_{state}"] = top[f"rho_{state}"]
        fam_rows.append(frow)
    fam = pd.DataFrame(fam_rows)
    if not fam.empty and "family_share_tumor" in fam.columns:
        fam = fam.sort_values("family_share_tumor", ascending=False).reset_index(drop=True)

    gdir = out_dir / gene
    gdir.mkdir(parents=True, exist_ok=True)
    reg.to_csv(gdir / "regulator_profile.tsv", sep="\t", index=False)
    fam.to_csv(gdir / "family_profile.tsv", sep="\t", index=False)
    _write_manifest(gdir, {"mode": "gene", "anchor": gene, "n_regulators": len(reg),
                           "states": list(bundles), "n_families": len(fam)})

    print(f"\n=== GENE deep dive: {gene} ({len(reg)} regulators, states={list(bundles)}) ===")
    acq_ok = healthy_anchor and reg["acquired_share_tumor"].notna().any()
    cols0 = ["miRNA", "family", "share_tumor"] + (["acquired_share_tumor"] if acq_ok else [])
    show = [c for c in cols0 + ["rho_tumor", "rho_nat", "rho_gtex",
                                "pct_tumor", "pct_nat", "pct_gtex"] if c in reg.columns]
    print("Top regulators by tumor pressure share"
          + (" (share_tumor=within-tumor, acquired=vs-healthy):" if acq_ok else ":"))
    print(reg.head(12)[show].to_string(index=False))
    if acq_ok:
        acq = reg.dropna(subset=["acquired_share_tumor"]).sort_values("acquired_share_tumor", ascending=False)
        print("\nTop regulators by ACQUIRED (vs-healthy) pressure share:")
        print(acq.head(10)[[c for c in ["miRNA", "family", "acquired_share_tumor", "share_tumor", "rho_tumor"] if c in acq.columns]].to_string(index=False))
    if not fam.empty:
        fshow = [c for c in ["family", "n_arms", "family_share_tumor", "carrier_tumor",
                             "carrier_rho_tumor", "carrier_rho_nat", "carrier_rho_gtex"] if c in fam.columns]
        print("\nTop families by tumor pressure share (within-family carrier + coupling):")
        print(fam.head(8)[fshow].to_string(index=False))
    return {"regulators": reg, "families": fam}


# --------------------------------------------------------------------------- #
# miRNA deep dive (mirror)
# --------------------------------------------------------------------------- #
def mirna_deep_dive(
    mirna: str, *, family: Optional[str] = None, out_dir: Path = OUT_DIR, max_targets: int = 120,
    healthy_anchor: bool = False,
) -> Dict[str, pd.DataFrame]:
    fam_map = _family_map()
    fam = family or _fam_of(mirna, fam_map)
    fam_arms = sorted({a for a, f in fam_map.items() if f == fam})
    hs = HallmarkSets.load()
    universe = sorted(hs.universe)

    edges = load_mirtar_edges(universe, resolve_arms=True)
    fam_edges = edges.loc[edges["miRNA"].isin(set(fam_arms) | {mirna})]
    if fam_edges.empty:
        raise SystemExit(f"[deep_dive] no Hallmark targets for {mirna} / family {fam}")
    m_targets = sorted(fam_edges.loc[fam_edges["miRNA"].eq(mirna), "gene"].unique())
    if not m_targets:
        raise SystemExit(f"[deep_dive] {mirna} has no Hallmark targets")
    # bound by evidence for tractable coupling
    ev_rank = (fam_edges.loc[fam_edges["miRNA"].eq(mirna)]
               .groupby("gene")["evidence_score"].max().sort_values(ascending=False))
    m_targets = list(ev_rank.head(max_targets).index)

    target_universe = sorted(set(fam_edges["gene"]) & set(m_targets) | set(m_targets))
    bundles = _state_bundles(target_universe)
    pcts = {s: _arm_pct(mir) for s, (mir, _r) in bundles.items()}
    shares = _shares_for_genes(m_targets, bundles, edges)
    acquired = _acquired_tumor_shares(m_targets, bundles["tumor"][0], edges) if healthy_anchor else pd.DataFrame()

    # ---- arm abundance trajectory (the anchor + its family) ----
    arm_rows: List[dict] = []
    for arm in sorted(set(fam_arms) | {mirna}):
        r = {"miRNA": arm, "family": fam, "is_anchor": arm == mirna,
             "n_hallmark_targets": int((fam_edges["miRNA"].eq(arm)).sum())}
        for state, (mir, _rna) in bundles.items():
            r[f"abund_{state}"] = float(mir.loc[arm].median()) if arm in mir.index else np.nan
            r[f"pct_{state}"] = float(pcts[state].get(arm, np.nan))
        arm_rows.append(r)
    arm_summary = pd.DataFrame(arm_rows).sort_values("pct_tumor", ascending=False, na_position="last")

    # ---- anchor target profile (gene-level: share + coupling per state) ----
    trows: List[dict] = []
    for gene in m_targets:
        ev = float(fam_edges.loc[fam_edges["miRNA"].eq(mirna) & fam_edges["gene"].eq(gene), "evidence_score"].max())
        row = {"miRNA": mirna, "gene": gene, "evidence_score": ev}
        for state, (mir, rna) in bundles.items():
            sh = shares[state]
            srow = sh.loc[sh["miRNA"].eq(mirna) & sh["gene"].eq(gene)] if not sh.empty else sh
            row[f"share_{state}"] = float(srow["global_abs_share"].iloc[0]) if len(srow) else np.nan
            rho, p, n = _raw_rho(mir, rna, mirna, gene)
            row[f"rho_{state}"] = rho
            row[f"rho_p_{state}"] = p
            row[f"n_{state}"] = n
        row["acquired_share_tumor"] = _lookup_share(acquired, gene, mirna)
        trows.append(row)
    target_profile = pd.DataFrame(trows).sort_values("share_tumor", ascending=False, na_position="last").reset_index(drop=True)

    # ---- within-family target resolution (targets shared by >1 family arm) ----
    res_rows: List[dict] = []
    shared_targets = (fam_edges.groupby("gene")["miRNA"].nunique())
    shared_targets = shared_targets[shared_targets > 1].index
    for gene in shared_targets:
        for arm in sorted(set(fam_edges.loc[fam_edges["gene"].eq(gene), "miRNA"])):
            row = {"gene": gene, "miRNA": arm, "is_anchor": arm == mirna}
            for state, (mir, rna) in bundles.items():
                sh = shares[state]
                srow = sh.loc[sh["miRNA"].eq(arm) & sh["gene"].eq(gene)] if not sh.empty else sh
                row[f"share_{state}"] = float(srow["global_abs_share"].iloc[0]) if len(srow) else np.nan
                rho, _p, _n = _raw_rho(mir, rna, arm, gene)
                row[f"rho_{state}"] = rho
            res_rows.append(row)
    resolution = pd.DataFrame(res_rows)

    mdir = out_dir / f"mirna_{mirna.replace('/', '_')}"
    mdir.mkdir(parents=True, exist_ok=True)
    arm_summary.to_csv(mdir / "family_arm_summary.tsv", sep="\t", index=False)
    target_profile.to_csv(mdir / "anchor_target_profile.tsv", sep="\t", index=False)
    resolution.to_csv(mdir / "family_target_resolution.tsv", sep="\t", index=False)
    _write_manifest(mdir, {"mode": "mirna", "anchor": mirna, "family": fam,
                           "family_arms": fam_arms, "n_targets": len(m_targets),
                           "states": list(bundles)})

    print(f"\n=== miRNA deep dive: {mirna}  (family '{fam}': {', '.join(fam_arms)}) ===")
    ashow = [c for c in ["miRNA", "is_anchor", "n_hallmark_targets", "pct_tumor", "pct_nat", "pct_gtex"] if c in arm_summary.columns]
    print("Family arm abundance trajectory (within-state percentile):")
    print(arm_summary[ashow].to_string(index=False))
    tshow = [c for c in ["gene", "share_tumor", "rho_tumor", "rho_nat", "rho_gtex"] if c in target_profile.columns]
    print(f"\nTop {mirna} targets by tumor pressure share:")
    print(target_profile.head(12)[tshow].to_string(index=False))
    return {"arm_summary": arm_summary, "target_profile": target_profile, "resolution": resolution}


def _write_manifest(d: Path, extra: dict) -> None:
    man = {"module": "mirna_hallmark.analyses.cross_state.cross_state_deep_dive",
           "generated_utc": datetime.now(timezone.utc).isoformat(),
           "attribution_expr_mode": ATTR_MODE, "min_coupling_n": MIN_COUPLING_N, **extra}
    (d / "method_manifest.json").write_text(json.dumps(man, indent=2), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", type=str, default=None, help="anchor gene for gene deep dive")
    ap.add_argument("--mirna", type=str, default=None, help="anchor miRNA arm for mirror deep dive")
    ap.add_argument("--family", type=str, default=None, help="override TargetScan family for --mirna")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--healthy-anchor", action="store_true",
                    help="also emit acquired-vs-healthy pressure shares (needs miTED reference)")
    args = ap.parse_args()
    C.ensure_output_dirs()
    if not args.gene and not args.mirna:
        args.gene = "PTEN"
    if args.gene:
        gene_deep_dive(args.gene, out_dir=args.out_dir, healthy_anchor=args.healthy_anchor)
    if args.mirna:
        mirna_deep_dive(args.mirna, family=args.family, out_dir=args.out_dir, healthy_anchor=args.healthy_anchor)


if __name__ == "__main__":
    main()
