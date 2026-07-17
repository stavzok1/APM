"""Orphan-edge discovery via CPTAC protein: miRNA→gene candidates absent from miRTarBase.

Follow-up to MH-34/MH-35. The curated `protein_resid` survivors (MH-35) were all *known*
edges. Here we ask the inverse: do **orphan** edges — miRNA→gene pairs with a sequence
(TargetScan) or CLIP (ENCORI) prior but **no miRTarBase functional-MTI curation** — show
protein-level anti-correlation across the CPTAC cohorts, nominating *new* candidate edges?

Why CPTAC adds credibility over the prior mRNA-only orphan work (MH-23/26):
  - **protein_resid** (protein below its own mRNA-predicted level) is a *translational*
    signature much harder to fake by indirect co-expression than mRNA anti-correlation.
  - **two-cohort** screen: prospective (CPTAC-2, primary) + TCGA-105 replication direction.

Orphan edge = (TargetScan |weighted context++| ≥ 0.25  OR  ENCORI clipExpNum ≥ 2) AND
**not** a miRTarBase high-evidence (functional-MTI) edge. `mirtar_any` flags pairs that have
*some* (non-functional-MTI) miRTarBase row vs fully uncurated.

Method (screen): partial Spearman of the orphan arm's expression vs the gene's `protein_z`
and `protein_resid`, **cov-residualized** via the SHARED `cptac_validation._covariates` block
(purity/CPE + CIN/HRD **plus the 8 Wu-major non-malignant lineages + malignant proliferation**),
across the cohort; BH-FDR over all tested edges. Candidates = prospective FDR-negative;
replication = TCGA-105 sign-concordant negative.

⚠ **MH-114: the cell-composition block is NOT optional here.** Until 2026-07-12 this lane residualized
on purity+CIN only, so it nominated orphans from a screen with **no composition block** while the model
it feeds has one. Bulk protein and bulk miRNA are both compartment-weighted averages, so a miRNA that is
merely *expressed in a different cell type* than the target protein produces an anti-correlation with no
cell-autonomous repression. Measured (MH-114, HE genes): ~60% of the apparent protein anti-correlation is
compartment arithmetic; **CAF-marker proteins collapse to ~0** (CALD1 −0.494 → +0.032, MMP2 −0.338 →
−0.108) while tumour-intrinsic ones retain about half (PTEN −0.290 → −0.151). The driver is the stromal
**MIX**, not tumour content (a block with tumour content removed reproduces it: ZEB1 −0.199 vs −0.172).

**These are NOMINATIONS for wet-lab (CLIP/luciferase), not validated edges** — protein
anti-correlation cannot prove a direct edge (co-expression/indirect paths, seed-family
ambiguity remain). The protein-layer analogue of `mirtar_gap_coupling_queue` (MH-26).

Outputs (`output/cptac_validation/orphan_discovery/`):
  - `orphan_edges_tested.tsv.gz`   every tested edge + both-cohort stats + priors
  - `orphan_candidates.tsv`        FDR-negative prospective candidates, ranked, with tags
  - `method_manifest.json`

Run::

    .venv/bin/python3 -m mirna_hallmark.analyses.cptac.cptac_orphan_discovery
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.encori_edges import load_collapsed_encori_pairs
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.eval.targetscan_orphan_coupling import (
    MIN_TS_WEIGHT,
    _he_pairs,
    _mirtar_lookup,
    build_orphan_edge_table,
)
from mirna_hallmark.eval.cptac_validation import (
    OUT_DIR,
    _covariates,                # SHARED composition-aware confounder block — do not rebuild it locally
    load_cptac_layers,
    load_prospective_clinical,
    load_prospective_mirna_arms,
)

DISCOVERY_DIR = OUT_DIR / "orphan_discovery"
MIN_CLIP = 2
MIN_N = 20
FDR = 0.1


# --------------------------------------------------------------------------- #
# Orphan edge set
# --------------------------------------------------------------------------- #
def build_orphan_edges(universe: Sequence[str]) -> pd.DataFrame:
    he = D.high_evidence_edges(D.load_hallmark_edges())
    he_pairs = _he_pairs(he)
    mirtar = _mirtar_lookup(C.MIRTAR_HALLMARK_SUMMARY)

    ts = build_orphan_edge_table(list(universe), he_pairs, min_ts=MIN_TS_WEIGHT, mirtar=mirtar)
    ts = ts[["miRNA", "gene", "ts_weight", "n_studies", "mirtar_any"]].copy()
    ts["in_targetscan"] = True

    enc = load_collapsed_encori_pairs(genes=list(universe))
    enc = enc[pd.to_numeric(enc["clipExpNum"], errors="coerce").fillna(0) >= MIN_CLIP].copy()
    enc["pair"] = list(zip(enc["miRNA"], enc["gene"]))
    enc = enc[~enc["pair"].isin(he_pairs)]
    enc = enc[["miRNA", "gene", "clipExpNum", "enc_depth"]].copy()
    enc["in_encori"] = True

    merged = ts.merge(enc, on=["miRNA", "gene"], how="outer")
    merged["in_targetscan"] = merged["in_targetscan"].fillna(False)
    merged["in_encori"] = merged["in_encori"].fillna(False)
    merged["mirtar_any"] = merged["mirtar_any"].fillna(False)
    merged["n_studies"] = merged["n_studies"].fillna(0).astype(int)
    merged["prior"] = np.where(
        merged["in_targetscan"] & merged["in_encori"], "targetscan+encori",
        np.where(merged["in_targetscan"], "targetscan_only", "encori_only"),
    )
    return merged


# --------------------------------------------------------------------------- #
# Fast cov-residualized partial Spearman screen
# --------------------------------------------------------------------------- #
def _residualize(M: pd.DataFrame, cov: pd.DataFrame) -> pd.DataFrame:
    """Row-wise OLS residual of M (feature×sample) on cov (sample×k), shared complete samples."""
    samples = [c for c in M.columns if c in cov.index]
    cv = cov.loc[samples].apply(pd.to_numeric, errors="coerce")
    keep = cv.dropna().index
    X = np.column_stack([np.ones(len(keep)), cv.loc[keep].values])
    H = X @ np.linalg.pinv(X.T @ X) @ X.T          # hat matrix (samples×samples)
    sub = M[list(keep)].astype(float)
    R = sub.values - sub.values @ H.T               # residualize each row
    out = pd.DataFrame(R, index=M.index, columns=list(keep))
    return out


def _rank_matrix(R: pd.DataFrame) -> pd.DataFrame:
    """Rank each row across samples (NaN-aware), for Spearman as Pearson-of-ranks."""
    vals = R.values.astype(float)
    out = np.full_like(vals, np.nan)
    for i in range(vals.shape[0]):
        row = vals[i]
        m = ~np.isnan(row)
        if m.sum() >= MIN_N:
            out[i, m] = rankdata(row[m])
    return pd.DataFrame(out, index=R.index, columns=R.columns)


def _edge_partial_spearman(
    arm_rank: pd.Series, gene_rank: pd.Series
) -> Tuple[float, int]:
    df = pd.concat([arm_rank.rename("a"), gene_rank.rename("g")], axis=1).dropna()
    if len(df) < MIN_N:
        return (np.nan, len(df))
    a = df["a"].values - df["a"].mean()
    g = df["g"].values - df["g"].mean()
    denom = np.sqrt((a * a).sum() * (g * g).sum())
    return (float((a * g).sum() / denom) if denom else np.nan, len(df))


def screen_cohort(
    edges: pd.DataFrame, arms: pd.DataFrame, layers: Dict[str, pd.DataFrame], cov: pd.DataFrame
) -> pd.DataFrame:
    """Return per-edge partial-Spearman rho (protein_z & protein_resid) for one cohort."""
    arm_rank = _rank_matrix(_residualize(arms, cov))
    prot_rank = _rank_matrix(_residualize(layers["protein_z"], cov))
    resid_rank = _rank_matrix(_residualize(layers["protein_resid"], cov))
    arm_set, prot_set = set(arm_rank.index), set(prot_rank.index)

    sub = edges[edges["miRNA"].isin(arm_set) & edges["gene"].isin(prot_set)]
    rows = []
    for mir, gene in zip(sub["miRNA"], sub["gene"]):
        pr, n = _edge_partial_spearman(arm_rank.loc[mir], prot_rank.loc[gene])
        rr, _ = _edge_partial_spearman(arm_rank.loc[mir], resid_rank.loc[gene])
        rows.append({"miRNA": mir, "gene": gene, "n": n, "protein_rho": pr, "resid_rho": rr})
    return pd.DataFrame(rows)


def _spearman_p(rho: float, n: int) -> float:
    if not np.isfinite(rho) or n < MIN_N or abs(rho) >= 1:
        return np.nan
    from scipy.stats import t as tdist
    tstat = rho * np.sqrt((n - 2) / (1 - rho * rho))
    return float(2 * tdist.sf(abs(tstat), n - 2))


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = DISCOVERY_DIR, batch_kind: str = "none") -> None:
    if batch_kind != "none":
        out_dir = out_dir / f"batch_{batch_kind}"
    out_dir.mkdir(parents=True, exist_ok=True)
    from mirna_hallmark import cptac_batch as B
    universe = list(HallmarkSets.load().universe)
    edges = build_orphan_edges(universe)
    print(f"[orphan] orphan edge set: {len(edges)} "
          f"(TS-only {int((edges.prior=='targetscan_only').sum())}, "
          f"ENCORI-only {int((edges.prior=='encori_only').sum())}, "
          f"both {int((edges.prior=='targetscan+encori').sum())})")

    # Prospective (primary).
    pros_arms = load_prospective_mirna_arms()
    pros_layers = load_cptac_layers("prospective")
    pros_clin = load_prospective_clinical().drop_duplicates("participant").set_index("participant")
    # MH-114: route through the SHARED composition-aware builder. This lane used to take `purity+cin`
    # straight off the clinical frame, so it screened orphans with NO cell-composition block while the
    # model it feeds has one — the MH-107 defect, in the lane that produces the orphan nominations.
    # Measured on the HE genes: ~60% of the apparent protein anti-correlation is compartment arithmetic
    # (stromal MIX, not tumour content), and CAF-marker proteins (CALD1, MMP2) go to ~0.
    pros_cov, pbc = B.augment_cov(_covariates(load_prospective_clinical(), "prospective"),
                                  "prospective", batch_kind)
    if pbc:
        print(f"[orphan] +batch ({batch_kind}) prospective: {len(pbc)} dummies")
    print(f"[orphan] prospective covariates ({pros_cov.shape[1]}): {list(pros_cov.columns)}")
    print("[orphan] screening prospective ...")
    pros = screen_cohort(edges, pros_arms, pros_layers, pros_cov)
    pros["protein_p"] = [_spearman_p(r, n) for r, n in zip(pros["protein_rho"], pros["n"])]
    pros["resid_p"] = [_spearman_p(r, n) for r, n in zip(pros["resid_rho"], pros["n"])]
    pros["protein_q"] = S.bh_fdr(pros["protein_p"].fillna(1.0).values)
    pros["resid_q"] = S.bh_fdr(pros["resid_p"].fillna(1.0).values)

    # TCGA-105 (replication direction).
    tcga_arms = D.load_mirna_arms()
    tcga_layers = load_cptac_layers("tcga105")
    tcga_clin = D.load_clinical_strata().drop_duplicates("participant").set_index("participant")
    tcga_cov, tbc = B.augment_cov(_covariates(D.load_clinical_strata(), "tcga105"), "tcga105", batch_kind)
    if tbc:
        print(f"[orphan] +batch ({batch_kind}) tcga105: {len(tbc)} dummies")
    print(f"[orphan] tcga105 covariates ({tcga_cov.shape[1]}): {list(tcga_cov.columns)}")
    print("[orphan] screening tcga105 (replication) ...")
    tcga = screen_cohort(edges, tcga_arms, tcga_layers, tcga_cov)[["miRNA", "gene", "protein_rho", "resid_rho"]]
    tcga = tcga.rename(columns={"protein_rho": "tcga_protein_rho", "resid_rho": "tcga_resid_rho"})

    out = (edges.merge(pros, on=["miRNA", "gene"], how="inner")
                .merge(tcga, on=["miRNA", "gene"], how="left"))

    # Tags.
    prot_sig = (out["protein_q"] < FDR) & (out["protein_rho"] < 0)
    resid_sig = (out["resid_q"] < FDR) & (out["resid_rho"] < 0)
    tcga_conc = out["tcga_protein_rho"] < 0
    out["protein_candidate"] = prot_sig
    out["translational_candidate"] = resid_sig
    out["tcga_replicated"] = tcga_conc
    out["strength_tag"] = np.select(
        [resid_sig & tcga_conc, prot_sig & tcga_conc, resid_sig, prot_sig],
        ["S_translational_replicated", "S_protein_replicated", "P_translational_prospective", "P_protein_prospective"],
        default="ns",
    )

    out = out.sort_values(["protein_candidate", "translational_candidate", "protein_rho"],
                          ascending=[False, False, True])
    out.to_csv(out_dir / "orphan_edges_tested.tsv.gz", sep="\t", index=False, compression="gzip")
    cand = out[prot_sig | resid_sig].copy()
    cand.to_csv(out_dir / "orphan_candidates.tsv", sep="\t", index=False)

    tags = out["strength_tag"].value_counts().to_dict()
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.analyses.cptac.cptac_orphan_discovery",
        "orphan_definition": f"(TargetScan ts_weight>={MIN_TS_WEIGHT} OR ENCORI clipExpNum>={MIN_CLIP}) AND not miRTarBase functional-MTI",
        "n_orphan_edges": int(len(edges)),
        "n_tested_prospective": int(len(out)),
        "n_protein_candidates": int(prot_sig.sum()),
        "n_translational_candidates": int(resid_sig.sum()),
        "n_translational_replicated": int((resid_sig & tcga_conc).sum()),
        "strength_tag_counts": tags,
        "batch_kind": batch_kind,
        "method": "cov-residualized partial Spearman screen (prospective primary, tcga105 replication direction); BH-FDR"
                  + (f"; +batch({batch_kind}) MS-plex/site dummies" if batch_kind != "none" else ""),
        "caveat": "NOMINATIONS for CLIP/luciferase, not validated edges; protein anti-corr != direct edge",
        "outputs": ["orphan_edges_tested.tsv.gz", "orphan_candidates.tsv"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"[orphan] wrote {out_dir}")
    print(f"  protein candidates: {int(prot_sig.sum())} | translational: {int(resid_sig.sum())} "
          f"| translational+replicated: {int((resid_sig & tcga_conc).sum())}")
    show = out[resid_sig | prot_sig].head(15)
    if not show.empty:
        print(show[["miRNA", "gene", "prior", "protein_rho", "protein_q", "resid_rho", "resid_q",
                    "tcga_protein_rho", "strength_tag"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=DISCOVERY_DIR)
    ap.add_argument("--batch-kind", default="none", choices=["none", "site", "plex", "auto"],
                    help="Add MS-plex/site batch dummies to the cov (writes to batch_<kind>/ subdir)")
    ap.add_argument("--min-purity", type=float, default=None, help="Accepted for runner compatibility (no-op)")
    args = ap.parse_args()
    run(out_dir=args.out_dir, batch_kind=args.batch_kind)


if __name__ == "__main__":
    main()
