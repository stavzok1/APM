"""Resolve the CPTAC protein-residual survivors: which miRNA, which subtype, vs literature.

Follow-up to `cptac_validation` (MH-34). The `protein_resid` layer isolates genes whose
**protein sits below its own mRNA-predicted level as miRNA pressure rises** — the genuine
*translational / co-translational* repression component (the part the cohort's RNA-seq does
NOT already capture). Those FDR survivors exist only in the **independent prospective cohort**
(44 genes; 0 in same-patient TCGA-105). This module asks, for each survivor:

1. **Which miRNA arm drives it?** Single-edge test (the Q1 refinement): for each
   high-evidence regulator of the gene, partial-Spearman the arm's *own expression* vs the
   gene's `protein_resid` across the prospective samples (adjust purity + CIN). The most
   negative, FDR-surviving arm is the candidate translational repressor.
2. **Is it subtype-specific?** Per-PAM50 raw Spearman of the driver arm vs `protein_resid`
   (descriptive; small n) — which lineage carries the coupling.
3. **What does the literature say?** Annotate each driver (arm→gene) edge with miRTarBase
   functional-MTI study count and the MH-31 breast-context PMID layer → validated / breast-
   specific / HE-only-prediction.

Cohort: **prospective** only (the survivors are prospective). Pressure/expression loaders are
reused from `cptac_validation`.

Outputs (`output/cptac_validation/prospective/resid_survivors/`):
  - `survivor_arm_candidates.tsv`  every (gene, candidate arm) coupling + literature columns
  - `survivor_drivers.tsv`         one row / survivor: top driver arm + subtype + lit flag
  - `survivor_subtype.tsv`         (gene, driver arm, PAM50) coupling
  - `method_manifest.json`

Run::

    .venv/bin/python3 -m mirna_hallmark.cptac_resid_survivors
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.eval.cptac_validation import (
    OUT_DIR,
    load_cptac_layers,
    load_prospective_clinical,
    load_prospective_mirna_arms,
)

SURVIVOR_DIR = OUT_DIR / "prospective" / "resid_survivors"
GENE_TBL = OUT_DIR / "prospective" / "gene_level_associations.tsv.gz"
BREAST_CTX = C.OUTPUT_ROOT / "tissue_reference/edge_breast_context/edge_breast_context_scores.tsv"
MIN_N = 15
MIN_SUBTYPE_N = 15


def survivor_genes() -> pd.DataFrame:
    """Survivor genes from the saved prospective gene-level protein_resid (q<0.1, rho<0)."""
    g = pd.read_csv(GENE_TBL, sep="\t")
    sub = g[(g["layer"] == "protein_resid") & (g["partial_q"] < 0.1) & (g["partial_rho"] < 0)]
    # Best (most negative) variant per gene + which variants flagged it.
    agg = (
        sub.sort_values("partial_rho")
        .groupby("gene_name")
        .agg(best_variant=("variant", "first"),
             best_resid_rho=("partial_rho", "first"),
             best_resid_q=("partial_q", "first"),
             n_variants=("variant", "nunique"),
             highev=("variant", lambda v: any(x.startswith("highev") for x in v)))
        .reset_index()
    )
    return agg.sort_values("best_resid_rho")


def _candidate_arms(gene: str, edges: pd.DataFrame, arm_index: set) -> List[str]:
    he = edges[(edges["gene"] == gene)]
    arms = sorted(set(he["miRNA"].astype(str)) & arm_index)
    return arms


def arm_couplings(
    survivors: pd.DataFrame,
    arms: pd.DataFrame,
    resid: pd.DataFrame,
    cov: pd.DataFrame,
    edges: pd.DataFrame,
    lit: pd.DataFrame,
) -> pd.DataFrame:
    """Per (gene, candidate arm): partial Spearman(arm expression, gene protein_resid)."""
    arm_index = set(arms.index.astype(str))
    he_lookup = edges.set_index(["miRNA", "gene"])["n_functional_mti_studies"].to_dict()
    lit_b = lit.set_index(["miRNA", "gene"]) if not lit.empty else None

    rows: List[dict] = []
    for _, srow in survivors.iterrows():
        gene = srow["gene_name"]
        if gene not in resid.index:
            continue
        y = resid.loc[gene]
        for arm in _candidate_arms(gene, edges, arm_index):
            x = arms.loc[arm]
            st = S.correlation_pair(y, x, cov)
            if st["n"] < MIN_N:
                continue
            n_fmti = he_lookup.get((arm, gene), np.nan)
            blit = lit_b.loc[(arm, gene)] if (lit_b is not None and (arm, gene) in lit_b.index) else None
            rows.append({
                "gene_name": gene,
                "miRNA": arm,
                "n": st["n"],
                "spearman_rho": _r(st["spearman_rho"]),
                "partial_rho": _r(st["partial_rho"]),
                "partial_p": st["partial_p"],
                "n_functional_mti_studies": n_fmti,
                "n_breast_pmids": int(blit["n_breast_pmids"]) if blit is not None else 0,
                "any_breast": bool(blit["any_breast"]) if blit is not None else False,
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["partial_q"] = S.bh_fdr(out["partial_p"].fillna(1.0).values)
    out["literature_status"] = np.where(
        out["any_breast"], "breast_validated",
        np.where(out["n_functional_mti_studies"].fillna(0) >= 1, "functional_mti_validated", "HE_prediction"),
    )
    return out.sort_values(["gene_name", "partial_rho"])


def subtype_coupling(
    drivers: pd.DataFrame, arms: pd.DataFrame, resid: pd.DataFrame, clinical: pd.DataFrame
) -> pd.DataFrame:
    pam = clinical.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"]
    rows: List[dict] = []
    for _, d in drivers.iterrows():
        gene, arm = d["gene_name"], d["driver_miRNA"]
        if gene not in resid.index or arm not in arms.index:
            continue
        y, x = resid.loc[gene], arms.loc[arm]
        shared = [c for c in y.index if c in x.index and c in pam.index]
        for sub in C.PAM50_ORDER:
            ss = [c for c in shared if pam.get(c) == sub]
            df = pd.concat([y[ss].rename("y"), x[ss].rename("x")], axis=1).dropna()
            if len(df) < MIN_SUBTYPE_N:
                continue
            from scipy.stats import spearmanr
            rho, p = spearmanr(df["x"], df["y"])
            rows.append({"gene_name": gene, "driver_miRNA": arm, "pam50": sub,
                         "n": len(df), "spearman_rho": _r(rho), "spearman_p": float(p)})
    return pd.DataFrame(rows)


def _r(v) -> float:
    return round(float(v), 4) if v is not None and np.isfinite(v) else np.nan


def run(*, out_dir: Path = SURVIVOR_DIR, batch_kind: str = "none") -> None:
    if batch_kind != "none":
        out_dir = out_dir / f"batch_{batch_kind}"
    out_dir.mkdir(parents=True, exist_ok=True)
    survivors = survivor_genes()
    print(f"[resid_survivors] {len(survivors)} protein_resid survivor genes (prospective)")

    arms = load_prospective_mirna_arms()
    layers = load_cptac_layers("prospective")
    resid = layers["protein_resid"]
    clinical = load_prospective_clinical()
    cov = clinical.drop_duplicates("participant").set_index("participant")[["purity", "cin"]].apply(
        pd.to_numeric, errors="coerce"
    )
    from mirna_hallmark import cptac_batch as B
    cov, bcols = B.augment_cov(cov, "prospective", batch_kind)
    if bcols:
        print(f"[resid_survivors] +batch ({batch_kind}): {len(bcols)} dummies -> {out_dir}")
    edges = D.high_evidence_edges(D.load_hallmark_edges())
    lit = pd.read_csv(BREAST_CTX, sep="\t") if BREAST_CTX.exists() else pd.DataFrame()

    cand = arm_couplings(survivors, arms, resid, cov, edges, lit)
    cand = cand.merge(survivors[["gene_name", "best_resid_rho", "highev"]], on="gene_name", how="left")

    # Driver = most negative FDR-surviving arm per gene (fallback: most negative arm).
    drivers_rows = []
    for gene, g in cand.groupby("gene_name"):
        sig = g[(g["partial_q"] < 0.1) & (g["partial_rho"] < 0)]
        pick = (sig if not sig.empty else g).sort_values("partial_rho").iloc[0]
        drivers_rows.append({
            "gene_name": gene,
            "driver_miRNA": pick["miRNA"],
            "driver_partial_rho": pick["partial_rho"],
            "driver_partial_q": pick["partial_q"],
            "driver_fdr_sig": bool(not sig.empty),
            "n_candidate_arms": int(len(g)),
            "literature_status": pick["literature_status"],
            "n_functional_mti_studies": pick["n_functional_mti_studies"],
            "n_breast_pmids": pick["n_breast_pmids"],
            "gene_resid_rho": pick.get("best_resid_rho", np.nan),
        })
    drivers = pd.DataFrame(drivers_rows).sort_values("driver_partial_rho")

    subtype = subtype_coupling(drivers, arms, resid, clinical)
    # Tag each driver with the PAM50 carrying the strongest (most negative) coupling.
    if not subtype.empty:
        top = subtype.sort_values("spearman_rho").groupby("gene_name").first()["pam50"]
        drivers["top_pam50"] = drivers["gene_name"].map(top)

    cand.to_csv(out_dir / "survivor_arm_candidates.tsv", sep="\t", index=False)
    drivers.to_csv(out_dir / "survivor_drivers.tsv", sep="\t", index=False)
    subtype.to_csv(out_dir / "survivor_subtype.tsv", sep="\t", index=False)

    n_sig = int(drivers["driver_fdr_sig"].sum())
    lit_counts = drivers["literature_status"].value_counts().to_dict()
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.cptac_resid_survivors",
        "cohort": "prospective (CPTAC-2)",
        "n_survivor_genes": int(len(survivors)),
        "n_with_fdr_driver_arm": n_sig,
        "driver_test": "partial Spearman(arm expression, gene protein_resid) adj purity+CIN"
                       + (f" + batch({batch_kind})" if batch_kind != "none" else "") + "; expect negative",
        "batch_kind": batch_kind,
        "literature_layers": "miRTarBase functional-MTI (n_functional_mti_studies) + MH-31 breast-context PMIDs",
        "literature_status_counts": lit_counts,
        "outputs": ["survivor_arm_candidates.tsv", "survivor_drivers.tsv", "survivor_subtype.tsv"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"[resid_survivors] wrote {out_dir}")
    print(f"  drivers with FDR arm: {n_sig}/{len(drivers)} | literature: {lit_counts}")
    print(drivers.head(12)[["gene_name", "driver_miRNA", "driver_partial_rho", "driver_partial_q",
                            "literature_status", "top_pam50"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=SURVIVOR_DIR)
    ap.add_argument("--batch-kind", default="none", choices=["none", "site", "plex", "auto"],
                    help="Add MS-plex/site batch dummies to purity+CIN (writes to batch_<kind>/ subdir)")
    ap.add_argument("--min-purity", type=float, default=None, help="Accepted for runner compatibility (no-op)")
    args = ap.parse_args()
    run(out_dir=args.out_dir, batch_kind=args.batch_kind)


if __name__ == "__main__":
    main()
