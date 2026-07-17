"""Held-out validation of the pressure-construction / denominator choice.

The spine fixes the expression construction (``softmax_z_logrpm``) and the
promiscuity denominator (``evidence_mass``) partly because they give the sharpest
realized coupling. Chosen and *validated* on the same cohort, that is circular:
the construction could be overfit to this cohort's coupling. This module breaks the
circularity with **patient cross-validation** -- the construction is selected on a
TRAIN fold and scored on a held-out TEST fold (with pressure scope-recomputed
*within* each split, so the test predictor never sees train samples). If the spine
ranks at/near the top **out-of-sample**, the choice generalizes; if a different
construction wins out-of-sample, the in-sample selection was overfit.

Score = gene-level coupling strength: the mean confounder-adjusted partial Spearman
(CPE+HRD+batch) of aggregate pressure vs each gene's own expression, over the
Hallmark universe (more negative = stronger realized repression), plus the fraction
negative and the count of strong-negative genes.

Run: ``.venv/bin/python3 -m mirna_hallmark.held_out_tuning [--folds K] [--smoke]``.
Outputs under ``output/held_out_tuning/``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import pressure_build as PB
from mirna_hallmark.coupling_inference import partial_spearman_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets

OUT_DIR = C.OUTPUT_ROOT / "held_out_tuning"

# (expr_mode, target_norm); the spine is first. Sweeps both the denominator (the
# §3.3 concern) and the expression construction.
CANDIDATES: List[Tuple[str, str]] = [
    ("softmax_z_logrpm", "evidence_mass"),    # SPINE
    ("softmax_z_logrpm", "none"),
    ("softmax_z_logrpm", "degree"),
    ("softmax_z_logrpm", "evidence_sum_log"),
    ("softmax_z_logrpm", "evidence_outer_log"),
    ("softmax_z", "evidence_mass"),
    ("softmax_logrpm", "evidence_mass"),
    ("z", "evidence_mass"),
]
SPINE = CANDIDATES[0]


def _covariates(clinical: pd.DataFrame, samples: pd.Index) -> Optional[pd.DataFrame]:
    clin = clinical.set_index("participant")
    cols = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    if not cols:
        return None
    return D.augment_tcga_batch(clin.reindex(samples)[cols]).reindex(samples)


def _score(
    genes, edges, mirna_cols, rna, clinical, expr_mode, target_norm,
) -> Dict[str, float]:
    """Gene-level coupling score on one sample subset (scope = that subset)."""
    gp = PB.compute_gene_pressure(
        genes, edges=edges, mirna=mirna_cols, expr_mode=expr_mode,
        target_norm=target_norm, resolve_arms=False,
    )
    g = [x for x in gp.index if x in rna.index]
    if not g:
        return {"mean_rho": np.nan, "frac_neg": np.nan, "n_strong_neg": 0, "n_genes": 0}
    cov = _covariates(clinical, gp.columns.intersection(rna.columns))
    rho = partial_spearman_matrix(gp.loc[g], rna.loc[g], covariates=cov).dropna()
    return {
        "mean_rho": float(rho.mean()),
        "frac_neg": float((rho < 0).mean()),
        "n_strong_neg": int((rho < -0.1).sum()),
        "n_genes": int(len(rho)),
    }


def run(*, out_dir: Path = OUT_DIR, folds: int = 5, smoke: bool = False) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    from sklearn.model_selection import KFold

    hs = HallmarkSets.load()
    genes = list(hs.universe)
    clinical = D.load_clinical_strata()
    mirna = D.load_mirna_arms()
    rna = D.load_rna()
    mirna = mirna[~mirna.index.duplicated(keep="first")]
    rna = rna[~rna.index.duplicated(keep="first")]
    cohort = sorted(set(clinical["participant"]) & set(mirna.columns) & set(rna.columns))
    mirna = mirna[cohort]; rna = rna[cohort]
    print(f"[heldout] cohort {len(cohort)} participants (tumour, Normal-like excluded)")

    edges = D.high_evidence_edges()[["miRNA", "gene", "evidence_score"]].drop_duplicates()
    edges = PB.resolve_pressure_edges(edges, mirna)
    edges = edges[edges["gene"].isin(set(rna.index)) & edges["miRNA"].isin(set(mirna.index))]
    edges = edges.drop_duplicates(subset=["miRNA", "gene"])

    cands = CANDIDATES[:3] if smoke else CANDIDATES
    n_splits = 2 if smoke else folds
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=0)
    arr = np.asarray(cohort)

    rows = []
    for fi, (tr, te) in enumerate(kf.split(arr)):
        tr_cols, te_cols = list(arr[tr]), list(arr[te])
        print(f"[heldout] fold {fi+1}/{n_splits}: {len(tr_cols)} train / {len(te_cols)} test")
        for (em, tn) in cands:
            s_tr = _score(genes, edges, mirna[tr_cols], rna.loc[:, tr_cols], clinical, em, tn)
            s_te = _score(genes, edges, mirna[te_cols], rna.loc[:, te_cols], clinical, em, tn)
            rows.append({
                "fold": fi, "expr_mode": em, "target_norm": tn,
                "config": f"{em}|{tn}", "is_spine": (em, tn) == SPINE,
                "train_mean_rho": s_tr["mean_rho"], "test_mean_rho": s_te["mean_rho"],
                "train_n_strong_neg": s_tr["n_strong_neg"], "test_n_strong_neg": s_te["n_strong_neg"],
                "test_frac_neg": s_te["frac_neg"],
            })
    per_fold = pd.DataFrame(rows)
    per_fold.to_csv(out_dir / "per_fold.tsv", sep="\t", index=False)

    # per-fold ranks (1 = strongest = most negative mean_rho), then aggregate
    per_fold["train_rank"] = per_fold.groupby("fold")["train_mean_rho"].rank()
    per_fold["test_rank"] = per_fold.groupby("fold")["test_mean_rho"].rank()
    agg = per_fold.groupby("config").agg(
        is_spine=("is_spine", "first"),
        train_mean_rho=("train_mean_rho", "mean"),
        test_mean_rho=("test_mean_rho", "mean"),
        test_rho_sd=("test_mean_rho", "std"),
        mean_train_rank=("train_rank", "mean"),
        mean_test_rank=("test_rank", "mean"),
        mean_test_strong_neg=("test_n_strong_neg", "mean"),
    ).reset_index().sort_values("test_mean_rho")
    agg["generalization_gap"] = agg["test_mean_rho"] - agg["train_mean_rho"]
    agg.to_csv(out_dir / "per_candidate.tsv", sep="\t", index=False)

    # ---- verdict ------------------------------------------------------------
    spine_row = agg.loc[agg["is_spine"]].iloc[0]
    best_test = float(agg["test_mean_rho"].min())
    spine_test = float(spine_row["test_mean_rho"])
    spine_test_rank = float(spine_row["mean_test_rank"])
    # does selecting on train predict the test ranking?
    tt_corr = float(agg[["train_mean_rho", "test_mean_rho"]].corr(method="spearman").iloc[0, 1])
    # selection stability: per fold, is the train-winner also the test-winner?
    stab = []
    for _, f in per_fold.groupby("fold"):
        stab.append(f.loc[f["train_mean_rho"].idxmin(), "config"]
                    == f.loc[f["test_mean_rho"].idxmin(), "config"])
    # Two DISTINCT questions, kept separate (conflating them is what made the old
    # binary verdict misleading):
    #  - generalizes / not overfit: does the train ranking reproduce out-of-sample,
    #    and is the spine within noise of the OOS best? (the circularity question)
    #  - coupling-optimal: is the spine the single best construction on coupling?
    selection_generalizes = bool(
        tt_corr > 0.9 and float(np.mean(stab)) >= 0.8 and (spine_test - best_test) <= 0.01
    )
    spine_coupling_optimal = bool(round(spine_test_rank) == 1)
    if selection_generalizes and spine_coupling_optimal:
        verdict = "generalizes AND coupling-optimal"
    elif selection_generalizes:
        verdict = ("selection generalizes (NOT overfit), but spine is not coupling-maximal "
                   f"(rank {spine_test_rank:.0f}/{len(cands)}, within {spine_test-best_test:.4f} of "
                   f"'{agg.iloc[0]['config']}') — consistent with the denominator trading a little "
                   "coupling for hub-attribution control")
    else:
        verdict = "selection does not cleanly generalize — revisit"
    summary = {
        "module": "mirna_hallmark.held_out_tuning",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_folds": n_splits, "n_candidates": len(cands), "n_cohort": len(cohort), "smoke": smoke,
        "spine_config": f"{SPINE[0]}|{SPINE[1]}",
        "spine_test_mean_rho": round(spine_test, 4),
        "best_oos_test_mean_rho": round(best_test, 4),
        "best_oos_config": str(agg.iloc[0]["config"]),
        "spine_mean_test_rank": round(spine_test_rank, 2),
        "train_test_rank_spearman": round(tt_corr, 3),
        "selection_stability_train_winner_is_test_winner": round(float(np.mean(stab)), 3),
        "spine_generalization_gap": round(float(spine_row["generalization_gap"]), 4),
        "selection_generalizes_not_overfit": selection_generalizes,
        "spine_is_coupling_optimal": spine_coupling_optimal,
        "verdict": verdict,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\n[heldout] out-of-sample candidate ranking (most-negative test_mean_rho = best):")
    print(agg[["config", "is_spine", "train_mean_rho", "test_mean_rho", "mean_test_rank",
               "generalization_gap"]].round(4).to_string(index=False))
    print(f"\n[heldout] verdict: spine choice "
          f"{'HOLDS out-of-sample' if holds else 'does NOT clearly hold (revisit)'}; wrote {out_dir}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--smoke", action="store_true", help="2 folds, 3 candidates")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, folds=args.folds, smoke=args.smoke)


if __name__ == "__main__":
    main()
