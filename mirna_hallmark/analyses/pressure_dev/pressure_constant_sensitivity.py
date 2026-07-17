"""Sensitivity of the coupling headline to the hand-tuned weighting constants.

The evidence numerator weights (reporter 3.0, protein 2.5, rna/perturbation 1.5,
binding 0.5, weak-discount 0.3) and the AGO-load weights (AGO2 2.0, …) are
transparent priors, not fitted values. This module perturbs them and reports
whether the realized-coupling headline moves -- so "we used 3.0" is backed by a
stability band rather than asserted.

Two sweeps, each scored by gene-level coupling (mean confounder-adjusted partial
Spearman of aggregate pressure vs each gene's expression, and the count of
strong-negative genes, ρ < -0.1) over the Hallmark universe, full cohort
(tumour, Normal-like excluded):

A. **Evidence weights** -- random log-normal jitter (±~25%) of all six weights, plus
   a systematic ±50% on each one in turn; edges re-scored, pressure rebuilt.
B. **AGO-load weights / co-limitation** -- AGO2 weight grid × effector co-limit on/off;
   the gate is recomputed and applied (gated coupling).

Run: ``.venv/bin/python3 -m mirna_hallmark.pressure_constant_sensitivity [--n-jitter N] [--smoke]``.
Outputs under ``output/pressure_constant_sensitivity/``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import pressure_build as PB
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.coupling_inference import partial_spearman_matrix
from mirna_hallmark.evidence_scoring import _COL, _num
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna

OUT_DIR = C.OUTPUT_ROOT / "pressure_constant_sensitivity"

BASE_WEIGHTS = {"reporter": 3.0, "protein": 2.5, "rna": 1.5, "perturbation": 1.5, "binding": 0.5}
BASE_WEAK = 0.3


def _score_confidence(df: pd.DataFrame, w: Dict[str, float], weak: float) -> pd.Series:
    """confidence_logclass evidence score with parametrized weights."""
    func = (
        w["reporter"] * np.log1p(_num(df, _COL["rf"]))
        + w["protein"] * np.log1p(_num(df, _COL["pf"]))
        + w["rna"] * np.log1p(_num(df, _COL["rna_f"]))
        + w["perturbation"] * np.log1p(_num(df, _COL["pert_f"]))
        + w["binding"] * np.log1p(_num(df, _COL["bf"]))
    )
    weak_term = (
        w["reporter"] * np.log1p(_num(df, _COL["rw"]))
        + w["protein"] * np.log1p(_num(df, _COL["pw"]))
        + w["rna"] * np.log1p(_num(df, _COL["rna_w"]))
        + w["perturbation"] * np.log1p(_num(df, _COL["pert_w"]))
        + w["binding"] * np.log1p(_num(df, _COL["bw"]))
    )
    return func + weak * weak_term


def _coupling_score(gp: pd.DataFrame, rna: pd.DataFrame, clinical: pd.DataFrame) -> Dict[str, float]:
    g = [x for x in gp.index if x in rna.index]
    if not g:
        return {"mean_rho": np.nan, "n_strong_neg": 0, "frac_neg": np.nan}
    clin = clinical.set_index("participant")
    cols = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    samples = gp.columns.intersection(rna.columns)
    cov = D.augment_tcga_batch(clin.reindex(samples)[cols]).reindex(samples) if cols else None
    rho = partial_spearman_matrix(gp.loc[g], rna.loc[g], covariates=cov).dropna()
    return {"mean_rho": float(rho.mean()), "n_strong_neg": int((rho < -0.1).sum()),
            "frac_neg": float((rho < 0).mean())}


def run(*, out_dir: Path = OUT_DIR, n_jitter: int = 12, smoke: bool = False) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    hs = HallmarkSets.load()
    genes = set(hs.universe)
    clinical = D.load_clinical_strata()
    mirna = D.load_mirna_arms(); rna = D.load_rna()
    mirna = mirna[~mirna.index.duplicated(keep="first")]
    rna = rna[~rna.index.duplicated(keep="first")]
    cohort = sorted(set(clinical["participant"]) & set(mirna.columns) & set(rna.columns))
    mirna, rna = mirna[cohort], rna[cohort]
    print(f"[sens] cohort {len(cohort)} participants")

    summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
    summary = summary[summary["gene"].isin(genes)
                      & summary["miRNA"].astype(str).str.startswith("hsa-", na=False)].copy()

    def _pressure_for_weights(w: Dict[str, float], weak: float) -> pd.DataFrame:
        s = summary.copy()
        s["evidence_score"] = _score_confidence(s, w, weak)
        s = s[pd.to_numeric(s["evidence_score"], errors="coerce").fillna(0) >= C.PRESSURE_MIN_EVIDENCE]
        e = s[["miRNA", "gene", "evidence_score"]].drop_duplicates()
        resolved, _ = resolve_edges_mirna(e, mirna, weight_col="evidence_score")
        return PB.compute_gene_pressure(list(genes), edges=resolved, mirna=mirna, resolve_arms=False)

    # ---- baseline -----------------------------------------------------------
    gp_base = _pressure_for_weights(BASE_WEIGHTS, BASE_WEAK)
    base = _coupling_score(gp_base, rna, clinical)
    print(f"[sens] baseline: mean_rho={base['mean_rho']:.4f} n_strong_neg={base['n_strong_neg']}")

    rows: List[dict] = []
    rows.append({"sweep": "baseline", "perturbation": "default", **base})

    # ---- A. evidence-weight jitter + systematic ±50% ------------------------
    n_j = 3 if smoke else n_jitter
    for j in range(n_j):
        w = {k: v * float(np.exp(rng.normal(0, 0.25))) for k, v in BASE_WEIGHTS.items()}
        weak = BASE_WEAK * float(np.exp(rng.normal(0, 0.25)))
        sc = _coupling_score(_pressure_for_weights(w, weak), rna, clinical)
        rows.append({"sweep": "evidence_jitter", "perturbation": f"jitter_{j}", **sc})
    if not smoke:
        for k in list(BASE_WEIGHTS) + ["weak"]:
            for mult, tag in ((0.5, "x0.5"), (1.5, "x1.5")):
                w = dict(BASE_WEIGHTS); weak = BASE_WEAK
                if k == "weak":
                    weak = BASE_WEAK * mult
                else:
                    w[k] = BASE_WEIGHTS[k] * mult
                sc = _coupling_score(_pressure_for_weights(w, weak), rna, clinical)
                rows.append({"sweep": "evidence_systematic", "perturbation": f"{k}_{tag}", **sc})

    # ---- B. AGO-load weights / co-limitation (gated coupling) ---------------
    ago2_grid = [2.0] if smoke else [1.0, 1.5, 2.0, 3.0]
    for a2 in ago2_grid:
        for colimit in ([True] if smoke else [True, False]):
            weights = (("AGO1", 1.0), ("AGO2", a2), ("AGO3", 1.0), ("AGO4", 0.5))
            params = replace(C.AGO_GATE, ago_load_weights=weights, effector_colimit=colimit)
            gate = compute_ago_gate(rna, params=params)["ago_gate"]
            shared = gp_base.columns.intersection(gate.index)
            gp_gated = gp_base[shared].mul(gate.reindex(shared), axis=1)
            sc = _coupling_score(gp_gated, rna, clinical)
            rows.append({"sweep": "ago_weights", "perturbation": f"AGO2={a2},colimit={colimit}", **sc})

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "sensitivity_scores.tsv", sep="\t", index=False)

    # ---- stability summary --------------------------------------------------
    def _band(sub: pd.DataFrame, col: str) -> Dict[str, float]:
        v = sub[col].dropna()
        return {"min": round(float(v.min()), 4), "max": round(float(v.max()), 4),
                "mean": round(float(v.mean()), 4), "cv": round(float(v.std() / abs(v.mean())), 4) if v.mean() else np.nan}

    ev = df[df["sweep"].isin(["evidence_jitter", "evidence_systematic"])]
    ago = df[df["sweep"] == "ago_weights"]
    base_sn = base["n_strong_neg"]
    summary_out = {
        "module": "mirna_hallmark.pressure_constant_sensitivity",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_cohort": len(cohort), "smoke": smoke,
        "baseline": {k: round(v, 4) if isinstance(v, float) else v for k, v in base.items()},
        "evidence_weights": {
            "mean_rho_band": _band(ev, "mean_rho"),
            "n_strong_neg_band": _band(ev, "n_strong_neg"),
            "n_strong_neg_pct_of_baseline": (
                [round(float(x / base_sn), 3) for x in (ev["n_strong_neg"].min(), ev["n_strong_neg"].max())]
                if base_sn else None
            ),
        },
        "ago_weights": {
            "mean_rho_band": _band(ago, "mean_rho"),
            "n_strong_neg_band": _band(ago, "n_strong_neg"),
        },
    }
    # headline robust if the coupling count stays within ±~15% across evidence perturbations
    sn_band = summary_out["evidence_weights"]["n_strong_neg_band"]
    summary_out["verdict_headline_robust_to_constants"] = bool(
        base_sn and (sn_band["min"] >= 0.8 * base_sn) and (sn_band["max"] <= 1.25 * base_sn)
    )
    (out_dir / "summary.json").write_text(json.dumps(summary_out, indent=2), encoding="utf-8")
    print(json.dumps(summary_out, indent=2))
    print(f"[sens] verdict: headline "
          f"{'ROBUST' if summary_out['verdict_headline_robust_to_constants'] else 'SENSITIVE (inspect)'}"
          f" to the tuned constants; wrote {out_dir}")
    return summary_out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-jitter", type=int, default=12)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, n_jitter=args.n_jitter, smoke=args.smoke)


if __name__ == "__main__":
    main()
