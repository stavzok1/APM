"""THE GENE LENS — per-gene competence, a priori. Can this gene's miRNA regulation be measured AT ALL?

MH-130 built this and it is load-bearing (`STATE_OF_PLAY.md` Axis 4 cites it): **17.6% of genes are NOT
MEASURABLE**, **48.1% have ONE seed family** (where β ≡ uniform exactly, so "does the learning help?" is a
NON-QUESTION), and the model's honest **domain of validity is 27% of the universe** (`n_fam ≥ 3` AND
`w_max > median`) — a rule that is **fit-free**: neither axis touches Y, β, or the decoy.

⚠⚠ **WHY THIS MODULE EXISTS (MH-144, 2026-07-17).** MH-130 ran as **scratchpad code that was never
committed** — MH-127 persisted its builders under `output/learned/mh127/code/`, MH-130 did not. So
`gene_atlas.tsv` / `gene_competence_map.tsv` were **UNREPRODUCIBLE**: the numbers above could be cited but
not regenerated, and `output/` is gitignored, so they existed on exactly one disk. That is a provenance
hole, not a filing detail — an artifact nobody can rebuild cannot be audited, corrected, or refreshed when
its inputs move.

⭐ **AND THE REBUILD FIXES A SECOND DEFECT.** MH-130's competence map carried 20 decoy-dependent columns
(`*__FAKE1` / `*__FAKE2`) computed on the decoy that **MH-135/137 later discredited** (6mer-contaminated,
purity ANTI-matched, and a 126× evidence hole that let CLIP-supported REAL edges serve as "decoys"). Those
columns are known-bad. This module therefore **does not reimplement them**: the ATLAS half is computed here,
and the gap columns are joined from **`eval/decoy_bench.py`**, the canonical decoy (MH-137: core gap
−0.0119; MH-139: −0.0129 [−0.0205,−0.0052] caliper-free, the two agreeing). One decoy, one home.

WHAT IS COMPUTED HERE (all decoy-independent, all a priori):
  * `ceiling_R2_oof_*` — 5-fold OOF R² of an **UNRESTRICTED OLS** on the gene's full real regulator design.
    This is the **upper bound for ANY linear weighting** of those regulators (real β, fake β, dose, uniform),
    so `ceiling ≤ 0` ⇒ **no weighting can explain anything out-of-fold** ⇒ the gene is not measurable, and no
    real-vs-fake contrast there is interpretable. C / z-scoring / dose floor are all **TRAIN-only** (fitting
    them on all samples leaks the test fold into the standardisation).
  * `n_fam`, `n_arms`, `w_max` — design width and curated evidence depth.
  * `apriori_class` — A_COMPETENT (`n_fam≥3` ∧ `w_max>median`) · B_PARTIAL · C_WEAK · D_UNDEFINED.

⚠ **The ceiling is DESIGN CAPACITY, not biology** (MH-130): spearman(ceiling, n_fam) = +0.564. A gene with a
low ceiling has too few/too collinear regulators to resolve, not necessarily no regulation. Target detection
is NOT the binding constraint (the detection floor fails only 7.9% of genes) — **regulator-design width is**.

CLI:  .venv/bin/python3 -m mirna_hallmark.learned.analyses.gene_atlas [--workers 8] [--limit N]
Out:  output/learned/gene_atlas.tsv           (the atlas — decoy-independent)
      output/learned/gene_competence_map.tsv  (atlas ⋈ decoy_bench gaps)
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

ROOT = Path(__file__).resolve().parents[2]
OUT_ATLAS = ROOT / "mirna_hallmark/output/learned/gene_atlas.tsv"
OUT_MAP = ROOT / "mirna_hallmark/output/learned/gene_competence_map.tsv"
DECOY_BENCH = ROOT / "mirna_hallmark/output/learned/decoy_bench.tsv"

N_SPLITS, SEED = 5, 0          # the folds MH-130's ladder used


def _oof_ceiling(y: np.ndarray, X: np.ndarray) -> tuple:
    """5-fold OOF R² + rho of an UNRESTRICTED OLS — the upper bound for any linear weighting of X.

    ⚠ Everything scale-dependent is fit on TRAIN ONLY. z-scoring on all samples would leak the held-out
    fold into the standardisation and inflate the ceiling — which is exactly the bound we are trying to
    measure honestly.
    """
    n, p = X.shape
    if n < N_SPLITS * 2 or p < 1:
        return np.nan, np.nan, np.nan, np.nan
    pred = np.full(n, np.nan)
    for tr, te in KFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED).split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0)
        keep = sd >= 0.1                                  # TRAIN-only dose floor
        if not keep.any():
            continue
        Z = (X[:, keep] - mu[keep]) / (sd[keep] + 1e-9)
        A = np.column_stack([np.ones(len(tr)), Z[tr]])
        try:
            coef, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        except np.linalg.LinAlgError:
            continue
        pred[te] = np.column_stack([np.ones(len(te)), Z[te]]) @ coef
    m = np.isfinite(pred)
    if m.sum() < 10 or np.std(pred[m]) < 1e-12:
        return np.nan, np.nan, np.nan, np.nan
    ss_res = float(((y[m] - pred[m]) ** 2).sum())
    ss_tot = float(((y[m] - y[m].mean()) ** 2).sum())
    r2_oof = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    from scipy.stats import spearmanr
    rho = float(spearmanr(pred[m], y[m]).statistic)
    # in-sample R² / adj-R² on the full design, for the in-vs-out gap
    mu, sd = X.mean(0), X.std(0)
    keep = sd >= 0.1
    r2_in = r2adj = np.nan
    if keep.any():
        Z = (X[:, keep] - mu[keep]) / (sd[keep] + 1e-9)
        A = np.column_stack([np.ones(n), Z])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ coef
        sst = float(((y - y.mean()) ** 2).sum())
        r2_in = 1.0 - float((res ** 2).sum()) / sst if sst > 0 else np.nan
        q = int(keep.sum())
        r2adj = 1.0 - (1.0 - r2_in) * (n - 1) / max(n - q - 1, 1) if r2_in == r2_in else np.nan
    return r2_oof, rho, r2_in, r2adj


def _one(gene: str) -> Optional[dict]:
    from mirna_hallmark.learned import data as LD, families as FAM
    from mirna_hallmark.learned.evidence import ledger as LG
    o = {"gene": gene}
    for tag, deconv in (("core", False), ("comp", True)):
        try:
            Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=deconv)
        except Exception:
            return None
        if X.shape[1] < 1:
            return None
        fam = FAM.family_of(pd.Index(X.columns))
        Xf, wf, _ = FAM.collapse_by_family(X, w, fam)
        # C-residualise the target, matching the canonical path (`AE._prep`): r = -resid(Y|C)
        from mirna_hallmark.learned.attribution_eb import _resid
        y = -_resid(Y.to_numpy(float), C.to_numpy(float))
        if tag == "core":
            o.update({"n_arms": int(X.shape[1]), "n_fam": int(Xf.shape[1]), "n_eff": int(len(Y)),
                      "w_max": float(np.nanmax(w.to_numpy(float))) if len(w) else np.nan,
                      "y_mean": float(np.mean(Y)), "y_sd": float(np.std(Y)),
                      "y_dropout": float((Y.to_numpy(float) <= 0).mean())})
        o[f"n_C_{tag}"] = int(C.shape[1])
        for lvl, M in (("fam", Xf), ("arm", X)):
            r2o, rho, r2i, r2a = _oof_ceiling(y, M.to_numpy(float))
            o[f"ceiling_R2_oof_{lvl}_{tag}"] = r2o
            o[f"ceiling_rho_oof_{lvl}_{tag}"] = rho
            o[f"ceiling_R2_{lvl}_{tag}"] = r2i
            o[f"ceiling_R2adj_{lvl}_{tag}"] = r2a
    return o


def build(genes=None, *, workers: int = 8, limit: Optional[int] = None) -> pd.DataFrame:
    from mirna_hallmark.learned.evidence import ledger as LG
    if genes is None:
        he = LG.pooled_he_edges()
        genes = sorted(he["gene"].unique())
    if limit:
        genes = genes[:limit]
    rows = []
    if workers and workers > 1:
        import multiprocessing as mp
        with mp.Pool(workers) as pool:
            for i, r in enumerate(pool.imap_unordered(_one, genes, chunksize=4)):
                if i % 100 == 0:
                    print(f"[gene_atlas] {i}/{len(genes)}", flush=True)
                if r:
                    rows.append(r)
    else:
        rows = [r for r in (_one(g) for g in genes) if r]
    A = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)

    # ---- the A-PRIORI competence class (fit-free: no Y, no β, no decoy) ----
    wmed = A["w_max"].median()
    A["apriori_class"] = np.select(
        [(A.n_fam >= 3) & (A.w_max > wmed),
         (A.n_fam >= 3) & (A.w_max <= wmed),
         (A.n_fam == 2)],
        ["A_COMPETENT", "C_WEAK", "B_PARTIAL"], default="D_UNDEFINED")
    A["measurable"] = A["ceiling_R2_oof_fam_core"] > 0
    return A


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    A = build(workers=a.workers, limit=a.limit)
    A.to_csv(OUT_ATLAS, sep="\t", index=False)

    # ---- the competence MAP = atlas ⋈ the CANONICAL decoy (never MH-130's discredited FAKE1/FAKE2) ----
    if DECOY_BENCH.exists():
        d = pd.read_csv(DECOY_BENCH, sep="\t")
        keep = [c for c in ["gene", "gap_core", "gap_deconv", "real_core", "dec_core", "d_dose"] if c in d.columns]
        M = A.merge(d[keep], on="gene", how="left")
    else:
        print(f"[gene_atlas] ⚠ {DECOY_BENCH} absent — emitting the atlas half only")
        M = A.copy()
    M.to_csv(OUT_MAP, sep="\t", index=False)

    n = len(A)
    print(f"\n=== THE GENE LENS — {n} genes ===")
    print(f"  NOT MEASURABLE (ceiling_R2_oof ≤ 0, family/core): "
          f"{int((A.ceiling_R2_oof_fam_core <= 0).sum())}/{n} = {(A.ceiling_R2_oof_fam_core <= 0).mean():.1%}")
    print(f"  ceiling ≤ 0.02                                  : {(A.ceiling_R2_oof_fam_core <= 0.02).mean():.1%}")
    print(f"  ONE seed family (β ≡ uniform ⇒ a NON-question)  : "
          f"{int((A.n_fam == 1).sum())}/{n} = {(A.n_fam == 1).mean():.1%}")
    print(f"  median ceiling                                  : {A.ceiling_R2_oof_fam_core.median():+.4f}")
    print(f"  ⚠ ceiling is DESIGN CAPACITY: spearman(ceiling, n_fam) = "
          f"{A[['ceiling_R2_oof_fam_core','n_fam']].corr(method='spearman').iloc[0,1]:+.3f}")
    print(f"\n  A-PRIORI CLASS (fit-free — no Y, no β, no decoy):")
    for k, v in A.apriori_class.value_counts().items():
        print(f"    {k:14s} {v:5d}  ({v/n:.1%})")
    print(f"\n  -> {OUT_ATLAS}\n  -> {OUT_MAP}")


if __name__ == "__main__":
    main()
