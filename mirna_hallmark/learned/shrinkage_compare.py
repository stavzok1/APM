"""Ridge-vs-lasso at full scale — does L2 shrinkage beat L1 selection on OOF coupling, and does any
advantage concentrate on COOPERATIVE-FAMILY (collinear) genes? (2026-07-07, chasing the eb_shrink lead.)

Motivation: the shipped model is a non-negative adaptive **lasso** (L1 selection). On ZEB1 — a tight
cooperative seed-family gene — Gaussian shrinkage beat it (`eb_shrink`, spike_slab). But that win
conflated two axes: **penalty type** (L1 vs L2) and **non-negativity** (the 'miRNAs repress' box). This
module disentangles them and sweeps genome-wide.

For every gene, on C-residualised, family-collapsed, adaptively-weighted predictors, run 5-fold OOF
coupling for a grid of estimators that differ ONLY in (l1_ratio, positive):
    nnlasso   l1=1.00 positive=True   ← the shipped model (L1, repress-box)
    nnen50    l1=0.50 positive=True   ← elastic-net midpoint
    nnridge   l1=0.02 positive=True   ← near-ridge, still non-negative
    uncridge  ridge   positive=False  ← unconstrained ridge (the eb_shrink flavour; negativity ALLOWED)
Regularisation strength is CV-selected per gene per estimator (inner CV on the training folds → no leakage
into the OOF coupling), so the comparison is fair on λ. Each gene also gets a COLLINEARITY score (mean
|Spearman| among its residualised family predictors) to test the interaction: Δρ(ridge−lasso) vs collinearity.

CLI:  .venv/bin/python3 -m mirna_hallmark.learned.shrinkage_compare            # hub panel
      .venv/bin/python3 -m mirna_hallmark.learned.shrinkage_compare --full     # genome-wide (background)
"""
from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM

warnings.filterwarnings("ignore", category=ConvergenceWarning)

HUB_PANEL = ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]
# All estimators share the SAME leak-free harness (train-only resid, adaptive w-scaling, raw-X aggregate),
# differing only in HOW the shrinkage is chosen and whether the repress-box (M≥0) is enforced:
#   nnlasso     L1 selection, λ by CV, M≥0                     ← the shipped model (baseline)
#   nnridge_cv  L2 shrinkage, λ by CV, M≥0
#   ebridge_nn  L2 shrinkage, λ LEARNED (half-normal EB, ν² sampled), M≥0   ← learned-τ², biology intact
#   ebridge_unc L2 shrinkage, λ LEARNED (MacKay type-II ML / evidence), free sign  ← learned-τ², eb_shrink flavour
#   blasso      LAPLACE prior, per-edge τ_j² + global λ² LEARNED, M≥0   ← Bayesian lasso: sparse+learned-shrink
#   ssridge     INCLUSION z~Bernoulli(π) + GAUSSIAN(ridge) slab + learned ν², M≥0 ← spike-slab-ridge
#   ssblasso    INCLUSION z~Bernoulli(π) + Laplace slab + learned λ², M≥0  ← spike-slab-lasso: selection+bayes
# → full 2×2: {dense=ebridge_nn/blasso, inclusion=ssridge/ssblasso} × {ridge slab, lasso slab}
CONFIGS = ["nnlasso", "nnridge_cv", "ebridge_nn", "blasso", "ssridge", "ssblasso"]


def _resid(V: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    return V - LinearRegression().fit(Cmat, V).predict(Cmat)


def _eb_ridge_unc(Xs: np.ndarray, y: np.ndarray, *, n_iter: int = 200) -> np.ndarray:
    """Unconstrained ridge with the shrinkage LEARNED by evidence maximization (MacKay type-II ML): iterate
    prior precision α (=1/τ²) and noise precision β from the effective d.o.f. γ. This is the learned-τ²
    mechanism in closed form (SVD) — no MCMC. Returns coef (free sign)."""
    U, s, Vt = np.linalg.svd(Xs, full_matrices=False)
    d2 = s ** 2
    Uty = U.T @ y
    n = Xs.shape[0]
    alpha = 1.0
    beta = 1.0 / (float(np.var(y)) + 1e-9)                       # noise precision
    m = np.zeros(Xs.shape[1])
    for _ in range(n_iter):
        denom = alpha + beta * d2
        coef_svd = beta * s * Uty / denom
        m = Vt.T @ coef_svd
        gamma = float(np.sum(beta * d2 / denom))                # effective number of parameters
        mm = float(m @ m)
        alpha_new = gamma / mm if mm > 1e-12 else alpha
        resid = y - Xs @ m
        rss = float(resid @ resid)
        beta_new = (n - gamma) / rss if rss > 1e-12 and n > gamma else beta
        if abs(alpha_new - alpha) < 1e-6 * alpha and abs(beta_new - beta) < 1e-6 * beta:
            alpha, beta = alpha_new, beta_new
            break
        alpha, beta = max(alpha_new, 1e-8), max(beta_new, 1e-8)
    return m


def _eb_ridge_nn(Xs: np.ndarray, y: np.ndarray, *, seed: int = 0) -> np.ndarray:
    """Non-negative ridge with the slab variance LEARNED (half-normal EB): the spike-and-slab Gibbs with
    inclusion forced ON for every edge (π≡1) and ν² sampled → a Bayesian NN-ridge, learned shrinkage,
    repress-box intact. Returns coef (≥0)."""
    from mirna_hallmark.learned import spike_slab as SS
    beta, _ = SS._gibbs_ss(Xs, y, np.ones(Xs.shape[1]), n_iter=1000, burn=300, seed=seed)
    return beta


def _fit_M(Xtr_res, ytr_res, w, method) -> np.ndarray:
    """Adaptive (column-scaled by w) fit on TRAIN-residualised (Xtr_res, ytr_res); return M on the raw-X
    scale (M = w·coef). ytr_res = −resid(Y|C) so POSITIVE coef = repression."""
    Xs = Xtr_res * w                                            # adaptive: scale columns by prior weight
    if method == "nnlasso":
        coef = ElasticNetCV(l1_ratio=[1.0], positive=True, cv=4, alphas=40, max_iter=20000).fit(Xs, ytr_res).coef_
    elif method == "nnridge_cv":
        coef = ElasticNetCV(l1_ratio=[0.02], positive=True, cv=4, alphas=40, max_iter=20000).fit(Xs, ytr_res).coef_
    elif method == "ebridge_unc":
        coef = _eb_ridge_unc(Xs, ytr_res)
    elif method == "ebridge_nn":
        coef = _eb_ridge_nn(Xs, ytr_res)
    elif method == "blasso":
        from mirna_hallmark.learned import spike_slab as SS
        coef, _ = SS._gibbs_blasso(Xs, ytr_res, n_iter=1000, burn=300, seed=0)
    elif method == "ssridge":
        from mirna_hallmark.learned import spike_slab as SS
        pi = SS.inclusion_prior(w)                             # inclusion + GAUSSIAN(half-normal) slab, learned ν²
        coef, _ = SS._gibbs_ss(Xs, ytr_res, pi, n_iter=1000, burn=300, seed=0)
    elif method == "ssblasso":
        from mirna_hallmark.learned import spike_slab as SS
        pi = SS.inclusion_prior(w)                             # evidence-graded inclusion from the prior weight
        coef, _, _ = SS._gibbs_ss_blasso(Xs, ytr_res, pi, n_iter=1000, burn=300, seed=0)
    else:
        raise ValueError(method)
    return w * coef                                            # unscale → coefficient on raw-X units


def _collinearity(Xr: np.ndarray) -> float:
    """Mean |Spearman| over predictor pairs (residualised X). High ⇒ cooperative/collinear predictor set."""
    p = Xr.shape[1]
    if p < 2:
        return 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        R = spearmanr(Xr).correlation
    if np.ndim(R) == 0 or R is None:
        return 0.0
    iu = np.triu_indices(p, 1)
    v = np.abs(np.asarray(R)[iu])
    v = v[np.isfinite(v)]
    return float(v.mean()) if v.size else 0.0


def coupling_row(gene: str, *, folds: int = 5, seed: int = 0, w_prior_source: str = "ledger") -> dict:
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source=w_prior_source)
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    Cmat = C.to_numpy(dtype=float)
    Yv = Y.to_numpy(dtype=float)
    Xv = X.to_numpy(dtype=float)
    yr = _resid(Yv, Cmat)                                       # full-data resid ONLY for the final scoring
    wv = np.clip(w.reindex(X.columns).fillna(0.0).to_numpy(dtype=float), 1e-3, None)
    n, p = Xv.shape
    row = {"gene": gene, "n": n, "n_pred": p, "collin": round(_collinearity(_resid(Xv, Cmat)), 3)}
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    preds = {k: np.full(n, np.nan) for k in CONFIGS}
    for tr, te in kf.split(Xv):
        Ctr = Cmat[tr]                                          # TRAIN-ONLY residualisation (no leak; matches fit_gene)
        ytr = -_resid(Yv[tr], Ctr)                             # repression → positive coefficient
        Xtr = _resid(Xv[tr], Ctr)
        for k in CONFIGS:
            M = _fit_M(Xtr, ytr, wv, k)
            preds[k][te] = Xv[te] @ M                          # aggregate on RAW held-out X (like mvp.oof_gate)
    for k in CONFIGS:
        rho = spearmanr(_resid(preds[k], Cmat), yr).correlation  # expect <0 (pressure ↓ target)
        row[f"rho_{k}"] = round(float(rho), 3) if rho == rho else np.nan
    row["d_ebnn_lasso"] = (round(row["rho_ebridge_nn"] - row["rho_nnlasso"], 3)
                           if pd.notna(row.get("rho_ebridge_nn")) and pd.notna(row.get("rho_nnlasso")) else np.nan)
    return row


def run(genes=None, *, folds: int = 5) -> pd.DataFrame:
    genes = genes or HUB_PANEL
    rows = []
    for g in genes:
        try:
            rows.append(coupling_row(g, folds=folds))
        except Exception as e:
            rows.append({"gene": g, "error": repr(e)[:70]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170):
        print(df.to_string(index=False))
    _summarize(df)
    return df


def _summarize(df: pd.DataFrame) -> None:
    d = df[df.get("rho_nnlasso").notna()] if "rho_nnlasso" in df else df.iloc[0:0]
    if not len(d):
        return
    print("\nmean OOF coupling (more negative = better):")
    for k in CONFIGS:
        c = f"rho_{k}"
        if c in d:
            print(f"  {k:12s} {d[c].mean():+.4f}   beats-nnlasso {int((d[c] < d['rho_nnlasso']).sum())}/{len(d)}")
    if "d_ebnn_lasso" in d and "collin" in d:                   # learned-τ² NN-ridge vs lasso, by collinearity
        dd = d[d["d_ebnn_lasso"].notna()]
        if len(dd) > 3:
            rho_int = spearmanr(dd["collin"], dd["d_ebnn_lasso"]).correlation  # Δ<0 = EB-ridge better
            print(f"\ninteraction: Spearman( collinearity , Δρ(EBridge_nn − lasso) ) = {rho_int:+.3f}")
            print("  (NEGATIVE ⇒ learned-τ² NN-ridge advantage GROWS with collinearity — cooperative-family test)")
            hi = dd[dd["collin"] >= dd["collin"].median()]
            lo = dd[dd["collin"] < dd["collin"].median()]
            print(f"  high-collin: mean Δρ(EBnn−lasso) {hi['d_ebnn_lasso'].mean():+.4f} (n={len(hi)}) | "
                  f"low-collin {lo['d_ebnn_lasso'].mean():+.4f} (n={len(lo)})")


def full(*, out: str = "mirna_hallmark/output/learned/shrinkage_compare.tsv",
         limit: int | None = None, progress: int = 50) -> pd.DataFrame:
    from pathlib import Path
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[shrinkage_compare] {i}/{len(genes)}", flush=True)
        try:
            rows.append(coupling_row(g))
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"\n=== GENOME-WIDE RIDGE vs LASSO ({len(df)} genes) → {out} ===")
    _summarize(df)
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Ridge vs lasso OOF-coupling sweep")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--genes", nargs="*", default=None)
    a = ap.parse_args()
    if a.full:
        full(limit=a.limit)
    else:
        run(a.genes)
