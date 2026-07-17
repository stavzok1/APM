"""Empirical-Bayes learned shrinkage TOWARD THE CURATED PRIOR (the untried variant, 2026-07-07).

Everywhere else the curated evidence prior `w` is an INVERSE PENALTY (ordering/selection), never a
location — because the program's headline empirical fact is "prior sets ordering, data sets magnitude"
(METHODS §2/§3). This module tests the opposite: make the prior a **location** the coefficients shrink
toward, and let the DATA decide how hard, via a learned between-effect variance τ².

Per gene, on C-residualised, z-scored predictors Xz and target r = −resid(Y|C):

    r = Xz·β + ε ,  ε ~ N(0, σ²)
    β_m ~ N(c·ŵ_m, τ²)             # SHRINK TOWARD THE PRIOR: ŵ = unit-normalised curated weight, c a
                                    #   learned scalar mapping prior→coefficient scale
    c ~ N(0, ω²) ,  σ², τ² ~ InvGamma

τ² is the adaptive knob: **small** ⇒ data agrees with the prior's magnitude pattern ⇒ shrink hard (trust
the prior); **large** ⇒ data disagrees ⇒ shrink little (β driven by data, prior discounted). c measures
how much the prior's *direction* is used at all (c≈0 ⇒ prior magnitude pattern carries no weight).
Conjugate Gibbs (all Gaussian/InvGamma) → fast. The diagnostic is (learned τ², learned c, OOF coupling)
vs the fixed-α adaptive lasso: does *learning* the anchor strength beat *guessing* α?

CLI: `.venv/bin/python3 -m mirna_hallmark.learned.analyses.eb_shrink`
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR

HUB_PANEL = ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]


def _resid(V: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    return V - LinearRegression().fit(Cmat, V).predict(Cmat)


def _gibbs_eb(Xz: np.ndarray, r: np.ndarray, w_unit: np.ndarray, *,
              n_iter: int = 1500, burn: int = 500, omega2: float = 100.0, seed: int = 0
              ) -> tuple[np.ndarray, float, float]:
    """Conjugate Gibbs for β ~ N(c·ŵ, τ²), learned c and τ². Returns (posterior-mean β, τ²_mean, c_mean)."""
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    XtX = Xz.T @ Xz
    Xtr = Xz.T @ r
    ww = float(w_unit @ w_unit) or 1.0
    beta = np.zeros(p)
    c = 0.0
    sig2 = float(np.var(r)) or 1.0
    tau2 = 1.0
    a0 = b0 = 1.0
    acc_b = np.zeros(p); acc_t = []; acc_c = []
    for it in range(n_iter):
        # β | ·   Gaussian: V = (XᵀX/σ² + I/τ²)⁻¹ , mean = V(Xᵀr/σ² + c·ŵ/τ²)
        V = np.linalg.inv(XtX / sig2 + np.eye(p) / tau2)
        m = V @ (Xtr / sig2 + c * w_unit / tau2)
        beta = m + np.linalg.cholesky(V) @ rng.standard_normal(p)
        # c | ·   Gaussian from β ~ N(c·ŵ, τ²): prec = ŵᵀŵ/τ² + 1/ω²
        prec_c = ww / tau2 + 1.0 / omega2
        mean_c = (float(w_unit @ beta) / tau2) / prec_c
        c = mean_c + rng.standard_normal() / np.sqrt(prec_c)
        # τ² | ·   InvGamma over (β − c·ŵ)
        d = beta - c * w_unit
        tau2 = 1.0 / rng.gamma(a0 + p / 2.0, 1.0 / (b0 + 0.5 * float(d @ d)))
        # σ² | ·
        res = r - Xz @ beta
        sig2 = 1.0 / rng.gamma(a0 + n / 2.0, 1.0 / (b0 + 0.5 * float(res @ res)))
        if it >= burn:
            acc_b += beta; acc_t.append(tau2); acc_c.append(c)
    return acc_b / len(acc_t), float(np.mean(acc_t)), float(np.mean(acc_c))


def fit_gene_eb(Y, X, C, w_prior, *, n_iter=1500, burn=500, seed=0):
    """β (Series), plus learned (τ², c). β is on z-scored scale → returned M is raw-abundance scale."""
    Cmat = C.to_numpy(dtype=float)
    r = -_resid(Y.to_numpy(dtype=float), Cmat)
    Xr = _resid(X.to_numpy(dtype=float), Cmat)
    sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    w = np.clip(w_prior.reindex(X.columns).fillna(0.0).to_numpy(dtype=float), 0, None)
    nrm = np.linalg.norm(w)
    w_unit = w / nrm if nrm > 0 else w
    beta, tau2, c = _gibbs_eb(Xz, r, w_unit, n_iter=n_iter, burn=burn, seed=seed)
    M = np.where(sd > 1e-9, beta / (sd + 1e-9), 0.0)
    return pd.Series(M, index=X.columns, name="M"), tau2, c


def oof_gate_eb(gene, *, folds=5, seed=0, alpha=0.005, n_iter=1200, burn=400, w_prior_source="ledger"):
    """OOF coupling for the EB-shrink-toward-prior estimator alongside the fixed-α lasso (same folds/C)."""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source=w_prior_source)
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    n = len(Y)
    oof_eb = np.full(n, np.nan); oof_l = np.full(n, np.nan); oof_a = np.full(n, np.nan)
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for fi, (tr, te) in enumerate(kf.split(X)):
        M_eb, _, _ = fit_gene_eb(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w, seed=seed + fi, n_iter=n_iter, burn=burn)
        M_l = LR.fit_gene(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w, alpha=alpha)
        oof_eb[te] = LR.aggregate(X.iloc[te], M_eb)
        oof_l[te] = LR.aggregate(X.iloc[te], M_l)
        oof_a[te] = X.iloc[te].to_numpy(dtype=float).mean(axis=1)
    Cm = C.to_numpy(dtype=float)
    yr = _resid(Y.to_numpy(dtype=float), Cm)
    rho_eb = spearmanr(_resid(oof_eb, Cm), yr).correlation
    rho_l = spearmanr(_resid(oof_l, Cm), yr).correlation
    rho_a = spearmanr(_resid(oof_a, Cm), yr).correlation
    _, tau2, c = fit_gene_eb(Y, X, C, w, n_iter=n_iter, burn=burn)     # full-data learned anchor strength
    # τ²/c² ratio: how big is the data-driven dispersion relative to the prior-anchored part (>>1 ⇒ prior discounted)
    disp = tau2 / (c * c) if abs(c) > 1e-9 else float("inf")
    return {"gene": gene, "n": n, "n_pred": X.shape[1],
            "rho_abund": round(float(rho_a), 3), "rho_lasso": round(float(rho_l), 3),
            "rho_eb": round(float(rho_eb), 3), "eb_vs_lasso": bool(rho_eb <= rho_l + 1e-6),
            "tau2": round(tau2, 3), "c": round(c, 3), "disp_tau2/c2": round(disp, 2)}


def run(genes=None, *, folds=5):
    genes = genes or HUB_PANEL
    rows = []
    for g in genes:
        try:
            rows.append(oof_gate_eb(g, folds=folds))
        except Exception as e:
            rows.append({"gene": g, "error": repr(e)[:80]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 160):
        print(df.to_string(index=False))
    if "eb_vs_lasso" in df:
        print(f"\nEB-shrink-toward-prior: matches/beats LASSO {int(df['eb_vs_lasso'].sum())}/{df['eb_vs_lasso'].notna().sum()}"
              f" | mean rho_eb {df['rho_eb'].mean():.3f} vs rho_lasso {df['rho_lasso'].mean():.3f}"
              f" | mean learned τ² {df['tau2'].mean():.2f} (large τ² ⇒ data discounts the prior magnitude)")
    return df


if __name__ == "__main__":
    run()
