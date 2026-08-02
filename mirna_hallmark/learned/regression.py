"""Gene-focused non-negative adaptive-lasso (Design §3a Phase-1 estimator, §Decision A alternative).

Model (one gene g):  Y = alpha + C·beta − X·M + eps ,  M_{g,m} >= 0.
- Confounders C are partialled out of Y and X first (the matrix-view identification, Design §2/§B):
  M is estimated in the residual space (I − P_C).
- Non-negativity encodes "miRNAs repress" (Design §Decision C sign/box).
- **Adaptive** penalty λ_{g,m} = λ / w_prior(g,m) via the standard feature-scaling reparametrisation
  (scale column m by w_m, plain lasso, unscale): strong curated edges shrink little, weak edges heavily.
  This is where the evidence prior enters as MAGNITUDE, not just support (Design §3b channel 1–2).

Unsaturated (linear) occupancy: contribution ≈ a·M. The saturating occupancy link + threshold gauge +
shared pool (Design §C/§G) replace this in occupancy.py at Phase 3.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression


def _residualize(V: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    """Return V minus its OLS fit on C (intercept included by LinearRegression)."""
    lr = LinearRegression().fit(Cmat, V)
    return V - lr.predict(Cmat)


def fit_gene(
    Y: pd.Series,
    X: pd.DataFrame,
    C: pd.DataFrame,
    w_prior: pd.Series,
    *,
    alpha: float = 0.01,
    max_iter: int = 50_000,
) -> pd.Series:
    """⛔ **RETIRED (MH-115) — the ADAPTIVE LASSO. Use `fit_gene_bayes` below for any ESTIMAND.**

    Fit M (Series[arm], >= 0) for one gene. Y,X,C aligned on the same participants.

    ⚠ This is still called in ~30 places and that is NOT all wrong — see `eval/_retired_lasso_audit.py`,
    which classifies every call site. It is legitimate as a **BASELINE** (the lasso IS the comparator in
    `bayes_parity`, `shuffled_null`, `spike_slab`, `attribution_eb`'s support test) and measured **INERT**
    as a **NUISANCE** (`he_agg` in `discovery`/`dossier`: swapping to Gibbs moved coupling by −0.0012,
    p=0.255, n=183). It is a **DEFECT only where the returned M is the reported quantity** — those sites
    were switched to `fit_gene_bayes` on 2026-08-01 (MH-184).
    ⛔ **Do not delete this function and do not bulk-replace its callers** — MH-142 records what happens
    when a load-bearing term is "fixed" without first checking its role."""
    Cmat = C.to_numpy(dtype=float)
    yr = _residualize(Y.to_numpy(dtype=float), Cmat)
    Xr = _residualize(X.to_numpy(dtype=float), Cmat)            # multi-output OLS residuals

    w = w_prior.reindex(X.columns).to_numpy(dtype=float)
    Xs = Xr * w                                                 # adaptive: scale columns by prior weight
    # Y_resid ≈ − X_resid·M ; fit (−yr) on Xs with M>=0, then unscale.
    lasso = Lasso(alpha=alpha, positive=True, max_iter=max_iter).fit(Xs, -yr)
    M = pd.Series(lasso.coef_ * w, index=X.columns, name="M")
    return M


def fit_gene_bayes(Y: pd.Series, X: pd.DataFrame, C: pd.DataFrame, w_prior: pd.Series, *,
                   n_iter: int = 1200, burn: int = 400, seed: int = 0, **_) -> pd.Series:
    """⭐ THE CANONICAL DROP-IN for `fit_gene` — learned-τ² dense Gibbs posterior mean, same signature/scale.

    Promoted here from `learned/eval/bayes_parity.py` on 2026-08-01 (MH-184): the replacement for the
    retired adaptive lasso lived inside an EVAL module, so every caller that wanted to stop using the
    lasso had to import from a test harness. It belongs next to the thing it replaces.

    ⚠ **SCALE CONVENTION — this matches `fit_gene`, NOT the canonical readouts path.** It residualises on
    C and scales columns by the prior weight `w` (exactly as the adaptive lasso does) so it is a genuine
    drop-in. `learned/readouts.py` instead z-scores the SEED-FAMILY columns and calls `_gibbs_posterior`
    with `pi = ones(p)`. ⇒ **swapping a call site here removes the RETIRED ESTIMATOR; it does not make the
    result equal to the shipped `beta` column.** Do not conflate the two (axiom 6).
    ⚠ Any persisted artifact produced by a switched call site is STALE until re-run.
    """
    from mirna_hallmark.learned import spike_slab as SS
    Cmat = C.to_numpy(float)
    yr = _residualize(Y.to_numpy(dtype=float), Cmat)
    Xr = _residualize(X.to_numpy(dtype=float), Cmat)
    w = np.nan_to_num(w_prior.reindex(X.columns).to_numpy(float), nan=0.0)
    w = np.clip(w, 1e-3, None)
    coef, _, _ = SS._gibbs_posterior(Xr * w, -yr, np.ones(Xr.shape[1]),
                                     n_iter=n_iter, burn=burn, seed=seed)
    return pd.Series(coef * w, index=X.columns, name="M")


def aggregate(X: pd.DataFrame, M: pd.Series) -> np.ndarray:
    """Predicted repression pressure ρ = X·M (higher ρ → more repression → lower Y)."""
    return X.to_numpy(dtype=float) @ M.reindex(X.columns).fillna(0.0).to_numpy(dtype=float)
