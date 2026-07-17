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
    """Fit M (Series[arm], >= 0) for one gene. Y,X,C aligned on the same participants."""
    Cmat = C.to_numpy(dtype=float)
    yr = _residualize(Y.to_numpy(dtype=float), Cmat)
    Xr = _residualize(X.to_numpy(dtype=float), Cmat)            # multi-output OLS residuals

    w = w_prior.reindex(X.columns).to_numpy(dtype=float)
    Xs = Xr * w                                                 # adaptive: scale columns by prior weight
    # Y_resid ≈ − X_resid·M ; fit (−yr) on Xs with M>=0, then unscale.
    lasso = Lasso(alpha=alpha, positive=True, max_iter=max_iter).fit(Xs, -yr)
    M = pd.Series(lasso.coef_ * w, index=X.columns, name="M")
    return M


def aggregate(X: pd.DataFrame, M: pd.Series) -> np.ndarray:
    """Predicted repression pressure ρ = X·M (higher ρ → more repression → lower Y)."""
    return X.to_numpy(dtype=float) @ M.reindex(X.columns).fillna(0.0).to_numpy(dtype=float)
