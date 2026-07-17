"""Occupancy link + abundance-threshold gauge + shared free-AGO pool (Design §Decision C/§G). STUB.

Replaces the linear a·M of regression.py at Phase 3. Key forms:
    occ(m,g,s) = a(m,s) / (a(m,s) + kappa(g,m)),  a = (X − x0)_+   # threshold-gated; x0 = Mullokandov floor
    phi(s)     = A_total(s) / (1 + sum occ)                         # shared pool → ceRNA + retires D(m)
    w_eff      = phi(s) · occ
The gauge is fixed by log2(RPM+1) abundance + the functional-copy threshold x0 (tuned out-of-fold).
"""
from __future__ import annotations
import numpy as np


def occupancy(a: np.ndarray, kappa: np.ndarray, x0: float = 0.0) -> np.ndarray:
    """Langmuir occupancy with a functional-threshold floor. a,kappa broadcastable."""
    ag = np.clip(a - x0, 0.0, None)
    return ag / (ag + kappa)


def shared_pool_scale(occ_sum: np.ndarray, a_total: np.ndarray) -> np.ndarray:
    """Per-sample free-AGO scaling φ (ceRNA / promiscuity emerges here; retires D(m))."""
    return a_total / (1.0 + occ_sum)


def transform(X, *, x0: float = 0.0, kappa=None):
    """Occupancy-transform the abundance matrix X (participant × arm), replacing linear a·w with the
    saturating occupancy a/(a+κ) above a functional-threshold floor x0 (Design §Decision C).

    - **Threshold gauge:** arms below x0 (the functional-copy floor, Mullokandov) → occ ≈ 0 (inert
      regardless of site) → they drop out. Breaks the a·w non-identifiability with a physical anchor.
    - **Saturation:** high-abundance arms saturate toward 1 → the aggregate stops being dominated by the
      single most-abundant arm; per-site affinity (the learned weight) regains leverage. This is what
      *retires D(m)* — promiscuous/abundant arms no longer swamp a gene by abundance alone.

    kappa: per-arm dissociation on the abundance scale. Default = per-arm median (so occ=0.5 at the arm's
    own median), a self-anchoring per-arm gauge.
    """
    import pandas as pd
    a = (X - x0).clip(lower=0.0)
    if kappa is None:
        kappa = a[a > 0].median().fillna(1.0)                  # per-arm median of the gated abundance
        kappa = kappa.replace(0.0, 1.0)
    occ = a / (a + kappa)
    return pd.DataFrame(occ, index=X.index, columns=X.columns)
