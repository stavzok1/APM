"""Small shared statistics helpers for the miRNA x Hallmark subproject.

Kept intentionally light; inferential tests follow the APM analysis-conduct
guardrails (report effect sizes, FDR-correct, confounder-adjust via
``analysis.utils.common.loaders.partial_spearman`` / ``partial_pearson``).
"""

from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np
import pandas as pd


def correlation_pair(
    y: pd.Series,
    x: pd.Series,
    covariates: pd.DataFrame | None = None,
    *,
    spearman_only: bool = False,
) -> Dict[str, object]:
    """Marginal + partial Spearman and Pearson for paired continuous vectors.

    Partial variants OLS-residualize on ``covariates`` then correlate residuals
    (same ladder as ``partial_spearman`` / ``partial_pearson`` in parent loaders).

    ``spearman_only=True`` skips the Pearson half (marginal + partial) for callers
    that never read it (e.g. cn_dosage_attribution) -- the Pearson keys stay NaN.
    Numerically identical for the Spearman outputs.
    """
    from scipy.stats import spearmanr
    from analysis.utils.common.loaders import residualize

    df = pd.concat([y.rename("y"), x.rename("x")], axis=1, join="inner").dropna()
    n = len(df)
    out: Dict[str, object] = {
        "n": n,
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "partial_rho": np.nan,
        "partial_p": np.nan,
        "partial_pearson_r": np.nan,
        "partial_pearson_p": np.nan,
    }
    if n < 15:
        return out

    sr, sp = spearmanr(df["x"], df["y"])
    out["spearman_rho"] = float(sr)
    out["spearman_p"] = float(sp)
    if not spearman_only:
        from scipy.stats import pearsonr
        pr, pp = pearsonr(df["x"], df["y"])
        out["pearson_r"] = float(pr)
        out["pearson_p"] = float(pp)

    if covariates is not None and len(covariates.columns):
        # OLS-residualize y and x on the covariates ONCE, then take both the Spearman
        # and Pearson of the residuals. Reproduces partial_spearman / partial_pearson
        # exactly (same ``residualize``, same dropna ladder, same <15 guards) but with
        # 2 OLS fits per call instead of 4 (each loader fn re-residualized internally).
        try:
            dfc = pd.concat(
                [df["y"].rename("y"), df["x"].rename("x"), covariates.reindex(df.index)],
                axis=1,
            ).dropna()
            if len(dfc) >= 15:
                cov_cols = dfc.drop(columns=["y", "x"])
                ry = residualize(dfc["y"], cov_cols)
                rx = residualize(dfc["x"], cov_cols)
                common = ry.index.intersection(rx.index)
                if len(common) >= 15:
                    ry_c, rx_c = ry.loc[common], rx.loc[common]
                    prho, ppv = spearmanr(ry_c, rx_c)
                    out["partial_rho"] = float(prho) if np.isfinite(prho) else np.nan
                    out["partial_p"] = float(ppv) if np.isfinite(ppv) else np.nan
                    if not spearman_only:
                        from scipy.stats import pearsonr
                        pr_r, pr_pv = pearsonr(ry_c, rx_c)
                        out["partial_pearson_r"] = float(pr_r) if np.isfinite(pr_r) else np.nan
                        out["partial_pearson_p"] = float(pr_pv) if np.isfinite(pr_pv) else np.nan
        except Exception:
            pass
    return out


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN-safe (NaNs pass through)."""
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    mask = ~np.isnan(p)
    m = int(mask.sum())
    if m == 0:
        return q
    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[mask] = out
    return q


def kruskal_across_strata(
    values: pd.Series,
    strata: pd.Series,
    *,
    min_per_group: int = 5,
) -> Tuple[float, float, int]:
    """Kruskal-Wallis across stratum levels for one numeric vector.

    Returns (H, p, n_groups). Requires >=2 groups each with >= min_per_group obs.
    """
    from scipy.stats import kruskal

    df = pd.DataFrame({"v": pd.to_numeric(values, errors="coerce"), "g": strata}).dropna()
    groups = [g["v"].values for _, g in df.groupby("g") if len(g) >= min_per_group]
    if len(groups) < 2:
        return np.nan, np.nan, len(groups)
    try:
        h, p = kruskal(*groups)
    except ValueError:
        return np.nan, np.nan, len(groups)
    return float(h), float(p), len(groups)


def hypergeom_enrichment(
    n_hit_in_set: int,
    n_set: int,
    n_hit_total: int,
    n_universe: int,
) -> Tuple[float, float]:
    """One-sided (enrichment) hypergeometric test.

    Probability of observing >= ``n_hit_in_set`` successes when drawing ``n_set``
    items from a universe of ``n_universe`` with ``n_hit_total`` successes.
    Returns (fold_enrichment, p_value).
    """
    from scipy.stats import hypergeom

    if n_set == 0 or n_universe == 0 or n_hit_total == 0:
        return np.nan, np.nan
    expected = n_set * n_hit_total / n_universe
    fold = (n_hit_in_set / expected) if expected > 0 else np.nan
    # P(X >= k) = sf(k-1)
    p = float(hypergeom.sf(n_hit_in_set - 1, n_universe, n_hit_total, n_set))
    return float(fold), p


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Per-row (per-gene/per-miRNA) z-score across columns (samples)."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0)
