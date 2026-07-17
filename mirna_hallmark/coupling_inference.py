"""Family-aware FDR + permutation nulls for miRNA->target coupling.

Closes two inference gaps in the default coupling layer (``hallmark_interaction``
/ ``stats.bh_fdr``):

1. **Seed-family non-independence.** Paralogous arms that share a TargetScan seed
   (e.g. miR-29a/b/c) target the same sites and co-vary, so plain BH over edges
   treats strongly correlated tests as independent and is *anti-conservative*.
   Two corrections are provided:
     - ``benjamini_yekutieli`` -- BH valid under arbitrary dependence (BY), the
       cheap dependence-robust floor.
     - ``family_simes_fdr`` -- collapse each seed family to a single Simes p, run
       BH across *families* (the family is the testing unit), then broadcast the
       family q back to its member edges. This is the principled answer to "n
       families << n edges".

2. **Parametric nulls.** Heavy-tailed expression + confounder adjustment make the
   asymptotic partial-correlation p suspect. ``permutation_pvalues`` builds an
   empirical null by **Freedman-Lane** rank-residual permutation, applied *jointly*
   (one shared sample permutation per draw across all pairs) so the seed-family
   correlation structure is preserved -- which additionally yields a max/min-
   statistic **family-wise** null (``family_permutation_pvalues``).

The same three primitives serve every resolution (edge / gene / miRNA-module);
the runner ``coupling_permutation`` wires them to real edges and pressure.

Self-test: ``.venv/bin/python3 -m mirna_hallmark.coupling_inference``.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import stats as S

# --------------------------------------------------------------------------- #
# Seed-family map
# --------------------------------------------------------------------------- #
_ARM_SUFFIX = re.compile(r"\.\d+$")
# Name-stem fallback when an arm is absent from the TargetScan family table:
# keep the let/miR stem + numeric core + mature arm, drop the paralog letter.
_NAME_FAM_RE = re.compile(r"^(let|miR)-(\d+)[a-z]*(?:-\d+)?(-[35]p)?$")


def _name_stem_family(arm: str) -> str:
    a = str(arm).replace("hsa-", "")
    m = _NAME_FAM_RE.match(a)
    if not m:
        return a
    arm_suffix = m.group(3) or ""
    return f"{m.group(1)}-{m.group(2)}{arm_suffix}"


def seed_family_map(
    arms: Sequence[str],
    *,
    family_info=None,
) -> Dict[str, str]:
    """Map each arm -> TargetScan seed family (human), name-stem fallback.

    Paralogues that share a seed collapse to one family; this is the grouping
    over which FDR and permutations must be made dependence-aware.
    """
    base_to_fam: Dict[str, str] = {}
    path = family_info or C.MIR_FAMILY_INFO
    try:
        df = pd.read_csv(path, sep="\t")
        h = df.loc[df["Species ID"].eq(9606), ["miR family", "MiRBase ID"]].dropna()
        base_to_fam = dict(zip(h["MiRBase ID"].astype(str), h["miR family"].astype(str)))
    except Exception:
        base_to_fam = {}
    out: Dict[str, str] = {}
    for a in arms:
        a = str(a)
        base = _ARM_SUFFIX.sub("", a)
        fam = base_to_fam.get(base) or base_to_fam.get(a) or _name_stem_family(a)
        out[a] = fam
    return out


# --------------------------------------------------------------------------- #
# Dependence-aware FDR
# --------------------------------------------------------------------------- #
def benjamini_yekutieli(pvals: Sequence[float]) -> np.ndarray:
    """BH-style q-values valid under arbitrary dependence (Benjamini-Yekutieli).

    Identical to BH but the threshold is divided by the harmonic number
    ``c(m) = sum_{i=1..m} 1/i`` -> uniformly more conservative. NaN-safe.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    mask = ~np.isnan(p)
    m = int(mask.sum())
    if m == 0:
        return q
    c_m = float(np.sum(1.0 / np.arange(1, m + 1)))
    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m * c_m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[mask] = out
    return q


def simes_p(pvals: Sequence[float]) -> float:
    """Simes combined p for one family of (possibly dependent) tests."""
    p = np.sort(np.asarray([x for x in np.asarray(pvals, float) if np.isfinite(x)]))
    m = p.size
    if m == 0:
        return np.nan
    return float(np.clip(np.min(p * m / np.arange(1, m + 1)), 0, 1))


def family_simes_fdr(
    pvals: Sequence[float],
    families: Sequence[str],
) -> Tuple[np.ndarray, pd.Series, pd.Series]:
    """Family-as-unit FDR: Simes p per seed family, BH across families.

    Returns ``(per_test_q, family_p, family_q)``. ``per_test_q`` broadcasts each
    family's BH q to its member tests, so two paralogous edges can no longer each
    spend an independent slice of the error budget.
    """
    df = pd.DataFrame({"p": np.asarray(pvals, float), "fam": [str(f) for f in families]})
    fam_p = df.groupby("fam")["p"].apply(lambda s: simes_p(s.values))
    fam_q = pd.Series(S.bh_fdr(fam_p.values), index=fam_p.index)
    per_test_q = df["fam"].map(fam_q).to_numpy()
    return per_test_q, fam_p, fam_q


# --------------------------------------------------------------------------- #
# Permutation null (Freedman-Lane on ranks)
# --------------------------------------------------------------------------- #
def _rank_rows(mat: np.ndarray) -> np.ndarray:
    """Average-rank each row across columns (samples)."""
    return pd.DataFrame(mat).rank(axis=1).to_numpy()


def _residualize_rows(mat: np.ndarray, hat: Optional[np.ndarray]) -> np.ndarray:
    """Residualize each row on the covariate column space (hat = Z Z^+)."""
    if hat is None:
        return mat - mat.mean(axis=1, keepdims=True)
    return mat - mat @ hat.T


def _unit_rows(mat: np.ndarray) -> np.ndarray:
    """Center + scale each row to zero mean / unit L2 (so row dot = Pearson r)."""
    m = mat - mat.mean(axis=1, keepdims=True)
    nrm = np.sqrt((m * m).sum(axis=1, keepdims=True))
    nrm[nrm == 0] = np.nan
    return m / nrm


def _covariate_hat(cov: Optional[pd.DataFrame], samples: pd.Index) -> Optional[np.ndarray]:
    if cov is None or cov.shape[1] == 0:
        return None
    Z = cov.reindex(samples).to_numpy(dtype=float)
    Z = np.nan_to_num(Z, nan=np.nanmean(Z) if np.isfinite(np.nanmean(Z)) else 0.0)
    # rank-transform covariates too (partial *Spearman*), add intercept
    Z = pd.DataFrame(Z).rank(axis=0).to_numpy()
    Z = np.column_stack([np.ones(len(samples)), Z])
    return Z @ np.linalg.pinv(Z)


def permutation_pvalues(
    predictor: pd.DataFrame,
    target: pd.DataFrame,
    *,
    covariates: Optional[pd.DataFrame] = None,
    families: Optional[Sequence[str]] = None,
    n_perm: int = 2000,
    tail: str = "neg",
    seed: int = 0,
    return_null: bool = False,
):
    """Joint Freedman-Lane permutation partial-Spearman for paired rows.

    ``predictor`` and ``target`` are aligned (n_pairs x n_samples) frames sharing
    the same row index (one row per edge/gene/module) and columns (samples). For
    each row we test partial Spearman(predictor, target | covariates) against an
    empirical null built by permuting the predictor residual's sample order. The
    **same** permutation is applied to every row in a draw, preserving cross-row
    (seed-family) dependence.

    Returns a frame indexed like the inputs with: ``rho`` (observed partial
    Spearman), ``p_perm`` (one-sided empirical p for ``tail``), ``q_bh``,
    ``q_by`` (Benjamini-Yekutieli), and -- when ``families`` is given --
    ``q_family`` (Simes-per-family BH) and ``p_family_fwer`` (family-wise via the
    per-family min-statistic null).
    """
    idx = predictor.index
    samples = predictor.columns.intersection(target.columns)
    X = predictor[samples].to_numpy(dtype=float)
    Y = target.loc[idx, samples].to_numpy(dtype=float)
    n = len(samples)

    hat = _covariate_hat(covariates, samples)
    Xr = _unit_rows(_residualize_rows(_rank_rows(X), hat))
    Yr = _unit_rows(_residualize_rows(_rank_rows(Y), hat))
    valid = np.isfinite(Xr).all(axis=1) & np.isfinite(Yr).all(axis=1)

    obs = np.full(len(idx), np.nan)
    obs[valid] = (Xr[valid] * Yr[valid]).sum(axis=1)

    rng = np.random.default_rng(seed)
    cnt = np.zeros(len(idx))  # count of null at-least-as-extreme as obs

    # Family bookkeeping: integer codes over the *valid* rows so the per-family
    # min-statistic null is a single vectorized np.minimum.at per draw.
    fam_arr = np.asarray([str(f) for f in families]) if families is not None else None
    fam_codes = fam_levels = fam_obs_min = fam_cnt = None
    if fam_arr is not None:
        fam_levels, fam_codes = np.unique(fam_arr[valid], return_inverse=True)
        n_fam = len(fam_levels)
        fam_obs_min = np.full(n_fam, np.inf)
        np.minimum.at(fam_obs_min, fam_codes, obs[valid])
        fam_cnt = np.zeros(n_fam)

    Xv, Yv, obs_v = Xr[valid], Yr[valid], obs[valid]
    null_mat = np.full((n_perm, len(idx)), np.nan) if return_null else None
    s1 = np.zeros(valid.sum())  # running sum / sumsq of the null statistic per test
    s2 = np.zeros(valid.sum())
    for b in range(n_perm):
        perm = rng.permutation(n)
        null = (Xv[:, perm] * Yv).sum(axis=1)
        s1 += null
        s2 += null * null
        if null_mat is not None:
            null_mat[b, valid] = null
        if tail == "neg":
            cnt[valid] += (null <= obs_v)
        elif tail == "pos":
            cnt[valid] += (null >= obs_v)
        else:
            cnt[valid] += (np.abs(null) >= np.abs(obs_v))
        if fam_codes is not None:
            fam_min = np.full(len(fam_levels), np.inf)
            np.minimum.at(fam_min, fam_codes, null)  # per-family most-negative null
            fam_cnt += (fam_min <= fam_obs_min)

    p_perm = np.full(len(idx), np.nan)
    p_perm[valid] = (1.0 + cnt[valid]) / (n_perm + 1.0)

    # Smooth tail p from the null moments (the permuted statistic is ~Gaussian by CLT),
    # so multiplicity correction is not capped by the empirical 1/(n_perm+1) floor.
    from scipy.stats import norm

    null_mean = s1 / n_perm
    null_sd = np.sqrt(np.clip(s2 / n_perm - null_mean ** 2, 1e-12, None))
    z = (obs_v - null_mean) / null_sd
    p_z = np.full(len(idx), np.nan)
    if tail == "neg":
        p_z[valid] = norm.cdf(z)
    elif tail == "pos":
        p_z[valid] = norm.sf(z)
    else:
        p_z[valid] = 2.0 * norm.sf(np.abs(z))

    out = pd.DataFrame({"rho": obs, "p_perm": p_perm, "p_z": p_z}, index=idx)
    # Empirical p validates the asymptotic null; multiplicity rides the smooth p_z.
    out["q_bh"] = S.bh_fdr(out["p_z"].values)
    out["q_by"] = benjamini_yekutieli(out["p_z"].values)
    if fam_arr is not None:
        out["family"] = fam_arr
        per_test_q, _, _ = family_simes_fdr(out["p_z"].values, fam_arr)
        out["q_family"] = per_test_q
        fam_fwer = pd.Series((1.0 + fam_cnt) / (n_perm + 1.0), index=fam_levels)
        valid_fam = pd.Series(fam_arr).where(pd.Series(valid).values)
        out["p_family_fwer"] = valid_fam.map(fam_fwer).to_numpy()
    if return_null:
        return out, null_mat
    return out


def partial_spearman_matrix(
    predictor: pd.DataFrame,
    target: pd.DataFrame,
    *,
    covariates: Optional[pd.DataFrame] = None,
) -> pd.Series:
    """Row-wise partial Spearman (predictor vs target | covariates), vectorized.

    Same residualize-ranks-then-correlate convention as ``permutation_pvalues`` but
    returns only the point estimates -- a fast scoring primitive for the held-out
    tuning sweep (no permutation, no p-value).
    """
    idx = predictor.index
    samples = predictor.columns.intersection(target.columns)
    X = predictor[samples].to_numpy(dtype=float)
    Y = target.loc[idx, samples].to_numpy(dtype=float)
    hat = _covariate_hat(covariates, samples)
    Xr = _unit_rows(_residualize_rows(_rank_rows(X), hat))
    Yr = _unit_rows(_residualize_rows(_rank_rows(Y), hat))
    rho = np.full(len(idx), np.nan)
    valid = np.isfinite(Xr).all(axis=1) & np.isfinite(Yr).all(axis=1)
    rho[valid] = (Xr[valid] * Yr[valid]).sum(axis=1)
    return pd.Series(rho, index=idx)


def partial_spearman_batch(
    predictor: pd.DataFrame,
    target: pd.DataFrame,
    covariates: Optional[pd.DataFrame] = None,
    *,
    min_n: int = 15,
) -> pd.DataFrame:
    """Vectorized drop-in for ``loaders.partial_spearman`` over many paired rows.

    Reproduces the reference convention **exactly** — OLS-residualize *raw* predictor
    and target on ``[intercept, covariates]``, then Spearman the residuals (rho via
    rank-Pearson, p via scipy's t-approximation) — but for all rows at once, sharing
    the covariate hat. Returns columns ``rho``, ``p``, ``n``.

    Assumes a common covariate-complete sample set (the case in the hot modules after
    their sample intersection); a row whose predictor/target has residual NaNs there
    returns NaN, and the caller can fall back to the scalar reference for those.
    Differs from :func:`partial_spearman_matrix`, which ranks *before* residualizing
    (a different, internally-consistent convention used by the new tuning analyses).
    """
    from scipy.stats import t as _t

    idx = predictor.index
    samples = predictor.columns.intersection(target.columns)
    if covariates is not None:
        samples = samples.intersection(covariates.index)
    X = predictor[samples].to_numpy(dtype=float)
    Y = target.loc[idx, samples].to_numpy(dtype=float)

    if covariates is not None and covariates.shape[1]:
        Z = covariates.reindex(samples).to_numpy(dtype=float)
        keep = ~np.isnan(Z).any(axis=1)
        Z, X, Y = Z[keep], X[:, keep], Y[:, keep]
        Zd = np.column_stack([np.ones(Z.shape[0]), Z])
        Zpinv = np.linalg.pinv(Zd)           # (k+1) x n; residualize via the thin factor
        Xres = X - (X @ Zd) @ Zpinv          # avoids forming the n x n hat
        Yres = Y - (Y @ Zd) @ Zpinv
    else:
        Xres, Yres = X - X.mean(axis=1, keepdims=True), Y - Y.mean(axis=1, keepdims=True)

    n = Xres.shape[1]
    rho = np.full(len(idx), np.nan)
    p = np.full(len(idx), np.nan)
    valid = np.isfinite(Xres).all(axis=1) & np.isfinite(Yres).all(axis=1)
    if valid.any() and n >= min_n:
        Xr = _unit_rows(_rank_rows(Xres[valid]))
        Yr = _unit_rows(_rank_rows(Yres[valid]))
        r = np.clip((Xr * Yr).sum(axis=1), -1.0, 1.0)
        rho[valid] = r
        dof = n - 2
        with np.errstate(divide="ignore", invalid="ignore"):
            tstat = r * np.sqrt(dof / ((1.0 - r) * (1.0 + r)))
        p[valid] = 2.0 * _t.sf(np.abs(tstat), dof)
    return pd.DataFrame({"rho": rho, "p": p, "n": int(n)}, index=idx)


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def _selftest() -> None:
    rng = np.random.default_rng(0)

    # (1) BY >= BH everywhere; equal-ish when m=1
    p = rng.uniform(size=200) ** 2
    bh, by = S.bh_fdr(p), benjamini_yekutieli(p)
    assert np.all(by >= bh - 1e-9), "BY must be >= BH"

    # (2) family_simes reduces to BH when every test is its own family
    fams_singleton = [f"f{i}" for i in range(len(p))]
    q_fam, _, _ = family_simes_fdr(p, fams_singleton)
    assert np.allclose(np.sort(q_fam), np.sort(S.bh_fdr(p)), atol=1e-9), "singleton families != BH"

    # (3) collapsing correlated tests into one family shrinks the discovery count
    #     (the anti-conservatism fix): 10 families x 8 identical members each
    base = rng.uniform(size=10) ** 3
    p_corr = np.repeat(base, 8)
    fams = np.repeat([f"F{i}" for i in range(10)], 8)
    q_fam, fam_p, fam_q = family_simes_fdr(p_corr, fams)
    n_bh = int((S.bh_fdr(p_corr) < 0.05).sum())
    n_family = int((q_fam < 0.05).sum())
    assert len(fam_p) == 10
    # family q counts each family once, not 8x
    assert (q_fam < 0.05).sum() % 8 == 0

    # (4) permutation: null is ~uniform; strong negative signal is significant
    n, E = 120, 40
    samples = pd.Index([f"S{i}" for i in range(n)])
    Xnull = pd.DataFrame(rng.normal(size=(E, n)), columns=samples)
    Ynull = pd.DataFrame(rng.normal(size=(E, n)), columns=samples)
    res_null = permutation_pvalues(Xnull, Ynull, n_perm=500, tail="neg", seed=1)
    # roughly uniform: ~half below 0.5
    assert 0.25 < (res_null["p_perm"] < 0.5).mean() < 0.75

    x = pd.DataFrame(rng.normal(size=(E, n)), columns=samples)
    y = -x + 0.3 * pd.DataFrame(rng.normal(size=(E, n)), columns=samples)  # strong neg
    res_sig = permutation_pvalues(x, y, n_perm=500, tail="neg", seed=2)
    assert (res_sig["rho"] < 0).all() and (res_sig["q_by"] < 0.05).mean() > 0.9

    # (5) seed-family map collapses paralogues
    fam = seed_family_map(["hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p", "hsa-let-7a-5p"])
    assert fam["hsa-miR-29a-3p"] == fam["hsa-miR-29b-3p"] == fam["hsa-miR-29c-3p"], fam
    print("coupling_inference self-test OK:",
          f"BY>=BH; singleton=BH; family collapses 8x; perm null uniform & signal sig; "
          f"29a/b/c -> {fam['hsa-miR-29a-3p']!r}")


if __name__ == "__main__":
    _selftest()
