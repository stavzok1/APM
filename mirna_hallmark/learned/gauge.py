"""The CROSS-COHORT GAUGE — shared by every channel whose evidence comes from a different cohort.

⚠⚠ **READ THE VERDICT FIRST (MH-102d, 2026-07-12): every cross-cohort channel tested here contributes
0.4–2.4% of the mRNA likelihood's own precision on β_T, and carries NO detectable per-edge payload.
The state channel was CANCELLED on this evidence. This module exists so the next channel is measured
BEFORE it is built — run `calibrate()` + `fit_gauge()` and read `info_ratio` and `tau`.**

The problem
-----------
A channel whose evidence lives in another cohort observes β on a different scale: bulk composition
dilutes it differently, the platform differs, n differs (⇒ the estimator shrinks differently), and the C
blocks cannot be column-identical (`CPE`/`target_cn` don't exist in healthy tissue). Feeding a raw
`β_source` in as `pihat` imports a systematic offset **that reads as biology** — under the naive additive
`Δ = β_T − β_H`, attenuation alone makes Δ>0 on nearly every edge, manufacturing "acquired wiring".

The model
---------
    β_source,e = a · β_target,e + c + δ_e ,     δ_e ~ N(0, τ²)

`a` = ONE global scale (nuisance: composition + platform + C differences). `c` = the NNLS non-negativity
floor. **δ is the payload** (acquired/lost wiring; protein-specific realization) — the global scale is a
nuisance we never claim. Channel term (Gaussian-conjugate ⇒ stays **J1/Gibbs**):

    pihat = (β̂_src − c)/a      s² = (τ² + max(se, se_min)²)/a²      members = [(col, 1.0)]

TWO CONVENTIONS THAT ARE NOT OPTIONAL
-------------------------------------
1. **E3 (`zscore_y=True`, the default)** — z-score the C-residualised Y within each cohort so β is a
   STANDARDIZED effect, comparable across cohorts by construction. Without it the gauge silently absorbs
   the cohort's Y-scale: median sd(resid Y) is **TCGA 0.237 vs NAT 0.600 / GTEx 0.612** (the normal tissues
   are ~2.5× more variable), and that entire factor was being reported as "composition attenuation"
   (a_raw/a_zY = 2.25 ≈ the sd ratio 2.58). CHANNEL_FUSION §E.
2. **`calibrate()` — the split-half control, MANDATORY.** `Gauge.usable` (and therefore
   `to_channel_terms`) refuses to emit anything until it passes.
   ⚠ **THIS IS ABOUT THE GAUGE, NOT ABOUT THE MODEL'S ESTIMATOR — an earlier framing here was wrong and is
   corrected.** Fed the **dense Gibbs posterior's** β/SD, this gauge returns **a = 4.1** (3.5e6 in one
   configuration) where the truth is 1.0. **Gibbs is NOT at fault — it is the BETTER estimator** (split-half
   reproducibility ρ=**0.822** vs bagged NNLS's 0.729). The fault is the errors-in-variables correction:
   it divides by `Var(β̂) − mean(se²)`, and for Gibbs `mean(se²)` is dominated by a **heavy tail of a few
   enormous posterior SDs** — `sqrt(mean(se²))=0.055` against a *typical* se of 0.015 — which collapses the
   denominator (reliability 0.17 vs NNLS's 0.72). `_MIN_RELIABILITY` now REFUSES that gauge instead of
   silently returning 4.1. **Use bagged NNLS for the GAUGE; keep Gibbs for the MODEL.**

Estimating `a`
--------------
Errors-in-variables ⇒ OLS is attenuated. `_eiv_slope` uses the attenuation-corrected (regression-
calibration) slope and returns the target's **reliability** = 1 − E[se²]/Var(β̂); below `_MIN_RELIABILITY`
the correction is not usable and the gauge is refused. Fit on **ALL** edges — never the both-nonzero
subset (that selects for agreement and biases `a` up; the β_T>0 / β_H=0 edges ARE the attenuation).
`a` is applied **out-of-fold over genes** (`a_for`), so the channel never uses a scale calibrated on its
own gene.

Falsification
-------------
`shuffled_gauge()` — the shared CROSS-COHORT SHUFFLE-NULL primitive (analogue of
`instrument.exclusion(shuffle_cn=)`): permute the source cohort's sample labels, refit β, refit the gauge;
**`a` must collapse**. It has already earned its keep twice — it caught a through-the-origin slope that
read the NNLS positivity floor as signal, and it currently **FAILS for NAT** (a_null = 0.076, CI excludes
0), which retires NAT's apparent per-edge payload.

CLI: `python -m mirna_hallmark.learned.gauge [--source gtex] [--n-genes 150]`
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned.attribution_eb import _prep, _bagged_nnls_meansd

_CACHE: dict = {}


# ============================ cohort matrices (the shared β producer) ============================
def cohort_matrices(cohort: str):
    """(X: arm×participant miRNA, Y: gene×participant mRNA) for a cohort, on ONE participant key."""
    key = f"mat_{cohort}"
    if key in _CACHE:
        return _CACHE[key]
    if cohort == "tcga":
        from mirna_hallmark.learned import data as LD
        d = LD._load()
        X, Y = d["X"], d["Y"]
    elif cohort == "nat":
        from mirna_hallmark.learned import states as STS
        X, Y = STS.state_matrices("11")
    elif cohort == "gtex":
        from mirna_hallmark.learned import state as ST
        X, Y = ST._gtex_mirna(), ST._gtex_mrna()
        X.columns = CF._to_participant("gtex", X.columns)
        Y.columns = CF._to_participant("gtex", Y.columns)
        X, Y = X.T.groupby(level=0).mean().T, Y.T.groupby(level=0).mean().T
    elif cohort == "cptac":
        from mirna_hallmark.eval import cptac_validation as CV
        X = CV.load_prospective_mirna_arms()
        Y = CF.expression("cptac")
    else:
        raise ValueError(cohort)
    X = X[~X.index.duplicated(keep="first")]
    Y = Y[~Y.index.duplicated(keep="first")]
    _CACHE[key] = (X, Y)
    return X, Y


def beta_table(cohort: str, genes: Optional[Sequence[str]] = None, *, n_boot: int = 30,
               min_n: int = 40, min_sd_rel: float = 0.25, permute: bool = False,
               zscore_y: bool = True, seed: int = 0, **c_kw) -> pd.DataFrame:
    """Per-(gene, family) β for a cohort, using `build_C` + the canonical bagged-NNLS estimator.

    **zscore_y=True is the E3 GAUGE CONVENTION (CHANNEL_FUSION §E) and the DEFAULT — do not turn it off for
    a cross-cohort gauge.** It z-scores the C-residualised Y within each cohort, so β is a STANDARDIZED
    effect ("SD of Y per SD of dose") and is comparable across cohorts BY CONSTRUCTION. Without it, β is on
    the raw-`r` scale and the gauge silently absorbs the cohort's Y-scale: median sd(residual Y) is
    **TCGA 0.237 vs NAT 0.600 / GTEx 0.612** — the normal tissues are ~2.5× more variable — and that whole
    factor was being read as "composition attenuation" (measured: a_raw/a_zY = 2.25 ≈ the sd ratio 2.58).
    ⚠ The model's own posterior runs on raw-`r`; this convention is for the GAUGE only.

    permute=True → SHUFFLE the miRNA sample labels within this cohort (breaking the miRNA↔mRNA pairing
    while preserving both marginals) = the cross-cohort null. β must collapse.

    ⚠ **`min_sd_rel` replaces an ABSOLUTE `min_sd=0.2`, which was a SCALE-DEPENDENT BUG** (found by the CPTAC
    session, 2026-07-12). Cohort median sd(Y): **TCGA 0.233 · NAT 0.595 · GTEx 0.612 · CPTAC 1.058** — the
    cohorts are on different scales, so a fixed 0.2 dropped **34% of TCGA genes but 0% of CPTAC**. Because
    `fit_gauge` merges on the SHARED (gene, family) keys, the gauge was therefore fit on "TCGA's high-variance
    third" — a biased subset selected by an absolute threshold on a cohort-specific scale. The filter is now
    RELATIVE to each cohort's own median gene sd (`sd < min_sd_rel · median_sd`), so it means the same thing
    everywhere: "this gene is flat *for this cohort*".
    Returns columns: gene · family · beta · se · n.
    """
    from mirna_hallmark.learned.evidence import ledger as LG
    X, Y = cohort_matrices(cohort)
    ed = LG.pooled_he_edges()
    genes = list(genes) if genes is not None else [g for g in ed["gene"].unique() if g in Y.index]

    dec = CF.deconv(cohort)
    parts = [p for p in Y.columns if p in X.columns and (dec is None or p in dec.index)]
    if len(parts) < min_n:
        raise ValueError(f"{cohort}: only {len(parts)} usable participants")
    C = CF.build_C(cohort, parts, **c_kw)
    Cd = np.c_[np.ones(len(parts)), C.to_numpy(float)] if C.shape[1] else np.ones((len(parts), 1))
    Cdf = pd.DataFrame(Cd, index=parts)

    from mirna_hallmark.learned import data as _LD
    _nm = _LD._arm_name_map(X)                           # built on the FULL cohort matrix (guide tiebreak = full abundance)

    Xp = X[parts]
    if permute:                                          # NULL: break the miRNA↔mRNA pairing
        rng = np.random.default_rng(seed)
        Xp = Xp.iloc[:, rng.permutation(len(parts))]
        Xp.columns = parts

    # SCALE-FREE variability floor on Y: relative to THIS cohort's own median gene sd (see the docstring).
    _sds = np.array([float(Y.loc[g, parts].astype(float).std()) for g in genes if g in Y.index])
    _sds = _sds[np.isfinite(_sds) & (_sds > 0)]
    min_sd = float(min_sd_rel * np.median(_sds)) if len(_sds) else 0.0
    # ...and the SAME on X. `attribution_eb._prep` zeroes any family with an ABSOLUTE sd < 0.1, but the
    # cohorts' miRNA dose scales differ badly (median sd(X_fam): TCGA 0.947 · NAT 0.682 · **GTEx 0.235** ·
    # CPTAC 0.861), so that floor zeroes **32% of GTEx's families vs 3% of TCGA's** — forcing β_H to zero
    # for a third of the source and biasing `a` DOWN. Same bug class as the Y floor, on the predictor side.
    # `_prep` itself is NOT changed (it is the model's TCGA-only path; touching it ripples everywhere).
    x_floor = float(0.1 * _x_scale(cohort, Xp, parts, genes, ed) / 0.947)

    rows = []
    for g in genes:
        if g not in Y.index:
            continue
        y = Y.loc[g, parts].astype(float)
        keep = y.notna()                                 # CPTAC's proteogenomic matrices are sparse —
        if int(keep.sum()) < min_n:                      # drop only the missing participants, per gene
            continue
        pk = [p for p, k in zip(parts, keep) if k]
        y = y.loc[pk]
        if y.std() < min_sd:                             # gene not variable here → β undefined, skip
            continue
        # edge-arm -> matrix-row resolution via `_arm_name_map` (the `.N`/case/suffixless-guide fix). A bare
        # `m in Xp.index` test SILENTLY DROPS `.N`-suffixed arms, and does so ASYMMETRICALLY across cohorts
        # (measured: TCGA loses 19 edge arms, CPTAC 58) — which would make a cross-cohort gauge compare
        # DIFFERENT arm sets. Pull by matrix name, relabel to the canonical edge name (as `assemble_gene` does).
        regs = [m for m in ed.loc[ed["gene"] == g, "miRNA"].unique()
                if (_nm.get(m) or _nm.get(str(m).lower()))]
        if len(regs) < 2:
            continue
        regs_x = [_nm.get(m) or _nm.get(str(m).lower()) or m for m in regs]
        Xg = Xp.loc[regs_x, pk].T.astype(float).fillna(0.0)
        Xg.columns = list(regs)
        Xf, _, _ = FAM.collapse_by_family(Xg, pd.Series(1.0, index=regs), FAM.family_of(pd.Index(regs)))
        if Xf.shape[1] < 1:
            continue
        yr, Xz, cols = _prep_scalefree(y, Xf, Cdf.loc[pk], x_floor)
        if zscore_y:                                     # E3: standardized, cross-cohort-comparable β
            s = float(np.std(yr))
            if s < 1e-9:
                continue
            yr = yr / s
        b, sd = _bagged_nnls_meansd(Xz, yr, n_boot=n_boot, seed=seed)
        for c, v, s in zip(cols, b, sd):
            rows.append({"gene": g, "family": c, "beta": float(v), "se": float(s), "n": len(pk)})
    return pd.DataFrame(rows)


def _x_scale(cohort, Xp, parts, genes, ed) -> float:
    """Median family-dose sd for this cohort — the scale the X floor must be relative to."""
    key = f"xscale_{cohort}"
    if key in _CACHE:
        return _CACHE[key]
    sds = []
    for g in list(genes)[:120]:
        regs = [m for m in ed.loc[ed["gene"] == g, "miRNA"].unique() if m in Xp.index]
        if len(regs) < 2:
            continue
        Xg = Xp.loc[regs, parts].T.astype(float).fillna(0.0)
        Xf, _, _ = FAM.collapse_by_family(Xg, pd.Series(1.0, index=regs), FAM.family_of(pd.Index(regs)))
        sds.extend(Xf.std(ddof=0).to_numpy())
    sds = [s for s in sds if np.isfinite(s) and s > 0]
    _CACHE[key] = float(np.median(sds)) if sds else 1.0
    return _CACHE[key]


def _prep_scalefree(y, Xf, Cd, x_floor: float):
    """`attribution_eb._prep`, but with the family-dose floor made RELATIVE to the cohort's own scale."""
    from sklearn.linear_model import LinearRegression
    Cm = Cd.to_numpy(float)
    yv = y.to_numpy(float)
    yr = -(yv - LinearRegression().fit(Cm, yv).predict(Cm))
    sd = Xf.std(ddof=0)
    Xz = ((Xf - Xf.mean()) / (sd + 1e-9)).fillna(0.0)
    Xz.loc[:, sd < x_floor] = 0.0                        # scale-free: same MEANING in every cohort
    return yr, Xz.to_numpy(float), list(Xf.columns)


# ============================ the gauge ============================
@dataclass
class Gauge:
    a: float                      # the global scale  β_source ≈ a · β_target + c
    a_lo: float
    a_hi: float                   # bootstrap 95% CI over edges
    c: float                      # intercept — absorbs the NNLS non-negativity floor
    tau: float                    # per-edge deviation SD (the PAYLOAD's scale)
    lam: float                    # error-variance ratio se_src²/se_tgt²
    n: int
    rho: float                    # Spearman(source, target) — diagnostic, not used in the fit
    reliability: float = 1.0      # 1 − E[se_tgt²]/Var(β̂_tgt); < _MIN_RELIABILITY ⇒ `a` is not estimable
    info_ratio: float = 0.0       # the channel's precision on β_T as a FRACTION of the mRNA likelihood's own
    calibrated: Optional[bool] = None   # did the split-half control pass? (None = not run)
    source: str = ""
    target: str = ""
    fold_a: dict = field(default_factory=dict)     # out-of-fold (a, c), keyed by fold id
    gene_fold: dict = field(default_factory=dict)  # gene -> fold id

    def __repr__(self):
        cal = {True: "cal✓", False: "CAL-FAIL", None: "uncal"}[self.calibrated]
        return (f"Gauge({self.source}→{self.target}: a={self.a:.3f} [{self.a_lo:.3f},{self.a_hi:.3f}] "
                f"c={self.c:+.4f} τ={self.tau:.4f} ρ={self.rho:.3f} rel={self.reliability:.2f} "
                f"info={self.info_ratio:.1%} n={self.n} {cal})")

    @property
    def dead(self) -> bool:
        """The channel carries no information about β_target (a=0 is not excluded by the CI)."""
        return self.a_lo <= 0.0

    @property
    def usable(self) -> bool:
        """Safe to emit channel terms from? Requires a live `a`, an estimable slope, and a PASSED
        split-half calibration. `to_channel_terms` refuses otherwise — silently returning a mis-scaled
        `pihat` is how this module nearly shipped a=4.1 as if it were 1.0."""
        return (not self.dead) and self.reliability >= _MIN_RELIABILITY and self.calibrated is True


#: Minimum RELIABILITY of the target β for the errors-in-variables correction to be usable.
#: reliability = 1 − E[se_x²]/Var(x)  — the fraction of β̂_target's variance that is real signal.
#: Below this the correction divides by ~0 and `a` explodes. THIS IS NOT HYPOTHETICAL: with the dense
#: Gibbs posterior's SDs (which are large — they carry prior uncertainty) the denominator collapsed and the
#: gauge returned a = 3,481,707 on one configuration, and **failed the split-half control at a = 4.1 when
#: the true answer is 1.0**. Guard it, and refuse to return a gauge that cannot be calibrated (`calibrate`).
_MIN_RELIABILITY = 0.25


def _eiv_slope(x, y, se_x) -> tuple[float, float, float]:
    """Attenuation-corrected (errors-in-variables) slope + intercept of `y ~ a·x + c`, and the target's
    RELIABILITY (the diagnostic that says whether the correction is even usable).

        a = Cov(x, y) / ( Var(x) − E[se_x²] )        c = mean(y) − a·mean(x)
        reliability = 1 − E[se_x²] / Var(x)          ∈ (0,1];  → 0 ⇒ the denominator vanishes ⇒ `a` blows up

    The denominator is x's TRUE variance (observed minus measurement-error variance), so `a` is not
    attenuated by the noise in β̂_target — the classic regression-calibration correction.

    ⚠ NOT through-the-origin TLS (first implementation, REJECTED): bagged **NNLS β is ≥0 with a POSITIVE
    MEAN even on pure noise**, so an origin-forced slope reads `E[x]·E[y] > 0` as signal — the SHUFFLED NULL
    DID NOT COLLAPSE (a_null=0.254 vs a_real=0.504, while the shuffled Spearman was correctly dead at
    0.084). Centering absorbs the NNLS positivity floor into the intercept `c`, and `a → 0` exactly when the
    covariance does.
    """
    vobs = float(np.var(x))
    verr = float(np.mean(se_x ** 2))
    rel = 1.0 - verr / max(vobs, 1e-12)                         # reliability of β̂_target
    vx = max(vobs - verr, 1e-9)
    a = float(np.cov(x, y, bias=True)[0, 1] / vx)
    c = float(np.mean(y) - a * np.mean(x))
    return a, c, float(rel)


#: Fraction of the target-β SD to floor τ at. This is where the **C-construction nuisance** lives — the
#: residual difference in what each cohort's confounder block removes (different reference: the measured
#: CPTAC-on-REF-B term, β ρ=0.986; different C columns: `CPE`/`target_cn` exist only in tumour). It is
#: NOT biology, but it inflates the source↔target discrepancy, so it must not be attributed to precision.
#: Without a floor the method-of-moments τ can clip to 0 (it does, for CPTAC) and the channel becomes
#: **overconfident** — s² would then rest entirely on the bootstrap se.
_TAU_FLOOR_FRAC = 0.15


def fit_gauge(src: pd.DataFrame, tgt: pd.DataFrame, *, n_folds: int = 5, n_boot: int = 300,
              seed: int = 0, source: str = "", target: str = "",
              tau_floor_frac: float = _TAU_FLOOR_FRAC,
              calibrated: Optional[bool] = None) -> Gauge:
    """Estimate the global scale `a` and the per-edge deviation SD `τ` from two β tables.

    τ is floored at `tau_floor_frac · a · sd(β_target)` — the C-construction/reference nuisance (see
    `_TAU_FLOOR_FRAC`). Pass `tau_floor_frac=0` to see the raw method-of-moments τ."""
    k = ["gene", "family"]
    m = src.merge(tgt, on=k, suffixes=("_s", "_t"))
    if m.empty:
        raise ValueError("no shared (gene, family) edges between the two cohorts")
    x, y = m["beta_t"].to_numpy(float), m["beta_s"].to_numpy(float)      # x = target, y = source
    se_t, se_s = m["se_t"].to_numpy(float), m["se_s"].to_numpy(float)
    lam = float(np.mean(se_s ** 2) / max(np.mean(se_t ** 2), 1e-12))
    a, c, rel = _eiv_slope(x, y, se_t)

    rng = np.random.default_rng(seed)                                    # bootstrap CI over edges
    boots = [_eiv_slope(x[i], y[i], se_t[i])[0]
             for i in (rng.integers(0, len(x), len(x)) for _ in range(n_boot))]
    a_lo, a_hi = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))

    # τ² by method of moments: var(resid) = τ² + se_src² + a²·se_tgt²
    resid = y - c - a * x
    tau2 = float(np.var(resid) - np.mean(se_s ** 2) - (a ** 2) * np.mean(se_t ** 2))
    tau = float(np.sqrt(max(tau2, 0.0)))
    tau = max(tau, tau_floor_frac * abs(a) * float(np.std(x)))            # C-nuisance floor (see above)

    genes = sorted(m["gene"].unique())                                   # OUT-OF-FOLD (a, c), over GENES
    gf = {g: i % n_folds for i, g in enumerate(rng.permutation(genes))}
    fold_a = {}
    for f in range(n_folds):
        keep = ~m["gene"].map(gf).eq(f).to_numpy()                       # fit on the OTHER folds
        fold_a[f] = _eiv_slope(x[keep], y[keep], se_t[keep])[:2] if keep.sum() > 10 else (a, c)

    # How much does this channel ACTUALLY add to β_T?  channel precision a²/(τ²+se_src²) as a fraction of
    # the mRNA likelihood's own 1/se_tgt².  This is the number that decides whether a channel is worth
    # building at all — and for every cross-cohort source tested it came out at ~1% (registry MH-102d).
    info = float(np.median(a ** 2 / (tau ** 2 + se_s ** 2)) / max(np.median(1.0 / np.maximum(se_t, 1e-9) ** 2), 1e-12))
    return Gauge(a=a, a_lo=a_lo, a_hi=a_hi, c=c, tau=tau, lam=lam, n=len(m),
                 rho=float(spearmanr(x, y).correlation), reliability=rel, info_ratio=info,
                 calibrated=calibrated, source=source, target=target, fold_a=fold_a, gene_fold=gf)


def a_for(gauge: Gauge, gene: str) -> tuple[float, float]:
    """The OUT-OF-FOLD (a, c) for a gene — the channel never uses a scale calibrated on its own gene."""
    f = gauge.gene_fold.get(gene)
    return gauge.fold_a.get(f, (gauge.a, gauge.c)) if f is not None else (gauge.a, gauge.c)


def se_floor(src: pd.DataFrame) -> float:
    """The per-edge SE floor: the median of the POSITIVE bootstrap SEs in the source table.

    ⚠ WHY THIS IS REQUIRED (a bug caught 2026-07-12 before wiring). Bagged NNLS pins many coefficients at
    the boundary β̂=0, and a coefficient that is 0 in *every* bootstrap resample gets `se = 0` — which is a
    BOUNDARY ARTIFACT, not precision. Combined with the intercept, such an edge produced
    `pihat = (0 − c)/a < 0` with `s² ≈ 0`: the channel would assert, at ~6σ, that β is NEGATIVE — for an
    edge the source merely had no evidence about. Under the half-normal (repression ≥ 0) gauge that
    actively crushes real edges toward zero. A boundary zero means "no evidence here", so it must carry
    TYPICAL uncertainty, not infinite confidence."""
    pos = src.loc[src["se"] > 0, "se"]
    return float(pos.median()) if len(pos) else 0.0


def to_channel_terms(gene: str, cols: Sequence[str], src: pd.DataFrame, gauge: Gauge,
                     *, min_a: float = 0.05, se_min: Optional[float] = None) -> list[dict]:
    """The `channels=` list for `spike_slab._gibbs_posterior` — one Gaussian term per family.

        pihat = (β̂_source − c) / a      s² = (τ² + max(se_src, se_min)²) / a²      members = [(col, 1.0)]

    `c` = the NNLS non-negativity floor (`_eiv_slope`); `se_min` = the boundary-zero SE floor
    (`se_floor`). Returns [] when the gauge is dead or `a` is too small to invert — the posterior then
    reverts to mRNA-only for this gene (axis P1, "drop the absent term")."""
    if not gauge.usable:                                 # dead / unreliable slope / uncalibrated
        return []
    a, c = a_for(gauge, gene)
    if a < min_a:
        return []
    if se_min is None:
        se_min = se_floor(src)
    ix = {c_: i for i, c_ in enumerate(cols)}
    terms = []
    for _, r in src[src["gene"] == gene].iterrows():
        i = ix.get(r["family"])
        if i is None:
            continue
        se = max(float(r["se"]), se_min)                 # boundary-zero ⇒ typical uncertainty, not 0
        s2 = (gauge.tau ** 2 + se ** 2) / (a ** 2)
        terms.append({"members": [(i, 1.0)], "pihat": (float(r["beta"]) - c) / a, "s2": float(s2)})
    return terms


def shuffled_gauge(source: str, target: str, tgt: pd.DataFrame, *, genes=None, seed: int = 0,
                   **kw) -> Gauge:
    """THE CROSS-COHORT SHUFFLE-NULL (shared primitive). Permute the source cohort's miRNA sample labels
    → refit β_source → refit the gauge. A real channel's `a` must COLLAPSE toward 0 here."""
    src0 = beta_table(source, genes, permute=True, seed=seed, **kw)
    return fit_gauge(src0, tgt, source=f"{source}[SHUFFLED]", target=target, seed=seed)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="gtex", choices=[c for c in CF.COHORTS if c != "tcga"])
    ap.add_argument("--target", default="tcga")
    ap.add_argument("--n-genes", type=int, default=150)
    args = ap.parse_args()

    from mirna_hallmark.learned.evidence import ledger as LG
    _, Yt = cohort_matrices(args.target)
    _, Ys = cohort_matrices(args.source)
    ed = LG.pooled_he_edges()
    genes = [g for g in ed["gene"].unique() if g in Yt.index and g in Ys.index][:args.n_genes]
    print(f"[gauge] {args.source} → {args.target} on {len(genes)} shared genes")

    tgt = beta_table(args.target, genes)
    src = beta_table(args.source, genes)
    g = fit_gauge(src, tgt, source=args.source, target=args.target)
    print(f"\n  REAL : {g}")
    print(f"         dead={g.dead}  (a_lo>0 ⇒ the source carries information about β_target)")
    print(f"         out-of-fold a: {[round(v[0],3) for v in g.fold_a.values()]}")

    gn = shuffled_gauge(args.source, args.target, tgt, genes=genes)
    print(f"  NULL : {gn}")
    print(f"         dead={gn.dead}   <-- MUST be True (a collapses) for the channel to be valid")
    print(f"\n  a_real/a_null = {g.a / gn.a:.2f}x" if abs(gn.a) > 1e-6 else "\n  null a ≈ 0 ✓")


# ============================ the CALIBRATION GATE (mandatory) ============================
def calibrate(target: str = "tcga", genes: Optional[Sequence[str]] = None, *, n_genes: int = 200,
              seed: int = 0, tol_a: float = 0.25, **bt_kw) -> Gauge:
    """SPLIT-HALF CONTROL — the gate every gauge must pass before it may emit channel terms.

    Fit β on two random halves of the SAME cohort: same C, same platform, same biology, same n.
    The truth is therefore **a = 1 and τ = 0**, exactly. Any deviation is the estimator + the slope
    machinery, not nature.

    ⚠ THIS IS NOT CEREMONY — it CAUGHT A REAL FAILURE (2026-07-12). With the dense **Gibbs posterior** as the
    β estimator, the split-half gauge returned **a = 4.1** (and 3,478,000 in one configuration) where the true
    answer is 1.0: the Gibbs posterior SDs are large enough that the errors-in-variables denominator
    `Var(β̂) − E[se²]` collapses toward zero and `a` explodes. Bagged NNLS passes (a = 1.10 raw-`r`, **1.04
    under the E3 z-Y convention**). Without this gate, a channel would have silently shipped a 4× mis-scaled
    `pihat` with no symptom. `Gauge.usable` requires `calibrated is True`.

    Passes iff |a − 1| ≤ tol_a AND τ ≤ its own null AND reliability ≥ _MIN_RELIABILITY.
    """
    from mirna_hallmark.learned.evidence import ledger as LG
    X, Y = cohort_matrices(target)
    dec = CF.deconv(target)
    parts = [p for p in Y.columns if p in X.columns and (dec is None or p in dec.index)]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(parts))
    hA = [parts[i] for i in perm[: len(parts) // 2]]
    hB = [parts[i] for i in perm[len(parts) // 2:]]
    if genes is None:
        ed = LG.pooled_he_edges()
        genes = [g for g in ed["gene"].unique() if g in Y.index][:n_genes]

    a_tab = _beta_on(target, genes, hA, seed=seed + 1, **bt_kw)
    b_tab = _beta_on(target, genes, hB, seed=seed + 2, **bt_kw)
    g = fit_gauge(a_tab, b_tab, source=f"{target}[half-A]", target=f"{target}[half-B]", seed=seed)
    ok = (abs(g.a - 1.0) <= tol_a) and (g.tau <= 1e-9 or g.tau <= _TAU_FLOOR_FRAC) \
        and (g.reliability >= _MIN_RELIABILITY)
    g.calibrated = bool(ok)
    return g


def _beta_on(cohort: str, genes, parts, *, seed: int = 0, **bt_kw) -> pd.DataFrame:
    """`beta_table` restricted to an explicit participant subset (used by `calibrate`)."""
    X, Y = cohort_matrices(cohort)
    from mirna_hallmark.learned.evidence import ledger as LG
    ed = LG.pooled_he_edges()
    C = CF.build_C(cohort, parts)
    Cdf = pd.DataFrame(np.c_[np.ones(len(parts)), C.to_numpy(float)], index=parts)
    zy = bt_kw.get("zscore_y", True)
    # same SCALE-FREE floor as beta_table (an absolute 0.2 is cohort-scale-dependent — see its docstring)
    _sds = np.array([float(Y.loc[g, parts].astype(float).std()) for g in genes if g in Y.index])
    _sds = _sds[np.isfinite(_sds) & (_sds > 0)]
    _floor = float(0.25 * np.median(_sds)) if len(_sds) else 0.0
    rows = []
    for g in genes:
        if g not in Y.index:
            continue
        y = Y.loc[g, parts].astype(float)
        k = y.notna()
        pk = [p for p, q in zip(parts, k) if q]
        if len(pk) < 40:
            continue
        y = y.loc[pk]
        if y.std() < _floor:
            continue
        regs = [m for m in ed.loc[ed["gene"] == g, "miRNA"].unique() if m in X.index]
        if len(regs) < 2:
            continue
        Xg = X.loc[regs, pk].T.astype(float).fillna(0.0)
        Xf, _, _ = FAM.collapse_by_family(Xg, pd.Series(1.0, index=regs), FAM.family_of(pd.Index(regs)))
        if Xf.shape[1] < 1:
            continue
        yr, Xz, cols = _prep(y, Xf, Cdf.loc[pk])
        if zy:
            s = float(np.std(yr))
            if s < 1e-9:
                continue
            yr = yr / s
        b, sd = _bagged_nnls_meansd(Xz, yr, n_boot=30, seed=seed)
        for c, v, s_ in zip(cols, b, sd):
            rows.append({"gene": g, "family": c, "beta": float(v), "se": float(s_), "n": len(pk)})
    return pd.DataFrame(rows)
