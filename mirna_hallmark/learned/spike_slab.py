"""Spike-and-slab inclusion model (Design §Decision D/E; LEARNED_MODEL_DESIGN_RESPONSE §910) — the
PARKED Bayesian alternative to the shipped adaptive lasso, now BUILT so it can be benchmarked head-to-head.

Where the adaptive lasso (regression.fit_gene) decides inclusion by soft L1 thresholding and lets the
evidence prior act as an INVERSE PENALTY, this puts an explicit, evidence-graded **inclusion indicator**
on every candidate edge and integrates over it:

    r_g = X_g β_g + ε_g ,   ε_g ~ N(0, σ_g²)          # r = −resid(Y|C): repression → positive
    β_m = z_m · θ_m                                    # coefficient = inclusion × magnitude
    z_m ~ Bernoulli(π_m)                               # π_m EVIDENCE-GRADED (monotone in the prior w_m)
    θ_m ~ N⁺(0, ν²)                                    # half-normal slab (θ ≥ 0: "miRNAs repress")
    σ_g², ν² ~ InvGamma

Estimation: componentwise SSVS Gibbs (Kuo–Mallick). Each edge's inclusion odds fold the evidence prior
π_m against the data marginal likelihood; the slab is a truncated-normal so the sign box (M ≥ 0) matches
the lasso. The point readout for the gate is the posterior mean coefficient M_m = E[z_m·θ_m] (inclusion
probability × conditional magnitude); PIP_m = E[z_m] is the posterior inclusion probability.

The evidence prior enters ONLY through π (inclusion) — the slab magnitude is evidence-neutral — so the
benchmark isolates exactly the design's claim: does an evidence-graded inclusion PRIOR beat the lasso's
implicit L1 selection on the same out-of-fold coupling gate?

CLI:  .venv/bin/python3 -m mirna_hallmark.learned.spike_slab            # hub panel head-to-head
      .venv/bin/python3 -m mirna_hallmark.learned.spike_slab --fdr      # genome-wide vs lasso gate
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.special import expit, log_ndtr, ndtr
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR

HUB_PANEL = ["PTEN", "GATA3", "ESR1", "ZEB1", "CDKN1A"]


def _resid(V: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    return V - LinearRegression().fit(Cmat, V).predict(Cmat)


def inclusion_prior(w: np.ndarray, *, p0: float = 0.3, kappa: float = 1.5,
                    lo: float = 0.02, hi: float = 0.98) -> np.ndarray:
    """Evidence-graded Bernoulli inclusion probability π_m, monotone in the prior weight w_m.
    π = expit( logit(p0) + κ·z ), z = standardized log1p(w) (centred on the median, robust-scaled).
    p0 = base inclusion rate at median evidence; κ = how sharply evidence moves inclusion. Clipped [lo,hi]."""
    s = np.log1p(np.clip(w, 0, None))
    med = np.median(s)
    mad = np.median(np.abs(s - med)) * 1.4826
    z = (s - med) / mad if mad > 1e-9 else np.zeros_like(s)
    pi = expit(np.log(p0 / (1 - p0)) + kappa * z)
    return np.clip(pi, lo, hi)


def _gibbs_ss(Xz: np.ndarray, r: np.ndarray, pi: np.ndarray, *,
              nu2: float = 1.0, slab_scale: np.ndarray | None = None, sample_nu2: bool = True,
              n_iter: int = 1500, burn: int = 500, seed: int = 0,
              a0: float = 1.0, b0: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Componentwise spike-and-slab Gibbs on standardized predictors Xz (n×p) and target r (n,).
    Half-normal slab N⁺(0, ν²·scale_m²), evidence-graded inclusion π. `slab_scale` (per-edge, from
    priors.slab_scale) loosens the slab for high-affinity/deep-evidence edges; None → uniform slab.
    Returns (posterior-mean β, PIP)."""
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    sc2 = np.ones(p) if slab_scale is None else np.asarray(slab_scale, dtype=float) ** 2  # per-edge scale²
    xtx = (Xz * Xz).sum(0)                                   # diag XᵀX (predictors standardized, ~n each)
    beta = np.zeros(p)
    z = (pi > 0.5).astype(float)
    theta = np.zeros(p)
    sig2 = float(np.var(r)) or 1.0                          # a0,b0 = InvGamma hyperprior on ν² and σ² (default vague)
    with np.errstate(divide="ignore"):                      # pi==0 (zero-variance arm) → −inf logit = never include
        logit_pi = np.where(pi > 0, np.log(pi) - np.log1p(-pi), -1e12)
    acc_beta = np.zeros(p)
    acc_z = np.zeros(p)
    keep = 0
    fitted = Xz @ beta                                       # running fit, kept in sync
    for it in range(n_iter):
        nu2_m = nu2 * sc2                                    # per-edge effective slab variance
        for m in range(p):
            r_m = r - fitted + beta[m] * Xz[:, m]            # partial residual excluding edge m
            A = xtx[m] / sig2
            B = float(Xz[:, m] @ r_m) / sig2
            s2 = 1.0 / (A + 1.0 / nu2_m[m])                  # posterior variance of θ_m | include
            s = np.sqrt(s2)
            mu = B * s2                                      # posterior mean of θ_m | include (untruncated)
            # log marginal-likelihood ratio for inclusion, half-normal slab truncated ≥0:
            #   L_m = 2·(s/ν)·exp(½(B·s)²)·Φ(μ/s)
            log_L = np.log(2.0) + np.log(s) - 0.5 * np.log(nu2_m[m]) + 0.5 * (B * s) ** 2 + log_ndtr(mu / s)
            log_odds = logit_pi[m] + log_L
            pin = expit(log_odds)                            # P(z_m = 1 | ·)
            if rng.random() < pin:
                z[m] = 1.0
                theta[m] = _rtnorm_pos(mu, s, rng)          # TN_{≥0}(mu, s²)
            else:
                z[m] = 0.0
                theta[m] = 0.0
            new_beta = z[m] * theta[m]
            fitted += (new_beta - beta[m]) * Xz[:, m]
            beta[m] = new_beta
        # σ² | ·
        resid = r - fitted
        sig2 = 1.0 / rng.gamma(a0 + n / 2.0, 1.0 / (b0 + 0.5 * float(resid @ resid)))
        # ν² | ·  (shared base scale, over the SCALE-standardized active slab values → keeps per-edge scale)
        if sample_nu2:
            act_idx = z > 0.5
            std_theta = theta[act_idx] / np.sqrt(sc2[act_idx])
            k = std_theta.size
            nu2 = 1.0 / rng.gamma(a0 + k / 2.0, 1.0 / (b0 + 0.5 * float(std_theta @ std_theta))) if k else nu2
            nu2 = float(np.clip(nu2, 1e-3, 1e3))
        if it >= burn:
            acc_beta += beta
            acc_z += z
            keep += 1
    return acc_beta / keep, acc_z / keep


def _gibbs_posterior(Xz: np.ndarray, r: np.ndarray, pi: np.ndarray, *,
                     n_iter: int = 2000, burn: int = 700, seed: int = 0,
                     a0: float = 1.0, b0: float = 1.0,
                     learn_pi0: bool = False, pi0_ab: tuple = (1.0, 1.0),
                     channels: list | None = None, nu: float | None = None,
                     return_samples: bool = False
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Half-normal spike-and-slab Gibbs → posterior (mean, SD, PIP) of β. π≡1 ⇒ learned-τ² ridge (all in);
    evidence-graded π ⇒ selection (un-informed edge → PIP→0).
    `learn_pi0=True` ⇒ **no hand-set inclusion prior**: a shared base rate π0 ~ Beta(pi0_ab) is LEARNED from
    the data (Beta–Bernoulli), replacing the input π's *values* for all variance-passing edges (the input π
    is used only as the variance-floor mask, π==0 → hard-excluded). The truly prior-free inclusion mode.
    `channels` = the exogenous-CHANNEL fusion (CHANNEL_FUSION_DESIGN §1): a list of Gaussian pseudo-observations of
    a *linear combination of β*, each `{'members': [(col, γ), ...], 'pihat': π̂_s, 's2': s²_π}` — the CN reduced form
    `π̂_s ~ N(Σ γ·β, s²_π)`. Each adds `A += γ_m²/s²_π`, `B += γ_m·e_s/s²_π` (`e_s = π̂_s − Σ_{other members} γ·β`) to
    column m's conditional — the exact `A_cn/B_cn` of §1. `channels=None` ⇒ current behaviour, bit-identical.
    `nu` = Student-t **observation** likelihood (CHANNEL_FUSION_DESIGN §AB / WHATS_NEXT ①): the mRNA residual is
    heavy-tailed (diagnostic: excess-kurt ~1.2 median, t-MLE df~7) from amplified-subset & near-floor tumours. `nu`
    turns the Gaussian `ε~N(0,σ²)` into `t_ν` via the scale-mixture `ε_i~N(0,σ²/λ_i), λ_i~Gamma(ν/2,ν/2)` — a per-
    observation latent precision `λ_i` down-weights outlier patients, STAYS Gibbs-conjugate (no HMC). `nu=None` ⇒
    Gaussian, **bit-identical** to §6b (no λ draws → RNG stream unchanged). Same augmentation applies to a CN channel
    term (set `t['nu']`); left Gaussian here by default."""
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    xtx = (Xz * Xz).sum(0)
    lam = np.ones(n)                                              # Student-t per-obs precision weights (nu=None ⇒ inert)
    active = pi > 0                                                # variance-floor / candidacy mask
    beta = np.zeros(p); z = (pi > 0.5).astype(float); theta = np.zeros(p)
    sig2 = float(np.var(r)) or 1.0; nu2 = 1.0; pi0 = 0.5
    ap, bp = pi0_ab
    with np.errstate(divide="ignore"):
        logit_pi = np.where(active, np.log(pi) - np.log1p(-pi), -1e12)
    s1 = np.zeros(p); s2 = np.zeros(p); sz = np.zeros(p); keep = 0
    bsamp: list = [] if return_samples else None                  # post-burn β draws (cross-family posterior COVARIANCE readout)
    fitted = Xz @ beta
    chan_by_col: list = [[] for _ in range(p)]                    # column m → [(term, γ_m), ...]  (CN channel §1)
    for t in (channels or []):
        for (col, gm) in t["members"]:
            if 0 <= col < p:
                chan_by_col[col].append((t, float(gm)))
    chan_eff_s2 = {id(t): float(t["s2"]) for t in (channels or [])}   # effective s²_π (Student-t inflates it on outlier segments; §AB)
    for it in range(n_iter):
        if learn_pi0:                                             # shared learned base rate for active edges
            logit_pi = np.where(active, np.log(pi0) - np.log1p(-pi0), -1e12)
        if nu is not None:                                        # Student-t: weight each obs by λ_i (held fixed over the β sweep)
            xw = Xz * lam[:, None]; xtxw = (Xz * xw).sum(0)
        else:
            xw = Xz; xtxw = xtx                                   # Gaussian: exact original arrays (bit-identical)
        for m in range(p):
            r_m = r - fitted + beta[m] * Xz[:, m]
            A = xtxw[m] / sig2
            B = float(xw[:, m] @ r_m) / sig2
            for (t, gm) in chan_by_col[m]:                        # CN channel: A += γ²/s²_π, B += γ·e_s/s²_π (§1)
                e_s = t["pihat"] - sum(g2 * beta[c2] for (c2, g2) in t["members"] if c2 != m)
                cs2 = chan_eff_s2[id(t)]
                A += gm * gm / cs2; B += gm * e_s / cs2
            sv = 1.0 / (A + 1.0 / nu2); s = np.sqrt(sv); mu = B * sv
            log_L = np.log(2.0) + np.log(s) - 0.5 * np.log(nu2) + 0.5 * (B * s) ** 2 + log_ndtr(mu / s)
            if rng.random() < expit(logit_pi[m] + log_L):
                z[m] = 1.0; theta[m] = _rtnorm_pos(mu, s, rng)
            else:
                z[m] = 0.0; theta[m] = 0.0
            nb = z[m] * theta[m]
            fitted += (nb - beta[m]) * Xz[:, m]; beta[m] = nb
        if learn_pi0:                                             # π0 | z ~ Beta(a + #on, b + #off) over active edges
            na = int(active.sum()); non = int(z[active].sum())
            pi0 = float(np.clip(rng.beta(ap + non, bp + na - non), 1e-3, 1 - 1e-3)) if na else pi0
        resid = r - fitted
        if nu is not None:                                       # weighted SSR (Σ λ_i ε_i²); then draw λ_i | ε,σ² ~ Gamma
            ssr = float((lam * resid * resid).sum())
            sig2 = 1.0 / rng.gamma(a0 + n / 2.0, 1.0 / (b0 + 0.5 * ssr))
            lam = rng.gamma((nu + 1.0) / 2.0, 2.0 / (nu + resid * resid / sig2))
        else:
            sig2 = 1.0 / rng.gamma(a0 + n / 2.0, 1.0 / (b0 + 0.5 * float(resid @ resid)))
        for t in (channels or []):                               # per-channel Student-t: latent scale on the CN pseudo-obs (t['nu'])
            nc = t.get("nu")
            if nc is None:
                continue
            e_full = t["pihat"] - sum(g2 * beta[c2] for (c2, g2) in t["members"])
            lam_c = rng.gamma((nc + 1.0) / 2.0, 2.0 / (nc + e_full * e_full / t["s2"]))
            chan_eff_s2[id(t)] = float(t["s2"]) / max(lam_c, 1e-6)   # outlier segment (large e) ⇒ small λ ⇒ inflated s² ⇒ down-weighted
        act = theta[z > 0.5]; k = act.size
        nu2 = 1.0 / rng.gamma(a0 + k / 2.0, 1.0 / (b0 + 0.5 * float(act @ act))) if k else nu2
        nu2 = float(np.clip(nu2, 1e-3, 1e3))
        if it >= burn:
            s1 += beta; s2 += beta * beta; sz += z; keep += 1
            if bsamp is not None:
                bsamp.append(beta.copy())
    mean = s1 / keep
    var = np.clip(s2 / keep - mean * mean, 0.0, None)
    if return_samples:
        return mean, np.sqrt(var), sz / keep, np.asarray(bsamp)
    return mean, np.sqrt(var), sz / keep


def _gibbs_blasso(Xz: np.ndarray, r: np.ndarray, *, n_iter: int = 1500, burn: int = 500, seed: int = 0,
                  r_lam: float = 1.0, delta_lam: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Non-negative **Bayesian lasso** (Park & Casella 2008) — learned-τ² with a LAPLACE prior, so shrinkage
    is per-edge and adaptive (heavy tails: strong edges escape, weak edges shrink hard) and the global
    strength λ is LEARNED, not guessed. Combines the lasso's sparse/adaptive character with the EB
    learned-shrinkage. Repress-box kept (β ≥ 0). Returns posterior (mean, sd) of β.

        r = Xz·β + ε ,  ε ~ N(0,σ²) ;  β_j ~ N⁺(0, τ_j²) ;  τ_j² ~ Exp(λ²/2) ;  λ² ~ Gamma(r_lam, δ_lam)

    Gibbs: β_j | · truncated-normal (componentwise); 1/τ_j² | · Inverse-Gaussian (numpy `wald`); λ² | ·
    Gamma (the learned global shrinkage); σ² | · InvGamma."""
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    xtx = (Xz * Xz).sum(0)
    beta = np.zeros(p); tau2 = np.ones(p)
    sig2 = float(np.var(r)) or 1.0; lam2 = 1.0
    s1 = np.zeros(p); s2 = np.zeros(p); keep = 0
    fitted = Xz @ beta
    for it in range(n_iter):
        for j in range(p):
            r_j = r - fitted + beta[j] * Xz[:, j]
            A = xtx[j] / sig2
            B = float(Xz[:, j] @ r_j) / sig2
            v = 1.0 / (A + 1.0 / tau2[j]); s = np.sqrt(v); mu = B * v
            nb = _rtnorm_pos(mu, s, rng)                       # β_j ≥ 0 (dense posterior, adaptively shrunk)
            fitted += (nb - beta[j]) * Xz[:, j]; beta[j] = nb
        # 1/τ_j² | · ~ InverseGaussian(mean=√(λ²)/β_j, shape=λ²); β_j→0 ⇒ tight prior (heavy shrinkage)
        bj = np.clip(beta, 1e-6, None)
        inv_tau2 = rng.wald(np.sqrt(lam2) / bj, lam2)
        tau2 = 1.0 / np.clip(inv_tau2, 1e-8, 1e8)
        # λ² | · ~ Gamma(r_lam + p, δ_lam + Στ_j²/2)  — the LEARNED global shrinkage
        lam2 = rng.gamma(r_lam + p, 1.0 / (delta_lam + 0.5 * float(tau2.sum())))
        resid = r - fitted
        sig2 = 1.0 / rng.gamma(1.0 + n / 2.0, 1.0 / (1.0 + 0.5 * float(resid @ resid)))
        if it >= burn:
            s1 += beta; s2 += beta * beta; keep += 1
    mean = s1 / keep
    var = np.clip(s2 / keep - mean * mean, 0.0, None)
    return mean, np.sqrt(var)


def _gibbs_ss_blasso(Xz: np.ndarray, r: np.ndarray, pi: np.ndarray, *, n_iter: int = 1500, burn: int = 500,
                     seed: int = 0, r_lam: float = 1.0, delta_lam: float = 1.0
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """**Spike-and-slab LASSO** (Ročková–George flavour): the inclusion indicator (genuine selection) AND a
    Bayesian-lasso slab (Laplace, per-edge τ_j², learned global λ²), non-negative. Combines all three ideas —
    selection (z), adaptive per-edge shrinkage (Laplace), learned strength (λ²). Returns (mean, sd, PIP).

        z_m ~ Bernoulli(π_m) ;  β_m | z_m=1 ~ N⁺(0, τ_m²) ;  τ_m² ~ Exp(λ²/2) ;  λ² ~ Gamma(r_lam, δ_lam)

    Gibbs: inclusion uses the marginal-L with the CURRENT per-edge slab variance τ_m² (conditional Gibbs);
    τ_m² updated (Inverse-Gaussian) for included edges, drawn from the Exp prior for excluded; λ² Gamma; σ²."""
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    xtx = (Xz * Xz).sum(0)
    beta = np.zeros(p); z = (pi > 0.5).astype(float); theta = np.zeros(p)
    tau2 = np.ones(p); sig2 = float(np.var(r)) or 1.0; lam2 = 1.0
    with np.errstate(divide="ignore"):
        logit_pi = np.where(pi > 0, np.log(pi) - np.log1p(-pi), -1e12)
    s1 = np.zeros(p); s2 = np.zeros(p); sz = np.zeros(p); keep = 0
    fitted = Xz @ beta
    for it in range(n_iter):
        for m in range(p):
            r_m = r - fitted + beta[m] * Xz[:, m]
            A = xtx[m] / sig2
            B = float(Xz[:, m] @ r_m) / sig2
            sv = 1.0 / (A + 1.0 / tau2[m]); s = np.sqrt(sv); mu = B * sv
            log_L = np.log(2.0) + np.log(s) - 0.5 * np.log(tau2[m]) + 0.5 * (B * s) ** 2 + log_ndtr(mu / s)
            if rng.random() < expit(logit_pi[m] + log_L):
                z[m] = 1.0; theta[m] = _rtnorm_pos(mu, s, rng)
            else:
                z[m] = 0.0; theta[m] = 0.0
            nb = z[m] * theta[m]
            fitted += (nb - beta[m]) * Xz[:, m]; beta[m] = nb
        inc = z > 0.5                                          # τ_m² | ·: InvGaussian for included, Exp prior else
        if inc.any():
            bj = np.clip(beta[inc], 1e-6, None)
            tau2[inc] = 1.0 / np.clip(rng.wald(np.sqrt(lam2) / bj, lam2), 1e-8, 1e8)
        if (~inc).any():
            tau2[~inc] = rng.exponential(2.0 / lam2, size=int((~inc).sum()))
        lam2 = rng.gamma(r_lam + p, 1.0 / (delta_lam + 0.5 * float(tau2.sum())))
        resid = r - fitted
        sig2 = 1.0 / rng.gamma(1.0 + n / 2.0, 1.0 / (1.0 + 0.5 * float(resid @ resid)))
        if it >= burn:
            s1 += beta; s2 += beta * beta; sz += z; keep += 1
    mean = s1 / keep
    var = np.clip(s2 / keep - mean * mean, 0.0, None)
    return mean, np.sqrt(var), sz / keep


def _rtnorm_pos(mu: float, s: float, rng) -> float:
    """One draw from N(mu, s²) truncated to [0, ∞) — Robert (1995), exact in the far-left tail.

    ⚠ THE BUG THIS REPLACES (MH-124, 2026-07-13). The old inverse-CDF form
    (`u = lo + (1-lo)·rand; u = clip(u, 1e-12, 1-1e-12); return mu + s·ndtri(u)`) SILENTLY BROKE ITS OWN
    SUPPORT. For `mu/s < -7.0345`, `ndtr(-mu/s)` saturates to **exactly 1.0** in float64 ⇒ `lo = 1.0` ⇒ `u`
    clips to `1-1e-12` ⇒ the "draw" becomes the DETERMINISTIC CONSTANT `mu + s·ndtri(1-1e-12)`
    = `mu + 7.0345·s`, which is **NEGATIVE** — in a model whose slab is HALF-NORMAL (β ≥ 0, repression ⇒
    positive weight). Measured: at mu/s = −7.05, 100% of draws negative and all IDENTICAL (1 unique value in
    200). It contaminated **3.15%** of persisted `readouts_edges.tsv` β (161/5117, 129 genes) and **3.09%** of
    `readouts_arm_edges.tsv` (179/5802) — and those impossible negatives were what made `share` (β_f/Σβ)
    explode to 43.7 in MH-119, which was mis-attributed to "anti-repressive biology / sign cancellation".

    Robert's algorithm is exact for every `mu/s`. Let `a = -mu/s` (the standardized truncation point):
      * `a <= 0`  → the mode is inside the support: plain rejection from N(0,1) until `z >= a`.
      * `a > 0`   → far-left-tail regime: EXPONENTIAL proposal with optimal rate
                    `α = (a + sqrt(a² + 4)) / 2`, accept w.p. `exp(-(z-α)²/2)`.
    As `mu/s → -∞` the correct draw → 0⁺ (never negative), which is what the model means by "this edge is off".
    """
    a = -mu / s
    if a <= 0.0:                                            # mode inside support — plain rejection
        for _ in range(200):
            z = rng.standard_normal()
            if z >= a:
                return float(mu + s * z)
        return float(max(mu, 0.0))                          # unreachable in practice
    alpha = 0.5 * (a + np.sqrt(a * a + 4.0))                # Robert's optimal exponential rate
    for _ in range(500):
        z = a + rng.exponential(1.0 / alpha)                # Exp(alpha) shifted to the truncation point
        if rng.random() <= np.exp(-0.5 * (z - alpha) ** 2):
            return float(mu + s * z)
    return float(mu + s * a)                                # = 0.0; the tail limit, still in support


def fit_gene_ss(Y: pd.Series, X: pd.DataFrame, C: pd.DataFrame, w_prior: pd.Series, *,
                n_iter: int = 1500, burn: int = 500, seed: int = 0,
                p0: float = 0.3, kappa: float = 1.5,
                gene: str | None = None, use_priors: bool = False, mu_gain: float = 1.0
                ) -> tuple[pd.Series, pd.Series]:
    """Spike-and-slab analogue of regression.fit_gene. Returns (M ≥ 0 posterior mean, PIP), Series[arm].
    `use_priors=True` (needs `gene`) draws π + the evidence-graded slab scale from `priors.edge_priors`
    (biochemical μ × evidence τ): high-affinity/deep-evidence edges get a looser slab, thin ones tighter."""
    Cmat = C.to_numpy(dtype=float)
    r = -_resid(Y.to_numpy(dtype=float), Cmat)             # repression → positive target
    Xr = _resid(X.to_numpy(dtype=float), Cmat)
    sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    w = w_prior.reindex(X.columns).fillna(0.0).to_numpy(dtype=float)
    slab = None
    if use_priors and gene is not None:
        from mirna_hallmark.learned import priors as PR
        pr = PR.edge_priors(gene, candidates=X.columns, w_evidence=w_prior, p0=p0, kappa=kappa)
        pi = pr["pi"].to_numpy(dtype=float)
        slab = PR.slab_scale(pr, mu_gain=mu_gain)
    else:
        pi = inclusion_prior(w, p0=p0, kappa=kappa)
    pi = np.where(sd > 1e-9, pi, 0.0)                       # zero-variance-in-fold arm → cannot enter
    beta, pip = _gibbs_ss(Xz, r, pi, slab_scale=slab, n_iter=n_iter, burn=burn, seed=seed)
    # β is on the z-scored scale; rescale to raw-abundance units so X·M matches the lasso's aggregate gauge.
    M = np.where(sd > 1e-9, beta / (sd + 1e-9), 0.0)
    return (pd.Series(M, index=X.columns, name="M"),
            pd.Series(pip, index=X.columns, name="PIP"))


def oof_gate_ss(gene: str, *, alpha: float = 0.005, folds: int = 5, seed: int = 0,
                family: bool = True, deconv: bool = False, w_prior_source: str = "ledger",
                n_iter: int = 1200, burn: int = 400, use_priors: bool = False, mu_gain: float = 1.0) -> dict:
    """OOF coupling gate for the spike-and-slab, RUN ALONGSIDE the lasso on the SAME folds/C/predictors so
    the two estimators are directly comparable. Returns both rho_model_ss and rho_model_lasso.
    `use_priors=True` swaps the uniform slab for the evidence-graded slab from priors.edge_priors."""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source=w_prior_source, deconv=deconv)
    if family:
        X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    n = len(Y)
    oof_ss = np.full(n, np.nan)
    oof_lasso = np.full(n, np.nan)
    oof_abund = np.full(n, np.nan)
    pip_sets: list = []
    lasso_sets: list = []
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for fi, (tr, te) in enumerate(kf.split(X)):
        M_ss, pip = fit_gene_ss(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w, seed=seed + fi,
                                n_iter=n_iter, burn=burn, gene=gene, use_priors=use_priors, mu_gain=mu_gain)
        M_l = LR.fit_gene(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w, alpha=alpha)
        oof_ss[te] = LR.aggregate(X.iloc[te], M_ss)
        oof_lasso[te] = LR.aggregate(X.iloc[te], M_l)
        oof_abund[te] = X.iloc[te].to_numpy(dtype=float).mean(axis=1)
        pip_sets.append(frozenset(pip[pip > 0.5].index))
        lasso_sets.append(frozenset(M_l[M_l > 0].index))

    Cmat = C.to_numpy(dtype=float)
    yr = _resid(Y.to_numpy(dtype=float), Cmat)
    rho_ss = spearmanr(_resid(oof_ss, Cmat), yr).correlation
    rho_lasso = spearmanr(_resid(oof_lasso, Cmat), yr).correlation
    rho_abund = spearmanr(_resid(oof_abund, Cmat), yr).correlation
    w_vec = w.reindex(X.columns).fillna(0.0).to_numpy(dtype=float)
    rho_curated = spearmanr(_resid(X.to_numpy(dtype=float) @ w_vec, Cmat), yr).correlation

    M_full, pip_full = fit_gene_ss(Y, X, C, w, n_iter=n_iter, burn=burn,
                                   gene=gene, use_priors=use_priors, mu_gain=mu_gain)
    top = M_full[M_full > 0].sort_values(ascending=False)
    return {
        "gene": gene, "n": n, "n_pred": X.shape[1],
        "nonzero_ss": int((pip_full > 0.5).sum()),
        "rho_abund": round(float(rho_abund), 3),
        "rho_curated": round(float(rho_curated), 3),
        "rho_lasso": round(float(rho_lasso), 3),
        "rho_ss": round(float(rho_ss), 3),
        "ss_vs_abund": bool(rho_ss < rho_abund),
        "ss_vs_lasso": bool(rho_ss <= rho_lasso + 1e-6),
        "stability_ss": round(_mean_jaccard(pip_sets), 2),
        "stability_lasso": round(_mean_jaccard(lasso_sets), 2),
        "top_ss": ", ".join(f"{m}={v:.2f}(p{pip_full[m]:.2f})" for m, v in top.head(4).items()),
    }


def _mean_jaccard(sets: list) -> float:
    js = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            js.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    return float(np.mean(js)) if js else float("nan")


def run(genes=None, *, folds: int = 5, family: bool = True, deconv: bool = False,
        use_priors: bool = False, mu_gain: float = 1.0) -> pd.DataFrame:
    genes = genes or HUB_PANEL
    rows = []
    for g in genes:
        try:
            rows.append(oof_gate_ss(g, folds=folds, family=family, deconv=deconv,
                                    use_priors=use_priors, mu_gain=mu_gain))
        except Exception as e:
            rows.append({"gene": g, "error": repr(e)[:80]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(df.to_string(index=False))
    if "ss_vs_lasso" in df:
        print(f"\nspike-and-slab: beats abundance {int(df['ss_vs_abund'].sum())}/{df['ss_vs_abund'].notna().sum()}"
              f" | matches/beats LASSO {int(df['ss_vs_lasso'].sum())}/{df['ss_vs_lasso'].notna().sum()}"
              f" | mean rho_ss {df['rho_ss'].mean():.3f} vs rho_lasso {df['rho_lasso'].mean():.3f}")
    return df


def gate_fdr(genes=None, *, family: bool = True, deconv: bool = False, w_prior_source: str = "ledger",
             out: str = "mirna_hallmark/output/learned/gate_fdr_spike_slab.tsv",
             progress: int = 50, limit: int | None = None) -> pd.DataFrame:
    """Genome-wide head-to-head: run the OOF gate for BOTH estimators on every HE gene, with the same
    one-sided partial-t p (df≈n−8) + BH/BY multiplicity as mvp.gate_fdr, for the spike-and-slab AND the
    lasso side by side. `limit` caps the gene count for a timed dry run."""
    from pathlib import Path
    from scipy import stats as _st
    from mirna_hallmark import stats as _S
    from mirna_hallmark.coupling_inference import benjamini_yekutieli
    from mirna_hallmark.learned.evidence import ledger as _LG
    genes = genes or sorted(set(_LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    rows = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[ss gate_fdr] {i}/{len(genes)}", flush=True)
        try:
            rows.append(oof_gate_ss(g, family=family, deconv=deconv, w_prior_source=w_prior_source))
        except Exception:
            pass
    df = pd.DataFrame(rows)
    df = df[df["rho_ss"].notna() & df["n"].notna()].copy()

    def _p_neg(rho, n):
        dfree = int(n) - 8
        if not (rho == rho) or dfree <= 1:
            return np.nan
        if rho >= 0:
            return 1.0
        t = rho * np.sqrt(dfree / max(1.0 - rho * rho, 1e-9))
        return float(_st.t.cdf(t, dfree))

    for tag, col in [("ss", "rho_ss"), ("lasso", "rho_lasso")]:
        df[f"p_{tag}"] = [_p_neg(getattr(r, col), r.n) for r in df.itertuples()]
        df[f"q_by_{tag}"] = benjamini_yekutieli(df[f"p_{tag}"].values)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    sig_ss = int((df["q_by_ss"] < 0.05).sum())
    sig_l = int((df["q_by_lasso"] < 0.05).sum())
    print(f"\n=== GENOME-WIDE GATE: SPIKE-AND-SLAB vs LASSO ({len(df)} genes) ===")
    print(f"FDR-significant coupling (q_BY<0.05 & rho<0):  spike-slab {sig_ss} ({100*sig_ss/max(len(df),1):.0f}%)"
          f"  |  lasso {sig_l} ({100*sig_l/max(len(df),1):.0f}%)")
    print(f"mean rho_model:  spike-slab {df['rho_ss'].mean():.3f}  |  lasso {df['rho_lasso'].mean():.3f}"
          f"  |  Δ(ss−lasso) {(df['rho_ss']-df['rho_lasso']).mean():+.3f}")
    print(f"spike-slab ≤ lasso coupling on {int((df['rho_ss']<=df['rho_lasso']+1e-6).sum())}/{len(df)} genes"
          f"  |  wrote {out}")
    return df


if __name__ == "__main__":
    import sys as _sys
    if "--fdr" in _sys.argv:
        lim = None
        if "--limit" in _sys.argv:
            lim = int(_sys.argv[_sys.argv.index("--limit") + 1])
        gate_fdr(limit=lim)
        _sys.exit(0)
    ap = argparse.ArgumentParser(description="Spike-and-slab inclusion model vs adaptive lasso")
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--deconv", action="store_true")
    ap.add_argument("--priors", action="store_true", help="use evidence-graded slab from priors.edge_priors")
    ap.add_argument("--mu-gain", type=float, default=1.0, help="strength of the biochemical slab-widening")
    a = ap.parse_args()
    run(a.genes, folds=a.folds, deconv=a.deconv, use_priors=a.priors, mu_gain=a.mu_gain)
