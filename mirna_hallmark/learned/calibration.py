"""How much can the model's NUMBERS be trusted? Two resampling diagnostics.

Both were forced by MH-102d's split-half control, and both answer questions the model cannot answer
about itself.

────────────────────────────────────────────────────────────────────────────────────────────────────
1. POSTERIOR CALIBRATION (`posterior_calibration`) — are the reported uncertainties honest?
────────────────────────────────────────────────────────────────────────────────────────────────────
Fit β on two INDEPENDENT halves of the cohort. Random halves have the same expected β, so
`sd(β̂_A − β̂_B)/√2` is an unbiased estimate of the estimator's TRUE sampling SD. A calibrated `se`
should match it.

    MEASURED (2026-07-12, TCGA, 150 genes, averaged over 3 random splits):
      estimator              reported se   true sampling sd   ratio   INFLATION FACTOR
      bagged NNLS              0.0056           0.0075        0.74×       1.35×
      Gibbs (Gaussian)         0.0142           0.0245        0.58×       1.73×
      Gibbs Student-t ν=7      0.0139           0.0205        0.68×       1.47×   ← helps, doesn't fix

**The reported widths are 35–70% too narrow.** MH-92's Student-t likelihood **partially** repairs it
(0.58 → 0.68) but leaves the posterior still 1.47× too confident ⇒ heavy residual tails are *part* of the
story, not all of it. The residue is **unmodeled between-participant heterogeneity** (the likelihood treats
patients as iid given X and C — they are not). Multiply reported widths by `inflation_factor` to be honest.

⚠ **THIS IS WHY SBC IS THE WRONG TOOL** (and the board item was mis-specified). Simulation-based
calibration asks "is the sampler correct *given the model*?" — it would **PASS**, and passing would be
misleading, because the fault is the MODEL, not the sampler. A correct sampler on a misspecified model
still yields narrow-but-wrong posteriors. The honest instrument is resampling, which is this function.

**WHO IS AFFECTED** (triage, so this is not over-read):
  • width claims → **MH-94 Bayesian-Shapley "honest width"**, **MH-98 δ-pooling confidence**. Understated.
    (Note the DIRECTION: wider intervals make MH-94's non-identifiability claim STRONGER, not weaker.)
  • inverse-variance FUSION WEIGHTS (MH-98) → **SAFE**: they depend on SE *ratios*, which cancel.
  • channel `s²` → immaterial (the CN channel carries ~0.7% weight).
  • OOF / permutation results (`prolif_verdict` 2×2, discovery FDR, coupling ρ) → **UNTOUCHED**: they never
    use the model's SEs. That covers most of the program's headline findings.

────────────────────────────────────────────────────────────────────────────────────────────────────
2. STABILITY SELECTION (`stability_selection`) — is a top-N discovery list reproducible?
────────────────────────────────────────────────────────────────────────────────────────────────────
Measured: the top-100 edge ranking reproduces only **63/100** between independent half-cohorts. That is a
half-cohort number (n≈520); Spearman–Brown puts the FULL-cohort reliability at `2r/(1+r)` = 2(0.728)/1.728
≈ **0.84**. So the picture is better than the raw number — but the churn is real and, until now,
**unquantified**, which is the actual problem.

Fix = Meinshausen–Bühlmann **stability selection**: bootstrap the participants, refit, and record how often
each edge lands in the top-N. Churn becomes a **per-edge confidence** that can be published alongside the
ranking. **A top-N list should not ship without it.**

CLI: `python -m mirna_hallmark.learned.calibration [--what calib|stability] [--n-genes 200]`
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import gauge as G

_CACHE: dict = {}


def _gene_designs(cohort: str, genes: Sequence[str], parts: Sequence[str]):
    """Per-gene (y, X_fam) on a fixed participant set — hoisted ONCE (the batch-loop axiom: the cost is
    loads and family collapse, not the numerics)."""
    key = (cohort, tuple(genes), tuple(parts))
    if key in _CACHE:
        return _CACHE[key]
    from mirna_hallmark.learned.evidence import ledger as LG
    X, Y = G.cohort_matrices(cohort)
    ed = LG.pooled_he_edges()
    # SCALE-FREE variability floor (an absolute 0.2 dropped 34% of TCGA genes on a cohort-specific scale —
    # same bug class as gauge.beta_table's, found by the CPTAC session 2026-07-12).
    _sds = np.array([float(Y.loc[g, parts].astype(float).std()) for g in genes if g in Y.index])
    _sds = _sds[np.isfinite(_sds) & (_sds > 0)]
    _floor = float(0.25 * np.median(_sds)) if len(_sds) else 0.0
    out = []
    for g in genes:
        if g not in Y.index:
            continue
        y = Y.loc[g, parts].astype(float)
        if y.std() < _floor:
            continue
                                                         # still drops ~34% of TCGA genes (median sd 0.233), so the
                                                         # calibration/stability numbers are computed on the variable
                                                         # two-thirds. Coverage caveat, NOT a bias across cohorts.
        regs = [m for m in ed.loc[ed["gene"] == g, "miRNA"].unique() if m in X.index]
        if len(regs) < 2:
            continue
        Xg = X.loc[regs, parts].T.astype(float).fillna(0.0)
        Xf, _, _ = FAM.collapse_by_family(Xg, pd.Series(1.0, index=regs), FAM.family_of(pd.Index(regs)))
        if Xf.shape[1] < 1:
            continue
        out.append((g, y.to_numpy(float), Xf.to_numpy(float), list(Xf.columns)))
    _CACHE[key] = out
    return out


def _fit_rows(designs, Cd: np.ndarray, idx: np.ndarray):
    """NNLS β for every gene on a participant subset `idx`. The C-residualiser is factorised ONCE per
    subset (shared by all genes) — the whole point of hoisting."""
    Q, _ = np.linalg.qr(Cd[idx])                                   # shared projector for this resample
    betas = []
    for g, y, Xf, cols in designs:
        yv = y[idx]
        yr = -(yv - Q @ (Q.T @ yv))                                # C-residualised, repression gauge
        Xs = Xf[idx]
        sd = Xs.std(0)
        Xz = np.where(sd > 1e-9, (Xs - Xs.mean(0)) / (sd + 1e-9), 0.0)
        Xz[:, sd < 0.1] = 0.0                                      # matches attribution_eb._prep
        b, _ = nnls(Xz, yr)
        for c, v in zip(cols, b):
            betas.append(((g, c), float(v)))
    return betas


# ───────────────────────────────── 1. posterior calibration ─────────────────────────────────
def posterior_calibration(cohort: str = "tcga", genes: Optional[Sequence[str]] = None, *,
                          n_genes: int = 150, n_splits: int = 3, seed: int = 0) -> pd.DataFrame:
    """Reported `se` vs the TRUE sampling SD (from independent half-cohorts), per estimator.
    ratio < 1 ⇒ the model claims more precision than it has."""
    from mirna_hallmark.learned import spike_slab as SS
    from mirna_hallmark.learned.attribution_eb import _prep, _bagged_nnls_meansd
    from mirna_hallmark.learned.evidence import ledger as LG
    X, Y = G.cohort_matrices(cohort)
    dec = CF.deconv(cohort)
    parts = [p for p in Y.columns if p in X.columns and (dec is None or p in dec.index)]
    if genes is None:
        ed = LG.pooled_he_edges()
        genes = [g for g in ed["gene"].unique() if g in Y.index][:n_genes]
    C = CF.build_C(cohort, parts)
    Cdf = pd.DataFrame(np.c_[np.ones(len(parts)), C.to_numpy(float)], index=parts)
    designs = _gene_designs(cohort, genes, parts)
    pos = {p: i for i, p in enumerate(parts)}

    def fit(sub, est, nu, seed_):
        rows = []
        for g, y, Xf, cols in designs:
            ii = np.array([pos[p] for p in sub])
            yr, Xz, _ = _prep(pd.Series(y[ii], index=sub),
                              pd.DataFrame(Xf[ii], index=sub, columns=cols), Cdf.loc[sub])
            if est == "gibbs":
                b, sd, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=1200, burn=400,
                                               seed=seed_, nu=nu)
            else:
                b, sd = _bagged_nnls_meansd(Xz, yr, n_boot=30, seed=seed_)
            for c, v, s in zip(cols, b, sd):
                rows.append({"edge": (g, c), "beta": float(v), "se": float(s)})
        return pd.DataFrame(rows)

    rng = np.random.default_rng(seed)
    out = []
    for lab, est, nu in (("bagged NNLS", "nnls", None), ("Gibbs (Gaussian)", "gibbs", None),
                         ("Gibbs Student-t nu=7", "gibbs", 7.0)):
        reps, trues = [], []
        for s in range(n_splits):                                  # average over several random splits
            perm = rng.permutation(len(parts))
            hA = [parts[i] for i in perm[: len(parts) // 2]]
            hB = [parts[i] for i in perm[len(parts) // 2:]]
            m = fit(hA, est, nu, 1).merge(fit(hB, est, nu, 2), on="edge", suffixes=("_a", "_b"))
            trues.append(float(np.std(m.beta_a - m.beta_b) / np.sqrt(2)))
            reps.append(float(np.mean([m.se_a.mean(), m.se_b.mean()])))
        rep, true = float(np.mean(reps)), float(np.mean(trues))
        out.append({"estimator": lab, "reported_se": rep, "true_sampling_sd": true,
                    "ratio": rep / true, "inflation_factor": true / rep,
                    "verdict": "CALIBRATED" if 0.9 <= rep / true <= 1.1 else
                               ("OVERCONFIDENT" if rep / true < 0.9 else "underconfident")})
    return pd.DataFrame(out).set_index("estimator")


# ───────────────────────────────── 2. stability selection ─────────────────────────────────
def stability_selection(cohort: str = "tcga", genes: Optional[Sequence[str]] = None, *,
                        n_genes: int = 300, n_boot: int = 40, top_n: Sequence[int] = (25, 50, 100),
                        seed: int = 0) -> pd.DataFrame:
    """Meinshausen–Bühlmann stability selection over the CROSS-GENE edge ranking.

    Bootstrap the participants → refit every gene → re-rank all edges → record top-N membership.
    Returns per edge: the full-cohort β and rank, plus `freq_topN` = the fraction of bootstrap resamples
    in which it made the top-N. **`freq` is the number to publish beside any top-N list** — an edge at
    rank 12 with freq 0.35 is a coin-flip, not a discovery."""
    from mirna_hallmark.learned.evidence import ledger as LG
    X, Y = G.cohort_matrices(cohort)
    dec = CF.deconv(cohort)
    parts = [p for p in Y.columns if p in X.columns and (dec is None or p in dec.index)]
    if genes is None:
        ed = LG.pooled_he_edges()
        genes = [g for g in ed["gene"].unique() if g in Y.index][:n_genes]
    C = CF.build_C(cohort, parts)
    Cd = np.c_[np.ones(len(parts)), C.to_numpy(float)]
    designs = _gene_designs(cohort, genes, parts)
    n = len(parts)

    full = dict(_fit_rows(designs, Cd, np.arange(n)))              # the full-cohort ranking
    order = sorted(full, key=lambda e: -full[e])
    rank = {e: i + 1 for i, e in enumerate(order)}

    rng = np.random.default_rng(seed)
    hits = {N: {e: 0 for e in full} for N in top_n}
    for b in range(n_boot):
        idx = rng.integers(0, n, n)                                # bootstrap participants
        bb = dict(_fit_rows(designs, Cd, idx))
        ob = sorted(bb, key=lambda e: -bb[e])
        for N in top_n:
            for e in ob[:N]:
                hits[N][e] += 1

    rows = []
    for e in full:
        r = {"gene": e[0], "family": e[1], "beta": round(full[e], 5), "rank": rank[e]}
        for N in top_n:
            r[f"freq_top{N}"] = round(hits[N][e] / n_boot, 3)
        rows.append(r)
    df = pd.DataFrame(rows).sort_values("rank").reset_index(drop=True)
    df["stable"] = df[f"freq_top{max(top_n)}"] >= 0.6              # the publishable core
    return df


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--what", default="both", choices=["calib", "stability", "both"])
    ap.add_argument("--n-genes", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=40)
    a = ap.parse_args()

    if a.what in ("calib", "both"):
        print("=" * 92)
        print("POSTERIOR CALIBRATION — are the reported uncertainties honest?")
        print("  (true sampling SD from INDEPENDENT half-cohorts; a calibrated se should MATCH it)")
        print("=" * 92)
        with pd.option_context("display.width", 140):
            print(posterior_calibration(n_genes=min(a.n_genes, 150)).round(4).to_string())
        print("\n  ratio < 1 ⇒ OVERCONFIDENT. Multiply reported widths by `inflation_factor` to be honest.")
        print("  ⚠ SBC would NOT catch this: the fault is the MODEL (unmodeled patient heterogeneity),")
        print("    not the sampler — and Student-t (MH-92) does NOT repair it.")

    if a.what in ("stability", "both"):
        print("\n" + "=" * 92)
        print("STABILITY SELECTION — is the top-N discovery ranking reproducible?")
        print("=" * 92)
        df = stability_selection(n_genes=a.n_genes, n_boot=a.n_boot)
        print(f"  {len(df)} edges over {df.gene.nunique()} genes, {a.n_boot} bootstraps\n")
        with pd.option_context("display.width", 140):
            print(df.head(20).to_string(index=False))
        for N in (25, 50, 100):
            top = df[df["rank"] <= N]
            print(f"\n  top-{N:<4} mean freq_top{N} = {top[f'freq_top{N}'].mean():.2f}   "
                  f"edges with freq ≥0.6: {(top[f'freq_top{N}'] >= 0.6).sum()}/{N}   "
                  f"(a freq of 0.5 = a coin-flip)")
        print(f"\n  ⇒ PUBLISH `freq` BESIDE ANY TOP-N LIST. {int((~df.head(100)['stable']).sum())}/100 of the "
              f"top-100 are NOT stable (freq_top100 < 0.6).")


# ───────────────────── 3. RE-DERIVED Shapley width (MH-94's "honest width", corrected) ─────────────────────
def shapley_resampled_width(gene: str, *, n_boot: int = 30, n_perm: int = 60, n_iter: int = 900,
                            seed: int = 0) -> pd.DataFrame:
    """The HONEST width of a Bayesian-Shapley identity share — re-derived, not rescaled.

    `attribution.bayes_shapley_identity` (MH-94) reports `identity_sd` = the spread of the Shapley share
    across POSTERIOR DRAWS. But the posterior itself is overconfident (`posterior_calibration`: Gibbs is
    0.58× the true sampling SD), so that width understates.

    ⚠ **AND IT CANNOT BE FIXED BY MULTIPLYING.** The calibration factor was measured on **β** — a linear
    coefficient — whereas a Shapley share is a **bounded, non-linear functional** of β on [0,1]. Naively
    scaling PTEN miR-141/200a's `0.77 ± 0.41` by the 1.73× β-factor gives ±0.71, i.e. [0.06, **1.48**] —
    past the bound. (An earlier registry row asserted exactly that ±0.71 and was RETRACTED.)

    The honest width is the **sampling SD of the share itself**: bootstrap the participants → refit the
    posterior → recompute the Shapley on the posterior mean → take the SD across resamples. Returns both,
    so the understatement is visible:
        identity_mean · posterior_sd (what MH-94 reports) · bootstrap_sd (honest) · ratio
    """
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD, families as FAM, states as ST, spike_slab as SS
    from mirna_hallmark.learned import attribution_eb as AE
    from mirna_hallmark.learned.attribution import shapley_identity

    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, _, _ = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    M = ST.canonical_M(gene, "01", arm_level=False)
    supp = [f for f in M.index if M.get(f, 0) > 0 and f in Xf.columns]
    if len(supp) < 2:
        print(f"{gene}: <2 nonzero families — identity undefined"); return pd.DataFrame()
    pi = np.clip(AE._evidence_pi(gene, supp), 0.0, 1.0)

    def _shares(idx, sd=False, seed_=0):
        Cm = C.to_numpy(float)[idx]
        yv = Y.to_numpy(float)[idx]
        yr = yv - LinearRegression().fit(Cm, yv).predict(Cm)
        Xs = Xf[supp].to_numpy(float)[idx]
        Xz = (Xs - Xs.mean(0)) / (Xs.std(0) + 1e-9)
        Xr = np.column_stack([c - LinearRegression().fit(Cm, c).predict(Cm) for c in Xz.T])
        if sd:                                            # posterior-draw spread (what MH-94 reports)
            _, _, _, samp = SS._gibbs_posterior(Xr, -yr, pi, n_iter=1400, burn=400, seed=seed_,
                                                return_samples=True)
            ii = np.linspace(0, len(samp) - 1, 60).astype(int)
            sh = np.array([(lambda p: p / p.sum() if p.sum() > 0 else p)(
                shapley_identity(Xr, samp[i], yr, n_perm=n_perm, seed=i).clip(min=0)) for i in ii])
            return sh.mean(0), sh.std(0)
        b, _, _ = SS._gibbs_posterior(Xr, -yr, pi, n_iter=n_iter, burn=n_iter // 3, seed=seed_)
        phi = shapley_identity(Xr, b, yr, n_perm=n_perm, seed=seed_).clip(min=0)
        return (phi / phi.sum() if phi.sum() > 0 else phi), None

    n = len(Y)
    mean_full, post_sd = _shares(np.arange(n), sd=True, seed_=seed)     # full data, posterior width
    rng = np.random.default_rng(seed)
    boots = np.array([_shares(rng.integers(0, n, n), seed_=b)[0] for b in range(n_boot)])
    boot_sd = boots.std(0)

    df = pd.DataFrame({"identity_mean": mean_full.round(3),
                       "posterior_sd (MH-94)": post_sd.round(3),
                       "bootstrap_sd (HONEST)": boot_sd.round(3),
                       "understatement": np.round(boot_sd / np.maximum(post_sd, 1e-6), 2)},
                      index=supp).sort_values("identity_mean", ascending=False)
    print(f"\n=== {gene} — Bayesian-Shapley identity: RE-DERIVED width ({n_boot} bootstraps) ===")
    print(df.to_string())
    print(f"\n  posterior width understates the true sampling width by "
          f"{np.median(boot_sd / np.maximum(post_sd, 1e-6)):.2f}× (median over families).")
    return df
