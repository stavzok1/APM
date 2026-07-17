"""Attribution: IDENTITY vs MAGNITUDE (Design §Decision I).

Two distinct objects, never conflated:
    MAGNITUDE / force  = the realized budget share  M(f)·X_fam  (states.budget_shift §5a) — "who exerts the
                         most realized pressure" (abundance-INCLUDED).
    IDENTITY / who     = the **Shapley/LMG** decomposition of the FIXED-weight aggregate's explained variance
                         R²(X·M, Y) across families — "who fairly OWNS the coupling", split under collinearity
                         (two collinear families share the credit, not double-count it). Abundance enters only
                         through the fitted aggregate, so identity ≠ realized magnitude.

For the unsaturated linear aggregate a·M, Shapley = LMG (average incremental R² over predictor orderings);
we sample orderings (exact = 2^p). **Family = the identified unit** (Design §F).

DISCIPLINE (repo): share ⊥ coupling — a high identity share is NEVER evidence of repression; a driver call
requires FDR-significant coupling (§6/§7). Identity says *who*, not *whether*.

CLI: `python -m mirna_hallmark.learned.attribution PTEN ESR1`
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd


def shapley_identity(Xr: np.ndarray, m: np.ndarray, yr: np.ndarray, *, n_perm: int = 400, seed: int = 0) -> np.ndarray:
    """Sampled Shapley (= LMG for the linear aggregate) credit of each predictor for the FIXED-weight
    aggregate's explained variance R²(Xr·m, yr). Xr: n×p (C-residualised, z-scored); m: fixed weights; yr:
    C-residualised target. Fair under collinearity — shared variance is split, not double-counted."""
    n, p = Xr.shape
    rng = np.random.default_rng(seed)

    # ⚡ INCREMENTAL ACCUMULATION (2026-07-17). The previous form rebuilt the aggregate FROM SCRATCH at every
    # step of every permutation — `pred = Xr[:, cols] @ m[cols]`, i.e. O(n·k) at depth k ⇒ O(n_perm · n · p²/2)
    # overall — and called `np.corrcoef(pred, yr)` at each of the n_perm·p steps, which re-derives yr's mean
    # and norm every time and allocates a 2×2 matrix. But along a permutation, admitting predictor j only ADDS
    # a fixed column: pred += Xr[:,j]·m[j]. So the aggregate is O(n) per step ⇒ O(n_perm · n · p) overall, a
    # ~p/2 speedup, and yr's centred form is hoisted out of the loop. MEASURED (n=1041, 400 perms):
    # **ESR1 p=15: 476 → 58 ms (8.2×) · PTEN p=64: 2405 → 247 ms (9.7×)**. VERIFIED OUTPUT-IDENTICAL to the
    # old form: max|Δφ| = 7.6e-17 (ESR1) / 6.9e-17 (PTEN) — float-summation noise only; same seed ⇒ same
    # permutations ⇒ same estimator, just evaluated without the redundant work.
    # ⚠ The realised speedup is ~8–10×, NOT the ~p/2 the flop count alone predicts (32× at p=64): the inner
    # step is now O(n)=1041 numpy ops, so per-step Python/dispatch overhead — not arithmetic — dominates.
    # Vectorising ACROSS permutations would attack that; not done (this is already 3% of a readouts run).
    yc = yr - yr.mean()
    y_norm = float(np.sqrt(yc @ yc))
    if y_norm < 1e-12:
        return np.zeros(p)
    W = Xr * m                                  # n×p — the weighted columns; admitting j is now `pred += W[:,j]`
    sqrt_n = np.sqrt(n)

    phi = np.zeros(p)
    for _ in range(n_perm):
        order = rng.permutation(p)
        pred = np.zeros(n)
        prev = 0.0
        for j in order:
            pred += W[:, j]                     # O(n), was O(n·k) rebuilt from scratch
            pc = pred - pred.mean()
            pn = float(np.sqrt(pc @ pc))
            # `np.std(pred) < 1e-12` in the old form; ||pc|| = std·sqrt(n), so this is the same test.
            now = 0.0 if pn < 1e-12 * sqrt_n else (float(pc @ yc) / (pn * y_norm)) ** 2
            phi[j] += now - prev
            prev = now
    return phi / n_perm


def _exonerate_between(Xf: pd.DataFrame, supp: list, yv: np.ndarray, Cm: np.ndarray,
                       *, rho_thr: float = 0.5, surv_thr: float = -0.1):
    """Between-family exoneration gate (ATTRIBUTION_PRIMITIVE §2). A family whose UNIQUE variation — residual on
    its abundance-collinear mates (|ρ|≥rho_thr) + C — does NOT anti-track Y (partial-Spearman ≥ surv_thr) is a
    collinear RIDER: it earns Shapley credit only through the shared variance it co-varies with. Drop it from the
    split. Families with no collinear mate cannot be riders → kept. Returns (survivors, exonerated)."""
    from scipy.stats import spearmanr
    from sklearn.linear_model import LinearRegression
    corr = Xf[supp].corr(method="spearman").abs()
    surv, exon = [], []
    for f in supp:
        mates = [o for o in supp if o != f and corr.loc[f, o] >= rho_thr]
        if not mates:
            surv.append(f); continue                            # nothing to ride on → not a rider
        Z = np.column_stack([Cm, Xf[mates].to_numpy(float)])
        xf = Xf[f].to_numpy(float)
        xr = xf - LinearRegression().fit(Z, xf).predict(Z)
        if float(np.std(xr)) < 1e-9:
            exon.append(f); continue                            # perfectly collinear with a mate → rider
        yr = yv - LinearRegression().fit(Z, yv).predict(Z)
        rho = spearmanr(xr, yr).correlation
        (surv if (rho == rho and rho < surv_thr) else exon).append(f)
    return surv, exon


def bayes_shapley_identity(gene: str, *, n_draw: int = 120, n_perm: int = 60, seed: int = 0) -> pd.DataFrame:
    """BAYESIAN Shapley identity (CHANNEL_FUSION §2.5 general reframe, MH-94): the same identity/ownership as
    `identity_vs_magnitude`, but the two binary/point pieces are made Bayesian, from ONE inclusion posterior:
      • continuous ENTRY — the binary `_exonerate_between` rider gate → the inclusion **PIP** (P(family enters); a
        collinear rider with no unique coupling gets β≈0 in most draws ⇒ ~0 Shapley automatically, no hard drop).
      • Shapley with WIDTH — the point Shapley on a fixed M → Shapley computed **per posterior draw** of β → a
        posterior mean ± sd of each family's identity share (vs the point estimate's hidden uncertainty).
    Channel-agnostic: any channel (CN, protein) that informs the posterior flows in here too. Returns per family:
    pip, identity_mean ± identity_sd (Bayesian), and the point identity_share for comparison."""
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD, families as FAM, states as ST, spike_slab as SS
    from mirna_hallmark.learned import attribution_eb as AE
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, _, _ = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    M = ST.canonical_M(gene, "01", arm_level=False)
    supp = [f for f in M.index if M.get(f, 0) > 0 and f in Xf.columns]
    if len(supp) < 2:
        print(f"{gene}: <2 nonzero families — identity undefined"); return pd.DataFrame()
    Cm = C.to_numpy(float)
    _resid = lambda v: v - LinearRegression().fit(Cm, v).predict(Cm)
    yr = _resid(Y.to_numpy(float))
    Xz = ((Xf[supp] - Xf[supp].mean()) / (Xf[supp].std(ddof=0) + 1e-9)).fillna(0.0)
    Xr = np.column_stack([_resid(Xz[f].to_numpy(float)) for f in supp])
    pi = np.clip(AE._evidence_pi(gene, supp), 0.0, 1.0)
    # ONE inclusion posterior over the FULL support (no pre-exoneration): β draws (sign-matched to canonical M, β≥0) + PIP
    _, _, pip, samp = SS._gibbs_posterior(Xr, -yr, pi, n_iter=1400, burn=400, seed=seed, return_samples=True)
    idx = np.linspace(0, len(samp) - 1, min(n_draw, len(samp))).astype(int)         # subsample draws
    shares = np.zeros((len(idx), len(supp)))
    for k, i in enumerate(idx):
        phi = shapley_identity(Xr, samp[i], yr, n_perm=n_perm, seed=i).clip(min=0)   # per-draw Shapley
        shares[k] = phi / phi.sum() if phi.sum() > 0 else phi
    df = pd.DataFrame({"pip": np.round(pip, 3),
                       "identity_mean": shares.mean(0).round(3), "identity_sd": shares.std(0).round(3),
                       "M": M[supp].round(3).to_numpy()}, index=supp).sort_values("identity_mean", ascending=False)
    print(f"\n=== {gene} — BAYESIAN Shapley identity (n_draw={len(idx)}): continuous-entry PIP + Shapley ± sd ===")
    print(df.to_string())
    print(f"  honest width: mean identity_sd = {df['identity_sd'].mean():.3f}; "
          f"widest = {df['identity_sd'].idxmax()} ({df.loc[df['identity_sd'].idxmax(),'identity_mean']}±{df['identity_sd'].max()})")
    return df


def identity_vs_magnitude(gene: str, *, n_perm: int = 400) -> pd.DataFrame:
    """The identity/magnitude plane for a gene (Decision I), at the seed-family unit.
    identity_share = Shapley credit for the coupling (collinearity-fair OWNERSHIP);
    magnitude_share = budget M·X_fam (realized FORCE). gap = identity − magnitude:
      + = owns more coupling than its realized budget (unique/abundance-light owner);
      − = high budget but shared/redundant coupling (collinear passenger)."""
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned import states as ST
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, _, _ = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    M = ST.canonical_M(gene, "01", arm_level=False)
    supp = [f for f in M.index if M.get(f, 0) > 0 and f in Xf.columns]
    if len(supp) < 2:
        print(f"{gene}: <2 nonzero families — identity undefined"); return pd.DataFrame()
    Cm = C.to_numpy(float)

    def _resid(v):
        return v - LinearRegression().fit(Cm, v).predict(Cm)

    yr = _resid(Y.to_numpy(float))
    Xz = ((Xf[supp] - Xf[supp].mean()) / (Xf[supp].std(ddof=0) + 1e-9)).fillna(0.0)
    from scipy.stats import spearmanr
    # EXONERATION GATE (ATTRIBUTION_PRIMITIVE §2): drop collinear riders BEFORE the Shapley — a family whose
    # unique variation doesn't anti-track Y would otherwise be handed shared-variance credit Shapley can't withhold.
    survivors, exon = _exonerate_between(Xf, supp, Y.to_numpy(float), Cm)
    if not survivors:
        survivors, exon = supp, []                              # guard: never blank a gene
    Xr_s = np.column_stack([_resid(Xz[f].to_numpy(float)) for f in survivors])
    phi_s = pd.Series(shapley_identity(Xr_s, M[survivors].to_numpy(float), yr, n_perm=n_perm), index=survivors).clip(lower=0)
    ident = pd.Series(0.0, index=supp)                          # exonerated riders → identity 0
    if phi_s.sum() > 0:
        ident[survivors] = phi_s / phi_s.sum()
    ident = ident.round(3)
    a = Xf[supp].mean()                                         # mean family abundance (level)
    mag = (M[supp] * a).clip(lower=0)                           # budget: force is force — NOT exonerated
    magsh = (mag / mag.sum()).round(3) if mag.sum() > 0 else mag
    Xr = np.column_stack([_resid(Xz[f].to_numpy(float)) for f in supp])
    ind = pd.Series([spearmanr(Xr[:, i], yr).correlation for i in range(len(supp))], index=supp).round(3)
    df = pd.DataFrame({"M": M[supp].round(3), "indiv_rho": ind, "identity_share": ident, "magnitude_share": magsh})
    df["gap"] = (df["identity_share"] - df["magnitude_share"]).round(3)
    df["exonerated"] = [f in exon for f in supp]
    df = df.sort_values("identity_share", ascending=False)
    tot_r2 = float(phi_s.sum())
    with pd.option_context("display.width", 150):
        print(f"\n=== {gene} — IDENTITY (Shapley, rider-exonerated) vs MAGNITUDE (budget), R²={tot_r2:.3f}, "
              f"{len(exon)}/{len(supp)} exonerated ===")
        print(df.to_string())
    if exon:
        print(f"  ⚑ exonerated collinear riders (budget>0 but no unique coupling): {exon}")
    print(f"  identity↔magnitude rank corr = {df['identity_share'].corr(df['magnitude_share']):+.2f}  "
          f"| quiet owner: {df['gap'].idxmax()} (+{df['gap'].max():.2f}); loud passenger: {df['gap'].idxmin()} ({df['gap'].min():.2f})")
    return df


if __name__ == "__main__":
    for g in ([a for a in sys.argv[1:] if not a.startswith("-")] or ["PTEN", "ESR1"]):
        identity_vs_magnitude(g)
