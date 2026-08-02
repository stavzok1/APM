"""Subtype-specific coupling test (user 2026-07-05 — fills the subtype gap with precursor stats).
Whole-cohort M gives the aggregate pressure ρ = X·M; then per gene:
  • per-PAM50-subtype coupling — partial-Spearman(pressure, target) within each subtype (descriptive);
  • **coupling SLOPE differs?** — nested rank-OLS interaction F-test: rank(target) ~ rank(pressure)*subtype + C
    vs the common-slope model (this is the formal "subtype-specific coupling" test), BH-FDR across genes;
  • **pressure LEVEL differs?** — `stats.kruskal_across_strata` on X·M across subtypes (the named precursor test).

CLI: `python -m mirna_hallmark.learned.subtype [GENES...]`
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy import stats as sst

from mirna_hallmark import stats as S
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import regression as LR
from mirna_hallmark.learned import states as ST

SUBS = ["LumA", "LumB", "Basal", "Her2"]


def subtype_coupling(genes, *, alpha: float = 0.005, deconv: bool = False) -> pd.DataFrame:
    pam = ST._pam50()
    rows = []
    for g in genes:
        try:
            Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=deconv)
            # ⭐ SWITCHED off the RETIRED adaptive lasso -> canonical Gibbs drop-in (MH-184, 2026-08-01).
            # This call's M is the REPORTED quantity (ESTIMAND class in `eval/_retired_lasso_audit.py`),
            # so the retired estimator was a real defect here. ⚠ Any persisted output is STALE until re-run.
            M = LR.fit_gene_bayes(Y, X, C, w, alpha=alpha)
        except Exception:
            continue
        press = X.to_numpy(float) @ M.reindex(X.columns).fillna(0.0).to_numpy()
        y = Y.to_numpy(float); Cm = C.to_numpy(float)
        st = np.array([pam.get(p) for p in Y.index], dtype=object)
        keep = np.isin(st, SUBS)
        if keep.sum() < 120:
            continue
        press, y, Cm, st = press[keep], y[keep], Cm[keep], st[keep]
        present = [s for s in SUBS if (st == s).sum() >= 25]
        if len(present) < 2:
            continue
        cpl = {s: float(sst.spearmanr(press[st == s], y[st == s]).correlation) for s in present}
        # interaction F-test on ranks (subtype-specific slope)
        pr, yr = sst.rankdata(press), sst.rankdata(y)
        Dm = np.column_stack([(st == s).astype(float) for s in present[1:]])   # subtype dummies (ref = present[0])
        n = len(y)
        base = np.column_stack([np.ones(n), Cm, pr, Dm])                       # common slope
        inter = pr[:, None] * Dm                                               # pressure × subtype
        full = np.column_stack([base, inter])

        def _rss(Xd):
            b, *_ = np.linalg.lstsq(Xd, yr, rcond=None)
            return float(((yr - Xd @ b) ** 2).sum())

        r0, r1 = _rss(base), _rss(full); q = inter.shape[1]; dfe = n - full.shape[1]
        F = ((r0 - r1) / q) / (r1 / dfe) if (r1 > 0 and dfe > 0 and q > 0) else np.nan
        p_int = float(sst.f.sf(F, q, dfe)) if F == F else np.nan
        _, p_kr, _ = S.kruskal_across_strata(pd.Series(press), pd.Series(st))  # pressure LEVEL differs?
        rows.append({"gene": g, **{s: round(cpl.get(s, np.nan), 2) for s in SUBS},
                     "F_int": round(F, 1) if F == F else np.nan,
                     "p_coupling_differs": p_int, "p_pressure_differs": p_kr})
    df = pd.DataFrame(rows)
    if len(df):
        df["q_coupling_differs"] = S.bh_fdr(df["p_coupling_differs"].values)
        df["q_pressure_differs"] = S.bh_fdr(df["p_pressure_differs"].values)
        df = df.sort_values("p_coupling_differs")
        with pd.option_context("display.width", 190):
            print("=== subtype-specific coupling (slope interaction F-test) + pressure level (Kruskal) ===")
            print(df.round(4).to_string(index=False))
        print(f"\ncoupling SLOPE differs across subtypes (q<0.1): {int((df['q_coupling_differs']<0.1).sum())}/{len(df)} | "
              f"pressure LEVEL differs (q<0.1): {int((df['q_pressure_differs']<0.1).sum())}/{len(df)}")
    return df


def subtype_wiring(gene, *, n_boot: int = 40, min_n: int = 60) -> pd.DataFrame:
    """Fit the CANONICAL bagged family M SEPARATELY WITHIN each PAM50 subtype → does the WIRING (M) itself
    differ by subtype? This is the direct test `subtype_coupling` cannot give (it uses one whole-cohort M, so
    a subtype-averaged wiring can look invariant). Small-n subtypes (Basal ~180, Her2 ~80) → bagged NNLS +
    family collapse for stability (same estimator as cross-state ΔM). Reports per-family M by subtype + the
    mean pairwise cross-subtype M correlation (high = conserved wiring; low = subtype-remodeled)."""
    import itertools
    from mirna_hallmark.learned import families as FAM
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    pam = ST._pam50()
    st = pd.Series([pam.get(p) for p in Y.index], index=Y.index)
    fam = FAM.family_of(pd.Index(X.columns))
    Ms = {}
    for s in SUBS:
        m = (st == s).to_numpy()
        if m.sum() < min_n:
            continue
        Xf, _, _ = FAM.collapse_by_family(X.loc[m], pd.Series(1.0, index=X.columns), fam)
        Ms[s] = ST._bagged_nnls(Y[m], Xf, C.loc[m].to_numpy(float), n_boot=n_boot)
    if len(Ms) < 2:
        print(f"{gene}: <2 subtypes with n>={min_n}"); return pd.DataFrame()
    shared = sorted(set.intersection(*[set(m.index) for m in Ms.values()]))
    df = pd.DataFrame({s: Ms[s].reindex(shared) for s in Ms}).round(3)
    corrs = [sst.spearmanr(df[a], df[b]).correlation for a, b in itertools.combinations(df.columns, 2)]
    mc = float(np.nanmean(corrs)) if corrs else np.nan
    # n-matched NOISE FLOOR: fit M on random cohort subsamples matched to each subtype's n (IGNORING subtype)
    # → cross-draw M corr = pure small-n estimation noise at these sizes. Real remodeling only if mc << floor.
    ns = {s: int((st == s).sum()) for s in Ms}
    rng = np.random.default_rng(0); allidx = np.arange(len(Y)); nulls = []
    for _ in range(3):
        Mr = {}
        for s, nn in ns.items():
            samp = rng.choice(allidx, nn, replace=False)
            Xf, _, _ = FAM.collapse_by_family(X.iloc[samp], pd.Series(1.0, index=X.columns), fam)
            Mr[s] = ST._bagged_nnls(Y.iloc[samp], Xf, C.iloc[samp].to_numpy(float), n_boot=n_boot).reindex(shared)
        nulls += [sst.spearmanr(Mr[a], Mr[b]).correlation for a, b in itertools.combinations(Mr, 2)]
    floor = float(np.nanmean(nulls)) if nulls else np.nan
    verdict = "REAL subtype remodeling" if (mc == mc and floor == floor and mc < floor - 0.1) else \
              "within small-n NOISE (not distinguishable from n-matched null)"
    df = df.sort_values(df.columns[0], ascending=False)
    print(f"\n=== {gene} — WITHIN-subtype canonical M ({', '.join(f'{s} n={ns[s]}' for s in Ms)}) ===")
    print(df[df.abs().sum(axis=1) > 0].head(14).to_string())
    print(f"cross-subtype M corr {mc:.2f}  vs  n-matched noise floor {floor:.2f}  →  {verdict}")
    return df


def subtype_wiring_pooled(gene, *, lam: float = 10.0, min_n: int = 60, n_boot: int = 40, n_null: int = 8) -> pd.DataFrame:
    """H — Bayesian state-nesting (subtype axis; Design §Decision H). Fit whole-cohort M_all (bagged NNLS) as
    the PRIOR MEAN, then per PAM50 subtype `M_s = M_all + δ_s` via ridge-toward-M_all:
        δ_s = argmin_δ ‖ (yr − Xz·M_all)[s] − Xz[s]·δ ‖² + λ‖δ‖²     (closed form)
    Δ is regularised so small-n noise shrinks toward 0 — the POOLED alternative to `subtype_wiring`, which
    differences two independent noisy fits (its low noise floor). Remodeling magnitude ‖δ_s‖₁ is compared to an
    n-matched null (random cohort subsamples), which is now ~0 by construction (ridge shrinks random δ), so a
    truly-remodeled subtype separates cleanly."""
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import families as FAM
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, _, _ = FAM.collapse_by_family(X, pd.Series(1.0, index=X.columns), fam)
    Cm = C.to_numpy(float)
    yr = -(Y.to_numpy(float) - LinearRegression().fit(Cm, Y.to_numpy(float)).predict(Cm))   # whole-cohort resid target
    sd = Xf.std(ddof=0)
    Xz = ((Xf - Xf.mean()) / (sd + 1e-9)).fillna(0.0)
    Xz.loc[:, sd < 0.1] = 0.0
    Xzv = Xz.to_numpy(float)
    m_all = ST._bagged_nnls(Y, Xf, Cm, n_boot=n_boot).reindex(Xf.columns).fillna(0.0).to_numpy()   # PRIOR MEAN
    p = Xzv.shape[1]; I = np.eye(p)
    resid_all = yr - Xzv @ m_all

    def _delta(idx):
        Xi = Xzv[idx]
        return np.linalg.solve(Xi.T @ Xi + lam * I, Xi.T @ resid_all[idx])

    pam = ST._pam50()
    st = np.array([pam.get(p_) for p_ in Y.index], dtype=object)
    rng = np.random.default_rng(0); allidx = np.arange(len(Y))
    cols = {"M_all": np.round(m_all, 3)}
    summ = []
    for s in SUBS:
        idx = np.where(st == s)[0]
        if len(idx) < min_n:
            continue
        d_s = _delta(idx)
        cols[f"M_{s}"] = np.round(m_all + d_s, 3)
        nulls = [np.abs(_delta(rng.choice(allidx, len(idx), replace=False))).sum() for _ in range(n_null)]
        summ.append({"subtype": s, "n": len(idx), "remodel_L1": round(float(np.abs(d_s).sum()), 3),
                     "null_L1": round(float(np.mean(nulls)), 3),
                     "verdict": "REAL remodeling" if np.abs(d_s).sum() > 1.5 * np.mean(nulls) else "within noise"})
    df = pd.DataFrame(cols, index=Xf.columns)
    df = df[df.abs().sum(axis=1) > 0].sort_values("M_all", ascending=False)
    print(f"\n=== {gene} — POOLED subtype wiring (ridge-toward-M_all, λ={lam}) ===")
    print(df.head(12).to_string())
    print(pd.DataFrame(summ).to_string(index=False))
    return pd.DataFrame(summ)


def _pcor(a, b, Cm):
    """partial corr of ranks a,b | C."""
    D = np.column_stack([np.ones(len(a)), Cm])
    ar = a - D @ np.linalg.lstsq(D, a, rcond=None)[0]; br = b - D @ np.linalg.lstsq(D, b, rcond=None)[0]
    return float(np.corrcoef(ar, br)[0, 1]) if ar.std() > 1e-9 and br.std() > 1e-9 else np.nan


def _het_ftest_p(rx, ry, Cm, Dm):
    """confound-mod-adjusted subtype-interaction F-test p: rank(y) ~ rx*subtype + C + rx×C (base) vs
    +rx×subtype (full) — the subtype slope-heterogeneity NET of pressure×C composition slope-modification."""
    n = len(ry); prC = rx[:, None] * Cm
    base = np.column_stack([np.ones(n), Cm, rx[:, None], Dm, prC]); full = np.column_stack([base, rx[:, None] * Dm])
    q = Dm.shape[1]; dfe = n - full.shape[1]
    if q <= 0 or dfe <= 0:
        return np.nan
    def rss(Xd):
        b, *_ = np.linalg.lstsq(Xd, ry, rcond=None); return float(((ry - Xd @ b) ** 2).sum())
    F = ((rss(base) - rss(full)) / q) / (rss(full) / dfe)
    return float(sst.f.sf(F, q, dfe)) if F == F and F > 0 else np.nan


def genomewide_edge_heterogeneity(*, n_null: int = 20, seed: int = 0):
    """PER-EDGE (arm→gene) HE coupling heterogeneity across PAM50 subtypes, genome-wide (MH-165).
    Estimand (apm-edge-question E2, realized coupling): RAW single-arm abundance predictor, partial
    Spearman(x_arm, target | C); repression = NEGATIVE ρ. F-test nets pressure×C slope-modification.
    ⭐ Calibrated against the MARGINAL-PRESERVING null (within-subtype y-permutation — preserves each
    subtype's variance/range, breaks only coupling; the honest floor. A label-shuffle is anti-conservative,
    MH-154 pattern). VERDICT (MH-165): heterogeneity REAL (16.1% vs 5.2% floor, 3.1×; best_sub does NOT
    track n → not a detectability artifact; ~half the edges sign-flip) — a DISTRIBUTIONAL shift, tag P,
    NOT validated per-edge lists (Her2/Basal n small). Surfacing (ρ_pool≈0 & subtype ρ<−0.15) is at chance
    (337 vs 287±35, z=+1.4) and family-incoherent; wiring (subtype_wiring_pooled) weak (~12% over placebo)."""
    from mirna_hallmark import config as CFG
    from mirna_hallmark.learned import families as FAM
    card = pd.read_csv(CFG.OUTPUT_ROOT / "learned" / "readouts_arm_edges.tsv", sep="\t")
    HE = {g: list(s["arm"]) for g, s in card.groupby("gene")}
    famof = FAM.family_of(pd.Index(card["arm"].unique()))
    pam = ST._pam50(); rng = np.random.default_rng(seed)
    genes, rows = [], []
    for g in sorted(HE):
        try:
            Y, X, C, _ = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)
        except Exception:
            continue
        if X is None or X.shape[1] == 0:
            continue
        idx = list(Y.index); st = np.array([pam.get(p) for p in idx], dtype=object)
        present = [s for s in SUBS if (st == s).sum() >= 25]
        if np.isin(st, SUBS).sum() < 120 or len(present) < 2:
            continue
        Cm = C.to_numpy(float); y = Y.to_numpy(float); ry = sst.rankdata(y)
        Dm = np.column_stack([(st == s).astype(float) for s in present[1:]])
        sub_idx = {s: np.where(st == s)[0] for s in present}
        arms = [(a, sst.rankdata(X[a].to_numpy(float))) for a in HE[g]
                if a in X.columns and np.std(X[a].to_numpy(float)) > 1e-9]
        for a, rx in arms:
            rho_s = {s: _pcor(rx[st == s], ry[st == s], Cm[st == s]) for s in present}
            best = min(rho_s, key=lambda s: (rho_s[s] if rho_s[s] == rho_s[s] else 9))
            rows.append({"gene": g, "arm": a, "family": famof.get(a), "rho_pool": _pcor(rx, ry, Cm),
                         **{f"rho_{s}": rho_s.get(s, np.nan) for s in SUBS}, "best_sub": best,
                         "rho_best": rho_s[best], "p_adj": _het_ftest_p(rx, ry, Cm, Dm)})
        genes.append((y, Cm, Dm, sub_idx, arms))
    R = pd.DataFrame(rows); R["q_adj"] = S.bh_fdr(R.p_adj.values)
    obs = float((R.p_adj < 0.05).mean())
    floor = []
    for _ in range(n_null):                                    # marginal-preserving null
        sig = tot = 0
        for y, Cm, Dm, sub_idx, arms in genes:
            yp = y.copy()
            for s, ii in sub_idx.items():
                yp[ii] = y[ii][rng.permutation(len(ii))]
            ryp = sst.rankdata(yp)
            for a, rx in arms:
                p = _het_ftest_p(rx, ryp, Cm, Dm)
                if p == p:
                    tot += 1; sig += (p < 0.05)
        floor.append(sig / tot)
    floor = np.array(floor)
    out = CFG.OUTPUT_ROOT / "learned" / "subtype_edge_heterogeneity.tsv"; R.to_csv(out, sep="\t", index=False)
    nflip = int(((R[[f"rho_{s}" for s in SUBS]].min(axis=1) < -0.05) & (R[[f"rho_{s}" for s in SUBS]].max(axis=1) > 0.05)).sum())
    print(f"HETEROGENEITY (marginal-preserving null): obs p<0.05 {obs:.1%} ({(R.p_adj<0.05).sum()}/{len(R)}) "
          f"vs floor {floor.mean():.1%}±{floor.std():.1%} = {obs/floor.mean():.1f}×  |  q_adj<0.1 {int((R.q_adj<0.1).sum())}  "
          f"|  sign-flip edges {nflip}  |  best_sub(sig) {R[R.q_adj<0.1].best_sub.value_counts().to_dict()}  ->  {out}")
    return R, obs, floor


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if "--genomewide" in sys.argv:
        genomewide_edge_heterogeneity()
    elif "--pooled" in sys.argv:
        for g in (args or ["PTEN", "RB1", "ZEB1", "CDKN1A"]):
            subtype_wiring_pooled(g)
    elif "--wiring" in sys.argv:
        for g in (args or ["RB1", "PTEN", "ZEB1", "CDKN1A"]):
            subtype_wiring(g)
    else:
        subtype_coupling(args or ["PTEN", "ESR1", "ZEB1", "GATA3", "BCL2", "CDKN1A", "MET", "MYC", "RB1", "EZH2"])
