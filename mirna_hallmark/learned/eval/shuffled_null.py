"""Bar 3 — shuffled-evidence null (Design §5 Bar 3 / §0 ref 3): does the prior do REAL work?

Permute the prior→edge assignment **within each gene's regulators** (preserves the marginal weight
distribution = the degree/strata-preserving shuffle) and refit. If the TRUE evidence→weight map couples
held-out target mRNA better than random reassignments, the prior carries signal beyond "any sparse
structure would do." This guards against the adaptive-lasso's coupling being purely data-driven (the lasso
fits coefficients regardless of the prior — so a shuffled prior still couples via the data; the question is
whether the REAL assignment couples *better*).

Metric = OOF `rho_model` (the Bar-1 gate coupling; more negative = better). Null = `n_shuffle` permutations
per gene. Permutation p = fraction of nulls at least as negative as real (small ⇒ prior does work). Also a
cohort z (real − null_mean)/null_sd, negative = real beats the null mean.

    python -m mirna_hallmark.learned.eval.shuffled_null --prior ledger --family
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import mvp as M
from mirna_hallmark.learned import regression as LR


def _safe_rho(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman ρ, nan if either side is constant (degenerate fold: all-zero M → constant aggregate)."""
    if np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return np.nan
    return spearmanr(a, b).correlation


def bar3_gene(gene: str, *, n_shuffle: int = 100, alpha: float = 0.005, folds: int = 5,
              w_prior_source: str = "ledger", family: bool = True, seed: int = 0) -> dict:
    """Assemble the gene ONCE, then evaluate OOF coupling for the real prior and `n_shuffle` in-memory
    permutations of it (same folds, same data) — the efficient Bar-3 null."""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source=w_prior_source)
    if family:
        X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    Cmat = C.to_numpy(float)
    yr = M._resid(Y.to_numpy(float), Cmat)
    splits = list(KFold(n_splits=folds, shuffle=True, random_state=seed).split(X))

    def oof_rho(wv: pd.Series) -> float:
        oof = np.full(len(Y), np.nan)
        for tr, te in splits:
            Mv = LR.fit_gene(Y.iloc[tr], X.iloc[tr], C.iloc[tr], wv, alpha=alpha)
            oof[te] = LR.aggregate(X.iloc[te], Mv)
        return _safe_rho(M._resid(oof, Cmat), yr)

    real_rho = oof_rho(w)
    rng = np.random.default_rng(seed)
    nulls = np.array([oof_rho(pd.Series(rng.permutation(w.to_numpy()), index=w.index))
                      for _ in range(n_shuffle)])
    nulls = nulls[np.isfinite(nulls)]
    if not len(nulls) or real_rho != real_rho:
        return {"gene": gene, "real_rho": real_rho, "n_null": len(nulls)}
    p = (np.sum(nulls <= real_rho) + 1) / (len(nulls) + 1)         # fraction of nulls ≥ as negative as real
    z = (real_rho - nulls.mean()) / (nulls.std() or np.nan)
    return {"gene": gene, "real_rho": round(real_rho, 4), "null_mean": round(float(nulls.mean()), 4),
            "null_sd": round(float(nulls.std()), 4), "delta": round(real_rho - float(nulls.mean()), 4),
            "z": round(float(z), 2), "p_perm": round(float(p), 4), "n_null": len(nulls),
            "prior_works": bool(p <= 0.05)}


def bar3_inclusion_gene(gene: str, *, n_shuffle: int = 100, alpha: float = 0.005, folds: int = 5,
                        w_prior_source: str = "ledger", seed: int = 0) -> dict:
    """Inclusion (π) null: does the curated HE arm-SET couple better than random arm-sets of the same
    size? Uniform weights on BOTH (isolates inclusion from magnitude). Random arms drawn from the full
    expressed-arm universe. Caveat: unmatched on abundance/variance, so a positive result is inclusion +
    some ascertainment; a null result is strong evidence inclusion adds nothing."""
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source=w_prior_source)
    real_arms = list(X.columns)
    k = len(real_arms)
    Xall = LD._load()["X"]                                        # arm × participant
    universe = [a for a in Xall.index if a not in set(real_arms)]
    Cmat = C.to_numpy(float)
    yr = M._resid(Y.to_numpy(float), Cmat)
    splits = list(KFold(n_splits=folds, shuffle=True, random_state=seed).split(X))

    def oof_rho_arms(arms) -> float:
        Xr = Xall.loc[arms, X.index].T.fillna(0.0)                # → participant × arm
        wv = pd.Series(1.0, index=arms)                           # uniform prior — isolate inclusion
        oof = np.full(len(Y), np.nan)
        for tr, te in splits:
            Mv = LR.fit_gene(Y.iloc[tr], Xr.iloc[tr], C.iloc[tr], wv, alpha=alpha)
            oof[te] = LR.aggregate(Xr.iloc[te], Mv)
        return _safe_rho(M._resid(oof, Cmat), yr)

    real_rho = oof_rho_arms(real_arms)
    rng = np.random.default_rng(seed)
    nulls = np.array([oof_rho_arms(list(rng.choice(universe, k, replace=False)))
                      for _ in range(n_shuffle)])
    nulls = nulls[np.isfinite(nulls)]
    if not len(nulls) or real_rho != real_rho:
        return {"gene": gene, "real_rho": real_rho, "n_null": len(nulls)}
    p = (np.sum(nulls <= real_rho) + 1) / (len(nulls) + 1)
    z = (real_rho - nulls.mean()) / (nulls.std() or np.nan)
    return {"gene": gene, "k_arms": k, "real_rho": round(real_rho, 4),
            "null_mean": round(float(nulls.mean()), 4), "delta": round(real_rho - float(nulls.mean()), 4),
            "z": round(float(z), 2), "p_perm": round(float(p), 4), "n_null": len(nulls),
            "prior_works": bool(p <= 0.05)}


def run(genes=None, *, n_shuffle: int = 100, alpha: float = 0.005, folds: int = 5,
        w_prior_source: str = "ledger", family: bool = True, test: str = "magnitude") -> pd.DataFrame:
    genes = genes or M.HUB_PANEL
    fn = bar3_inclusion_gene if test == "inclusion" else bar3_gene
    kw = {} if test == "inclusion" else {"family": family}
    print(f"[Bar 3] {'INCLUSION (π): HE arm-set vs random arm-sets' if test=='inclusion' else 'MAGNITUDE (μ): shuffle prior→edge weights'}")
    rows = [fn(g, n_shuffle=n_shuffle, alpha=alpha, folds=folds, w_prior_source=w_prior_source, **kw)
            for g in genes]
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170):
        print(df.to_string(index=False))
    if "p_perm" in df:
        # combined evidence: how many genes show the prior beating its null, and Fisher-combined p
        works = int(df["prior_works"].sum())
        valid = df["p_perm"].notna()
        chi = -2 * np.log(df.loc[valid, "p_perm"].clip(lower=1e-6)).sum()
        from scipy.stats import chi2
        fisher_p = chi2.sf(chi, 2 * int(valid.sum()))
        print(f"\nBar 3: prior beats shuffled null (p≤0.05) on {works}/{int(valid.sum())} genes | "
              f"Fisher-combined p = {fisher_p:.2e} | mean Δρ(real−null) = {df['delta'].mean():+.4f} "
              f"(negative = real couples better)")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Bar 3 — shuffled-evidence null")
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--n-shuffle", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.005)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--prior", default="ledger",
                    choices=["evidence_score", "ledger", "ledger_mrna", "scanmir", "fused"])
    ap.add_argument("--family", action="store_true")
    ap.add_argument("--test", choices=["magnitude", "inclusion"], default="magnitude",
                    help="magnitude = shuffle prior weights (μ); inclusion = HE arm-set vs random sets (π)")
    a = ap.parse_args()
    run(a.genes, n_shuffle=a.n_shuffle, alpha=a.alpha, folds=a.folds,
        w_prior_source=a.prior, family=a.family, test=a.test)
