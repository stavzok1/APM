"""Phase 3 — program-wise hierarchical Bayes (Design §Decision A primary, §6 Phase 3). Cross-gene pooling
+ posterior uncertainty, the two payoffs Bar 3 says the prior should deliver (inclusion works, magnitude is
weak/data-driven ⇒ borrow strength and quantify it, don't hand-weight).

Model, per program (a set of genes sharing regulators), on C-residualised, z-scored abundance:
    y_g = X_g β_g + ε_g ,   ε_g ~ N(0, σ_g²)
    β(m,g) ~ N(μ_m, τ²)            # miRNA-level pooling ACROSS the genes m targets (borrow strength)
    μ_m ~ N(0, ω²) ,  τ²,σ_g² ~ InvGamma       # μ_m = miRNA's typical per-target weight
Gibbs sampler → posterior mean + **sd** of β(m,g) (the identifiability object, Decision F) and the shared
μ_m. **Inclusion** = the curated HE candidate set (Bar 3: what carries coupling); **magnitude** = pooled data.

The gate (Design §6 Phase 3): does pooling improve OOF coupling, *especially on weak genes* (CDKN1A),
vs the same model with no cross-gene sharing (μ_m≡0 = independent ridge)? CLI: `python -m mirna_hallmark.learned.hierarchical`
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM

PROGRAM = ["PTEN", "CDKN1A", "CDKN1B", "BCL2", "RB1", "E2F1", "CCND1", "MYC"]  # cell-cycle/p53, shared regulators


def _prep(genes, *, family=True, deconv=False):
    """Per gene: C-residualised y and z-scored abundance X (shared gauge). Returns aligned arrays + the
    global regulator index (union across genes; families if family=True). deconv=True puts CIBERSORTx
    non-malignant fractions in C → cell-intrinsic edges (drops stroma, e.g. miR-29→collagen)."""
    data, reg = {}, {}
    for g in genes:
        try:
            Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=deconv)
        except Exception:
            continue
        if family:
            X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
        Cm = C.to_numpy(float)
        yr = Y.to_numpy(float) - LinearRegression().fit(Cm, Y.to_numpy(float)).predict(Cm)
        Xr = X.to_numpy(float)
        Xr = Xr - LinearRegression().fit(Cm, Xr).predict(Cm)                 # residualise X on C too
        Xz = (Xr - Xr.mean(0)) / (Xr.std(0) + 1e-9)                          # shared abundance gauge
        data[g] = {"y": (yr - yr.mean()) / (yr.std() + 1e-9), "X": Xz, "cols": list(X.columns)}
        for c in X.columns:
            reg.setdefault(c, len(reg))
    return data, reg


def _gibbs(data, reg, *, pool=True, n_iter=700, burn=300, seed=0):
    """Hierarchical normal Gibbs. pool=True shares μ_m across genes; pool=False fixes μ_m≡0 (independent
    ridge = the no-pooling control). Returns per-gene posterior-mean β and (pooled) per-miRNA μ + β sd."""
    rng = np.random.default_rng(seed)
    genes = list(data)
    beta = {g: np.zeros(len(data[g]["cols"])) for g in genes}
    mu = np.zeros(len(reg)); tau2 = 1.0; sig2 = {g: 1.0 for g in genes}
    a0 = b0 = 1.0
    acc = {g: [] for g in genes}; acc_mu = []
    XtX = {g: data[g]["X"].T @ data[g]["X"] for g in genes}
    Xty = {g: data[g]["X"].T @ data[g]["y"] for g in genes}
    for it in range(n_iter):
        # 1. per-gene β | ·  (Gaussian conjugate ridge toward μ_{reg})
        for g in genes:
            idx = np.array([reg[c] for c in data[g]["cols"]])
            p = len(idx)
            V = np.linalg.inv(XtX[g] / sig2[g] + np.eye(p) / tau2)
            m = V @ (Xty[g] / sig2[g] + (mu[idx] if pool else 0.0) / tau2)
            L = np.linalg.cholesky(V)
            beta[g] = m + L @ rng.standard_normal(p)
        # 2. μ_m | ·  (only when pooling)
        if pool:
            sums = np.zeros(len(reg)); cnts = np.zeros(len(reg))
            for g in genes:
                idx = np.array([reg[c] for c in data[g]["cols"]])
                np.add.at(sums, idx, beta[g]); np.add.at(cnts, idx, 1.0)
            prec = cnts / tau2 + 1.0 / 100.0
            mean = (sums / tau2) / prec
            mu = mean + rng.standard_normal(len(reg)) / np.sqrt(prec)
        # 3. τ² | ·
        ss = 0.0; k = 0
        for g in genes:
            idx = np.array([reg[c] for c in data[g]["cols"]])
            d = beta[g] - (mu[idx] if pool else 0.0); ss += float(d @ d); k += len(idx)
        tau2 = 1.0 / rng.gamma(a0 + k / 2, 1.0 / (b0 + ss / 2))
        # 4. σ_g² | ·
        for g in genes:
            r = data[g]["y"] - data[g]["X"] @ beta[g]
            sig2[g] = 1.0 / rng.gamma(a0 + len(r) / 2, 1.0 / (b0 + float(r @ r) / 2))
        if it >= burn:
            for g in genes:
                acc[g].append(beta[g].copy())
            acc_mu.append(mu.copy())
    post = {g: {"mean": np.mean(acc[g], 0), "sd": np.std(acc[g], 0), "cols": data[g]["cols"],
                "draws": np.asarray(acc[g])} for g in genes}  # draws → posterior covariance (Decision F: the identifiability object)
    return post, (np.mean(acc_mu, 0) if pool and acc_mu else mu)


def oof(genes=None, *, folds=5, family=True, seed=0, n_sub=None) -> pd.DataFrame:
    """OOF predictive coupling per gene, POOLED (hierarchical) vs UNPOOLED (independent ridge). `n_sub`
    subsamples each gene to that many patients (pooling's value is the SMALL-n regime)."""
    genes = genes or PROGRAM
    data0, reg0 = _prep(genes, family=family)
    genes = list(data0)
    if n_sub:
        rng = np.random.default_rng(seed)
        for g in genes:
            take = rng.permutation(len(data0[g]["y"]))[:n_sub]
            data0[g] = {"y": data0[g]["y"][take], "X": data0[g]["X"][take], "cols": data0[g]["cols"]}
    rows = {g: {"pool": np.full(len(data0[g]["y"]), np.nan), "solo": np.full(len(data0[g]["y"]), np.nan)}
            for g in genes}
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    idxsplit = {g: list(kf.split(data0[g]["y"])) for g in genes}
    for f in range(folds):
        train = {g: {"y": data0[g]["y"][idxsplit[g][f][0]], "X": data0[g]["X"][idxsplit[g][f][0]],
                     "cols": data0[g]["cols"]} for g in genes}
        for pool, key in [(True, "pool"), (False, "solo")]:
            post, _ = _gibbs(train, reg0, pool=pool, seed=seed + f)
            for g in genes:
                te = idxsplit[g][f][1]
                rows[g][key][te] = data0[g]["X"][te] @ post[g]["mean"]
    out = []
    for g in genes:
        y = data0[g]["y"]
        rp = spearmanr(rows[g]["pool"], y).correlation
        rs = spearmanr(rows[g]["solo"], y).correlation
        out.append({"gene": g, "n": len(y), "n_reg": len(data0[g]["cols"]),
                    "rho_pooled": round(rp, 3), "rho_solo": round(rs, 3),
                    "delta_pool": round(rp - rs, 3), "pool_helps": bool(rp > rs)})
    df = pd.DataFrame(out)
    with pd.option_context("display.width", 160):
        print(df.to_string(index=False))
    print(f"\nPhase 3 gate: pooling helps {int(df['pool_helps'].sum())}/{len(df)} genes | "
          f"mean Δρ(pooled−solo) = {df['delta_pool'].mean():+.3f} | "
          f"weak genes (|rho_solo|<0.15): Δ = {df.loc[df['rho_solo'].abs()<0.15,'delta_pool'].mean():+.3f}")
    return df


def uncertainty(genes=None, *, family=True, deconv=False) -> pd.DataFrame:
    """Posterior mean ± sd of β for the program's regulators — the identifiability object (Decision F)."""
    genes = genes or PROGRAM
    data, reg = _prep(genes, family=family, deconv=deconv)
    post, mu = _gibbs(data, reg, pool=True, n_iter=1200, burn=500)
    rows = []
    for g in post:
        for c, m, s in zip(post[g]["cols"], post[g]["mean"], post[g]["sd"]):
            rows.append({"gene": g, "arm": c, "beta": round(float(m), 3), "sd": round(float(s), 3),
                         "z": round(float(m / s), 2) if s else np.nan})
    df = pd.DataFrame(rows)
    strong = df[df["z"].abs() > 2].sort_values("beta")
    print(f"identified edges (|β/sd|>2): {len(strong)}/{len(df)}")
    print(strong.head(12).to_string(index=False))
    return df


if __name__ == "__main__":
    oof()
