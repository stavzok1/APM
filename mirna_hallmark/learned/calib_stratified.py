"""IS THE POSTERIOR'S OVERCONFIDENCE DRIVEN BY SUBTYPE MIX? — random vs SUBTYPE-STRATIFIED half-splits.

    .venv/bin/python3 -m mirna_hallmark.learned.calib_stratified [--seeds 4] [--n-genes 200] [--arm random|strat]

WHY (2026-08-01). `calibration.posterior_calibration` compares the reported posterior SD against the TRUE
sampling SD measured from two INDEPENDENT half-cohorts, and finds the model OVERCONFIDENT. Measured across
4 seeds at 200 genes / 3 splits: **bagged NNLS 0.745 ± 0.009 · Gibbs Gaussian 0.815 ± 0.084 · Student-t
ν=7 0.849 ± 0.051** (ratio < 1 ⇒ overconfident; CALIBRATED = 0.9–1.1).

⭐ **THE DIAGNOSTIC QUESTION THIS ANSWERS.** The harness's own footer blames "unmodeled patient
heterogeneity", and there is independent evidence for it: **MH-165 measured that coupling GENUINELY
DIFFERS across PAM50 subtypes** (a real distributional shift). If β truly varies by subtype, then two
RANDOM halves differ partly because they drew different subtype MIXES — and a single-β model cannot
report that as uncertainty. **Prediction (pre-registered): if subtype mix is the driver, balancing the
subtype composition of the two halves should SHRINK the true sampling SD and push the ratio toward 1.**

  * ratio rises materially under stratification ⇒ the overconfidence is **subtype-mix heterogeneity** ⇒
    the fix is MODELLING (a subtype random effect on β), and a single global inflation constant is wrong.
  * ratio unchanged ⇒ the heterogeneity is **finer-grained than subtype** ⇒ modelling it is not tractable
    here and the honest fix is REPORTING (carry a calibrated width).

⚠ **STRATIFYING THE SPLIT IS NOT CONDITIONING ON SUBTYPE.** MH-102f forbids PAM50 in the confounder block
`C` — it is computed FROM the mRNA matrix and 27 of its 50 classifier genes are our own targets, so
conditioning on it conditions on a function of the outcome. Using it to BALANCE two halves does not enter
any regression: the estimate is unchanged, only WHICH patients land in which half. Legitimate, and the
distinction matters — do not let this become a precedent for putting PAM50 in `C`.

⚠ Underpowered by construction: the Gibbs ratio's own seed-SD is ~0.084, so a difference smaller than
~0.1 is not resolvable at 4 seeds. Report the spread, not a point estimate.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats


def _subtype_labels(parts):
    """PAM50 (or nearest available) per participant — used ONLY to balance the split."""
    from mirna_hallmark import data_loaders as D
    cl = D.load_clinical_strata()
    col = next((c for c in ("PAM50_final", "PAM50", "tnbc_subtype_4") if c in cl.columns), None)
    if col is None:
        return None
    key = "participant" if "participant" in cl.columns else cl.columns[0]
    s = cl.drop_duplicates(key).set_index(key)[col].reindex(parts)
    return s.fillna("NA").astype(str)


def _split(parts, labels, rng, stratified: bool):
    """Two halves — random, or with each subtype's members split ~50/50 (balanced MIX)."""
    if not stratified or labels is None:
        perm = rng.permutation(len(parts))
        return [parts[i] for i in perm[:len(parts)//2]], [parts[i] for i in perm[len(parts)//2:]]
    A, B = [], []
    for lv in labels.unique():                       # balance WITHIN each subtype => same mix in both
        mem = [p for p in parts if labels.get(p) == lv]
        idx = rng.permutation(len(mem))
        A += [mem[i] for i in idx[:len(mem)//2]]
        B += [mem[i] for i in idx[len(mem)//2:]]
    return A, B


def run(n_genes: int = 200, seeds: int = 4, arm: str = "both") -> pd.DataFrame:
    from mirna_hallmark.learned import confounders as CF
    from mirna_hallmark.learned import gauge as G
    from mirna_hallmark.learned import spike_slab as SS
    from mirna_hallmark.learned.attribution_eb import _prep
    from mirna_hallmark.learned.calibration import _gene_designs
    from mirna_hallmark.learned.evidence import ledger as LG

    X, Y = G.cohort_matrices("tcga")
    dec = CF.deconv("tcga")
    parts = [p for p in Y.columns if p in X.columns and (dec is None or p in dec.index)]
    ed = LG.pooled_he_edges()
    genes = [g for g in ed["gene"].unique() if g in Y.index][:n_genes]
    C = CF.build_C("tcga", parts)
    Cdf = pd.DataFrame(np.c_[np.ones(len(parts)), C.to_numpy(float)], index=parts)
    designs = _gene_designs("tcga", genes, parts)
    pos = {p: i for i, p in enumerate(parts)}
    labels = _subtype_labels(parts)
    print(f"[calib_stratified] {len(parts)} participants · {len(designs)} designs · "
          f"subtype labels: {'YES ' + str(dict(labels.value_counts())) if labels is not None else 'NONE'}")

    def fit(sub, seed_):
        rows = []
        for g, y, Xf, cols in designs:
            ii = np.array([pos[p] for p in sub])
            yr, Xz, _ = _prep(pd.Series(y[ii], index=sub),
                              pd.DataFrame(Xf[ii], index=sub, columns=cols), Cdf.loc[sub])
            b, sd, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=1200, burn=400, seed=seed_)
            rows += [{"edge": (g, c), "beta": float(v), "se": float(s)} for c, v, s in zip(cols, b, sd)]
        return pd.DataFrame(rows)

    arms = ("random", "stratified") if arm == "both" else (arm,)
    out = []
    for a in arms:
        for s in range(seeds):
            rng = np.random.default_rng(1000 + s)
            hA, hB = _split(parts, labels, rng, a == "stratified")
            m = fit(hA, 1).merge(fit(hB, 2), on="edge", suffixes=("_a", "_b"))
            true = float(np.std(m.beta_a - m.beta_b) / np.sqrt(2))
            rep = float(np.mean([m.se_a.mean(), m.se_b.mean()]))
            out.append({"arm": a, "seed": s, "reported_se": rep, "true_sd": true, "ratio": rep / true})
            print(f"  {a:11s} seed {s}: ratio {rep/true:.4f}  (reported {rep:.5f}, true {true:.5f})",
                  flush=True)
    R = pd.DataFrame(out)
    print("\n=== random vs SUBTYPE-STRATIFIED half-splits (Gibbs Gaussian) ===")
    for a, s in R.groupby("arm"):
        print(f"  {a:11s} n={len(s)}  ratio {s.ratio.mean():.4f} ± {s.ratio.std(ddof=1):.4f}   "
              f"true_sd {s.true_sd.mean():.5f}")
    if R.arm.nunique() == 2:
        x = R[R.arm == "stratified"].ratio.to_numpy(); y = R[R.arm == "random"].ratio.to_numpy()
        d = x.mean() - y.mean()
        se = np.sqrt(x.var(ddof=1)/len(x) + y.var(ddof=1)/len(y))
        print(f"\n  Δ(stratified − random) = {d:+.4f} ± {se:.4f}   "
              f"(Welch p={stats.ttest_ind(x, y, equal_var=False).pvalue:.3f})")
        print("  ⇒ a materially POSITIVE Δ means subtype MIX drove the excess spread (⇒ model it);")
        print("    ~zero means the heterogeneity is finer-grained (⇒ report a calibrated width instead).")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-genes", type=int, default=200)
    ap.add_argument("--arm", default="both", choices=["both", "random", "stratified"])
    a = ap.parse_args()
    run(a.n_genes, a.seeds, a.arm)
