"""THE WEIGHT GAIN, OUT-OF-FOLD — does the learned β still beat sum-abundance when it never saw the sample?

    .venv/bin/python3 -m mirna_hallmark.learned.weight_gain_oof [--workers 8] [--genes A,B]

⭐ WHY THIS EXISTS (user challenge, 2026-08-01). MH-177/MH-178 measured the β-weighted vs sum-abundance
gain using the CARD's `beta` — and `readouts.py` fits that with `_gibbs_posterior` on the **FULL sample,
with no fold splitting** (verified in code). Correlating `Σβ·X` with `Y` on the SAME patients β was fit on
is an **IN-SAMPLE** statistic, and in-sample fit rises MECHANICALLY with the number of free parameters.
Since the reported gain grew steeply with design width (`n_fam` 2 → 10+ : +0.023 → +0.153), the width
gradient could be pure parameter-count optimism rather than model quality. **That is not a caveat, it is a
confound, and it is testable.**

⚠ The OOF machinery already existed — `eval/decoy_bench.oof_budget` folds properly — it simply was never
used for this question. This module reuses that exact protocol.

THE PROTOCOL (identical to `decoy_bench.oof_budget`, extended to score TWO aggregates per fold):
    for each of 5 folds:
        fit on TRAIN only:  C's coefficients · the z-scoring mean/sd · the variance floor · the Gibbs β
        score on TEST:      agg  = Z_test · β        (the learned weighting)
                            abund = Z_test · 1       (the UNWEIGHTED sum — same arms, same z-scoring)
    rho = Spearman(pooled OOF prediction, pooled C-residualised Y)
    gain = rho_abund − rho_agg          (+ = the weighting helps, OUT OF FOLD)

⭐ The EXACT internal null survives the fold split: with one seed family β is uniform, so `agg` and
`abund` are the same vector up to a positive scale ⇒ Spearman identical ⇒ gain ≡ 0. It is re-checked here
rather than assumed, because the fold machinery is new code around it.

⚠ Family collapse is applied to BOTH arms exactly as `decoy_bench` does (a TRUE RPM pool, never a mean),
so the comparison is weighted-vs-unweighted over the SAME family design — the only difference is β.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "weight_gain_oof.tsv"
N_FOLDS = 5
GIBBS_ITER, GIBBS_BURN = 200, 80          # the decoy_bench setting, validated there (max|Δρ| <= 0.0013)


def _one(gene: str):
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.eval import decoy_bench as DB
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import spike_slab as SS

    try:
        Yv, X, Cm, _ = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return None
    if X.shape[1] < 1:
        return None
    Xf = DB.pool_family(list(X.columns), Yv.index)          # SAME family pool as the decoy bench
    n = len(Yv)
    if n < 100:
        return None
    y, Xv, Cmat = Yv.to_numpy(float), Xf.to_numpy(float), Cm.to_numpy(float)
    fold = np.random.default_rng(0).permutation(np.arange(n) % N_FOLDS)
    A, B, Yo = [], [], []
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        if te.sum() < 20 or tr.sum() < 50:
            continue
        lc = LinearRegression().fit(Cmat[tr], y[tr])            # C on TRAIN
        mu, sd = Xv[tr].mean(0), Xv[tr].std(0)                  # z-params from TRAIN
        Ztr = (Xv[tr] - mu) / (sd + 1e-9); Ztr[:, sd < 0.1] = 0.0
        Zte = (Xv[te] - mu) / (sd + 1e-9); Zte[:, sd < 0.1] = 0.0
        try:
            b, _, _ = SS._gibbs_posterior(Ztr, -(y[tr] - lc.predict(Cmat[tr])), np.ones(Ztr.shape[1]),
                                          n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
        except Exception:
            return None
        # ⭐ STANDARDISE EACH FOLD'S PREDICTION BEFORE POOLING. β is REFIT PER FOLD, so pooling raw
        # predictions concatenates differently-scaled pieces — and Spearman over the pooled vector is
        # NOT invariant to that. Measured without this: the EXACT internal null (one family ⇒ agg = b·abund)
        # FAILED at max|gain| = 9.2e-02 on 563 one-family genes. Per-fold z-scoring restores it.
        def _z(v):
            sdv = np.std(v)
            return (v - np.mean(v)) / sdv if sdv > 1e-12 else v * 0.0
        A.append(_z(Zte @ b))                                   # learned weighting
        B.append(_z(Zte @ np.ones(Zte.shape[1])))               # UNWEIGHTED sum, same design
        Yo.append(y[te] - lc.predict(Cmat[te]))
    if not A:
        return None
    a, bb, yo = np.concatenate(A), np.concatenate(B), np.concatenate(Yo)
    ra = stats.spearmanr(a, yo).correlation if np.std(a) > 1e-9 else np.nan
    rb = stats.spearmanr(bb, yo).correlation if np.std(bb) > 1e-9 else np.nan
    return {"gene": gene, "n": n, "n_fam_oof": Xf.shape[1],
            "oof_rho_agg": ra, "oof_rho_abund": rb,
            "oof_gain": (rb - ra) if np.isfinite(ra) and np.isfinite(rb) else np.nan}


def run(genes=None, workers: int = 8) -> pd.DataFrame:
    from multiprocessing import Pool

    from mirna_hallmark.eval import decoy_bench as DB
    ctx = DB._ctx()
    genes = genes or sorted({g for g in ctx["he"].gene.unique() if g in ctx["Y"].index})
    print(f"[weight_gain_oof] {len(genes):,} genes · {N_FOLDS}-fold · {workers} workers")
    with Pool(workers) as p:
        R = pd.DataFrame([r for r in p.imap_unordered(_one, genes, chunksize=4) if r])
    R.to_csv(DEST, sep="\t", index=False)
    print(f"-> {DEST}")
    return R


def report(R: pd.DataFrame) -> None:
    W = pd.read_csv(OUT / "weight_gain_profile.tsv", sep="\t")[
        ["gene", "beta_const", "beta_conc", "n_fam", "gain_tcga_rna"]]
    d = R.merge(W, on="gene", how="left")
    v = d[~d.beta_const.fillna(True)].dropna(subset=["oof_gain"])
    c = d[d.beta_const.fillna(False)].dropna(subset=["oof_gain"])

    print(f"\n=== OUT-OF-FOLD WEIGHT GAIN — {len(d):,} genes ===\n")
    print("--- internal null must survive the fold split (uniform β ⇒ agg ∝ abund) ---")
    print(f"  constant-β genes n={len(c):4d}  max|gain| {c.oof_gain.abs().max():.3e}  "
          f"{'✅ EXACT' if c.oof_gain.abs().max() < 1e-9 else '⚠ NON-ZERO'}")
    print(f"\n--- variable-β genes, n={len(v):,} ---")
    print(f"  OOF   rho_agg {v.oof_rho_agg.mean():+.4f}   rho_abund {v.oof_rho_abund.mean():+.4f}   "
          f"gain {v.oof_gain.mean():+.4f}  helps {(v.oof_gain > 0).mean():.0%}  "
          f"p={stats.wilcoxon(v.oof_gain, alternative='greater').pvalue:.2g}")
    print(f"  IN-SAMPLE (MH-177)                                    "
          f"gain {v.gain_tcga_rna.mean():+.4f}  helps {(v.gain_tcga_rna > 0).mean():.0%}")
    sh = 1 - v.oof_gain.mean() / v.gain_tcga_rna.mean() if v.gain_tcga_rna.mean() else np.nan
    print(f"  ⇒ OOF retains {v.oof_gain.mean()/v.gain_tcga_rna.mean():.2f} of the in-sample gain "
          f"({sh:.0%} was optimism)")

    print("\n--- ⭐ THE TEST THAT MATTERS: does the WIDTH gradient survive out of fold? ---")
    print(f"  {'stratum':12s} {'n':>5s} {'OOF agg':>9s} {'OOF abund':>10s} {'OOF gain':>9s} "
          f"{'IN-SAMPLE gain':>15s} {'retained':>9s}")
    for lab, sel in [("n_fam 2", v.n_fam == 2), ("n_fam 3-4", v.n_fam.between(3, 4)),
                     ("n_fam 5-9", v.n_fam.between(5, 9)), ("n_fam 10+", v.n_fam >= 10)]:
        s = v[sel]
        if len(s) < 15:
            continue
        ret = s.oof_gain.mean() / s.gain_tcga_rna.mean() if s.gain_tcga_rna.mean() else np.nan
        print(f"  {lab:12s} {len(s):5d} {s.oof_rho_agg.mean():+9.4f} {s.oof_rho_abund.mean():+10.4f} "
              f"{s.oof_gain.mean():+9.4f} {s.gain_tcga_rna.mean():+15.4f} {ret:9.2f}")
    r = stats.spearmanr(v.n_fam, v.oof_gain)
    print(f"\n  spearman(OOF gain, n_fam)   = {r.correlation:+.4f}  p={r.pvalue:.2g}")
    r2 = stats.spearmanr(v.beta_conc, v.oof_gain)
    print(f"  spearman(OOF gain, beta_conc) = {r2.correlation:+.4f}  p={r2.pvalue:.2g}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--genes")
    a = ap.parse_args()
    report(run(a.genes.split(",") if a.genes else None, a.workers))
