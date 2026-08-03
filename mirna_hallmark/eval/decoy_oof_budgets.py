"""OUT-OF-FOLD real & decoy BUDGETS on the CURRENT matched pairs — the artifact the decoy arc lacked.

    .venv/bin/python3 -m mirna_hallmark.eval.decoy_oof_budgets [--workers 6]

⭐⭐ WHY THIS EXISTS — A STALE-ARTIFACT DEFECT, CAUGHT BY THE USER ASKING *"why are there so few n?"*.
`decoy_bench` persists a per-gene **ρ** but never the per-sample **budget**, so anything wanting to score a
decoy against another layer (protein), or within a stratum (PAM50), had to fall back on
`decoy_family_betas.tsv`. That table is dated **2026-07-13** — three days older than
`decoy_full_pairs.tsv` (07-16) and three weeks older than the current `decoy_bench.tsv` (08-02).

**MEASURED CONSEQUENCE: only 63 of its 3,558 (gene, family) decoy cells — 1.8% — correspond to the CURRENT
Hungarian-matched partners, and it covers just 516 of 1,395 genes.** So the "few n" was never a power limit;
it was the footprint of a superseded artifact. ⚠ And the 07-13 decoys predate MH-137's three fixes (the
126× evidence hole, the Poisson-gated 6mer, the dose caliper), whose recorded effect is
**−0.0306 (6mer-contaminated) → −0.0119 (fixed)** — i.e. **pre-fix decoys are too WEAK and flatter real**.
Any result built on them is biased toward the real arm. (This is axiom 2: a shared artifact consumed without
checking its provenance.)

WHAT THIS FIXES, BOTH DEFECTS AT ONCE
------------------------------------
1. **The CURRENT pairs** — `decoy_full_pairs.tsv` (1,395 genes), i.e. site-free by TargetScan context++ ∧
   scanMiR ∧ Poisson-gated 6mer ∧ the evidence union, cluster-excluded, sparse-residual decorrelated, and
   Hungarian-matched on signed loadings + dose + variance under a dose caliper.
2. **OUT-OF-FOLD budgets for BOTH arms** — β refit inside each training fold and applied to held-out
   patients, for the real design *and* the decoy design, with **identical** folds, C, prep and estimator.
   That removes the in-sample inflation that made an earlier in-sample gap read −0.043 where
   `decoy_bench`'s strict-OOF figure is −0.012 to −0.015.

⚠ **PER-FOLD STANDARDISATION BEFORE POOLING (MH-181).** `Z·β` is an arbitrary-scale index and β is refit per
fold, so raw concatenation glues differently-scaled pieces together and the pooled statistic is not
scale-invariant. Both arms are z-scored within each fold.

⭐ The output is a per-sample budget for each arm, which is what lets the SAME decoy control be reused
against mRNA, against RPPA protein, and inside a PAM50 stratum — the three places MH-203 needed it.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "decoy_oof_budgets.parquet"
#: ⭐ drop-in replacement for the STALE `decoy_family_betas.tsv` (2026-07-13) — same schema, current pairs.
BETA_DEST = OUT / "decoy_family_betas_oof.tsv"
PAIRS = OUT / "decoy_full_pairs.tsv"
N_FOLDS = 5
GIBBS_ITER, GIBBS_BURN = 200, 80
_MEM: dict = {}


def _pool(arms, index, Xall):
    """TRUE RPM family pool over the given arms — `log2(1 + Σ(2^x − 1))`, never a mean."""
    from mirna_hallmark.coupling_inference import seed_family_map
    fam = _MEM.get("fammap")
    if fam is None:
        fam = seed_family_map(pd.Index(Xall.index)); _MEM["fammap"] = fam
    by: dict = {}
    for a in arms:
        by.setdefault(str(fam.get(a)), []).append(a)
    cols = {}
    for f, mem in by.items():
        if f == "None":
            continue
        lin = np.nan_to_num(np.power(2.0, Xall.loc[mem, index].to_numpy(float)) - 1.0, nan=0.0)
        cols[f] = np.log2(1.0 + np.clip(lin, 0, None).sum(axis=0))
    return pd.DataFrame(cols, index=index) if cols else None


def _one(gene: str):
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.learned import attribution_eb as AE
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned import spike_slab as SS

    pairs = _MEM["pairs"]
    sub = pairs[pairs.gene == gene]
    if sub.empty:
        return None
    try:
        Y, X, Cm, w = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return None
    if len(Y) < 200:
        return None
    Xall = LD._load()["X"]
    Xr, _, _ = FAM.collapse_by_family(X, w, FAM.family_of(pd.Index(X.columns)))
    Xd = _pool([a for a in sub.fake_arm if a in Xall.index], Y.index, Xall)
    if Xd is None or Xd.shape[1] < 1 or Xr.shape[1] < 1:
        return None
    y, Cmat, n = Y.to_numpy(float), Cm.to_numpy(float), len(Y)
    fold = np.random.default_rng(0).permutation(np.arange(n) % N_FOLDS)
    out = {"real": np.full(n, np.nan), "decoy": np.full(n, np.nan)}
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        if te.sum() < 20 or tr.sum() < 100:
            continue
        lc = LinearRegression().fit(Cmat[tr], y[tr])
        yr = -(y[tr] - lc.predict(Cmat[tr]))
        for lab, M in (("real", Xr), ("decoy", Xd)):
            v = M.to_numpy(float)
            mu, sd = v[tr].mean(0), v[tr].std(0)
            Ztr = (v[tr] - mu) / (sd + 1e-9); Ztr[:, sd < 0.1] = 0.0
            Zte = (v[te] - mu) / (sd + 1e-9); Zte[:, sd < 0.1] = 0.0
            try:
                b, _, _ = SS._gibbs_posterior(Ztr, yr, np.ones(Ztr.shape[1]),
                                              n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
            except Exception:
                return None
            m = Zte @ b
            s = np.std(m)
            out[lab][te] = (m - np.mean(m)) / s if s > 1e-12 else 0.0   # ⚠ per-fold z (MH-181)
    if not np.isfinite(out["real"]).any() or not np.isfinite(out["decoy"]).any():
        return None
    # ⭐ ALSO emit a FULL-DATA decoy β per family — the drop-in replacement for the stale
    # `decoy_family_betas.tsv`. Full-data (not per-fold) because a TRANSPORTED decoy must be a frozen
    # weight vector applied in another cohort, exactly like the real family-card β it is compared against.
    yr_all = -(y - LinearRegression().fit(Cmat, y).predict(Cmat))
    v = Xd.to_numpy(float)
    mu, sd = v.mean(0), v.std(0)
    Z = (v - mu) / (sd + 1e-9); Z[:, sd < 0.1] = 0.0
    try:
        bfull, sfull, _ = SS._gibbs_posterior(Z, yr_all, np.ones(Z.shape[1]),
                                              n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
    except Exception:
        bfull = sfull = np.full(Xd.shape[1], np.nan)
    betas = pd.DataFrame({"gene": gene, "unit": [str(c) for c in Xd.columns], "decoy": True,
                          "beta": bfull, "beta_sd": sfull})
    return (pd.DataFrame({"gene": gene, "participant": list(Y.index),
                          "budget_real": out["real"], "budget_decoy": out["decoy"],
                          "n_fam_real": Xr.shape[1], "n_fam_decoy": Xd.shape[1]}), betas)


def run(workers: int = 6, genes=None) -> pd.DataFrame:
    from multiprocessing import get_context

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.learned import data as LD

    pairs = pd.read_csv(PAIRS, sep="\t")
    _MEM["pairs"] = pairs
    genes = genes or sorted(pairs.gene.unique())
    yidx = set(LD._load()["Y"].index)
    D.load_cnv_target_genes(sorted(set(genes) & yidx))      # ONE batched warm (axiom 3a)
    _pool(list(pairs.fake_arm.unique()[:2]), list(LD._load()["X"].columns)[:5], LD._load()["X"])  # prime fammap
    print(f"[decoy_oof] {len(genes):,} genes on the CURRENT pairs · {workers} workers", flush=True)
    with get_context("fork").Pool(workers) as p:
        res = [r for r in p.imap_unordered(_one, genes, chunksize=4) if r is not None]
    T = pd.concat([r[0] for r in res], ignore_index=True)
    Bt = pd.concat([r[1] for r in res], ignore_index=True)
    T.to_parquet(DEST)
    Bt.to_csv(BETA_DEST, sep="\t", index=False)
    print(f"-> {DEST}   {T.gene.nunique():,} genes × {T.participant.nunique():,} participants")
    print(f"-> {BETA_DEST}   {len(Bt):,} decoy (gene,family) cells over {Bt.gene.nunique():,} genes"
          f"   [replaces the STALE decoy_family_betas.tsv: 3,558 cells / 516 genes]")
    return T


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    run(a.workers)
