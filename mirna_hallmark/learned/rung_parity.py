"""RUNG PARITY — the FAMILY-rung coupling the card never had, and a WITHIN-CELL arm identity.

    .venv/bin/python3 -m mirna_hallmark.learned.rung_parity [--workers 8]

Second delivery on the both-rungs directive (first was `arm_rung.py`, MH-186). Two gaps it closes:

⭐ **(1) THE CARD MIXES RUNGS AND NEVER SAID SO.** `beta` is estimated on seed-FAMILY columns and
broadcast to arms, but `card._edge_coupling` is a **single-ARM** partial Spearman. So on any multi-arm
edge the card asserts a FAMILY coupling weight next to an ARM coupling correlation, and nothing marks the
mismatch. This emits `coupling_fam` — the same partial Spearman on the family's **true RPM pool**
(`log2(1+Σ(2^x−1))`, never a mean) — so both rungs of the coupling readout sit side by side.

⭐ **(2) WITHIN-CELL ARM IDENTITY.** `identity` is Shapley/LMG credit per FAMILY. Nothing answers *within*
a family: given that this seed matters, WHICH ARM carries it? `identity_arm` is an exact Shapley over the
cell's arm columns on R²(arms, y|C) — exact because cells are small (max 7 arms ⇒ ≤128 subsets).

⚠⚠ **THE COHERENCE CHECK IS THE POINT, NOT THE SPLIT.** Same-seed arms share a binding site, so a
within-cell split is the most collinear question in the whole program and a confident per-arm credit
would be suspicious. MH-186 already produced an independent verdict on the same cells
(`arm_resolvable`: separable AND pays out-of-fold, true for 28.7%). **If the two readouts are measuring
anything real, identity concentration must be HIGHER where the arms are resolvable and near-EQUAL where
they are not.** That is an internal-consistency test between two independently-computed quantities, and
it is cheaper and sharper than inventing a permutation null for a collinear design.

PRE-REGISTERED (axiom 1a): concentration (max within-cell share) is higher in `arm_resolvable` cells than
in the rest. If it is NOT, at least one of the two readouts is noise and neither should be shipped.

⛔ Family β and the arm-level `coupling_*` both remain the CANONICAL card defaults. These columns exist so
the rung is VISIBLE and checkable, not to replace either.
"""
from __future__ import annotations

import argparse
import os
from itertools import combinations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "rung_parity.tsv"
MAX_ARMS = 7          # exact Shapley: 2^7 = 128 subsets, trivial


def _r2(y, Z):
    if Z.shape[1] == 0:
        return 0.0
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    r = y - Z @ b
    sst = float(((y - y.mean()) ** 2).sum())
    return max(0.0, 1.0 - float((r ** 2).sum()) / sst) if sst > 0 else 0.0


def _shapley(y, Z):
    """Exact Shapley over columns of Z on R² — the LMG decomposition, collinearity-fair."""
    p = Z.shape[1]
    idx = list(range(p))
    from math import factorial
    phi = np.zeros(p)
    for i in idx:
        rest = [j for j in idx if j != i]
        for k in range(len(rest) + 1):
            wgt = factorial(k) * factorial(p - k - 1) / factorial(p)
            for S in combinations(rest, k):
                phi[i] += wgt * (_r2(y, Z[:, list(S) + [i]]) - _r2(y, Z[:, list(S)]))
    return np.maximum(phi, 0.0)


def _one(gene: str):
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM

    try:
        Y, X, Cm, w = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return None
    if X.shape[1] < 1:
        return None
    fam = FAM.family_of(pd.Index(X.columns))
    lab = pd.Series([str(fam.get(a)) for a in X.columns], index=X.columns)
    Cmat = Cm.to_numpy(float)
    lc = LinearRegression().fit(Cmat, Y.to_numpy(float))
    yr = Y.to_numpy(float) - lc.predict(Cmat)

    rows = []
    for cell, members in lab.groupby(lab).groups.items():
        members = list(members)
        # ── (1) FAMILY-rung coupling: the cell's TRUE RPM pool vs the target
        lin = np.nan_to_num(np.power(2.0, X[members].to_numpy(float)) - 1.0, nan=0.0)
        pool = np.log2(1.0 + lin.sum(axis=1))
        lp = LinearRegression().fit(Cmat, pool)
        cf = stats.spearmanr(pool - lp.predict(Cmat), yr).correlation
        # ── (2) WITHIN-CELL arm identity (multi-arm cells only)
        share = {}
        if 1 < len(members) <= MAX_ARMS:
            Z = X[members].to_numpy(float)
            Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
            Zr = Z - Cmat @ np.linalg.lstsq(Cmat, Z, rcond=None)[0]
            phi = _shapley(yr, Zr)
            tot = phi.sum()
            share = {m: (float(v / tot) if tot > 1e-12 else np.nan) for m, v in zip(members, phi)}
        for m in members:
            rows.append({"gene": gene, "arm": m, "n_arm_cell": len(members),
                         "coupling_fam": float(cf) if cf == cf else np.nan,
                         "identity_arm": share.get(m, np.nan)})
    return rows


def run(genes=None, workers: int = 8) -> pd.DataFrame:
    from multiprocessing import Pool

    from mirna_hallmark.learned.evidence import ledger as LG
    genes = genes or sorted(LG.pooled_he_edges().gene.unique())
    print(f"[rung_parity] {len(genes):,} genes · {workers} workers")
    with Pool(workers) as p:
        res = [r for r in p.imap_unordered(_one, genes, chunksize=4) if r]
    R = pd.DataFrame([x for rs in res for x in rs])
    R.to_csv(DEST, sep="\t", index=False)
    print(f"-> {DEST}")
    return R


def report(R: pd.DataFrame) -> None:
    card = pd.read_csv(OUT / "realization/edge_card.tsv", sep="\t", low_memory=False,
                       usecols=lambda c: c in ("gene", "arm", "coupling_tum", "arm_resolvable",
                                               "arm_sep_z", "n_arm_in_cell"))
    d = R.merge(card, on=["gene", "arm"], how="left")

    print("\n=== (1) THE COUPLING RUNG: arm-level (card default) vs FAMILY pool ===")
    m = d.dropna(subset=["coupling_fam", "coupling_tum"])
    multi = m[m.n_arm_cell > 1]
    print(f"  all edges   n={len(m):,}   arm {m.coupling_tum.mean():+.4f}  family {m.coupling_fam.mean():+.4f}"
          f"   spearman {stats.spearmanr(m.coupling_tum, m.coupling_fam).correlation:.4f}")
    print(f"  MULTI-arm   n={len(multi):,}   arm {multi.coupling_tum.mean():+.4f}  "
          f"family {multi.coupling_fam.mean():+.4f}"
          f"   spearman {stats.spearmanr(multi.coupling_tum, multi.coupling_fam).correlation:.4f}")
    print("  ⇒ where they diverge, the card was asserting a FAMILY weight beside an ARM correlation.")

    print("\n=== (2) WITHIN-CELL ARM IDENTITY — the pre-registered coherence check ===")
    q = d.dropna(subset=["identity_arm"])
    cell = q.groupby(["gene", "arm"]).first().reset_index()
    top = q.groupby(["gene"]).identity_arm.max()
    res = q.dropna(subset=["arm_resolvable"])
    if len(res):
        a = res[res.arm_resolvable == True].groupby("gene").identity_arm.max()
        b = res[res.arm_resolvable == False].groupby("gene").identity_arm.max()
        if len(a) > 5 and len(b) > 5:
            u = stats.mannwhitneyu(a, b, alternative="greater")
            print(f"  max within-cell share — RESOLVABLE cells {a.mean():.3f} (n={len(a)})  "
                  f"vs NOT {b.mean():.3f} (n={len(b)})   MWU p={u.pvalue:.3g}")
            print(f"  {'✅ COHERENT' if u.pvalue < 0.05 else '⚠ NOT COHERENT'} — "
                  "concentration should be higher exactly where MH-186 says the arms are separable.")
    print(f"  overall max share {top.mean():.3f} (equal split for a 2-arm cell = 0.500)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    report(run(workers=a.workers))
