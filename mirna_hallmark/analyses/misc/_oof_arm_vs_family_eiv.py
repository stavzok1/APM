"""MH-ADHOC: does the ARM-design's out-of-sample advantage over the FAMILY-collapsed design
survive in MEASUREMENT-HOMOGENEOUS families?  (unequal biology vs errors-in-variables reweighting)

5-fold OOF over patients; per gene: partial Spearman(pred, Y | C) on held-out predictions,
delta = rho_arm - rho_family  (negative = arm better, repression is negative).
"""
from __future__ import annotations
import sys, warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from scipy.stats import spearmanr, rankdata

warnings.filterwarnings("ignore")

from mirna_hallmark.learned import data as LD, families as FAM, regression as REG

SEED = 0
NFOLD = 5


def _resid(a: np.ndarray, C: np.ndarray) -> np.ndarray:
    D = np.column_stack([np.ones(len(a)), C])
    beta, *_ = np.linalg.lstsq(D, a, rcond=None)
    return a - D @ beta


def partial_spearman(pred: np.ndarray, y: np.ndarray, C: np.ndarray):
    ok = np.isfinite(pred) & np.isfinite(y)
    if ok.sum() < 30 or np.nanstd(pred[ok]) == 0:
        return np.nan
    Cr = np.column_stack([rankdata(C[ok, j]) for j in range(C.shape[1])])
    keep = Cr.std(0) > 0
    Cr = Cr[:, keep]
    rp = _resid(rankdata(pred[ok]).astype(float), Cr)
    ry = _resid(rankdata(y[ok]).astype(float), Cr)
    if rp.std() == 0 or ry.std() == 0:
        return np.nan
    return float(spearmanr(rp, ry).statistic)


def run_gene(gene: str, det: pd.Series, nmap: dict):
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    fam = FAM.family_of(X.columns)
    Xf, wf, members = FAM.collapse_by_family(X, w, fam)
    multi = {f: m for f, m in members.items() if len(m) > 1}
    if not multi:
        return None
    Cn = C.apply(pd.to_numeric, errors="coerce")
    Cn = Cn.loc[:, Cn.notna().all() & (Cn.std() > 0)]
    Cm = Cn.to_numpy(float)
    yv = Y.to_numpy(float)
    n = len(Y)
    pa = np.full(n, np.nan)
    pf = np.full(n, np.nan)
    kf = KFold(n_splits=NFOLD, shuffle=True, random_state=SEED)
    for tr, te in kf.split(np.arange(n)):
        # ⭐ SWITCHED off the RETIRED adaptive lasso -> canonical Gibbs drop-in (MH-184, 2026-08-01).
        # This call's M is the REPORTED quantity (ESTIMAND class in `eval/_retired_lasso_audit.py`),
        # so the retired estimator was a real defect here. ⚠ Any persisted output is STALE until re-run.
        Ma = REG.fit_gene_bayes(Y.iloc[tr], X.iloc[tr], Cn.iloc[tr], w)
        Mf = REG.fit_gene_bayes(Y.iloc[tr], Xf.iloc[tr], Cn.iloc[tr], wf)
        pa[te] = REG.aggregate(X.iloc[te], Ma)
        pf[te] = REG.aggregate(Xf.iloc[te], Mf)
    ra = partial_spearman(pa, yv, Cm)
    rf = partial_spearman(pf, yv, Cm)

    # --- measurement-quality homogeneity within the multi-arm families ---
    det_min, det_spread, expr_spread, var_spread = 1.0, 0.0, 0.0, 0.0
    for f, arms in multi.items():
        xs = [nmap.get(a) or nmap.get(str(a).lower()) or a for a in arms]
        dv = det.reindex(xs).to_numpy(float)
        det_min = min(det_min, float(np.nanmin(dv)))
        det_spread = max(det_spread, float(np.nanmax(dv) - np.nanmin(dv)))
        ex = X[arms].mean(0).to_numpy(float)
        expr_spread = max(expr_spread, float(ex.max() - ex.min()))
        vr = X[arms].std(0).to_numpy(float)
        var_spread = max(var_spread, float(vr.max() - vr.min()))
    return dict(gene=gene, n=n, n_arms=X.shape[1], n_fam=Xf.shape[1],
                n_multi_fam=len(multi), max_fam_size=max(len(m) for m in multi.values()),
                rho_arm=ra, rho_fam=rf, delta=(ra - rf) if np.isfinite(ra) and np.isfinite(rf) else np.nan,
                det_min=det_min, det_spread=det_spread, expr_spread=expr_spread, var_spread=var_spread)


def main():
    r = pd.read_csv("mirna_hallmark/output/learned/readouts_arm_edges.tsv", sep="\t")
    genes = sorted(r.loc[r.family_size > 1, "gene"].unique())
    print(f"[genes with >=1 multi-arm family in readouts] {len(genes)}", flush=True)
    Xall = LD._load()["X"]
    det = pd.Series(np.isfinite(Xall.to_numpy(float)).mean(1), index=Xall.index)
    nmap = LD._arm_name_map(Xall)
    rows = []
    for i, g in enumerate(genes):
        try:
            out = run_gene(g, det, nmap)
        except Exception as e:
            print(f"  SKIP {g}: {type(e).__name__}: {e}", flush=True)
            continue
        if out:
            rows.append(out)
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(genes)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv("/sci/labs/michall/stavzok/APM/mirna_hallmark/output/learned/_oof_arm_vs_family_eiv.tsv",
              sep="\t", index=False)
    print(df.shape)


if __name__ == "__main__":
    main()
