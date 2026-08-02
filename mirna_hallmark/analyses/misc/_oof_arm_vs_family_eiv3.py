"""Part 3 — the decisive discriminator.

In an abundance-skewed family the pool log2(1+sum RPM) ~= the DOMINANT member, so a noise-optimal
estimator has nothing left to down-weight (the minor member is already ~absent from the pool).
EIV can only ATTENUATE a noisy member's weight; it can never INFLATE it. So:

  DROP design = family pool with the minor member(s) DELETED (minor weight forced to 0 = the EIV limit).

  H_EIV      : rho_drop <= rho_fam  and  rho_arm ~= rho_drop   (arm wins by discarding the noisy minor)
  H_biology  : rho_arm  <  rho_drop                            (the minor member ADDS real information)

Also emits corr(X_fam, X_dominant) and the fitted arm weight on the minor member.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
from mirna_hallmark.learned import data as LD, families as FAM, regression as REG
from mirna_hallmark.analyses.misc._oof_arm_vs_family_eiv import partial_spearman

SEED, NFOLD = 0, 5


def run(gene, det, nmap):
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    fam = FAM.family_of(X.columns)
    Xf, wf, members = FAM.collapse_by_family(X, w, fam)
    multi = {f: m for f, m in members.items() if len(m) > 1}
    if not multi:
        return None
    C = C.apply(pd.to_numeric, errors="coerce")
    C = C.loc[:, C.notna().all() & (C.std() > 0)]
    lin = pd.DataFrame(np.power(2.0, X.to_numpy(float)) - 1.0, index=X.index, columns=X.columns)

    # dominant member per family; DROP design = dominant-only columns (minor deleted)
    dom, minors = {}, []
    for f, arms in members.items():
        d = lin[arms].mean(0).idxmax()
        dom[f] = d
        minors += [a for a in arms if a != d]
    Xd = X[[dom[f] for f in members]].copy()
    Xd.columns = list(members)
    wd = pd.Series({f: float(w.reindex(members[f]).max()) for f in members})

    # how close is the pool to the dominant member alone?
    fam_dom_r, minor_share, det_min, det_sp, expr_sp = [], [], 1.0, 0.0, 0.0
    for f, arms in multi.items():
        fam_dom_r.append(float(np.corrcoef(Xf[f], X[dom[f]])[0, 1]))
        tot = lin[arms].mean(0).to_numpy(float)
        minor_share.append(float(1.0 - tot.max() / max(tot.sum(), 1e-9)))
        xs = [nmap.get(a) or nmap.get(str(a).lower()) or a for a in arms]
        dv = det.reindex(xs).to_numpy(float)
        det_min = min(det_min, float(np.nanmin(dv)))
        det_sp = max(det_sp, float(np.nanmax(dv) - np.nanmin(dv)))
        ex = X[arms].mean(0).to_numpy(float)
        expr_sp = max(expr_sp, float(ex.max() - ex.min()))

    n = len(Y)
    pa, pf, pd_ = (np.full(n, np.nan) for _ in range(3))
    for tr, te in KFold(n_splits=NFOLD, shuffle=True, random_state=SEED).split(np.arange(n)):
        # ⭐ SWITCHED off the RETIRED adaptive lasso -> canonical Gibbs drop-in (MH-184, 2026-08-01).
        # This call's M is the REPORTED quantity (ESTIMAND class in `eval/_retired_lasso_audit.py`),
        # so the retired estimator was a real defect here. ⚠ Any persisted output is STALE until re-run.
        Ma = REG.fit_gene_bayes(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w)
        Mf = REG.fit_gene_bayes(Y.iloc[tr], Xf.iloc[tr], C.iloc[tr], wf)
        Md = REG.fit_gene_bayes(Y.iloc[tr], Xd.iloc[tr], C.iloc[tr], wd)
        pa[te] = REG.aggregate(X.iloc[te], Ma)
        pf[te] = REG.aggregate(Xf.iloc[te], Mf)
        pd_[te] = REG.aggregate(Xd.iloc[te], Md)
    Cm = C.to_numpy(float)
    yv = Y.to_numpy(float)
    ra = partial_spearman(pa, yv, Cm)
    rf = partial_spearman(pf, yv, Cm)
    rd = partial_spearman(pd_, yv, Cm)

    # full-data arm fit: does the model actually LOAD the minor member?
    Ma_full = REG.fit_gene_bayes(Y, X, C, w)
    mw = Ma_full.reindex([a for a in minors if a in Ma_full.index]).fillna(0.0)
    dw = Ma_full.reindex([dom[f] for f in multi]).fillna(0.0)
    minor_sel = float((mw > 0).mean()) if len(mw) else np.nan
    # per-molecule inflation: minor's fitted weight vs dominant's, within multi-arm families
    infl = []
    for f, arms in multi.items():
        d = dom[f]
        for a in arms:
            if a == d:
                continue
            infl.append(float(Ma_full.get(a, 0.0) - Ma_full.get(d, 0.0)))
    return dict(gene=gene, n=n, n_multi_fam=len(multi),
                rho_arm=ra, rho_fam=rf, rho_drop=rd,
                d_arm_fam=ra - rf, d_drop_fam=rd - rf, d_arm_drop=ra - rd,
                fam_dom_r=float(np.nanmean(fam_dom_r)), minor_share=float(np.nanmax(minor_share)),
                det_min=det_min, det_spread=det_sp, expr_spread=expr_sp,
                minor_selected=minor_sel, minor_w_mean=float(mw.mean()) if len(mw) else np.nan,
                dom_w_mean=float(dw.mean()) if len(dw) else np.nan,
                minor_minus_dom_w=float(np.nanmean(infl)) if infl else np.nan)


def main():
    r = pd.read_csv("mirna_hallmark/output/learned/readouts_arm_edges.tsv", sep="\t")
    genes = sorted(r.loc[r.family_size > 1, "gene"].unique())
    Xall = LD._load()["X"]
    det = pd.Series(np.isfinite(Xall.to_numpy(float)).mean(1), index=Xall.index)
    nmap = LD._arm_name_map(Xall)
    rows = []
    for i, g in enumerate(genes):
        try:
            o = run(g, det, nmap)
        except Exception as e:
            print(f"SKIP {g}: {e}", flush=True); continue
        if o:
            rows.append(o)
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(genes)}", flush=True)
    pd.DataFrame(rows).to_csv(
        "/sci/labs/michall/stavzok/APM/mirna_hallmark/output/learned/_oof_arm_vs_family_eiv3.tsv",
        sep="\t", index=False)
    print("done", len(rows))


if __name__ == "__main__":
    main()
