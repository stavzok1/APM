"""Part 2: richer covariates (opportunity + measurement heterogeneity) + two calibration SIMULATIONS.

sim1 (doctrine TRUE, no EIV): Y_sim = -(X_fam @ M_fam) + C.g + permuted resid, generated from the
      MEASURED family pool. Family design is then EXACTLY correct; free arm weights can only add variance.
      -> expected delta_sim >= 0.  Measures estimator/evaluation bias.

sim2 (doctrine TRUE, EIV LIVE): the TRUE arms are the measured arms; Y_sim generated from the pool of the
      TRUE arms with EQUAL per-molecule weight; then the *observed* minor member of each multi-arm family is
      corrupted with extra multiplicative noise. Both designs see the corrupted data. -> how negative can
      noise-optimal reweighting alone drive delta?
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from scipy.stats import rankdata, spearmanr

warnings.filterwarnings("ignore")
from mirna_hallmark.learned import data as LD, families as FAM, regression as REG
from mirna_hallmark.analyses.misc._oof_arm_vs_family_eiv import partial_spearman, _resid

SEED, NFOLD = 0, 5


def oof_delta(Y, X, Xf, C, w, wf):
    n = len(Y)
    pa, pf = np.full(n, np.nan), np.full(n, np.nan)
    for tr, te in KFold(n_splits=NFOLD, shuffle=True, random_state=SEED).split(np.arange(n)):
        # ⭐ SWITCHED off the RETIRED adaptive lasso -> canonical Gibbs drop-in (MH-184, 2026-08-01).
        # This call's M is the REPORTED quantity (ESTIMAND class in `eval/_retired_lasso_audit.py`),
        # so the retired estimator was a real defect here. ⚠ Any persisted output is STALE until re-run.
        Ma = REG.fit_gene_bayes(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w)
        Mf = REG.fit_gene_bayes(Y.iloc[tr], Xf.iloc[tr], C.iloc[tr], wf)
        pa[te] = REG.aggregate(X.iloc[te], Ma)
        pf[te] = REG.aggregate(Xf.iloc[te], Mf)
    Cm = C.to_numpy(float)
    ra = partial_spearman(pa, Y.to_numpy(float), Cm)
    rf = partial_spearman(pf, Y.to_numpy(float), Cm)
    return ra, rf


def pool(lin_df, members):
    return pd.DataFrame({f: np.log2(1.0 + lin_df[m].sum(axis=1)) for f, m in members.items()},
                        index=lin_df.index)


def run(gene, det, nmap, rng):
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    fam = FAM.family_of(X.columns)
    Xf, wf, members = FAM.collapse_by_family(X, w, fam)
    multi = {f: m for f, m in members.items() if len(m) > 1}
    if not multi:
        return None
    C = C.apply(pd.to_numeric, errors="coerce")
    C = C.loc[:, C.notna().all() & (C.std() > 0)]

    # ---- covariates -------------------------------------------------------
    extra_dof = sum(len(m) - 1 for m in multi.values())
    pair_r, minor_share, det_min, det_sp, expr_sp = [], [], 1.0, 0.0, 0.0
    lin = pd.DataFrame(np.power(2.0, X.to_numpy(float)) - 1.0, index=X.index, columns=X.columns)
    for f, arms in multi.items():
        xs = [nmap.get(a) or nmap.get(str(a).lower()) or a for a in arms]
        dv = det.reindex(xs).to_numpy(float)
        det_min = min(det_min, float(np.nanmin(dv)))
        det_sp = max(det_sp, float(np.nanmax(dv) - np.nanmin(dv)))
        ex = X[arms].mean(0).to_numpy(float)
        expr_sp = max(expr_sp, float(ex.max() - ex.min()))
        cm = X[arms].corr().to_numpy(float)
        pair_r.append(float(np.nanmean(cm[np.triu_indices(len(arms), 1)])))
        tot = lin[arms].mean(0).to_numpy(float)
        minor_share.append(float(1.0 - tot.max() / max(tot.sum(), 1e-9)))   # linear-RPM share NOT in the top member

    ra, rf = oof_delta(Y, X, Xf, C, w, wf)

    # ---- sim1: DGP = measured family pool (doctrine true, no EIV) ----------
    Mf_full = REG.fit_gene_bayes(Y, Xf, C, wf)
    fit = -REG.aggregate(Xf, Mf_full)
    Cm = C.to_numpy(float)
    D = np.column_stack([np.ones(len(Y)), Cm])
    b, *_ = np.linalg.lstsq(D, Y.to_numpy(float) - fit, rcond=None)
    mu = fit + D @ b
    resid = Y.to_numpy(float) - mu
    d1 = []
    for _ in range(3):
        ys = pd.Series(mu + rng.permutation(resid), index=Y.index)
        a1, f1 = oof_delta(ys, X, Xf, C, w, wf)
        d1.append(a1 - f1)

    # ---- sim2: DGP = pool of TRUE arms, minor member then CORRUPTED (EIV) --
    d2 = {}
    for sd in (0.5, 1.0, 2.0):
        dd = []
        for _ in range(3):
            Xc = X.copy()
            linc = lin.copy()
            for f, arms in multi.items():
                tot = lin[arms].mean(0)
                minor = [a for a in arms if a != tot.idxmax()]
                for a in minor:                      # multiplicative log-normal noise on the minor member
                    linc[a] = lin[a].to_numpy(float) * np.exp(rng.normal(0, sd, len(lin)))
            for a in Xc.columns:
                Xc[a] = np.log2(1.0 + linc[a])
            Xf_c = pool(linc, members)               # observed pool (built from corrupted arms)
            # truth uses the CLEAN pool + equal per-molecule weights (= the doctrine)
            a2, f2 = oof_delta(pd.Series(mu + rng.permutation(resid), index=Y.index),
                               Xc, Xf_c, C, w, wf)
            dd.append(a2 - f2)
        d2[sd] = float(np.nanmean(dd))

    return dict(gene=gene, n=len(Y), n_arms=X.shape[1], n_multi_fam=len(multi), extra_dof=extra_dof,
                rho_arm=ra, rho_fam=rf, delta=ra - rf,
                det_min=det_min, det_spread=det_sp, expr_spread=expr_sp,
                pair_r=float(np.nanmean(pair_r)), minor_share=float(np.nanmax(minor_share)),
                delta_sim1=float(np.nanmean(d1)),
                delta_sim2_sd05=d2[0.5], delta_sim2_sd10=d2[1.0], delta_sim2_sd20=d2[2.0])


def main():
    r = pd.read_csv("mirna_hallmark/output/learned/readouts_arm_edges.tsv", sep="\t")
    genes = sorted(r.loc[r.family_size > 1, "gene"].unique())
    Xall = LD._load()["X"]
    det = pd.Series(np.isfinite(Xall.to_numpy(float)).mean(1), index=Xall.index)
    nmap = LD._arm_name_map(Xall)
    rng = np.random.default_rng(1)
    rows = []
    for i, g in enumerate(genes):
        try:
            o = run(g, det, nmap, rng)
        except Exception as e:
            print(f"SKIP {g}: {e}", flush=True); continue
        if o:
            rows.append(o)
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(genes)}", flush=True)
    pd.DataFrame(rows).to_csv(
        "/sci/labs/michall/stavzok/APM/mirna_hallmark/output/learned/_oof_arm_vs_family_eiv2.tsv",
        sep="\t", index=False)
    print("done", len(rows))


if __name__ == "__main__":
    main()
