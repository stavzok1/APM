"""MH-118 — DOES THE §8 SEED-FAMILY COLLAPSE THROW AWAY REAL SIGNAL? A head-to-head, OOF.

    .venv/bin/python3 -m mirna_hallmark.eval.within_family_showdown

THE §8 RATIONALE (`learned/families.py`): *"Same-seed arms are near-collinear ⇒ plain lasso keeps one member
arbitrarily, UNSTABLY across folds."* So the model collapses same-seed arms into ONE family predictor, and
the latent beta is per FAMILY. That collapse is load-bearing for all four readouts.

WHAT PROVOKED THIS (MH-117). Within a family, a direct ARM-level OOF beta predicted the arm's realized
held-out coupling at Spearman −0.72 (−0.48 after conditioning on abundance/detection/variance), and it
tracked AGO loading (the biological guide label, p=1.2e-03) — while the family beta is CONSTANT within a
family and dose-delivery share had NO within-family power (p=0.34, wrong sign).

⚠ BUT THAT WAS THE WRONG ESTIMAND. Correlating beta_arm(train) with the arm's MARGINAL held-out rho is close
to tautological: an arm the fit weights is an arm that correlates with Y. It would come out negative even if
the collapse were entirely correct. **The estimand is OOF PREDICTION of the gene's mRNA**, and §8's actual
claim is about **STABILITY**. So test exactly those two things:

  (1) PREDICTION  : OOF partial rho(M_g, Y_g | C) on held-out patients, for M built from
                    ARM-level beta vs FAMILY-collapsed beta vs a DETECTION-GATED arm-level beta.
  (2) STABILITY   : across folds, does the arm-level fit pick the SAME member? (§8's claim.)
                    Measured as fold-to-fold top-member agreement + the Spearman of beta across folds,
                    for arms vs families.

Only genes with at least one MULTI-ARM family are used — they are the only place the collapse can matter.

Output: `output/learned/within_family_showdown.tsv`.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import KFold

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR

OUT = Path("mirna_hallmark/output/learned")
FOLDS = 5
DET_MIN = 0.80          # detection-gate: keep arms measured in >=80% of patients (MH-117: detection rate is
                        # the ONE artifact proxy the arm-level beta tracks, +0.305 p=4e-3)


def _partial(y, x, Cm):
    m = np.isfinite(y) & np.isfinite(x) & np.isfinite(Cm).all(1)
    if m.sum() < 30:
        return np.nan
    A = np.column_stack([np.ones(m.sum()), Cm[m]])
    ry, rx = stats.rankdata(y[m]), stats.rankdata(x[m])
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    if ey.std() < 1e-9 or ex.std() < 1e-9:
        return np.nan
    return float((ey - ey.mean()) @ (ex - ex.mean()) / (len(ey) * ey.std() * ex.std()))


def run(folds: int = FOLDS, seed: int = 0, limit: int | None = None) -> pd.DataFrame:
    e = D.load_hallmark_edges()
    genes = sorted(e[e.high_evidence].gene.unique())
    if limit:
        genes = genes[:limit]
    Xall = LD._load()["X"]
    det_all = pd.Series(np.isfinite(Xall.to_numpy(float)).mean(1), index=Xall.index)

    rows, stab = [], []
    for i, g in enumerate(genes):
        if i % 200 == 0:
            print(f"[wf] {i}/{len(genes)}", flush=True)
        try:
            Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger")
        except Exception:
            continue
        if X is None or X.shape[1] < 2 or len(Y) < 100:
            continue
        fam = FAM.family_of(list(X.columns))
        fam = fam.reindex(X.columns)
        sizes = fam.value_counts()
        if (sizes > 1).sum() == 0:                     # no multi-arm family -> collapse cannot matter
            continue
        Xf, wf, _members = FAM.collapse_by_family(X, w, fam)

        det = det_all.reindex(X.columns).fillna(0.0)
        keep = [c for c in X.columns if det.get(c, 0) >= DET_MIN]
        Xg = X[keep] if len(keep) >= 2 else None

        yv = Y.to_numpy(float)
        Cm = C.to_numpy(float)
        n = len(yv)
        preds = {k: np.full(n, np.nan) for k in ["arm", "family", "arm_detgated"]}
        betas = {"arm": [], "family": []}
        kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
        for tr, te in kf.split(np.arange(n)):
            for key, Xd, wd in [("arm", X, w), ("family", Xf, wf),
                                ("arm_detgated", Xg, (w.reindex(keep) if Xg is not None else None))]:
                if Xd is None or Xd.shape[1] == 0:
                    continue
                try:
                    b = LR.fit_gene(Y.iloc[tr], Xd.iloc[tr], C.iloc[tr], wd)
                except Exception:
                    continue
                bv = b.reindex(Xd.columns).fillna(0.0).to_numpy(float)
                preds[key][te] = np.nan_to_num(Xd.iloc[te].to_numpy(float)) @ bv
                if key in betas:
                    betas[key].append(pd.Series(bv, index=Xd.columns))

        r = {"gene": g, "n_arms": X.shape[1], "n_fams": Xf.shape[1],
             "n_multi": int((sizes > 1).sum()), "n_kept": len(keep)}
        for key in preds:
            r[f"oof_{key}"] = _partial(yv, preds[key], Cm)
        rows.append(r)

        # --- STABILITY: does the fit pick the same member across folds?
        for key in ["arm", "family"]:
            B = betas[key]
            if len(B) < 2:
                continue
            M = pd.concat(B, axis=1)
            agree = np.mean([M.iloc[:, a].abs().idxmax() == M.iloc[:, b].abs().idxmax()
                             for a, b in combinations(range(M.shape[1]), 2)])
            rr = [stats.spearmanr(M.iloc[:, a], M.iloc[:, b]).correlation
                  for a, b in combinations(range(M.shape[1]), 2)]
            stab.append({"gene": g, "level": key, "top_agree": float(agree),
                         "beta_rank_corr": float(np.nanmean(rr))})

    R = pd.DataFrame(rows)
    S = pd.DataFrame(stab)
    R.to_csv(OUT / "within_family_showdown.tsv", sep="\t", index=False)
    S.to_csv(OUT / "within_family_stability.tsv", sep="\t", index=False)

    print(f"\n=== (1) OOF PREDICTION of held-out mRNA — {len(R)} genes with a multi-arm family ===")
    for k in ["family", "arm", "arm_detgated"]:
        v = R[f"oof_{k}"].dropna()
        print(f"   {k:14s} mean OOF partial rho = {v.mean():+.4f}   (median {v.median():+.4f}, n={len(v)})")
    d = R.dropna(subset=["oof_arm", "oof_family"])
    dd = d.oof_arm - d.oof_family
    print(f"\n   ARM minus FAMILY: mean Δ = {dd.mean():+.4f}  (negative = arm predicts BETTER)")
    print(f"     wins: arm {int((dd < 0).sum())} / family {int((dd > 0).sum())}   "
          f"Wilcoxon p = {stats.wilcoxon(d.oof_arm, d.oof_family).pvalue:.2e}")
    d2 = R.dropna(subset=["oof_arm_detgated", "oof_family"])
    if len(d2):
        dd2 = d2.oof_arm_detgated - d2.oof_family
        print(f"   ARM(det-gated) minus FAMILY: mean Δ = {dd2.mean():+.4f}   "
              f"Wilcoxon p = {stats.wilcoxon(d2.oof_arm_detgated, d2.oof_family).pvalue:.2e}")

    print(f"\n=== (2) STABILITY across folds (§8's actual claim) ===")
    for k in ["arm", "family"]:
        s = S[S.level == k]
        print(f"   {k:8s}: top-member agreement {s.top_agree.mean():.1%}   "
              f"beta rank-corr across folds {s.beta_rank_corr.mean():+.3f}   (n={len(s)} genes)")
    return R


if __name__ == "__main__":
    run()
