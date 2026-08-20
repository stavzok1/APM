"""THE ARM RUNG — an arm-level β beside the canonical family β, plus an IDENTIFIABILITY verdict.

    .venv/bin/python3 -m mirna_hallmark.learned.arm_rung [--workers 8] [--genes A,B]

STANDING DIRECTIVE (user, 2026-08-01): build BOTH an arm-level and a family-level model and propagate
BOTH to the cards. Today `readouts.py` estimates β on z-scored **SEED-FAMILY** columns and **broadcasts it
equally** to member arms — confirmed on the card: exactly one distinct β within **90.6%** of (gene,family)
cells. There is no arm-level β at all; the arm column is an inherited copy.

⭐ SCOPE — measured first, because it decides the design. The rung is only a DISTINCT question where a
family contributes >1 arm to a gene's design:
    467 / 4,964 (gene,family) cells have >1 arm   =  9.4%
    1,147 / 5,648 EDGES sit inside them           = 20.3%   <- the arm rung's real footprint
    273 / 1,420 genes touched                     = 19.2%   (median 2 arms/cell, max 7)
So this is a **20%-of-edges question, not a model rebuild**. The other 90.6% of cells are single-arm,
where the arm β *is* the family β by construction — fitting an "arm model" there renames a number.

⚠ AND THE ASSUMPTION IS DOING REAL WORK: within a multi-arm cell the member arms sit a **median 2.37 log2
apart in abundance, 59.7% of cells >4× apart**. Equal β across a 4-fold abundance gap is a substantive
choice, not a formality.

WHAT THIS EMITS (and deliberately does NOT)
-------------------------------------------
Per arm, within multi-arm cells only:  `beta_arm`, `sd_arm`, `z_arm`
Per cell:  `arm_dbeta`  (max−min β across member arms)
           `arm_sep_z`  (that difference in units of its own pooled posterior SD)
           `oof_drho`   (held-out ρ of the ARM design minus the FAMILY design — does splitting PAY?)
           `arm_resolvable`  = the arms are separable AND splitting them helps out-of-fold
⛔ **It does NOT replace the family β, and the card keeps family as the canonical default.** Within a
family the arms share a binding site, so an arm-level fit is asking *which COPY of this seed carries the
signal* — near-collinear by construction. That collinearity is the documented reason family collapse
exists and the reason the identifiability ceiling (median posterior SD/|β| ≈ 0.8) is what it is. The
honest deliverable is therefore **a verdict about whether the arms can be told apart at all**, not a
confident per-arm β.

⚠ JUDGED OUT-OF-FOLD, NOT IN-SAMPLE (MH-179). An arm design has strictly MORE columns than its family
collapse, so it will *always* look better in-sample. Only a held-out comparison can say whether the split
buys anything. Per-fold predictions are STANDARDISED before pooling (MH-181): β is refit per fold and
`Z·β` is an arbitrary-scale index, so pooling raw folds is not scale-invariant.
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
DEST = OUT / "arm_rung.tsv"
N_FOLDS = 5
GIBBS_ITER, GIBBS_BURN = 400, 150
SEP_Z = 2.0          # |Δβ| must exceed this many pooled SDs for the arms to be called separable


def _z(v):
    s = np.std(v)
    return (v - np.mean(v)) / s if s > 1e-12 else v * 0.0


def _one(gene: str):
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned import spike_slab as SS

    try:
        Y, X, Cm, w = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return None
    if X.shape[1] < 2:
        return None
    fam = FAM.family_of(pd.Index(X.columns))
    lab = pd.Series([str(fam.get(a)) for a in X.columns], index=X.columns)
    if lab.value_counts().max() < 2:          # no multi-arm cell -> the rung is not a distinct question
        return None
    Xf, wf, _ = FAM.collapse_by_family(X, w, fam)
    y, Cmat = Y.to_numpy(float), Cm.to_numpy(float)
    n = len(y)
    if n < 100:
        return None

    def _fit(mat):
        v = mat.to_numpy(float)
        mu, sd = v.mean(0), v.std(0)
        Z = (v - mu) / (sd + 1e-9); Z[:, sd < 0.1] = 0.0
        lc = LinearRegression().fit(Cmat, y)
        b, s, _ = SS._gibbs_posterior(Z, -(y - lc.predict(Cmat)), np.ones(Z.shape[1]),
                                      n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
        return b, s

    try:
        b_arm, s_arm = _fit(X)
    except Exception:
        return None

    # ── OOF: does SPLITTING a family into its arms buy anything held-out? (MH-179/181)
    fold = np.random.default_rng(0).permutation(np.arange(n) % N_FOLDS)
    A, F, Yo = [], [], []
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        if te.sum() < 20 or tr.sum() < 50:
            continue
        lc = LinearRegression().fit(Cmat[tr], y[tr])
        yr = -(y[tr] - lc.predict(Cmat[tr]))
        preds = []
        for mat in (X, Xf):
            v = mat.to_numpy(float)
            mu, sd = v[tr].mean(0), v[tr].std(0)
            Ztr = (v[tr] - mu) / (sd + 1e-9); Ztr[:, sd < 0.1] = 0.0
            Zte = (v[te] - mu) / (sd + 1e-9); Zte[:, sd < 0.1] = 0.0
            try:
                bb, _, _ = SS._gibbs_posterior(Ztr, yr, np.ones(Ztr.shape[1]),
                                               n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
            except Exception:
                return None
            preds.append(_z(Zte @ bb))        # ⚠ standardise per fold before pooling (MH-181)
        A.append(preds[0]); F.append(preds[1])
        Yo.append(y[te] - lc.predict(Cmat[te]))
    if not A:
        return None
    a, f, yo = np.concatenate(A), np.concatenate(F), np.concatenate(Yo)
    r_arm = stats.spearmanr(a, yo).correlation if np.std(a) > 1e-9 else np.nan
    r_fam = stats.spearmanr(f, yo).correlation if np.std(f) > 1e-9 else np.nan
    dr = (r_arm - r_fam) if np.isfinite(r_arm) and np.isfinite(r_fam) else np.nan

    rows = []
    B = pd.Series(b_arm, index=X.columns); S = pd.Series(s_arm, index=X.columns)
    for cell, members in lab.groupby(lab).groups.items():
        members = list(members)
        if len(members) < 2:
            continue
        bs, ss = B[members].to_numpy(float), S[members].to_numpy(float)
        dbeta = float(bs.max() - bs.min())
        pooled = float(np.sqrt((ss ** 2).sum() / len(ss)))
        sep = dbeta / pooled if pooled > 1e-12 else np.nan
        for m in members:
            rows.append({"gene": gene, "arm": m, "seed_family": cell, "n_arm_in_cell": len(members),
                         "beta_arm": float(B[m]), "sd_arm": float(S[m]),
                         "z_arm": float(B[m] / S[m]) if S[m] > 1e-12 else np.nan,
                         "arm_dbeta": dbeta, "arm_sep_z": sep,
                         "oof_rho_arm": r_arm, "oof_rho_fam": r_fam, "oof_drho": dr,
                         # ⭐ BOTH conditions: distinguishable AND the split pays out-of-fold
                         "cell_arms_resolvable": bool(np.isfinite(sep) and sep > SEP_Z
                                                and np.isfinite(dr) and dr < 0)})
    return rows


def run(genes=None, workers: int = 8) -> pd.DataFrame:
    from multiprocessing import Pool

    from mirna_hallmark.learned.evidence import ledger as LG
    genes = genes or sorted(LG.pooled_he_edges().gene.unique())
    print(f"[arm_rung] scanning {len(genes):,} genes for multi-arm cells · {workers} workers")
    with Pool(workers) as p:
        res = [r for r in p.imap_unordered(_one, genes, chunksize=4) if r]
    R = pd.DataFrame([x for rs in res for x in rs])
    R.to_csv(DEST, sep="\t", index=False)
    print(f"-> {DEST}")
    return R


def report(R: pd.DataFrame) -> None:
    cells = R.drop_duplicates(["gene", "seed_family"])
    print(f"\n=== THE ARM RUNG — {len(R):,} arm-edges in {len(cells):,} multi-arm cells "
          f"({R.gene.nunique():,} genes) ===\n")
    print("--- can the same-seed arms be TOLD APART? (|Δβ| vs its pooled posterior SD) ---")
    print(f"  arm_sep_z: median {cells.arm_sep_z.median():.2f}  "
          f"> {SEP_Z} in {(cells.arm_sep_z > SEP_Z).mean():.1%} of cells")
    print("\n--- does SPLITTING the family into arms PAY out-of-fold? ---")
    c = cells.dropna(subset=["oof_drho"])
    print(f"  OOF ρ  arm {c.oof_rho_arm.mean():+.4f}   family {c.oof_rho_fam.mean():+.4f}   "
          f"Δ {c.oof_drho.mean():+.4f}")
    print(f"  the arm design WINS out-of-fold in {(c.oof_drho < 0).mean():.1%} of genes  "
          f"(Wilcoxon p={stats.wilcoxon(c.oof_drho).pvalue:.3g})")
    print(f"\n--- ⭐ VERDICT: `arm_resolvable` (separable AND pays OOF) ---")
    print(f"  {int(cells.arm_resolvable.sum()):,} / {len(cells):,} cells = "
          f"{cells.arm_resolvable.mean():.1%}   covering "
          f"{int(R[R.arm_resolvable].shape[0]):,} of {len(R):,} arm-edges")
    print("\n  ⚠ family β remains the CANONICAL default on the card; these columns FLAG where the")
    print("    broadcast is an assumption you can check, they do not replace it.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--genes")
    a = ap.parse_args()
    report(run(a.genes.split(",") if a.genes else None, a.workers))
