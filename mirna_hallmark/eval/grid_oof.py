"""MH-117 — the NON-CIRCULAR TCGA cell: OUT-OF-FOLD (OOF) learned weights, all 1,041 patients.

    .venv/bin/python3 -m mirna_hallmark.eval.grid_oof

THE PROBLEM (user 2026-07-13). In `grid_full` the TCGA cells score a beta that was FITTED ON TCGA mRNA —
against TCGA mRNA. The edge-level -0.712 is therefore an IN-SAMPLE upper bound, not evidence. CPTAC is
genuinely out-of-sample (beta never saw those patients), which is why its numbers are so much smaller.

THE FIX. K-fold over PATIENTS. For each fold, fit beta on the training patients only, then score the
HELD-OUT patients with it. Every one of the 1,041 TCGA patients gets a prediction from a model that never
saw them, and the full n is retained (no subsampling).

  GENE level : M_g(s) = X[s] @ beta_train(fold containing s)  for every s  ->  partial rho(M_g, Y_g | C) on all n
  EDGE level : per fold, beta from TRAIN; the edge's realized rho from the TEST patients ONLY;
               Spearman(beta_train, rho_test) per fold, averaged over folds. Nothing is scored on the
               patients that fitted it.

ARMS: abundance (uniform w=1) vs learned (OOF beta). The heuristic evidence weight is fixed a priori (it is
not fitted), so it has no in-sample problem and is carried over from `grid_full` unchanged.

Output: `output/learned/grid_oof.tsv`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import KFold

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import regression as LR

OUT = Path("mirna_hallmark/output/learned")
FOLDS = 5
MIN_N = 30


def _partial(y: np.ndarray, x: np.ndarray, Cm: np.ndarray) -> float:
    m = np.isfinite(y) & np.isfinite(x) & np.isfinite(Cm).all(1)
    if m.sum() < MIN_N:
        return np.nan
    A = np.column_stack([np.ones(m.sum()), Cm[m]])
    ry, rx = stats.rankdata(y[m]), stats.rankdata(x[m])
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    ey, ex = ry - A @ by, rx - A @ bx
    if ey.std() < 1e-9 or ex.std() < 1e-9:
        return np.nan
    return float((ey - ey.mean()) @ (ex - ex.mean()) / (len(ey) * ey.std() * ex.std()))


def run(universe: str = "he", folds: int = FOLDS, seed: int = 0) -> pd.DataFrame:
    e = D.load_hallmark_edges()
    if universe == "he":
        genes = sorted(e[e.high_evidence].gene.unique())
        kw = dict(w_prior_source="ledger")
    else:
        dos = pd.read_csv(OUT / "discovery_dossier.tsv", sep="\t")
        genes = sorted(dos.gene.unique())
        kw = dict(w_prior_source="ledger", orphans=True, orphan_source="targetscan")

    gene_rows, edge_rows = [], []
    for i, g in enumerate(genes):
        if i % 200 == 0:
            print(f"[oof:{universe}] {i}/{len(genes)}", flush=True)
        try:
            Y, X, C, w = LD.assemble_gene(g, **kw)
        except Exception:
            continue
        if X is None or X.shape[1] == 0 or len(Y) < 100:
            continue
        yv = Y.to_numpy(float)
        Xv = X.to_numpy(float)
        Cm = C.to_numpy(float)
        n = len(yv)
        kf = KFold(n_splits=folds, shuffle=True, random_state=seed)

        oof_learned = np.full(n, np.nan)
        oof_abund = np.full(n, np.nan)
        fold_edge = []                                     # (arm, beta_train, rho_test) per fold
        for tr, te in kf.split(np.arange(n)):
            try:
                b = LR.fit_gene(Y.iloc[tr], X.iloc[tr], C.iloc[tr], w)
            except Exception:
                continue
            bv = b.reindex(X.columns).fillna(0.0).to_numpy(float)
            oof_learned[te] = np.nan_to_num(Xv[te]) @ bv          # scored on patients the fit never saw
            oof_abund[te] = np.nan_to_num(Xv[te]) @ np.ones(len(bv))
            for j, arm in enumerate(X.columns):                   # edge level: rho on the TEST patients only
                r = _partial(yv[te], Xv[te, j], Cm[te])
                if np.isfinite(r) and np.isfinite(bv[j]):
                    fold_edge.append((arm, float(bv[j]), r, float(np.nanmedian(Xv[te, j]))))

        for arm_name, sc in [("learned", oof_learned), ("abundance", oof_abund)]:
            r = _partial(yv, sc, Cm)
            if np.isfinite(r):
                gene_rows.append({"universe": universe, "arm": arm_name, "gene": g, "rho": r})
        for arm, bj, r, ab in fold_edge:
            edge_rows.append({"universe": universe, "gene": g, "miRNA": arm,
                              "beta_train": bj, "rho_test": r, "abundance": ab})

    G = pd.DataFrame(gene_rows)
    E = pd.DataFrame(edge_rows)
    return G, E


def main():
    rows = []
    allG, allE = [], []
    for uni in ["he", "orphan"]:
        G, E = run(uni)
        allG.append(G); allE.append(E)
        for arm in ["abundance", "learned"]:
            v = G[G.arm == arm].rho.dropna()
            rows.append(dict(unit="gene", universe=uni.upper(), arm=arm, n=len(v),
                             stat=v.mean(), p=stats.ttest_1samp(v, 0).pvalue))
        for arm, col in [("abundance", "abundance"), ("learned", "beta_train")]:
            v = E.dropna(subset=[col, "rho_test"])
            r = stats.spearmanr(v[col], v.rho_test)
            rows.append(dict(unit="edge", universe=uni.upper(), arm=arm, n=len(v),
                             stat=r.correlation, p=r.pvalue))
    pd.concat(allG).to_csv(OUT / "grid_oof_genes.tsv", sep="\t", index=False)
    pd.concat(allE).to_csv(OUT / "grid_oof_edges.tsv", sep="\t", index=False)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "grid_oof.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\n=== TCGA mRNA, OUT-OF-FOLD (non-circular), all patients ===")
    print(s.to_string(index=False))
    return s


if __name__ == "__main__":
    main()
