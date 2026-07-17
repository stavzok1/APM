"""MH-116 — THE FULL GRID: layer x cohort x unit x universe, with THREE aggregation arms.

    .venv/bin/python3 -m mirna_hallmark.eval.grid_full

AXES (user 2026-07-13)
  layer    : mRNA | protein
  cohort   : TCGA (n~1041) | CPTAC-2 prospective (n~101)   [+ BOTH CPTAC mRNA tables: linkedomics | star]
  unit     : gene (aggregate) | edge (per miRNA-arm x gene)
  universe : HE (curated high-evidence) | ORPHAN (miRTarBase-absent -> NEVER curated)

ARMS (the benchmark; only the WEIGHT differs)
  abundance : uniform w=1 (gene) / the arm's own median expression (edge)  -- "is it just how much miRNA?"
  heuristic : the production evidence score / `compute_gene_pressure`      -- EXISTS ONLY FOR HE (orphans are
              uncurated BY DEFINITION, so the heuristic cannot even be formed there -- that is the point)
  learned   : beta from `regression.fit_gene` on TCGA mRNA (the canonical core-C fit), for the HE design and
              for the ORPHAN-augmented design respectively

WHY THE ORPHAN UNIVERSE IS THE **EXISTENCE** CLAIM (user 2026-07-13). On HE edges one can object that beta
"knows which edges are curated". (In fact curation is CONSTANT within HE, so it cannot drive a WITHIN-HE
correlation -- the objection is weaker than it looks.) But the objection dies completely on ORPHAN edges:
they carry NO curation at all, so if beta -- fitted on TCGA mRNA -- still predicts an edge's realized
coupling in CPTAC PROTEIN (different patients, different molecular layer, no curation), the model is
capturing biology, not curation. That is an EXISTENCE claim in a universe the curation cannot reach.

Output: `output/learned/grid_full.tsv` (one row per cell x arm).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import data_loaders as D
from mirna_hallmark.eval.cptac_validation import (_covariates, load_cptac_layers, load_prospective_clinical)
from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import cptac_data as CD
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import regression as LR

OUT = Path("mirna_hallmark/output/learned")
MIN_N = 30


# --------------------------------------------------------------------------- #
# LEARNED beta, fitted on TCGA mRNA (the canonical core-C fit), per universe
# --------------------------------------------------------------------------- #
def fit_beta(universe: str, genes: list[str]) -> pd.DataFrame:
    """Per-(gene, arm) learned weight from the TCGA fit. universe='he' | 'orphan'."""
    rows = []
    for i, g in enumerate(genes):
        if i % 250 == 0:
            print(f"[beta:{universe}] {i}/{len(genes)}", flush=True)
        try:
            if universe == "he":
                Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger")
            else:
                Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True,
                                              orphan_source="targetscan")
        except Exception:
            continue
        if X is None or X.shape[1] == 0:
            continue
        try:
            b = LR.fit_gene(Y, X, C, w)
        except Exception:
            continue
        for arm, v in b.items():
            rows.append({"gene": g, "miRNA": arm, "beta": float(v)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# per-EDGE realized coupling in a given cohort/layer, vectorized
# --------------------------------------------------------------------------- #
def edge_rho(edges: pd.DataFrame, X: pd.DataFrame, Y: pd.DataFrame, C: pd.DataFrame) -> np.ndarray:
    s = [c for c in X.columns if c in Y.columns and c in C.index]
    cv = C.loc[s].apply(pd.to_numeric, errors="coerce")
    keep = list(cv.dropna().index)
    A = np.column_stack([np.ones(len(keep)), cv.loc[keep].to_numpy(float)])
    Xm, Ym = X[keep], Y[keep]
    xi = {m: i for i, m in enumerate(Xm.index)}
    yi = {g: i for i, g in enumerate(Ym.index)}
    XA, YA = Xm.to_numpy(float), Ym.to_numpy(float)
    out = np.full(len(edges), np.nan)
    for k, (m, g) in enumerate(zip(edges.miRNA, edges.gene)):
        if m not in xi or g not in yi:
            continue
        xv, yv = XA[xi[m]], YA[yi[g]]
        msk = ~np.isnan(xv) & ~np.isnan(yv)
        if msk.sum() < MIN_N:
            continue
        Ai = A[msk]
        rx, ry = stats.rankdata(xv[msk]), stats.rankdata(yv[msk])
        bx, *_ = np.linalg.lstsq(Ai, rx, rcond=None)
        by, *_ = np.linalg.lstsq(Ai, ry, rcond=None)
        ex, ey = rx - Ai @ bx, ry - Ai @ by
        sx, sy = ex.std(), ey.std()
        if sx < 1e-9 or sy < 1e-9:
            continue
        out[k] = float((ex - ex.mean()) @ (ey - ey.mean()) / (len(ex) * sx * sy))
    return out


def gene_rho(score: pd.DataFrame, Y: pd.DataFrame, C: pd.DataFrame) -> pd.DataFrame:
    from mirna_hallmark.eval.coupling_grid import _score
    return _score(Y, score, C, {})


def main() -> pd.DataFrame:
    e = D.load_hallmark_edges()
    HE = e[e.high_evidence][["miRNA", "gene"]].drop_duplicates()
    ev = e[e.high_evidence].groupby(["miRNA", "gene"])["evidence_score"].max()

    dos = pd.read_csv(OUT / "discovery_dossier.tsv", sep="\t")
    ORPH = dos.rename(columns={"arm": "miRNA"})[["miRNA", "gene"]].drop_duplicates()
    curated = set(zip(e.miRNA, e.gene))
    ORPH = ORPH[[(m, g) not in curated for m, g in zip(ORPH.miRNA, ORPH.gene)]]   # truly uncurated
    print(f"[grid] HE edges {len(HE)} | ORPHAN edges {len(ORPH)}")

    # --- learned beta, TCGA-fitted, one per universe
    bHE = fit_beta("he", sorted(HE.gene.unique()))
    bOR = fit_beta("orphan", sorted(ORPH.gene.unique()))
    bHE.to_csv(OUT / "grid_beta_he.tsv", sep="\t", index=False)
    bOR.to_csv(OUT / "grid_beta_orphan.tsv", sep="\t", index=False)
    print(f"[grid] beta: HE {len(bHE)} | ORPHAN {len(bOR)}")

    # --- cohorts / layers
    d = LD._load()
    Xt, Yt = d["X"], d["Y"]
    pt = [p for p in Xt.columns if p in Yt.columns and p in CF.deconv("tcga").index]
    Xc = CD.arms()
    clin = load_prospective_clinical()
    lay = load_cptac_layers("prospective")
    CELLS = [
        ("mRNA",    "TCGA",  Xt[pt], Yt[pt],           lambda c: CF.build_C("tcga", pt, composition=c)),
        ("mRNA-LO", "CPTAC", Xc,     CD.mrna("linkedomics"), lambda c: _covariates(clin, "prospective", composition=c)),
        ("mRNA-STAR", "CPTAC", Xc,   CD.mrna("star"),  lambda c: _covariates(clin, "prospective", composition=c)),
        ("protein", "CPTAC", Xc,     lay["protein_z"], lambda c: _covariates(clin, "prospective", composition=c)),
    ]

    rows = []
    for layer, cohort, X, Y, Cf in CELLS:
        abund = X.median(axis=1)                       # the arm's own abundance
        for uni, E, beta in [("HE", HE, bHE), ("ORPHAN", ORPH, bOR)]:
            ed = E.merge(beta, on=["miRNA", "gene"], how="left")
            ed["evidence"] = [ev.get((m, g), np.nan) for m, g in zip(ed.miRNA, ed.gene)]
            ed["abundance"] = ed.miRNA.map(abund)
            for cname, comp in [("core", False), ("deconv", True)]:
                C = Cf(comp)
                ed[f"rho_{cname}"] = edge_rho(ed, X, Y, C)
                sub = ed.dropna(subset=[f"rho_{cname}"])
                # ---- EDGE level: which weight predicts the realized coupling?
                for arm in ["abundance", "heuristic", "learned"]:
                    col = {"abundance": "abundance", "heuristic": "evidence", "learned": "beta"}[arm]
                    v = sub.dropna(subset=[col])
                    if len(v) < 50:
                        rows.append(dict(layer=layer, cohort=cohort, unit="edge", universe=uni, C=cname,
                                         arm=arm, n=len(v), stat=np.nan, p=np.nan))
                        continue
                    r = stats.spearmanr(v[col], v[f"rho_{cname}"])
                    rows.append(dict(layer=layer, cohort=cohort, unit="edge", universe=uni, C=cname, arm=arm,
                                     n=len(v), stat=r.correlation, p=r.pvalue))
                # ---- GENE level: aggregate the edges, then correlate with the outcome
                S = {}
                for arm in ["abundance", "learned"]:
                    w = ed.dropna(subset=["beta"]).copy()
                    if arm == "abundance":
                        w["w"] = 1.0
                    else:
                        w["w"] = w["beta"]
                    mats = {}
                    for g, s2 in w.groupby("gene"):
                        arms_ = [a for a in s2.miRNA if a in X.index]
                        if not arms_:
                            continue
                        ww = s2.set_index("miRNA")["w"].reindex(arms_).to_numpy(float)
                        mats[g] = ww @ np.nan_to_num(X.loc[arms_].to_numpy(float))
                    if mats:
                        S[arm] = pd.DataFrame(mats, index=X.columns).T
                for arm, Sm in S.items():
                    gr = gene_rho(Sm, Y, C).dropna(subset=["rho"])
                    if len(gr) < 20:
                        continue
                    rows.append(dict(layer=layer, cohort=cohort, unit="gene", universe=uni, C=cname, arm=arm,
                                     n=len(gr), stat=gr.rho.mean(),
                                     p=stats.ttest_1samp(gr.rho, 0).pvalue))
            print(f"[grid] done {layer}/{cohort}/{uni}", flush=True)

    g = pd.DataFrame(rows)
    g.to_csv(OUT / "grid_full.tsv", sep="\t", index=False)
    pd.set_option("display.width", 240)
    print("\n" + g.to_string(index=False))
    return g


if __name__ == "__main__":
    main()
