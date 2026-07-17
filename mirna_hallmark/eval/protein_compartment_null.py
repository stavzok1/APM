"""MH-114 — the SHUFFLED NULL that proves the CPTAC bulk-protein miRNA signal is compartment ARITHMETIC.

    .venv/bin/python3 -m mirna_hallmark.eval.protein_compartment_null

THE PROBLEM. Bulk protein and bulk miRNA are both compartment-weighted averages over the cell types in the
slice. So a miRNA that merely *lives in a different cell type* than a target protein will anti-correlate with
it across patients — with NO cell-autonomous repression anywhere. Conditioning on cell composition removes
this, but conditioning could equally be an OVER-CONTROL that destroys real biology (a miRNA acting through a
cell-STATE shift PRODUCES a composition change: `LEARNED_MODEL_METHODS §1`, axiom 2a). Cross-cohort
reproducibility CANNOT arbitrate: both CPTAC cohorts are bulk breast tumour with the same CAF-rich stroma, so
a shared confound REPLICATES (measured: adjustment *lowers* prospective<->TCGA-105 concordance, +0.263 ->
+0.123 — because the ARTIFACT is what was reproducing).

THE ARBITER. Two facts that point in opposite directions:
  * compartment ARITHMETIC is BLIND to 3'UTR sites   — a fake edge is inflated exactly like a real one;
  * true REPRESSION is BLIND to compartment          — it does not care where the miRNA lives.
So score REAL curated HE edges against DEGREE-PRESERVING SHUFFLED non-edges, STRATIFIED by
    orientation = sign(corr(miRNA, CAF)) x sign(corr(protein, CAF))
  OPPOSITE-compartment: bulk mixing manufactures a NEGATIVE rho (the artifact ADDS to apparent repression)
  SAME-compartment    : bulk mixing manufactures a POSITIVE rho (the artifact MASKS repression)

  artifact  => a LARGE orientation gradient that is THE SAME in shuffled non-edges, killed by adjustment
  biology   => no gradient in shuffled non-edges; the real-vs-shuffled gap is compartment-BLIND

⚠ THE TRAP THIS MODULE EXISTS TO AVOID: a shuffled null compared on the MEAN is USELESS here. A random pair
is equally likely to be same- or opposite-compartment, so the artifact is SIGN-SYMMETRIC and CANCELS in the
mean (measured: shuffled mean rho = +0.0001). **A shuffled null must be STRATIFIED BY THE CONFOUND'S OWN
AXIS.** An unstratified mean-based null said "adjustment removes biology" — the opposite of the truth.

MEASURED (prospective, 4,130 HE edges vs 10,032 shuffled non-edges):
    REAL     OPPOSITE -0.0735 | SAME +0.0537  => gradient -0.1272  (p=1.3e-96)
    SHUFFLED OPPOSITE -0.0636 | SAME +0.0654  => gradient -0.1290  (p=8.5e-236)
  The gradients are IDENTICAL: pairs that CANNOT repress produce the same bulk anti-correlation as curated
  ones. Composition adjustment removes 90%/89% of it. ==> the unadjusted signal IS compartment arithmetic.
  TRUE effect = the compartment-BLIND real-minus-shuffled contrast: -0.0099 (opposite) / -0.0117 (same)
  ~= -0.011, the SAME in both orientations, SURVIVING adjustment at ~-0.0065 (MW p=2.0e-03).

Outputs (`output/learned/`): `protein_compartment_null_edges.tsv`, `protein_compartment_null_summary.tsv`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cptac import cptac_orphan_discovery as OD
from mirna_hallmark.eval.cptac_validation import (
    _covariates,
    load_cptac_layers,
    load_prospective_clinical,
    load_prospective_mirna_arms,
)
from mirna_hallmark.learned import confounders as CF

OUT = Path("mirna_hallmark/output/learned")
N_SHUFFLE = 3


def _fast_rho(edges: pd.DataFrame, arms: pd.DataFrame, prot: pd.DataFrame,
              cov: pd.DataFrame) -> np.ndarray:
    """Per-edge partial Spearman, vectorized. Bit-equivalent to the per-edge estimator in
    `cptac_orphan_discovery` (verified max|diff| = 5.6e-17) but ~1000x faster: the arm matrix has NO NaN, so
    an edge's valid-sample mask depends only on the GENE -> group by gene and do one matmul per gene."""
    ar = OD._rank_matrix(OD._residualize(arms, cov))
    pr = OD._rank_matrix(OD._residualize(prot, cov))
    cols = [c for c in ar.columns if c in pr.columns]
    A, P = ar[cols].to_numpy(float), pr[cols].to_numpy(float)
    ai = {m: i for i, m in enumerate(ar.index)}
    pi = {g: i for i, g in enumerate(pr.index)}
    pos = {k: i for i, k in enumerate(edges.index)}
    out = np.full(len(edges), np.nan)
    for g, idx in edges.groupby("gene").groups.items():
        if g not in pi:
            continue
        pv = P[pi[g]]
        msk = ~np.isnan(pv)
        rows = [(k, ai[m]) for k, m in zip(idx, edges.loc[idx, "miRNA"]) if m in ai]
        if not rows:
            continue
        ks, mi = zip(*rows)
        sub, y = A[list(mi)][:, msk], pv[msk]
        good = ~np.isnan(sub).any(0)
        sub, y = sub[:, good], y[good]
        if y.size < OD.MIN_N:
            continue
        yc = y - y.mean()
        xc = sub - sub.mean(1, keepdims=True)
        den = np.sqrt((xc ** 2).sum(1) * (yc ** 2).sum())
        r = np.where(den > 0, (xc @ yc) / np.where(den > 0, den, 1), np.nan)
        for k, val in zip(ks, r):
            out[pos[k]] = val
    return out


def _caf_loading(M: pd.DataFrame, caf: pd.Series) -> pd.Series:
    """Spearman(row, CAF fraction) — each feature's position on the stromal<->epithelial axis."""
    s = [c for c in M.columns if c in caf.index]
    x = caf.loc[s].to_numpy(float)
    V = M[s].to_numpy(float)
    out = np.full(V.shape[0], np.nan)
    for i in range(V.shape[0]):
        m = ~np.isnan(V[i]) & ~np.isnan(x)
        if m.sum() >= 30:
            out[i] = stats.spearmanr(V[i][m], x[m]).correlation
    return pd.Series(out, index=M.index)


def run(seed: int = 0) -> pd.DataFrame:
    arms = load_prospective_mirna_arms()
    layers = load_cptac_layers("prospective")
    prot = layers["protein_z"]
    clin = load_prospective_clinical()

    caf = CF.deconv("cptac")["CAFs"]
    caf_mir, caf_gene = _caf_loading(arms, caf), _caf_loading(prot, caf)

    e = D.load_hallmark_edges()
    E = e[e.high_evidence][["miRNA", "gene"]].drop_duplicates()
    E = E[E.miRNA.isin(arms.index) & E.gene.isin(prot.index)].reset_index(drop=True)
    curated = set(zip(e.miRNA, e.gene))          # exclude ANY curated pair from the null, not just HE
    rng = np.random.default_rng(seed)
    sh = []
    for _ in range(N_SHUFFLE):                    # degree-preserving: permute the gene column
        g = rng.permutation(E.gene.to_numpy())
        s = pd.DataFrame({"miRNA": E.miRNA.to_numpy(), "gene": g})
        sh.append(s[[(m, gg) not in curated for m, gg in zip(s.miRNA, s.gene)]])
    S = pd.concat(sh).drop_duplicates().reset_index(drop=True)
    print(f"[null] HE edges {len(E)} | shuffled non-edges {len(S)}")

    d = pd.concat([E.assign(real=1), S.assign(real=0)], ignore_index=True)
    d["caf_mir"] = d.miRNA.map(caf_mir)
    d["caf_gene"] = d.gene.map(caf_gene)
    d["orientation"] = np.where(d.caf_mir * d.caf_gene < 0, "OPPOSITE", "SAME")
    for lbl, comp in [("rho_base", False), ("rho_adj", True)]:
        d[lbl] = _fast_rho(d, arms, prot, _covariates(clin, "prospective", composition=comp))

    rows = []
    for real in (1, 0):
        for o in ("OPPOSITE", "SAME"):
            k = d[(d.real == real) & (d.orientation == o)]
            rows.append({"group": "REAL_HE" if real else "SHUFFLED", "orientation": o, "n": int(k.rho_base.notna().sum()),
                         "rho_base": k.rho_base.mean(), "rho_adj": k.rho_adj.mean()})
    summ = pd.DataFrame(rows)

    OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT / "protein_compartment_null_edges.tsv", sep="\t", index=False)
    summ.to_csv(OUT / "protein_compartment_null_summary.tsv", sep="\t", index=False)

    print("\n" + summ.to_string(index=False))
    gr = lambda g, c: (summ[(summ.group == g) & (summ.orientation == "OPPOSITE")][c].iloc[0]
                       - summ[(summ.group == g) & (summ.orientation == "SAME")][c].iloc[0])
    print(f"\n  GRADIENT (opposite - same):  REAL {gr('REAL_HE','rho_base'):+.4f} -> {gr('REAL_HE','rho_adj'):+.4f}"
          f"   SHUFFLED {gr('SHUFFLED','rho_base'):+.4f} -> {gr('SHUFFLED','rho_adj'):+.4f}")
    print("  ⇒ the shuffled gradient MATCHES the real one: pairs that CANNOT repress show the same bulk "
          "anti-correlation ⇒ COMPARTMENT ARITHMETIC.")
    for lbl in ("rho_base", "rho_adj"):
        a, b = d[d.real == 1][lbl].dropna(), d[d.real == 0][lbl].dropna()
        p = stats.mannwhitneyu(a, b, alternative="less").pvalue
        print(f"  TRUE effect ({lbl}): real {a.mean():+.4f} vs shuffled {b.mean():+.4f} "
              f"(Δ={a.mean()-b.mean():+.4f}, MW p={p:.1e})")
    return d


if __name__ == "__main__":
    run()
