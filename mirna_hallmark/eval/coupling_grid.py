"""MH-115 — THE COUPLING GRID: layer x cohort x unit x universe, and — critically — HEURISTIC vs LEARNED
gene aggregation.

    .venv/bin/python3 -m mirna_hallmark.eval.coupling_grid

WHY THIS EXISTS (user, 2026-07-13). Every gene-level protein/mRNA validation in this subproject aggregates
miRNA -> gene with `pressure_build.compute_gene_pressure` — the **evidence-weighted abundance-share
HEURISTIC** (softmax_z share x cohort-z, target_norm=evidence_mass, sum). **That is NOT the learned model.**
So `cptac_validation`'s "independent protein validation" validates the PRESSURE HEURISTIC, not the §6b
posterior. The learned model's own gene-level repression field is

    M_g(s) = sum_f  beta_{g,f} * Xfam_f(s)          (beta from `learned/readouts_edges.tsv`, TCGA-fitted)

which is the quantity Bar-5 transports. This module scores BOTH aggregations, on BOTH layers, in BOTH
cohorts, under BOTH confounder blocks — so the grid has no missing cells and the two aggregations are
compared head-to-head on identical samples.

⚠ The learned beta is fitted on TCGA mRNA. Scoring it on TCGA mRNA is IN-SAMPLE (reported for reference and
marked `in_sample`); the honest comparisons are CPTAC mRNA and CPTAC protein, where beta is transported.

Outputs (`output/learned/`): `coupling_grid_genes.tsv`, `coupling_grid_summary.tsv`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark.eval.cptac_validation import (
    _covariates, build_pressure_matrices, load_cptac_layers, load_prospective_clinical,
    load_prospective_mirna_arms,
)
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import cptac_data as CD
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark import data_loaders as D

OUT = Path("mirna_hallmark/output/learned")
BETA = OUT / "readouts_edges.tsv"
MIN_N = 30


# --------------------------------------------------------------------------- #
# family dose matrices (families x samples), one per cohort
# --------------------------------------------------------------------------- #
def _fam_doses(X: pd.DataFrame) -> pd.DataFrame:
    """Family-collapsed dose matrix. §8: the latent beta is per SEED FAMILY, so the design column a gene's
    beta multiplies is the family dose, not the arm.

    ⚠ MH-125 FIX. This used `groupby(...).mean()` over member arms — but the CANONICAL collapse
    (`families.collapse_by_family`) is a TRUE RPM POOL, `log2(1 + Σ(2^x − 1))`. A mean is a DIFFERENT
    QUANTITY, and it is not the design the beta was fitted against. Also `.fillna(0)` to match
    `data.assemble_gene` (an undetected arm is RPM 0, not missing): the TCGA arm matrix is 66% NaN, and the
    old raw matmul propagated that into ~28% NaN gene scores, silently dropping 141 genes and every sample
    that was missing ANY one of a gene's families.
    """
    fam = FAM.family_of(list(X.index))
    fam = fam[fam.notna()]
    Xa = X.loc[[a for a in fam.index if a in X.index]]
    lin = np.power(2.0, Xa.to_numpy(float)) - 1.0                        # log2(RPM+1) → RPM
    lin = np.nan_to_num(lin, nan=0.0, posinf=0.0, neginf=0.0)            # undetected ⇒ RPM 0 (assemble_gene)
    pooled = pd.DataFrame(lin, index=Xa.index, columns=Xa.columns) \
        .groupby(fam.reindex(Xa.index).to_numpy()).sum()                 # POOL the members' RPM
    return np.log2(1.0 + pooled)                                         # back to log2 space


def _fam_z(Xf: pd.DataFrame) -> pd.DataFrame:
    """The scale the beta was FITTED on (`attribution_eb._prep`): each family dose z-scored across samples,
    with the §5 VARIANCE FLOOR (sd < 0.1 ⇒ column zeroed).

    ⚠ MH-125 FIX. `_learned_M` multiplied a beta fitted on the Z-SCORED design straight into RAW log2-RPM
    doses — which silently re-weights every beta_f by that family's σ_f. The weights and the design must live
    on the same scale.
    """
    sd = Xf.std(axis=1)
    Z = Xf.sub(Xf.mean(axis=1), axis=0).div(sd + 1e-9, axis=0)
    Z.loc[sd < 0.1] = 0.0
    return Z.fillna(0.0)


def _learned_M(beta: pd.DataFrame, Xf: pd.DataFrame, col: str = "beta",
               uniform: bool = False, raw_scale: bool = False) -> pd.DataFrame:
    """M_g(s) = Σ_f w_{g,f} · Xfam_f(s) over the gene's targeting seed families, on the **z-scored** design
    (the scale beta was fitted on — see `_fam_z`).

    THE TWO ABUNDANCE BENCHMARKS (`uniform=True`) — report BOTH; the choice is NOT innocent, and it is what
    flipped the MH-115 "abundance beats the heuristic" sub-claim:
      * `uniform=True, raw_scale=False` — uniform weights on the **z-scored** design. This is the MATCHED
        null: same genes, same families, same design, ONLY the weights differ (uniform vs beta). It answers
        "does the model know anything beyond *which* families, weighted equally?"
      * `uniform=True, raw_scale=True`  — uniform weights on the **raw** doses = TOTAL MIRNA ABUNDANCE on the
        gene. It answers the cruder "does the model know anything beyond *how much miRNA is around*?"
    """
    Xu = Xf if raw_scale else _fam_z(Xf)
    b = beta[beta.family.isin(Xu.index)]
    rows, idx = [], []
    for g, sub in b.groupby("gene"):
        w = sub.set_index("family")[col]
        w = w[w.notna()]
        if w.empty:
            continue
        wv = np.ones(len(w)) if uniform else w.to_numpy()
        rows.append(wv @ Xu.loc[w.index].to_numpy(float))
        idx.append(g)
    return pd.DataFrame(rows, index=idx, columns=Xu.columns)


def _score(Y: pd.DataFrame, S: pd.DataFrame, C: pd.DataFrame, tag: dict) -> pd.DataFrame:
    """Partial Spearman(score_g, outcome_g | C) for every gene, VECTORIZED.

    The obvious implementation — a `pd.concat` + `lstsq` per gene — costs ~18k pandas concats over the full
    grid (12 configs x 1.5k genes) and dominates the runtime. Instead: take the common sample set ONCE, rank
    and residualize each gene's outcome and score row on C with a single hat-matrix multiply, then read the
    correlation off a dot product. Same numbers (verified against the per-gene path, max|diff| < 1e-12).
    NaNs sit only in the outcome (scores are dense), so a gene's valid mask depends on the outcome row alone.
    """
    s = [c for c in Y.columns if c in S.columns and c in C.index]
    cv = C.loc[s].apply(pd.to_numeric, errors="coerce")
    keep = cv.dropna().index
    if len(keep) < MIN_N:
        return pd.DataFrame()
    A = np.column_stack([np.ones(len(keep)), cv.loc[keep].to_numpy(float)])
    H = A @ np.linalg.pinv(A.T @ A) @ A.T                     # hat matrix, samples x samples
    genes = [g for g in S.index if g in Y.index]
    Ym = Y.loc[genes, keep].to_numpy(float)
    Sm = S.loc[genes, keep].to_numpy(float)
    n_s = len(keep)
    rho = np.full(len(genes), np.nan)
    nn = np.zeros(len(genes), int)
    for i in range(len(genes)):
        m = ~np.isnan(Ym[i]) & ~np.isnan(Sm[i])
        k = int(m.sum())
        nn[i] = k
        if k < MIN_N:
            continue
        if k == n_s:                                          # complete row: reuse the prebuilt hat matrix
            ry = stats.rankdata(Ym[i]); rx = stats.rankdata(Sm[i])
            ey = ry - H @ ry; ex = rx - H @ rx
        else:                                                 # partial row: refit on its own sample subset
            Ai = A[m]
            ry = stats.rankdata(Ym[i][m]); rx = stats.rankdata(Sm[i][m])
            by, *_ = np.linalg.lstsq(Ai, ry, rcond=None)
            bx, *_ = np.linalg.lstsq(Ai, rx, rcond=None)
            ey = ry - Ai @ by; ex = rx - Ai @ bx
        sy, sx = ey.std(), ex.std()
        if sy < 1e-9 or sx < 1e-9:
            continue
        rho[i] = float((ey - ey.mean()) @ (ex - ex.mean()) / (len(ey) * sy * sx))
    df = pd.DataFrame({**tag, "gene": genes, "rho": rho, "n": nn})
    with np.errstate(divide="ignore", invalid="ignore"):      # two-sided p from the correlation + df
        dfree = np.maximum(df.n - A.shape[1] - 1, 1)
        t = df.rho * np.sqrt(dfree / np.maximum(1 - df.rho ** 2, 1e-12))
        df["p"] = 2 * stats.t.sf(np.abs(t), dfree)
    df.loc[df.rho.isna(), "p"] = np.nan
    return df


def run() -> pd.DataFrame:
    beta = pd.read_csv(BETA, sep="\t")
    hs = HallmarkSets.load()
    res = []

    # ---------------- TCGA (mRNA) -------------------------------------------
    d = LD._load()
    Xt, Yt = d["X"], d["Y"]
    caf_idx = CF.deconv("tcga").index
    pt = [p for p in Xt.columns if p in Yt.columns and p in caf_idx]
    Xt, Yt = Xt[pt], Yt[pt]
    Xft = _fam_doses(Xt)
    Mt = _learned_M(beta, Xft)                                   # LEARNED aggregation (z-scored design)
    At = _learned_M(beta, Xft, uniform=True)                     # ABUNDANCE, MATCHED null (uniform on z)
    Rt = _learned_M(beta, Xft, uniform=True, raw_scale=True)     # ABUNDANCE, RAW doses (total miRNA around)
    Pt = build_pressure_matrices(hs)["highev|gated"]             # HEURISTIC aggregation
    Pt = Pt[[c for c in Pt.columns if c in pt]]
    print(f"[grid] TCGA n={len(pt)} | learned genes {Mt.shape[0]} | heuristic genes {Pt.shape[0]}")
    for cname, comp in [("core", False), ("deconv", True)]:
        C = CF.build_C("tcga", pt, composition=comp)
        for agg, S in [("learned", Mt), ("heuristic", Pt), ("abundance", At), ("abundance_raw", Rt)]:
            res.append(_score(Yt, S, C, dict(layer="mRNA", cohort="TCGA", C=cname, agg=agg,
                                             in_sample=(agg == "learned"))))

    # ---------------- CPTAC (mRNA + protein) --------------------------------
    Xc = CD.arms()
    Mc_rna = CD.mrna()
    Pc = CD.protein()
    clin = load_prospective_clinical()
    Xfc = _fam_doses(Xc)
    Mc = _learned_M(beta, Xfc)                                   # LEARNED beta TRANSPORTED to CPTAC
    Ac = _learned_M(beta, Xfc, uniform=True)                     # ABUNDANCE, MATCHED null (uniform on z)
    Rc = _learned_M(beta, Xfc, uniform=True, raw_scale=True)     # ABUNDANCE, RAW doses
    Pc_h = build_pressure_matrices(hs, mirna=load_prospective_mirna_arms(),
                                   rna_for_gate=None)["highev|gated"] if False else None
    lay = load_cptac_layers("prospective")
    from mirna_hallmark.eval.cptac_validation import get_cohort_config, load_cct
    Pc_h = build_pressure_matrices(hs, mirna=load_prospective_mirna_arms(),
                                   rna_for_gate=load_cct(get_cohort_config("pancan122").rna_cct))["highev|gated"]
    for layer, Y in [("mRNA", Mc_rna), ("protein", lay["protein_z"])]:
        parts = [p for p in Y.columns if p in Xfc.columns]
        print(f"[grid] CPTAC {layer}: n={len(parts)}")
        for cname, comp in [("core", False), ("deconv", True)]:
            C = _covariates(clin, "prospective", composition=comp)
            for agg, S in [("learned", Mc), ("heuristic", Pc_h), ("abundance", Ac), ("abundance_raw", Rc)]:
                Ss = S[[c for c in S.columns if c in parts]]
                res.append(_score(Y[parts], Ss, C, dict(layer=layer, cohort="CPTAC", C=cname, agg=agg,
                                                        in_sample=False)))

    g = pd.concat([r for r in res if r is not None and len(r)], ignore_index=True)

    # --- PER-GENE COMPOSITION CHARACTERIZATION (user 2026-07-13): stratify every cell by the gene's own
    # retention class. Shapley of a linear aggregate is exactly beta_f => additive => gene_retention =
    # sum(beta_deconv)/sum(beta) (MH-111). PREDICTION: cell_intrinsic genes should validate at protein;
    # composition_explained genes should NOT (their bulk signal is compartment arithmetic).
    b = beta.dropna(subset=["beta"])
    alloc = (b.groupby("gene").apply(lambda s: s.beta_deconv.sum() / s.beta.sum()
                                     if abs(s.beta.sum()) > 1e-9 else np.nan)
             .rename("gene_retention").reset_index())
    alloc["comp_class"] = pd.cut(alloc.gene_retention, [-np.inf, 0.4, 0.7, np.inf],
                                 labels=["composition_explained", "partial", "cell_intrinsic"])
    g = g.merge(alloc, on="gene", how="left")

    OUT.mkdir(parents=True, exist_ok=True)
    g.to_csv(OUT / "coupling_grid_genes.tsv", sep="\t", index=False)

    from mirna_hallmark import stats as S_

    def _summ(df, keys):
        rows = []
        for k, s in df.groupby(keys, observed=True):
            v = s.rho.dropna()
            if len(v) < 10:
                continue
            q = s.dropna(subset=["p"]).copy()
            q["q"] = S_.bh_fdr(q["p"].fillna(1).to_numpy())
            rows.append({**dict(zip(keys, k if isinstance(k, tuple) else (k,))), "n_genes": len(v),
                         "mean_rho": v.mean(), "median_rho": v.median(), "frac_neg": (v < 0).mean(),
                         "fdr_neg": int(((q.q < 0.05) & (q.rho < 0)).sum()),
                         "p_vs_zero": stats.ttest_1samp(v, 0).pvalue})
        return pd.DataFrame(rows)

    summ = _summ(g, ["layer", "cohort", "C", "agg", "in_sample"]).sort_values(["layer", "cohort", "C", "agg"])
    strat = _summ(g[g.comp_class.notna()], ["layer", "cohort", "C", "agg", "comp_class"])
    summ.to_csv(OUT / "coupling_grid_summary.tsv", sep="\t", index=False)
    strat.to_csv(OUT / "coupling_grid_by_composition.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n=== THE GRID ===\n" + summ.to_string(index=False))
    print("\n=== STRATIFIED BY THE GENE'S COMPOSITION CLASS ===\n"
          + strat.sort_values(["layer", "cohort", "agg", "C", "comp_class"]).to_string(index=False))
    return g


if __name__ == "__main__":
    run()
