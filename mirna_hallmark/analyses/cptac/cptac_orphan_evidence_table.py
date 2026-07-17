"""Consolidated evidence table for the prediction-only orphan candidate edges (MH-36/39).

For each genuinely-novel edge (TargetScan-sequence-only: no CLIP, no miRTarBase) that showed a
translational protein signal and replicated, assemble one row with:

  * **TargetScan weighted score** of the novel prediction (``ts_weight`` = Σ|weighted context++|).
  * **Confounder-adjusted mRNA anti-correlation** (arm vs target mRNA, partial Spearman) in the
    three RNA cohorts: TCGA-BRCA (CPE+HRD), CPTAC-prospective (purity+CIN), CPTAC-105 (CPE+HRD).
  * **CPTAC protein** anti-correlation (``protein_z``) and the slope-free **translational**
    residual (``protein_resid``), prospective + TCGA-105, from `cptac_orphan_discovery`.
  * **Composition-robust** flag from the confounder pilot (`cptac_orphan_confound_pilot`).
  * **Target-gene context**: is the gene a curated miRTarBase **high-evidence** target at all, how
    many *other* experimentally-validated miRNAs hit it, and the strongest few (the "company it
    keeps" — distinguishes uncurated targets from densely-regulated hubs).

Output: `output/cptac_validation/orphan_discovery/prediction_only_evidence_table.tsv`

Run::

    .venv/bin/python3 -m mirna_hallmark.cptac_orphan_evidence_table
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cptac.cptac_orphan_discovery import (
    DISCOVERY_DIR,
    MIN_N,
    _edge_partial_spearman,
    _rank_matrix,
    _residualize,
    _spearman_p,
)
from mirna_hallmark.analyses.cptac.cptac_orphan_confound_pilot import (
    ENDO_MARKERS,
    PILOT_DIR,
    _estimate_fractions,
    _marker_score,
    _pam50_dummies,
)
from mirna_hallmark.eval.cptac_validation import (
    load_cptac_layers,
    load_prospective_clinical,
    load_prospective_mirna_arms,
)

N_TOP_HE = 6
CPTAC_CNA = C.REPO_ROOT / "data" / "CPTAC" / "HS_CPTAC_BRCA_2018_CNA.cct"


def _mrna_partial(edges: pd.DataFrame, arms: pd.DataFrame, rna: pd.DataFrame,
                  cov: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Partial Spearman(arm, target mRNA | cov) for each edge in one cohort."""
    arms_sub = arms.loc[arms.index.intersection(set(edges["miRNA"]))]
    rna_sub = rna.loc[rna.index.intersection(set(edges["gene"]))]
    arm_rank = _rank_matrix(_residualize(arms_sub, cov))
    rna_rank = _rank_matrix(_residualize(rna_sub, cov))
    rows = []
    for mir, gene in zip(edges["miRNA"], edges["gene"]):
        if mir in arm_rank.index and gene in rna_rank.index:
            rho, n = _edge_partial_spearman(arm_rank.loc[mir], rna_rank.loc[gene])
            rows.append({"miRNA": mir, "gene": gene, f"{tag}_mrna_rho": round(rho, 3) if np.isfinite(rho) else np.nan,
                         f"{tag}_mrna_p": _spearman_p(rho, n), f"{tag}_mrna_n": n})
        else:
            rows.append({"miRNA": mir, "gene": gene, f"{tag}_mrna_rho": np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Widest-confounder partial: base + PAM50 + ESTIMATE Immune/Stromal + endo + TARGET CN
# --------------------------------------------------------------------------- #
def _zrows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1, ddof=0).replace(0, np.nan), axis=0)


def _static_cov(base_cov: pd.DataFrame, rna_z: pd.DataFrame, clin: pd.DataFrame,
                estimate: pd.DataFrame, target_genes: set) -> pd.DataFrame:
    """base + endothelial marker + ESTIMATE Immune/Stromal (z) + PAM50 dummies. Gene-independent."""
    samples = list(base_cov.index)
    comp = pd.DataFrame(index=samples)
    comp["endo_score"] = _marker_score(rna_z, ENDO_MARKERS, target_genes).reindex(samples)
    comp["immune_estimate"] = estimate["ImmuneScore"].reindex(samples)
    comp["stromal_estimate"] = estimate["StromalScore"].reindex(samples)
    pam = _pam50_dummies(clin, samples)
    return pd.concat([base_cov, comp, pam], axis=1)


def _cptac_target_cn(genes: set) -> pd.DataFrame:
    cna = pd.read_csv(CPTAC_CNA, sep="\t", index_col=0)
    return cna.loc[cna.index.intersection(genes)]


def _resid_vec(v: pd.Series, cov: pd.DataFrame) -> pd.Series | None:
    """OLS residual of v on [1|cov], over samples complete in both."""
    v = pd.to_numeric(v, errors="coerce").dropna()
    cv = cov.reindex(v.index).apply(pd.to_numeric, errors="coerce")
    keep = cv.dropna().index
    if len(keep) < MIN_N:
        return None
    vv = v.loc[keep].values.astype(float)
    X = np.column_stack([np.ones(len(keep)), cv.loc[keep].values.astype(float)])
    beta, *_ = np.linalg.lstsq(X, vv, rcond=None)
    return pd.Series(vv - X @ beta, index=keep)


def _partial_vec(x: pd.Series, y: pd.Series, cov: pd.DataFrame) -> tuple[float, int]:
    """Partial Spearman(x, y | cov): residualize raw, rank residuals, correlate."""
    rx, ry = _resid_vec(x, cov), _resid_vec(y, cov)
    if rx is None or ry is None:
        return (np.nan, 0)
    common = rx.index.intersection(ry.index)
    if len(common) < MIN_N:
        return (np.nan, len(common))
    a = rankdata(rx.loc[common].values); b = rankdata(ry.loc[common].values)
    a = a - a.mean(); b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return (float((a * b).sum() / den) if den else np.nan, len(common))


def _wide_cohort(edges: pd.DataFrame, arm_mat: pd.DataFrame, resp_mats: dict,
                 static_cov: pd.DataFrame, cn_mat: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Per edge: widest-confounder partial (static_cov + target-gene CN) for each response layer."""
    rows = []
    for mir, gene in zip(edges["miRNA"], edges["gene"]):
        rec = {"miRNA": mir, "gene": gene}
        cov = static_cov.copy()
        cov["target_cn"] = cn_mat.loc[gene] if gene in cn_mat.index else np.nan
        if mir not in arm_mat.index:
            rows.append(rec); continue
        for layer, mat in resp_mats.items():
            if gene in mat.index:
                rho, n = _partial_vec(arm_mat.loc[mir], mat.loc[gene], cov)
                rec[f"{prefix}_{layer}_rho_wide"] = round(rho, 3) if np.isfinite(rho) else np.nan
                rec[f"{prefix}_{layer}_p_wide"] = _spearman_p(rho, n)
                rec[f"{prefix}_wide_n"] = n
        rows.append(rec)
    return pd.DataFrame(rows)


def _arm_guide_passenger(arms: Sequence[str]) -> pd.DataFrame:
    """miRGeneDB guide/passenger call per arm (MH-32 axis C), with `.N` multi-copy fallback."""
    import re
    from mirna_hallmark.analyses.edge_panels.edge_prior_refinement import mirgenedb_arm_status
    st = mirgenedb_arm_status(do_fetch=False).drop_duplicates("miRNA").set_index("miRNA")
    rows = []
    for arm in arms:
        cls, src = "absent/neutral", "-"
        for cand in (arm, re.sub(r"\.\d+$", "", arm)):  # exact, then strip ".1" multi-copy
            if cand in st.index:
                cls, src = st.loc[cand, "arm_class"], st.loc[cand, "arm_class_source"]
                break
        rows.append({"miRNA": arm, "arm_class": cls, "arm_class_source": src})
    return pd.DataFrame(rows)


def _he_target_context() -> pd.DataFrame:
    """Per target gene: HE-target flag, # of *other* experimentally-validated miRNAs, top ones."""
    he = D.high_evidence_edges(D.load_hallmark_edges())
    he = he[he["high_evidence"].astype(bool)].copy()
    # dedup (miRNA,gene), keep max evidence_score
    he = (he.sort_values("evidence_score", ascending=False)
            .drop_duplicates(["miRNA", "gene"]))
    rows = []
    for gene, sub in he.groupby("gene"):
        mirnas = list(dict.fromkeys(sub["miRNA"]))  # ordered by evidence (already sorted)
        rows.append({"gene": gene, "gene_is_HE_target": True,
                     "n_HE_mirnas": len(mirnas),
                     "top_HE_mirnas": ", ".join(mirnas[:N_TOP_HE])})
    return pd.DataFrame(rows)


def build(*, out_dir: Path = DISCOVERY_DIR) -> pd.DataFrame:
    cand = pd.read_csv(out_dir / "orphan_candidates.tsv", sep="\t")
    edges = cand[(cand["prior"] == "targetscan_only")
                 & cand["translational_candidate"].astype(bool)
                 & cand["tcga_replicated"].astype(bool)].copy()
    cols = ["miRNA", "gene", "ts_weight",
            "protein_rho", "protein_q", "resid_rho", "resid_q",
            "tcga_protein_rho", "tcga_resid_rho"]
    base = edges[cols].rename(columns={
        "protein_rho": "cptac_pros_protein_rho", "protein_q": "cptac_pros_protein_q",
        "resid_rho": "cptac_pros_resid_rho", "resid_q": "cptac_pros_resid_q",
        "tcga_protein_rho": "cptac_105_protein_rho", "tcga_resid_rho": "cptac_105_resid_rho"})

    # --- mRNA partials across the three RNA cohorts ---
    tcga_arms = D.load_mirna_arms()
    tcga_rna = D.load_rna()
    tcga_clin = D.load_clinical_strata().drop_duplicates("participant").set_index("participant")
    tcga_cov = tcga_clin[[c for c in C.CONFOUNDER_NUMERIC if c in tcga_clin.columns]]
    tcga_mrna = _mrna_partial(edges, tcga_arms, tcga_rna, tcga_cov, "tcga")

    pros_arms = load_prospective_mirna_arms()
    pros_layers = load_cptac_layers("prospective")
    pros_clin = load_prospective_clinical().drop_duplicates("participant").set_index("participant")
    pros_mrna = _mrna_partial(edges, pros_arms, pros_layers["rna_z"], pros_clin[["purity", "cin"]], "cptac_pros")

    t105_arms = D.load_mirna_arms()
    t105_layers = load_cptac_layers("tcga105")
    t105_mrna = _mrna_partial(edges, t105_arms, t105_layers["rna_z"], tcga_cov, "cptac_105")

    # --- WIDEST confounder: base + PAM50 + ESTIMATE Immune/Stromal + endo + TARGET CN ---
    target_genes = set(edges["gene"])
    tcga_cn = D.load_cnv_target_genes(sorted(target_genes))
    tcga_estimate = _zrows(pd.read_csv(
        C.OUTPUT_ROOT / "tissue_reference" / "estimate" / "tumor_estimate_scores.tsv",
        sep="\t", index_col=0).T).T  # participant x {Immune,Stromal} z
    tcga_rnaz_endo = _zrows(tcga_rna.loc[tcga_rna.index.intersection(ENDO_MARKERS)])
    tcga_static = _static_cov(tcga_cov, tcga_rnaz_endo, tcga_clin, tcga_estimate, target_genes)
    tcga_wide = _wide_cohort(edges, tcga_arms, {"mrna": tcga_rna}, tcga_static, tcga_cn, "tcga")

    pros_estimate = _estimate_fractions("prospective", PILOT_DIR)
    pros_static = _static_cov(pros_clin[["purity", "cin"]], pros_layers["rna_z"], pros_clin,
                              pros_estimate, target_genes)
    pros_cn = _cptac_target_cn(target_genes)
    pros_wide = _wide_cohort(edges, pros_arms,
                             {"mrna": pros_layers["rna_z"], "protein": pros_layers["protein_z"],
                              "resid": pros_layers["protein_resid"]},
                             pros_static, pros_cn, "cptac_pros")

    t105_estimate = _estimate_fractions("tcga105", PILOT_DIR)
    t105_static = _static_cov(tcga_cov, t105_layers["rna_z"], tcga_clin, t105_estimate, target_genes)
    t105_wide = _wide_cohort(edges, t105_arms,
                             {"mrna": t105_layers["rna_z"], "protein": t105_layers["protein_z"],
                              "resid": t105_layers["protein_resid"]},
                             t105_static, tcga_cn, "cptac_105")

    # --- composition-robust flag from the pilot (marker mode) ---
    pilot_f = PILOT_DIR / "pilot_edge_confound_comparison.tsv"
    if pilot_f.exists():
        pilot = pd.read_csv(pilot_f, sep="\t")[
            ["miRNA", "gene", "pros_comp_surv", "pros_comp_resid_rho"]].rename(
            columns={"pros_comp_surv": "composition_robust",
                     "pros_comp_resid_rho": "cptac_pros_resid_rho_compadj"})
    else:
        pilot = pd.DataFrame(columns=["miRNA", "gene"])

    # --- target-gene HE context + arm guide/passenger (miRGeneDB) ---
    he_ctx = _he_target_context()
    arm_gp = _arm_guide_passenger(sorted(edges["miRNA"].unique()))

    wide_cols = (tcga_wide.merge(pros_wide, on=["miRNA", "gene"], how="outer")
                          .merge(t105_wide, on=["miRNA", "gene"], how="outer"))
    wide_cols = wide_cols[["miRNA", "gene"] + [c for c in wide_cols.columns if c.endswith("_wide")]]
    # BH-FDR over the pilot set on the prospective translational (resid) widest p-values.
    from mirna_hallmark import stats as _S
    wide_cols["cptac_pros_resid_q_wide"] = _S.bh_fdr(
        wide_cols["cptac_pros_resid_p_wide"].fillna(1.0).values)
    wide_cols["survives_widest"] = (wide_cols["cptac_pros_resid_rho_wide"] < 0) & \
        (wide_cols["cptac_pros_resid_q_wide"] < 0.1)

    out = (base
           .merge(tcga_mrna[["miRNA", "gene", "tcga_mrna_rho", "tcga_mrna_p"]], on=["miRNA", "gene"], how="left")
           .merge(pros_mrna[["miRNA", "gene", "cptac_pros_mrna_rho"]], on=["miRNA", "gene"], how="left")
           .merge(t105_mrna[["miRNA", "gene", "cptac_105_mrna_rho"]], on=["miRNA", "gene"], how="left")
           .merge(wide_cols, on=["miRNA", "gene"], how="left")
           .merge(pilot, on=["miRNA", "gene"], how="left")
           .merge(arm_gp, on="miRNA", how="left")
           .merge(he_ctx, on="gene", how="left"))
    out["gene_is_HE_target"] = out["gene_is_HE_target"] == True  # noqa: E712 (NaN->False, no downcast warn)
    out["n_HE_mirnas"] = out["n_HE_mirnas"].fillna(0).astype(int)
    out["top_HE_mirnas"] = out["top_HE_mirnas"].fillna("")

    order = ["miRNA", "gene", "arm_class", "arm_class_source", "ts_weight",
             "composition_robust", "survives_widest",
             # base-adjusted (CPE+HRD / purity+CIN only)
             "tcga_mrna_rho", "cptac_pros_mrna_rho", "cptac_105_mrna_rho",
             "cptac_pros_protein_rho", "cptac_pros_resid_rho",
             "cptac_105_protein_rho", "cptac_105_resid_rho",
             # WIDEST-confounder (base + PAM50 + ESTIMATE Immune/Stromal + endo + target CN)
             "tcga_mrna_rho_wide", "tcga_mrna_p_wide",
             "cptac_pros_mrna_rho_wide", "cptac_105_mrna_rho_wide",
             "cptac_pros_protein_rho_wide", "cptac_pros_resid_rho_wide", "cptac_pros_resid_q_wide",
             "cptac_105_protein_rho_wide", "cptac_105_resid_rho_wide",
             "gene_is_HE_target", "n_HE_mirnas", "top_HE_mirnas"]
    out = out[[c for c in order if c in out.columns]].sort_values(
        "cptac_pros_resid_rho_wide")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=DISCOVERY_DIR)
    args = ap.parse_args()
    tbl = build(out_dir=args.out_dir)
    path = args.out_dir / "prediction_only_evidence_table.tsv"
    tbl.to_csv(path, sep="\t", index=False)
    pd.set_option("display.width", 240, "display.max_columns", 40)
    show = ["miRNA", "gene", "ts_weight",
            "cptac_pros_resid_rho", "cptac_pros_resid_rho_wide" if "cptac_pros_resid_rho_wide" in tbl else "cptac_pros_resid_rho",
            "cptac_pros_protein_rho_wide", "tcga_mrna_rho_wide", "cptac_105_resid_rho_wide",
            "n_HE_mirnas", "top_HE_mirnas"]
    print(tbl[[c for c in dict.fromkeys(show) if c in tbl.columns]].to_string(index=False))
    print(f"\n[evidence-table] wrote {path}  ({len(tbl)} edges)")


if __name__ == "__main__":
    main()
