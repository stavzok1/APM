"""Validate lost miRNA->gene couplings: is a broken coupling a genuine loss of *realized*
repression, or an artifact / a coupling masked by a competing driver?

Confound-safe lines of evidence, none of which reads raw expression *direction* (that would
re-introduce the level/platform confound the partial-Spearman design removes):

1. **Platform-matched (NAT->tumor) confidence flag.** NAT and tumor are the same TCGA platform and
   both ESTIMATE-adjusted, so a partial-rho significantly negative in NAT and ~0 in tumor is not a
   cross-platform / composition artifact (`nat_anchored_confident`).
2. **Module control (arm-wide vs target-specific).** Joining the Stage-2 target-module composite
   partial-rho: if the arm's WHOLE module decouples too (composite collapses NAT->tumor) the arm
   loses realized control broadly (`arm_wide_failure`); if the module stays coupled while this one
   edge breaks, it is `target_specific_escape`.
3. **Confounder conditioning (CN / meth / ATAC).** Re-fit the tumor edge partial-rho with the TARGET
   GENE's copy number (ASCAT3), promoter methylation, or ATAC accessibility added to the covariate
   set. If the repression RE-EMERGES (markedly more negative once the factor is controlled), the
   de-coupling is `{cnv,meth,atac}_masked` -- that factor was hiding an otherwise-intact miRNA effect.
4. **Competition control (`p_others`).** The linear analog of the softmax competition: condition the
   tumor partial on the gene's **competitive field** `P_others = Σ_{m'∈R(g), m'≠m} c(m',g,s)` (its
   total OTHER incoming realized pressure). This is the only term that makes a rank-coupling
   edge-specific (raw abundance is arm-global; `edge_w`/z/logRPM are rank-invisible or rank-identical
   — see `coupling_predictor_comparison`). If repression RE-EMERGES once the field is held constant a
   co-regulator was masking m (`competition_masked`: m is the real, hidden regulator, NOT a true
   loss); if the residual coupling COLLAPSES toward 0 the apparent effect was the field, not m
   (`competition_explained`).
5. **GTEx-primary healthy anchor.** The "was repressive in healthy" reference is GTEx (true healthy,
   consistent with the §11 acquired axis which rejects NAT for field-effect); NAT is a same-platform
   *confirmation* tier, not an independent trigger. `healthy_anchor_class ∈ {gtex_anchored, nat_only,
   neither}` reframes the (permissive GTEx-OR-NAT) lost set: `nat_only` edges are field-effect-suspect.

Output: `output/tissue_reference/decoupling_validation/decoupling_validation.tsv`.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cross_state.cross_state_coupling import EPI_MARKERS, Q_ALPHA, RHO_FLOOR, _state_covariates
from mirna_hallmark.analyses.cross_state.cross_state_deep_dive import _state_bundles
from mirna_hallmark.family_normal_reference import _participant
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.TISSUE_REFERENCE_DIR / "decoupling_validation"
STATE_DIR = C.TISSUE_REFERENCE_DIR / "mirna_state_class"
COREP = C.TISSUE_REFERENCE_DIR / "mirna_comovement" / "target_corepression.tsv"
MIN_N = 25
CNV_REEMERGE = 0.10  # tumor adj_rho must drop by >= this once CN is controlled to call cnv_masked


CNV_REEMERGE = 0.10  # tumor adj_rho must drop by >= this once CN is controlled to call cnv_masked
RHO_COUPLED = -0.10
RHO_DECOUPLED = -0.05


def _coregulator_context(edge_all: pd.DataFrame, sub: pd.DataFrame) -> pd.DataFrame:
    """For each lost edge, summarise share/uniqueness among co-regulators of the same gene."""
    ea = edge_all.copy()
    ea["c_tumor"] = pd.to_numeric(ea.get("c_tumor"), errors="coerce").fillna(0.0)
    for c in ["adj_rho_tumor", "adj_rho_nat", "edge_share_tumor", "n_targets_mirna"]:
        if c in ea.columns:
            ea[c] = pd.to_numeric(ea[c], errors="coerce")
    gm = ea.groupby("gene")["c_tumor"].transform("sum").replace(0, np.nan)
    ea["gene_regulator_share"] = ea["c_tumor"] / gm
    ea["regulator_share_rank"] = ea.groupby("gene")["c_tumor"].rank(ascending=False, method="min")
    ea["n_regulators"] = ea.groupby("gene")["miRNA"].transform("nunique")
    ea["regulator_share_rank_pct"] = ea["regulator_share_rank"] / ea["n_regulators"]
    # STRUCTURAL regulator rank (edge_w·sm, abundance spent once) — the de-double-counted
    # "is the focal arm the gene's preferred regulator" axis. The "competitors carry the
    # brake" call gates on this so a merely-abundant focal arm is not auto-called dominant.
    if "gene_struct_share_rank_tumor" in ea.columns:
        ea["struct_share_rank"] = pd.to_numeric(ea["gene_struct_share_rank_tumor"], errors="coerce")
    elif "gene_struct_share_tumor" in ea.columns:
        ss = pd.to_numeric(ea["gene_struct_share_tumor"], errors="coerce")
        ea["struct_share_rank"] = ss.groupby(ea["gene"]).rank(ascending=False, method="min")
    else:
        ea["struct_share_rank"] = np.nan  # pre-regen edge table: gate falls back to abundance rank
    ea["struct_share_rank_pct"] = ea["struct_share_rank"] / ea["n_regulators"]
    ea["coupled_tumor"] = ea["adj_rho_tumor"].fillna(0) < RHO_COUPLED
    ea["decoupled_tumor"] = ea["adj_rho_tumor"].fillna(0) > RHO_DECOUPLED

    focal = ea.set_index(["miRNA", "gene"])
    rows: List[dict] = []
    for _, r in sub.iterrows():
        arm, g = r["miRNA"], r["gene"]
        out: dict = {}
        if (arm, g) in focal.index:
            f = focal.loc[(arm, g)]
            struct_rank_pct = f.get("struct_share_rank_pct")
            # prefer the structural axis; fall back to abundance share only if absent
            high_gate = struct_rank_pct if pd.notna(struct_rank_pct) \
                else f.get("regulator_share_rank_pct", 1)
            out.update({
                "gene_regulator_share": float(f.get("gene_regulator_share", np.nan)),
                "regulator_share_rank": int(f.get("regulator_share_rank", np.nan))
                if pd.notna(f.get("regulator_share_rank")) else np.nan,
                "regulator_share_rank_pct": float(f.get("regulator_share_rank_pct", np.nan)),
                "gene_struct_share_tumor": float(f.get("gene_struct_share_tumor", np.nan)),
                "struct_share_rank_pct": float(struct_rank_pct)
                if pd.notna(struct_rank_pct) else np.nan,
                "arm_edge_share_tumor": float(f.get("edge_share_tumor", np.nan)),
                "arm_n_targets": int(f.get("n_targets_mirna", np.nan))
                if pd.notna(f.get("n_targets_mirna")) else np.nan,
                "focal_high_regulator_share": bool(high_gate <= 0.25),
            })
        comp = ea.loc[ea["gene"].eq(g) & ~ea["miRNA"].eq(arm)]
        out["n_competitors"] = int(len(comp))
        out["n_competitors_coupled_tumor"] = int(comp["coupled_tumor"].sum()) if len(comp) else 0
        out["n_competitors_decoupled_tumor"] = int(comp["decoupled_tumor"].sum()) if len(comp) else 0
        if len(comp):
            top = comp.sort_values("c_tumor", ascending=False).iloc[0]
            out["top_competitor_mirna"] = str(top["miRNA"])
            out["top_competitor_rho_tumor"] = float(top.get("adj_rho_tumor", np.nan))
            out["top_competitor_gene_share"] = float(top.get("gene_regulator_share", np.nan))
            out["top_competitor_coupled"] = bool(top.get("coupled_tumor", False))
            # are competitors still doing the repression while focal lost?
            out["competitors_carry_coupling"] = bool(out["n_competitors_coupled_tumor"] > 0
                                                     and out.get("focal_high_regulator_share") is False)
        else:
            out.update({"top_competitor_mirna": "", "top_competitor_rho_tumor": np.nan,
                        "top_competitor_gene_share": np.nan, "top_competitor_coupled": False,
                        "competitors_carry_coupling": False})
        rows.append({"miRNA": arm, "gene": g, **out})
    return pd.DataFrame(rows)


def _partial(y: pd.Series, x: pd.Series, cov: pd.DataFrame):
    from analysis.utils.common.loaders import partial_spearman
    shared = sorted(set(x.dropna().index) & set(y.dropna().index) & set(cov.dropna().index))
    if len(shared) < MIN_N or x.loc[shared].nunique() < 2 or y.loc[shared].nunique() < 2:
        return np.nan, np.nan, len(shared)
    r, p, n = partial_spearman(y.loc[shared], x.loc[shared], cov.reindex(shared))
    return float(r), float(p), int(n)


def run(*, out_dir: Path = OUT_DIR, classes=("lost", "nat_decoupled"),
        cov_source: str = "estimate_epi") -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    edge = pd.read_csv(STATE_DIR / "mirna_gene_edge_class.tsv", sep="\t")
    sub = edge[edge["joint_edge_class"].isin(classes)].copy()
    for c in ["adj_rho_gtex", "adj_rho_nat", "adj_rho_tumor",
              "adj_q_gtex", "adj_q_nat", "adj_q_tumor"]:
        if c in sub.columns:
            sub[c] = pd.to_numeric(sub[c], errors="coerce")
    print(f"[decouple] {len(sub)} edges in classes {classes}")

    # tumor matrices + covariates + target CNV (participant-keyed)
    hs = HallmarkSets.load()
    genes = sorted(set(sub["gene"]))
    # Include EPI_MARKERS so the estimate_epi epithelial covariate can be computed:
    # `_state_covariates(..., estimate_epi)` builds `epi = _metagene(rna, EPI_MARKERS)`,
    # and `_partial` drops any sample with a NaN covariate — if the epithelial markers are
    # absent from this (restricted) RNA bundle the `epi` column is all-NaN and every
    # conditioning refit (CN/meth/ATAC) silently returns n=0.
    bundles = _state_bundles(sorted(set(genes) | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                                    | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))
                                    | set(EPI_MARKERS)))
    mir_t, rna_t = bundles["tumor"]
    cov_t = _state_covariates("tumor", rna_t, hs, cov_source)
    cnv = D.load_cnv_target_genes(genes)  # gene x participant integer CN
    # map tumor columns -> participant for CN lookup
    col_part = {c: _participant(c) for c in rna_t.columns}
    # per-sample methylation (Hallmark genes) + ATAC (direct + subtype-representative fallback)
    from mirna_hallmark.analyses.misc.escape_mechanism import load_atac_per_sample, load_meth_per_sample
    meth_mat = load_meth_per_sample(genes)            # gene x participant beta
    atac_val, atac_imp = load_atac_per_sample(genes)  # gene x participant + imputed mask
    print(f"[decouple] meth genes={0 if meth_mat.empty else meth_mat.shape[0]} "
          f"atac genes={0 if atac_val.empty else atac_val.shape[0]}")

    # Competitive field P_others(g,s) = the gene's TOTAL OTHER incoming realized pressure
    # (Σ over R(g)\{m} of c(m',g,s), softmax_logrpm). Built per-sample from the edge pressure map
    # over ALL regulators of the decoupling genes, so block (6) can hold it constant. This is the
    # linear, interpretable form of the softmax competition the predictor cannot encode for a rank
    # test (the softmax is its nonlinear rank-1 proxy).
    from mirna_hallmark.pressure_build import load_mirtar_edges
    from mirna_hallmark.pressure_engine import compute_edge_pressure_map
    edges_all = load_mirtar_edges(genes, resolve_arms=True)
    pmap = compute_edge_pressure_map(
        edges_all, mir_t, genes=genes,
        expr_mode=C.PRESSURE_ATTRIBUTION_EXPR_MODE,  # type: ignore[arg-type]
        target_norm=C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
    )
    gene_total: Dict[str, pd.Series] = {}
    for (a_, g_), s_ in pmap.items():
        gene_total[g_] = s_.copy() if g_ not in gene_total else gene_total[g_].add(s_, fill_value=0.0)
    print(f"[decouple] competitive field: {len(gene_total)} genes, {len(pmap)} edge-pressure series")

    # module-control: arm-level composite partial-rho per state
    mod = pd.DataFrame()
    if COREP.exists():
        mod = pd.read_csv(COREP, sep="\t").set_index("miRNA")
    # gene-level realized net-repression (gene-side mirror): is this gene under genuine net miRNA
    # repression, and did that survive into tumor? Lets us read a lost edge as edge-local escape
    # (gene still net-repressed by other regulators) vs gene-wide de-repression.
    gcr = pd.DataFrame()
    gcr_path = COREP.parent / "gene_corepression.tsv"
    if gcr_path.exists():
        gcr = pd.read_csv(gcr_path, sep="\t").set_index("gene")

    rows: List[dict] = []
    coreg = _coregulator_context(edge, sub)
    coreg_idx = coreg.set_index(["miRNA", "gene"])
    for _, e in sub.iterrows():
        arm, g = e["miRNA"], e["gene"]
        row = {"miRNA": arm, "gene": g, "joint_edge_class": e["joint_edge_class"],
               "adj_rho_nat": e.get("adj_rho_nat"), "adj_rho_tumor": e.get("adj_rho_tumor"),
               "adj_rho_gtex": e.get("adj_rho_gtex"),
               "n_regulators_gene": e.get("n_regulators_gene"),
               "edge_share_tumor": e.get("edge_share_tumor"),
               "gene_share_tumor": e.get("gene_share_tumor"),
               "gene_share_gtex": e.get("gene_share_gtex"),
               "gene_share_rank_tumor": e.get("gene_share_rank_tumor")}
        # gene-side realized net-repression context (whole-regulator view of the target)
        if not gcr.empty and g in gcr.index:
            gr = gcr.loc[g]
            row["gene_rho_pressure_tumor"] = gr.get("rho_gene_pressure_tumor")
            row["gene_rho_pressure_nat"] = gr.get("rho_gene_pressure_nat")
            row["gene_repression_class"] = gr.get("gene_repression_class")
            row["gene_net_repressed_tumor"] = bool(gr.get("gene_net_repressed_tumor", False))
            # is the edge's de-coupling local (gene still net-repressed by others) or gene-wide?
            row["escape_scope"] = ("edge_local_escape" if row["gene_net_repressed_tumor"]
                                   else "gene_wide_derepression"
                                   if str(gr.get("gene_repression_class")) == "lost_repression"
                                   else "gene_not_repressed")
        else:
            row["gene_rho_pressure_tumor"] = row["gene_rho_pressure_nat"] = np.nan
            row["gene_repression_class"] = "na"
            row["gene_net_repressed_tumor"] = False
            row["escape_scope"] = "unknown"
        # (1) platform-matched confidence
        row["nat_anchored_confident"] = bool(
            (e.get("adj_rho_nat", 0) < -0.10) and (e.get("adj_q_nat", 1) < 0.10)
            and (e.get("adj_rho_tumor", 0) > -0.05))
        # (1b) GTEx-primary healthy anchor: true-healthy repression is the biological call; NAT
        # confirms (same platform) but does not independently establish a "loss" (field-effect).
        rg, qg = e.get("adj_rho_gtex"), e.get("adj_q_gtex")
        rn, qn = e.get("adj_rho_nat"), e.get("adj_q_nat")
        row["gtex_repressive"] = bool(pd.notna(rg) and rg < -RHO_FLOOR and pd.notna(qg) and qg < Q_ALPHA)
        row["nat_repressive"] = bool(pd.notna(rn) and rn < -RHO_FLOOR and pd.notna(qn) and qn < Q_ALPHA)
        row["healthy_anchor_class"] = ("gtex_anchored" if row["gtex_repressive"]
                                       else "nat_only" if row["nat_repressive"]
                                       else "neither")
        # (2) module control
        if arm in mod.index:
            cn = float(mod.loc[arm].get("rho_composite_nat", np.nan))
            ct = float(mod.loc[arm].get("rho_composite_tumor", np.nan))
            row["module_rho_nat"], row["module_rho_tumor"] = cn, ct
            if np.isfinite(cn) and np.isfinite(ct):
                # module still repressed in tumor (<-0.1) but edge lost = target-specific escape;
                # module also collapses (tumor weaker, >-0.05) = arm-wide realized failure
                row["decoupling_mechanism"] = ("target_specific_escape" if ct < -0.10
                                               else "arm_wide_failure" if ct > -0.05
                                               else "intermediate")
            else:
                row["decoupling_mechanism"] = "unknown"
        else:
            row["module_rho_nat"] = row["module_rho_tumor"] = np.nan
            row["decoupling_mechanism"] = "unknown"
        # (3) competing-regulator (target CN) conditioning in tumor
        if arm in mir_t.index and g in rna_t.index and g in cnv.index:
            x = pd.to_numeric(mir_t.loc[arm], errors="coerce")
            y = pd.to_numeric(rna_t.loc[g], errors="coerce")
            cn_g = pd.Series({c: cnv.loc[g].get(col_part.get(c)) for c in rna_t.columns},
                             dtype=float)
            base_r, _bp, _bn = _partial(y, x, cov_t)
            cov_cn = cov_t.copy()
            cov_cn["target_cn"] = cn_g.reindex(cov_cn.index)
            cnv_r, _cp, cnv_n = _partial(y, x, cov_cn)
            row["adj_rho_tumor_estimate"] = base_r
            row["adj_rho_tumor_plus_cn"] = cnv_r
            row["cn_reemergence"] = (base_r - cnv_r) if (np.isfinite(base_r) and np.isfinite(cnv_r)) else np.nan
            row["n_cn"] = cnv_n
            row["cnv_masked"] = bool(np.isfinite(row["cn_reemergence"])
                                     and row["cn_reemergence"] >= CNV_REEMERGE
                                     and cnv_r < -0.10)
        else:
            row["adj_rho_tumor_estimate"] = row["adj_rho_tumor_plus_cn"] = np.nan
            row["cn_reemergence"] = np.nan
            row["cnv_masked"] = False
        # (4) methylation conditioning (target promoter beta as extra covariate)
        if (not meth_mat.empty and arm in mir_t.index and g in rna_t.index
                and g in meth_mat.index):
            x = pd.to_numeric(mir_t.loc[arm], errors="coerce")
            y = pd.to_numeric(rna_t.loc[g], errors="coerce")
            mrow = meth_mat.loc[g]
            meth_g = pd.Series({c: mrow.get(col_part.get(c)) for c in rna_t.columns}, dtype=float)
            base_r, _bp, _bn = _partial(y, x, cov_t)
            cov_m = cov_t.copy()
            cov_m["target_meth"] = meth_g.reindex(cov_m.index)
            meth_r, _mp, meth_n = _partial(y, x, cov_m)
            row["adj_rho_tumor_plus_meth"] = meth_r
            row["meth_reemergence"] = (base_r - meth_r) if (np.isfinite(base_r) and np.isfinite(meth_r)) else np.nan
            row["n_meth"] = meth_n
            row["meth_masked"] = bool(np.isfinite(row["meth_reemergence"])
                                      and row["meth_reemergence"] >= CNV_REEMERGE and meth_r < -0.10)
        else:
            row["adj_rho_tumor_plus_meth"] = np.nan
            row["meth_reemergence"] = np.nan
            row["meth_masked"] = False
            row["n_meth"] = 0
        # (5) ATAC conditioning (direct where measured, subtype-representative otherwise)
        if (not atac_val.empty and arm in mir_t.index and g in rna_t.index
                and g in atac_val.index):
            x = pd.to_numeric(mir_t.loc[arm], errors="coerce")
            y = pd.to_numeric(rna_t.loc[g], errors="coerce")
            arow = atac_val.loc[g]
            irow = atac_imp.loc[g] if g in atac_imp.index else None
            atac_g = pd.Series({c: arow.get(col_part.get(c)) for c in rna_t.columns}, dtype=float)
            base_r, _bp, _bn = _partial(y, x, cov_t)
            cov_a = cov_t.copy()
            cov_a["target_atac"] = atac_g.reindex(cov_a.index)
            atac_r, _ap, atac_n = _partial(y, x, cov_a)
            row["adj_rho_tumor_plus_atac"] = atac_r
            row["atac_reemergence"] = (base_r - atac_r) if (np.isfinite(base_r) and np.isfinite(atac_r)) else np.nan
            row["n_atac_cond"] = atac_n
            # fraction of conditioned samples whose ATAC is subtype-imputed (interpretation caveat)
            if irow is not None:
                imp_frac = float(pd.Series({c: bool(irow.get(col_part.get(c), True))
                                            for c in rna_t.columns}).mean())
            else:
                imp_frac = np.nan
            row["atac_frac_imputed"] = imp_frac
            row["atac_masked"] = bool(np.isfinite(row["atac_reemergence"])
                                      and row["atac_reemergence"] >= CNV_REEMERGE and atac_r < -0.10)
        else:
            row["adj_rho_tumor_plus_atac"] = np.nan
            row["atac_reemergence"] = np.nan
            row["atac_masked"] = False
            row["n_atac_cond"] = 0
            row["atac_frac_imputed"] = np.nan
        # (5b) STRICT ATAC proof-of-concept: ONLY the directly-assayed tumors (n≈74, no subtype
        # imputation). Asks what *measured* promoter accessibility explains in exactly those cases —
        # base partial-ρ and the +ATAC partial-ρ are BOTH recomputed on the direct subset so the
        # comparison is apples-to-apples (a smaller, noisier, but un-imputed test).
        row.update({"n_atac_direct": 0, "adj_rho_tumor_direct": np.nan,
                    "adj_rho_tumor_direct_plus_atac": np.nan,
                    "atac_direct_reemergence": np.nan, "atac_direct_masked": False})
        if (not atac_val.empty and arm in mir_t.index and g in rna_t.index
                and g in atac_val.index and g in atac_imp.index):
            irow = atac_imp.loc[g]
            direct_cols = [c for c in rna_t.columns
                           if not bool(irow.get(col_part.get(c), True))]
            if len(direct_cols) >= MIN_N:
                xs = pd.to_numeric(mir_t.loc[arm], errors="coerce").reindex(direct_cols)
                ys = pd.to_numeric(rna_t.loc[g], errors="coerce").reindex(direct_cols)
                arow = atac_val.loc[g]
                atac_gs = pd.Series({c: arow.get(col_part.get(c)) for c in direct_cols}, dtype=float)
                cov_ds = cov_t.reindex(direct_cols)
                base_d, _bp, base_dn = _partial(ys, xs, cov_ds)
                cov_da = cov_ds.copy()
                cov_da["target_atac"] = atac_gs.reindex(cov_da.index)
                atac_d, _adp, _adn = _partial(ys, xs, cov_da)
                row["n_atac_direct"] = base_dn
                row["adj_rho_tumor_direct"] = base_d
                row["adj_rho_tumor_direct_plus_atac"] = atac_d
                row["atac_direct_reemergence"] = ((base_d - atac_d)
                                                  if (np.isfinite(base_d) and np.isfinite(atac_d)) else np.nan)
                row["atac_direct_masked"] = bool(np.isfinite(row["atac_direct_reemergence"])
                                                 and row["atac_direct_reemergence"] >= CNV_REEMERGE
                                                 and atac_d < -0.10)
            else:
                row["n_atac_direct"] = len(direct_cols)
        # (6) COMPETITION control: condition the tumor partial on the competitive field P_others.
        #     re-emergence => a co-regulator masked m (competition_masked: not a true loss);
        #     collapse toward 0 => the residual coupling was the field, not m (competition_explained).
        row.update({"adj_rho_tumor_plus_compet": np.nan, "compet_reemergence": np.nan,
                    "n_compet": 0, "competition_masked": False, "competition_explained": False,
                    "p_others_present": False})
        if arm in mir_t.index and g in rna_t.index and g in gene_total and (arm, g) in pmap:
            x = pd.to_numeric(mir_t.loc[arm], errors="coerce")
            y = pd.to_numeric(rna_t.loc[g], errors="coerce")
            p_others = (gene_total[g] - pmap[(arm, g)])
            cov_c = cov_t.copy()
            cov_c["p_others"] = p_others.reindex(cov_c.index)
            # base on the SAME samples as the +p_others fit so Δ is apples-to-apples
            valid = cov_c.dropna().index
            base_r, _bp, _bn = _partial(y, x, cov_t.loc[valid])
            comp_r, _cp, comp_n = _partial(y, x, cov_c)
            row["adj_rho_tumor_plus_compet"] = comp_r
            row["n_compet"] = comp_n
            row["p_others_present"] = True
            if np.isfinite(base_r) and np.isfinite(comp_r):
                row["compet_reemergence"] = base_r - comp_r
                row["competition_masked"] = bool(row["compet_reemergence"] >= CNV_REEMERGE
                                                 and comp_r < -0.10)
                row["competition_explained"] = bool(row["compet_reemergence"] <= -CNV_REEMERGE)
        if (arm, g) in coreg_idx.index:
            row.update(coreg_idx.loc[(arm, g)].to_dict())
        rows.append(row)

    df = pd.DataFrame(rows)
    sort_cols = [c for c in ["competition_masked", "gtex_repressive", "nat_anchored_confident",
                             "compet_reemergence", "cn_reemergence"] if c in df.columns]
    df = df.sort_values(sort_cols, ascending=[False] * len(sort_cols)).reset_index(drop=True)
    df.to_csv(out_dir / "decoupling_validation.tsv", sep="\t", index=False)

    mech = df["decoupling_mechanism"].value_counts().to_dict()
    anchor = df["healthy_anchor_class"].value_counts().to_dict() if "healthy_anchor_class" in df.columns else {}
    manifest = {"module": "mirna_hallmark.decoupling_validation",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "classes": list(classes), "cov_source": cov_source, "n_edges": int(len(df)),
                "n_nat_confident": int(df["nat_anchored_confident"].sum()),
                "n_cnv_masked": int(df["cnv_masked"].sum()),
                "n_meth_masked": int(df["meth_masked"].sum()) if "meth_masked" in df.columns else 0,
                "n_atac_masked": int(df["atac_masked"].sum()) if "atac_masked" in df.columns else 0,
                "n_atac_direct_tested": int((df.get("n_atac_direct", pd.Series(0)) >= MIN_N).sum()),
                "n_atac_direct_masked": int(df["atac_direct_masked"].sum()) if "atac_direct_masked" in df.columns else 0,
                "n_competition_tested": int(df["p_others_present"].sum()) if "p_others_present" in df.columns else 0,
                "n_competition_masked": int(df["competition_masked"].sum()) if "competition_masked" in df.columns else 0,
                "n_competition_explained": int(df["competition_explained"].sum()) if "competition_explained" in df.columns else 0,
                "healthy_anchor_class_counts": anchor,
                "n_gtex_anchored": int(df["gtex_repressive"].sum()) if "gtex_repressive" in df.columns else 0,
                "mechanism_counts": mech}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[decouple] NAT-anchored confident: {int(df['nat_anchored_confident'].sum())}/{len(df)}")
    print(f"[decouple] mechanism: {mech}")
    print(f"[decouple] CNV-masked: {int(df['cnv_masked'].sum())} | "
          f"meth-masked: {int(df['meth_masked'].sum()) if 'meth_masked' in df.columns else 0} "
          f"(n with meth: {int((df.get('n_meth', pd.Series(0)) > 0).sum())}) | "
          f"ATAC-masked: {int(df['atac_masked'].sum()) if 'atac_masked' in df.columns else 0} "
          f"(n with atac: {int((df.get('n_atac_cond', pd.Series(0)) > 0).sum())})")
    if "atac_direct_masked" in df.columns:
        nd = int((df["n_atac_direct"] >= MIN_N).sum())
        print(f"[decouple] STRICT direct-ATAC (n≈74, no imputation): {nd} edges testable, "
              f"{int(df['atac_direct_masked'].sum())} atac_direct_masked; "
              f"median n_direct={int(df.loc[df['n_atac_direct']>=MIN_N,'n_atac_direct'].median()) if nd else 0}")
    if "competitors_carry_coupling" in df.columns:
        print(f"[decouple] focal lost but competitors still coupled: "
              f"{int(df['competitors_carry_coupling'].sum())}/{len(df)}")
    if "p_others_present" in df.columns:
        nt = int(df["p_others_present"].sum())
        print(f"[decouple] COMPETITION control (P_others): {nt} edges tested | "
              f"competition_masked {int(df['competition_masked'].sum())} (co-regulator hid m → not a true loss) | "
              f"competition_explained {int(df['competition_explained'].sum())} (residual was the field, not m)")
    if "healthy_anchor_class" in df.columns:
        print(f"[decouple] GTEx-primary anchor: {anchor} "
              f"(nat_only = field-effect-suspect, not GTEx-confirmed)")
        print(f"[decouple] focal high regulator share (top quartile on gene): "
              f"{int(df.get('focal_high_regulator_share', pd.Series(False)).sum())}/{len(df)}")
    if "escape_scope" in df.columns:
        print(f"[decouple] escape scope: {df['escape_scope'].value_counts().to_dict()}")
    conf = df[df["nat_anchored_confident"]].head(12)
    if not conf.empty:
        print("\n[decouple] top NAT-anchored confident lost edges:")
        cols = ["miRNA", "gene", "adj_rho_nat", "adj_rho_tumor", "module_rho_tumor",
                "decoupling_mechanism", "cn_reemergence", "cnv_masked"]
        print(conf[cols].round(3).to_string(index=False))
    return df


def _fast_pspear(Y: np.ndarray, X: np.ndarray, C: np.ndarray) -> float:
    """OLS-residualize Y and X on design C (intercept added), then Spearman of residuals.
    Bit-matches `analysis.utils.common.loaders.partial_spearman`, ~5x faster (no pandas)."""
    from scipy.stats import spearmanr
    A = np.column_stack([np.ones(len(Y)), C])
    ry = Y - A @ np.linalg.lstsq(A, Y, rcond=None)[0]
    rx = X - A @ np.linalg.lstsq(A, X, rcond=None)[0]
    r = spearmanr(ry, rx).statistic
    return float(r) if np.isfinite(r) else np.nan


def _boot_condition(Y: np.ndarray, X: np.ndarray, C0: np.ndarray, V: np.ndarray,
                    n_boot: int, rng) -> dict:
    """Bootstrap CI on the re-emergence Δ = base_ρ − (ρ | +V) for one conditioning variable V
    (V already aligned & NaN-clean with Y/X/C0). Returns point Δ, 95% CI, and the fraction of
    resamples in which the masked rule (Δ≥0.10 ∧ ρ_plus<−0.10) holds (robustness)."""
    C1 = np.column_stack([C0, V])
    base_r = _fast_pspear(Y, X, C0)
    plus_r = _fast_pspear(Y, X, C1)
    ds = np.empty(n_boot); maskb = 0
    m = len(Y)
    for b in range(n_boot):
        idx = rng.integers(0, m, m)
        br = _fast_pspear(Y[idx], X[idx], C0[idx]); pr = _fast_pspear(Y[idx], X[idx], C1[idx])
        ds[b] = br - pr
        if (br - pr) >= CNV_REEMERGE and pr < -0.10:
            maskb += 1
    d = base_r - plus_r
    return {"reemergence": d, "plus_rho": plus_r, "n": m,
            "ci_lo": float(np.percentile(ds, 2.5)), "ci_hi": float(np.percentile(ds, 97.5)),
            "masked_robust_frac": maskb / n_boot,
            "masked_robust": bool(d >= CNV_REEMERGE and plus_r < -0.10 and maskb / n_boot >= 0.90)}


def evaluate(*, out_dir: Path = OUT_DIR, classes=("lost", "nat_decoupled"),
             cov_source: str = "estimate_epi", n_boot: int = 200, n_perm: int = 200,
             seed: int = 0) -> pd.DataFrame:
    """PERFORMANCE evaluation of the new competition-controlled decoupling modeling.

    For each lost edge it answers three questions the point estimate cannot:
      1. **Effect-size + bootstrap CI** on the competition re-emergence Δ = base_ρ − (ρ | P_others)
         and on the conditioned ρ (participant resampling). A `competition_masked` call is *robust*
         only if it holds in ≥90% of bootstraps.
      2. **Permutation null** — shuffle P_others vs the samples and recompute Δ; the null masked
         rate calibrates whether the observed masking is signal or noise (the key 'does it perform').
      3. **Loss-call stability** — fraction of bootstraps in which the tumour coupling stays flat
         (the decoupling premise).
    Output: `decoupling_validation_eval.tsv` + a printed performance summary.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    edge = pd.read_csv(STATE_DIR / "mirna_gene_edge_class.tsv", sep="\t")
    sub = edge[edge["joint_edge_class"].isin(classes)].copy()
    genes = sorted(set(sub["gene"]))
    hs = HallmarkSets.load()
    bundles = _state_bundles(sorted(set(genes) | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                                    | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))
                                    | set(EPI_MARKERS)))
    mir_t, rna_t = bundles["tumor"]
    cov_t = _state_covariates("tumor", rna_t, hs, cov_source)
    from mirna_hallmark.pressure_build import load_mirtar_edges
    from mirna_hallmark.pressure_engine import compute_edge_pressure_map
    edges_all = load_mirtar_edges(genes, resolve_arms=True)
    pmap = compute_edge_pressure_map(
        edges_all, mir_t, genes=genes,
        expr_mode=C.PRESSURE_ATTRIBUTION_EXPR_MODE, target_norm=C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
    )
    gene_total: Dict[str, pd.Series] = {}
    for (a_, g_), s_ in pmap.items():
        gene_total[g_] = s_.copy() if g_ not in gene_total else gene_total[g_].add(s_, fill_value=0.0)
    # confounder matrices (same as run) so every masking verdict gets a bootstrap-robust flag
    cnv = D.load_cnv_target_genes(genes)
    col_part = {c: _participant(c) for c in rna_t.columns}
    from mirna_hallmark.analyses.misc.escape_mechanism import load_atac_per_sample, load_meth_per_sample
    meth_mat = load_meth_per_sample(genes)
    atac_val, _atac_imp = load_atac_per_sample(genes)

    print(f"[eval] {len(sub)} edges | n_boot={n_boot} n_perm={n_perm}")
    rows: List[dict] = []
    for _, e in sub.iterrows():
        arm, g = e["miRNA"], e["gene"]
        rec = {"miRNA": arm, "gene": g, "joint_edge_class": e["joint_edge_class"]}
        if not (arm in mir_t.index and g in rna_t.index and g in gene_total and (arm, g) in pmap):
            rows.append({**rec, "n": 0}); continue
        x = pd.to_numeric(mir_t.loc[arm], errors="coerce")
        y = pd.to_numeric(rna_t.loc[g], errors="coerce")
        p_others = gene_total[g] - pmap[(arm, g)]
        dfm = pd.concat([y.rename("y"), x.rename("x"), cov_t,
                         p_others.rename("p_others")], axis=1).dropna()
        m = len(dfm)
        if m < MIN_N:
            rows.append({**rec, "n": m}); continue
        Y = dfm["y"].values; X = dfm["x"].values
        C0 = dfm[[c for c in cov_t.columns]].values
        P = dfm["p_others"].values
        C1 = np.column_stack([C0, P])
        base_r = _fast_pspear(Y, X, C0)
        comp_r = _fast_pspear(Y, X, C1)
        d_obs = base_r - comp_r
        # bootstrap (participant resampling)
        ds = np.empty(n_boot); crs = np.empty(n_boot); maskb = 0; tumflat = 0
        for b in range(n_boot):
            idx = rng.integers(0, m, m)
            br = _fast_pspear(Y[idx], X[idx], C0[idx]); cr = _fast_pspear(Y[idx], X[idx], C1[idx])
            ds[b] = br - cr; crs[b] = cr
            if (br - cr) >= CNV_REEMERGE and cr < -0.10:
                maskb += 1
            if br > -0.05:
                tumflat += 1
        # permutation null (break P_others <-> sample alignment)
        dperm = np.empty(n_perm); maskp = 0
        for pq in range(n_perm):
            Cp = np.column_stack([C0, P[rng.permutation(m)]])
            cr = _fast_pspear(Y, X, Cp); dperm[pq] = base_r - cr
            if (base_r - cr) >= CNV_REEMERGE and cr < -0.10:
                maskp += 1
        masked_obs = bool(d_obs >= CNV_REEMERGE and comp_r < -0.10)
        dperm_p = float(np.mean(np.abs(dperm) >= abs(d_obs)))
        # bootstrap CI + robustness for the CN / meth / ATAC confounders (uniform with competition)
        cond_out: dict = {}
        cond_series = {}
        if g in cnv.index:
            cond_series["cn"] = pd.Series({c: cnv.loc[g].get(col_part.get(c)) for c in rna_t.columns}, dtype=float)
        if not meth_mat.empty and g in meth_mat.index:
            cond_series["meth"] = pd.Series({c: meth_mat.loc[g].get(col_part.get(c)) for c in rna_t.columns}, dtype=float)
        if not atac_val.empty and g in atac_val.index:
            cond_series["atac"] = pd.Series({c: atac_val.loc[g].get(col_part.get(c)) for c in rna_t.columns}, dtype=float)
        for vname, vser in cond_series.items():
            dfv = pd.concat([y.rename("y"), x.rename("x"), cov_t, vser.rename("v")], axis=1).dropna()
            if len(dfv) < MIN_N:
                cond_out[f"{vname}_masked_robust"] = False
                continue
            r = _boot_condition(dfv["y"].values, dfv["x"].values,
                                dfv[[c for c in cov_t.columns]].values, dfv["v"].values, n_boot, rng)
            cond_out[f"{vname}_reemergence"] = r["reemergence"]
            cond_out[f"{vname}_ci_lo"] = r["ci_lo"]
            cond_out[f"{vname}_ci_hi"] = r["ci_hi"]
            cond_out[f"{vname}_masked_robust"] = r["masked_robust"]
        rows.append({**rec, "n": m, "base_rho": base_r, "comp_rho": comp_r, **cond_out,
                     "compet_reemergence": d_obs,
                     "d_ci_lo": float(np.percentile(ds, 2.5)), "d_ci_hi": float(np.percentile(ds, 97.5)),
                     "comp_rho_ci_lo": float(np.percentile(crs, 2.5)),
                     "comp_rho_ci_hi": float(np.percentile(crs, 97.5)),
                     "masked_boot_frac": maskb / n_boot, "masked_null_frac": maskp / n_perm,
                     "dperm_p": dperm_p, "tumor_flat_boot_frac": tumflat / n_boot,
                     "competition_masked": masked_obs,
                     "competition_masked_robust": bool(masked_obs and maskb / n_boot >= 0.90
                                                       and dperm_p < 0.05)})
    ev = pd.DataFrame(rows)
    ev.to_csv(out_dir / "decoupling_validation_eval.tsv", sep="\t", index=False)
    tested = ev[ev["n"] >= MIN_N].copy()
    n_masked = int(tested["competition_masked"].sum()) if "competition_masked" in tested else 0
    n_robust = int(tested["competition_masked_robust"].sum()) if "competition_masked_robust" in tested else 0
    exp_null = float(tested["masked_null_frac"].sum()) if "masked_null_frac" in tested else np.nan
    summ = {
        "module": "mirna_hallmark.decoupling_validation.evaluate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_boot": n_boot, "n_perm": n_perm, "n_tested": int(len(tested)),
        "competition_masked_observed": n_masked,
        "competition_masked_robust": n_robust,
        "competition_masked_expected_null": exp_null,
        "median_abs_reemergence": float(tested["compet_reemergence"].abs().median()),
        "median_ci_width": float((tested["d_ci_hi"] - tested["d_ci_lo"]).median()),
        "median_dperm_p": float(tested["dperm_p"].median()),
        "median_tumor_flat_boot_frac": float(tested["tumor_flat_boot_frac"].median()),
        "frac_edges_dperm_sig": float((tested["dperm_p"] < 0.05).mean()),
    }
    for vname in ("cn", "meth", "atac"):
        col = f"{vname}_masked_robust"
        if col in tested.columns:
            summ[f"{vname}_masked_robust"] = int(tested[col].astype("boolean").fillna(False).sum())
    (out_dir / "eval_method_manifest.json").write_text(json.dumps(summ, indent=2), encoding="utf-8")
    print("\n========== DECOUPLING MODELING — PERFORMANCE ==========")
    print(f"  competition_masked: observed {n_masked} | bootstrap-robust {n_robust} | "
          f"expected-under-null {exp_null:.2f}  <- calibration")
    print(f"  edges with significant competition effect (perm p<0.05): "
          f"{int((tested['dperm_p']<0.05).sum())}/{len(tested)} ({summ['frac_edges_dperm_sig']:.0%})")
    print(f"  median |Δ re-emergence| {summ['median_abs_reemergence']:.4f} | "
          f"median 95%CI width {summ['median_ci_width']:.3f} | median perm-p {summ['median_dperm_p']:.2f}")
    print(f"  loss-call stability: tumour stays flat in median "
          f"{summ['median_tumor_flat_boot_frac']:.0%} of bootstraps")
    rob = {v: summ.get(f"{v}_masked_robust") for v in ("cn", "meth", "atac")
           if f"{v}_masked_robust" in summ}
    if rob:
        print(f"  bootstrap-ROBUST confounder masks (CN/meth/ATAC): {rob} "
              f"(of point-estimate cnv/meth/atac flags in the main table)")
    interp = ("competition signal ABOVE null" if n_masked > exp_null + 1
              else "competition signal AT/BELOW null -> lost set is genuinely NOT competition-confounded")
    print(f"  VERDICT: {interp}")
    print("=======================================================")
    return ev


def _boot_repressive_stability(y: pd.Series, x: pd.Series, cov: pd.DataFrame,
                               n_boot: int, rng) -> Tuple[float, float, int]:
    """Bootstrap stability of the 'repressive' call: fraction of participant-resamples in which the
    composition-adjusted partial ρ stays < −0.10. Returns (point ρ, repressive-call stability, n)."""
    df = pd.concat([y.rename("y"), x.rename("x"), cov], axis=1).dropna()
    if len(df) < MIN_N:
        return np.nan, np.nan, len(df)
    Y = df["y"].values; X = df["x"].values; Cm = df[[c for c in cov.columns]].values
    m = len(df)
    base = _fast_pspear(Y, X, Cm)
    cnt = 0
    for _ in range(n_boot):
        idx = rng.integers(0, m, m)
        if _fast_pspear(Y[idx], X[idx], Cm[idx]) < -0.10:
            cnt += 1
    return base, cnt / n_boot, m


def go_up_stability(*, out_dir: Path = OUT_DIR, classes=("lost", "nat_decoupled"),
                    cov_source: str = "estimate_epi", n_boot: int = 200, seed: int = 0) -> pd.DataFrame:
    """GO UP: does the healthy(GTEx) repression call — the foundation of any decoupling — become
    more bootstrap-STABLE when read at the arm-module level than per single edge? (Effect sizes do
    not magnify with aggregation, so stability is the only lever; this measures it head-to-head.)
    Output: `decoupling_go_up_stability.tsv` + a printed edge-vs-arm comparison."""
    out_dir.mkdir(parents=True, exist_ok=True)
    from mirna_hallmark.analyses.misc.mirna_comovement import _composite_repression
    from mirna_hallmark.pressure_build import load_mirtar_edges
    rng = np.random.default_rng(seed)
    edge = pd.read_csv(STATE_DIR / "mirna_gene_edge_class.tsv", sep="\t")
    sub = edge[edge["joint_edge_class"].isin(classes)].copy()
    lost_pairs = list(map(tuple, sub[["miRNA", "gene"]].drop_duplicates().values))
    arms = sorted(set(sub["miRNA"]))
    hs = HallmarkSets.load()
    edges_full = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    targets = edges_full.groupby("miRNA")["gene"].apply(lambda s: sorted(set(s)))
    genes_needed = sorted(set().union(*[set(targets.get(a, [])) for a in arms])
                          | set(g for _, g in lost_pairs) | set(EPI_MARKERS))
    bundles = _state_bundles(genes_needed)
    states = [s for s in ("gtex", "tumor") if s in bundles]
    cov = {s: _state_covariates(s, bundles[s][1], hs, cov_source) for s in states}

    rows: List[dict] = []
    print(f"[go-up] {len(lost_pairs)} lost edges, {len(arms)} arms | n_boot={n_boot}")
    for arm, g in lost_pairs:
        rec = {"level": "edge", "unit": f"{arm}|{g}", "miRNA": arm, "gene": g}
        for s in states:
            mir, rna = bundles[s]
            if arm in mir.index and g in rna.index:
                r, stab, n = _boot_repressive_stability(rna.loc[g], mir.loc[arm], cov[s], n_boot, rng)
                rec[f"rho_{s}"] = r; rec[f"repr_stab_{s}"] = stab; rec[f"n_{s}"] = n
        rows.append(rec)
    for arm in arms:
        rec = {"level": "arm_module", "unit": arm, "miRNA": arm, "gene": ""}
        for s in states:
            mir, rna = bundles[s]
            comp = _composite_repression(rna, targets.get(arm, []))
            if not comp.empty and arm in mir.index:
                r, stab, n = _boot_repressive_stability(comp, mir.loc[arm], cov[s], n_boot, rng)
                rec[f"rho_{s}"] = r; rec[f"repr_stab_{s}"] = stab; rec[f"n_{s}"] = n
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "decoupling_go_up_stability.tsv", sep="\t", index=False)

    print("\n========== GO UP: healthy(GTEx) repression-call STABILITY ==========")
    for lvl in ("edge", "arm_module"):
        d = df[df["level"] == lvl]
        if "rho_gtex" not in d.columns:
            continue
        rep = d[pd.to_numeric(d["rho_gtex"], errors="coerce") < -0.10]   # units repressive in healthy
        stab = pd.to_numeric(rep.get("repr_stab_gtex"), errors="coerce").dropna()
        allrep = pd.to_numeric(d.get("rho_gtex"), errors="coerce")
        print(f"  {lvl:11s}: n={len(d):4d} | healthy-repressive {int((allrep<-0.10).sum()):4d} "
              f"({(allrep<-0.10).mean():.0%}) | median repressive-call stability "
              f"{stab.median() if len(stab) else float('nan'):.0%} "
              f"| frac stable≥0.9: {(stab>=0.9).mean() if len(stab) else float('nan'):.0%}")
    print("  (higher stability at arm_module => aggregation denoises the call => GO UP justified)")
    print("===================================================================")
    return df


def functional_edge_decoupling(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """FUNCTIONAL-EDGE selection: of the 'lost' edges, keep only those whose healthy(GTEx) repression
    is both STRONG and bootstrap-STABLE (where the go-up test showed the signal actually lives), then
    call decoupling only on those. Reads `decoupling_go_up_stability.tsv` (edge rows, written by
    `go_up_stability`). The gold set = strong+stable healthy repression that goes genuinely flat in
    tumour, with a stable flat call — the only edges the data can support as real decoupling events.
    Output: `functional_edge_decoupling.tsv`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    src = out_dir / "decoupling_go_up_stability.tsv"
    if not src.exists():
        raise FileNotFoundError(f"run --go-up first to produce {src}")
    d = pd.read_csv(src, sep="\t")
    e = d[d["level"] == "edge"].copy()
    for c in ["rho_gtex", "repr_stab_gtex", "rho_tumor", "repr_stab_tumor"]:
        if c in e.columns:
            e[c] = pd.to_numeric(e[c], errors="coerce")
    # selection thresholds (documented, not magic): strong+stable healthy repression
    STRONG, STABLE, FLAT, FLAT_STAB = -0.15, 0.90, -0.05, 0.10
    e["functional_healthy"] = (e["rho_gtex"] <= STRONG) & (e["repr_stab_gtex"] >= STABLE)
    e["tumor_flat"] = e["rho_tumor"] > FLAT
    e["tumor_flat_stable"] = e["repr_stab_tumor"] <= FLAT_STAB
    e["gold_decoupling"] = e["functional_healthy"] & e["tumor_flat"] & e["tumor_flat_stable"]
    e = e.sort_values(["gold_decoupling", "functional_healthy", "rho_gtex"],
                      ascending=[False, False, True])
    e.to_csv(out_dir / "functional_edge_decoupling.tsv", sep="\t", index=False)
    n_total = len(e)
    n_func = int(e["functional_healthy"].sum())
    n_gold = int(e["gold_decoupling"].sum())
    print("\n========== FUNCTIONAL-EDGE DECOUPLING (where the signal lives) ==========")
    print(f"  of {n_total} 'lost' edges: {n_func} have STRONG+STABLE healthy repression "
          f"(ρ_gtex≤{STRONG}, stability≥{STABLE:.0%})")
    print(f"  of those, {n_gold} are GOLD decoupling (also genuinely + stably flat in tumour)")
    print(f"  => the trustworthy decoupling set is {n_gold}/{n_total} edges "
          f"({n_gold/max(n_total,1):.0%}); the rest are point-estimate 'losses' the data can't support")
    if n_gold:
        cols = [c for c in ["miRNA", "gene", "rho_gtex", "repr_stab_gtex", "rho_tumor",
                            "repr_stab_tumor"] if c in e.columns]
        print("\n  GOLD decoupling edges:")
        print(e[e["gold_decoupling"]][cols].round(3).to_string(index=False))
    print("========================================================================")
    return e


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate lost miRNA->gene couplings")
    ap.add_argument("--classes", nargs="+", default=["lost", "nat_decoupled"])
    ap.add_argument("--cov-source", default="estimate_epi",
                    choices=["estimate", "estimate_epi", "metagene"],
                    help="composition covariates for tumor partial-rho refits (default: estimate_epi)")
    ap.add_argument("--evaluate", action="store_true",
                    help="run the performance evaluation (bootstrap CI + permutation null + stability)")
    ap.add_argument("--go-up", action="store_true",
                    help="edge-vs-arm-module healthy-repression call stability comparison")
    ap.add_argument("--functional-edges", action="store_true",
                    help="select strong+stable healthy-repression edges -> gold decoupling set")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=200)
    args = ap.parse_args()
    if args.functional_edges:
        functional_edge_decoupling()
    elif args.go_up:
        go_up_stability(classes=tuple(args.classes), cov_source=args.cov_source, n_boot=args.n_boot)
    elif args.evaluate:
        evaluate(classes=tuple(args.classes), cov_source=args.cov_source,
                 n_boot=args.n_boot, n_perm=args.n_perm)
    else:
        run(classes=tuple(args.classes), cov_source=args.cov_source)


if __name__ == "__main__":
    main()
