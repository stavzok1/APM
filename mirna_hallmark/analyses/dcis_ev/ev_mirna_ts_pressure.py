"""TargetScan + miRTar pressure bridge for EV miRNA cohort arms.

Shared spine used by per-cohort screen modules (GSE270497, GSE255660, …):
  1. Top EV-discriminative arms (AUC or Mann-Whitney delta)
  2. TCGA ``target_combined_anticorr`` (miRTar HE)
  3. TargetScan top predicted targets + cohort partial anticorr (orphan vs curated)

Run GSE270497 bridge:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_gse270497_ts
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.dcis_ev.ev_mirna_replication import _anticorr_summary_for_arms
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.robustness_checks import _load_targetscan_weights, _partial_ladder, _proliferation_proxies
from pipeline.miRNA_exp.loaders import _arm_label_variants

TS_MIN_WEIGHT = 0.25
TS_TOP_GENES_PER_ARM = 20


def arm_variants(arm: str) -> Set[str]:
    out: Set[str] = {str(arm).strip()}
    for v in _arm_label_variants(arm):
        out.add(v)
        if v.startswith("hsa-miR-"):
            base = v.rsplit("-", 1)[0] if v.endswith(("-5p", "-3p")) else v
            out.add(base)
            out.add(base.replace("hsa-miR-", "hsa-mir-"))
    return out


def load_targetscan_for_arms(
    arms: Sequence[str],
    genes: Sequence[str],
    *,
    min_weight: float = TS_MIN_WEIGHT,
    top_per_arm: int = TS_TOP_GENES_PER_ARM,
) -> pd.DataFrame:
    """Map TargetScan rows to canonical hsa-miR-* arm labels."""
    ts = _load_targetscan_weights(list(genes))
    if ts.empty:
        return ts
    variant_to_arm: Dict[str, str] = {}
    for arm in arms:
        for v in arm_variants(arm):
            variant_to_arm[v] = arm
            variant_to_arm[v.replace("hsa-miR-", "hsa-mir-")] = arm
    ts["arm"] = ts["miRNA"].map(variant_to_arm)
    ts = ts.dropna(subset=["arm"])
    ts = ts.loc[ts["ts_weight"] >= min_weight]
    ts = ts.sort_values(["arm", "ts_weight"], ascending=[True, False])
    return ts.groupby("arm", as_index=False).head(top_per_arm)


def load_targetscan_all_for_arms(
    arms: Sequence[str],
    genes: Sequence[str],
    *,
    min_weight: float = TS_MIN_WEIGHT,
) -> pd.DataFrame:
    """All TS pairs above weight floor (no per-arm cap)."""
    ts = load_targetscan_for_arms(
        arms, genes, min_weight=min_weight, top_per_arm=10_000_000
    )
    return ts


def targetscan_cohort_anticorr(
    arms: Sequence[str],
    ts_edges: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-arm TS target coupling + per-edge detail."""
    if ts_edges.empty:
        return pd.DataFrame(), pd.DataFrame()

    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    clinical = D.load_clinical_strata()
    clin = clinical.set_index("participant")
    hs = HallmarkSets.load()
    proxies = _proliferation_proxies(rna, hs)
    edges = D.load_hallmark_edges()
    he = edges.loc[edges["high_evidence"]].drop_duplicates(["miRNA", "gene"])
    he_pairs = set(zip(he["miRNA"], he["gene"]))

    genes = sorted(set(ts_edges["gene"]))
    cnv = D.load_cnv_target_genes(genes)

    edge_rows: List[dict] = []
    for _, edge in ts_edges.iterrows():
        arm, gene = edge["arm"], edge["gene"]
        if arm not in mirna.index or gene not in rna.index:
            continue
        cols = mirna.columns.intersection(rna.columns).intersection(clin.index)
        if len(cols) < 40:
            continue
        cn_row = cnv.loc[gene] if gene in cnv.index else None
        ladder = _partial_ladder(
            rna.loc[gene],
            mirna.loc[arm],
            clin,
            pd.Index(sorted(cols)),
            proxies,
            target_cn=cn_row,
        )
        edge_rows.append(
            {
                "arm": arm,
                "gene": gene,
                "ts_weight": float(edge["ts_weight"]),
                "mirtar_he": (arm, gene) in he_pairs,
                "ts_orphan": (arm, gene) not in he_pairs,
                **ladder,
            }
        )
    detail = pd.DataFrame(edge_rows)
    if detail.empty:
        return pd.DataFrame(), detail

    rho_col = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in detail.columns else "rho_CPE_HRD"
    p_col = "p_CPE_HRD_CN" if "p_CPE_HRD_CN" in detail.columns else "p_CPE_HRD"
    detail["neg_sig_cohort"] = (detail[rho_col] < 0) & (detail[p_col] < 0.05)

    summary_rows: List[dict] = []
    for arm, sub in detail.groupby("arm"):
        neg = sub.loc[sub["neg_sig_cohort"]]
        orphan_neg = neg.loc[neg["ts_orphan"]]
        top_neg = neg.reindex(neg[rho_col].abs().sort_values(ascending=False).index).head(5)
        summary_rows.append(
            {
                "arm": arm,
                "n_ts_top_targets": len(sub),
                "n_ts_neg_sig_cohort": len(neg),
                "n_ts_orphan_neg_sig": len(orphan_neg),
                "median_ts_partial_rho": round(float(neg[rho_col].median()), 4) if len(neg) else np.nan,
                "top_ts_neg_sig_genes": ";".join(top_neg["gene"].astype(str).tolist()),
                "top_ts_orphan_neg_sig": ";".join(orphan_neg["gene"].astype(str).head(5).tolist()),
            }
        )
    return pd.DataFrame(summary_rows), detail


def build_ev_pressure_bridge(
    ev_table: pd.DataFrame,
    arms: Sequence[str],
    *,
    contrast_col: str = "contrast",
    auc_col: str = "discriminative_auc",
    delta_col: str = "delta_case_minus_control",
    direction_col: str = "direction",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Join miRTar anticorr + TargetScan summary for a set of EV arms."""
    arm_list = sorted(set(arms))
    mirtar = _anticorr_summary_for_arms(arm_list)
    if not mirtar.empty:
        mirtar = mirtar.rename(
            columns={
                "n_anticorr_gene_x_pam50": "mirtar_n_anticorr_pam50",
                "n_ecm_anticorr_hits": "mirtar_n_ecm_anticorr",
                "median_abs_partial_rho": "mirtar_median_abs_rho",
            }
        )

    rna = D.load_rna()
    ts_edges = load_targetscan_for_arms(arm_list, rna.index)
    ts_summary, ts_detail = targetscan_cohort_anticorr(arm_list, ts_edges)

    sub = ev_table.loc[ev_table["arm"].isin(arm_list)].copy()
    if contrast_col in sub.columns:
        wide_auc = sub.pivot_table(
            index="arm",
            columns=contrast_col,
            values=[auc_col, delta_col, direction_col],
            aggfunc="first",
        )
        wide_auc.columns = [f"{a}_{b}" for a, b in wide_auc.columns]
        wide_auc = wide_auc.reset_index()
        bridge = wide_auc
    else:
        bridge = sub.drop_duplicates("arm")[["arm"]].copy()
        for c in (auc_col, delta_col, direction_col):
            if c in sub.columns:
                bridge = bridge.merge(sub[["arm", c]].drop_duplicates("arm"), on="arm", how="left")

    bridge = bridge.merge(mirtar, on="arm", how="left")
    if not ts_summary.empty:
        bridge = bridge.merge(ts_summary, on="arm", how="left")

    if auc_col in sub.columns:
        bridge["best_discriminative_auc"] = (
            sub.groupby("arm")[auc_col].max().reindex(bridge["arm"]).values
        )
    bridge["pressure_coherence"] = (
        (bridge.get("mirtar_n_anticorr_pam50", 0).fillna(0) > 0)
        | (bridge.get("n_ts_neg_sig_cohort", 0).fillna(0) >= 3)
    )
    bridge["ts_orphan_story"] = bridge.get("n_ts_orphan_neg_sig", 0).fillna(0) >= 2
    bridge["discovery_score"] = (
        bridge.get("best_discriminative_auc", bridge.get(auc_col, 0.5)).fillna(0.5) * 2
        + np.log1p(bridge.get("mirtar_n_anticorr_pam50", 0).fillna(0))
        + np.log1p(bridge.get("n_ts_orphan_neg_sig", 0).fillna(0))
    )
    return bridge.sort_values("discovery_score", ascending=False), ts_detail


def coherent_story_table(
    bridge: pd.DataFrame,
    ts_detail: pd.DataFrame,
    ev_table: pd.DataFrame,
    *,
    contrast_col: str = "contrast",
    auc_col: str = "discriminative_auc",
    head: int = 40,
) -> pd.DataFrame:
    if bridge.empty:
        return bridge
    rho_col = "rho_CPE_HRD_CN" if not ts_detail.empty and "rho_CPE_HRD_CN" in ts_detail.columns else "rho_CPE_HRD"
    rows: List[dict] = []
    for _, r in bridge.head(head).iterrows():
        arm = r["arm"]
        sub = ts_detail[(ts_detail["arm"] == arm) & ts_detail["neg_sig_cohort"]].copy() if not ts_detail.empty else pd.DataFrame()
        orphans = sub.loc[sub["ts_orphan"]].sort_values(rho_col).head(5) if not sub.empty else sub
        best = ev_table.loc[ev_table["arm"] == arm].sort_values(auc_col, ascending=False).head(1)
        rows.append(
            {
                "arm": arm,
                "best_contrast": best[contrast_col].iloc[0] if contrast_col in best.columns and len(best) else "",
                "best_discriminative_auc": r.get("best_discriminative_auc"),
                "mirtar_n_anticorr_pam50": r.get("mirtar_n_anticorr_pam50"),
                "n_ts_orphan_neg_sig": r.get("n_ts_orphan_neg_sig"),
                "top_ts_orphan_targets": ";".join(orphans["gene"].astype(str).tolist()) if len(orphans) else "",
                "ts_orphan_story": bool(r.get("ts_orphan_story")),
                "pressure_coherence": bool(r.get("pressure_coherence")),
                "note": (
                    "TS-predicted orphan coupling"
                    if r.get("n_ts_orphan_neg_sig", 0) >= 2 and (r.get("mirtar_n_anticorr_pam50") or 0) < 3
                    else "miRTar pressure-supported"
                    if (r.get("mirtar_n_anticorr_pam50") or 0) >= 5
                    else "EV-only / weak pressure"
                ),
            }
        )
    return pd.DataFrame(rows)
