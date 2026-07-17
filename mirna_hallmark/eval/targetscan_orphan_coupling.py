"""TargetScan-orphan coupling: repression not in the miRTarBase HE graph.

Tests whether sequence-predicted miRNA→target pairs **excluded from the
high-evidence miRTarBase graph** show repression-consistent coupling in
expression data — at **gene** resolution (hub targets) and **Hallmark program**
resolution (TS-only orphan pressure aggregated over member genes).

Orphan definition (per miRNA, gene pair):
- TargetScan ``|weighted context++|`` ≥ ``MIN_TS_WEIGHT``
- **Not** a miRTarBase high-evidence edge
- Annotated: ``mirtar_any`` (any miRTarBase row), ``n_studies``, ``evidence_score``

Program pressure uses top ``TOP_ORPHAN_ARMS_PER_GENE_PRESSURE`` orphan arms per gene
(log1p(ts_weight) × z(miRNA), AGO-gated), then mean over Hallmark members.

Specificity: partial Spearman within each PAM50 subtype (cohort + LumA/LumB/Her2/Basal);
summary tables report neg-sig counts and Basal-concentration flags.

Outputs (``output/robustness/targetscan_orphan/``):
- ``orphan_gene_coupling_by_scope.tsv``
- ``orphan_gene_specificity_summary.tsv``
- ``orphan_gene_aggregate_by_scope.tsv`` — per hub × scope orphan-arm counts + survival
- ``mirtar_gap_coupling_queue.tsv`` — edges absent from miRTarBase HE that still couple in-tumor (curation gap, NOT literature novelty; ``literature_status`` spot-checks)
- ``orphan_mirna_cluster_summary.tsv`` — seed-family buckets + Basal coupling coherence
- ``orphan_mirna_coupling_matrix_basal.tsv`` — miRNA × hub ρ (Basal, CN+prolif) for heatmaps
- ``orphan_mirna_coupling_by_subtype.tsv`` — miRNA × PAM50 median ρ across hubs
- ``orphan_hub_gene_pressure_per_sample.tsv.gz`` — aggregate orphan pressure per hub gene
- ``orphan_mirna_hallmark_coupling_by_pam50.tsv`` — orphan-panel miRNA × Hallmark × PAM50
- ``orphan_mirna_hallmark_coupling_matrix_basal.tsv`` — key Hallmarks @ Basal (wide)
- ``mirtar_mirna_hallmark_coupling_by_pam50.tsv`` — miRTar HE hub-route miRNAs × Hallmark
- ``mirtar_mirna_hallmark_coupling_matrix_basal.tsv`` — key Hallmarks @ Basal (wide)
- ``orphan_hallmark_coupling_by_pam50.tsv``
- ``orphan_hallmark_specificity_summary.tsv``
- ``orphan_vs_mirtar_coupling_comparison.tsv``
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.robustness_checks import (
    AIM1_DIR,
    HUB_ROUTES,
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    ROBUSTNESS_DIR,
    _f,
    _hallmark_cn_matrix,
    _hallmark_expr_matrix,
    _load_targetscan_weights,
    _pam50_scope_iter,
    _partial_ladder,
    _proliferation_proxies,
    _r,
    _zscore_within,
)

ORPHAN_DIR = ROBUSTNESS_DIR / "targetscan_orphan"

# Sequence-strength floor (context++ sum per pair).
MIN_TS_WEIGHT = 0.25
# Gene-level screen: top orphan arms per hub target by TargetScan weight.
TOP_ORPHAN_ARMS_PER_HUB_GENE = 15
# Program pressure: cap orphan arms per gene to keep aggregation stable.
TOP_ORPHAN_ARMS_PER_GENE_PRESSURE = 5


def _he_pairs(edges: pd.DataFrame) -> Set[Tuple[str, str]]:
    he = edges.loc[edges["high_evidence"]].drop_duplicates(["miRNA", "gene"])
    return set(zip(he["miRNA"], he["gene"]))


def _mirtar_lookup(summary_path: Path) -> pd.DataFrame:
    """(miRNA, gene) -> n_studies, evidence_score from Hallmark miRTarBase summary."""
    if not summary_path.exists():
        return pd.DataFrame(columns=["miRNA", "gene", "n_studies", "evidence_score"])
    s = pd.read_csv(summary_path, sep=",", usecols=lambda c: c in {"miRNA", "gene", "n_studies", "evidence_score"})
    return s.drop_duplicates(["miRNA", "gene"])


def build_orphan_edge_table(
    genes: Sequence[str],
    he_pairs: Set[Tuple[str, str]],
    *,
    min_ts: float = MIN_TS_WEIGHT,
    top_n_per_gene: Optional[int] = None,
    mirtar: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """TargetScan pairs for ``genes`` that are not high-evidence miRTarBase edges."""
    ts = _load_targetscan_weights(list(genes))
    if ts.empty:
        return ts
    ts = ts.loc[ts["ts_weight"] >= min_ts].copy()
    ts["is_he_edge"] = [pair in he_pairs for pair in zip(ts["miRNA"], ts["gene"])]
    ts = ts.loc[~ts["is_he_edge"]]
    if mirtar is not None and not mirtar.empty:
        ts = ts.merge(mirtar, on=["miRNA", "gene"], how="left")
        ts["n_studies"] = ts["n_studies"].fillna(0).astype(int)
        ts["evidence_score"] = ts["evidence_score"].fillna(0).astype(float)
        ts["mirtar_any"] = ts["n_studies"] > 0
    else:
        ts["n_studies"] = 0
        ts["evidence_score"] = 0.0
        ts["mirtar_any"] = False
    ts = ts.sort_values(["gene", "ts_weight"], ascending=[True, False])
    if top_n_per_gene is not None:
        ts = ts.groupby("gene", as_index=False).head(int(top_n_per_gene))
    return ts.reset_index(drop=True)


def compute_orphan_gene_pressure(
    orphan_edges: pd.DataFrame,
    mirna: pd.DataFrame,
    gate: pd.Series,
) -> pd.DataFrame:
    """Gene × sample orphan pressure; uses spine weighting + TS abundance floor."""
    from mirna_hallmark.pressure_build import compute_gene_pressure

    if orphan_edges.empty:
        return pd.DataFrame()
    edges = (
        orphan_edges[["miRNA", "gene", "ts_weight"]]
        .rename(columns={"ts_weight": "evidence_score"})
        .drop_duplicates(["miRNA", "gene"])
    )
    gp = compute_gene_pressure(
        sorted(edges["gene"].unique()),
        edges=edges,
        mirna=mirna,
        abundance_floor=C.PRESSURE_TS_ABUNDANCE_FLOOR,
    )
    shared = gp.columns.intersection(gate.index)
    return gp[shared].mul(gate.reindex(shared), axis=1)


def orphan_gene_coupling_by_scope(
    rna: pd.DataFrame,
    mirna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    orphan_edges: pd.DataFrame,
    cnv: pd.DataFrame,
    *,
    prolif_key: str = "e2f_g2m",
) -> pd.DataFrame:
    """Per (hub gene, orphan miRNA, scope): partial Spearman ladder."""
    clin = clinical.set_index("participant")
    pro_all = proxies[prolif_key]
    rows: List[dict] = []

    for _, edge in orphan_edges.iterrows():
        target, m = edge["gene"], edge["miRNA"]
        if target not in rna.index or m not in mirna.index:
            continue
        y_all = rna.loc[target]
        x_all = mirna.loc[m]
        cn_row = cnv.loc[target] if target in cnv.index else None
        for scope, samples in _pam50_scope_iter(clinical):
            cols = x_all.index.intersection(y_all.index)
            if samples is not None:
                cols = cols.intersection(samples)
            cols = pd.Index(sorted(cols)).intersection(clin.index)
            if len(cols) < 20:
                continue
            ladder = _partial_ladder(y_all, x_all, clin, cols, proxies, target_cn=cn_row)
            rows.append({
                "target": target,
                "miRNA": m,
                "scope": scope,
                "ts_weight": round(float(edge["ts_weight"]), 4),
                "mirtar_any": bool(edge.get("mirtar_any", False)),
                "n_studies": int(edge.get("n_studies", 0)),
                "mirtar_evidence_score": float(edge.get("evidence_score", 0)),
                "n": int(len(cols)),
                **ladder,
            })
    return pd.DataFrame(rows)


def orphan_gene_specificity_summary(detail: pd.DataFrame) -> pd.DataFrame:
    """Per (target, miRNA): neg-sig counts by scope + Basal-concentration flag."""
    if detail.empty:
        return detail
    rho_col, p_col = "rho_e2f_g2m_CN", "p_e2f_g2m_CN"
    if rho_col not in detail.columns:
        rho_col, p_col = "rho_e2f_g2m", "p_e2f_g2m"
    g = detail.copy()
    surv_col = "survives_e2f_g2m_CN" if "survives_e2f_g2m_CN" in g.columns else "survives_e2f_g2m"
    g["neg_sig"] = (g[rho_col] < 0) & (g[p_col] < 0.05)
    g["survives"] = g[surv_col].astype(bool)
    rows: List[dict] = []
    for (target, mir), sub in g.groupby(["target", "miRNA"]):
        by_scope = sub.set_index("scope")
        basal_neg = int(by_scope.loc["Basal", "neg_sig"]) if "Basal" in by_scope.index else 0
        other_neg = sum(
            int(by_scope.loc[s, "neg_sig"])
            for s in PAM50_SUBTYPES if s != "Basal" and s in by_scope.index
        )
        rows.append({
            "target": target,
            "miRNA": mir,
            "ts_weight": sub["ts_weight"].iloc[0],
            "mirtar_any": bool(sub["mirtar_any"].iloc[0]),
            "n_studies": int(sub["n_studies"].iloc[0]),
            "neg_sig_Basal": basal_neg,
            "neg_sig_other_subtypes": other_neg,
            "survives_Basal_CN": bool(by_scope.loc["Basal", "survives"]) if "Basal" in by_scope.index else False,
            "survives_cohort_CN": bool(by_scope.loc["cohort", "survives"]) if "cohort" in by_scope.index else False,
            "basal_concentrated": bool(basal_neg and basal_neg > other_neg),
            "median_rho_Basal": _r(by_scope.loc["Basal", rho_col]) if "Basal" in by_scope.index else np.nan,
        })
    out = pd.DataFrame(rows).sort_values(["basal_concentrated", "neg_sig_Basal", "ts_weight"],
                                         ascending=[False, False, False])
    return out.reset_index(drop=True)


def _orphan_seed_family(mirna: str) -> str:
    """Coarse seed-family bucket for orphan clustering (name-based)."""
    m = mirna.lower()
    if "let-7" in m:
        return "let-7"
    if "mir-98" in m:
        return "let-7/98"
    for tag in ("mir-130", "mir-301", "mir-454", "mir-520", "mir-148", "mir-107", "mir-499"):
        if tag in m:
            return tag.replace("mir-", "miR-")
    return "other"


def orphan_gene_aggregate_by_scope(detail: pd.DataFrame) -> pd.DataFrame:
    """Per (hub target, scope): orphan-arm counts, neg-sig, CN survival."""
    if detail.empty:
        return detail
    rho_col, p_col = "rho_e2f_g2m_CN", "p_e2f_g2m_CN"
    surv_col = "survives_e2f_g2m_CN"
    if rho_col not in detail.columns:
        rho_col, p_col = "rho_e2f_g2m", "p_e2f_g2m"
        surv_col = "survives_e2f_g2m"
    rows: List[dict] = []
    for (target, scope), sub in detail.groupby(["target", "scope"]):
        neg = (sub[rho_col] < 0) & (sub[p_col] < 0.05)
        surv = sub[surv_col].astype(bool)
        best = sub.loc[neg].sort_values(rho_col).head(1) if neg.any() else sub.sort_values(rho_col).head(1)
        rows.append({
            "target": target,
            "scope": scope,
            "n_orphan_arms": int(len(sub)),
            "n_neg_sig": int(neg.sum()),
            "n_survives_CN": int(surv.sum()),
            "median_rho": _r(sub[rho_col].median()),
            "best_arm": best["miRNA"].iloc[0] if not best.empty else "",
            "best_rho": _r(best[rho_col].iloc[0]) if not best.empty else np.nan,
        })
    return pd.DataFrame(rows).sort_values(["scope", "n_survives_CN", "n_neg_sig"],
                                          ascending=[True, False, False]).reset_index(drop=True)


# Spot-checked PubMed status (2026-06) for the strongest curation-gap edges.
# IMPORTANT: absence from miRTarBase HE (``n_studies == 0``) is a *curation gap*,
# NOT literature novelty — several of these edges are luciferase-validated in other
# tissues. Genuine novelty requires per-edge manual literature review.
_LIT_SPOT_CHECK: Dict[Tuple[str, str], str] = {
    # --- published direct edges (luciferase elsewhere; miRTar-HE curation gap) ---
    ("FOXO1", "miR-183"): "published: luciferase (esophageal/endometrial; breast debated, miR-96/182 validated)",
    ("FOXO1", "miR-135a"): "published: luciferase (HCC, bladder)",
    ("CDKN1A", "miR-301a"): "published: luciferase p21 (oral/prostate)",
    ("CDKN1A", "miR-499"): "published: luciferase p21 (cardiomyocyte)",
    ("CDKN1A", "miR-130b"): "published: miR-130b->CDKN1A direct (cervical); p21 modulation prostate/CRC",
    ("TGFBR2", "miR-301a"): "published: luciferase (colorectal)",
    ("IRF1", "miR-383"): "published: luciferase (testicular EC, cholangiocarcinoma)",
    ("IRF1", "miR-17"): "published: IRF1 a confirmed miR-17~92 cluster target (review/luciferase)",
    ("IRF1", "miR-124"): "published: luciferase IRF1 3'UTR (SLE/CD4+ T cell)",
    ("PTEN", "miR-148b"): "published: luciferase (bladder, cardiomyocyte)",
    ("PTEN", "miR-454"): "published: luciferase (NSCLC)",
    ("VIM", "miR-30"): "published: luciferase miR-30c->VIM (breast/lung EMT); same seed family",
    ("BCL2L11", "miR-181b"): "published: luciferase miR-181 family->BIM (astrocyte)",
    # --- likely indirect (no direct 3'UTR edge) ---
    ("VIM", "miR-124"): "reported breast association, likely indirect (EMT), not direct 3'UTR",
    # --- genuine direct-edge gaps (no luciferase validation found) ---
    ("CDKN1A", "miR-454"): "NO direct p21 validation found; miR-454 validated on PTEN/E2F6/BTG1 — candidate",
    ("CDKN1A", "miR-301b"): "NO direct p21 validation found; miR-301b validated on NR3C2(breast)/HOXB1/TAp63 — candidate",
    ("IRF1", "miR-24"): "not directly validated (network/db only) — candidate",
    ("IRF1", "miR-519d"): "not directly validated (network/db only) — candidate",
}


def _lit_spot_check(target: str, mirna: str) -> str:
    m = mirna.lower()
    for (t, tok), note in _LIT_SPOT_CHECK.items():
        if t == target and tok.lower() in m:
            return note
    return ""


def mirtar_gap_coupling_queue(detail: pd.DataFrame) -> pd.DataFrame:
    """Edges absent from miRTarBase HE that still show in-tumor coupling.

    Uses miRTarBase HE membership as a *negative* filter on the TargetScan
    predictions: keep only (target, miRNA) pairs with ``mirtar_any == False`` /
    ``n_studies == 0`` — sequence-predicted but absent from any curated miRTarBase
    record — that *also* show CN+proliferation-adjusted repression coupling
    (neg-sig, ideally CN-surviving) in >=1 scope.

    queue = (TS-predicted) AND (zero miRTar HE studies) AND (observational coupling)

    CAVEAT (verified): miRTarBase-HE absence is a **curation gap**, NOT literature
    novelty. Spot-checks (``literature_status`` column, ``_LIT_SPOT_CHECK``) show the
    strongest edges (miR-183→FOXO1, miR-301a→CDKN1A/TGFBR2, miR-383→IRF1,
    miR-148b→PTEN) are luciferase-validated elsewhere. Each row is therefore a
    candidate for **breast-context / in-tumor confirmation**, not a new edge claim;
    true novelty needs per-edge manual literature review.
    """
    if detail.empty:
        return detail
    rho_col, p_col = "rho_e2f_g2m_CN", "p_e2f_g2m_CN"
    surv_col = "survives_e2f_g2m_CN"
    if rho_col not in detail.columns:
        rho_col, p_col = "rho_e2f_g2m", "p_e2f_g2m"
        surv_col = "survives_e2f_g2m"
    g = detail.loc[~detail["mirtar_any"].astype(bool)].copy()
    if g.empty:
        return g
    g["neg_sig"] = (g[rho_col] < 0) & (g[p_col] < 0.05)
    g["survives"] = g[surv_col].astype(bool)
    rows: List[dict] = []
    for (target, mir), sub in g.groupby(["target", "miRNA"]):
        surviving = sub.loc[sub["survives"]]
        coupling = surviving if not surviving.empty else sub.loc[sub["neg_sig"]]
        if coupling.empty:
            continue
        best = coupling.sort_values(rho_col).iloc[0]
        scopes = sorted(coupling["scope"].tolist())
        rows.append({
            "target": target,
            "miRNA": mir,
            "seed_family": _orphan_seed_family(mir),
            "ts_weight": _r(sub["ts_weight"].iloc[0]),
            "n_studies_mirtar_he": int(sub["n_studies"].iloc[0]),
            "literature_status": _lit_spot_check(target, mir),
            "best_scope": best["scope"],
            "best_rho": _r(best[rho_col]),
            "best_p": _r(best[p_col]),
            "n_scopes_coupled": int(len(scopes)),
            "scopes_coupled": ";".join(scopes),
            "survives_CN_any": bool(sub["survives"].any()),
            "basal_coupled": bool((coupling["scope"] == "Basal").any()),
        })
    if not rows:
        return pd.DataFrame(rows)
    out = pd.DataFrame(rows).sort_values(
        ["survives_CN_any", "best_rho", "ts_weight"],
        ascending=[False, True, False],
    )
    return out.reset_index(drop=True)


def orphan_mirna_cluster_summary(
    detail: pd.DataFrame,
    gene_spec: pd.DataFrame,
) -> pd.DataFrame:
    """Per orphan miRNA: seed family, hub-target breadth, Basal coupling + cluster flags."""
    if detail.empty:
        return detail
    rho_col = "rho_e2f_g2m_CN" if "rho_e2f_g2m_CN" in detail.columns else "rho_e2f_g2m"
    p_col = "p_e2f_g2m_CN" if "p_e2f_g2m_CN" in detail.columns else "p_e2f_g2m"
    basal = detail.loc[detail["scope"] == "Basal"].copy()
    basal["neg_sig"] = (basal[rho_col] < 0) & (basal[p_col] < 0.05)

    spec = gene_spec.set_index(["target", "miRNA"]) if not gene_spec.empty else None
    rows: List[dict] = []
    for mir, sub in basal.groupby("miRNA"):
        fam = _orphan_seed_family(mir)
        n_targets = sub["target"].nunique()
        n_neg = int(sub["neg_sig"].sum())
        targets_neg = sorted(sub.loc[sub["neg_sig"], "target"].unique())
        basal_conc = False
        if spec is not None:
            idx = [(t, mir) for t in sub["target"].unique()]
            flags = [spec.loc[i, "basal_concentrated"] for i in idx if i in spec.index]
            basal_conc = any(bool(x) for x in flags)
        rows.append({
            "miRNA": mir,
            "seed_family": fam,
            "n_hub_targets": int(n_targets),
            "n_basal_neg_sig": n_neg,
            "basal_neg_targets": ";".join(targets_neg),
            "median_rho_basal": _r(sub[rho_col].median()),
            "basal_concentrated_any": basal_conc,
            "zero_mirtar_studies": int(sub["n_studies"].iloc[0]) == 0,
        })
    out = pd.DataFrame(rows).sort_values(
        ["seed_family", "n_basal_neg_sig", "n_hub_targets"],
        ascending=[True, False, False],
    )

    # Family-level roll-up (same table, family rows appended with miRNA=family summary)
    fam_rows: List[dict] = []
    for fam, sub in out.groupby("seed_family"):
        fam_rows.append({
            "miRNA": f"[family:{fam}]",
            "seed_family": fam,
            "n_hub_targets": int(sub["n_hub_targets"].sum()),
            "n_basal_neg_sig": int(sub["n_basal_neg_sig"].sum()),
            "basal_neg_targets": "",
            "median_rho_basal": _r(sub["median_rho_basal"].median()),
            "basal_concentrated_any": bool(sub["basal_concentrated_any"].any()),
            "zero_mirtar_studies": bool(sub["zero_mirtar_studies"].all()),
            "n_arms_in_family": int(len(sub)),
        })
    return pd.concat([out, pd.DataFrame(fam_rows)], ignore_index=True)


def orphan_mirna_coupling_matrix(detail: pd.DataFrame, *, scope: str = "Basal") -> pd.DataFrame:
    """miRNA × hub target partial ρ for one scope (wide, for clustered heatmaps)."""
    if detail.empty:
        return pd.DataFrame()
    rho_col = "rho_e2f_g2m_CN" if "rho_e2f_g2m_CN" in detail.columns else "rho_e2f_g2m"
    sub = detail.loc[detail["scope"] == scope, ["miRNA", "target", rho_col]]
    wide = sub.pivot(index="miRNA", columns="target", values=rho_col)
    return wide.sort_index()


def orphan_mirna_coupling_by_subtype(detail: pd.DataFrame) -> pd.DataFrame:
    """miRNA × PAM50 subtype: median ρ across hub targets (CN+prolif ladder)."""
    if detail.empty:
        return detail
    rho_col = "rho_e2f_g2m_CN" if "rho_e2f_g2m_CN" in detail.columns else "rho_e2f_g2m"
    sub = detail.loc[detail["scope"].isin(PAM50_SUBTYPES)]
    agg = (
        sub.groupby(["miRNA", "scope"])[rho_col]
        .median()
        .unstack("scope")
        .reindex(columns=list(PAM50_SUBTYPES))
    )
    agg.index.name = "miRNA"
    agg["seed_family"] = [_orphan_seed_family(m) for m in agg.index]
    return agg.reset_index()


def mirna_hallmark_coupling_by_pam50(
    rna: pd.DataFrame,
    mirna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    cnv: pd.DataFrame,
    mirna_panel: Sequence[str],
    *,
    edge_class: str,
    prolif_key: str = "e2f_g2m",
) -> pd.DataFrame:
    """Per (miRNA, Hallmark, PAM50): miRNA expr vs program expr (CN+prolif+within-subtype z)."""
    from analysis.utils.common.loaders import partial_spearman

    if not mirna_panel:
        return pd.DataFrame()

    clin = clinical.set_index("participant")
    he = _hallmark_expr_matrix(rna, hs)
    hc = _hallmark_cn_matrix(cnv, hs)
    pro_all = proxies[prolif_key]
    panel = sorted(set(mirna_panel) & set(mirna.index))

    rows: List[dict] = []
    for m in panel:
        x_all = mirna.loc[m]
        for sub in PAM50_SUBTYPES:
            sub_samples = set(clin.index[clin["PAM50_final"].eq(sub)])
            for hset in he.index:
                cols = pd.Index(sorted(
                    x_all.index.intersection(he.columns).intersection(sub_samples).intersection(clin.index)
                ))
                if len(cols) < 20:
                    continue
                x = pd.to_numeric(x_all.reindex(cols), errors="coerce")
                e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
                if x.nunique(dropna=True) < 5 or e.nunique(dropna=True) < 5:
                    continue
                cn = pd.to_numeric(hc.loc[hset].reindex(cols), errors="coerce") if hset in hc.index else None
                cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
                hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
                pro = pd.to_numeric(pro_all.reindex(cols), errors="coerce")
                ez, xz = _zscore_within(e), _zscore_within(x)
                cov_base = pd.concat([cpe, hrd, pro], axis=1)
                rho_wcn, p_wcn, _ = partial_spearman(
                    ez, xz,
                    pd.concat([cpe, hrd, pro, cn], axis=1) if cn is not None else cov_base,
                )
                rows.append({
                    "miRNA": m,
                    "hallmark_set": hset,
                    "subtype": sub,
                    "edge_class": edge_class,
                    "n": int(len(cols)),
                    "key_hallmark": hset in KEY_HALLMARKS,
                    "rho_prolif_cn_wsd_adj": _r(rho_wcn),
                    "p_prolif_cn_wsd_adj": _f(p_wcn),
                })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    parts = []
    for _, g in out.groupby("subtype"):
        g = g.copy()
        valid = g["p_prolif_cn_wsd_adj"].notna()
        if valid.any():
            g.loc[valid, "q_prolif_cn_wsd_adj"] = S.bh_fdr(g.loc[valid, "p_prolif_cn_wsd_adj"].values)
        parts.append(g)
    return pd.concat(parts, ignore_index=True)


def mirna_hallmark_coupling_matrix(
    detail: pd.DataFrame,
    *,
    scope: str = "Basal",
    key_only: bool = True,
) -> pd.DataFrame:
    """miRNA × Hallmark wide ρ matrix for heatmaps."""
    if detail.empty:
        return pd.DataFrame()
    sub = detail.loc[detail["subtype"] == scope]
    if key_only and "key_hallmark" in sub.columns:
        sub = sub.loc[sub["key_hallmark"]]
    wide = sub.pivot(index="miRNA", columns="hallmark_set", values="rho_prolif_cn_wsd_adj")
    return wide.sort_index()


def orphan_hallmark_coupling_by_subtype(
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    orphan_pressure: pd.DataFrame,
    cnv: pd.DataFrame,
    *,
    prolif_key: str = "e2f_g2m",
) -> pd.DataFrame:
    """Hallmark × PAM50: orphan TS-only pressure vs program expression."""
    from analysis.utils.common.loaders import partial_spearman

    if orphan_pressure.empty:
        return pd.DataFrame()

    clin = clinical.set_index("participant")
    he = _hallmark_expr_matrix(rna, hs)
    hc = _hallmark_cn_matrix(cnv, hs)
    pro_all = proxies[prolif_key]

    rows: List[dict] = []
    for sub in PAM50_SUBTYPES:
        sub_samples = set(clin.index[clin["PAM50_final"].eq(sub)])
        for hset in orphan_pressure.index.intersection(he.index):
            cols = pd.Index(sorted(
                orphan_pressure.columns.intersection(he.columns)
                .intersection(sub_samples).intersection(clin.index)
            ))
            if len(cols) < 20:
                continue
            p = pd.to_numeric(orphan_pressure.loc[hset].reindex(cols), errors="coerce")
            if p.nunique(dropna=True) < 5:
                continue
            e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
            cn = pd.to_numeric(hc.loc[hset].reindex(cols), errors="coerce") if hset in hc.index else None
            cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
            hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
            pro = pd.to_numeric(pro_all.reindex(cols), errors="coerce")
            ez, pz = _zscore_within(e), _zscore_within(p)
            cov_base = pd.concat([cpe, hrd, pro], axis=1)
            rho_adj, p_adj, _ = partial_spearman(e, p, cov_base)
            rho_wcn, p_wcn, _ = partial_spearman(
                ez, pz,
                pd.concat([cpe, hrd, pro, cn], axis=1) if cn is not None else cov_base,
            )
            rows.append({
                "hallmark_set": hset,
                "subtype": sub,
                "n": int(len(cols)),
                "key_hallmark": hset in KEY_HALLMARKS,
                "rho_prolif_adj": _r(rho_adj),
                "p_prolif_adj": _f(p_adj),
                "rho_prolif_cn_wsd_adj": _r(rho_wcn),
                "p_prolif_cn_wsd_adj": _f(p_wcn),
            })
    out = pd.DataFrame(rows)
    for pcol, qcol in (("p_prolif_adj", "q_prolif_adj"), ("p_prolif_cn_wsd_adj", "q_prolif_cn_wsd_adj")):
        if pcol not in out.columns:
            continue
        parts = []
        for _, g in out.groupby("subtype"):
            g = g.copy()
            valid = g[pcol].notna()
            if valid.any():
                g.loc[valid, qcol] = S.bh_fdr(g.loc[valid, pcol].values)
            parts.append(g)
        out = pd.concat(parts, ignore_index=True)
    return out


def orphan_hallmark_specificity_summary(detail: pd.DataFrame) -> pd.DataFrame:
    """Per subtype: neg-sig Hallmark counts for orphan program pressure."""
    if detail.empty:
        return detail
    rho_col = "rho_prolif_cn_wsd_adj"
    q_col = "q_prolif_cn_wsd_adj"
    g = detail.assign(neg_sig=(detail[rho_col] < 0) & (detail[q_col] < 0.10))
    rows: List[dict] = []
    for subtype, sub in g.groupby("subtype"):
        rows.append({
            "subtype": subtype,
            "n_hallmarks": sub["hallmark_set"].nunique(),
            "n_neg_sig": int(sub["neg_sig"].sum()),
            "n_key_neg_sig": int((sub["neg_sig"] & sub["key_hallmark"]).sum()),
            "median_rho": round(float(sub[rho_col].median()), 4),
        })
    summ = pd.DataFrame(rows)
    summ["median_rho"] = summ["median_rho"].round(4)
    basal_row = summ.loc[summ["subtype"] == "Basal"]
    if not basal_row.empty:
        b = int(basal_row["n_neg_sig"].iloc[0])
        others = summ.loc[summ["subtype"] != "Basal", "n_neg_sig"].max()
        summ["basal_concentrated_vs_others"] = summ["subtype"].eq("Basal") & (b > others)
    return summ


def pair_subtype_membership(
    detail: pd.DataFrame,
    *,
    edge_class: str,
    rho_col: str = "rho_e2f_g2m_CN",
    p_col: str = "p_e2f_g2m_CN",
    surv_col: str = "survives_e2f_g2m_CN",
) -> pd.DataFrame:
    """Long-form (target, miRNA, scope) flags for subtype overlap / UpSet plots."""
    if detail.empty:
        return detail
    sub = detail.loc[detail["scope"].isin(PAM50_SUBTYPES)].copy()
    sub["neg_sig"] = (sub[rho_col] < 0) & (sub[p_col] < 0.05)
    sub["survives_CN"] = sub[surv_col].astype(bool)
    sub["neg_sig_survives"] = sub["neg_sig"] & sub["survives_CN"]
    sub["edge_class"] = edge_class
    sub["seed_family"] = sub["miRNA"].map(_orphan_seed_family)
    cols = [
        "target", "miRNA", "scope", "edge_class", "seed_family",
        "neg_sig", "survives_CN", "neg_sig_survives", rho_col, p_col,
    ]
    if "ts_weight" in sub.columns:
        cols.insert(5, "ts_weight")
    if "high_evidence_edge" in sub.columns:
        cols.insert(5, "high_evidence_edge")
    return sub[cols].rename(columns={rho_col: "rho", p_col: "p"})


def mirna_subtype_membership_matrix(membership: pd.DataFrame) -> pd.DataFrame:
    """miRNA × PAM50: any neg_sig+CN-surviving hub edge in that subtype."""
    if membership.empty:
        return membership
    rows: List[dict] = []
    for mir, g in membership.groupby("miRNA"):
        row = {"miRNA": mir, "seed_family": g["seed_family"].iloc[0], "edge_class": g["edge_class"].iloc[0]}
        for scope in PAM50_SUBTYPES:
            gs = g.loc[g["scope"] == scope]
            row[f"{scope}_neg_sig_survives"] = bool(gs["neg_sig_survives"].any())
            row[f"{scope}_n_targets"] = int(gs.loc[gs["neg_sig_survives"], "target"].nunique())
        active = [s for s in PAM50_SUBTYPES if row[f"{s}_neg_sig_survives"]]
        row["n_subtypes_active"] = len(active)
        row["subtypes_active"] = "+".join(active)
        row["basal_only"] = active == ["Basal"]
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["n_subtypes_active", "Basal_n_targets"],
        ascending=[False, False],
        na_position="last",
    )


def pair_subtype_membership_matrix(
    membership: pd.DataFrame,
    *,
    min_subtypes: int = 1,
) -> pd.DataFrame:
    """(target, miRNA) × PAM50 binary matrix for heatmaps."""
    if membership.empty:
        return membership
    pairs = membership[["target", "miRNA", "seed_family", "edge_class"]].drop_duplicates()
    wide = membership.pivot_table(
        index=["target", "miRNA"],
        columns="scope",
        values="neg_sig_survives",
        aggfunc="max",
        fill_value=False,
    ).reindex(columns=list(PAM50_SUBTYPES))
    wide = wide.join(pairs.set_index(["target", "miRNA"])[["seed_family", "edge_class"]])
    wide["n_subtypes"] = wide[list(PAM50_SUBTYPES)].sum(axis=1).astype(int)
    return wide.loc[wide["n_subtypes"] >= min_subtypes].sort_values("n_subtypes", ascending=False)


def subtype_overlap_counts(membership: pd.DataFrame, *, item: str = "miRNA") -> pd.DataFrame:
    """Count miRNAs (or edges) per active-subtype pattern, e.g. Basal+LumB."""
    if membership.empty:
        return membership
    if item == "edge":
        key_cols = ["target", "miRNA"]
        sub = membership.loc[membership["neg_sig_survives"]].drop_duplicates(key_cols + ["scope"])
        patterns: Dict[str, Set[Tuple[str, str]]] = {}
        for (t, m), g in sub.groupby(key_cols):
            active = tuple(s for s in PAM50_SUBTYPES if g.loc[g["scope"] == s, "neg_sig_survives"].any())
            if not active:
                continue
            patterns.setdefault("+".join(active), set()).add((t, m))
        rows = [{"pattern": k, "n_edges": len(v), "n_mirnas": len({m for _, m in v})} for k, v in patterns.items()]
    else:
        mat = mirna_subtype_membership_matrix(membership)
        patterns: Dict[str, int] = {}
        for r in mat.itertuples():
            if r.n_subtypes_active == 0:
                continue
            patterns[r.subtypes_active] = patterns.get(r.subtypes_active, 0) + 1
        rows = [{"pattern": k, "n_mirnas": v, "n_edges": np.nan} for k, v in patterns.items()]
    out = pd.DataFrame(rows).sort_values("n_mirnas" if item == "miRNA" else "n_edges", ascending=False)
    return out.reset_index(drop=True)


def subtype_pairwise_jaccard(membership: pd.DataFrame, *, item: str = "miRNA") -> pd.DataFrame:
    """4×4 Jaccard overlap of neg_sig+CN-surviving items between subtypes."""
    sets: Dict[str, Set] = {}
    for scope in PAM50_SUBTYPES:
        g = membership.loc[(membership["scope"] == scope) & membership["neg_sig_survives"]]
        if item == "edge":
            sets[scope] = set(zip(g["target"], g["miRNA"]))
        else:
            sets[scope] = set(g["miRNA"].unique())
    rows = []
    for a in PAM50_SUBTYPES:
        for b in PAM50_SUBTYPES:
            sa, sb = sets[a], sets[b]
            union = len(sa | sb)
            rows.append({
                "subtype_a": a,
                "subtype_b": b,
                "n_intersection": len(sa & sb),
                "n_union": union,
                "jaccard": round(len(sa & sb) / union, 4) if union else np.nan,
            })
    return pd.DataFrame(rows)


def orphan_vs_mirtar_comparison(
    orphan_hm: pd.DataFrame,
    mirtar_hm: pd.DataFrame,
    *,
    entities: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Side-by-side subtype coupling for key Hallmarks (orphan TS pressure vs miRTar HE)."""
    if orphan_hm.empty or mirtar_hm.empty:
        return pd.DataFrame()
    keys = set(KEY_HALLMARKS)
    if entities:
        keys = set(entities) & keys
    o = orphan_hm.loc[orphan_hm["hallmark_set"].isin(keys)].copy()
    m = mirtar_hm.loc[mirtar_hm["hallmark_set"].isin(keys)].copy()
    merged = o.merge(
        m,
        on=["hallmark_set", "subtype"],
        how="outer",
        suffixes=("_orphan", "_mirtar"),
    )
    merged["orphan_neg_sig"] = (
        (merged["rho_prolif_cn_wsd_adj_orphan"] < 0)
        & (merged["q_prolif_cn_wsd_adj_orphan"] < 0.10)
    )
    merged["mirtar_neg_sig"] = (
        (merged["rho_prolif_cn_wsd_adj_mirtar"] < 0)
        & (merged["q_prolif_cn_wsd_adj_mirtar"] < 0.10)
    )
    merged["orphan_only_neg_sig"] = merged["orphan_neg_sig"] & ~merged["mirtar_neg_sig"]
    merged["mirtar_only_neg_sig"] = merged["mirtar_neg_sig"] & ~merged["orphan_neg_sig"]
    return merged.sort_values(["hallmark_set", "subtype"]).reset_index(drop=True)


def run(*, out_dir: Path = ORPHAN_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    edges = D.load_hallmark_edges()
    he_pairs = _he_pairs(edges)
    mirtar = _mirtar_lookup(C.MIRTAR_HALLMARK_SUMMARY)
    proxies = _proliferation_proxies(rna, hs)

    hub_genes = sorted(HUB_ROUTES.keys())
    cnv = D.load_cnv_target_genes(hub_genes + sorted({g for gs in hs.sets.values() for g in gs}))

    from mirna_hallmark.ago_gate import compute_ago_gate
    from mirna_hallmark.config import AgoGateParams

    gate_df = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(include_tnrc6=False, gate_min=C.AGO_GATE.gate_min,
                             gate_k=C.AGO_GATE.gate_k, gate_midpoint=C.AGO_GATE.gate_midpoint),
    )
    gate = gate_df["ago_gate"]

    # ---- Gene level (hub targets) ---- #
    print("[targetscan_orphan] gene-level orphan edges (hub targets) ...")
    orphan_hub = build_orphan_edge_table(
        hub_genes, he_pairs, min_ts=MIN_TS_WEIGHT,
        top_n_per_gene=TOP_ORPHAN_ARMS_PER_HUB_GENE, mirtar=mirtar,
    )
    print(f"[targetscan_orphan]   {len(orphan_hub):,} orphan (target, miRNA) pairs on {len(hub_genes)} hub genes")

    print("[targetscan_orphan] gene-level coupling by scope ...")
    gene_detail = orphan_gene_coupling_by_scope(rna, mirna, clinical, proxies, orphan_hub, cnv)
    gene_detail.to_csv(out_dir / "orphan_gene_coupling_by_scope.tsv", sep="\t", index=False)

    gene_summ = orphan_gene_specificity_summary(gene_detail)
    gene_summ.to_csv(out_dir / "orphan_gene_specificity_summary.tsv", sep="\t", index=False)

    gene_agg = orphan_gene_aggregate_by_scope(gene_detail)
    gene_agg.to_csv(out_dir / "orphan_gene_aggregate_by_scope.tsv", sep="\t", index=False)

    gap_q = mirtar_gap_coupling_queue(gene_detail)
    gap_q.to_csv(out_dir / "mirtar_gap_coupling_queue.tsv", sep="\t", index=False)
    print(f"[targetscan_orphan]   {len(gap_q)} miRTar-HE-gap coupling edges queued "
          "(curation gap, NOT literature novelty)")

    mir_cluster = orphan_mirna_cluster_summary(gene_detail, gene_summ)
    mir_cluster.to_csv(out_dir / "orphan_mirna_cluster_summary.tsv", sep="\t", index=False)

    orphan_mirna_coupling_matrix(gene_detail, scope="Basal").to_csv(
        out_dir / "orphan_mirna_coupling_matrix_basal.tsv", sep="\t"
    )
    orphan_mirna_coupling_by_subtype(gene_detail).to_csv(
        out_dir / "orphan_mirna_coupling_by_subtype.tsv", sep="\t", index=False
    )

    print("[targetscan_orphan] subtype edge / miRNA overlap tables ...")
    orphan_mem = pair_subtype_membership(gene_detail, edge_class="orphan_ts")
    orphan_mem.to_csv(out_dir / "orphan_pair_subtype_membership.tsv", sep="\t", index=False)
    mirna_subtype_membership_matrix(orphan_mem).to_csv(
        out_dir / "orphan_mirna_subtype_membership.tsv", sep="\t", index=False
    )
    pair_subtype_membership_matrix(orphan_mem).to_csv(
        out_dir / "orphan_pair_subtype_membership_matrix.tsv", sep="\t"
    )
    subtype_overlap_counts(orphan_mem, item="miRNA").to_csv(
        out_dir / "orphan_mirna_subtype_pattern_counts.tsv", sep="\t", index=False
    )
    subtype_overlap_counts(orphan_mem, item="edge").to_csv(
        out_dir / "orphan_edge_subtype_pattern_counts.tsv", sep="\t", index=False
    )
    subtype_pairwise_jaccard(orphan_mem, item="miRNA").to_csv(
        out_dir / "orphan_mirna_subtype_jaccard.tsv", sep="\t", index=False
    )
    subtype_pairwise_jaccard(orphan_mem, item="edge").to_csv(
        out_dir / "orphan_edge_subtype_jaccard.tsv", sep="\t", index=False
    )

    routes_path = AIM1_DIR / "hub_route_partial_corr.tsv"
    if routes_path.exists():
        routes = pd.read_csv(routes_path, sep="\t")
        mirtar_mem = pair_subtype_membership(routes, edge_class="mirtar_he")
        mirtar_mem.to_csv(out_dir / "mirtar_pair_subtype_membership.tsv", sep="\t", index=False)
        mirna_subtype_membership_matrix(mirtar_mem).to_csv(
            out_dir / "mirtar_mirna_subtype_membership.tsv", sep="\t", index=False
        )
        pair_subtype_membership_matrix(mirtar_mem).to_csv(
            out_dir / "mirtar_pair_subtype_membership_matrix.tsv", sep="\t"
        )
        subtype_overlap_counts(mirtar_mem, item="miRNA").to_csv(
            out_dir / "mirtar_mirna_subtype_pattern_counts.tsv", sep="\t", index=False
        )
        subtype_pairwise_jaccard(mirtar_mem, item="miRNA").to_csv(
            out_dir / "mirtar_mirna_subtype_jaccard.tsv", sep="\t", index=False
        )

    orphan_panel = sorted(
        gene_summ.loc[gene_summ["neg_sig_Basal"] > 0, "miRNA"].unique()
        if not gene_summ.empty else []
    )
    mirtar_panel = sorted({m for ms in HUB_ROUTES.values() for m in ms})

    print("[targetscan_orphan] miRNA x Hallmark program coupling (orphan panel) ...")
    orphan_prog = mirna_hallmark_coupling_by_pam50(
        rna, mirna, clinical, proxies, hs, cnv, orphan_panel, edge_class="orphan_ts",
    )
    orphan_prog.to_csv(out_dir / "orphan_mirna_hallmark_coupling_by_pam50.tsv", sep="\t", index=False)
    mirna_hallmark_coupling_matrix(orphan_prog, scope="Basal").to_csv(
        out_dir / "orphan_mirna_hallmark_coupling_matrix_basal.tsv", sep="\t"
    )

    print("[targetscan_orphan] miRNA x Hallmark program coupling (miRTar HE hub panel) ...")
    mirtar_prog = mirna_hallmark_coupling_by_pam50(
        rna, mirna, clinical, proxies, hs, cnv, mirtar_panel, edge_class="mirtar_he",
    )
    mirtar_prog.to_csv(out_dir / "mirtar_mirna_hallmark_coupling_by_pam50.tsv", sep="\t", index=False)
    mirna_hallmark_coupling_matrix(mirtar_prog, scope="Basal").to_csv(
        out_dir / "mirtar_mirna_hallmark_coupling_matrix_basal.tsv", sep="\t"
    )

    gp_hub = compute_orphan_gene_pressure(orphan_hub, mirna, gate)
    if not gp_hub.empty:
        gp_hub.to_csv(out_dir / "orphan_hub_gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

    if not gene_summ.empty:
        hits = gene_summ.loc[gene_summ["survives_Basal_CN"]].head(8)
        print("[targetscan_orphan]   Basal CN-surviving orphan pairs (top):")
        for r in hits.itertuples():
            print(f"[targetscan_orphan]     {r.target} <- {r.miRNA} ts={r.ts_weight} "
                  f"studies={r.n_studies} basal_conc={r.basal_concentrated}")

    # ---- Program level ---- #
    print("[targetscan_orphan] building Hallmark orphan pressure ...")
    orphan_prog_edges = build_orphan_edge_table(
        hs.universe, he_pairs, min_ts=MIN_TS_WEIGHT,
        top_n_per_gene=TOP_ORPHAN_ARMS_PER_GENE_PRESSURE, mirtar=mirtar,
    )
    gp_orphan = compute_orphan_gene_pressure(orphan_prog_edges, mirna, gate)
    from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
    hp_orphan = hallmark_pressure_matrix(gp_orphan, hs)
    hp_orphan.to_csv(out_dir / "orphan_hallmark_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

    print("[targetscan_orphan] program-level coupling by PAM50 ...")
    hm_orphan = orphan_hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp_orphan, cnv)
    hm_orphan.to_csv(out_dir / "orphan_hallmark_coupling_by_pam50.tsv", sep="\t", index=False)

    hm_spec = orphan_hallmark_specificity_summary(hm_orphan)
    hm_spec.to_csv(out_dir / "orphan_hallmark_specificity_summary.tsv", sep="\t", index=False)
    if not hm_spec.empty:
        print("[targetscan_orphan]   orphan program neg-sig by subtype:")
        for r in hm_spec.itertuples():
            print(f"[targetscan_orphan]     {r.subtype}: {int(r.n_neg_sig)}/{int(r.n_hallmarks)} "
                  f"(key {int(r.n_key_neg_sig)}) median_rho={r.median_rho:+.3f}")

    mirtar_hm = pd.read_csv(AIM1_DIR / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t")
    cmp_tbl = orphan_vs_mirtar_comparison(hm_orphan, mirtar_hm)
    cmp_tbl.to_csv(out_dir / "orphan_vs_mirtar_coupling_comparison.tsv", sep="\t", index=False)

    (out_dir / "method_manifest.json").write_text(json.dumps({
        "module": "mirna_hallmark.eval.targetscan_orphan_coupling",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "Do TargetScan-predicted miRNA→target pairs excluded from miRTarBase HE "
                    "show repression-consistent coupling at gene and Hallmark program resolution?",
        "orphan_definition": {
            "min_ts_weight": MIN_TS_WEIGHT,
            "exclude": "miRTarBase high_evidence edge",
            "hub_gene_top_arms": TOP_ORPHAN_ARMS_PER_HUB_GENE,
            "program_top_arms_per_gene": TOP_ORPHAN_ARMS_PER_GENE_PRESSURE,
        },
        "pressure": "Σ log1p(ts_weight)×z(miRNA), AGO-gated; Hallmark = mean over members",
        "coupling": "partial Spearman | CPE+HRD+E2F/G2M; program also +CN+within-subtype z",
        "specificity": "neg_sig counts per PAM50; basal_concentrated if Basal > other subtypes",
        "n_orphan_hub_pairs": int(len(orphan_hub)),
        "n_orphan_program_edges": int(len(orphan_prog_edges)),
        "aggregate_outputs": [
            "orphan_gene_aggregate_by_scope.tsv",
            "mirtar_gap_coupling_queue.tsv",
            "orphan_mirna_cluster_summary.tsv",
            "orphan_mirna_coupling_matrix_basal.tsv",
            "orphan_mirna_coupling_by_subtype.tsv",
            "orphan_hub_gene_pressure_per_sample.tsv.gz",
            "orphan_mirna_hallmark_coupling_by_pam50.tsv",
            "orphan_mirna_hallmark_coupling_matrix_basal.tsv",
            "mirtar_mirna_hallmark_coupling_by_pam50.tsv",
            "mirtar_mirna_hallmark_coupling_matrix_basal.tsv",
        ],
    }, indent=2), encoding="utf-8")

    print(f"[targetscan_orphan] done. Outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=ORPHAN_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
