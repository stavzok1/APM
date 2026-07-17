"""Robustness checks for the basal miRNA-hub claims (grant Aims 1-2).

This module pre-empts the two strongest reviewer objections to the
``mirna_hallmark`` basal-hub story:

**Aim 1 - proliferation confound.** The hub *is* the proliferation machinery:
miR-17-5p / miR-20a-5p are the MYC/E2F-driven miR-17~92 cluster (OncomiR-1) and
miR-106b-5p its paralog miR-106b~25. Basal is the most proliferative PAM50
subtype, so elevated hub-miRNA AND depressed PTEN/FOXO1/p21/BIM could both be
downstream of a proliferative program without the miRNAs doing the repressing.
We therefore add a **proliferation covariate** (mean per-sample z of the
HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT metagene; MKI67 reported as a
sensitivity check) to the confounder ladder and ask whether the
miRNA -> target negative coupling survives, cohort-wide and within Basal.

**Aim 2 - curation-bias circularity.** An evidence-weighted (study-count)
miRTarBase graph upweights the most-studied oncomiRs (miR-17~92, miR-21, miR-155),
so "the hub is the famous oncomiRs" is partly baked into the prior. We recompute
each miRNA's basal regulatory load under weightings that do **not** depend on
study count -- a **binary / degree** weighting and a **sequence-based TargetScan
weighted context++** weighting -- and test whether the hub stays on top.

Outputs (``output/robustness/``):
- ``aim1_proliferation/hub_route_partial_corr.tsv``  per (target, miRNA, scope):
  raw / CPE+HRD / CPE+HRD+target_CN / +proliferation partial Spearman of target expr vs miRNA.
  Scopes: cohort + each PAM50 subtype (``LumA``, ``LumB``, ``Her2``, ``Basal``).
- ``aim1_proliferation/gene_pressure_route_partial_corr.tsv`` per hub target: target
  expr vs **combined** evidence-weighted miRNA pressure (+ CN ladder).
- ``aim1_proliferation/hallmark_coupling_prolif.tsv`` per Hallmark, cohort+each PAM50:
  pressure<->target-expression partial rho under the confounder ladder.
- ``aim2_circularity/mirna_basal_load_by_weighting.tsv`` per miRNA basal load +
  rank under evidence / binary / TargetScan weighting.
- ``aim2_circularity/hub_rank_summary.tsv`` hub-miRNA ranks + rank-correlation.
- ``aim2_circularity/per_target_top_regulators.tsv`` top predicted regulators of
  each hub target under TargetScan (sequence) weighting.
- ``method_manifest.json`` in each.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.hallmark_sets import HallmarkSets


ROBUSTNESS_DIR = C.OUTPUT_ROOT / "robustness"
AIM1_DIR = ROBUSTNESS_DIR / "aim1_proliferation"
AIM2_DIR = ROBUSTNESS_DIR / "aim2_circularity"

TARGETSCAN_CONTEXT = C.REPO_ROOT / "data" / "miRNA" / "Predicted_Targets_Context_Scores.default_predictions.txt"

# Proliferation metagene = union of these Hallmark programs (canonical cell-cycle).
PROLIFERATION_SETS = ("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
# Gene sets for PC1 proliferation axis (cell-cycle Hallmarks; union used as gene panel).
PROLIF_PC1_SETS = (
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
)
# Gene sets used to build an E2F/MYC-INDEPENDENT proliferation proxy (over-adjustment
# control): start from structural mitotic / cell-cycle programs and strip every gene
# that is also an E2F or MYC direct-target readout (those share the hub's upstream
# driver). MITOTIC_SPINDLE + G2M capture mitosis without being an E2F-target list.
ORTHO_PROLIF_BASE_SETS = ("HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE")
ORTHO_PROLIF_STRIP_SETS = (
    "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
)

# The compact basal hub under scrutiny.
HUB_MIRNAS = ("hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p", "hsa-miR-23a-3p")

# Top basal-load miRNA families for the per-family basal-specificity decomposition.
# Mature arms grouped by polycistron / seed family. We test whether each family's
# pressure->repression coupling is *basal-restricted* (like the lead cluster) or
# pan-subtype, distinguishing "basal-enriched expression" (the load metric) from
# "basal-specific coupling" (the Fig. 1A property). Order = descending basal load.
BASAL_MIRNA_FAMILIES: Dict[str, List[str]] = {
    "miR-17~92/106b~25/106a~363": [
        "hsa-miR-17-5p", "hsa-miR-18a-5p", "hsa-miR-19a-3p", "hsa-miR-19b-3p",
        "hsa-miR-20a-5p", "hsa-miR-92a-3p", "hsa-miR-17-3p", "hsa-miR-20a-3p",
        "hsa-miR-106b-5p", "hsa-miR-93-5p", "hsa-miR-25-3p", "hsa-miR-106b-3p",
        "hsa-miR-106a-5p", "hsa-miR-18b-5p", "hsa-miR-20b-5p", "hsa-miR-363-3p",
        "hsa-miR-19b-1-5p", "hsa-miR-92a-1-5p",
    ],
    "miR-23~27~24": [
        "hsa-miR-23a-3p", "hsa-miR-27a-3p", "hsa-miR-24-3p", "hsa-miR-23b-3p",
        "hsa-miR-27b-3p", "hsa-miR-24-2-5p", "hsa-miR-24-1-5p", "hsa-miR-23a-5p",
        "hsa-miR-27a-5p",
    ],
    "miR-221/222": [
        "hsa-miR-221-3p", "hsa-miR-222-3p", "hsa-miR-221-5p", "hsa-miR-222-5p",
    ],
}

# Hub target routes (target gene -> the hub miRNAs we test as its regulators).
# Restricted at runtime to (target, miRNA) pairs that actually carry an edge.
HUB_ROUTES: Dict[str, List[str]] = {
    "PTEN": ["hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p"],
    "FOXO1": ["hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p"],
    "CDKN1A": ["hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p"],
    "BCL2L11": ["hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p"],
    "IRF1": ["hsa-miR-23a-3p", "hsa-miR-106b-5p"],
    "VIM": ["hsa-miR-17-5p", "hsa-miR-20a-5p"],
    "SDC4": ["hsa-miR-17-5p", "hsa-miR-20a-5p"],
    "TGFBR2": ["hsa-miR-17-5p", "hsa-miR-20a-5p"],
}

# Programs to headline in the Hallmark-level coupling robustness.
KEY_HALLMARKS = (
    "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
)


# --------------------------------------------------------------------------- #
# Shared inputs
# --------------------------------------------------------------------------- #
def _metagene(rna: pd.DataFrame, genes) -> pd.Series:
    present = [g for g in genes if g in rna.index]
    z = S.zscore_rows(rna.loc[present])
    return z.mean(axis=0)


def _prolif_pc1(rna: pd.DataFrame, hs: HallmarkSets, ref: pd.Series) -> pd.Series:
    """PC1 over z-scored cell-cycle Hallmark genes (samples × genes SVD)."""
    genes = sorted({g for s in PROLIF_PC1_SETS for g in hs.sets.get(s, [])})
    present = [g for g in genes if g in rna.index]
    if len(present) < 10:
        return pd.Series(dtype=float, name="prolif_pc1")
    z = S.zscore_rows(rna.loc[present]).T  # samples × genes
    x = z.values.astype(float)
    x = x - x.mean(axis=0, keepdims=True)
    if x.shape[0] < 3:
        return pd.Series(dtype=float, name="prolif_pc1")
    _u, _s, vt = np.linalg.svd(x, full_matrices=False)
    pc1 = x @ vt[0]
    out = pd.Series(pc1, index=z.index, name="prolif_pc1")
    # Orient so PC1 aligns with e2f_g2m (positive = proliferative).
    aligned = ref.reindex(out.index)
    if aligned.notna().sum() >= 10 and out.corr(aligned) < 0:
        out = -out
    return out


def _proliferation_proxies(rna: pd.DataFrame, hs: HallmarkSets) -> Dict[str, pd.Series]:
    """Three proliferation proxies, deliberately differing in E2F/MYC dependence.

    - ``e2f_g2m``  : mean z of E2F_TARGETS + G2M_CHECKPOINT (the original, and the
      one a reviewer attacks: the hub miR-17~92 / miR-106b~25 are E2F/MYC targets,
      so partialling this can over-adjust by closing a shared-E2F path).
    - ``mki67``    : MKI67 alone (the standard single-marker clinical proxy).
    - ``ortho_noE2F_MYC`` : mean z of G2M + MITOTIC_SPINDLE genes **after removing
      every E2F_TARGETS / MYC_TARGETS member** — proliferation read out by genes
      that are *largely* E2F/MYC-independent (G2M genes are still partly
      E2F-influenced, so MKI67-alone is the genuinely independent proxy; this one
      is the metagene-scale cross-check). The over-adjustment objection only
      survives if the strengthening is ABSENT under the independent proxies.
    - ``pc1`` : PC1 of z-scored union of E2F + G2M + MITOTIC_SPINDLE + MYC Hallmarks.
    """
    proxies: Dict[str, pd.Series] = {}
    proxies["e2f_g2m"] = _metagene(
        rna, sorted({g for s in PROLIFERATION_SETS for g in hs.sets.get(s, [])})
    ).rename("prolif_e2f_g2m")
    if "MKI67" in rna.index:
        proxies["mki67"] = S.zscore_rows(rna.loc[["MKI67"]]).iloc[0].rename("prolif_mki67")
    strip = {g for s in ORTHO_PROLIF_STRIP_SETS for g in hs.sets.get(s, [])}
    base = {g for s in ORTHO_PROLIF_BASE_SETS for g in hs.sets.get(s, [])}
    ortho_genes = sorted(base - strip)
    proxies["ortho_noE2F_MYC"] = _metagene(rna, ortho_genes).rename("prolif_ortho_noE2F_MYC")
    proxies["pc1"] = _prolif_pc1(rna, hs, proxies["e2f_g2m"])
    print(f"[robustness]   proliferation proxies: e2f_g2m, mki67, "
          f"ortho_noE2F_MYC (n={len(ortho_genes)} genes after stripping E2F/MYC), pc1")
    return proxies


PAM50_SUBTYPES = ("LumA", "LumB", "Her2", "Basal")


def _basal_samples(clinical: pd.DataFrame) -> List[str]:
    pam = clinical.set_index("participant")["PAM50_final"]
    return pam.index[pam.eq("Basal")].tolist()


def _pam50_scope_iter(clinical: pd.DataFrame) -> List[Tuple[str, Optional[set]]]:
    """Analysis scopes for cohort + each PAM50 subtype (sample sets for filtering)."""
    clin = clinical.set_index("participant")
    pam = clin["PAM50_final"]
    scopes: List[Tuple[str, Optional[set]]] = [("cohort", None)]
    for sub in PAM50_SUBTYPES:
        scopes.append((sub, set(pam.index[pam.eq(sub)])))
    return scopes


def _zscore_within(s: pd.Series) -> pd.Series:
    """Z-score a vector over the provided sample index (dynamic-range control)."""
    v = pd.to_numeric(s, errors="coerce")
    sd = v.std()
    if pd.isna(sd) or sd == 0:
        return v * 0.0
    return (v - v.mean()) / sd


def _partial_batch(y, x, cov):
    """partial_spearman with TCGA sequencing-batch covariates appended to ``cov``
    (config.TCGA_BATCH_KIND; no-op when 'none'). Single chokepoint so every robustness
    partial gets the same batch adjustment as the rest of the spine."""
    from analysis.utils.common.loaders import partial_spearman as _ps_raw

    if cov is not None:
        cov = D.augment_tcga_batch(cov)
    return _ps_raw(y, x, cov)


def _build_cov(
    clin: pd.DataFrame,
    cols: pd.Index,
    *,
    target_cn: Optional[pd.Series] = None,
    extra: Optional[pd.Series] = None,
) -> pd.DataFrame:
    parts = [
        pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce").rename("CPE"),
        pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce").rename("HRD"),
    ]
    if target_cn is not None:
        parts.append(pd.to_numeric(target_cn.reindex(cols), errors="coerce").rename("target_CN"))
    if extra is not None:
        parts.append(pd.to_numeric(extra.reindex(cols), errors="coerce").rename(extra.name or "extra"))
    return pd.concat(parts, axis=1)


def _hallmark_cn_matrix(cnv: pd.DataFrame, hs: HallmarkSets) -> pd.DataFrame:
    """Hallmark x sample = mean ASCAT3 integer CN over present member genes."""
    rows = {}
    for hset, members in hs.sets.items():
        present = [g for g in members if g in cnv.index]
        if present:
            rows[hset] = cnv.loc[present].mean(axis=0)
    return pd.DataFrame(rows).T


def _partial_ladder(
    y: pd.Series,
    x: pd.Series,
    clin: pd.DataFrame,
    cols: pd.Index,
    proxies: Dict[str, pd.Series],
    *,
    target_cn: Optional[pd.Series] = None,
    zscore_within_scope: bool = False,
) -> dict:
    """Partial-Spearman ladder: CPE+HRD, +CN, +prolif (+CN), optional within-scope z."""
    from analysis.utils.common.loaders import partial_spearman
    from scipy.stats import spearmanr

    yv = pd.to_numeric(y.reindex(cols), errors="coerce")
    xv = pd.to_numeric(x.reindex(cols), errors="coerce")
    if zscore_within_scope:
        yv = _zscore_within(yv)
        xv = _zscore_within(xv)

    raw_rho, raw_p = spearmanr(*_drop2(yv, xv))
    cov_ch = _build_cov(clin, cols)
    rho_ch, p_ch, _ = _partial_batch(yv, xv, cov_ch)
    row = {
        "raw_rho": round(float(raw_rho), 4), "raw_p": float(raw_p),
        "rho_CPE_HRD": _r(rho_ch), "p_CPE_HRD": _f(p_ch),
    }
    if target_cn is not None:
        cov_cn = _build_cov(clin, cols, target_cn=target_cn)
        rho_cn, p_cn, _ = _partial_batch(yv, xv, cov_cn)
        row["rho_CPE_HRD_CN"] = _r(rho_cn)
        row["p_CPE_HRD_CN"] = _f(p_cn)
    for name, pr in proxies.items():
        pro = pr.rename(name)
        cov_p = _build_cov(clin, cols, extra=pro)
        rho_p, p_p, _ = _partial_batch(yv, xv, cov_p)
        row[f"rho_{name}"] = _r(rho_p)
        row[f"p_{name}"] = _f(p_p)
        row[f"survives_{name}"] = bool(pd.notna(rho_p) and rho_p < 0 and pd.notna(p_p) and p_p < 0.05)
        if target_cn is not None:
            cov_pc = _build_cov(clin, cols, target_cn=target_cn, extra=pro)
            rho_pc, p_pc, _ = _partial_batch(yv, xv, cov_pc)
            row[f"rho_{name}_CN"] = _r(rho_pc)
            row[f"p_{name}_CN"] = _f(p_pc)
            row[f"survives_{name}_CN"] = bool(
                pd.notna(rho_pc) and rho_pc < 0 and pd.notna(p_pc) and p_pc < 0.05
            )
    return row


# --------------------------------------------------------------------------- #
# AIM 1 - proliferation-adjusted attribution
# --------------------------------------------------------------------------- #
def _export_mirna_coupling_mats(
    routes: pd.DataFrame,
    out_dir: Path,
    *,
    prefix: str,
    rho_col: str = "rho_e2f_g2m_CN",
) -> None:
    """Write miRNA×hub and miRNA×PAM50 wide tables for notebook heatmaps."""
    if routes.empty or rho_col not in routes.columns:
        return
    basal = routes.loc[routes["scope"] == "Basal"]
    if not basal.empty:
        wide = basal.pivot(index="miRNA", columns="target", values=rho_col).sort_index()
        wide.to_csv(out_dir / f"{prefix}_mirna_coupling_matrix_basal.tsv", sep="\t")
    sub = routes.loc[routes["scope"].isin(PAM50_SUBTYPES)]
    if not sub.empty:
        agg = (
            sub.groupby(["miRNA", "scope"])[rho_col]
            .median()
            .unstack("scope")
            .reindex(columns=list(PAM50_SUBTYPES))
        )
        agg.reset_index().to_csv(out_dir / f"{prefix}_mirna_coupling_by_subtype.tsv", sep="\t", index=False)


def hub_route_partial_corr(
    rna: pd.DataFrame,
    mirna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    edges: pd.DataFrame,
    cnv: pd.DataFrame,
) -> pd.DataFrame:
    """Per (target, hub-miRNA, scope): partial Spearman under the confounder ladder.

    Tests target_expr vs miRNA_expr (Spearman; negative = repression-consistent)
    under: raw, +CPE+HRD, +CPE+HRD+target_CN, then +prolif (+CN). Cohort + each PAM50.
    """
    clin = clinical.set_index("participant")
    edge_pairs = set(zip(edges["miRNA"], edges["gene"]))
    he_pairs = set(zip(edges.loc[edges["high_evidence"], "miRNA"],
                       edges.loc[edges["high_evidence"], "gene"]))

    rows: List[dict] = []
    for target, mirs in HUB_ROUTES.items():
        if target not in rna.index:
            continue
        y_all = rna.loc[target]
        cn_row = cnv.loc[target] if target in cnv.index else None
        for m in mirs:
            if m not in mirna.index or (m, target) not in edge_pairs:
                continue
            x_all = mirna.loc[m]
            for scope, samples in _pam50_scope_iter(clinical):
                cols = x_all.index.intersection(y_all.index)
                if samples is not None:
                    cols = cols.intersection(samples)
                cols = pd.Index(sorted(cols)).intersection(clin.index)
                if len(cols) < 20:
                    continue
                ladder = _partial_ladder(
                    y_all, x_all, clin, cols, proxies, target_cn=cn_row,
                )
                rows.append({
                    "target": target, "miRNA": m, "scope": scope,
                    "high_evidence_edge": (m, target) in he_pairs, "n": int(len(cols)),
                    **ladder,
                })
    out = pd.DataFrame(rows)
    return out.sort_values(["target", "scope", "miRNA"]).reset_index(drop=True)


def gene_pressure_route_partial_corr(
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    gene_pressure_gated: pd.DataFrame,
    cnv: pd.DataFrame,
) -> pd.DataFrame:
    """Per hub target: target RNA vs **combined** evidence-weighted miRNA pressure.

    Unlike ``hub_route_partial_corr`` (one miRNA arm at a time), this uses the
    full per-gene pressure score (Σ z(miRNA)×log1p(evidence) over all high-evidence
    regulators, AGO-gated) — the natural test for "does aggregate miRNA repression
    track lower expression at this node?"
    """
    clin = clinical.set_index("participant")
    targets = sorted(set(HUB_ROUTES) & set(gene_pressure_gated.index) & set(rna.index))

    rows: List[dict] = []
    for target in targets:
        y_all = rna.loc[target]
        x_all = gene_pressure_gated.loc[target]
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
                "target": target, "scope": scope, "n": int(len(cols)),
                "predictor": "combined_gene_pressure", **ladder,
            })
    return pd.DataFrame(rows).sort_values(["target", "scope"]).reset_index(drop=True)


def _hallmark_expr_matrix(rna: pd.DataFrame, hs: HallmarkSets) -> pd.DataFrame:
    he_rows = {}
    for hset, members in hs.sets.items():
        present = [g for g in members if g in rna.index]
        if present:
            he_rows[hset] = rna.loc[present].mean(axis=0)
    return pd.DataFrame(he_rows).T


def hallmark_coupling_prolif(
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    hp_gated: pd.DataFrame,
) -> pd.DataFrame:
    """Per-Hallmark pressure<->target-expression partial rho under the ladder.

    Reuses the cached gated hallmark pressure matrix; expected rho<0 = pressure
    couples to repression. We compare CPE+HRD against CPE+HRD + each proliferation
    proxy, cohort and within each PAM50 subtype. The ``ortho_noE2F_MYC`` proxy is the
    over-adjustment control.
    """
    from analysis.utils.common.loaders import partial_spearman

    clin = clinical.set_index("participant")
    he = _hallmark_expr_matrix(rna, hs)

    rows: List[dict] = []
    for scope, samples in _pam50_scope_iter(clinical):
        for hset in hp_gated.index.intersection(he.index):
            cols = hp_gated.columns.intersection(he.columns)
            if samples is not None:
                cols = cols.intersection(samples)
            cols = pd.Index(sorted(cols)).intersection(clin.index)
            if len(cols) < 20:
                continue
            p = pd.to_numeric(hp_gated.loc[hset].reindex(cols), errors="coerce")
            e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
            cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
            hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
            base_cov = pd.concat([cpe, hrd], axis=1)
            rho_ch, p_ch, _ = _partial_batch(e, p, base_cov)
            row = {
                "hallmark_set": hset, "scope": scope, "n": int(len(cols)),
                "key_hallmark": hset in KEY_HALLMARKS,
                "rho_CPE_HRD": _r(rho_ch), "p_CPE_HRD": _f(p_ch),
            }
            for name, pr in proxies.items():
                pro = pd.to_numeric(pr.reindex(cols), errors="coerce")
                rho_p, p_p, _ = _partial_batch(e, p, pd.concat([base_cov, pro], axis=1))
                row[f"rho_{name}"] = _r(rho_p)
                row[f"p_{name}"] = _f(p_p)
                row[f"survives_{name}"] = bool(pd.notna(rho_p) and rho_p < 0 and pd.notna(p_p) and p_p < 0.05)
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["scope", "rho_ortho_noE2F_MYC"]).reset_index(drop=True)


def hallmark_coupling_by_subtype(
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    hp_gated: pd.DataFrame,
    cnv: pd.DataFrame,
    prolif_key: str = "e2f_g2m",
) -> pd.DataFrame:
    """Per (Hallmark, PAM50 subtype): proliferation-adjusted pressure↔repression rho.

    Primary readout: partial Spearman | CPE + HRD + proliferation.
    Sensitivity columns add (i) mean Hallmark-member **target CN** and (ii) **within-
    subtype z-scoring** of pressure and expression before partialing (dynamic-range
    control so Basal is not favored by larger internal variance alone).
    """
    from analysis.utils.common.loaders import partial_spearman

    clin = clinical.set_index("participant")
    pam = clin["PAM50_final"]
    he = _hallmark_expr_matrix(rna, hs)
    hc = _hallmark_cn_matrix(cnv, hs)
    pro_all = proxies[prolif_key]

    rows: List[dict] = []
    for sub in PAM50_SUBTYPES:
        sub_samples = set(pam.index[pam.eq(sub)])
        for hset in hp_gated.index.intersection(he.index):
            cols = pd.Index(sorted(
                hp_gated.columns.intersection(he.columns)
                .intersection(sub_samples).intersection(clin.index)
            ))
            if len(cols) < 20:
                rows.append({"hallmark_set": hset, "subtype": sub, "n": int(len(cols)),
                             "rho_CPE_HRD": np.nan, "rho_prolif_adj": np.nan, "p_prolif_adj": np.nan})
                continue
            p = pd.to_numeric(hp_gated.loc[hset].reindex(cols), errors="coerce")
            e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
            cn = pd.to_numeric(hc.loc[hset].reindex(cols), errors="coerce") if hset in hc.index else None
            cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
            hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
            pro = pd.to_numeric(pro_all.reindex(cols), errors="coerce")
            rho_ch, _p_ch, _ = _partial_batch(e, p, pd.concat([cpe, hrd], axis=1))
            rho_adj, p_adj, _ = _partial_batch(e, p, pd.concat([cpe, hrd, pro], axis=1))
            row = {
                "hallmark_set": hset, "subtype": sub, "n": int(len(cols)),
                "key_hallmark": hset in KEY_HALLMARKS,
                "rho_CPE_HRD": _r(rho_ch),
                "rho_prolif_adj": _r(rho_adj), "p_prolif_adj": _f(p_adj),
            }
            if cn is not None:
                rho_cn, p_cn, _ = _partial_batch(e, p, pd.concat([cpe, hrd, pro, cn], axis=1))
                row["rho_prolif_cn_adj"] = _r(rho_cn)
                row["p_prolif_cn_adj"] = _f(p_cn)
            ez, pz = _zscore_within(e), _zscore_within(p)
            rho_wsd, p_wsd, _ = _partial_batch(
                ez, pz, pd.concat([cpe, hrd, pro], axis=1)
            )
            row["rho_prolif_wsd_adj"] = _r(rho_wsd)
            row["p_prolif_wsd_adj"] = _f(p_wsd)
            if cn is not None:
                rho_wcn, p_wcn, _ = _partial_batch(
                    ez, pz, pd.concat([cpe, hrd, pro, cn], axis=1)
                )
                row["rho_prolif_cn_wsd_adj"] = _r(rho_wcn)
                row["p_prolif_cn_wsd_adj"] = _f(p_wcn)
            rows.append(row)
    out = pd.DataFrame(rows)
    for pcol, qcol in (
        ("p_prolif_adj", "q_prolif_adj"),
        ("p_prolif_cn_adj", "q_prolif_cn_adj"),
        ("p_prolif_wsd_adj", "q_prolif_wsd_adj"),
        ("p_prolif_cn_wsd_adj", "q_prolif_cn_wsd_adj"),
    ):
        if pcol not in out.columns:
            continue
        parts = []
        for sub, g in out.groupby("subtype"):
            g = g.copy()
            valid = g[pcol].notna()
            if valid.any():
                g.loc[valid, qcol] = S.bh_fdr(g.loc[valid, pcol].values)
            parts.append(g)
        out = pd.concat(parts, ignore_index=True)
    return out


def _base_name(m: str) -> str:
    return m.replace("hsa-", "")


def family_load_share(load_tbl: pd.DataFrame) -> pd.DataFrame:
    """Per top-basal family: share of total positive basal load (evidence + sequence).

    Makes the grant's "compact hub" claim traceable: what fraction of the total
    basal repressive load each polycistron / canonical-basal family carries, under
    both evidence and TargetScan (sequence) weighting, plus how many of the top-10
    miRNAs by load it supplies.
    """
    load_tbl = load_tbl.copy()
    load_tbl["base"] = load_tbl["miRNA"].map(_base_name)
    fam_lookup = {_base_name(m): name for name, mirs in BASAL_MIRNA_FAMILIES.items() for m in mirs}
    tot_ev = load_tbl.loc[load_tbl["load_evidence"] > 0, "load_evidence"].sum()
    tot_ts = load_tbl.loc[load_tbl["load_targetscan"] > 0, "load_targetscan"].sum()
    top10 = set(load_tbl.nlargest(10, "load_evidence")["base"])
    top20 = set(load_tbl.nlargest(20, "load_evidence")["base"])
    rows: List[dict] = []
    for name, mirs in BASAL_MIRNA_FAMILIES.items():
        bases = {_base_name(m) for m in mirs}
        sub = load_tbl[load_tbl["base"].isin(bases)]
        ev = sub.loc[sub["load_evidence"] > 0, "load_evidence"].sum()
        ts = sub.loc[sub["load_targetscan"] > 0, "load_targetscan"].sum()
        rows.append({
            "family": name, "n_arms_present": int((sub["load_evidence"] != 0).sum()),
            "load_evidence": round(float(ev), 2),
            "share_evidence_pct": round(100 * ev / tot_ev, 1) if tot_ev else np.nan,
            "load_targetscan": round(float(ts), 2),
            "share_targetscan_pct": round(100 * ts / tot_ts, 1) if tot_ts else np.nan,
            "n_in_top10_by_load": len(bases & top10),
            "n_in_top20_by_load": len(bases & top20),
        })
    out = pd.DataFrame(rows).sort_values("load_evidence", ascending=False).reset_index(drop=True)
    out.attrs["total_positive_load_evidence"] = round(float(tot_ev), 1)
    out.attrs["total_positive_load_targetscan"] = round(float(tot_ts), 1)
    return out


def _family_hallmark_pressure(
    family_mirnas: Sequence[str], hs: HallmarkSets, gate: pd.Series
) -> pd.DataFrame:
    """Gated Hallmark x sample pressure restricted to one miRNA family's edges.

    Identical construction to the main pipeline pressure (compute_pressure =
    evidence-weighted Σ z(miRNA), then AGO/RISC gate, then mean over a Hallmark's
    member genes) but with edges filtered to ``family_mirnas`` only.
    """
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges

    edges = load_mirtar_edges(hs.universe)
    fam = set(family_mirnas)
    edges = edges.loc[edges["miRNA"].isin(fam)]
    if edges.empty:
        return pd.DataFrame()
    mirna = D.load_mirna_arms()
    gp = compute_gene_pressure(hs.universe, edges=edges, mirna=mirna)
    shared = gp.columns.intersection(gate.index)
    gp = gp[shared].mul(gate.reindex(shared), axis=1)
    return hallmark_pressure_matrix(gp, hs)


def family_coupling_by_subtype(
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    gate: pd.Series,
    cnv: pd.DataFrame,
    prolif_key: str = "e2f_g2m",
) -> pd.DataFrame:
    """Per (family, PAM50 subtype, Hallmark): proliferation-adjusted pressure↔repression rho.

    Builds a *family-restricted* gated Hallmark pressure and asks, within each
    subtype, whether that family's pressure couples to target repression
    (partial Spearman | CPE+HRD+proliferation; rho<0 = repression). FDR within each
    (family, subtype). This separates **basal-enriched expression** (the load metric)
    from **basal-specific coupling**: a family whose negative, FDR-significant
    coupling concentrates in Basal is basal-specific like the lead cluster; one that
    couples across subtypes (or not at all) is not.
    """
    from analysis.utils.common.loaders import partial_spearman

    clin = clinical.set_index("participant")
    pam = clin["PAM50_final"]
    he = _hallmark_expr_matrix(rna, hs)
    hc = _hallmark_cn_matrix(cnv, hs)
    pro_all = proxies[prolif_key]

    rows: List[dict] = []
    for fam_name, mirs in BASAL_MIRNA_FAMILIES.items():
        hp_fam = _family_hallmark_pressure(mirs, hs, gate)
        if hp_fam.empty:
            continue
        for sub in PAM50_SUBTYPES:
            sub_samples = set(pam.index[pam.eq(sub)])
            fam_rows: List[dict] = []
            for hset in hp_fam.index.intersection(he.index):
                cols = pd.Index(sorted(
                    hp_fam.columns.intersection(he.columns)
                    .intersection(sub_samples).intersection(clin.index)
                ))
                if len(cols) < 20:
                    continue
                p = pd.to_numeric(hp_fam.loc[hset].reindex(cols), errors="coerce")
                # skip Hallmarks the family does not actually pressure (all-zero/constant)
                if p.nunique(dropna=True) < 5:
                    continue
                e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
                cn = pd.to_numeric(hc.loc[hset].reindex(cols), errors="coerce") if hset in hc.index else None
                cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
                hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
                pro = pd.to_numeric(pro_all.reindex(cols), errors="coerce")
                ez, pz = _zscore_within(e), _zscore_within(p)
                cov_base = pd.concat([cpe, hrd, pro], axis=1)
                rho_adj, p_adj, _ = _partial_batch(e, p, cov_base)
                rho_wcn, p_wcn, _ = _partial_batch(
                    ez, pz,
                    pd.concat([cpe, hrd, pro, cn], axis=1) if cn is not None else cov_base,
                )
                fam_rows.append({
                    "family": fam_name, "subtype": sub, "hallmark_set": hset,
                    "key_hallmark": hset in KEY_HALLMARKS, "n": int(len(cols)),
                    "rho_prolif_adj": _r(rho_adj), "p_prolif_adj": _f(p_adj),
                    "rho_prolif_cn_wsd_adj": _r(rho_wcn), "p_prolif_cn_wsd_adj": _f(p_wcn),
                })
            if fam_rows:
                g = pd.DataFrame(fam_rows)
                valid = g["p_prolif_adj"].notna()
                g.loc[valid, "q_prolif_adj"] = S.bh_fdr(g.loc[valid, "p_prolif_adj"].values)
                if "p_prolif_cn_wsd_adj" in g.columns:
                    v2 = g["p_prolif_cn_wsd_adj"].notna()
                    g.loc[v2, "q_prolif_cn_wsd_adj"] = S.bh_fdr(g.loc[v2, "p_prolif_cn_wsd_adj"].values)
                rows.extend(g.to_dict("records"))
    return pd.DataFrame(rows)


def family_specificity_summary(detail: pd.DataFrame) -> pd.DataFrame:
    """Collapse family_coupling_by_subtype to one row per (family, subtype).

    Reports, over all pressured Hallmarks: n negative & FDR-significant, median rho,
    and a basal-specificity flag = (Basal neg-sig count) vs the max over other
    subtypes — the quantity that says whether the family is basal-*specific* in
    coupling, not merely basal-enriched in expression.
    """
    if detail.empty:
        return detail
    rho_col = "rho_prolif_cn_wsd_adj" if "rho_prolif_cn_wsd_adj" in detail.columns else "rho_prolif_adj"
    q_col = "q_prolif_cn_wsd_adj" if "q_prolif_cn_wsd_adj" in detail.columns else "q_prolif_adj"
    g = (detail.assign(neg_sig=(detail[rho_col] < 0) & (detail[q_col] < 0.10))
         .groupby(["family", "subtype"])
         .agg(n_hallmarks=("hallmark_set", "nunique"),
              n_neg_sig=("neg_sig", "sum"),
              median_rho=(rho_col, "median"))
         .reset_index())
    g["median_rho"] = g["median_rho"].round(4)
    # order subtypes
    order = {s: i for i, s in enumerate(PAM50_SUBTYPES)}
    g = g.sort_values(["family", "subtype"], key=lambda c: c.map(order) if c.name == "subtype" else c)
    return g.reset_index(drop=True)


def prolif_sign_structure(
    rna: pd.DataFrame,
    mirna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    hp_gated: pd.DataFrame,
) -> pd.DataFrame:
    """Show the suppressor sign structure that makes adjustment strengthen the coupling.

    Reports, within Basal, Spearman(proliferation, "pressure/predictor") and
    Spearman(proliferation, target-expression) at BOTH resolutions:

    - ``level="program"`` (the 8 key Hallmarks): pressure = gated Hallmark pressure,
      target = mean member expression.
    - ``level="gene"`` (the hub routes PTEN/CDKN1A/TGFBR2/VIM/IRF1/BCL2L11): the
      "pressure/predictor" axis is the *miRNA arm's expression* and the target axis
      is the *single target gene's expression* — confirming the gene-resolved robust
      claim individually, since PTEN/p21 need not track proliferation the way a
      growth-program score does.

    A *positive* association with BOTH axes means proliferation injects a positive
    component into the predictor↔target correlation, so removing it deepens the
    negative residual — genuine suppressor confounding, not a manufactured coupling.
    """
    from scipy.stats import spearmanr

    clin = clinical.set_index("participant")
    basal = set(_basal_samples(clinical))
    he = _hallmark_expr_matrix(rna, hs)

    def _sign_row(base: dict, pred: pd.Series, tgt: pd.Series) -> dict:
        cols = pred.index.intersection(tgt.index).intersection(basal).intersection(clin.index)
        cols = pd.Index(sorted(cols))
        if len(cols) < 20:
            return {}
        p = pd.to_numeric(pred.reindex(cols), errors="coerce")
        e = pd.to_numeric(tgt.reindex(cols), errors="coerce")
        row = dict(base, n=int(len(cols)))
        for name, pr in proxies.items():
            pro = pd.to_numeric(pr.reindex(cols), errors="coerce")
            r_pp, _ = spearmanr(*_drop2(pro, p))
            r_pe, _ = spearmanr(*_drop2(pro, e))
            row[f"rho_prolif_pressure__{name}"] = _r(r_pp)
            row[f"rho_prolif_targetexpr__{name}"] = _r(r_pe)
        return row

    rows: List[dict] = []
    for hset in [h for h in KEY_HALLMARKS if h in hp_gated.index and h in he.index]:
        r = _sign_row({"level": "program", "entity": hset, "partner": "", "scope": "basal"},
                      hp_gated.loc[hset], he.loc[hset])
        if r:
            rows.append(r)
    for target, mirs in HUB_ROUTES.items():
        if target not in rna.index:
            continue
        for m in mirs:
            if m not in mirna.index:
                continue
            r = _sign_row({"level": "gene", "entity": target, "partner": m, "scope": "basal"},
                          mirna.loc[m], rna.loc[target])
            if r:
                rows.append(r)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# AIM 2 - curation-bias circularity (alternative edge weightings)
# --------------------------------------------------------------------------- #
def _load_targetscan_weights(genes: Sequence[str]) -> pd.DataFrame:
    """(miRNA, gene) -> sequence-based weight = Σ |weighted context++ score|.

    Human (`hsa-`) rows only; restricted to the Hallmark gene universe. Study-count
    independent (pure sequence/site model), so it cannot carry miRTarBase curation
    bias.
    """
    gene_set = set(genes)
    usecols = ["Gene Symbol", "miRNA", "weighted context++ score"]
    chunks = []
    for chunk in pd.read_csv(
        TARGETSCAN_CONTEXT, sep="\t", usecols=lambda c: c in usecols, chunksize=2_000_000
    ):
        chunk = chunk.loc[chunk["miRNA"].astype(str).str.startswith("hsa-")]
        chunk = chunk.loc[chunk["Gene Symbol"].isin(gene_set)]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        return pd.DataFrame(columns=["miRNA", "gene", "ts_weight"])
    ts = pd.concat(chunks, ignore_index=True)
    ts["w"] = pd.to_numeric(ts["weighted context++ score"], errors="coerce").abs()
    ts = ts.dropna(subset=["w"])
    ts = ts.groupby(["miRNA", "Gene Symbol"], as_index=False)["w"].sum()
    return ts.rename(columns={"Gene Symbol": "gene", "w": "ts_weight"})


def mirna_basal_load_by_weighting(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    clinical: pd.DataFrame,
    hs: HallmarkSets,
) -> pd.DataFrame:
    """Per-miRNA basal regulatory load under evidence / binary / TargetScan weights.

    load(m) = structural_weight(m) x basal_mean_z(m). The expression term is shared;
    only the structural weight changes (this is where curation bias lives), so the
    comparison isolates the curation question.
    """
    he = edges.loc[edges["high_evidence"]].dropna(subset=["miRNA", "gene"]).copy()
    he = he.drop_duplicates(["miRNA", "gene"])  # dedup multi-Hallmark context rows

    # evidence-weighted (current) and binary structural weights
    he["w_evidence"] = np.log1p(pd.to_numeric(he["evidence_score"], errors="coerce").fillna(0.0))
    he["w_binary"] = 1.0
    struct = he.groupby("miRNA").agg(
        sw_evidence=("w_evidence", "sum"),
        sw_binary=("w_binary", "sum"),
        n_edges=("gene", "nunique"),
    )

    # sequence-based TargetScan weight over the SAME high-evidence (miRNA,gene) pairs
    ts = _load_targetscan_weights(hs.universe)
    he_ts = he.merge(ts, on=["miRNA", "gene"], how="left")
    he_ts["ts_weight"] = he_ts["ts_weight"].fillna(0.0)
    struct["sw_targetscan"] = he_ts.groupby("miRNA")["ts_weight"].sum()
    struct["sw_targetscan"] = struct["sw_targetscan"].fillna(0.0)
    struct["ts_covered_edges"] = he_ts.assign(c=he_ts["ts_weight"] > 0).groupby("miRNA")["c"].sum()

    # basal mean z of each miRNA
    z = S.zscore_rows(mirna).fillna(0.0)
    basal = [s for s in _basal_samples(clinical) if s in z.columns]
    basal_z = z[basal].mean(axis=1) if basal else pd.Series(0.0, index=z.index)
    struct = struct.loc[struct.index.intersection(z.index)].copy()
    struct["basal_mean_z"] = struct.index.map(basal_z)

    for w in ("evidence", "binary", "targetscan"):
        struct[f"load_{w}"] = (struct[f"sw_{w}"] * struct["basal_mean_z"]).round(4)
        # rank by realized load (higher = more basal repressive load)
        struct[f"rank_{w}"] = struct[f"load_{w}"].rank(ascending=False, method="min").astype(int)

    struct = struct.reset_index().rename(columns={"index": "miRNA"})
    struct["is_hub"] = struct["miRNA"].isin(HUB_MIRNAS)
    return struct.sort_values("load_evidence", ascending=False).reset_index(drop=True)


def hub_rank_summary(load_tbl: pd.DataFrame) -> pd.DataFrame:
    """Hub-miRNA ranks under each weighting + cross-weighting rank correlations."""
    from scipy.stats import spearmanr

    rows = []
    n = len(load_tbl)
    for _, r in load_tbl.loc[load_tbl["is_hub"]].iterrows():
        rows.append({
            "miRNA": r["miRNA"], "n_miRNAs_ranked": n, "n_edges": int(r["n_edges"]),
            "ts_covered_edges": int(r.get("ts_covered_edges", 0)),
            "rank_evidence": int(r["rank_evidence"]),
            "rank_binary": int(r["rank_binary"]),
            "rank_targetscan": int(r["rank_targetscan"]),
            "load_evidence": r["load_evidence"],
            "load_binary": r["load_binary"],
            "load_targetscan": r["load_targetscan"],
        })
    hub = pd.DataFrame(rows)
    # rank-correlation across weightings (over all miRNAs with positive basal load)
    pos = load_tbl.loc[load_tbl["load_evidence"] > 0]
    corr = {}
    for a, b in (("evidence", "binary"), ("evidence", "targetscan"), ("binary", "targetscan")):
        if len(pos) > 10:
            rho, _ = spearmanr(pos[f"load_{a}"], pos[f"load_{b}"])
            corr[f"rho_{a}_vs_{b}"] = round(float(rho), 4)
    hub.attrs["rank_correlations"] = corr
    return hub


def per_target_top_regulators(edges: pd.DataFrame, hs: HallmarkSets, *, top_n: int = 6) -> pd.DataFrame:
    """Top predicted regulators of each hub target under TargetScan (sequence) weight."""
    ts = _load_targetscan_weights(list(HUB_ROUTES.keys()))
    rows = []
    he = edges.loc[edges["high_evidence"]].drop_duplicates(["miRNA", "gene"])
    he_pairs = set(zip(he["miRNA"], he["gene"]))
    for target in HUB_ROUTES:
        sub = ts.loc[ts["gene"] == target].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("ts_weight", ascending=False).head(top_n)
        for _, r in sub.iterrows():
            rows.append({
                "target": target, "miRNA": r["miRNA"],
                "ts_weight_abs_wcs": round(float(r["ts_weight"]), 4),
                "is_hub": r["miRNA"] in HUB_MIRNAS,
                "has_mirtarbase_he_edge": (r["miRNA"], target) in he_pairs,
            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _drop2(a: pd.Series, b: pd.Series):
    df = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    return df["a"].values, df["b"].values


def _r(v) -> float:
    return round(float(v), 4) if pd.notna(v) else np.nan


def _f(v) -> float:
    return float(v) if pd.notna(v) else np.nan


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = ROBUSTNESS_DIR) -> None:
    AIM1_DIR.mkdir(parents=True, exist_ok=True)
    AIM2_DIR.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    edges = D.load_hallmark_edges()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted(set(HUB_ROUTES) | {g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    print(f"[robustness] proliferation proxies over {rna.shape[1]} samples; "
          f"Basal n={len(_basal_samples(clinical))}; CNV for {cnv.shape[0]} genes")

    # ---- Aim 1 ---- #
    print("[robustness] Aim 1: hub-route partial correlations (+ target CN) ...")
    routes = hub_route_partial_corr(rna, mirna, clinical, proxies, edges, cnv)
    routes.to_csv(AIM1_DIR / "hub_route_partial_corr.tsv", sep="\t", index=False)
    _export_mirna_coupling_mats(routes, AIM1_DIR, prefix="mirtar")

    print("[robustness] Aim 1: combined gene-pressure routes (+ target CN) ...")
    from mirna_hallmark.hallmark_interaction import compute_gene_pressure
    from mirna_hallmark.ago_gate import compute_ago_gate
    from mirna_hallmark.config import AgoGateParams
    gp = compute_gene_pressure(hs.universe)
    gate_df = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(include_tnrc6=False, gate_min=C.AGO_GATE.gate_min,
                             gate_k=C.AGO_GATE.gate_k, gate_midpoint=C.AGO_GATE.gate_midpoint),
    )
    gate = gate_df["ago_gate"]
    shared = gp.columns.intersection(gate.index)
    gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)
    hub_gp = gp_gated.loc[sorted(HUB_ROUTES.keys())]
    hub_gp.to_csv(AIM1_DIR / "mirtar_hub_gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")
    gene_routes = gene_pressure_route_partial_corr(rna, clinical, proxies, gp_gated, cnv)
    gene_routes.to_csv(AIM1_DIR / "gene_pressure_route_partial_corr.tsv", sep="\t", index=False)

    print("[robustness] Aim 1: Hallmark coupling under the confounder ladder ...")
    hp_gated = pd.read_csv(
        C.INTERACTION_DIR / "hallmark_pressure_per_sample.tsv.gz", sep="\t", index_col=0
    )
    coupling = hallmark_coupling_prolif(rna, clinical, proxies, hs, hp_gated)
    coupling.to_csv(AIM1_DIR / "hallmark_coupling_prolif.tsv", sep="\t", index=False)

    sign = prolif_sign_structure(rna, mirna, clinical, proxies, hs, hp_gated)
    sign.to_csv(AIM1_DIR / "prolif_sign_structure.tsv", sep="\t", index=False)

    print("[robustness] Aim 1: proliferation-adjusted coupling by PAM50 subtype "
          "(+ Hallmark CN + within-subtype z-score) ...")
    by_sub = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp_gated, cnv)
    by_sub.to_csv(AIM1_DIR / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)

    # survival counts per proxy (within Basal), for the manifest
    cb = coupling[coupling["scope"] == "Basal"]
    survive_by_proxy = {n: f"{int(cb[f'survives_{n}'].sum())}/{len(cb)}" for n in proxies}
    survive_key_by_proxy = {
        n: f"{int(cb.loc[cb['key_hallmark'], f'survives_{n}'].sum())}/{int(cb['key_hallmark'].sum())}"
        for n in proxies
    }
    rb = routes[routes["scope"] == "Basal"]
    route_survive_by_proxy = {n: f"{int(rb[f'survives_{n}'].sum())}/{len(rb)}" for n in proxies}
    (AIM1_DIR / "method_manifest.json").write_text(json.dumps({
        "module": "mirna_hallmark.robustness_checks (Aim 1)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "Does the basal miRNA->target repression survive proliferation adjustment, "
                    "including an E2F/MYC-INDEPENDENT proliferation proxy (over-adjustment control)?",
        "proliferation_proxies": {
            "e2f_g2m": "mean z of E2F_TARGETS + G2M_CHECKPOINT (E2F/MYC-dependent; the attackable one)",
            "mki67": "MKI67 single-gene z (standard clinical proliferation marker)",
            "ortho_noE2F_MYC": "mean z of G2M_CHECKPOINT + MITOTIC_SPINDLE genes after removing all "
                               "E2F_TARGETS / MYC_TARGETS members (largely E2F/MYC-independent metagene; "
                               "MKI67-alone is the genuinely independent proxy)",
        },
        "over_adjustment_logic": "If strengthening appears ONLY under e2f_g2m it is a shared-E2F-variance "
                                 "artifact; if it also holds under ortho_noE2F_MYC the objection is closed.",
        "confounder_ladder": ["raw", "CPE+HRD", "CPE+HRD+<each proxy>"],
        "partial_method": "analysis.utils.common.loaders.partial_spearman (OLS residualize)",
        "basal_hallmarks_neg_p05_by_proxy": survive_by_proxy,
        "basal_key_hallmarks_neg_p05_by_proxy": survive_key_by_proxy,
        "basal_hub_routes_neg_p05_by_proxy": route_survive_by_proxy,
        "sign_structure_note": "prolif_sign_structure.tsv shows Spearman(prolif,pressure) and "
                               "Spearman(prolif,target-expr); positive-positive => genuine suppressor confound.",
    }, indent=2), encoding="utf-8")

    # ---- Aim 2 ---- #
    print("[robustness] Aim 2: per-miRNA basal load under evidence/binary/TargetScan ...")
    load_tbl = mirna_basal_load_by_weighting(edges, mirna, clinical, hs)
    load_tbl.to_csv(AIM2_DIR / "mirna_basal_load_by_weighting.tsv", sep="\t", index=False)

    hub = hub_rank_summary(load_tbl)
    hub.to_csv(AIM2_DIR / "hub_rank_summary.tsv", sep="\t", index=False)

    print("[robustness] Aim 2: per-target top regulators under TargetScan ...")
    per_target = per_target_top_regulators(edges, hs)
    per_target.to_csv(AIM2_DIR / "per_target_top_regulators.tsv", sep="\t", index=False)

    (AIM2_DIR / "method_manifest.json").write_text(json.dumps({
        "module": "mirna_hallmark.robustness_checks (Aim 2)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "Is the basal hub an artifact of miRTarBase study-count weighting?",
        "weightings": {
            "evidence": "Σ log1p(evidence_score) over high-evidence edges (current, curation-prone)",
            "binary": "count of high-evidence target edges (degree; equal weight per edge)",
            "targetscan": "Σ |weighted context++ score| over the same edges (sequence-based, study-count-independent)",
        },
        "targetscan_source": str(TARGETSCAN_CONTEXT),
        "hub": list(HUB_MIRNAS),
        "rank_correlations": hub.attrs.get("rank_correlations", {}),
    }, indent=2), encoding="utf-8")

    # ---- Aim 1b: per-family basal-specificity decomposition ---- #
    print("[robustness] Aim 1b: per-family load share + within-subtype coupling ...")
    fam_load = family_load_share(load_tbl)
    fam_load.to_csv(AIM2_DIR / "family_load_share.tsv", sep="\t", index=False)

    fam_detail = family_coupling_by_subtype(rna, clinical, proxies, hs, gate, cnv)
    fam_detail.to_csv(AIM1_DIR / "family_coupling_by_pam50.tsv", sep="\t", index=False)
    fam_summary = family_specificity_summary(fam_detail)
    fam_summary.to_csv(AIM1_DIR / "family_specificity_summary.tsv", sep="\t", index=False)
    if not fam_summary.empty:
        print("[robustness]   per-family (subtype: n_neg_sig / n_hallmarks, median rho):")
        for fam, g in fam_summary.groupby("family"):
            cells = ", ".join(
                f"{r.subtype} {int(r.n_neg_sig)}/{int(r.n_hallmarks)} ({r.median_rho:+.2f})"
                for r in g.itertuples()
            )
            print(f"[robustness]     {fam}: {cells}")
    (AIM1_DIR / "family_specificity_manifest.json").write_text(json.dumps({
        "module": "mirna_hallmark.robustness_checks (Aim 1b: per-family basal-specificity)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "For each top basal-load miRNA family, is its pressure->repression "
                    "coupling basal-SPECIFIC (concentrated in Basal) or merely the family "
                    "being basal-ENRICHED in expression (the load metric)?",
        "families": {k: v for k, v in BASAL_MIRNA_FAMILIES.items()},
        "load_share_total_evidence": fam_load.attrs.get("total_positive_load_evidence"),
        "load_share_total_targetscan": fam_load.attrs.get("total_positive_load_targetscan"),
        "pressure": "family-restricted compute_pressure (evidence-weighted Σ z), AGO-gated; "
                    "Hallmark = mean over member genes",
        "coupling": "partial Spearman(target expr, family pressure | CPE+HRD+E2F/G2M prolif) "
                    "within each PAM50 subtype; rho<0 = repression; FDR within (family, subtype) at q<0.10",
        "specificity_readout": "family is basal-specific if negative & FDR-significant coupling "
                               "concentrates in Basal vs LumA/LumB/Her2",
    }, indent=2), encoding="utf-8")

    print(f"[robustness] done. Outputs under {out_dir}")

    print("[robustness] TargetScan orphan coupling (gene + program, specificity) ...")
    from mirna_hallmark.eval.targetscan_orphan_coupling import run as run_orphan
    run_orphan()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=ROBUSTNESS_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
