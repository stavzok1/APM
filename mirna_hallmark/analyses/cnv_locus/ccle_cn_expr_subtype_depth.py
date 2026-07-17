"""Precursor-level TCGA↔CCLE CN→expression concordance + breast PAM50 / TAD line depth.

Fair comparison: aggregate mirBase arms that share one CCLE NanoString probe
(e.g. hsa-miR-30d-5p + hsa-miR-30d-3p → ``hsa-miR-30d``) before Spearman.

Breast depth uses ``pipeline.cell_line_subtype_map`` (TADs-aligned PAM50 labels)
on DepMap ``StrippedCellLineName`` for CCLE breast lines, including the 12 lines
that overlap Kim/TAD reference panel (T47D, HCC70, …).
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
from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import (
    MatchMode,
    PanelMode,
    RhoKind,
    _arm_cnv_from_loci,
    _ccle_covariates_for_stratum,
    _correlation_pair,
    _resolve_cn_profiles,
    _tcga_confounder_covariates,
    build_ccle_arm_expression_map,
    build_sample_locus_cnv,
    resolve_panel_arms,
)
from pipeline.biosample_names import normalize_cell_line_label
from pipeline.cell_line_subtype_map import (
    PAM50_GROUP_TO_CELL_LINES,
    PAM50_TUMOR_GROUPS,
    subtype_for_key,
)


def audit_tad_breast_line_coverage(
    panel_arms: Sequence[str],
    *,
    match_mode: MatchMode = "alias",
) -> pd.DataFrame:
    """Why each TAD breast line is included or excluded from CN×expr joins."""
    with open(C.CCLE_MIRNA_GCT) as f:
        f.readline()
        f.readline()
        gct_cols = set(f.readline().strip().split("\t")[2:])

    prof_map = _resolve_cn_profiles()
    default = pd.read_csv(C.CCLE_DEFAULT_PROFILES)
    default_dna = default.loc[default["ProfileType"].astype(str).str.lower() == "dna"]
    default_dna = default_dna.drop_duplicates("ModelID").set_index("ModelID")["ProfileID"]

    op = pd.read_csv(C.CCLE_OMICS_PROFILES)
    seg_p = set(pd.read_csv(C.CCLE_CN_SEGMENTS, usecols=["ProfileID"])["ProfileID"])

    arm_to_row, _ = build_ccle_arm_expression_map(panel_arms, match_mode=match_mode)
    from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import _ccle_expr_matrix

    expr_cols = set(_ccle_expr_matrix(arm_to_row).columns) if arm_to_row else set()

    meta = _ccle_model_subtype_table()
    tad = meta.loc[meta["tad_reference_line"] & (meta["OncotreeLineage"] == "Breast")].copy()
    rows: List[dict] = []
    for _, r in tad.iterrows():
        mid = r["ModelID"]
        def_pid = default_dna.get(mid)
        use_pid = prof_map.get(mid)
        wes = op.loc[
            (op["ModelID"] == mid) & (op["Datatype"].astype(str).str.lower() == "wes"),
            "ProfileID",
        ].tolist()
        reasons: List[str] = []
        if r["CCLEName"] not in gct_cols:
            reasons.append("no_CCLE_miRNA_GCT_column")
        if pd.isna(def_pid):
            reasons.append("no_default_DNA_profile")
        elif def_pid not in seg_p:
            reasons.append("default_WGS_not_in_PureCN_segments")
        if use_pid is None:
            reasons.append("no_WES/WGS_profile_with_segments")
        elif use_pid != def_pid:
            reasons.append(f"using_fallback_profile_{use_pid}")
        if mid not in expr_cols:
            reasons.append("expr_matrix_missing")
        rows.append({
            "line": r["line"],
            "ModelID": mid,
            "CCLEName": r["CCLEName"],
            "pam50_group": r["pam50_group"],
            "expr_in_gct": r["CCLEName"] in gct_cols,
            "default_dna_profile": def_pid,
            "cn_profile_used": use_pid,
            "wes_profiles_in_segments": ";".join(p for p in wes if p in seg_p),
            "in_expr_matrix": mid in expr_cols,
            "in_cn_cache": mid in prof_map.index,
            "join_ok": mid in expr_cols and mid in prof_map.index,
            "exclude_reason": ";".join(reasons) if reasons else "",
        })
    return pd.DataFrame(rows)


def _precursor_groups(map_notes: pd.DataFrame) -> Dict[str, List[str]]:
    """ccle_description (precursor probe) -> mirBase arms."""
    out: Dict[str, List[str]] = {}
    for _, r in map_notes.iterrows():
        prec = str(r["ccle_description"])
        out.setdefault(prec, []).append(str(r["arm"]))
    return out


def _aggregate_wide(
    wide: pd.DataFrame,
    groups: Dict[str, List[str]],
    *,
    how: str = "median",
) -> pd.DataFrame:
    """Collapse arm rows to precursor probe rows (median across arms per column)."""
    rows = []
    for prec, arms in groups.items():
        present = [a for a in arms if a in wide.index]
        if not present:
            continue
        sub = wide.loc[present]
        if how == "mean":
            rows.append((prec, sub.mean(axis=0)))
        else:
            rows.append((prec, sub.median(axis=0)))
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame({k: v for k, v in rows}).T


def _correlation_table(
    cn_wide: pd.DataFrame,
    expr_wide: pd.DataFrame,
    *,
    strata: Dict[str, Set[str]],
    min_n: int,
    entity_col: str = "precursor",
    rho_kind: RhoKind = "partial",
    tcga_strata: Optional[Set[str]] = None,
    min_pair_n: Optional[int] = None,
) -> pd.DataFrame:
    """Per-entity marginal or partial Spearman(CN, expr) by stratum."""
    pair_n = min_n if min_pair_n is None else min_pair_n
    rows: List[dict] = []
    entities = sorted(set(cn_wide.index) & set(expr_wide.index))
    for ent in entities:
        cn = cn_wide.loc[ent].astype(float)
        ex = expr_wide.loc[ent].astype(float)
        for stratum, ids in strata.items():
            mids = sorted(set(cn.index) & set(ex.index) & ids)
            if len(mids) < min_n:
                continue
            c = cn.reindex(mids)
            e = ex.reindex(mids)
            ok = c.notna() & e.notna()
            if ok.sum() < min_n:
                continue
            mids_ok = [m for m in mids if ok.get(m, False)]
            if stratum.startswith("tcga_"):
                cov = _tcga_confounder_covariates(mids_ok) if rho_kind == "partial" else None
            else:
                cov = (
                    _ccle_covariates_for_stratum(stratum, mids_ok)
                    if rho_kind == "partial"
                    else None
                )
            rho, p, n = _correlation_pair(
                c.loc[mids_ok], e.loc[mids_ok], cov, rho_kind=rho_kind, min_pair_n=pair_n,
            )
            row = {
                entity_col: ent,
                "stratum": stratum,
                "n": int(n),
                "rho_kind": rho_kind,
            }
            if rho_kind == "partial":
                row["partial_rho"] = round(float(rho), 4) if pd.notna(rho) else np.nan
                row["partial_p"] = float(p) if pd.notna(p) else np.nan
            else:
                row["spearman_rho"] = round(float(rho), 4) if pd.notna(rho) else np.nan
                row["spearman_p"] = float(p) if pd.notna(p) else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    pcol = "partial_p" if rho_kind == "partial" else "spearman_p"
    qcol = "partial_q" if rho_kind == "partial" else "spearman_q"
    sort_col = "partial_rho" if rho_kind == "partial" else "spearman_rho"
    for stratum in df["stratum"].unique():
        m = df["stratum"] == stratum
        valid = df.loc[m, pcol].notna()
        if valid.any():
            df.loc[m & valid, qcol] = S.bh_fdr(df.loc[m & valid, pcol].values)
    return df.sort_values(["stratum", sort_col], ascending=[True, False])


def _tcga_arm_cn_expr(panel_arms: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """TCGA arm-level CN and expression (participant columns)."""
    from mirna_hallmark import data_loaders as D

    long_path = C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz"
    if not long_path.exists():
        raise FileNotFoundError(f"Run mirna_locus_cnv first: {long_path}")
    long = pd.read_csv(long_path, sep="\t", low_memory=False)
    arm_long = _arm_cnv_from_loci(long, panel_arms)
    cn = (
        arm_long.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack("participant")
    )
    expr = D.load_mirna_arms()
    expr = expr.loc[[a for a in panel_arms if a in expr.index]]
    return cn, expr


def _ccle_model_subtype_table() -> pd.DataFrame:
    """ModelID × PAM50 group from TAD-aligned cell_line_subtype_map."""
    model = pd.read_csv(C.CCLE_MODEL, usecols=["ModelID", "CCLEName", "StrippedCellLineName", "OncotreeLineage"])
    rows = []
    for _, r in model.iterrows():
        token = str(r["StrippedCellLineName"] or r["CCLEName"].rsplit("_", 1)[0])
        line = normalize_cell_line_label(token)
        rec = subtype_for_key(line)
        rows.append({
            "ModelID": r["ModelID"],
            "CCLEName": r["CCLEName"],
            "line": line,
            "OncotreeLineage": r["OncotreeLineage"],
            "pam50_group": rec.pam50_group,
            "er_status": rec.er_status,
            "her2_status": rec.her2_status,
            "tnbc": rec.tnbc,
            "subtype_source": rec.source,
        })
    df = pd.DataFrame(rows)
    tad_lines = {normalize_cell_line_label(x) for xs in PAM50_GROUP_TO_CELL_LINES.values() for x in xs}
    df["tad_reference_line"] = df["line"].isin(tad_lines)
    return df


def _breast_strata(min_n: int) -> Dict[str, Set[str]]:
    """PAM50 strata for breast CCLE lines with known subtype."""
    meta = _ccle_model_subtype_table()
    br = meta.loc[meta["OncotreeLineage"] == "Breast"]
    strata: Dict[str, Set[str]] = {"breast_all": set(br["ModelID"])}
    for pam in list(PAM50_TUMOR_GROUPS) + ["Normal_like"]:
        mids = set(br.loc[br["pam50_group"] == pam, "ModelID"])
        if len(mids) >= min_n:
            strata[f"breast_{pam}"] = mids
    tad = set(br.loc[br["tad_reference_line"], "ModelID"])
    if len(tad) >= min_n:
        strata["breast_tad_panel"] = tad
    return strata


def _line_level_correlation(
    cn_wide: pd.DataFrame,
    expr_wide: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    entities: Sequence[str],
    min_lines: int = 4,
) -> pd.DataFrame:
    """Across cell lines (one median CN/expr per line): Spearman within PAM50 group."""
    from scipy.stats import spearmanr

    br = meta.loc[(meta["OncotreeLineage"] == "Breast") & meta["line"].notna()]
    rows: List[dict] = []
    for ent in entities:
        if ent not in cn_wide.index or ent not in expr_wide.index:
            continue
        for pam in list(PAM50_TUMOR_GROUPS) + ["Normal_like"]:
            sub = br.loc[br["pam50_group"] == pam]
            pts = []
            for _, r in sub.iterrows():
                mid = r["ModelID"]
                if mid not in cn_wide.columns or mid not in expr_wide.columns:
                    continue
                cn = cn_wide.loc[ent, mid]
                ex = expr_wide.loc[ent, mid]
                if pd.notna(cn) and pd.notna(ex):
                    pts.append((r["line"], float(cn), float(ex)))
            if len(pts) < min_lines:
                continue
            cn_v = [p[1] for p in pts]
            ex_v = [p[2] for p in pts]
            rho, p = spearmanr(cn_v, ex_v)
            rows.append({
                "entity": ent,
                "pam50_group": pam,
                "n_lines": len(pts),
                "lines": ";".join(p[0] for p in pts),
                "spearman_rho": round(float(rho), 4),
                "spearman_p": float(p),
            })
    return pd.DataFrame(rows)


def run_precursor_tcga_ccle(
    *,
    panel: PanelMode = "focus_8q",
    match_mode: MatchMode = "alias",
    out_dir: Optional[Path] = None,
    min_n: int = 15,
    rho_kind: RhoKind = "partial",
) -> pd.DataFrame:
    if out_dir is None:
        base = C.CCLE_CN_CONCORDANCE_DIR / "precursor_tcga_ccle"
        out_dir = base / panel if panel != "focus_8q" else base
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_arms = resolve_panel_arms(panel, rho_kind=rho_kind)
    _, map_notes = build_ccle_arm_expression_map(panel_arms, match_mode=match_mode)
    map_notes.to_csv(out_dir / "arm_to_ccle_probe_map.tsv", sep="\t", index=False)
    groups = _precursor_groups(map_notes)
    mapped_arms = sorted({a for arms in groups.values() for a in arms})

    # TCGA
    tcga_cn, tcga_expr = _tcga_arm_cn_expr(mapped_arms)
    tcga_cn_p = _aggregate_wide(tcga_cn, groups)
    tcga_expr_p = _aggregate_wide(tcga_expr, groups)
    tcga_strata = {"all": set(tcga_cn_p.columns)}
    clin = D.load_clinical_strata()
    for pam in PAM50_TUMOR_GROUPS:
        mids = set(clin.loc[clin["PAM50_final"] == pam, "participant"])
        mids &= set(tcga_cn_p.columns)
        if len(mids) >= min_n:
            tcga_strata[f"tcga_{pam}"] = mids

    tcga_conc = _correlation_table(
        tcga_cn_p, tcga_expr_p, strata=tcga_strata, min_n=min_n, rho_kind=rho_kind,
    )
    tcga_conc.to_csv(out_dir / "tcga_precursor_concordance.tsv", sep="\t", index=False)

    # CCLE (reuse cached locus CN)
    from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import _ccle_expr_matrix

    arm_to_row, _ = build_ccle_arm_expression_map(mapped_arms, match_mode=match_mode)
    expr = _ccle_expr_matrix(arm_to_row)
    locus_long = build_sample_locus_cnv(
        cache_path=C.CCLE_CN_CONCORDANCE_DIR / "tables" / "sample_locus_cn.tsv.gz",
        force=False,
        refresh_missing=True,
    )
    locus_long = locus_long.loc[locus_long["participant"].isin(expr.columns)]
    arm_cn = _arm_cnv_from_loci(locus_long, mapped_arms)
    ccle_cn = (
        arm_cn.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack("participant")
    )
    ccle_cn_p = _aggregate_wide(ccle_cn, groups)
    ccle_expr_p = _aggregate_wide(expr, groups)

    ccle_strata = {"all": set(ccle_cn_p.columns)}
    ccle_strata.update(_breast_strata(min_n=min(8, min_n)))
    ccle_conc = _correlation_table(
        ccle_cn_p, ccle_expr_p, strata=ccle_strata, min_n=min(8, min_n), rho_kind=rho_kind,
    )
    ccle_conc.to_csv(out_dir / "ccle_precursor_concordance.tsv", sep="\t", index=False)

    rho_col = "partial_rho" if rho_kind == "partial" else "spearman_rho"
    q_col = "partial_q" if rho_kind == "partial" else "spearman_q"

    # Side-by-side compare (precursor entities)
    if not tcga_conc.empty and not ccle_conc.empty:
        t = tcga_conc.loc[tcga_conc["stratum"] == "all"].rename(columns={
            "n": "tcga_n", rho_col: "tcga_rho", q_col: "tcga_q",
        })
        c = ccle_conc.loc[ccle_conc["stratum"] == "all"].rename(columns={
            "n": "ccle_n", rho_col: "ccle_rho", q_col: "ccle_q",
        })
        cmp = t.merge(c, on="precursor", how="inner")
        cmp["rho_delta"] = (cmp["ccle_rho"] - cmp["tcga_rho"]).round(4)
        br = ccle_conc.loc[ccle_conc["stratum"] == "breast_tad_panel"].rename(columns={
            rho_col: "ccle_breast_tad_rho",
            "partial_p" if rho_kind == "partial" else "spearman_p": "ccle_breast_tad_p",
        })
        if not br.empty:
            cmp = cmp.merge(
                br[["precursor", "ccle_breast_tad_rho", "n"]].rename(columns={"n": "ccle_breast_tad_n"}),
                on="precursor",
                how="left",
            )
        cmp.to_csv(out_dir / "precursor_tcga_ccle_compare.tsv", sep="\t", index=False)

    sig_def = (
        "partial_q<0.05 & partial_rho>0 (CPE+HRD)"
        if panel == "tcga_sig" and rho_kind == "partial"
        else "partial_q<0.05 & partial_rho>0"
        if rho_kind == "partial"
        else "spearman_q<0.05 & spearman_rho>0"
    )
    manifest = {
        "module": "mirna_hallmark.ccle_cn_expr_subtype_depth",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "panel": panel,
        "match_mode": match_mode,
        "rho_kind": rho_kind,
        "panel_definition": (
            f"focus_8q: TCGA focus amplicon arms ({sig_def} not used for panel pick)"
            if panel == "focus_8q"
            else f"tcga_sig: TCGA CN→expr {sig_def} (NOT CCLE intersection)"
            if panel == "tcga_sig"
            else "high_evidence: Hallmark high-evidence arms"
        ),
        "tcga_partial_covariates": list(C.CONFOUNDER_NUMERIC),
        "ccle_partial_covariates": "lineage dummies (all); PAM50 dummies (breast)",
        "n_panel_arms": len(panel_arms),
        "n_ccle_expr_mapped_arms": len(arm_to_row),
        "n_precursor_probes": len(groups),
        "n_tcga_precursor_tests": int(len(tcga_conc)),
        "n_ccle_precursor_tests": int(len(ccle_conc)),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[precursor] wrote under {out_dir}")
    return tcga_conc


def run_breast_tad_line_depth(
    *,
    panel: PanelMode = "focus_8q",
    match_mode: MatchMode = "alias",
    out_dir: Optional[Path] = None,
    rho_kind: RhoKind = "partial",
) -> pd.DataFrame:
    """Per-line CN/expr for TAD reference breast lines + line-level PAM50 correlations."""
    if out_dir is None:
        base = C.CCLE_CN_CONCORDANCE_DIR / "breast_tad_lines"
        out_dir = base / panel if panel != "focus_8q" else base
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_arms = resolve_panel_arms(panel, rho_kind=rho_kind)
    arm_to_row, map_notes = build_ccle_arm_expression_map(panel_arms, match_mode=match_mode)
    groups = _precursor_groups(map_notes)
    mapped_arms = sorted(arm_to_row.keys())

    from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import _ccle_expr_matrix

    expr = _ccle_expr_matrix(arm_to_row)
    locus_long = build_sample_locus_cnv(
        cache_path=C.CCLE_CN_CONCORDANCE_DIR / "tables" / "sample_locus_cn.tsv.gz",
        force=False,
        refresh_missing=True,
    )
    arm_cn = _arm_cnv_from_loci(locus_long, mapped_arms)
    cn = (
        arm_cn.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack("participant")
    )
    cn_p = _aggregate_wide(cn, groups)
    expr_p = _aggregate_wide(expr, groups)

    meta = _ccle_model_subtype_table()
    meta.to_csv(out_dir / "ccle_breast_subtype_map.tsv", sep="\t", index=False)
    audit_tad_breast_line_coverage(panel_arms, match_mode=match_mode).to_csv(
        out_dir / "tad_line_coverage_audit.tsv", sep="\t", index=False
    )

    tad = meta.loc[meta["tad_reference_line"] & (meta["OncotreeLineage"] == "Breast")]
    long_rows: List[dict] = []
    for _, line in tad.iterrows():
        mid = line["ModelID"]
        for ent in cn_p.index:
            if mid in cn_p.columns and mid in expr_p.columns:
                long_rows.append({
                    "line": line["line"],
                    "CCLEName": line["CCLEName"],
                    "ModelID": mid,
                    "pam50_group": line["pam50_group"],
                    "precursor": ent,
                    "copy_number": cn_p.loc[ent, mid],
                    "log2_expr": expr_p.loc[ent, mid],
                })
    per_line = pd.DataFrame(long_rows)
    per_line.to_csv(out_dir / "tad_line_precursor_cn_expr.tsv", sep="\t", index=False)

    # PAM50-group concordance at sample level (all breast lines with known subtype)
    br_strata = _breast_strata(min_n=4)
    br_conc = _correlation_table(
        cn_p, expr_p, strata=br_strata, min_n=4, entity_col="precursor", rho_kind=rho_kind,
    )
    br_conc.to_csv(out_dir / "breast_pam50_precursor_concordance.tsv", sep="\t", index=False)

    # Line-level: correlate across TAD lines within each PAM50 bucket
    line_corr = _line_level_correlation(cn_p, expr_p, meta, entities=cn_p.index, min_lines=3)
    line_corr.to_csv(out_dir / "breast_line_level_pam50_correlation.tsv", sep="\t", index=False)

    # Arm-level breast PAM50 for focus entities
    arm_br = _correlation_table(
        cn, expr, strata=br_strata, min_n=4, entity_col="arm", rho_kind=rho_kind,
    )
    arm_br.to_csv(out_dir / "breast_pam50_arm_concordance.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.ccle_cn_expr_subtype_depth",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "panel": panel,
        "rho_kind": rho_kind,
        "n_tad_breast_lines": int(len(tad)),
        "subtype_map": "pipeline.cell_line_subtype_map (TADs-aligned)",
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[breast_tad] wrote under {out_dir} ({len(tad)} TAD reference lines)")
    return per_line


def _tad_breast_table(meta: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    meta = meta if meta is not None else _ccle_model_subtype_table()
    return meta.loc[
        meta["tad_reference_line"] & (meta["OncotreeLineage"] == "Breast")
    ].copy()


def _tad_panel_strata(tad: pd.DataFrame, *, min_n: int = 3) -> Dict[str, Set[str]]:
    """ModelID sets for TAD breast lines only."""
    strata: Dict[str, Set[str]] = {"tad_panel_11": set(tad["ModelID"])}
    for pam in list(PAM50_TUMOR_GROUPS) + ["Normal_like"]:
        mids = set(tad.loc[tad["pam50_group"] == pam, "ModelID"])
        if len(mids) >= min_n:
            strata[f"tad_{pam}"] = mids
        elif len(mids) == 1:
            strata[f"tad_{pam}_n1"] = mids
    return strata


def _within_line_arm_correlation(
    cn_wide: pd.DataFrame,
    expr_wide: pd.DataFrame,
    tad: pd.DataFrame,
    *,
    min_arms: int = 30,
) -> pd.DataFrame:
    """Per TAD line: Spearman(CN, expr) across high-evidence arms (one point per arm)."""
    from scipy.stats import spearmanr

    rows: List[dict] = []
    for _, r in tad.iterrows():
        mid = r["ModelID"]
        if mid not in cn_wide.columns or mid not in expr_wide.columns:
            rows.append({
                "line": r["line"],
                "ModelID": mid,
                "pam50_group": r["pam50_group"],
                "n_arms": 0,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "note": "missing_cn_or_expr_matrix",
            })
            continue
        cn = cn_wide[mid].astype(float)
        ex = expr_wide[mid].astype(float)
        ok = cn.notna() & ex.notna()
        if ok.sum() < min_arms:
            rows.append({
                "line": r["line"],
                "ModelID": mid,
                "pam50_group": r["pam50_group"],
                "n_arms": int(ok.sum()),
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "note": f"fewer_than_{min_arms}_arms",
            })
            continue
        rho, p = spearmanr(cn[ok], ex[ok])
        rows.append({
            "line": r["line"],
            "ModelID": mid,
            "pam50_group": r["pam50_group"],
            "n_arms": int(ok.sum()),
            "spearman_rho": round(float(rho), 4),
            "spearman_p": float(p),
            "note": "",
        })
    df = pd.DataFrame(rows)
    if not df.empty and df["spearman_p"].notna().any():
        df.loc[df["spearman_p"].notna(), "spearman_q"] = S.bh_fdr(
            df.loc[df["spearman_p"].notna(), "spearman_p"].values
        )
    return df.sort_values("spearman_rho", ascending=False, na_position="last")


def run_tad_panel_high_evidence_concordance(
    *,
    panel: PanelMode = "high_evidence",
    match_mode: MatchMode = "alias",
    out_dir: Optional[Path] = None,
    rho_kind: RhoKind = "marginal",
    min_arms_per_line: int = 30,
    min_lines_panel: int = 3,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """TAD 11-line panel: panel-wide per-arm CN×expr + per-line across-arm CN×expr.

    Panel-wide: for each mapped arm, Spearman(CN, expr) across the 11 TAD lines
    (and within LumA / Basal subsets when n ≥ 3). Uses **marginal** Spearman by
    default — partial with PAM50 dummies is not identified at n = 11.

    Per-line: for each TAD line, Spearman(CN, expr) across all mapped panel arms
    (one point per miRNA arm within that line).
    """
    from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import (
        _ccle_expr_matrix,
        _ensure_partial_q,
    )

    if out_dir is None:
        base = C.CCLE_CN_CONCORDANCE_DIR / "breast_tad_lines"
        out_dir = base / panel if panel != "focus_8q" else base / "high_evidence_tad"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_arms = resolve_panel_arms(panel, rho_kind="partial" if panel == "tcga_sig" else rho_kind)
    arm_to_row, map_notes = build_ccle_arm_expression_map(panel_arms, match_mode=match_mode)
    map_notes.to_csv(out_dir / "arm_to_ccle_probe_map.tsv", sep="\t", index=False)
    mapped_arms = sorted(arm_to_row.keys())
    print(
        f"[tad_he] panel={panel} arms={len(panel_arms)} "
        f"ccle_mapped={len(mapped_arms)} match={match_mode}"
    )
    if not mapped_arms:
        raise ValueError("No arms mapped to CCLE expression")

    expr = _ccle_expr_matrix(arm_to_row)
    locus_long = build_sample_locus_cnv(
        cache_path=C.CCLE_CN_CONCORDANCE_DIR / "tables" / "sample_locus_cn.tsv.gz",
        force=False,
        refresh_missing=True,
    )
    arm_cn = _arm_cnv_from_loci(locus_long, mapped_arms)
    cn = (
        arm_cn.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack("participant")
    )

    meta = _ccle_model_subtype_table()
    tad = _tad_breast_table(meta)
    join_ok = tad.loc[tad["ModelID"].isin(cn.columns) & tad["ModelID"].isin(expr.columns)]
    print(f"[tad_he] TAD breast lines with CN+expr: {len(join_ok)} / {len(tad)}")

    tad_mids = sorted(join_ok["ModelID"])
    cn_t = cn.loc[[a for a in mapped_arms if a in cn.index], tad_mids]
    ex_t = expr.loc[[a for a in mapped_arms if a in expr.index], tad_mids]

    # Long per-line × arm table
    long_rows: List[dict] = []
    for _, r in join_ok.iterrows():
        mid = r["ModelID"]
        for arm in cn_t.index:
            if arm not in ex_t.index:
                continue
            cn_v, ex_v = cn_t.loc[arm, mid], ex_t.loc[arm, mid]
            if pd.notna(cn_v) and pd.notna(ex_v):
                long_rows.append({
                    "line": r["line"],
                    "ModelID": mid,
                    "pam50_group": r["pam50_group"],
                    "arm": arm,
                    "copy_number": float(cn_v),
                    "log2_expr": float(ex_v),
                })
    per_line_arm = pd.DataFrame(long_rows)
    per_line_arm.to_csv(out_dir / "per_line_arm_cn_expr.tsv.gz", sep="\t", index=False, compression="gzip")

    # Panel-wide: each arm across TAD lines
    tad_strata = _tad_panel_strata(join_ok, min_n=3)
    panel_conc = _correlation_table(
        cn_t,
        ex_t,
        strata=tad_strata,
        min_n=min_lines_panel,
        entity_col="arm",
        rho_kind="marginal",
    )
    panel_conc.to_csv(out_dir / "tad_panel_arm_concordance.tsv", sep="\t", index=False)

    # Per-line: across all arms within line
    per_line = _within_line_arm_correlation(
        cn_t, ex_t, join_ok, min_arms=min_arms_per_line,
    )
    per_line.to_csv(out_dir / "per_line_across_arms_concordance.tsv", sep="\t", index=False)

    # TCGA partial compare for panel-wide tad_panel_11 stratum
    tad_only = panel_conc.loc[panel_conc["stratum"] == "tad_panel_11"].copy()
    rename_map = {
        "spearman_rho": "ccle_tad_rho",
        "n": "ccle_tad_n",
    }
    if "spearman_q" in tad_only.columns:
        rename_map["spearman_q"] = "ccle_tad_q"
    tad_only = tad_only.rename(columns=rename_map)
    tcga = _ensure_partial_q(
        pd.read_csv(C.MIRNA_LOCUS_CNV_DIR / "mirna_cnv_expr_concordance.tsv", sep="\t")
    )
    cmp = tcga[["arm", "partial_rho", "partial_q", "n"]].rename(columns={
        "partial_rho": "tcga_partial_rho",
        "partial_q": "tcga_partial_q",
        "n": "tcga_n",
    })
    cmp_cols = ["arm", "ccle_tad_n", "ccle_tad_rho"]
    if "ccle_tad_q" in tad_only.columns:
        cmp_cols.append("ccle_tad_q")
    cmp = cmp.merge(tad_only[cmp_cols], on="arm", how="inner")
    cmp["rho_delta_tad_minus_tcga"] = (cmp["ccle_tad_rho"] - cmp["tcga_partial_rho"]).round(4)
    cmp = cmp.sort_values("tcga_partial_rho", ascending=False)
    cmp.to_csv(out_dir / "tcga_partial_vs_tad_panel_compare.tsv", sep="\t", index=False)

    # Summary counts
    n11 = panel_conc.loc[panel_conc["stratum"] == "tad_panel_11"]
    sig11 = n11.loc[(n11["spearman_q"] < C.FDR_ALPHA) & (n11["spearman_rho"] > 0)]
    manifest = {
        "module": "mirna_hallmark.ccle_cn_expr_subtype_depth.run_tad_panel_high_evidence_concordance",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "panel": panel,
        "match_mode": match_mode,
        "n_panel_arms": len(panel_arms),
        "n_ccle_mapped_arms": len(mapped_arms),
        "n_tad_lines_joined": int(len(join_ok)),
        "panel_wide_test": "Spearman(CN, expr) per arm across TAD lines (marginal)",
        "per_line_test": "Spearman(CN, expr) per line across mapped arms",
        "min_arms_per_line": min_arms_per_line,
        "min_lines_panel_stratum": min_lines_panel,
        "tad_panel_11_pos_fdr": int(len(sig11)),
        "tad_panel_11_tested": int(len(n11)),
        "per_line_pos_fdr": int(
            len(per_line.loc[(per_line["spearman_q"] < C.FDR_ALPHA) & (per_line["spearman_rho"] > 0)])
        ) if "spearman_q" in per_line.columns else 0,
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(
        f"[tad_he] panel-wide tad_11: {len(sig11)}/{len(n11)} FDR+ pos; "
        f"per-line rows: {len(per_line)} → {out_dir}"
    )
    return panel_conc, per_line


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", choices=("focus_8q", "tcga_sig", "high_evidence"), default="focus_8q")
    ap.add_argument("--match-mode", choices=("alias", "exact", "precursor_stem"), default="alias")
    ap.add_argument(
        "--rho-kind",
        choices=("partial", "marginal"),
        default="partial",
        help="Partial (default): TCGA adjusts CPE+HRD; CCLE adjusts lineage or PAM50",
    )
    ap.add_argument("--skip-precursor", action="store_true")
    ap.add_argument("--skip-breast", action="store_true")
    ap.add_argument(
        "--tad-panel-he",
        action="store_true",
        help="TAD 11-line panel-wide + per-line concordance for --panel (default high_evidence)",
    )
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.tad_panel_he:
        run_tad_panel_high_evidence_concordance(
            panel=args.panel,
            match_mode=args.match_mode,
            rho_kind="marginal",
        )
    if not args.skip_precursor:
        run_precursor_tcga_ccle(
            panel=args.panel, match_mode=args.match_mode, rho_kind=args.rho_kind,
        )
    if not args.skip_breast:
        run_breast_tad_line_depth(
            panel=args.panel, match_mode=args.match_mode, rho_kind=args.rho_kind,
        )


if __name__ == "__main__":
    main()
