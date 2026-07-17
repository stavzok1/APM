"""CCLE cell-line replication: miRNA locus CN vs mature-arm expression concordance.

Uses DepMap PureCN absolute segments at pre-miRNA loci (hg38 overlap, mirroring
``mirna_locus_cnv``) and CCLE NanoString miRNA expression (log2(50 + nSolver)).

Panel modes (``--panel``):
- ``high_evidence`` — Hallmark high-evidence arms with CCLE expression
- ``focus_8q`` — TCGA focus amplicon arms (miR-151a/30b/30d/… from CNV maps)
- ``tcga_sig`` — TCGA partial CN→expr concordance (partial_q<0.05, partial_rho>0; CPE+HRD adjusted)

Match modes (``--match-mode``):
- ``alias`` (default) — mirBase→CCLE alias map (151a→151-3p, precursor stems, …)
- ``exact`` — only exact Description matches (stricter; misses many TCGA stars)
- ``precursor_stem`` — legacy stem-only mapping

Outputs under ``output/ccle_mirna_cn_concordance/`` (+ ``--panel`` subfolder when not
``high_evidence``).
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S

PanelMode = Literal["high_evidence", "focus_8q", "tcga_sig"]
MatchMode = Literal["exact", "alias", "precursor_stem"]
RhoKind = Literal["marginal", "partial"]


def _ensure_partial_q(tcga: pd.DataFrame) -> pd.DataFrame:
    if "partial_q" in tcga.columns and tcga["partial_q"].notna().any():
        return tcga
    out = tcga.copy()
    out["partial_q"] = np.nan
    m = out["partial_p"].notna()
    if m.any():
        out.loc[m, "partial_q"] = S.bh_fdr(out.loc[m, "partial_p"].values)
    return out


def _tcga_concordance_sig_arms(rho_kind: RhoKind = "partial") -> List[str]:
    """Arms significant in TCGA CN→expr concordance (marginal or partial Spearman)."""
    tcga_path = C.MIRNA_LOCUS_CNV_DIR / "mirna_cnv_expr_concordance.tsv"
    tcga = pd.read_csv(tcga_path, sep="\t")
    if rho_kind == "partial":
        tcga = _ensure_partial_q(tcga)
        sig = tcga.loc[
            tcga["partial_p"].notna()
            & (tcga["partial_q"] < C.FDR_ALPHA)
            & (tcga["partial_rho"] > 0)
        ]
    else:
        sig = tcga.loc[(tcga["spearman_q"] < C.FDR_ALPHA) & (tcga["spearman_rho"] > 0)]
    return sorted(sig["arm"].dropna().unique())


def _ccle_lineage_covariates(model_ids: Sequence[str]) -> pd.DataFrame:
    """OncotreeLineage one-hot for CCLE partial Spearman."""
    model = _load_model_lineage().set_index("ModelID")
    lineages = model.reindex(model_ids)["OncotreeLineage"].fillna("Unknown").astype(str)
    dummies = pd.get_dummies(lineages, prefix="lineage", dtype=float)
    return dummies.loc[list(model_ids)]


def _ccle_pam50_covariates(model_ids: Sequence[str]) -> pd.DataFrame:
    """PAM50 group one-hot for breast CCLE lines."""
    from pipeline.biosample_names import normalize_cell_line_label
    from pipeline.cell_line_subtype_map import subtype_for_key

    model = pd.read_csv(
        C.CCLE_MODEL, usecols=["ModelID", "StrippedCellLineName", "CCLEName"]
    ).set_index("ModelID")
    groups = []
    for mid in model_ids:
        if mid not in model.index:
            groups.append("Unknown")
            continue
        r = model.loc[mid]
        token = str(r["StrippedCellLineName"] or r["CCLEName"].rsplit("_", 1)[0])
        groups.append(subtype_for_key(normalize_cell_line_label(token)).pam50_group)
    dummies = pd.get_dummies(pd.Series(groups, index=list(model_ids)), prefix="pam50", dtype=float)
    return dummies


def _ccle_covariates_for_stratum(stratum: str, model_ids: Sequence[str]) -> pd.DataFrame:
    if stratum.startswith("breast"):
        return _ccle_pam50_covariates(model_ids)
    return _ccle_lineage_covariates(model_ids)


def _tcga_confounder_covariates(participants: Sequence[str]) -> pd.DataFrame:
    clin = D.load_clinical_strata().set_index("participant")
    cols = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    return clin.reindex(participants)[cols].astype(float)


def _correlation_pair(
    cn: pd.Series,
    expr: pd.Series,
    cov: Optional[pd.DataFrame],
    *,
    rho_kind: RhoKind,
    min_pair_n: int = 15,
) -> Tuple[float, float, int]:
    """Marginal or partial Spearman(CN, expr); returns (rho, p, n)."""
    from scipy.stats import spearmanr

    merged = pd.concat([cn.rename("cn"), expr.rename("expr")], axis=1)
    if cov is not None and not cov.empty and rho_kind == "partial":
        merged = merged.join(cov, how="inner")
    merged = merged.dropna()
    if len(merged) < min_pair_n:
        return np.nan, np.nan, len(merged)
    if rho_kind == "partial" and cov is not None and not cov.empty:
        from analysis.utils.common.loaders import partial_spearman

        cov_cols = merged.drop(columns=["cn", "expr"])
        if cov_cols.shape[1] == 0:
            rho, p = spearmanr(merged["cn"], merged["expr"])
            return float(rho), float(p), len(merged)
        rho, p, n = partial_spearman(merged["expr"], merged["cn"], cov_cols)
        return float(rho), float(p), n
    rho, p = spearmanr(merged["cn"], merged["expr"])
    return float(rho), float(p), len(merged)


def _arm_stem(arm: str) -> str:
    return re.sub(r"-[35]p(\*)?$", "", str(arm).strip())


def _ccle_alias_candidates(arm: str) -> List[str]:
    """Ordered CCLE Description candidates for a mirBase mature arm."""
    arm = str(arm).strip()
    out: List[str] = [arm]
    if arm in C.CCLE_COMBINED_MIRNA_PROBES:
        out.append(C.CCLE_COMBINED_MIRNA_PROBES[arm])
    # mirBase letter before -3p/-5p: hsa-miR-151a-3p → hsa-miR-151-3p
    m = re.match(r"^(hsa-miR-\d+)[a-z]-([35]p(\*)?)$", arm)
    if m:
        out.append(f"{m.group(1)}-{m.group(2)}")
    stem = _arm_stem(arm)
    if stem not in out:
        out.append(stem)
    return out


def resolve_panel_arms(panel: PanelMode, *, rho_kind: RhoKind = "partial") -> List[str]:
    """Arms to test for a given panel mode."""
    if panel == "focus_8q":
        p = C.CCLE_FOCUS_8Q_ARMS_PATH
        if not p.exists():
            raise FileNotFoundError(f"Focus arm table missing: {p}")
        return sorted(pd.read_csv(p, sep="\t")["arm_label"].dropna().unique())

    if panel == "tcga_sig":
        return _tcga_concordance_sig_arms(rho_kind)

    he = D.high_evidence_edges()["miRNA"].drop_duplicates().tolist()
    return sorted(he)


def build_ccle_arm_expression_map(
    arms: Sequence[str],
    *,
    match_mode: MatchMode = "alias",
) -> Tuple[Dict[str, str], pd.DataFrame]:
    """Map mature arm -> CCLE GCT row id (nmiR*.1 index)."""
    if not C.CCLE_MIRNA_GCT.exists():
        raise FileNotFoundError(f"CCLE miRNA GCT missing: {C.CCLE_MIRNA_GCT}")
    gct = pd.read_csv(C.CCLE_MIRNA_GCT, sep="\t", skiprows=2, index_col=0)
    desc = gct["Description"].astype(str)
    desc_to_row: Dict[str, str] = {}
    for row_id, d in desc.items():
        if d and d != "nan" and d not in desc_to_row:
            desc_to_row[d] = str(row_id)

    arm_to_row: Dict[str, str] = {}
    notes: List[dict] = []
    for arm in arms:
        if match_mode == "exact":
            candidates = [arm]
        elif match_mode == "precursor_stem":
            candidates = [arm, _arm_stem(arm)]
        else:
            candidates = _ccle_alias_candidates(arm)
        hit_desc = next((c for c in candidates if c in desc_to_row), None)
        if hit_desc is None:
            continue
        arm_to_row[arm] = desc_to_row[hit_desc]
        notes.append({
            "arm": arm,
            "ccle_description": hit_desc,
            "ccle_row": desc_to_row[hit_desc],
            "match": (
                "exact" if hit_desc == arm
                else "combined_probe" if hit_desc in C.CCLE_COMBINED_MIRNA_PROBES.values()
                else "alias" if match_mode == "alias"
                else "precursor_stem"
            ),
        })
    return arm_to_row, pd.DataFrame(notes)


def _load_model_lineage() -> pd.DataFrame:
    if not C.CCLE_MODEL.exists():
        raise FileNotFoundError(f"DepMap Model.csv missing: {C.CCLE_MODEL}")
    model = pd.read_csv(C.CCLE_MODEL, usecols=["ModelID", "CCLEName", "OncotreeLineage"])
    model = model.drop_duplicates("ModelID")
    model["lineage_slug"] = (
        model["OncotreeLineage"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True)
    )
    return model


def _stratum_model_sets(min_n: int) -> Dict[str, Set[str]]:
    """ModelID sets per stratum label (all + configured lineages with enough lines)."""
    model = _load_model_lineage()
    strata: Dict[str, Set[str]] = {"all": set(model["ModelID"])}
    for lineage in C.CCLE_LINEAGE_STRATA:
        mids = set(model.loc[model["OncotreeLineage"] == lineage, "ModelID"])
        if len(mids) >= min_n:
            slug = lineage.lower().replace("/", "_").replace(" ", "_")
            strata[slug] = mids
    other = set(model["ModelID"])
    for k, v in list(strata.items()):
        if k != "all":
            other -= v
    if len(other) >= min_n:
        strata["other_lineages"] = other
    return strata


def _ccle_expr_matrix(arm_to_row: Dict[str, str]) -> pd.DataFrame:
    """arm x ModelID log2(50 + nSolver) for CCLE lines with DepMap ModelID."""
    gct = pd.read_csv(C.CCLE_MIRNA_GCT, sep="\t", skiprows=2, index_col=0)
    model = _load_model_lineage()
    ccle_to_model = model.set_index("CCLEName")["ModelID"].to_dict()
    sample_cols = [c for c in gct.columns if c != "Description" and c in ccle_to_model]
    if not sample_cols:
        raise ValueError("No CCLE miRNA columns matched Model.CCLEName")

    mat = pd.DataFrame(
        index=sorted(arm_to_row),
        columns=sorted({ccle_to_model[c] for c in sample_cols}),
        dtype=float,
    )
    for arm, row_id in arm_to_row.items():
        if row_id not in gct.index:
            continue
        row = pd.to_numeric(gct.loc[row_id, sample_cols], errors="coerce")
        row.index = [ccle_to_model[c] for c in sample_cols]
        mat.loc[arm] = np.log2(C.CCLE_MIRNA_EXPR_PSEUDOCOUNT + row.values)
    return mat


def _default_dna_profiles() -> pd.Series:
    if not C.CCLE_DEFAULT_PROFILES.exists():
        raise FileNotFoundError(f"Missing {C.CCLE_DEFAULT_PROFILES}")
    d = pd.read_csv(C.CCLE_DEFAULT_PROFILES)
    dna = d.loc[d["ProfileType"].astype(str).str.lower() == "dna", ["ModelID", "ProfileID"]]
    return dna.drop_duplicates("ModelID").set_index("ModelID")["ProfileID"]


def _resolve_cn_profiles() -> pd.Series:
    """ModelID -> ProfileID present in ``OmicsAbsoluteCNSegmentsProfile.csv``.

    DepMap's default DNA profile is often WGS; PureCN segment files may only
    ship WES for that model (e.g. SKBR3 PR-miPDjx vs PR-44ugos). Prefer default
    when it has segments, else first WES/WGS profile with segments.
    """
    seg_profiles = set(pd.read_csv(C.CCLE_CN_SEGMENTS, usecols=["ProfileID"])["ProfileID"])
    default = _default_dna_profiles()

    op = pd.read_csv(C.CCLE_OMICS_PROFILES)
    dtype = op["Datatype"].astype(str).str.lower()
    dna_op = op.loc[dtype.isin(["wes", "wgs", "copy_number"])].drop_duplicates(["ModelID", "ProfileID"])

    out: Dict[str, str] = {}
    models = sorted(set(default.index) | set(dna_op["ModelID"]))
    for model_id in models:
        candidates: List[str] = []
        if model_id in default.index:
            candidates.append(str(default[model_id]))
        wes = dna_op.loc[
            (dna_op["ModelID"] == model_id) & (dna_op["Datatype"].astype(str).str.lower() == "wes"),
            "ProfileID",
        ]
        wgs = dna_op.loc[
            (dna_op["ModelID"] == model_id) & (dna_op["Datatype"].astype(str).str.lower() == "wgs"),
            "ProfileID",
        ]
        for pid in pd.concat([wes, wgs]).astype(str):
            if pid not in candidates:
                candidates.append(pid)
        hit = next((p for p in candidates if p in seg_profiles), None)
        if hit:
            out[model_id] = hit
    return pd.Series(out)


def _depmap_segments_to_ascat(seg: pd.DataFrame) -> pd.DataFrame:
    from pipeline.CNV.features import _classify_cn_state

    out = seg.copy()
    out["chrom"] = out["Chromosome"].astype(str).map(
        lambda c: f"chr{c}" if not str(c).startswith("chr") else str(c)
    )
    out["start"] = pd.to_numeric(out["Start"], errors="coerce")
    out["end"] = pd.to_numeric(out["End"], errors="coerce")
    out["cn_total"] = pd.to_numeric(out["SegmentAbsoluteCN"], errors="coerce")
    out["loh_flag"] = 0
    out["cn_state"] = out["cn_total"].map(_classify_cn_state)
    return out.dropna(subset=["start", "end"])


def build_sample_locus_cnv(
    *,
    cache_path: Path,
    force: bool = False,
    refresh_missing: bool = False,
    max_models: Optional[int] = None,
) -> pd.DataFrame:
    """Long ModelID × locus CN from DepMap absolute segments (cached)."""
    if cache_path.exists() and not force:
        existing = pd.read_csv(cache_path, sep="\t", low_memory=False)
        if not refresh_missing:
            print(f"[ccle_mirna_cn] reusing cached locus CN: {cache_path}")
            return existing
        prof_map = _resolve_cn_profiles()
        have = set(existing["participant"].unique())
        missing = [m for m in prof_map.index if m not in have]
        if not missing:
            print(f"[ccle_mirna_cn] reusing cached locus CN: {cache_path}")
            return existing
        print(f"[ccle_mirna_cn] appending {len(missing)} models missing from cache …")
        extra = _build_locus_cnv_for_models(missing, prof_map=prof_map)
        combined = pd.concat([existing, extra], ignore_index=True)
        combined.to_csv(cache_path, sep="\t", index=False, compression="gzip")
        print(f"[ccle_mirna_cn] updated cache ({len(combined):,} rows) -> {cache_path}")
        return combined

    from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import (
        assign_mirna_locus_cn_from_segments,
        build_entity_catalog,
        load_mirna_family_map,
        _mirna_locus_catalog,
    )

    if not C.CCLE_CN_SEGMENTS.exists():
        raise FileNotFoundError(f"DepMap CN segments missing: {C.CCLE_CN_SEGMENTS}")

    prof_map = _resolve_cn_profiles()
    models = list(prof_map.index)
    if max_models:
        models = models[: max_models]
    long = _build_locus_cnv_for_models(models, prof_map=prof_map)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    print(f"[ccle_mirna_cn] cached locus CN ({len(long):,} rows) -> {cache_path}")
    return long


def _build_locus_cnv_for_models(
    models: Sequence[str],
    *,
    prof_map: pd.Series,
) -> pd.DataFrame:
    from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import (
        assign_mirna_locus_cn_from_segments,
        build_entity_catalog,
        load_mirna_family_map,
        _mirna_locus_catalog,
    )

    family_map = load_mirna_family_map(C.MIR_FAMILY_INFO)
    catalog = build_entity_catalog(
        loci_path=C.MIRNA_PRECURSOR_LOCI,
        arms_path=C.MIRNA_MATURE_LOCI,
        family_map=family_map,
        panel_arm_ids=None,
    )
    locus_cat = _mirna_locus_catalog(catalog)

    seg_all = pd.read_csv(C.CCLE_CN_SEGMENTS)
    parts: List[pd.DataFrame] = []
    n = len(models)
    for i, model_id in enumerate(models, start=1):
        pid = prof_map[model_id]
        seg = _depmap_segments_to_ascat(seg_all.loc[seg_all["ProfileID"] == pid])
        parts.append(
            assign_mirna_locus_cn_from_segments(
                seg, locus_cat, sample_vial=model_id, participant=model_id
            )
        )
        if i % 200 == 0 or i == n:
            print(f"[ccle_mirna_cn] locus CN: {i}/{n} cell lines …")
    return pd.concat(parts, ignore_index=True)


def _arm_cnv_from_loci(locus_long: pd.DataFrame, arms: Sequence[str]) -> pd.DataFrame:
    from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import (
        assemble_mirna_cnv_long_table,
        build_entity_catalog,
        load_mirna_family_map,
    )

    family_map = load_mirna_family_map(C.MIR_FAMILY_INFO)
    catalog = build_entity_catalog(
        loci_path=C.MIRNA_PRECURSOR_LOCI,
        arms_path=C.MIRNA_MATURE_LOCI,
        family_map=family_map,
        panel_arm_ids=list(arms),
    )
    full = assemble_mirna_cnv_long_table(locus_long, catalog, loci_path=C.MIRNA_PRECURSOR_LOCI)
    return full.loc[full["entity_level"] == "arm"].copy()


def concordance_by_stratum(
    arm_cn: pd.DataFrame,
    expr: pd.DataFrame,
    *,
    min_n: int = C.CCLE_CONCORDANCE_MIN_N,
    strata: Optional[Dict[str, Set[str]]] = None,
    rho_kind: RhoKind = "partial",
) -> pd.DataFrame:
    """Per-arm Spearman(CN, log2 expr) across CCLE strata (marginal or partial)."""
    strata = strata or _stratum_model_sets(min_n)

    cn_wide = (
        arm_cn.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack("participant")
    )

    rows: List[dict] = []
    for arm in expr.index:
        if arm not in cn_wide.index:
            continue
        e = expr.loc[arm].dropna()
        for stratum, model_ids in strata.items():
            mids = sorted(set(e.index) & set(cn_wide.columns) & model_ids)
            if len(mids) < min_n:
                continue
            cn = cn_wide.loc[arm, mids].astype(float)
            ex = e.reindex(mids).astype(float)
            ok = cn.notna() & ex.notna()
            if ok.sum() < min_n:
                continue
            mids_ok = [m for m in mids if ok.get(m, False)]
            cov = (
                _ccle_covariates_for_stratum(stratum, mids_ok)
                if rho_kind == "partial"
                else None
            )
            rho, p, n = _correlation_pair(
                cn.loc[mids_ok], ex.loc[mids_ok], cov, rho_kind=rho_kind,
            )
            row = {
                "arm": arm,
                "stratum": stratum,
                "n": int(n),
                "rho_kind": rho_kind,
                "spearman_rho": round(float(rho), 4) if pd.notna(rho) else np.nan,
                "spearman_p": float(p) if pd.notna(p) else np.nan,
            }
            if rho_kind == "partial":
                row["partial_rho"] = row["spearman_rho"]
                row["partial_p"] = row["spearman_p"]
            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    pcol = "partial_p" if rho_kind == "partial" else "spearman_p"
    qcol = "partial_q" if rho_kind == "partial" else "spearman_q"
    for stratum in df["stratum"].unique():
        mask = df["stratum"] == stratum
        valid = df.loc[mask, pcol].notna()
        if valid.any():
            df.loc[mask & valid, qcol] = S.bh_fdr(df.loc[mask & valid, pcol].values)
    sort_col = "partial_rho" if rho_kind == "partial" else "spearman_rho"
    return df.sort_values(["stratum", sort_col], ascending=[True, False])


def compare_tcga_ccle(
    conc_ccle: pd.DataFrame,
    tcga_path: Optional[Path] = None,
    *,
    lineage_cols: Optional[Sequence[str]] = None,
    rho_kind: RhoKind = "partial",
) -> pd.DataFrame:
    tcga_path = Path(tcga_path or C.MIRNA_LOCUS_CNV_DIR / "mirna_cnv_expr_concordance.tsv")
    if not tcga_path.exists():
        return pd.DataFrame()
    tcga = pd.read_csv(tcga_path, sep="\t")
    tcga = _ensure_partial_q(tcga)
    if rho_kind == "partial":
        tcga = tcga.rename(columns={
            "n": "tcga_n",
            "partial_rho": "tcga_rho",
            "partial_q": "tcga_q",
        })
    else:
        tcga = tcga.rename(columns={
            "n": "tcga_n",
            "spearman_rho": "tcga_rho",
            "spearman_q": "tcga_q",
        })
    out = tcga[["arm", "tcga_n", "tcga_rho", "tcga_q"]].copy()
    lineage_cols = lineage_cols or ["all", "breast"]
    rho_col = "partial_rho" if rho_kind == "partial" else "spearman_rho"
    q_col = "partial_q" if rho_kind == "partial" else "spearman_q"
    if rho_col not in conc_ccle.columns:
        rho_col, q_col = "spearman_rho", "spearman_q"
    for stratum in lineage_cols:
        sub = conc_ccle.loc[conc_ccle["stratum"] == stratum]
        if sub.empty:
            continue
        sfx = stratum if stratum != "all" else "all"
        merged = sub.rename(columns={
            "n": f"ccle_{sfx}_n",
            rho_col: f"ccle_{sfx}_rho",
            q_col: f"ccle_{sfx}_q",
        })
        out = out.merge(
            merged[["arm", f"ccle_{sfx}_n", f"ccle_{sfx}_rho", f"ccle_{sfx}_q"]],
            on="arm",
            how="left",
        )
    if "ccle_all_rho" in out.columns:
        out["rho_delta_ccle_minus_tcga"] = (out["ccle_all_rho"] - out["tcga_rho"]).round(4)
        out["sign_agree_all"] = np.sign(out["tcga_rho"]) == np.sign(out["ccle_all_rho"])
    if "ccle_breast_rho" in out.columns:
        out["sign_agree_breast"] = np.sign(out["tcga_rho"]) == np.sign(out["ccle_breast_rho"])
    return out.sort_values("tcga_rho", ascending=False)


def _print_concordance_summary(conc: pd.DataFrame, *, rho_kind: RhoKind = "partial") -> None:
    if conc.empty:
        return
    rho_col = "partial_rho" if rho_kind == "partial" else "spearman_rho"
    q_col = "partial_q" if rho_kind == "partial" else "spearman_q"
    for stratum in conc["stratum"].unique():
        sub = conc.loc[
            (conc["stratum"] == stratum)
            & (conc[q_col] < C.FDR_ALPHA)
            & (conc[rho_col] > 0)
        ]
        n_tested = int((conc["stratum"] == stratum).sum())
        top = conc.loc[conc["stratum"] == stratum].sort_values(rho_col, ascending=False).head(3)
        tops = ", ".join(f"{r.arm} ρ={getattr(r, rho_col):.3f}" for r in top.itertuples())
        print(
            f"[ccle_mirna_cn] {stratum} ({rho_kind}): q<{C.FDR_ALPHA} pos "
            f"{len(sub)}/{n_tested}; top: {tops}"
        )


def run(
    *,
    out_dir: Path = C.CCLE_CN_CONCORDANCE_DIR,
    force_cnv: bool = False,
    max_models: Optional[int] = None,
    panel: PanelMode = "high_evidence",
    match_mode: MatchMode = "alias",
    rho_kind: RhoKind = "partial",
) -> pd.DataFrame:
    if panel != "high_evidence":
        out_dir = out_dir / panel
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)

    locus_cache = C.CCLE_CN_CONCORDANCE_DIR / "tables" / "sample_locus_cn.tsv.gz"
    panel_arms = resolve_panel_arms(panel, rho_kind=rho_kind)
    arm_to_row, map_notes = build_ccle_arm_expression_map(panel_arms, match_mode=match_mode)
    map_notes.to_csv(out_dir / "ccle_arm_expression_map.tsv", sep="\t", index=False)
    mapped_arms = sorted(arm_to_row)
    print(
        f"[ccle_mirna_cn] panel={panel} match={match_mode} "
        f"arms={len(panel_arms)} mapped={len(mapped_arms)}"
    )
    if not mapped_arms:
        raise ValueError("No arms mapped to CCLE expression; check match_mode / panel")

    expr = _ccle_expr_matrix(arm_to_row)
    locus_long = build_sample_locus_cnv(
        cache_path=locus_cache,
        force=force_cnv,
        max_models=max_models,
    )
    locus_long = locus_long.loc[locus_long["participant"].isin(expr.columns)].copy()
    arm_cn = _arm_cnv_from_loci(locus_long, mapped_arms)
    arm_cn.to_csv(out_dir / "tables" / "sample_arm_cn.tsv.gz", sep="\t", index=False, compression="gzip")
    cn_arms = set(arm_cn["entity_label"].unique())
    missing_cn = sorted(set(mapped_arms) - cn_arms)
    if missing_cn:
        print(f"[ccle_mirna_cn] warning: {len(missing_cn)} mapped arms lack CN (e.g. {missing_cn[:5]})")

    test_arms = sorted(cn_arms & set(mapped_arms))
    conc = concordance_by_stratum(arm_cn, expr.loc[test_arms], rho_kind=rho_kind)
    conc.to_csv(out_dir / "ccle_mirna_cn_expr_concordance.tsv", sep="\t", index=False)

    compare = compare_tcga_ccle(
        conc,
        lineage_cols=["all", "breast", "lung", "ovary_fallopian_tube"],
        rho_kind=rho_kind,
    )
    if not compare.empty:
        compare.to_csv(out_dir / "tcga_ccle_concordance_compare.tsv", sep="\t", index=False)

    _print_concordance_summary(conc, rho_kind=rho_kind)

    manifest = {
        "module": "mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "panel": panel,
        "match_mode": match_mode,
        "cnv_source": "DepMap OmicsAbsoluteCNSegmentsProfile (PureCN absolute), locus max-overlap",
        "expression_source": str(C.CCLE_MIRNA_GCT),
        "expression_transform": f"log2({C.CCLE_MIRNA_EXPR_PSEUDOCOUNT} + nSolver)",
        "n_panel_arms": len(panel_arms),
        "n_ccle_mapped_arms": len(mapped_arms),
        "n_arms_with_cn": len(test_arms),
        "n_cell_lines_expr": int(expr.shape[1]),
        "lineage_strata": list(_stratum_model_sets(C.CCLE_CONCORDANCE_MIN_N).keys()),
        "concordance_min_n": C.CCLE_CONCORDANCE_MIN_N,
        "n_concordance_rows": int(len(conc)),
        "rho_kind": rho_kind,
        "tcga_partial_covariates": list(C.CONFOUNDER_NUMERIC),
        "ccle_partial_covariates": (
            "OncotreeLineage dummies (all strata); PAM50 dummies (breast strata)"
        ),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[ccle_mirna_cn] wrote outputs under {out_dir}")
    return conc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.CCLE_CN_CONCORDANCE_DIR)
    ap.add_argument("--force-cnv", action="store_true", help="Re-extract locus CN from DepMap segments")
    ap.add_argument("--max-models", type=int, default=None, help="Debug: limit cell lines")
    ap.add_argument(
        "--panel",
        choices=("high_evidence", "focus_8q", "tcga_sig"),
        default="high_evidence",
        help="miRNA arm panel (default: high_evidence Hallmark arms)",
    )
    ap.add_argument(
        "--match-mode",
        choices=("alias", "exact", "precursor_stem"),
        default="alias",
        help="CCLE GCT Description matching (default: alias incl. 151a→151-3p)",
    )
    ap.add_argument(
        "--rho-kind",
        choices=("partial", "marginal"),
        default="partial",
        help="Spearman on raw pairs (marginal) or OLS-residualized (partial; default)",
    )
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(
        out_dir=args.out_dir,
        force_cnv=args.force_cnv,
        max_models=args.max_models,
        panel=args.panel,
        match_mode=args.match_mode,
        rho_kind=args.rho_kind,
    )


if __name__ == "__main__":
    main()
