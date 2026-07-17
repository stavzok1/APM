"""SV DUP/DEL overlap at miRNA hairpin loci and dense clusters.

Scans the SV pipeline strict set (``07_final_sv_with_fimo`` by default) for
``SVTYPE`` in {DEL, DUP} with **direct overlap** of a pre-miRNA locus
(``mir_hits`` ``region_hit == overlaps`` or ``overlap_bp > 0``). Joins to the
cached ASCAT3 segment CNV table to flag CNV-only, SV-only, and both layers.

Outputs under ``output/mirna_locus_cnv/sv_overlap/``:
  - ``mirna_locus_sv_hits.tsv.gz`` — long (participant × SV × locus)
  - ``mirna_locus_sv_recurrence.tsv`` — per-locus cohort recurrence
  - ``mirna_cluster_sv_recurrence.tsv`` — dense-cluster rollup
  - ``mirna_locus_cnv_sv_layers.tsv.gz`` — participant × locus layer class
  - ``mirna_locus_sv_both_patterns.tsv`` / ``mirna_locus_sv_only_patterns.tsv`` — layer QC summaries
  - ``mirna_locus_sv_both_concordance_by_cn_state.tsv`` — DUP/DEL direction vs CNV state
  - ``review_queues/*.tsv.gz`` — IGV review queues (focal sv_only, discordant both, mixed both)
  - ``chr19_megacluster/`` — megacluster locus + participant dossier
  - ``method_manifest.json``

Future work (not implemented here): adjudicate SV vs CNV when both exist,
partial overlaps, and BND neojunction at miRNA loci.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps import (
    CHR19_MEGACLUSTER_ID,
    _loci_in_cluster_block,
    build_dense_locus_clusters,
)

_SV_COLS = [
    "participant",
    "sample_vial",
    "chrom",
    "pos",
    "END",
    "SVLEN",
    "id",
    "SVTYPE",
    "mir_hits",
    "tumor_alt",
    "TUMOR_PR_ref",
    "TUMOR_SR_ref",
]
_VAF_LOOKUP_COLS = ["participant", "sample_vial", "id", "tumor_alt", "TUMOR_PR_ref", "TUMOR_SR_ref"]
_CNV_ALTERED = {"gain", "amp", "loss", "deep_del"}
_NAN_RE = re.compile(r"\bnan\b", re.IGNORECASE)


def _parse_listish(v: object) -> List[dict]:
    if v is None:
        return []
    try:
        if isinstance(v, float) and np.isnan(v):
            return []
    except (TypeError, ValueError):
        pass
    if isinstance(v, (list, tuple)):
        return [x for x in v if isinstance(x, dict)]
    if isinstance(v, str):
        s = v.strip()
        if not s or s == "[]" or s.lower() == "nan":
            return []
        try:
            parsed = ast.literal_eval(_NAN_RE.sub("None", s))
        except (ValueError, SyntaxError):
            return []
        if isinstance(parsed, list):
            return [x for x in parsed if isinstance(x, dict)]
    return []


def _mir_overlap_hits(mir_hits: object) -> List[dict]:
    out: List[dict] = []
    for h in _parse_listish(mir_hits):
        overlap_bp = int(pd.to_numeric(h.get("overlap_bp"), errors="coerce") or 0)
        region = str(h.get("region_hit") or "")
        if overlap_bp <= 0 and region != "overlaps":
            continue
        lid = str(h.get("pre_gene_id") or h.get("locus_id") or "").strip()
        if not lid:
            continue
        out.append(
            {
                "locus_id": lid,
                "locus_name": str(h.get("pre_gene_name") or h.get("gene_name") or ""),
                "mirbase_mature_id": str(h.get("mirbase_mature_id") or ""),
                "overlap_bp": overlap_bp,
                "overlap_percent": round(float(pd.to_numeric(h.get("overlap_percent"), errors="coerce") or 0), 2),
                "region_hit": region,
            }
        )
    return out


def _sv_span_bp(row: pd.Series) -> float:
    svlen = pd.to_numeric(row.get("SVLEN"), errors="coerce")
    if pd.notna(svlen) and abs(float(svlen)) > 0:
        return abs(float(svlen))
    end = pd.to_numeric(row.get("END"), errors="coerce")
    pos = pd.to_numeric(row.get("pos"), errors="coerce")
    if pd.notna(end) and pd.notna(pos):
        return abs(float(end) - float(pos))
    return np.nan


def _tumor_vaf_from_manta(row: object) -> float:
    if isinstance(row, pd.Series):
        d = row
    elif hasattr(row, "_asdict"):
        d = row._asdict()
    else:
        d = row
    alt = float(pd.to_numeric(d.get("tumor_alt"), errors="coerce") or 0)
    ref = float(pd.to_numeric(d.get("TUMOR_PR_ref"), errors="coerce") or 0) + float(
        pd.to_numeric(d.get("TUMOR_SR_ref"), errors="coerce") or 0
    )
    tot = alt + ref
    if tot <= 0:
        return float("nan")
    return float(min(1.0, alt / tot))


def _cn_minor_estimate(copy_number: object, loh_flag: object, cn_state: object) -> float:
    if str(loh_flag).lower() in ("true", "1", "yes") or str(cn_state) in ("loss", "deep_del"):
        return 0.0
    cn = pd.to_numeric(copy_number, errors="coerce")
    if np.isfinite(cn) and cn > 0:
        return float(max(0.0, cn / 2.0))
    return 1.0


def enrich_sv_hits_with_vaf(
    hits: pd.DataFrame,
    *,
    layers_path: Optional[Path] = None,
    sv_dir: Optional[Path] = None,
    save_path: Optional[Path] = None,
) -> pd.DataFrame:
    """Attach Manta read VAF + CPE-adjusted CCF and locus-ASCAT expected-VAF ratio."""
    from analysis.variants.sv.ascat_vaf import expected_vaf_from_ascat, purity_adjusted_vaf

    if hits.empty:
        return hits

    out = hits.copy()
    sv_dir = Path(sv_dir or C.SV_FINAL_DIR)
    layers_path = Path(layers_path or C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz")

    if "tumor_vaf" not in out.columns or out["tumor_vaf"].isna().mean() > 0.5:
        wanted = out.groupby("participant")["sv_id"].apply(lambda s: set(map(str, s))).to_dict()
        vaf_rows: List[dict] = []
        for csv_path in sorted(sv_dir.glob("*_strict_sv_set.csv")):
            head = pd.read_csv(csv_path, nrows=0)
            use = [c for c in _VAF_LOOKUP_COLS if c in head.columns]
            if "id" not in use:
                continue
            df = pd.read_csv(csv_path, usecols=use, low_memory=False)
            if df.empty:
                continue
            part = str(df["participant"].iloc[0]) if "participant" in df.columns else ""
            if not part:
                stem = csv_path.name.split("_strict_sv_set.csv")[0]
                part = "-".join(stem.split("-")[:3]) if stem.startswith("TCGA-") else stem
            ids = wanted.get(part)
            if not ids:
                continue
            sub = df.loc[df["id"].astype(str).isin(ids)].drop_duplicates("id")
            for row in sub.itertuples(index=False):
                d = row._asdict()
                vaf = _tumor_vaf_from_manta(d)
                vaf_rows.append(
                    {
                        "participant": part,
                        "sv_id": str(d.get("id")),
                        "tumor_vaf": round(vaf, 6) if np.isfinite(vaf) else np.nan,
                        "tumor_alt_support": pd.to_numeric(d.get("tumor_alt"), errors="coerce"),
                        "tumor_pr_ref": pd.to_numeric(d.get("TUMOR_PR_ref"), errors="coerce"),
                        "tumor_sr_ref": pd.to_numeric(d.get("TUMOR_SR_ref"), errors="coerce"),
                    }
                )
        if vaf_rows:
            vaf_df = pd.DataFrame(vaf_rows).drop_duplicates(["participant", "sv_id"])
            drop_cols = [c for c in ("tumor_vaf", "tumor_pr_ref", "tumor_sr_ref") if c in out.columns]
            out = out.drop(columns=drop_cols, errors="ignore")
            out = out.merge(vaf_df, on=["participant", "sv_id"], how="left", suffixes=("", "_vaf"))
            if "tumor_alt_support_vaf" in out.columns:
                out["tumor_alt_support"] = out["tumor_alt_support"].fillna(out["tumor_alt_support_vaf"])
                out = out.drop(columns=["tumor_alt_support_vaf"])
    elif "tumor_alt_support" in out.columns and "tumor_vaf" in out.columns:
        miss = out["tumor_vaf"].isna()
        out.loc[miss, "tumor_vaf"] = [
            _tumor_vaf_from_manta(out.loc[idx]) if pd.notna(out.loc[idx, "tumor_alt_support"]) else np.nan
            for idx in out.index[miss]
        ]

    if layers_path.exists():
        lay = pd.read_csv(
            layers_path,
            sep="\t",
            usecols=["participant", "entity_id", "CPE", "copy_number", "loh_flag", "cn_state"],
            low_memory=False,
        )
        lay = lay.rename(columns={"entity_id": "locus_id"}).drop_duplicates(["participant", "locus_id"])
        out = out.merge(lay, on=["participant", "locus_id"], how="left", suffixes=("", "_lay"))

    cpe = pd.to_numeric(out.get("CPE"), errors="coerce")
    vaf = pd.to_numeric(out.get("tumor_vaf"), errors="coerce")
    out["tumor_vaf_adj"] = [
        purity_adjusted_vaf(float(a), float(p), method="ccf")
        if np.isfinite(a) and np.isfinite(p)
        else np.nan
        for a, p in zip(vaf, cpe)
    ]
    cn = pd.to_numeric(out.get("copy_number"), errors="coerce")
    cn_minor = [
        _cn_minor_estimate(c, lf, cs)
        for c, lf, cs in zip(cn, out.get("loh_flag", pd.Series([np.nan] * len(out))), out.get("cn_state", pd.Series([""] * len(out))))
    ]
    out["tumor_vaf_expected_locus"] = [
        expected_vaf_from_ascat(float(p), float(ct), float(cm))
        if np.isfinite(p) and np.isfinite(ct)
        else np.nan
        for p, ct, cm in zip(cpe, cn, cn_minor)
    ]
    vadj = pd.to_numeric(out["tumor_vaf_adj"], errors="coerce")
    exp = pd.to_numeric(out["tumor_vaf_expected_locus"], errors="coerce")
    out["tumor_vaf_ratio_adj"] = vadj / exp.replace(0, np.nan)

    if save_path:
        out.to_csv(save_path, sep="\t", index=False, compression="gzip")
    return out


def _scan_one_sample(csv_path: str) -> List[dict]:
    path = Path(csv_path)
    head = pd.read_csv(path, nrows=0)
    use = [c for c in _SV_COLS if c in head.columns]
    df = pd.read_csv(path, usecols=use, low_memory=False)
    if df.empty:
        return []
    df = df.loc[df["SVTYPE"].astype(str).isin(["DEL", "DUP"])].copy()
    rows: List[dict] = []
    for r in df.itertuples(index=False):
        row = r._asdict()
        hits = _mir_overlap_hits(row.get("mir_hits"))
        if not hits:
            continue
        span = _sv_span_bp(pd.Series(row))
        vaf = _tumor_vaf_from_manta(row)
        for h in hits:
            rows.append(
                {
                    "participant": str(row.get("participant") or ""),
                    "sample_vial": str(row.get("sample_vial") or ""),
                    "sv_id": str(row.get("id") or ""),
                    "svtype": str(row.get("SVTYPE") or ""),
                    "sv_chrom": str(row.get("chrom") or ""),
                    "sv_start": int(pd.to_numeric(row.get("pos"), errors="coerce") or 0),
                    "sv_end": int(pd.to_numeric(row.get("END"), errors="coerce") or 0)
                    if pd.notna(row.get("END"))
                    else np.nan,
                    "sv_span_bp": round(span, 0) if np.isfinite(span) else np.nan,
                    "tumor_alt_support": pd.to_numeric(row.get("tumor_alt"), errors="coerce"),
                    "tumor_vaf": round(vaf, 6) if np.isfinite(vaf) else np.nan,
                    "tumor_pr_ref": pd.to_numeric(row.get("TUMOR_PR_ref"), errors="coerce"),
                    "tumor_sr_ref": pd.to_numeric(row.get("TUMOR_SR_ref"), errors="coerce"),
                    **h,
                }
            )
    return rows


def _assign_cluster_ids(locus_map: pd.DataFrame, clusters: pd.DataFrame) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for _, cl in clusters.iterrows():
        cid = f"{cl['chrom']}:{cl['cluster_start_mb']:.3f}-{cl['cluster_end_mb']:.3f}"
        block = _loci_in_cluster_block(locus_map, cl)
        for lid in block["locus_id"].astype(str):
            mapping[lid] = cid
    return mapping


def build_locus_recurrence(hits: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return pd.DataFrame()
    rows: List[dict] = []
    for lid, grp in hits.groupby("locus_id"):
        parts = grp.groupby("participant")["svtype"].apply(lambda s: set(s.astype(str))).to_dict()
        n_del = sum(1 for s in parts.values() if "DEL" in s)
        n_dup = sum(1 for s in parts.values() if "DUP" in s)
        meta = grp.iloc[0]
        rows.append(
            {
                "locus_id": lid,
                "locus_name": meta.get("locus_name", ""),
                "cluster_id": meta.get("cluster_id", ""),
                "n_participants_any": len(parts),
                "n_participants_DEL": n_del,
                "n_participants_DUP": n_dup,
                "n_SV_events": int(len(grp.drop_duplicates(["participant", "sv_id", "locus_id"]))),
                "median_sv_span_bp": round(
                    float(pd.to_numeric(grp["sv_span_bp"], errors="coerce").median()), 0
                ),
                "median_overlap_bp": round(
                    float(pd.to_numeric(grp["overlap_bp"], errors="coerce").median()), 0
                ),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values("n_participants_any", ascending=False)


def build_cluster_recurrence(hits: pd.DataFrame) -> pd.DataFrame:
    sub = hits.loc[hits["cluster_id"].astype(str).str.len() > 0].copy()
    if sub.empty:
        return pd.DataFrame()
    rows: List[dict] = []
    for cid, grp in sub.groupby("cluster_id"):
        loci = grp["locus_id"].nunique()
        parts = grp.groupby("participant")["svtype"].apply(lambda s: set(s.astype(str))).to_dict()
        n_del = sum(1 for s in parts.values() if "DEL" in s)
        n_dup = sum(1 for s in parts.values() if "DUP" in s)
        rows.append(
            {
                "cluster_id": cid,
                "chrom": cid.split(":")[0],
                "n_hairpin_loci_in_cluster": int(loci),
                "n_participants_any": len(parts),
                "n_participants_DEL": n_del,
                "n_participants_DUP": n_dup,
                "n_SV_events": int(len(grp.drop_duplicates(["participant", "sv_id", "locus_id"]))),
            }
        )
    return pd.DataFrame(rows).sort_values("n_participants_any", ascending=False)


def build_cnv_sv_layers(hits: pd.DataFrame, cnv_long: pd.DataFrame, cluster_map: Dict[str, str]) -> pd.DataFrame:
    """Participant × locus dosage layer: cnv_only | sv_only | both | neither."""
    loc = cnv_long.loc[cnv_long["entity_level"] == "locus"].drop_duplicates(
        subset=["participant", "entity_id"], keep="first"
    ).copy()
    loc["locus_id"] = loc["entity_id"].astype(str)
    loc["cluster_id"] = loc["locus_id"].map(cluster_map).fillna("")

    if hits.empty:
        loc["has_sv_DEL"] = False
        loc["has_sv_DUP"] = False
        loc["has_sv_dup_del"] = False
    else:
        sv_flags = (
            hits.assign(
                is_del=hits["svtype"].astype(str) == "DEL",
                is_dup=hits["svtype"].astype(str) == "DUP",
            )
            .groupby(["participant", "locus_id"], as_index=False)
            .agg(has_sv_DEL=("is_del", "any"), has_sv_DUP=("is_dup", "any"))
        )
        loc = loc.merge(sv_flags, on=["participant", "locus_id"], how="left")
        loc["has_sv_DEL"] = loc["has_sv_DEL"].eq(True).fillna(False)
        loc["has_sv_DUP"] = loc["has_sv_DUP"].eq(True).fillna(False)
    loc["has_sv_dup_del"] = loc["has_sv_DEL"] | loc["has_sv_DUP"]
    loc["cnv_altered"] = loc["cn_state"].astype(str).isin(_CNV_ALTERED)

    def _layer(row) -> str:
        cnv = bool(row["cnv_altered"])
        sv = bool(row["has_sv_dup_del"])
        if cnv and sv:
            return "both"
        if cnv:
            return "cnv_only"
        if sv:
            return "sv_only"
        return "neither"

    loc["dosage_layer"] = loc.apply(_layer, axis=1)
    return loc.rename(columns={"entity_label": "locus_name"})


def _expected_sv_type(cn_state: str) -> Optional[str]:
    st = str(cn_state)
    if st in ("amp", "gain"):
        return "DUP"
    if st in ("loss", "deep_del"):
        return "DEL"
    return None


def _segment_length_mb(seg_label: object) -> float:
    if not seg_label or pd.isna(seg_label) or ":" not in str(seg_label):
        return np.nan
    try:
        _, coords = str(seg_label).split(":", 1)
        start_s, end_s = coords.split("-", 1)
        return (int(end_s) - int(start_s)) / 1e6
    except (ValueError, TypeError):
        return np.nan


def _both_participant_locus_table(layers: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    both = layers.loc[layers["dosage_layer"] == "both", ["participant", "locus_id", "cn_state", "PAM50_final"]]
    merged = both.merge(hits, on=["participant", "locus_id"], how="inner")
    if merged.empty:
        return pd.DataFrame()

    rows: List[dict] = []
    for (participant, locus_id), grp in merged.groupby(["participant", "locus_id"]):
        cn = str(grp["cn_state"].iloc[0])
        exp = _expected_sv_type(cn)
        svtypes = grp["svtype"].astype(str).value_counts()
        has_dup = int(svtypes.get("DUP", 0)) > 0
        has_del = int(svtypes.get("DEL", 0)) > 0
        conc = (exp == "DUP" and has_dup) or (exp == "DEL" and has_del)
        discord_only = (exp == "DUP" and has_del and not has_dup) or (exp == "DEL" and has_dup and not has_del)
        span_mb = pd.to_numeric(grp["sv_span_bp"], errors="coerce") / 1e6
        rows.append(
            {
                "participant": participant,
                "locus_id": locus_id,
                "locus_name": str(grp["locus_name"].iloc[0]),
                "cluster_id": str(grp["cluster_id"].iloc[0]) if "cluster_id" in grp.columns else "",
                "PAM50_final": str(grp["PAM50_final"].iloc[0]),
                "cn_state": cn,
                "expected_svtype": exp or "",
                "n_sv_events": int(grp["sv_id"].nunique()),
                "has_SV_DUP": has_dup,
                "has_SV_DEL": has_del,
                "any_concordant_sv": conc,
                "discordant_sv_only": discord_only,
                "mixed_dup_and_del": has_dup and has_del,
                "median_sv_span_mb": round(float(span_mb.median()), 3) if span_mb.notna().any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _sv_only_participant_locus_table(layers: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    svo = layers.loc[layers["dosage_layer"] == "sv_only", ["participant", "locus_id", "PAM50_final", "cluster_id"]]
    merged = svo.merge(hits, on=["participant", "locus_id"], how="inner")
    if merged.empty:
        return pd.DataFrame()

    rows: List[dict] = []
    for (participant, locus_id), grp in merged.groupby(["participant", "locus_id"]):
        svtypes = grp["svtype"].astype(str).value_counts()
        span_mb = pd.to_numeric(grp["sv_span_bp"], errors="coerce") / 1e6
        rows.append(
            {
                "participant": participant,
                "locus_id": locus_id,
                "locus_name": str(grp["locus_name"].iloc[0]),
                "cluster_id": str(grp["cluster_id"].iloc[0] if "cluster_id" in grp.columns else ""),
                "PAM50_final": str(grp["PAM50_final"].iloc[0]),
                "n_sv_events": int(grp["sv_id"].nunique()),
                "has_SV_DUP": int(svtypes.get("DUP", 0)) > 0,
                "has_SV_DEL": int(svtypes.get("DEL", 0)) > 0,
                "median_sv_span_mb": round(float(span_mb.median()), 3) if span_mb.notna().any() else np.nan,
                "focal_sv_le_1mb": bool(span_mb.median() <= 1.0) if span_mb.notna().any() else False,
                "median_tumor_alt": round(
                    float(pd.to_numeric(grp["tumor_alt_support"], errors="coerce").median()), 2
                ),
            }
        )
    return pd.DataFrame(rows)


def _summarize_both_patterns(pl: pd.DataFrame, *, stratum: str) -> dict:
    if pl.empty:
        return {"stratum": stratum, "n_participant_locus": 0}
    return {
        "stratum": stratum,
        "n_participant_locus": int(len(pl)),
        "pct_any_concordant_sv": round(100.0 * float(pl["any_concordant_sv"].mean()), 2),
        "pct_discordant_sv_only": round(100.0 * float(pl["discordant_sv_only"].mean()), 2),
        "pct_mixed_dup_and_del": round(100.0 * float(pl["mixed_dup_and_del"].mean()), 2),
        "median_sv_span_mb": round(float(pl["median_sv_span_mb"].median()), 2),
        "median_n_sv_events": round(float(pl["n_sv_events"].median()), 2),
    }


def _summarize_sv_only_patterns(pl: pd.DataFrame, *, stratum: str) -> dict:
    if pl.empty:
        return {"stratum": stratum, "n_participant_locus": 0}
    return {
        "stratum": stratum,
        "n_participant_locus": int(len(pl)),
        "pct_focal_sv_le_1mb": round(100.0 * float(pl["focal_sv_le_1mb"].mean()), 2),
        "pct_sv_DEL_only": round(100.0 * float((pl["has_SV_DEL"] & ~pl["has_SV_DUP"]).mean()), 2),
        "pct_sv_DUP_only": round(100.0 * float((pl["has_SV_DUP"] & ~pl["has_SV_DEL"]).mean()), 2),
        "pct_mixed_dup_and_del": round(100.0 * float((pl["has_SV_DUP"] & pl["has_SV_DEL"]).mean()), 2),
        "median_sv_span_mb": round(float(pl["median_sv_span_mb"].median()), 2),
        "median_tumor_alt": round(float(pl["median_tumor_alt"].median()), 2),
        "median_n_sv_events": round(float(pl["n_sv_events"].median()), 2),
    }


def build_layer_pattern_summaries(
    layers: pd.DataFrame,
    hits: pd.DataFrame,
    locus_map: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (both_summary, sv_only_summary, sv_only_top_loci, both_by_cn_state)."""
    both_pl = _both_participant_locus_table(layers, hits)
    svo_pl = _sv_only_participant_locus_table(layers, hits)

    both_rows = [_summarize_both_patterns(both_pl, stratum="cohort")]
    for cn, sub in both_pl.groupby("cn_state"):
        both_rows.append(_summarize_both_patterns(sub, stratum=f"cn_state={cn}"))
    for pam, sub in both_pl.groupby("PAM50_final"):
        if str(pam) and str(pam) != "nan":
            both_rows.append(_summarize_both_patterns(sub, stratum=f"PAM50={pam}"))
    both_summary = pd.DataFrame(both_rows)

    svo_rows = [_summarize_sv_only_patterns(svo_pl, stratum="cohort")]
    for pam, sub in svo_pl.groupby("PAM50_final"):
        if str(pam) and str(pam) != "nan":
            svo_rows.append(_summarize_sv_only_patterns(sub, stratum=f"PAM50={pam}"))
    svo_summary = pd.DataFrame(svo_rows)

    top_svo = pd.DataFrame()
    if not svo_pl.empty:
        top = (
            svo_pl.groupby(["locus_id", "locus_name"], as_index=False)
            .agg(
                n_participant_locus=("participant", "nunique"),
                median_sv_span_mb=("median_sv_span_mb", "median"),
                pct_focal_le_1mb=("focal_sv_le_1mb", "mean"),
            )
            .sort_values("n_participant_locus", ascending=False)
        )
        meta = locus_map.set_index("locus_id")[["chrom", "delta_vs_diploid"]]
        top_svo = top.merge(meta, left_on="locus_id", right_index=True, how="left")
        top_svo["pct_focal_le_1mb"] = (100.0 * top_svo["pct_focal_le_1mb"]).round(2)

    cn_detail = pd.DataFrame()
    if not both_pl.empty:
        cn_detail = (
            both_pl.groupby("cn_state", as_index=False)
            .agg(
                n_participant_locus=("participant", "size"),
                pct_any_concordant_sv=("any_concordant_sv", lambda s: round(100 * s.mean(), 2)),
                pct_discordant_sv_only=("discordant_sv_only", lambda s: round(100 * s.mean(), 2)),
                pct_mixed_dup_and_del=("mixed_dup_and_del", lambda s: round(100 * s.mean(), 2)),
                median_sv_span_mb=("median_sv_span_mb", "median"),
            )
            .sort_values("n_participant_locus", ascending=False)
        )

    return both_summary, svo_summary, top_svo, cn_detail


def _sv_interval_label(chrom: object, start: object, end: object) -> str:
    if pd.isna(chrom) or pd.isna(start):
        return ""
    c = str(chrom)
    s = int(pd.to_numeric(start, errors="coerce") or 0)
    e = int(pd.to_numeric(end, errors="coerce") or s)
    if e < s:
        s, e = e, s
    return f"{c}:{s}-{e}"


def _igv_window(chrom: object, start: object, end: object, pad_mb: float) -> str:
    if pd.isna(chrom) or pd.isna(start) or pd.isna(end):
        return ""
    pad = int(max(0.25, float(pad_mb)) * 1e6)
    s = max(0, int(start) - pad)
    e = int(end) + pad
    return f"{chrom}:{s}-{e}"


def _enrich_review_rows(
    queue: pd.DataFrame,
    *,
    locus_map: pd.DataFrame,
    layers: pd.DataFrame,
    hits: pd.DataFrame,
    review_queue: str,
) -> pd.DataFrame:
    if queue.empty:
        return queue

    meta = locus_map[
        ["locus_id", "locus_name", "chrom", "start", "end", "delta_vs_diploid"]
    ].drop_duplicates("locus_id")
    lay = layers.copy()
    lay["locus_id"] = lay["entity_id"].astype(str)
    lay = lay.drop_duplicates(["participant", "locus_id"])
    lay_cols = [
        "participant",
        "locus_id",
        "sample_vial",
        "copy_number",
        "cn_state",
        "overlap_segment",
        "dosage_layer",
    ]
    if "PAM50_final" in lay.columns:
        lay_cols.append("PAM50_final")
    if "CPE" in lay.columns:
        lay_cols.append("CPE")

    def _sv_intervals(g: pd.DataFrame) -> str:
        parts = []
        for row in g.drop_duplicates("sv_id").itertuples(index=False):
            parts.append(_sv_interval_label(row.sv_chrom, row.sv_start, row.sv_end))
        return ";".join(p for p in parts if p)[:800]

    sv_agg = (
        hits.groupby(["participant", "locus_id"], as_index=False)
        .agg(
            sample_vial=("sample_vial", lambda s: next((x for x in s.astype(str) if x and x != "nan"), "")),
            sv_ids=("sv_id", lambda s: ";".join(sorted(set(map(str, s)))[:500])),
            svtypes=("svtype", lambda s: ",".join(sorted(set(map(str, s))))),
            median_sv_span_mb=(
                "sv_span_bp",
                lambda s: round(float(pd.to_numeric(s, errors="coerce").median()) / 1e6, 3),
            ),
            median_tumor_alt=("tumor_alt_support", "median"),
            median_tumor_vaf=("tumor_vaf", "median"),
            median_tumor_vaf_adj=("tumor_vaf_adj", "median"),
            median_vaf_ratio_adj=("tumor_vaf_ratio_adj", "median"),
        )
    )
    sv_intervals = (
        hits.groupby(["participant", "locus_id"], group_keys=False)
        .apply(_sv_intervals, include_groups=False)
        .reset_index(name="sv_intervals")
    )
    sv_agg = sv_agg.merge(sv_intervals, on=["participant", "locus_id"], how="left")

    core = queue[["participant", "locus_id"]].drop_duplicates()
    for col in [
        "cluster_id",
        "cn_state",
        "expected_svtype",
        "n_sv_events",
        "has_SV_DUP",
        "has_SV_DEL",
        "any_concordant_sv",
        "discordant_sv_only",
        "mixed_dup_and_del",
        "focal_sv_le_1mb",
        "median_sv_span_mb",
        "median_tumor_alt",
    ]:
        if col in queue.columns:
            core = core.merge(queue[["participant", "locus_id", col]].drop_duplicates(), on=["participant", "locus_id"], how="left")

    out = core.merge(meta, on="locus_id", how="left")
    out = out.merge(lay[lay_cols], on=["participant", "locus_id"], how="left", suffixes=("", "_lay"))
    if "PAM50_final_lay" in out.columns:
        out["PAM50_final"] = out["PAM50_final"].fillna(out["PAM50_final_lay"])
        out = out.drop(columns=["PAM50_final_lay"])
    out = out.merge(sv_agg, on=["participant", "locus_id"], how="left", suffixes=("", "_hit"))
    for col in ("median_sv_span_mb", "median_tumor_alt", "sample_vial"):
        hit_col = f"{col}_hit"
        if hit_col in out.columns:
            out[col] = out[col].fillna(out[hit_col]) if col in out.columns else out[hit_col]
            out = out.drop(columns=[hit_col])
    for col in ("median_tumor_vaf", "median_tumor_vaf_adj", "median_vaf_ratio_adj"):
        hit_col = f"{col}_hit"
        if hit_col in out.columns:
            out[col] = out[col].fillna(out[hit_col]) if col in out.columns else out[hit_col]
            out = out.drop(columns=[hit_col])

    out.insert(0, "review_queue", review_queue)
    pad_mb = pd.to_numeric(out["median_sv_span_mb"], errors="coerce").fillna(0.5).clip(lower=0.5)
    out["igv_locus"] = out.apply(
        lambda r: _sv_interval_label(r["chrom"], r["start"], r["end"]),
        axis=1,
    )
    out["igv_window"] = [
        _igv_window(r.chrom, r.start, r.end, p)
        for r, p in zip(out.itertuples(index=False), pad_mb)
    ]
    out["gdc_case_query"] = out["participant"].astype(str)

    preferred = [
        "review_queue",
        "participant",
        "sample_vial",
        "PAM50_final",
        "locus_id",
        "locus_name",
        "chrom",
        "start",
        "end",
        "igv_locus",
        "igv_window",
        "copy_number",
        "cn_state",
        "overlap_segment",
        "dosage_layer",
        "delta_vs_diploid",
        "sv_ids",
        "svtypes",
        "sv_intervals",
        "median_sv_span_mb",
        "median_tumor_alt",
        "median_tumor_vaf",
        "median_tumor_vaf_adj",
        "median_vaf_ratio_adj",
        "CPE",
        "n_sv_events",
        "has_SV_DUP",
        "has_SV_DEL",
        "any_concordant_sv",
        "discordant_sv_only",
        "mixed_dup_and_del",
        "focal_sv_le_1mb",
        "cluster_id",
        "gdc_case_query",
    ]
    cols = [c for c in preferred if c in out.columns] + [c for c in out.columns if c not in preferred]
    out = out[cols]
    sort_cols = ["review_queue", "PAM50_final", "participant", "locus_id"]
    sort_cols = [c for c in sort_cols if c in out.columns]
    return out.sort_values(sort_cols)


def build_review_queues(
    both_pl: pd.DataFrame,
    svo_pl: pd.DataFrame,
    layers: pd.DataFrame,
    hits: pd.DataFrame,
    locus_map: pd.DataFrame,
    out_dir: Path,
) -> Dict[str, Path]:
    """Write three IGV-oriented review queues."""
    queue_dir = out_dir / "review_queues"
    queue_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}

    specs = [
        ("sv_only_focal_le_1mb", svo_pl.loc[svo_pl["focal_sv_le_1mb"]].copy()),
        ("both_discordant_sv_only", both_pl.loc[both_pl["discordant_sv_only"]].copy()),
        ("both_mixed_dup_del", both_pl.loc[both_pl["mixed_dup_and_del"]].copy()),
    ]
    summary_rows: List[dict] = []
    for name, sub in specs:
        enriched = _enrich_review_rows(sub, locus_map=locus_map, layers=layers, hits=hits, review_queue=name)
        path = queue_dir / f"review_queue_{name}.tsv.gz"
        enriched.to_csv(path, sep="\t", index=False, compression="gzip")
        paths[name] = path
        summary_rows.append(
            {
                "review_queue": name,
                "n_rows": int(len(enriched)),
                "n_participants": int(enriched["participant"].nunique()) if not enriched.empty else 0,
                "n_loci": int(enriched["locus_id"].nunique()) if not enriched.empty else 0,
            }
        )

    summary_path = queue_dir / "review_queue_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    paths["summary"] = summary_path
    paths.update(export_igv_review_bundle(queue_dir, locus_map, hits))
    return paths


def _bed_line(chrom: str, start1: int, end1: int, name: str, rgb: str = "0,102,204") -> str:
    bed0 = max(0, int(start1) - 1)
    return f"{chrom}\t{bed0}\t{int(end1)}\t{name}\t0\t+\t{bed0}\t{int(end1)}\t{rgb}\n"


def export_igv_review_bundle(
    queue_dir: Path,
    locus_map: pd.DataFrame,
    hits: pd.DataFrame,
) -> Dict[str, Path]:
    """Write IGV helper BED + session guide + README (review is one sample at a time)."""
    igv_dir = queue_dir / "igv"
    igv_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}

    bed_path = igv_dir / "mirna_hairpin_loci.bed"
    bed_lines = []
    for row in locus_map.itertuples(index=False):
        bed_lines.append(
            _bed_line(
                str(row.chrom),
                int(row.start),
                int(row.end),
                f"{row.locus_id}|{row.locus_name}",
            )
        )
    bed_path.write_text("".join(bed_lines), encoding="utf-8")
    paths["igv_mirna_bed"] = bed_path

    session_rows: List[pd.DataFrame] = []
    for path in sorted(queue_dir.glob("review_queue_*.tsv.gz")):
        session_rows.append(pd.read_csv(path, sep="\t", low_memory=False))
    session = pd.concat(session_rows, ignore_index=True) if session_rows else pd.DataFrame()
    if not session.empty and "igv_window" in session.columns:
        session["ucsc_url"] = session["igv_window"].astype(str).map(
            lambda w: f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={w}" if w and w != "nan" else ""
        )
    session_path = igv_dir / "igv_session_guide.tsv"
    if not session.empty:
        session.to_csv(session_path, sep="\t", index=False)
    else:
        session_path.write_text("review_queue\tparticipant\tlocus_id\n", encoding="utf-8")
    paths["igv_session_guide"] = session_path

    sv_bed_path = igv_dir / "review_queue_sv_intervals.bed"
    if not session.empty and not hits.empty:
        sub_hits = hits.merge(
            session[["participant", "locus_id"]].drop_duplicates(),
            on=["participant", "locus_id"],
            how="inner",
        )
        sv_lines = []
        for row in sub_hits.drop_duplicates(["participant", "sv_id"]).itertuples(index=False):
            rgb = "220,50,50" if str(row.svtype) == "DEL" else "50,120,220"
            name = f"{row.participant}|{row.svtype}|{row.sv_id}|{row.locus_id}"
            end = row.sv_end if pd.notna(row.sv_end) else row.sv_start
            sv_lines.append(_bed_line(str(row.sv_chrom), int(row.sv_start), int(end), name, rgb=rgb))
        sv_bed_path.write_text("".join(sv_lines), encoding="utf-8")
        paths["igv_sv_bed"] = sv_bed_path

    readme = """Review queues — miRNA locus SV/CNV
====================================

Each row is **one participant × one hairpin locus**.

**No BAM/CRAM? Start here instead**
----------------------------------
1. **Static schematics (recommended)** — ``../case_plots/index.html``
   Gray = ASCAT3 segment; red = hairpin; blue/red = DUP/DEL SV. Generated by:
   ``python -m mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_case_plots``
2. **UCSC Genome Browser** — open ``ucsc_url`` in ``igv_session_guide.tsv`` (hg38,
   no login). Optionally add custom tracks: ``mirna_hairpin_loci.bed``,
   ``review_queue_sv_intervals.bed`` (Track → Custom track → paste BED).
3. **Tables** — ``igv_session_guide.tsv`` columns ``overlap_segment``, ``sv_intervals``,
   ``cn_state``, ``dosage_layer`` describe the call without reads.

**IGV with BAM/CRAM (optional)**
--------------------------------
Workflow (repeat for each row you want to inspect):

1. **Pick a row** from ``igv_session_guide.tsv`` (sort by ``review_queue``, ``PAM50_final``,
   or ``median_sv_span_mb``). Note ``participant``, ``sample_vial``, ``igv_window``.

2. **Load the tumor BAM/CRAM** for that participant in IGV:
   - GDC Data Portal → search case ID = ``participant`` (e.g. TCGA-3C-AALI).
   - Open the **primary tumor** aliquot; ``sample_vial`` in the table is the
     barcode suffix (…-01A) when present in the SV layer.
   - Load aligned reads (BAM/CRAM) for that aliquot.

3. **Load copy-number context** (optional but recommended):
   - ASCAT3 allelic segments or SNP-array copy ratio for the same sample, if available locally.
   - The table column ``overlap_segment`` is the ASCAT3 segment covering this hairpin.

4. **Load SV evidence**:
   - ``review_queue_sv_intervals.bed`` — DEL (red) / DUP (blue) intervals for all queued cases.
     Filter mentally to the current ``participant`` (name prefix) or subset the BED.
   - Alternatively load the per-sample ``*_strict_sv_set.csv`` from the SV pipeline
     (``data/SV/pipeline_output/07_final_sv_with_fimo/``) and search ``sv_ids`` from the row.

5. **Load miRNA annotation**:
   - ``mirna_hairpin_loci.bed`` — all 506 MirGeneDB hairpins (blue).

6. **Navigate**: paste ``igv_window`` into IGV (e.g. ``chr22:45612752-46613765``). ``igv_locus`` is the
   tight hairpin interval.

7. **What to check**:
   - ``sv_only_focal_le_1mb``: CN=2 segment but focal DEL/DUP overlapping hairpin — subclonal SV?
   - ``both_discordant_sv_only``: CNV amp/loss but SV direction mismatches (e.g. DEL on amp).
   - ``both_mixed_dup_del``: both DEL and DUP calls at the same locus on an altered segment.

Cohort-level bedGraphs (all samples) live under ``analysis/output/sv_recurrence/`` via
``python -m analysis.variants.sv.export_igv_tracks`` — use those for recurrence context,
then drill into individual participants from this guide.
"""
    readme_path = igv_dir / "README_IGV.txt"
    readme_path.write_text(readme, encoding="utf-8")
    paths["igv_readme"] = readme_path
    return paths


def build_chr19_megacluster_dossier(
    layers: pd.DataFrame,
    hits: pd.DataFrame,
    locus_map: pd.DataFrame,
    clusters: pd.DataFrame,
    both_pl: pd.DataFrame,
    svo_pl: pd.DataFrame,
    out_dir: Path,
    *,
    cluster_id: str = CHR19_MEGACLUSTER_ID,
) -> Dict[str, Path]:
    """chr19:53.7 megacluster — locus table, participant layers, summary."""
    dossier_dir = out_dir / "chr19_megacluster"
    dossier_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {}

    cl_row = clusters.loc[
        clusters.apply(
            lambda r: f"{r['chrom']}:{r['cluster_start_mb']:.3f}-{r['cluster_end_mb']:.3f}",
            axis=1,
        )
        == cluster_id
    ]
    if cl_row.empty:
        return paths
    cl_row = cl_row.iloc[0]
    block = _loci_in_cluster_block(locus_map, cl_row)
    loci = block["locus_id"].astype(str).tolist()

    lay = layers.copy()
    lay["locus_id"] = lay["entity_id"].astype(str)
    cl_layers = lay.loc[lay["locus_id"].isin(loci)].copy()

    locus_rows: List[dict] = []
    for lid in loci:
        sub = cl_layers.loc[cl_layers["locus_id"] == lid]
        meta = block.loc[block["locus_id"].astype(str) == lid].iloc[0]
        locus_rows.append(
            {
                "locus_id": lid,
                "locus_name": meta["locus_name"],
                "start": int(meta["start"]),
                "end": int(meta["end"]),
                "mid_mb": round(float(meta["midpoint"]) / 1e6, 4),
                "cohort_mean_delta": meta["delta_vs_diploid"],
                "pct_cnv_altered": round(100 * sub["cnv_altered"].mean(), 2) if len(sub) else np.nan,
                "pct_sv_only": round(100 * (sub["dosage_layer"] == "sv_only").mean(), 2) if len(sub) else np.nan,
                "pct_both": round(100 * (sub["dosage_layer"] == "both").mean(), 2) if len(sub) else np.nan,
                "n_sv_participants": int(
                    hits.loc[hits["locus_id"].astype(str) == lid, "participant"].nunique()
                )
                if not hits.empty
                else 0,
            }
        )
    locus_dossier = pd.DataFrame(locus_rows).sort_values("mid_mb")
    paths["loci"] = dossier_dir / "chr19_megacluster_loci.tsv"
    locus_dossier.to_csv(paths["loci"], sep="\t", index=False)

    part = cl_layers.merge(
        both_pl[["participant", "locus_id", "any_concordant_sv", "discordant_sv_only", "mixed_dup_and_del"]],
        on=["participant", "locus_id"],
        how="left",
    ).merge(
        svo_pl[["participant", "locus_id", "focal_sv_le_1mb", "median_sv_span_mb"]],
        on=["participant", "locus_id"],
        how="left",
    )
    paths["participant_locus"] = dossier_dir / "chr19_megacluster_participant_locus.tsv.gz"
    part.to_csv(paths["participant_locus"], sep="\t", index=False, compression="gzip")

    layer_counts = cl_layers["dosage_layer"].value_counts().to_dict()
    summary = {
        "cluster_id": cluster_id,
        "cluster_span_kb": float(cl_row["cluster_span_kb"]),
        "n_hairpin_loci": int(len(loci)),
        "n_participant_locus_rows": int(len(cl_layers)),
        "layer_counts": layer_counts,
        "n_participants_sv_any": int(hits.loc[hits["locus_id"].isin(loci), "participant"].nunique())
        if not hits.empty
        else 0,
        "cohort_mean_locus_delta": round(float(block["delta_vs_diploid"].mean()), 4),
    }
    paths["summary_json"] = dossier_dir / "chr19_megacluster_summary.json"
    paths["summary_json"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return paths


def write_review_queues_from_cache(out_dir: Optional[Path] = None) -> Dict[str, Path]:
    """Rebuild review queues + chr19 dossier from cached SV overlap tables (no SV rescan)."""
    out_dir = Path(out_dir or C.MIRNA_LOCUS_SV_DIR)
    locus_map = pd.read_csv(C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_genome_map.tsv", sep="\t")
    clusters = build_dense_locus_clusters(locus_map)
    layers = pd.read_csv(out_dir / "mirna_locus_cnv_sv_layers.tsv.gz", sep="\t", low_memory=False)
    hits_path = out_dir / "mirna_locus_sv_hits.tsv.gz"
    hits = pd.read_csv(hits_path, sep="\t", low_memory=False)
    hits = enrich_sv_hits_with_vaf(
        hits,
        layers_path=C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz",
        sv_dir=C.SV_FINAL_DIR,
        save_path=hits_path,
    )
    both_pl = pd.read_csv(out_dir / "mirna_locus_sv_both_participant_locus.tsv.gz", sep="\t")
    svo_pl = pd.read_csv(out_dir / "mirna_locus_sv_only_participant_locus.tsv.gz", sep="\t")

    paths: Dict[str, Path] = {}
    paths.update(build_review_queues(both_pl, svo_pl, layers, hits, locus_map, out_dir))
    paths.update(
        build_chr19_megacluster_dossier(
            layers, hits, locus_map, clusters, both_pl, svo_pl, out_dir
        )
    )
    from mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_case_plots import write_case_plot_gallery

    paths.update(write_case_plot_gallery(queue_dir=out_dir / "review_queues", hits_path=out_dir / "mirna_locus_sv_hits.tsv.gz"))
    from mirna_hallmark.analyses.cnv_locus.mirna_locus_chr19_megacluster_gallery import write_chr19_megacluster_gallery

    paths.update(write_chr19_megacluster_gallery(sv_dir=out_dir))
    print(f"[mirna_locus_sv_overlap] review queues + chr19 dossier -> {out_dir}")
    return paths


def scan_sv_mirna_overlaps(
    *,
    sv_dir: Path,
    locus_map: pd.DataFrame,
    clusters: pd.DataFrame,
    workers: int = 8,
    max_samples: Optional[int] = None,
) -> pd.DataFrame:
    files = sorted(sv_dir.glob("*_strict_sv_set.csv"))
    if max_samples:
        files = files[: max_samples]
    if not files:
        raise FileNotFoundError(f"No SV strict-set CSVs under {sv_dir}")

    cluster_map = _assign_cluster_ids(locus_map, clusters)
    parts: List[pd.DataFrame] = []
    n = len(files)
    if workers <= 1:
        for i, f in enumerate(files):
            rows = _scan_one_sample(str(f))
            if rows:
                parts.append(pd.DataFrame(rows))
            if (i + 1) % 100 == 0:
                print(f"  SV×miRNA: {i + 1}/{n} samples…")
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_scan_one_sample, str(f)): f for f in files}
            for i, fut in enumerate(as_completed(futs), start=1):
                rows = fut.result()
                if rows:
                    parts.append(pd.DataFrame(rows))
                if i % 100 == 0:
                    print(f"  SV×miRNA: {i}/{n} samples…")

    if not parts:
        return pd.DataFrame()
    hits = pd.concat(parts, ignore_index=True)
    hits["cluster_id"] = hits["locus_id"].astype(str).map(cluster_map).fillna("")
    return hits


def write_sv_overlap_tables(
    *,
    sv_dir: Optional[Path] = None,
    cnv_long_path: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    workers: int = 8,
    max_samples: Optional[int] = None,
) -> Dict[str, Path]:
    sv_dir = Path(sv_dir or C.SV_FINAL_DIR)
    out_dir = Path(out_dir or C.MIRNA_LOCUS_SV_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    cnv_path = Path(cnv_long_path or C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz")

    locus_map = pd.read_csv(C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_genome_map.tsv", sep="\t")
    clusters = build_dense_locus_clusters(locus_map)

    print(f"[mirna_locus_sv_overlap] scanning {sv_dir} …")
    hits = scan_sv_mirna_overlaps(
        sv_dir=sv_dir,
        locus_map=locus_map,
        clusters=clusters,
        workers=workers,
        max_samples=max_samples,
    )

    paths: Dict[str, Path] = {}
    hits_path = out_dir / "mirna_locus_sv_hits.tsv.gz"
    hits = enrich_sv_hits_with_vaf(
        hits,
        layers_path=cnv_path if cnv_path.exists() else None,
        sv_dir=sv_dir,
        save_path=hits_path,
    )
    paths["hits"] = hits_path

    locus_rec = build_locus_recurrence(hits)
    cluster_rec = build_cluster_recurrence(hits)
    paths["locus_recurrence"] = out_dir / "mirna_locus_sv_recurrence.tsv"
    paths["cluster_recurrence"] = out_dir / "mirna_cluster_sv_recurrence.tsv"
    locus_rec.to_csv(paths["locus_recurrence"], sep="\t", index=False)
    cluster_rec.to_csv(paths["cluster_recurrence"], sep="\t", index=False)

    layers = pd.DataFrame()
    both_pl = pd.DataFrame()
    svo_pl = pd.DataFrame()
    if cnv_path.exists():
        cnv_long = pd.read_csv(cnv_path, sep="\t", low_memory=False)
        cluster_map = _assign_cluster_ids(locus_map, clusters)
        layers = build_cnv_sv_layers(hits, cnv_long, cluster_map)
        layers_path = out_dir / "mirna_locus_cnv_sv_layers.tsv.gz"
        layers.to_csv(layers_path, sep="\t", index=False, compression="gzip")
        paths["layers"] = layers_path

        both_pl = _both_participant_locus_table(layers, hits)
        svo_pl = _sv_only_participant_locus_table(layers, hits)
        both_summary, svo_summary, svo_top, both_cn = build_layer_pattern_summaries(
            layers, hits, locus_map
        )
        paths["both_patterns"] = out_dir / "mirna_locus_sv_both_patterns.tsv"
        paths["sv_only_patterns"] = out_dir / "mirna_locus_sv_only_patterns.tsv"
        paths["sv_only_top_loci"] = out_dir / "mirna_locus_sv_only_top_loci.tsv"
        paths["both_concordance_by_cn_state"] = out_dir / "mirna_locus_sv_both_concordance_by_cn_state.tsv"
        both_summary.to_csv(paths["both_patterns"], sep="\t", index=False)
        svo_summary.to_csv(paths["sv_only_patterns"], sep="\t", index=False)
        svo_top.to_csv(paths["sv_only_top_loci"], sep="\t", index=False)
        both_cn.to_csv(paths["both_concordance_by_cn_state"], sep="\t", index=False)
        if not both_pl.empty:
            paths["both_participant_locus"] = out_dir / "mirna_locus_sv_both_participant_locus.tsv.gz"
            both_pl.to_csv(paths["both_participant_locus"], sep="\t", index=False, compression="gzip")
        if not svo_pl.empty:
            paths["sv_only_participant_locus"] = out_dir / "mirna_locus_sv_only_participant_locus.tsv.gz"
            svo_pl.to_csv(paths["sv_only_participant_locus"], sep="\t", index=False, compression="gzip")

        paths.update(build_review_queues(both_pl, svo_pl, layers, hits, locus_map, out_dir))
        paths.update(
            build_chr19_megacluster_dossier(
                layers, hits, locus_map, clusters, both_pl, svo_pl, out_dir
            )
        )

    pattern_notes: Dict[str, object] = {}
    if not layers.empty and not hits.empty:
        if not both_pl.empty:
            pattern_notes["both_pct_any_concordant_sv"] = round(100 * both_pl["any_concordant_sv"].mean(), 2)
            pattern_notes["both_pct_mixed_dup_and_del"] = round(100 * both_pl["mixed_dup_and_del"].mean(), 2)
        if not svo_pl.empty:
            pattern_notes["sv_only_pct_focal_le_1mb"] = round(100 * svo_pl["focal_sv_le_1mb"].mean(), 2)
            pattern_notes["sv_only_median_sv_span_mb"] = round(float(svo_pl["median_sv_span_mb"].median()), 2)

    manifest = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sv_dir": str(sv_dir),
        "cnv_long": str(cnv_path),
        "n_sv_files_scanned": len(list(sv_dir.glob("*_strict_sv_set.csv"))) if not max_samples else max_samples,
        "n_hit_rows": int(len(hits)),
        "n_participants_with_hits": int(hits["participant"].nunique()) if not hits.empty else 0,
        "n_loci_with_hits": int(hits["locus_id"].nunique()) if not hits.empty else 0,
        "overlap_rule": "mir_hits region_hit=overlaps OR overlap_bp>0; SVTYPE in DEL,DUP",
        "dosage_layers": {
            "cnv_only": "ASCAT3 cn_state altered, no overlapping SV DEL/DUP at locus",
            "sv_only": "overlapping SV DEL/DUP, ASCAT3 not altered",
            "both": "CNV altered and SV DEL/DUP overlap",
            "neither": "diploid CNV and no SV overlap",
        },
        "layer_row_counts": (
            layers["dosage_layer"].value_counts().to_dict() if not layers.empty else {}
        ),
        "layer_pattern_notes": pattern_notes,
    }
    manifest_path = out_dir / "method_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    paths["manifest"] = manifest_path

    print(
        f"[mirna_locus_sv_overlap] hits: {len(hits):,}; loci: {manifest['n_loci_with_hits']}; "
        f"participants: {manifest['n_participants_with_hits']}"
    )
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sv-dir", type=Path, default=C.SV_FINAL_DIR)
    ap.add_argument("--cnv-long", type=Path, default=C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz")
    ap.add_argument("--out-dir", type=Path, default=C.MIRNA_LOCUS_SV_DIR)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-samples", type=int, default=None)
    ap.add_argument(
        "--queues-only",
        action="store_true",
        help="Rebuild review queues + chr19 dossier + case schematics from cached tables (skip SV scan)",
    )
    ap.add_argument(
        "--case-plots-only",
        action="store_true",
        help="Regenerate static case schematic gallery only (no BAM required)",
    )
    ap.add_argument(
        "--chr19-gallery-only",
        action="store_true",
        help="Regenerate chr19 megacluster participant HTML gallery only",
    )
    ap.add_argument(
        "--enrich-vaf-only",
        action="store_true",
        help="Refresh tumor_vaf / CCF / ratio columns on cached hits (no full SV rescan)",
    )
    args = ap.parse_args()
    if args.enrich_vaf_only:
        write_review_queues_from_cache(out_dir=args.out_dir)
        return
    if args.chr19_gallery_only:
        from mirna_hallmark.analyses.cnv_locus.mirna_locus_chr19_megacluster_gallery import write_chr19_megacluster_gallery

        write_chr19_megacluster_gallery(sv_dir=args.out_dir)
        return
    if args.case_plots_only:
        from mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_case_plots import write_case_plot_gallery

        write_case_plot_gallery(
            queue_dir=args.out_dir / "review_queues",
            hits_path=args.out_dir / "mirna_locus_sv_hits.tsv.gz",
        )
        return
    if args.queues_only:
        write_review_queues_from_cache(out_dir=args.out_dir)
        return
    write_sv_overlap_tables(
        sv_dir=args.sv_dir,
        cnv_long_path=args.cnv_long,
        out_dir=args.out_dir,
        workers=args.workers,
        max_samples=args.max_samples,
    )


if __name__ == "__main__":
    main()
