"""Genome-wide miRNA CNV maps: locus-level and MIMAT paralog decomposition.

Healthy reference:
  - **Locus:** diploid CN = 2  →  Δ = observed − 2
  - **MIMAT (multi-paralog):** cohort expression weights ``w_norm`` per (locus, mimat);
    weighted healthy CN = Σ (2 × w_norm) = 2; per-locus contribution = 2 × w_norm;
    weighted observed CN = Σ (w_norm × CN_locus); Δ = weighted − 2

Outputs under ``output/mirna_locus_cnv/maps/``:
  - ``mirna_cnv_locus_genome_map.tsv`` (+ per-PAM50 columns)
  - ``mirna_cnv_mimat_paralog_map.tsv`` (long: mimat × locus with weights)
  - ``mirna_cnv_mimat_summary.tsv`` (one row per MIMAT)
  - ``mirna_cnv_chromosome_summary.tsv`` / ``mirna_cnv_chromosome_arm_summary.tsv``
    / ``mirna_cnv_locus_cluster_summary.tsv`` (landscape counts + sizes)
  - ``mirna_cnv_cluster_segment_focality.tsv`` (dense-cluster vs ASCAT3 segment span)
  - ``mirna_cnv_interesting_cases.tsv``  -- flagged loci/MIMATs with reasons
  - ``mirna_cnv_locus_genome_map.png`` / ``mirna_cnv_mimat_paralog_genome_map.png``
  - ``mirna_cnv_locus_genome_map_by_pam50.png`` / ``mirna_cnv_mimat_paralog_genome_map_by_pam50.png``
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

DIPLOID_CN = 2.0
PAM50_ORDER = ["Normal", "LumA", "LumB", "Basal", "Her2"]
TUMOR_PAM50 = ["LumA", "LumB", "Basal", "Her2"]
CHR_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]

# hg38 chromosome and cytoband-arm sizes (Mb) for landscape summaries.
CHR_SIZE_MB: Dict[str, float] = {
    "chr1": 248.96,
    "chr2": 242.19,
    "chr3": 198.22,
    "chr4": 190.21,
    "chr5": 181.54,
    "chr6": 170.80,
    "chr7": 159.34,
    "chr8": 145.08,
    "chr9": 138.39,
    "chr10": 133.79,
    "chr11": 135.09,
    "chr12": 133.27,
    "chr13": 114.36,
    "chr14": 107.04,
    "chr15": 101.99,
    "chr16": 90.34,
    "chr17": 83.19,
    "chr18": 80.37,
    "chr19": 58.62,
    "chr20": 64.44,
    "chr21": 46.71,
    "chr22": 50.82,
    "chrX": 156.04,
    "chrY": 57.22,
}
CENTROMERE_MB: Dict[str, float] = {
    "chr1": 121.5,
    "chr2": 92.4,
    "chr3": 91.0,
    "chr4": 50.4,
    "chr5": 48.4,
    "chr6": 60.0,
    "chr7": 60.2,
    "chr8": 45.2,
    "chr9": 49.0,
    "chr10": 40.2,
    "chr11": 53.7,
    "chr12": 36.0,
    "chr13": 17.7,
    "chr14": 17.2,
    "chr15": 19.0,
    "chr16": 36.4,
    "chr17": 25.2,
    "chr18": 18.5,
    "chr19": 26.2,
    "chr20": 26.4,
    "chr21": 10.9,
    "chr22": 14.5,
    "chrX": 58.1,
    "chrY": 10.2,
}

# Curated arms for explicit biology flags (oncogenic / amplicon-neighbor context).
CHR19_MEGACLUSTER_ID = "chr19:53.667-53.762"

CURATED_ARMS = {
    "hsa-miR-21-5p": "chr17/MYC-pathway neighbor; LumB–Her2 amp",
    "hsa-miR-21-3p": "chr17/MYC-pathway neighbor; LumB–Her2 amp",
    "hsa-miR-151a-3p": "chr8 amplicon (paralog MI0000809)",
    "hsa-miR-151a-5p": "chr8 amplicon (paralog MI0000809)",
    "hsa-miR-17-5p": "miR-17~92 cluster",
    "hsa-miR-200c-3p": "miR-200 family / EMT axis",
    "hsa-miR-205-5p": "chr1 amp; luminal-enriched gain",
    "hsa-miR-135b-5p": "chr1 amp hotspot",
    "hsa-let-7a-5p": "multi-paralog let-7 (chr9/11/22)",
}


def _load_locus_catalog() -> pd.DataFrame:
    df = pd.read_csv(C.MIRNA_PRECURSOR_LOCI)
    df = df.rename(columns={"gene_id": "locus_id", "gene_name": "locus_name"})
    df["chrom"] = df["chrom"].astype(str).str.strip()
    df.loc[~df["chrom"].str.startswith("chr"), "chrom"] = "chr" + df["chrom"]
    df["midpoint"] = (pd.to_numeric(df["start"], errors="coerce") + pd.to_numeric(df["end"], errors="coerce")) / 2.0
    return df


def _load_weight_reference() -> pd.DataFrame:
    from pipeline.miRNA_exp.locus_mimat_expression import load_locus_mimat_weight_reference

    w = load_locus_mimat_weight_reference()
    w["weight_ref"] = pd.to_numeric(w["weight_ref"], errors="coerce").fillna(0.0)
    totals = w.groupby("mimat")["weight_ref"].transform("sum")
    w["w_norm"] = np.where(totals > 0, w["weight_ref"] / totals, 0.0)
    w["healthy_cn_contribution"] = DIPLOID_CN * w["w_norm"]
    return w


def _load_arm_labels() -> pd.DataFrame:
    arms = pd.read_csv(C.MIRNA_MATURE_LOCI)
    id_col = "mature_accession" if "mature_accession" in arms.columns else "mimat"
    label_col = "mirbase_mature_id" if "mirbase_mature_id" in arms.columns else "gene_name"
    out = arms[[id_col, label_col]].drop_duplicates(subset=[id_col], keep="first")
    return out.rename(columns={id_col: "mimat", label_col: "arm_label"})


def _mean_locus_cn(
    long: pd.DataFrame,
    *,
    stratum_col: Optional[str] = None,
    stratum_val: Optional[str] = None,
) -> pd.Series:
    sub = long.loc[long["entity_level"] == "locus"].copy()
    if stratum_col and stratum_val is not None:
        sub = sub.loc[sub[stratum_col].astype(str) == str(stratum_val)]
    sub = sub.drop_duplicates(subset=["participant", "entity_id"], keep="first")
    sub["copy_number"] = pd.to_numeric(sub["copy_number"], errors="coerce")
    return sub.groupby("entity_id")["copy_number"].mean()


def _genome_x(chrom: str, pos: float, offsets: Dict[str, float]) -> float:
    if chrom not in offsets or pd.isna(pos):
        return np.nan
    return offsets[chrom] + float(pos)


def _cytoband_arm(chrom: str, midpoint: float) -> str:
    """Return cytoband arm label (e.g. chr8q) from hg38 midpoint."""
    cen = CENTROMERE_MB.get(str(chrom))
    if cen is None or pd.isna(midpoint):
        return ""
    pos_mb = float(midpoint) / 1e6
    return f"{chrom}{'p' if pos_mb < cen else 'q'}"


def _arm_size_mb(cytoband_arm: str) -> float:
    chrom = cytoband_arm[:-1] if cytoband_arm.endswith(("p", "q")) else cytoband_arm
    if chrom not in CHR_SIZE_MB or chrom not in CENTROMERE_MB:
        return np.nan
    cen = CENTROMERE_MB[chrom]
    total = CHR_SIZE_MB[chrom]
    if cytoband_arm.endswith("p"):
        return cen
    if cytoband_arm.endswith("q"):
        return total - cen
    return total


def build_chromosome_cn_summary(
    locus_map: pd.DataFrame,
    paralog_map: pd.DataFrame,
) -> pd.DataFrame:
    """One row per chromosome: locus counts, mature-arm counts, mean CN delta."""
    mimat_by_chrom = paralog_map.groupby("chrom")["mimat"].nunique()
    rows: List[dict] = []
    for chrom in CHR_ORDER:
        sub = locus_map.loc[locus_map["chrom"] == chrom]
        if sub.empty:
            continue
        delta = pd.to_numeric(sub["delta_vs_diploid"], errors="coerce")
        rows.append(
            {
                "chrom": chrom,
                "chr_size_mb": round(CHR_SIZE_MB.get(chrom, np.nan), 2),
                "n_hairpin_loci": int(len(sub)),
                "n_mature_arms": int(mimat_by_chrom.get(chrom, 0)),
                "mean_locus_delta_vs_diploid": round(float(delta.mean()), 4),
                "n_loci_delta_gt_1": int((delta > 1.0).sum()),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values("mean_locus_delta_vs_diploid", ascending=False, na_position="last")


def build_chromosome_arm_cn_summary(
    locus_map: pd.DataFrame,
    paralog_map: pd.DataFrame,
) -> pd.DataFrame:
    """One row per cytoband arm (chrNp/chrNq): counts, arm size, mean CN delta."""
    loc = locus_map.copy()
    loc["cytoband_arm"] = [
        _cytoband_arm(r.chrom, r.midpoint) for r in loc.itertuples(index=False)
    ]
    para = paralog_map.copy()
    para["cytoband_arm"] = [
        _cytoband_arm(r.chrom, r.midpoint) for r in para.itertuples(index=False)
    ]
    mimat_by_arm = para.groupby("cytoband_arm")["mimat"].nunique()
    rows: List[dict] = []
    for arm in sorted(set(loc["cytoband_arm"]) - {""}):
        sub = loc.loc[loc["cytoband_arm"] == arm]
        delta = pd.to_numeric(sub["delta_vs_diploid"], errors="coerce")
        chrom = arm[:-1] if arm.endswith(("p", "q")) else arm
        rows.append(
            {
                "cytoband_arm": arm,
                "chrom": chrom,
                "arm_size_mb": round(_arm_size_mb(arm), 2),
                "n_hairpin_loci": int(len(sub)),
                "n_mature_arms": int(mimat_by_arm.get(arm, 0)),
                "mean_locus_delta_vs_diploid": round(float(delta.mean()), 4),
                "n_loci_delta_gt_1": int((delta > 1.0).sum()),
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values("mean_locus_delta_vs_diploid", ascending=False, na_position="last")


def build_dense_locus_clusters(
    locus_map: pd.DataFrame,
    *,
    window_bp: int = 100_000,
    min_loci: int = 3,
) -> pd.DataFrame:
    """Hairpin clusters within ``window_bp`` on the same chromosome."""
    rows: List[dict] = []
    for chrom, grp in locus_map.sort_values(["chrom", "start"]).groupby("chrom"):
        g = grp.sort_values("start").reset_index(drop=True)
        i = 0
        while i < len(g):
            j = i
            while j + 1 < len(g) and g.loc[j + 1, "start"] - g.loc[i, "start"] <= window_bp:
                j += 1
            if j - i + 1 >= min_loci:
                block = g.iloc[i : j + 1]
                delta = pd.to_numeric(block["delta_vs_diploid"], errors="coerce")
                rows.append(
                    {
                        "chrom": chrom,
                        "cluster_start_mb": round(block["start"].min() / 1e6, 3),
                        "cluster_end_mb": round(block["end"].max() / 1e6, 3),
                        "cluster_span_kb": round((block["end"].max() - block["start"].min()) / 1e3, 2),
                        "n_hairpin_loci": int(len(block)),
                        "n_mature_arms": int(
                            block["mature_accessions"]
                            .astype(str)
                            .str.split(",")
                            .explode()
                            .nunique()
                        ),
                        "mean_locus_delta_vs_diploid": round(float(delta.mean()), 4),
                        "example_loci": "; ".join(block["locus_name"].head(4).astype(str).tolist()),
                    }
                )
            i = j + 1
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("mean_locus_delta_vs_diploid", ascending=False, na_position="last")


def _parse_overlap_segment(seg_label: str) -> Optional[Tuple[str, int, int]]:
    if not seg_label or pd.isna(seg_label) or ":" not in str(seg_label):
        return None
    try:
        chrom, coords = str(seg_label).split(":", 1)
        start_s, end_s = coords.split("-", 1)
        return chrom, int(start_s), int(end_s)
    except (ValueError, TypeError):
        return None


def _segment_length_mb(seg_label: str) -> float:
    parsed = _parse_overlap_segment(seg_label)
    if parsed is None:
        return np.nan
    _, start, end = parsed
    return (end - start) / 1e6


def _loci_in_cluster_block(locus_map: pd.DataFrame, cluster_row: pd.Series) -> pd.DataFrame:
    chrom = cluster_row["chrom"]
    start_mb = float(cluster_row["cluster_start_mb"])
    end_mb = float(cluster_row["cluster_end_mb"])
    return locus_map.loc[
        (locus_map["chrom"] == chrom)
        & (locus_map["start"] / 1e6 >= start_mb - 1e-3)
        & (locus_map["end"] / 1e6 <= end_mb + 1e-3)
    ].copy()


def build_cluster_segment_focality(
    long: pd.DataFrame,
    locus_map: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    *,
    focal_overhang_mb: float = 5.0,
    focal_segment_mb: float = 10.0,
    focal_cluster_multiplier: float = 10.0,
) -> pd.DataFrame:
    """Per dense cluster: are ASCAT3 amp/del segments focal around the cluster window?

    Uses participant-level ``overlap_segment`` labels from segment-based miRNA CNV.
    A cluster is *segment-coherent* when all hairpins share one segment; *focal*
    when that segment's span (or overhang beyond the cluster interval) is small
    relative to the cluster genomic span.
    """
    if cluster_summary.empty:
        return pd.DataFrame()

    loc = (
        long.loc[long["entity_level"] == "locus"]
        .drop_duplicates(subset=["participant", "entity_id"], keep="first")
        .copy()
    )
    loc["entity_id"] = loc["entity_id"].astype(str)
    loc["copy_number"] = pd.to_numeric(loc["copy_number"], errors="coerce")
    altered_states = {"gain", "amp", "loss", "deep_del"}

    rows: List[dict] = []
    for _, cl in cluster_summary.iterrows():
        block = _loci_in_cluster_block(locus_map, cl)
        if block.empty:
            continue
        lids = block["locus_id"].astype(str).tolist()
        cluster_start = int(block["start"].min())
        cluster_end = int(block["end"].max())
        cluster_span_mb = (cluster_end - cluster_start) / 1e6
        cluster_id = f"{cl['chrom']}:{cl['cluster_start_mb']:.3f}-{cl['cluster_end_mb']:.3f}"
        focal_ratio_mb = max(focal_cluster_multiplier * cluster_span_mb, 0.5)

        sub = loc.loc[loc["entity_id"].isin(lids)]
        per_participant: List[dict] = []
        for participant, grp in sub.groupby("participant"):
            seg_labels = grp["overlap_segment"].dropna().astype(str)
            seg_labels = seg_labels[seg_labels.str.contains(":", regex=False)]
            if seg_labels.empty:
                continue
            seg_label = seg_labels.mode().iloc[0]
            parsed = _parse_overlap_segment(seg_label)
            if parsed is None:
                continue
            _, seg_start, seg_end = parsed
            seg_mb = (seg_end - seg_start) / 1e6
            overhang_mb = max(0, cluster_start - seg_start) / 1e6 + max(0, seg_end - cluster_end) / 1e6
            per_participant.append(
                {
                    "same_segment": seg_labels.nunique() == 1,
                    "seg_mb": seg_mb,
                    "overhang_mb": overhang_mb,
                    "contains_cluster": seg_start <= cluster_start and seg_end >= cluster_end,
                    "focal_overhang": overhang_mb <= focal_overhang_mb,
                    "focal_segment": seg_mb <= focal_segment_mb,
                    "focal_vs_cluster": seg_mb <= focal_ratio_mb,
                    "altered": grp["cn_state"].astype(str).isin(altered_states).any(),
                    "amp_or_gain": (grp["copy_number"] >= 3).any()
                    or grp["cn_state"].astype(str).isin({"gain", "amp"}).any(),
                }
            )

        pdf = pd.DataFrame(per_participant)
        if pdf.empty:
            continue
        altered = pdf.loc[pdf["altered"]]
        amp = pdf.loc[pdf["amp_or_gain"]]

        def _pct(series: pd.Series) -> float:
            return round(100.0 * float(series.mean()), 2) if len(series) else np.nan

        rows.append(
            {
                "cluster_id": cluster_id,
                "chrom": cl["chrom"],
                "cluster_start_mb": cl["cluster_start_mb"],
                "cluster_end_mb": cl["cluster_end_mb"],
                "cluster_span_kb": round(cluster_span_mb * 1e3, 2),
                "n_hairpin_loci": int(len(block)),
                "mean_locus_delta_vs_diploid": cl["mean_locus_delta_vs_diploid"],
                "example_loci": cl.get("example_loci", ""),
                "n_participants": int(len(pdf)),
                "pct_loci_same_segment": _pct(pdf["same_segment"]),
                "median_segment_mb": round(float(pdf["seg_mb"].median()), 2),
                "median_overhang_mb": round(float(pdf["overhang_mb"].median()), 2),
                "pct_segment_contains_cluster": _pct(pdf["contains_cluster"]),
                "pct_focal_overhang_le_5mb": _pct(pdf["focal_overhang"]),
                "pct_focal_segment_le_10mb": _pct(pdf["focal_segment"]),
                "pct_focal_segment_le_10x_cluster": _pct(pdf["focal_vs_cluster"]),
                "pct_focal_overhang_amp": _pct(amp["focal_overhang"]) if len(amp) else np.nan,
                "pct_focal_10x_cluster_amp": _pct(amp["focal_vs_cluster"]) if len(amp) else np.nan,
                "n_amp_or_gain_participants": int(len(amp)),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("mean_locus_delta_vs_diploid", ascending=False, na_position="last")


def _chrom_offsets(catalog: pd.DataFrame, *, gap: float = 5e6) -> Dict[str, float]:
    offsets: Dict[str, float] = {}
    x = 0.0
    for chrom in CHR_ORDER:
        if chrom not in set(catalog["chrom"]):
            continue
        offsets[chrom] = x
        sub = catalog.loc[catalog["chrom"] == chrom]
        x += float(sub["end"].max()) + gap
    return offsets


def build_locus_genome_map(
    long: pd.DataFrame,
    *,
    catalog: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    catalog = _load_locus_catalog() if catalog is None else catalog
    rows: List[dict] = []
    cohort_cn = _mean_locus_cn(long)
    pam50_cols = {}
    for pam in PAM50_ORDER:
        if "PAM50_final" not in long.columns:
            break
        pam50_cols[pam] = _mean_locus_cn(long, stratum_col="PAM50_final", stratum_val=pam)

    for _, loc in catalog.iterrows():
        lid = str(loc["locus_id"])
        cn = float(cohort_cn.get(lid, np.nan))
        row = {
            "locus_id": lid,
            "locus_name": loc["locus_name"],
            "chrom": loc["chrom"],
            "start": int(loc["start"]),
            "end": int(loc["end"]),
            "midpoint": loc["midpoint"],
            "healthy_cn": DIPLOID_CN,
            "mean_copy_number": round(cn, 4) if np.isfinite(cn) else np.nan,
            "delta_vs_diploid": round(cn - DIPLOID_CN, 4) if np.isfinite(cn) else np.nan,
            "mature_accessions": loc.get("mature_accessions", ""),
        }
        for pam, ser in pam50_cols.items():
            pcn = float(ser.get(lid, np.nan))
            row[f"mean_cn_{pam}"] = round(pcn, 4) if np.isfinite(pcn) else np.nan
            row[f"delta_vs_diploid_{pam}"] = round(pcn - DIPLOID_CN, 4) if np.isfinite(pcn) else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    offsets = _chrom_offsets(catalog)
    out["genome_x"] = [_genome_x(r.chrom, r.midpoint, offsets) for r in out.itertuples()]
    return out.sort_values(["chrom", "start"]).reset_index(drop=True)


def build_mimat_paralog_maps(
    long: pd.DataFrame,
    *,
    catalog: Optional[pd.DataFrame] = None,
    weight_ref: Optional[pd.DataFrame] = None,
    arm_labels: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (paralog long table, MIMAT summary table)."""
    catalog = _load_locus_catalog() if catalog is None else catalog
    weight_ref = _load_weight_reference() if weight_ref is None else weight_ref
    arm_labels = _load_arm_labels() if arm_labels is None else arm_labels
    loc_meta = catalog.set_index("locus_id")
    cohort_locus_cn = _mean_locus_cn(long)
    pam50_locus_cn = {
        pam: _mean_locus_cn(long, stratum_col="PAM50_final", stratum_val=pam)
        for pam in PAM50_ORDER
        if "PAM50_final" in long.columns
    }

    paralog_rows: List[dict] = []
    summary_rows: List[dict] = []

    for mimat, grp in weight_ref.groupby("mimat"):
        mimat = str(mimat)
        parts: List[dict] = []
        for _, wr in grp.iterrows():
            lid = str(wr["locus_id"])
            if lid not in loc_meta.index:
                continue
            meta = loc_meta.loc[lid]
            cn = float(cohort_locus_cn.get(lid, np.nan))
            w = float(wr["w_norm"])
            healthy_contrib = float(wr["healthy_cn_contribution"])
            part = {
                "mimat": mimat,
                "locus_id": lid,
                "locus_name": meta["locus_name"],
                "chrom": meta["chrom"],
                "start": int(meta["start"]),
                "end": int(meta["end"]),
                "midpoint": meta["midpoint"],
                "w_norm": round(w, 6),
                "weight_ref_rpm": round(float(wr["weight_ref"]), 4),
                "healthy_cn_contribution": round(healthy_contrib, 6),
                "mean_locus_cn": round(cn, 4) if np.isfinite(cn) else np.nan,
                "delta_locus_vs_diploid": round(cn - DIPLOID_CN, 4) if np.isfinite(cn) else np.nan,
                "weighted_delta_contribution": round(w * (cn - DIPLOID_CN), 4) if np.isfinite(cn) else np.nan,
            }
            for pam, ser in pam50_locus_cn.items():
                pcn = float(ser.get(lid, np.nan))
                part[f"mean_locus_cn_{pam}"] = round(pcn, 4) if np.isfinite(pcn) else np.nan
                part[f"delta_locus_vs_diploid_{pam}"] = round(pcn - DIPLOID_CN, 4) if np.isfinite(pcn) else np.nan
            parts.append(part)
            paralog_rows.append(part)

        if not parts:
            continue
        pdf = pd.DataFrame(parts)
        w = pdf["w_norm"].values
        cn = pd.to_numeric(pdf["mean_locus_cn"], errors="coerce").values
        mask = np.isfinite(cn)
        weighted_cn = float(np.average(cn[mask], weights=w[mask])) if mask.any() and w[mask].sum() > 0 else np.nan
        median_cn = float(np.nanmedian(cn)) if mask.any() else np.nan
        healthy_ref = float(pdf["healthy_cn_contribution"].sum())  # = 2 when weights sum to 1

        srow = {
            "mimat": mimat,
            "n_paralog_loci": len(pdf),
            "healthy_ref_cn_weighted": round(healthy_ref, 4),
            "mean_cn_weighted": round(weighted_cn, 4) if np.isfinite(weighted_cn) else np.nan,
            "mean_cn_median": round(median_cn, 4) if np.isfinite(median_cn) else np.nan,
            "delta_vs_healthy_weighted": round(weighted_cn - healthy_ref, 4) if np.isfinite(weighted_cn) else np.nan,
            "delta_vs_diploid_median": round(median_cn - DIPLOID_CN, 4) if np.isfinite(median_cn) else np.nan,
        }
        for pam in pam50_locus_cn:
            pcn = pd.to_numeric(pdf[f"mean_locus_cn_{pam}"], errors="coerce").values
            m = np.isfinite(pcn)
            wcn = float(np.average(pcn[m], weights=w[m])) if m.any() and w[m].sum() > 0 else np.nan
            srow[f"mean_cn_weighted_{pam}"] = round(wcn, 4) if np.isfinite(wcn) else np.nan
            srow[f"delta_vs_healthy_weighted_{pam}"] = round(wcn - healthy_ref, 4) if np.isfinite(wcn) else np.nan
        summary_rows.append(srow)

    paralog = pd.DataFrame(paralog_rows)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = summary.drop_duplicates(subset=["mimat"], keep="first")
        summary = summary.merge(arm_labels.drop_duplicates(subset=["mimat"], keep="first"), on="mimat", how="left")
        summary = summary.sort_values("delta_vs_healthy_weighted", ascending=False, na_position="last")

    offsets = _chrom_offsets(catalog)
    if not paralog.empty:
        paralog["genome_x"] = [_genome_x(r.chrom, r.midpoint, offsets) for r in paralog.itertuples()]
        mimat_delta = summary.drop_duplicates("mimat").set_index("mimat")["delta_vs_healthy_weighted"]
        paralog["mimat_delta_vs_healthy_weighted"] = paralog["mimat"].map(mimat_delta)
        for pam in TUMOR_PAM50 + ["Normal"]:
            col = f"delta_vs_healthy_weighted_{pam}"
            if col in summary.columns:
                paralog[f"mimat_delta_vs_healthy_{pam}"] = paralog["mimat"].map(
                    summary.drop_duplicates("mimat").set_index("mimat")[col]
                )
    return paralog, summary


def _attach_paralog_pam50_deltas(paralog_map: pd.DataFrame, mimat_summary: pd.DataFrame) -> pd.DataFrame:
    """Merge per-PAM50 parent MIMAT deltas onto paralog rows for faceted plots."""
    if paralog_map.empty or mimat_summary.empty:
        return paralog_map
    delta_cols = [c for c in mimat_summary.columns if c.startswith("delta_vs_healthy_weighted_")]
    if not delta_cols:
        return paralog_map
    meta = mimat_summary[["mimat"] + delta_cols].drop_duplicates("mimat")
    rename = {c: c.replace("delta_vs_healthy_weighted_", "mimat_delta_vs_healthy_") for c in delta_cols}
    meta = meta.rename(columns=rename)
    out = paralog_map.merge(meta, on="mimat", how="left", suffixes=("", "_dup"))
    return out


def flag_interesting_cases(
    locus_map: pd.DataFrame,
    paralog_map: pd.DataFrame,
    mimat_summary: pd.DataFrame,
) -> pd.DataFrame:
    """Flag loci/MIMATs with extreme, subtype-skewed, or paralog-divergent CN."""
    rows: List[dict] = []

    def _add(entity_level, entity_id, label, chrom, flag_type, pam50, metric, value, detail, priority):
        rows.append(
            {
                "entity_level": entity_level,
                "entity_id": entity_id,
                "entity_label": label,
                "chrom": chrom,
                "flag_type": flag_type,
                "pam50": pam50,
                "metric": metric,
                "value": round(float(value), 4) if pd.notna(value) else np.nan,
                "detail": detail,
                "priority": priority,
            }
        )

    # --- MIMAT-level flags ---
    for _, r in mimat_summary.iterrows():
        mimat = str(r["mimat"])
        label = str(r.get("arm_label") or mimat)
        n_para = int(r.get("n_paralog_loci") or 1)
        d_norm = r.get("delta_vs_healthy_weighted_Normal", np.nan)
        d_wt = r.get("delta_vs_healthy_weighted", np.nan)
        d_med = r.get("delta_vs_diploid_median", np.nan)

        if label in CURATED_ARMS:
            best_pam, best_d = None, np.nan
            for pam in TUMOR_PAM50:
                d = r.get(f"delta_vs_healthy_weighted_{pam}", np.nan)
                if pd.notna(d) and (not pd.notna(best_d) or d > best_d):
                    best_pam, best_d = pam, d
            _add(
                "arm", mimat, label, "",
                "curated_amp_context", best_pam or "cohort", "delta_vs_healthy_weighted",
                best_d if pd.notna(best_d) else d_wt,
                CURATED_ARMS[label], 1,
            )

        for pam in TUMOR_PAM50:
            col = f"delta_vs_healthy_weighted_{pam}"
            d = r.get(col, np.nan)
            if pd.notna(d) and d >= 2.0:
                _add(
                    "arm", mimat, label, "",
                    "extreme_weighted_amp", pam, col, d,
                    f"Weighted MIMAT CN Δ≥+2 vs healthy (Σ2·w) in {pam}", 1,
                )

        if pd.notna(d_norm):
            for pam in TUMOR_PAM50:
                d = r.get(f"delta_vs_healthy_weighted_{pam}", np.nan)
                if pd.notna(d) and (d - d_norm) >= 1.5:
                    _add(
                        "arm", mimat, label, "",
                        "subtype_enriched_vs_normal", pam, f"delta_minus_normal_{pam}",
                        d - d_norm,
                        f"{pam} weighted Δ exceeds Normal by ≥1.5 ({d:.2f} vs Normal {d_norm:.2f})", 2,
                    )

        tumor_ds = [r.get(f"delta_vs_healthy_weighted_{p}", np.nan) for p in TUMOR_PAM50]
        tumor_ds = [float(x) for x in tumor_ds if pd.notna(x)]
        if len(tumor_ds) >= 2 and (max(tumor_ds) - min(tumor_ds)) >= 1.5:
            _add(
                "arm", mimat, label, "",
                "subtype_dispersion", "across_tumors", "delta_range_tumor_pam50",
                max(tumor_ds) - min(tumor_ds),
                f"Tumor PAM50 spread {min(tumor_ds):.2f}–{max(tumor_ds):.2f}", 2,
            )

        if n_para > 1 and pd.notna(d_wt) and pd.notna(d_med) and abs(d_wt - d_med) >= 0.4:
            _add(
                "arm", mimat, label, "",
                "paralog_weight_vs_median", "cohort", "weighted_minus_median_delta",
                d_wt - d_med,
                f"Multi-paralog: weighted Δ {d_wt:.2f} vs median-paralog Δ {d_med:.2f}", 2,
            )

    # --- Paralog spread within MIMAT ---
    for mimat, grp in paralog_map.groupby("mimat"):
        if len(grp) < 2:
            continue
        label = mimat_summary.loc[mimat_summary["mimat"] == mimat, "arm_label"]
        label = str(label.iloc[0]) if len(label) else str(mimat)
        spread = grp["delta_locus_vs_diploid"].max() - grp["delta_locus_vs_diploid"].min()
        if pd.notna(spread) and spread >= 1.0:
            locs = grp.sort_values("delta_locus_vs_diploid", ascending=False)
            hi = locs.iloc[0]
            lo = locs.iloc[-1]
            _add(
                "arm", mimat, label, "",
                "paralog_locus_spread", "cohort", "locus_delta_spread", spread,
                f"Paralog CN spread {spread:.2f}: {hi['locus_id']} Δ={hi['delta_locus_vs_diploid']:.2f} "
                f"vs {lo['locus_id']} Δ={lo['delta_locus_vs_diploid']:.2f}", 1,
            )
        for pam in TUMOR_PAM50:
            col = f"delta_locus_vs_diploid_{pam}"
            if col not in grp.columns:
                continue
            pam_spread = grp[col].max() - grp[col].min()
            if pd.notna(pam_spread) and pam_spread >= 1.5:
                _add(
                    "arm", mimat, label, "",
                    "paralog_locus_spread", pam, f"locus_delta_spread_{pam}", pam_spread,
                    f"{pam}: paralog loci differ by Δ={pam_spread:.2f}", 2,
                )

    # --- Locus-level flags ---
    for _, loc in locus_map.iterrows():
        lid = str(loc["locus_id"])
        label = str(loc["locus_name"])
        chrom = str(loc["chrom"])
        for pam in TUMOR_PAM50:
            col = f"delta_vs_diploid_{pam}"
            d = loc.get(col, np.nan)
            if pd.notna(d) and d >= 2.0:
                _add(
                    "locus", lid, label, chrom,
                    "locus_extreme_amp", pam, col, d,
                    f"Hairpin locus mean CN Δ≥+2 vs diploid in {pam}", 1,
                )
        if chrom == "chr8":
            d8 = max(loc.get(f"delta_vs_diploid_{p}", np.nan) for p in TUMOR_PAM50)
            if pd.notna(d8) and d8 >= 1.5:
                _add(
                    "locus", lid, label, chrom,
                    "chr8_amp_neighborhood", "tumor_max", "max_tumor_delta", d8,
                    "chr8 amplicon neighborhood (≥+1.5 Δ in some tumor subtype)", 2,
                )

    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out = out.sort_values(["priority", "value"], ascending=[True, False], na_position="last")
    out = out.drop_duplicates(
        subset=["entity_level", "entity_id", "flag_type", "pam50", "metric"], keep="first"
    )
    return out.reset_index(drop=True)


def _plot_manhattan(
    df: pd.DataFrame,
    *,
    y_col: str,
    title: str,
    out_path: Path,
    size_col: Optional[str] = None,
    alpha: float = 0.85,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot = df.dropna(subset=["genome_x", y_col]).copy()
    if plot.empty:
        return

    chrom_palette = plt.cm.tab20.colors
    chroms = [c for c in CHR_ORDER if c in set(plot["chrom"])]
    color_map = {c: chrom_palette[i % len(chrom_palette)] for i, c in enumerate(chroms)}

    fig, ax = plt.subplots(figsize=(14, 4.5), dpi=150)
    sizes = None
    if size_col and size_col in plot.columns:
        s = pd.to_numeric(plot[size_col], errors="coerce").fillna(0.05)
        sizes = 12 + 80 * (s / max(s.max(), 1e-9))
    for chrom in chroms:
        sub = plot.loc[plot["chrom"] == chrom]
        ax.scatter(
            sub["genome_x"],
            sub[y_col],
            c=[color_map[chrom]],
            s=sizes.loc[sub.index] if sizes is not None else 14,
            alpha=alpha,
            linewidths=0,
            label=chrom,
        )
    ax.axhline(0, color="0.35", lw=0.8, ls="--")
    ax.set_xlabel("Genome position (hg38)")
    ax.set_ylabel("Δ copy number vs healthy reference")
    ax.set_title(title)
    ax.set_xticks([])
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[c], markersize=6, label=c) for c in chroms[:12]]
    ax.legend(handles=handles, ncol=6, fontsize=7, loc="upper right", frameon=False)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_pam50_faceted(
    df: pd.DataFrame,
    *,
    y_cols: Sequence[str],
    panel_labels: Sequence[str],
    title: str,
    out_path: Path,
    size_col: Optional[str] = None,
    alpha: float = 0.75,
    highlight: Optional[pd.DataFrame] = None,
) -> None:
    """Stacked PAM50 panels sharing genome x-axis."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(y_cols)
    if n == 0:
        return

    chrom_palette = plt.cm.tab20.colors
    chroms = [c for c in CHR_ORDER if c in set(df.get("chrom", pd.Series(dtype=str)))]
    color_map = {c: chrom_palette[i % len(chrom_palette)] for i, c in enumerate(chroms)}

    fig, axes = plt.subplots(n, 1, figsize=(14, 2.6 * n), dpi=150, sharex=True)
    if n == 1:
        axes = [axes]

    highlight_ids = set()
    highlight_mimat = set()
    if highlight is not None and not highlight.empty:
        highlight_ids = set(highlight.loc[highlight["entity_level"] == "locus", "entity_id"].astype(str))
        highlight_mimat = set(highlight.loc[highlight["entity_level"] == "arm", "entity_id"].astype(str))

    for ax, y_col, plab in zip(axes, y_cols, panel_labels):
        plot = df.dropna(subset=["genome_x", y_col]).copy()
        sizes = None
        if size_col and size_col in plot.columns:
            s = pd.to_numeric(plot[size_col], errors="coerce").fillna(0.05)
            sizes = 10 + 70 * (s / max(s.max(), 1e-9))
        for chrom in chroms:
            sub = plot.loc[plot["chrom"] == chrom]
            if sub.empty:
                continue
            ax.scatter(
                sub["genome_x"], sub[y_col], c=[color_map[chrom]],
                s=sizes.loc[sub.index] if sizes is not None else 12,
                alpha=alpha, linewidths=0,
            )
        hi = pd.DataFrame()
        if "mimat" in plot.columns and highlight_mimat:
            hi = plot.loc[plot["mimat"].astype(str).isin(highlight_mimat)]
        elif "locus_id" in plot.columns and highlight_ids:
            hi = plot.loc[plot["locus_id"].astype(str).isin(highlight_ids)]
        if not hi.empty:
                ax.scatter(
                    hi["genome_x"], hi[y_col], facecolors="none", edgecolors="black",
                    s=55, linewidths=0.9, alpha=0.9, zorder=5,
                )
        ax.axhline(0, color="0.35", lw=0.7, ls="--")
        ax.set_ylabel("Δ CN")
        ax.set_title(plab, fontsize=10, loc="left")
        ax.set_xticks([])

    axes[-1].set_xlabel("Genome position (hg38)")
    fig.suptitle(title, fontsize=12, y=1.002)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_genome_maps(
    long: pd.DataFrame,
    out_dir: Path,
    *,
    write_figures: bool = True,
) -> Dict[str, Path]:
    maps_dir = out_dir / "maps"
    maps_dir.mkdir(parents=True, exist_ok=True)

    locus_map = build_locus_genome_map(long)
    paralog_map, mimat_summary = build_mimat_paralog_maps(long)

    paths = {
        "locus_tsv": maps_dir / "mirna_cnv_locus_genome_map.tsv",
        "mimat_paralog_tsv": maps_dir / "mirna_cnv_mimat_paralog_map.tsv",
        "mimat_summary_tsv": maps_dir / "mirna_cnv_mimat_summary.tsv",
    }
    locus_map.to_csv(paths["locus_tsv"], sep="\t", index=False)
    paralog_map.to_csv(paths["mimat_paralog_tsv"], sep="\t", index=False)
    mimat_summary.to_csv(paths["mimat_summary_tsv"], sep="\t", index=False)

    chrom_summary = build_chromosome_cn_summary(locus_map, paralog_map)
    arm_summary = build_chromosome_arm_cn_summary(locus_map, paralog_map)
    cluster_summary = build_dense_locus_clusters(locus_map)
    cluster_focality = build_cluster_segment_focality(long, locus_map, cluster_summary)
    paths["chromosome_summary_tsv"] = maps_dir / "mirna_cnv_chromosome_summary.tsv"
    paths["chromosome_arm_summary_tsv"] = maps_dir / "mirna_cnv_chromosome_arm_summary.tsv"
    paths["locus_cluster_summary_tsv"] = maps_dir / "mirna_cnv_locus_cluster_summary.tsv"
    paths["cluster_segment_focality_tsv"] = maps_dir / "mirna_cnv_cluster_segment_focality.tsv"
    chrom_summary.to_csv(paths["chromosome_summary_tsv"], sep="\t", index=False)
    arm_summary.to_csv(paths["chromosome_arm_summary_tsv"], sep="\t", index=False)
    cluster_summary.to_csv(paths["locus_cluster_summary_tsv"], sep="\t", index=False)
    cluster_focality.to_csv(paths["cluster_segment_focality_tsv"], sep="\t", index=False)

    flags = flag_interesting_cases(locus_map, paralog_map, mimat_summary)
    flags_path = maps_dir / "mirna_cnv_interesting_cases.tsv"
    flags.to_csv(flags_path, sep="\t", index=False)
    paths["interesting_cases_tsv"] = flags_path

    # Highlight priority-1 entity ids on PAM50 figures
    hi_entities = flags.loc[flags["priority"] == 1] if not flags.empty else pd.DataFrame()

    if write_figures:
        locus_fig = maps_dir / "mirna_cnv_locus_genome_map.png"
        mimat_fig = maps_dir / "mirna_cnv_mimat_paralog_genome_map.png"
        _plot_manhattan(
            locus_map,
            y_col="delta_vs_diploid",
            title="miRNA hairpin loci — Δ CN vs diploid (CN=2)",
            out_path=locus_fig,
        )
        _plot_manhattan(
            paralog_map,
            y_col="mimat_delta_vs_healthy_weighted",
            title="MIMAT paralog loci — parent MIMAT Δ vs weighted healthy (Σ 2·w_locus)",
            out_path=mimat_fig,
            size_col="w_norm",
            alpha=0.55,
        )
        paths["locus_fig"] = locus_fig
        paths["mimat_fig"] = mimat_fig

        pam_present = [p for p in PAM50_ORDER if f"delta_vs_diploid_{p}" in locus_map.columns]
        if pam_present:
            locus_pam_fig = maps_dir / "mirna_cnv_locus_genome_map_by_pam50.png"
            _plot_pam50_faceted(
                locus_map,
                y_cols=[f"delta_vs_diploid_{p}" for p in pam_present],
                panel_labels=[f"PAM50 {p} — locus Δ vs diploid (CN=2)" for p in pam_present],
                title="miRNA locus CN by PAM50 (Δ vs diploid)",
                out_path=locus_pam_fig,
                highlight=hi_entities.loc[hi_entities["entity_level"] == "locus"] if not hi_entities.empty else None,
            )
            paths["locus_pam50_fig"] = locus_pam_fig

        paralog_map = _attach_paralog_pam50_deltas(paralog_map, mimat_summary)
        mimat_pam_cols = [f"mimat_delta_vs_healthy_{p}" for p in PAM50_ORDER if f"mimat_delta_vs_healthy_{p}" in paralog_map.columns]
        if mimat_pam_cols:
            mimat_pam_fig = maps_dir / "mirna_cnv_mimat_paralog_genome_map_by_pam50.png"
            _plot_pam50_faceted(
                paralog_map,
                y_cols=mimat_pam_cols,
                panel_labels=[f"PAM50 {c.replace('mimat_delta_vs_healthy_', '')} — MIMAT weighted Δ" for c in mimat_pam_cols],
                title="MIMAT paralog CN by PAM50 (parent weighted Δ vs Σ2·w_locus)",
                out_path=mimat_pam_fig,
                size_col="w_norm",
                alpha=0.5,
                highlight=hi_entities.loc[hi_entities["entity_level"] == "arm"] if not hi_entities.empty else None,
            )
            paths["mimat_pam50_fig"] = mimat_pam_fig

        if not flags.empty:
            summary_json = {
                "n_flags": int(len(flags)),
                "n_priority1": int((flags["priority"] == 1).sum()),
                "flag_types": flags["flag_type"].value_counts().to_dict(),
                "top_priority1": flags.loc[flags["priority"] == 1].head(20)[
                    ["entity_label", "flag_type", "pam50", "value", "detail"]
                ].to_dict(orient="records"),
            }
            summary_path = maps_dir / "mirna_cnv_interesting_cases_summary.json"
            summary_path.write_text(json.dumps(summary_json, indent=2), encoding="utf-8")
            paths["interesting_cases_json"] = summary_path

    n_flags = len(flags) if not flags.empty else 0
    n_p1 = int((flags["priority"] == 1).sum()) if not flags.empty else 0
    print(f"[mirna_cnv_genome_maps] locus rows: {len(locus_map):,}; mimat summary: {len(mimat_summary):,}; "
          f"paralog rows: {len(paralog_map):,}; interesting flags: {n_flags} ({n_p1} priority-1)")

    from mirna_hallmark.analyses.cnv_locus.mirna_cnv_subtype_depth import write_subtype_depth

    paths.update(write_subtype_depth(long, out_dir))
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--long-cnv", type=Path, default=C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz")
    ap.add_argument("--out-dir", type=Path, default=C.MIRNA_LOCUS_CNV_DIR)
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()
    long = pd.read_csv(args.long_cnv, sep="\t", low_memory=False)
    write_genome_maps(long, args.out_dir, write_figures=not args.no_figures)


if __name__ == "__main__":
    main()
