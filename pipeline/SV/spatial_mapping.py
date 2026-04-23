"""
Spatial mapping of SVs to genes, regulatory elements, and lncRNAs.

Functions for:
- Computing signed distances between intervals
- Classifying SV-gene spatial relationships
- Mapping SVs to genes, elements, and lncRNAs
"""

from typing import Tuple, Dict, Any, Optional, List

import numpy as np
import pandas as pd

from .vcf_loader import add_bnd_remote_coords


# Shared empty frame used as a fallback for unknown chromosomes in chrom indexes.
_EMPTY_DF = pd.DataFrame()


def _chrom_index(df: Optional[pd.DataFrame], col: str = "chrom") -> Dict[Any, pd.DataFrame]:
    """Group ``df`` by ``col`` once so per-SV code can do O(1) chrom lookup.

    The values are views into ``df``; callers must treat them as read-only.
    """
    if df is None or df.empty or col not in df.columns:
        return {}
    return {key: grp for key, grp in df.groupby(col, sort=False)}


# =============================================================================
# DISTANCE CALCULATIONS
# =============================================================================

def compute_signed_distance_with_overlap(
    sv_start: int,
    sv_end: int,
    target_start: int,
    target_end: int,
    strand: str = "+",
) -> Tuple[int, Optional[int], Optional[int]]:
    """
    Compute signed distance between SV and target interval.
    
    Returns:
        Tuple of (signed_distance, overlap_start, overlap_end)
        - signed_distance = 0 if overlapping
        - overlap_start/end = genomic coordinates of overlap, else None
    """
    # Check overlap
    if sv_end >= target_start and sv_start <= target_end:
        overlap_start = max(sv_start, target_start)
        overlap_end = min(sv_end, target_end)
        return 0, overlap_start, overlap_end

    # No overlap
    overlap_start = None
    overlap_end = None

    # SV completely before target
    if sv_end < target_start:
        d = target_start - sv_end
        if strand == '+':
            return -d, overlap_start, overlap_end
        else:
            return +d, overlap_start, overlap_end

    # SV completely after target
    if sv_start > target_end:
        d = sv_start - target_end
        if strand == '+':
            return +d, overlap_start, overlap_end
        else:
            return -d, overlap_start, overlap_end

    return 0, overlap_start, overlap_end


# =============================================================================
# SV ↔ miRNA mapping (pre-miRNA loci with mature-arm metadata)
# =============================================================================

def _map_span_to_mirnas(
    df: pd.DataFrame,
    mirnas: pd.DataFrame,
    svtype: str,
    window: int,
) -> pd.DataFrame:
    """
    Map span SVs (DEL/DUP/INV) to nearby/overlapping miRNA loci.

    Expects mirnas with columns: chrom,start,end,strand,gene_name,gene_id
    and optionally: mature_names,mature_accessions.
    """
    out = df.copy()
    mask = out.get("SVTYPE") == svtype
    if mask is None or not mask.any():
        return out

    sub = out.loc[mask].copy()
    mirnas = mirnas.copy()
    mirnas["chrom"] = mirnas["chrom"].astype(str)

    hits_col: List[list] = []
    for _, r in sub.iterrows():
        chrom = str(r.get("chrom"))
        sv_start = int(r.get("pos"))
        sv_end = int(r.get("END")) if pd.notna(r.get("END")) else sv_start

        mc = mirnas.loc[mirnas["chrom"] == chrom]
        if mc.empty:
            hits_col.append([])
            continue

        cand = mc.loc[(mc["end"] >= sv_start - window) & (mc["start"] <= sv_end + window)]
        if cand.empty:
            hits_col.append([])
            continue

        hits = []
        for _, m in cand.iterrows():
            ms, me = int(m["start"]), int(m["end"])
            strand = str(m.get("strand", "+"))
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                sv_start, sv_end, ms, me, strand=strand
            )
            if abs(signed_dist) > window:
                continue
            overlap_bp = 0
            overlap_percent = 0.0
            if overlap_start is not None and overlap_end is not None and overlap_end >= overlap_start:
                overlap_bp = overlap_end - overlap_start + 1
                overlap_percent = (overlap_bp / (me - ms + 1)) * 100
            hit = {
                "gene_name": m.get("gene_name"),
                "gene_id": m.get("gene_id"),
                "strand": strand,
                "signed_dist": int(signed_dist),
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                "overlap_bp": int(overlap_bp),
                "overlap_percent": float(overlap_percent),
                "region_hit": ("overlaps" if signed_dist == 0 else ("proximal" if abs(signed_dist) <= 5000 else "distal")),
                "hit_side": "span",
            }
            for extra in (
                "mature_names",
                "mature_accessions",
                "mirbase_mature_id",
                "pre_gene_name",
                "pre_gene_id",
            ):
                if extra in m and pd.notna(m.get(extra)):
                    hit[extra] = m.get(extra)
            hits.append(hit)
        hits_col.append(hits)

    sub["mir_hits"] = hits_col
    out.loc[sub.index, "mir_hits"] = sub["mir_hits"]
    return out


def _map_ins_to_mirnas(df: pd.DataFrame, mirnas: pd.DataFrame, window: int) -> pd.DataFrame:
    out = df.copy()
    mask = out.get("SVTYPE") == "INS"
    if mask is None or not mask.any():
        return out
    sub = out.loc[mask].copy()
    mirnas = mirnas.copy()
    mirnas["chrom"] = mirnas["chrom"].astype(str)

    hits_col: List[list] = []
    for _, r in sub.iterrows():
        chrom = str(r.get("chrom"))
        pos = int(r.get("pos"))
        mc = mirnas.loc[mirnas["chrom"] == chrom]
        if mc.empty:
            hits_col.append([])
            continue
        cand = mc.loc[(mc["end"] >= pos - window) & (mc["start"] <= pos + window)]
        hits = []
        for _, m in cand.iterrows():
            ms, me = int(m["start"]), int(m["end"])
            strand = str(m.get("strand", "+"))
            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(pos, pos, ms, me, strand=strand)
            if abs(signed_dist) > window:
                continue
            overlap_bp = 0
            overlap_percent = 0.0
            if overlap_start is not None and overlap_end is not None and overlap_end >= overlap_start:
                overlap_bp = overlap_end - overlap_start + 1
                overlap_percent = (overlap_bp / (me - ms + 1)) * 100
            hit = {
                "gene_name": m.get("gene_name"),
                "gene_id": m.get("gene_id"),
                "strand": strand,
                "signed_dist": int(signed_dist),
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                "overlap_bp": int(overlap_bp),
                "overlap_percent": float(overlap_percent),
                "region_hit": ("overlaps" if signed_dist == 0 else ("proximal" if abs(signed_dist) <= 5000 else "distal")),
                "hit_side": "point",
            }
            for extra in (
                "mature_names",
                "mature_accessions",
                "mirbase_mature_id",
                "pre_gene_name",
                "pre_gene_id",
            ):
                if extra in m and pd.notna(m.get(extra)):
                    hit[extra] = m.get(extra)
            hits.append(hit)
        hits_col.append(hits)
    sub["mir_hits"] = hits_col
    out.loc[sub.index, "mir_hits"] = sub["mir_hits"]
    return out


def _map_bnds_to_mirnas(df: pd.DataFrame, mirnas: pd.DataFrame, window: int) -> pd.DataFrame:
    out = df.copy()
    mask = out.get("SVTYPE") == "BND"
    if mask is None or not mask.any():
        return out
    sub = add_bnd_remote_coords(out.loc[mask].copy())
    mirnas = mirnas.copy()
    mirnas["chrom"] = mirnas["chrom"].astype(str)

    hits_col: List[list] = []
    for _, r in sub.iterrows():
        chrom1 = str(r.get("chrom"))
        pos1 = int(r.get("pos"))
        chrom2 = r.get("bnd_remote_chrom")
        pos2 = r.get("bnd_remote_pos")
        hits: List[dict] = []

        def process_side(chrom: str, pos: int, side: str) -> None:
            mc = mirnas.loc[mirnas["chrom"] == str(chrom)]
            if mc.empty:
                return
            cand = mc.loc[(mc["end"] >= pos - window) & (mc["start"] <= pos + window)]
            for _, m in cand.iterrows():
                ms, me = int(m["start"]), int(m["end"])
                strand = str(m.get("strand", "+"))
                signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(pos, pos, ms, me, strand=strand)
                if abs(signed_dist) > window:
                    continue
                overlap_bp = 0
                overlap_percent = 0.0
                if overlap_start is not None and overlap_end is not None and overlap_end >= overlap_start:
                    overlap_bp = overlap_end - overlap_start + 1
                    overlap_percent = (overlap_bp / (me - ms + 1)) * 100
                hit = {
                    "gene_name": m.get("gene_name"),
                    "gene_id": m.get("gene_id"),
                    "strand": strand,
                    "signed_dist": int(signed_dist),
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_bp": int(overlap_bp),
                    "overlap_percent": float(overlap_percent),
                    "region_hit": ("overlaps" if signed_dist == 0 else ("proximal" if abs(signed_dist) <= 5000 else "distal")),
                    "hit_side": side,
                }
                for extra in (
                    "mature_names",
                    "mature_accessions",
                    "mirbase_mature_id",
                    "pre_gene_name",
                    "pre_gene_id",
                ):
                    if extra in m and pd.notna(m.get(extra)):
                        hit[extra] = m.get(extra)
                hits.append(hit)

        process_side(chrom1, pos1, "bp1")
        if isinstance(chrom2, str) and pd.notna(pos2):
            process_side(str(chrom2), int(pos2), "bp2")
        hits_col.append(hits)

    sub["mir_hits"] = hits_col
    out.loc[sub.index, "mir_hits"] = sub["mir_hits"]
    return out


# =============================================================================
# SV-GENE CLASSIFICATION
# =============================================================================

from typing import Any, Dict
import pandas as pd

def classify_span_sv_gene_hit(
    sv_start: int,
    sv_end: int,
    gene_start: int,
    gene_end: int,
    strand: str,
    signed_dist: int,
    gene_features: pd.DataFrame,
    promoter_up: int = 2000,
    promoter_down: int = 500,
    proximal_window: int = 5000,
) -> Dict[str, Any]:
    """
    Fixed v2:
    - initializes hits_cds
    - fixes MANE logic (no "last row wins" bug)
    - makes intron_only consistent: intron_only only if gene-body overlap AND no CDS/exon/UTR hit
    - keeps start/stop codons correctly separated
    """

    # --- gene body overlap ---
    inter_body_start = max(sv_start, gene_start)
    inter_body_end = min(sv_end, gene_end)
    if inter_body_end >= inter_body_start:
        overlap_bp = inter_body_end - inter_body_start + 1
        overlap_percent = (overlap_bp / (gene_end - gene_start + 1)) * 100
        gene_body_flag = 1
    else:
        overlap_bp = 0
        overlap_percent = 0.0
        gene_body_flag = 0

    # --- promoter region (strand-aware) ---
    if strand == "+":
        tss = gene_start
        prom_start = tss - promoter_up
        prom_end = tss + promoter_down
    else:
        tss = gene_end
        prom_start = tss - promoter_down
        prom_end = tss + promoter_up

    inter_prom_start = max(sv_start, prom_start)
    inter_prom_end = min(sv_end, prom_end)
    promoter_flag = 1 if inter_prom_end >= inter_prom_start else 0
    promoter_len_bp = int(prom_end - prom_start + 1)
    promoter_overlap_bp = int(inter_prom_end - inter_prom_start + 1) if promoter_flag else 0
    promoter_overlap_frac = float(promoter_overlap_bp / promoter_len_bp) if promoter_len_bp > 0 else 0.0

    # --- overlap features (robust + avoids iterrows logic bugs) ---
    if gene_features is None or gene_features.empty:
        ov = gene_features
    else:
        gf = gene_features.copy()

        # Ensure numeric starts/ends (avoid surprises from strings)
        gf["start"] = pd.to_numeric(gf["start"], errors="coerce")
        gf["end"] = pd.to_numeric(gf["end"], errors="coerce")

        # Keep only rows with valid coords
        gf = gf.dropna(subset=["start", "end"])

        ov = gf.loc[~((sv_end < gf["start"]) | (sv_start > gf["end"]))].copy()

    # Default flags
    hits_exon = False
    hits_utr = False
    hits_cds = False
    hits_start_codon = False
    hits_stop_codon = False

    transcript_id = None
    transcript_type = None
    mane_transcript_flag = 0

    mane_cds_hit = False  # "SV overlaps a CDS that belongs to a MANE transcript"

    if ov is not None and not ov.empty:
        # Normalize missing columns
        if "feature" not in ov.columns:
            ov["feature"] = None
        if "is_MANE" not in ov.columns:
            ov["is_MANE"] = False

        # Feature hits
        feats = ov["feature"].astype(str)

        hits_exon = (feats == "exon").any()
        hits_utr = (feats == "UTR").any()
        hits_cds = (feats == "CDS").any()
        hits_start_codon = (feats == "start_codon").any()
        hits_stop_codon = (feats == "stop_codon").any()

        # Record the specific overlapped feature intervals.
        # If the features table contains stable exon IDs / exon numbers from the original GTF
        # conversion, we emit them. Otherwise we fall back to deterministic identifiers based on
        # transcript_id + feature + coordinates.
        def _feature_interval_id(r) -> str:
            tid = r.get("transcript_id")
            tid = "" if tid is None or (isinstance(tid, float) and pd.isna(tid)) else str(tid)
            f = str(r.get("feature", "") or "")
            s = r.get("start")
            e = r.get("end")
            try:
                s_i = int(s) if pd.notna(s) else -1
            except Exception:
                s_i = -1
            try:
                e_i = int(e) if pd.notna(e) else -1
            except Exception:
                e_i = -1
            return f"{tid}:{f}:{s_i}-{e_i}"

        def _feature_exon_id(r) -> Optional[str]:
            exid = r.get("exon_id")
            if exid is not None and not (isinstance(exid, float) and pd.isna(exid)):
                return str(exid)
            # Optional: exon_number as a weaker identifier (transcript-specific)
            exn = r.get("exon_number")
            tid = r.get("transcript_id")
            if exn is not None and not (isinstance(exn, float) and pd.isna(exn)) and tid is not None:
                try:
                    return f"{str(tid)}:exon_number:{int(float(exn))}"
                except Exception:
                    return None
            return None

        exon_interval_ids = [
            _feature_interval_id(r) for _, r in ov.loc[feats == "exon"].iterrows()
        ]
        exon_ids = [
            x for x in (_feature_exon_id(r) for _, r in ov.loc[feats == "exon"].iterrows()) if x
        ]
        cds_interval_ids = [
            _feature_interval_id(r) for _, r in ov.loc[feats == "CDS"].iterrows()
        ]
        utr_interval_ids = [
            _feature_interval_id(r) for _, r in ov.loc[feats == "UTR"].iterrows()
        ]
        start_codon_interval_ids = [
            _feature_interval_id(r) for _, r in ov.loc[feats == "start_codon"].iterrows()
        ]
        stop_codon_interval_ids = [
            _feature_interval_id(r) for _, r in ov.loc[feats == "stop_codon"].iterrows()
        ]

        # Determine transcript (prefer MANE transcript if present among overlapped transcript rows)
        transcript_rows = ov.loc[feats == "transcript"].copy()
        if not transcript_rows.empty:
            # If is_MANE is present on transcript rows, prefer it
            if transcript_rows["is_MANE"].astype(bool).any():
                trow = transcript_rows.loc[transcript_rows["is_MANE"].astype(bool)].iloc[0]
                mane_transcript_flag = 1
            else:
                trow = transcript_rows.iloc[0]

            transcript_id = trow.get("transcript_id")
            transcript_type = trow.get("transcript_type")

        # MANE_CDS: robust to is_MANE being stored either on CDS rows or transcript rows
        # Case A: is_MANE annotated per-CDS row
        cds_rows = ov.loc[feats == "CDS"].copy()
        if not cds_rows.empty and cds_rows["is_MANE"].astype(bool).any():
            mane_cds_hit = True
        else:
            # Case B: is_MANE annotated on transcript rows; if any overlapped MANE transcript AND any overlapped CDS,
            # treat as MANE_CDS (best-effort without transcript_id mapping)
            if mane_transcript_flag == 1 and hits_cds:
                mane_cds_hit = True

    cds_flag = 1 if hits_cds else 0
    exon_flag = 1 if hits_exon else 0
    utr_flag = 1 if hits_utr else 0
    start_codon_flag = 1 if hits_start_codon else 0
    stop_codon_flag = 1 if hits_stop_codon else 0

    mane_cds_flag = 1 if mane_cds_hit else 0

    # intron_only only if within gene body and no exon/CDS/UTR hit
    intron_only_flag = 1 if (gene_body_flag and not (cds_flag or exon_flag or utr_flag)) else 0

    # --- classify region_hit (coarse) ---
    region_hit_parts = []
    if promoter_flag:
        region_hit_parts.append("promoter")

    if overlap_bp > 0:
        if mane_cds_flag:
            region_hit_parts.append("MANE_CDS")
        elif cds_flag:
            region_hit_parts.append("CDS")
        elif exon_flag:
            region_hit_parts.append("exon")
        elif utr_flag:
            region_hit_parts.append("UTR")
        else:
            region_hit_parts.append("intron")
    else:
        if not promoter_flag:
            if abs(signed_dist) <= proximal_window:
                region_hit_parts.append("upstream_5kb" if signed_dist < 0 else "downstream_5kb")
            else:
                region_hit_parts.append("intergenic")

    region_hit = ",".join(region_hit_parts) + ("," if region_hit_parts else "")

    upstream_5kb_flag = int("upstream_5kb" in region_hit)
    downstream_5kb_flag = int("downstream_5kb" in region_hit)

    return {
        "overlap_bp": overlap_bp,
        "overlap_percent": overlap_percent,
        "promoter_flag": promoter_flag,
        # Promoter overlap is computed against the strand-aware promoter window
        # [TSS - promoter_up, TSS + promoter_down] (or swapped for '-' strand).
        "promoter_len_bp": promoter_len_bp,
        "promoter_overlap_bp": promoter_overlap_bp,
        "promoter_overlap_frac": promoter_overlap_frac,
        "gene_body_flag": gene_body_flag,
        "mane_cds_flag": mane_cds_flag,
        "cds_flag": cds_flag,
        "utr_flag": utr_flag,
        "exon_flag": exon_flag,
        "mane_transcript_flag": mane_transcript_flag,
        "intron_only_flag": intron_only_flag,
        "upstream_5kb_flag": upstream_5kb_flag,
        "downstream_5kb_flag": downstream_5kb_flag,
        "region_hit": region_hit,
        "hit_side": "span",
        "stop_codon_flag": stop_codon_flag,
        "start_codon_flag": start_codon_flag,
        "transcript_id": transcript_id,
        "transcript_type": transcript_type,
        # Feature-interval attribution (deterministic IDs; see above).
        "exon_interval_ids": exon_interval_ids if ov is not None and not ov.empty else [],
        "exon_ids": exon_ids if ov is not None and not ov.empty else [],
        "cds_interval_ids": cds_interval_ids if ov is not None and not ov.empty else [],
        "utr_interval_ids": utr_interval_ids if ov is not None and not ov.empty else [],
        "start_codon_interval_ids": start_codon_interval_ids if ov is not None and not ov.empty else [],
        "stop_codon_interval_ids": stop_codon_interval_ids if ov is not None and not ov.empty else [],
    }

def classify_span_sv_elem_hit(
    sv_start: int,
    sv_end: int,
    elem_start: int,
    elem_end: int,
    signed_dist: int,
    proximal_window: int = 5000,
) -> Dict[str, Any]:
    """
    Classify SV-element spatial relationship.
    
    Args:
        sv_start: SV start position
        sv_end: SV end position
        elem_start: Element start position
        elem_end: Element end position
        signed_dist: Signed distance from SV to element
        proximal_window: Proximal window size
    
    Returns:
        Dict with region classification flags
    """
    # Element overlap
    inter_start = max(sv_start, elem_start)
    inter_end = min(sv_end, elem_end)
    if inter_end >= inter_start:
        overlap_bp = inter_end - inter_start + 1
        overlap_percent = (overlap_bp / (elem_end - elem_start + 1)) * 100
        overlaps_flag = 1
    else:
        overlap_bp = 0
        overlap_percent = 0
        overlaps_flag = 0

    # Classify coarse region_hit
    if overlaps_flag:
        region_hit = "overlaps"
    else:
        if abs(signed_dist) <= proximal_window:
            if signed_dist < 0:
                region_hit = "proximal_upstream"
            elif signed_dist > 0:
                region_hit = "proximal_downstream"
            else:
                region_hit = "proximal"
        else:
            region_hit = "distal"

    proximal_flag = int(region_hit.startswith("proximal"))
    distal_flag = int(region_hit == "distal")

    return {
        "overlap_bp": overlap_bp,
        "overlap_percent": overlap_percent,
        "overlaps_flag": overlaps_flag,
        "region_hit": region_hit,
        "proximal_flag": proximal_flag,
        "distal_flag": distal_flag,
        "hit_side": "span",
    }


# =============================================================================
# SV-GENE MAPPING (DEL, DUP, INS, BND)
# =============================================================================

def _map_span_to_genes(
    df: pd.DataFrame,
    genes: pd.DataFrame,
    genes_only: pd.DataFrame,
    svtype: str,
    is_lncrnas: bool,
    window: int,
    genes_only_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
    genes_features_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Internal function to map span-type SVs (DEL, DUP) to genes.
    """
    df = df.copy()
    svs = df[df["SVTYPE"] == svtype].copy()
    hits_col_name = "lncRNA_hits" if is_lncrnas else "gene_hits"

    gene_hits_col = []

    for idx, sv in svs.iterrows():
        chrom = sv["chrom"]
        sv_start = int(sv["pos"])

        if "END" in sv and not pd.isna(sv["END"]):
            sv_end = int(sv["END"])
        elif "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            sv_end = sv_start + abs(int(sv["SVLEN"]))
        else:
            sv_end = sv_start

        g_chr = (
            genes_only_by_chrom.get(chrom, _EMPTY_DF)
            if genes_only_by_chrom is not None
            else genes_only[genes_only["chrom"] == chrom]
        )
        if g_chr.empty:
            gene_hits_col.append([])
            continue
        g_chr = g_chr[
            (g_chr["end"] >= sv_start - window) &
            (g_chr["start"] <= sv_end + window)
        ]

        if g_chr.empty:
            gene_hits_col.append([])
            continue

        if genes_features_by_chrom is not None:
            feats_chr = genes_features_by_chrom.get(chrom, _EMPTY_DF)
            genes_features_chr = feats_chr[feats_chr["gene_name"].isin(g_chr["gene_name"])]
        else:
            genes_features_chr = genes[genes["gene_name"].isin(g_chr["gene_name"])]
        connections = []

        for _, g in g_chr.iterrows():
            gene_start = int(g["start"])
            gene_end = int(g["end"])
            strand = g.get("strand", "+")
            gene_name = g["gene_name"]
            gene_id = g.get("gene_id", None)

            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                sv_start, sv_end, gene_start, gene_end, strand
            )

            if abs(signed_dist) > window:
                continue

            gene_feats = genes_features_chr[genes_features_chr["gene_name"] == gene_name]

            region_info = classify_span_sv_gene_hit(
                sv_start=sv_start,
                sv_end=sv_end,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                signed_dist=signed_dist,
                gene_features=gene_feats,
            )

            gene_hit = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(gene_hit)

        gene_hits_col.append(connections)

    svs[hits_col_name] = gene_hits_col
    df.loc[svs.index, hits_col_name] = svs[hits_col_name]

    return df


def _map_ins_to_genes(
    df: pd.DataFrame,
    genes: pd.DataFrame,
    genes_only: pd.DataFrame,
    is_lncrnas: bool,
    window: int,
    genes_only_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
    genes_features_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Map INS SVs to genes (point-based).
    """
    df = df.copy()
    ins_df = df[df["SVTYPE"] == "INS"].copy()
    hits_col_name = "lncRNA_hits" if is_lncrnas else "gene_hits"

    gene_hits_col = []

    for idx, sv in ins_df.iterrows():
        chrom = sv["chrom"]
        ins_pos = int(sv["pos"])

        g_chr = (
            genes_only_by_chrom.get(chrom, _EMPTY_DF)
            if genes_only_by_chrom is not None
            else genes_only[genes_only["chrom"] == chrom]
        )
        if g_chr.empty:
            gene_hits_col.append([])
            continue
        g_chr = g_chr[
            (g_chr["start"] <= ins_pos + window) &
            (g_chr["end"] >= ins_pos - window)
        ]

        if g_chr.empty:
            gene_hits_col.append([])
            continue

        if genes_features_by_chrom is not None:
            feats_chr = genes_features_by_chrom.get(chrom, _EMPTY_DF)
            genes_features_chr = feats_chr[feats_chr["gene_name"].isin(g_chr["gene_name"])]
        else:
            genes_features_chr = genes[genes["gene_name"].isin(g_chr["gene_name"])]
        connections = []

        for _, g in g_chr.iterrows():
            gene_start = int(g["start"])
            gene_end = int(g["end"])
            strand = g.get("strand", "+")
            gene_name = g["gene_name"]
            gene_id = g.get("gene_id", None)

            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                ins_pos, ins_pos, gene_start, gene_end, strand
            )

            if abs(signed_dist) > window:
                continue

            gene_feats = genes_features_chr[genes_features_chr["gene_name"] == gene_name]

            region_info = classify_span_sv_gene_hit(
                sv_start=ins_pos,
                sv_end=ins_pos,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                signed_dist=signed_dist,
                gene_features=gene_feats,
            )
            region_info["hit_side"] = "point"

            gene_hit = {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "strand": strand,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(gene_hit)

        gene_hits_col.append(connections)

    ins_df[hits_col_name] = gene_hits_col
    df.loc[ins_df.index, hits_col_name] = ins_df[hits_col_name]

    return df


def _map_bnds_to_genes(
    df: pd.DataFrame,
    genes: pd.DataFrame,
    genes_only: pd.DataFrame,
    is_lncrnas: bool,
    window: int,
    genes_only_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
    genes_features_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Map BND SVs to genes (both breakpoints).
    """
    df = df.copy()
    bnds = df[df["SVTYPE"] == "BND"].copy()
    hits_col_name = "lncRNA_hits" if is_lncrnas else "gene_hits"

    gene_hits_col = []

    for idx, sv in bnds.iterrows():
        chrom1 = sv["chrom"]
        pos1 = int(sv["pos"])

        chrom2 = sv.get("bnd_remote_chrom", None)
        pos2 = sv.get("bnd_remote_pos", None)

        if pd.isna(chrom2) or chrom2 is None or pd.isna(pos2) or pos2 is None:
            chrom2 = None
            pos2 = None
        else:
            chrom2 = str(chrom2)
            pos2 = int(pos2)

        gene_hits = []

        def process_side(bk_chrom, bk_pos, side_label, bp_index, mate_chrom, mate_pos):
            nonlocal gene_hits

            if bk_chrom is None:
                return

            g_side = (
                genes_only_by_chrom.get(bk_chrom, _EMPTY_DF)
                if genes_only_by_chrom is not None
                else genes_only[genes_only["chrom"] == bk_chrom]
            )
            if g_side.empty:
                return
            g_side = g_side[
                (g_side["start"] <= bk_pos + window) &
                (g_side["end"] >= bk_pos - window)
            ]

            if g_side.empty:
                return

            if genes_features_by_chrom is not None:
                feats_side = genes_features_by_chrom.get(bk_chrom, _EMPTY_DF)
                genes_features_side = feats_side[feats_side["gene_name"].isin(g_side["gene_name"])]
            else:
                genes_features_side = genes[genes["gene_name"].isin(g_side["gene_name"])]

            for _, g in g_side.iterrows():
                gene_start = int(g["start"])
                gene_end = int(g["end"])
                strand = g.get("strand", "+")
                gene_name = g["gene_name"]
                gene_id = g.get("gene_id", None)

                signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                    bk_pos, bk_pos, gene_start, gene_end, strand
                )

                if abs(signed_dist) > window:
                    continue

                gene_feats = genes_features_side[genes_features_side["gene_name"] == gene_name]

                region_info = classify_span_sv_gene_hit(
                    sv_start=bk_pos,
                    sv_end=bk_pos,
                    gene_start=gene_start,
                    gene_end=gene_end,
                    strand=strand,
                    signed_dist=signed_dist,
                    gene_features=gene_feats,
                )

                region_info["hit_side"] = side_label
                region_info["bp_index"] = bp_index
                region_info["bp_chrom"] = bk_chrom
                region_info["bp_pos"] = bk_pos
                region_info["mate_chrom"] = mate_chrom
                region_info["mate_pos"] = mate_pos

                gene_hit = {
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "strand": strand,
                    "signed_dist": signed_dist,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    **region_info,
                }

                gene_hits.append(gene_hit)

        process_side(chrom1, pos1, "bp1", 1, chrom2, pos2)
        process_side(chrom2, pos2, "bp2", 2, chrom1, pos1)

        gene_hits_col.append(gene_hits)

    bnds[hits_col_name] = gene_hits_col
    df.loc[bnds.index, hits_col_name] = bnds[hits_col_name]

    return df


# =============================================================================
# SV-ELEMENT MAPPING
# =============================================================================

def _map_span_to_elements(
    df: pd.DataFrame,
    elements: pd.DataFrame,
    svtype: str,
    window: int,
    elements_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Map span-type SVs (DEL, DUP) to regulatory elements.
    """
    df = df.copy()
    svs = df[df["SVTYPE"] == svtype].copy()

    elem_hits_col = []

    for idx, sv in svs.iterrows():
        chrom = sv["chrom"]
        sv_start = int(sv["pos"])

        if "END" in sv and not pd.isna(sv["END"]):
            sv_end = int(sv["END"])
        elif "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            sv_end = sv_start + abs(int(sv["SVLEN"]))
        else:
            sv_end = sv_start

        e_chr = (
            elements_by_chrom.get(chrom, _EMPTY_DF)
            if elements_by_chrom is not None
            else elements[elements["chrom"] == chrom]
        )
        if e_chr.empty:
            elem_hits_col.append([])
            continue
        e_chr = e_chr[
            (e_chr["end"] >= sv_start - window) &
            (e_chr["start"] <= sv_end + window)
        ]

        if e_chr.empty:
            elem_hits_col.append([])
            continue

        connections = []

        for _, e in e_chr.iterrows():
            elem_start = int(e["start"])
            elem_end = int(e["end"])
            elem_id = e["elem_id"]
            elem_type = e.get("type", None)

            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                sv_start, sv_end, elem_start, elem_end, strand="+"
            )

            if abs(signed_dist) > window:
                continue

            region_info = classify_span_sv_elem_hit(
                sv_start=sv_start,
                sv_end=sv_end,
                elem_start=elem_start,
                elem_end=elem_end,
                signed_dist=signed_dist,
            )

            elem_hit = {
                "elem_id": elem_id,
                "elem_type": elem_type,
                "chrom": chrom,
                "elem_start": elem_start,
                "elem_end": elem_end,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(elem_hit)

        elem_hits_col.append(connections)

    svs["elem_hits"] = elem_hits_col
    df.loc[svs.index, "elem_hits"] = svs["elem_hits"]

    return df


def _map_ins_to_elements(
    df: pd.DataFrame,
    elements: pd.DataFrame,
    window: int,
    elements_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Map INS SVs to regulatory elements.
    """
    df = df.copy()
    ins_df = df[df["SVTYPE"] == "INS"].copy()

    elem_hits_col = []

    for idx, sv in ins_df.iterrows():
        chrom = sv["chrom"]
        ins_pos = int(sv["pos"])

        e_chr = (
            elements_by_chrom.get(chrom, _EMPTY_DF)
            if elements_by_chrom is not None
            else elements[elements["chrom"] == chrom]
        )
        if e_chr.empty:
            elem_hits_col.append([])
            continue
        e_chr = e_chr[
            (e_chr["start"] <= ins_pos + window) &
            (e_chr["end"] >= ins_pos - window)
        ]

        if e_chr.empty:
            elem_hits_col.append([])
            continue

        connections = []

        for _, e in e_chr.iterrows():
            elem_start = int(e["start"])
            elem_end = int(e["end"])
            elem_id = e["elem_id"]
            elem_type = e.get("type", None)

            signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                ins_pos, ins_pos, elem_start, elem_end, strand="+"
            )

            if abs(signed_dist) > window:
                continue

            region_info = classify_span_sv_elem_hit(
                sv_start=ins_pos,
                sv_end=ins_pos,
                elem_start=elem_start,
                elem_end=elem_end,
                signed_dist=signed_dist,
            )
            region_info["hit_side"] = "point"

            elem_hit = {
                "elem_id": elem_id,
                "elem_type": elem_type,
                "chrom": chrom,
                "elem_start": elem_start,
                "elem_end": elem_end,
                "signed_dist": signed_dist,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                **region_info,
            }

            connections.append(elem_hit)

        elem_hits_col.append(connections)

    ins_df["elem_hits"] = elem_hits_col
    df.loc[ins_df.index, "elem_hits"] = ins_df["elem_hits"]

    return df


def _map_bnds_to_elements(
    df: pd.DataFrame,
    elements: pd.DataFrame,
    window: int,
    elements_by_chrom: Optional[Dict[Any, pd.DataFrame]] = None,
) -> pd.DataFrame:
    """
    Map BND SVs to regulatory elements.
    """
    df = df.copy()
    bnds = df[df["SVTYPE"] == "BND"].copy()

    elem_hits_col = []

    for idx, sv in bnds.iterrows():
        chrom1 = sv["chrom"]
        pos1 = int(sv["pos"])

        chrom2 = sv.get("bnd_remote_chrom", None)
        pos2 = sv.get("bnd_remote_pos", None)

        if pd.isna(chrom2) or chrom2 is None or pd.isna(pos2) or pos2 is None:
            chrom2 = None
            pos2 = None
        else:
            chrom2 = str(chrom2)
            pos2 = int(pos2)

        elem_hits = []

        def process_side(bk_chrom, bk_pos, side_label, bp_index, mate_chrom, mate_pos):
            nonlocal elem_hits

            if bk_chrom is None:
                return

            e_side = (
                elements_by_chrom.get(bk_chrom, _EMPTY_DF)
                if elements_by_chrom is not None
                else elements[elements["chrom"] == bk_chrom]
            )
            if e_side.empty:
                return
            e_side = e_side[
                (e_side["start"] <= bk_pos + window) &
                (e_side["end"] >= bk_pos - window)
            ]

            if e_side.empty:
                return

            for _, e in e_side.iterrows():
                elem_start = int(e["start"])
                elem_end = int(e["end"])
                elem_id = e["elem_id"]
                elem_type = e.get("type", None)

                signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
                    bk_pos, bk_pos, elem_start, elem_end, strand="+"
                )

                if abs(signed_dist) > window:
                    continue

                region_info = classify_span_sv_elem_hit(
                    sv_start=bk_pos,
                    sv_end=bk_pos,
                    elem_start=elem_start,
                    elem_end=elem_end,
                    signed_dist=signed_dist,
                )

                region_info["hit_side"] = side_label
                region_info["bp_index"] = bp_index
                region_info["bp_chrom"] = bk_chrom
                region_info["bp_pos"] = bk_pos
                region_info["mate_chrom"] = mate_chrom
                region_info["mate_pos"] = mate_pos

                elem_hit = {
                    "elem_id": elem_id,
                    "elem_type": elem_type,
                    "chrom": bk_chrom,
                    "elem_start": elem_start,
                    "elem_end": elem_end,
                    "signed_dist": signed_dist,
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    **region_info,
                }

                elem_hits.append(elem_hit)

        process_side(chrom1, pos1, "bp1", 1, chrom2, pos2)
        process_side(chrom2, pos2, "bp2", 2, chrom1, pos1)

        elem_hits_col.append(elem_hits)

    bnds["elem_hits"] = elem_hits_col
    df.loc[bnds.index, "elem_hits"] = bnds["elem_hits"]

    return df


# =============================================================================
# PUBLIC API
# =============================================================================

def map_svs_to_genes(
    sv_df: pd.DataFrame,
    genes: pd.DataFrame,
    genes_only: pd.DataFrame,
    window: int = 1_000_000,
) -> pd.DataFrame:
    """
    Map all SVs to genes.
    
    Args:
        sv_df: SV DataFrame
        genes: Full gene features DataFrame
        genes_only: Gene-level only DataFrame (feature == "gene")
        window: Maximum distance to consider
    
    Returns:
        DataFrame with gene_hits column added
    """
    sv_df = sv_df.copy()

    if "gene_hits" not in sv_df.columns:
        sv_df["gene_hits"] = [[] for _ in range(len(sv_df))]

    sv_df = add_bnd_remote_coords(sv_df)

    genes_only_by_chrom = _chrom_index(genes_only)
    genes_features_by_chrom = _chrom_index(genes)

    sv_df = _map_span_to_genes(
        sv_df, genes, genes_only, "DEL", False, window,
        genes_only_by_chrom=genes_only_by_chrom,
        genes_features_by_chrom=genes_features_by_chrom,
    )
    sv_df = _map_span_to_genes(
        sv_df, genes, genes_only, "DUP", False, window,
        genes_only_by_chrom=genes_only_by_chrom,
        genes_features_by_chrom=genes_features_by_chrom,
    )
    sv_df = _map_ins_to_genes(
        sv_df, genes, genes_only, False, window,
        genes_only_by_chrom=genes_only_by_chrom,
        genes_features_by_chrom=genes_features_by_chrom,
    )
    sv_df = _map_bnds_to_genes(
        sv_df, genes, genes_only, False, window,
        genes_only_by_chrom=genes_only_by_chrom,
        genes_features_by_chrom=genes_features_by_chrom,
    )

    return sv_df


def map_svs_to_elements(
    sv_df: pd.DataFrame,
    elements: pd.DataFrame,
    window: int = 1_000_000,
) -> pd.DataFrame:
    """
    Map all SVs to regulatory elements.
    
    Args:
        sv_df: SV DataFrame
        elements: Regulatory elements DataFrame with elem_id column
        window: Maximum distance to consider
    
    Returns:
        DataFrame with elem_hits column added
    """
    sv_df = sv_df.copy()

    if "elem_hits" not in sv_df.columns:
        sv_df["elem_hits"] = [[] for _ in range(len(sv_df))]

    elements_by_chrom = _chrom_index(elements)

    sv_df = _map_span_to_elements(sv_df, elements, "DEL", window, elements_by_chrom=elements_by_chrom)
    sv_df = _map_span_to_elements(sv_df, elements, "DUP", window, elements_by_chrom=elements_by_chrom)
    sv_df = _map_ins_to_elements(sv_df, elements, window, elements_by_chrom=elements_by_chrom)
    sv_df = _map_bnds_to_elements(sv_df, elements, window, elements_by_chrom=elements_by_chrom)

    return sv_df


def map_svs_to_lncrnas(
    sv_df: pd.DataFrame,
    lncrnas: pd.DataFrame,
    lncrnas_only: pd.DataFrame,
    window: int = 1_000_000,
) -> pd.DataFrame:
    """
    Map all SVs to lncRNAs.
    
    Args:
        sv_df: SV DataFrame
        lncrnas: Full lncRNA features DataFrame
        lncrnas_only: lncRNA-level only DataFrame
        window: Maximum distance to consider
    
    Returns:
        DataFrame with lncRNA_hits column added
    """
    sv_df = sv_df.copy()

    if "lncRNA_hits" not in sv_df.columns:
        sv_df["lncRNA_hits"] = [[] for _ in range(len(sv_df))]

    sv_df = add_bnd_remote_coords(sv_df)

    lncrnas_only_by_chrom = _chrom_index(lncrnas_only)
    lncrnas_features_by_chrom = _chrom_index(lncrnas)

    sv_df = _map_span_to_genes(
        sv_df, lncrnas, lncrnas_only, "DEL", True, window,
        genes_only_by_chrom=lncrnas_only_by_chrom,
        genes_features_by_chrom=lncrnas_features_by_chrom,
    )
    sv_df = _map_span_to_genes(
        sv_df, lncrnas, lncrnas_only, "DUP", True, window,
        genes_only_by_chrom=lncrnas_only_by_chrom,
        genes_features_by_chrom=lncrnas_features_by_chrom,
    )
    sv_df = _map_ins_to_genes(
        sv_df, lncrnas, lncrnas_only, True, window,
        genes_only_by_chrom=lncrnas_only_by_chrom,
        genes_features_by_chrom=lncrnas_features_by_chrom,
    )
    sv_df = _map_bnds_to_genes(
        sv_df, lncrnas, lncrnas_only, True, window,
        genes_only_by_chrom=lncrnas_only_by_chrom,
        genes_features_by_chrom=lncrnas_features_by_chrom,
    )

    return sv_df


def map_svs_to_mirnas(
    sv_df: pd.DataFrame,
    mirnas: pd.DataFrame,
    window: int = 1_000_000,
) -> pd.DataFrame:
    """
    Map all SVs to miRNA loci (pre-miRNA coordinates), carrying mature-arm metadata when present.

    Adds column: `mir_hits` (list[dict]) similar to CNV mir_hits.
    """
    sv_df = sv_df.copy()
    if "mir_hits" not in sv_df.columns:
        sv_df["mir_hits"] = [[] for _ in range(len(sv_df))]

    sv_df = add_bnd_remote_coords(sv_df)
    sv_df = _map_span_to_mirnas(sv_df, mirnas, "DEL", window)
    sv_df = _map_span_to_mirnas(sv_df, mirnas, "DUP", window)
    # Other span types may exist (INV); map if present.
    sv_df = _map_span_to_mirnas(sv_df, mirnas, "INV", window)
    sv_df = _map_ins_to_mirnas(sv_df, mirnas, window)
    sv_df = _map_bnds_to_mirnas(sv_df, mirnas, window)
    return sv_df
