"""
Annotate feature DataFrames with ATAC peak links.

Supports:
- Interval features (genes, cCREs, lncRNAs): uses start/end
- SV features: handles DEL/DUP as intervals, INS/BND as points
- SNV features: uses pos as point

Adds atac_peak_links column with list of dicts containing:
- Peak info: peak_id, chrom, start, end, center, score, percentGC, annotation
- Spatial info: signed_dist, overlap_bp, upstream_5kb_flag, downstream_5kb_flag
"""

import bisect
from typing import List, Dict, Any, Optional, Union, Literal

import numpy as np
import pandas as pd

from ..utils import normalize_chrom, compute_interval_overlap


# =============================================================================
# FEATURE TYPE DEFINITIONS
# =============================================================================

FeatureKind = Literal["interval", "gene", "ccre", "lncrna", "sv", "snv", "point"]


# =============================================================================
# FEATURE INTERVAL EXTRACTION
# =============================================================================

def _extract_feature_intervals(
    row: pd.Series,
    kind: FeatureKind,
) -> List[Dict[str, Any]]:
    """
    Extract interval(s) from a feature row based on feature type.
    
    Returns list of dicts with:
        - chrom: str
        - start: int
        - end: int
        - center: int (midpoint for distance calculations)
        - label: str (feature type)
    
    For point features (SNV, INS, BND), start == end == center.
    """
    chrom = row.get("chrom")
    if chrom is None or (isinstance(chrom, float) and pd.isna(chrom)):
        return []
    
    intervals = []
    
    # Standard interval features
    if kind in {"interval", "gene", "ccre", "lncrna"}:
        start = row.get("start")
        end = row.get("end")
        
        if pd.isna(start) or pd.isna(end):
            return []
        
        start, end = int(start), int(end)
        if end < start:
            start, end = end, start
        
        center = (start + end) // 2
        intervals.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "center": center,
            "label": kind,
        })
        return intervals
    
    # SNV: point feature at pos
    if kind == "snv":
        pos = row.get("pos")
        if pd.isna(pos):
            return []
        
        pos = int(pos)
        intervals.append({
            "chrom": chrom,
            "start": pos,
            "end": pos,
            "center": pos,
            "label": "snv",
        })
        return intervals
    
    # Point feature (generic)
    if kind == "point":
        pos = row.get("pos", row.get("position", row.get("start")))
        if pd.isna(pos):
            return []
        
        pos = int(pos)
        intervals.append({
            "chrom": chrom,
            "start": pos,
            "end": pos,
            "center": pos,
            "label": "point",
        })
        return intervals
    
    # SV: depends on SVTYPE
    if kind == "sv":
        svtype = row.get("SVTYPE", "")
        pos = row.get("pos")
        end = row.get("END")
        svlen = row.get("SVLEN")
        
        if pd.isna(pos):
            return []
        
        pos = int(pos)
        
        # DEL/DUP: proper interval
        if svtype in {"DEL", "DUP"}:
            if pd.notna(end):
                end = int(end)
            elif pd.notna(svlen):
                end = pos + abs(int(svlen))
            else:
                end = pos
            
            if end < pos:
                pos, end = end, pos
            
            center = (pos + end) // 2
            intervals.append({
                "chrom": chrom,
                "start": pos,
                "end": end,
                "center": center,
                "label": f"sv_{svtype}",
            })
            return intervals
        
        # INS/BND/other: point at pos
        intervals.append({
            "chrom": chrom,
            "start": pos,
            "end": pos,
            "center": pos,
            "label": f"sv_{svtype}" if svtype else "sv_unknown",
        })
        return intervals
    
    raise ValueError(f"Unknown feature kind: {kind}")


# =============================================================================
# PEAK LINK COMPUTATION
# =============================================================================

def _compute_peak_link(
    feature_start: int,
    feature_end: int,
    feature_center: int,
    peak_start: int,
    peak_end: int,
    peak_center: int,
    peak_row: pd.Series,
    upstream_threshold: int = 5000,
    downstream_threshold: int = 5000,
    strand: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compute link between a feature and a peak.
    
    Args:
        feature_start, feature_end, feature_center: Feature coordinates
        peak_start, peak_end, peak_center: Peak coordinates
        peak_row: Peak DataFrame row (for extracting metadata)
        upstream_threshold: Distance threshold for upstream flag (bp)
        downstream_threshold: Distance threshold for downstream flag (bp)
        strand: Feature strand (for directional distance, None = unsigned)
    
    Returns:
        Dict with peak info and spatial relationship
    """
    # Compute overlap
    overlap = compute_interval_overlap(feature_start, feature_end, peak_start, peak_end)
    
    # Compute signed distance from feature center to peak center
    # Positive = peak is downstream (right), Negative = peak is upstream (left)
    raw_dist = peak_center - feature_center
    
    # If strand is negative, flip the sign convention
    if strand == "-":
        signed_dist = -raw_dist
    else:
        signed_dist = raw_dist
    
    # Absolute distance
    abs_dist = abs(signed_dist)
    
    # Upstream/downstream flags
    # Upstream: peak is before the feature (considering strand)
    # Downstream: peak is after the feature (considering strand)
    if overlap["overlaps"]:
        # Overlapping peaks are neither strictly upstream nor downstream
        upstream_flag = False
        downstream_flag = False
    else:
        if strand == "-":
            # For negative strand, upstream is positive direction
            upstream_flag = raw_dist > 0 and raw_dist <= upstream_threshold
            downstream_flag = raw_dist < 0 and abs(raw_dist) <= downstream_threshold
        else:
            # For positive strand (or no strand), upstream is negative direction
            upstream_flag = raw_dist < 0 and abs(raw_dist) <= upstream_threshold
            downstream_flag = raw_dist > 0 and raw_dist <= downstream_threshold
    
    return {
        # Peak identifiers and metadata
        "peak_id": peak_row.get("peak_id"),
        "chrom": peak_row.get("chrom"),
        "start": int(peak_start),
        "end": int(peak_end),
        "center": int(peak_center),
        "score": peak_row.get("score"),
        "percentGC": peak_row.get("percentGC"),
        "annotation": peak_row.get("annotation"),
        # Spatial relationship
        "signed_dist": int(signed_dist),
        "abs_dist": int(abs_dist),
        "overlap_bp": overlap["overlap_bp"],
        "overlaps": overlap["overlaps"],
        "upstream_5kb_flag": upstream_flag,
        "downstream_5kb_flag": downstream_flag,
    }


# =============================================================================
# MAIN ANNOTATION FUNCTION
# =============================================================================

def annotate_df_with_peaks(
    df: pd.DataFrame,
    peaks: pd.DataFrame,
    kind: FeatureKind,
    window_bp: int = 100_000,
    out_col: str = "atac_peak_links",
    strand_col: Optional[str] = "strand",
    upstream_threshold: int = 5000,
    downstream_threshold: int = 5000,
    bin_size: int = 100_000,
    copy_df: bool = True,
) -> pd.DataFrame:
    """
    Annotate a feature DataFrame with nearby ATAC peaks.
    
    Args:
        df: Feature DataFrame (genes, cCREs, SVs, SNVs, etc.)
        peaks: ATAC peaks DataFrame with peak_id, chrom, start, end, center
        kind: Feature type - "interval", "gene", "ccre", "lncrna", "sv", "snv", "point"
        window_bp: Maximum distance from feature center to peak center
        out_col: Output column name for peak links
        strand_col: Column name for strand info (for directional distance)
        upstream_threshold: Distance threshold for upstream_5kb_flag
        downstream_threshold: Distance threshold for downstream_5kb_flag
        bin_size: Bin size for coarse prefiltering
        copy_df: If True (default), work on a copy of ``df``. Set False only when the caller
            owns the frame and no longer needs the pre-annotation object identity (saves
            RAM for very large tables such as ``elem_focus``).
    
    Returns:
        DataFrame with atac_peak_links column added
    
    Each row's atac_peak_links is a list of dicts sorted by abs_dist:
        [
            {
                "peak_id": "chr6:29000000-29001000",
                "chrom": "chr6",
                "start": 29000000,
                "end": 29001000,
                "center": 29000500,
                "score": 5.2,
                "percentGC": 0.45,
                "annotation": "Promoter",
                "signed_dist": -5000,
                "abs_dist": 5000,
                "overlap_bp": 0,
                "overlaps": False,
                "upstream_5kb_flag": True,
                "downstream_5kb_flag": False,
            },
            ...
        ]
    """
    if copy_df:
        df = df.copy()
    
    # Initialize output column
    df[out_col] = [[] for _ in range(len(df))]
    
    # Normalize chromosomes
    df["chrom"] = normalize_chrom(df["chrom"])
    peaks = peaks.copy()
    peaks["chrom"] = normalize_chrom(peaks["chrom"])
    
    # Ensure peaks have center column
    if "center" not in peaks.columns:
        peaks["center"] = ((peaks["start"] + peaks["end"]) // 2).astype(int)
    
    # Check for strand column
    has_strand = strand_col and strand_col in df.columns
    
    # Prepare peaks with bins for prefiltering
    peaks["start"] = peaks["start"].astype(int)
    peaks["end"] = peaks["end"].astype(int)
    peaks["center"] = peaks["center"].astype(int)
    
    # Process per chromosome
    for chrom in df["chrom"].dropna().unique():
        feature_mask = df["chrom"] == chrom
        if not feature_mask.any():
            continue
        
        peaks_chrom = peaks[peaks["chrom"] == chrom]
        if peaks_chrom.empty:
            continue
        
        # Build peak lookup sorted by center for window queries (avoids O(n_peaks) per feature)
        peak_data = [
            (int(row["start"]), int(row["end"]), int(row["center"]), row)
            for _, row in peaks_chrom.iterrows()
        ]
        peak_data.sort(key=lambda t: t[2])
        peak_centers = [t[2] for t in peak_data]

        # Process each feature on this chromosome
        for idx in df.index[feature_mask]:
            row = df.loc[idx]
            
            # Get strand if available
            strand = row.get(strand_col) if has_strand else None
            if pd.isna(strand):
                strand = None
            
            # Extract feature intervals
            intervals = _extract_feature_intervals(row, kind)
            if not intervals:
                continue
            
            # For each interval, find nearby peaks
            all_links = []
            for itv in intervals:
                f_start = itv["start"]
                f_end = itv["end"]
                f_center = itv["center"]

                lo_c = f_center - window_bp
                hi_c = f_center + window_bp
                j0 = bisect.bisect_left(peak_centers, lo_c)
                j1 = bisect.bisect_right(peak_centers, hi_c)

                for j in range(j0, j1):
                    p_start, p_end, p_center, p_row = peak_data[j]

                    link = _compute_peak_link(
                        f_start, f_end, f_center,
                        p_start, p_end, p_center,
                        p_row,
                        upstream_threshold=upstream_threshold,
                        downstream_threshold=downstream_threshold,
                        strand=strand,
                    )
                    all_links.append(link)
            
            # Sort by absolute distance and deduplicate by peak_id
            seen_peaks = set()
            unique_links = []
            for link in sorted(all_links, key=lambda x: x["abs_dist"]):
                if link["peak_id"] not in seen_peaks:
                    seen_peaks.add(link["peak_id"])
                    unique_links.append(link)
            
            df.at[idx, out_col] = unique_links
    
    # Print summary
    n_with_peaks = df[out_col].apply(lambda x: len(x) > 0).sum()
    total_links = df[out_col].apply(len).sum()
    print(f"Annotated {n_with_peaks}/{len(df)} features with ATAC peaks")
    print(f"  Total peak links: {total_links}")
    print(f"  Mean peaks per feature: {total_links / len(df):.1f}")
    
    return df


# =============================================================================
# CONVENIENCE WRAPPERS
# =============================================================================

def annotate_genes_with_peaks(
    genes: pd.DataFrame,
    peaks: pd.DataFrame,
    window_bp: int = 100_000,
    **kwargs,
) -> pd.DataFrame:
    """Annotate genes with nearby ATAC peaks."""
    return annotate_df_with_peaks(
        genes, peaks,
        kind="gene",
        window_bp=window_bp,
        **kwargs,
    )


def annotate_ccres_with_peaks(
    ccres: pd.DataFrame,
    peaks: pd.DataFrame,
    window_bp: int = 100_000,
    **kwargs,
) -> pd.DataFrame:
    """Annotate cCREs with nearby ATAC peaks."""
    return annotate_df_with_peaks(
        ccres, peaks,
        kind="ccre",
        window_bp=window_bp,
        **kwargs,
    )


def annotate_svs_with_peaks(
    svs: pd.DataFrame,
    peaks: pd.DataFrame,
    window_bp: int = 100_000,
    **kwargs,
) -> pd.DataFrame:
    """
    Annotate SVs with nearby ATAC peaks.
    
    Handles:
    - DEL/DUP: uses full interval
    - INS/BND: uses pos as point
    """
    return annotate_df_with_peaks(
        svs, peaks,
        kind="sv",
        window_bp=window_bp,
        strand_col=None,  # SVs don't have strand
        **kwargs,
    )


def annotate_snvs_with_peaks(
    snvs: pd.DataFrame,
    peaks: pd.DataFrame,
    window_bp: int = 100_000,
    **kwargs,
) -> pd.DataFrame:
    """Annotate SNVs with nearby ATAC peaks."""
    return annotate_df_with_peaks(
        snvs, peaks,
        kind="snv",
        window_bp=window_bp,
        strand_col=None,  # SNVs don't have strand
        **kwargs,
    )


# =============================================================================
# QUERY HELPERS
# =============================================================================

def get_features_with_overlapping_peaks(
    df: pd.DataFrame,
    peak_col: str = "atac_peak_links",
) -> pd.DataFrame:
    """Get features that have at least one overlapping peak."""
    if peak_col not in df.columns:
        return df.iloc[0:0].copy()
    
    def has_overlap(links):
        if not isinstance(links, list):
            return False
        return any(link.get("overlaps", False) for link in links)
    
    mask = df[peak_col].apply(has_overlap)
    return df[mask].copy()


def get_features_with_upstream_peaks(
    df: pd.DataFrame,
    peak_col: str = "atac_peak_links",
) -> pd.DataFrame:
    """Get features that have at least one upstream peak (within threshold)."""
    if peak_col not in df.columns:
        return df.iloc[0:0].copy()
    
    def has_upstream(links):
        if not isinstance(links, list):
            return False
        return any(link.get("upstream_5kb_flag", False) for link in links)
    
    mask = df[peak_col].apply(has_upstream)
    return df[mask].copy()


def get_features_with_downstream_peaks(
    df: pd.DataFrame,
    peak_col: str = "atac_peak_links",
) -> pd.DataFrame:
    """Get features that have at least one downstream peak (within threshold)."""
    if peak_col not in df.columns:
        return df.iloc[0:0].copy()
    
    def has_downstream(links):
        if not isinstance(links, list):
            return False
        return any(link.get("downstream_5kb_flag", False) for link in links)
    
    mask = df[peak_col].apply(has_downstream)
    return df[mask].copy()


def get_closest_peak(
    links: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Get the closest peak from a list of peak links."""
    if not links:
        return None
    return min(links, key=lambda x: x.get("abs_dist", float("inf")))


def summarize_peak_links(
    df: pd.DataFrame,
    peak_col: str = "atac_peak_links",
) -> pd.DataFrame:
    """
    Create summary statistics for peak links per feature.
    
    Returns DataFrame with columns:
    - n_peaks: number of linked peaks
    - n_overlapping: number of overlapping peaks
    - n_upstream: peaks with upstream flag
    - n_downstream: peaks with downstream flag
    - min_dist: minimum absolute distance
    - closest_peak_id: ID of closest peak
    - closest_peak_score: score of closest peak
    """
    if peak_col not in df.columns:
        return pd.DataFrame()
    
    def summarize_links(links):
        if not isinstance(links, list) or len(links) == 0:
            return {
                "n_peaks": 0,
                "n_overlapping": 0,
                "n_upstream": 0,
                "n_downstream": 0,
                "min_dist": None,
                "closest_peak_id": None,
                "closest_peak_score": None,
            }
        
        n_overlapping = sum(1 for l in links if l.get("overlaps", False))
        n_upstream = sum(1 for l in links if l.get("upstream_5kb_flag", False))
        n_downstream = sum(1 for l in links if l.get("downstream_5kb_flag", False))
        
        closest = get_closest_peak(links)
        
        return {
            "n_peaks": len(links),
            "n_overlapping": n_overlapping,
            "n_upstream": n_upstream,
            "n_downstream": n_downstream,
            "min_dist": closest.get("abs_dist") if closest else None,
            "closest_peak_id": closest.get("peak_id") if closest else None,
            "closest_peak_score": closest.get("score") if closest else None,
        }
    
    summary = df[peak_col].apply(summarize_links).apply(pd.Series)
    return summary


def flatten_peak_links(
    df: pd.DataFrame,
    id_col: str,
    peak_col: str = "atac_peak_links",
) -> pd.DataFrame:
    """
    Flatten peak links into a long-form DataFrame.
    
    Args:
        df: DataFrame with peak links
        id_col: Column to use as feature identifier
        peak_col: Column containing peak links
    
    Returns:
        DataFrame with one row per feature-peak pair
    """
    if peak_col not in df.columns:
        return pd.DataFrame()
    
    rows = []
    for _, row in df.iterrows():
        feature_id = row[id_col]
        links = row[peak_col]
        
        if not isinstance(links, list):
            continue
        
        for link in links:
            rows.append({
                id_col: feature_id,
                **link,
            })
    
    return pd.DataFrame(rows)
