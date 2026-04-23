"""
ATAC peak to cCRE matching.

Provides:
- Overlap detection between peaks and cCREs
- Proximity/distance calculation for non-overlapping pairs
- Combined ccre_links list building
"""

from typing import List, Dict, Any, Optional

import numpy as np
import pandas as pd

from ..utils import compute_interval_overlap
from ..config import THRESHOLDS


# =============================================================================
# DISTANCE CALCULATIONS
# =============================================================================

def _min_interval_distance(
    a_start: np.ndarray,
    a_end: np.ndarray,
    b_start: np.ndarray,
    b_end: np.ndarray,
) -> np.ndarray:
    """
    Minimum genomic distance between intervals.
    Returns 0 if intervals overlap.
    """
    left_gap = b_start - a_end
    right_gap = a_start - b_end
    return np.maximum(0, np.maximum(left_gap, right_gap))


# =============================================================================
# CORE MATCHING FUNCTION
# =============================================================================

def match_peaks_to_ccres(
    peaks: pd.DataFrame,
    ccres: pd.DataFrame,
    max_distance: int = 0,
    bin_size: int = 100_000,
) -> pd.DataFrame:
    """
    Match ATAC peaks to cCREs by overlap and proximity.
    
    Args:
        peaks: ATAC peaks DataFrame with peak_id, chrom, start, end
        ccres: cCRE DataFrame with cCRE_id, ENCODE_id, chrom, start, end, raw_type
        max_distance: Maximum distance for proximity matching (0 = overlap only)
        bin_size: Bin size for coarse prefiltering
    
    Returns:
        DataFrame with all peak-cCRE pairs (overlapping or within max_distance)
    """
    # Prepare peaks
    p = peaks[["peak_id", "chrom", "start", "end"]].copy()
    p["start"] = p["start"].astype(np.int64)
    p["end"] = p["end"].astype(np.int64)
    p = p.rename(columns={"start": "peak_start", "end": "peak_end"})

    # Prepare cCREs
    ccre_cols = ["cCRE_id", "ENCODE_id", "chrom", "start", "end"]
    if "raw_type" in ccres.columns:
        ccre_cols.append("raw_type")
    
    c = ccres[ccre_cols].copy()
    c["start"] = c["start"].astype(np.int64)
    c["end"] = c["end"].astype(np.int64)
    c = c.rename(columns={"start": "ccre_start", "end": "ccre_end"})

    # Fast path: overlap-only matching avoids list-of-bins + explode (much lower RAM).
    if int(max_distance) == 0:
        p["bin_start"] = (p["peak_start"] // bin_size).astype(np.int64)
        p["bin_end"] = (p["peak_end"] // bin_size).astype(np.int64)

        pb = p[["peak_id", "chrom", "peak_start", "peak_end", "bin_start", "bin_end"]].copy()
        pb = pb.rename(columns={"bin_start": "bins"})
        pb = pb[["peak_id", "chrom", "peak_start", "peak_end", "bins"]]

        # Rare case: interval crosses bin boundary (duplicate into second bin)
        cross = p["bin_end"] != p["bin_start"]
        if cross.any():
            pb2 = p.loc[cross, ["peak_id", "chrom", "peak_start", "peak_end", "bin_end"]].copy()
            pb2 = pb2.rename(columns={"bin_end": "bins"})
            pb = pd.concat([pb, pb2], ignore_index=True)

        c["bin_start"] = (c["ccre_start"] // bin_size).astype(np.int64)
        c["bin_end"] = (c["ccre_end"] // bin_size).astype(np.int64)

        cb = c[["cCRE_id", "ENCODE_id", "chrom", "ccre_start", "ccre_end", "bin_start", "bin_end"]].copy()
        cb = cb.rename(columns={"bin_start": "bins"})
        cb = cb[["cCRE_id", "ENCODE_id", "chrom", "ccre_start", "ccre_end", "bins"]]

        cross_c = c["bin_end"] != c["bin_start"]
        if cross_c.any():
            cb2 = c.loc[cross_c, ["cCRE_id", "ENCODE_id", "chrom", "ccre_start", "ccre_end", "bin_end"]].copy()
            cb2 = cb2.rename(columns={"bin_end": "bins"})
            cb = pd.concat([cb, cb2], ignore_index=True)

        if "raw_type" in c.columns:
            raw = c[["cCRE_id", "raw_type"]].drop_duplicates("cCRE_id")
            cb = cb.merge(raw, on="cCRE_id", how="left")

        pref = pb.merge(
            cb,
            on=["chrom", "bins"],
            how="inner",
        )

        # Exact overlap filter (distance is always 0 in overlap-only mode)
        pref = pref[
            (pref["ccre_end"] >= pref["peak_start"])
            & (pref["ccre_start"] <= pref["peak_end"])
        ].copy()

        pref["distance"] = np.int64(0)
    else:
        # General path: overlap + proximity window
        p["win_start"] = (p["peak_start"] - max_distance).clip(lower=0).astype(np.int64)
        p["win_end"] = (p["peak_end"] + max_distance).astype(np.int64)

        p["bin_start"] = (p["win_start"] // bin_size).astype(np.int64)
        p["bin_end"] = (p["win_end"] // bin_size).astype(np.int64)
        p["bins"] = [list(range(s, e + 1)) for s, e in zip(p["bin_start"], p["bin_end"])]

        c["bin_start"] = (c["ccre_start"] // bin_size).astype(np.int64)
        c["bin_end"] = (c["ccre_end"] // bin_size).astype(np.int64)
        c["bins"] = [list(range(s, e + 1)) for s, e in zip(c["bin_start"], c["bin_end"])]

        pb = p.explode("bins")
        cb = c.explode("bins")

        ccre_merge_cols = ["cCRE_id", "ENCODE_id", "chrom", "ccre_start", "ccre_end", "bins"]
        if "raw_type" in c.columns:
            ccre_merge_cols.insert(2, "raw_type")

        pref = (
            pb[["peak_id", "chrom", "peak_start", "peak_end", "win_start", "win_end", "bins"]]
            .merge(
                cb[ccre_merge_cols],
                on=["chrom", "bins"],
                how="inner",
            )
        )

        pref = pref[
            (pref["ccre_end"] >= pref["win_start"])
            & (pref["ccre_start"] <= pref["win_end"])
        ].copy()

        pref["distance"] = _min_interval_distance(
            pref["peak_start"].to_numpy(),
            pref["peak_end"].to_numpy(),
            pref["ccre_start"].to_numpy(),
            pref["ccre_end"].to_numpy(),
        ).astype(np.int64)

        pref = pref[pref["distance"] <= max_distance].copy()
    
    # Compute overlaps
    def compute_overlap_row(row):
        ov = compute_interval_overlap(
            int(row["peak_start"]), int(row["peak_end"]),
            int(row["ccre_start"]), int(row["ccre_end"]),
        )
        # Rename for cCRE context
        return {
            "overlaps": ov["overlaps"],
            "overlap_bp": ov["overlap_bp"],
            "overlap_interval": ov["overlap_interval"],
            "overlap_frac_of_peak": ov["overlap_frac_of_a"],
            "overlap_frac_of_ccre": ov["overlap_frac_of_b"],
        }
    
    print(f"Computing overlaps for {len(pref)} peak-cCRE pairs...")
    pref["overlap"] = pref.apply(compute_overlap_row, axis=1)
    
    # Select final columns
    final_cols = [
        "peak_id", "cCRE_id", "ENCODE_id",
        "chrom", "peak_start", "peak_end", "ccre_start", "ccre_end",
        "distance", "overlap",
    ]
    if "raw_type" in pref.columns:
        final_cols.insert(3, "raw_type")
    
    pair_df = pref[final_cols].copy()
    
    # Drop duplicates from bin explosion
    pair_df = pair_df.drop_duplicates(
        subset=["peak_id", "cCRE_id"]
    ).reset_index(drop=True)
    
    n_overlapping = (pair_df["distance"] == 0).sum()
    print(f"Matched {len(pair_df)} peak-cCRE pairs within {max_distance}bp")
    print(f"  Overlapping: {n_overlapping}")
    print(f"  Proximal only: {len(pair_df) - n_overlapping}")
    
    return pair_df


# =============================================================================
# cCRE LINKS LIST BUILDING
# =============================================================================

def build_ccre_links(
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build ccre_links list column from peak-cCRE pairs.
    
    Args:
        pair_df: Output from match_peaks_to_ccres
    
    Returns:
        DataFrame with [peak_id, ccre_links] where ccre_links is list of dicts
    """
    def build_links_for_peak(df_sub: pd.DataFrame) -> List[Dict[str, Any]]:
        links = []
        # Sort by distance so overlapping cCREs come first
        df_sub = df_sub.sort_values("distance")
        
        for _, row in df_sub.iterrows():
            links.append({
                "cCRE_id": row["cCRE_id"],
                "ENCODE_id": row["ENCODE_id"],
                "raw_type": row.get("raw_type"),
                "distance": int(row["distance"]),
                "overlap": row["overlap"],
            })
        return links
    
    ccre_links_series = (
        pair_df
        .groupby("peak_id", sort=False)
        .apply(build_links_for_peak)
        .rename("ccre_links")
    )
    
    return ccre_links_series.reset_index()


# =============================================================================
# AGGREGATION HELPERS
# =============================================================================

def aggregate_ccres_per_peak(
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create peak-level summary of linked cCREs.
    """
    def agg_peak(df_sub):
        n_total = len(df_sub)
        n_overlapping = (df_sub["distance"] == 0).sum()
        min_dist = df_sub["distance"].min()
        
        # Count by raw_type
        if "raw_type" in df_sub.columns:
            type_counts = df_sub["raw_type"].value_counts().to_dict()
        else:
            type_counts = {}
        
        return pd.Series({
            "n_ccres_total": n_total,
            "n_ccres_overlapping": n_overlapping,
            "ccre_types": type_counts,
            "min_distance": int(min_dist),
        })
    
    return pair_df.groupby("peak_id").apply(agg_peak).reset_index()


def get_peaks_overlapping_ccre_type(
    pair_df: pd.DataFrame,
    ccre_type: str,
    require_overlap: bool = True,
) -> pd.DataFrame:
    """
    Get peaks that overlap or are near a specific cCRE type.
    """
    if "raw_type" not in pair_df.columns:
        raise ValueError("pair_df must have 'raw_type' column")
    
    # Handle raw_type which may contain multiple types
    mask = pair_df["raw_type"].str.contains(ccre_type, na=False, regex=False)
    
    if require_overlap:
        mask &= pair_df["distance"] == 0
    
    return pair_df[mask].sort_values(["peak_id", "distance"]).reset_index(drop=True)


def summarize_ccre_overlap_by_type(
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarize peak-cCRE overlaps by cCRE type.
    """
    def extract_primary_type(raw_type):
        if pd.isna(raw_type):
            return "unknown"
        return str(raw_type).split(",")[0].strip()
    
    df = pair_df.copy()
    df["primary_type"] = df["raw_type"].apply(extract_primary_type)
    
    def agg_type(df_sub):
        overlapping = df_sub[df_sub["distance"] == 0]
        mean_overlap = overlapping["overlap"].apply(
            lambda x: x.get("overlap_bp", 0) if isinstance(x, dict) else 0
        ).mean() if len(overlapping) > 0 else 0
        
        return pd.Series({
            "n_peaks": df_sub["peak_id"].nunique(),
            "n_pairs": len(df_sub),
            "n_overlapping": len(overlapping),
            "mean_overlap_bp": round(mean_overlap, 1),
        })
    
    return df.groupby("primary_type").apply(agg_type).reset_index()


# =============================================================================
# BIDIRECTIONAL LOOKUP HELPERS
# =============================================================================

def get_ccres_for_peak(
    pair_df: pd.DataFrame,
    peak_id: str,
    max_distance: Optional[int] = None,
    require_overlap: bool = False,
) -> pd.DataFrame:
    """Get all cCREs linked to a specific peak."""
    df = pair_df[pair_df["peak_id"] == peak_id].copy()
    
    if max_distance is not None:
        df = df[df["distance"] <= max_distance]
    
    if require_overlap:
        df = df[df["distance"] == 0]
    
    return df.sort_values("distance").reset_index(drop=True)


def get_peaks_for_ccre(
    pair_df: pd.DataFrame,
    ccre_id: str,
    id_col: str = "cCRE_id",
    max_distance: Optional[int] = None,
    require_overlap: bool = False,
) -> pd.DataFrame:
    """Get all peaks linked to a specific cCRE."""
    df = pair_df[pair_df[id_col] == ccre_id].copy()
    
    if max_distance is not None:
        df = df[df["distance"] <= max_distance]
    
    if require_overlap:
        df = df[df["distance"] == 0]
    
    return df.sort_values("distance").reset_index(drop=True)
