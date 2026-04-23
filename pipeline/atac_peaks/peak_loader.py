"""
ATAC peak loading and normalization.

Functions for:
- Loading ATAC peaks from CSV/TSV
- Generating composite peak IDs
- Creating peak_id ↔ original_name mapping table
"""

from pathlib import Path
from typing import Union, Optional, Tuple, List

import pandas as pd
import numpy as np

from ..utils import normalize_chrom


# =============================================================================
# PEAK ID GENERATION
# =============================================================================

def generate_peak_id(chrom: str, start: int, end: int) -> str:
    """
    Generate composite peak ID from coordinates.
    
    Format: {chrom}:{start}-{end}
    Example: chr19:6012-90651
    """
    return f"{chrom}:{start}-{end}"


def generate_peak_ids(df: pd.DataFrame) -> pd.Series:
    """
    Generate peak IDs for entire DataFrame.
    
    Args:
        df: DataFrame with chrom, start, end columns
    
    Returns:
        Series of peak IDs
    """
    return (
        df["chrom"].astype(str) + ":" +
        df["start"].astype(str) + "-" +
        df["end"].astype(str)
    )


# =============================================================================
# LOADING
# =============================================================================

def load_atac_peaks(
    path: Union[str, Path],
    sep: Optional[str] = None,
    chrom_col: str = "seqnames",
    start_col: str = "start",
    end_col: str = "end",
    name_col: str = "name",
) -> pd.DataFrame:
    """
    Load ATAC peaks from CSV/TSV file.
    
    Args:
        path: Path to peaks file
        sep: Separator (auto-detected if None)
        chrom_col: Column name for chromosome
        start_col: Column name for start position
        end_col: Column name for end position
        name_col: Column name for original peak name
    
    Returns:
        DataFrame with normalized columns:
        - peak_id: composite ID (chr:start-end)
        - chrom, start, end: coordinates
        - original_name: original peak name from file
        - score, annotation, percentGC: if present in input
    
    Example:
        >>> peaks = load_atac_peaks("peaks.csv")
        >>> peaks.columns
        ['peak_id', 'chrom', 'start', 'end', 'original_name', 'score', 'annotation', 'percentGC']
    """
    path = Path(path)
    
    # Auto-detect separator
    if sep is None:
        with open(path, 'r') as f:
            first_line = f.readline()
        sep = '\t' if '\t' in first_line else ','
    
    df = pd.read_csv(path, sep=sep)
    
    # Rename core columns
    rename_map = {}
    if chrom_col in df.columns and chrom_col != "chrom":
        rename_map[chrom_col] = "chrom"
    if name_col in df.columns and name_col != "original_name":
        rename_map[name_col] = "original_name"
    
    if rename_map:
        df = df.rename(columns=rename_map)
    
    # Normalize chromosome
    df["chrom"] = normalize_chrom(df["chrom"])
    
    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    # Sort by coordinates to ensure consistent ordering
    df = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    
    # Generate peak IDs
    df["peak_id"] = generate_peak_ids(df)
    
    # Check for duplicate peak IDs (shouldn't happen with coordinate-based IDs)
    n_dup = df["peak_id"].duplicated().sum()
    if n_dup > 0:
        print(f"Warning: {n_dup} duplicate peak IDs found (same coordinates)")
    
    # Reorder columns
    core_cols = ["peak_id", "chrom", "start", "end", "original_name"]
    other_cols = [c for c in df.columns if c not in core_cols]
    df = df[core_cols + other_cols]
    
    # Compute center position
    df["center"] = ((df["start"] + df["end"]) // 2).astype("Int64")
    
    # Compute peak length
    df["length"] = (df["end"] - df["start"]).astype("Int64")
    
    print(f"Loaded {len(df)} ATAC peaks from {path}")
    print(f"  Chromosomes: {df['chrom'].nunique()}")
    print(f"  Mean peak length: {df['length'].mean():.0f} bp")
    
    return df


# =============================================================================
# MAPPING TABLE
# =============================================================================

def create_peak_id_mapping(
    peaks: pd.DataFrame,
    output_path: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Create mapping table between peak_id and original_name.
    
    Args:
        peaks: ATAC peaks DataFrame
        output_path: Optional path to save mapping CSV
    
    Returns:
        DataFrame with [peak_id, original_name, chrom, start, end]
    """
    mapping = peaks[["peak_id", "original_name", "chrom", "start", "end"]].copy()
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        mapping.to_csv(output_path, index=False)
        print(f"Saved peak ID mapping to {output_path}")
    
    return mapping


# =============================================================================
# FILTERING HELPERS
# =============================================================================

def filter_peaks_by_chromosomes(
    peaks: pd.DataFrame,
    include_chroms: Optional[List[str]] = None,
    exclude_chroms: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Filter peaks by chromosome.
    
    Args:
        peaks: ATAC peaks DataFrame
        include_chroms: List of chromosomes to keep (None = keep all)
        exclude_chroms: List of chromosomes to remove
    
    Returns:
        Filtered DataFrame
    """
    df = peaks.copy()
    
    if include_chroms:
        include_set = {normalize_chrom(c) for c in include_chroms}
        df = df[df["chrom"].isin(include_set)]
    
    if exclude_chroms:
        exclude_set = {normalize_chrom(c) for c in exclude_chroms}
        df = df[~df["chrom"].isin(exclude_set)]
    
    return df


def filter_peaks_by_score(
    peaks: pd.DataFrame,
    min_score: Optional[float] = None,
    max_score: Optional[float] = None,
    score_col: str = "score",
) -> pd.DataFrame:
    """
    Filter peaks by score threshold.
    
    Args:
        peaks: ATAC peaks DataFrame
        min_score: Minimum score (inclusive)
        max_score: Maximum score (inclusive)
        score_col: Column name for score
    
    Returns:
        Filtered DataFrame
    """
    if score_col not in peaks.columns:
        raise ValueError(f"Score column '{score_col}' not found in DataFrame")
    
    df = peaks.copy()
    
    if min_score is not None:
        df = df[df[score_col] >= min_score]
    
    if max_score is not None:
        df = df[df[score_col] <= max_score]
    
    return df


def filter_peaks_by_length(
    peaks: pd.DataFrame,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
) -> pd.DataFrame:
    """
    Filter peaks by length.
    
    Args:
        peaks: ATAC peaks DataFrame
        min_length: Minimum peak length in bp
        max_length: Maximum peak length in bp
    
    Returns:
        Filtered DataFrame
    """
    df = peaks.copy()
    
    if "length" not in df.columns:
        df["length"] = df["end"] - df["start"]
    
    if min_length is not None:
        df = df[df["length"] >= min_length]
    
    if max_length is not None:
        df = df[df["length"] <= max_length]
    
    return df


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def summarize_peaks(peaks: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for peaks per chromosome.
    
    Returns DataFrame with per-chromosome stats:
    - n_peaks: number of peaks
    - total_bp: total covered bp
    - mean_length: mean peak length
    - median_length: median peak length
    - mean_score: mean score (if present)
    """
    df = peaks.copy()
    
    if "length" not in df.columns:
        df["length"] = df["end"] - df["start"]
    
    agg_dict = {
        "peak_id": "count",
        "length": ["sum", "mean", "median"],
    }
    
    if "score" in df.columns:
        agg_dict["score"] = "mean"
    
    summary = df.groupby("chrom").agg(agg_dict)
    summary.columns = ["_".join(col).strip("_") for col in summary.columns]
    summary = summary.rename(columns={
        "peak_id_count": "n_peaks",
        "length_sum": "total_bp",
        "length_mean": "mean_length",
        "length_median": "median_length",
        "score_mean": "mean_score",
    })
    
    return summary.reset_index()
