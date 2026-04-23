"""
Shared utility functions used across the pipeline.

Includes:
- Chromosome normalization
- Distance calculations
- Binning helpers
- Parsing helpers
- Type conversions
"""

import ast
import re
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd


# =============================================================================
# CHROMOSOME HANDLING
# =============================================================================

def normalize_chrom(s: pd.Series) -> pd.Series:
    """
    Normalize chromosome names to 'chr' prefix format.
    
    Examples:
        '1' -> 'chr1'
        'chr1' -> 'chr1'
        1 -> 'chr1'
    """
    s = s.astype(str)
    # Add 'chr' prefix if not present
    return s.str.replace(r"^(?!chr)", "chr", regex=True)


def harmonize_chrom_column(df: pd.DataFrame, chrom_col: str = "chrom") -> Tuple[pd.DataFrame, bool]:
    """
    Harmonize chromosome column in a DataFrame.
    
    Returns:
        Tuple of (modified DataFrame, whether changes were made)
    """
    df = df.copy()
    was_changed = False
    
    # Handle 'seqname' -> 'chrom' rename
    if "seqname" in df.columns and chrom_col not in df.columns:
        df.rename(columns={"seqname": chrom_col}, inplace=True)
        was_changed = True
    
    if chrom_col in df.columns:
        original = df[chrom_col].copy()
        df[chrom_col] = normalize_chrom(df[chrom_col])
        if not df[chrom_col].equals(original):
            was_changed = True
    
    return df, was_changed


# =============================================================================
# DISTANCE CALCULATIONS
# =============================================================================

def min_interval_distance(
    a_start: Union[int, np.ndarray],
    a_end: Union[int, np.ndarray],
    b_start: Union[int, np.ndarray],
    b_end: Union[int, np.ndarray],
) -> Union[int, np.ndarray]:
    """
    Minimum genomic distance between intervals [a_start, a_end] and [b_start, b_end].
    Returns 0 if intervals overlap.
    
    Works with scalars or numpy arrays.
    """
    left_gap = b_start - a_end
    right_gap = a_start - b_end
    return np.maximum(0, np.maximum(left_gap, right_gap))


def tss_to_interval_distance(
    tss: Union[int, np.ndarray],
    interval_start: Union[int, np.ndarray],
    interval_end: Union[int, np.ndarray],
) -> Union[int, np.ndarray]:
    """
    Distance from a TSS point to an interval.
    Returns 0 if TSS is inside the interval.
    """
    inside = (tss >= interval_start) & (tss <= interval_end)
    dist_to_start = np.abs(tss - interval_start)
    dist_to_end = np.abs(tss - interval_end)
    return np.where(inside, 0, np.minimum(dist_to_start, dist_to_end))


# =============================================================================
# GENERAL INTERVAL OVERLAP (shared utility)
# =============================================================================

def compute_interval_overlap(
    a_start: int,
    a_end: int,
    b_start: int,
    b_end: int,
) -> Dict[str, Any]:
    """
    Compute overlap between two intervals.
    
    Args:
        a_start, a_end: First interval (e.g., peak)
        b_start, b_end: Second interval (e.g., gene body, cCRE)
    
    Returns:
        Dict with overlap statistics:
        - overlaps: bool
        - overlap_bp: int
        - overlap_interval: [start, end] or None
        - overlap_frac_of_a: float
        - overlap_frac_of_b: float
    """
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    
    if overlap_start >= overlap_end:
        return {
            "overlaps": False,
            "overlap_bp": 0,
            "overlap_interval": None,
            "overlap_frac_of_a": 0.0,
            "overlap_frac_of_b": 0.0,
        }
    
    overlap_bp = overlap_end - overlap_start
    a_len = max(1, a_end - a_start)
    b_len = max(1, b_end - b_start)
    
    return {
        "overlaps": True,
        "overlap_bp": overlap_bp,
        "overlap_interval": [overlap_start, overlap_end],
        "overlap_frac_of_a": round(overlap_bp / a_len, 6),
        "overlap_frac_of_b": round(overlap_bp / b_len, 6),
    }

# =============================================================================
# BINNING FOR FAST PREFILTER
# =============================================================================

def add_bins(
    df: pd.DataFrame,
    start_col: str,
    end_col: str,
    bin_size: int,
) -> pd.DataFrame:
    """
    Add coarse binning columns for fast prefilter joins.
    
    Adds columns:
        - bin_start: start position bin index
        - bin_end: end position bin index
        - bins: list of all bins the interval spans
    """
    out = df.copy()
    out["bin_start"] = (out[start_col] // bin_size).astype(np.int64)
    out["bin_end"] = (out[end_col] // bin_size).astype(np.int64)
    out["bins"] = [
        list(range(s, e + 1))
        for s, e in zip(out["bin_start"], out["bin_end"])
    ]
    return out


# =============================================================================
# TIER/CATEGORY HELPERS
# =============================================================================

def make_cut_edges_for_labels(
    base_edges: List[int],
    include_zero_bin: bool = True,
) -> List[float]:
    """
    Make edges so that len(labels) == len(edges) - 1.
    
    With base_edges=[0, 100k, 250k, 500k, 1M] and 4 labels,
    we need edges=[-0.1, 100k, 250k, 500k, 1M].
    """
    if include_zero_bin:
        return [-0.1] + base_edges[1:]
    return base_edges


def assign_distance_tier(
    distances: pd.Series,
    tier_edges: List[int],
    tier_labels: List[str],
) -> pd.Categorical:
    """
    Assign distance tier labels to a series of distances.
    """
    bins = make_cut_edges_for_labels(tier_edges, include_zero_bin=True)
    return pd.cut(
        distances,
        bins=bins,
        labels=tier_labels,
        include_lowest=True,
        right=True,
        ordered=True,
    )


# =============================================================================
# PARSING HELPERS
# =============================================================================

def safe_parse_list(x) -> List:
    """
    Safely parse a value that might be a list, string representation of list, or empty.
    """
    if pd.isna(x) or x in ["", "[]"]:
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        try:
            return ast.literal_eval(x)
        except (ValueError, SyntaxError):
            # Fall back: split comma-separated strings
            return [g.strip() for g in x.split(",") if g.strip()]
    return []


def extract_primary_class(s: str) -> str:
    """
    Extract primary cCRE class from comma-separated classification string.
    
    Example:
        'pELS,CTCF-bound' -> 'pELS'
        'CTCF-only' -> 'CTCF-only'
    """
    if not isinstance(s, str):
        return ""
    parts = [p.strip() for p in s.split(",") if p.strip()]
    return parts[0] if parts else ""


def parse_hichip_anchor2(value: str, source_path: str = "") -> Tuple[str, int, int, float]:
    """
    Parse LoopCatalog anchor2 format: chr2:start2-end2,score
    
    Returns:
        Tuple of (chr2, start2, end2, score)
    """
    m = re.match(r"(chr[^:]+):(\d+)-(\d+),(.+)", str(value))
    if not m:
        raise ValueError(f"Unparseable anchor2 format in {source_path}: {value}")
    
    c2, s2, e2, v = m.groups()
    v = v.strip()
    score = float("inf") if v.lower() == "inf" else float(v)
    return c2, int(s2), int(e2), score


# =============================================================================
# TYPE CONVERSIONS
# =============================================================================

def as_categorical(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """Convert specified columns to categorical dtype."""
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = out[c].astype("category")
    return out


def safe_float(val, default: Optional[float] = None) -> Optional[float]:
    """Safely convert to float, returning default if NaN or conversion fails."""
    if pd.isna(val):
        return default
    try:
        result = float(val)
        return result if np.isfinite(result) else default
    except (ValueError, TypeError):
        return default


def safe_int(val, default: Optional[int] = None) -> Optional[int]:
    """Safely convert to int, returning default if NaN or conversion fails."""
    if pd.isna(val):
        return default
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


# =============================================================================
# DATAFRAME HELPERS
# =============================================================================

def ensure_columns_exist(
    df: pd.DataFrame,
    required_cols: List[str],
    source_name: str = "DataFrame",
) -> None:
    """Raise ValueError if any required columns are missing."""
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"{source_name} missing required columns: {missing}")


def fill_numeric_columns(df: pd.DataFrame, cols: List[str], fill_value: int = 0) -> pd.DataFrame:
    """Fill NaN in numeric columns and convert to int."""
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(fill_value).astype(int)
    return df


def fill_string_columns(df: pd.DataFrame, cols: List[str], fill_value: str = "") -> pd.DataFrame:
    """Fill NaN in string columns."""
    df = df.copy()
    for c in cols:
        if c in df.columns:
            df[c] = df[c].astype(object).fillna(fill_value)
    return df


# =============================================================================
# PROMOTER COORDINATE HELPERS
# =============================================================================

def compute_tss(
    start: pd.Series,
    end: pd.Series,
    strand: pd.Series,
) -> pd.Series:
    """Compute TSS position based on strand."""
    return np.where(strand.astype(str) == "+", start, end)


def compute_promoter_coords(
    start: pd.Series,
    end: pd.Series,
    strand: pd.Series,
    upstream_bp: int = 2000,
    downstream_bp: int = 500,
) -> Tuple[pd.Series, pd.Series]:
    """
    Compute promoter start and end coordinates.
    
    Returns:
        Tuple of (prom_start, prom_end) Series
    """
    is_plus = strand.astype(str) == "+"
    
    prom_start = np.where(
        is_plus,
        (start - upstream_bp).clip(lower=0),
        end - downstream_bp,
    )
    
    prom_end = np.where(
        is_plus,
        start + downstream_bp,
        end + upstream_bp,
    )
    
    return pd.Series(prom_start, index=start.index), pd.Series(prom_end, index=end.index)


# =============================================================================
# STRENGTH CLASSIFICATION
# =============================================================================

def classify_strength(is_strong: bool, is_weak: bool) -> str:
    """Classify link strength based on boolean flags."""
    if is_strong:
        return "strong"
    return "weak" if is_weak else "none"


def classify_strength_vectorized(is_strong: pd.Series, is_weak: pd.Series) -> pd.Series:
    """Vectorized strength classification."""
    return pd.Series([
        classify_strength(s, w) for s, w in zip(is_strong, is_weak)
    ], index=is_strong.index)
