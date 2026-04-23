"""
SV filtering functions for strict and lenient variant sets.

Provides quality-based filtering for somatic SV calls from Manta.
"""

from typing import List, Optional

import pandas as pd


# =============================================================================
# FILTERING FUNCTIONS
# =============================================================================

def get_strict_sv_set(
    df: pd.DataFrame,
    allowed_types: Optional[List[str]] = None,
    min_tumor_sr_alt: int = 2,
    min_tumor_alt: int = 8,
    max_normal_alt: int = 1,
    min_somatic_score: int = 25,
) -> pd.DataFrame:
    """
    Apply strict filtering to SV DataFrame.
    
    Strict criteria:
    - PASS filter
    - Allowed SV types (DEL, DUP, INV, INS, BND)
    - Minimum split-read support in tumor
    - Minimum total alt support in tumor
    - Maximum alt support in normal
    - Minimum somatic score
    
    Args:
        df: SV DataFrame from load_manta_sv_vcf
        allowed_types: List of allowed SVTYPE values
        min_tumor_sr_alt: Minimum tumor split-read alt count
        min_tumor_alt: Minimum tumor total alt count
        max_normal_alt: Maximum normal alt count
        min_somatic_score: Minimum SOMATICSCORE
    
    Returns:
        Filtered DataFrame
    """
    if allowed_types is None:
        allowed_types = ["DEL", "DUP", "INV", "INS", "BND"]

    has_svlen = "SVLEN" in df.columns
    has_score = "SOMATICSCORE" in df.columns

    q = (
        (df["filter"] == "PASS") &
        (df["SVTYPE"].isin(allowed_types)) &
        (df["tumor_sr_alt"] >= min_tumor_sr_alt) &
        (df["tumor_alt"] >= min_tumor_alt) &
        (df["normal_alt"] <= max_normal_alt)
    )

    if has_score:
        q &= df["SOMATICSCORE"].fillna(0) >= min_somatic_score

    return df[q].copy()


def get_lenient_sv_set(
    df: pd.DataFrame,
    allowed_types: Optional[List[str]] = None,
    min_tumor_sr_alt: int = 1,
    min_tumor_pr_alt: int = 5,
    max_normal_alt: int = 1,
    min_somatic_score: int = 15,
) -> pd.DataFrame:
    """
    Apply lenient filtering to SV DataFrame.
    
    Lenient criteria:
    - PASS filter
    - Allowed SV types
    - Either split-read OR paired-read support
    - Maximum alt support in normal
    - Lower somatic score threshold
    
    Args:
        df: SV DataFrame from load_manta_sv_vcf
        allowed_types: List of allowed SVTYPE values
        min_tumor_sr_alt: Minimum tumor split-read alt count
        min_tumor_pr_alt: Minimum tumor paired-read alt count (alternative)
        max_normal_alt: Maximum normal alt count
        min_somatic_score: Minimum SOMATICSCORE
    
    Returns:
        Filtered DataFrame
    """
    if allowed_types is None:
        allowed_types = ["DEL", "DUP", "INV", "INS", "BND"]

    has_svlen = "SVLEN" in df.columns
    has_score = "SOMATICSCORE" in df.columns

    q = (
        (df["filter"] == "PASS") &
        (df["SVTYPE"].isin(allowed_types)) &
        ((df["tumor_sr_alt"] >= min_tumor_sr_alt) | (df["tumor_pr_alt"] >= min_tumor_pr_alt)) &
        (df["normal_alt"] <= max_normal_alt)
    )

    if has_score:
        q &= df["SOMATICSCORE"].fillna(0) >= min_somatic_score

    return df[q].copy()


def filter_by_svlen(
    df: pd.DataFrame,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
) -> pd.DataFrame:
    """
    Filter SVs by length.
    
    Args:
        df: SV DataFrame
        min_len: Minimum absolute SVLEN
        max_len: Maximum absolute SVLEN
    
    Returns:
        Filtered DataFrame
    """
    if "SVLEN" not in df.columns:
        return df.copy()
    
    mask = pd.Series([True] * len(df), index=df.index)
    abs_len = df["SVLEN"].abs()
    
    if min_len is not None:
        mask &= abs_len >= min_len
    if max_len is not None:
        mask &= abs_len <= max_len
    
    return df[mask].copy()


def filter_by_region(
    df: pd.DataFrame,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> pd.DataFrame:
    """
    Filter SVs by genomic region.
    
    Args:
        df: SV DataFrame
        chrom: Chromosome to filter to
        start: Start position (inclusive)
        end: End position (inclusive)
    
    Returns:
        Filtered DataFrame
    """
    mask = pd.Series([True] * len(df), index=df.index)
    
    if chrom is not None:
        mask &= df["chrom"] == chrom
    if start is not None:
        mask &= df["pos"] >= start
    if end is not None:
        mask &= df["pos"] <= end
    
    return df[mask].copy()


def get_sv_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get summary statistics for SV DataFrame.
    
    Returns DataFrame with counts per SVTYPE.
    """
    if df.empty:
        return pd.DataFrame(columns=["SVTYPE", "count"])
    
    summary = df.groupby("SVTYPE").size().reset_index(name="count")
    summary = summary.sort_values("count", ascending=False)
    
    return summary
