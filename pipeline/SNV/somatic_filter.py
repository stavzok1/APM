"""
High-confidence somatic variant filtering.

Provides configurable filtering for Mutect2 SNV/indel calls based on:
- VAF thresholds (tumor/normal)
- Read depth requirements
- TLOD scores
- Population allele frequency (POPAF)
- FILTER status
"""

from typing import Optional

import numpy as np
import pandas as pd


def high_conf_somatic_mask(
    df: pd.DataFrame,
    min_tumor_vaf: Optional[float] = None,
    max_normal_vaf: Optional[float] = None,
    min_tlod: Optional[float] = None,
    min_popaf: Optional[float] = None,
    min_tumor_dp: Optional[int] = None,
    min_normal_dp: Optional[int] = None,
    use_popaf: Optional[bool] = None,
    require_tlod: Optional[bool] = None,
    require_pass: Optional[bool] = None,
    tumor_dp_col: str = "TUMOR_DP",
    normal_dp_col: str = "NORMAL_DP",
) -> pd.Series:
    """
    Generate boolean mask for high-confidence somatic variants.
    
    Args:
        df: DataFrame with SNV data (from load_mutect_snv_vcf)
        min_tumor_vaf: Minimum tumor VAF (default from config: 0.05 = 5%)
        max_normal_vaf: Maximum normal VAF (default from config: 0.02 = 2%)
        min_tlod: Minimum tumor log-odds score (default from config: 6.0)
        min_popaf: Minimum -log10(population AF) (default from config: 3.0 → AF ≤ 1e-3)
        min_tumor_dp: Minimum tumor read depth (default from config: 20)
        min_normal_dp: Minimum normal read depth (default from config: 10)
        use_popaf: Whether to filter on population AF (default from config: True)
        require_tlod: If True, TLOD must be present and pass threshold.
                      If False, missing TLOD is allowed (default from config: False)
        require_pass: Require FILTER == "PASS" (default from config: True)
        tumor_dp_col: Column name for tumor depth
        normal_dp_col: Column name for normal depth
    
    Returns:
        Boolean Series indicating high-confidence somatic variants
    
    Notes:
        - POPAF from Mutect2 is stored as -log10(population AF)
        - Higher POPAF = rarer variant = more likely somatic
        - Missing values in optional fields (TLOD, POPAF, depth) are allowed
          unless require_* flags are set
        - All None parameters pull defaults from pipeline.config.THRESHOLDS
    """
    # Import config and apply defaults for any None values
    try:
        from ..config import THRESHOLDS
        if min_tumor_vaf is None:
            min_tumor_vaf = THRESHOLDS.snv_min_tumor_vaf
        if max_normal_vaf is None:
            max_normal_vaf = THRESHOLDS.snv_max_normal_vaf
        if min_tlod is None:
            min_tlod = THRESHOLDS.snv_min_tlod
        if min_popaf is None:
            min_popaf = THRESHOLDS.snv_min_popaf
        if min_tumor_dp is None:
            min_tumor_dp = THRESHOLDS.snv_min_tumor_dp
        if min_normal_dp is None:
            min_normal_dp = THRESHOLDS.snv_min_normal_dp
        if use_popaf is None:
            use_popaf = THRESHOLDS.snv_use_popaf
        if require_tlod is None:
            require_tlod = THRESHOLDS.snv_require_tlod
        if require_pass is None:
            require_pass = THRESHOLDS.snv_require_pass
    except ImportError:
        # Fallback defaults if config not available (standalone usage)
        if min_tumor_vaf is None:
            min_tumor_vaf = 0.05
        if max_normal_vaf is None:
            max_normal_vaf = 0.02
        if min_tlod is None:
            min_tlod = 6.0
        if min_popaf is None:
            min_popaf = 3.0
        if min_tumor_dp is None:
            min_tumor_dp = 20
        if min_normal_dp is None:
            min_normal_dp = 10
        if use_popaf is None:
            use_popaf = True
        if require_tlod is None:
            require_tlod = False
        if require_pass is None:
            require_pass = True
    
    # FILTER status
    df_filter = df.get("filter", pd.Series("", index=df.index))
    if require_pass:
        pass_mask = df_filter == "PASS"
    else:
        pass_mask = pd.Series(True, index=df.index)
    
    # VAF thresholds
    tumor_vaf = df.get("tumor_vaf", pd.Series(np.nan, index=df.index))
    tumor_vaf_mask = tumor_vaf >= min_tumor_vaf
    
    normal_vaf = df.get("normal_vaf", pd.Series(np.nan, index=df.index))
    # Allow NaN normal_vaf (tumor-only samples)
    normal_vaf_mask = (normal_vaf <= max_normal_vaf) | normal_vaf.isna()
    
    # Read depth thresholds
    tumor_dp = pd.to_numeric(df.get(tumor_dp_col, np.nan), errors="coerce")
    tumor_dp_mask = (tumor_dp >= min_tumor_dp) | tumor_dp.isna()
    
    normal_dp = pd.to_numeric(df.get(normal_dp_col, np.nan), errors="coerce")
    normal_dp_mask = (normal_dp >= min_normal_dp) | normal_dp.isna()
    
    # TLOD threshold
    tlod = pd.to_numeric(df.get("TLOD", np.nan), errors="coerce")
    if require_tlod:
        tlod_mask = tlod >= min_tlod
    else:
        tlod_mask = (tlod >= min_tlod) | tlod.isna()
    
    # Population AF threshold
    if use_popaf:
        popaf = pd.to_numeric(df.get("POPAF", np.nan), errors="coerce")
        popaf_mask = (popaf >= min_popaf) | popaf.isna()
    else:
        popaf_mask = pd.Series(True, index=df.index)
    
    # Combine all masks
    mask = (
        pass_mask
        & tumor_vaf_mask
        & normal_vaf_mask
        & tumor_dp_mask
        & normal_dp_mask
        & tlod_mask
        & popaf_mask
    )
    
    return mask


def apply_somatic_filter(
    df: pd.DataFrame,
    **filter_kwargs,
) -> pd.DataFrame:
    """
    Apply high-confidence somatic filter and return filtered DataFrame.
    
    Args:
        df: DataFrame with SNV data
        **filter_kwargs: Arguments passed to high_conf_somatic_mask
    
    Returns:
        Filtered DataFrame containing only high-confidence somatic variants
    """
    mask = high_conf_somatic_mask(df, **filter_kwargs)
    return df[mask].copy()


def summarize_filtering(
    df_before: pd.DataFrame,
    df_after: pd.DataFrame,
) -> dict:
    """
    Summarize the effect of somatic filtering.
    
    Returns:
        Dictionary with filtering statistics
    """
    n_before = len(df_before)
    n_after = len(df_after)
    
    summary = {
        "n_variants_before": n_before,
        "n_variants_after": n_after,
        "n_filtered": n_before - n_after,
        "retention_rate": n_after / n_before if n_before > 0 else 0.0,
    }
    
    # Breakdown by filter reason if possible
    if "filter" in df_before.columns:
        filter_counts = df_before["filter"].value_counts().to_dict()
        summary["filter_breakdown"] = filter_counts
    
    return summary