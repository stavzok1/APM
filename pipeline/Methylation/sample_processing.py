"""
Per-sample methylation data processing.

Handles loading, validation, and enrichment of sample-level beta values.

Main functions:
    load_sample_beta: Load sample beta values from file
    validate_sample_beta: Validate beta value DataFrame
    enrich_sample_with_annotations: Add reference annotations to sample data
    compute_m_values: Convert beta to M-values
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, List, Tuple, Union

import pandas as pd
import numpy as np

from .meth_schemas import (
    METHYLATION_COLUMNS as MC,
    validate_sample_beta_df,
)


# =============================================================================
# SAMPLE LOADING
# =============================================================================

def load_sample_beta(
    path: Path,
    probe_col: str = "probeID",
    beta_col: str = "beta",
    sep: str = "\t",
    **kwargs,
) -> pd.DataFrame:
    """
    Load sample beta values from file.
    
    Expected format: two-column file with probe ID and beta value.
    
    Args:
        path: Path to sample beta file
        probe_col: Column name or index for probe ID (can be int for positional)
        beta_col: Column name or index for beta value
        sep: Column separator
        **kwargs: Additional arguments passed to pd.read_csv
    
    Returns:
        DataFrame with probeID and beta columns
    """
    path = Path(path)
    
    # Load file
    df = pd.read_csv(path, sep=sep, header=None, **kwargs)
    
    # Handle column specification
    # if isinstance(probe_col, int):
    #     probe_col = df.columns[probe_col]
    # if isinstance(beta_col, int):
    #     beta_col = df.columns[beta_col]
    
    probe_col = df.columns[0]
    beta_col = df.columns[1]
    
    # Rename to standard names
    if probe_col != "probeID":
        df = df.rename(columns={probe_col: "probeID"})
    if beta_col != "beta":
        df = df.rename(columns={beta_col: "beta"})
    
    # Ensure we have required columns
    if "probeID" not in df.columns or "beta" not in df.columns:
        raise ValueError(
            f"Could not find required columns. "
            f"Available: {df.columns.tolist()}"
        )
    
    # Convert beta to numeric
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    
    return df[["probeID", "beta"]].copy()


def load_sample_beta_matrix(
    path: Path,
    probe_col: str = "probeID",
    sample_cols: Optional[List[str]] = None,
    sep: str = "\t",
) -> pd.DataFrame:
    """
    Load sample beta matrix (probes x samples).
    
    For cases where multiple samples are in one file.
    
    Args:
        path: Path to matrix file
        probe_col: Column name for probe IDs
        sample_cols: List of sample column names (None = all non-probe columns)
        sep: Column separator
    
    Returns:
        DataFrame with probeID as index and samples as columns
    """
    df = pd.read_csv(path, sep=sep)
    
    if probe_col not in df.columns:
        raise ValueError(f"Probe column '{probe_col}' not found")
    
    df = df.set_index(probe_col)
    
    if sample_cols is not None:
        df = df[sample_cols]
    
    # Convert all to numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    
    return df


# =============================================================================
# VALIDATION
# =============================================================================

def validate_sample_beta(
    df: pd.DataFrame,
    strict: bool = True,
) -> Tuple[bool, List[str]]:
    """
    Validate sample beta DataFrame.
    
    Checks:
        - Required columns present
        - Beta values in [0, 1] range
        - No duplicate probe IDs
        - Reasonable number of valid values
    
    Args:
        df: Sample beta DataFrame
        strict: If True, treat warnings as errors
    
    Returns:
        Tuple of (is_valid, list of issues)
    """
    issues = validate_sample_beta_df(df)
    
    # Additional checks
    if "probeID" in df.columns:
        n_dup = df["probeID"].duplicated().sum()
        if n_dup > 0:
            issues.append(f"Duplicate probe IDs: {n_dup}")
    
    if "beta" in df.columns:
        n_valid = df["beta"].notna().sum()
        n_total = len(df)
        pct_valid = 100 * n_valid / n_total if n_total > 0 else 0
        
        if pct_valid < 50:
            msg = f"Low valid beta percentage: {pct_valid:.1f}%"
            if strict:
                issues.append(msg)
            else:
                print(f"  Warning: {msg}")
    
    is_valid = len(issues) == 0
    return is_valid, issues


# =============================================================================
# M-VALUE CONVERSION
# =============================================================================

def compute_m_values(
    beta: Union[pd.Series, np.ndarray],
    offset: float = 0.001,
) -> Union[pd.Series, np.ndarray]:
    """
    Convert beta values to M-values.
    
    M = log2(beta / (1 - beta))
    
    M-values are more appropriate for statistical analysis as they
    are approximately normally distributed.
    
    Args:
        beta: Beta values (0-1 range)
        offset: Small value to avoid log(0) and log(inf)
    
    Returns:
        M-values
    """
    # Clip beta to avoid log issues
    beta_clipped = np.clip(beta, offset, 1 - offset)
    
    m_values = np.log2(beta_clipped / (1 - beta_clipped))
    
    if isinstance(beta, pd.Series):
        return pd.Series(m_values, index=beta.index)
    return m_values


def compute_beta_from_m(
    m_values: Union[pd.Series, np.ndarray],
) -> Union[pd.Series, np.ndarray]:
    """
    Convert M-values back to beta values.
    
    beta = 2^M / (2^M + 1)
    """
    two_pow_m = np.power(2, m_values)
    beta = two_pow_m / (two_pow_m + 1)
    
    if isinstance(m_values, pd.Series):
        return pd.Series(beta, index=m_values.index)
    return beta


# =============================================================================
# SAMPLE ENRICHMENT
# =============================================================================

def enrich_sample_with_annotations(
    sample_df: pd.DataFrame,
    probe_reference: pd.DataFrame,
    columns_to_add: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Add reference annotations to sample beta data.
    
    Args:
        sample_df: Sample DataFrame with probeID and beta
        probe_reference: Annotated probe reference
        columns_to_add: Specific columns to add (None = key columns)
    
    Returns:
        Enriched sample DataFrame
    """
    if columns_to_add is None:
        columns_to_add = [
            "chrom", "start", "end",
            "in_CGI", "CGI_context",
            "in_promoter", "promoter_genes",
            "in_gene_body", "gene_body_genes",
            "overlapping_ccres", "ccre_types",
            "in_lncrna_promoter", "lncrna_promoter_genes",
        ]
    
    # Filter to columns that exist in reference
    cols_available = [c for c in columns_to_add if c in probe_reference.columns]
    
    # Merge
    enriched = sample_df.merge(
        probe_reference[["probeID"] + cols_available],
        on="probeID",
        how="left",
    )
    
    return enriched


# =============================================================================
# SAMPLE PROCESSING PIPELINE
# =============================================================================

def process_sample_methylation(
    sample_path: Path,
    probe_reference: pd.DataFrame,
    sample_id: Optional[str] = None,
    compute_m: bool = True,
    enrich: bool = True,
    validate: bool = True,
) -> pd.DataFrame:
    """
    Complete processing pipeline for a single sample.
    
    Steps:
        1. Load beta values
        2. Validate data
        3. Compute M-values (optional)
        4. Enrich with reference annotations (optional)
    
    Args:
        sample_path: Path to sample beta file
        probe_reference: Annotated probe reference
        sample_id: Sample identifier (for logging)
        compute_m: Whether to compute M-values
        enrich: Whether to add reference annotations
        validate: Whether to validate data
    
    Returns:
        Processed sample DataFrame
    """
    sample_id = sample_id or Path(sample_path).stem
    print(f"Processing sample: {sample_id}")
    
    # Load
    sample_df = load_sample_beta(sample_path)
    print(f"  Loaded {len(sample_df)} probes")
    
    # Validate
    if validate:
        is_valid, issues = validate_sample_beta(sample_df, strict=False)
        if issues:
            for issue in issues:
                print(f"  Warning: {issue}")
    
    # Compute M-values
    if compute_m:
        sample_df["m_value"] = compute_m_values(sample_df["beta"])
    
    # Enrich with annotations
    if enrich:
        sample_df = enrich_sample_with_annotations(sample_df, probe_reference)
        print(f"  Enriched with annotations")
    
    # Add sample ID
    sample_df["sample_id"] = sample_id
    
    # Summary stats
    n_valid = sample_df["beta"].notna().sum()
    mean_beta = sample_df["beta"].mean()
    print(f"  Valid probes: {n_valid}")
    print(f"  Mean beta: {mean_beta:.3f}")
    
    return sample_df


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def process_sample_batch(
    sample_paths: List[Path],
    probe_reference: pd.DataFrame,
    sample_ids: Optional[List[str]] = None,
    output_dir: Optional[Path] = None,
    compute_m: bool = True,
    enrich: bool = True,
) -> List[pd.DataFrame]:
    """
    Process multiple samples.
    
    Args:
        sample_paths: List of paths to sample files
        probe_reference: Annotated probe reference
        sample_ids: Optional list of sample IDs
        output_dir: Optional directory to save processed files
        compute_m: Whether to compute M-values
        enrich: Whether to add reference annotations
    
    Returns:
        List of processed sample DataFrames
    """
    if sample_ids is None:
        sample_ids = [Path(p).stem for p in sample_paths]
    
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    
    for path, sid in zip(sample_paths, sample_ids):
        try:
            sample_df = process_sample_methylation(
                path,
                probe_reference,
                sample_id=sid,
                compute_m=compute_m,
                enrich=enrich,
            )
            
            if output_dir:
                out_path = output_dir / f"{sid}_probes.parquet"
                sample_df.to_parquet(out_path, index=False)
                print(f"  Saved to {out_path}")
            
            results.append(sample_df)
            
        except Exception as e:
            print(f"  Error processing {sid}: {e}")
            continue
    
    print(f"\nProcessed {len(results)}/{len(sample_paths)} samples successfully")
    
    return results


# =============================================================================
# QUALITY CONTROL
# =============================================================================

def compute_sample_qc_metrics(
    sample_df: pd.DataFrame,
) -> dict:
    """
    Compute QC metrics for a sample.
    
    Returns dict with metrics like:
        - n_probes: total probes
        - n_valid: probes with valid beta
        - pct_valid: percentage valid
        - mean_beta: mean beta value
        - median_beta: median beta value
        - std_beta: standard deviation
        - pct_hypermeth: fraction beta > 0.7
        - pct_hypometh: fraction beta < 0.3
        - pct_cgi_probes: fraction of probes in CGI
    """
    beta = sample_df["beta"].dropna()
    print("******entered compute_sample_qc_metrics******")
   
    
    metrics = {
        "n_probes": len(sample_df),
        "n_valid": len(beta),
        "pct_valid": 100 * len(beta) / len(sample_df) if len(sample_df) > 0 else 0,
        "mean_beta": float(beta.mean()) if len(beta) > 0 else None,
        "median_beta": float(beta.median()) if len(beta) > 0 else None,
        "std_beta": float(beta.std()) if len(beta) > 0 else None,
        "pct_hypermeth": 100 * (beta > 0.7).mean() if len(beta) > 0 else None,
        "pct_hypometh": 100 * (beta < 0.3).mean() if len(beta) > 0 else None,
    }
    

    # CGI metrics if available
    if "in_CGI" in sample_df.columns:
        cgi_mask = sample_df["in_CGI"].fillna(False).astype(bool)
        cgi_beta = sample_df.loc[cgi_mask, "beta"].dropna()
        
        metrics["n_cgi_probes"] = int(cgi_mask.sum())
        metrics["pct_cgi_probes"] = 100 * cgi_mask.mean()
        metrics["mean_cgi_beta"] = float(cgi_beta.mean()) if len(cgi_beta) > 0 else None
        print("******exited compute_sample_qc_metrics******")

    
    return metrics


def filter_low_quality_probes(
    sample_df: pd.DataFrame,
    detection_pval_col: Optional[str] = None,
    detection_pval_threshold: float = 0.01,
    min_beta: float = 0.0,
    max_beta: float = 1.0,
) -> pd.DataFrame:
    """
    Filter probes based on quality criteria.
    
    Args:
        sample_df: Sample DataFrame
        detection_pval_col: Column with detection p-values (optional)
        detection_pval_threshold: Max detection p-value to keep
        min_beta: Minimum valid beta
        max_beta: Maximum valid beta
    
    Returns:
        Filtered DataFrame
    """
    mask = pd.Series(True, index=sample_df.index)
    
    # Beta range filter
    if "beta" in sample_df.columns:
        beta = sample_df["beta"]
        mask &= (beta >= min_beta) & (beta <= max_beta) | beta.isna()
    
    # Detection p-value filter
    if detection_pval_col and detection_pval_col in sample_df.columns:
        pval = sample_df[detection_pval_col]
        mask &= (pval <= detection_pval_threshold) | pval.isna()
    
    n_filtered = (~mask).sum()
    if n_filtered > 0:
        print(f"  Filtered {n_filtered} low-quality probes")
    
    return sample_df[mask].copy()
