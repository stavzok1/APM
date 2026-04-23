"""
Methylation aggregation to gene/lncRNA/cCRE level.

Aggregates probe-level beta values to feature-level summaries with
comprehensive metrics for promoter, gene body, and regulatory element
methylation.

Main functions:
    aggregate_to_genes: Gene-level promoter and body methylation
    aggregate_to_lncrnas: lncRNA-level methylation
    aggregate_to_ccres: cCRE-level methylation
    compute_promoter_methylation_score: Composite methylation score
"""

from __future__ import annotations

from typing import List, Dict, Any, Optional, Callable

import pandas as pd
import numpy as np

from .meth_schemas import (
    empty_gene_methylation_entry,
    empty_lncrna_methylation_entry,
    empty_ccre_methylation_entry,
    CGI_CONTEXTS,
)

from ..config import THRESHOLDS, PATHS, OUTPUT_SUBDIRS



# =============================================================================
# THRESHOLDS
# =============================================================================

HYPERMETH_THRESHOLD = THRESHOLDS.meth_hypermeth_threshold
HYPOMETH_THRESHOLD = THRESHOLDS.meth_hypometh_threshold
INTERMEDIATE_LOW = THRESHOLDS.meth_intermediate_low
INTERMEDIATE_HIGH = THRESHOLDS.meth_intermediate_high


# =============================================================================
# GENE-LEVEL AGGREGATION
# =============================================================================

def aggregate_to_genes(
    sample_probes: pd.DataFrame,
    probe_reference: pd.DataFrame,
    gene_panel: List[str],
    include_gene_body: bool = True,
    hypermeth_threshold: float = HYPERMETH_THRESHOLD,
    hypometh_threshold: float = HYPOMETH_THRESHOLD,
) -> pd.DataFrame:
    """
    Aggregate probe beta values to gene-level metrics.
    
    Computes metrics for:
        - Promoter region (defined by promoter_genes annotation)
        - Gene body (if include_gene_body=True)
        - CGI-specific promoter methylation
        - Shore-specific promoter methylation
    
    Args:
        sample_probes: DataFrame with probeID and beta columns
        probe_reference: Annotated probe reference with promoter_genes column
        gene_panel: List of genes to aggregate
        include_gene_body: Whether to include gene body metrics
        hypermeth_threshold: Beta threshold for hypermethylation
        hypometh_threshold: Beta threshold for hypomethylation
    
    Returns:
        DataFrame with one row per gene and aggregation columns
    """
    # Merge beta values with reference
    merged = sample_probes[["probeID", "beta"]].merge(
        probe_reference[[
            "probeID", "promoter_genes", "gene_body_genes",
            "in_CGI", "CGI_context"
        ]],
        on="probeID",
        how="inner",
    )
    
    # Drop probes with no beta value
    merged = merged.dropna(subset=["beta"])
    
    if merged.empty:
        return _empty_gene_agg_df(gene_panel)
    
    # === PROMOTER AGGREGATION ===
    # Explode promoter_genes to get one row per (probe, gene)
    prom_df = merged[merged["promoter_genes"].apply(len) > 0].copy()
    prom_df = prom_df.explode("promoter_genes")
    prom_df = prom_df[prom_df["promoter_genes"].isin(gene_panel)]
    prom_df = prom_df.rename(columns={"promoter_genes": "gene_name"})
    
    promoter_agg = _aggregate_region(
        prom_df,
        group_col="gene_name",
        prefix="promoter",
        hypermeth_threshold=hypermeth_threshold,
        hypometh_threshold=hypometh_threshold,
        include_cgi_context=True,
    )
    
    # === GENE BODY AGGREGATION ===
    if include_gene_body:
        body_df = merged[merged["gene_body_genes"].apply(len) > 0].copy()
        body_df = body_df.explode("gene_body_genes")
        body_df = body_df[body_df["gene_body_genes"].isin(gene_panel)]
        body_df = body_df.rename(columns={"gene_body_genes": "gene_name"})
        
        body_agg = _aggregate_region(
            body_df,
            group_col="gene_name",
            prefix="gene_body",
            hypermeth_threshold=hypermeth_threshold,
            hypometh_threshold=hypometh_threshold,
            include_cgi_context=False,
        )
    else:
        body_agg = pd.DataFrame({"gene_name": []})
    
    # === COMBINE ===
    # Start with full gene panel
    result = pd.DataFrame({"gene_name": gene_panel})
    
    # Merge promoter metrics
    result = result.merge(promoter_agg, on="gene_name", how="left")
    
    # Merge body metrics
    if include_gene_body and not body_agg.empty:
        result = result.merge(body_agg, on="gene_name", how="left")
    
    # Fill NaN counts with 0
    count_cols = [c for c in result.columns if c.endswith("_n_probes")]
    for col in count_cols:
        result[col] = result[col].fillna(0).astype(int)
    
    return result


def _aggregate_region(
    df: pd.DataFrame,
    group_col: str,
    prefix: str,
    hypermeth_threshold: float,
    hypometh_threshold: float,
    include_cgi_context: bool = False,
) -> pd.DataFrame:
    """
    Aggregate methylation for a genomic region.
    
    Generic helper used by gene, lncRNA, and cCRE aggregation.
    """
    if df.empty:
        return pd.DataFrame({group_col: []})
    
    # Basic aggregation
    agg = df.groupby(group_col).agg(
        beta_mean=("beta", "mean"),
        beta_median=("beta", "median"),
        beta_std=("beta", "std"),
        beta_min=("beta", "min"),
        beta_max=("beta", "max"),
        n_probes=("beta", "count"),
    ).reset_index()
    
    # Rename with prefix
    agg = agg.rename(columns={
        "beta_mean": f"{prefix}_beta_mean",
        "beta_median": f"{prefix}_beta_median",
        "beta_std": f"{prefix}_beta_std",
        "beta_min": f"{prefix}_beta_min",
        "beta_max": f"{prefix}_beta_max",
        "n_probes": f"{prefix}_n_probes",
    })
    
    # Compute range
    agg[f"{prefix}_beta_range"] = agg[f"{prefix}_beta_max"] - agg[f"{prefix}_beta_min"]
    
    # Compute hyper/hypo fractions
    hyper_counts = df[df["beta"] > hypermeth_threshold].groupby(group_col).size()
    hypo_counts = df[df["beta"] < hypometh_threshold].groupby(group_col).size()
    
    agg[f"{prefix}_frac_hypermeth"] = (
        agg[group_col].map(hyper_counts).fillna(0) / agg[f"{prefix}_n_probes"]
    )
    agg[f"{prefix}_frac_hypometh"] = (
        agg[group_col].map(hypo_counts).fillna(0) / agg[f"{prefix}_n_probes"]
    )
    
    # CGI-specific metrics
    if include_cgi_context and "in_CGI" in df.columns:
        # CGI probes
        cgi_df = df[df["in_CGI"] == True]
        if not cgi_df.empty:
            cgi_means = cgi_df.groupby(group_col)["beta"].mean()
            cgi_counts = cgi_df.groupby(group_col).size()
            agg[f"{prefix}_CGI_beta_mean"] = agg[group_col].map(cgi_means)
            agg[f"{prefix}_CGI_n_probes"] = agg[group_col].map(cgi_counts).fillna(0).astype(int)
        else:
            agg[f"{prefix}_CGI_beta_mean"] = np.nan
            agg[f"{prefix}_CGI_n_probes"] = 0
        
        # Shore probes
        if "CGI_context" in df.columns:
            shore_df = df[df["CGI_context"].isin(CGI_CONTEXTS.SHORES)]
            if not shore_df.empty:
                shore_means = shore_df.groupby(group_col)["beta"].mean()
                agg[f"{prefix}_shore_beta_mean"] = agg[group_col].map(shore_means)
            else:
                agg[f"{prefix}_shore_beta_mean"] = np.nan
        else:
            agg[f"{prefix}_shore_beta_mean"] = np.nan

    return agg


def _empty_gene_agg_df(gene_panel: List[str]) -> pd.DataFrame:
    """Return empty gene aggregation DataFrame with correct schema."""
    return pd.DataFrame({
        "gene_name": gene_panel,
        "promoter_beta_mean": np.nan,
        "promoter_beta_median": np.nan,
        "promoter_beta_std": np.nan,
        "promoter_beta_min": np.nan,
        "promoter_beta_max": np.nan,
        "promoter_beta_range": np.nan,
        "promoter_n_probes": 0,
        "promoter_frac_hypermeth": np.nan,
        "promoter_frac_hypometh": np.nan,
        "promoter_CGI_beta_mean": np.nan,
        "promoter_shore_beta_mean": np.nan,
        "gene_body_beta_mean": np.nan,
        "gene_body_n_probes": 0,
    })


# =============================================================================
# LNCRNA-LEVEL AGGREGATION
# =============================================================================

def _lncrna_agg_output_columns() -> List[str]:
    """Column names produced by ``_aggregate_region`` for lncRNA promoters (CGI context on)."""
    return [
        "lncrna_name",
        "promoter_beta_mean",
        "promoter_beta_median",
        "promoter_beta_std",
        "promoter_beta_min",
        "promoter_beta_max",
        "promoter_beta_range",
        "promoter_n_probes",
        "promoter_frac_hypermeth",
        "promoter_frac_hypometh",
        "promoter_CGI_beta_mean",
        "promoter_CGI_n_probes",
        "promoter_shore_beta_mean",
    ]


def aggregate_to_lncrnas(
    sample_probes: pd.DataFrame,
    probe_reference: pd.DataFrame,
    lncrna_panel: Optional[List[str]] = None,
    hypermeth_threshold: float = HYPERMETH_THRESHOLD,
    hypometh_threshold: float = HYPOMETH_THRESHOLD,
) -> pd.DataFrame:
    """
    Aggregate probe beta values to lncRNA-level metrics.
    
    Same structure as gene aggregation but for lncRNAs.
    
    Args:
        sample_probes: DataFrame with probeID and beta
        probe_reference: Annotated probe reference with lncrna_promoter_genes
        lncrna_panel: List of lncRNA names to aggregate
        hypermeth_threshold: Beta threshold for hypermethylation
        hypometh_threshold: Beta threshold for hypomethylation
    
    Returns:
        DataFrame with one row per lncRNA
    """
    # Check if lncRNA annotations exist
    if "lncrna_promoter_genes" not in probe_reference.columns:
        return _empty_lncrna_agg_df(lncrna_panel)
    
    # Merge beta values with reference
    merged = sample_probes[["probeID", "beta"]].merge(
        probe_reference[["probeID", "lncrna_promoter_genes", "in_CGI", "CGI_context"]],
        on="probeID",
        how="inner",
    )
    
    merged = merged.dropna(subset=["beta"])
    
    if merged.empty:
        return _empty_lncrna_agg_df(lncrna_panel)
    
    # Explode lncRNA promoter genes
    prom_df = merged[merged["lncrna_promoter_genes"].apply(len) > 0].copy()
    prom_df = prom_df.explode("lncrna_promoter_genes")
    if lncrna_panel is not None:
        prom_df = prom_df[prom_df["lncrna_promoter_genes"].isin(lncrna_panel)]
    prom_df = prom_df.rename(columns={"lncrna_promoter_genes": "lncrna_name"})
    
    if prom_df.empty:
        return _empty_lncrna_agg_df(lncrna_panel)
    
    # Aggregate
    agg = _aggregate_region(
        prom_df,
        group_col="lncrna_name",
        prefix="promoter",
        hypermeth_threshold=hypermeth_threshold,
        hypometh_threshold=hypometh_threshold,
        include_cgi_context=True,
    )
    
    # Ensure all panel lncRNAs present (skip merge when panel is None or empty: use observed lncRNAs only)
    if lncrna_panel:
        result = pd.DataFrame({"lncrna_name": list(lncrna_panel)})
        result = result.merge(agg, on="lncrna_name", how="left")
    else:
        result = agg

    # Fill NaN counts
    count_cols = [c for c in result.columns if c.endswith("_n_probes")]
    for col in count_cols:
        result[col] = result[col].fillna(0).astype(int)
    
    return result


def _empty_lncrna_agg_df(lncrna_panel: Optional[List[str]]) -> pd.DataFrame:
    """
    Return an lncRNA aggregation frame with no data rows.

    When ``lncrna_panel`` is None or empty, returns zero rows but the same columns
    as a real aggregation so downstream code does not fail. When the panel is
    non-empty, returns one NaN-filled row per name (matches gene empty helper).
    """
    cols = _lncrna_agg_output_columns()
    if not lncrna_panel:
        return pd.DataFrame({c: [] for c in cols})

    n = len(lncrna_panel)
    return pd.DataFrame(
        {
            "lncrna_name": list(lncrna_panel),
            "promoter_beta_mean": [np.nan] * n,
            "promoter_beta_median": [np.nan] * n,
            "promoter_beta_std": [np.nan] * n,
            "promoter_beta_min": [np.nan] * n,
            "promoter_beta_max": [np.nan] * n,
            "promoter_beta_range": [np.nan] * n,
            "promoter_n_probes": [0] * n,
            "promoter_frac_hypermeth": [np.nan] * n,
            "promoter_frac_hypometh": [np.nan] * n,
            "promoter_CGI_beta_mean": [np.nan] * n,
            "promoter_CGI_n_probes": [0] * n,
            "promoter_shore_beta_mean": [np.nan] * n,
        }
    )


# =============================================================================
# CCRE-LEVEL AGGREGATION
# =============================================================================

def aggregate_to_ccres(
    sample_probes: pd.DataFrame,
    probe_reference: pd.DataFrame,
    ccre_ids: Optional[List[str]] = None,
    hypermeth_threshold: float = HYPERMETH_THRESHOLD,
    hypometh_threshold: float = HYPOMETH_THRESHOLD,
) -> pd.DataFrame:
    """
    Aggregate probe beta values to cCRE-level metrics.
    
    Args:
        sample_probes: DataFrame with probeID and beta
        probe_reference: Annotated probe reference with overlapping_ccres
        ccre_ids: Optional list of cCRE IDs to include (None = all)
        hypermeth_threshold: Beta threshold for hypermethylation
        hypometh_threshold: Beta threshold for hypomethylation
    
    Returns:
        DataFrame with one row per cCRE with methylation metrics
    """
    # Check if cCRE annotations exist
    if "overlapping_ccres" not in probe_reference.columns:
        return pd.DataFrame(columns=[
            "cCRE_id", "ccre_beta_mean", "ccre_beta_median",
            "ccre_n_probes", "ccre_frac_hypermeth", "ccre_frac_hypometh"
        ])
    
    # Merge beta values with reference
    cols_to_use = ["probeID", "overlapping_ccres", "in_CGI"]
    if "ccre_types" in probe_reference.columns:
        cols_to_use.append("ccre_types")
    
    merged = sample_probes[["probeID", "beta"]].merge(
        probe_reference[cols_to_use],
        on="probeID",
        how="inner",
    )
    
    merged = merged.dropna(subset=["beta"])
    
    if merged.empty:
        return pd.DataFrame(columns=[
            "cCRE_id", "ccre_beta_mean", "ccre_beta_median",
            "ccre_n_probes", "ccre_frac_hypermeth", "ccre_frac_hypometh"
        ])
    
    # Explode to one row per (probe, cCRE)
    ccre_df = merged[merged["overlapping_ccres"].apply(len) > 0].copy()
    ccre_df = ccre_df.explode("overlapping_ccres")
    ccre_df = ccre_df.rename(columns={"overlapping_ccres": "cCRE_id"})
    
    # Filter to specified cCREs if provided
    if ccre_ids is not None:
        ccre_df = ccre_df[ccre_df["cCRE_id"].isin(ccre_ids)]
    
    if ccre_df.empty:
        return pd.DataFrame(columns=[
            "cCRE_id", "ccre_beta_mean", "ccre_beta_median",
            "ccre_n_probes", "ccre_frac_hypermeth", "ccre_frac_hypometh"
        ])
    
    # Aggregate
    agg = _aggregate_region(
        ccre_df,
        group_col="cCRE_id",
        prefix="ccre",
        hypermeth_threshold=hypermeth_threshold,
        hypometh_threshold=hypometh_threshold,
        include_cgi_context=False,
    )
    
    # Add CGI overlap info per cCRE
    # A cCRE overlaps CGI if ANY of its probes are in CGI
    cgi_flags = ccre_df.groupby("cCRE_id")["in_CGI"].any()
    agg["ccre_CGI_overlap"] = agg["cCRE_id"].map(cgi_flags).fillna(False)
    
    # CGI-specific mean for cCREs that have CGI probes
    cgi_probes = ccre_df[ccre_df["in_CGI"] == True]
    if not cgi_probes.empty:
        cgi_means = cgi_probes.groupby("cCRE_id")["beta"].mean()
        cgi_counts = cgi_probes.groupby("cCRE_id").size()
        agg["ccre_CGI_beta_mean"] = agg["cCRE_id"].map(cgi_means)
        agg["ccre_CGI_n_probes"] = agg["cCRE_id"].map(cgi_counts).fillna(0).astype(int)
    else:
        agg["ccre_CGI_beta_mean"] = np.nan
        agg["ccre_CGI_n_probes"] = 0
    
    return agg



# =============================================================================
# ATAC-LEVEL AGGREGATION
# =============================================================================

def aggregate_to_atac(
    sample_probes: pd.DataFrame,
    probe_reference: pd.DataFrame,
    atac_peak_ids: Optional[List[str]] = None,
    hypermeth_threshold: float = HYPERMETH_THRESHOLD,
    hypometh_threshold: float = HYPOMETH_THRESHOLD,
) -> pd.DataFrame:
    """
    Aggregate probe beta values to ATAC-peak-level metrics.
    
    Args:
        sample_probes: DataFrame with probeID and beta
        probe_reference: Annotated probe reference with overlapping_atac_peaks
        atac_peak_ids: Optional list of ATAC peak IDs to include (None = all)
        hypermeth_threshold: Beta threshold for hypermethylation
        hypometh_threshold: Beta threshold for hypomethylation
    
    Returns:
        DataFrame with one row per ATAC peak with methylation metrics:
            - atac_peak_id: Peak identifier (chr:start-end format)
            - atac_beta_mean: Mean beta across probes in peak
            - atac_beta_median: Median beta across probes
            - atac_beta_std: Standard deviation of beta
            - atac_beta_min: Minimum beta value
            - atac_beta_max: Maximum beta value
            - atac_n_probes: Number of probes in peak
            - atac_frac_hypermeth: Fraction of probes above hypermeth threshold
            - atac_frac_hypometh: Fraction of probes below hypometh threshold
            - atac_CGI_overlap: Whether any probe in peak overlaps CGI
            - atac_CGI_beta_mean: Mean beta of CGI-overlapping probes (if any)
            - atac_CGI_n_probes: Count of CGI-overlapping probes
    
    Notes:
        - ATAC peaks are identified by composite IDs like 'chr22:41939978-41940479'
        - Probes can map to multiple peaks (exploded for aggregation)
        - CGI context is preserved when available in probe_reference
    """
    output_columns = [
        "atac_peak_id", "atac_beta_mean", "atac_beta_median",
        "atac_beta_std", "atac_beta_min", "atac_beta_max",
        "atac_n_probes", "atac_frac_hypermeth", "atac_frac_hypometh",
        "atac_CGI_overlap", "atac_CGI_beta_mean", "atac_CGI_n_probes",
    ]
    
    # Check if ATAC annotations exist
    if "overlapping_atac_peaks" not in probe_reference.columns:
        return pd.DataFrame(columns=output_columns)
    
    # Merge beta values with reference
    cols_to_use = ["probeID", "overlapping_atac_peaks"]
    if "in_CGI" in probe_reference.columns:
        cols_to_use.append("in_CGI")
    
    merged = sample_probes[["probeID", "beta"]].merge(
        probe_reference[cols_to_use],
        on="probeID",
        how="inner",
    )
    
    merged = merged.dropna(subset=["beta"])
    
    if merged.empty:
        return pd.DataFrame(columns=output_columns)
    
    # Handle overlapping_atac_peaks column - could be list or string representation
    def parse_atac_peaks(val):
        """Parse ATAC peak list from various formats."""
        if val is None:
            return []
        if isinstance(val, list):
            return val
        if isinstance(val, str):
            # Handle string representation of list
            if val.startswith("[") and val.endswith("]"):
                try:
                    import ast
                    return ast.literal_eval(val)
                except (ValueError, SyntaxError):
                    pass
            # Handle empty string
            if val.strip() == "" or val == "[]":
                return []
            # Handle comma-separated string
            return [p.strip() for p in val.split(",") if p.strip()]
        return []
    
    merged["overlapping_atac_peaks"] = merged["overlapping_atac_peaks"].apply(parse_atac_peaks)
    
    # Filter to rows with at least one ATAC peak
    atac_df = merged[merged["overlapping_atac_peaks"].apply(len) > 0].copy()
    
    if atac_df.empty:
        return pd.DataFrame(columns=output_columns)
    
    # Explode to one row per (probe, ATAC peak)
    atac_df = atac_df.explode("overlapping_atac_peaks")
    atac_df = atac_df.rename(columns={"overlapping_atac_peaks": "atac_peak_id"})
    
    # Filter to specified ATAC peaks if provided
    if atac_peak_ids is not None:
        atac_df = atac_df[atac_df["atac_peak_id"].isin(atac_peak_ids)]
    
    if atac_df.empty:
        return pd.DataFrame(columns=output_columns)
    
    # Aggregate core metrics
    agg = _aggregate_region(
        atac_df,
        group_col="atac_peak_id",
        prefix="atac",
        hypermeth_threshold=hypermeth_threshold,
        hypometh_threshold=hypometh_threshold,
        include_cgi_context=False,
    )
    
    # Add CGI overlap info per ATAC peak (if CGI info available)
    if "in_CGI" in atac_df.columns:
        # A peak overlaps CGI if ANY of its probes are in CGI
        cgi_flags = atac_df.groupby("atac_peak_id")["in_CGI"].any()
        agg["atac_CGI_overlap"] = agg["atac_peak_id"].map(cgi_flags).fillna(False)
        
        # CGI-specific mean for peaks that have CGI probes
        cgi_probes = atac_df[atac_df["in_CGI"] == True]
        if not cgi_probes.empty:
            cgi_means = cgi_probes.groupby("atac_peak_id")["beta"].mean()
            cgi_counts = cgi_probes.groupby("atac_peak_id").size()
            agg["atac_CGI_beta_mean"] = agg["atac_peak_id"].map(cgi_means)
            agg["atac_CGI_n_probes"] = agg["atac_peak_id"].map(cgi_counts).fillna(0).astype(int)
        else:
            agg["atac_CGI_beta_mean"] = np.nan
            agg["atac_CGI_n_probes"] = 0
    else:
        agg["atac_CGI_overlap"] = False
        agg["atac_CGI_beta_mean"] = np.nan
        agg["atac_CGI_n_probes"] = 0
    
    return agg

# =============================================================================
# COMPOSITE METHYLATION SCORES
# =============================================================================

def compute_promoter_methylation_score(
    gene_agg: pd.DataFrame,
    method: str = "mean",
    weights: Optional[Dict[str, float]] = None,
) -> pd.Series:
    """
    Compute a composite promoter methylation score.
    
    Args:
        gene_agg: Gene-level aggregation DataFrame
        method: Score computation method:
            - "mean": Simple mean beta
            - "cgi_weighted": Weight CGI probes higher
            - "custom": Use provided weights
        weights: Custom weights for "custom" method
    
    Returns:
        Series with methylation scores indexed by gene_name
    """
    if method == "mean":
        return gene_agg.set_index("gene_name")["promoter_beta_mean"]
    
    elif method == "cgi_weighted":
        # Weight CGI probes 2x compared to non-CGI
        cgi_mean = gene_agg["promoter_CGI_beta_mean"]
        overall_mean = gene_agg["promoter_beta_mean"]
        
        # Where CGI data exists, use weighted average
        has_cgi = cgi_mean.notna()
        score = overall_mean.copy()
        score[has_cgi] = 0.67 * cgi_mean[has_cgi] + 0.33 * overall_mean[has_cgi]
        
        return pd.Series(score.values, index=gene_agg["gene_name"])
    
    elif method == "custom":
        if weights is None:
            raise ValueError("weights required for custom method")
        
        # Compute weighted sum of specified columns
        score = pd.Series(0.0, index=gene_agg.index)
        for col, weight in weights.items():
            if col in gene_agg.columns:
                score += weight * gene_agg[col].fillna(0)
        
        return pd.Series(score.values, index=gene_agg["gene_name"])
    
    else:
        raise ValueError(f"Unknown method: {method}")


def classify_methylation_status(
    beta_mean: pd.Series,
    hypermeth_threshold: float = HYPERMETH_THRESHOLD,
    hypometh_threshold: float = HYPOMETH_THRESHOLD,
) -> pd.Series:
    """
    Classify methylation status based on mean beta.
    
    Returns:
        Series with values: "hypermethylated", "hypomethylated", "intermediate", "NA"
    """
    status = pd.Series("intermediate", index=beta_mean.index)
    status[beta_mean > hypermeth_threshold] = "hypermethylated"
    status[beta_mean < hypometh_threshold] = "hypomethylated"
    status[beta_mean.isna()] = "NA"
    
    return status


# =============================================================================
# DIFFERENTIAL METHYLATION HELPERS
# =============================================================================

def compute_delta_methylation(
    sample_agg: pd.DataFrame,
    reference_agg: pd.DataFrame,
    feature_col: str = "gene_name",
    beta_col: str = "promoter_beta_mean",
) -> pd.DataFrame:
    """
    Compute difference in methylation between sample and reference.
    
    Useful for tumor vs normal comparisons.
    
    Args:
        sample_agg: Sample aggregation (e.g., tumor)
        reference_agg: Reference aggregation (e.g., normal)
        feature_col: Column with feature IDs (gene_name, cCRE_id, etc.)
        beta_col: Column with beta values to compare
    
    Returns:
        DataFrame with delta methylation values
    """
    merged = sample_agg[[feature_col, beta_col]].merge(
        reference_agg[[feature_col, beta_col]],
        on=feature_col,
        suffixes=("_sample", "_reference"),
    )
    
    merged["delta_beta"] = (
        merged[f"{beta_col}_sample"] - merged[f"{beta_col}_reference"]
    )
    
    # Classify change direction
    merged["change_direction"] = "stable"
    merged.loc[merged["delta_beta"] > 0.1, "change_direction"] = "hypermethylated"
    merged.loc[merged["delta_beta"] < -0.1, "change_direction"] = "hypomethylated"
    
    return merged


# =============================================================================
# BATCH AGGREGATION
# =============================================================================

def aggregate_sample_batch_to_genes(
    sample_dfs: List[pd.DataFrame],
    sample_ids: List[str],
    probe_reference: pd.DataFrame,
    gene_panel: List[str],
) -> pd.DataFrame:
    """
    Aggregate multiple samples to gene-level, returning a genes x samples matrix.
    
    Args:
        sample_dfs: List of sample DataFrames with probeID and beta
        sample_ids: List of sample identifiers
        probe_reference: Annotated probe reference
        gene_panel: List of genes to aggregate
    
    Returns:
        DataFrame with genes as rows, samples as columns, promoter_beta_mean as values
    """
    results = {}
    
    for sample_df, sample_id in zip(sample_dfs, sample_ids):
        gene_agg = aggregate_to_genes(
            sample_df,
            probe_reference,
            gene_panel,
            include_gene_body=False,
        )
        results[sample_id] = gene_agg.set_index("gene_name")["promoter_beta_mean"]
    
    matrix = pd.DataFrame(results)
    matrix.index.name = "gene_name"
    
    return matrix


def aggregate_sample_batch_to_ccres(
    sample_dfs: List[pd.DataFrame],
    sample_ids: List[str],
    probe_reference: pd.DataFrame,
    ccre_ids: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Aggregate multiple samples to cCRE-level, returning a cCREs x samples matrix.
    
    Args:
        sample_dfs: List of sample DataFrames
        sample_ids: List of sample identifiers
        probe_reference: Annotated probe reference
        ccre_ids: Optional list of cCRE IDs to include
    
    Returns:
        DataFrame with cCREs as rows, samples as columns
    """
    results = {}
    
    for sample_df, sample_id in zip(sample_dfs, sample_ids):
        ccre_agg = aggregate_to_ccres(
            sample_df,
            probe_reference,
            ccre_ids=ccre_ids,
        )
        if not ccre_agg.empty:
            results[sample_id] = ccre_agg.set_index("cCRE_id")["ccre_beta_mean"]
    
    if not results:
        return pd.DataFrame()
    
    matrix = pd.DataFrame(results)
    matrix.index.name = "cCRE_id"
    
    return matrix
