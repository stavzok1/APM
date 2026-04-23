"""
Gene and lncRNA loading from GENCODE annotation.

Functions for:
- Loading gene annotations from GTF-derived CSV
- Filtering to specific gene panels
- Adding promoter coordinates
- Creating BED files
"""

import os
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd

from ..utils import (
    harmonize_chrom_column,
    compute_promoter_coords,
    compute_tss,
)


# =============================================================================
# LOADING
# =============================================================================

def load_genes(path: Path) -> pd.DataFrame:
    """
    Load gene annotations from GENCODE GTF-derived CSV.
    
    Expected columns: chrom/seqname, start, end, strand, gene_name, gene_id, feature, gene_type
    """
    df = pd.read_csv(path)
    df, _ = harmonize_chrom_column(df)
    
    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    return df


def load_lncrnas(path: Path) -> pd.DataFrame:
    """
    Load lncRNA annotations from CSV.
    Same structure as gene annotations but filtered to lncRNAs.
    """
    df = pd.read_csv(path)
    df, _ = harmonize_chrom_column(df)
    
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    return df


# =============================================================================
# FILTERING
# =============================================================================

def filter_genes_by_names(
    genes: pd.DataFrame,
    gene_names: List[str],
    feature_type: str = "gene",
) -> pd.DataFrame:
    """
    Filter genes to those in the specified name list.
    
    Args:
        genes: Full gene DataFrame
        gene_names: List of gene symbols to keep
        feature_type: GTF feature type to filter to (default: "gene")
    
    Returns:
        Filtered DataFrame
    """
    mask = genes["gene_name"].isin(gene_names)
    
    if "feature" in genes.columns and feature_type:
        mask &= genes["feature"] == feature_type
    
    return genes[mask].copy()


def filter_lncrnas(
    genes: pd.DataFrame,
    feature_type: str = "gene",
) -> pd.DataFrame:
    """
    Filter to lncRNA genes only.
    
    Args:
        genes: Full gene DataFrame
        feature_type: GTF feature type to filter to (default: "gene")
    
    Returns:
        DataFrame with only lncRNA genes
    """
    mask = genes["gene_type"] == "lncRNA"
    
    if "feature" in genes.columns and feature_type:
        mask &= genes["feature"] == feature_type
    
    return genes[mask].copy()


# =============================================================================
# PROMOTER COORDINATES
# =============================================================================

def add_promoter_columns(
    genes: pd.DataFrame,
    upstream_bp: int = 2000,
    downstream_bp: int = 500,
) -> pd.DataFrame:
    """
    Add promoter coordinate columns (prom_start, prom_end) to gene DataFrame.
    
    Promoter is defined relative to TSS:
    - For + strand: [TSS - upstream, TSS + downstream]
    - For - strand: [TSS - downstream, TSS + upstream]
    
    Also adds 'tss' column.
    """
    genes = genes.copy()
    
    # Compute TSS
    genes["tss"] = compute_tss(genes["start"], genes["end"], genes["strand"])
    
    # Compute promoter coordinates
    prom_start, prom_end = compute_promoter_coords(
        genes["start"],
        genes["end"],
        genes["strand"],
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
    )
    
    genes["prom_start"] = prom_start.astype("Int64")
    genes["prom_end"] = prom_end.astype("Int64")
    
    return genes


def add_tss_window(
    genes: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    """
    Add TSS-centered window columns (win_start, win_end) to gene DataFrame.
    """
    genes = genes.copy()
    
    if "tss" not in genes.columns:
        genes["tss"] = compute_tss(genes["start"], genes["end"], genes["strand"])
    
    genes["win_start"] = (genes["tss"] - window_bp).clip(lower=0).astype("Int64")
    genes["win_end"] = (genes["tss"] + window_bp).astype("Int64")
    
    return genes


# =============================================================================
# OUTPUT FORMATS
# =============================================================================

def create_genes_bed(
    genes: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Create BED file from gene DataFrame.
    
    BED format is 0-based, half-open: [start, end)
    """
    df_bed = genes[["chrom", "start", "end", "gene_name", "strand"]].copy()
    
    # Convert to 0-based
    df_bed["start"] = df_bed["start"] - 1
    
    # Ensure no header, tab-separated
    df_bed.to_csv(output_path, sep="\t", header=False, index=False)
    print(f"Created BED file: {output_path}")


def save_gene_tables(
    genes: pd.DataFrame,
    primary_genes: List[str],
    lncrna_names: List[str],
    output_dir: Path,
) -> None:
    """
    Save filtered gene tables to CSV.
    
    Creates:
    - primary_genes_all_features.csv: All GTF features for primary genes
    - primary_genes_only.csv: Gene-level entries only for primary genes
    - lncRNAs_genes_all_features.csv: All features for matched lncRNAs
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Primary genes - all features
    primary_all = genes[genes["gene_name"].isin(primary_genes)]
    primary_all.to_csv(output_dir / "primary_genes_all_features.csv", index=False)
    
    # Primary genes - gene level only
    primary_gene = primary_all[primary_all.get("feature", "gene") == "gene"]
    primary_gene.to_csv(output_dir / "primary_genes_only.csv", index=False)
    
    # lncRNAs
    lncrna_all = genes[genes["gene_name"].isin(lncrna_names)]
    lncrna_all.to_csv(output_dir / "lncRNAs_genes_all_features.csv", index=False)
    
    print(f"Saved gene tables to {output_dir}")


# =============================================================================
# HARMONIZATION ACROSS MULTIPLE DATAFRAMES
# =============================================================================

def harmonize_multiple_dfs(
    dfs: List[pd.DataFrame],
    paths: Optional[List[Path]] = None,
    save_if_changed: bool = False,
) -> Tuple[List[pd.DataFrame], List[bool]]:
    """
    Harmonize chromosome columns across multiple DataFrames.
    
    Args:
        dfs: List of DataFrames to harmonize
        paths: Optional list of paths for saving modified DataFrames
        save_if_changed: Whether to save back to paths if changes were made
    
    Returns:
        Tuple of (harmonized DataFrames, list of was_changed flags)
    """
    results = []
    changes = []
    
    for i, df in enumerate(dfs):
        df_harmonized, was_changed = harmonize_chrom_column(df)
        results.append(df_harmonized)
        changes.append(was_changed)
        
        if save_if_changed and was_changed and paths and i < len(paths):
            df_harmonized.to_csv(paths[i], index=False)
            print(f"Saved harmonized DataFrame to {paths[i]}")
    
    return results, changes
