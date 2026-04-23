"""
lncRNA-to-gene proximity matching.

Identifies lncRNAs within a specified window of protein-coding genes,
computing exact distances and creating bidirectional mappings.
"""

import os
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from ..utils import min_interval_distance


# =============================================================================
# CORE MATCHING LOGIC
# =============================================================================

def match_lncrnas_to_genes(
    genes: pd.DataFrame,
    lncrnas: pd.DataFrame,
    window_bp: int = 1_000_000,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Find lncRNAs within a window of each gene.
    
    Uses optimized per-chromosome processing with searchsorted for efficiency.
    
    Args:
        genes: Gene DataFrame with chrom, start, end, strand, gene_name
        lncrnas: lncRNA DataFrame with same columns
        window_bp: Window size in base pairs around each gene
    
    Returns:
        Tuple of:
        - pairs_df: All gene-lncRNA pairs with distances
        - genes_with_lists: Genes with lncRNAs_within_{window} column
        - lncrnas_with_lists: lncRNAs with genes_within_{window} column
    """
    # Defensive copies & minimal selection
    g = genes[["chrom", "start", "end", "strand", "gene_name"]].copy()
    l = lncrnas[["chrom", "start", "end", "strand", "gene_name"]].copy()
    
    # Normalize types
    for df in (g, l):
        df["start"] = df["start"].astype(np.int64)
        df["end"] = df["end"].astype(np.int64)
        df["chrom"] = df["chrom"].astype("category")
    
    # Precompute expanded gene windows
    g["w_start"] = (g["start"] - window_bp).clip(lower=0)
    g["w_end"] = g["end"] + window_bp
    
    pairs = []
    
    # Process per chromosome
    chroms = np.intersect1d(
        g["chrom"].cat.categories.astype(str),
        l["chrom"].cat.categories.astype(str),
    )
    
    for c in chroms:
        g_c = g[g["chrom"].astype(str) == c].sort_values("w_start").reset_index(drop=True)
        l_c = l[l["chrom"].astype(str) == c].sort_values("start").reset_index(drop=True)
        
        if g_c.empty or l_c.empty:
            continue
        
        # Arrays for lncRNA positions
        l_start = l_c["start"].to_numpy()
        l_end = l_c["end"].to_numpy()
        
        # Gene window positions
        w_start = g_c["w_start"].to_numpy()
        w_end = g_c["w_end"].to_numpy()
        
        # searchsorted for efficient candidate finding
        idx_end = np.searchsorted(l_start, w_end, side="right")
        
        for i in range(g_c.shape[0]):
            upto = idx_end[i]
            if upto == 0:
                continue
            
            # Candidates by start position
            cand_idx = np.arange(upto, dtype=np.int64)
            
            # Filter by end >= w_start[i]
            mask = l_end[cand_idx] >= w_start[i]
            if not np.any(mask):
                continue
            
            sel = cand_idx[mask]
            if sel.size == 0:
                continue
            
            # Compute minimal distance between gene interval and lncRNA interval
            dist = min_interval_distance(
                g_c.at[i, "start"], g_c.at[i, "end"],
                l_c["start"].to_numpy()[sel], l_c["end"].to_numpy()[sel],
            )
            
            # Build result rows
            sub = pd.DataFrame({
                "chrom": c,
                "gene_name": g_c.at[i, "gene_name"],
                "gene_strand": g_c.at[i, "strand"],
                "gene_start": g_c.at[i, "start"],
                "gene_end": g_c.at[i, "end"],
                "lncRNA_name": l_c["gene_name"].to_numpy()[sel],
                "lncRNA_start": l_c["start"].to_numpy()[sel],
                "lncRNA_end": l_c["end"].to_numpy()[sel],
                "lncRNA_strand": l_c["strand"].to_numpy()[sel],
                "min_distance_bp": dist,
            })
            pairs.append(sub)
    
    # Combine results
    if len(pairs) == 0:
        pairs_df = pd.DataFrame(columns=[
            "chrom", "gene_name", "gene_start", "gene_end",
            "lncRNA_name", "lncRNA_start", "lncRNA_end", "min_distance_bp",
        ])
    else:
        pairs_df = pd.concat(pairs, ignore_index=True)
    
    # Deduplicate
    pairs_df = pairs_df.drop_duplicates(
        ["chrom", "gene_name", "lncRNA_name", "gene_start", "gene_end", "lncRNA_start", "lncRNA_end"]
    ).reset_index(drop=True)
    
    # Build per-gene aggregation
    genes_with_lists = _aggregate_lncrnas_per_gene(genes, pairs_df, window_bp)
    
    # Build per-lncRNA aggregation
    lncrnas_with_lists = _aggregate_genes_per_lncrna(lncrnas, pairs_df, window_bp)
    
    return pairs_df, genes_with_lists, lncrnas_with_lists


def _aggregate_lncrnas_per_gene(
    genes: pd.DataFrame,
    pairs_df: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    """Add list of lncRNAs within window to each gene."""
    col_name = f"lncRNAs_within_{window_bp // 1000}kb" if window_bp >= 1000 else f"lncRNAs_within_{window_bp}bp"
    
    lnc_per_gene = (
        pairs_df.groupby("gene_name", sort=False)["lncRNA_name"]
        .agg(lambda s: sorted(pd.unique(s)))
        .rename(col_name)
        .reset_index()
    )
    
    genes_out = genes.copy()
    genes_out = genes_out.merge(lnc_per_gene, on="gene_name", how="left")
    genes_out[col_name] = genes_out[col_name].apply(
        lambda x: x if isinstance(x, list) else []
    )
    
    return genes_out


def _aggregate_genes_per_lncrna(
    lncrnas: pd.DataFrame,
    pairs_df: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    """Add list of genes within window to each lncRNA, filtering to those with matches."""
    col_name = f"genes_within_{window_bp // 1000}kb" if window_bp >= 1000 else f"genes_within_{window_bp}bp"
    
    genes_per_lnc = (
        pairs_df.groupby("lncRNA_name", sort=False)["gene_name"]
        .agg(lambda s: sorted(pd.unique(s)))
        .rename(col_name)
        .reset_index()
    )
    
    lncrnas_out = lncrnas.copy()
    lncrnas_out = lncrnas_out.merge(
        genes_per_lnc,
        left_on="gene_name",
        right_on="lncRNA_name",
        how="left",
    )
    
    # Clean up merge artifacts
    if "lncRNA_name" in lncrnas_out.columns:
        lncrnas_out.drop(columns=["lncRNA_name"], inplace=True)
    
    lncrnas_out[col_name] = lncrnas_out[col_name].apply(
        lambda x: x if isinstance(x, list) else []
    )
    
    # Add count column
    lncrnas_out["num_genes"] = lncrnas_out[col_name].apply(len)
    
    # Filter to lncRNAs with at least one gene match
    lncrnas_out = lncrnas_out[lncrnas_out["num_genes"] > 0].copy()
    
    return lncrnas_out


# =============================================================================
# OUTPUT
# =============================================================================

def save_lncrna_matching(
    output_dir: Path,
    pairs_df: pd.DataFrame,
    genes_with_lists: pd.DataFrame,
    lncrnas_with_lists: pd.DataFrame,
    window_bp: int,
) -> None:
    """
    Save lncRNA matching results to CSV files.
    
    Creates:
    - genes_lncRNAs_{window}bp_distances.csv: All pairs with distances
    - genes_with_lncRNAs_{window}bp.csv: Gene-centric view
    - lncRNAs_with_genes_{window}bp.csv: lncRNA-centric view
    - lncRNAs_names_with_genes_{window}bp_names.csv: Just lncRNA names
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine column name based on window
    window_str = f"{window_bp // 1000}kb" if window_bp >= 1000 else f"{window_bp}bp"
    lnc_col = f"lncRNAs_within_{window_str}"
    gene_col = f"genes_within_{window_str}"
    
    # Pairs table
    pairs_df.to_csv(
        output_dir / f"genes_lncRNAs_{window_bp}bp_distances.csv",
        index=False,
    )
    
    # Genes with lncRNA lists
    gene_cols = ["chrom", "start", "end", "gene_name", "gene_id"]
    if lnc_col in genes_with_lists.columns:
        gene_cols.append(lnc_col)
    genes_subset = genes_with_lists[[c for c in gene_cols if c in genes_with_lists.columns]]
    genes_subset.to_csv(
        output_dir / f"genes_with_lncRNAs_{window_bp}bp.csv",
        index=False,
    )
    
    # lncRNAs with gene lists
    lnc_cols = ["chrom", "start", "end", "gene_name", "gene_id", "strand", "num_genes"]
    if gene_col in lncrnas_with_lists.columns:
        lnc_cols.insert(-1, gene_col)
    lncrnas_subset = lncrnas_with_lists[[c for c in lnc_cols if c in lncrnas_with_lists.columns]]
    lncrnas_subset.to_csv(
        output_dir / f"lncRNAs_with_genes_{window_bp}bp.csv",
        index=False,
    )
    
    # Just names
    lncrnas_with_lists["gene_name"].to_csv(
        output_dir / f"lncRNAs_names_with_genes_{window_bp}bp_names.csv",
        index=False,
    )
    
    print(f"Saved lncRNA matching results to {output_dir}")


def get_matched_lncrna_names(lncrnas_with_lists: pd.DataFrame) -> List[str]:
    """Extract list of lncRNA names that have gene matches."""
    return lncrnas_with_lists["gene_name"].tolist()
