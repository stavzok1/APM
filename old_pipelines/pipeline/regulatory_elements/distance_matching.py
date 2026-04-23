"""
cCRE to gene distance matching.

Matches regulatory elements to genes based on TSS distance,
organizing results into distance tiers.
"""

import os
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from ..utils import (
    add_bins,
    assign_distance_tier,
    tss_to_interval_distance,
    as_categorical,
    compute_tss,
)
from ..config import THRESHOLDS


# =============================================================================
# CORE MATCHING
# =============================================================================

def match_ccres_to_genes(
    genes: pd.DataFrame,
    ccres: pd.DataFrame,
    window_bp: int = 1_000_000,
    bin_size: int = 100_000,
    tier_edges: List[int] = None,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Match cCREs to genes based on distance from TSS.
    
    Uses binning for efficient prefiltering, then computes exact distances.
    
    Args:
        genes: Gene DataFrame with chrom, start, end, strand, gene_name
        ccres: cCRE DataFrame with chrom, start, end, ENCODE_id, type
        window_bp: Maximum distance from TSS to consider
        bin_size: Bin size for coarse prefiltering
        tier_edges: Distance tier boundaries in bp
        tier_labels: Labels for each tier
    
    Returns:
        DataFrame with all gene-cCRE pairs within window, with distances and tiers
    """
    # Use defaults from config if not provided
    if tier_edges is None:
        tier_edges = THRESHOLDS.tier_edges_bp
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels
    
    # Prepare genes
    g = genes[["chrom", "start", "end", "strand", "gene_name"]].copy()
    g["start"] = g["start"].astype(np.int64)
    g["end"] = g["end"].astype(np.int64)
    g = as_categorical(g, ["chrom", "gene_name", "strand"])
    
    # Compute TSS and window
    g["tss"] = compute_tss(g["start"], g["end"], g["strand"]).astype(np.int64)
    g["win_start"] = (g["tss"] - window_bp).clip(lower=0).astype(np.int64)
    g["win_end"] = (g["tss"] + window_bp).astype(np.int64)
    
    # Prepare cCREs
    c = ccres[["chrom", "start", "end", "ENCODE_id", "cCRE_id", "type", "raw_type"]].copy()
    c["start"] = c["start"].astype(np.int64)
    c["end"] = c["end"].astype(np.int64)
    c["center"] = ((c["start"] + c["end"]) // 2).astype(np.int64)
    c = as_categorical(c, ["chrom", "type"])
    
    # Add bins for prefiltering
    g_bins = add_bins(g, "win_start", "win_end", bin_size)
    c_bins = add_bins(c, "start", "end", bin_size)
    
    # Explode bins and join
    gb = g_bins.explode("bins")
    cb = c_bins.explode("bins")
    
    pref = (
        gb[["gene_name", "chrom", "win_start", "win_end", "tss", "bins"]]
        .merge(
            cb[["chrom", "start", "end", "center", "cCRE_id", "ENCODE_id", "type", "raw_type", "bins"]],
            on=["chrom", "bins"],
            how="inner",
            copy=False,
        )
    )
    
    # Tighten with real window overlap
    pref = pref[
        (pref["end"] >= pref["win_start"]) & 
        (pref["start"] <= pref["win_end"])
    ].copy()
    
    # Compute exact distance from TSS to element interval
    pref["dist_to_tss"] = tss_to_interval_distance(
        pref["tss"].to_numpy(),
        pref["start"].to_numpy(),
        pref["end"].to_numpy(),
    ).astype(np.int64)
    
    # Keep within window by exact distance
    pref = pref[pref["dist_to_tss"] <= window_bp].copy()
    
    # Assign distance tiers
    pref["tier"] = assign_distance_tier(pref["dist_to_tss"], tier_edges, tier_labels)
    
    # Select final columns
    pair_df = pref[[
        "gene_name", "cCRE_id", "ENCODE_id", "type", "raw_type",
        "chrom", "start", "end", "center", "tss", "dist_to_tss", "tier",
    ]].copy()
    
    # Drop duplicates from bin explosion
    pair_df = pair_df.drop_duplicates(
        subset=["gene_name", "cCRE_id", "ENCODE_id"]
    ).reset_index(drop=True)
    
    print(f"Matched {len(pair_df)} gene-cCRE pairs within {window_bp}bp")
    
    return pair_df


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def save_distance_matching(
    output_dir: Path,
    pair_df: pd.DataFrame,
    prefix: str = "",
) -> None:
    """
    Save gene-cCRE distance matching results.
    
    Args:
        output_dir: Directory to save files
        pair_df: DataFrame from match_ccres_to_genes
        prefix: Optional prefix for filenames (e.g., "lncRNA_")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filename = f"{prefix}gene_to_elements.csv" if prefix else "gene_to_elements.csv"
    pair_df.to_csv(output_dir / filename, index=False)
    
    print(f"Saved distance matching to {output_dir / filename}")



def save_all_matching_outputs(
    output_dir: Path,
    pair_df: pd.DataFrame,
    elem_focus: pd.DataFrame,
    gene_summary: pd.DataFrame,
    label: str = "coding",
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # pairs
    save_distance_matching(output_dir, pair_df, prefix=f"{label}_")

    # tables
    elem_focus.to_csv(output_dir / f"{label}_element_focus.csv", index=False)
    gene_summary.to_csv(output_dir / f"{label}_gene_summary.csv", index=False)

    print(f"Saved all matching outputs for {label} to {output_dir}")


    


# =============================================================================
# AGGREGATION BY GENE
# =============================================================================

def aggregate_ccres_per_gene(
    pair_df: pd.DataFrame,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Create gene-level summary of cCRE counts and IDs per type/tier.
    
    Returns wide DataFrame with columns like:
    - "pELS | 0-100kb | count"
    - "pELS | 0-100kb | ids"
    """
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels
    
    # Counts per (gene, type, tier)
    cnt = (
        pair_df
        .groupby(["gene_name", "type", "tier"], observed=True)["cCRE_id"]
        .nunique()
        .rename("count")
        .reset_index()
    )
    
    cnt_wide = (
        cnt
        .assign(col=lambda d: d["type"].astype(str) + " | " + d["tier"].astype(str) + " | count")
        .pivot(index="gene_name", columns="col", values="count")
        .reset_index()
    )
    
    # IDs per (gene, type, tier)
    ids = (
        pair_df
        .groupby(["gene_name", "type", "tier"], observed=True)["cCRE_id"]
        .apply(lambda s: ",".join(sorted(pd.unique(s.astype(str)))))
        .rename("ids")
        .reset_index()
    )
    
    ids_wide = (
        ids
        .assign(col=lambda d: d["type"].astype(str) + " | " + d["tier"].astype(str) + " | ids")
        .pivot(index="gene_name", columns="col", values="ids")
        .reset_index()
    )
    
    # Merge counts and IDs
    gene_summary = cnt_wide.merge(ids_wide, on="gene_name", how="outer")
    
    # Fill NaN appropriately
    num_cols = [c for c in gene_summary.columns if c.endswith(" | count")]
    id_cols = [c for c in gene_summary.columns if c.endswith(" | ids")]
    
    for c in num_cols:
        gene_summary[c] = pd.to_numeric(gene_summary[c], errors="coerce").fillna(0).astype(int)
    for c in id_cols:
        gene_summary[c] = gene_summary[c].astype(object).fillna("")
    
    return gene_summary


# =============================================================================
# AGGREGATION BY ELEMENT
# =============================================================================

def aggregate_genes_per_ccre(
    pair_df: pd.DataFrame,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Create element-level summary with genes organized by distance tier.
    
    Returns DataFrame with columns:
    - cCRE_id, ENCODE_id, type, raw_type, chrom, start, end
    - One column per tier containing comma-separated gene names
    - min_dist_to_any_gene
    - genes_by_exact_dist (format: "gene1:123,gene2:456")
    """
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels
    
    # Genes per tier
    agg_genes = (
        pair_df
        .groupby(["cCRE_id", "ENCODE_id", "type", "raw_type", "chrom", "start", "end", "tier"], observed=True)["gene_name"]
        .apply(lambda s: ",".join(sorted(pd.unique(s.astype(str)))))
        .reset_index()
    )
    
    elem_focus = (
        agg_genes
        .pivot(
            index=["cCRE_id", "ENCODE_id", "type", "raw_type", "chrom", "start", "end"],
            columns="tier",
            values="gene_name",
        )
        .reset_index()
    )
    
    # Ensure tier columns are object type before fillna
    tier_cols = [c for c in elem_focus.columns if c in tier_labels]
    for c in tier_cols:
        elem_focus[c] = elem_focus[c].astype(object).fillna("")
    
    # Min distance to any gene
    min_dist = (
        pair_df
        .groupby("cCRE_id", observed=True)["dist_to_tss"]
        .min()
        .rename("min_dist_to_any_gene")
        .reset_index()
    )
    elem_focus = elem_focus.merge(min_dist, on="cCRE_id", how="left")
    
    # Exact "gene:dist" string
    exact_dist = (
        pair_df
        .sort_values(["cCRE_id", "dist_to_tss", "gene_name"])
        .groupby("cCRE_id", observed=True)
        .apply(lambda d: ",".join(f"{g}:{int(dd)}" for g, dd in zip(d["gene_name"], d["dist_to_tss"])))
        .rename("genes_by_exact_dist")
        .reset_index()
    )
    elem_focus = elem_focus.merge(exact_dist, on="cCRE_id", how="left")
    
    return elem_focus


def build_distance_matrix(pair_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build gene × cCRE distance matrix.
    
    Returns pivot table with genes as rows, cCREs as columns,
    minimum distance as values.
    """
    matrix = (
        pair_df
        .groupby(["gene_name", "cCRE_id"], observed=True)["dist_to_tss"]
        .min()
        .unstack("cCRE_id")
        .astype("float64")
    )
    
    # Order columns by global min distance
    elem_order = pair_df.groupby("cCRE_id", observed=True)["dist_to_tss"].min().sort_values().index
    matrix = matrix.reindex(columns=elem_order)
    
    return matrix


