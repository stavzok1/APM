"""
ATAC peak to gene matching.

Provides:
- TSS distance matching with distance tiers
- Gene body overlap detection
- Combined gene_links dict building
"""

from typing import List, Dict, Any, Optional

import numpy as np
import pandas as pd

from ..schemas import (
    empty_atac_gene_link_entry,
    empty_atac_body_overlap_entry,
    
)

from ..utils import (
    compute_interval_overlap,
)
from ..config import THRESHOLDS


# =============================================================================
# DISTANCE CALCULATIONS
# =============================================================================

def _tss_to_interval_distance(
    tss: np.ndarray,
    interval_start: np.ndarray,
    interval_end: np.ndarray,
) -> np.ndarray:
    """
    Distance from TSS point to an interval.
    Returns 0 if TSS is inside the interval.
    """
    inside = (tss >= interval_start) & (tss <= interval_end)
    dist_to_start = np.abs(tss - interval_start)
    dist_to_end = np.abs(tss - interval_end)
    return np.where(inside, 0, np.minimum(dist_to_start, dist_to_end))


def _compute_tss(start: pd.Series, end: pd.Series, strand: pd.Series) -> pd.Series:
    """Compute TSS position based on strand."""
    return np.where(strand.astype(str) == "+", start, end)


# =============================================================================
# DISTANCE TIER ASSIGNMENT
# =============================================================================

def _assign_distance_tier(
    distances: pd.Series,
    tier_edges: List[int],
    tier_labels: List[str],
) -> pd.Categorical:
    """Assign distance tier labels to a series of distances."""
    bins = [-0.1] + tier_edges[1:]
    return pd.cut(
        distances,
        bins=bins,
        labels=tier_labels,
        include_lowest=True,
        right=True,
        ordered=True,
    )


# =============================================================================
# GENE BODY OVERLAP
# =============================================================================

def _classify_overlap_type(
    peak_start: int,
    peak_end: int,
    gene_start: int,
    gene_end: int,
    tss: int,
    strand: str,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 500,
) -> Optional[str]:
    """
    Classify the type of overlap between peak and gene.
    
    Returns one of:
    - "promoter": overlaps TSS +/- window
    - "gene_body": overlaps gene body but not promoter
    - None: no overlap
    """
    # Check any overlap first
    if peak_end <= gene_start or peak_start >= gene_end:
        return None
    
    # Define promoter region
    if strand == "+":
        prom_start = max(0, tss - promoter_upstream)
        prom_end = tss + promoter_downstream
    else:
        prom_start = tss - promoter_downstream
        prom_end = tss + promoter_upstream
    
    # Check promoter overlap
    if not (peak_end <= prom_start or peak_start >= prom_end):
        return "promoter"
    
    return "gene_body"


def _compute_body_overlap(
    peak_start: int,
    peak_end: int,
    gene_start: int,
    gene_end: int,
    tss: int,
    strand: str,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 500,
) -> Dict[str, Any]:
    """
    Compute detailed gene body overlap for a peak-gene pair.
    """
    # Get basic overlap stats
    overlap = compute_interval_overlap(peak_start, peak_end, gene_start, gene_end)
    
    # Rename for gene context
    result = {
        "overlaps": overlap["overlaps"],
        "overlap_bp": overlap["overlap_bp"],
        "overlap_interval": overlap["overlap_interval"],
        "overlap_frac_of_peak": overlap["overlap_frac_of_a"],
        "overlap_frac_of_gene": overlap["overlap_frac_of_b"],
        "overlap_type": None,
    }
    
    if overlap["overlaps"]:
        result["overlap_type"] = _classify_overlap_type(
            peak_start, peak_end,
            gene_start, gene_end,
            tss, strand,
            promoter_upstream, promoter_downstream,
        )
    
    return result


# =============================================================================
# CORE MATCHING FUNCTION
# =============================================================================

def match_peaks_to_genes(
    peaks: pd.DataFrame,
    genes: pd.DataFrame,
    window_bp: int = 1_000_000,
    bin_size: int = 100_000,
    tier_edges: Optional[List[int]] = None,
    tier_labels: Optional[List[str]] = None,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 500,
) -> pd.DataFrame:
    """
    Match ATAC peaks to genes by TSS distance and compute body overlaps.
    
    Args:
        peaks: ATAC peaks DataFrame with peak_id, chrom, start, end
        genes: Gene DataFrame with gene_name, gene_id, chrom, start, end, strand
        window_bp: Maximum distance from TSS to consider
        bin_size: Bin size for coarse prefiltering
        tier_edges: Distance tier boundaries
        tier_labels: Labels for tiers
        promoter_upstream: bp upstream of TSS for promoter definition
        promoter_downstream: bp downstream of TSS for promoter definition
    
    Returns:
        DataFrame with all peak-gene pairs within window
    """
    tier_edges = tier_edges or THRESHOLDS.tier_edges_bp
    tier_labels = tier_labels or THRESHOLDS.tier_labels
    
    # Prepare genes
    gene_cols = ["chrom", "start", "end", "strand", "gene_name"]
    if "gene_id" in genes.columns:
        gene_cols.append("gene_id")
    if "gene_type" in genes.columns:
        gene_cols.append("gene_type")
    
    g = genes[gene_cols].copy()
    g["start"] = g["start"].astype(np.int64)
    g["end"] = g["end"].astype(np.int64)
    g["tss"] = _compute_tss(g["start"], g["end"], g["strand"]).astype(np.int64)
    g["win_start"] = (g["tss"] - window_bp).clip(lower=0).astype(np.int64)
    g["win_end"] = (g["tss"] + window_bp).astype(np.int64)
    
    # Prepare peaks
    p = peaks[["peak_id", "chrom", "start", "end"]].copy()
    p["start"] = p["start"].astype(np.int64)
    p["end"] = p["end"].astype(np.int64)
    p = p.rename(columns={"start": "peak_start", "end": "peak_end"})
    
    # Add bins for prefiltering
    g["bin_start"] = (g["win_start"] // bin_size).astype(np.int64)
    g["bin_end"] = (g["win_end"] // bin_size).astype(np.int64)
    g["bins"] = [list(range(s, e + 1)) for s, e in zip(g["bin_start"], g["bin_end"])]
    
    p["bin_start"] = (p["peak_start"] // bin_size).astype(np.int64)
    p["bin_end"] = (p["peak_end"] // bin_size).astype(np.int64)
    p["bins"] = [list(range(s, e + 1)) for s, e in zip(p["bin_start"], p["bin_end"])]
    
    # Explode bins and merge
    gb = g.explode("bins")
    pb = p.explode("bins")
    
    merge_cols = ["gene_name", "chrom", "start", "end", "strand", "tss", "win_start", "win_end", "bins"]
    if "gene_id" in g.columns:
        merge_cols.insert(1, "gene_id")
    if "gene_type" in g.columns:
        merge_cols.insert(2 if "gene_id" in g.columns else 1, "gene_type")
    
    pref = (
        gb[merge_cols]
        .merge(
            pb[["peak_id", "chrom", "peak_start", "peak_end", "bins"]],
            on=["chrom", "bins"],
            how="inner",
        )
    )
    
    # Filter to real overlaps with window
    pref = pref[
        (pref["peak_end"] >= pref["win_start"]) & 
        (pref["peak_start"] <= pref["win_end"])
    ].copy()
    
    # Compute exact TSS distance
    pref["dist_to_tss"] = _tss_to_interval_distance(
        pref["tss"].to_numpy(),
        pref["peak_start"].to_numpy(),
        pref["peak_end"].to_numpy(),
    ).astype(np.int64)
    
    # Filter by distance
    pref = pref[pref["dist_to_tss"] <= window_bp].copy()
    
    # Assign distance tiers
    pref["tier"] = _assign_distance_tier(pref["dist_to_tss"], tier_edges, tier_labels)
    
    # Compute body overlaps
    def compute_overlap_row(row):
        return _compute_body_overlap(
            int(row["peak_start"]), int(row["peak_end"]),
            int(row["start"]), int(row["end"]),
            int(row["tss"]), str(row["strand"]),
            promoter_upstream, promoter_downstream,
        )
    
    print(f"Computing body overlaps for {len(pref)} peak-gene pairs...")
    pref["body_overlap"] = pref.apply(compute_overlap_row, axis=1)
    pref["body_overlap"].fillna(empty_atac_body_overlap_entry, inplace=True)
    
    # Select final columns
    final_cols = [
        "peak_id", "gene_name",
        "chrom", "peak_start", "peak_end",
        "start", "end", "strand", "tss",
        "dist_to_tss", "tier", "body_overlap",
    ]
    if "gene_id" in pref.columns:
        final_cols.insert(2, "gene_id")
    if "gene_type" in pref.columns:
        final_cols.insert(3 if "gene_id" in pref.columns else 2, "gene_type")
    
    pair_df = pref[final_cols].copy()
    
    # Rename gene coords
    pair_df = pair_df.rename(columns={
        "start": "gene_start",
        "end": "gene_end",
    })
    
    # Drop duplicates from bin explosion
    pair_df = pair_df.drop_duplicates(
        subset=["peak_id", "gene_name"]
    ).reset_index(drop=True)
    
    print(f"Matched {len(pair_df)} peak-gene pairs within {window_bp}bp")
    n_overlaps = pair_df["body_overlap"].apply(lambda x: x.get("overlaps", False)).sum()
    print(f"  Including {n_overlaps} pairs with body overlap")
    
    return pair_df


# =============================================================================
# GENE LINKS DICT BUILDING
# =============================================================================

def build_gene_links(
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build gene_links dict column from peak-gene pairs.
    
    Args:
        pair_df: Output from match_peaks_to_genes
    
    Returns:
        DataFrame with [peak_id, gene_links] where gene_links is nested dict
    """
    def build_links_for_peak(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        links = {}
        for _, row in df_sub.iterrows():
            links[row["gene_name"]] = {
                "gene_id": row.get("gene_id"),
                "gene_type": row.get("gene_type"),
                "dist_to_tss": int(row["dist_to_tss"]),
                "tier": str(row["tier"]),
                "tss_position": int(row["tss"]),
                "strand": str(row["strand"]),
                "body_overlap": row["body_overlap"],
            }
        return links
    
    gene_links_series = (
        pair_df
        .groupby("peak_id", sort=False)
        .apply(build_links_for_peak)
        .rename("gene_links")
    )
    gene_links_series.fillna(empty_atac_gene_link_entry, inplace=True)
    return gene_links_series.reset_index()


# =============================================================================
# AGGREGATION HELPERS
# =============================================================================

def aggregate_genes_per_peak(
    pair_df: pd.DataFrame,
    tier_labels: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Create peak-level summary with genes per distance tier.
    """
    tier_labels = tier_labels or THRESHOLDS.tier_labels
    
    def agg_peak(df_sub):
        genes_by_tier = {tier: [] for tier in tier_labels}
        all_genes = []
        overlapping = []
        
        for _, row in df_sub.iterrows():
            gene = row["gene_name"]
            tier = str(row["tier"])
            all_genes.append(gene)
            
            if tier in genes_by_tier:
                genes_by_tier[tier].append(gene)
            
            if row["body_overlap"].get("overlaps", False):
                overlapping.append(gene)
        
        return pd.Series({
            "linked_genes": sorted(set(all_genes)),
            "genes_by_tier": {k: sorted(set(v)) for k, v in genes_by_tier.items()},
            "n_genes_total": len(set(all_genes)),
            "n_genes_overlapping": len(set(overlapping)),
        })
    
    return pair_df.groupby("peak_id").apply(agg_peak).reset_index()


def get_peaks_for_gene(
    pair_df: pd.DataFrame,
    gene_name: str,
    max_dist: Optional[int] = None,
    require_overlap: bool = False,
) -> pd.DataFrame:
    """
    Get all peaks linked to a specific gene.
    """
    df = pair_df[pair_df["gene_name"] == gene_name].copy()
    
    if max_dist is not None:
        df = df[df["dist_to_tss"] <= max_dist]
    
    if require_overlap:
        df = df[df["body_overlap"].apply(lambda x: x.get("overlaps", False))]
    
    return df.sort_values("dist_to_tss").reset_index(drop=True)
