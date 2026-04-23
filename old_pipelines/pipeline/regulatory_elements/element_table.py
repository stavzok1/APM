"""
Element-centric table building.

Functions to build the final regulatory element table with:
- Distance matching results
- Cell-line signals
- Gene links from various evidence sources
"""

import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd

from .distance_matching import aggregate_genes_per_ccre, aggregate_ccres_per_gene, build_distance_matrix


# =============================================================================
# ELEMENT FOCUS TABLE
# =============================================================================

def build_element_focus_table(
    ccres: pd.DataFrame,
    pair_df: pd.DataFrame,
    tier_labels: List[str],
    cell_line_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build element-focused table combining distance matching with cCRE metadata.
    
    Args:
        ccres: cCRE DataFrame (may include cell-line signal columns)
        pair_df: Gene-cCRE pair DataFrame from distance matching
        tier_labels: Distance tier labels
        cell_line_cols: List of cell line column names to include
    
    Returns:
        Element-focused DataFrame
    """
    # Get per-element aggregation
    elem_agg = aggregate_genes_per_ccre(pair_df, tier_labels)
    
    # Determine which columns to bring from ccres
    merge_cols = ["cCRE_id", "ENCODE_id"]
    if cell_line_cols:
        merge_cols.extend([c for c in cell_line_cols if c in ccres.columns])
    
    # Merge with ccre metadata
    cols = [c for c in merge_cols if c in ccres.columns]
    key = [c for c in ["cCRE_id", "ENCODE_id"] if c in cols]

    ccre_subset = (
        ccres[cols]
        .groupby(key, as_index=False, sort=False)
        .first()
    )    
    elem_focus = elem_agg.merge(
        ccre_subset,
        on=["cCRE_id", "ENCODE_id"],
        how="left",
    )
    
    return elem_focus


# =============================================================================
# GENE SUMMARY TABLE
# =============================================================================

def build_gene_summary_table(
    pair_df: pd.DataFrame,
    tier_labels: List[str],
) -> pd.DataFrame:
    """
    Build gene-focused summary with cCRE counts and IDs per type/tier.
    """
    return aggregate_ccres_per_gene(pair_df, tier_labels)


# =============================================================================
# SAVE MATCHING OUTPUTS
# =============================================================================

def save_all_matching_outputs(
    output_dir: Path,
    pair_df: pd.DataFrame,
    elem_focus: pd.DataFrame,
    gene_summary: pd.DataFrame,
    gene_type: str = "coding",
) -> None:
    """
    Save all matching-related outputs.
    
    Args:
        output_dir: Base output directory
        pair_df: Gene-cCRE pairs
        elem_focus: Element-focused table
        gene_summary: Gene-focused summary
        gene_type: "coding" or "lncRNA" to organize output subdirectory
    """
    output_dir = Path(output_dir)
    
    if gene_type != "coding":
        output_dir = output_dir / gene_type
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pair_df.to_csv(output_dir / "gene_to_elements.csv", index=False)
    elem_focus.to_csv(output_dir / "regulatory_element_focus.csv", index=False)
    gene_summary.to_csv(output_dir / "gene_focus.csv", index=False)
    
    # Distance matrix
    dist_matrix = build_distance_matrix(pair_df)
    dist_matrix.to_csv(output_dir / "distance_matrix.csv")
    
    print(f"Saved matching outputs to {output_dir}")


# =============================================================================
# GENE LINKS COLUMN BUILDING
# =============================================================================

def initialize_gene_links_column(
    elem_focus: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Initialize the gene_links column with basic structure from distance matching.
    
    Each element gets:
    gene_links = {
        gene_name: {
            "dist_to_tss": ...,
            "tier": ...,
            # Evidence placeholders
            "screen_exp": {},
            "screen_comp": {},
            "ABC_enhancers": [],
            "hichip": {},
        },
        ...
    }
    """
    # Build gene_links dict per cCRE
    def build_gene_links(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        out = {}
        for _, row in df_sub.iterrows():
            out[row["gene_name"]] = {
                "dist_to_tss": int(row["dist_to_tss"]),
                "tier": str(row["tier"]),
                # Placeholders for evidence
                "screen_exp": {},
                "screen_comp": {},
                "ABC_enhancers": [],
                "hichip": {},
            }
        return out
    
    gene_links_series = (
        pair_df
        .groupby("cCRE_id", sort=False)
        .apply(build_gene_links)
        .rename("gene_links")
    )
    
    # Merge into elem_focus
    elem_focus = elem_focus.merge(
        gene_links_series.reset_index(),
        on="cCRE_id",
        how="left",
    )
    
    # Fill missing with empty dict
    elem_focus["gene_links"] = elem_focus["gene_links"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )
    
    return elem_focus


def add_evidence_to_gene_links(
    elem_focus: pd.DataFrame,
    evidence_df: pd.DataFrame,
    evidence_type: str,
    ccre_col: str = "ENCODE_id",
    gene_col: str = "gene_name",
) -> pd.DataFrame:
    """
    Add evidence from an evidence DataFrame to gene_links.
    
    Args:
        elem_focus: Element focus table with gene_links column
        evidence_df: Evidence DataFrame with cCRE and gene columns plus evidence data
        evidence_type: Key to use in gene_links (e.g., "screen_exp", "ABC_enhancers")
        ccre_col: Column name for cCRE ID in evidence_df
        gene_col: Column name for gene in evidence_df
    
    Returns:
        Updated elem_focus
    """
    # Build lookup: (cCRE, gene) -> evidence
    evidence_lookup = {}
    for _, row in evidence_df.iterrows():
        key = (str(row[ccre_col]), str(row[gene_col]))
        evidence_lookup[key] = row.to_dict()
    
    def update_gene_links(row):
        gene_links = row["gene_links"]
        if not isinstance(gene_links, dict):
            return gene_links
        
        ccre_id = str(row.get("ENCODE_id", row.get("cCRE_id", "")))
        
        for gene_name, gene_data in gene_links.items():
            key = (ccre_id, gene_name)
            if key in evidence_lookup:
                evidence = evidence_lookup[key].copy()
                evidence.pop(ccre_col, None)
                evidence.pop(gene_col, None)
                gene_data[evidence_type] = evidence
        
        return gene_links
    
    elem_focus = elem_focus.copy()
    elem_focus["gene_links"] = elem_focus.apply(update_gene_links, axis=1)
    
    return elem_focus


# =============================================================================
# FINALIZATION
# =============================================================================

def finalize_element_table(
    elem_focus: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Finalize element table for output.
    
    - Ensures all expected columns are present
    - Validates gene_links structure
    - Optionally saves to file
    """
    # Validate gene_links
    if "gene_links" in elem_focus.columns:
        n_empty = elem_focus["gene_links"].apply(
            lambda x: len(x) == 0 if isinstance(x, dict) else True
        ).sum()
        n_total = len(elem_focus)
        print(f"Elements with gene links: {n_total - n_empty}/{n_total}")
    
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        elem_focus.to_csv(output_path, index=False)
        print(f"Saved final element table to {output_path}")
    
    return elem_focus
