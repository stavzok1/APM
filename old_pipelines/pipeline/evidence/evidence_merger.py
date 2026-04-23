"""
Evidence merger - combines all evidence sources into unified gene_links.

Orchestrates the merging of:
- SCREEN experimental links
- SCREEN computational links  
- ABC model predictions
- HiChIP loop data
"""

from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd

from ..schemas import (
    ensure_screen_block,
    empty_screen_block,
)


# =============================================================================
# MERGE LINK SOURCES
# =============================================================================

def merge_screen_links(
    links_exp: pd.DataFrame,
    links_comp: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge experimental and computational SCREEN links.
    
    Both should have been collapsed to nested structure with
    'screen_exp' or 'screen_comp' columns.
    """
    print("Merging SCREEN links...")
    
    # Ensure string types for merge keys
    merge_cols = ["ENCODE_id", "gene_id", "gene_name", "gene_type"]
    for df in (links_exp, links_comp):
        for col in merge_cols:
            if col in df.columns:
                df[col] = df[col].astype(str)
    
    # Outer merge to keep all links
    merged = links_exp.merge(
        links_comp,
        on=["ENCODE_id", "gene_id", "gene_name", "gene_type"],
        how="outer",
        suffixes=("_exp", "_comp"),
    )
    
    # Handle region column if present
    if "region" in merged.columns:
        pass  # Keep as is
    elif "region_comp" in merged.columns:
        merged["region"] = merged["region_comp"]
        merged.drop(columns=["region_comp"], inplace=True, errors="ignore")
    
    print(f"  Merged {len(merged)} total links")
    return merged


def merge_abc_links(
    screen_merged: pd.DataFrame,
    abc_collapsed: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge ABC predictions with SCREEN links.
    
    abc_collapsed should have columns: ENCODE_id, gene_name, ABC_enhancers
    """
    print("Merging ABC links...")
    
    if abc_collapsed.empty:
        screen_merged["ABC_enhancers"] = [[] for _ in range(len(screen_merged))]
        return screen_merged
    
    merged = screen_merged.merge(
        abc_collapsed[["ENCODE_id", "gene_name", "ABC_enhancers"]],
        on=["ENCODE_id", "gene_name"],
        how="left",
    )
    
    # Fill missing ABC with empty list
    merged["ABC_enhancers"] = merged["ABC_enhancers"].apply(
        lambda x: x if isinstance(x, list) else []
    )
    
    n_with_abc = (merged["ABC_enhancers"].apply(len) > 0).sum()
    print(f"  {n_with_abc}/{len(merged)} links have ABC predictions")
    
    return merged


# =============================================================================
# BUILD UNIFIED GENE LINKS
# =============================================================================

def build_unified_gene_links(
    merged_links: pd.DataFrame,
    screen_exp_biosamples: List[str],
    screen_exp_assays: List[str],
    screen_comp_biosamples: List[str],
    screen_comp_assays: List[str],
) -> pd.DataFrame:
    """
    Build final element-level table with nested gene_links column.
    
    Structure per element:
    {
        gene_name: {
            "gene_id": ...,
            "gene_type": ...,
            "screen_exp": {per_biosample, conservation_global, conservation_breast},
            "screen_comp": {per_biosample, conservation_global, conservation_breast},
            "ABC_enhancers": [...],
            "hichip": {},
        },
        ...
    }
    """
    print("Building unified gene_links...")
    
    def build_gene_links_for_element(df_element: pd.DataFrame) -> Dict[str, Dict]:
        out = {}
        
        for gene_name, df_gene in df_element.groupby("gene_name", sort=False):
            base = df_gene.iloc[0].to_dict()
            
            # Remove outer keys
            base.pop("ENCODE_id", None)
            base.pop("gene_name", None)
            
            # Ensure screen blocks are complete
            base["screen_exp"] = ensure_screen_block(
                base.get("screen_exp"),
                screen_exp_biosamples,
                screen_exp_assays,
            )
            base["screen_comp"] = ensure_screen_block(
                base.get("screen_comp"),
                screen_comp_biosamples,
                screen_comp_assays,
            )
            
            # Ensure ABC_enhancers exists
            abc_lists = [x for x in df_gene.get("ABC_enhancers", []) if isinstance(x, list)]
            base["ABC_enhancers"] = [e for lst in abc_lists for e in lst] if abc_lists else []
            
            # Initialize empty hichip (will be populated later)
            base["hichip"] = {}
            
            out[gene_name] = base
        
        return out
    
    gene_links_series = (
        merged_links
        .groupby("ENCODE_id", sort=False)
        .apply(build_gene_links_for_element)
        .rename("gene_links")
    )
    
    result = gene_links_series.reset_index()
    print(f"  Built {len(result)} element records with gene_links")
    
    return result


# =============================================================================
# FULL MERGE PIPELINE
# =============================================================================

def merge_all_evidence(
    links_exp: pd.DataFrame,
    links_comp: pd.DataFrame,
    abc_collapsed: pd.DataFrame,
    screen_exp_biosamples: List[str],
    screen_exp_assays: List[str],
    screen_comp_biosamples: List[str],
    screen_comp_assays: List[str],
) -> pd.DataFrame:
    """
    Complete evidence merging pipeline.
    
    Returns element-level table with gene_links column containing
    all evidence types.
    """
    # Step 1: Merge SCREEN sources
    merged = merge_screen_links(links_exp, links_comp)
    
    # Step 2: Merge ABC
    merged = merge_abc_links(merged, abc_collapsed)
    
    # Step 3: Build unified gene_links structure
    elements = build_unified_gene_links(
        merged,
        screen_exp_biosamples,
        screen_exp_assays,
        screen_comp_biosamples,
        screen_comp_assays,
    )
    elements.to_csv(r"C:\Users\stavz\Desktop\masters\APM\test1\elements_merged.csv", index=False)
    return elements


# =============================================================================
# ATTACH TO ELEMENT TABLE
# =============================================================================

def attach_gene_links_to_elements(
    elem_focus: pd.DataFrame,
    gene_links_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Attach gene_links column to existing element focus table.
    
    Merges on ENCODE_id.
    """
    print("Attaching gene_links to element table...")
    
    result = elem_focus.merge(
        gene_links_df[["ENCODE_id", "gene_links"]],
        on="ENCODE_id",
        how="left",
    )
    
    # Fill missing gene_links with empty dict
    result["gene_links"] = result["gene_links"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )
    
    n_with_links = (result["gene_links"].apply(len) > 0).sum()
    print(f"  {n_with_links}/{len(result)} elements have gene_links")
    
    return result


# =============================================================================
# VALIDATION
# =============================================================================

def validate_gene_links(elem_focus: pd.DataFrame) -> Dict[str, Any]:
    """
    Validate gene_links structure and return summary stats.
    """
    stats = {
        "total_elements": len(elem_focus),
        "elements_with_links": 0,
        "total_links": 0,
        "links_with_screen_exp": 0,
        "links_with_screen_comp": 0,
        "links_with_abc": 0,
        "links_with_hichip": 0,
    }
    
    if "gene_links" not in elem_focus.columns:
        return stats
    
    for _, row in elem_focus.iterrows():
        gl = row["gene_links"]
        if not isinstance(gl, dict) or len(gl) == 0:
            continue
        
        stats["elements_with_links"] += 1
        
        for gene_name, gdata in gl.items():
            stats["total_links"] += 1
            
            if gdata.get("screen_exp", {}).get("per_biosample"):
                stats["links_with_screen_exp"] += 1
            
            if gdata.get("screen_comp", {}).get("per_biosample"):
                stats["links_with_screen_comp"] += 1
            
            if gdata.get("ABC_enhancers"):
                stats["links_with_abc"] += 1
            
            if gdata.get("hichip"):
                stats["links_with_hichip"] += 1
    
    return stats


def print_evidence_summary(elem_focus: pd.DataFrame) -> None:
    """Print summary of evidence coverage."""
    stats = validate_gene_links(elem_focus)
    
    print("\n" + "=" * 50)
    print("EVIDENCE SUMMARY")
    print("=" * 50)
    print(f"Total elements: {stats['total_elements']}")
    print(f"Elements with gene links: {stats['elements_with_links']}")
    print(f"Total gene-element links: {stats['total_links']}")
    print(f"  - With SCREEN experimental: {stats['links_with_screen_exp']}")
    print(f"  - With SCREEN computational: {stats['links_with_screen_comp']}")
    print(f"  - With ABC predictions: {stats['links_with_abc']}")
    print(f"  - With HiChIP loops: {stats['links_with_hichip']}")
    print("=" * 50)
