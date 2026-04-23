"""
cCRE overlap matching for SNVs/indels.

Matches point mutations and small indels to overlapping cCREs
to identify regulatory element disruptions.
"""

from typing import List, Dict, Any, Optional

import pandas as pd


def match_snvs_to_ccres(
    snv_df: pd.DataFrame,
    elements: pd.DataFrame,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
) -> pd.DataFrame:
    """
    Match SNVs to overlapping cCREs.
    
    For each SNV, finds all cCREs that contain the variant position,
    storing results in a new 'cCRE_hits' column.
    
    Args:
        snv_df: DataFrame with SNV data (must have chrom, pos columns)
        elements: cCRE DataFrame with chrom, start, end, cCRE_id, raw_type, 
                  and optionally genes_by_exact_dist
        chrom_col: Name of chromosome column in snv_df
        pos_col: Name of position column in snv_df
    
    Returns:
        DataFrame with 'cCRE_hits' column containing list of dicts:
        [
            {
                "cCRE_id": str,
                "elem_type": str,
                "chrom": str,
                "elem_start": int,
                "elem_end": int,
                "genes_by_exact_dist": str or None,
            },
            ...
        ]
    """
    df = snv_df.copy()
    
    # Pre-index elements by chromosome for faster lookup
    elements_by_chrom = {
        chrom: group for chrom, group in elements.groupby("chrom")
    }
    
    ccre_hits_col = []
    
    for idx, snv in df.iterrows():
        chrom = snv.get(chrom_col)
        pos = snv.get(pos_col)
        
        # Handle missing data
        if pd.isna(chrom) or pd.isna(pos):
            ccre_hits_col.append([])
            continue
        
        pos = int(pos)
        
        # Get elements on same chromosome
        e_chr = elements_by_chrom.get(chrom)
        if e_chr is None or e_chr.empty:
            ccre_hits_col.append([])
            continue
        
        # Find overlapping elements (pos within [start, end])
        overlapping = e_chr[
            (e_chr["start"] <= pos) & (e_chr["end"] >= pos)
        ]
        
        if overlapping.empty:
            ccre_hits_col.append([])
            continue
        
        # Build hit records
        hits = []
        for _, elem in overlapping.iterrows():
            hit = {
                "cCRE_id": elem.get("cCRE_id"),
                "elem_type": elem.get("raw_type"),
                "chrom": chrom,
                "elem_start": int(elem["start"]),
                "elem_end": int(elem["end"]),
                "genes_by_exact_dist": elem.get("genes_by_exact_dist"),
            }
            hits.append(hit)
        
        ccre_hits_col.append(hits)
    
    df["cCRE_hits"] = ccre_hits_col
    return df


def summarize_ccre_hits(snv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary statistics for cCRE hits across all SNVs.
    
    Args:
        snv_df: DataFrame with 'cCRE_hits' column
    
    Returns:
        DataFrame with per-SNV summary:
        - n_ccre_hits: Number of cCREs overlapping
        - ccre_types: Comma-separated unique cCRE types
        - ccre_ids: Comma-separated cCRE IDs
        - linked_genes: Comma-separated genes linked to hit cCREs
    """
    def summarize_row(row):
        hits = row.get("cCRE_hits", [])
        if not hits:
            return pd.Series({
                "n_ccre_hits": 0,
                "ccre_types": "",
                "ccre_ids": "",
                "linked_genes": "",
            })
        
        types = sorted(set(h.get("elem_type", "") for h in hits if h.get("elem_type")))
        ids = sorted(set(h.get("cCRE_id", "") for h in hits if h.get("cCRE_id")))
        
        # Extract genes from genes_by_exact_dist strings
        genes = set()
        for h in hits:
            genes_str = h.get("genes_by_exact_dist", "")
            if genes_str and isinstance(genes_str, str):
                # Format is "gene1:dist1,gene2:dist2,..."
                for gene_dist in genes_str.split(","):
                    if ":" in gene_dist:
                        gene = gene_dist.split(":")[0].strip()
                        if gene:
                            genes.add(gene)
        
        return pd.Series({
            "n_ccre_hits": len(hits),
            "ccre_types": ",".join(types),
            "ccre_ids": ",".join(ids),
            "linked_genes": ",".join(sorted(genes)),
        })
    
    summary = snv_df.apply(summarize_row, axis=1)
    return pd.concat([snv_df, summary], axis=1)


def filter_snvs_by_ccre_type(
    snv_df: pd.DataFrame,
    ccre_types: Optional[List[str]] = None,
    require_hit: bool = True,
) -> pd.DataFrame:
    """
    Filter SNVs based on cCRE hit properties.
    
    Args:
        snv_df: DataFrame with 'cCRE_hits' column
        ccre_types: If provided, only keep SNVs hitting these cCRE types
                    (e.g., ["pELS", "dELS", "PLS"])
        require_hit: If True, filter out SNVs with no cCRE hits
    
    Returns:
        Filtered DataFrame
    """
    def check_row(row):
        hits = row.get("cCRE_hits", [])
        
        if require_hit and not hits:
            return False
        
        if ccre_types is not None:
            hit_types = {h.get("elem_type", "").split(",")[0] for h in hits}
            if not hit_types.intersection(set(ccre_types)):
                return False
        
        return True
    
    mask = snv_df.apply(check_row, axis=1)
    return snv_df[mask].copy()


def get_snvs_in_enhancers(snv_df: pd.DataFrame) -> pd.DataFrame:
    """Convenience filter for SNVs in enhancer elements (pELS, dELS)."""
    return filter_snvs_by_ccre_type(snv_df, ccre_types=["pELS", "dELS"])


def get_snvs_in_promoters(snv_df: pd.DataFrame) -> pd.DataFrame:
    """Convenience filter for SNVs in promoter-like elements (PLS)."""
    return filter_snvs_by_ccre_type(snv_df, ccre_types=["PLS"])


def get_snvs_in_ctcf(snv_df: pd.DataFrame) -> pd.DataFrame:
    """Convenience filter for SNVs in CTCF-bound elements."""
    return filter_snvs_by_ccre_type(snv_df, ccre_types=["CTCF-only", "CTCF-bound"])
