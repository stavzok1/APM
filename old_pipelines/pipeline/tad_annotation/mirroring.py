"""
Mirror feature annotations into TAD domain tables.

After annotating features (genes, cCREs, lncRNAs) with their TAD context,
this module mirrors those annotations back into the TAD domains table,
creating columns like gene_hits, cCRE_hits, lncRNA_hits.

This enables querying from either direction:
    - "Which TAD contains gene X?" → annotator.py
    - "What genes are in TAD Y?" → mirroring.py

Main function:
    mirror_hits_into_domains: Add feature hits to tad_domains DataFrame
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional

import pandas as pd


# =============================================================================
# UTILITIES
# =============================================================================

def _ensure_dict_col(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Ensure column exists and contains dicts (not NaN/None)."""
    if col not in df.columns:
        df[col] = [{} for _ in range(len(df))]
    else:
        df[col] = df[col].apply(
            lambda x: {} if x is None or (isinstance(x, float) and pd.isna(x)) else x
        )
    return df


# =============================================================================
# HIT EXTRACTION
# =============================================================================

def _gene_hit_from_row(row: pd.Series, payload: dict) -> Dict[str, Any]:
    """Extract gene hit record from a feature row."""
    return {
        "gene_id": row.get("gene_id", None),
        "gene_name": row.get("gene_name", None),
        "gene_type": row.get("gene_type", row.get("type", None)),
        "chrom": row.get("chrom", None),
        "start": int(row.get("start")) if pd.notna(row.get("start")) else None,
        "end": int(row.get("end")) if pd.notna(row.get("end")) else None,
        "strand": row.get("strand", None),
        "rel": payload.get("primary", {}).get("rel") if payload.get("primary") else None,
    }


def _lncrna_hit_from_row(row: pd.Series, payload: dict) -> Dict[str, Any]:
    """Extract lncRNA hit record from a feature row."""
    return {
        "gene_id": row.get("gene_id", None),
        "gene_name": row.get("gene_name", None),
        "gene_type": row.get("gene_type", row.get("type", None)),
        "chrom": row.get("chrom", None),
        "start": int(row.get("start")) if pd.notna(row.get("start")) else None,
        "end": int(row.get("end")) if pd.notna(row.get("end")) else None,
        "strand": row.get("strand", None),
        "rel": payload.get("primary", {}).get("rel") if payload.get("primary") else None,
    }


def _ccre_hit_from_row(row: pd.Series, payload: dict) -> Dict[str, Any]:
    """Extract cCRE hit record from a feature row."""
    return {
        "cCRE_id": row.get("cCRE_id", row.get("ccre_id", None)),
        "chrom": row.get("chrom", None),
        "start": int(row.get("start")) if pd.notna(row.get("start")) else None,
        "end": int(row.get("end")) if pd.notna(row.get("end")) else None,
        "rel": payload.get("primary", {}).get("rel") if payload.get("primary") else None,
    }


# =============================================================================
# MAIN MIRRORING FUNCTION
# =============================================================================

def mirror_hits_into_domains(
    tad_domains: pd.DataFrame,
    features_df: pd.DataFrame,
    *,
    biosample: str,
    feature_kind: str,
    tad_col: str = "TAD_domains",
    hit_col: Optional[str] = None,
    mode: str = "primary",
    max_hits_per_domain: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror feature→domain annotations into tad_domains as a dict keyed by biosample.
    
    After annotating features with TAD context via annotate_df_with_tads(),
    this function reads those annotations and writes them back into the
    tad_domains table, enabling queries like "what genes are in this domain?"
    
    Args:
        tad_domains: TAD domains DataFrame with domain_id column
        features_df: Features DataFrame (genes, lncRNAs, cCREs) with TAD_domains column
        biosample: Cell line / sample name to process
        feature_kind: "gene", "lncrna", or "ccre"
        tad_col: Column name in features_df containing TAD annotations
        hit_col: Output column name in tad_domains (auto-derived if None)
        mode: "primary" (only primary domain) or "all" (all overlapping domains)
        max_hits_per_domain: Limit hits per domain per biosample (None = no limit)
    
    Returns:
        tad_domains DataFrame with hit_col containing:
            tad_domains.loc[domain_row, hit_col][biosample] = [hit, hit, ...]
    
    Example:
        >>> tad_domains = mirror_hits_into_domains(
        ...     tad_domains,
        ...     genes_df,
        ...     biosample="Kim_T47D",
        ...     feature_kind="gene",
        ...     mode="primary",
        ... )
        >>> # Now tad_domains["gene_hits"]["Kim_T47D"] contains gene records
    """
    if feature_kind not in {"gene", "lncrna", "ccre"}:
        raise ValueError("feature_kind must be 'gene', 'lncrna', or 'ccre'")
    
    # Auto-derive hit column name
    if hit_col is None:
        hit_col = {
            "gene": "gene_hits",
            "lncrna": "lncRNA_hits",
            "ccre": "cCRE_hits",
        }[feature_kind]
    
    # Select hit extractor
    hit_extractor = {
        "gene": _gene_hit_from_row,
        "lncrna": _lncrna_hit_from_row,
        "ccre": _ccre_hit_from_row,
    }[feature_kind]
    
    td = tad_domains.copy()
    td = _ensure_dict_col(td, hit_col)
    
    if "domain_id" not in td.columns:
        raise ValueError("tad_domains must include 'domain_id' column")
    
    # Build domain_id → row index mapping
    dom_to_idx = pd.Series(td.index.values, index=td["domain_id"]).to_dict()
    
    # Process each feature row
    for _, row in features_df.iterrows():
        ann = row.get(tad_col, None)
        if not isinstance(ann, dict):
            continue
        
        cl = ann.get(biosample, None)
        if cl is None:
            continue
        
        # Handle both single-interval and multi-interval payloads
        payloads = cl.get("intervals", None)
        if payloads is None:
            payloads = [cl]
        
        for p in payloads:
            domains_dict = p.get("domains", {}) or {}
            primary = p.get("primary", None)
            
            # Select which domains to annotate based on mode
            if mode == "primary":
                if not (primary and primary.get("domain_id") is not None):
                    continue
                domain_ids = [primary["domain_id"]]
            else:  # "all"
                domain_ids = list(domains_dict.keys())
                if not domain_ids:
                    continue
            
            hit = hit_extractor(row, p)
            
            for did in domain_ids:
                idx = dom_to_idx.get(did)
                if idx is None:
                    continue
                
                cur = td.at[idx, hit_col]
                bucket = cur.get(biosample)
                if bucket is None:
                    bucket = []
                    cur[biosample] = bucket
                
                # Apply max hits limit if specified
                if max_hits_per_domain is None or len(bucket) < max_hits_per_domain:
                    bucket.append(hit)
                
                td.at[idx, hit_col] = cur
    
    return td


# =============================================================================
# CONVENIENCE WRAPPERS
# =============================================================================

def mirror_genes_into_domains(
    tad_domains: pd.DataFrame,
    genes_df: pd.DataFrame,
    biosample: str,
    mode: str = "primary",
    max_hits_per_domain: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror gene annotations into TAD domains table.
    
    Args:
        tad_domains: TAD domains DataFrame
        genes_df: Genes DataFrame with TAD_domains column
        biosample: Cell line name
        mode: "primary" or "all"
        max_hits_per_domain: Optional limit
    
    Returns:
        tad_domains with gene_hits column populated
    """
    return mirror_hits_into_domains(
        tad_domains,
        genes_df,
        biosample=biosample,
        feature_kind="gene",
        mode=mode,
        max_hits_per_domain=max_hits_per_domain,
    )


def mirror_lncrnas_into_domains(
    tad_domains: pd.DataFrame,
    lncrnas_df: pd.DataFrame,
    biosample: str,
    mode: str = "primary",
    max_hits_per_domain: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror lncRNA annotations into TAD domains table.
    
    Args:
        tad_domains: TAD domains DataFrame
        lncrnas_df: lncRNAs DataFrame with TAD_domains column
        biosample: Cell line name
        mode: "primary" or "all"
        max_hits_per_domain: Optional limit
    
    Returns:
        tad_domains with lncRNA_hits column populated
    """
    return mirror_hits_into_domains(
        tad_domains,
        lncrnas_df,
        biosample=biosample,
        feature_kind="lncrna",
        mode=mode,
        max_hits_per_domain=max_hits_per_domain,
    )


def mirror_ccres_into_domains(
    tad_domains: pd.DataFrame,
    ccre_df: pd.DataFrame,
    biosample: str,
    mode: str = "primary",
    max_hits_per_domain: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror cCRE annotations into TAD domains table.
    
    Args:
        tad_domains: TAD domains DataFrame
        ccre_df: cCREs DataFrame with TAD_domains column
        biosample: Cell line name
        mode: "primary" or "all"
        max_hits_per_domain: Optional limit
    
    Returns:
        tad_domains with cCRE_hits column populated
    """
    return mirror_hits_into_domains(
        tad_domains,
        ccre_df,
        biosample=biosample,
        feature_kind="ccre",
        mode=mode,
        max_hits_per_domain=max_hits_per_domain,
    )


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def mirror_all_features_into_domains(
    tad_domains: pd.DataFrame,
    genes_df: Optional[pd.DataFrame] = None,
    lncrnas_df: Optional[pd.DataFrame] = None,
    ccre_df: Optional[pd.DataFrame] = None,
    biosample: str = None,
    mode: str = "primary",
) -> pd.DataFrame:
    """
    Mirror all feature types into TAD domains in one call.
    
    Args:
        tad_domains: TAD domains DataFrame
        genes_df: Optional genes DataFrame with TAD_domains
        lncrnas_df: Optional lncRNAs DataFrame with TAD_domains
        ccre_df: Optional cCREs DataFrame with TAD_domains
        biosample: Cell line name
        mode: "primary" or "all"
    
    Returns:
        tad_domains with all hit columns populated
    """
    if biosample is None:
        raise ValueError("biosample must be specified")
    
    result = tad_domains
    
    if genes_df is not None:
        result = mirror_genes_into_domains(result, genes_df, biosample, mode)
    
    if lncrnas_df is not None:
        result = mirror_lncrnas_into_domains(result, lncrnas_df, biosample, mode)
    
    if ccre_df is not None:
        result = mirror_ccres_into_domains(result, ccre_df, biosample, mode)
    
    return result


# =============================================================================
# QUERY HELPERS
# =============================================================================

def get_domain_gene_count(
    tad_domains: pd.DataFrame,
    biosample: str,
    hit_col: str = "gene_hits",
) -> pd.Series:
    """
    Get gene count per domain for a biosample.
    
    Returns:
        Series indexed by domain_id with gene counts
    """
    def count_hits(d):
        if not isinstance(d, dict):
            return 0
        bucket = d.get(biosample, [])
        return len(bucket) if bucket else 0
    
    return tad_domains.set_index("domain_id")[hit_col].apply(count_hits)


def get_domains_containing_gene(
    tad_domains: pd.DataFrame,
    gene_name: str,
    biosample: str,
    hit_col: str = "gene_hits",
) -> List[str]:
    """
    Find all domains containing a specific gene.
    
    Returns:
        List of domain_ids
    """
    result = []
    for _, row in tad_domains.iterrows():
        hits_dict = row.get(hit_col, {})
        if not isinstance(hits_dict, dict):
            continue
        bucket = hits_dict.get(biosample, [])
        if not bucket:
            continue
        for hit in bucket:
            if hit.get("gene_name") == gene_name:
                result.append(row["domain_id"])
                break
    return result


def get_genes_in_domain(
    tad_domains: pd.DataFrame,
    domain_id: str,
    biosample: str,
    hit_col: str = "gene_hits",
) -> List[Dict[str, Any]]:
    """
    Get all genes in a specific domain.
    
    Returns:
        List of gene hit dicts
    """
    row = tad_domains[tad_domains["domain_id"] == domain_id]
    if row.empty:
        return []
    
    hits_dict = row.iloc[0].get(hit_col, {})
    if not isinstance(hits_dict, dict):
        return []
    
    return hits_dict.get(biosample, [])
