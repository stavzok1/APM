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

def _gene_hit_from_row(row: pd.Series, payload: Optional[dict]) -> Dict[str, Any]:
    """Extract gene hit record from a feature row."""
    payload = payload or {}
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


def _lncrna_hit_from_row(row: pd.Series, payload: Optional[dict]) -> Dict[str, Any]:
    """Extract lncRNA hit record from a feature row."""
    payload = payload or {}
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


def _ccre_hit_from_row(row: pd.Series, payload: Optional[dict]) -> Dict[str, Any]:
    """Extract cCRE hit record from a feature row."""
    payload = payload or {}
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




def mirror_hits_into_boundaries_by_overlap(
    boundaries: pd.DataFrame,
    features_df: pd.DataFrame,
    *,
    biosample: str,
    feature_kind: str,
    hit_col: Optional[str] = None,
    max_hits_per_boundary: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror feature hits into boundaries based only on direct interval overlap.

    A feature is mirrored into a boundary only if:
        feature.start <= boundary.end and feature.end >= boundary.start

    Args:
        boundaries: Boundary DataFrame with columns:
                    chrom, start, end, boundary_id
        features_df: Feature DataFrame with columns:
                     chrom, start, end, ...
        biosample: Cell line / sample name to process
        feature_kind: "gene", "lncrna", or "ccre"
        hit_col: Output column name in boundaries (auto-derived if None)
        max_hits_per_boundary: Limit hits per boundary per biosample

    Returns:
        boundaries DataFrame with hit_col containing:
            boundaries.loc[row, hit_col][biosample] = [hit, hit, ...]
    """
    if feature_kind not in {"gene", "lncrna", "ccre"}:
        raise ValueError("feature_kind must be 'gene', 'lncrna', or 'ccre'")

    if hit_col is None:
        hit_col = {
            "gene": "gene_hits",
            "lncrna": "lncRNA_hits",
            "ccre": "cCRE_hits",
        }[feature_kind]

    hit_extractor = {
        "gene": _gene_hit_from_row,
        "lncrna": _lncrna_hit_from_row,
        "ccre": _ccre_hit_from_row,
    }[feature_kind]

    bd = boundaries.copy()
    bd = _ensure_dict_col(bd, hit_col)

    # Some sources use `name` instead of `boundary_id`.
    if "boundary_id" not in bd.columns and "name" in bd.columns:
        bd["boundary_id"] = bd["name"]

    required_boundary_cols = {"chrom", "start", "end", "boundary_id"}
    missing_boundary = required_boundary_cols - set(bd.columns)
    if missing_boundary:
        raise ValueError(
            f"boundaries missing required columns: {sorted(missing_boundary)}"
        )

    required_feature_cols = {"chrom", "start", "end"}
    missing_feature = required_feature_cols - set(features_df.columns)
    if missing_feature:
        raise ValueError(
            f"features_df missing required columns: {sorted(missing_feature)}"
        )

    # Normalize types for overlap checks and indexing.
    bd["boundary_id"] = bd["boundary_id"].astype(str)
    bd["chrom"] = bd["chrom"].astype(str)
    bd["start"] = pd.to_numeric(bd["start"], errors="coerce").astype("Int64")
    bd["end"] = pd.to_numeric(bd["end"], errors="coerce").astype("Int64")

    # Pre-index boundaries by chromosome
    boundaries_by_chrom = {chrom: group for chrom, group in bd.groupby("chrom", sort=False)}

    # boundary_id -> dataframe index
    bnd_to_idx = pd.Series(bd.index.values, index=bd["boundary_id"]).to_dict()

    for _, row in features_df.iterrows():
        chrom = row.get("chrom")
        f_start = row.get("start")
        f_end = row.get("end")

        if pd.isna(chrom) or pd.isna(f_start) or pd.isna(f_end):
            continue

        f_start = int(f_start)
        f_end = int(f_end)

        if f_start > f_end:
            f_start, f_end = f_end, f_start

        b_chr = boundaries_by_chrom.get(str(chrom))
        if b_chr is None or b_chr.empty:
            continue

        # direct overlap only
        overlapping = b_chr[
            (b_chr["start"] <= f_end) &
            (b_chr["end"] >= f_start)
        ]

        if overlapping.empty:
            continue

        # Boundary overlap mirroring does not have per-row TAD payloads; pass None safely.
        hit = hit_extractor(row, None)

        for _, b in overlapping.iterrows():
            bid = b["boundary_id"]
            idx = bnd_to_idx.get(bid)
            if idx is None:
                continue

            cur = bd.at[idx, hit_col]
            bucket = cur.get(biosample)
            if bucket is None:
                bucket = []
                cur[biosample] = bucket

            if max_hits_per_boundary is None or len(bucket) < max_hits_per_boundary:
                bucket.append(hit)

            bd.at[idx, hit_col] = cur

    return bd


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





def mirror_ccres_into_boundaries(
    boundaries: pd.DataFrame,
    ccre_df: pd.DataFrame,
    biosample: str,
    mode: str = "primary",
    max_hits_per_boundary: Optional[int] = None,
) -> pd.DataFrame:
    """
    Mirror cCRE annotations into boundaries table based on direct genomic overlap.

    Args:
        boundaries: Boundary DataFrame with columns:
                    chrom, start, end, boundary_id
        ccre_df: cCRE DataFrame with columns:
                 chrom, start, end, ...
        biosample: Cell line name
        mode: Kept for API compatibility. Ignored for overlap-based mirroring.
        max_hits_per_boundary: Optional limit

    Returns:
        boundaries with cCRE_hits column populated
    """
    return mirror_hits_into_boundaries_by_overlap(
        boundaries,
        ccre_df,
        biosample=biosample,
        feature_kind="ccre",
        max_hits_per_boundary=max_hits_per_boundary,
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
