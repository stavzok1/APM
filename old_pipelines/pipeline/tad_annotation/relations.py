"""
TAD domain relation and boundary proximity calculations.

Pure geometry functions for determining how features relate to TAD domains
and boundaries. No DataFrame manipulation - these are building blocks.

Functions:
    get_element_tad_relations: Classify overlap type vs all overlapping domains
    pick_primary_domain_id: Select single "primary" domain for a feature
    boundary_overlap_and_dist: Compute overlap/distance to a boundary interval
"""

from __future__ import annotations

from typing import Dict, Any, Optional, Tuple

import pandas as pd


# =============================================================================
# DOMAIN RELATIONS
# =============================================================================

def get_element_tad_relations(
    start: int,
    end: int,
    tad_domains: pd.DataFrame,
) -> Dict[str, str]:
    """
    For a feature interval [start, end], return relation label vs every overlapping domain.
    
    Args:
        start: Feature start position
        end: Feature end position  
        tad_domains: DataFrame with columns [domain_id, start, end] for ONE chromosome
    
    Returns:
        Dict mapping domain_id → relation label:
            - "contains": domain fully contains the feature
            - "contained": feature fully contains the domain
            - "overlap_left": feature extends past domain's left boundary
            - "overlap_right": feature extends past domain's right boundary
            - "overlap": generic overlap (shouldn't happen with above logic)
    
    Example:
        >>> rels = get_element_tad_relations(1000, 2000, tad_chrom)
        >>> # {'domain_001': 'contains', 'domain_002': 'overlap_right'}
    """
    rel = {}
    
    # Filter to overlapping domains only
    ov = tad_domains[(tad_domains["start"] <= end) & (tad_domains["end"] >= start)]
    
    for _, d in ov.iterrows():
        a, b = int(d["start"]), int(d["end"])
        did = d["domain_id"]
        
        if a <= start and end <= b:
            rel[did] = "contains"
        elif start <= a and b <= end:
            rel[did] = "contained"
        elif a <= start <= b < end:
            rel[did] = "overlap_left"
        elif start < a <= end <= b:
            rel[did] = "overlap_right"
        else:
            rel[did] = "overlap"
    
    return rel


def pick_primary_domain_id(
    start: int,
    end: int,
    tad_domains_chrom: pd.DataFrame,
) -> Optional[str]:
    """
    Choose a single "primary" domain for the feature.
    
    Selection logic:
        1. Prefer domain that contains the feature midpoint (stable, intuitive)
        2. Fallback: nearest domain by distance to the feature interval (handles gaps)
    
    Args:
        start: Feature start position
        end: Feature end position
        tad_domains_chrom: DataFrame with [domain_id, start, end] for ONE chromosome
    
    Returns:
        domain_id of primary domain, or None if no domains on chromosome
    
    Example:
        >>> primary = pick_primary_domain_id(1500, 2500, tad_chrom)
        >>> # 'domain_003'
    """
    if tad_domains_chrom.empty:
        return None
    
    mid = (start + end) // 2
    
    # First try: domain containing midpoint
    hit = tad_domains_chrom[
        (tad_domains_chrom["start"] <= mid) & (tad_domains_chrom["end"] >= mid)
    ]
    if not hit.empty:
        return hit.iloc[0]["domain_id"]
    
    # Fallback: nearest domain by distance to [start, end]
    # distance = 0 if overlaps; else gap to closest edge
    a = tad_domains_chrom["start"].astype(int).to_numpy()
    b = tad_domains_chrom["end"].astype(int).to_numpy()
    
    # dist to interval [start, end]:
    # - if domain ends before start => start - domain_end
    # - if domain starts after end => domain_start - end
    # - else overlap => 0
    dist = (start - b).clip(min=0) + (a - end).clip(min=0)
    i = int(dist.argmin())
    
    return tad_domains_chrom.iloc[i]["domain_id"]


# =============================================================================
# BOUNDARY PROXIMITY
# =============================================================================

def boundary_overlap_and_dist(
    feat_start: int,
    feat_end: int,
    b_start: int,
    b_end: int,
    b_pos: int,
) -> Dict[str, Any]:
    """
    Compute overlap and distance between a feature and a boundary.
    
    Args:
        feat_start: Feature interval start
        feat_end: Feature interval end
        b_start: Boundary interval start
        b_end: Boundary interval end
        b_pos: Boundary anchor position (typically midpoint or peak)
    
    Returns:
        Dict with:
            - overlap: bool, whether feature intersects boundary interval
            - dist_bp: int, distance from boundary anchor to feature (0 if inside)
    
    Example:
        >>> boundary_overlap_and_dist(1000, 2000, 1500, 1600, 1550)
        >>> # {'overlap': True, 'dist_bp': 0}
    """
    # Do intervals overlap?
    overlap = (b_start <= feat_end) and (b_end >= feat_start)
    
    # Distance to anchor from feature interval (0 if anchor inside feature)
    if b_pos < feat_start:
        dist_bp = feat_start - b_pos
    elif b_pos > feat_end:
        dist_bp = b_pos - feat_end
    else:
        dist_bp = 0
    
    return {"overlap": bool(overlap), "dist_bp": int(dist_bp)}


# =============================================================================
# NORMALIZED POSITION
# =============================================================================

def compute_normalized_position(
    feat_start: int,
    feat_end: int,
    domain_start: int,
    domain_end: int,
) -> Dict[str, float]:
    """
    Compute normalized position of feature within domain.
    
    Uses feature midpoint for stable positioning.
    
    Args:
        feat_start: Feature start
        feat_end: Feature end
        domain_start: Domain start
        domain_end: Domain end
    
    Returns:
        Dict with:
            - frac_from_left: [0,1] position from left boundary
            - frac_from_right: [0,1] position from right boundary  
            - frac_to_nearest_boundary: [0,0.5] distance to nearest boundary
    
    Note:
        Values outside [0,1] indicate feature extends beyond domain.
    """
    mid = (feat_start + feat_end) // 2
    d_len = max(1, domain_end - domain_start)
    
    frac_from_left = (mid - domain_start) / d_len
    frac_from_right = (domain_end - mid) / d_len
    frac_to_nearest = min(frac_from_left, frac_from_right)
    
    return {
        "frac_from_left": float(frac_from_left),
        "frac_from_right": float(frac_from_right),
        "frac_to_nearest_boundary": float(frac_to_nearest),
    }
