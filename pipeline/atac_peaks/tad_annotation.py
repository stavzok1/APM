"""
TAD annotation for ATAC peaks.

Provides:
- TAD domain annotation (using existing pipeline annotator)
- TAD boundary overlap detection (new functionality)
"""

from pathlib import Path
from typing import Dict, Any, List, Optional, Union

import pandas as pd
import numpy as np

from ..schemas import (
    empty_tad_boundary_overlaps_biosample,
)
from ..utils import normalize_chrom, compute_interval_overlap


# =============================================================================
# TAD DOMAIN ANNOTATION
# =============================================================================

def annotate_peaks_with_tads(
    peaks: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
    out_col: str = "TAD_domains",
) -> pd.DataFrame:
    """
    Annotate ATAC peaks with TAD domain context.
    
    This is a wrapper that uses the main pipeline's annotator.
    
    Args:
        peaks: ATAC peaks DataFrame
        tad_domains: TAD domains DataFrame
        domain_flanks: Domain flanks DataFrame
        boundaries: Boundaries DataFrame
        biosample: Biosample name
        out_col: Output column name
    
    Returns:
        Peaks DataFrame with TAD_domains column added
    """
    from ..tad_annotation.annotator import annotate_df_with_tads
    
    return annotate_df_with_tads(
        peaks,
        kind="interval",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
        out_col=out_col,
    )


# =============================================================================
# TAD BOUNDARY OVERLAP DETECTION
# =============================================================================

def annotate_peaks_with_boundary_overlaps(
    peaks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
    out_col: str = "TAD_boundary_overlaps",
    boundary_strength_col: Optional[str] = "strength",
) -> pd.DataFrame:
    """
    Annotate ATAC peaks with TAD boundary overlaps.
    
    Adds df[out_col][biosample] = {
        "overlaps_boundary": bool,
        "n_boundaries": int,
        "boundaries": [
            {
                "boundary_id": str,
                "overlap_bp": int,
                "overlap_interval": [start, end],
                "boundary_strength": str or None,
            },
            ...
        ]
    }
    
    Args:
        peaks: ATAC peaks DataFrame
        boundaries: TAD boundaries DataFrame with [boundary_id, chrom, start, end]
        biosample: Biosample name (key in nested dict)
        out_col: Output column name
        boundary_strength_col: Column name for boundary strength (optional)
    
    Returns:
        Peaks DataFrame with boundary overlap annotations
    """
    df = peaks.copy()
    
    # Ensure output column exists
    if out_col not in df.columns:
        df[out_col] = [{} for _ in range(len(df))]
    else:
        df[out_col] = df[out_col].apply(
            lambda x: {} if x is None or (isinstance(x, float) and pd.isna(x)) else x
        )
    
    # Normalize chromosomes
    boundaries = boundaries.copy()
    boundaries["chrom"] = normalize_chrom(boundaries["chrom"])
    df["chrom"] = normalize_chrom(df["chrom"])
    
    # Determine boundary ID column
    bid_col = "boundary_id" if "boundary_id" in boundaries.columns else (
        "name" if "name" in boundaries.columns else None
    )
    if bid_col is None:
        raise ValueError("boundaries must contain 'boundary_id' or 'name' column")
    
    # Process per chromosome for efficiency
    for chrom in df["chrom"].dropna().unique():
        peak_mask = df["chrom"] == chrom
        if not peak_mask.any():
            continue
        
        bounds_chrom = boundaries[boundaries["chrom"] == chrom]
        if bounds_chrom.empty:
            # No boundaries on this chromosome
            for idx in df.index[peak_mask]:
                current = df.at[idx, out_col]
                current[biosample] = empty_tad_boundary_overlaps_biosample()
                df.at[idx, out_col] = current
            continue
        
        # Build boundary lookup
        bounds_list = bounds_chrom.to_dict("records")
        
        def find_boundary_overlaps(row) -> Dict[str, Any]:
            peak_start = int(row["start"])
            peak_end = int(row["end"])
            
            overlapping = []
            for b in bounds_list:
                b_start = int(b["start"])
                b_end = int(b["end"])
                
                ov = compute_interval_overlap(peak_start, peak_end, b_start, b_end)
                
                if ov["overlaps"]:
                    entry = {
                        "boundary_id": b[bid_col],
                        "overlap_bp": ov["overlap_bp"],
                        "overlap_interval": ov["overlap_interval"],
                        "boundary_strength": b.get(boundary_strength_col) if boundary_strength_col else None,
                    }
                    overlapping.append(entry)
            
            return {
                "overlaps_boundary": len(overlapping) > 0,
                "n_boundaries": len(overlapping),
                "boundaries": overlapping,
            }
        
        # Apply to masked rows
        payloads = df.loc[peak_mask].apply(find_boundary_overlaps, axis=1)
        
        # Merge into existing dict
        for idx, payload in zip(df.index[peak_mask], payloads):
            current = df.at[idx, out_col]
            current[biosample] = payload
            df.at[idx, out_col] = current
    
    return df


# =============================================================================
# MULTI-BIOSAMPLE ANNOTATION
# =============================================================================

def annotate_peaks_with_all_tad_sources(
    peaks: pd.DataFrame,
    processed_dir: Union[str, Path],
    biosamples: Optional[List[str]] = None,
    verbose: bool = True,
    emit_summary: bool = True,
) -> pd.DataFrame:
    """
    Annotate peaks with TAD domains and boundary overlaps from multiple sources.
    
    Args:
        peaks: ATAC peaks DataFrame
        processed_dir: Path to processed TAD data directory
        biosamples: List of biosamples to use (None = all available)
        verbose: Print per-biosample progress
        emit_summary: Print final coverage counts (set False when a caller batches rows)
    
    Returns:
        Annotated peaks DataFrame
    """
    processed_dir = Path(processed_dir)
    
    # Discover available biosamples
    if biosamples is None:
        biosamples = [d.name for d in processed_dir.iterdir() 
                      if d.is_dir() and (d / "domains.tsv").exists()]
    
    if verbose:
        print(f"Annotating peaks with {len(biosamples)} TAD sources...")
    
    df = peaks.copy()
    
    for biosample in biosamples:
        biosample_dir = processed_dir / biosample
        
        # Load TAD data files
        domains_path = biosample_dir / "domains.tsv"
        flanks_path = biosample_dir / "flanks.tsv"
        bounds_path = biosample_dir / "boundaries.tsv"
        
        if not all(p.exists() for p in [domains_path, flanks_path, bounds_path]):
            if verbose:
                print(f"  Skipping {biosample}: missing files")
            continue
        
        domains = pd.read_csv(domains_path, sep="\t")
        flanks = pd.read_csv(flanks_path, sep="\t")
        bounds = pd.read_csv(bounds_path, sep="\t")
        
        if verbose:
            print(f"  Processing {biosample}...")
        
        # TAD domain annotation
        df = annotate_peaks_with_tads(
            df,
            tad_domains=domains,
            domain_flanks=flanks,
            boundaries=bounds,
            biosample=biosample,
        )
        
        # TAD boundary overlap annotation
        df = annotate_peaks_with_boundary_overlaps(
            df,
            boundaries=bounds,
            biosample=biosample,
        )
    
    if emit_summary and verbose:
        # Summary
        n_with_tad = df["TAD_domains"].apply(
            lambda x: len(x) > 0 if isinstance(x, dict) else False
        ).sum()
        n_boundary_overlap = df["TAD_boundary_overlaps"].apply(
            lambda x: any(v.get("overlaps_boundary", False) for v in x.values()) 
            if isinstance(x, dict) else False
        ).sum()
        print(f"  Peaks with TAD annotations: {n_with_tad}/{len(df)}")
        print(f"  Peaks overlapping boundaries: {n_boundary_overlap}/{len(df)}")
    
    return df


# =============================================================================
# QUERY HELPERS
# =============================================================================

def get_peaks_at_boundaries(
    peaks: pd.DataFrame,
    biosample: Optional[str] = None,
    boundary_col: str = "TAD_boundary_overlaps",
) -> pd.DataFrame:
    """
    Get peaks that overlap TAD boundaries.
    
    Args:
        peaks: Annotated peaks DataFrame
        biosample: Specific biosample to query (None = any biosample)
        boundary_col: Column containing boundary annotations
    
    Returns:
        Filtered DataFrame of peaks at boundaries
    """
    if boundary_col not in peaks.columns:
        raise ValueError(f"Column '{boundary_col}' not found")
    
    def check_boundary_overlap(x):
        if not isinstance(x, dict):
            return False
        
        if biosample:
            return x.get(biosample, {}).get("overlaps_boundary", False)
        else:
            return any(v.get("overlaps_boundary", False) for v in x.values())
    
    mask = peaks[boundary_col].apply(check_boundary_overlap)
    return peaks[mask].copy()


def summarize_boundary_overlaps(
    peaks: pd.DataFrame,
    boundary_col: str = "TAD_boundary_overlaps",
) -> pd.DataFrame:
    """
    Summarize boundary overlap statistics per biosample.
    """
    if boundary_col not in peaks.columns:
        raise ValueError(f"Column '{boundary_col}' not found")
    
    # Collect all biosamples
    all_biosamples = set()
    for x in peaks[boundary_col]:
        if isinstance(x, dict):
            all_biosamples.update(x.keys())
    
    results = []
    for biosample in sorted(all_biosamples):
        n_at_boundary = 0
        total_overlaps = 0
        overlap_bps = []
        
        for x in peaks[boundary_col]:
            if not isinstance(x, dict):
                continue
            bs_data = x.get(biosample, {})
            if bs_data.get("overlaps_boundary", False):
                n_at_boundary += 1
                bounds = bs_data.get("boundaries", [])
                total_overlaps += len(bounds)
                for b in bounds:
                    overlap_bps.append(b.get("overlap_bp", 0))
        
        results.append({
            "biosample": biosample,
            "n_peaks_at_boundary": n_at_boundary,
            "n_total_boundary_overlaps": total_overlaps,
            "mean_overlap_bp": np.mean(overlap_bps) if overlap_bps else 0,
        })
    
    return pd.DataFrame(results)
