"""
HiChIP loop data processing.

Complete pipeline: load loops → normalize formats → compute overlaps → build nested dicts
"""

import os
import re
import glob
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    import pyranges as pr
    HAS_PYRANGES = True
except ImportError:
    HAS_PYRANGES = False
    print("Warning: pyranges not installed. HiChIP integration will be limited.")

from ..utils import safe_float, safe_int, parse_hichip_anchor2


# =============================================================================
# LOOP FILE LOADING & NORMALIZATION
# =============================================================================

def load_loops_unified(path: Path) -> pd.DataFrame:
    """
    Normalize different HiChIP loop formats into a standard DataFrame:
        chr1, start1, end1, chr2, start2, end2, counts, score, n_reps
    
    Handles:
    - BEDPE with header (HiChIPper format)
    - BEDPE without header
    - LoopCatalog 4-column BED format
    """
    path = Path(path)
    ext = path.suffix.lower()

    if ext == ".bedpe":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep="\t", header=None)

    ncols = df.shape[1]

    # BEDPE with header
    if {"chr1", "start1", "end1", "chr2", "start2", "end2"}.issubset(df.columns):
        return _normalize_bedpe_with_header(df, path)

    # BEDPE without header (9 columns)
    if ncols == 9:
        df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "score", "n_reps"]
        return _finalize_loop_types(df, path)

    # BEDPE without header (8 columns, no score)
    if ncols == 8:
        df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "n_reps"]
        df["score"] = np.nan
        return _finalize_loop_types(df, path)

    # LoopCatalog 4-column BED
    if ncols == 4:
        return _parse_loopcatalog_bed(df, path)

    raise ValueError(f"Unknown loop format with {ncols} columns in {path}")


def _normalize_bedpe_with_header(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Normalize BEDPE with header (HMEC/HiChIPper format)."""
    coord_cols = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
    missing = [c for c in coord_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing coordinate columns in {path}: {missing}")

    if "PETs" in df.columns:
        counts = df["PETs"].astype(float)
    elif "weight" in df.columns:
        counts = df["weight"].astype(float)
    else:
        counts = pd.Series([np.nan] * len(df))

    out = df.copy()
    out["counts"] = counts
    
    if "n_reps" not in out.columns:
        out["n_reps"] = 1
    if "score" not in out.columns:
        out["score"] = np.nan

    return _finalize_loop_types(out, path)


def _parse_loopcatalog_bed(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Parse LoopCatalog 4-column .bed format."""
    df.columns = ["chr1", "start1", "end1", "anchor2"]
    
    parsed = df["anchor2"].apply(lambda x: parse_hichip_anchor2(x, str(path)))
    df["chr2"] = parsed.str[0]
    df["start2"] = parsed.str[1].astype(int)
    df["end2"] = parsed.str[2].astype(int)
    df["score"] = parsed.str[3].astype(float)
    df["counts"] = np.nan
    df["n_reps"] = 1
    
    return _finalize_loop_types(df, path)


def _finalize_loop_types(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Ensure final types and column presence for loop data."""
    out = df.copy()
    
    for col in ["counts", "score"]:
        if col not in out.columns:
            out[col] = np.nan
    if "n_reps" not in out.columns:
        out["n_reps"] = 1

    out["start1"] = out["start1"].astype(int)
    out["end1"] = out["end1"].astype(int)
    out["start2"] = out["start2"].astype(int)
    out["end2"] = out["end2"].astype(int)
    out["counts"] = out["counts"].astype(float)
    out["score"] = out["score"].astype(float)
    out["n_reps"] = out["n_reps"].astype(int)

    return out[["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "score", "n_reps"]]


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def find_hichip_loops_file(
    hichip_dir: Path,
    cell_type: str,
    experiment: str = "H3K27ac",
) -> Optional[Path]:
    """Find the loops file for a given cell type and experiment."""
    cell_type_dir = Path(hichip_dir) / cell_type
    
    if not cell_type_dir.exists():
        return None
    
    # Try common naming patterns
    candidates = [
        cell_type_dir / f"{cell_type}_{experiment}_loops_hg38.bedpe",
        cell_type_dir / f"{cell_type}_{experiment}_loops_hg38.bed",
    ]
    
    for c in candidates:
        if c.exists():
            return c

    # Fall back to glob
    pattern = str(cell_type_dir / "*_loops_hg38.bed*")
    matches = glob.glob(pattern)
    if matches:
        return Path(matches[0])

    return None


# =============================================================================
# ANCHOR BUILDING
# =============================================================================

def build_anchors(loops: pd.DataFrame) -> pd.DataFrame:
    """Build anchor DataFrame from loops (one row per anchor A/B)."""
    if "loop_id" not in loops.columns:
        loops = loops.copy()
        loops["loop_id"] = np.arange(len(loops))

    anchor_a = loops[["chr1", "start1", "end1", "loop_id", "counts", "score", "n_reps"]].copy()
    anchor_a.rename(columns={"chr1": "chr", "start1": "start", "end1": "end"}, inplace=True)
    anchor_a["anchor"] = "A"

    anchor_b = loops[["chr2", "start2", "end2", "loop_id", "counts", "score", "n_reps"]].copy()
    anchor_b.rename(columns={"chr2": "chr", "start2": "start", "end2": "end"}, inplace=True)
    anchor_b["anchor"] = "B"

    return pd.concat([anchor_a, anchor_b], ignore_index=True)


# =============================================================================
# OVERLAP COMPUTATION (requires pyranges)
# =============================================================================

def compute_hichip_overlaps(
    loops: pd.DataFrame,
    ccres: pd.DataFrame,
    genes: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute anchor overlaps with cCREs and promoters.
    
    Requires pyranges to be installed.
    
    Returns:
        tuple: (anchor_metrics, overlap_ccre, overlap_prom)
    """
    if not HAS_PYRANGES:
        raise ImportError("pyranges is required for HiChIP overlap computation")
    
    loops = loops.copy()
    if "loop_id" not in loops.columns:
        loops["loop_id"] = np.arange(len(loops))
    
    anchors = build_anchors(loops)
    
    # Separate metrics (contains NaN/inf that PyRanges can't handle)
    anchor_metrics = anchors[["loop_id", "anchor", "counts", "score", "n_reps"]].drop_duplicates()
    anchors_for_pr = anchors[["chr", "start", "end", "loop_id", "anchor"]].copy()
    
    # Prepare promoters
    promoters = genes[["chrom", "prom_start", "prom_end", "gene_id", "gene_name"]].copy()
    promoters.rename(columns={"prom_start": "start", "prom_end": "end"}, inplace=True)
    
    # Prepare cCREs
    ccre_for_pr = ccres[["chrom", "start", "end", "cCRE_id"]].copy()
    
    # Convert to PyRanges format
    anchors_pr = anchors_for_pr.rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"})
    promoters_pr = promoters.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    ccre_pr = ccre_for_pr.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    
    # Compute overlaps
    gr_anchors = pr.PyRanges(anchors_pr)
    gr_promoters = pr.PyRanges(promoters_pr)
    gr_ccres = pr.PyRanges(ccre_pr)
    
    overlap_prom = gr_anchors.join(gr_promoters).df
    overlap_ccre = gr_anchors.join(gr_ccres).df
    
    # Rejoin metrics
    overlap_prom = overlap_prom.merge(anchor_metrics, on=["loop_id", "anchor"], how="left")
    overlap_ccre = overlap_ccre.merge(anchor_metrics, on=["loop_id", "anchor"], how="left")
    
    return anchor_metrics, overlap_ccre, overlap_prom


# =============================================================================
# DICT BUILDERS
# =============================================================================

def build_ccre_hichip_dict(
    overlap_ccre: pd.DataFrame,
    loops: pd.DataFrame,
) -> Dict[str, Dict[str, Any]]:
    """
    Build per-cCRE HiChIP summary.
    
    Returns:
        {cCRE_id: {n_loops, max_counts, max_score, loops: [...]}}
    """
    if "loop_id" not in loops.columns:
        loops = loops.copy()
        loops["loop_id"] = np.arange(len(loops))
    
    # Build loop coordinates lookup
    loop_coords = (
        loops.set_index("loop_id")
        .apply(lambda r: {
            "anchor_a": {"chr": r["chr1"], "start": int(r["start1"]), "end": int(r["end1"])},
            "anchor_b": {"chr": r["chr2"], "start": int(r["start2"]), "end": int(r["end2"])},
        }, axis=1)
        .to_dict()
    )
    
    # Build loop -> cCREs mapping
    loop_to_ccres = (
        overlap_ccre.groupby("loop_id")["cCRE_id"]
        .apply(lambda s: set(s.astype(str)))
        .to_dict()
    )
    
    per_ccre = {}
    
    for ccre_id, df_cc in overlap_ccre.groupby("cCRE_id"):
        seen_loops = set()
        loops_list = []
        max_counts = None
        max_score = None
        ccre_id_str = str(ccre_id)
        
        for _, row in df_cc.iterrows():
            lid = row["loop_id"]
            if lid in seen_loops:
                continue
            seen_loops.add(lid)
            
            cnt = safe_float(row["counts"])
            scr = safe_float(row["score"])
            
            if cnt is not None:
                max_counts = cnt if max_counts is None else max(max_counts, cnt)
            if scr is not None:
                max_score = scr if max_score is None else max(max_score, scr)
            
            base = {
                "loop_id": int(lid),
                "counts": cnt,
                "score": scr,
                "n_reps": int(row["n_reps"]),
                "partner_cCREs": sorted(loop_to_ccres.get(lid, set()) - {ccre_id_str}),
            }
            
            coords = loop_coords.get(lid)
            if coords:
                base.update(coords)
            
            loops_list.append(base)
        
        if loops_list:
            per_ccre[ccre_id_str] = {
                "n_loops": len(loops_list),
                "max_counts": max_counts,
                "max_score": max_score,
                "loops": loops_list,
            }
    
    return per_ccre


def build_ccre_gene_hichip_map(
    overlap_ccre: pd.DataFrame,
    overlap_prom: pd.DataFrame,
) -> Dict[str, Dict[str, List[Dict]]]:
    """
    Build cCRE -> gene -> [loops] mapping for opposite-anchor pairs.
    
    Returns:
        {cCRE_id: {gene_name: [{loop_id, counts, score, ...}, ...]}}
    """
    cc = overlap_ccre.rename(columns={"anchor": "anchor_ccre"})
    cp = overlap_prom.rename(columns={"anchor": "anchor_prom"})
    
    merged = cc.merge(cp, on="loop_id", suffixes=("_ccre", "_prom"))
    # Keep only opposite-anchor pairs (cCRE on A, gene on B or vice versa)
    merged_ep = merged[merged["anchor_ccre"] != merged["anchor_prom"]]
    
    cg = defaultdict(lambda: defaultdict(list))
    for _, r in merged_ep.iterrows():
        cg[r["cCRE_id"]][r["gene_name"]].append({
            "loop_id": r["loop_id"],
            "counts": safe_float(r["counts_ccre"]),
            "score": safe_float(r["score_ccre"]),
            "n_reps": r["n_reps_ccre"],
            "anchor_ccre": r["anchor_ccre"],
            "anchor_prom": r["anchor_prom"],
        })
    
    return dict(cg)


# =============================================================================
# PUBLIC API
# =============================================================================

def build_hichip_links(
    hichip_dir: Path,
    cell_types: List[str],
    ccres: pd.DataFrame,
    genes: pd.DataFrame,
    experiment: str = "H3K27ac",
) -> Dict[str, Dict[str, Any]]:
    """
    Build HiChIP links for multiple cell types.
    
    Returns:
        {cell_type: {cCRE_id: {n_loops, max_counts, max_score, loops: [...]}}}
    """
    if not HAS_PYRANGES:
        print("Warning: pyranges not installed. Skipping HiChIP integration.")
        return {}
    
    print("Building HiChIP links...")
    all_hichip = {}
    
    for cell_type in cell_types:
        loops_path = find_hichip_loops_file(hichip_dir, cell_type, experiment)
        
        if loops_path is None:
            print(f"  {cell_type}: No loops file found, skipping")
            continue
        
        print(f"  {cell_type}: Loading {loops_path}")
        
        try:
            loops = load_loops_unified(loops_path)
            loops["loop_id"] = np.arange(len(loops))
            
            _, overlap_ccre, overlap_prom = compute_hichip_overlaps(loops, ccres, genes)
            
            per_ccre = build_ccre_hichip_dict(overlap_ccre, loops)
            ccre_gene_map = build_ccre_gene_hichip_map(overlap_ccre, overlap_prom)
            
            all_hichip[cell_type] = {
                "per_ccre": per_ccre,
                "ccre_gene_map": ccre_gene_map,
            }
            
            print(f"    Found {len(per_ccre)} cCREs with loops")
            
        except Exception as e:
            print(f"  {cell_type}: Error processing - {e}")
            continue
    
    return all_hichip


def integrate_hichip_to_element_table(
    elem_focus: pd.DataFrame,
    hichip_data: Dict[str, Dict[str, Any]],
) -> pd.DataFrame:
    """
    Add HiChIP data to element table.
    
    Adds per-cell-type HiChIP columns and updates gene_links.
    """
    elem_focus = elem_focus.copy()
    
    for cell_type, ct_data in hichip_data.items():
        per_ccre = ct_data.get("per_ccre", {})
        ccre_gene_map = ct_data.get("ccre_gene_map", {})
        
        # Initialize cell type column if needed
        if cell_type not in elem_focus.columns:
            elem_focus[cell_type] = [{} for _ in range(len(elem_focus))]
        
        # Add hichip to cell type column
        def add_hichip_to_cell_col(row):
            cell_dict = row[cell_type] if isinstance(row[cell_type], dict) else {}
            ccre_id = str(row.get("cCRE_id", row.get("ENCODE_id", "")))
            cell_dict["hichip"] = per_ccre.get(ccre_id, {
                "n_loops": 0, "max_counts": None, "max_score": None, "loops": []
            })
            return cell_dict
        
        elem_focus[cell_type] = elem_focus.apply(add_hichip_to_cell_col, axis=1)
        
        # Update gene_links with HiChIP evidence
        def update_gene_links_hichip(row):
            gene_links = row.get("gene_links", {})
            if not isinstance(gene_links, dict):
                return gene_links
            
            ccre_id = str(row.get("cCRE_id", row.get("ENCODE_id", "")))
            per_gene_loops = ccre_gene_map.get(ccre_id, {})
            
            for gene_name, gdict in gene_links.items():
                hichip = gdict.setdefault("hichip", {})
                loops_for_pair = per_gene_loops.get(gene_name, [])
                
                if not loops_for_pair:
                    continue
                
                counts_vals = [l.get("counts") for l in loops_for_pair if l.get("counts") is not None]
                score_vals = [l.get("score") for l in loops_for_pair if l.get("score") is not None]
                
                hichip[cell_type] = {
                    "n_loops": len(loops_for_pair),
                    "max_counts": max(counts_vals) if counts_vals else None,
                    "max_score": max(score_vals) if score_vals else None,
                    "loops": loops_for_pair,
                }
            
            return gene_links
        
        if "gene_links" in elem_focus.columns:
            elem_focus["gene_links"] = elem_focus.apply(update_gene_links_hichip, axis=1)
    
    return elem_focus
