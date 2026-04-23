"""
BED interval generation for SV FASTA extraction.

Functions for:
- Creating BED intervals from processed SVs
- Generating flank regions around breakpoints
- Creating overlap regions with regulatory elements
"""

import ast
import os
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd


# =============================================================================
# BED INTERVAL BUILDING
# =============================================================================

def build_sv_flanks_and_overlaps_bed(
    sv_df: pd.DataFrame,
    output_path: Optional[Path] = None,
    flank: int = 150,
) -> pd.DataFrame:
    """
    Build BED intervals for SV flanks and regulatory element overlaps.
    
    Creates intervals for:
    - Flanking regions around breakpoints (for motif scanning)
    - Overlap regions with regulatory elements (for motif scanning)
    
    Args:
        sv_df: Processed SV DataFrame with elem_hits column
        output_path: Optional path to write BED file
        flank: Flank size in bp around breakpoints
    
    Returns:
        DataFrame with BED intervals (chrom, start, end, name)
    """
    mask = sv_df["SVTYPE"].isin(["DEL", "DUP", "INS", "BND"]) & (sv_df["elem_hits"] != "[]")
    df = sv_df[mask].copy()
    
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name"])

    # Normalize chromosome names
    if "chr" in str(df["chrom"].iloc[0]):
        df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)

    records = []

    for idx, row in df.iterrows():
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        svt = str(row["SVTYPE"])
        sv_id = str(row.get("id", idx))

        end_val = row.get("END", np.nan)
        has_end = not pd.isna(end_val)
        end = int(end_val) if has_end else None

        # ---- Flanks ----
        if svt == "DEL":
            # Left flank
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank_left"))

            # Right flank
            if end is not None:
                s2 = max(e1 + 1, end - flank)
                e2 = end + flank
                records.append((chrom, s2 - 1, e2, f"{sv_id}|{svt}_flank_right"))

        elif svt == "INS":
            # Left flank
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank_left"))

            # Right flank (insertion point + 1)
            center2 = pos + 1
            s2 = max(e1 + 1, center2 - flank)
            e2 = center2 + flank
            records.append((chrom, s2 - 1, e2, f"{sv_id}|{svt}_flank_right"))

        elif svt == "BND":
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank"))

        elif svt == "DUP":
            # Left flank
            s1 = max(1, pos - flank)
            e1 = pos + flank
            records.append((chrom, s1 - 1, e1, f"{sv_id}|{svt}_flank_left"))

            # Right flank
            if end is not None:
                s2 = max(e1 + 1, end - flank)
                e2 = end + flank
                records.append((chrom, s2 - 1, e2, f"{sv_id}|{svt}_flank_right"))

        # ---- Element overlaps ----
        try:
            hits = ast.literal_eval(row["elem_hits"]) if isinstance(row["elem_hits"], str) else row["elem_hits"]
        except Exception:
            continue

        if not isinstance(hits, list):
            continue

        for h in hits:
            if not isinstance(h, dict):
                continue
            if not h.get("overlaps_flag", 0):
                continue

            e_chrom = str(h["chrom"]).replace("chr", "")
            o_start = int(h["overlap_start"])
            o_end = int(h["overlap_end"])

            elem_id = str(h["elem_id"])
            elem_type = str(h.get("elem_type", ""))
            hit_side = str(h.get("hit_side", ""))

            name = f"{sv_id}|{svt}|elem:{elem_id}|type:{elem_type}|hit_side:{hit_side}"
            records.append((e_chrom, o_start - 1, o_end, name))

    bed = pd.DataFrame.from_records(records, columns=["chrom", "start", "end", "name"])
    # Identical BED lines (e.g. duplicate elem_hits entries) would multiply bedtools
    # intersect output in 06; one genomic interval × one name should appear once.
    if not bed.empty:
        bed = bed.drop_duplicates(ignore_index=True)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        bed.to_csv(output_path, sep="\t", header=False, index=False)
        print(f"Wrote BED file: {output_path}")

    return bed


def create_bed_from_sv_csv(
    csv_path: Path,
    output_dir: Path,
    flank: int = 150,
) -> pd.DataFrame:
    """
    Create BED intervals from a processed SV CSV file.
    
    Args:
        csv_path: Path to SV CSV file
        output_dir: Directory for output BED file
        flank: Flank size in bp
    
    Returns:
        BED DataFrame
    """
    sv_df = pd.read_csv(csv_path)
    file_name = csv_path.stem
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_bed = output_dir / f"{file_name}_flanks_and_overlaps.bed"
    
    return build_sv_flanks_and_overlaps_bed(sv_df, out_bed, flank)


def create_beds_from_directory(
    csv_dir: Path,
    output_dir: Path,
    flank: int = 150,
) -> List[Path]:
    """
    Create BED files for all CSV files in a directory.
    
    Args:
        csv_dir: Directory containing SV CSV files
        output_dir: Directory for output BED files
        flank: Flank size in bp
    
    Returns:
        List of created BED file paths
    """
    csv_dir = Path(csv_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    created_files = []
    
    for fname in os.listdir(csv_dir):
        if fname.endswith(".csv") or fname.endswith(".tsv"):
            fpath = csv_dir / fname
            print(f"Processing: {fpath}")
            
            bed = create_bed_from_sv_csv(fpath, output_dir, flank)
            
            if not bed.empty:
                out_path = output_dir / f"{fpath.stem}_flanks_and_overlaps.bed"
                created_files.append(out_path)
    
    return created_files


# =============================================================================
# BED UTILITIES
# =============================================================================

def merge_bed_files(
    bed_paths: List[Path],
    output_path: Path,
) -> pd.DataFrame:
    """
    Merge multiple BED files into one.
    
    Args:
        bed_paths: List of BED file paths
        output_path: Output path for merged BED
    
    Returns:
        Merged BED DataFrame
    """
    dfs = []
    for path in bed_paths:
        if path.exists():
            df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "name"])
            dfs.append(df)
    
    if not dfs:
        return pd.DataFrame(columns=["chrom", "start", "end", "name"])
    
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.drop_duplicates()
    merged = merged.sort_values(["chrom", "start", "end"])
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", header=False, index=False)
    
    return merged


def filter_bed_by_region(
    bed_df: pd.DataFrame,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> pd.DataFrame:
    """
    Filter BED intervals by genomic region.
    """
    mask = pd.Series([True] * len(bed_df), index=bed_df.index)
    
    if chrom is not None:
        mask &= bed_df["chrom"] == chrom.replace("chr", "")
    if start is not None:
        mask &= bed_df["end"] >= start
    if end is not None:
        mask &= bed_df["start"] <= end
    
    return bed_df[mask].copy()


def summarize_bed(bed_df: pd.DataFrame) -> Dict:
    """
    Get summary statistics for BED intervals.
    """
    if bed_df.empty:
        return {"n_intervals": 0, "total_bp": 0}
    
    lengths = bed_df["end"] - bed_df["start"]
    
    return {
        "n_intervals": len(bed_df),
        "total_bp": int(lengths.sum()),
        "mean_length": float(lengths.mean()),
        "median_length": float(lengths.median()),
        "min_length": int(lengths.min()),
        "max_length": int(lengths.max()),
        "n_chromosomes": bed_df["chrom"].nunique(),
    }
