"""
Motif scanning integration with FIMO.

Functions for:
- Extracting selected motifs from MEME databases
- Parsing FIMO output
- Recombining motif hits with SV tables
"""

import ast
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import List, Set, Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# MOTIF EXTRACTION
# =============================================================================

def extract_selected_motifs(
    input_meme: Path,
    output_meme: Path,
    target_symbols: Optional[List[str]] = None,
) -> Tuple[int, List[str]]:
    """
    Extract motifs matching target TF symbols from a MEME file.
    
    Args:
        input_meme: Path to input MEME file
        output_meme: Path to output MEME file
        target_symbols: List of TF symbols to match (case-insensitive substring)
    
    Returns:
        Tuple of (n_selected, missing_symbols)
    """
    if target_symbols is None:
        target_symbols = [
            "STAT1", "STAT2", "STAT3",
            "IRF1", "IRF2", "IRF3",
            "NFKB1", "NFKB2", "RELA", "p65", "RELB", "REL",
            "FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3",
            "CTCF", "CTCFL",
            "NLRC5", "MYC",
        ]
    
    target_patterns = [re.compile(sym, re.IGNORECASE) for sym in target_symbols]
    
    def motif_matches_any(motif_name: str) -> bool:
        return any(pat.search(motif_name) for pat in target_patterns)
    
    input_meme = Path(input_meme)
    output_meme = Path(output_meme)
    
    in_text = input_meme.read_text()
    
    # Split header vs motifs
    parts = in_text.split("\nMOTIF ")
    if len(parts) == 1:
        raise ValueError("Could not find any 'MOTIF ' blocks in the input MEME file.")
    
    header = parts[0]
    motif_blocks_raw = parts[1:]
    
    selected_blocks = []
    found_symbols = set()
    
    for block_raw in motif_blocks_raw:
        block = "MOTIF " + block_raw
        first_line = block.splitlines()[0]
        tokens = first_line.split()
        
        if len(tokens) < 2:
            continue
        
        motif_name = tokens[1]
        
        if motif_matches_any(motif_name):
            selected_blocks.append(block)
            for sym in target_symbols:
                if re.search(sym, motif_name, re.IGNORECASE):
                    found_symbols.add(sym)
    
    # Build output
    out_lines = [header.rstrip(), ""]
    for block in selected_blocks:
        out_lines.append(block.rstrip())
        out_lines.append("")
    
    out_text = "\n".join(out_lines).rstrip() + "\n"
    
    output_meme.parent.mkdir(parents=True, exist_ok=True)
    output_meme.write_text(out_text)
    
    missing = [sym for sym in target_symbols if sym not in found_symbols]
    
    print(f"Selected {len(selected_blocks)} motif blocks into {output_meme}")
    if missing:
        print(f"No motif found for: {', '.join(missing)}")
    
    return len(selected_blocks), missing


# =============================================================================
# FIMO OUTPUT PARSING
# =============================================================================

def parse_fimo_tsv_to_bed(
    fimo_tsv_path: Path,
    output_bed_path: Path,
) -> pd.DataFrame:
    """
    Convert FIMO TSV output to BED format.
    
    Args:
        fimo_tsv_path: Path to FIMO TSV file
        output_bed_path: Path to output BED file
    
    Returns:
        BED DataFrame
    """
    fimo_tsv_path = Path(fimo_tsv_path)
    output_bed_path = Path(output_bed_path)
    
    df = pd.read_csv(fimo_tsv_path, sep="\t")
    
    # Coerce start/stop to numeric
    df["start_num"] = pd.to_numeric(df["start"], errors="coerce")
    df["stop_num"] = pd.to_numeric(df["stop"], errors="coerce")
    df = df.dropna(subset=["start_num", "stop_num"]).copy()
    
    # Clean chromosome names
    chrom = df["sequence_name"].astype(str)
    chrom = chrom.str.replace(r"\.0$", "", regex=True)
    chrom = chrom.str.replace(r"^chr", "", regex=True)
    
    bed = pd.DataFrame({
        "chrom": chrom,
        "start": (df["start_num"].astype(int) - 1).clip(lower=0),
        "end": df["stop_num"].astype(int),
        "motif_id": df["motif_id"],
        "p_value": df["p-value"],
        "q_value": df.get("q-value", np.nan),
        "score": df.get("score", np.nan),
        "strand": df["strand"],
    })
    
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    bed.to_csv(output_bed_path, sep="\t", header=False, index=False)
    
    return bed


# =============================================================================
# MOTIF RECOMBINATION
# =============================================================================

def _parse_elem_hits_cell(x) -> List[Dict]:
    """Parse elem_hits cell into list of dicts."""
    if isinstance(x, list):
        return x
    if pd.isna(x) or x in ("", "[]"):
        return []
    try:
        return ast.literal_eval(x)
    except Exception:
        return []


def _parse_flank_side(name_tail: str) -> str:
    """Extract flank side from name tail."""
    parts = name_tail.split("_")
    if len(parts) >= 3 and parts[-2] == "flank":
        return parts[-1]
    return parts[-1]


def _extract_elem_id_from_name(name: str) -> Optional[str]:
    """Extract elem_id from BED name."""
    for part in name.split("|"):
        if part.startswith("elem:"):
            return part[len("elem:"):]
    return None


def recombine_sv_fimo(
    merged_bed_path: Path,
    sv_csv_path: Path,
    output_csv_path: Path,
) -> None:
    """
    Recombine FIMO motif hits into processed SV CSV.
    
    Args:
        merged_bed_path: Path to merged FIMO-SV intersection BED
        sv_csv_path: Path to processed SV CSV
        output_csv_path: Path for output CSV with motif hits
    """
    merged_bed_path = Path(merged_bed_path)
    sv_csv_path = Path(sv_csv_path)
    output_csv_path = Path(output_csv_path)
    
    if not merged_bed_path.is_file():
        print(f"[WARN] Merged BED not found: {merged_bed_path}")
        return
    if not sv_csv_path.is_file():
        print(f"[WARN] SV CSV not found: {sv_csv_path}")
        return
    
    # Load SV CSV
    s = pd.read_csv(sv_csv_path)
    
    if "id" not in s.columns or "elem_hits" not in s.columns or "pos" not in s.columns:
        print(f"[ERROR] CSV missing required columns: {sv_csv_path}")
        return
    
    # Prepare flank_motif_hits column
    if "flank_motif_hits" not in s.columns:
        s["flank_motif_hits"] = [[] for _ in range(len(s))]
    else:
        s["flank_motif_hits"] = s["flank_motif_hits"].apply(
            lambda x: [] if (pd.isna(x) or x in ("", "[]")) else (
                x if isinstance(x, list) else ast.literal_eval(x)
            )
        )
    
    # Parse elem_hits
    s["elem_hits_parsed"] = s["elem_hits"].apply(_parse_elem_hits_cell)
    
    # Build SV id index
    id_to_indices = defaultdict(list)
    for idx, sv_id in s["id"].items():
        id_to_indices[str(sv_id)].append(idx)
    
    # Load merged BED
    if os.path.getsize(merged_bed_path) == 0:
        print(f"[INFO] Empty BED file: {merged_bed_path}")
        s.drop(columns=["elem_hits_parsed"], inplace=True)
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        s.to_csv(output_csv_path, index=False)
        return
    
    mcols = [
        "sv_chrom", "sv_start", "sv_end", "sv_name",
        "motif_chrom", "motif_start", "motif_end",
        "motif_id", "p_value", "q_value", "score", "strand",
    ]
    merged = pd.read_csv(
        merged_bed_path, sep="\t", header=None, names=mcols,
        dtype={"sv_chrom": str, "motif_chrom": str}
    )
    
    if merged.empty:
        s.drop(columns=["elem_hits_parsed"], inplace=True)
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        s.to_csv(output_csv_path, index=False)
        return
    
    # Separate flank vs regulatory hits
    is_flank = merged["sv_name"].str.contains("flank", na=False)
    flank_rows = merged[is_flank].copy()
    reg_rows = merged[~is_flank].copy()
    
    # Collect flank motif hits
    flank_hits_by_sv = defaultdict(list)
    for _, row in flank_rows.iterrows():
        name = row["sv_name"]
        sv_id = name.split("|")[0]
        tail = name.split("|", 1)[1] if "|" in name else ""
        flank_side = _parse_flank_side(tail)
        
        motif_start = int(row["motif_start"])
        motif_end = int(row["motif_end"])
        motif_id = str(row["motif_id"])
        
        try:
            p_val = float(row["p_value"])
        except Exception:
            p_val = np.nan
        try:
            q_val = float(row["q_value"])
        except Exception:
            q_val = np.nan
        try:
            score = float(row["score"])
        except Exception:
            score = np.nan
        
        strand = str(row["strand"])
        tf_short = motif_id.split(".", 1)[0] if "." in motif_id else motif_id
        
        hit_dict = {
            "flank_side": flank_side,
            "start": motif_start,
            "end": motif_end,
            "TF": tf_short,
            "motif_id": motif_id,
            "score": score,
            "p_value": p_val,
            "q_value": q_val,
            "strand": strand,
        }
        flank_hits_by_sv[sv_id].append(hit_dict)
    
    # Collect regulatory motif hits
    elem_motif_hits = defaultdict(lambda: defaultdict(list))
    for _, row in reg_rows.iterrows():
        name = row["sv_name"]
        sv_id = name.split("|")[0]
        elem_id = _extract_elem_id_from_name(name)
        if elem_id is None:
            continue
        
        motif_start = int(row["motif_start"])
        motif_end = int(row["motif_end"])
        motif_id = str(row["motif_id"])
        
        try:
            p_val = float(row["p_value"])
        except Exception:
            p_val = np.nan
        try:
            q_val = float(row["q_value"])
        except Exception:
            q_val = np.nan
        try:
            score = float(row["score"])
        except Exception:
            score = np.nan
        
        strand = str(row["strand"])
        tf_short = motif_id.split(".", 1)[0] if "." in motif_id else motif_id
        
        hit_dict = {
            "start": motif_start,
            "end": motif_end,
            "TF": tf_short,
            "motif_id": motif_id,
            "score": score,
            "p_value": p_val,
            "q_value": q_val,
            "strand": strand,
        }
        elem_motif_hits[sv_id][elem_id].append(hit_dict)
    
    # Attach flank hits
    for sv_id, hits in flank_hits_by_sv.items():
        idxs = id_to_indices.get(sv_id, [])
        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None
            
            existing = s.at[idx, "flank_motif_hits"]
            if not isinstance(existing, list):
                existing = []
            
            hits_with_dist = []
            for h in hits:
                h_copy = dict(h)
                if bp_pos is not None:
                    center = (h_copy["start"] + h_copy["end"]) // 2
                    h_copy["distance_to_pos"] = center - bp_pos
                else:
                    h_copy["distance_to_pos"] = None
                hits_with_dist.append(h_copy)
            
            existing.extend(hits_with_dist)
            s.at[idx, "flank_motif_hits"] = existing
    
    # Attach regulatory motif hits
    for sv_id, elem_dict in elem_motif_hits.items():
        idxs = id_to_indices.get(sv_id, [])
        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None
            
            elem_list = s.at[idx, "elem_hits_parsed"]
            if not isinstance(elem_list, list):
                continue
            
            for elem in elem_list:
                if not isinstance(elem, dict):
                    continue
                eid = elem.get("elem_id")
                if not eid:
                    continue
                
                hits_for_elem = elem_dict.get(eid)
                if not hits_for_elem:
                    elem["motif_hits"] = []
                    continue
                
                existing = elem.get("motif_hits", [])
                if not isinstance(existing, list):
                    existing = []
                
                hits_with_dist = []
                for h in hits_for_elem:
                    h_copy = dict(h)
                    if bp_pos is not None:
                        center = (h_copy["start"] + h_copy["end"]) // 2
                        h_copy["distance_to_pos"] = center - bp_pos
                    else:
                        h_copy["distance_to_pos"] = None
                    hits_with_dist.append(h_copy)
                
                existing.extend(hits_with_dist)
                elem["motif_hits"] = existing
    
    # Finalize
    s["elem_hits"] = s["elem_hits_parsed"].apply(repr)
    s.drop(columns=["elem_hits_parsed"], inplace=True)
    
    s["flank_motif_hits"] = s["flank_motif_hits"].apply(
        lambda x: "[]" if (x is None or (isinstance(x, float) and np.isnan(x)) or x == []) else repr(x)
    )
    
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(output_csv_path, index=False)
    print(f"[INFO] Wrote: {output_csv_path}")


def recombine_all_sv_fimo(
    merged_bed_dir: Path,
    sv_csv_dir: Path,
    output_csv_dir: Path,
) -> List[Path]:
    """
    Recombine FIMO hits for all SVs in a directory.
    
    Args:
        merged_bed_dir: Directory with merged FIMO-SV BED files
        sv_csv_dir: Directory with processed SV CSVs
        output_csv_dir: Directory for output CSVs
    
    Returns:
        List of created output file paths
    """
    merged_bed_dir = Path(merged_bed_dir)
    sv_csv_dir = Path(sv_csv_dir)
    output_csv_dir = Path(output_csv_dir)
    output_csv_dir.mkdir(parents=True, exist_ok=True)
    
    created_files = []
    
    for fname in os.listdir(merged_bed_dir):
        if not fname.endswith("_fimo_merged.bed"):
            continue
        
        merged_bed_path = merged_bed_dir / fname
        
        # Extract base name
        base = fname[:-len("_fimo_merged.bed")]
        if base.endswith("_flanks_and_overlaps"):
            base_csv = base[:-len("_flanks_and_overlaps")]
        else:
            base_csv = base
        
        csv_name = base_csv + ".csv"
        sv_csv_path = sv_csv_dir / csv_name
        out_csv_path = output_csv_dir / csv_name
        
        recombine_sv_fimo(merged_bed_path, sv_csv_path, out_csv_path)
        
        if out_csv_path.exists():
            created_files.append(out_csv_path)
    
    return created_files


# =============================================================================
# MOTIF SUMMARY
# =============================================================================

def summarize_motif_hits(sv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize motif hits across all SVs.
    
    Returns DataFrame with TF hit counts.
    """
    tf_counts = defaultdict(int)
    
    # Count flank motif hits
    if "flank_motif_hits" in sv_df.columns:
        for hits in sv_df["flank_motif_hits"]:
            if isinstance(hits, str):
                try:
                    hits = ast.literal_eval(hits)
                except Exception:
                    continue
            if isinstance(hits, list):
                for h in hits:
                    if isinstance(h, dict):
                        tf = h.get("TF", "unknown")
                        tf_counts[tf] += 1
    
    # Count elem motif hits
    if "elem_hits" in sv_df.columns:
        for elem_hits in sv_df["elem_hits"]:
            if isinstance(elem_hits, str):
                try:
                    elem_hits = ast.literal_eval(elem_hits)
                except Exception:
                    continue
            if isinstance(elem_hits, list):
                for elem in elem_hits:
                    if isinstance(elem, dict):
                        motif_hits = elem.get("motif_hits", [])
                        if isinstance(motif_hits, list):
                            for h in motif_hits:
                                if isinstance(h, dict):
                                    tf = h.get("TF", "unknown")
                                    tf_counts[tf] += 1
    
    if not tf_counts:
        return pd.DataFrame(columns=["TF", "count"])
    
    return pd.DataFrame([
        {"TF": tf, "count": count}
        for tf, count in sorted(tf_counts.items(), key=lambda x: -x[1])
    ])
