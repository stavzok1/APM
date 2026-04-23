"""
cCRE loading and cell-line signal integration.

Functions for:
- Loading ENCODE cCREs from CSV
- Integrating per-cell-line signal tracks (H3K27ac, CTCF, etc.)
"""

import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import numpy as np
import pandas as pd

from ..biosample_names import canonical_ccre_signal_column_name
from ..utils import harmonize_chrom_column, extract_primary_class


# =============================================================================
# LOADING
# =============================================================================

def load_ccres(path: Path) -> pd.DataFrame:
    """
    Load cCREs from ENCODE SCREEN CSV.
    
    Normalizes chromosome names and extracts primary element class.
    """
    df = pd.read_csv(path)
    df, _ = harmonize_chrom_column(df)
    
    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    # Extract primary class from type column
    if "type" in df.columns:
        df["raw_type"] = df["type"].astype(str)
        df["type"] = df["raw_type"].apply(extract_primary_class)
    
    # Compute center
    df["center"] = ((df["start"] + df["end"]) // 2).astype("Int64")
    
    return df


# =============================================================================
# CELL-LINE SIGNAL INTEGRATION
# =============================================================================

def add_cell_line_signals(
    ccres: pd.DataFrame,
    cell_line_dir: Path,
    cell_line: str,
    wanted_signals: List[str],
    ccre_id_col_in_bed: int = 3,
) -> pd.DataFrame:
    """
    Add cell-line-specific signal values to cCRE table.
    
    For each element, creates a dict in column `cell_line`:
    {
        "in_{cell_line}": bool,
        "H3K27ac": float or None,
        "CTCF": float or None,
        ...
    }
    
    Args:
        ccres: cCRE DataFrame with ENCODE_id column
        cell_line_dir: Directory containing BED and signal files for cell line
        cell_line: Name of cell line (used as column name)
        wanted_signals: List of signal names to extract
        ccre_id_col_in_bed: Column index (0-based) of cCRE ID in BED file
    
    Returns:
        DataFrame with new cell_line column containing nested dict for cCREs
        present in the cell line; non-present rows are set to None to avoid
        building massive per-row dicts for the entire genome.
    """
    df = ccres  # mutate in-place to avoid 1M-row copies
    cell_line_dir = Path(cell_line_dir)
    
    # 1. Find and load the BED file to get cCREs present in this cell line
    bed_files = [f for f in os.listdir(cell_line_dir) if f.endswith(".bed")]
    if len(bed_files) != 1:
        raise ValueError(
            f"Expected exactly one .bed file in {cell_line_dir}, found {len(bed_files)}"
        )
    
    bed_path = cell_line_dir / bed_files[0]
    bed = pd.read_csv(bed_path, sep="\t", header=None)
    ccre_ids_in_cell_line = set(bed.iloc[:, ccre_id_col_in_bed].astype(str).tolist())

    encode_ids = df["ENCODE_id"].astype(str)
    in_mask = encode_ids.isin(ccre_ids_in_cell_line).to_numpy()

    # 2. Load signal maps (ENCODE_id -> value) for wanted signals
    signal_maps: Dict[str, pd.Series] = {}
    for file in os.listdir(cell_line_dir):
        base, ext = os.path.splitext(file)
        if base in wanted_signals:
            signal_path = cell_line_dir / file
            signal_df = pd.read_csv(signal_path, sep="\t", header=None)
            signal_map = signal_df.set_index(signal_df.columns[0])[signal_df.columns[1]]
            signal_maps[base] = pd.to_numeric(signal_map, errors="coerce")

    # 3. Build dicts only for present cCREs; store None for non-present
    in_key = f"in_{cell_line}"
    idx_present = np.where(in_mask)[0]
    out = [None] * len(df)

    if len(idx_present) > 0:
        # Precompute signal vectors for present IDs only (keeps memory bounded)
        present_ids = encode_ids.iloc[idx_present]
        sig_vectors: Dict[str, np.ndarray] = {}
        for sig_name in wanted_signals:
            if sig_name in signal_maps:
                vals = present_ids.map(signal_maps[sig_name]).to_numpy()
                sig_vectors[sig_name] = vals
            else:
                sig_vectors[sig_name] = np.full(len(idx_present), np.nan, dtype="float64")

        for j, i in enumerate(idx_present):
            d: Dict[str, Any] = {in_key: True}
            for sig_name in wanted_signals:
                v = sig_vectors[sig_name][j]
                d[sig_name] = None if (v is None or (isinstance(v, float) and pd.isna(v)) or pd.isna(v)) else float(v)
            out[i] = d

    df[cell_line] = out
    
    return df


def add_multiple_cell_line_signals(
    ccres: pd.DataFrame,
    cell_lines_base_dir: Path,
    cell_lines: List[str],
    wanted_signals: List[str],
    ccre_id_col_in_bed: int = 3,
) -> pd.DataFrame:
    """
    Add signals from multiple cell lines.
    
    Args:
        ccres: cCRE DataFrame
        cell_lines_base_dir: Base directory containing cell line subdirectories
        cell_lines: List of cell line names to process
        wanted_signals: Signal names to extract
        ccre_id_col_in_bed: Column index of cCRE ID in BED files
    
    Returns:
        DataFrame with one column per cell line containing signal dicts
    """
    cell_lines_base_dir = Path(cell_lines_base_dir)
    
    for cell_line_disk in cell_lines:
        cell_line_dir = cell_lines_base_dir / cell_line_disk
        column_key = canonical_ccre_signal_column_name(str(cell_line_disk))
        if cell_line_dir.exists():
            print(f"Processing cell line signals: disk={cell_line_disk!r} → column={column_key!r}")
            ccres = add_cell_line_signals(
                ccres,
                cell_line_dir,
                column_key,
                wanted_signals,
                ccre_id_col_in_bed,
            )
        else:
            print(f"Warning: Cell line directory not found: {cell_line_dir}")
    
    return ccres


# =============================================================================
# SIGNAL EXTRACTION HELPERS
# =============================================================================

def extract_signal_value(
    ccres: pd.DataFrame,
    cell_line: str,
    signal_name: str,
) -> pd.Series:
    """
    Extract a specific signal value from the nested cell line dict.
    
    Useful for downstream analysis or filtering.
    """
    if cell_line not in ccres.columns:
        return pd.Series([None] * len(ccres), index=ccres.index)
    
    return ccres[cell_line].apply(
        lambda d: d.get(signal_name) if isinstance(d, dict) else None
    )


def filter_ccres_by_signal(
    ccres: pd.DataFrame,
    cell_line: str,
    signal_name: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
    require_present: bool = True,
) -> pd.DataFrame:
    """
    Filter cCREs based on signal values in a specific cell line.
    
    Args:
        ccres: cCRE DataFrame with cell line columns
        cell_line: Cell line to filter on
        signal_name: Signal to filter (e.g., "H3K27ac")
        min_value: Minimum signal value (inclusive)
        max_value: Maximum signal value (inclusive)
        require_present: If True, only keep cCREs present in cell line
    
    Returns:
        Filtered DataFrame
    """
    if cell_line not in ccres.columns:
        raise ValueError(f"Cell line '{cell_line}' not found in DataFrame")
    
    mask = pd.Series([True] * len(ccres), index=ccres.index)
    
    if require_present:
        in_key = f"in_{cell_line}"
        mask &= ccres[cell_line].apply(
            lambda d: d.get(in_key, False) if isinstance(d, dict) else False
        )
    
    signal_values = extract_signal_value(ccres, cell_line, signal_name)
    
    if min_value is not None:
        mask &= signal_values >= min_value
    
    if max_value is not None:
        mask &= signal_values <= max_value
    
    return ccres[mask].copy()
