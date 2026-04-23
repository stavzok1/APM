"""
SV intersection with unified ChIP-seq peaks.

Adds a `chip_hits` column to processed SV tables, analogous to `elem_hits`
but for TF binding evidence. Each SV gets a list of dicts, one per
ChIP peak it disrupts (overlaps or comes within `window` bp of).

ChIP integration is SV-as-disruptor only: we record peaks whose binding
could be lost / fused / brought into novel proximity by the SV. No FIMO
motif validation is performed at this stage.

Per-SVTYPE handling mirrors `map_svs_to_elements`:
  - DEL / DUP / INS: span-based intersection, hit_side="span"
  - BND: two-sided breakpoint intersection, hit_side="bp1"/"bp2"
"""

from pathlib import Path
from typing import Optional, List

import numpy as np
import pandas as pd

from ..SV.spatial_mapping import compute_signed_distance_with_overlap
from ..biosample_names import normalize_cell_line_label


# =============================================================================
# CHIP-HIT BUILDER
# =============================================================================

def _build_chip_hit(
    peak_row: pd.Series,
    sv_start: int,
    sv_end: int,
    chrom: str,
    hit_side: str,
    overlap_window: int,
) -> Optional[dict]:
    """
    Build a single chip_hit dict for a peak that passed the window filter.

    Returns None if the peak is outside `overlap_window`.
    """
    peak_start = int(peak_row["start"])
    peak_end = int(peak_row["end"])

    signed_dist, overlap_start, overlap_end = compute_signed_distance_with_overlap(
        sv_start, sv_end, peak_start, peak_end, strand="+"
    )

    if abs(signed_dist) > overlap_window:
        return None

    overlaps = overlap_start is not None and overlap_end is not None
    overlap_bp = (overlap_end - overlap_start) if overlaps else 0

    return {
        "tf": peak_row["tf"],
        "cell_type": normalize_cell_line_label(str(peak_row["cell_type"])),
        "source": peak_row["source"],
        "sample_id": (
            None
            if ("sample_id" not in peak_row.index or pd.isna(peak_row["sample_id"]))
            else str(peak_row["sample_id"])
        ),
        "score_norm": (
            None if pd.isna(peak_row["score_norm"])
            else float(peak_row["score_norm"])
        ),
        "chrom": chrom,
        "peak_start": peak_start,
        "peak_end": peak_end,
        "signed_dist": signed_dist,
        "overlap_start": overlap_start,
        "overlap_end": overlap_end,
        "overlap_bp": overlap_bp,
        "overlaps_flag": int(overlaps),
        "hit_side": hit_side,
    }


# =============================================================================
# PER-SVTYPE MAPPING
# =============================================================================

def _map_span_svs_to_chip(
    df: pd.DataFrame,
    chip: pd.DataFrame,
    svtype: str,
    window: int,
) -> pd.DataFrame:
    """Map DEL / DUP / INS (span-based) to ChIP peaks."""
    sub = df[df["SVTYPE"] == svtype].copy()
    if sub.empty:
        return df

    chip_hits_col: List[list] = []

    for _, sv in sub.iterrows():
        chrom = sv["chrom"]
        sv_start = int(sv["pos"])

        # End fallback for SVTYPEs that may lack END (e.g., INS as point)
        if "END" in sv and not pd.isna(sv["END"]):
            sv_end = int(sv["END"])
        elif "SVLEN" in sv and not pd.isna(sv["SVLEN"]):
            sv_end = sv_start + abs(int(sv["SVLEN"]))
        else:
            sv_end = sv_start

        # Coarse pre-filter: same chrom, within window
        c_chr = chip[chip["chrom"] == chrom]
        c_chr = c_chr[
            (c_chr["end"] >= sv_start - window) &
            (c_chr["start"] <= sv_end + window)
        ]

        if c_chr.empty:
            chip_hits_col.append([])
            continue

        hits = []
        for _, peak in c_chr.iterrows():
            hit = _build_chip_hit(peak, sv_start, sv_end, chrom,
                                  hit_side="span", overlap_window=window)
            if hit is not None:
                hits.append(hit)

        chip_hits_col.append(hits)

    sub["chip_hits"] = chip_hits_col
    df.loc[sub.index, "chip_hits"] = sub["chip_hits"]
    return df


def _map_bnds_to_chip(
    df: pd.DataFrame,
    chip: pd.DataFrame,
    window: int,
) -> pd.DataFrame:
    """Map BND two-sided breakpoints to ChIP peaks."""
    bnds = df[df["SVTYPE"] == "BND"].copy()
    if bnds.empty:
        return df

    chip_hits_col: List[list] = []

    for _, sv in bnds.iterrows():
        chrom1 = sv["chrom"]
        pos1 = int(sv["pos"])
        chrom2 = sv.get("bnd_remote_chrom", None)
        pos2 = sv.get("bnd_remote_pos", None)
        pos2 = int(pos2) if pd.notna(pos2) else None

        hits = []

        def _process_side(bk_chrom, bk_pos, side_label, bp_index,
                          mate_chrom, mate_pos):
            if bk_chrom is None or bk_pos is None or pd.isna(bk_chrom):
                return

            c_chr = chip[chip["chrom"] == bk_chrom]
            c_chr = c_chr[
                (c_chr["end"] >= bk_pos - window) &
                (c_chr["start"] <= bk_pos + window)
            ]
            if c_chr.empty:
                return

            for _, peak in c_chr.iterrows():
                hit = _build_chip_hit(
                    peak, bk_pos, bk_pos, bk_chrom,
                    hit_side=side_label, overlap_window=window,
                )
                if hit is None:
                    continue
                hit["bp_index"] = bp_index
                hit["bp_chrom"] = bk_chrom
                hit["bp_pos"] = bk_pos
                hit["mate_chrom"] = mate_chrom
                hit["mate_pos"] = mate_pos
                hits.append(hit)

        _process_side(chrom1, pos1, "bp1", 1, chrom2, pos2)
        _process_side(chrom2, pos2, "bp2", 2, chrom1, pos1)

        chip_hits_col.append(hits)

    bnds["chip_hits"] = chip_hits_col
    df.loc[bnds.index, "chip_hits"] = bnds["chip_hits"]
    return df


# =============================================================================
# PUBLIC API
# =============================================================================

def map_svs_to_chip(
    sv_df: pd.DataFrame,
    chip: pd.DataFrame,
    window: int = 1_000_000,
    cell_type_whitelist: Optional[List[str]] = None,
    tf_whitelist: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Add a `chip_hits` column to an SV DataFrame.

    Args:
        sv_df: Processed SV DataFrame (output of process_single_vcf).
               Must contain SVTYPE, chrom, pos, and for BND:
               bnd_remote_chrom / bnd_remote_pos.
        chip: Unified ChIP peak table from chip_loader.load_unified_chip.
              Required columns: chrom, start, end, tf, cell_type,
              source, score_norm.
        window: Maximum distance (bp) between SV span/breakpoint and a
                ChIP peak to count as a hit.
        cell_type_whitelist: Optional list of cell types to retain.
                             None = keep all.
        tf_whitelist: Optional list of TFs to retain. None = keep all.

    Returns:
        sv_df with an added `chip_hits` column (list of dicts per SV).
    """
    df = sv_df.copy()

    # Filter ChIP table once up-front (cheaper than per-SV)
    chip_f = chip.copy()
    if "cell_type" in chip_f.columns:
        chip_f["cell_type"] = chip_f["cell_type"].map(normalize_cell_line_label)
    if cell_type_whitelist is not None:
        canon_wl = {normalize_cell_line_label(x) for x in cell_type_whitelist}
        chip_f = chip_f[chip_f["cell_type"].isin(canon_wl)]
    if tf_whitelist is not None:
        chip_f = chip_f[chip_f["tf"].isin(tf_whitelist)]

    # Initialize column with empty lists (object dtype)
    df["chip_hits"] = [[] for _ in range(len(df))]

    if chip_f.empty:
        print("  [WARN] No ChIP peaks remain after filtering")
        return df

    for svtype in ("DEL", "DUP", "INS"):
        df = _map_span_svs_to_chip(df, chip_f, svtype, window)

    df = _map_bnds_to_chip(df, chip_f, window)

    return df


# =============================================================================
# BATCH ANNOTATION
# =============================================================================

def annotate_sv_csvs_with_chip(
    sv_csv_dir: Path,
    chip: pd.DataFrame,
    output_dir: Path,
    window: int = 1_000_000,
    cell_type_whitelist: Optional[List[str]] = None,
    tf_whitelist: Optional[List[str]] = None,
) -> List[Path]:
    """
    Add chip_hits column to every SV CSV in a directory.

    Intended to consume `08_neojunction_enriched/` and write to
    `09_chip_enriched/`.

    Args:
        sv_csv_dir: Directory with SV CSVs.
        chip: Unified ChIP peak table.
        output_dir: Destination directory.
        window: Distance window (bp).
        cell_type_whitelist: Optional cell-type filter.
        tf_whitelist: Optional TF filter.

    Returns:
        List of written CSV paths.
    """
    sv_csv_dir = Path(sv_csv_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for csv_file in sorted(sv_csv_dir.glob("*.csv")):
        sv_df = pd.read_csv(csv_file)
        sv_df = map_svs_to_chip(
            sv_df, chip, window=window,
            cell_type_whitelist=cell_type_whitelist,
            tf_whitelist=tf_whitelist,
        )

        # Serialize chip_hits for CSV (parquet preserves natively)
        sv_df["chip_hits"] = sv_df["chip_hits"].apply(
            lambda x: "[]" if (x is None or x == []) else repr(x)
        )

        out_path = output_dir / csv_file.name
        sv_df.to_csv(out_path, index=False)
        written.append(out_path)

        n_with_hits = (sv_df["chip_hits"] != "[]").sum()
        print(f"  Wrote: {out_path.name}  ({n_with_hits}/{len(sv_df)} SVs with ChIP hits)")

    return written
