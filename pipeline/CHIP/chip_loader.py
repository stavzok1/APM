"""
ChIP-seq peak loading and unification.

Functions for:
- Loading ENCODE ChIP-seq BED files (chrom start end name score strand
  signalValue pValue qValue peak tf cell_type)
- Loading ChIP-Atlas **narrow 6-column** BED (``chrom start end experiment tf cell_type``),
  optionally with a header row. UCSC / gffTags multi-column exports must be converted first
  using ``scripts/chip/preprocess_chip_atlas_beds.py`` (see ``pipeline.CHIP.chip_atlas_convert``).
- Unifying both into a single concatenated peak table with schema:
    chrom, start, end, tf, cell_type, source, score_norm, sample_id, cell_subtype

Both sources are assumed to be HG38. Duplicates across sources are kept
(distinguished by `source` column) for downstream filtering flexibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ..config import BIOSAMPLES
from ..utils import harmonize_chrom_column

from .chip_atlas_convert import needs_ucsc_conversion
from .chip_cell_line_policy import apply_chip_atlas_cell_line_policy

# ChIP-Atlas legacy BED sometimes leaves the cell-line column blank for an SRX
# (no co-occurring non-empty rows in the same file to backfill). Use a stable
# label so peaks remain usable in groupbys / nested summaries.
CHIP_ATLAS_UNSPECIFIED_CELL_LINE = "unspecified_cell_line"


# =============================================================================
# UNIFIED SCHEMA
# =============================================================================

UNIFIED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "tf",
    "cell_type",
    "source",
    "score_norm",
    # Identifier for the specific ChIP sample/replicate (derived from filename stem).
    # Examples: CTCF_MCF7_1, CTCF_MCF7_2, CEBPE_MCF7, CTCF_mammary_1, ...
    "sample_id",
    # Subtype / tissue-state label (ChIP-Atlas rows after ``apply_chip_atlas_cell_line_policy``).
    "cell_subtype",
]

ENCODE_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "peak", "tf", "cell_type",
]


def load_encode_bed(path: Path) -> pd.DataFrame:
    """
    Load a single ENCODE ChIP-seq BED file.

    pValue is dropped (constant -1.0 sentinel in this dataset).
    score_norm is set from signalValue.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=ENCODE_COLUMNS,
                     low_memory=False)
    df, _ = harmonize_chrom_column(df)

    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df["score_norm"] = pd.to_numeric(df["signalValue"], errors="coerce")

    df["source"] = "ENCODE"
    df["sample_id"] = path.stem
    df["cell_subtype"] = ""
    return df[UNIFIED_COLUMNS]


def load_chip_atlas_bed(path: Path) -> pd.DataFrame:
    """
    Load a single **narrow** ChIP-Atlas BED (six data columns, tab-separated).

    Optional header row ``chrom\\tstart\\tend\\texperiment\\ttf\\tcell_type`` is skipped.

    UCSC / gffTags exports must be converted first::

        .venv/bin/python3 scripts/chip/preprocess_chip_atlas_beds.py

    Peaks are assumed pre-filtered at MACS2 -10*log10(Q) >= 50 where applicable.
    No signal column is provided; score_norm is set to NaN.
    """
    path = Path(path)
    if needs_ucsc_conversion(path):
        raise ValueError(
            f"{path} still contains UCSC/gffTags rows. Run:\n"
            f"  python3 scripts/chip/preprocess_chip_atlas_beds.py "
            f"--chip-atlas-dir {path.parent}\n"
            "then rebuild the unified parquet."
        )

    expected = ["chrom", "start", "end", "experiment", "tf", "cell_type"]
    try:
        df = pd.read_csv(path, sep="\t", header=None, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    if df.empty:
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    if str(df.iloc[0, 0]).strip().lower() == "chrom":
        df = df.iloc[1:].reset_index(drop=True)
    if df.shape[1] < 6:
        raise ValueError(
            f"{path}: expected at least 6 tab columns in narrow ChIP-Atlas BED, got {df.shape[1]}"
        )
    if df.shape[1] > 6:
        df = df.iloc[:, :6].copy()
    else:
        df = df.copy()
    df.columns = expected

    df, _ = harmonize_chrom_column(df)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df["tf"] = df["tf"].astype(str).str.strip()
    df["cell_type"] = df["cell_type"].fillna("").astype(str).str.strip()
    df["cell_type"] = df["cell_type"].replace("", CHIP_ATLAS_UNSPECIFIED_CELL_LINE)

    df["score_norm"] = pd.NA
    df["source"] = "CHIP_ATLAS"
    df["sample_id"] = path.stem
    df["cell_subtype"] = ""
    return df[UNIFIED_COLUMNS]


# =============================================================================
# DIRECTORY-LEVEL LOADING
# =============================================================================

def load_source_directory(
    source_dir: Path,
    source: str,
    pattern: str = "*.bed",
) -> pd.DataFrame:
    """
    Load all BED files from a source directory and concatenate.

    Args:
        source_dir: Directory of BED files.
        source: Either "ENCODE" or "CHIP_ATLAS".
        pattern: Glob pattern for BED files.

    Returns:
        Concatenated unified-schema DataFrame.
    """
    source_dir = Path(source_dir)
    if source == "ENCODE":
        loader = load_encode_bed
    elif source == "CHIP_ATLAS":
        loader = load_chip_atlas_bed
    else:
        raise ValueError(f"Unknown source: {source}")

    bed_files = sorted(source_dir.glob(pattern))
    if not bed_files:
        print(f"  [WARN] No BED files found in {source_dir}")
        return pd.DataFrame(columns=UNIFIED_COLUMNS)

    parts = []
    for bed in bed_files:
        try:
            parts.append(loader(bed))
            print(f"  Loaded {source}: {bed.name} ({len(parts[-1]):,} peaks)")
        except Exception as e:
            print(f"  [ERROR] Failed to load {bed}: {e}")

    return pd.concat(parts, ignore_index=True) if parts else \
        pd.DataFrame(columns=UNIFIED_COLUMNS)


# =============================================================================
# UNIFIED LOADER
# =============================================================================

def load_unified_chip(
    working_dir: Path,
    encode_subdir: str = "ENCODE",
    chip_atlas_subdir: str = "CHIP_ATLAS",
    output_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Load and concatenate ChIP-seq peaks from both ENCODE and ChIP-Atlas.

    ChIP-Atlas rows are filtered to **breast / mammary–concordant** cell contexts
    (see ``BiosampleConfig.chip_brca_celltypes`` and ``chip_atlas_mammary_tissue_keywords``)
    and receive ``cell_subtype`` labels from ``chip_cell_line_subtype_map``.

    Both sources are assumed HG38. Duplicates (same TF + cell type appearing
    in both sources) are kept and distinguished by the `source` column.

    Args:
        working_dir: Parent directory containing the two source subdirectories.
        encode_subdir: Name of ENCODE subdirectory.
        chip_atlas_subdir: Name of ChIP-Atlas subdirectory.
        output_path: If provided, save unified table to this path
                     (parquet if .parquet extension, else CSV).

    Returns:
        Unified DataFrame with columns in ``UNIFIED_COLUMNS``.
    """
    working_dir = Path(working_dir)

    print(f"Loading ENCODE peaks from {working_dir / encode_subdir} ...")
    encode = load_source_directory(working_dir / encode_subdir, source="ENCODE")

    print(f"Loading ChIP-Atlas peaks from {working_dir / chip_atlas_subdir} ...")
    chip_atlas = load_source_directory(
        working_dir / chip_atlas_subdir, source="CHIP_ATLAS"
    )
    chip_atlas = apply_chip_atlas_cell_line_policy(chip_atlas, BIOSAMPLES)

    unified = pd.concat([encode, chip_atlas], ignore_index=True)

    # Drop rows with bad coordinates
    n_pre = len(unified)
    unified = unified.dropna(subset=["chrom", "start", "end"])
    n_dropped = n_pre - len(unified)
    if n_dropped:
        print(f"  Dropped {n_dropped:,} rows with missing coordinates")

    # Normalize TF / cell_type strings (strip whitespace, keep case)
    for col in ("tf", "cell_type"):
        unified[col] = unified[col].astype(str).str.strip()

    print(f"\nUnified ChIP table: {len(unified):,} peaks")
    print(f"  By source:\n{unified['source'].value_counts().to_string()}")
    print(f"  Unique TFs: {unified['tf'].nunique()}")
    print(f"  Unique cell types: {unified['cell_type'].nunique()}")

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix == ".parquet":
            unified.to_parquet(output_path, index=False)
        else:
            unified.to_csv(output_path, index=False)
        print(f"  Saved unified table to: {output_path}")

    return unified
