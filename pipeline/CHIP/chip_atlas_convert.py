"""
Convert ChIP-Atlas **UCSC / gffTags** BED exports to **narrow 6-column** BED.

Narrow format (tab-separated, optional header row)::

    chrom  start  end  experiment  tf  cell_type

``experiment`` is taken from ``ID=…`` in the gffTags column; ``tf`` from the
filename stem (same convention as the unified ChIP table); ``cell_type`` from
gffTags keys / ``Name=TF (@ cell)`` (URL-decoded), matching the logic previously
embedded in ``chip_loader.load_chip_atlas_bed``.

Run ``scripts/chip/preprocess_chip_atlas_beds.py`` on ``data/CHIP/CHIP_ATLAS`` before
``scripts/chip/build_unified_chip_peaks.py`` when atlas downloads are UCSC-style.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List, Optional, TextIO
from urllib.parse import unquote

from ..utils import harmonize_chrom_column


def looks_like_chromosome_token(chrom: str) -> bool:
    s = str(chrom).strip()
    if not s or s.lower() in ("chrom", "#chrom"):
        return False
    return s.startswith("chr")


def chip_atlas_gfftags_cell_line(meta: str) -> str:
    """Cell / biosample label from ChIP-Atlas ``gffTags`` (BED column 4)."""
    if not meta or not isinstance(meta, str):
        return ""
    m = meta.strip()
    for pat in (
        r"cell%20line=([^;]+)",
        r"(?:^|;)cell=([^;]+)",
        r"ArrayExpress-CellType=([^;]+)",
        r"<br>source_name=([^;]+)",
        r"cell%20type=([^;]+)",
        r"(?:^|;)source_name=([^;]+)",
        r"cell_line=([^;]+)",
        r"cell_type=([^;]+)",
    ):
        mm = re.search(pat, m, flags=re.IGNORECASE)
        if mm:
            return unquote(mm.group(1).strip())
    mm = re.search(r"Name=([^;]+)", m)
    if mm:
        name = unquote(mm.group(1).strip())
        if " (@ " in name:
            return name.split(" (@ ", 1)[1].strip().rstrip(")")
    return ""


def extract_experiment_from_gfftags(meta: str) -> str:
    mm = re.search(r"(?:^|;)ID=([^;]+)", str(meta).strip())
    if mm:
        return unquote(mm.group(1).strip())
    mm = re.search(r"^ID=([^;]+)", str(meta).strip())
    if mm:
        return unquote(mm.group(1).strip())
    return ""


def row_is_ucsc_gfftags_bed(parts: List[str]) -> bool:
    if len(parts) < 6:
        return False
    meta = parts[3]
    if "Name=" in meta:
        return True
    if meta.startswith("ID="):
        return True
    return False


def gfftags_row_to_narrow(
    parts: List[str],
    tf_from_filename: str,
) -> Optional[tuple[str, str, str, str, str, str]]:
    """Return (chrom, start, end, experiment, tf, cell_type) or None if invalid."""
    chrom = parts[0].strip()
    if not looks_like_chromosome_token(chrom):
        return None
    meta = parts[3]
    exp = extract_experiment_from_gfftags(meta)
    cell = chip_atlas_gfftags_cell_line(meta)
    return (chrom, parts[1].strip(), parts[2].strip(), exp, tf_from_filename, cell)


def iter_narrow_rows_from_atlas_bed(
    path: Path,
    tf_from_filename: Optional[str] = None,
) -> Iterable[tuple[str, str, str, str, str, str]]:
    """
    Yield narrow 6-tuples from a ChIP-Atlas BED (either pre-narrow or UCSC/gffTags).

    Skips ``track`` / ``browser`` / ``#`` lines and header rows.
    """
    path = Path(path)
    tf_stem = tf_from_filename if tf_from_filename is not None else path.stem
    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line or line.startswith(("track ", "browser ", "#")):
                continue
            parts = [p.strip() for p in line.split("\t")]
            while parts and parts[-1] == "":
                parts.pop()
            if len(parts) < 6:
                continue
            chrom = parts[0]
            if not looks_like_chromosome_token(chrom):
                continue
            if str(chrom).lower() == "chrom" and str(parts[1]).lower() == "start":
                continue
            if row_is_ucsc_gfftags_bed(parts):
                row = gfftags_row_to_narrow(parts, tf_stem)
                if row is not None:
                    yield row
            else:
                # Legacy narrow: chrom start end experiment tf cell_type
                if len(parts) != 6:
                    continue
                yield (
                    chrom,
                    parts[1],
                    parts[2],
                    parts[3],
                    parts[4],
                    parts[5],
                )


def write_narrow_chip_atlas_bed(
    in_path: Path,
    out: TextIO,
    *,
    tf_from_filename: Optional[str] = None,
    write_header: bool = True,
) -> int:
    """Write narrow 6-column BED to open text stream. Returns number of data rows."""
    n = 0
    if write_header:
        out.write("chrom\tstart\tend\texperiment\ttf\tcell_type\n")
    for row in iter_narrow_rows_from_atlas_bed(in_path, tf_from_filename=tf_from_filename):
        out.write("\t".join(row) + "\n")
        n += 1
    return n


def needs_ucsc_conversion(path: Path, max_scan_lines: int = 12_000) -> bool:
    """True if file appears to contain at least one UCSC/gffTags row."""
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for i, raw in enumerate(fh):
            if i > max_scan_lines:
                break
            line = raw.rstrip("\r\n")
            if not line or line.startswith(("track ", "browser ", "#")):
                continue
            parts = [p.strip() for p in line.split("\t")]
            while parts and parts[-1] == "":
                parts.pop()
            if len(parts) < 6:
                continue
            if looks_like_chromosome_token(parts[0]) and row_is_ucsc_gfftags_bed(parts):
                return True
    return False


def convert_chip_atlas_file_to_narrow(
    in_path: Path,
    out_path: Path,
    *,
    tf_from_filename: Optional[str] = None,
) -> int:
    """Convert in_path to narrow BED at out_path. Returns row count."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        return write_narrow_chip_atlas_bed(
            in_path, out, tf_from_filename=tf_from_filename, write_header=True
        )
