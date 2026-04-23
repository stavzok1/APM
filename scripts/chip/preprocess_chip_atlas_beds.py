#!/usr/bin/env python3
"""
Convert ChIP-Atlas UCSC / gffTags BED files under ``CHIP_ATLAS/`` to narrow 6-column BED.

Writes ``chrom``, ``start``, ``end``, ``experiment`` (from ``ID=``), ``tf`` (filename stem),
and ``cell_type`` (from gffTags / ``Name=``), with a header row — the format expected by
``pipeline.CHIP.chip_loader.load_chip_atlas_bed`` after preprocessing.

Examples::

    .venv/bin/python3 scripts/chip/preprocess_chip_atlas_beds.py --dry-run
    .venv/bin/python3 scripts/chip/preprocess_chip_atlas_beds.py --chip-atlas-dir data/CHIP/CHIP_ATLAS
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.CHIP.chip_atlas_convert import (  # noqa: E402
    convert_chip_atlas_file_to_narrow,
    needs_ucsc_conversion,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--chip-atlas-dir",
        type=Path,
        default=ROOT / "data" / "CHIP" / "CHIP_ATLAS",
    )
    ap.add_argument(
        "--backup-suffix",
        type=str,
        default=".pre_narrow.bak",
        help="When replacing in place, copy original to name+suffix first.",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    d = Path(args.chip_atlas_dir)
    if not d.is_dir():
        raise SystemExit(f"Not a directory: {d}")

    beds = sorted(d.glob("*.bed"))
    if not beds:
        print(f"No *.bed under {d}")
        return

    converted = 0
    skipped = 0
    for bed in beds:
        if not needs_ucsc_conversion(bed):
            print(f"  skip (already narrow): {bed.name}")
            skipped += 1
            continue
        n = 0
        tmp = bed.with_suffix(bed.suffix + ".narrow_tmp")
        if args.dry_run:
            print(f"  would convert: {bed.name}")
            converted += 1
            continue
        n = convert_chip_atlas_file_to_narrow(bed, tmp)
        bak = bed.with_name(bed.name + args.backup_suffix)
        shutil.copy2(bed, bak)
        shutil.move(str(tmp), str(bed))
        print(f"  converted {bed.name}: {n:,} rows (backup {bak.name})")
        converted += 1

    print(f"\nDone. converted={converted} skipped={skipped} dry_run={args.dry_run}")


if __name__ == "__main__":
    main()
