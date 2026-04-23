#!/usr/bin/env python3
"""
Build (or rebuild) the cached unified ChIP peak table at `PATHS.chip_unified`.

This reads all BED files under `PATHS.chip_dir/ENCODE` and `PATHS.chip_dir/CHIP_ATLAS`
and writes a single parquet with a `sample_id` column (derived from each BED filename stem),
so multiple TF×cell-type samples remain distinguishable (e.g. CTCF_MCF7_1, CTCF_MCF7_2).

**ChIP-Atlas UCSC exports:** run ``scripts/chip/preprocess_chip_atlas_beds.py`` first so atlas
BEDs are narrow 6-column tables (see ``pipeline.CHIP.chip_atlas_convert``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd  # noqa: E402

from pipeline.config import PATHS  # noqa: E402
from pipeline.CHIP.chip_loader import load_unified_chip  # noqa: E402


def _print_parquet_summary(path: Path) -> None:
    if not path.is_file():
        print(f"Missing parquet: {path}")
        return
    df = pd.read_parquet(path)
    print("path:", path)
    print("rows:", f"{len(df):,}")
    print("columns:", list(df.columns))
    print("\nby source:")
    print(df["source"].value_counts().to_string())
    print("\nunique tf:", df["tf"].nunique())
    print(df["tf"].value_counts().head(15).to_string())
    print("\nunique cell_type:", df["cell_type"].nunique())
    n_empty_ct = int((df["cell_type"].astype(str).str.strip() == "").sum())
    print("empty cell_type rows:", n_empty_ct)
    print(df["cell_type"].value_counts().to_string())
    if "sample_id" in df.columns:
        print("\nunique sample_id:", df["sample_id"].nunique())
        sids = sorted({str(x) for x in df["sample_id"].dropna().astype(str).unique()})
        print("sample_id values:", sids)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build unified ChIP peaks parquet (ENCODE + ChIP-Atlas).")
    ap.add_argument("--chip-dir", type=str, default=str(PATHS.chip_dir), help="CHIP root containing ENCODE/ and CHIP_ATLAS/")
    ap.add_argument("--out", type=str, default=str(PATHS.chip_unified), help="Output parquet path")
    ap.add_argument(
        "--inspect-only",
        action="store_true",
        help="Only print a summary of the existing --out parquet (no rebuild).",
    )
    args = ap.parse_args()

    out = Path(args.out)
    if args.inspect_only:
        _print_parquet_summary(out)
        return

    chip_dir = Path(args.chip_dir)

    df = load_unified_chip(working_dir=chip_dir, output_path=out)
    print(f"\nDone. Unified peaks: {len(df):,}")
    if "sample_id" in df.columns:
        print(f"Unique sample_id: {df['sample_id'].nunique()}")
    _print_parquet_summary(out)


if __name__ == "__main__":
    main()

