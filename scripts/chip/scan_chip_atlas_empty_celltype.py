#!/usr/bin/env python3
"""One-off: per BED file, count rows where ChIP-Atlas cell_type parsing returned empty."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pipeline.CHIP.chip_loader import load_chip_atlas_bed  # noqa: E402


def main() -> None:
    d = ROOT / "data" / "CHIP" / "CHIP_ATLAS"
    rows = []
    for bed in sorted(d.glob("*.bed")):
        df = load_chip_atlas_bed(bed)
        n = int((df["cell_type"].astype(str).str.strip() == "").sum())
        if n:
            rows.append((n, bed.name, len(df)))
    rows.sort(reverse=True)
    for n, name, tot in rows:
        print(f"{name}\t{n}\tof\t{tot}\t({100 * n / tot:.2f}%)")


if __name__ == "__main__":
    main()
