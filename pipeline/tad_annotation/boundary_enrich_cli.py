from __future__ import annotations

import argparse
from pathlib import Path

from pipeline.config import PATHS

from .boundary_enrichment import build_boundaries_enriched_table


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a single TAD boundary enriched table (cCRE/ATAC/probes + adjacent-domain genes).")
    ap.add_argument("--processed-dir", type=str, default="", help="TAD processed dir (default: PATHS.tads_processed).")
    ap.add_argument("--out", type=str, default="", help="Output parquet path.")
    ap.add_argument("--biosamples", type=str, default="", help="Comma-separated biosamples (default: all discovered).")
    args = ap.parse_args()

    processed_dir = Path(args.processed_dir) if args.processed_dir else Path(PATHS.tads_processed)
    out_path = Path(args.out) if args.out else (processed_dir / "boundaries_enriched.parquet")

    bios = None
    if args.biosamples.strip():
        bios = [x.strip() for x in args.biosamples.split(",") if x.strip()]

    df = build_boundaries_enriched_table(processed_dir=processed_dir, biosamples=bios)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    print(f"Wrote: {out_path}  (rows={len(df)}, cols={df.shape[1]})")


if __name__ == "__main__":
    main()

