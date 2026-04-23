#!/usr/bin/env python3
"""
Build a full multi-feature Parquet from the GTF-derived CSV.

Goal:
- Avoid runtime dependence on the ~GB-scale CSV (`data/gencode.v49.annotation.gtf.csv`)
- Preserve exon-level metadata (`exon_id`, `exon_number`, `protein_id`, etc.)
- Enable fast predicate-pushdown reads for panel subsets (SV mapping, gene tables).

Example:

  python scripts/gencode/build_gencode_annotation_parquet.py \\
    --from-csv data/gencode.v49.annotation.gtf.csv \\
    --out data/gencode.v49.annotation.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def _normalize_chunk_types(df: pd.DataFrame) -> pd.DataFrame:
    """
    Force a stable schema across CSV chunks.

    The upstream CSV has mixed-type columns in practice (e.g. artif_dupl sometimes numeric, sometimes string).
    For a robust parquet, we:
    - keep core coordinates as integers
    - keep is_mane_select as boolean
    - store all other columns as strings (nullable)
    """
    df = df.copy()
    for c in ("start", "end", "level"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

    if "is_mane_select" in df.columns:
        # allow 0/1, True/False, empty
        s = df["is_mane_select"]
        if s.dtype != bool:
            df["is_mane_select"] = s.map(
                lambda x: False
                if x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() in ("", "0", "False", "false", "no", "NO")
                else True
                if str(x).strip() in ("1", "True", "true", "yes", "YES")
                else bool(x)
            ).astype(bool)

    # Everything else: string (keeps exon_id/exon_number stable even if numeric in some chunks)
    for c in df.columns:
        if c in ("start", "end", "level", "is_mane_select"):
            continue
        df[c] = df[c].astype("string")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-csv", type=Path, required=True, help="Input GTF-derived CSV (all features).")
    ap.add_argument("--out", type=Path, required=True, help="Output parquet path.")
    ap.add_argument("--chunksize", type=int, default=400_000, help="CSV chunk size.")
    ap.add_argument(
        "--columns",
        type=str,
        default="",
        help="Comma-separated columns to keep (default: keep all columns).",
    )
    args = ap.parse_args()

    csv_path: Path = args.from_csv
    out: Path = args.out
    if not csv_path.is_file():
        sys.exit(f"Not a file: {csv_path}")
    out.parent.mkdir(parents=True, exist_ok=True)

    cols: Optional[list[str]] = None
    if args.columns.strip():
        cols = [c.strip() for c in args.columns.split(",") if c.strip()]

    writer: Optional[pq.ParquetWriter] = None
    schema: Optional[pa.Schema] = None
    for chunk in pd.read_csv(csv_path, chunksize=int(args.chunksize), low_memory=False):
        if cols is not None:
            missing = [c for c in cols if c not in chunk.columns]
            if missing:
                raise ValueError(f"Missing requested columns in CSV: {missing[:10]}")
            chunk = chunk[cols].copy()

        chunk = _normalize_chunk_types(chunk)
        table = pa.Table.from_pandas(chunk, preserve_index=False)
        if writer is None:
            schema = table.schema
            writer = pq.ParquetWriter(out, schema, compression="zstd")
        elif schema is not None:
            # ensure consistent schema (pyarrow enforces this, but error messages are clearer this way)
            if table.schema != schema:
                table = table.cast(schema)
        writer.write_table(table)

    if writer is not None:
        writer.close()

    print(f"Wrote parquet: {out}")


if __name__ == "__main__":
    main()

