#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


COLS = [
    "chrom",
    "start",
    "end",
    "peak_id",
    "strand",
    "rbp",
    "assay",
    "cell_tissue",
    "source_accession",
    "score",
]


def _iter_chunks(
    path: Path,
    chunksize: int,
) -> Iterable[pd.DataFrame]:
    # POSTAR3.txt appears to be a BED-like 10-column TSV with no header.
    # We stream it to parquet to avoid holding it in memory.
    for chunk in pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=COLS,
        chunksize=int(chunksize),
        low_memory=False,
        dtype={
            "chrom": "string",
            "start": "int64",
            "end": "int64",
            "peak_id": "string",
            "strand": "string",
            "rbp": "string",
            "assay": "string",
            "cell_tissue": "string",
            "source_accession": "string",
            "score": "float64",
        },
    ):
        yield chunk


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert POSTAR3.txt to parquet (streaming).")
    ap.add_argument("--in", dest="inp", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    ap.add_argument("--rbp", type=str, default="", help="Comma-separated RBPs to keep (optional).")
    ap.add_argument("--cell", type=str, default="", help="Substring filter on cell_tissue (optional).")
    args = ap.parse_args()

    rbps: Optional[set[str]] = None
    if args.rbp.strip():
        rbps = {x.strip() for x in args.rbp.split(",") if x.strip()}

    args.out.parent.mkdir(parents=True, exist_ok=True)

    writer: Optional[pq.ParquetWriter] = None
    schema: Optional[pa.Schema] = None
    for chunk in _iter_chunks(args.inp, args.chunksize):
        if rbps is not None:
            chunk = chunk[chunk["rbp"].isin(rbps)]
        if args.cell.strip():
            chunk = chunk[chunk["cell_tissue"].astype(str).str.contains(args.cell.strip(), na=False)]
        if chunk.empty:
            continue
        table = pa.Table.from_pandas(chunk, preserve_index=False)
        if writer is None:
            schema = table.schema
            writer = pq.ParquetWriter(args.out, schema, compression="zstd")
        elif schema is not None and table.schema != schema:
            table = table.cast(schema)
        writer.write_table(table)

    if writer is not None:
        writer.close()

    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()

