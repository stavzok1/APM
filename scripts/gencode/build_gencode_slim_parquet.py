#!/usr/bin/env python3
"""
Build a **gene-row-only** Parquet from GENCODE (``feature == gene``).

**Naming note:** in this repo, ``data/gencode.v49.slim.parquet`` has historically been a
**full** GENCODE extract in Parquet (~7.7M rows: exon, CDS, transcript, … plus **~79k**
``gene`` rows). ``load_genes()`` therefore loads the whole table unless you replace it
with a **true** gene-only file produced by this script. After replacement,
``load_genes()`` is lightweight and ``load_gencode_multifeature_subset()`` still pulls
exons/transcripts from ``gencode.v49.annotation.gtf.csv`` when needed.

The regulatory pipeline uses ``PATHS.gencode_gtf_pq`` for ``load_genes()`` joins; the
full multi-feature table stays in ``gencode.v49.annotation.gtf.csv`` (or your
``*.annotation.gtf.gz``) for chunked multifeature loads.

Inputs (pick one):

1. **CSV** (recommended if you already converted GTF → CSV, as in ``data/gencode.v49.annotation.gtf.csv``):
   chunked read, filter ``feature == "gene"``, write Parquet.

2. **GTF.GZ** (streaming, low memory): scan lines, keep rows where column 3 is ``gene``,
   parse ``gene_id``, ``gene_name``, ``gene_type`` from the attributes column.

Examples::

  python scripts/gencode/build_gencode_slim_parquet.py \\
    --from-csv data/gencode.v49.annotation.gtf.csv \\
    --out data/gencode.v49.slim.parquet

  python scripts/gencode/build_gencode_slim_parquet.py \\
    --from-gtf-gz /path/to/gencode.v49.annotation.gtf.gz \\
    --out data/gencode.v49.slim.parquet

After rebuilding, update ``pipeline/config.py`` ``PathConfig.gencode_gtf_pq`` if you change
the filename, and clear any downstream caches that assume frozen coordinates.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional

import pandas as pd

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


def _attr_get(attrs: str, key: str) -> str:
    m = re.search(rf'{re.escape(key)}\s+"([^"]*)"', attrs)
    return m.group(1) if m else ""


def _iter_gtf_genes_gz(path: Path) -> Iterator[Dict[str, object]]:
    """Yield dicts for GTF lines with feature == gene."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            yield {
                "seqname": parts[0],
                "source": parts[1],
                "feature": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5],
                "strand": parts[6],
                "frame": parts[7],
                "gene_id": _attr_get(attrs, "gene_id"),
                "gene_name": _attr_get(attrs, "gene_name"),
                "gene_type": _attr_get(attrs, "gene_type"),
                "transcript_id": "",
                "transcript_name": "",
                "transcript_type": "",
                "tag_list": "",
                "tag": _attr_get(attrs, "tag"),
                "is_MANE": False,
            }


def _slim_from_csv(csv_path: Path, chunksize: int) -> pd.DataFrame:
    parts: List[pd.DataFrame] = []
    for chunk in pd.read_csv(csv_path, chunksize=chunksize, low_memory=False):
        if "feature" not in chunk.columns:
            raise ValueError(f"{csv_path}: expected a 'feature' column")
        g = chunk.loc[chunk["feature"].astype(str) == "gene"].copy()
        if not g.empty:
            parts.append(g)
    if not parts:
        return pd.DataFrame()
    out = pd.concat(parts, ignore_index=True)
    return _harmonize_slim_columns(out)


def _harmonize_slim_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Match columns expected by ``load_genes`` / existing slim parquet."""
    df = df.copy()
    if "tag_list" not in df.columns and "tag" in df.columns:
        df["tag_list"] = df["tag"].astype(str).replace("nan", "")
    if "tag_list" not in df.columns:
        df["tag_list"] = ""
    if "tag" not in df.columns:
        df["tag"] = ""
    for c in ("transcript_id", "transcript_name", "transcript_type"):
        if c not in df.columns:
            df[c] = pd.NA
    if "is_mane_select" in df.columns and "is_MANE" not in df.columns:
        df["is_MANE"] = df["is_mane_select"].fillna(False).astype(bool)
    elif "is_MANE" not in df.columns:
        df["is_MANE"] = False
    keep = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "gene_id",
        "gene_name",
        "gene_type",
        "transcript_id",
        "transcript_name",
        "transcript_type",
        "tag_list",
        "tag",
        "is_MANE",
    ]
    for c in keep:
        if c not in df.columns:
            df[c] = pd.NA
    return df[keep]


def _slim_from_gtf_gz(path: Path, batch: int) -> pd.DataFrame:
    batch_rows: List[Dict[str, object]] = []
    parts: List[pd.DataFrame] = []
    for row in _iter_gtf_genes_gz(path):
        batch_rows.append(row)
        if len(batch_rows) >= batch:
            parts.append(pd.DataFrame(batch_rows))
            batch_rows = []
    if batch_rows:
        parts.append(pd.DataFrame(batch_rows))
    if not parts:
        return pd.DataFrame()
    out = pd.concat(parts, ignore_index=True)
    return _harmonize_slim_columns(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-csv", type=Path, help="GENCODE table as CSV (any row layout with a `feature` column).")
    ap.add_argument("--from-gtf-gz", type=Path, help="GENCODE GTF or GTF.GZ (streams `gene` rows only).")
    ap.add_argument("--out", type=Path, required=True, help="Output Parquet path.")
    ap.add_argument("--chunksize", type=int, default=400_000, help="CSV chunk size.")
    ap.add_argument("--gtf-batch", type=int, default=50_000, help="GTF rows per concat batch.")
    args = ap.parse_args()

    if bool(args.from_csv) == bool(args.from_gtf_gz):
        ap.error("Specify exactly one of --from-csv or --from-gtf-gz")

    if args.from_csv:
        if not args.from_csv.is_file():
            sys.exit(f"Not a file: {args.from_csv}")
        slim = _slim_from_csv(args.from_csv, args.chunksize)
    else:
        if not args.from_gtf_gz.is_file():
            sys.exit(f"Not a file: {args.from_gtf_gz}")
        slim = _slim_from_gtf_gz(args.from_gtf_gz, args.gtf_batch)

    if slim.empty:
        sys.exit("No gene rows found — check input format.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    slim.to_parquet(args.out, index=False)
    print(f"Wrote {args.out}  ({len(slim)} gene rows)")


if __name__ == "__main__":
    main()
