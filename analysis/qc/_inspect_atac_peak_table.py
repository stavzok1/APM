#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    import pyarrow.parquet as pq

    p = Path("data/atac_peaks/atac_peaks_annotated.parquet")
    if not p.exists():
        raise SystemExit(f"missing: {p}")

    pf = pq.ParquetFile(p)
    cols = pf.schema.names
    print("n_cols:", len(cols))
    print("first_cols:", cols[:40])

    # Load a tiny slice with nested fields decoded.
    from pipeline.atac_peaks import load_atac_peaks_annotated

    want = [c for c in ("peak_id", "chrom", "start", "end", "gene_links", "genes_by_tier", "linked_genes") if c in cols]
    df = load_atac_peaks_annotated(p, columns=want)
    print("loaded:", df.shape)
    if df.empty:
        return
    r0 = df.iloc[0].to_dict()
    print("row0 keys:", list(r0.keys()))
    for k in ("peak_id", "chrom", "start", "end"):
        if k in r0:
            print(k, "=", r0.get(k))
    gl = r0.get("gene_links")
    print("gene_links type:", type(gl))
    if isinstance(gl, dict):
        print("gene_links n_genes:", len(gl))
        if gl:
            g0 = next(iter(gl.keys()))
            print("first gene:", g0)
            print("first gene val keys:", list(gl[g0].keys()) if isinstance(gl[g0], dict) else type(gl[g0]))


if __name__ == "__main__":
    main()

