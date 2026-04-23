#!/usr/bin/env python3
from __future__ import annotations

import pyarrow.dataset as ds

from pipeline.config import PATHS


def main() -> None:
    p = str(PATHS.gencode_gtf_full_pq)
    d = ds.dataset(p, format="parquet")
    cols = d.schema.names
    print("path:", p)
    print("ncols:", len(cols))
    print("gene-ish columns:", [c for c in cols if "gene" in c.lower()][:120])

    t = d.to_table(columns=cols, filter=(ds.field("feature") == "exon")).slice(0, 1)
    row = t.to_pandas().iloc[0].to_dict()
    print("sample_exon_row_keys:", sorted(row.keys()))


if __name__ == "__main__":
    main()
