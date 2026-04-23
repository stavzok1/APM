#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    from pipeline.config import PATHS
    from pipeline.regulatory_elements import load_regulatory_element_focus

    cols = ["cCRE_id", "gene_links", "atac_peak_links"]
    df = load_regulatory_element_focus(PATHS.regulatory_elements_table_with_evidence_parquet, columns=cols, decode_nested=True)
    if df.empty:
        print("empty")
        return
    # Find first row with some ATAC peak links (if any).
    row = None
    for _, r in df.iterrows():
        apl = r.get("atac_peak_links")
        if isinstance(apl, (list, dict)) and len(apl) > 0:
            row = r.to_dict()
            break
    if row is None:
        row = df.iloc[0].to_dict()
    # Print shallow keys and a short JSON preview for nested fields
    print("cCRE_id:", row.get("cCRE_id"))
    gl = row.get("gene_links")
    apl = row.get("atac_peak_links")
    print("gene_links type:", type(gl), "n_genes:", (len(gl) if isinstance(gl, dict) else None))
    print("atac_peak_links type:", type(apl))
    if isinstance(apl, dict):
        print("atac_peak_links keys (first 20):", list(apl.keys())[:20])
        # show one entry
        if apl:
            k0 = next(iter(apl.keys()))
            v0 = apl[k0]
            print("first key:", k0)
            try:
                print("first val json:", json.dumps(v0, default=str)[:800])
            except Exception:
                print("first val repr:", repr(v0)[:800])
    else:
        print("atac_peak_links repr:", repr(apl)[:1200])


if __name__ == "__main__":
    main()

