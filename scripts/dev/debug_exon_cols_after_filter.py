#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd
import pyarrow.dataset as ds

from pipeline.config import PATHS, PRIMARY_GENES, TIER1_LNCRNA_GENES
from pipeline.lncRNA_interactions.encori import build_encori_lncrna_target_list


def main() -> None:
    p = str(PATHS.gencode_gtf_full_pq)
    d = ds.dataset(p, format="parquet")
    cols = set(d.schema.names)
    chrom_col = "chrom" if "chrom" in cols else "seqname"

    base = [chrom_col, "start", "end", "strand", "gene_id", "feature"]
    if "gene_name" in cols:
        base.insert(4, "gene_name")
    elif "gene" in cols:
        base.insert(4, "gene")
    else:
        raise RuntimeError("no gene symbol column")

    sel = build_encori_lncrna_target_list(
        panel_lncrnas=TIER1_LNCRNA_GENES,
        lncrna_proximity_pairs_csv=PATHS.working_dir / "lncRNA_matching" / "genes_lncRNAs_1000000bp_distances.csv",
        primary_genes=PRIMARY_GENES,
        n_extra_close_lncrnas=20,
    )
    wanted = set(sel)

    ex = d.to_table(columns=base, filter=(ds.field("feature") == "exon")).to_pandas()
    if "gene_name" in ex.columns:
        sym = ex["gene_name"].astype(str)
    else:
        sym = ex["gene"].astype(str)
    ex = ex[sym.isin(wanted)].copy()
    if "gene_name" not in ex.columns and "gene" in ex.columns:
        ex["gene_name"] = ex["gene"].astype(str)

    print("cols:", ex.columns.tolist())
    print(ex.head(3).to_string(index=False))


if __name__ == "__main__":
    main()
