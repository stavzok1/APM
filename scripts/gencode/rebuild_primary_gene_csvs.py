#!/usr/bin/env python3
"""
Rebuild gene panel CSVs written by ``save_gene_tables`` (including
``primary_genes_all_features.csv``) without running the full regulatory pipeline.

Uses ``PATHS`` output locations (default under ``PATHS.working_dir``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pipeline.config import (
    CNV_GENES,
    PATHS,
    PIPELINE_GENE_PANEL,
    THRESHOLDS,
    use_legacy_lncrna_intervals_csv_input,
)
from pipeline.genes import (
    add_promoter_columns,
    filter_genes_by_names,
    filter_lncrnas,
    load_genes,
    load_lncrnas,
    match_lncrnas_to_genes,
)
from pipeline.genes.gene_loader import (
    harmonize_multiple_dfs,
    lncrna_gene_intervals_from_annotation,
    save_gene_tables,
)
from pipeline.genes.lncrna_matching import get_matched_lncrna_names


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-panel",
        type=Path,
        default=None,
        help="Optional path to a newline-separated gene symbol list (default: PIPELINE_GENE_PANEL)",
    )
    args = parser.parse_args()

    if args.gene_panel is not None:
        text = args.gene_panel.read_text(encoding="utf-8")
        gene_panel = [ln.strip() for ln in text.splitlines() if ln.strip()]
    else:
        gene_panel = list(PIPELINE_GENE_PANEL)

    genes_all = load_genes(PATHS.gencode_gtf_pq)
    if use_legacy_lncrna_intervals_csv_input():
        lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
    else:
        lncrnas_all = lncrna_gene_intervals_from_annotation(genes_all)
    [genes_all, lncrnas_all], _ = harmonize_multiple_dfs([genes_all, lncrnas_all])
    genes_all_harmonized = genes_all

    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(
        genes,
        THRESHOLDS.promoter_upstream_bp,
        THRESHOLDS.promoter_downstream_bp,
    )
    cnv_genes = filter_genes_by_names(genes_all, CNV_GENES)
    cnv_genes = add_promoter_columns(
        cnv_genes,
        THRESHOLDS.promoter_upstream_bp,
        THRESHOLDS.promoter_downstream_bp,
    )
    lncrnas = filter_lncrnas(lncrnas_all)
    lncrnas = add_promoter_columns(
        lncrnas,
        THRESHOLDS.promoter_upstream_bp,
        THRESHOLDS.promoter_downstream_bp,
    )

    _pairs_df, _genes_with_lnc, lncrnas_with_genes = match_lncrnas_to_genes(
        genes, lncrnas, THRESHOLDS.lncrna_window_bp
    )
    matched_lncrna_names = get_matched_lncrna_names(lncrnas_with_genes)

    out_dir = Path(PATHS.working_dir)
    save_gene_tables(
        genes,
        cnv_genes,
        gene_panel,
        matched_lncrna_names,
        out_dir,
        genes_all_harmonized=genes_all_harmonized,
    )
    print(f"Wrote: {PATHS.genes_all_features}")
    print(f"Wrote: {PATHS.genes_only}")
    print(f"Wrote: {PATHS.lncrnas_all_features}")
    print(f"Wrote: {PATHS.cnv_genes}")
    print(f"Wrote: {PATHS.tier2_medium_genes_only}")
    print(f"Wrote: {PATHS.tier2_medium_genes_all_features}")
    print(f"Wrote: {PATHS.tier3_cnv_only_genes_only}")
    print(f"Wrote: {PATHS.tier3_cnv_only_genes_all_features}")
    print(f"Wrote: {PATHS.tier4_readout_genes_only}")
    print(f"Wrote: {PATHS.tier4_readout_genes_all_features}")


if __name__ == "__main__":
    main()
