from __future__ import annotations

import argparse

from pipeline.lncRNA_interactions.postar3_summary import (
    Postar3SummaryPaths,
    build_postar3_overlaps_for_selected_lncrnas,
)


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize POSTAR3 RBP peaks overlapping selected lncRNA loci.")
    ap.add_argument("--n-extra-close-lncrnas", type=int, default=20)
    ap.add_argument("--chunk-rows", type=int, default=1_500_000)
    ap.add_argument("--region-mode", type=str, default="gene", choices=["gene", "exons", "introns", "promoter"])
    ap.add_argument("--promoter-upstream-bp", type=int, default=2000)
    ap.add_argument("--promoter-downstream-bp", type=int, default=500)
    args = ap.parse_args()

    paths = Postar3SummaryPaths()
    overlaps, rbp_summary, lncrna_summary = build_postar3_overlaps_for_selected_lncrnas(
        n_extra_close_lncrnas=int(args.n_extra_close_lncrnas),
        paths=paths,
        chunk_rows=int(args.chunk_rows),
        region_mode=str(args.region_mode),
        promoter_upstream_bp=int(args.promoter_upstream_bp),
        promoter_downstream_bp=int(args.promoter_downstream_bp),
    )
    print(f"overlaps: {paths.out_overlap_parquet} ({overlaps.shape[0]} rows)")
    print(f"rbp_summary: {paths.out_rbp_summary_parquet} ({rbp_summary.shape[0]} rows)")
    print(f"lncrna_summary: {paths.out_lncrna_summary_parquet} ({lncrna_summary.shape[0]} rows)")


if __name__ == "__main__":
    main()

