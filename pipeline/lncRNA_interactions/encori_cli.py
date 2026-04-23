from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pipeline.config import PATHS
from pipeline.lncRNA_interactions.encori import (
    EncoriQuery,
    build_encori_lncrna_target_list,
    fetch_encori_miRNA_targets,
    fetch_encori_rbp_targets,
)
from pipeline.config import PRIMARY_GENES, TIER1_LNCRNA_GENES


def main() -> None:
    ap = argparse.ArgumentParser(description="Download ENCORI interactions (cached TSVs).")
    ap.add_argument("--cache-dir", type=Path, default=Path(PATHS.working_dir) / "external_cache")
    ap.add_argument("--assembly", type=str, default="hg38")
    ap.add_argument("--clip-exp-num", type=int, default=2)
    ap.add_argument("--pancancer-num", type=int, default=0)
    ap.add_argument("--cell-type", type=str, default="all")
    ap.add_argument("--gene-type", type=str, default="lncRNA")

    ap.add_argument(
        "--lncrna-proximity-pairs",
        type=Path,
        default=Path(PATHS.working_dir) / "lncRNA_matching" / "genes_lncRNAs_1000000bp_distances.csv",
        help="Proximity pairs (primary gene ↔ lncRNA) to pick extra close lncRNAs.",
    )
    ap.add_argument("--n-extra-close-lncrnas", type=int, default=20)

    ap.add_argument("--rbp-list", type=str, default="", help="Comma-separated RBPs (recommended).")
    ap.add_argument("--out-rbp", type=Path, default=Path(PATHS.working_dir) / "lncRNA_interactions" / "encori_rbp_targets.parquet")

    ap.add_argument("--do-mirna", action="store_true", help="Also fetch miRNA targets for modelled lncRNAs.")
    ap.add_argument("--out-mirna", type=Path, default=Path(PATHS.working_dir) / "lncRNA_interactions" / "encori_mirna_targets.parquet")

    args = ap.parse_args()

    q = EncoriQuery(
        assembly=args.assembly,
        gene_type=args.gene_type,
        cell_type=args.cell_type,
        clip_exp_num=int(args.clip_exp_num),
        pancancer_num=int(args.pancancer_num),
    )

    if args.rbp_list.strip():
        rbps = [x.strip() for x in args.rbp_list.split(",") if x.strip()]
        df_rbp = fetch_encori_rbp_targets(rbps=rbps, query=q, cache_dir=args.cache_dir)
        # Preserve a stable schema even if ENCORI returns no valid rows.
        if df_rbp.empty:
            df_rbp = pd.DataFrame(
                columns=[
                    "RBP",
                    "geneID",
                    "geneName",
                    "geneType",
                    "clusterNum",
                    "totalClipExpNum",
                    "totalClipSiteNum",
                    "clusterID",
                    "chromosome",
                    "narrowStart",
                    "narrowEnd",
                    "broadStart",
                    "broadEnd",
                    "strand",
                    "clipExpNum",
                    "HepG2(shRNA)",
                    "K562(shRNA)",
                    "HepG2(CRISPR)",
                    "K562(CRISPR)",
                    "pancancerNum",
                    "cellline/tissue",
                    "__encori_url",
                ]
            )
        args.out_rbp.parent.mkdir(parents=True, exist_ok=True)
        df_rbp.to_parquet(args.out_rbp, index=False)
        print(f"Wrote: {args.out_rbp} ({df_rbp.shape[0]} rows)")

    if args.do_mirna:
        lncrnas = build_encori_lncrna_target_list(
            panel_lncrnas=TIER1_LNCRNA_GENES,
            lncrna_proximity_pairs_csv=args.lncrna_proximity_pairs,
            primary_genes=PRIMARY_GENES,
            n_extra_close_lncrnas=int(args.n_extra_close_lncrnas),
        )
        df_m = fetch_encori_miRNA_targets(targets=lncrnas, query=q, cache_dir=args.cache_dir)
        args.out_mirna.parent.mkdir(parents=True, exist_ok=True)
        df_m.to_parquet(args.out_mirna, index=False)
        print(f"Wrote: {args.out_mirna} ({df_m.shape[0]} rows; targets={len(lncrnas)})")


if __name__ == "__main__":
    main()

