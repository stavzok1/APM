from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from pipeline.config import PATHS, PRIMARY_GENES, TIER1_LNCRNA_GENES
from pipeline.lncRNA_interactions.encori import (
    EncoriQuery,
    build_encori_lncrna_target_list,
    fetch_encori_miRNA_targets_with_symbol_fallback,
    fetch_encori_rbp_targets_for_lncrnas,
    fetch_encori_rbp_targets_for_lncrnas_with_fallback,
)
from pipeline.lncRNA_interactions.postar3_summary import (
    Postar3SummaryPaths,
    build_postar3_overlaps_for_selected_lncrnas,
)


@dataclass(frozen=True)
class LncRnaInteractionsPaths:
    out_dir: Path = Path(PATHS.working_dir) / "lncRNA_interactions"
    encori_mirna_targets_pq: Path = out_dir / "encori_mirna_targets.parquet"
    encori_mirna_targets_diag_csv: Path = out_dir / "encori_mirna_targets_diagnostics.csv"
    encori_rbp_targets_pq: Path = out_dir / "encori_rbp_targets.parquet"
    encori_rbp_targets_diag_csv: Path = out_dir / "encori_rbp_targets_diagnostics.csv"
    postar3_rbp_summary_pq: Path = out_dir / "postar3_rbp_summary.parquet"
    postar3_lncrna_summary_pq: Path = out_dir / "postar3_lncrna_summary.parquet"
    postar3_overlaps_pq: Path = out_dir / "postar3_overlaps.parquet"
    selected_lncrnas_txt: Path = out_dir / "selected_lncrnas.txt"
    recommended_rbps_csv: Path = out_dir / "recommended_rbps_from_postar3.csv"


def build_all_lncRNA_interactions(
    *,
    n_extra_close_lncrnas: int = 20,
    encori_clip_exp_num: int = 2,
    encori_rbp_clip_exp_num: int = 2,
    encori_rbp_top_n: int = 12,
    cache_dir: Path | None = None,
    postar3_chunk_rows: int = 1_500_000,
    postar3_region_mode: str = "gene",
    out_paths: LncRnaInteractionsPaths = LncRnaInteractionsPaths(),
) -> None:
    out_paths.out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = cache_dir or (Path(PATHS.working_dir) / "external_cache")

    selected_lncrnas = build_encori_lncrna_target_list(
        panel_lncrnas=TIER1_LNCRNA_GENES,
        lncrna_proximity_pairs_csv=Path(PATHS.working_dir) / "lncRNA_matching" / "genes_lncRNAs_1000000bp_distances.csv",
        primary_genes=PRIMARY_GENES,
        n_extra_close_lncrnas=int(n_extra_close_lncrnas),
    )
    out_paths.selected_lncrnas_txt.write_text("\n".join(selected_lncrnas) + "\n", encoding="utf-8")

    # 1) ENCORI miRNA→lncRNA (arm-level miRNAs, MiRBase MIMAT ids)
    q = EncoriQuery(clip_exp_num=int(encori_clip_exp_num))
    df_m, df_diag = fetch_encori_miRNA_targets_with_symbol_fallback(
        targets=selected_lncrnas,
        query=q,
        cache_dir=cache_dir,
    )
    df_m.to_parquet(out_paths.encori_mirna_targets_pq, index=False)
    df_diag.to_csv(out_paths.encori_mirna_targets_diag_csv, index=False)

    # 2) POSTAR3 overlap-derived RBP summaries
    ppaths = Postar3SummaryPaths(
        out_overlap_parquet=out_paths.postar3_overlaps_pq,
        out_rbp_summary_parquet=out_paths.postar3_rbp_summary_pq,
        out_lncrna_summary_parquet=out_paths.postar3_lncrna_summary_pq,
    )
    overlaps, rbp_summary, lncrna_summary = build_postar3_overlaps_for_selected_lncrnas(
        n_extra_close_lncrnas=int(n_extra_close_lncrnas),
        paths=ppaths,
        chunk_rows=int(postar3_chunk_rows),
        region_mode=str(postar3_region_mode),
    )

    # 3) Recommended RBPs list (derived from POSTAR3 for your lncRNA set)
    #    Keep a compact table you can use for “follow-up” (ENCORI/eCLIP/interpretation).
    rec = rbp_summary[["rbp", "n_peaks", "n_assays", "n_cell_tissue"]].copy()
    rec.to_csv(out_paths.recommended_rbps_csv, index=False)

    # 4) ENCORI RBPTarget for a small curated RBP set derived from POSTAR3 (practical subset of ENCORI universe)
    top = (
        rbp_summary["rbp"]
        .astype(str)
        .loc[lambda s: s.str.len() > 1]
        .loc[lambda s: ~s.str.contains("parameter", case=False, na=False)]
        .head(int(encori_rbp_top_n))
        .tolist()
    )
    q_rbp = EncoriQuery(
        assembly=q.assembly,
        gene_type="lncRNA",
        cell_type="all",
        clip_exp_num=int(encori_rbp_clip_exp_num),
        pancancer_num=int(q.pancancer_num),
    )
    df_rbp, df_rbp_diag = fetch_encori_rbp_targets_for_lncrnas_with_fallback(
        rbps=top,
        lncrna_targets=selected_lncrnas,
        query=q_rbp,
        cache_dir=cache_dir,
    )
    df_rbp.to_parquet(out_paths.encori_rbp_targets_pq, index=False)
    df_rbp_diag.to_csv(out_paths.encori_rbp_targets_diag_csv, index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build unified lncRNA interaction tables (ENCORI + POSTAR3).")
    ap.add_argument("--n-extra-close-lncrnas", type=int, default=20)
    ap.add_argument("--encori-clip-exp-num", type=int, default=2)
    ap.add_argument("--encori-rbp-clip-exp-num", type=int, default=2)
    ap.add_argument("--encori-rbp-top-n", type=int, default=12)
    ap.add_argument("--postar3-region-mode", type=str, default="gene", choices=["gene", "exons", "introns", "promoter"])
    ap.add_argument("--postar3-chunk-rows", type=int, default=1_500_000)
    ap.add_argument("--cache-dir", type=Path, default=Path(PATHS.working_dir) / "external_cache")
    args = ap.parse_args()

    build_all_lncRNA_interactions(
        n_extra_close_lncrnas=int(args.n_extra_close_lncrnas),
        encori_clip_exp_num=int(args.encori_clip_exp_num),
        encori_rbp_clip_exp_num=int(args.encori_rbp_clip_exp_num),
        encori_rbp_top_n=int(args.encori_rbp_top_n),
        postar3_region_mode=str(args.postar3_region_mode),
        postar3_chunk_rows=int(args.postar3_chunk_rows),
        cache_dir=args.cache_dir,
    )
    print("Wrote unified lncRNA_interactions outputs under data/lncRNA_interactions/")


if __name__ == "__main__":
    main()

