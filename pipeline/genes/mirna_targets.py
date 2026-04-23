"""
miRNA target prediction processing from TargetScan.

Loads TargetScan predictions, filters to genes of interest,
and identifies top miRNA regulators per gene.
"""

import os
from pathlib import Path
from typing import List, Optional

import pandas as pd


# =============================================================================
# TARGETSCAN PROCESSING
# =============================================================================

def get_mirna_targets(
    genes: pd.DataFrame,
    gtf: pd.DataFrame,
    targetscan_path: Path,
    top_n: int = 200,
    weight_threshold: float = -0.2,
    output_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Extract miRNA targets for genes of interest from TargetScan predictions.
    """
    print(f"Loading TargetScan predictions from {targetscan_path}...")
    targetscan_df = pd.read_csv(targetscan_path, sep="\t")

    if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    ):
        from pipeline.genes.symbol_normalization import (
            apply_symbol_mapping_series,
            default_symbol_mapping,
        )

        if "Gene Symbol" in targetscan_df.columns:
            targetscan_df["Gene Symbol"] = apply_symbol_mapping_series(
                targetscan_df["Gene Symbol"], default_symbol_mapping()
            )

    # Extract gene IDs without version numbers for matching
    gene_ids_no_version = genes["gene_id"].str.split(".").str[0].unique()

    # Filter to human genes in our panel
    df = targetscan_df[
        (targetscan_df["Gene Tax ID"] == 9606) &
        (targetscan_df["Gene ID"].str.split(".").str[0].isin(gene_ids_no_version))
    ].copy()

    print(f"Found {len(df)} TargetScan rows for {len(gene_ids_no_version)} genes")

    if df.empty:
        print("Warning: No TargetScan predictions found for the gene panel")
        return pd.DataFrame()

    # Keep only rows whose transcript exists in GTF
    df = df[df["Transcript ID"].isin(gtf["transcript_id"])].copy()

    print(f"Found {len(df)} TargetScan rows with transcripts in GTF")

    # --- NEW: build per-site dicts ---
    # (keeping original column name "weighted context++ score")
    df["site_dict"] = df.apply(
        lambda r: {
            "transcript_id": r["Transcript ID"],
            "UTR_start": int(r["UTR_start"]) if pd.notna(r["UTR_start"]) else None,
            "UTR_end": int(r["UTR end"]) if pd.notna(r["UTR end"]) else None,
            "weighted_context_score": float(r["weighted context++ score"])
                if pd.notna(r["weighted context++ score"]) else None,
        },
        axis=1,
    )

    # Aggregate per (gene, miRNA)
    # - best score = min
    # - sites = list of dicts
    gene_mirna_best = (
        df.groupby(["Gene ID", "Gene Symbol", "miRNA"], as_index=False)
        .agg(
            weighted_context_score=("weighted context++ score", "min"),
            sites=("site_dict", list),
        )
    )

    # Rename ID/Symbol columns to your standard names
    gene_mirna_best = gene_mirna_best.rename(columns={
        "Gene ID": "gene_id",
        "Gene Symbol": "gene_symbol",
        "miRNA": "miRNA",
    })

    # Add number of binding sites
    gene_mirna_best["num_sites"] = gene_mirna_best["sites"].apply(len)

    # Get top N miRNAs per gene by score
    top_hits = (
        gene_mirna_best
        .sort_values(["gene_id", "weighted_context_score"])
        .groupby("gene_id")
        .head(top_n)
    )

    # Apply threshold for strong repression
    top_hits = top_hits[top_hits["weighted_context_score"] <= weight_threshold].copy()

    print(f"Found {len(top_hits)} miRNA-gene pairs meeting threshold")

    # Save results if output directory provided
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        gene_mirna_best.to_csv(
            output_dir / "all_miRNAs_targetscan.csv",
            index=False,
        )
        top_hits.to_csv(
            output_dir / f"top_{top_n}_threshold_{weight_threshold}_targetscan.csv",
            index=False,
        )
        print(f"Saved miRNA results to {output_dir}")

    return top_hits


def summarize_mirnas_per_gene(mirna_hits: pd.DataFrame) -> pd.DataFrame:
    """
    Create gene-level summary of miRNA targeting.
    
    Returns DataFrame with:
    - gene_id
    - n_mirnas: Number of distinct miRNAs targeting this gene
    - top_mirna: Best scoring miRNA
    - top_score: Score of best miRNA
    - mirna_list: List of all targeting miRNAs
    """
    if mirna_hits.empty:
        return pd.DataFrame(columns=[
            "gene_id", "n_mirnas", "top_mirna", "top_score", "mirna_list"
        ])
    
    summary = (
        mirna_hits
        .sort_values(["gene_id", "weighted_context_score"])
        .groupby("gene_id")
        .agg({
            "miRNA": ["count", "first", lambda x: list(x)],
            "weighted_context_score": "first",
        })
    )
    
    summary.columns = ["n_mirnas", "top_mirna", "mirna_list", "top_score"]
    summary = summary.reset_index()
    
    return summary


def get_mirna_gene_matrix(mirna_hits: pd.DataFrame) -> pd.DataFrame:
    """
    Create gene × miRNA matrix of weighted context scores.
    
    Returns pivot table with genes as rows, miRNAs as columns,
    and scores as values (NaN where no prediction).
    """
    if mirna_hits.empty:
        return pd.DataFrame()
    
    matrix = mirna_hits.pivot_table(
        index="gene_id",
        columns="miRNA",
        values="weighted_context_score",
        aggfunc="min",  # Best score if multiple
    )
    
    return matrix
