"""
miRNA overlap matching for SNVs/indels.

Matches point mutations and small indels to overlapping miRNA loci
to identify variants that fall within annotated miRNA genes.
"""

from typing import Dict

import pandas as pd


def match_snvs_to_mirnas(
    snv_df: pd.DataFrame,
    mirna_df: pd.DataFrame,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
) -> pd.DataFrame:
    """
    Match SNVs to overlapping miRNA loci.

    For each SNV, finds all miRNA records that contain the variant position,
    storing results in a new 'miRNA_hits' column.

    Args:
        snv_df: DataFrame with SNV data (must have chrom, pos columns)
        mirna_df: miRNA annotation DataFrame with columns:
                  chrom, start, end, strand, gene_name, gene_id
        chrom_col: Name of chromosome column in snv_df
        pos_col: Name of position column in snv_df

    Returns:
        DataFrame with 'miRNA_hits' column containing list of dicts:
        [
            {
                "gene_name": str,
                "gene_id": str,
                "strand": str,
                "chrom": str,
                "mirna_start": int,
                "mirna_end": int,
            },
            ...
        ]
    """
    df = snv_df.copy()

    # Pre-index miRNA annotations by chromosome for faster lookup
    mirnas_by_chrom: Dict[str, pd.DataFrame] = {
        chrom: group for chrom, group in mirna_df.groupby("chrom")
    }

    mirna_hits_col = []

    for _, snv in df.iterrows():
        chrom = snv.get(chrom_col)
        pos = snv.get(pos_col)

        # Handle missing data
        if pd.isna(chrom) or pd.isna(pos):
            mirna_hits_col.append([])
            continue

        pos = int(pos)

        # Get miRNAs on same chromosome
        m_chr = mirnas_by_chrom.get(chrom)
        if m_chr is None or m_chr.empty:
            mirna_hits_col.append([])
            continue

        # Find overlapping miRNAs (pos within [start, end])
        overlapping = m_chr[
            (m_chr["start"] <= pos) & (m_chr["end"] >= pos)
        ]

        if overlapping.empty:
            mirna_hits_col.append([])
            continue

        # Build hit records
        hits = []
        for _, mirna in overlapping.iterrows():
            hit = {
                "gene_name": mirna.get("gene_name"),
                "gene_id": mirna.get("gene_id"),
                "strand": mirna.get("strand"),
                "chrom": chrom,
                "mirna_start": int(mirna["start"]),
                "mirna_end": int(mirna["end"]),
            }
            # Optional extra miRNA metadata (depends on which loci table is used).
            # - precursor table: mature_names/mature_accessions
            # - mature-arm table: mirbase_mature_id + pre_gene_* context
            for extra in (
                "mature_names",
                "mature_accessions",
                "mirbase_mature_id",
                "pre_gene_name",
                "pre_gene_id",
            ):
                if extra in mirna and pd.notna(mirna.get(extra)):
                    hit[extra] = mirna.get(extra)
            hits.append(hit)

        mirna_hits_col.append(hits)

    # Standardize to snake_case to match other "hits" columns in pipeline
    df["mirna_hits"] = mirna_hits_col
    # Backwards-compat alias (some downstream code may expect this exact spelling)
    df["miRNA_hits"] = df["mirna_hits"]
    return df


def summarize_mirna_hits(snv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary statistics for miRNA hits across all SNVs.

    Args:
        snv_df: DataFrame with 'miRNA_hits' column

    Returns:
        DataFrame with per-SNV summary:
        - n_mirna_hits: Number of miRNA loci overlapping
        - mirna_gene_names: Comma-separated unique gene names
        - mirna_gene_ids: Comma-separated unique gene IDs
        - mirna_strands: Comma-separated unique strands
    """
    def summarize_row(row):
        hits = row.get("miRNA_hits", [])
        if not hits:
            return pd.Series({
                "n_mirna_hits": 0,
                "mirna_gene_names": "",
                "mirna_gene_ids": "",
                "mirna_strands": "",
            })

        gene_names = sorted(
            set(h.get("gene_name", "") for h in hits if h.get("gene_name"))
        )
        gene_ids = sorted(
            set(h.get("gene_id", "") for h in hits if h.get("gene_id"))
        )
        strands = sorted(
            set(h.get("strand", "") for h in hits if h.get("strand"))
        )

        return pd.Series({
            "n_mirna_hits": len(hits),
            "mirna_gene_names": ",".join(gene_names),
            "mirna_gene_ids": ",".join(gene_ids),
            "mirna_strands": ",".join(strands),
        })

    summary = snv_df.apply(summarize_row, axis=1)
    return pd.concat([snv_df, summary], axis=1)