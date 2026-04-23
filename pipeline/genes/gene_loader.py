"""
Gene and lncRNA loading from GENCODE annotation.

Functions for:
- Loading gene annotations from GTF-derived CSV
- Filtering to specific gene panels
- Adding promoter coordinates
- Creating BED files
"""

import os
from pathlib import Path
from typing import Collection, List, Tuple, Optional, Set, Iterable

import numpy as np
import pandas as pd

from ..config import PATHS

from ..utils import (
    harmonize_chrom_column,
    compute_promoter_coords,
    compute_tss,
)


# =============================================================================
# LOADING
# =============================================================================

def _is_gene_only_slim_gencode(df: pd.DataFrame) -> bool:
    """True if ``df`` looks like a gene-only extract (no transcript/exon/CDS rows)."""
    if df is None or df.empty:
        return True
    if "feature" not in df.columns:
        return True
    uniq = set(df["feature"].dropna().astype(str).unique())
    return len(uniq) == 0 or uniq == {"gene"}


def load_gencode_multifeature_subset(gene_names: Collection[str]) -> pd.DataFrame:
    """
    Load GENCODE rows for every GTF ``feature`` type (gene, transcript, exon, CDS, …),
    restricted to ``gene_names`` (symbols in ``gene_name``).

    The default ``PATHS.gencode_gtf_pq`` is a **slim** gene-only Parquet for fast joins; this
    function then reads ``{working_dir}/gencode.v49.annotation.gtf.csv`` in **chunks** so the
    full annotation is never loaded into memory. If the Parquet already contains multiple
    ``feature`` values, it is filtered in-memory instead (expected to be rare / small).
    """
    names: Set[str] = {str(x) for x in gene_names if x is not None and str(x)}
    if not names:
        return pd.DataFrame()

    pq = Path(PATHS.gencode_gtf_pq)
    pq_full = Path(getattr(PATHS, "gencode_gtf_full_pq", ""))
    csv_path = Path(PATHS.working_dir) / "gencode.v49.annotation.gtf.csv"

    # Preferred runtime source: full multi-feature parquet (keeps exon_id/exon_number, etc.)
    if pq_full and pq_full.exists() and pq_full.suffix == ".parquet":
        try:
            # Predicate pushdown keeps this fast even if the parquet is large.
            df_full = pd.read_parquet(
                pq_full,
                filters=[("gene_name", "in", sorted(names))],
            )
            if df_full is not None and not df_full.empty:
                df_full, _ = harmonize_chrom_column(df_full)
                df_full["start"] = pd.to_numeric(df_full["start"], errors="coerce").astype("Int64")
                df_full["end"] = pd.to_numeric(df_full["end"], errors="coerce").astype("Int64")
                return df_full
        except Exception as e:
            print(f"[WARN] Could not read multifeature subset from {pq_full}: {e}")

    df: Optional[pd.DataFrame] = None
    if pq.exists() and pq.suffix == ".parquet":
        try:
            df = pd.read_parquet(pq)
        except Exception as e:
            print(f"[WARN] Could not read gene parquet {pq}: {e}")
            df = None

    if df is not None and not _is_gene_only_slim_gencode(df):
        df = df[df["gene_name"].isin(names)].copy()
        df, _ = harmonize_chrom_column(df)
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
        return df

    if df is not None:
        print(
            f"[WARN] Full GENCODE CSV not found ({csv_path}); "
            "multifeature gene exports may miss exon_id/exon_number metadata."
        )
        df = df[df["gene_name"].isin(names)].copy()
        df, _ = harmonize_chrom_column(df)
        df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
        return df
    # Last resort: CSV (large). Keep for back-compat, but not the intended runtime path.
    if csv_path.exists():
        print(f"[WARN] Loading multifeature GENCODE subset from CSV (large): {csv_path}")
        parts: List[pd.DataFrame] = []
        for chunk in pd.read_csv(csv_path, chunksize=400_000, low_memory=False):
            hit = chunk.loc[chunk["gene_name"].isin(names)]
            if not hit.empty:
                parts.append(hit)
        if not parts:
            return pd.DataFrame()
        out = pd.concat(parts, ignore_index=True)
        out, _ = harmonize_chrom_column(out)
        out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
        out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")
        return out

    raise FileNotFoundError(
        f"No multifeature GENCODE source: {pq_full} missing; parquet unusable; CSV missing: {csv_path}"
    )


def load_genes(path: Path) -> pd.DataFrame:
    """
    Load gene annotations from GENCODE GTF-derived Parquet (preferred),
    falling back to CSV if parquet support or file is unavailable.
    
    Expected columns: chrom/seqname, start, end, strand, gene_name, gene_id, feature, gene_type
    """
    path = Path(path)
    df = None
    if path.exists() and path.suffix == ".parquet":
        try:
            df = pd.read_parquet(path)
        except Exception as e:
            print(f"[WARN] Failed to read parquet genes at {path}: {e}")
            df = None

    if df is None:
        # Fallback: use CSV present on disk in this repo
        csv_fallback = getattr(PATHS, "working_dir", Path("."))
        csv_fallback = Path(csv_fallback) / "gencode.v49.annotation.gtf.csv"
        if not csv_fallback.exists():
            raise FileNotFoundError(
                f"Could not load genes from parquet ({path}) and CSV fallback not found: {csv_fallback}"
            )
        print(f"[WARN] Falling back to CSV genes: {csv_fallback}")
        df = pd.read_csv(csv_fallback)

    df, _ = harmonize_chrom_column(df)
    
    # Ensure numeric coordinates
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    
    return df


def load_lncrnas(path: Path) -> pd.DataFrame:
    """
    Load lncRNA annotations from CSV.

    Expects at least ``chrom``, ``start``, ``end``, ``strand``, ``gene_name`` (lncRNA id/symbol).
    GENCODE-style tables may include ``gene_type`` / ``feature``; gene-centric panel files
    (e.g. ``lncRNAs_with_genes_1000000bp.csv``) often omit ``gene_type`` because every row
    is already an lncRNA locus.
    """
    df = pd.read_csv(path)
    df, _ = harmonize_chrom_column(df)
    
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")

    if "gene_name" in df.columns:
        from .symbol_normalization import normalize_annotation_gene_names

        df = normalize_annotation_gene_names(df, ("gene_name",))

    return df


# =============================================================================
# FILTERING
# =============================================================================

def filter_genes_by_names(
    genes: pd.DataFrame,
    gene_names: List[str],
    feature_type: str = "gene",
) -> pd.DataFrame:
    """
    Filter genes to those in the specified name list.
    
    Args:
        genes: Full gene DataFrame
        gene_names: List of gene symbols to keep
        feature_type: GTF feature type to filter to (default: "gene")
    
    Returns:
        Filtered DataFrame
    """
    mask = genes["gene_name"].isin(gene_names)
    
    if "feature" in genes.columns and feature_type:
        mask &= genes["feature"] == feature_type
    
    return genes[mask].copy()


def lncrna_gene_intervals_from_annotation(genes_all: pd.DataFrame) -> pd.DataFrame:
    """
    All lncRNA ``feature == 'gene'`` intervals from a GENCODE-style ``genes_all`` frame.

    ``genes_all`` must be the same table passed to ``filter_genes_by_names`` for protein-coding
    anchors (typically ``load_genes(PATHS.gencode_gtf_pq)``): it must include a ``gene_type``
    column so lncRNAs are not confused with other biotypes.

    Applies the same symbol normalization as ``load_lncrnas``.
    """
    if "gene_type" not in genes_all.columns:
        raise ValueError(
            "lncrna_gene_intervals_from_annotation requires a `gene_type` column on the "
            "annotation (use a full GENCODE gene-level extract, e.g. from "
            "scripts/gencode/build_gencode_slim_parquet.py)."
        )
    lnc = filter_lncrnas(genes_all)
    if lnc.empty:
        n_lnc = int((genes_all["gene_type"].astype(str) == "lncRNA").sum())
        if n_lnc == 0:
            print(
                "[WARN] lncrna_gene_intervals_from_annotation: no rows with gene_type == "
                "'lncRNA' — lncRNA–gene matching will be empty."
            )
        else:
            print(
                f"[WARN] lncrna_gene_intervals_from_annotation: 0 rows after feature filter "
                f"(raw lncRNA gene rows in table: {n_lnc})."
            )
    from .symbol_normalization import normalize_annotation_gene_names

    return normalize_annotation_gene_names(lnc.copy(), ("gene_name",))


def filter_lncrnas(
    genes: pd.DataFrame,
    feature_type: str = "gene",
) -> pd.DataFrame:
    """
    Filter to lncRNA genes only.
    
    Args:
        genes: Full gene DataFrame
        feature_type: GTF feature type to filter to (default: "gene")
    
    Returns:
        DataFrame with only lncRNA genes
    """
    if "gene_type" in genes.columns:
        mask = genes["gene_type"] == "lncRNA"
    else:
        # Gene-centric lncRNA panel CSVs (no GENCODE ``gene_type`` column): keep all rows.
        mask = pd.Series(True, index=genes.index)
    
    if "feature" in genes.columns and feature_type:
        mask &= genes["feature"] == feature_type
    
    return genes[mask].copy()


# =============================================================================
# PROMOTER COORDINATES
# =============================================================================

def add_promoter_columns(
    genes: pd.DataFrame,
    upstream_bp: int = 2000,
    downstream_bp: int = 500,
) -> pd.DataFrame:
    """
    Add promoter coordinate columns (prom_start, prom_end) to gene DataFrame.
    
    Promoter is defined relative to TSS:
    - For + strand: [TSS - upstream, TSS + downstream]
    - For - strand: [TSS - downstream, TSS + upstream]
    
    Also adds 'tss' column.
    """
    genes = genes.copy()
    
    # Compute TSS
    genes["tss"] = compute_tss(genes["start"], genes["end"], genes["strand"])
    
    # Compute promoter coordinates
    prom_start, prom_end = compute_promoter_coords(
        genes["start"],
        genes["end"],
        genes["strand"],
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
    )
    
    genes["prom_start"] = prom_start.astype("Int64")
    genes["prom_end"] = prom_end.astype("Int64")
    
    return genes


def add_tss_window(
    genes: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    """
    Add TSS-centered window columns (win_start, win_end) to gene DataFrame.
    """
    genes = genes.copy()
    
    if "tss" not in genes.columns:
        genes["tss"] = compute_tss(genes["start"], genes["end"], genes["strand"])
    
    genes["win_start"] = (genes["tss"] - window_bp).clip(lower=0).astype("Int64")
    genes["win_end"] = (genes["tss"] + window_bp).astype("Int64")
    
    return genes


# =============================================================================
# OUTPUT FORMATS
# =============================================================================

def create_genes_bed(
    genes: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Create BED file from gene DataFrame.
    
    BED format is 0-based, half-open: [start, end)
    """
    df_bed = genes[["chrom", "start", "end", "gene_name", "strand"]].copy()
    
    # Convert to 0-based
    df_bed["start"] = df_bed["start"] - 1
    
    # Ensure no header, tab-separated
    df_bed.to_csv(output_path, sep="\t", header=False, index=False)
    print(f"Created BED file: {output_path}")


def _write_tier_gene_tables(
    gff: pd.DataFrame,
    genes_all_harmonized: pd.DataFrame,
    tier_symbols: Iterable[str],
    genes_only_path: Path,
    all_features_path: Path,
    label: str,
) -> None:
    """Gene-level + multifeature CSVs for one tier list (Tier 2 / 3 / 4)."""
    from ..config import THRESHOLDS

    names = [str(x) for x in tier_symbols if x is not None and str(x)]
    if not names:
        return
    sub = filter_genes_by_names(genes_all_harmonized, names)
    if sub.empty:
        print(f"[WARN] save_gene_tables: no GENCODE `gene` rows for {label} (requested {len(names)} symbols)")
        sub_out = sub
    else:
        sub_out = add_promoter_columns(
            sub,
            upstream_bp=THRESHOLDS.promoter_upstream_bp,
            downstream_bp=THRESHOLDS.promoter_downstream_bp,
        )
    sub_out.to_csv(genes_only_path, index=False)
    tier_all = gff[gff["gene_name"].isin(names)].copy()
    tier_all.to_csv(all_features_path, index=False)
    print(
        f"  {label}: {genes_only_path.name} ({len(sub_out)} gene rows), "
        f"{all_features_path.name} ({len(tier_all)} multifeature rows)"
    )


def save_gene_tables(
    genes: pd.DataFrame,
    cnv_genes: pd.DataFrame,
    gene_panel: List[str],
    lncrna_names: List[str],
    output_dir: Path,
    genes_all_harmonized: Optional[pd.DataFrame] = None,
    gencode_genes_full: Optional[pd.DataFrame] = None,
) -> None:
    """
    Save filtered gene tables to CSV.

    ``PATHS.genes_*`` rows follow ``gene_panel`` (typically ``PIPELINE_GENE_PANEL`` =
    ``FULL_INTEGRATION_GENES``). Legacy ``primary_*`` filenames are not limited to the 66
    ``PRIMARY_GENES``. Tier 2–4 tables require ``genes_all_harmonized``.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from ..config import (
        TIER2_MEDIUM_GENES,
        TIER3_CNV_ONLY_GENES,
        TIER4_READOUT_GENES,
    )

    need_names: Set[str] = {str(x) for x in gene_panel if x is not None and str(x)}
    need_names |= {str(x) for x in lncrna_names if x is not None and str(x)}
    if "gene_name" in cnv_genes.columns:
        need_names |= set(cnv_genes["gene_name"].dropna().astype(str))
    need_names |= set(TIER2_MEDIUM_GENES)
    need_names |= set(TIER3_CNV_ONLY_GENES)
    need_names |= set(TIER4_READOUT_GENES)

    gff = gencode_genes_full
    if gff is None or _is_gene_only_slim_gencode(gff):
        gff = load_gencode_multifeature_subset(need_names)
    else:
        gff = gff[gff["gene_name"].isin(need_names)].copy()

    primary_all = gff[gff["gene_name"].isin(gene_panel)].copy()
    primary_all.to_csv(PATHS.genes_all_features, index=False)

    primary_gene = genes[genes["gene_name"].isin(gene_panel)].copy()
    primary_gene.to_csv(PATHS.genes_only, index=False)

    cnv_gene = cnv_genes[cnv_genes.get("feature", "gene") == "gene"]
    cnv_gene.to_csv(PATHS.cnv_genes, index=False)

    lncrna_all = gff[gff["gene_name"].isin(lncrna_names)].copy()
    lncrna_all.to_csv(PATHS.lncrnas_all_features, index=False)

    if genes_all_harmonized is not None and not genes_all_harmonized.empty:
        _write_tier_gene_tables(
            gff,
            genes_all_harmonized,
            TIER2_MEDIUM_GENES,
            PATHS.tier2_medium_genes_only,
            PATHS.tier2_medium_genes_all_features,
            "Tier 2 medium-depth genes",
        )
        _write_tier_gene_tables(
            gff,
            genes_all_harmonized,
            TIER3_CNV_ONLY_GENES,
            PATHS.tier3_cnv_only_genes_only,
            PATHS.tier3_cnv_only_genes_all_features,
            "Tier 3 CNV-only genes",
        )
        _write_tier_gene_tables(
            gff,
            genes_all_harmonized,
            TIER4_READOUT_GENES,
            PATHS.tier4_readout_genes_only,
            PATHS.tier4_readout_genes_all_features,
            "Tier 4 readout genes",
        )
    else:
        print(
            "[WARN] save_gene_tables: genes_all_harmonized not passed — "
            "skipping tier2/tier3/tier4 gene CSV exports."
        )

    print(f"Saved gene tables to {output_dir}")


# =============================================================================
# HARMONIZATION ACROSS MULTIPLE DATAFRAMES
# =============================================================================

def harmonize_multiple_dfs(
    dfs: List[pd.DataFrame],
    paths: Optional[List[Path]] = None,
    save_if_changed: bool = False,
) -> Tuple[List[pd.DataFrame], List[bool]]:
    """
    Harmonize chromosome columns across multiple DataFrames.
    
    Args:
        dfs: List of DataFrames to harmonize
        paths: Optional list of paths for saving modified DataFrames
        save_if_changed: Whether to save back to paths if changes were made
    
    Returns:
        Tuple of (harmonized DataFrames, list of was_changed flags)
    """
    results = []
    changes = []
    
    for i, df in enumerate(dfs):
        df_harmonized, was_changed = harmonize_chrom_column(df)
        results.append(df_harmonized)
        changes.append(was_changed)
        
        if save_if_changed and was_changed and paths and i < len(paths):
            df_harmonized.to_csv(paths[i], index=False)
            print(f"Saved harmonized DataFrame to {paths[i]}")
    
    return results, changes
