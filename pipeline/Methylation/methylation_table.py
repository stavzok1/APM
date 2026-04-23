"""
Main methylation pipeline orchestrator.

Coordinates all methylation module components to build:
- Annotated probe reference table
- Per-sample methylation tables (probe and aggregated levels)
- Cohort-level methylation matrices

This is the primary entry point for methylation integration.
"""

from __future__ import annotations

import gc
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import pandas as pd
import numpy as np

from .probe_loader import (
    load_probe_reference,
    annotate_probes_full,
    annotate_probes_with_tads,
    resolve_methylation_ccre_id_column,
)
from .sample_processing import (
    load_sample_beta,
    process_sample_methylation,
    compute_sample_qc_metrics,
    enrich_sample_with_annotations,
)
from .aggregation import (
    aggregate_to_genes,
    aggregate_to_lncrnas,
    aggregate_to_ccres,
    aggregate_to_atac,
    aggregate_sample_batch_to_genes,
    aggregate_sample_batch_to_ccres,
)
from ..config import THRESHOLDS, PATHS, PIPELINE_GENE_PANEL
from ..genes.symbol_normalization import normalize_annotation_gene_names
from ..utils import harmonize_chrom_column
from ..sample_ids import add_tcga_id_columns_inplace

from ..tad_annotation import load_tad_source, discover_tad_sources


from .meth_schemas import METHYLATION_COLUMNS as MC

# Columns used by enrich_sample_with_annotations + aggregation. When reusing a
# saved probe reference, we load only these (via Parquet column projection) so
# huge optional columns (e.g. TAD_domains) do not OOM on read.
METHYLATION_PROBE_REFERENCE_LOAD_COLUMNS: Tuple[str, ...] = (
    "probeID",
    "chrom",
    "start",
    "end",
    "strand",
    "center",
    "gene_list",
    "distToTSS_parsed",
    "in_CGI",
    "CGI_context",
    "in_promoter",
    "promoter_genes",
    "promoter_gene_ids",
    "in_gene_body",
    "gene_body_genes",
    "in_cds",
    "cds_genes",
    "overlapping_ccres",
    "ccre_types",
    "n_overlapping_ccres",
    "in_lncrna_promoter",
    "lncrna_promoter_genes",
    "overlapping_atac_peaks",
    "n_overlapping_atac",
)


def load_existing_methylation_probe_reference(path: Path) -> pd.DataFrame:
    """
    Load a saved annotated probe reference without pulling unused heavy columns.

    By default skips columns not needed for the methylation pipeline (notably
    ``TAD_domains``). Set ``APM_METH_LOAD_FULL_PROBE_REFERENCE=1`` to read every
    column (high RAM if TAD was stored in full).
    """
    path = Path(path)
    full = os.environ.get("APM_METH_LOAD_FULL_PROBE_REFERENCE", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    if path.suffix.lower() == ".parquet":
        if full:
            df = pd.read_parquet(path)
            print(f"  Loaded {len(df)} probes (all Parquet columns) from {path}")
            return df
        try:
            import pyarrow.parquet as pq

            schema_cols = set(pq.ParquetFile(path).schema_arrow.names)
        except Exception:
            schema_cols = set(pd.read_parquet(path, engine="pyarrow").columns)
        want = [c for c in METHYLATION_PROBE_REFERENCE_LOAD_COLUMNS if c in schema_cols]
        if "probeID" not in want:
            df = pd.read_parquet(path)
            print(f"  Loaded {len(df)} probes from {path} (full; unexpected schema)")
            return df
        skip = sorted(schema_cols - set(want))
        if skip:
            preview = skip[:4]
            more = "…" if len(skip) > len(preview) else ""
            print(
                f"  Loading {len(want)}/{len(schema_cols)} columns "
                f"(skipping heavy/unused: {preview}{more}) — "
                f"APM_METH_LOAD_FULL_PROBE_REFERENCE=1 loads all."
            )
        df = pd.read_parquet(path, columns=want)
        print(f"  Loaded {len(df)} probes from {path}")
        return df

    # CSV / TSV (no column projection in pandas if we cannot parse header cheaply)
    if full:
        df = pd.read_csv(path, low_memory=False)
        print(f"  Loaded {len(df)} probes (all columns) from {path}")
        return df
    try:
        hdr = pd.read_csv(path, nrows=0, low_memory=False)
    except Exception:
        df = pd.read_csv(path, low_memory=False)
        print(f"  Loaded {len(df)} probes from {path}")
        return df
    names = set(hdr.columns)
    want = [c for c in METHYLATION_PROBE_REFERENCE_LOAD_COLUMNS if c in names]
    if "probeID" not in want:
        df = pd.read_csv(path, low_memory=False)
        print(f"  Loaded {len(df)} probes from {path} (full; no probeID in subset)")
        return df
    df = pd.read_csv(path, usecols=want, low_memory=False)
    if len(want) < len(names):
        print(
            f"  Loading {len(want)}/{len(names)} CSV columns "
            f"(APM_METH_LOAD_FULL_PROBE_REFERENCE unset)…"
        )
    print(f"  Loaded {len(df)} probes from {path}")
    return df


def _write_probe_reference_parquet(probes: pd.DataFrame, output_path: Path) -> None:
    """
    Write probe reference Parquet in row chunks to cap PyArrow conversion RAM.

    Env ``APM_METH_PROBE_PARQUET_CHUNK`` (default ``50000``): rows per slice.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        probes.to_parquet(output_path, index=False)
        return

    chunk = int(os.environ.get("APM_METH_PROBE_PARQUET_CHUNK", "50000"))
    chunk = max(5_000, chunk)
    n = len(probes)
    writer = None
    try:
        for start in range(0, n, chunk):
            sub = probes.iloc[start : start + chunk]
            table = pa.Table.from_pandas(sub, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    str(output_path), table.schema, compression="zstd"
                )
            writer.write_table(table)
            del table, sub
            gc.collect()
    finally:
        if writer is not None:
            writer.close()


# =============================================================================
# PROBE REFERENCE TABLE BUILDER
# =============================================================================

def build_probe_reference_table(
    probe_reference_path: Path,
    genes: pd.DataFrame,
    genes_cds: pd.DataFrame,
    ccres: pd.DataFrame,
    lncrnas: Optional[pd.DataFrame] = None,
    atac_peaks: Optional[pd.DataFrame] = None,
    gene_panel: Optional[List[str]] = None,
    lncrna_names: Optional[List[str]] = None,
    tad_sources_dir: Optional[Path] = None,
    include_tads: bool = True,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
    output_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Build the complete annotated probe reference table.
    
    This is the master probe annotation that gets built once and reused
    across all samples.
    
    Args:
        probe_reference_path: Path to GDC probe reference file
        genes: Gene DataFrame (from GENCODE)
        genes_cds: Gene CDS DataFrame (from GENCODE)
        ccres: cCRE DataFrame (from ENCODE)
        lncrnas: Optional lncRNA DataFrame
        atac_peaks: Optional ATAC peaks DataFrame
        gene_panel: List of genes for promoter annotation
        lncrna_names: List of lncRNAs for annotation
        tad_sources_dir: Directory containing TAD sources
        upstream_bp: Promoter upstream of TSS
        downstream_bp: Promoter downstream of TSS
        output_path: Optional path to save the annotated reference
    
    Returns:
        Fully annotated probe reference DataFrame
    """
    print("=" * 70)
    print("BUILDING PROBE REFERENCE TABLE")
    print("=" * 70)
    
    # Load raw probe reference
    probes = load_probe_reference(probe_reference_path)
    
    # Apply all annotations
    probes = annotate_probes_full(
        probes,
        genes=genes,
        genes_cds=genes_cds,
        ccres=ccres,
        lncrnas=lncrnas,
        atac_peaks=atac_peaks,
        gene_panel=gene_panel,
        lncrna_names=lncrna_names,
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
    )

    probes = probes[(probes["in_promoter"] == True) | (probes["in_gene_body"] == True) |
     (probes["n_overlapping_ccres"] > 0 ) | (probes["in_lncrna_promoter"] == True) | (probes["n_overlapping_atac"] > 0 )]
    
    # TAD annotation (per biosample)
    if include_tads and tad_sources_dir:
        tad_sources_dict = discover_tad_sources(tad_sources_dir)
        print("\nAnnotating with TAD domains...")
        for biosample, tad_source_paths in tad_sources_dict.items():
            domains, boundaries, flanks = load_tad_source(tad_source_paths)
            probes = annotate_probes_with_tads(
                probes,
                tad_domains=domains,
                domain_flanks=flanks,
                boundaries=boundaries,
                biosample=biosample,
            )
            del domains, boundaries, flanks
            gc.collect()  # drop TAD tables before the next biosample
    
    # Save if output path provided (parquet preserves list/dict columns; CSV is lossy for nested cells)
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix.lower() == ".parquet":
            print(f"\nWriting probe reference Parquet (chunked): {output_path} ...")
            _write_probe_reference_parquet(probes, output_path)
            print(f"Saved probe reference to: {output_path}")
        else:
            probes.to_csv(output_path, index=False)
            print(f"\nSaved probe reference to: {output_path}")
    
    print("\n" + "=" * 70)
    print("PROBE REFERENCE COMPLETE")
    print("=" * 70)
    
    return probes


# =============================================================================
# PER-SAMPLE METHYLATION TABLES
# =============================================================================

def build_sample_methylation_tables(
    sample_path: Path,
    probe_reference: pd.DataFrame,
    sample_id: str,
    gene_panel: List[str],
    lncrna_panel: Optional[List[str]] = None,
    ccre_ids: Optional[List[str]] = None,
    output_dir: Optional[Path] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Build all methylation tables for a single sample.
    
    Produces:
        - probe_table: Per-probe beta values with annotations
        - gene_table: Gene-level aggregation
        - lncrna_table: lncRNA-level aggregation (if panel provided)
        - ccre_table: cCRE-level aggregation
    
    Args:
        sample_path: Path to sample beta file
        probe_reference: Annotated probe reference
        sample_id: Sample identifier
        gene_panel: List of genes for aggregation
        lncrna_panel: Optional list of lncRNAs
        ccre_ids: Optional list of cCRE IDs (None = all with probes)
        output_dir: Optional directory to save tables
    
    Returns:
        Dict with keys: "probes", "genes", "lncrnas", "ccres"
    """
    print(f"\n{'='*60}")
    print(f"Processing sample: {sample_id}")
    print(f"{'='*60}")
    
    # Load and process sample
    sample_df = process_sample_methylation(
        sample_path,
        probe_reference,
        sample_id=sample_id,
        compute_m=True,
        enrich=True,
    )
    
    # QC metrics
    qc = compute_sample_qc_metrics(sample_df)
    print(f"\nQC Metrics:")
    print(f"  Valid probes: {qc['n_valid']} ({qc['pct_valid']:.1f}%)")
    print(f"  Mean beta: {qc['mean_beta']:.3f}" if qc['mean_beta'] else "  Mean beta: N/A")
    
    results = {"probes": sample_df}
    
    # Gene aggregation
    print("\nAggregating to gene level...")
    gene_agg = aggregate_to_genes(
        sample_df,
        probe_reference,
        gene_panel,
        include_gene_body=True,
    )
    gene_agg["sample_id"] = sample_id
    results["genes"] = gene_agg
    
    n_genes_with_data = (gene_agg["promoter_n_probes"] > 0).sum()
    print(f"  Genes with promoter data: {n_genes_with_data}/{len(gene_panel)}")
    
    # lncRNA aggregation

    print("Aggregating to lncRNA level...")
    lncrna_agg = aggregate_to_lncrnas(
        sample_df,
        probe_reference,
        lncrna_panel,
    )
    lncrna_agg["sample_id"] = sample_id
    results["lncrnas"] = lncrna_agg
    
    n_lnc_with_data = int((lncrna_agg["promoter_n_probes"] > 0).sum()) if "promoter_n_probes" in lncrna_agg.columns else 0
    if lncrna_panel:
        print(f"  lncRNAs with promoter data: {n_lnc_with_data}/{len(lncrna_panel)}")
    else:
        n_lnc_total = len(lncrna_agg)
        print(f"  lncRNAs with promoter data: {n_lnc_with_data}/{n_lnc_total} (no fixed panel; observed lncRNAs only)")
    
    # cCRE aggregation
    print("Aggregating to cCRE level...")
    ccre_agg = aggregate_to_ccres(
        sample_df,
        probe_reference,
        ccre_ids=ccre_ids,
    )
    ccre_agg["sample_id"] = sample_id
    results["ccres"] = ccre_agg

    atac_agg = aggregate_to_atac(
        sample_df,
        probe_reference,
    )
    atac_agg["sample_id"] = sample_id
    results["atac"] = atac_agg
    
    print(f"  cCREs with data: {len(ccre_agg)}")
    
    # Save outputs
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Probes (parquet for nested columns)
        # probes_path = output_dir / f"{sample_id}_probes.parquet"
        # sample_df.to_parquet(probes_path, index=False)

        add_tcga_id_columns_inplace(sample_df, raw_id=sample_id)
        sample_df.to_csv(output_dir / f"{sample_id}_probes.csv", index=False)
        
        # Aggregations (CSV is fine, no nested columns)
        gene_path = output_dir / f"{sample_id}_gene_meth.csv"
        add_tcga_id_columns_inplace(gene_agg, raw_id=sample_id)
        gene_agg.to_csv(gene_path, index=False)
        
        if "lncrnas" in results:
            print("saving lncRNA level...")
            lncrna_path = output_dir / f"{sample_id}_lncrna_meth.csv"
            add_tcga_id_columns_inplace(lncrna_agg, raw_id=sample_id)
            lncrna_agg.to_csv(lncrna_path, index=False)
        
        ccre_path = output_dir / f"{sample_id}_ccre_meth.csv"
        add_tcga_id_columns_inplace(ccre_agg, raw_id=sample_id)
        ccre_agg.to_csv(ccre_path, index=False)
        print("saving lncRNA level...")
        atac_path = output_dir / f"{sample_id}_atac_meth.csv"
        add_tcga_id_columns_inplace(atac_agg, raw_id=sample_id)
        atac_agg.to_csv(atac_path, index=False)
        
        # QC metrics
        qc_path = output_dir / f"{sample_id}_qc.json"
        with open(qc_path, "w") as f:
            json.dump(qc, f, indent=2)
        
        print(f"\nSaved to: {output_dir}")
    
    return results


# =============================================================================
# COHORT-LEVEL MATRICES
# =============================================================================

# def build_cohort_matrices(
#     sample_paths: List[Path],
#     sample_ids: List[str],
#     probe_reference: pd.DataFrame,
#     gene_panel: List[str],
#     lncrna_panel: Optional[List[str]] = None,
#     ccre_ids: Optional[List[str]] = None,
#     output_dir: Optional[Path] = None,
# ) -> Dict[str, pd.DataFrame]:
#     """
#     Build cohort-level methylation matrices (features x samples).
    
#     Produces:
#         - gene_matrix: Genes x Samples (promoter beta mean)
#         - lncrna_matrix: lncRNAs x Samples (if panel provided)
#         - ccre_matrix: cCREs x Samples
    
#     Args:
#         sample_paths: List of paths to sample beta files
#         sample_ids: List of sample identifiers
#         probe_reference: Annotated probe reference
#         gene_panel: List of genes
#         lncrna_panel: Optional list of lncRNAs
#         ccre_ids: Optional list of cCRE IDs
#         output_dir: Optional directory to save matrices
    
#     Returns:
#         Dict with keys: "genes", "lncrnas", "ccres"
#     """
#     print("\n" + "=" * 70)
#     print("BUILDING COHORT METHYLATION MATRICES")
#     print("=" * 70)
#     print(f"Samples: {len(sample_paths)}")
#     print(f"Genes: {len(gene_panel)}")
    
#     # Process all samples and collect aggregations
#     gene_results = {}
#     lncrna_results = {}
#     ccre_results = {}
#     atac_results = {}
#     qc_metrics = []
    
#     for i, (sample_path, sample_id) in enumerate(zip(sample_paths, sample_ids)):
#         print(f"\n[{i+1}/{len(sample_paths)}] Processing {sample_id}...")
        
#         try:
#             # Load sample
#             sample_df = load_sample_beta(sample_path)
#             sample_df = enrich_sample_with_annotations(sample_df, probe_reference)

#             print("sample_df.columns", sample_df.columns)
            
#             # QC
#             qc = compute_sample_qc_metrics(sample_df)
#             qc["sample_id"] = sample_id
#             qc_metrics.append(qc)
            
#             # Gene aggregation
#             gene_agg = aggregate_to_genes(
#                 sample_df,
#                 probe_reference,
#                 gene_panel,
#                 include_gene_body=False,
#             )
#             # print("gene_agg.columns", gene_agg.columns)

#             gene_results[sample_id] = gene_agg.set_index("gene_name")["promoter_beta_mean"]
            
#             # lncRNA aggregation
       
#             lncrna_agg = aggregate_to_lncrnas(
#                 sample_df,
#                 probe_reference,
#                 lncrna_panel,
#             )
#             # print("lncrna_agg.columns", lncrna_agg.columns)

#             lncrna_results[sample_id] = lncrna_agg.set_index("lncrna_name")["promoter_beta_mean"]
            
#             # cCRE aggregation
#             ccre_agg = aggregate_to_ccres(
#                 sample_df,
#                 probe_reference,
#                 ccre_ids=ccre_ids,
#             )
#             print("ccre_agg.columns", ccre_agg.columns)
#             atac_agg = aggregate_to_atac(
#                 sample_df,
#                 probe_reference,
#             )
#             print("atac_agg.columns", atac_agg.columns)
#             atac_results[sample_id] = atac_agg.set_index("atac_peak_id")["atac_beta_mean"]


#             if not ccre_agg.empty:
#                 ccre_results[sample_id] = ccre_agg.set_index("cCRE_id")["ccre_beta_mean"]
            
#             # Memory management
#             del sample_df
#             gc.collect()
            
#         except Exception as e:
#             print(f"  Error: {e}")
#             continue
    
#     # Build matrices
#     results = {}
    
#     # Gene matrix
#     if gene_results:
#         gene_matrix = pd.DataFrame(gene_results)
#         gene_matrix.index.name = "gene_name"
#         results["genes"] = gene_matrix
#         print(f"\nGene matrix: {gene_matrix.shape}")
    
#     # lncRNA matrix
#     if lncrna_results:
#         lncrna_matrix = pd.DataFrame(lncrna_results)
#         lncrna_matrix.index.name = "lncrna_name"
#         results["lncrnas"] = lncrna_matrix
#         print(f"lncRNA matrix: {lncrna_matrix.shape}")
    
#     # cCRE matrix
#     if ccre_results:
#         ccre_matrix = pd.DataFrame(ccre_results)
#         ccre_matrix.index.name = "cCRE_id"
#         results["ccres"] = ccre_matrix
#         print(f"cCRE matrix: {ccre_matrix.shape}")
    

#     # ATAC matrix
#     if atac_results:
#         atac_matrix = pd.DataFrame(atac_results)
#         atac_matrix.index.name = "atac_peak_id"
#         results["atac"] = atac_matrix
#         print(f"ATAC matrix: {atac_matrix.shape}")

#     # QC summary
#     if qc_metrics:
#         qc_df = pd.DataFrame(qc_metrics)
#         results["qc_summary"] = qc_df
    
#     # Save outputs
#     if output_dir:
#         output_dir = Path(output_dir)
#         output_dir.mkdir(parents=True, exist_ok=True)
        
#         if "genes" in results:
#             results["genes"].to_csv(output_dir / "gene_meth_matrix.csv")
#             results["genes"].to_parquet(output_dir / "gene_meth_matrix.parquet")
        
#         if "lncrnas" in results:
#             results["lncrnas"].to_csv(output_dir / "lncrna_meth_matrix.csv")
        
#         if "ccres" in results:
#             results["ccres"].to_csv(output_dir / "ccre_meth_matrix.csv")
#             results["ccres"].to_parquet(output_dir / "ccre_meth_matrix.parquet")
        
#         if "atac" in results:
#             results["atac"].to_csv(output_dir / "atac_meth_matrix.csv")
#             results["atac"].to_parquet(output_dir / "atac_meth_matrix.parquet")
        
#         if "qc_summary" in results:
#             results["qc_summary"].to_csv(output_dir / "sample_qc_summary.csv", index=False)
        
#         print(f"\nSaved cohort matrices to: {output_dir}")
    
#     return results


def build_cohort_matrices(
    sample_paths: List[Path],
    sample_ids: List[str],
    probe_reference: pd.DataFrame,
    gene_panel: List[str],
    lncrna_panel: Optional[List[str]] = None,
    ccre_ids: Optional[List[str]] = None,
    output_dir: Optional[Path] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Build cohort-level methylation matrices (features x samples) + per-feature meta.

    Produces:
        - genes: Genes x Samples (promoter_beta_mean) + meta column
        - lncrnas: lncRNAs x Samples (promoter_beta_mean) + meta column
        - ccres: cCREs x Samples (ccre_beta_mean) + meta column
        - atac: ATAC peaks x Samples (atac_beta_mean) + meta column
        - qc_summary: per-sample QC table

    For each matrix, adds a column "meta" where each row contains:
        {sample_id: {"n_probes": <int or None>}, ...}
    """
    print("\n" + "=" * 70)
    print("BUILDING COHORT METHYLATION MATRICES")
    print("=" * 70)
    print(f"Samples: {len(sample_paths)}")
    print(f"Genes: {len(gene_panel)}")

    # beta values (feature -> sample -> value)
    gene_results = {}
    lncrna_results = {}
    ccre_results = {}
    atac_results = {}

    # n_probes (feature -> sample -> n_probes)
    gene_nprobes = {}
    lncrna_nprobes = {}
    ccre_nprobes = {}
    atac_nprobes = {}

    qc_metrics = []

    for i, (sample_path, sample_id) in enumerate(zip(sample_paths, sample_ids)):
        print(f"\n[{i+1}/{len(sample_paths)}] Processing {sample_id}...")

        try:
            # Load + annotate
            sample_df = load_sample_beta(sample_path)
            sample_df = enrich_sample_with_annotations(sample_df, probe_reference)

            # QC
            qc = compute_sample_qc_metrics(sample_df)
            qc["sample_id"] = sample_id
            qc_metrics.append(qc)

            # ----------------
            # Gene aggregation
            # ----------------
            gene_agg = aggregate_to_genes(
                sample_df,
                probe_reference,
                gene_panel,
                include_gene_body=False,
            )
            if not gene_agg.empty:
                gi = gene_agg.set_index("gene_name")
                gene_results[sample_id] = gi["promoter_beta_mean"]
                if "promoter_n_probes" in gi.columns:
                    gene_nprobes[sample_id] = gi["promoter_n_probes"]

            # ------------------
            # lncRNA aggregation
            # ------------------
            lncrna_agg = aggregate_to_lncrnas(
                sample_df,
                probe_reference,
                lncrna_panel,
            )
            if not lncrna_agg.empty:
                li = lncrna_agg.set_index("lncrna_name")
                lncrna_results[sample_id] = li["promoter_beta_mean"]
                if "promoter_n_probes" in li.columns:
                    lncrna_nprobes[sample_id] = li["promoter_n_probes"]

            # ----------------
            # cCRE aggregation
            # ----------------
            ccre_agg = aggregate_to_ccres(
                sample_df,
                probe_reference,
                ccre_ids=ccre_ids,
            )
            if not ccre_agg.empty:
                ci = ccre_agg.set_index("cCRE_id")
                ccre_results[sample_id] = ci["ccre_beta_mean"]
                if "ccre_n_probes" in ci.columns:
                    ccre_nprobes[sample_id] = ci["ccre_n_probes"]

            # ----------------
            # ATAC aggregation
            # ----------------
            atac_agg = aggregate_to_atac(
                sample_df,
                probe_reference,
            )
            if not atac_agg.empty:
                ai = atac_agg.set_index("atac_peak_id")
                atac_results[sample_id] = ai["atac_beta_mean"]
                if "atac_n_probes" in ai.columns:
                    atac_nprobes[sample_id] = ai["atac_n_probes"]

            # Memory management
            del sample_df
            gc.collect()

        except Exception as e:
            print(f"  Error: {e}")
            continue

    def _attach_meta(beta_mat: pd.DataFrame, nprobe_mat: Optional[pd.DataFrame]) -> pd.DataFrame:
        """
        Add a 'meta' column where each row is:
            {sample_id: {"n_probes": <value or None>}, ...}
        """
        if nprobe_mat is None or nprobe_mat.empty:
            # still add meta (all None) for consistency
            samples = list(beta_mat.columns)
            meta_col = [{sid: {"n_probes": None} for sid in samples} for _ in range(beta_mat.shape[0])]
            out = beta_mat.copy()
            out["meta"] = meta_col
            return out

        # align indices/columns
        nprobe_mat = nprobe_mat.reindex(index=beta_mat.index, columns=beta_mat.columns)

        samples = list(beta_mat.columns)

        # build row-wise dicts efficiently
        meta_col = []
        for _, row in nprobe_mat.iterrows():
            d = {}
            for sid in samples:
                v = row.get(sid)
                if pd.isna(v):
                    d[sid] = {"n_probes": None}
                else:
                    # keep as int when possible
                    try:
                        d[sid] = {"n_probes": int(v)}
                    except Exception:
                        d[sid] = {"n_probes": v}
            meta_col.append(d)

        out = beta_mat.copy()
        out["meta"] = meta_col
        return out

    results: Dict[str, pd.DataFrame] = {}

    # ----------------
    # Gene matrix
    # ----------------
    if gene_results:
        gene_matrix = pd.DataFrame(gene_results)
        gene_matrix.index.name = "gene_name"
        gene_nprobe_matrix = pd.DataFrame(gene_nprobes) if gene_nprobes else pd.DataFrame(index=gene_matrix.index)
        gene_nprobe_matrix.index.name = "gene_name"

        gene_matrix = _attach_meta(gene_matrix, gene_nprobe_matrix)
        results["genes"] = gene_matrix
        print(f"\nGene matrix (+meta): {gene_matrix.shape}")

    # ----------------
    # lncRNA matrix
    # ----------------
    if lncrna_results:
        lncrna_matrix = pd.DataFrame(lncrna_results)
        lncrna_matrix.index.name = "lncrna_name"
        lncrna_nprobe_matrix = pd.DataFrame(lncrna_nprobes) if lncrna_nprobes else pd.DataFrame(index=lncrna_matrix.index)
        lncrna_nprobe_matrix.index.name = "lncrna_name"

        lncrna_matrix = _attach_meta(lncrna_matrix, lncrna_nprobe_matrix)
        results["lncrnas"] = lncrna_matrix
        print(f"lncRNA matrix (+meta): {lncrna_matrix.shape}")

    # ----------------
    # cCRE matrix
    # ----------------
    if ccre_results:
        ccre_matrix = pd.DataFrame(ccre_results)
        ccre_matrix.index.name = "cCRE_id"
        ccre_nprobe_matrix = pd.DataFrame(ccre_nprobes) if ccre_nprobes else pd.DataFrame(index=ccre_matrix.index)
        ccre_nprobe_matrix.index.name = "cCRE_id"

        ccre_matrix = _attach_meta(ccre_matrix, ccre_nprobe_matrix)
        results["ccres"] = ccre_matrix
        print(f"cCRE matrix (+meta): {ccre_matrix.shape}")

    # ----------------
    # ATAC matrix
    # ----------------
    if atac_results:
        atac_matrix = pd.DataFrame(atac_results)
        atac_matrix.index.name = "atac_peak_id"
        atac_nprobe_matrix = pd.DataFrame(atac_nprobes) if atac_nprobes else pd.DataFrame(index=atac_matrix.index)
        atac_nprobe_matrix.index.name = "atac_peak_id"

        atac_matrix = _attach_meta(atac_matrix, atac_nprobe_matrix)
        results["atac"] = atac_matrix
        print(f"ATAC matrix (+meta): {atac_matrix.shape}")

    # QC summary
    if qc_metrics:
        qc_df = pd.DataFrame(qc_metrics)
        results["qc_summary"] = qc_df

    # Save outputs
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if "genes" in results:
            results["genes"].to_csv(output_dir / "gene_meth_matrix.csv")
            results["genes"].to_parquet(output_dir / "gene_meth_matrix.parquet")

        if "lncrnas" in results:
            results["lncrnas"].to_csv(output_dir / "lncrna_meth_matrix.csv")
            results["lncrnas"].to_parquet(output_dir / "lncrna_meth_matrix.parquet")

        if "ccres" in results:
            results["ccres"].to_csv(output_dir / "ccre_meth_matrix.csv")
            results["ccres"].to_parquet(output_dir / "ccre_meth_matrix.parquet")

        if "atac" in results:
            results["atac"].to_csv(output_dir / "atac_meth_matrix.csv")
            results["atac"].to_parquet(output_dir / "atac_meth_matrix.parquet")

        if "qc_summary" in results:
            results["qc_summary"].to_csv(output_dir / "sample_qc_summary.csv", index=False)

        print(f"\nSaved cohort matrices to: {output_dir}")

    return results


# =============================================================================
# MAIN PIPELINE RUNNER
# =============================================================================

def run_methylation_pipeline(
    probe_reference_path: Path = PATHS.methylation_probe_reference,
    sample_manifest_path: Path = PATHS.methylation_sample_manifest,  # columns: Sample ID, File Name
    sample_beta_dir: Path = PATHS.methylation_samples_dir,
    genes: pd.DataFrame = None,
    genes_cds: pd.DataFrame = None,
    ccres: pd.DataFrame = None,
    gene_panel: Optional[List[str]] = None,
    working_dir: Path = PATHS.methylation_output_dir,
    lncrnas: Optional[pd.DataFrame] = None,
    lncrna_panel: Optional[List[str]] = None,
    atac_peaks: Optional[pd.DataFrame] = None,
    tad_sources_dir: Optional[Path] = PATHS.tads_processed,
    include_tads: bool = True,
    ccre_ids: Optional[List[str]] = None,
    upstream_bp: int = THRESHOLDS.meth_promoter_upstream_bp,
    downstream_bp: int = THRESHOLDS.meth_promoter_downstream_bp,
    build_reference: bool = True,
    build_per_sample: bool = True,
    build_cohort: bool = True,
) -> Dict[str, Any]:
    """
    Run the complete methylation pipeline.
    
    This is the main entry point that orchestrates:
        1. Building the annotated probe reference
        2. Processing each sample
        3. Building cohort matrices
    
    Args:
        probe_reference_path: Path to GDC probe reference
        sample_manifest_path: Path to sample manifest file
        sample_beta_dir: Directory containing sample beta files
        genes: Gene DataFrame
        ccres: cCRE DataFrame
        gene_panel: List of panel genes
        working_dir: Base output directory
        lncrnas: Optional lncRNA DataFrame
        lncrna_panel: Optional list of lncRNAs
        atac_peaks: Optional ATAC peaks DataFrame
        tad_sources_dir: Optional Directory containing TAD sources
        ccre_ids: Optional list of cCRE IDs to include
        upstream_bp: Promoter upstream definition
        downstream_bp: Promoter downstream definition
        build_reference: If True, always rebuild the annotated table. If False, load
            ``{working_dir}/reference/probe_annotations.parquet`` when present, else
            ``probe_annotations.csv`` (legacy); if neither exists, builds once.
        build_per_sample: Whether to build per-sample tables
        build_cohort: Whether to build cohort matrices
    
    Returns:
        Dict with results including paths to outputs
    """
    print("\n" + "=" * 70)
    print("METHYLATION PIPELINE")
    print("=" * 70)
    
    working_dir = Path(working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)

    results = {"working_dir": str(working_dir)}
    panel_from_caller = gene_panel is not None

    # === STEP 0: Load data ===
    if genes is None:
        genes = pd.read_csv(PATHS.genes_only)

    genes = normalize_annotation_gene_names(genes, ("gene_name",))

    if gene_panel is None:
        gene_panel = list(PIPELINE_GENE_PANEL)
    if not panel_from_caller and "gene_name" in genes.columns:
        present = set(genes["gene_name"].astype(str))
        gene_panel = [g for g in gene_panel if g in present]
    if genes_cds is None:
        genes_cds = pd.read_csv(PATHS.genes_all_features)
        genes_cds = genes_cds[genes_cds["feature"] == "CDS"]

    if ccres is None:
        from ..regulatory_elements import load_regulatory_element_focus
        ccres = load_regulatory_element_focus(PATHS.regulatory_elements_table)
    if lncrnas is None:
        centric = PATHS.lncrnas_genes_centric
        legacy = PATHS.lncrnas_all_features
        if centric.is_file():
            lncrnas = pd.read_csv(centric, low_memory=False)
        elif legacy.is_file():
            lncrnas = pd.read_csv(legacy, low_memory=False)
            if "feature" in lncrnas.columns:
                lncrnas = lncrnas[lncrnas["feature"] == "gene"]
        else:
            raise FileNotFoundError(
                f"No lncRNA table found. Expected gene-centric {centric} "
                f"or legacy GENCODE export {legacy}."
            )
        lncrnas, _ = harmonize_chrom_column(lncrnas)
    if lncrnas is not None and "gene_name" in lncrnas.columns:
        lncrnas = normalize_annotation_gene_names(lncrnas, ("gene_name",))
    if atac_peaks is None:
        atac_peaks = pd.read_csv(PATHS.atac_peaks_mapping_csv)
        atac_peaks, _ = harmonize_chrom_column(atac_peaks)
    
    # === STEP 1: Probe Reference ===
    reference_dir = working_dir / "reference"
    reference_parquet = reference_dir / "probe_annotations.parquet"
    reference_csv = reference_dir / "probe_annotations.csv"
    existing_ref: Optional[Path] = None
    if reference_parquet.is_file():
        existing_ref = reference_parquet
    elif reference_csv.is_file():
        existing_ref = reference_csv

    if build_reference or existing_ref is None:
        print("\n" + "-" * 50)
        print("STEP 1: Building probe reference")
        print("-" * 50)
        
        probe_reference = build_probe_reference_table(
            probe_reference_path=probe_reference_path,
            genes=genes,
            genes_cds=genes_cds,
            ccres=ccres,
            lncrnas=lncrnas,
            atac_peaks=atac_peaks,
            gene_panel=gene_panel,
            lncrna_names=lncrna_panel,
            tad_sources_dir=tad_sources_dir,
            include_tads=include_tads,
            upstream_bp=upstream_bp,
            downstream_bp=downstream_bp,
            output_path=reference_parquet,
        )
        results["probe_reference_path"] = str(reference_parquet)
    else:
        print("\n" + "-" * 50)
        print("STEP 1: Loading existing probe reference")
        print("-" * 50)
        probe_reference = load_existing_methylation_probe_reference(existing_ref)
        results["probe_reference_path"] = str(existing_ref)
    
    # === STEP 2: Per-Sample Processing ===
    if build_per_sample:
        print("\n" + "-" * 50)
        print("STEP 2: Processing individual samples")
        print("-" * 50)
        
        per_sample_dir = working_dir / "per_sample"
        per_sample_dir.mkdir(parents=True, exist_ok=True)
        
        if sample_manifest_path.exists():
            sample_manifest = pd.read_csv(sample_manifest_path, sep="\t")
        else:
            raise FileNotFoundError(f"Sample manifest file not found: {sample_manifest_path}")
        
        for _, row in sample_manifest.iterrows():
            sample_id = row["Sample ID"]
            beta_path = sample_beta_dir / row["File Name"]
            
            if not beta_path.exists():
                print(f"  Warning: File not found for {sample_id}: {beta_path}")
                continue
            
            sample_out_dir = per_sample_dir / sample_id
            
            try:
                build_sample_methylation_tables(
                    sample_path=beta_path,
                    probe_reference=probe_reference,
                    sample_id=sample_id,
                    gene_panel=gene_panel,
                    lncrna_panel=lncrna_panel,
                    ccre_ids=ccre_ids,
                    output_dir=sample_out_dir,
                )
            except Exception as e:
                print(f"  Error processing {sample_id}: {e}")
            
            gc.collect()
        
        results["per_sample_dir"] = str(per_sample_dir)
    
    # === STEP 3: Cohort Matrices ===
    if build_cohort:
        print("\n" + "-" * 50)
        print("STEP 3: Building cohort matrices")
        print("-" * 50)
        
        cohort_dir = working_dir / "cohort"
        
        # Filter to existing files
        valid_samples = sample_manifest[
            sample_manifest["File Name"].apply(lambda x: (sample_beta_dir / x).exists())
        ]
        
        cohort_results = build_cohort_matrices(
            sample_paths=[sample_beta_dir / p for p in valid_samples["File Name"]],
            sample_ids=valid_samples["Sample ID"].tolist(),
            probe_reference=probe_reference,
            gene_panel=gene_panel,
            lncrna_panel=lncrna_panel,
            ccre_ids=ccre_ids,
            output_dir=cohort_dir,
        )
        
        results["cohort_dir"] = str(cohort_dir)
        results["cohort_stats"] = {
            "n_samples": len(valid_samples),
            "gene_matrix_shape": list(cohort_results.get("genes", pd.DataFrame()).shape),
            "ccre_matrix_shape": list(cohort_results.get("ccres", pd.DataFrame()).shape),
        }
    
    # === Summary ===
    print("\n" + "=" * 70)
    print("METHYLATION PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Working directory: {working_dir}")
    
    # Save pipeline summary
    summary_path = working_dir / "pipeline_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    
    return results


# =============================================================================
# INTEGRATION HELPERS
# =============================================================================

def integrate_methylation_with_element_table(
    element_table: pd.DataFrame,
    sample_ccre_meth: pd.DataFrame,
    ccre_id_col: str = "cCRE_id",
) -> pd.DataFrame:
    """
    Add sample-level cCRE methylation to the regulatory element table.
    
    This allows methylation to be combined with other cCRE evidence
    (SCREEN, ABC, HiChIP, etc.).
    
    Args:
        element_table: Regulatory element focus table
        sample_ccre_meth: cCRE methylation from aggregate_to_ccres
        ccre_id_col: Column name for cCRE ID in element table
    
    Returns:
        Element table with methylation columns added
    """
    et2, merge_key = resolve_methylation_ccre_id_column(
        element_table,
        ccre_id_col=ccre_id_col,
        encode_id_col="ENCODE_id",
    )

    # Merge
    meth_cols = [
        "cCRE_id", "ccre_beta_mean", "ccre_beta_median",
        "ccre_n_probes", "ccre_frac_hypermeth", "ccre_frac_hypometh",
        "ccre_CGI_overlap", "ccre_CGI_beta_mean",
    ]
    meth_cols = [c for c in meth_cols if c in sample_ccre_meth.columns]
    
    # Rename to avoid collision
    sample_id = sample_ccre_meth["sample_id"].iloc[0] if "sample_id" in sample_ccre_meth.columns else "sample"
    rename_map = {c: f"meth_{c}" for c in meth_cols if c != "cCRE_id"}
    
    meth_to_merge = sample_ccre_meth[meth_cols].rename(columns=rename_map)
    
    result = et2.merge(
        meth_to_merge,
        left_on=merge_key,
        right_on="cCRE_id",
        how="left",
    )
    
    # Clean up duplicate cCRE_id column if created
    if "cCRE_id_y" in result.columns:
        result = result.drop(columns=["cCRE_id_y"])
        result = result.rename(columns={"cCRE_id_x": "cCRE_id"})
    
    return result


def integrate_methylation_with_gene_table(
    gene_table: pd.DataFrame,
    sample_gene_meth: pd.DataFrame,
    gene_col: str = "gene_name",
) -> pd.DataFrame:
    """
    Add sample-level gene methylation to a gene-focused table.
    
    Args:
        gene_table: Gene-focused table
        sample_gene_meth: Gene methylation from aggregate_to_genes
        gene_col: Column name for gene in gene_table
    
    Returns:
        Gene table with methylation columns added
    """
    meth_cols = [
        "gene_name",
        "promoter_beta_mean", "promoter_beta_median", "promoter_beta_std",
        "promoter_n_probes", "promoter_frac_hypermeth", "promoter_frac_hypometh",
        "promoter_CGI_beta_mean", "promoter_shore_beta_mean",
        "gene_body_beta_mean", "gene_body_n_probes",
    ]
    meth_cols = [c for c in meth_cols if c in sample_gene_meth.columns]
    
    # Rename to add prefix
    rename_map = {c: f"meth_{c}" for c in meth_cols if c != "gene_name"}
    meth_to_merge = sample_gene_meth[meth_cols].rename(columns=rename_map)
    
    result = gene_table.merge(
        meth_to_merge,
        left_on=gene_col,
        right_on="gene_name",
        how="left",
    )
    
    # Clean up
    if "gene_name_y" in result.columns:
        result = result.drop(columns=["gene_name_y"])
        result = result.rename(columns={"gene_name_x": "gene_name"})
    
    return result
