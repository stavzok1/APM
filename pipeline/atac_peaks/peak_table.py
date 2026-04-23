"""
ATAC peak table building and output.

Main entry point: build_atac_peak_table()

Orchestrates:
- Peak loading and ID generation
- Gene matching (TSS distance + body overlap)
- cCRE matching (overlap + proximity)
- TAD domain and boundary annotations
- Final table assembly and saving
"""

import gc
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Queue
from threading import Thread
from typing import List, Dict, Any, Optional, Sequence, Union

import pandas as pd

from ..config import THRESHOLDS, PATHS, OUTPUT_SUBDIRS
from ..genes.symbol_normalization import normalize_annotation_gene_names

from .peak_loader import (
    load_atac_peaks,
    create_peak_id_mapping,
    filter_peaks_by_chromosomes,
)
from .gene_matching import (
    match_peaks_to_genes,
    build_gene_links,
    aggregate_genes_per_peak,
)
from .ccre_matching import (
    match_peaks_to_ccres,
    build_ccre_links,
    aggregate_ccres_per_peak,
)
from .tad_annotation import (
    annotate_peaks_with_all_tad_sources,
)


def _pyarrow_from_pandas_kwargs() -> Dict[str, Any]:
    """Optional ``n_threads`` for ``pa.Table.from_pandas`` (env ``APM_ATAC_PARQUET_N_THREADS``)."""
    raw = os.environ.get("APM_ATAC_PARQUET_N_THREADS", "").strip()
    if raw.isdigit():
        return {"n_threads": max(1, int(raw))}
    return {}


def _table_from_pandas(pa_module: Any, df: pd.DataFrame) -> Any:
    kw = _pyarrow_from_pandas_kwargs()
    try:
        return pa_module.Table.from_pandas(df, preserve_index=False, **kw)
    except TypeError:
        return pa_module.Table.from_pandas(df, preserve_index=False)


def _atac_stream_workers() -> int:
    """Worker threads for parallel JSON-prep of ATAC write-chunks.

    Env ``APM_ATAC_STREAM_WORKERS`` overrides; default ``min(4, cpu_count - 1)``.
    """
    raw = os.environ.get("APM_ATAC_STREAM_WORKERS", "").strip()
    if raw.isdigit():
        return max(1, int(raw))
    return max(1, min(4, (os.cpu_count() or 2) - 1))


def _effective_atac_tad_chunk_rows(explicit: Optional[int]) -> int:
    """
    Rows per chunk for Step 10 TAD. ``<= 0`` means legacy single-pass (no chunking).
    """
    if explicit is not None:
        return int(explicit)
    raw = os.environ.get("APM_ATAC_TAD_CHUNK_ROWS", "").strip()
    if raw.isdigit():
        return int(raw)
    return int(getattr(THRESHOLDS, "atac_tad_chunk_rows", 20_000))


def _effective_atac_tad_stream_write_rows(explicit: Optional[int]) -> int:
    if explicit is not None:
        return max(256, int(explicit))
    raw = os.environ.get("APM_ATAC_TAD_STREAM_WRITE_ROWS", "").strip()
    if raw.isdigit():
        return max(256, int(raw))
    return max(256, int(getattr(THRESHOLDS, "atac_tad_stream_write_rows", 2_048)))


def _annotate_peaks_with_tad_sources_maybe_chunked(
    peak_table: pd.DataFrame,
    tad_processed_dir: Path,
    tad_biosamples: Optional[List[str]],
    tad_chunk_rows: Optional[int],
    verbose: bool,
) -> pd.DataFrame:
    """
    Same semantics as ``annotate_peaks_with_all_tad_sources`` on the full table, but
    processes disjoint row slices when *chunk_rows* > 0.

    Important: we **always** use this loop when chunking is enabled (even if
    ``len(peak_table) <= chunk_rows``), so we never take the legacy "whole table at
    once" path for medium-sized tables — that path was easy to OOM-kill (~10–15+ GB)
    once every biosample had merged nested ``TAD_*`` dicts onto every peak row.
    """
    chunk_rows = _effective_atac_tad_chunk_rows(tad_chunk_rows)
    n = len(peak_table)
    if chunk_rows <= 0:
        return annotate_peaks_with_all_tad_sources(
            peak_table,
            processed_dir=tad_processed_dir,
            biosamples=tad_biosamples,
            verbose=verbose,
            emit_summary=True,
        )
    if n == 0:
        return peak_table
    if verbose:
        n_chunks = (n + chunk_rows - 1) // chunk_rows
        print(
            f"  [ATAC TAD] row-chunked mode: {chunk_rows} peaks/chunk, {n} peaks, {n_chunks} chunk(s) "
            f"(all biosamples per chunk; then concat — env APM_ATAC_TAD_CHUNK_ROWS / "
            f"THRESHOLDS.atac_tad_chunk_rows)"
        )
    parts: List[pd.DataFrame] = []
    for start in range(0, n, chunk_rows):
        sub = peak_table.iloc[start : start + chunk_rows].copy()
        sub = annotate_peaks_with_all_tad_sources(
            sub,
            processed_dir=tad_processed_dir,
            biosamples=tad_biosamples,
            verbose=verbose and start == 0,
            emit_summary=False,
        )
        if verbose and start > 0:
            end = start + len(sub)
            print(f"  [ATAC TAD] finished chunk rows [{start}, {end})")
        parts.append(sub)
        gc.collect()
    try:
        out = pd.concat(parts, ignore_index=True, copy=False)
    except TypeError:
        out = pd.concat(parts, ignore_index=True)
    parts.clear()
    gc.collect()
    if verbose:
        n_with_tad = out["TAD_domains"].apply(
            lambda x: len(x) > 0 if isinstance(x, dict) else False
        ).sum()
        n_boundary_overlap = out["TAD_boundary_overlaps"].apply(
            lambda x: any(v.get("overlaps_boundary", False) for v in x.values())
            if isinstance(x, dict) else False
        ).sum()
        print(f"  Peaks with TAD annotations: {n_with_tad}/{len(out)}")
        print(f"  Peaks overlapping boundaries: {n_boundary_overlap}/{len(out)}")
    return out


_ATAC_SLIM_RETURN_COLS: Sequence[str] = (
    "peak_id",
    "chrom",
    "start",
    "end",
    "center",
    "length",
    "original_name",
    "score",
    "annotation",
    "percentGC",
    "linked_genes",
    "linked_lncrnas",
    "n_ccres_total",
    "n_genes_total",
    "n_lncrnas_total",
    "n_ccres_overlapping",
    "n_genes_overlapping",
    "n_lncrnas_overlapping",
)

# Columns reloaded into RAM after ``tad_parquet_stream`` (nested dict/list cols stay on parquet only).
ATAC_SLIM_MEMORY_COLUMNS: Sequence[str] = _ATAC_SLIM_RETURN_COLS


def load_atac_peaks_annotated(
    path: Union[str, Path],
    columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    Load the saved ATAC peak table from parquet or CSV.

    After ``build_atac_peak_table(..., tad_parquet_stream=...)``, the pipeline keeps only
    ``ATAC_SLIM_MEMORY_COLUMNS`` in RAM; use this loader (optionally with ``columns``)
    when you need ``gene_links``, ``TAD_*``, etc. from disk.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    suf = path.suffix.lower()
    if suf == ".parquet":
        cols = list(columns) if columns is not None else None
        return pd.read_parquet(path, columns=cols)
    if suf == ".csv":
        return pd.read_csv(path, usecols=list(columns) if columns is not None else None, low_memory=False)
    raise ValueError(f"Unsupported ATAC table format: {path}")


def _annotate_peaks_with_tad_sources_stream_parquet(
    peak_table: pd.DataFrame,
    tad_processed_dir: Path,
    tad_biosamples: Optional[List[str]],
    tad_chunk_rows: Optional[int],
    parquet_path: Path,
    parquet_write_rows: Optional[int],
    verbose: bool,
) -> pd.DataFrame:
    """
    Row-chunk TAD annotation, appending each annotated slice to *parquet_path* without
    ``pd.concat`` of the full table (avoids a second in-RAM copy of every peak row).

    Returns a **slim** DataFrame (coordinates + overlap counts) read back from the
    parquet file; full nested columns exist only on disk.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    chunk_rows = _effective_atac_tad_chunk_rows(tad_chunk_rows)
    if chunk_rows <= 0:
        raise ValueError(
            "tad_parquet_stream requires positive TAD row chunking "
            "(set THRESHOLDS.atac_tad_chunk_rows or APM_ATAC_TAD_CHUNK_ROWS; not 0)."
        )

    parquet_path = Path(parquet_path)
    parquet_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(peak_table)
    if n == 0:
        if verbose:
            print("  [ATAC TAD] no peaks after assembly; skipping stream-to-parquet")
        return peak_table.copy()

    if parquet_path.exists():
        parquet_path.unlink()

    wrows = _effective_atac_tad_stream_write_rows(parquet_write_rows)
    writer: Optional[pq.ParquetWriter] = None
    n_with_tad = 0
    n_boundary_overlap = 0

    if verbose:
        print(
            f"  [ATAC TAD] streaming TAD → {parquet_path} "
            f"(tad_row_chunk={chunk_rows}, parquet_slice={wrows}; no full-table concat)"
        )

    # ------------------------------------------------------------------
    # Pipeline: a background writer thread pulls tables from a bounded queue
    # while the main thread runs TAD annotation; Arrow-table preparation for
    # each write-slice is parallelized with a small ThreadPool. Pure Python
    # JSON encoding releases the GIL enough for threads to help here.
    # ------------------------------------------------------------------
    write_q: "Queue[Optional[Any]]" = Queue(maxsize=4)
    writer_err: Dict[str, Exception] = {}

    def _writer_loop() -> None:
        nonlocal writer
        try:
            while True:
                item = write_q.get()
                if item is None:
                    return
                if writer is None:
                    writer = pq.ParquetWriter(
                        str(parquet_path),
                        item.schema,
                        compression="zstd",
                        use_dictionary=True,
                    )
                else:
                    if item.schema != writer.schema:
                        item = item.cast(writer.schema, safe=False)
                writer.write_table(item)
        except Exception as exc:  # pragma: no cover — surfaced below
            writer_err["err"] = exc

    writer_thread = Thread(target=_writer_loop, name="atac-tad-parquet-writer")
    writer_thread.start()

    def _prep_slice(w0: int, sub_df: pd.DataFrame) -> Any:
        wchunk = sub_df.iloc[w0 : w0 + wrows].copy()
        _prepare_atac_chunk_for_parquet_inplace(wchunk)
        return _table_from_pandas(pa, wchunk)

    try:
        prep_workers = _atac_stream_workers()
        prep_pool = ThreadPoolExecutor(max_workers=prep_workers) if prep_workers > 1 else None
        try:
            for start in range(0, n, chunk_rows):
                if writer_err:
                    raise writer_err["err"]
                sub = peak_table.iloc[start : start + chunk_rows].copy()
                sub = annotate_peaks_with_all_tad_sources(
                    sub,
                    processed_dir=tad_processed_dir,
                    biosamples=tad_biosamples,
                    verbose=verbose and start == 0,
                    emit_summary=False,
                )
                n_with_tad += int(
                    sub["TAD_domains"]
                    .map(lambda x: len(x) > 0 if isinstance(x, dict) else False)
                    .sum()
                )
                n_boundary_overlap += int(
                    sub["TAD_boundary_overlaps"]
                    .map(
                        lambda x: any(v.get("overlaps_boundary", False) for v in x.values())
                        if isinstance(x, dict)
                        else False
                    )
                    .sum()
                )
                if verbose and start > 0:
                    end = start + len(sub)
                    print(f"  [ATAC TAD] streamed TAD rows [{start}, {end})")

                slice_starts = list(range(0, len(sub), wrows))
                if prep_pool is not None and len(slice_starts) > 1:
                    tables = list(prep_pool.map(lambda w0: _prep_slice(w0, sub), slice_starts))
                else:
                    tables = [_prep_slice(w0, sub) for w0 in slice_starts]
                for t in tables:
                    write_q.put(t)
                del sub, tables
                gc.collect()
        finally:
            if prep_pool is not None:
                prep_pool.shutdown(wait=True)
    finally:
        write_q.put(None)
        writer_thread.join()
        if writer is not None:
            writer.close()
        if writer_err:
            raise writer_err["err"]

    schema = pq.read_schema(str(parquet_path))
    slim = [c for c in _ATAC_SLIM_RETURN_COLS if c in schema.names]
    out = pd.read_parquet(parquet_path, columns=slim if slim else None)
    gc.collect()

    if verbose:
        print(
            f"  Peaks with TAD annotations: {n_with_tad}/{n}\n"
            f"  Peaks overlapping boundaries: {n_boundary_overlap}/{n}\n"
            f"  [ATAC TAD] reloaded slim in-memory frame ({len(out)} rows, {len(slim)} cols) "
            f"for downstream steps; full table remains in {parquet_path.name}"
        )
    return out


# Columns with dict/list cell values: JSON-stringify before parquet so PyArrow does not
# materialize huge nested arrays (can trigger OOM right after the peak table is built).
_ATAC_PARQUET_JSON_COLS = (
    "gene_links",
    "genes_by_tier",
    "lncrna_links",
    "lncrnas_by_tier",
    "ccre_links",
    "ccre_types",
    "TAD_domains",
    "TAD_boundary_overlaps",
)


def _prepare_atac_chunk_for_parquet_inplace(chunk: pd.DataFrame) -> None:
    """Mutate *chunk* so nested columns are JSON strings (or CSV-style lists)."""
    jd = json.dumps

    def _json_cell(x: Any) -> Any:
        return jd(x, default=str) if isinstance(x, (dict, list)) else x

    for col in _ATAC_PARQUET_JSON_COLS:
        if col not in chunk.columns:
            continue
        chunk[col] = [_json_cell(x) for x in chunk[col].tolist()]
    if "linked_genes" in chunk.columns:
        chunk["linked_genes"] = [
            ",".join(x) if isinstance(x, list) else x for x in chunk["linked_genes"].tolist()
        ]
    if "linked_lncrnas" in chunk.columns:
        chunk["linked_lncrnas"] = [
            ",".join(x) if isinstance(x, list) else x for x in chunk["linked_lncrnas"].tolist()
        ]


def _write_atac_peak_table_parquet_chunked(
    peak_table: pd.DataFrame,
    path: Path,
    chunksize: int = 4096,
) -> None:
    """
    Write peak_table to parquet in row chunks with nested columns JSON-encoded.

    Avoids a single PyArrow conversion of very large nested object columns, which often
    spikes RAM and triggers the OOM killer after ``PEAK TABLE COMPLETE``.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    path = Path(path)
    n = len(peak_table)
    if n == 0:
        peak_table.to_parquet(path, index=False)
        return

    cs = max(1, int(chunksize))
    writer: Optional[pq.ParquetWriter] = None
    write_q: "Queue[Optional[Any]]" = Queue(maxsize=4)
    writer_err: Dict[str, Exception] = {}

    def _writer_loop() -> None:
        nonlocal writer
        try:
            while True:
                item = write_q.get()
                if item is None:
                    return
                if writer is None:
                    writer = pq.ParquetWriter(
                        str(path),
                        item.schema,
                        compression="zstd",
                        use_dictionary=True,
                    )
                else:
                    if item.schema != writer.schema:
                        item = item.cast(writer.schema, safe=False)
                writer.write_table(item)
        except Exception as exc:
            writer_err["err"] = exc

    writer_thread = Thread(target=_writer_loop, name="atac-parquet-writer")
    writer_thread.start()

    def _prep(start: int) -> Any:
        chunk = peak_table.iloc[start : start + cs].copy()
        _prepare_atac_chunk_for_parquet_inplace(chunk)
        return _table_from_pandas(pa, chunk)

    try:
        prep_workers = _atac_stream_workers()
        if prep_workers > 1 and n > cs:
            with ThreadPoolExecutor(max_workers=prep_workers) as pool:
                for tbl in pool.map(_prep, range(0, n, cs)):
                    if writer_err:
                        raise writer_err["err"]
                    write_q.put(tbl)
        else:
            for start in range(0, n, cs):
                if writer_err:
                    raise writer_err["err"]
                write_q.put(_prep(start))
    finally:
        write_q.put(None)
        writer_thread.join()
        if writer is not None:
            writer.close()
        if writer_err:
            raise writer_err["err"]


# =============================================================================
# MAIN TABLE BUILDER
# =============================================================================

def build_atac_peak_table(
    peaks: Union[str, Path, pd.DataFrame],
    genes: pd.DataFrame,
    ccres: pd.DataFrame,
    lncrnas: Optional[pd.DataFrame] = None,
    gene_panel: Optional[List[str]] = None,
    gene_window_bp: Optional[int] = None,
    lncrna_window_bp: Optional[int] = None,
    ccre_max_distance: Optional[int] = None,
    tier_edges: Optional[List[int]] = None,
    tier_labels: Optional[List[str]] = None,
    tad_processed_dir: Optional[Union[str, Path]] = None,
    tad_biosamples: Optional[List[str]] = None,
    tad_chunk_rows: Optional[int] = None,
    tad_parquet_stream: Optional[Union[str, Path]] = None,
    tad_stream_write_rows: Optional[int] = None,
    exclude_chroms: Optional[List[str]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Build complete ATAC peak annotation table.
    
    Args:
        peaks: Path to ATAC peaks file OR pre-loaded DataFrame
        genes: Gene annotations DataFrame
        ccres: cCRE annotations DataFrame
        lncrnas: Optional lncRNA annotations DataFrame
        gene_panel: Optional list of genes to filter to
        gene_window_bp: Maximum distance from TSS for gene matching
        ccre_max_distance: Maximum distance for cCRE matching (0 = overlap only)
        tier_edges: Distance tier boundaries
        tier_labels: Labels for tiers
        tad_processed_dir: Path to processed TAD data (None = skip TAD annotation)
        tad_biosamples: Specific TAD biosamples to use (None = all)
        tad_chunk_rows: Peaks per RAM-bounded TAD chunk (None → ``THRESHOLDS.atac_tad_chunk_rows`` / env).
            Use ``0`` via env ``APM_ATAC_TAD_CHUNK_ROWS=0`` for legacy single-pass.
        tad_parquet_stream: If set with ``tad_processed_dir``, TAD-annotated row chunks are appended
            directly to this parquet path (no ``pd.concat`` of the full table). The returned
            DataFrame is a **slim** subset read back from that file; nested columns exist only on disk.
        tad_stream_write_rows: Rows per PyArrow slice when streaming (None → thresholds / env).
        exclude_chroms: Chromosomes to exclude (e.g., ["chrM", "chrY"])
        verbose: Print progress
    
    Returns:
        Peak annotation table. With ``tad_parquet_stream``, nested gene/cCRE/lncRNA/TAD
        columns are written only to that parquet file; the returned frame is a slim subset
        for downstream RAM (see ``_ATAC_SLIM_RETURN_COLS``).
    """
    # Use config defaults
    gene_window_bp = gene_window_bp or THRESHOLDS.ccre_window_bp
    lncrna_window_bp = lncrna_window_bp or THRESHOLDS.ccre_window_bp
    ccre_max_distance = ccre_max_distance or THRESHOLDS.atac_ccre_max_distance
    tier_edges = tier_edges or THRESHOLDS.tier_edges_bp
    tier_labels = tier_labels or THRESHOLDS.tier_labels
    
    if verbose:
        print("\n" + "=" * 60)
        print("BUILDING ATAC PEAK TABLE")
        print("=" * 60)
    
    # =========================================================================
    # STEP 1: Load peaks
    # =========================================================================
    if verbose:
        print("\n" + "-" * 40)
        print("Loading ATAC peaks")
        print("-" * 40)
    
    if isinstance(peaks, (str, Path)):
        peaks = load_atac_peaks(peaks)
    else:
        peaks = peaks.copy()
        # Ensure peak_id exists
        if "peak_id" not in peaks.columns:
            from .peak_loader import generate_peak_ids
            peaks["peak_id"] = generate_peak_ids(peaks)
    
    # Filter chromosomes if specified
    if exclude_chroms:
        peaks = filter_peaks_by_chromosomes(peaks, exclude_chroms=exclude_chroms)
        if verbose:
            print(f"  After filtering: {len(peaks)} peaks")
    
    # =========================================================================
    # STEP 2: Canonicalize symbols, then filter genes to panel (if specified)
    # =========================================================================
    genes_work = normalize_annotation_gene_names(genes.copy(), ("gene_name",))
    if gene_panel:
        genes_filtered = genes_work[genes_work["gene_name"].isin(gene_panel)].copy()
        if verbose:
            print(f"\n  Filtered to {len(genes_filtered)} genes in panel")
    else:
        genes_filtered = genes_work

    if lncrnas is not None:
        lncrnas = normalize_annotation_gene_names(lncrnas.copy(), ("gene_name",))
    
    # =========================================================================
    # STEP 3: Gene matching
    # =========================================================================
    if verbose:
        print("\n" + "-" * 40)
        print("Matching peaks to genes")
        print("-" * 40)
    
    peak_gene_pairs = match_peaks_to_genes(
        peaks,
        genes_filtered,
        window_bp=gene_window_bp,
        tier_edges=tier_edges,
        tier_labels=tier_labels,
        promoter_upstream=THRESHOLDS.promoter_upstream_bp,
        promoter_downstream=THRESHOLDS.promoter_downstream_bp,
    )
    
    # Build gene_links dict
    gene_links_df = build_gene_links(peak_gene_pairs)
    
    # Build aggregated gene info
    gene_agg = aggregate_genes_per_peak(peak_gene_pairs, tier_labels)
    
    if lncrnas is not None:
        if verbose:
            print("\n" + "-" * 40)
            print("Matching peaks to lncRNAs")
            print("-" * 40)
        
        peak_lncrna_pairs = match_peaks_to_genes(
            peaks,
            lncrnas,
            window_bp=lncrna_window_bp,
            tier_edges=tier_edges,
            tier_labels=tier_labels,
            promoter_upstream=THRESHOLDS.promoter_upstream_bp,
            promoter_downstream=THRESHOLDS.promoter_downstream_bp,
        )
        
        # Build lncrna_links dict
        lncrna_links_df = build_gene_links(peak_lncrna_pairs)
        lncrna_links_df.rename(columns={"gene_links": "lncrna_links"}, inplace=True)
        
        # Build aggregated lncrna info
        lncrna_agg = aggregate_genes_per_peak(peak_lncrna_pairs, tier_labels)
        lncrna_agg.rename(columns={"linked_genes": "linked_lncrnas", "genes_by_tier": "lncrnas_by_tier", "n_genes_total": "n_lncrnas_total", "n_genes_overlapping": "n_lncrnas_overlapping"}, inplace=True)
    else:
        lncrna_links_df = pd.DataFrame()
        lncrna_agg = pd.DataFrame()


    # =========================================================================
    # STEP 4: cCRE matching
    # =========================================================================
    if verbose:
        print("\n" + "-" * 40)
        print("Matching peaks to cCREs")
        print("-" * 40)
    
    peak_ccre_pairs = match_peaks_to_ccres(
        peaks,
        ccres,
        max_distance=ccre_max_distance,
    )
    
    # Build ccre_links list
    ccre_links_df = build_ccre_links(peak_ccre_pairs)
    
    # Build aggregated cCRE info
    ccre_agg = aggregate_ccres_per_peak(peak_ccre_pairs)
    
    # =========================================================================
    # STEP 5: Assemble peak table
    # =========================================================================
    if verbose:
        print("\n" + "-" * 40)
        print("Assembling peak table")
        print("-" * 40)
    
    # Start with base peaks
    peak_table = peaks.copy()
    
    # Merge gene links
    peak_table = peak_table.merge(
        gene_links_df,
        on="peak_id",
        how="left",
    )
    peak_table["gene_links"] = peak_table["gene_links"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )

    # Merge lncrna links
    if not lncrna_links_df.empty:
        peak_table = peak_table.merge(
            lncrna_links_df,
            on="peak_id",
            how="left",
        )
        peak_table["lncrna_links"] = peak_table["lncrna_links"].apply(
            lambda x: x if isinstance(x, dict) else {}
        )
    else:
        peak_table["lncrna_links"] = [{}] * len(peak_table)
    
    # Merge gene aggregates
    peak_table = peak_table.merge(
        gene_agg[["peak_id", "linked_genes", "genes_by_tier", "n_genes_total", "n_genes_overlapping"]],
        on="peak_id",
        how="left",
    )

    # Merge lncrna aggregates
    if not lncrna_agg.empty:
        peak_table = peak_table.merge(
            lncrna_agg[["peak_id", "linked_lncrnas", "lncrnas_by_tier", "n_lncrnas_total", "n_lncrnas_overlapping"]],
            on="peak_id",
            how="left",
        )
    else:
        peak_table["linked_lncrnas"] = [{}] * len(peak_table)
        peak_table["lncrnas_by_tier"] = [{}] * len(peak_table)
        peak_table["n_lncrnas_total"] = 0
        peak_table["n_lncrnas_overlapping"] = 0
    
    # Fill missing gene aggregates
    peak_table["linked_genes"] = peak_table["linked_genes"].apply(
        lambda x: x if isinstance(x, list) else []
    )
    peak_table["genes_by_tier"] = peak_table["genes_by_tier"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )
    peak_table["n_genes_total"] = peak_table["n_genes_total"].fillna(0).astype(int)
    peak_table["n_genes_overlapping"] = peak_table["n_genes_overlapping"].fillna(0).astype(int)

    # Fill missing lncrna aggregates
    peak_table["linked_lncrnas"] = peak_table["linked_lncrnas"].apply(
        lambda x: x if isinstance(x, list) else []
    )
    peak_table["lncrnas_by_tier"] = peak_table["lncrnas_by_tier"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )
    peak_table["n_lncrnas_total"] = peak_table["n_lncrnas_total"].fillna(0).astype(int)
    peak_table["n_lncrnas_overlapping"] = peak_table["n_lncrnas_overlapping"].fillna(0).astype(int)
    
    # Merge cCRE links
    peak_table = peak_table.merge(
        ccre_links_df,
        on="peak_id",
        how="left",
    )
    peak_table["ccre_links"] = peak_table["ccre_links"].apply(
        lambda x: x if isinstance(x, list) else []
    )
    
    # Merge cCRE aggregates
    if not ccre_agg.empty:
        peak_table = peak_table.merge(
            ccre_agg[["peak_id", "n_ccres_total", "n_ccres_overlapping", "ccre_types"]],
            on="peak_id",
            how="left",
        )
        peak_table["n_ccres_total"] = peak_table["n_ccres_total"].fillna(0).astype(int)
        peak_table["n_ccres_overlapping"] = peak_table["n_ccres_overlapping"].fillna(0).astype(int)
        peak_table["ccre_types"] = peak_table["ccre_types"].apply(
            lambda x: x if isinstance(x, dict) else {}
        )
    else:
        peak_table["n_ccres_total"] = 0
        peak_table["n_ccres_overlapping"] = 0
        peak_table["ccre_types"] = [{}] * len(peak_table)
    
    if verbose:
        print(f"  Assembled table with {len(peak_table)} peaks")
        print(f"  Peaks with gene links: {(peak_table['n_genes_total'] > 0).sum()}")
        print(f"  Peaks with cCRE overlaps: {(peak_table['n_ccres_overlapping'] > 0).sum()}")
    
    peak_table = peak_table[(peak_table["n_ccres_total"] > 0) | (peak_table["n_genes_total"] > 0) | (peak_table["n_lncrnas_total"] > 0)]
    # =========================================================================
    # STEP 6: TAD annotation (if requested)
    # =========================================================================
    if tad_processed_dir:
        if verbose:
            print("\n" + "-" * 40)
            print("Adding TAD annotations")
            print("-" * 40)
        
        tad_processed_dir = Path(tad_processed_dir)
        stream_path = Path(tad_parquet_stream) if tad_parquet_stream else None
        if stream_path is not None:
            stream_path.parent.mkdir(parents=True, exist_ok=True)
            if stream_path.exists():
                stream_path.unlink()

        if stream_path is not None and len(peak_table) > 0:
            peak_table = _annotate_peaks_with_tad_sources_stream_parquet(
                peak_table,
                tad_processed_dir,
                tad_biosamples,
                tad_chunk_rows,
                stream_path,
                tad_stream_write_rows,
                verbose,
            )
        else:
            peak_table = _annotate_peaks_with_tad_sources_maybe_chunked(
                peak_table,
                tad_processed_dir,
                tad_biosamples,
                tad_chunk_rows,
                verbose,
            )
    
    # =========================================================================
    # Final column ordering
    # =========================================================================
    core_cols = [
        "peak_id", "chrom", "start", "end", "center", "length",
        "original_name", "score", "annotation", "percentGC",
    ]
    
    gene_cols = [
        "linked_genes", "genes_by_tier", "n_genes_total", "n_genes_overlapping",
        "gene_links",
    ]
    
    lncrna_cols = [
        "linked_lncrnas", "lncrnas_by_tier", "n_lncrnas_total", "n_lncrnas_overlapping",
        "lncrna_links",
    ]
    
    ccre_cols = [
        "n_ccres_total", "n_ccres_overlapping", "ccre_types",
        "ccre_links",
    ]
    
    tad_cols = ["TAD_domains", "TAD_boundary_overlaps"]
    
    # Build final column order
    final_cols = []
    for c in core_cols + gene_cols + lncrna_cols + ccre_cols + tad_cols:
        if c in peak_table.columns:
            final_cols.append(c)
    
    # Add any remaining columns
    remaining = [c for c in peak_table.columns if c not in final_cols]
    final_cols.extend(remaining)
    
    peak_table = peak_table[final_cols]

    gc.collect()
    if verbose:
        print("\n" + "=" * 60)
        print("PEAK TABLE COMPLETE")
        print("=" * 60)
        print(f"  Total peaks: {len(peak_table)}")
        print(f"  Columns: {len(peak_table.columns)}")
    
    return peak_table


# =============================================================================
# SAVE OUTPUTS
# =============================================================================

def save_atac_outputs(
    peak_table: pd.DataFrame,
    output_dir: Union[str, Path],
    peak_gene_pairs: Optional[pd.DataFrame] = None,
    peak_lncrna_pairs: Optional[pd.DataFrame] = None,
    peak_ccre_pairs: Optional[pd.DataFrame] = None,
    save_parquet: bool = False,
    save_csv: bool = True,
    parquet_chunksize: int = 2048,
    reuse_main_parquet_if: Optional[Union[str, Path]] = None,
) -> Dict[str, Path]:
    """
    Save ATAC peak table and related outputs.
    
    Args:
        peak_table: Main peak annotation table
        output_dir: Output directory
        peak_gene_pairs: Optional gene matching pairs DataFrame
        peak_lncrna_pairs: Optional lncRNA matching pairs DataFrame
        peak_ccre_pairs: Optional cCRE matching pairs DataFrame
        save_parquet: Save parquet format. Large tables are written in row chunks; nested
            dict/list columns are stored as JSON strings to limit peak RAM during write.
        save_csv: Save CSV format (nested dicts as JSON strings)
        parquet_chunksize: Row batch size when ``save_parquet`` is True (ignored for CSV).
        reuse_main_parquet_if: If this path resolves to the same file as ``atac_peaks_annotated.parquet``
            in *output_dir* and the file already exists (e.g. written during TAD streaming), skip rewriting.
    
    Returns:
        Dict mapping output names to file paths
    """
    # Add derived scalar columns for scanning (pushdown-friendly).
    write_scan = os.environ.get("APM_WRITE_SCANNING_COLUMNS", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    )
    if write_scan:
        try:
            from pipeline.scanning_columns import derive_atac_peak_scanning_columns

            peak_table = derive_atac_peak_scanning_columns(peak_table)
        except Exception as e:
            raise RuntimeError(f"Failed to derive scanning columns for ATAC peaks: {e}") from e

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = {}
    
    # Main peak table
    if save_parquet:
        path = output_dir / "atac_peaks_annotated.parquet"
        n_peaks = len(peak_table)
        reuse = (
            Path(reuse_main_parquet_if).resolve()
            if reuse_main_parquet_if
            else None
        )
        if reuse is not None and reuse == path.resolve() and path.exists():
            saved_files["peak_table_parquet"] = path
            print(f"ATAC main parquet already on disk (TAD stream): {path}")
        elif n_peaks == 0:
            peak_table.to_parquet(path, index=False)
            saved_files["peak_table_parquet"] = path
            print(f"Saved: {path}")
        else:
            # Never use ``peak_table.copy(); to_parquet()`` for large tables: it briefly
            # **doubles** RSS and pushes PyArrow to convert the entire nested object graph at
            # once — a common OOM right after Step 10 TAD (RAM can sit ~80% then jump as the
            # copy + Arrow build land). Always stream row chunks with JSON-encoded nested cols.
            cs = max(256, min(int(parquet_chunksize), n_peaks))
            print(
                f"Writing {n_peaks} ATAC peaks to parquet in row chunks of {cs} "
                f"(nested columns JSON-encoded per chunk to cap PyArrow RAM)..."
            )
            gc.collect()
            _write_atac_peak_table_parquet_chunked(peak_table, path, chunksize=cs)
            saved_files["peak_table_parquet"] = path
            print(f"Saved: {path}")
    
    if save_csv:
        # Convert nested structures to JSON for CSV
        csv_table = peak_table.copy()
        for col in ["gene_links", "genes_by_tier", "lncrna_links", "lncrnas_by_tier", "ccre_links", "ccre_types", 
                    "TAD_domains", "TAD_boundary_overlaps"]:
            if col in csv_table.columns:
                csv_table[col] = csv_table[col].apply(
                    lambda x: json.dumps(x) if isinstance(x, (dict, list)) else x
                )
        
        # Convert linked_genes list to comma-separated
        if "linked_genes" in csv_table.columns:
            csv_table["linked_genes"] = csv_table["linked_genes"].apply(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )
        if "linked_lncrnas" in csv_table.columns:
            csv_table["linked_lncrnas"] = csv_table["linked_lncrnas"].apply(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )
        
        path = output_dir / "atac_peaks_annotated.csv"
        csv_table.to_csv(path, index=False)
        saved_files["peak_table_csv"] = path
        print(f"Saved: {path}")
    
    # Peak ID mapping
    mapping = peak_table[["peak_id", "original_name", "chrom", "start", "end"]].copy()
    path = output_dir / "peak_id_mapping.csv"
    mapping.to_csv(path, index=False)
    saved_files["peak_id_mapping"] = path
    print(f"Saved: {path}")
    
    # Optional pair tables
    if peak_gene_pairs is not None:
        pairs_csv = peak_gene_pairs.copy()
        if "body_overlap" in pairs_csv.columns:
            pairs_csv["body_overlap"] = pairs_csv["body_overlap"].apply(json.dumps)
        
        path = output_dir / "peak_gene_pairs.csv"
        pairs_csv.to_csv(path, index=False)
        saved_files["peak_gene_pairs"] = path
        print(f"Saved: {path}")

    if peak_lncrna_pairs is not None:
        pairs_csv = peak_lncrna_pairs.copy()
        if "body_overlap" in pairs_csv.columns:
            pairs_csv["body_overlap"] = pairs_csv["body_overlap"].apply(json.dumps)
        
        path = output_dir / "peak_lncrna_pairs.csv"
        pairs_csv.to_csv(path, index=False)
        saved_files["peak_lncrna_pairs"] = path
        print(f"Saved: {path}")
    
    if peak_ccre_pairs is not None:
        pairs_csv = peak_ccre_pairs.copy()
        if "overlap" in pairs_csv.columns:
            pairs_csv["overlap"] = pairs_csv["overlap"].apply(json.dumps)
        
        path = output_dir / "peak_ccre_pairs.csv"
        pairs_csv.to_csv(path, index=False)
        saved_files["peak_ccre_pairs"] = path
        print(f"Saved: {path}")
    
    # Summary statistics
    summary = {
        "total_peaks": len(peak_table),
        "peaks_with_gene_links": int((peak_table["n_genes_total"] > 0).sum()),
        "peaks_with_gene_overlap": int((peak_table["n_genes_overlapping"] > 0).sum()),
        "peaks_with_ccre_overlap": int((peak_table["n_ccres_overlapping"] > 0).sum()),
    }
    
    if "TAD_domains" in peak_table.columns:
        summary["peaks_with_tad_annotation"] = int(peak_table["TAD_domains"].apply(
            lambda x: len(x) > 0 if isinstance(x, dict) else False
        ).sum())
    
    if "TAD_boundary_overlaps" in peak_table.columns:
        summary["peaks_at_tad_boundaries"] = int(peak_table["TAD_boundary_overlaps"].apply(
            lambda x: any(v.get("overlaps_boundary", False) for v in x.values())
            if isinstance(x, dict) else False
        ).sum())
    
    path = output_dir / "summary.json"
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
    saved_files["summary"] = path
    print(f"Saved: {path}")
    
    return saved_files


# =============================================================================
# CONVENIENCE: QUICK ANALYSIS FUNCTIONS
# =============================================================================

def get_peaks_near_gene(
    peak_table: pd.DataFrame,
    gene_name: str,
    max_distance: Optional[int] = None,
) -> pd.DataFrame:
    """
    Get all peaks linked to a specific gene (requires ``gene_links`` on *peak_table*).

    After TAD streaming, load from parquet with ``load_atac_peaks_annotated(path, columns=[...])``.
    """
    mask = peak_table["gene_links"].apply(
        lambda x: gene_name in x if isinstance(x, dict) else False
    )
    
    result = peak_table[mask].copy()
    
    if max_distance is not None:
        result = result[result["gene_links"].apply(
            lambda x: x.get(gene_name, {}).get("dist_to_tss", float("inf")) <= max_distance
        )]
    
    return result


def get_peaks_overlapping_gene(
    peak_table: pd.DataFrame,
    gene_name: str,
) -> pd.DataFrame:
    """Get peaks that overlap a specific gene's body (needs ``gene_links``; see ``load_atac_peaks_annotated``)."""
    mask = peak_table["gene_links"].apply(
        lambda x: x.get(gene_name, {}).get("body_overlap", {}).get("overlaps", False)
        if isinstance(x, dict) else False
    )
    return peak_table[mask].copy()


def get_peaks_at_promoters(
    peak_table: pd.DataFrame,
    gene_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Get peaks that overlap gene promoters (needs ``gene_links``; see ``load_atac_peaks_annotated``)."""
    def check_promoter_overlap(x):
        if not isinstance(x, dict):
            return False
        for gene, data in x.items():
            if gene_names and gene not in gene_names:
                continue
            if data.get("body_overlap", {}).get("overlap_type") == "promoter":
                return True
        return False
    
    mask = peak_table["gene_links"].apply(check_promoter_overlap)
    return peak_table[mask].copy()
