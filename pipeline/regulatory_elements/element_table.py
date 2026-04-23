"""
Element-centric table building.

Functions to build the final regulatory element table with:
- Distance matching results
- Cell-line signals
- Gene links from various evidence sources
"""

import gc
import json
import os
from pathlib import Path
from typing import List, Dict, Any, Optional, Sequence, Tuple, Union

import pandas as pd

from .distance_matching import aggregate_genes_per_ccre, aggregate_ccres_per_gene, build_distance_matrix


# =============================================================================
# ELEMENT FOCUS TABLE
# =============================================================================

def build_element_focus_table(
    ccres: pd.DataFrame,
    pair_df: pd.DataFrame,
    tier_labels: List[str],
    cell_line_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build element-focused table combining distance matching with cCRE metadata.
    
    Args:
        ccres: cCRE DataFrame (may include cell-line signal columns)
        pair_df: Gene-cCRE pair DataFrame from distance matching
        tier_labels: Distance tier labels
        cell_line_cols: List of cell line column names to include
    
    Returns:
        Element-focused DataFrame
    """
    # Get per-element aggregation
    elem_agg = aggregate_genes_per_ccre(pair_df, tier_labels)
    
    # Determine which columns to bring from ccres
    merge_cols = ["cCRE_id", "ENCODE_id"]
    if cell_line_cols:
        merge_cols.extend([c for c in cell_line_cols if c in ccres.columns])
    
    # Merge with ccre metadata
    cols = [c for c in merge_cols if c in ccres.columns]
    key = [c for c in ["cCRE_id", "ENCODE_id"] if c in cols]

    ccre_subset = (
        ccres[cols]
        .groupby(key, as_index=False, sort=False)
        .first()
    )    
    elem_focus = elem_agg.merge(
        ccre_subset,
        on=["cCRE_id", "ENCODE_id"],
        how="left",
    )
    
    return elem_focus


# =============================================================================
# GENE SUMMARY TABLE
# =============================================================================

def build_gene_summary_table(
    pair_df: pd.DataFrame,
    tier_labels: List[str],
) -> pd.DataFrame:
    """
    Build gene-focused summary with cCRE counts and IDs per type/tier.
    """
    return aggregate_ccres_per_gene(pair_df, tier_labels)


# =============================================================================
# SAVE MATCHING OUTPUTS
# =============================================================================

def save_all_matching_outputs(
    output_dir: Path,
    pair_df: pd.DataFrame,
    elem_focus: pd.DataFrame,
    gene_summary: pd.DataFrame,
    gene_type: str = "coding",
) -> None:
    """
    Save all matching-related outputs.
    
    Args:
        output_dir: Base output directory
        pair_df: Gene-cCRE pairs
        elem_focus: Element-focused table
        gene_summary: Gene-focused summary
        gene_type: "coding" or "lncRNA" to organize output subdirectory
    """
    output_dir = Path(output_dir)
    
    if gene_type != "coding":
        output_dir = output_dir / gene_type
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pair_df.to_csv(output_dir / "gene_to_elements.csv", index=False)
    elem_focus.to_csv(output_dir / "regulatory_element_focus.csv", index=False)
    gene_summary.to_csv(output_dir / "gene_focus.csv", index=False)
    
    # Distance matrix
    dist_matrix = build_distance_matrix(pair_df)
    dist_matrix.to_csv(output_dir / "distance_matrix.csv")
    
    print(f"Saved matching outputs to {output_dir}")


# =============================================================================
# GENE LINKS COLUMN BUILDING
# =============================================================================

def initialize_gene_links_column(
    elem_focus: pd.DataFrame,
    pair_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Initialize the gene_links column with basic structure from distance matching.
    
    Each element gets:
    gene_links = {
        gene_name: {
            "dist_to_tss": ...,
            "tier": ...,
            # Evidence placeholders
            "screen_exp": {},
            "screen_comp": {},
            "ABC_enhancers": [],
            "hichip": {},
        },
        ...
    }
    """
    # Build gene_links dict per cCRE
    def build_gene_links(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        out = {}
        for _, row in df_sub.iterrows():
            out[row["gene_name"]] = {
                "dist_to_tss": int(row["dist_to_tss"]),
                "tier": str(row["tier"]),
                # Placeholders for evidence
                "screen_exp": {},
                "screen_comp": {},
                "ABC_enhancers": [],
                "hichip": {},
            }
        return out
    
    gene_links_series = (
        pair_df
        .groupby("cCRE_id", sort=False)
        .apply(build_gene_links)
        .rename("gene_links")
    )
    
    # Merge into elem_focus
    elem_focus = elem_focus.merge(
        gene_links_series.reset_index(),
        on="cCRE_id",
        how="left",
    )
    
    # Fill missing with empty dict
    elem_focus["gene_links"] = elem_focus["gene_links"].apply(
        lambda x: x if isinstance(x, dict) else {}
    )
    
    return elem_focus


def add_evidence_to_gene_links(
    elem_focus: pd.DataFrame,
    evidence_df: pd.DataFrame,
    evidence_type: str,
    ccre_col: str = "ENCODE_id",
    gene_col: str = "gene_name",
) -> pd.DataFrame:
    """
    Add evidence from an evidence DataFrame to gene_links.
    
    Args:
        elem_focus: Element focus table with gene_links column
        evidence_df: Evidence DataFrame with cCRE and gene columns plus evidence data
        evidence_type: Key to use in gene_links (e.g., "screen_exp", "ABC_enhancers")
        ccre_col: Column name for cCRE ID in evidence_df
        gene_col: Column name for gene in evidence_df
    
    Returns:
        Updated elem_focus
    """
    # Build lookup: (cCRE, gene) -> evidence
    evidence_lookup = {}
    for _, row in evidence_df.iterrows():
        key = (str(row[ccre_col]), str(row[gene_col]))
        evidence_lookup[key] = row.to_dict()
    
    def update_gene_links(row):
        gene_links = row["gene_links"]
        if not isinstance(gene_links, dict):
            return gene_links
        
        ccre_id = str(row.get("ENCODE_id", row.get("cCRE_id", "")))
        
        for gene_name, gene_data in gene_links.items():
            key = (ccre_id, gene_name)
            if key in evidence_lookup:
                evidence = evidence_lookup[key].copy()
                evidence.pop(ccre_col, None)
                evidence.pop(gene_col, None)
                gene_data[evidence_type] = evidence
        
        return gene_links
    
    elem_focus = elem_focus.copy()
    elem_focus["gene_links"] = elem_focus.apply(update_gene_links, axis=1)
    
    return elem_focus


# =============================================================================
# LARGE EXPORT (final element table) — parquet-first to cap RAM at STEP 12
# =============================================================================


_ELEM_FOCUS_NESTED_COLS = (
    "gene_links",
    "TAD_domains",
    "TAD_boundary_overlaps",
    "atac_peaks",
    "atac_peak_links",
    "chip_hits",
    "ABC_enhancers",
    "hichip",
    "screen_exp",
    "screen_comp",
)


# =============================================================================
# UNIFIED LOADER (parquet preferred, CSV fallback; JSON-decoded nested columns)
# =============================================================================


def _ensure_elem_id_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Guarantee an ``elem_id`` column for SV/spatial consumers.

    STEP 12 / element-focus tables typically use ``cCRE_id`` (and ``ENCODE_id``);
    older paths used ``elem_id`` directly.
    """
    if "elem_id" in df.columns:
        return df
    out = df.copy()
    if "cCRE_id" in out.columns:
        base = out["cCRE_id"]
        if "ENCODE_id" in out.columns:
            base = base.fillna(out["ENCODE_id"])
        out["elem_id"] = base.astype(str)
    elif "ENCODE_id" in out.columns:
        out["elem_id"] = out["ENCODE_id"].astype(str)
    else:
        raise ValueError(
            "Regulatory element table needs elem_id, cCRE_id, or ENCODE_id; "
            f"found columns: {list(out.columns)}"
        )
    return out


def _json_decode_nested_inplace(df: pd.DataFrame) -> None:
    """Decode JSON strings in known-nested columns back into dict/list in-place.

    Accepts both the historical ``object`` dtype and pandas ``string[python]`` /
    ``string[pyarrow]`` dtypes that parquet readers now return, since the
    on-disk values are still JSON-encoded strings either way.
    """
    for col in _ELEM_FOCUS_NESTED_COLS:
        if col not in df.columns:
            continue
        ser = df[col]
        # Skip true numeric/bool columns; accept object and pandas string dtypes.
        if not (
            pd.api.types.is_object_dtype(ser) or pd.api.types.is_string_dtype(ser)
        ):
            continue
        sample = ser.dropna().head(1)
        if sample.empty:
            continue
        v0 = sample.iloc[0]
        if isinstance(v0, str):
            first = v0.lstrip()
            if first.startswith("{") or first.startswith("["):
                df[col] = [
                    json.loads(x) if isinstance(x, str) and x else x for x in ser.tolist()
                ]


# Module-level cache: within a single Python process, SNV/SV/CNV (and sometimes Methylation)
# all call ``load_regulatory_element_focus`` with the same path. Parsing the table can be
# tens of seconds for evidence-enriched tables, so we cache the decoded DataFrame keyed by
# ``(resolved_file, mtime_ns, size, columns, decode_nested)`` and return a shallow copy.
_ELEM_FOCUS_CACHE: Dict[Tuple[Any, ...], pd.DataFrame] = {}


def clear_regulatory_element_focus_cache() -> None:
    """Drop the in-process cache used by ``load_regulatory_element_focus``."""
    _ELEM_FOCUS_CACHE.clear()


def load_regulatory_element_focus(
    path: Optional[Union[str, Path]] = None,
    columns: Optional[Sequence[str]] = None,
    decode_nested: bool = True,
    prefer_evidence: bool = True,
    use_cache: bool = True,
) -> pd.DataFrame:
    """
    Load the regulatory element focus table.

    Resolution order (when *path* is ``None`` or its sibling exists):
      1. ``<stem>_with_evidence.parquet``  (STEP 12 primary output)
      2. ``<stem>_with_evidence.csv``      (legacy / opt-in)
      3. *path* itself (e.g. pre-STEP 12 ``regulatory_element_focus.csv`` from STEP 3)

    Nested JSON-encoded columns (``gene_links``, ``TAD_domains``, ``atac_peaks``, …)
    are decoded back into ``dict``/``list`` when *decode_nested* is True.

    When *use_cache* is True (default) the decoded DataFrame is cached per (file, mtime,
    size, columns, decode_nested) so repeated calls in the same process return a shallow
    copy rather than re-parsing the file. Call ``clear_regulatory_element_focus_cache``
    to invalidate.
    """
    candidates: List[Path] = []

    if path is not None:
        p = Path(path)
        candidates.append(p.with_suffix(".parquet"))
        if prefer_evidence:
            stem = p.stem
            if "_with_evidence" not in stem:
                candidates.append(
                    p.with_name(f"{stem}_with_evidence").with_suffix(".parquet")
                )
                candidates.append(
                    p.with_name(f"{stem}_with_evidence").with_suffix(".csv")
                )
        candidates.append(p)

    chosen: Optional[Path] = None
    for c in candidates:
        if c and c.is_file():
            chosen = c
            break

    if chosen is None:
        raise FileNotFoundError(
            f"regulatory_element_focus not found. Tried: {[str(c) for c in candidates]}"
        )

    cols = list(columns) if columns is not None else None

    cache_key: Optional[Tuple[Any, ...]] = None
    if use_cache:
        try:
            st = chosen.stat()
            cache_key = (
                str(chosen.resolve()),
                st.st_mtime_ns,
                st.st_size,
                tuple(cols) if cols is not None else None,
                bool(decode_nested),
            )
        except OSError:
            cache_key = None
        if cache_key is not None:
            cached = _ELEM_FOCUS_CACHE.get(cache_key)
            if cached is not None:
                return cached.copy(deep=False)

    suffix = chosen.suffix.lower()
    if suffix == ".parquet":
        df = pd.read_parquet(chosen, columns=cols)
    else:
        df = pd.read_csv(chosen, usecols=cols, low_memory=False)

    if decode_nested:
        _json_decode_nested_inplace(df)

    df = _ensure_elem_id_column(df)

    if cache_key is not None:
        _ELEM_FOCUS_CACHE[cache_key] = df

    return df.copy(deep=False) if cache_key is not None else df


def _elem_focus_jsonify_nested_inplace(chunk: pd.DataFrame) -> None:
    """Mutate *chunk* so dict/list cells in known-nested columns become JSON strings."""
    jd = json.dumps
    for col in chunk.columns:
        ser = chunk[col]
        if ser.dtype != object:
            continue
        head = ser.dropna()
        if head.empty:
            continue
        v0 = head.iloc[0]
        if not isinstance(v0, (dict, list)):
            if col not in _ELEM_FOCUS_NESTED_COLS:
                continue
        chunk[col] = [
            jd(x, default=str) if isinstance(x, (dict, list)) else x
            for x in ser.tolist()
        ]


def _elem_focus_chunk_jsonify_objects(chunk: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of *chunk* with dict/list cells JSON-encoded for CSV."""
    out = chunk.copy()
    _elem_focus_jsonify_nested_inplace(out)
    return out


def save_regulatory_element_focus_evidence_parquet(
    elem_focus: pd.DataFrame,
    path: Union[str, Path],
    chunksize: Optional[int] = None,
    verbose: bool = True,
) -> Path:
    """
    Write ``elem_focus`` to parquet in row chunks with nested columns JSON-encoded.

    Designed for STEP 12: a single ``DataFrame.to_parquet`` on 70k rows × large nested
    ``gene_links`` / ``TAD_domains`` / ``atac_peaks`` can double RSS and OOM. This path
    keeps peak RAM proportional to *chunksize*.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    # Add derived scalar columns for scanning (pushdown-friendly).
    write_scan = os.environ.get("APM_WRITE_SCANNING_COLUMNS", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    )
    if write_scan:
        try:
            from pipeline.scanning_columns import derive_elem_focus_scanning_columns

            elem_focus = derive_elem_focus_scanning_columns(elem_focus)
        except Exception as e:
            raise RuntimeError(f"Failed to derive scanning columns for elem_focus parquet: {e}") from e

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    n = len(elem_focus)
    if n == 0:
        elem_focus.to_parquet(path, index=False)
        if verbose:
            print(f"  Saved empty table: {path}")
        return path

    raw = os.environ.get("APM_ELEM_FOCUS_PARQUET_CHUNK_ROWS", "").strip()
    if chunksize is None:
        cs = int(raw) if raw.isdigit() else 512
    else:
        cs = int(chunksize)
    cs = max(64, min(cs, n))

    if verbose:
        print(
            f"  Writing regulatory element table to parquet in row chunks of {cs} rows "
            f"({n} rows; nested dict/list columns JSON-encoded per chunk)..."
        )

    writer: Optional[pq.ParquetWriter] = None
    try:
        for start in range(0, n, cs):
            sub = elem_focus.iloc[start : start + cs].copy()
            _elem_focus_jsonify_nested_inplace(sub)
            table = pa.Table.from_pandas(sub, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    str(path),
                    table.schema,
                    compression="zstd",
                    use_dictionary=True,
                )
            else:
                if table.schema != writer.schema:
                    table = table.cast(writer.schema, safe=False)
            writer.write_table(table)
            del sub, table
            gc.collect()
    finally:
        if writer is not None:
            writer.close()

    if verbose:
        print(f"  Saved: {path}")
    return path


def save_regulatory_element_focus_evidence_csv(
    elem_focus: pd.DataFrame,
    path: Union[str, Path],
    chunksize: Optional[int] = None,
    verbose: bool = True,
) -> None:
    """
    Write ``elem_focus`` to CSV in row chunks with nested columns JSON-encoded per chunk.

    A single ``DataFrame.to_csv`` on very wide nested ``gene_links`` / evidence columns can
    spike RAM (PyArrow/pandas string materialization) and trigger OOM during STEP 12.
    """
    write_scan = os.environ.get("APM_WRITE_SCANNING_COLUMNS", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    )
    if write_scan:
        try:
            from pipeline.scanning_columns import derive_elem_focus_scanning_columns

            elem_focus = derive_elem_focus_scanning_columns(elem_focus)
        except Exception as e:
            raise RuntimeError(f"Failed to derive scanning columns for elem_focus CSV: {e}") from e

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    n = len(elem_focus)
    if n == 0:
        elem_focus.to_csv(path, index=False)
        if verbose:
            print(f"  Saved empty table: {path}")
        return

    raw = os.environ.get("APM_ELEM_FOCUS_CSV_CHUNK_ROWS", "").strip()
    if chunksize is None:
        cs = int(raw) if raw.isdigit() else 256
    else:
        cs = int(chunksize)
    cs = max(64, min(cs, n))

    if verbose:
        print(
            f"  Writing regulatory element table in CSV chunks of {cs} rows "
            f"({n} rows; nested dict/list columns JSON-encoded per chunk)..."
        )

    for start in range(0, n, cs):
        sub = elem_focus.iloc[start : start + cs].copy()
        sub = _elem_focus_chunk_jsonify_objects(sub)
        sub.to_csv(path, mode="a", header=(start == 0), index=False)
        del sub
        gc.collect()

    if verbose:
        print(f"  Saved: {path}")


# =============================================================================
# FINALIZATION
# =============================================================================

def finalize_element_table(
    elem_focus: pd.DataFrame,
    output_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Finalize element table for output.
    
    - Ensures all expected columns are present
    - Validates gene_links structure
    - Optionally saves to file
    """
    # Validate gene_links
    if "gene_links" in elem_focus.columns:
        n_empty = elem_focus["gene_links"].apply(
            lambda x: len(x) == 0 if isinstance(x, dict) else True
        ).sum()
        n_total = len(elem_focus)
        print(f"Elements with gene links: {n_total - n_empty}/{n_total}")
    
    if output_path:
        save_regulatory_element_focus_evidence_csv(
            elem_focus, output_path, verbose=True
        )
    
    return elem_focus
