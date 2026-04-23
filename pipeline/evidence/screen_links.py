"""
ENCODE SCREEN enhancer-gene link processing.

Handles both:
- Experimental (3D chromatin) links
- Computational (ABC, etc.) links

Complete pipeline: load → filter → aggregate → strength classify → build nested dicts
"""

import os
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Set, Dict, Any, Optional, Tuple, Iterable

import numpy as np
import pandas as pd


def _screen_parallel_workers(n_tasks: int) -> int:
    """Worker count for per-biosample / per-conservation SCREEN aggregation.

    Env ``APM_SCREEN_WORKERS`` overrides; default is ``os.cpu_count() - 1`` capped at *n_tasks*.
    """
    raw = os.environ.get("APM_SCREEN_WORKERS", "").strip()
    if raw.isdigit():
        n = int(raw)
    else:
        n = max(1, (os.cpu_count() or 2) - 1)
    return max(1, min(n, max(1, n_tasks)))

from ..utils import classify_strength, safe_float
from ..biosample_names import canonicalize_screen_biosample
from ..schemas import (
    empty_assay_entry,
    empty_biosample_assays,
    empty_conservation_entry,
    empty_conservation_block,
    empty_screen_block,
    ensure_screen_block,
)


# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

COLS_SCREEN_EXP = [
    "cCRE ID", "Gene ID", "Common Gene Name", "Gene Type",
    "Assay Type", "Experiment ID", "Biosample", "Score", "P-value",
]

COLS_SCREEN_COMP = [
    "cCRE ID", "Gene ID", "Common Gene Name", "Gene Type",
    "Assay Type", "Region", "Experiment ID", "Biosample", "Score",
]

RENAME_SCREEN = {
    "cCRE ID": "ENCODE_id",
    "Gene ID": "gene_id",
    "Common Gene Name": "gene_name",
    "Gene Type": "gene_type",
    "Assay Type": "assay_type",
    "Experiment ID": "experiment_id",
    "Biosample": "biosample",
    "Score": "score",
    "P-value": "p_value",
    "Region": "region",
}

LINK_COLS = ["ENCODE_id", "gene_id", "gene_name"]


# =============================================================================
# STREAMING & FILTERING
# =============================================================================

def _stream_and_filter_genes(
    file_handle,
    gene_set: Set[str],
    columns: List[str],
    chunksize: int,
    sep: str = "\t",
    compression: Optional[str] = None,
    add_dummy_pvalue: bool = False,
) -> pd.DataFrame:
    """Stream a tabular file and retain only rows matching the gene set."""
    chunks = []
    reader = pd.read_csv(
        file_handle,
        sep=sep,
        header=None,
        names=columns,
        chunksize=chunksize,
        low_memory=False,
        compression=compression,
    )

    for i, chunk in enumerate(reader):
        chunk["Common Gene Name"] = chunk["Common Gene Name"].astype(str).str.strip()
        chunk["Score"] = pd.to_numeric(chunk["Score"], errors="coerce")

        
        if add_dummy_pvalue:
            chunk["P-value"] = np.nan

        mask = chunk["Common Gene Name"].isin(gene_set)
        if mask.any():
            chunks.append(chunk.loc[mask])

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1} chunks...")

    if not chunks:
        raise ValueError("No rows found for selected genes.")

    try:
        df = pd.concat(chunks, ignore_index=True, copy=False)
    except TypeError:
        df = pd.concat(chunks, ignore_index=True)
    print(f"  Total rows extracted: {len(df)}")
    return df


# =============================================================================
# STRENGTH CLASSIFICATION
# =============================================================================

def _apply_strength_classification(
    df: pd.DataFrame,
    weak_quantile: float,
    strong_quantile: float,
    intact_pvalue_threshold: Optional[float] = None,
    quantiles_out_path: Optional[Path] = None,
    provenance: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """Compute per-assay quantiles and assign strength labels."""
    quantiles = (
        df.groupby("assay_type")["score"]
        .quantile([weak_quantile, strong_quantile])
        .unstack()
        .rename(columns={weak_quantile: "q_weak", strong_quantile: "q_strong"})
    )
    print("  Per-assay score quantiles:")
    for line in quantiles.to_string().split('\n'):
        print(f"    {line}")

    df = df.merge(quantiles, left_on="assay_type", right_index=True, how="left")

    df["is_strong"] = df["score"] >= df["q_strong"]
    df["is_weak"] = (df["score"] >= df["q_weak"]) & ~df["is_strong"]

    if intact_pvalue_threshold is not None:
        mask = df["assay_type"] == "Intact-HiC"
        df.loc[mask, "is_strong"] &= df.loc[mask, "p_value"] <= intact_pvalue_threshold

    df["strength"] = [classify_strength(s, w) for s, w in zip(df["is_strong"], df["is_weak"])]

    if quantiles_out_path is not None:
        try:
            quantiles_out_path = Path(quantiles_out_path)
            quantiles_out_path.parent.mkdir(parents=True, exist_ok=True)
            q = quantiles.reset_index().copy()
            q["weak_quantile"] = float(weak_quantile)
            q["strong_quantile"] = float(strong_quantile)
            q["intact_pvalue_threshold"] = float(intact_pvalue_threshold) if intact_pvalue_threshold is not None else np.nan
            q.to_csv(quantiles_out_path, index=False)

            if provenance is not None:
                import json

                meta_path = quantiles_out_path.with_suffix(quantiles_out_path.suffix + ".meta.json")
                meta = dict(provenance)
                meta.update(
                    {
                        "weak_quantile": float(weak_quantile),
                        "strong_quantile": float(strong_quantile),
                        "intact_pvalue_threshold": float(intact_pvalue_threshold) if intact_pvalue_threshold is not None else None,
                        "n_rows_scored": int(len(df)),
                        "assay_types": sorted([str(x) for x in quantiles.index.tolist()]),
                    }
                )
                meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        except Exception:
            # Never block pipeline execution on metadata persistence.
            pass
    return df


# =============================================================================
# AGGREGATION
# =============================================================================

def _aggregate_per_link_biosample_assay(
    df: pd.DataFrame,
    include_region: bool = False,
) -> pd.DataFrame:
    """Aggregate SCREEN data to one row per (link, biosample, assay_type)."""
    group_cols = LINK_COLS + ["biosample", "assay_type"]

    agg_dict = {
        "score": ("score", "max"),
        "p_value": ("p_value", "min"),
        "any_strong": ("is_strong", "max"),
        "any_weak": ("is_weak", "max"),
        "gene_type": ("gene_type", "first"),
    }
    if include_region and "region" in df.columns:
        agg_dict["region"] = ("region", "first")

    result = df.groupby(group_cols).agg(**agg_dict).reset_index()
    result["strength"] = [
        classify_strength(s, w) for s, w in zip(result["any_strong"], result["any_weak"])
    ]
    return result


# =============================================================================
# DICT BUILDERS
# =============================================================================

def _assay_dict_from_group_rows(
    df_sub: pd.DataFrame,
    assay_types: List[str],
) -> Dict[str, Dict[str, Any]]:
    """One nested dict per (link, biosample) — matches legacy ``iterrows`` last-wins per assay."""
    out: Dict[str, Dict[str, Any]] = {a: empty_assay_entry() for a in assay_types}
    if df_sub.empty:
        return out
    at_arr = df_sub["assay_type"].to_numpy(copy=False)
    sc_arr = df_sub["score"].to_numpy(copy=False)
    pv_arr = df_sub["p_value"].to_numpy(copy=False)
    st_arr = df_sub["strength"].to_numpy(copy=False)
    idx_by_assay: Dict[Any, int] = {}
    for i in range(len(df_sub)):
        a = at_arr[i]
        idx_by_assay[a] = i
        idx_by_assay[str(a)] = i
    for assay in assay_types:
        j = idx_by_assay.get(assay)
        if j is None:
            j = idx_by_assay.get(str(assay))
        if j is None:
            continue
        out[assay] = {
            "score": safe_float(sc_arr[j]),
            "p_value": safe_float(pv_arr[j]),
            "strength": st_arr[j],
        }
    return out


def _build_biosample_columns(
    df_agg: pd.DataFrame,
    df_links: pd.DataFrame,
    biosamples: List[str],
    assay_types: List[str],
) -> pd.DataFrame:
    """Add one column per biosample with nested {assay: {score, p_value, strength}}."""

    default = empty_biosample_assays(assay_types)
    need = set(biosamples)
    by_bio = {
        b: g
        for b, g in df_agg.groupby("biosample", sort=False)
        if b in need
    }
    n_links = len(df_links)

    def _build_bio_frame(bio: str) -> Tuple[str, Optional[pd.DataFrame]]:
        df_bio = by_bio.get(bio)
        if df_bio is None or df_bio.empty:
            return bio, None
        rows_out: List[Dict[str, Any]] = []
        for key, g in df_bio.groupby(LINK_COLS, sort=False):
            d = _assay_dict_from_group_rows(g, assay_types)
            if isinstance(key, tuple):
                rec = {LINK_COLS[j]: key[j] for j in range(len(LINK_COLS))}
            else:
                rec = {LINK_COLS[0]: key}
            rec[bio] = d
            rows_out.append(rec)
        return bio, pd.DataFrame(rows_out)

    n_workers = _screen_parallel_workers(len(biosamples))
    if n_workers > 1 and len(biosamples) > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            built = list(ex.map(_build_bio_frame, biosamples))
    else:
        built = [_build_bio_frame(b) for b in biosamples]

    for bio, bio_df in built:
        if bio_df is None:
            df_links[bio] = [default.copy() for _ in range(n_links)]
            continue
        df_links = df_links.merge(bio_df, on=LINK_COLS, how="left")
        df_links[bio] = df_links[bio].map(
            lambda x: x if isinstance(x, dict) else default.copy()
        )

    return df_links


def _build_conservation_column(
    df_agg: pd.DataFrame,
    biosample_set: Set[str],
    assay_types: List[str],
    column_name: str,
) -> Optional[pd.DataFrame]:
    """Build conservation summary: {assay_type: {n_biosamples, n_strong, n_weak}}."""
    filtered = df_agg[df_agg["biosample"].isin(biosample_set)]

    if filtered.empty:
        return None

    cons = (
        filtered.groupby(LINK_COLS + ["assay_type"])
        .agg(
            n_biosamples=("biosample", "nunique"),
            n_strong=("any_strong", "sum"),
            n_weak=("any_weak", "sum"),
        )
        .reset_index()
    )

    def conservation_dict_for_link(df_sub: pd.DataFrame) -> Dict[str, Dict[str, int]]:
        result: Dict[str, Dict[str, int]] = {}
        at_arr = df_sub["assay_type"].to_numpy(copy=False)
        nb = df_sub["n_biosamples"].to_numpy(copy=False)
        ns = df_sub["n_strong"].to_numpy(copy=False)
        nw = df_sub["n_weak"].to_numpy(copy=False)
        for i in range(len(df_sub)):
            a = at_arr[i]
            result[a] = {
                "n_biosamples": int(nb[i]),
                "n_strong": int(ns[i]),
                "n_weak": int(nw[i]),
            }
        for assay in assay_types:
            if assay not in result:
                result[assay] = empty_conservation_entry()
        return result

    rows: List[Dict[str, Any]] = []
    for key, df_sub in cons.groupby(LINK_COLS, sort=False):
        d = conservation_dict_for_link(df_sub)
        if isinstance(key, tuple):
            rec = {LINK_COLS[j]: key[j] for j in range(len(LINK_COLS))}
        else:
            rec = {LINK_COLS[0]: key}
        rec[column_name] = d
        rows.append(rec)
    return pd.DataFrame(rows)


def _attach_conservation_columns(
    df_links: pd.DataFrame,
    df_agg: pd.DataFrame,
    assay_types: List[str],
    global_biosamples: Optional[Set[str]],
    breast_biosamples: Set[str],
) -> pd.DataFrame:
    """Attach conservation_global and conservation_breast columns (parallel)."""
    default_cons = empty_conservation_block(assay_types)

    global_subset = df_agg if global_biosamples is None else \
        df_agg[df_agg["biosample"].isin(global_biosamples)]
    global_bs = set(global_subset["biosample"]) if global_biosamples is None else global_biosamples

    def _global() -> Optional[pd.DataFrame]:
        return _build_conservation_column(
            global_subset, global_bs, assay_types, "conservation_global"
        )

    def _breast() -> Optional[pd.DataFrame]:
        return _build_conservation_column(
            df_agg, breast_biosamples, assay_types, "conservation_breast"
        )

    if _screen_parallel_workers(2) > 1:
        with ThreadPoolExecutor(max_workers=2) as ex:
            fut_g = ex.submit(_global)
            fut_b = ex.submit(_breast)
            cons_global = fut_g.result()
            cons_breast = fut_b.result()
    else:
        cons_global = _global()
        cons_breast = _breast()

    if cons_global is not None:
        df_links = df_links.merge(cons_global, on=LINK_COLS, how="left")
    df_links["conservation_global"] = df_links.get(
        "conservation_global", pd.Series([None] * len(df_links))
    ).map(lambda x: x if isinstance(x, dict) else default_cons.copy())

    if cons_breast is not None:
        df_links = df_links.merge(cons_breast, on=LINK_COLS, how="left")
    df_links["conservation_breast"] = df_links.get(
        "conservation_breast", pd.Series([None] * len(df_links))
    ).map(lambda x: x if isinstance(x, dict) else default_cons.copy())

    return df_links


# =============================================================================
# PUBLIC API: BUILD SCREEN LINKS
# =============================================================================

def build_screen_exp_links(
    zip_path: Path,
    inner_file: str,
    gene_list: Iterable[str],
    panel_biosamples: List[str],
    breast_biosamples: Optional[List[str]] = None,
    global_biosamples: Optional[List[str]] = None,
    weak_quantile: float = 0.90,
    strong_quantile: float = 0.95,
    intact_pvalue_threshold: Optional[float] = None,
    chunksize: int = 300_000,
    quantiles_out_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build enhancer-gene links from ENCODE SCREEN 3D-Chromatin (experimental) data.
    
    Returns:
        Tuple of (links DataFrame, list of assay types found)
    """
    print("Building SCREEN experimental links...")

    gene_set = set(gene_list)
    breast_set = {canonicalize_screen_biosample(b) for b in (breast_biosamples or panel_biosamples)}
    global_set = (
        {canonicalize_screen_biosample(b) for b in global_biosamples}
        if global_biosamples
        else None
    )

    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_file) as f:
            df = _stream_and_filter_genes(f, gene_set, COLS_SCREEN_EXP, chunksize)

    df = df.rename(columns=RENAME_SCREEN)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["biosample"] = df["biosample"].map(canonicalize_screen_biosample)

    assay_types = sorted(df["assay_type"].dropna().unique())
    print(f"  Assay types found: {assay_types}")

    df = _apply_strength_classification(
        df,
        weak_quantile,
        strong_quantile,
        intact_pvalue_threshold,
        quantiles_out_path=quantiles_out_path,
        provenance={
            "source": "SCREEN_experimental",
            "zip_path": str(zip_path),
            "inner_file": str(inner_file),
            "gene_set_size": int(len(gene_set)),
        },
    )
    df_agg = _aggregate_per_link_biosample_assay(df)
    
    links = df_agg[LINK_COLS + ["gene_type"]].drop_duplicates().reset_index(drop=True)
    links = _build_biosample_columns(df_agg, links, panel_biosamples, assay_types)
    links = _attach_conservation_columns(links, df_agg, assay_types, global_set, breast_set)

    print(f"  Built {len(links)} links")
    return links, assay_types


def build_screen_comp_links(
    gz_path: Path,
    gene_list: Iterable[str],
    panel_biosamples: List[str],
    breast_biosamples: Optional[List[str]] = None,
    global_biosamples: Optional[List[str]] = None,
    weak_quantile: float = 0.90,
    strong_quantile: float = 0.95,
    chunksize: int = 300_000,
    quantiles_out_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build enhancer-gene links from ENCODE SCREEN Computational-Methods data.
    
    Returns:
        Tuple of (links DataFrame, list of assay types found)
    """
    print("Building SCREEN computational links...")

    gene_set = set(gene_list)
    breast_set = {canonicalize_screen_biosample(b) for b in (breast_biosamples or panel_biosamples)}
    global_set = (
        {canonicalize_screen_biosample(b) for b in global_biosamples}
        if global_biosamples
        else None
    )

    df = _stream_and_filter_genes(
        str(gz_path), gene_set, COLS_SCREEN_COMP, chunksize,
        compression="gzip", add_dummy_pvalue=True
    )

    df = df.rename(columns=RENAME_SCREEN)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["biosample"] = df["biosample"].map(canonicalize_screen_biosample)

    assay_types = sorted(df["assay_type"].dropna().unique())
    print(f"  Assay types found: {assay_types}")

    df = _apply_strength_classification(
        df,
        weak_quantile,
        strong_quantile,
        quantiles_out_path=quantiles_out_path,
        provenance={
            "source": "SCREEN_computational",
            "gz_path": str(gz_path),
            "gene_set_size": int(len(gene_set)),
        },
    )
    df_agg = _aggregate_per_link_biosample_assay(df, include_region=True)
    
    links = df_agg[LINK_COLS + ["gene_type"]].drop_duplicates().reset_index(drop=True)

    # Add region column
    region_series = df_agg.groupby(LINK_COLS)["region"].first().reset_index()
    links = links.merge(region_series, on=LINK_COLS, how="left")

    links = _build_biosample_columns(df_agg, links, panel_biosamples, assay_types)
    links = _attach_conservation_columns(links, df_agg, assay_types, global_set, breast_set)

    print(f"  Built {len(links)} links")
    return links, assay_types


# =============================================================================
# COLLAPSE TO NESTED STRUCTURE
# =============================================================================

def collapse_screen_to_nested(
    df_links: pd.DataFrame,
    biosamples: List[str],
    assay_types: List[str],
    target_column: str,
) -> pd.DataFrame:
    """
    Collapse biosample columns into a single nested dict column.
    
    Structure:
    {
        "per_biosample": {biosample: {assay: {score, p_value, strength}, ...}, ...},
        "conservation_global": {...},
        "conservation_breast": {...},
    }
    """
    bio_arrays = {b: df_links[b].to_numpy(copy=False) for b in biosamples if b in df_links.columns}
    has_cg = "conservation_global" in df_links.columns
    has_cb = "conservation_breast" in df_links.columns
    cg_arr = df_links["conservation_global"].to_numpy(copy=False) if has_cg else None
    cb_arr = df_links["conservation_breast"].to_numpy(copy=False) if has_cb else None

    n = len(df_links)
    nested: List[Dict[str, Any]] = []
    for i in range(n):
        per_bio: Dict[str, Any] = {}
        for bio in biosamples:
            if bio in bio_arrays:
                v = bio_arrays[bio][i]
                per_bio[bio] = v if isinstance(v, dict) else empty_biosample_assays(assay_types)
            else:
                per_bio[bio] = empty_biosample_assays(assay_types)

        cons_global = cg_arr[i] if has_cg else None
        cons_breast = cb_arr[i] if has_cb else None

        nested.append(
            {
                "per_biosample": per_bio,
                "conservation_global": cons_global
                if isinstance(cons_global, dict)
                else empty_conservation_block(assay_types),
                "conservation_breast": cons_breast
                if isinstance(cons_breast, dict)
                else empty_conservation_block(assay_types),
            }
        )

    df_links[target_column] = nested

    drop_cols = [c for c in biosamples + ["conservation_global", "conservation_breast"] if c in df_links.columns]
    return df_links.drop(columns=drop_cols)
