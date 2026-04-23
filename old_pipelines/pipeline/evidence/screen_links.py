"""
ENCODE SCREEN enhancer-gene link processing.

Handles both:
- Experimental (3D chromatin) links
- Computational (ABC, etc.) links

Complete pipeline: load → filter → aggregate → strength classify → build nested dicts
"""

import zipfile
from pathlib import Path
from typing import List, Set, Dict, Any, Optional, Tuple, Iterable

import numpy as np
import pandas as pd

from ..utils import classify_strength, safe_float
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

def _build_biosample_columns(
    df_agg: pd.DataFrame,
    df_links: pd.DataFrame,
    biosamples: List[str],
    assay_types: List[str],
) -> pd.DataFrame:
    """Add one column per biosample with nested {assay: {score, p_value, strength}}."""

    def make_assay_dict(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        rows_by_assay = {r["assay_type"]: r for _, r in df_sub.iterrows()}
        out = {}
        for assay in assay_types:
            if assay in rows_by_assay:
                r = rows_by_assay[assay]
                out[assay] = {
                    "score": safe_float(r["score"]),
                    "p_value": safe_float(r["p_value"]),
                    "strength": r["strength"],
                }
            else:
                out[assay] = empty_assay_entry()
        return out

    default = empty_biosample_assays(assay_types)

    for bio in biosamples:
        df_bio = df_agg[df_agg["biosample"] == bio]
        if df_bio.empty:
            df_links[bio] = [default.copy() for _ in range(len(df_links))]
        else:
            bio_series = df_bio.groupby(LINK_COLS).apply(make_assay_dict)
            df_links = df_links.merge(bio_series.rename(bio).reset_index(), on=LINK_COLS, how="left")
            df_links[bio] = df_links[bio].apply(lambda x: x if isinstance(x, dict) else default.copy())

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

    def to_dict(df_sub: pd.DataFrame) -> Dict[str, Dict[str, int]]:
        result = {r["assay_type"]: {
            "n_biosamples": int(r["n_biosamples"]),
            "n_strong": int(r["n_strong"]),
            "n_weak": int(r["n_weak"]),
        } for _, r in df_sub.iterrows()}
        
        for assay in assay_types:
            if assay not in result:
                result[assay] = empty_conservation_entry()
        return result

    return cons.groupby(LINK_COLS).apply(to_dict).rename(column_name).reset_index()


def _attach_conservation_columns(
    df_links: pd.DataFrame,
    df_agg: pd.DataFrame,
    assay_types: List[str],
    global_biosamples: Optional[Set[str]],
    breast_biosamples: Set[str],
) -> pd.DataFrame:
    """Attach conservation_global and conservation_breast columns."""
    default_cons = empty_conservation_block(assay_types)

    # Global conservation
    global_subset = df_agg if global_biosamples is None else \
        df_agg[df_agg["biosample"].isin(global_biosamples)]
    
    cons_global = _build_conservation_column(
        global_subset, set(global_subset["biosample"]), assay_types, "conservation_global"
    )
    if cons_global is not None:
        df_links = df_links.merge(cons_global, on=LINK_COLS, how="left")
    
    df_links["conservation_global"] = df_links.get("conservation_global", pd.Series([None] * len(df_links))).apply(
        lambda x: x if isinstance(x, dict) else default_cons.copy()
    )

    # Breast conservation
    cons_breast = _build_conservation_column(
        df_agg, breast_biosamples, assay_types, "conservation_breast"
    )
    if cons_breast is not None:
        df_links = df_links.merge(cons_breast, on=LINK_COLS, how="left")
    
    df_links["conservation_breast"] = df_links.get("conservation_breast", pd.Series([None] * len(df_links))).apply(
        lambda x: x if isinstance(x, dict) else default_cons.copy()
    )

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
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build enhancer-gene links from ENCODE SCREEN 3D-Chromatin (experimental) data.
    
    Returns:
        Tuple of (links DataFrame, list of assay types found)
    """
    print("Building SCREEN experimental links...")

    gene_set = set(gene_list)
    breast_set = set(breast_biosamples or panel_biosamples)
    global_set = set(global_biosamples) if global_biosamples else None

    with zipfile.ZipFile(zip_path) as z:
        with z.open(inner_file) as f:
            df = _stream_and_filter_genes(f, gene_set, COLS_SCREEN_EXP, chunksize)

    df = df.rename(columns=RENAME_SCREEN)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    assay_types = sorted(df["assay_type"].dropna().unique())
    print(f"  Assay types found: {assay_types}")

    df = _apply_strength_classification(df, weak_quantile, strong_quantile, intact_pvalue_threshold)
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
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Build enhancer-gene links from ENCODE SCREEN Computational-Methods data.
    
    Returns:
        Tuple of (links DataFrame, list of assay types found)
    """
    print("Building SCREEN computational links...")

    gene_set = set(gene_list)
    breast_set = set(breast_biosamples or panel_biosamples)
    global_set = set(global_biosamples) if global_biosamples else None

    df = _stream_and_filter_genes(
        str(gz_path), gene_set, COLS_SCREEN_COMP, chunksize,
        compression="gzip", add_dummy_pvalue=True
    )

    df = df.rename(columns=RENAME_SCREEN)
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    assay_types = sorted(df["assay_type"].dropna().unique())
    print(f"  Assay types found: {assay_types}")

    df = _apply_strength_classification(df, weak_quantile, strong_quantile)
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
    def row_to_nested(row) -> Dict[str, Any]:
        per_bio = {}
        for bio in biosamples:
            if bio in row.index and isinstance(row[bio], dict):
                per_bio[bio] = row[bio]
            else:
                per_bio[bio] = empty_biosample_assays(assay_types)

        cons_global = row.get("conservation_global")
        cons_breast = row.get("conservation_breast")

        return {
            "per_biosample": per_bio,
            "conservation_global": cons_global if isinstance(cons_global, dict) else empty_conservation_block(assay_types),
            "conservation_breast": cons_breast if isinstance(cons_breast, dict) else empty_conservation_block(assay_types),
        }

    df_links[target_column] = df_links.apply(row_to_nested, axis=1)

    drop_cols = [c for c in biosamples + ["conservation_global", "conservation_breast"] if c in df_links.columns]
    return df_links.drop(columns=drop_cols)
