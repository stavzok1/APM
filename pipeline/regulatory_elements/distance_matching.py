"""
cCRE to gene distance matching.

Matches regulatory elements to genes based on TSS distance,
organizing results into distance tiers.
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from ..utils import (
    assign_distance_tier,
    compute_tss,
)
from ..config import THRESHOLDS


# =============================================================================
# CORE MATCHING
# =============================================================================

def _per_chrom_match(
    gc_: pd.DataFrame,
    cc_: pd.DataFrame,
    window_bp: int,
) -> pd.DataFrame:
    """
    Match genes ↔ cCREs on a single chromosome via sorted-array searchsorted.

    For each gene with TSS window ``[win_start, win_end]`` we select the sorted cCREs whose
    ``start <= win_end`` (right edge) and then filter ``end >= win_start``. This replaces the
    global ``explode``-then-``merge`` path (100M-row intermediates) with O(n_genes · log n_cCREs)
    index lookups and a single boolean mask per gene.
    """
    if gc_.empty or cc_.empty:
        return pd.DataFrame()

    cc_sorted = cc_.sort_values("start", kind="mergesort")
    cc_start_arr = cc_sorted["start"].to_numpy(dtype=np.int64, copy=False)
    cc_end_arr = cc_sorted["end"].to_numpy(dtype=np.int64, copy=False)
    cc_center_arr = cc_sorted["center"].to_numpy(dtype=np.int64, copy=False)
    cc_id_arr = cc_sorted["cCRE_id"].astype(object).to_numpy(copy=False)
    cc_enc_arr = cc_sorted["ENCODE_id"].astype(object).to_numpy(copy=False)
    cc_type_arr = cc_sorted["type"].astype(object).to_numpy(copy=False)
    cc_raw_arr = cc_sorted["raw_type"].astype(object).to_numpy(copy=False)
    cc_chrom_arr = cc_sorted["chrom"].astype(object).to_numpy(copy=False)

    gene_names = gc_["gene_name"].astype(object).to_numpy(copy=False)
    tss_arr = gc_["tss"].to_numpy(dtype=np.int64, copy=False)
    ws_arr = gc_["win_start"].to_numpy(dtype=np.int64, copy=False)
    we_arr = gc_["win_end"].to_numpy(dtype=np.int64, copy=False)

    # Right-edge cutoff via searchsorted (1 log query per gene).
    right_idx = np.searchsorted(cc_start_arr, we_arr, side="right")

    gene_parts: List[np.ndarray] = []
    cidx_parts: List[np.ndarray] = []
    dist_parts: List[np.ndarray] = []
    tss_parts: List[np.ndarray] = []

    for i in range(len(gc_)):
        r = int(right_idx[i])
        if r == 0:
            continue
        ws = int(ws_arr[i])
        # Candidates: cCREs with start <= win_end AND end >= win_start.
        ends_sub = cc_end_arr[:r]
        mask = ends_sub >= ws
        if not mask.any():
            continue
        idx = np.nonzero(mask)[0]
        tss_i = int(tss_arr[i])
        starts_sub = cc_start_arr[idx]
        ends_hit = ends_sub[idx]
        # Exact TSS-to-interval distance (vectorized numpy).
        inside = (tss_i >= starts_sub) & (tss_i <= ends_hit)
        dist = np.where(
            inside,
            0,
            np.minimum(np.abs(tss_i - starts_sub), np.abs(tss_i - ends_hit)),
        ).astype(np.int64)
        keep = dist <= window_bp
        if not keep.any():
            continue
        idx = idx[keep]
        dist = dist[keep]

        gene_parts.append(np.full(idx.shape[0], gene_names[i], dtype=object))
        tss_parts.append(np.full(idx.shape[0], tss_i, dtype=np.int64))
        cidx_parts.append(idx)
        dist_parts.append(dist)

    if not gene_parts:
        return pd.DataFrame()

    gene_flat = np.concatenate(gene_parts)
    cidx_flat = np.concatenate(cidx_parts)
    dist_flat = np.concatenate(dist_parts)
    tss_flat = np.concatenate(tss_parts)

    return pd.DataFrame(
        {
            "gene_name": gene_flat,
            "cCRE_id": cc_id_arr[cidx_flat],
            "ENCODE_id": cc_enc_arr[cidx_flat],
            "type": cc_type_arr[cidx_flat],
            "raw_type": cc_raw_arr[cidx_flat],
            "chrom": cc_chrom_arr[cidx_flat],
            "start": cc_start_arr[cidx_flat],
            "end": cc_end_arr[cidx_flat],
            "center": cc_center_arr[cidx_flat],
            "tss": tss_flat,
            "dist_to_tss": dist_flat,
        }
    )


def match_ccres_to_genes(
    genes: pd.DataFrame,
    ccres: pd.DataFrame,
    window_bp: int = 1_000_000,
    bin_size: int = 100_000,
    tier_edges: List[int] = None,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Match cCREs to genes based on distance from TSS.

    Per-chromosome sorted-array range join (``numpy.searchsorted``) replaces the global
    bin-explode + merge path: intermediate row blow-up is avoided entirely, so runtime is
    proportional to the number of kept pairs (plus a tiny ``log(n_cCREs)`` overhead per gene).
    ``bin_size`` is accepted for API compatibility and currently unused.
    """
    # Use defaults from config if not provided
    if tier_edges is None:
        tier_edges = THRESHOLDS.tier_edges_bp
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels

    # Prepare genes (cheap: int casts + TSS/window)
    g = genes[["chrom", "start", "end", "strand", "gene_name"]].copy()
    g["start"] = g["start"].astype(np.int64)
    g["end"] = g["end"].astype(np.int64)
    g["tss"] = compute_tss(g["start"], g["end"], g["strand"]).astype(np.int64)
    g["win_start"] = (g["tss"] - window_bp).clip(lower=0).astype(np.int64)
    g["win_end"] = (g["tss"] + window_bp).astype(np.int64)

    # Prepare cCREs
    c = ccres[["chrom", "start", "end", "ENCODE_id", "cCRE_id", "type", "raw_type"]].copy()
    c["start"] = c["start"].astype(np.int64)
    c["end"] = c["end"].astype(np.int64)
    c["center"] = ((c["start"] + c["end"]) // 2).astype(np.int64)

    # Cast chrom to str for robust grouping across categorical inputs.
    g["chrom"] = g["chrom"].astype(str)
    c["chrom"] = c["chrom"].astype(str)

    # Per-chromosome range join (very cheap vs. global explode/merge).
    parts: List[pd.DataFrame] = []
    g_by_chrom = {ch: sub for ch, sub in g.groupby("chrom", sort=False)}
    for chrom, cc_ in c.groupby("chrom", sort=False):
        gc_ = g_by_chrom.get(chrom)
        if gc_ is None or gc_.empty:
            continue
        df_part = _per_chrom_match(gc_, cc_, window_bp)
        if not df_part.empty:
            parts.append(df_part)

    if not parts:
        pair_df = pd.DataFrame(
            columns=[
                "gene_name", "cCRE_id", "ENCODE_id", "type", "raw_type",
                "chrom", "start", "end", "center", "tss", "dist_to_tss", "tier",
            ]
        )
        print(f"Matched 0 gene-cCRE pairs within {window_bp}bp")
        return pair_df

    try:
        pair_df = pd.concat(parts, ignore_index=True, copy=False)
    except TypeError:
        pair_df = pd.concat(parts, ignore_index=True)
    parts.clear()

    # Assign distance tiers (vectorized pd.cut).
    pair_df["tier"] = assign_distance_tier(pair_df["dist_to_tss"], tier_edges, tier_labels)

    # Per-chromosome builder can emit duplicates only if (gene, cCRE) truly overlaps
    # multiple times; in practice never, but keep the safety net (fast: hash-based).
    pair_df = pair_df.drop_duplicates(
        subset=["gene_name", "cCRE_id", "ENCODE_id"]
    ).reset_index(drop=True)

    print(f"Matched {len(pair_df)} gene-cCRE pairs within {window_bp}bp")

    return pair_df


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def save_distance_matching(
    output_dir: Path,
    pair_df: pd.DataFrame,
    prefix: str = "",
) -> None:
    """
    Save gene-cCRE distance matching results.
    
    Args:
        output_dir: Directory to save files
        pair_df: DataFrame from match_ccres_to_genes
        prefix: Optional prefix for filenames (e.g., "lncRNA_")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filename = f"{prefix}gene_to_elements.csv" if prefix else "gene_to_elements.csv"
    pair_df.to_csv(output_dir / filename, index=False)
    
    print(f"Saved distance matching to {output_dir / filename}")



def save_all_matching_outputs(
    output_dir: Path,
    pair_df: pd.DataFrame,
    elem_focus: pd.DataFrame,
    gene_summary: pd.DataFrame,
    label: str = "coding",
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # pairs
    save_distance_matching(output_dir, pair_df, prefix=f"{label}_")

    # tables
    elem_focus.to_csv(output_dir / f"{label}_element_focus.csv", index=False)
    gene_summary.to_csv(output_dir / f"{label}_gene_summary.csv", index=False)

    print(f"Saved all matching outputs for {label} to {output_dir}")


    


# =============================================================================
# AGGREGATION BY GENE
# =============================================================================

def _groupby_join_unique_sorted(
    df: pd.DataFrame,
    group_cols: List[str],
    value_col: str,
) -> pd.DataFrame:
    """
    Fast equivalent of ``groupby(...).apply(lambda s: ','.join(sorted(unique(s))))``.

    One global sort + one ``drop_duplicates`` + one ``groupby.agg(','.join)`` — avoids the
    Python-per-group lambda overhead of the previous implementation.
    """
    if df.empty:
        cols = group_cols + ["ids"]
        return pd.DataFrame(columns=cols)
    work = df[group_cols + [value_col]].drop_duplicates()
    work = work.sort_values(group_cols + [value_col], kind="mergesort")
    work[value_col] = work[value_col].astype(str)
    out = (
        work.groupby(group_cols, observed=True, sort=False)[value_col]
        .agg(",".join)
        .rename("ids")
        .reset_index()
    )
    return out


def aggregate_ccres_per_gene(
    pair_df: pd.DataFrame,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Create gene-level summary of cCRE counts and IDs per type/tier.
    
    Returns wide DataFrame with columns like:
    - "pELS | 0-100kb | count"
    - "pELS | 0-100kb | ids"
    """
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels

    # Counts per (gene, type, tier)
    cnt = (
        pair_df
        .groupby(["gene_name", "type", "tier"], observed=True)["cCRE_id"]
        .nunique()
        .rename("count")
        .reset_index()
    )

    cnt_wide = (
        cnt
        .assign(col=lambda d: d["type"].astype(str) + " | " + d["tier"].astype(str) + " | count")
        .pivot(index="gene_name", columns="col", values="count")
        .reset_index()
    )

    # IDs per (gene, type, tier) — fast path via pre-sorted drop_duplicates + agg ','.join.
    ids = _groupby_join_unique_sorted(
        pair_df, ["gene_name", "type", "tier"], "cCRE_id"
    )

    ids_wide = (
        ids
        .assign(col=lambda d: d["type"].astype(str) + " | " + d["tier"].astype(str) + " | ids")
        .pivot(index="gene_name", columns="col", values="ids")
        .reset_index()
    )

    # Merge counts and IDs
    gene_summary = cnt_wide.merge(ids_wide, on="gene_name", how="outer")

    # Fill NaN appropriately
    num_cols = [c for c in gene_summary.columns if c.endswith(" | count")]
    id_cols = [c for c in gene_summary.columns if c.endswith(" | ids")]

    for c in num_cols:
        gene_summary[c] = pd.to_numeric(gene_summary[c], errors="coerce").fillna(0).astype(int)
    for c in id_cols:
        gene_summary[c] = gene_summary[c].astype(object).fillna("")

    return gene_summary


# =============================================================================
# AGGREGATION BY ELEMENT
# =============================================================================

def aggregate_genes_per_ccre(
    pair_df: pd.DataFrame,
    tier_labels: List[str] = None,
) -> pd.DataFrame:
    """
    Create element-level summary with genes organized by distance tier.
    
    Returns DataFrame with columns:
    - cCRE_id, ENCODE_id, type, raw_type, chrom, start, end
    - One column per tier containing comma-separated gene names
    - min_dist_to_any_gene
    - genes_by_exact_dist (format: "gene1:123,gene2:456")
    """
    if tier_labels is None:
        tier_labels = THRESHOLDS.tier_labels
    
    # Genes per tier — fast path via pre-sorted unique rows + agg ','.join.
    idx_cols = ["cCRE_id", "ENCODE_id", "type", "raw_type", "chrom", "start", "end", "tier"]
    agg_genes = _groupby_join_unique_sorted(pair_df, idx_cols, "gene_name")
    agg_genes = agg_genes.rename(columns={"ids": "gene_name"})

    elem_focus = (
        agg_genes
        .pivot(
            index=["cCRE_id", "ENCODE_id", "type", "raw_type", "chrom", "start", "end"],
            columns="tier",
            values="gene_name",
        )
        .reset_index()
    )

    # Ensure tier columns are object type before fillna
    tier_cols = [c for c in elem_focus.columns if c in tier_labels]
    for c in tier_cols:
        elem_focus[c] = elem_focus[c].astype(object).fillna("")

    # Min distance to any gene
    min_dist = (
        pair_df
        .groupby("cCRE_id", observed=True)["dist_to_tss"]
        .min()
        .rename("min_dist_to_any_gene")
        .reset_index()
    )
    elem_focus = elem_focus.merge(min_dist, on="cCRE_id", how="left")

    # Exact "gene:dist" string — build a sorted ``gene:dist`` token column once and then use
    # a single vectorized ``groupby.agg(','.join)`` to avoid the per-group Python lambda.
    sort_df = pair_df.sort_values(
        ["cCRE_id", "dist_to_tss", "gene_name"], kind="mergesort"
    )
    tokens = (
        sort_df["gene_name"].astype(str).to_numpy()
        + ":"
        + sort_df["dist_to_tss"].astype(np.int64).astype(str).to_numpy()
    )
    tokens_df = pd.DataFrame(
        {"cCRE_id": sort_df["cCRE_id"].to_numpy(copy=False), "_tok": tokens}
    )
    exact_dist = (
        tokens_df.groupby("cCRE_id", observed=True, sort=False)["_tok"]
        .agg(",".join)
        .rename("genes_by_exact_dist")
        .reset_index()
    )
    elem_focus = elem_focus.merge(exact_dist, on="cCRE_id", how="left")

    return elem_focus


def build_distance_matrix(pair_df: pd.DataFrame) -> pd.DataFrame:
    """
    Build gene × cCRE distance matrix.
    
    Returns pivot table with genes as rows, cCREs as columns,
    minimum distance as values.
    """
    matrix = (
        pair_df
        .groupby(["gene_name", "cCRE_id"], observed=True)["dist_to_tss"]
        .min()
        .unstack("cCRE_id")
        .astype("float64")
    )
    
    # Order columns by global min distance
    elem_order = pair_df.groupby("cCRE_id", observed=True)["dist_to_tss"].min().sort_values().index
    matrix = matrix.reindex(columns=elem_order)
    
    return matrix


