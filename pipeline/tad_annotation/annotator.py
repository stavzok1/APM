"""
TAD domain annotation for features (genes, cCREs, SVs, lncRNAs).

Adds TAD_domains column to feature DataFrames with per-biosample payloads
containing domain relations, primary domain, boundary proximity, and
normalized positions.

Main function:
    annotate_df_with_tads: Annotate any feature DataFrame with TAD context

The annotation is designed to be biosample-aware, storing results keyed by
cell line name so multiple TAD sources can coexist.
"""

from __future__ import annotations

import gc
import os
from typing import Any, Callable, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .relations import (
    get_element_tad_relations,
    pick_primary_domain_id,
    boundary_overlap_and_dist,
    compute_normalized_position,
)


# =============================================================================
# UTILITIES
# =============================================================================

def ensure_chr_prefix(chrom: Union[pd.Series, str]) -> Union[pd.Series, str]:
    """Ensure chromosome names have 'chr' prefix."""
    if isinstance(chrom, str):
        return chrom if chrom.startswith("chr") else f"chr{chrom}"
    chrom = chrom.astype(str)
    return chrom.where(chrom.str.startswith("chr"), "chr" + chrom)


def _as_int(x) -> Optional[int]:
    """Safely convert to int, returning None for invalid values."""
    if x is None:
        return None
    if isinstance(x, float) and pd.isna(x):
        return None
    try:
        return int(x)
    except Exception:
        return None


def _ensure_dict_col(df: pd.DataFrame, col: str) -> pd.DataFrame:
    """Ensure column exists and contains dicts (not NaN/None)."""
    if col not in df.columns:
        df[col] = [{} for _ in range(len(df))]
    else:
        df[col] = df[col].apply(
            lambda x: {} if x is None or (isinstance(x, float) and pd.isna(x)) else x
        )
    return df


def slim_tad_payload_for_probe(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compact TAD payload for methylation probes (450K/EPIC-scale row counts).

    Drops the full per-domain relation map and large boundary blobs while
    keeping primary domain identity, coarse geometry, and normalized position
    within the domain—enough for typical probe-level QC / stratification.

    Used for methylation probes when slimming is on (default; disable with
    ``APM_METH_TAD_SLIM=0`` or ``slim_tad_payload=False``).
    """
    if not isinstance(payload, dict):
        return {}
    if "intervals" in payload:
        ints = payload["intervals"]
        if isinstance(ints, list):
            return {
                "intervals": [
                    slim_tad_payload_for_probe(x)
                    for x in ints
                    if isinstance(x, dict)
                ]
            }
        return {}

    doms = payload.get("domains")
    feat = payload.get("feature") or {}
    out: Dict[str, Any] = {}
    if isinstance(doms, dict):
        out["n_overlapping_domains"] = len(doms)
    out["feature"] = {
        k: feat.get(k)
        for k in ("chrom", "start", "end", "mid")
        if isinstance(feat, dict) and k in feat
    }

    primary = payload.get("primary")
    if not primary:
        out["primary"] = None
        return out

    norm = primary.get("normalized") or {}
    boundaries = primary.get("boundaries") or {}
    nearest = boundaries.get("nearest") if isinstance(boundaries, dict) else None
    slim_primary: Dict[str, Any] = {
        "domain_id": primary.get("domain_id"),
        "rel": primary.get("rel"),
        "normalized": {
            k: norm.get(k)
            for k in ("frac_from_left", "frac_from_right", "frac_to_nearest_boundary")
            if isinstance(norm, dict) and k in norm
        },
        "nearest_boundary": None
        if not isinstance(nearest, dict)
        else {
            "side": nearest.get("side"),
            "dist_bp": nearest.get("dist_bp"),
        },
    }
    dom = primary.get("domain")
    if isinstance(dom, dict):
        slim_primary["domain"] = {
            k: dom.get(k)
            for k in ("chrom", "start", "end", "len")
            if k in dom and dom.get(k) is not None
        }
    out["primary"] = slim_primary
    return out


def _element_tad_relations_subset(
    start: int,
    end: int,
    d_starts: np.ndarray,
    d_ends: np.ndarray,
    d_ids: np.ndarray,
    idxs: np.ndarray,
) -> Dict[str, str]:
    """Same labels as :func:`get_element_tad_relations` but only for ``idxs`` rows."""
    rel: Dict[str, str] = {}
    for j in idxs:
        a, b = int(d_starts[j]), int(d_ends[j])
        did = d_ids[j]
        key = str(did)
        if a <= start and end <= b:
            rel[key] = "contains"
        elif start <= a and b <= end:
            rel[key] = "contained"
        elif a <= start <= b < end:
            rel[key] = "overlap_left"
        elif start < a <= end <= b:
            rel[key] = "overlap_right"
        else:
            rel[key] = "overlap"
    return rel


# =============================================================================
# FEATURE INTERVAL EXTRACTION
# =============================================================================

def feature_intervals_from_row(
    row: pd.Series,
    *,
    kind: str,
    ins_mode: str = "point",
    ins_span_policy: str = "right",
    ins_fallback_bp: int = 1,
    bnd_mode: str = "point",
) -> List[Dict[str, Any]]:
    """
    Convert a DataFrame row into interval(s) for TAD annotation.
    
    Args:
        row: Single row from feature DataFrame
        kind: Feature type - "gene", "ccre", "lncrna", "sv", "bed", "interval"
        ins_mode: For INS SVs - "point" or "span"
        ins_span_policy: If ins_mode="span" - "right" or "sym"
        ins_fallback_bp: Default span if SVLEN missing
        bnd_mode: For BND SVs - currently only "point"
    
    Returns:
        List of interval dicts, each with:
            - chrom: str
            - start: int
            - end: int
            - label: str (feature type)
            - meta: dict (additional info)
    
    Notes:
        - Genes/cCREs/lncRNAs: single interval from [start, end]
        - SVs: depends on SVTYPE
            - DEL/DUP: proper interval
            - INS: point or span based on ins_mode
            - BND: point at breakend
    """
    chrom = row.get("chrom")
    if chrom is None or (isinstance(chrom, float) and pd.isna(chrom)):
        return []
    
    out: List[Dict[str, Any]] = []
    
    # Standard interval features
    if kind in {"bed", "ccre", "gene", "lncrna", "interval"}:
        s = _as_int(row.get("start"))
        e = _as_int(row.get("end"))
        if s is None or e is None:
            return []
        if e < s:
            s, e = e, s
        
        out.append({
            "chrom": chrom,
            "start": s,
            "end": e,
            "label": kind,
            "meta": {
                "strand": row.get("strand", None),
                "id": row.get("gene_name", row.get("gene_id", row.get("ccre_id", row.get("cCRE_id", None)))),
            }
        })
        return out
    
    # Structural variants - handle different SVTYPES
    if kind == "sv":
        svt = row.get("SVTYPE")
        pos = _as_int(row.get("pos"))
        end = _as_int(row.get("END"))
        svlen = _as_int(row.get("SVLEN"))
        
        if pos is None:
            return []
        
        # DEL/DUP: proper interval
        if svt in {"DEL", "DUP"}:
            if end is None:
                # Fallback: construct END from SVLEN if available
                if svlen is not None:
                    end = pos + abs(svlen)
                else:
                    end = pos
            s, e = (pos, end) if end >= pos else (end, pos)
            out.append({
                "chrom": chrom,
                "start": s,
                "end": e,
                "label": f"sv_{svt}",
                "meta": {"SVTYPE": svt, "SVLEN": svlen, "pos": pos, "END": end},
            })
            return out
        
        # INS: point or span
        if svt == "INS":
            L = abs(svlen) if svlen is not None else ins_fallback_bp
            if ins_mode == "point":
                out.append({
                    "chrom": chrom,
                    "start": pos,
                    "end": pos,
                    "label": "sv_INS",
                    "meta": {"SVTYPE": svt, "SVLEN": svlen, "pos": pos},
                })
            else:
                if ins_span_policy == "sym":
                    half = max(0, L // 2)
                    s, e = pos - half, pos + (L - half)
                else:  # "right"
                    s, e = pos, pos + L
                if s < 0:
                    s = 0
                out.append({
                    "chrom": chrom,
                    "start": int(s),
                    "end": int(e),
                    "label": "sv_INS",
                    "meta": {"SVTYPE": svt, "SVLEN": svlen, "pos": pos, "span_policy": ins_span_policy},
                })
            return out
        
        # BND: point at breakend
        if svt == "BND":
            out.append({
                "chrom": chrom,
                "start": pos,
                "end": pos,
                "label": "sv_BND",
                "meta": {"SVTYPE": svt, "pos": pos},
            })
            return out
        
        # Unknown SVTYPE: treat as point at pos
        out.append({
            "chrom": chrom,
            "start": pos,
            "end": pos,
            "label": "sv_unknown",
            "meta": {"SVTYPE": svt, "pos": pos, "END": end, "SVLEN": svlen},
        })
        return out
    
    raise ValueError(f"Unknown kind={kind!r}")


# =============================================================================
# MAIN ANNOTATOR
# =============================================================================

def annotate_df_with_tads(
    df: pd.DataFrame,
    *,
    kind: str,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
    out_col: str = "TAD_domains",
    ins_mode: str = "point",
    ins_span_policy: str = "right",
    payload_postprocess: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None,
    use_vectorized_interval: bool = False,
) -> pd.DataFrame:
    """
    Annotate a feature DataFrame with TAD domain context.
    
    Adds df[out_col][i][biosample] = payload for each row.
    
    Args:
        df: Feature DataFrame (genes, cCREs, SVs, lncRNAs)
        kind: Feature type - "gene", "ccre", "lncrna", "sv"
        tad_domains: DataFrame with [domain_id, chrom, start, end]
        domain_flanks: DataFrame with [domain_id, left_boundary_id, right_boundary_id, 
                                       left_strength, right_strength, n_flanking_boundaries]
        boundaries: DataFrame with [boundary_id, chrom, start, end, pos]
        biosample: Cell line / sample name (key in the nested dict)
        out_col: Column name for TAD annotation
        ins_mode: For SVs - "point" or "span" for insertions
        ins_span_policy: If ins_mode="span" - "right" or "sym"
    
    Returns:
        DataFrame with TAD_domains column containing per-biosample payloads.
    
    Payload structure (single interval):
        {
            "domains": {domain_id: relation, ...},
            "primary": {
                "domain_id": str,
                "rel": str,
                "domain": {chrom, start, end, len},
                "feature": {chrom, start, end, mid, ...},
                "boundaries": {left: {...}, right: {...}, nearest: {...}},
                "normalized": {frac_from_left, frac_from_right, frac_to_nearest_boundary}
            },
            "feature": {chrom, start, end, ...}
        }
    
    For multi-interval features (some SVs):
        {"intervals": [payload, payload, ...]}
    
    Example:
        >>> genes_df = annotate_df_with_tads(
        ...     genes_df,
        ...     kind="gene",
        ...     tad_domains=domains,
        ...     domain_flanks=flanks,
        ...     boundaries=bounds,
        ...     biosample="Kim_T47D",
        ... )

    Environment:
        ``APM_TAD_ANNOTATE_CHUNK`` (default ``40000``): max rows per chromosome
        batch when applying row-wise payloads. Lower this (e.g. ``20000``) if
        annotating very wide tables (many columns) still risks OOM; ``0`` disables
        chunking (single apply per chromosome).

        Methylation probes use ``slim_tad_payload_for_probe`` by default (see
        ``annotate_probes_with_tads``); set ``APM_METH_TAD_SLIM=0`` for full payloads.

        ``use_vectorized_interval=True`` (probe path): interval–domain overlap via
        NumPy broadcasting instead of per-row DataFrame filters (faster).
    """
    
    # Ensure output column exists and is dict-per-row
    df = _ensure_dict_col(df, out_col)
    
    # Normalize chromosome naming across all inputs
    tad_domains["chrom"] = ensure_chr_prefix(tad_domains["chrom"])
    boundaries["chrom"] = ensure_chr_prefix(boundaries["chrom"])
    df["chrom"] = ensure_chr_prefix(df["chrom"])

    tad_by_chrom = {c: g for c, g in tad_domains.groupby("chrom", sort=False)}
    
    # Build fast lookups
    flanks_by_domain = domain_flanks.set_index("domain_id", drop=False).to_dict(orient="index")
    
    # Handle different boundary ID column names
    bid_col = "boundary_id" if "boundary_id" in boundaries.columns else (
        "name" if "name" in boundaries.columns else None
    )
    if bid_col is None:
        raise ValueError("boundaries must contain 'boundary_id' or 'name' column")
    bounds_by_id = boundaries.set_index(bid_col, drop=False).to_dict(orient="index")
    
    # Helper: build payload for a single interval (``tad_chrom`` = domains on this chromosome only;
    # hoisted once per chromosome — previously this filter ran inside every row call).
    def build_payload_for_interval(
        chrom: str,
        s: int,
        e: int,
        meta: Dict[str, Any],
        *,
        tad_chrom: pd.DataFrame,
        dom_ix: Optional[pd.DataFrame] = None,
        rels_precomputed: Optional[Dict[str, str]] = None,
        primary_id_precomputed: Optional[Any] = None,
    ) -> Dict[str, Any]:
        if tad_chrom.empty:
            return {"domains": {}, "primary": None, "feature": {"chrom": chrom, "start": s, "end": e, **meta}}

        if rels_precomputed is not None:
            rels = rels_precomputed
        else:
            rels = get_element_tad_relations(s, e, tad_chrom)
        if primary_id_precomputed is not None:
            primary_id = primary_id_precomputed
        else:
            primary_id = pick_primary_domain_id(s, e, tad_chrom)

        primary = None
        if primary_id is not None:
            if dom_ix is not None and primary_id in dom_ix.index:
                drow = dom_ix.loc[primary_id]
                if isinstance(drow, pd.DataFrame):
                    drow = drow.iloc[0]
            else:
                drow = tad_chrom[tad_chrom["domain_id"] == primary_id].iloc[0]
            d_start, d_end = int(drow["start"]), int(drow["end"])
            d_len = max(1, d_end - d_start)

            mid = (s + e) // 2
            primary_rel = rels.get(primary_id)
            if primary_rel is None:
                primary_rel = rels.get(str(primary_id), "none")
            
            flank = flanks_by_domain.get(primary_id)
            left_b = right_b = nearest = None
            
            if flank is not None:
                left_id = flank.get("left_boundary_id")
                right_id = flank.get("right_boundary_id")
                
                # Build left boundary info
                if left_id in bounds_by_id:
                    b = bounds_by_id[left_id]
                    b_start = int(b["start"])
                    b_end = int(b["end"])
                    b_pos = int(b.get("pos", (b_start + b_end) // 2))
                    left_b = {
                        "boundary_id": left_id,
                        "start": b_start,
                        "end": b_end,
                        "pos": b_pos,
                        "strength": flank.get("left_strength"),
                        **boundary_overlap_and_dist(s, e, b_start, b_end, b_pos),
                    }
                
                # Build right boundary info
                if right_id in bounds_by_id:
                    b = bounds_by_id[right_id]
                    b_start = int(b["start"])
                    b_end = int(b["end"])
                    b_pos = int(b.get("pos", (b_start + b_end) // 2))
                    right_b = {
                        "boundary_id": right_id,
                        "start": b_start,
                        "end": b_end,
                        "pos": b_pos,
                        "strength": flank.get("right_strength"),
                        **boundary_overlap_and_dist(s, e, b_start, b_end, b_pos),
                    }
                
                # Determine nearest boundary
                candidates = []
                if left_b is not None:
                    candidates.append(("left", left_b["dist_bp"]))
                if right_b is not None:
                    candidates.append(("right", right_b["dist_bp"]))
                if candidates:
                    side, dist_bp = min(candidates, key=lambda t: t[1])
                    nearest = {"side": side, "dist_bp": int(dist_bp)}
            
            # Compute normalized position
            normalized = compute_normalized_position(s, e, d_start, d_end)
            
            primary = {
                "domain_id": primary_id,
                "rel": primary_rel,
                "domain": {"chrom": chrom, "start": d_start, "end": d_end, "len": int(d_len)},
                "feature": {"chrom": chrom, "start": s, "end": e, "mid": int(mid), **meta},
                "boundaries": {"left": left_b, "right": right_b, "nearest": nearest},
                "normalized": normalized,
            }
        
        return {"domains": rels, "primary": primary, "feature": {"chrom": chrom, "start": s, "end": e, **meta}}
    
    # Per-chromosome: build payloads then merge into df[out_col][row][biosample].
    # IMPORTANT: merge in-place into each row's dict. Using ``{**old, biosample: payload}``
    # allocates a new dict for every row on every biosample (e.g. methylation probes ×
    # 20+ TAD maps) and can drive RSS into OOM range.
    chunk_raw = os.environ.get("APM_TAD_ANNOTATE_CHUNK", "40000").strip()
    try:
        row_chunk = max(0, int(chunk_raw))
    except ValueError:
        row_chunk = 40_000

    for chrom in df["chrom"].dropna().unique():
        mask = df["chrom"] == chrom
        if not mask.any():
            continue

        tad_chrom = tad_by_chrom.get(chrom, tad_domains.iloc[0:0])
        dom_ix: Optional[pd.DataFrame] = None
        if not tad_chrom.empty and "domain_id" in tad_chrom.columns:
            dom_ix = tad_chrom.set_index("domain_id", drop=False)

        chrom_str = str(chrom)

        def _assign_payload(idx: Any, payload: Dict[str, Any]) -> None:
            if payload_postprocess is not None:
                payload = payload_postprocess(payload)
            cell = df.at[idx, out_col]
            if not isinstance(cell, dict):
                cell = {}
                df.at[idx, out_col] = cell
            cell[biosample] = payload

        def build_row_payload(row: pd.Series) -> Dict[str, Any]:
            intervals = feature_intervals_from_row(
                row,
                kind=kind,
                ins_mode=ins_mode,
                ins_span_policy=ins_span_policy,
            )
            if not intervals:
                return {"domains": {}, "primary": None}

            per_int = []
            for itv in intervals:
                s, e = int(itv["start"]), int(itv["end"])
                meta = {"label": itv.get("label")}
                meta.update(itv.get("meta", {}) or {})
                per_int.append(
                    build_payload_for_interval(
                        chrom_str, s, e, meta, tad_chrom=tad_chrom, dom_ix=dom_ix
                    )
                )

            if len(per_int) == 1:
                return per_int[0]
            return {"intervals": per_int}

        sub_idx = df.index[mask]
        use_vec = bool(use_vectorized_interval and kind == "interval")

        def _flush_chunk_vectorized(chunk_index: pd.Index) -> None:
            chunk_df = df.loc[chunk_index]
            s_ser = pd.to_numeric(chunk_df["start"], errors="coerce")
            e_ser = pd.to_numeric(chunk_df["end"], errors="coerce")
            valid = s_ser.notna() & e_ser.notna()
            s_np = s_ser.to_numpy(dtype=np.float64, copy=False)
            e_np = e_ser.to_numpy(dtype=np.float64, copy=False)
            s_i = np.nan_to_num(s_np, nan=0.0).astype(np.int64)
            e_i = np.nan_to_num(e_np, nan=0.0).astype(np.int64)
            valid_np = valid.to_numpy()

            if tad_chrom.empty:
                ov = np.zeros((len(chunk_df), 0), dtype=bool)
                d_starts = np.array([], dtype=np.int64)
                d_ends = np.array([], dtype=np.int64)
                d_ids = np.array([], dtype=object)
            else:
                d_starts = tad_chrom["start"].to_numpy(dtype=np.int64, copy=False)
                d_ends = tad_chrom["end"].to_numpy(dtype=np.int64, copy=False)
                d_ids = tad_chrom["domain_id"].to_numpy(copy=False)
                ov = (d_ends[None, :] >= s_i[:, None]) & (d_starts[None, :] <= e_i[:, None])

            for row_k, idx in enumerate(chunk_df.index.tolist()):
                if not valid_np[row_k]:
                    _assign_payload(idx, {"domains": {}, "primary": None})
                    continue
                s, e_row = int(s_i[row_k]), int(e_i[row_k])
                idxs = np.flatnonzero(ov[row_k])
                if idxs.size == 0:
                    payload = build_payload_for_interval(
                        chrom_str,
                        s,
                        e_row,
                        {"label": "interval"},
                        tad_chrom=tad_chrom,
                        dom_ix=dom_ix,
                    )
                else:
                    rels = _element_tad_relations_subset(
                        s, e_row, d_starts, d_ends, d_ids, idxs
                    )
                    sub = tad_chrom.iloc[idxs]
                    primary_id = pick_primary_domain_id(s, e_row, sub)
                    if primary_id is None:
                        primary_id = pick_primary_domain_id(s, e_row, tad_chrom)
                    payload = build_payload_for_interval(
                        chrom_str,
                        s,
                        e_row,
                        {"label": "interval"},
                        tad_chrom=tad_chrom,
                        dom_ix=dom_ix,
                        rels_precomputed=rels,
                        primary_id_precomputed=primary_id,
                    )
                _assign_payload(idx, payload)
            del chunk_df, ov

        def _flush_chunk_apply(chunk_index: pd.Index) -> None:
            chunk_df = df.loc[chunk_index]
            pay = chunk_df.apply(build_row_payload, axis=1)
            for idx, payload in pay.items():
                _assign_payload(idx, payload)
            del pay, chunk_df

        flush = _flush_chunk_vectorized if use_vec else _flush_chunk_apply
        if row_chunk <= 0 or len(sub_idx) <= row_chunk:
            flush(sub_idx)
        else:
            for start in range(0, len(sub_idx), row_chunk):
                flush(sub_idx[start : start + row_chunk])
                gc.collect()

        gc.collect()

    return df


# =============================================================================
# CONVENIENCE WRAPPERS
# =============================================================================

def annotate_genes_with_tads(
    genes_df: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
) -> pd.DataFrame:
    """Convenience wrapper for gene annotation."""
    return annotate_df_with_tads(
        genes_df,
        kind="gene",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
    )


def annotate_ccres_with_tads(
    ccre_df: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
) -> pd.DataFrame:
    """Convenience wrapper for cCRE annotation."""
    return annotate_df_with_tads(
        ccre_df,
        kind="ccre",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
    )


def annotate_lncrnas_with_tads(
    lncrna_df: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
) -> pd.DataFrame:
    """Convenience wrapper for lncRNA annotation."""
    return annotate_df_with_tads(
        lncrna_df,
        kind="lncrna",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
    )


def annotate_svs_with_tads(
    sv_df: pd.DataFrame,
    tad_domains: pd.DataFrame,
    domain_flanks: pd.DataFrame,
    boundaries: pd.DataFrame,
    biosample: str,
    ins_mode: str = "point",
    ins_span_policy: str = "right",
) -> pd.DataFrame:
    """Convenience wrapper for SV annotation."""
    return annotate_df_with_tads(
        sv_df,
        kind="sv",
        tad_domains=tad_domains,
        domain_flanks=domain_flanks,
        boundaries=boundaries,
        biosample=biosample,
        ins_mode=ins_mode,
        ins_span_policy=ins_span_policy,
    )
