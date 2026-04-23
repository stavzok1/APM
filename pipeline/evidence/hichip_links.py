"""
HiChIP loop data processing.

Complete pipeline: load loops → normalize formats → compute overlaps → build nested dicts
"""

import os
import re
import glob
import gc
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Iterable
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    import pyranges as pr
    HAS_PYRANGES = True
except ImportError:
    HAS_PYRANGES = False
    print("Warning: pyranges not installed. HiChIP integration will be limited.")

from ..biosample_names import canonical_hichip_output_key
from ..utils import safe_float, safe_int, parse_hichip_anchor2


# =============================================================================
# LOOP FILE LOADING & NORMALIZATION
# =============================================================================

def load_loops_unified(path: Path) -> pd.DataFrame:
    """
    Normalize different HiChIP loop formats into a standard DataFrame:
        chr1, start1, end1, chr2, start2, end2, counts, score, n_reps
    
    Handles:
    - BEDPE with header (HiChIPper format)
    - BEDPE without header
    - LoopCatalog 4-column BED format
    """
    path = Path(path)
    ext = path.suffix.lower()

    if ext == ".bedpe":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep="\t", header=None)

    ncols = df.shape[1]

    # BEDPE with header
    if {"chr1", "start1", "end1", "chr2", "start2", "end2"}.issubset(df.columns):
        return _normalize_bedpe_with_header(df, path)

    # BEDPE without header (9 columns)
    if ncols == 9:
        df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "score", "n_reps"]
        return _finalize_loop_types(df, path)

    # BEDPE without header (8 columns, no score)
    if ncols == 8:
        df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "n_reps"]
        df["score"] = np.nan
        return _finalize_loop_types(df, path)

    # LoopCatalog 4-column BED
    if ncols == 4:
        return _parse_loopcatalog_bed(df, path)

    raise ValueError(f"Unknown loop format with {ncols} columns in {path}")


def _normalize_bedpe_with_header(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Normalize BEDPE with header (HMEC/HiChIPper format)."""
    coord_cols = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
    missing = [c for c in coord_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing coordinate columns in {path}: {missing}")

    if "PETs" in df.columns:
        counts = df["PETs"].astype(float)
    elif "weight" in df.columns:
        counts = df["weight"].astype(float)
    else:
        counts = pd.Series([np.nan] * len(df))

    out = df.copy()
    out["counts"] = counts
    
    if "n_reps" not in out.columns:
        out["n_reps"] = 1
    if "score" not in out.columns:
        out["score"] = np.nan

    return _finalize_loop_types(out, path)


def _parse_loopcatalog_bed(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Parse LoopCatalog 4-column .bed format."""
    df.columns = ["chr1", "start1", "end1", "anchor2"]
    
    parsed = df["anchor2"].apply(lambda x: parse_hichip_anchor2(x, str(path)))
    df["chr2"] = parsed.str[0]
    df["start2"] = parsed.str[1].astype(int)
    df["end2"] = parsed.str[2].astype(int)
    df["score"] = parsed.str[3].astype(float)
    df["counts"] = np.nan
    df["n_reps"] = 1
    
    return _finalize_loop_types(df, path)


def _finalize_loop_types(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Ensure final types and column presence for loop data."""
    out = df.copy()
    
    for col in ["counts", "score"]:
        if col not in out.columns:
            out[col] = np.nan
    if "n_reps" not in out.columns:
        out["n_reps"] = 1

    out["start1"] = out["start1"].astype(int)
    out["end1"] = out["end1"].astype(int)
    out["start2"] = out["start2"].astype(int)
    out["end2"] = out["end2"].astype(int)
    out["counts"] = out["counts"].astype(float)
    out["score"] = out["score"].astype(float)
    out["n_reps"] = out["n_reps"].astype(int)

    return out[["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "score", "n_reps"]]


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def find_hichip_loops_file(
    hichip_dir: Path,
    cell_type: str,
    experiment: str = "H3K27ac",
) -> Optional[Path]:
    """Find the loops file for a given cell type and experiment."""
    cell_type_dir = Path(hichip_dir) / cell_type
    
    if not cell_type_dir.exists():
        return None
    
    # Try common naming patterns
    candidates = [
        cell_type_dir / f"{cell_type}_{experiment}_loops_hg38.bedpe",
        cell_type_dir / f"{cell_type}_{experiment}_loops_hg38.bed",
    ]
    
    for c in candidates:
        if c.exists():
            return c

    # Fall back to glob
    pattern = str(cell_type_dir / "*_loops_hg38.bed*")
    matches = glob.glob(pattern)
    if matches:
        return Path(matches[0])

    return None


# =============================================================================
# ANCHOR BUILDING
# =============================================================================

def build_anchors(loops: pd.DataFrame) -> pd.DataFrame:
    """Build anchor DataFrame from loops (one row per anchor A/B)."""
    if "loop_id" not in loops.columns:
        loops = loops.copy()
        loops["loop_id"] = np.arange(len(loops))

    anchor_a = loops[["chr1", "start1", "end1", "loop_id", "counts", "score", "n_reps"]].copy()
    anchor_a.rename(columns={"chr1": "chr", "start1": "start", "end1": "end"}, inplace=True)
    anchor_a["anchor"] = "A"

    anchor_b = loops[["chr2", "start2", "end2", "loop_id", "counts", "score", "n_reps"]].copy()
    anchor_b.rename(columns={"chr2": "chr", "start2": "start", "end2": "end"}, inplace=True)
    anchor_b["anchor"] = "B"

    return pd.concat([anchor_a, anchor_b], ignore_index=True)


# =============================================================================
# OVERLAP COMPUTATION (requires pyranges)
# =============================================================================

def compute_hichip_overlaps(
    loops: pd.DataFrame,
    ccres: pd.DataFrame,
    genes: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute anchor overlaps with cCREs and promoters.
    
    Requires pyranges to be installed.
    
    Returns:
        tuple: (anchor_metrics, overlap_ccre, overlap_prom)
    """
    if not HAS_PYRANGES:
        raise ImportError("pyranges is required for HiChIP overlap computation")
    
    loops = loops.copy()
    if "loop_id" not in loops.columns:
        loops["loop_id"] = np.arange(len(loops))
    
    anchors = build_anchors(loops)
    
    # Separate metrics (contains NaN/inf that PyRanges can't handle)
    anchor_metrics = anchors[["loop_id", "anchor", "counts", "score", "n_reps"]].drop_duplicates()
    anchors_for_pr = anchors[["chr", "start", "end", "loop_id", "anchor"]].copy()
    
    # Prepare promoters
    promoters = genes[["chrom", "prom_start", "prom_end", "gene_id", "gene_name"]].copy()
    promoters.rename(columns={"prom_start": "start", "prom_end": "end"}, inplace=True)
    
    # Prepare cCREs
    ccre_for_pr = ccres[["chrom", "start", "end", "cCRE_id"]].copy()
    
    # Convert to PyRanges format
    anchors_pr = anchors_for_pr.rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"})
    promoters_pr = promoters.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    ccre_pr = ccre_for_pr.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    
    # Compute overlaps
    gr_anchors = pr.PyRanges(anchors_pr)
    gr_promoters = pr.PyRanges(promoters_pr)
    gr_ccres = pr.PyRanges(ccre_pr)
    
    overlap_prom = gr_anchors.join(gr_promoters).df
    overlap_ccre = gr_anchors.join(gr_ccres).df
    
    # Rejoin metrics
    overlap_prom = overlap_prom.merge(anchor_metrics, on=["loop_id", "anchor"], how="left")
    overlap_ccre = overlap_ccre.merge(anchor_metrics, on=["loop_id", "anchor"], how="left")
    
    return anchor_metrics, overlap_ccre, overlap_prom


# =============================================================================
# DICT BUILDERS
# =============================================================================

def build_ccre_hichip_dict(
    overlap_ccre: pd.DataFrame,
    loops: pd.DataFrame,
    *,
    max_loops_per_ccre: Optional[int] = None,
    max_partner_ccres_per_loop: Optional[int] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Build per-cCRE HiChIP summary.
    
    Returns:
        {cCRE_id: {n_loops, max_counts, max_score, loops: [...]}}
    """
    if "loop_id" not in loops.columns:
        loops = loops.copy()
        loops["loop_id"] = np.arange(len(loops))
    
    # Build loop coordinates lookup
    loop_coords = (
        loops.set_index("loop_id")
        .apply(lambda r: {
            "anchor_a": {"chr": r["chr1"], "start": int(r["start1"]), "end": int(r["end1"])},
            "anchor_b": {"chr": r["chr2"], "start": int(r["start2"]), "end": int(r["end2"])},
        }, axis=1)
        .to_dict()
    )
    
    # Build loop -> cCREs mapping
    loop_to_ccres = (
        overlap_ccre.groupby("loop_id")["cCRE_id"]
        .apply(lambda s: set(s.astype(str)))
        .to_dict()
    )
    
    per_ccre = {}
    
    for ccre_id, df_cc in overlap_ccre.groupby("cCRE_id"):
        seen_loops = set()
        loops_list = []
        max_counts = None
        max_score = None
        ccre_id_str = str(ccre_id)
        
        for _, row in df_cc.iterrows():
            lid = row["loop_id"]
            if lid in seen_loops:
                continue
            seen_loops.add(lid)
            
            cnt = safe_float(row["counts"])
            scr = safe_float(row["score"])
            
            if cnt is not None:
                max_counts = cnt if max_counts is None else max(max_counts, cnt)
            if scr is not None:
                max_score = scr if max_score is None else max(max_score, scr)
            
            partners = sorted(loop_to_ccres.get(lid, set()) - {ccre_id_str})
            if max_partner_ccres_per_loop is not None and int(max_partner_ccres_per_loop) > 0:
                partners = partners[: int(max_partner_ccres_per_loop)]

            base = {
                "loop_id": int(lid),
                "counts": cnt,
                "score": scr,
                "n_reps": int(row["n_reps"]),
                "partner_cCREs": partners,
            }
            
            coords = loop_coords.get(lid)
            if coords:
                base.update(coords)
            
            loops_list.append(base)

        if max_loops_per_ccre is not None and int(max_loops_per_ccre) > 0 and len(loops_list) > int(max_loops_per_ccre):
            # Prefer higher-confidence loops when trimming.
            def _rank_key(d: Dict[str, Any]) -> Tuple[float, float, int]:
                s = d.get("score")
                c = d.get("counts")
                s_v = float(s) if s is not None and s == s else float("-inf")
                c_v = float(c) if c is not None and c == c else float("-inf")
                return (s_v, c_v, int(d.get("loop_id", -1)))

            loops_list = sorted(loops_list, key=_rank_key, reverse=True)[: int(max_loops_per_ccre)]
            # Recompute aggregates on the retained subset (keeps max_* consistent with stored loops)
            max_counts = None
            max_score = None
            for d in loops_list:
                cnt = d.get("counts")
                scr = d.get("score")
                if cnt is not None and cnt == cnt:
                    max_counts = float(cnt) if max_counts is None else max(max_counts, float(cnt))
                if scr is not None and scr == scr:
                    max_score = float(scr) if max_score is None else max(max_score, float(scr))

        if loops_list:
            per_ccre[ccre_id_str] = {
                "n_loops": len(loops_list),
                "max_counts": max_counts,
                "max_score": max_score,
                "loops": loops_list,
            }
    
    return per_ccre


def build_ccre_gene_hichip_map(
    overlap_ccre: pd.DataFrame,
    overlap_prom: pd.DataFrame,
    *,
    max_loops_per_gene: Optional[int] = None,
) -> Dict[str, Dict[str, List[Dict]]]:
    """
    Build cCRE -> gene -> [loops] mapping for opposite-anchor pairs.
    
    Returns:
        {cCRE_id: {gene_name: [{loop_id, counts, score, ...}, ...]}}
    """
    cc = overlap_ccre.rename(columns={"anchor": "anchor_ccre"})
    cp = overlap_prom.rename(columns={"anchor": "anchor_prom"})
    
    merged = cc.merge(cp, on="loop_id", suffixes=("_ccre", "_prom"))
    # Keep only opposite-anchor pairs (cCRE on A, gene on B or vice versa)
    merged_ep = merged[merged["anchor_ccre"] != merged["anchor_prom"]]
    
    cg = defaultdict(lambda: defaultdict(list))
    for _, r in merged_ep.iterrows():
        cg[r["cCRE_id"]][r["gene_name"]].append(
            {
                "loop_id": r["loop_id"],
                "counts": safe_float(r["counts_ccre"]),
                "score": safe_float(r["score_ccre"]),
                "n_reps": r["n_reps_ccre"],
                "anchor_ccre": r["anchor_ccre"],
                "anchor_prom": r["anchor_prom"],
            }
        )

    if max_loops_per_gene is not None and int(max_loops_per_gene) > 0:
        for ccre_id, gmap in list(cg.items()):
            for gname, lst in list(gmap.items()):
                if len(lst) <= int(max_loops_per_gene):
                    continue

                def _rank(d: Dict[str, Any]) -> Tuple[float, float, int]:
                    s = d.get("score")
                    c = d.get("counts")
                    s_v = float(s) if s is not None and s == s else float("-inf")
                    c_v = float(c) if c is not None and c == c else float("-inf")
                    return (s_v, c_v, int(d.get("loop_id", -1)))

                gmap[gname] = sorted(lst, key=_rank, reverse=True)[: int(max_loops_per_gene)]

    return dict(cg)


def _hichip_limits_from_env() -> Dict[str, Any]:
    """
    Memory-control knobs for HiChIP link materialization.

    By default (APM_HICHIP_SAFE_DEFAULTS=1), conservative caps apply unless you set
    explicit APM_HICHIP_MAX_* overrides. Set APM_HICHIP_SAFE_DEFAULTS=0 to disable
    those defaults (not recommended on memory-limited hosts).

    Env vars:
    - APM_HICHIP_LOOP_CHUNK_ROWS: chunk size when splitting huge loop files (default 250000)
    - APM_HICHIP_MAX_LOOPS_PER_CCRE: cap loops stored per cCRE in per_ccre dict
    - APM_HICHIP_MAX_LOOPS_PER_GENE: cap loops stored per (cCRE,gene) in ccre_gene_map
    - APM_HICHIP_MAX_PARTNER_CCRES: cap partner_cCREs list length per loop record
    - APM_HICHIP_SAFE_DEFAULTS: if "1" (default), apply conservative defaults when the explicit caps above are unset
    - APM_HICHIP_STREAM_MIN_BYTES: if file size exceeds this, stream loop files in chunks (default 50_000_000)
    """
    def _i(name: str) -> Optional[int]:
        v = os.environ.get(name, "").strip()
        if not v:
            return None
        try:
            x = int(v)
            return x if x > 0 else None
        except Exception:
            return None

    safe = os.environ.get("APM_HICHIP_SAFE_DEFAULTS", "1").strip().lower() not in ("0", "false", "no")

    max_ccre = _i("APM_HICHIP_MAX_LOOPS_PER_CCRE")
    max_gene = _i("APM_HICHIP_MAX_LOOPS_PER_GENE")
    max_part = _i("APM_HICHIP_MAX_PARTNER_CCRES")

    if safe:
        # These defaults are intentionally conservative: they preserve the strongest loop records
        # while preventing multi-million-row explosions in RAM during dict materialization.
        if max_ccre is None:
            max_ccre = 200
        if max_gene is None:
            max_gene = 50
        if max_part is None:
            max_part = 200

    return {
        "loop_chunk_rows": _i("APM_HICHIP_LOOP_CHUNK_ROWS") or 250_000,
        "max_loops_per_ccre": max_ccre,
        "max_loops_per_gene": max_gene,
        "max_partner_ccres": max_part,
        "stream_min_bytes": _i("APM_HICHIP_STREAM_MIN_BYTES") or 50_000_000,
        "force_stream": os.environ.get("APM_HICHIP_FORCE_STREAM", "").strip().lower() in ("1", "true", "yes"),
    }


def _bedpe_reader_kwargs(path: Path) -> Dict[str, Any]:
    """
    BEDPE files may be headered (HiChIPper) or headerless (9/8 columns). Choose a stable reader mode.
    """
    p = Path(path)
    if p.suffix.lower() != ".bedpe":
        return {"sep": "\t", "header": None}
    peek = pd.read_csv(p, sep="\t", nrows=5)
    cols = set(peek.columns.astype(str))
    if {"chr1", "start1", "end1", "chr2", "start2", "end2"}.issubset(cols):
        return {"sep": "\t", "header": 0}
    return {"sep": "\t", "header": None}


def _merge_ccre_hichip_dict(dst: Dict[str, Dict[str, Any]], src: Dict[str, Dict[str, Any]], *, max_loops_per_ccre: Optional[int]) -> None:
    for ccre_id, rec in src.items():
        if ccre_id not in dst:
            dst[ccre_id] = rec
            continue
        cur = dst[ccre_id]
        nxt_loops = (cur.get("loops") or []) + (rec.get("loops") or [])
        # de-dup by loop_id (first wins)
        seen = set()
        dedup = []
        for d in nxt_loops:
            lid = int(d.get("loop_id"))
            if lid in seen:
                continue
            seen.add(lid)
            dedup.append(d)

        if max_loops_per_ccre is not None and int(max_loops_per_ccre) > 0 and len(dedup) > int(max_loops_per_ccre):

            def _rank(d: Dict[str, Any]) -> Tuple[float, float, int]:
                s = d.get("score")
                c = d.get("counts")
                s_v = float(s) if s is not None and s == s else float("-inf")
                c_v = float(c) if c is not None and c == c else float("-inf")
                return (s_v, c_v, int(d.get("loop_id", -1)))

            dedup = sorted(dedup, key=_rank, reverse=True)[: int(max_loops_per_ccre)]

        max_counts = cur.get("max_counts")
        max_score = cur.get("max_score")
        for d in dedup:
            cnt = d.get("counts")
            scr = d.get("score")
            if cnt is not None and cnt == cnt:
                max_counts = float(cnt) if max_counts is None else max(float(max_counts), float(cnt))
            if scr is not None and scr == scr:
                max_score = float(scr) if max_score is None else max(float(max_score), float(scr))

        dst[ccre_id] = {"n_loops": len(dedup), "max_counts": max_counts, "max_score": max_score, "loops": dedup}


def _merge_ccre_gene_map(
    dst: Dict[str, Dict[str, List[Dict]]],
    src: Dict[str, Dict[str, List[Dict]]],
    *,
    max_loops_per_gene: Optional[int],
) -> None:
    for ccre_id, gmap in src.items():
        d0 = dst.setdefault(ccre_id, {})
        for gname, lst in gmap.items():
            d0.setdefault(gname, [])
            d0[gname].extend(lst)
            if max_loops_per_gene is not None and int(max_loops_per_gene) > 0 and len(d0[gname]) > int(max_loops_per_gene):

                def _rank(d: Dict[str, Any]) -> Tuple[float, float, int]:
                    s = d.get("score")
                    c = d.get("counts")
                    s_v = float(s) if s is not None and s == s else float("-inf")
                    c_v = float(c) if c is not None and c == c else float("-inf")
                    return (s_v, c_v, int(d.get("loop_id", -1)))

                d0[gname] = sorted(d0[gname], key=_rank, reverse=True)[: int(max_loops_per_gene)]


def _process_one_cell_type_hichip(
    *,
    loops: pd.DataFrame,
    ccres: pd.DataFrame,
    genes: pd.DataFrame,
    max_loops_per_ccre: Optional[int],
    max_loops_per_gene: Optional[int],
    max_partner_ccres: Optional[int],
) -> Dict[str, Any]:
    _, overlap_ccre, overlap_prom = compute_hichip_overlaps(loops, ccres, genes)
    per_ccre = build_ccre_hichip_dict(
        overlap_ccre,
        loops,
        max_loops_per_ccre=max_loops_per_ccre,
        max_partner_ccres_per_loop=max_partner_ccres,
    )
    ccre_gene_map = build_ccre_gene_hichip_map(overlap_ccre, overlap_prom, max_loops_per_gene=max_loops_per_gene)
    return {"per_ccre": per_ccre, "ccre_gene_map": ccre_gene_map}


def _iter_loop_chunks(loops_path: Path, chunk_rows: int) -> Iterable[pd.DataFrame]:
    """
    Stream a loops file in chunks without loading the full table at once.
    """
    path = Path(loops_path)
    ext = path.suffix.lower()
    if ext == ".bedpe":
        kw = _bedpe_reader_kwargs(path)
        reader = pd.read_csv(path, chunksize=int(chunk_rows), **kw)
    else:
        reader = pd.read_csv(path, sep="\t", header=None, chunksize=int(chunk_rows))

    offset = 0
    for chunk in reader:
        df = chunk
        # Reuse the same normalization logic as load_loops_unified, but on a chunk.
        # We mimic load_loops_unified by writing to a temp-ish path object for error messages only.
        if ext == ".bedpe" and {"chr1", "start1", "end1", "chr2", "start2", "end2"}.issubset(df.columns):
            df = _normalize_bedpe_with_header(df, path)
        else:
            ncols = df.shape[1]
            if ncols == 9:
                df = df.copy()
                df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "score", "n_reps"]
                df = _finalize_loop_types(df, path)
            elif ncols == 8:
                df = df.copy()
                df.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "counts", "n_reps"]
                df["score"] = np.nan
                df = _finalize_loop_types(df, path)
            elif ncols == 4:
                df = _parse_loopcatalog_bed(df, path)
            else:
                raise ValueError(f"Unknown loop format with {ncols} columns in {path} (chunked read)")

        df = df.copy()
        df["loop_id"] = np.arange(offset, offset + len(df), dtype=int)
        offset += len(df)
        yield df


# =============================================================================
# PUBLIC API
# =============================================================================

def build_hichip_links(
    hichip_dir: Path,
    cell_types: List[str],
    ccres: pd.DataFrame,
    genes: pd.DataFrame,
    experiment: str = "H3K27ac",
    *,
    loop_chunk_rows: Optional[int] = None,
    max_loops_per_ccre: Optional[int] = None,
    max_loops_per_gene: Optional[int] = None,
    max_partner_ccres: Optional[int] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Build HiChIP links for multiple cell types.
    
    Returns:
        {canonical_cell_type: {cCRE_id: {n_loops, max_counts, max_score, loops: [...]}}}
        Keys are from ``canonical_hichip_output_key`` (usually the same as the on-disk folder name).
    """
    if not HAS_PYRANGES:
        print("Warning: pyranges not installed. Skipping HiChIP integration.")
        return {}
    
    print("Building HiChIP links...")
    all_hichip = {}

    env = _hichip_limits_from_env()
    loop_chunk_rows = int(loop_chunk_rows or env["loop_chunk_rows"] or 250_000)
    max_loops_per_ccre = max_loops_per_ccre if max_loops_per_ccre is not None else env["max_loops_per_ccre"]
    max_loops_per_gene = max_loops_per_gene if max_loops_per_gene is not None else env["max_loops_per_gene"]
    max_partner_ccres = max_partner_ccres if max_partner_ccres is not None else env["max_partner_ccres"]

    stream_min_bytes = int(env.get("stream_min_bytes") or 50_000_000)
    force_stream = bool(env.get("force_stream"))

    print(
        "  HiChIP memory controls: "
        f"loop_chunk_rows={loop_chunk_rows}, "
        f"stream_if_bytes>={stream_min_bytes}, "
        f"force_stream={force_stream}, "
        f"always_stream_loopcatalog_bed=True, "
        f"max_loops_per_ccre={max_loops_per_ccre}, "
        f"max_loops_per_gene={max_loops_per_gene}, "
        f"max_partner_ccres={max_partner_ccres}"
    )
    
    for cell_type in cell_types:
        loops_path = find_hichip_loops_file(hichip_dir, cell_type, experiment)
        out_key = canonical_hichip_output_key(cell_type)
        
        if loops_path is None:
            print(f"  {cell_type}: No loops file found, skipping")
            continue
        
        print(f"  {cell_type}: Loading {loops_path} (output key: {out_key!r})")
        
        try:
            # Estimate whether we should stream. If unknown, load once and fall back.
            # For very large files, chunked overlap avoids holding the full loops table + overlaps at peak.
            per_ccre_acc: Dict[str, Dict[str, Any]] = {}
            ccre_gene_acc: Dict[str, Dict[str, List[Dict]]] = {}

            ext_l = loops_path.suffix.lower()
            try:
                sz = int(loops_path.stat().st_size)
            except Exception:
                sz = 0
            # LoopCatalog `.bed` inputs are often far larger in rows than BEDPE for similar bytes;
            # loading them whole is the common OOM path, so always chunk-read `.bed` loops here.
            use_stream = bool(force_stream or sz >= stream_min_bytes or ext_l == ".bed")

            if use_stream:
                print(f"    Streaming loops in chunks of {loop_chunk_rows} rows (memory-safe mode)")
                for loops in _iter_loop_chunks(loops_path, loop_chunk_rows):
                    part = _process_one_cell_type_hichip(
                        loops=loops,
                        ccres=ccres,
                        genes=genes,
                        max_loops_per_ccre=max_loops_per_ccre,
                        max_loops_per_gene=max_loops_per_gene,
                        max_partner_ccres=max_partner_ccres,
                    )
                    _merge_ccre_hichip_dict(per_ccre_acc, part["per_ccre"], max_loops_per_ccre=max_loops_per_ccre)
                    _merge_ccre_gene_map(ccre_gene_acc, part["ccre_gene_map"], max_loops_per_gene=max_loops_per_gene)
                all_hichip[out_key] = {"per_ccre": per_ccre_acc, "ccre_gene_map": ccre_gene_acc}
                print(f"    Found {len(per_ccre_acc)} cCREs with loops (chunked)")
            else:
                loops = load_loops_unified(loops_path)
                loops["loop_id"] = np.arange(len(loops))
                part = _process_one_cell_type_hichip(
                    loops=loops,
                    ccres=ccres,
                    genes=genes,
                    max_loops_per_ccre=max_loops_per_ccre,
                    max_loops_per_gene=max_loops_per_gene,
                    max_partner_ccres=max_partner_ccres,
                )
                all_hichip[out_key] = part
                print(f"    Found {len(part.get('per_ccre', {}))} cCREs with loops")
                del loops

            gc.collect()

        except Exception as e:
            print(f"  {cell_type}: Error processing - {e}")
            continue
    
    return all_hichip


def integrate_hichip_to_element_table(
    elem_focus: pd.DataFrame,
    hichip_data: Dict[str, Dict[str, Any]],
) -> pd.DataFrame:
    """
    Add HiChIP data to element table.
    
    Adds per-cell-type HiChIP columns and updates gene_links.
    """
    elem_focus = elem_focus.copy()
    
    for cell_type, ct_data in hichip_data.items():
        out_key = canonical_hichip_output_key(cell_type)
        per_ccre = ct_data.get("per_ccre", {})
        ccre_gene_map = ct_data.get("ccre_gene_map", {})
        
        # Initialize cell type column if needed
        if out_key not in elem_focus.columns:
            elem_focus[out_key] = [{} for _ in range(len(elem_focus))]
        
        # Add hichip to cell type column
        def add_hichip_to_cell_col(row):
            cell_dict = row[out_key] if isinstance(row[out_key], dict) else {}
            ccre_id = str(row.get("cCRE_id", row.get("ENCODE_id", "")))
            cell_dict["hichip"] = per_ccre.get(ccre_id, {
                "n_loops": 0, "max_counts": None, "max_score": None, "loops": []
            })
            return cell_dict
        
        elem_focus[out_key] = elem_focus.apply(add_hichip_to_cell_col, axis=1)
        
        # Update gene_links with HiChIP evidence
        def update_gene_links_hichip(row):
            gene_links = row.get("gene_links", {})
            if not isinstance(gene_links, dict):
                return gene_links
            
            ccre_id = str(row.get("cCRE_id", row.get("ENCODE_id", "")))
            per_gene_loops = ccre_gene_map.get(ccre_id, {})
            
            for gene_name, gdict in gene_links.items():
                hichip = gdict.setdefault("hichip", {})
                loops_for_pair = per_gene_loops.get(gene_name, [])
                
                if not loops_for_pair:
                    continue
                
                counts_vals = [l.get("counts") for l in loops_for_pair if l.get("counts") is not None]
                score_vals = [l.get("score") for l in loops_for_pair if l.get("score") is not None]
                
                hichip[out_key] = {
                    "n_loops": len(loops_for_pair),
                    "max_counts": max(counts_vals) if counts_vals else None,
                    "max_score": max(score_vals) if score_vals else None,
                    "loops": loops_for_pair,
                }
            
            return gene_links
        
        if "gene_links" in elem_focus.columns:
            elem_focus["gene_links"] = elem_focus.apply(update_gene_links_hichip, axis=1)
    
    return elem_focus
