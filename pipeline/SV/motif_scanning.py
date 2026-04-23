"""
Motif scanning integration with FIMO.

Functions for:
- Extracting selected motifs from MEME databases
- Parsing FIMO output
- Recombining motif hits with SV tables
- BND neojunction motif annotation (gained cCRE motifs from rearranged adjacency)
"""

import ast
import gc
import heapq
import os
import re
import shutil
import traceback
from io import StringIO
from collections import defaultdict
from pathlib import Path
from typing import List, Set, Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline.config import THRESHOLDS


# =============================================================================
# MOTIF EXTRACTION
# =============================================================================

def extract_selected_motifs(
    input_meme: Path,
    output_meme: Path,
    target_symbols: Optional[List[str]] = None,
) -> Tuple[int, List[str]]:
    """
    Extract motifs matching target TF symbols from a MEME file.
    
    Args:
        input_meme: Path to input MEME file
        output_meme: Path to output MEME file
        target_symbols: List of TF symbols to match (case-insensitive substring)
    
    Returns:
        Tuple of (n_selected, missing_symbols)
    """
    if target_symbols is None:
        from pipeline.config import SV_TARGET_TF_SYMBOLS

        target_symbols = list(SV_TARGET_TF_SYMBOLS)
    
    target_patterns = [re.compile(sym, re.IGNORECASE) for sym in target_symbols]
    
    def motif_matches_any(motif_name: str) -> bool:
        return any(pat.search(motif_name) for pat in target_patterns)
    
    input_meme = Path(input_meme)
    output_meme = Path(output_meme)
    
    in_text = input_meme.read_text()
    
    # Split header vs motifs
    parts = in_text.split("\nMOTIF ")
    if len(parts) == 1:
        raise ValueError("Could not find any 'MOTIF ' blocks in the input MEME file.")
    
    header = parts[0]
    motif_blocks_raw = parts[1:]
    
    selected_blocks = []
    found_symbols = set()
    
    for block_raw in motif_blocks_raw:
        block = "MOTIF " + block_raw
        first_line = block.splitlines()[0]
        tokens = first_line.split()
        
        if len(tokens) < 2:
            continue
        
        motif_name = tokens[1]
        
        if motif_matches_any(motif_name):
            selected_blocks.append(block)
            for sym in target_symbols:
                if re.search(sym, motif_name, re.IGNORECASE):
                    found_symbols.add(sym)
    
    # Build output
    out_lines = [header.rstrip(), ""]
    for block in selected_blocks:
        out_lines.append(block.rstrip())
        out_lines.append("")
    
    out_text = "\n".join(out_lines).rstrip() + "\n"
    
    output_meme.parent.mkdir(parents=True, exist_ok=True)
    output_meme.write_text(out_text)
    
    missing = [sym for sym in target_symbols if sym not in found_symbols]
    
    print(f"Selected {len(selected_blocks)} motif blocks into {output_meme}")
    if missing:
        print(f"No motif found for: {', '.join(missing)}")
    
    return len(selected_blocks), missing


# =============================================================================
# FIMO OUTPUT PARSING
# =============================================================================

def _fimo_field_key(name: str) -> str:
    s = str(name).strip().lstrip("\ufeff").lower()
    return s.replace(" ", "_").replace("-", "_")


def _parse_fimo_header_tokens(line: str) -> Optional[List[str]]:
    """
    If ``line`` is a FIMO TSV header row, return stripped column names; else None.

    Accepts MEME 5.5+ banner-prefixed headers like ``# motif_id<TAB>...`` and plain headers.
    """
    s = line.strip().lstrip("\ufeff")
    if not s or "\t" not in s:
        return None
    # Header rows may start with '#'; strip one leading '#' and whitespace for parsing only
    s_body = s[1:].lstrip() if s.startswith("#") else s
    toks = [t.strip() for t in s_body.split("\t")]
    if not any(toks):
        return None
    keys = {_fimo_field_key(t) for t in toks if t}
    if "start" in keys or "begin" in keys:
        return toks
    if "motif_id" in keys and ("sequence_name" in keys or "sequence" in keys):
        return toks
    return None


def _load_fimo_tsv_as_dataframe(fimo_tsv_path: Path) -> pd.DataFrame:
    """
    Read FIMO tab output into a DataFrame.

    Handles: UTF-8 BOM, leading ``#`` banner lines (FIMO 5.5.x), ``# motif_id<TAB>…`` headers,
    no-hit outputs, and minor column-name variants across MEME suite versions.
    """
    path = Path(fimo_tsv_path)

    def _read_csv_file(**kwargs: Any) -> pd.DataFrame:
        try:
            return pd.read_csv(path, sep="\t", encoding="utf-8-sig", **kwargs)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
    except OSError:
        return pd.DataFrame()

    lines = text.splitlines()
    header_idx: Optional[int] = None
    header_tokens: Optional[List[str]] = None
    for i, line in enumerate(lines[:8000]):
        toks = _parse_fimo_header_tokens(line)
        if toks:
            header_idx = i
            header_tokens = toks
            break

    if header_idx is not None and header_tokens:
        body = "\n".join(lines[header_idx + 1 :])
        if not body.strip():
            df = pd.DataFrame(columns=header_tokens)
        else:
            read_kw: Dict[str, Any] = dict(
                sep="\t",
                names=header_tokens,
                engine="python",
                comment="#",
            )
            # pandas >= 1.3
            read_kw["on_bad_lines"] = "skip"
            try:
                df = pd.read_csv(StringIO(body), **read_kw)
            except TypeError:
                read_kw.pop("on_bad_lines", None)
                df = pd.read_csv(StringIO(body), **read_kw)
    else:
        df = _read_csv_file(comment="#", engine="python")
        if "start" not in df.columns and len(df.columns) == 1:
            lone = str(df.columns[0])
            if lone.lstrip("\ufeff").lstrip().startswith("#") and "\t" not in lone:
                # Banner-only parse (e.g. BOM + single-column junk); treat as empty table
                df = pd.DataFrame()

    if df.empty and len(df.columns) == 0:
        return df

    # Strip BOM / whitespace from column labels
    df = df.rename(columns={c: str(c).strip().lstrip("\ufeff") for c in df.columns})

    # Map common aliases to names expected below
    lower_map: Dict[str, str] = {}
    for c in df.columns:
        lower_map[_fimo_field_key(c)] = c

    def pick_col(*candidates: str) -> Optional[str]:
        for raw in candidates:
            k = _fimo_field_key(raw)
            if k in lower_map:
                return lower_map[k]
        return None

    rename_map: Dict[str, str] = {}
    seq_c = pick_col("sequence_name", "chrom", "chr", "seq_name", "sequence")
    start_c = pick_col("start", "begin")
    stop_c = pick_col("stop", "end")
    strand_c = pick_col("strand")
    motif_c = pick_col("motif_id", "motif_name", "motif")
    p_c = pick_col("p-value", "p_value", "pvalue")
    q_c = pick_col("q-value", "q_value", "qvalue")
    score_c = pick_col("score")

    if seq_c and _fimo_field_key(seq_c) != "sequence_name":
        rename_map[seq_c] = "sequence_name"
    if start_c and start_c != "start":
        rename_map[start_c] = "start"
    if stop_c and stop_c != "stop":
        rename_map[stop_c] = "stop"
    if strand_c and strand_c != "strand":
        rename_map[strand_c] = "strand"
    if motif_c and motif_c != "motif_id":
        rename_map[motif_c] = "motif_id"
    if p_c and p_c != "p-value":
        rename_map[p_c] = "p-value"
    if q_c and q_c != "q-value":
        rename_map[q_c] = "q-value"
    if score_c and score_c != "score":
        rename_map[score_c] = "score"

    if rename_map:
        df = df.rename(columns=rename_map)

    return df


def parse_fimo_tsv_to_bed(
    fimo_tsv_path: Path,
    output_bed_path: Path,
) -> pd.DataFrame:
    """
    Convert FIMO TSV output to BED format.
    
    Args:
        fimo_tsv_path: Path to FIMO TSV file
        output_bed_path: Path to output BED file
    
    Returns:
        BED DataFrame
    """
    fimo_tsv_path = Path(fimo_tsv_path)
    output_bed_path = Path(output_bed_path)

    df = _load_fimo_tsv_as_dataframe(fimo_tsv_path)
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)

    if "start" not in df.columns and len(df) > 0:
        raise ValueError(
            f"FIMO TSV {fimo_tsv_path} has data rows but no start/begin column. "
            f"Columns: {list(df.columns)}"
        )

    if df.empty or "start" not in df.columns:
        bed = pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "motif_id",
                "p_value",
                "q_value",
                "score",
                "strand",
            ]
        )
        bed.to_csv(output_bed_path, sep="\t", header=False, index=False)
        return bed

    if "stop" not in df.columns:
        raise ValueError(
            f"FIMO TSV {fimo_tsv_path} has 'start' but no 'stop'/'end' column. "
            f"Columns: {list(df.columns)}"
        )
    if "sequence_name" not in df.columns:
        raise ValueError(
            f"FIMO TSV {fimo_tsv_path} missing sequence/chrom column. "
            f"Columns: {list(df.columns)}"
        )

    # Coerce start/stop to numeric
    df["start_num"] = pd.to_numeric(df["start"], errors="coerce")
    df["stop_num"] = pd.to_numeric(df["stop"], errors="coerce")
    df = df.dropna(subset=["start_num", "stop_num"]).copy()

    if df.empty:
        bed = pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "motif_id",
                "p_value",
                "q_value",
                "score",
                "strand",
            ]
        )
        bed.to_csv(output_bed_path, sep="\t", header=False, index=False)
        return bed

    # Clean chromosome names
    chrom = df["sequence_name"].astype(str)
    chrom = chrom.str.replace(r"\.0$", "", regex=True)
    chrom = chrom.str.replace(r"^chr", "", regex=True)

    motif_col = (
        df["motif_id"]
        if "motif_id" in df.columns
        else pd.Series(".", index=df.index, dtype=object)
    )
    strand_col = (
        df["strand"]
        if "strand" in df.columns
        else pd.Series("+", index=df.index, dtype=object)
    )
    p_col = (
        df["p-value"]
        if "p-value" in df.columns
        else pd.Series(np.nan, index=df.index, dtype=float)
    )

    bed = pd.DataFrame({
        "chrom": chrom,
        "start": (df["start_num"].astype(int) - 1).clip(lower=0),
        "end": df["stop_num"].astype(int),
        "motif_id": motif_col,
        "p_value": p_col,
        "q_value": df["q-value"] if "q-value" in df.columns else np.nan,
        "score": df["score"] if "score" in df.columns else np.nan,
        "strand": strand_col,
    })
    bed = bed.drop_duplicates(ignore_index=True)

    bed.to_csv(output_bed_path, sep="\t", header=False, index=False)

    return bed


# =============================================================================
# MOTIF RECOMBINATION
# =============================================================================

def _parse_elem_hits_cell(x) -> List[Dict]:
    """Parse elem_hits cell into list of dicts."""
    if isinstance(x, list):
        return x
    if pd.isna(x) or x in ("", "[]"):
        return []
    try:
        return ast.literal_eval(x)
    except Exception:
        return []


def _parse_flank_side(name_tail: str) -> str:
    """Extract flank side from name tail."""
    parts = name_tail.split("_")
    if len(parts) >= 3 and parts[-2] == "flank":
        return parts[-1]
    return parts[-1]


def _extract_elem_id_from_name(name: str) -> Optional[str]:
    """Extract elem_id from BED name."""
    for part in name.split("|"):
        if part.startswith("elem:"):
            return part[len("elem:"):]
    return None


def _fimo_recombine_chunk_rows() -> int:
    """Rows per ``read_csv`` chunk when scanning merged FIMO–SV BEDs (``THRESHOLDS.sv_fimo_recombine_chunk_rows``)."""
    return max(25_000, min(2_000_000, int(THRESHOLDS.sv_fimo_recombine_chunk_rows)))


def _fimo_max_flank_hits_per_tf() -> int:
    """Max motif hits per (SV id, TF) for flank intervals (``THRESHOLDS.sv_fimo_max_flank_hits_per_tf``)."""
    return min(10_000, max(20, int(THRESHOLDS.sv_fimo_max_flank_hits_per_tf)))


def _fimo_max_elem_hits_per_tf() -> int:
    """Max motif hits per (SV id, elem_id, TF) (``THRESHOLDS.sv_fimo_max_elem_hits_per_tf``)."""
    return min(10_000, max(20, int(THRESHOLDS.sv_fimo_max_elem_hits_per_tf)))


def _fimo_recombine_stream_if_merged_bed_bytes() -> int:
    """Merged BED size at which STEP 5 streams the SV CSV (``THRESHOLDS.sv_fimo_recombine_stream_if_merged_bed_bytes``)."""
    return max(0, int(THRESHOLDS.sv_fimo_recombine_stream_if_merged_bed_bytes))


def _fimo_recombine_sv_csv_chunk_rows() -> int:
    """SV rows per ``read_csv`` chunk when streaming STEP 5 (``THRESHOLDS.sv_fimo_recombine_sv_csv_chunk_rows``)."""
    return max(100, min(50_000, int(THRESHOLDS.sv_fimo_recombine_sv_csv_chunk_rows)))


def _prepare_sv_df_for_motif_attach(df: pd.DataFrame) -> None:
    if "flank_motif_hits" not in df.columns:
        df["flank_motif_hits"] = [[] for _ in range(len(df))]
    else:
        df["flank_motif_hits"] = df["flank_motif_hits"].apply(
            lambda x: [] if (pd.isna(x) or x in ("", "[]")) else (
                x if isinstance(x, list) else ast.literal_eval(x)
            )
        )
    df["elem_hits_parsed"] = df["elem_hits"].apply(_parse_elem_hits_cell)


def _build_sv_id_to_indices_from_df(df: pd.DataFrame) -> Dict[str, List[int]]:
    id_to_indices: Dict[str, List[int]] = defaultdict(list)
    for idx, sv_id in df["id"].items():
        id_to_indices[str(sv_id)].append(int(idx))
    return id_to_indices


def _build_sv_id_to_indices_ids_only(sv_csv_path: Path) -> Dict[str, List[int]]:
    ids = pd.read_csv(sv_csv_path, usecols=["id"], dtype={"id": str})
    id_to_indices: Dict[str, List[int]] = defaultdict(list)
    for idx, sv_id in enumerate(ids["id"].astype(str)):
        id_to_indices[sv_id].append(idx)
    return id_to_indices


def _attach_flank_and_elem_motif_hits(
    s: pd.DataFrame,
    id_to_indices: Dict[str, List[int]],
    flank_hits_by_sv: Dict[str, List[Dict[str, Any]]],
    elem_motif_hits: Dict[str, Dict[str, List[Dict[str, Any]]]],
) -> None:
    for sv_id, hits in flank_hits_by_sv.items():
        idxs = id_to_indices.get(sv_id, [])
        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None

            existing = s.at[idx, "flank_motif_hits"]
            if not isinstance(existing, list):
                existing = []

            hits_with_dist: List[Dict[str, Any]] = []
            for h in hits:
                h_copy = dict(h)
                if bp_pos is not None:
                    center = (h_copy["start"] + h_copy["end"]) // 2
                    h_copy["distance_to_pos"] = center - bp_pos
                else:
                    h_copy["distance_to_pos"] = None
                hits_with_dist.append(h_copy)

            existing.extend(hits_with_dist)
            s.at[idx, "flank_motif_hits"] = existing

    for sv_id, elem_dict in elem_motif_hits.items():
        idxs = id_to_indices.get(sv_id, [])
        for idx in idxs:
            try:
                bp_pos = int(s.at[idx, "pos"])
            except Exception:
                bp_pos = None

            elem_list = s.at[idx, "elem_hits_parsed"]
            if not isinstance(elem_list, list):
                continue

            for elem in elem_list:
                if not isinstance(elem, dict):
                    continue
                eid = elem.get("elem_id")
                if not eid:
                    continue

                hits_for_elem = elem_dict.get(eid)
                if not hits_for_elem:
                    elem["motif_hits"] = []
                    continue

                existing = elem.get("motif_hits", [])
                if not isinstance(existing, list):
                    existing = []

                hits_with_dist: List[Dict[str, Any]] = []
                for h in hits_for_elem:
                    h_copy = dict(h)
                    if bp_pos is not None:
                        center = (h_copy["start"] + h_copy["end"]) // 2
                        h_copy["distance_to_pos"] = center - bp_pos
                    else:
                        h_copy["distance_to_pos"] = None
                    hits_with_dist.append(h_copy)

                existing.extend(hits_with_dist)
                elem["motif_hits"] = existing


def _finalize_sv_motif_columns_for_csv(s: pd.DataFrame) -> None:
    s["elem_hits"] = s["elem_hits_parsed"].apply(repr)
    s.drop(columns=["elem_hits_parsed"], inplace=True)
    s["flank_motif_hits"] = s["flank_motif_hits"].apply(
        lambda x: "[]" if (x is None or (isinstance(x, float) and np.isnan(x)) or x == []) else repr(x)
    )


class _PValueHitHeap:
    """
    Keep up to ``max_hits`` motif records with the lowest p-values within one
    grouping key (e.g. one TF for flank or one TF for a regulatory element).

    bedtools can emit tens of millions of intersection lines; storing every dict
    in RAM is what blows memory. Heaps are keyed per TF separately so hits from
    different TFs do not compete for the same slots.
    """

    __slots__ = ("max_hits", "_h", "_seq")

    def __init__(self, max_hits: int) -> None:
        self.max_hits = max(1, max_hits)
        self._h: List[Tuple[float, int, Dict[str, Any]]] = []
        self._seq = 0

    def add(self, hit: Dict[str, Any]) -> None:
        p_raw = hit.get("p_value", float("nan"))
        try:
            p = float(p_raw)
        except Exception:
            p = float("nan")
        if p != p:  # NaN
            p = 1.0
        neg_p = -p
        self._seq += 1
        entry = (neg_p, self._seq, hit)
        if len(self._h) < self.max_hits:
            heapq.heappush(self._h, entry)
        elif neg_p > self._h[0][0]:
            heapq.heapreplace(self._h, entry)

    def to_sorted_hits(self) -> List[Dict[str, Any]]:
        hits = [t[2] for t in self._h]
        hits.sort(key=lambda d: float(d.get("p_value") or 1.0))
        return hits


def _safe_float_cell(x: Any) -> float:
    try:
        v = float(x)
        if v != v:  # NaN
            return float("nan")
        return v
    except Exception:
        return float("nan")


def _accumulate_merged_bed_chunk(
    chunk: pd.DataFrame,
    flank_heaps: Dict[str, Dict[str, _PValueHitHeap]],
    elem_heaps: Dict[str, Dict[str, Dict[str, _PValueHitHeap]]],
    k_flank_tf: int,
    k_elem_tf: int,
) -> None:
    """One pass over a merged bedtools chunk (no full-file copy)."""
    for tup in chunk.itertuples(index=False, name=None):
        (
            _sv_chrom,
            _sv_start,
            _sv_end,
            sv_name,
            _motif_chrom,
            motif_start,
            motif_end,
            motif_id,
            p_value,
            q_value,
            score,
            strand,
        ) = tup
        name = str(sv_name)
        sv_id = name.split("|")[0]
        motif_start_i = int(motif_start)
        motif_end_i = int(motif_end)
        motif_id_s = str(motif_id)
        p_val = _safe_float_cell(p_value)
        q_val = _safe_float_cell(q_value)
        score_f = _safe_float_cell(score)
        strand_s = str(strand)
        tf_short = motif_id_s.split(".", 1)[0] if "." in motif_id_s else motif_id_s
        hit_core = {
            "start": motif_start_i,
            "end": motif_end_i,
            "TF": tf_short,
            "motif_id": motif_id_s,
            "score": score_f,
            "p_value": p_val,
            "q_value": q_val,
            "strand": strand_s,
        }
        if "flank" in name:
            tail = name.split("|", 1)[1] if "|" in name else ""
            flank_side = _parse_flank_side(tail)
            hit_dict = {
                "flank_side": flank_side,
                **hit_core,
            }
            tf_map = flank_heaps[sv_id]
            if tf_short not in tf_map:
                tf_map[tf_short] = _PValueHitHeap(k_flank_tf)
            tf_map[tf_short].add(hit_dict)
        else:
            elem_id = _extract_elem_id_from_name(name)
            if elem_id is None:
                continue
            tf_m = elem_heaps[sv_id][elem_id]
            if tf_short not in tf_m:
                tf_m[tf_short] = _PValueHitHeap(k_elem_tf)
            tf_m[tf_short].add(dict(hit_core))


def recombine_sv_fimo(
    merged_bed_path: Path,
    sv_csv_path: Path,
    output_csv_path: Path,
) -> None:
    """
    Recombine FIMO motif hits into processed SV CSV.

    Merged BEDs from ``bedtools intersect -wa -wb`` can be very large (every overlap
    of each SV/flank/element interval with each FIMO hit, across all motifs that pass
    FIMO's ``--thresh``). Rows are read in chunks of ``THRESHOLDS.sv_fimo_recombine_chunk_rows``.
    Retained hits are **capped per TF** using ``THRESHOLDS.sv_fimo_max_flank_hits_per_tf`` and
    ``THRESHOLDS.sv_fimo_max_elem_hits_per_tf`` (best p-value within each TF bucket). Python
    never holds tens of millions of dicts in RAM.

    When the merged BED is large (``THRESHOLDS.sv_fimo_recombine_stream_if_merged_bed_bytes``),
    the SV CSV is read and written in ``THRESHOLDS.sv_fimo_recombine_sv_csv_chunk_rows`` chunks so
    the full table with serialized ``elem_hits`` / ``flank_motif_hits`` is never held at once.

    Args:
        merged_bed_path: Path to merged FIMO-SV intersection BED
        sv_csv_path: Path to processed SV CSV
        output_csv_path: Path for output CSV with motif hits
    """
    merged_bed_path = Path(merged_bed_path)
    sv_csv_path = Path(sv_csv_path)
    output_csv_path = Path(output_csv_path)

    if not merged_bed_path.is_file():
        print(f"[WARN] Merged BED not found: {merged_bed_path}")
        return
    if not sv_csv_path.is_file():
        print(f"[WARN] SV CSV not found: {sv_csv_path}")
        return

    hdr = pd.read_csv(sv_csv_path, nrows=0)
    for col in ("id", "elem_hits", "pos"):
        if col not in hdr.columns:
            print(f"[ERROR] CSV missing required column {col!r}: {sv_csv_path}")
            return

    merged_size = os.path.getsize(merged_bed_path)
    stream_bytes = _fimo_recombine_stream_if_merged_bed_bytes()
    stream_out = stream_bytes > 0 and merged_size >= stream_bytes

    if merged_size == 0:
        print(f"[INFO] Empty BED file: {merged_bed_path}")
        s_empty = pd.read_csv(sv_csv_path)
        _prepare_sv_df_for_motif_attach(s_empty)
        s_empty.drop(columns=["elem_hits_parsed"], inplace=True)
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        s_empty.to_csv(output_csv_path, index=False)
        return

    if stream_out:
        print(
            f"[INFO] {merged_bed_path.name}: merged BED {merged_size / (1024 * 1024):.1f} MiB >= "
            f"{stream_bytes / (1024 * 1024):.1f} MiB; streaming SV CSV in chunks of "
            f"{_fimo_recombine_sv_csv_chunk_rows()} rows "
            "(THRESHOLDS.sv_fimo_recombine_stream_if_merged_bed_bytes / "
            "sv_fimo_recombine_sv_csv_chunk_rows)"
        )
        id_to_indices = _build_sv_id_to_indices_ids_only(sv_csv_path)
        s: Optional[pd.DataFrame] = None
    else:
        s = pd.read_csv(sv_csv_path)
        _prepare_sv_df_for_motif_attach(s)
        id_to_indices = _build_sv_id_to_indices_from_df(s)

    mcols = [
        "sv_chrom", "sv_start", "sv_end", "sv_name",
        "motif_chrom", "motif_start", "motif_end",
        "motif_id", "p_value", "q_value", "score", "strand",
    ]
    str_dtype = {c: "str" for c in ("sv_chrom", "sv_name", "motif_chrom", "motif_id", "strand")}
    chunk_rows = _fimo_recombine_chunk_rows()
    k_flank_tf = _fimo_max_flank_hits_per_tf()
    k_elem_tf = _fimo_max_elem_hits_per_tf()

    flank_heaps: Dict[str, Dict[str, _PValueHitHeap]] = defaultdict(dict)
    elem_heaps: Dict[str, Dict[str, Dict[str, _PValueHitHeap]]] = defaultdict(
        lambda: defaultdict(dict)
    )

    read_kw: Dict[str, Any] = dict(
        sep="\t",
        header=None,
        names=mcols,
        dtype=str_dtype,
        chunksize=chunk_rows,
    )
    try:
        reader = pd.read_csv(merged_bed_path, **read_kw, engine="c")
    except Exception:
        reader = pd.read_csv(merged_bed_path, **read_kw, engine="python")

    n_merged = 0
    n_chunk = 0
    for merged in reader:
        n_chunk += 1
        n_merged += len(merged)
        if not merged.empty:
            _accumulate_merged_bed_chunk(merged, flank_heaps, elem_heaps, k_flank_tf, k_elem_tf)
        del merged
        if n_chunk % 8 == 0:
            gc.collect()

    if n_merged == 0:
        if stream_out:
            s0 = pd.read_csv(sv_csv_path)
            _prepare_sv_df_for_motif_attach(s0)
            s0.drop(columns=["elem_hits_parsed"], inplace=True)
            output_csv_path.parent.mkdir(parents=True, exist_ok=True)
            s0.to_csv(output_csv_path, index=False)
        else:
            assert s is not None
            s.drop(columns=["elem_hits_parsed"], inplace=True)
            output_csv_path.parent.mkdir(parents=True, exist_ok=True)
            s.to_csv(output_csv_path, index=False)
        return

    if n_merged > chunk_rows:
        print(
            f"[INFO] {merged_bed_path.name}: processed {n_merged} bedtools intersection rows "
            f"in chunks of {chunk_rows} (THRESHOLDS.sv_fimo_recombine_chunk_rows)"
        )
    if n_merged > 500_000:
        print(
            f"[INFO] {merged_bed_path.name}: retaining at most {k_flank_tf} flank hits per (SV,TF) and "
            f"{k_elem_tf} per (SV,elem,TF) by best p-value (THRESHOLDS.sv_fimo_max_flank_hits_per_tf / "
            f"sv_fimo_max_elem_hits_per_tf)"
        )

    def _sort_hits(lst: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        lst.sort(key=lambda d: float(d.get("p_value") or 1.0))
        return lst

    flank_hits_by_sv: Dict[str, List[Dict[str, Any]]] = {}
    for sv_id, tf_map in flank_heaps.items():
        merged_hits: List[Dict[str, Any]] = []
        for _tf, h in tf_map.items():
            merged_hits.extend(h.to_sorted_hits())
        flank_hits_by_sv[sv_id] = _sort_hits(merged_hits)

    elem_motif_hits: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for sv_id, em in elem_heaps.items():
        elem_motif_hits[sv_id] = {}
        for elem_id, tf_m in em.items():
            merged_e: List[Dict[str, Any]] = []
            for _tf, h in tf_m.items():
                merged_e.extend(h.to_sorted_hits())
            elem_motif_hits[sv_id][elem_id] = _sort_hits(merged_e)
    del flank_heaps, elem_heaps
    gc.collect()

    output_csv_path.parent.mkdir(parents=True, exist_ok=True)

    if stream_out:
        sv_chunk = _fimo_recombine_sv_csv_chunk_rows()
        row_base = 0
        first_write = True
        n_sv_chunk = 0
        for chunk in pd.read_csv(sv_csv_path, chunksize=sv_chunk):
            n_sv_chunk += 1
            chunk = chunk.copy()
            chunk.index = pd.RangeIndex(row_base, row_base + len(chunk))
            _prepare_sv_df_for_motif_attach(chunk)
            _attach_flank_and_elem_motif_hits(chunk, id_to_indices, flank_hits_by_sv, elem_motif_hits)
            _finalize_sv_motif_columns_for_csv(chunk)
            chunk.to_csv(
                output_csv_path,
                mode="w" if first_write else "a",
                header=first_write,
                index=False,
            )
            first_write = False
            row_base += len(chunk)
            del chunk
            if n_sv_chunk % 8 == 0:
                gc.collect()
    else:
        assert s is not None
        _attach_flank_and_elem_motif_hits(s, id_to_indices, flank_hits_by_sv, elem_motif_hits)
        _finalize_sv_motif_columns_for_csv(s)
        s.to_csv(output_csv_path, index=False)

    print(f"[INFO] Wrote: {output_csv_path}")


_FLANKS_SV_BED_SUFFIX = "_flanks_and_overlaps"


def _sv_csv_stem_from_flanks_bed_stem(bed_stem: str) -> str:
    """Map ``{csv_stem}_flanks_and_overlaps`` BED/FASTA/FIMO base to processed SV CSV stem."""
    if bed_stem.endswith(_FLANKS_SV_BED_SUFFIX):
        return bed_stem[: -len(_FLANKS_SV_BED_SUFFIX)]
    return bed_stem


def _sv_fimo_recombine_skip_csv_stems() -> frozenset:
    raw = THRESHOLDS.sv_fimo_recombine_skip_csv_basenames
    if not raw:
        return frozenset()
    return frozenset(str(x).strip() for x in raw if str(x).strip())


def recombine_all_sv_fimo(
    merged_bed_dir: Path,
    sv_csv_dir: Path,
    output_csv_dir: Path,
) -> List[Path]:
    """
    Recombine FIMO hits for all SVs in a directory.
    
    Args:
        merged_bed_dir: Directory with merged FIMO-SV BED files
        sv_csv_dir: Directory with processed SV CSVs
        output_csv_dir: Directory for output CSVs
    
    Returns:
        List of created output file paths
    """
    merged_bed_dir = Path(merged_bed_dir)
    sv_csv_dir = Path(sv_csv_dir)
    output_csv_dir = Path(output_csv_dir)
    output_csv_dir.mkdir(parents=True, exist_ok=True)
    
    created_files: List[Path] = []
    skip_stems = _sv_fimo_recombine_skip_csv_stems()

    for csv_stem in skip_stems:
        csv_name = f"{csv_stem}.csv"
        sv_csv_path = sv_csv_dir / csv_name
        out_csv_path = output_csv_dir / csv_name
        if not sv_csv_path.is_file():
            print(
                f"[WARN] Skip-list stem {csv_stem!r}: missing processed SV CSV {sv_csv_path}; "
                "cannot copy to final output.",
                flush=True,
            )
            continue
        shutil.copy2(sv_csv_path, out_csv_path)
        created_files.append(out_csv_path)
        print(
            f"[INFO] Skipped FIMO recombine for {csv_stem!r}: copied processed SV to {out_csv_path.name} "
            "(THRESHOLDS.sv_fimo_recombine_skip_csv_basenames). Existing merged BED is ignored.",
            flush=True,
        )

    merged_fnames = sorted(
        f for f in os.listdir(merged_bed_dir) if f.endswith("_fimo_merged.bed")
    )
    work: List[Tuple[str, Path, Path, Path]] = []
    for fname in merged_fnames:
        merged_bed_path = merged_bed_dir / fname
        base = fname[: -len("_fimo_merged.bed")]
        if base.endswith(_FLANKS_SV_BED_SUFFIX):
            base_csv = base[: -len(_FLANKS_SV_BED_SUFFIX)]
        else:
            base_csv = base
        if base_csv in skip_stems:
            continue
        csv_name = base_csv + ".csv"
        work.append((fname, merged_bed_path, sv_csv_dir / csv_name, output_csv_dir / csv_name))

    n_work = len(work)
    for i, (fname, merged_bed_path, sv_csv_path, out_csv_path) in enumerate(work, 1):
        sz_mb = merged_bed_path.stat().st_size / (1024 * 1024)
        print(
            f"[STEP 5] ({i}/{n_work}) {fname} (~{sz_mb:.1f} MiB merged) → {out_csv_path.name} ...",
            flush=True,
        )
        try:
            recombine_sv_fimo(merged_bed_path, sv_csv_path, out_csv_path)
        except Exception as exc:
            print(f"[ERROR] recombine_sv_fimo failed for {fname}: {exc}", flush=True)
            traceback.print_exc()
            raise
        gc.collect()

        if out_csv_path.exists():
            created_files.append(out_csv_path)

    return created_files


# =============================================================================
# MOTIF SUMMARY
# =============================================================================

def summarize_motif_hits(sv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize motif hits across all SVs.
    
    Returns DataFrame with TF hit counts.
    """
    tf_counts = defaultdict(int)
    
    # Count flank motif hits
    if "flank_motif_hits" in sv_df.columns:
        for hits in sv_df["flank_motif_hits"]:
            if isinstance(hits, str):
                try:
                    hits = ast.literal_eval(hits)
                except Exception:
                    continue
            if isinstance(hits, list):
                for h in hits:
                    if isinstance(h, dict):
                        tf = h.get("TF", "unknown")
                        tf_counts[tf] += 1
    
    # Count elem motif hits
    if "elem_hits" in sv_df.columns:
        for elem_hits in sv_df["elem_hits"]:
            if isinstance(elem_hits, str):
                try:
                    elem_hits = ast.literal_eval(elem_hits)
                except Exception:
                    continue
            if isinstance(elem_hits, list):
                for elem in elem_hits:
                    if isinstance(elem, dict):
                        motif_hits = elem.get("motif_hits", [])
                        if isinstance(motif_hits, list):
                            for h in motif_hits:
                                if isinstance(h, dict):
                                    tf = h.get("TF", "unknown")
                                    tf_counts[tf] += 1
    
    if not tf_counts:
        return pd.DataFrame(columns=["TF", "count"])
    
    return pd.DataFrame([
        {"TF": tf, "count": count}
        for tf, count in sorted(tf_counts.items(), key=lambda x: -x[1])
    ])


# =============================================================================
# ALL-cCRE FIMO LOADER
# =============================================================================

def load_all_ccre_fimo(
    fimo_tsv_path: Path,
    ccre_table: pd.DataFrame,
) -> Dict[str, List[Dict]]:
    """
    Load genome-wide FIMO results and assign motif hits to cCREs.

    In this FIMO run the ``sequence_name`` column is the **chromosome**
    (e.g. ``3``, ``22``, ``X``) and ``start``/``stop`` are **genomic
    coordinates** of the motif match.  Each motif hit is assigned to
    every cCRE it overlaps (a hit may fall in zero or more cCREs).

    The cCRE table must have columns: ``cCRE_id``, ``chrom``, ``start``,
    ``end`` (and optionally ``ENCODE_id``, ``type``).

    Args:
        fimo_tsv_path: Path to ``all_ccre_fimo.tsv``
        ccre_table:    Full cCRE coordinate table

    Returns:
        Dict mapping ``cCRE_id`` to a list of motif-hit dicts::

            {
                "motif_id":  str,
                "TF":        str,   # short name (before first '.')
                "start":     int,   # genomic start of motif hit
                "end":       int,   # genomic end   of motif hit
                "strand":    str,
                "score":     float,
                "p_value":   float,
                "q_value":   float,
            }
    """
    fimo_tsv_path = Path(fimo_tsv_path)

    print(f"[INFO] Reading FIMO TSV: {fimo_tsv_path}")
    df = pd.read_csv(fimo_tsv_path, sep="\t", comment="#")

    # Drop summary / malformed rows
    df["start_num"] = pd.to_numeric(df["start"], errors="coerce")
    df["stop_num"] = pd.to_numeric(df["stop"], errors="coerce")
    df = df.dropna(subset=["start_num", "stop_num"]).copy()
    df["start_num"] = df["start_num"].astype(int)
    df["stop_num"] = df["stop_num"].astype(int)

    # Normalise chromosome: sequence_name is the chrom
    df["chrom"] = (
        df["sequence_name"]
        .astype(str)
        .str.strip()
        .str.replace("chr", "", regex=False)
    )

    print(f"[INFO] {len(df)} motif hits across {df['chrom'].nunique()} chromosomes")

    # ── Prepare cCRE table ──
    ccre = ccre_table[["cCRE_id", "chrom", "start", "end"]].copy()
    ccre["chrom"] = ccre["chrom"].astype(str).str.replace("chr", "", regex=False)
    ccre["start"] = pd.to_numeric(ccre["start"], errors="coerce").astype(int)
    ccre["end"] = pd.to_numeric(ccre["end"], errors="coerce").astype(int)

    # ── Intersect: assign each motif hit to overlapping cCREs ──
    # Process per-chromosome to keep memory manageable
    lookup: Dict[str, List[Dict]] = defaultdict(list)
    n_assigned = 0

    for chrom_val, fimo_chr in df.groupby("chrom"):
        ccre_chr = ccre[ccre["chrom"] == chrom_val]
        if ccre_chr.empty:
            continue

        # Vectorised overlap: for each FIMO hit, find overlapping cCREs
        # motif overlaps cCRE if motif.end >= ccre.start AND motif.start <= ccre.end
        fimo_starts = fimo_chr["start_num"].values
        fimo_ends = fimo_chr["stop_num"].values
        ccre_starts = ccre_chr["start"].values
        ccre_ends = ccre_chr["end"].values
        ccre_ids = ccre_chr["cCRE_id"].values

        for f_idx in range(len(fimo_chr)):
            f_row = fimo_chr.iloc[f_idx]
            f_start = fimo_starts[f_idx]
            f_end = fimo_ends[f_idx]

            # Boolean mask over cCREs on this chrom
            overlaps = (f_end >= ccre_starts) & (f_start <= ccre_ends)
            if not overlaps.any():
                continue

            motif_id = str(f_row["motif_id"])
            tf_short = motif_id.split(".", 1)[0] if "." in motif_id else motif_id

            try:
                p_val = float(f_row["p-value"])
            except Exception:
                p_val = np.nan
            try:
                q_val = float(f_row.get("q-value", np.nan))
            except Exception:
                q_val = np.nan
            try:
                score = float(f_row.get("score", np.nan))
            except Exception:
                score = np.nan

            hit = {
                "motif_id": motif_id,
                "TF": tf_short,
                "start": f_start,
                "end": f_end,
                "strand": str(f_row["strand"]),
                "score": score,
                "p_value": p_val,
                "q_value": q_val,
            }

            for cid in ccre_ids[overlaps]:
                lookup[cid].append(hit)
                n_assigned += 1

    print(f"[INFO] Assigned {n_assigned} motif-cCRE associations "
          f"across {len(lookup)} cCREs")

    return dict(lookup)


def _build_ccre_chrom_index(
    ccre_table: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """
    Pre-index cCRE table by chromosome for fast spatial queries.

    Accepts columns named either ``cCRE_id`` or ``elem_id`` (the former
    is the canonical name in the full cCRE table; the latter may appear
    in the SV pipeline's filtered focus set).  Whichever is present is
    aliased to ``_ccre_key`` internally, but the original column name is
    preserved for downstream use.

    Chromosome names are normalised (strip 'chr' prefix) and the
    resulting per-chromosome DataFrames are sorted by start position.
    """
    ccre = ccre_table.copy()

    # Normalise the id column name — keep original, add alias
    if "cCRE_id" in ccre.columns:
        ccre["_ccre_key"] = ccre["cCRE_id"]
    elif "elem_id" in ccre.columns:
        ccre["_ccre_key"] = ccre["elem_id"]
    else:
        raise ValueError(
            "cCRE table must have a 'cCRE_id' or 'elem_id' column. "
            f"Found columns: {list(ccre.columns)}"
        )

    ccre["chrom"] = ccre["chrom"].astype(str).str.replace("chr", "", regex=False)
    ccre["start"] = pd.to_numeric(ccre["start"], errors="coerce")
    ccre["end"] = pd.to_numeric(ccre["end"], errors="coerce")
    ccre = ccre.dropna(subset=["start", "end"])
    ccre["start"] = ccre["start"].astype(int)
    ccre["end"] = ccre["end"].astype(int)
    ccre = ccre.sort_values(["chrom", "start"])

    return {chrom: grp for chrom, grp in ccre.groupby("chrom")}


# =============================================================================
# BND ORIENTATION HELPERS
# =============================================================================

def _parse_bnd_orientations(row: pd.Series) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract (orientation_self, orientation_remote) from a BND row.

    Tries, in order:
    1. Dedicated columns ``orientation_self`` / ``orientation_remote``
       (present if previously parsed).
    2. The ``orig_alt`` column containing the Manta BreakEnd(...) repr.
    3. The ``alt`` column in VCF BND notation  (e.g. ``N]chr5:123]``).

    Returns (None, None) if orientation cannot be determined.
    """
    # --- attempt 1: explicit columns ---
    o_self = row.get("orientation_self")
    o_remote = row.get("orientation_remote")
    if (
        isinstance(o_self, str) and o_self in "+-"
        and isinstance(o_remote, str) and o_remote in "+-"
    ):
        return o_self, o_remote

    # --- attempt 2: BreakEnd repr in orig_alt ---
    orig_alt = row.get("orig_alt", None)
    if isinstance(orig_alt, str) and "BreakEnd(" in orig_alt:
        m = re.search(
            r"BreakEnd\('[^']+',\s*[0-9]+,\s*'([+-])',\s*'([+-])'",
            orig_alt,
        )
        if m:
            return m.group(1), m.group(2)

    # --- attempt 3: VCF BND alt field ---
    alt = row.get("alt", None)
    if isinstance(alt, str):
        # VCF BND encoding rules  (spec §5.4):
        #   t[p[  =>  self +, remote -       (N extends right, mate extends right)
        #   t]p]  =>  self +, remote +       (N extends right, mate extends left)
        #   ]p]t  =>  self -, remote +       (N extends left,  mate extends left)
        #   [p[t  =>  self -, remote -       (N extends left,  mate extends right)
        if re.match(r"^[A-Za-z]\[", alt):
            return "+", "-"
        elif re.match(r"^[A-Za-z]\]", alt):
            return "+", "+"
        elif alt.startswith("]"):
            return "-", "+"
        elif alt.startswith("["):
            return "-", "-"

    return None, None


def _get_neojunction_window(
    mate_chrom: str,
    mate_pos: int,
    orientation_remote: str,
    window: int,
) -> Tuple[str, int, int]:
    """
    Return the genomic interval on the mate's side that is physically
    juxtaposed to the local breakpoint after the rearrangement.

    The orientation of the mate breakpoint tells us which direction the
    joined sequence extends from the mate position:

    - ``orientation_remote == "+"``  →  the *upstream* (left) side of
      the mate is joined.  The sequence that becomes adjacent is
      ``[mate_pos - window, mate_pos]``.

    - ``orientation_remote == "-"``  →  the *downstream* (right) side
      of the mate is joined.  The gained interval is
      ``[mate_pos, mate_pos + window]``.

    Args:
        mate_chrom:         Chromosome of the mate breakpoint.
        mate_pos:           Position of the mate breakpoint.
        orientation_remote: ``"+"`` or ``"-"``.
        window:             How far (bp) from the mate to search.

    Returns:
        ``(chrom, start, end)`` of the gained genomic interval.
    """
    mate_chrom = str(mate_chrom).replace("chr", "")
    if orientation_remote == "+":
        return (mate_chrom, max(1, mate_pos - window), mate_pos)
    else:
        return (mate_chrom, mate_pos, mate_pos + window)


# =============================================================================
# NEOJUNCTION MOTIF ANNOTATION
# =============================================================================

def annotate_bnd_neojunction_motifs(
    sv_df: pd.DataFrame,
    ccre_fimo_lookup: Dict[str, List[Dict]],
    ccre_table: pd.DataFrame,
    window: int = 500_000,
) -> pd.DataFrame:
    """
    For each BND, identify cCREs in the **neojunction window** — the
    genomic region near the mate breakpoint that becomes physically
    adjacent to the local breakpoint after the rearrangement — and
    attach their pre-scanned motif content from *ccre_fimo_lookup*.

    This captures the regulatory potential *gained* by the rearrangement:
    a cCRE that was previously far away (or on another chromosome) may
    now sit within enhancer-range of an APM gene.

    The result is stored in a new column ``neojunction_motif_hits``,
    which is a list of dicts per BND row::

        [
            {
                "gained_at":           "bp1" | "bp2",
                "cCRE_id":             str,
                "elem_type":           str | None,
                "elem_chrom":          str,
                "elem_start":          int,
                "elem_end":            int,
                "distance_from_mate":  int,   # bp from mate breakpoint
                                              #   to cCRE midpoint
                "mate_chrom":          str,
                "mate_pos":            int,
                "local_chrom":         str,
                "local_pos":           int,
                "orientation_self":    str,
                "orientation_remote":  str,
                "motif_hits":          [  ... motif dicts from FIMO ... ],
                "n_motif_hits":        int,
                "tf_summary":          str,   # comma-separated unique TFs
            },
            ...
        ]

    Non-BND rows get an empty list.

    Args:
        sv_df:             Processed SV DataFrame (must contain BND rows
                           with ``bnd_remote_chrom``, ``bnd_remote_pos``,
                           and either orientation columns or parseable alt).
        ccre_fimo_lookup:  Output of :func:`load_all_ccre_fimo` — keyed
                           by ``cCRE_id``.
        ccre_table:        Full cCRE coordinate table (``cCRE_id`` or
                           ``elem_id``, ``chrom``, ``start``, ``end``;
                           optionally ``type``).
        window:            How far (bp) from the mate breakpoint to search
                           for gained cCREs.  Default 500 kb.

    Returns:
        Copy of *sv_df* with ``neojunction_motif_hits`` column added.
    """
    sv_df = sv_df.copy()
    sv_df["neojunction_motif_hits"] = [[] for _ in range(len(sv_df))]

    bnds = sv_df[sv_df["SVTYPE"] == "BND"].copy()
    if bnds.empty:
        return sv_df

    # Pre-index cCREs by chromosome
    ccre_by_chrom = _build_ccre_chrom_index(ccre_table)

    n_annotated = 0
    n_gained_ccres = 0

    for idx, row in bnds.iterrows():
        bp1_chrom = str(row["chrom"]).replace("chr", "")
        bp1_pos = int(row["pos"])

        bp2_chrom = row.get("bnd_remote_chrom")
        bp2_pos = row.get("bnd_remote_pos")

        if pd.isna(bp2_chrom) or bp2_chrom is None or pd.isna(bp2_pos) or bp2_pos is None:
            continue

        bp2_chrom = str(bp2_chrom).replace("chr", "")
        bp2_pos = int(bp2_pos)

        orient_self, orient_remote = _parse_bnd_orientations(row)
        if orient_self is None or orient_remote is None:
            # Cannot determine directionality — skip
            continue

        all_hits = []

        # --- cCREs gained at bp1 (from bp2's side) ---
        win_chrom, win_start, win_end = _get_neojunction_window(
            bp2_chrom, bp2_pos, orient_remote, window,
        )
        gained_at_bp1 = _query_ccres_in_window(
            ccre_by_chrom, win_chrom, win_start, win_end,
        )
        for _, ccre_row in gained_at_bp1.iterrows():
            ccre_id = str(ccre_row["_ccre_key"])
            motifs = ccre_fimo_lookup.get(ccre_id, [])

            mid = (int(ccre_row["start"]) + int(ccre_row["end"])) // 2
            dist = abs(mid - bp2_pos)

            tfs = sorted({h["TF"] for h in motifs})

            all_hits.append({
                "gained_at": "bp1",
                "cCRE_id": ccre_id,
                "elem_type": ccre_row.get("type", None),
                "elem_chrom": win_chrom,
                "elem_start": int(ccre_row["start"]),
                "elem_end": int(ccre_row["end"]),
                "distance_from_mate": dist,
                "mate_chrom": bp2_chrom,
                "mate_pos": bp2_pos,
                "local_chrom": bp1_chrom,
                "local_pos": bp1_pos,
                "orientation_self": orient_self,
                "orientation_remote": orient_remote,
                "motif_hits": motifs,
                "n_motif_hits": len(motifs),
                "tf_summary": ",".join(tfs),
            })

        # --- cCREs gained at bp2 (from bp1's side) ---
        #
        # Symmetric: now the "mate" from bp2's perspective is bp1,
        # and the orientation that matters is orientation_self (which
        # describes bp1's joining direction).
        win_chrom2, win_start2, win_end2 = _get_neojunction_window(
            bp1_chrom, bp1_pos, orient_self, window,
        )
        gained_at_bp2 = _query_ccres_in_window(
            ccre_by_chrom, win_chrom2, win_start2, win_end2,
        )
        for _, ccre_row in gained_at_bp2.iterrows():
            ccre_id = str(ccre_row["_ccre_key"])
            motifs = ccre_fimo_lookup.get(ccre_id, [])

            mid = (int(ccre_row["start"]) + int(ccre_row["end"])) // 2
            dist = abs(mid - bp1_pos)

            tfs = sorted({h["TF"] for h in motifs})

            all_hits.append({
                "gained_at": "bp2",
                "cCRE_id": ccre_id,
                "elem_type": ccre_row.get("type", None),
                "elem_chrom": win_chrom2,
                "elem_start": int(ccre_row["start"]),
                "elem_end": int(ccre_row["end"]),
                "distance_from_mate": dist,
                "mate_chrom": bp1_chrom,
                "mate_pos": bp1_pos,
                "local_chrom": bp2_chrom,
                "local_pos": bp2_pos,
                "orientation_self": orient_remote,
                "orientation_remote": orient_self,
                "motif_hits": motifs,
                "n_motif_hits": len(motifs),
                "tf_summary": ",".join(tfs),
            })

        if all_hits:
            n_annotated += 1
            n_gained_ccres += len(all_hits)

        sv_df.at[idx, "neojunction_motif_hits"] = all_hits

    print(f"[INFO] Neojunction motif annotation: "
          f"{n_annotated}/{len(bnds)} BNDs annotated, "
          f"{n_gained_ccres} gained cCRE associations total")

    return sv_df


def _query_ccres_in_window(
    ccre_by_chrom: Dict[str, pd.DataFrame],
    chrom: str,
    start: int,
    end: int,
) -> pd.DataFrame:
    """
    Return cCREs that overlap the interval ``[start, end]`` on *chrom*.

    Uses the pre-indexed dict from :func:`_build_ccre_chrom_index`.
    A cCRE overlaps the window if ``cCRE.end >= start AND cCRE.start <= end``.
    """
    chrom = str(chrom).replace("chr", "")
    chr_df = ccre_by_chrom.get(chrom)
    if chr_df is None or chr_df.empty:
        return pd.DataFrame()

    mask = (chr_df["end"] >= start) & (chr_df["start"] <= end)
    return chr_df[mask]


# =============================================================================
# NEOJUNCTION — BATCH / DIRECTORY HELPERS
# =============================================================================

def annotate_neojunction_motifs_from_directory(
    sv_csv_dir: Path,
    output_dir: Path,
    all_ccre_fimo_path: Path,
    ccre_table_path: Path,
    window: int = 500_000,
) -> List[Path]:
    """
    Apply neojunction motif annotation to all SV CSVs in a directory.

    This is intended to be called as a post-processing step after the
    main SV pipeline (after ``recombine_all_sv_fimo`` or on the final
    CSVs).

    Args:
        sv_csv_dir:         Directory with processed SV CSVs.
        output_dir:         Directory for enriched output CSVs.
        all_ccre_fimo_path: Path to ``all_ccre_fimo.tsv``.
        ccre_table_path:    Path to the full cCRE coordinate CSV
                            (must have ``cCRE_id``, ``chrom``, ``start``,
                            ``end`` columns).
        window:             Neojunction search window in bp.

    Returns:
        List of created output file paths.
    """
    sv_csv_dir = Path(sv_csv_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load cCRE table first (needed by both the FIMO loader and the annotator)
    print("[INFO] Loading cCRE coordinate table...")
    ccre_table = pd.read_csv(ccre_table_path)

    # Load FIMO and intersect with cCRE coordinates — once, reused across samples
    print("[INFO] Loading all-cCRE FIMO results...")
    ccre_fimo_lookup = load_all_ccre_fimo(all_ccre_fimo_path, ccre_table)

    created = []

    for fname in sorted(os.listdir(sv_csv_dir)):
        if not fname.endswith(".csv"):
            continue

        csv_path = sv_csv_dir / fname
        sv_df = pd.read_csv(csv_path)

        n_bnd = (sv_df["SVTYPE"] == "BND").sum() if "SVTYPE" in sv_df.columns else 0
        if n_bnd == 0:
            # Nothing to annotate — copy through unchanged
            out_path = output_dir / fname
            sv_df.to_csv(out_path, index=False)
            created.append(out_path)
            continue

        print(f"\n  Processing {fname} ({n_bnd} BNDs)...")
        sv_df = annotate_bnd_neojunction_motifs(
            sv_df, ccre_fimo_lookup, ccre_table, window=window,
        )

        # Serialise the new column for CSV storage
        sv_df["neojunction_motif_hits"] = sv_df["neojunction_motif_hits"].apply(
            lambda x: "[]" if (x is None or x == []) else repr(x)
        )

        out_path = output_dir / fname
        sv_df.to_csv(out_path, index=False)
        created.append(out_path)
        print(f"  Wrote: {out_path}")

    return created


def summarize_neojunction_motifs(sv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise neojunction motif hits across all BNDs in a DataFrame.

    Returns a DataFrame with columns:
        TF, n_ccres (how many gained cCREs carry this TF),
        n_total_hits (total motif instances), n_bnds (how many BNDs
        contribute).
    """
    tf_stats: Dict[str, Dict[str, int]] = defaultdict(
        lambda: {"n_ccres": 0, "n_total_hits": 0, "bnd_ids": set()}
    )

    if "neojunction_motif_hits" not in sv_df.columns:
        return pd.DataFrame(columns=["TF", "n_ccres", "n_total_hits", "n_bnds"])

    for row_idx, cell in sv_df["neojunction_motif_hits"].items():
        if isinstance(cell, str):
            try:
                cell = ast.literal_eval(cell)
            except Exception:
                continue
        if not isinstance(cell, list):
            continue

        sv_id = str(sv_df.at[row_idx, "id"]) if "id" in sv_df.columns else str(row_idx)

        for gained in cell:
            if not isinstance(gained, dict):
                continue
            motifs = gained.get("motif_hits", [])
            if not isinstance(motifs, list):
                continue

            tfs_in_ccre = set()
            for h in motifs:
                if isinstance(h, dict):
                    tf = h.get("TF", "unknown")
                    tfs_in_ccre.add(tf)
                    tf_stats[tf]["n_total_hits"] += 1
                    tf_stats[tf]["bnd_ids"].add(sv_id)

            for tf in tfs_in_ccre:
                tf_stats[tf]["n_ccres"] += 1

    if not tf_stats:
        return pd.DataFrame(columns=["TF", "n_ccres", "n_total_hits", "n_bnds"])

    rows = []
    for tf, stats in sorted(tf_stats.items(), key=lambda x: -x[1]["n_total_hits"]):
        rows.append({
            "TF": tf,
            "n_ccres": stats["n_ccres"],
            "n_total_hits": stats["n_total_hits"],
            "n_bnds": len(stats["bnd_ids"]),
        })

    return pd.DataFrame(rows)
