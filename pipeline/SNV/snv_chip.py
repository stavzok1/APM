"""
ChIP-seq overlap for SNVs against ``PATHS.chip_unified`` (ENCODE + ChIP-Atlas).

**Strict overlap:** a peak ``[start, end)`` (0-based BED) overlaps the variant iff it
covers the VCF ``POS`` base (1-based), i.e. ``start <= POS - 1`` and ``end >= POS``.
No flanking window.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..biosample_names import normalize_cell_line_label
from ..utils import harmonize_chrom_column


_CHIP_READ_COLUMNS = (
    "chrom",
    "start",
    "end",
    "tf",
    "cell_type",
    "source",
    "score_norm",
    "sample_id",
    "cell_subtype",
)


def _snv_chip_stratum_from_peak_row(row: pd.Series) -> str:
    """Coarse biosample bucket for aggregation (subtype when available)."""
    st = str(row.get("cell_subtype", "") or "").strip()
    if st:
        return f"subtype:{st}"
    raw = str(row.get("cell_type", "") or "").lower()
    sid = str(row.get("sample_id", "") or "")
    if "tissue" in raw or "TISSUE" in sid.upper():
        return "context:tissue"
    if any(k in raw for k in ("hmec", "normal-like", "normal ", "blood", "lymphocyte")):
        return "context:normal_like"
    if "mcf10a" in raw or raw.strip() in ("mcf 10a", "mcf10a"):
        return "context:mcf10a_like"
    return "context:cell_line"


def _peak_row_to_hit(row: pd.Series, pos_1based: int, chrom: str) -> Dict[str, Any]:
    ps = int(row["start"])
    pe = int(row["end"])
    pos0 = pos_1based - 1
    sn = row.get("sample_id")
    sc = row.get("score_norm")
    csub = row.get("cell_subtype", "")
    return {
        "tf": str(row.get("tf", "") or ""),
        "cell_type": normalize_cell_line_label(str(row.get("cell_type", "") or "")),
        "cell_subtype": "" if pd.isna(csub) else str(csub).strip(),
        "source": str(row.get("source", "") or ""),
        "sample_id": None if pd.isna(sn) else str(sn),
        "score_norm": None if pd.isna(sc) else float(sc),
        "chrom": chrom,
        "peak_start": ps,
        "peak_end": pe,
        "variant_pos": int(pos_1based),
        "overlap_bp": 1,
        "overlap_start_0based": pos0,
        "overlap_end_0based_exclusive": pos0 + 1,
        "stratum": _snv_chip_stratum_from_peak_row(row),
    }


def _aggregate_hits(hits: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Per (tf, source, stratum): peak counts, score summaries, distinct ChIP ``sample_id``s.
    """
    # (tf, source, stratum) -> lists
    scores: Dict[Tuple[str, str, str], List[float]] = defaultdict(list)
    samples: Dict[Tuple[str, str, str], set] = defaultdict(set)
    counts: Dict[Tuple[str, str, str], int] = defaultdict(int)

    for h in hits:
        tf = str(h.get("tf") or "")
        src = str(h.get("source") or "")
        st = str(h.get("stratum") or "")
        key = (tf, src, st)
        counts[key] += 1
        sc = h.get("score_norm")
        if sc is not None and isinstance(sc, (int, float)) and not (isinstance(sc, float) and np.isnan(sc)):
            scores[key].append(float(sc))
        sid = h.get("sample_id")
        if sid:
            samples[key].add(str(sid))

    rows: List[Dict[str, Any]] = []
    for key, n in sorted(counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1], x[0][2])):
        tf, src, st = key
        scs = scores[key]
        rows.append(
            {
                "tf": tf,
                "source": src,
                "stratum": st,
                "n_peaks": int(n),
                "mean_score_norm": float(sum(scs) / len(scs)) if scs else None,
                "max_score_norm": float(max(scs)) if scs else None,
                "n_chip_samples": len(samples[key]),
                "chip_sample_ids": sorted(samples[key]),
            }
        )
    return {"by_tf_source_stratum": rows}


def annotate_snvs_with_chip(
    df: pd.DataFrame,
    chip_unified: Optional[Union[str, Path]],
    *,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
) -> pd.DataFrame:
    """
    Add ``snv_chip_hits`` (list of overlapping peaks, strict POS overlap) and
    ``snv_chip_aggregate`` (summaries keyed by tf × source × stratum).

    If ``chip_unified`` is missing or the table is empty, both columns are empty
    structures per row.
    """
    out = df.copy()
    n = len(out)
    if n == 0:
        return out

    if chip_unified is None or not Path(chip_unified).is_file():
        out["snv_chip_hits"] = [[] for _ in range(n)]
        out["snv_chip_aggregate"] = [{"by_tf_source_stratum": []} for _ in range(n)]
        return out

    p = Path(chip_unified)
    try:
        import pyarrow.parquet as pq

        names = set(pq.read_schema(p).names)
        usecols = [c for c in _CHIP_READ_COLUMNS if c in names]
        chip = pd.read_parquet(p, columns=usecols)
    except Exception:
        chip = pd.read_parquet(p)
        chip = chip.loc[:, [c for c in _CHIP_READ_COLUMNS if c in chip.columns]].copy()
    if "cell_subtype" not in chip.columns:
        chip["cell_subtype"] = ""
    chip, _ = harmonize_chrom_column(chip, "chrom")

    tmp = out[[chrom_col]].rename(columns={chrom_col: "chrom"})
    tmp, _ = harmonize_chrom_column(tmp, "chrom")
    out["_chip_chrom"] = tmp["chrom"]

    by_chrom: Dict[str, pd.DataFrame] = {c: g for c, g in chip.groupby("chrom", sort=False)}

    hits_col: List[List[Dict[str, Any]]] = []
    agg_col: List[Dict[str, Any]] = []

    for idx, row in out.iterrows():
        chrom = str(row["_chip_chrom"])
        pos = int(row[pos_col])
        peaks = by_chrom.get(chrom)
        if peaks is None or peaks.empty:
            hits_col.append([])
            agg_col.append({"by_tf_source_stratum": []})
            continue
        m = peaks[(peaks["start"] <= pos - 1) & (peaks["end"] >= pos)]
        if m.empty:
            hits_col.append([])
            agg_col.append({"by_tf_source_stratum": []})
            continue
        hits = [_peak_row_to_hit(r, pos, chrom) for _, r in m.iterrows()]
        hits_col.append(hits)
        agg_col.append(_aggregate_hits(hits))

    out["snv_chip_hits"] = hits_col
    out["snv_chip_aggregate"] = agg_col
    out.drop(columns=["_chip_chrom"], inplace=True)
    return out
