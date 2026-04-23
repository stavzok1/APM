"""
HLA-region LOH from ASCAT3 allelic segment files (GDC export).

Gene-level ASCAT tables sometimes carry min_copy_number == max_copy_number for HLA
loci even when allelic imbalance exists; overlapping ``*.ascat3.allelic_specific.seg.txt``
segments expose ``Major_Copy_Number`` / ``Minor_Copy_Number`` and recover
copy-neutral LOH (total CN 2, minor 0).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from pipeline.CNV.loader import extract_sample_id_from_annotations


def build_tumor_vial_to_ascat_seg_path(*, cnv_seg_dir: Path, ann_path: Path) -> Dict[str, Path]:
    """Map TCGA tumor vial id (e.g. ``TCGA-3C-AALI-01A``) → allelic ASCAT segment file path."""
    ann = pd.read_csv(ann_path, sep="\t", low_memory=False)
    if "File Name" not in ann.columns:
        return {}
    out: Dict[str, Path] = {}
    for fn in ann["File Name"].astype(str).tolist():
        if "ascat3.allelic_specific.seg" not in fn.lower():
            continue
        p = cnv_seg_dir / Path(fn).name
        if not p.is_file():
            continue
        try:
            vial = extract_sample_id_from_annotations(ann, fn)
        except Exception:
            continue
        vial = str(vial).strip()
        if vial:
            out[vial] = p
    return out


def load_ascat_allelic_seg(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cmap = {c.lower(): c for c in df.columns}
    need = {"chromosome", "start", "end", "major_copy_number", "minor_copy_number", "copy_number"}
    missing = need - set(cmap.keys())
    if missing:
        raise ValueError(f"{path}: missing columns {missing}; have {list(df.columns)}")
    out = pd.DataFrame(
        {
            "chrom": df[cmap["chromosome"]].astype(str).str.strip(),
            "start": pd.to_numeric(df[cmap["start"]], errors="coerce").astype("Int64"),
            "end": pd.to_numeric(df[cmap["end"]], errors="coerce").astype("Int64"),
            "cn_total": pd.to_numeric(df[cmap["copy_number"]], errors="coerce"),
            "cn_major": pd.to_numeric(df[cmap["major_copy_number"]], errors="coerce"),
            "cn_minor": pd.to_numeric(df[cmap["minor_copy_number"]], errors="coerce"),
        }
    )
    out["chrom"] = out["chrom"].map(lambda x: x if str(x).lower().startswith("chr") else f"chr{x}")
    return out


def _overlap_len(a0: int, a1: int, b0: int, b1: int) -> int:
    s = max(int(a0), int(b0))
    e = min(int(a1), int(b1))
    return max(0, e - s + 1)


def best_segment_for_interval(
    segs: pd.DataFrame,
    chrom: str,
    g0: int,
    g1: int,
) -> Tuple[int, Optional[int], Optional[int], Optional[int], Optional[float]]:
    """
    Pick the overlapping segment with largest overlap_bp on ``chrom``.

    Returns:
        (overlap_bp, cn_minor, cn_major, cn_total, overlap_frac_of_gene)
    """
    ch = str(chrom).strip()
    if not ch.lower().startswith("chr"):
        ch = f"chr{ch}"
    m = segs["chrom"].astype(str) == ch
    sub = segs.loc[m].copy()
    if sub.empty or g1 < g0:
        return 0, None, None, None, None
    gene_len = int(g1) - int(g0) + 1
    best_ov = 0
    best_row = None
    for _, r in sub.iterrows():
        s = int(r["start"]) if pd.notna(r["start"]) else None
        e = int(r["end"]) if pd.notna(r["end"]) else None
        if s is None or e is None:
            continue
        ov = _overlap_len(g0, g1, s, e)
        if ov > best_ov:
            best_ov = ov
            best_row = r
    if best_row is None or best_ov <= 0:
        return 0, None, None, None, None
    frac = float(best_ov) / float(max(1, gene_len))
    return (
        int(best_ov),
        int(best_row["cn_minor"]) if pd.notna(best_row["cn_minor"]) else None,
        int(best_row["cn_major"]) if pd.notna(best_row["cn_major"]) else None,
        int(best_row["cn_total"]) if pd.notna(best_row["cn_total"]) else None,
        float(frac),
    )


def hla_loh_from_segments_for_genes(
    segs: pd.DataFrame,
    gene_intervals: Dict[str, Tuple[str, int, int]],
    *,
    min_overlap_frac: float = 0.30,
) -> Dict[str, object]:
    """
    ``gene_intervals``: symbol -> (chrom, start, end) in bp, inclusive.

    LOH call: best-overlap segment has ``cn_minor == 0`` and overlap covers at least
    ``min_overlap_frac`` of the gene interval length.
    """
    out: Dict[str, object] = {}
    any_loh = False
    for sym, (ch, g0, g1) in gene_intervals.items():
        ov, mn, mj, tot, frac = best_segment_for_interval(segs, ch, int(g0), int(g1))
        loh = bool(
            mn == 0
            and ov > 0
            and frac is not None
            and frac >= float(min_overlap_frac)
        )
        if loh:
            any_loh = True
        out[f"{sym}_seg_overlap_bp"] = int(ov)
        out[f"{sym}_seg_minor"] = mn
        out[f"{sym}_seg_major"] = mj
        out[f"{sym}_seg_total_cn"] = tot
        out[f"{sym}_seg_overlap_frac"] = frac
        out[f"{sym}_seg_loh"] = bool(loh)
    out["hla_loh_seg_any"] = bool(any_loh)
    return out


def load_hla_gene_intervals_from_gene_table(path: Path) -> Dict[str, Tuple[str, int, int]]:
    """Read HLA-A/B/C + B2M intervals from a per-sample ASCAT gene-level table."""
    head = pd.read_csv(path, nrows=0)
    usecols = [c for c in ("gene_name", "chromosome", "start", "end") if c in head.columns]
    if len(usecols) < 4:
        return {}
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df["gene_name"] = df["gene_name"].astype(str)
    want = {"HLA-A", "HLA-B", "HLA-C", "B2M"}
    sub = df[df["gene_name"].isin(want)].drop_duplicates(subset=["gene_name"], keep="first")
    out: Dict[str, Tuple[str, int, int]] = {}
    for sym in sorted(want):
        row = sub.loc[sub["gene_name"] == sym]
        if row.empty:
            continue
        r = row.iloc[0]
        ch = str(r["chromosome"]).strip()
        if ch and not ch.lower().startswith("chr"):
            ch = f"chr{ch}"
        g0 = int(pd.to_numeric(r["start"], errors="coerce"))
        g1 = int(pd.to_numeric(r["end"], errors="coerce"))
        out[sym] = (ch, g0, g1)
    return out


def load_hla_types_by_aliquot(path: Path) -> pd.DataFrame:
    p = Path(path)
    if not p.is_file():
        return pd.DataFrame()
    df = pd.read_csv(p, low_memory=False)
    if "aliquot_id" not in df.columns:
        return pd.DataFrame()
    return df


def match_hla_types_row(hla_df: pd.DataFrame, tumor_vial: str) -> Optional[pd.Series]:
    if hla_df.empty or not tumor_vial:
        return None
    v = str(tumor_vial).strip()
    # Longest-prefix match: tumor vial is a prefix of GDC RNA/DNA aliquot barcodes.
    m = hla_df[hla_df["aliquot_id"].astype(str).str.startswith(v)]
    if m.empty:
        return None
    return m.iloc[0]
