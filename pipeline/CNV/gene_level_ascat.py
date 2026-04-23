"""
GDC / ASCAT3 gene-level copy-number TSV → per-sample gene tables.

Input files match names like ``*ascat3.gene_level_copy_number*.tsv`` with columns
``gene_id``, ``gene_name``, ``chromosome``, ``start``, ``end``,
``copy_number``, ``min_copy_number``, ``max_copy_number`` (numeric columns may
be empty for some genes).

Outputs mirror the segment-derived gene tables location (vial id stem) with
extra columns: discrete CN state, LOH flags from allele counts, optional
promoter coordinates from a gene annotation table, and a coarse regulatory
window minimum total copy number.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from .features import _classify_cn_state
from ..genes.gene_aliases import strip_ensembl_version
from ..genes.symbol_normalization import normalize_annotation_gene_names

# Filename substrings (case-insensitive) for GDC ASCAT3 gene-level exports.
ASCAT_GENE_LEVEL_MARKERS = (
    "ascat3.gene_level_copy_number",
    "gene_level_copy_number",
)

# Minimum columns; allele totals may be absent on sparse rows.
_REQUIRED_CORE = {"gene_id", "gene_name", "chromosome", "start", "end"}


def is_ascat_gene_level_filename(filename: str) -> bool:
    fn = filename.lower()
    return any(m in fn for m in ASCAT_GENE_LEVEL_MARKERS)


def _normalize_chr(c: Any) -> str:
    s = str(c).strip()
    if not s or s.lower() in ("na", "nan", "none"):
        return ""
    if s.lower().startswith("chr"):
        return s
    return f"chr{s}"


def _lower_colmap(df: pd.DataFrame) -> Dict[str, str]:
    return {c.lower(): c for c in df.columns}


def load_ascat_gene_level_tsv(path: str, sep: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, low_memory=False)
    cmap = _lower_colmap(df)
    missing = _REQUIRED_CORE - set(cmap.keys())
    if missing:
        raise ValueError(f"{path}: missing required columns {missing}; have {list(df.columns)}")
    out = pd.DataFrame(
        {
            "gene_id": df[cmap["gene_id"]].astype(str),
            "gene_name": df[cmap["gene_name"]].astype(str),
            "chromosome": df[cmap["chromosome"]].map(_normalize_chr),
            "start": pd.to_numeric(df[cmap["start"]], errors="coerce"),
            "end": pd.to_numeric(df[cmap["end"]], errors="coerce"),
        }
    )
    for opt in ("copy_number", "min_copy_number", "max_copy_number"):
        if opt in cmap:
            out[opt] = pd.to_numeric(df[cmap[opt]], errors="coerce")
        else:
            out[opt] = np.nan
    out = normalize_annotation_gene_names(out, ("gene_name",))
    return out


def _window_min_copy_number(
    starts: np.ndarray,
    ends: np.ndarray,
    cns: np.ndarray,
    window_bp: int,
) -> np.ndarray:
    """For each row on one chromosome, min copy_number over intervals overlapping expanded span."""
    n = len(starts)
    out = np.full(n, np.nan, dtype=float)
    if n == 0:
        return out
    order = np.argsort(starts)
    s = starts[order]
    e = ends[order]
    cn = cns[order]
    for ii, i in enumerate(order):
        qs = float(s[ii]) - window_bp
        qe = float(e[ii]) + window_bp
        j_lo = 0
        while j_lo < n and float(e[j_lo]) < qs:
            j_lo += 1
        best = np.nan
        j2 = j_lo
        while j2 < n and float(s[j2]) <= qe:
            if float(e[j2]) >= qs:
                v = cn[j2]
                if np.isfinite(v):
                    best = v if not np.isfinite(best) else min(best, v)
            j2 += 1
        out[i] = best
    return out


def _regulatory_window_min_cn(df: pd.DataFrame, window_bp: int) -> pd.Series:
    chroms = df["chromosome"].astype(str)
    starts = pd.to_numeric(df["start"], errors="coerce").to_numpy(dtype=float)
    ends = pd.to_numeric(df["end"], errors="coerce").to_numpy(dtype=float)
    cns = pd.to_numeric(df["copy_number"], errors="coerce").to_numpy(dtype=float)
    out = np.full(len(df), np.nan, dtype=float)
    for ch in chroms.unique():
        if not ch:
            continue
        m = (chroms == ch).to_numpy()
        if not m.any():
            continue
        out[m] = _window_min_copy_number(starts[m], ends[m], cns[m], window_bp)
    return pd.Series(out, index=df.index, dtype=float)


def _promoter_frame_from_genes(genes_df: pd.DataFrame) -> pd.DataFrame:
    """One row per (chrom, gene_name) with TSS / promoter interval when available."""
    g = genes_df.copy()
    if "feature" in g.columns:
        g = g[g["feature"].astype(str) == "gene"]
    chrom_col = "chrom" if "chrom" in g.columns else (
        "chromosome" if "chromosome" in g.columns else None
    )
    if chrom_col is None:
        return pd.DataFrame()
    g["_chrom"] = g[chrom_col].astype(str).map(_normalize_chr)
    g["_gene_upper"] = g["gene_name"].astype(str).str.upper()
    cols = ["_chrom", "_gene_upper"]
    for c in ("tss", "prom_start", "prom_end", "strand", "gene_id"):
        if c in g.columns:
            cols.append(c)
    sub = g[cols].drop_duplicates(subset=["_chrom", "_gene_upper"], keep="first")
    if "gene_id" in sub.columns:
        sub = sub.assign(_ens_base=sub["gene_id"].map(strip_ensembl_version))
    return sub


def enrich_ascat_gene_level(
    df: pd.DataFrame,
    *,
    genes_lookup: Optional[pd.DataFrame] = None,
    regulatory_window_bp: int = 250_000,
) -> pd.DataFrame:
    """
    Add ``cn_state``, LOH / allele columns, optional promoter fields, window min CN.
    """
    out = df.copy()
    out["cn_state"] = out["copy_number"].map(_classify_cn_state)

    mn = pd.to_numeric(out["min_copy_number"], errors="coerce")
    mx = pd.to_numeric(out["max_copy_number"], errors="coerce")
    minor = np.minimum(mn.to_numpy(dtype=float), mx.to_numpy(dtype=float))
    major = np.maximum(mn.to_numpy(dtype=float), mx.to_numpy(dtype=float))
    out["cn_minor"] = minor
    out["cn_major"] = major
    out["allele_delta"] = major - minor
    both = np.isfinite(minor) & np.isfinite(major)
    out["loh_flag"] = ((minor == 0.0) & (major >= 1.0) & both).astype(int)

    # Diploid total CN with one allele lost to zero (segment-style LOH on gene body).
    tot = pd.to_numeric(out["copy_number"], errors="coerce")
    out["loh_only"] = (
        (tot >= 1.5)
        & (tot <= 2.5)
        & (out["loh_flag"] == 1)
    ).astype(int)

    if regulatory_window_bp and regulatory_window_bp > 0:
        out["regulatory_window_min_cn"] = _regulatory_window_min_cn(out, int(regulatory_window_bp))
        out["regulatory_window_bp"] = int(regulatory_window_bp)
    else:
        out["regulatory_window_min_cn"] = np.nan
        out["regulatory_window_bp"] = 0

    if genes_lookup is not None and not genes_lookup.empty:
        pf = _promoter_frame_from_genes(genes_lookup)
        if not pf.empty:
            out["_chrom_key"] = out["chromosome"].astype(str)
            out["_gene_upper"] = out["gene_name"].astype(str).str.upper()
            out["_ens_base"] = out["gene_id"].map(strip_ensembl_version)
            pf2 = pf.rename(columns={"_chrom": "_chrom_key"})
            prom_cols = [c for c in ("tss", "prom_start", "prom_end", "strand") if c in pf2.columns]
            if prom_cols:
                m1 = pf2[["_chrom_key", "_gene_upper"] + prom_cols].drop_duplicates(
                    subset=["_chrom_key", "_gene_upper"], keep="first"
                )
                out = out.merge(m1, on=["_chrom_key", "_gene_upper"], how="left")
                if "_ens_base" in pf2.columns and out[prom_cols[0]].isna().any():
                    m2 = pf2.drop_duplicates(subset=["_ens_base"], keep="first")[
                        ["_ens_base"] + prom_cols
                    ]
                    out = out.merge(m2, on="_ens_base", how="left", suffixes=("", "_ensfb"))
                    for c in prom_cols:
                        alt = f"{c}_ensfb"
                        if alt in out.columns:
                            out[c] = out[c].where(out[c].notna(), out[alt])
                            out.drop(columns=[alt], inplace=True, errors="ignore")
            out.drop(columns=["_chrom_key", "_gene_upper", "_ens_base"], inplace=True, errors="ignore")

    drop_cols = [c for c in out.columns if str(c).startswith("_")]
    out = out.drop(columns=drop_cols, errors="ignore")
    return out


def build_ascat_gene_sample_table(
    raw: pd.DataFrame,
    sample_id: str,
    *,
    genes_lookup: Optional[pd.DataFrame] = None,
    regulatory_window_bp: int = 250_000,
) -> pd.DataFrame:
    """Return a wide per-gene table with ``sample_id`` first."""
    enriched = enrich_ascat_gene_level(
        raw,
        genes_lookup=genes_lookup,
        regulatory_window_bp=regulatory_window_bp,
    )
    enriched.insert(0, "sample_id", sample_id)
    return enriched
