"""
Annotate cCREs with ChIP-seq peak "hits" from a unified peak table.

The unified table is produced by `pipeline.CHIP.chip_loader.load_unified_chip`
and cached at `PATHS.chip_unified` with schema:
    chrom, start, end, tf, cell_type, source, score_norm, sample_id

This module overlaps ChIP peaks with cCRE intervals and produces a compact
per-cCRE summary suitable for storing as a complex column in the
regulatory element focus table.

**Cell-type filtering:** unlike SV ``map_svs_to_chip`` (which can apply
``cell_type_whitelist`` / ``tf_whitelist`` to the unified peak table), cCRE
``chip_hits`` here use **all** peaks in ``chip_unified``—only chromosome
restriction and coordinate validity apply.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from ..biosample_names import normalize_cell_line_label


def _compute_overlap_bp(a_start: np.ndarray, a_end: np.ndarray, b_start: np.ndarray, b_end: np.ndarray) -> np.ndarray:
    """Compute base-pair overlap for interval pairs (0 if no overlap)."""
    left = np.maximum(a_start, b_start)
    right = np.minimum(a_end, b_end)
    return np.maximum(0, right - left).astype(np.int64)


def build_chip_hits_for_ccres(
    ccres: pd.DataFrame,
    chip_unified: Union[str, Path, pd.DataFrame],
    *,
    bin_size: int = 50_000,
    chrom_col: str = "chrom",
    ccre_id_col: str = "cCRE_id",
    encode_id_col: str = "ENCODE_id",
) -> pd.DataFrame:
    """
    Build compact hierarchical `chip_hits` for each cCRE by overlapping ChIP peaks.

    Args:
        ccres: cCRE DataFrame with at least: chrom, start, end, cCRE_id/ENCODE_id
        chip_unified: parquet path or pre-loaded DataFrame with columns:
            chrom, start, end, tf, cell_type, source, score_norm, sample_id
        bin_size: coarse genomic bin size for reducing candidate pairs
        chrom_col: chromosome column name in both tables
        ccre_id_col: cCRE id column name
        encode_id_col: ENCODE id column name

    Returns:
        DataFrame with columns:
            cCRE_id, chip_hits
        where chip_hits is a nested dict keyed by TF → cell_type → source:

            chip_hits = {
              "CTCF": {
                "MCF7": {
                  "ENCODE": {
                    "n_peaks": 3,
                    "n_samples": 2,
                    "sample_ids": ["CTCF_MCF7_1", "CTCF_MCF7_2"],
                    "max_overlap_bp": 180,
                    "max_score_norm": 12.3,
                    "mean_score_norm": 7.1,
                    "per_sample": {
                      "CTCF_MCF7_1": {"n_peaks": 2, "max_overlap_bp": 180, "max_score_norm": 12.3, "mean_score_norm": 8.1},
                      "CTCF_MCF7_2": {"n_peaks": 1, "max_overlap_bp": 120, "max_score_norm": 9.0, "mean_score_norm": 9.0}
                    }
                  },
                  "CHIP_ATLAS": { ... }
                },
                ...
              },
              ...
            }
    """
    # --- load chip peaks ---
    if isinstance(chip_unified, (str, Path)):
        chip_path = Path(chip_unified)
        cols = [
            chrom_col,
            "start",
            "end",
            "tf",
            "cell_type",
            "source",
            "score_norm",
            "sample_id",
        ]
        try:
            chip = pd.read_parquet(chip_path, columns=cols)
        except Exception:
            # Back-compat: older cached parquet may lack `sample_id`.
            try:
                chip = pd.read_parquet(
                    chip_path,
                    columns=[
                        chrom_col,
                        "start",
                        "end",
                        "tf",
                        "cell_type",
                        "source",
                        "score_norm",
                    ],
                )
                chip["sample_id"] = pd.NA
            except Exception:
                chip = pd.read_parquet(chip_path)
                for c in cols:
                    if c not in chip.columns:
                        chip[c] = pd.NA
                chip = chip[cols].copy()
    else:
        chip = chip_unified.copy()
        if "sample_id" not in chip.columns:
            chip["sample_id"] = pd.NA

    # --- minimal input normalization ---
    c = ccres[[ccre_id_col, encode_id_col, chrom_col, "start", "end"]].copy()
    c["start"] = pd.to_numeric(c["start"], errors="coerce").astype("Int64")
    c["end"] = pd.to_numeric(c["end"], errors="coerce").astype("Int64")
    c = c.dropna(subset=[chrom_col, "start", "end", ccre_id_col])

    chip["start"] = pd.to_numeric(chip["start"], errors="coerce").astype("Int64")
    chip["end"] = pd.to_numeric(chip["end"], errors="coerce").astype("Int64")
    chip = chip.dropna(subset=[chrom_col, "start", "end", "tf", "cell_type", "source"])
    chip["cell_type"] = chip["cell_type"].map(normalize_cell_line_label)
    chip["sample_id"] = chip["sample_id"].astype(str).str.strip()

    # Filter chip peaks to chromosomes present in the cCRE set we are annotating
    chroms = set(c[chrom_col].dropna().unique().tolist())
    chip = chip[chip[chrom_col].isin(chroms)].copy()

    # --- coarse bin join to create candidate overlaps ---
    # cCRE bins
    c["ccre_start"] = c["start"].astype(np.int64)
    c["ccre_end"] = c["end"].astype(np.int64)
    c["bin_start"] = (c["ccre_start"] // bin_size).astype(np.int64)
    c["bin_end"] = (c["ccre_end"] // bin_size).astype(np.int64)

    cb = c[[ccre_id_col, chrom_col, "ccre_start", "ccre_end", "bin_start", "bin_end"]].copy()
    cb = cb.rename(columns={"bin_start": "bins"})
    cb = cb[[ccre_id_col, chrom_col, "ccre_start", "ccre_end", "bins"]]
    cross = c["bin_end"] != c["bin_start"]
    if cross.any():
        cb2 = c.loc[cross, [ccre_id_col, chrom_col, "ccre_start", "ccre_end", "bin_end"]].copy()
        cb2 = cb2.rename(columns={"bin_end": "bins"})
        cb = pd.concat([cb, cb2], ignore_index=True)

    # ChIP bins
    chip = chip.rename(columns={"start": "chip_start", "end": "chip_end"})
    chip["chip_start"] = chip["chip_start"].astype(np.int64)
    chip["chip_end"] = chip["chip_end"].astype(np.int64)
    chip["bin_start"] = (chip["chip_start"] // bin_size).astype(np.int64)
    chip["bin_end"] = (chip["chip_end"] // bin_size).astype(np.int64)

    pb = chip[[chrom_col, "chip_start", "chip_end", "tf", "cell_type", "source", "score_norm", "sample_id", "bin_start", "bin_end"]].copy()
    pb = pb.rename(columns={"bin_start": "bins"})
    pb = pb[[chrom_col, "chip_start", "chip_end", "tf", "cell_type", "source", "score_norm", "sample_id", "bins"]]
    cross_p = chip["bin_end"] != chip["bin_start"]
    if cross_p.any():
        pb2 = chip.loc[cross_p, [chrom_col, "chip_start", "chip_end", "tf", "cell_type", "source", "score_norm", "sample_id", "bin_end"]].copy()
        pb2 = pb2.rename(columns={"bin_end": "bins"})
        pb = pd.concat([pb, pb2], ignore_index=True)

    # Merge candidates and filter to true overlaps
    pref = cb.merge(pb, on=[chrom_col, "bins"], how="inner")
    pref = pref[
        (pref["chip_end"] >= pref["ccre_start"]) & (pref["chip_start"] <= pref["ccre_end"])
    ].copy()

    if pref.empty:
        ids = pd.Series(c[ccre_id_col].unique())
        return pd.DataFrame({ccre_id_col: ids, "chip_hits": [{} for _ in range(len(ids))]})

    # Exact overlap bp
    pref["overlap_bp"] = _compute_overlap_bp(
        pref["ccre_start"].to_numpy(),
        pref["ccre_end"].to_numpy(),
        pref["chip_start"].to_numpy(),
        pref["chip_end"].to_numpy(),
    )
    pref = pref[pref["overlap_bp"] > 0].copy()

    if pref.empty:
        ids = pd.Series(c[ccre_id_col].unique())
        return pd.DataFrame({ccre_id_col: ids, "chip_hits": [{} for _ in range(len(ids))]})

    # Per-sample summaries first (so we can record replicate/sample coverage).
    per_sample = (
        pref.groupby([ccre_id_col, "tf", "cell_type", "source", "sample_id"], sort=False)
        .agg(
            n_peaks=("overlap_bp", "size"),
            max_overlap_bp=("overlap_bp", "max"),
            max_score_norm=("score_norm", "max"),
            mean_score_norm=("score_norm", "mean"),
        )
        .reset_index()
    )

    # Aggregate to compact per-(cCRE, tf, cell_type, source) summaries (across all samples/replicates).
    agg = (
        pref.groupby([ccre_id_col, "tf", "cell_type", "source"], sort=False)
        .agg(
            n_peaks=("overlap_bp", "size"),
            n_samples=("sample_id", pd.Series.nunique),
            max_overlap_bp=("overlap_bp", "max"),
            max_score_norm=("score_norm", "max"),
            mean_score_norm=("score_norm", "mean"),
        )
        .reset_index()
    )

    # sample_id inventory per (cCRE, tf, cell_type, source)
    sample_ids = (
        per_sample.groupby([ccre_id_col, "tf", "cell_type", "source"], sort=False)["sample_id"]
        .apply(lambda s: sorted({str(x) for x in s.tolist() if str(x)}))
        .rename("sample_ids")
        .reset_index()
    )

    # Build nested per_sample dict per (cCRE, tf, cell_type, source)
    def _per_sample_dict(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        out: Dict[str, Dict[str, Any]] = {}
        for r in df_sub.itertuples(index=False):
            sid = str(r.sample_id)
            out[sid] = {
                "n_peaks": int(r.n_peaks),
                "max_overlap_bp": int(r.max_overlap_bp),
                "max_score_norm": None if pd.isna(r.max_score_norm) else float(r.max_score_norm),
                "mean_score_norm": None if pd.isna(r.mean_score_norm) else float(r.mean_score_norm),
            }
        return out

    per_sample_dict = (
        per_sample.groupby([ccre_id_col, "tf", "cell_type", "source"], sort=False)
        .apply(_per_sample_dict)
        .rename("per_sample")
        .reset_index()
    )

    agg = agg.merge(sample_ids, on=[ccre_id_col, "tf", "cell_type", "source"], how="left")
    agg = agg.merge(per_sample_dict, on=[ccre_id_col, "tf", "cell_type", "source"], how="left")
    agg["sample_ids"] = agg["sample_ids"].apply(lambda x: x if isinstance(x, list) else [])
    agg["per_sample"] = agg["per_sample"].apply(lambda x: x if isinstance(x, dict) else {})

    def _to_hierarchy(df_sub: pd.DataFrame) -> Dict[str, Dict[str, Dict[str, Dict[str, Any]]]]:
        """
        Convert per-cCRE aggregated rows into:
            {tf: {cell_type: {source: summary_dict}}}
        """
        out: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = {}
        for r in df_sub.itertuples(index=False):
            tf = str(r.tf)
            cell_type = normalize_cell_line_label(str(r.cell_type))
            source = str(r.source)
            out.setdefault(tf, {}).setdefault(cell_type, {})[source] = {
                "n_peaks": int(r.n_peaks),
                "n_samples": int(r.n_samples),
                "sample_ids": list(r.sample_ids) if isinstance(r.sample_ids, list) else [],
                "max_overlap_bp": int(r.max_overlap_bp),
                "max_score_norm": None if pd.isna(r.max_score_norm) else float(r.max_score_norm),
                "mean_score_norm": None if pd.isna(r.mean_score_norm) else float(r.mean_score_norm),
                "per_sample": r.per_sample if isinstance(r.per_sample, dict) else {},
            }
        return out

    chip_hits = (
        agg.groupby(ccre_id_col, sort=False)
        .apply(_to_hierarchy)
        .rename("chip_hits")
        .reset_index()
    )

    return chip_hits


def integrate_chip_hits_to_element_table(
    elem_focus: pd.DataFrame,
    ccres_full: pd.DataFrame,
    chip_unified: Union[str, Path, pd.DataFrame],
    *,
    ccre_id_col: str = "cCRE_id",
) -> pd.DataFrame:
    """
    Add `chip_hits` column to element focus table.

    elem_focus lacks genomic coordinates; we use `ccres_full` (with chrom/start/end)
    to compute overlaps and then merge back onto elem_focus by cCRE_id.
    """
    ef = elem_focus.copy()
    if ccre_id_col not in ef.columns:
        raise ValueError(f"elem_focus missing required column: {ccre_id_col}")

    wanted = set(ef[ccre_id_col].dropna().unique().tolist())
    cc = ccres_full[ccres_full[ccre_id_col].isin(wanted)].copy()

    hits_df = build_chip_hits_for_ccres(cc, chip_unified, ccre_id_col=ccre_id_col)
    ef = ef.merge(hits_df, on=ccre_id_col, how="left")
    ef["chip_hits"] = ef["chip_hits"].apply(lambda x: x if isinstance(x, dict) else {})
    return ef

