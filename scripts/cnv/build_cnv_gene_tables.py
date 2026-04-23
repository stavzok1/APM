#!/usr/bin/env python3
"""
Build gene-centric CNV tables from per-segment annotated CNV CSVs.

Inputs are the outputs of `pipeline/CNV/runner.py` (per sample: *_cnv_annotated.csv),
which include `cn_total`, `cn_state`, `loh_flag`, and `gene_hits`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Allow running via `python3 scripts/...` without requiring installation.
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from pipeline.CNV.gene_summary import GeneCnvRules, summarize_gene_cnv  # noqa: E402


def _infer_sample_id(path: Path) -> str:
    # Expected: TCGA-XX-XXXX-01A_cnv_annotated.csv
    name = path.name
    if name.endswith("_cnv_annotated.csv"):
        return name[: -len("_cnv_annotated.csv")]
    return path.stem


def build_one(in_csv: Path, out_csv: Path, *, rules: GeneCnvRules) -> pd.DataFrame:
    df = pd.read_csv(in_csv)
    g = summarize_gene_cnv(df, rules=rules)
    g.insert(0, "sample_id", _infer_sample_id(in_csv))
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(out_csv, index=False)
    return g


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", type=str, required=True, help="Directory with *_cnv_annotated.csv files")
    ap.add_argument("--out-dir", type=str, required=True, help="Output directory for gene tables")
    ap.add_argument(
        "--min-body-overlap-pct",
        type=float,
        default=None,
        help="Minimum gene body overlap percent for calling loss/gain/amp. Default: pipeline THRESHOLDS.cnv_gene_min_body_overlap_pct",
    )
    ap.add_argument(
        "--combined",
        type=str,
        default="cnv_gene_calls_all_samples.csv",
        help="Filename for combined all-samples table (written under --out-dir). Set empty to skip.",
    )
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    min_pct = args.min_body_overlap_pct
    if min_pct is None:
        try:
            from pipeline.config import THRESHOLDS

            min_pct = float(getattr(THRESHOLDS, "cnv_gene_min_body_overlap_pct", 30.0))
        except Exception:
            min_pct = 30.0
    rules = GeneCnvRules(min_body_overlap_pct=float(min_pct))

    paths = sorted(in_dir.glob("*_cnv_annotated.csv"))
    if not paths:
        raise SystemExit(f"No '*_cnv_annotated.csv' files found under {in_dir}")

    combined = []
    for p in paths:
        sid = _infer_sample_id(p)
        out_csv = out_dir / f"{sid}_cnv_gene_calls.csv"
        print(f"[CNV gene] {sid} -> {out_csv}")
        g = build_one(p, out_csv, rules=rules)
        combined.append(g)

    if args.combined:
        all_df = pd.concat(combined, ignore_index=True) if combined else pd.DataFrame()
        out_all = out_dir / args.combined
        all_df.to_csv(out_all, index=False)
        print(f"[CNV gene] wrote combined -> {out_all} ({len(all_df)} rows)")


if __name__ == "__main__":
    main()

