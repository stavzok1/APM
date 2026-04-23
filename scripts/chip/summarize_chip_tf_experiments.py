#!/usr/bin/env python3
"""
Print TF × source experiment counts from unified ChIP peaks.

``experiment`` = one per distinct ``sample_id`` (BED filename stem). If the ``tf``
column is empty (common for some ENCODE rows where the 11th BED field is blank),
TF is inferred as the first ``_``-separated token of ``sample_id`` (e.g.
``CTCF_MCF7_1`` → ``CTCF``; ``EP300_MCF7`` → ``EP300``).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from pipeline.config import PATHS


def _tf_raw(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def _bad_tf_mask(tf: pd.Series) -> pd.Series:
    t = _tf_raw(tf)
    return t.eq("") | t.str.lower().isin(("nan", "none"))


def _infer_tf_from_sample_id(sample_id: pd.Series) -> pd.Series:
    return sample_id.astype(str).str.split("_").str[0].str.strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--chip-unified",
        type=Path,
        default=PATHS.chip_unified,
        help="unified_chip_peaks.parquet path",
    )
    ap.add_argument(
        "--list-samples",
        action="store_true",
        help="Print distinct sample_id per reporting-TF for each source (verbose).",
    )
    args = ap.parse_args()
    p = args.chip_unified
    if not p.is_file():
        print(f"Missing {p}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_parquet(p, columns=["tf", "source", "sample_id"])
    raw = _tf_raw(df["tf"])
    bad = _bad_tf_mask(df["tf"])
    inferred = _infer_tf_from_sample_id(df["sample_id"])
    df["tf_report"] = np.where(bad, inferred, raw)
    n_bad = int(bad.sum())
    if n_bad:
        print(
            f"[note] {n_bad:,} peak rows have empty/invalid `tf`; "
            f"using first token of sample_id for grouping (`tf_report`).",
            file=sys.stderr,
        )

    g = df.groupby(["tf_report", "source"])["sample_id"].nunique().unstack(fill_value=0)
    for col in ("ENCODE", "CHIP_ATLAS"):
        if col not in g.columns:
            g[col] = 0
    g = g[["ENCODE", "CHIP_ATLAS"]]
    g["total"] = g["ENCODE"] + g["CHIP_ATLAS"]
    g = g.sort_values("total", ascending=False)

    print(f"File: {p}")
    print("Experiment = distinct sample_id per source BED; grouped by tf_report (see stderr note).")
    print()
    print(g.to_string())
    print()
    print("Totals (distinct sample_id, any TF):")
    for src in ("ENCODE", "CHIP_ATLAS"):
        n = df.loc[df["source"] == src, "sample_id"].nunique()
        print(f"  {src}: {n}")

    if args.list_samples:
        for src in ("ENCODE", "CHIP_ATLAS"):
            sub = df[df["source"] == src]
            print()
            print(f"=== {src} sample_id by tf_report ===")
            for tf in sorted(sub["tf_report"].unique()):
                ids = sorted(sub.loc[sub["tf_report"] == tf, "sample_id"].astype(str).unique())
                print(f"  {tf}: {len(ids)} -> {ids}")


if __name__ == "__main__":
    main()
