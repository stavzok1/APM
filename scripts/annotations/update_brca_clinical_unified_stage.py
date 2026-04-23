#!/usr/bin/env python3
"""Add pathologic_stage_collapsed to BRCA_clinical_immune_unified.tsv (idempotent if column exists)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.clinical_labels import add_pathologic_stage_collapsed_column  # noqa: E402
from pipeline.config import PATHS  # noqa: E402


def main() -> None:
    p = PATHS.brca_clinical_immune_unified
    df = pd.read_csv(p, sep="\t", low_memory=False)
    col = "pathologic_stage_collapsed"
    if col not in df.columns:
        df[col] = add_pathologic_stage_collapsed_column(df, source_col="pathologic_stage")
    else:
        df[col] = add_pathologic_stage_collapsed_column(df, source_col="pathologic_stage")
    df.to_csv(p, sep="\t", index=False)
    print(f"Updated {p} with column {col} ({df[col].astype(str).replace('', pd.NA).notna().sum()} non-empty values)")


if __name__ == "__main__":
    main()
