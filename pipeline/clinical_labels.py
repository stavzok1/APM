"""Helpers for TCGA BRCA clinical labels used in annotations and downstream analysis."""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd


def collapse_pathologic_stage(val: object) -> str:
    """
    Map detailed AJCC-style strings (e.g. Stage IIA, Stage IIIB) to coarse Roman stages I–IV, plus X for unknown / Stage X.

    Examples:
        Stage IIA -> II
        Stage IA -> I
        Stage IIIB -> III
        Stage IV -> IV
        Stage X -> X
    """
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    s = str(val).strip()
    if not s or s.lower() in ("nan", "none", "na", "[not available]", "[discrepancy]", "not reported"):
        return ""
    u = s.upper().replace(" ", "")
    if "STAGEX" in u or re.search(r"STAGE\s*X\b", s.upper()):
        return "X"
    # Strip leading STAGE for suffix checks
    rest = u
    if rest.startswith("STAGE"):
        rest = rest[5:]
    if not rest:
        return ""
    if rest.startswith("IV") or "STAGEIV" in u:
        return "IV"
    if rest.startswith("III") or "STAGEIII" in u:
        return "III"
    if rest.startswith("II"):
        return "II"
    if rest.startswith("I"):
        return "I"
    return ""


def add_pathologic_stage_collapsed_column(df: pd.DataFrame, *, source_col: str = "pathologic_stage") -> pd.Series:
    """Return a series aligned to df with collapsed stages (empty string if unknown)."""
    if source_col not in df.columns:
        return pd.Series([""] * len(df), index=df.index)
    return df[source_col].map(collapse_pathologic_stage)
