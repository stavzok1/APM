"""
Pipeline-wide gene symbol normalization (HGNC aliases + UCSC + legacy fixes).

Use the same mapping everywhere external datasets may use alternate symbols
while ``pipeline/config.py`` uses GENCODE/HGNC-approved names.
"""

from __future__ import annotations

import os
from functools import lru_cache
from typing import Collection, Dict, List, Mapping, Optional, Tuple

import pandas as pd

from .gene_aliases import strip_ensembl_version
from .panel_alias_registry import get_gene_flat_alias_map


def _max_lines_cache_key(max_lines: Optional[int]) -> int:
    """lru_cache needs hashable key; reserve -1 for None (full scan)."""
    return -1 if max_lines is None else int(max_lines)


@lru_cache(maxsize=8)
def _panel_symbol_mapping_cached(max_lines_key: int) -> tuple:
    """Thin cache over ``panel_alias_registry.get_gene_flat_alias_map``."""
    max_lines: Optional[int] = None if max_lines_key < 0 else max_lines_key
    return (get_gene_flat_alias_map(max_lines),)


def get_panel_symbol_mapping(max_lines: Optional[int] = None) -> Dict[str, str]:
    """
    Merged old_symbol -> canonical_symbol map for the full panel seed.

    Results are cached per ``max_lines`` (``None`` = scan entire HGNC table).
    """
    return dict(_panel_symbol_mapping_cached(_max_lines_cache_key(max_lines))[0])


def default_symbol_mapping() -> Dict[str, str]:
    """Full HGNC scan (slow once per distinct process + cache key)."""
    return get_panel_symbol_mapping(None)


def symbol_mapping_for_tests(max_lines: int) -> Dict[str, str]:
    """Small HGNC scan for fast unit/smoke tests (uses registry cache key)."""
    return dict(get_gene_flat_alias_map(max_lines))


def apply_symbol_mapping_series(series: pd.Series, mapping: Mapping[str, str]) -> pd.Series:
    """Replace values using ``mapping`` (unknown symbols unchanged; NaN preserved)."""
    if not len(mapping):
        return series

    def _one(x: object) -> object:
        if pd.isna(x):
            return x
        s = str(x).strip()
        return mapping.get(s, s)

    return series.map(_one)


def normalize_annotation_gene_names(
    df: pd.DataFrame,
    columns: Tuple[str, ...] = ("gene_name",),
    *,
    mapping: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Copy-on-write: return a dataframe with HGNC-style remaps applied to symbol columns.

    Used for GENCODE-style gene / lncRNA annotation tables (CNV, SV, ATAC, methylation).
    Disabled when ``APM_USE_GENE_SYMBOL_MAPPING`` is ``0`` / ``false`` / ``no``.
    """
    if df is None or getattr(df, "empty", True):
        return df
    if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() in (
        "0",
        "false",
        "no",
    ):
        return df
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return df
    out = df.copy()
    m = mapping if mapping is not None else default_symbol_mapping()
    for c in cols:
        out[c] = apply_symbol_mapping_series(out[c], m)
    return out


def resolve_symbol_to_panel(
    symbol: Optional[str],
    panel: Collection[str],
    mapping: Mapping[str, str],
) -> Optional[str]:
    """
    If ``symbol`` is in ``panel``, or ``mapping[symbol]`` is in ``panel``,
    return the canonical panel symbol; else None.
    """
    if symbol is None or (isinstance(symbol, float) and pd.isna(symbol)):
        return None
    s = str(symbol).strip()
    if not s or s == "nan":
        return None
    panel_set = set(panel)
    if s in panel_set:
        return s
    mapped = mapping.get(s)
    if mapped and mapped in panel_set:
        return mapped
    if s.startswith("ENS"):
        base = strip_ensembl_version(s)
        mid = mapping.get(s) or mapping.get(base)
        if mid and mid in panel_set:
            return mid
    return None


def normalize_gene_name_list(tokens: List[str], mapping: Mapping[str, str]) -> List[str]:
    """Apply mapping token-wise; drop empties; dedupe while preserving order."""
    out: List[str] = []
    seen: set[str] = set()
    for t in tokens:
        if t is None or (isinstance(t, float) and pd.isna(t)):
            continue
        u = str(t).strip()
        if not u:
            continue
        v = mapping.get(u, u)
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out
