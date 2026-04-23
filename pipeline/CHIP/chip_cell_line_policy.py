"""
ChIP-Atlas cell-line filtering and subtype labels (breast / mammary concordance).

Applied in ``load_unified_chip`` to **CHIP_ATLAS** rows only (ENCODE unchanged).
Uses ``BIOSAMPLES.chip_brca_celltypes`` plus tissue keywords and optional subtype map.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..biosample_names import normalize_cell_line_label

if TYPE_CHECKING:
    from ..config import BiosampleConfig


def _cell_matches_mammary_panel(norm: str, raw: str, biosamples: "BiosampleConfig") -> bool:
    canon = {normalize_cell_line_label(x) for x in biosamples.chip_brca_celltypes}
    if norm in canon:
        return True
    ru = str(raw).lower()
    if "mammary" in ru or "breast" in ru:
        return True
    for kw in biosamples.chip_atlas_mammary_tissue_keywords:
        if kw.lower() in ru:
            return True
    return False


def attach_chip_cell_subtype(cell_series: pd.Series, biosamples: "BiosampleConfig") -> pd.Series:
    """Map normalized cell line label to subtype/state string (empty if unknown)."""
    out = []
    mapping = {normalize_cell_line_label(k): v for k, v in biosamples.chip_cell_line_subtype_map.items()}
    for raw in cell_series.astype(str):
        n = normalize_cell_line_label(raw)
        out.append(mapping.get(n, ""))
    return pd.Series(out, index=cell_series.index, dtype=object)


def filter_chip_atlas_by_mammary_policy(
    df: pd.DataFrame,
    biosamples: "BiosampleConfig",
) -> pd.DataFrame:
    """
    Keep only rows whose ``cell_type`` matches BRCA panel / mammary / breast tissue keywords.
    """
    if df.empty or "cell_type" not in df.columns:
        return df
    raw = df["cell_type"].astype(str)
    norm = raw.map(normalize_cell_line_label)
    mask = pd.Series(
        [_cell_matches_mammary_panel(n, r, biosamples) for n, r in zip(norm, raw)],
        index=df.index,
    )
    return df.loc[mask].copy()


def apply_chip_atlas_cell_line_policy(
    chip_atlas: pd.DataFrame,
    biosamples: "BiosampleConfig",
) -> pd.DataFrame:
    """Filter + ``cell_subtype`` column for ChIP-Atlas unified rows."""
    if chip_atlas.empty:
        chip_atlas = chip_atlas.copy()
        if "cell_subtype" not in chip_atlas.columns:
            chip_atlas["cell_subtype"] = pd.Series(dtype=object)
        return chip_atlas
    df = chip_atlas.copy()
    n0 = len(df)
    df = filter_chip_atlas_by_mammary_policy(df, biosamples)
    n1 = len(df)
    print(
        f"  ChIP-Atlas mammary/tissue filter: retained {n1:,}/{n0:,} peaks "
        f"({100.0 * n1 / max(n0, 1):.1f}%)"
    )
    df["cell_subtype"] = attach_chip_cell_subtype(df["cell_type"], biosamples)
    return df
