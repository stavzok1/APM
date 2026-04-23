"""
Canonical labels for cell lines / biosamples used as dict keys and columns.

ENCODE, SCREEN, ChIP-atlas, and local folder names use different spellings
(MCF-7 vs MCF7, T-47D vs T47D). Pipeline outputs (SCREEN `per_biosample`, ABC
`ABC_full`, HiChIP columns, cCRE signal columns, ChIP `chip_hits`) use
**canonical** keys from this module.

**Disk vs data:**
- Paths (cCRE signals, HiChIP loops) still use `BiosampleConfig` directory names
  (`ccre_signal_cell_lines`, `hichip_panel`).
- SCREEN `Biosample` and ABC `CellType` values in downloaded files are mapped
  via `canonicalize_screen_biosample` / `canonicalize_abc_cell_type` before
  aggregation and nested dict construction.
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Set


def _compact_key(name: str) -> str:
    """Lowercase alphanumerics only (strip punctuation / spaces / underscores)."""
    return re.sub(r"[^a-z0-9]", "", str(name).strip().lower())


# Maps compact_key -> canonical display token used in pipeline outputs
_CELL_LINE_CANONICAL: Dict[str, str] = {
    _compact_key("MCF-7"): "MCF7",
    _compact_key("MCF7"): "MCF7",
    _compact_key("T-47D"): "T47D",
    _compact_key("T47D"): "T47D",
    _compact_key("SK-BR-3"): "SKBR3",
    _compact_key("SKBR3"): "SKBR3",
    _compact_key("MCF-10A"): "MCF10A",
    _compact_key("MCF10A"): "MCF10A",
    _compact_key("MCF_10A"): "MCF10A",
    _compact_key("MDA-MB-231"): "MDA-MB-231",
    _compact_key("MDA-MB-468"): "MDA-MB-468",
    _compact_key("HCC1954"): "HCC1954",
    _compact_key("ZR-75-1"): "ZR751",
    _compact_key("ZR751"): "ZR751",
    _compact_key("HMEC"): "HMEC",
    _compact_key("HMEC1"): "HMEC1",
    _compact_key("HMEC2"): "HMEC2",
    _compact_key("breast_tissue"): "breast_tissue",
    _compact_key("breast tissue"): "breast_tissue",
}


def normalize_cell_line_label(name: str) -> str:
    """
    Return a canonical cell-line label for use as nested-dict keys / columns.

    Unknown names are returned stripped (no aggressive guess).
    """
    if name is None or (isinstance(name, float) and str(name) == "nan"):
        return ""
    s = str(name).strip()
    if not s:
        return s
    k = _compact_key(s)
    return _CELL_LINE_CANONICAL.get(k, s)


def normalize_cell_line_labels(names) -> list:
    """Normalize an iterable of labels (dedupe preserves order)."""
    out: list = []
    seen = set()
    for n in names:
        c = normalize_cell_line_label(str(n))
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def slugify_biosample_key(name: str) -> str:
    """Stable ASCII key for tissue / long ENCODE biosample labels (lowercase snake)."""
    t = re.sub(r"[^a-zA-Z0-9]+", "_", str(name).strip().lower())
    t = re.sub(r"_+", "_", t).strip("_")
    return t or "unknown"


# Exact strings as they appear in ENCODE SCREEN "Biosample" → canonical output / join keys
SCREEN_UPSTREAM_TO_CANONICAL: Dict[str, str] = {
    "MCF-7": "MCF7",
    "MCF_10A": "MCF10A",
    "Breast Mammary Tissue": "breast_mammary_tissue",
    "breast_epithelium_tissue_female_adult_(53_years)": "breast_epithelium_tissue_female_adult_53_years",
    "mammary_epithelial_cell_female_adult_(19_years)": "mammary_epithelial_cell_female_adult_19_years",
    "T47D": "T47D",
    "mammary_epithelial_cell": "mammary_epithelial_cell",
    "breast_epithelium": "breast_epithelium",
}


def canonicalize_screen_biosample(raw: str) -> str:
    """Map SCREEN file biosample text to a canonical key (columns, per_biosample, conservation)."""
    if raw is None or (isinstance(raw, float) and str(raw) == "nan"):
        return ""
    s = str(raw).strip()
    if not s:
        return s
    if s in SCREEN_UPSTREAM_TO_CANONICAL:
        return SCREEN_UPSTREAM_TO_CANONICAL[s]
    ck = _compact_key(s)
    if ck in _CELL_LINE_CANONICAL:
        return _CELL_LINE_CANONICAL[ck]
    return slugify_biosample_key(s)


# Exact ABC file "CellType" values used in BiosampleConfig → nested ABC_full keys
ABC_UPSTREAM_TO_CANONICAL: Dict[str, str] = {
    "MCF-7-ENCODE": "MCF7_ENCODE",
    "MCF10A-Ji2017": "MCF10A_Ji2017",
    "MCF10A_treated_with_TAM24hr-Ji2017": "MCF10A_TAM24hr_Ji2017",
    "MDA-MB-231": "MDA-MB-231",
    "mammary_epithelial_cell-Roadmap": "mammary_epithelial_cell_Roadmap",
    "breast_epithelium-ENCODE": "breast_epithelium_ENCODE",
}


def canonicalize_abc_cell_type(raw: str) -> str:
    """Map ABC CellType strings to canonical keys inside ABC_full dicts."""
    if raw is None or (isinstance(raw, float) and str(raw) == "nan"):
        return ""
    s = str(raw).strip()
    if not s:
        return s
    if s in ABC_UPSTREAM_TO_CANONICAL:
        return ABC_UPSTREAM_TO_CANONICAL[s]
    ck = _compact_key(s)
    if ck in _CELL_LINE_CANONICAL:
        return _CELL_LINE_CANONICAL[ck]
    return slugify_biosample_key(s)


def canonical_hichip_output_key(name: str) -> str:
    """HiChIP folder name (disk) → key/column name in elem_focus / gene_links (usually unchanged)."""
    if name is None or (isinstance(name, float) and str(name) == "nan"):
        return ""
    s = str(name).strip()
    if not s:
        return s
    n = normalize_cell_line_label(s)
    return n if n else s


def canonical_ccre_signal_column_name(disk_subdir: str) -> str:
    """cCRE signal subdirectory name on disk → DataFrame column name for that cell line."""
    s = str(disk_subdir).strip()
    if not s:
        return s
    n = normalize_cell_line_label(s)
    if n:
        return n
    return slugify_biosample_key(s)


def ccre_signal_output_columns(cell_lines_on_disk: Iterable[str]) -> List[str]:
    """Ordered unique output column names for cCRE cell-line dict columns (from disk dir list)."""
    out: List[str] = []
    seen: Set[str] = set()
    for d in cell_lines_on_disk:
        c = canonical_ccre_signal_column_name(str(d))
        if c and c not in seen:
            seen.add(c)
            out.append(c)
    return out
