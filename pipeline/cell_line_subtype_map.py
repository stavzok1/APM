"""
Biology-based breast model labels for cell lines and biosample keys.

PAM50-aligned **tumor** groups used here: ``LumA``, ``LumB``, ``HER2``, ``Basal``
(TADs uses ``Her2`` spelling; we normalize to ``HER2`` for the four-class label).

**Non-tumor** buckets:
- ``Normal_like`` — benign / non-tumorigenic cultures (e.g. HMEC, MCF10A) per TADs metadata.
- ``Healthy_tissue`` — primary or organoid breast epithelium labels from ENCODE / Roadmap-style
  biosample strings (SCREEN tissues, ``breast_epithelium_*``, etc.).

Source of truth for tumor lines: ``../TADs/config.py`` ``CELL_LINE_METADATA`` (copied and extended
for lines that appear in this pipeline but were not listed there).

Lookup accepts **canonical** pipeline keys (``MCF7``, ``MCF7_ENCODE``, ``breast_mammary_tissue``, …)
as produced by ``biosample_names`` / ``BiosampleConfig``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .biosample_names import normalize_cell_line_label


# ---------------------------------------------------------------------------
# Subtype vocabulary (tumor four-class + non-tumor)
# ---------------------------------------------------------------------------

PAM50_TUMOR_GROUPS = ("LumA", "LumB", "HER2", "Basal")
NON_TUMOR_GROUPS = ("Normal_like", "Healthy_tissue")
UNKNOWN_GROUP = "Unknown"


@dataclass(frozen=True)
class SubtypeRecord:
    """Resolved biology labels for one biosample / cell-line key."""

    input_key: str
    pam50_group: str
    er_status: str
    her2_status: str
    tnbc: Optional[bool]
    source: str  # "tads_cell_line" | "apm_extension" | "tissue" | "abc_composite" | "unknown"


# Normalized from TADs ``CellLineInfo`` (pam50_subtype Her2 -> HER2).
_CELL_LINES: Dict[str, SubtypeRecord] = {
    "HMEC": SubtypeRecord(
        "HMEC", "Normal_like", "N/A", "N/A", False, "tads_cell_line",
    ),
    "MCF10A": SubtypeRecord(
        "MCF10A", "Normal_like", "ER-", "HER2-", False, "tads_cell_line",
    ),
    "MCF7": SubtypeRecord(
        "MCF7", "LumA", "ER+", "HER2-", False, "tads_cell_line",
    ),
    "T47D": SubtypeRecord(
        "T47D", "LumA", "ER+", "HER2-", False, "tads_cell_line",
    ),
    "BT474": SubtypeRecord(
        "BT474", "LumB", "ER+", "HER2+", False, "tads_cell_line",
    ),
    "SKBR3": SubtypeRecord(
        "SKBR3", "HER2", "ER-", "HER2+", False, "tads_cell_line",
    ),
    "ZR751": SubtypeRecord(
        "ZR751", "LumA", "ER+", "HER2-", False, "tads_cell_line",
    ),
    "ZR7530": SubtypeRecord(
        "ZR7530", "LumA", "ER+", "HER2-", False, "tads_cell_line",
    ),
    "HCC1954": SubtypeRecord(
        "HCC1954", "Basal", "ER-", "HER2+", False, "tads_cell_line",
    ),
    "HCC70": SubtypeRecord(
        "HCC70", "Basal", "ER-", "HER2-", True, "tads_cell_line",
    ),
    "BT549": SubtypeRecord(
        "BT549", "Basal", "ER-", "HER2-", True, "tads_cell_line",
    ),
    "MDA-MB-231": SubtypeRecord(
        "MDA-MB-231", "Basal", "ER-", "HER2-", True, "tads_cell_line",
    ),
    # Present in APM ChIP whitelist / cCRE but not in TADs CELL_LINE_METADATA
    "MDA-MB-468": SubtypeRecord(
        "MDA-MB-468", "Basal", "ER-", "HER2-", True, "apm_extension",
    ),
    "HMEC1": SubtypeRecord(
        "HMEC1", "Normal_like", "N/A", "N/A", False, "apm_extension",
    ),
    "HMEC2": SubtypeRecord(
        "HMEC2", "Normal_like", "N/A", "N/A", False, "apm_extension",
    ),
}

_TISSUE_OR_NONMALIGNANT_KEY: Dict[str, SubtypeRecord] = {
    # SCREEN / ENCODE-style canonical biosample keys (see biosample_names)
    "breast_mammary_tissue": SubtypeRecord(
        "breast_mammary_tissue", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
    "breast_epithelium_tissue_female_adult_53_years": SubtypeRecord(
        "breast_epithelium_tissue_female_adult_53_years",
        "Healthy_tissue",
        "N/A",
        "N/A",
        False,
        "tissue",
    ),
    "mammary_epithelial_cell_female_adult_19_years": SubtypeRecord(
        "mammary_epithelial_cell_female_adult_19_years",
        "Healthy_tissue",
        "N/A",
        "N/A",
        False,
        "tissue",
    ),
    "mammary_epithelial_cell": SubtypeRecord(
        "mammary_epithelial_cell", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
    "breast_epithelium": SubtypeRecord(
        "breast_epithelium", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
    "mammary_epithelial_cell_Roadmap": SubtypeRecord(
        "mammary_epithelial_cell_Roadmap", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
    "breast_epithelium_ENCODE": SubtypeRecord(
        "breast_epithelium_ENCODE", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
    "breast_tissue": SubtypeRecord(
        "breast_tissue", "Healthy_tissue", "N/A", "N/A", False, "tissue",
    ),
}

# ABC / composite keys where the cell-line token is not a simple prefix of the key
_ABC_COMPOSITE_TO_LINE: Dict[str, str] = {
    "MCF7_ENCODE": "MCF7",
    "MCF10A_Ji2017": "MCF10A",
    "MCF10A_TAM24hr_Ji2017": "MCF10A",
    "MDA-MB-231": "MDA-MB-231",
    "mammary_epithelial_cell_Roadmap": "mammary_epithelial_cell_Roadmap",
    "breast_epithelium_ENCODE": "breast_epithelium_ENCODE",
}

# HiChIP folder names not in TADs table — conservative assignments where well known
_HICHIP_EXTRA: Dict[str, SubtypeRecord] = {
    # Luminal progenitor / basal-enriched culture names used in some HiChIP panels
    "B80T5": SubtypeRecord("B80T5", "Basal", "N/A", "N/A", None, "apm_extension"),
    "K5plusK19plus": SubtypeRecord(
        "K5plusK19plus", "Normal_like", "N/A", "N/A", False, "apm_extension",
    ),
}

# Inverted index: PAM50 group -> cell lines (from TADs PAM50_TO_CELL_LINES, HER2 key normalized)
PAM50_GROUP_TO_CELL_LINES: Dict[str, List[str]] = {
    "Basal": ["HCC70", "BT549", "MDA-MB-231", "HCC1954", "MDA-MB-468", "B80T5"],
    "HER2": ["SKBR3", "HCC1954"],
    "LumA": ["MCF7", "T47D", "ZR751", "ZR7530"],
    "LumB": ["BT474"],
    "Normal_like": ["MCF10A", "HMEC", "HMEC1", "HMEC2", "K5plusK19plus"],
}


def _unknown_record(key: str) -> SubtypeRecord:
    return SubtypeRecord(
        key, UNKNOWN_GROUP, "N/A", "N/A", None, "unknown",
    )


def _resolve_cell_line_token(token: str) -> Optional[SubtypeRecord]:
    t = str(token).strip()
    if not t:
        return None
    if t in _CELL_LINES:
        return _CELL_LINES[t]
    norm = normalize_cell_line_label(t)
    if norm and norm in _CELL_LINES:
        return _CELL_LINES[norm]
    if norm and norm in _HICHIP_EXTRA:
        return _HICHIP_EXTRA[norm]
    if t in _HICHIP_EXTRA:
        return _HICHIP_EXTRA[t]
    return None


def _abc_like_key_to_line(key: str) -> str:
    """Map canonical ABC / composite keys to a cell-line token for ``_CELL_LINES`` lookup."""
    k = str(key).strip()
    if k in _ABC_COMPOSITE_TO_LINE:
        return _ABC_COMPOSITE_TO_LINE[k]
    m = re.match(
        r"^([A-Za-z0-9][A-Za-z0-9+\-]*?)_(ENCODE|Ji2017|Roadmap|TAM24hr.*)$",
        k,
    )
    if m:
        return m.group(1)
    return k


def subtype_for_key(key: str) -> SubtypeRecord:
    """
    Return PAM50 / tissue group for a canonical biosample or cell-line key.

    Accepts keys such as ``MCF7``, ``MCF7_ENCODE``, ``breast_mammary_tissue``,
    ``HMEC1``, HiChIP folder names, etc.
    """
    k = str(key).strip()
    if not k:
        return _unknown_record(k)

    if k in _TISSUE_OR_NONMALIGNANT_KEY:
        return _TISSUE_OR_NONMALIGNANT_KEY[k]

    if k in _HICHIP_EXTRA:
        return _HICHIP_EXTRA[k]

    line_guess = _abc_like_key_to_line(k)
    hit = _resolve_cell_line_token(line_guess)
    if hit:
        src = "abc_composite" if k in _ABC_COMPOSITE_TO_LINE else hit.source
        return SubtypeRecord(
            k,
            hit.pam50_group,
            hit.er_status,
            hit.her2_status,
            hit.tnbc,
            src,
        )

    hit = _resolve_cell_line_token(k)
    if hit:
        return SubtypeRecord(
            k, hit.pam50_group, hit.er_status, hit.her2_status, hit.tnbc, hit.source,
        )

    return _unknown_record(k)


def subtype_table_for_keys(keys: List[str]) -> List[Tuple[str, SubtypeRecord]]:
    """Ordered (key, record) pairs for reporting."""
    return [(k, subtype_for_key(k)) for k in keys]
