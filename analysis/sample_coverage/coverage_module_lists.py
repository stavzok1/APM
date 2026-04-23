"""
Shared module name lists for sample coverage + clinical stratification.

Kept in one place to avoid circular imports between
``analysis/sample_module_coverage.py`` and ``analysis/clinical_omics_stratification.py``.
"""

from __future__ import annotations

# Column order for presence matrices (vial-level omics first)
VIAL_MODULES = [
    "SNV",
    "SV",
    "CNV",
    "Methylation",
    "RPPA",
    "miRNA",
    "RNA",
    "ATAC",
    "HLA",
]

PARTICIPANT_EXTRA = [
    "Immune_advanced",
    "Immune_thornsson",
    "Clinical_unified",
]

# Participant-level assay modules (not keyed at sample_vial in available inputs)
PARTICIPANT_OMICS_ONLY = ["HiCHIP"]

PARTICIPANT_MODULES = VIAL_MODULES + PARTICIPANT_OMICS_ONLY + PARTICIPANT_EXTRA

# Omics = assay modules keyed at sample_vial; metadata = participant-level annotations only.
OMICS_MODULES = list(VIAL_MODULES) + list(PARTICIPANT_OMICS_ONLY)
METADATA_MODULES = list(PARTICIPANT_EXTRA)
