"""
QC package (pipeline-owned).

Design:
- Each QC module exposes:
  - compute(ctx, ...) -> (per_sample_df, per_gene_df, summary_dict)
  - assert_cohort(summary_dict, ...) -> list[Finding]

QC compute functions should be safe to run on partially-available modalities:
missing inputs produce empty outputs + warnings in the summary, not crashes.
"""

from .base import Finding, Level

"""
Generic QC utilities used by analysis sanity checks.

Keep this package modality-agnostic (e.g. NaN/empty checks). Expression-derived
scores live under `pipeline/RNA_exp/`.
"""

from .suspicious import (
    summarize_suspicious_columns,
    summarize_suspicious_nested_cells,
)

__all__ = [
    # suspicious / structural QC
    "summarize_suspicious_columns",
    "summarize_suspicious_nested_cells",
]

