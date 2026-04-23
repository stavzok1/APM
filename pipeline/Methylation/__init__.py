"""
Methylation Module

Integrates DNA methylation data (450K/EPIC array) with the regulatory pipeline.

Provides:
- Probe annotation with genomic context (genes, cCREs, TADs, ATAC peaks)
- Per-sample beta value processing
- Aggregation to gene/lncRNA promoters and cCRE regions
- Cohort-level matrix construction

Main entry points:
    - build_probe_reference_table(): Build annotated probe reference
    - process_sample_methylation(): Process single sample
    - build_cohort_matrices(): Build gene/cCRE methylation matrices
"""

from .probe_loader import (
    load_probe_reference,
    annotate_probes_with_promoters,
    annotate_probes_with_gene_bodies,
    annotate_probes_with_ccres,
    annotate_probes_with_lncrnas,
)

from .sample_processing import (
    load_sample_beta,
    process_sample_methylation,
    validate_sample_beta,
)

from .aggregation import (
    aggregate_to_genes,
    aggregate_to_lncrnas,
    aggregate_to_ccres,
    compute_promoter_methylation_score,
    aggregate_to_atac,
)

from .methylation_table import (
    build_probe_reference_table,
    build_sample_methylation_tables,
    build_cohort_matrices,
    run_methylation_pipeline,
)

from .meth_schemas import (
    empty_probe_annotation,
    empty_gene_methylation_entry,
    empty_lncrna_methylation_entry,
    empty_ccre_methylation_entry,
    empty_sample_probe_entry,
    METHYLATION_COLUMNS,
)

__version__ = "0.1.0"

__all__ = [
    # Probe loading
    "load_probe_reference",
    "annotate_probes_with_promoters",
    "annotate_probes_with_gene_bodies",
    "annotate_probes_with_ccres",
    "annotate_probes_with_lncrnas",
    # Sample processing
    "load_sample_beta",
    "process_sample_methylation",
    "validate_sample_beta",
    # Aggregation
    "aggregate_to_genes",
    "aggregate_to_lncrnas",
    "aggregate_to_ccres",
    "aggregate_to_atac",
    "compute_promoter_methylation_score",
    # Main orchestrators
    "build_probe_reference_table",
    "build_sample_methylation_tables",
    "build_cohort_matrices",
    "run_methylation_pipeline",
    # Schemas
    "empty_probe_annotation",
    "empty_gene_methylation_entry",
    "empty_lncrna_methylation_entry",
    "empty_ccre_methylation_entry",
    "empty_sample_probe_entry",
    "METHYLATION_COLUMNS",
]
