"""
CNV Pipeline Module
====================
Copy-number variant annotation pipeline for the APM regulatory analysis.

Sub-modules:
    loader, features, gene_hits, elem_hits, geometry, runner
"""

from .runner import process_cnv_directory, annotate_single_sample
from .loader import load_cnv_file, extract_sample_id_from_annotations
from .features import add_basic_cnv_features
from .gene_hits import (
    annotate_cnv_with_gene_hits,
    annotate_cnv_with_lncrna_hits,
    annotate_cnv_with_mirna_hits,
)
from .elem_hits import annotate_cnv_with_elem_hits
from .gene_level_ascat import (
    build_ascat_gene_sample_table,
    enrich_ascat_gene_level,
    is_ascat_gene_level_filename,
    load_ascat_gene_level_tsv,
)

__all__ = [
    "process_cnv_directory", "annotate_single_sample",
    "load_cnv_file", "extract_sample_id_from_annotations",
    "add_basic_cnv_features",
    "annotate_cnv_with_gene_hits", "annotate_cnv_with_lncrna_hits",
    "annotate_cnv_with_mirna_hits", "annotate_cnv_with_elem_hits",
    "build_ascat_gene_sample_table", "enrich_ascat_gene_level",
    "is_ascat_gene_level_filename", "load_ascat_gene_level_tsv",
]
