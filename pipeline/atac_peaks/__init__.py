"""
ATAC peaks processing module.

Provides functionality for:
- Loading and normalizing TCGA ATAC peaks
- Distance matching to genes (TSS + gene body overlap)
- Overlap/proximity matching to cCREs
- TAD domain and boundary annotations

Main entry points:
    - load_atac_peaks(): Load and normalize peaks from CSV
    - build_atac_peak_table(): Complete pipeline for peak annotation
    - match_peaks_to_genes(): Gene distance/overlap matching
    - match_peaks_to_ccres(): cCRE overlap/proximity matching
"""

from .peak_loader import (
    load_atac_peaks,
    generate_peak_id,
    create_peak_id_mapping,
    filter_peaks_by_chromosomes,
    filter_peaks_by_score,
    filter_peaks_by_length,
)

from .gene_matching import (
    match_peaks_to_genes,
    build_gene_links,
    aggregate_genes_per_peak,
    get_peaks_for_gene,
)

from .ccre_matching import (
    match_peaks_to_ccres,
    build_ccre_links,
    aggregate_ccres_per_peak,
    get_ccres_for_peak,
    get_peaks_for_ccre,
)

from .tad_annotation import (
    annotate_peaks_with_tads,
    annotate_peaks_with_boundary_overlaps,
    annotate_peaks_with_all_tad_sources,
    get_peaks_at_boundaries,
    summarize_boundary_overlaps,
)

from .peak_table import (
    ATAC_SLIM_MEMORY_COLUMNS,
    build_atac_peak_table,
    load_atac_peaks_annotated,
    save_atac_outputs,
    get_peaks_near_gene,
    get_peaks_overlapping_gene,
    get_peaks_at_promoters,
)

from .annotate_df_with_peaks import (
    annotate_df_with_peaks,
    annotate_genes_with_peaks,
    annotate_ccres_with_peaks,
    annotate_svs_with_peaks,
    annotate_snvs_with_peaks,
    get_features_with_overlapping_peaks,
    get_features_with_upstream_peaks,
    get_features_with_downstream_peaks,
    get_closest_peak,
    summarize_peak_links,
    flatten_peak_links,
)

__all__ = [
    # Loading
    "load_atac_peaks",
    "generate_peak_id",
    "create_peak_id_mapping",
    "filter_peaks_by_chromosomes",
    "filter_peaks_by_score",
    "filter_peaks_by_length",
    # Gene matching
    "match_peaks_to_genes",
    "build_gene_links",
    "aggregate_genes_per_peak",
    "get_peaks_for_gene",
    # cCRE matching
    "match_peaks_to_ccres",
    "build_ccre_links",
    "aggregate_ccres_per_peak",
    "get_ccres_for_peak",
    "get_peaks_for_ccre",
    # TAD annotation
    "annotate_peaks_with_tads",
    "annotate_peaks_with_boundary_overlaps",
    "annotate_peaks_with_all_tad_sources",
    "get_peaks_at_boundaries",
    "summarize_boundary_overlaps",
    # Table building
    "ATAC_SLIM_MEMORY_COLUMNS",
    "build_atac_peak_table",
    "load_atac_peaks_annotated",
    "save_atac_outputs",
    "get_peaks_near_gene",
    "get_peaks_overlapping_gene",
    "get_peaks_at_promoters",
      # Feature annotation with peaks
    "annotate_df_with_peaks",
    "annotate_genes_with_peaks",
    "annotate_ccres_with_peaks",
    "annotate_svs_with_peaks",
    "annotate_snvs_with_peaks",
    "get_features_with_overlapping_peaks",
    "get_features_with_upstream_peaks",
    "get_features_with_downstream_peaks",
    "get_closest_peak",
    "summarize_peak_links",
    "flatten_peak_links",
]
