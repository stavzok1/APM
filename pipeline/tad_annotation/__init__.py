"""
TAD annotation module for the regulatory elements pipeline.

This module annotates features (genes, cCREs, lncRNAs, SVs) with their
TAD domain context and mirrors those annotations back into TAD domain tables.

Architecture:
    - relations.py: Pure geometry functions (overlap, primary selection, distances)
    - annotator.py: Feature → TAD annotation (adds TAD_domains column to features)
    - mirroring.py: TAD → feature annotation (adds gene_hits/cCRE_hits to domains)
    - tad_config.py: Biosample registry, PAM50 metadata, source discovery
    - loader.py: Multi-biosample annotation orchestration

Typical workflow (single biosample):
    >>> from pipeline.tad_annotation import annotate_genes_with_tads
    >>> genes_df = annotate_genes_with_tads(
    ...     genes_df, domains, flanks, bounds, biosample="Kim_T47D"
    ... )

Typical workflow (all biosamples):
    >>> from pipeline.tad_annotation import annotate_with_all_tad_sources
    >>> genes_df, ccre_df, lncrnas_df = annotate_with_all_tad_sources(
    ...     genes_df, ccre_df, lncrnas_df,
    ...     processed_dir=Path("/home/stavz/masters/gdc/TADs/processed"),
    ... )
"""

# Relations (pure geometry)
from .relations import (
    get_element_tad_relations,
    pick_primary_domain_id,
    boundary_overlap_and_dist,
    compute_normalized_position,
)

# Annotators (feature → TAD)
from .annotator import (
    annotate_df_with_tads,
    annotate_genes_with_tads,
    annotate_ccres_with_tads,
    annotate_lncrnas_with_tads,
    annotate_svs_with_tads,
    feature_intervals_from_row,
)

# Mirroring (TAD → feature)
from .mirroring import (
    mirror_hits_into_domains,
    mirror_genes_into_domains,
    mirror_lncrnas_into_domains,
    mirror_ccres_into_domains,
    mirror_all_features_into_domains,
    # Query helpers
    get_domain_gene_count,
    get_domains_containing_gene,
    get_genes_in_domain,
)

# Config and metadata
from .tad_config import (
    TADBiosample,
    TADSourcePaths,
    TAD_BIOSAMPLE_REGISTRY,
    PAM50_GROUPS,
    TISSUE_GROUPS,
    PAM50_REFERENCE_MAP,
    discover_tad_sources,
    get_biosamples_by_pam50,
    get_biosamples_by_tissue_type,
    get_biosamples_by_study,
    get_best_tad_source_for_pam50,
)

# Multi-biosample loading and annotation
from .loader import (
    load_tad_file,
    load_tad_source,
    annotate_with_all_tad_sources,
    annotate_svs_with_all_tad_sources,
    mirror_and_save_all_domains,
    get_pam50_matched_sources,
    annotate_with_pam50_matched,
)

from .boundary_enrichment import build_boundaries_enriched_table


__all__ = [
    # Relations
    "get_element_tad_relations",
    "pick_primary_domain_id",
    "boundary_overlap_and_dist",
    "compute_normalized_position",
    # Annotators
    "annotate_df_with_tads",
    "annotate_genes_with_tads",
    "annotate_ccres_with_tads",
    "annotate_lncrnas_with_tads",
    "annotate_svs_with_tads",
    "feature_intervals_from_row",
    # Mirroring
    "mirror_hits_into_domains",
    "mirror_genes_into_domains",
    "mirror_lncrnas_into_domains",
    "mirror_ccres_into_domains",
    "mirror_all_features_into_domains",
    # Query helpers
    "get_domain_gene_count",
    "get_domains_containing_gene",
    "get_genes_in_domain",
    # Config
    "TADBiosample",
    "TADSourcePaths",
    "TAD_BIOSAMPLE_REGISTRY",
    "PAM50_GROUPS",
    "TISSUE_GROUPS",
    "PAM50_REFERENCE_MAP",
    "discover_tad_sources",
    "get_biosamples_by_pam50",
    "get_biosamples_by_tissue_type",
    "get_biosamples_by_study",
    "get_best_tad_source_for_pam50",
    # Loader
    "load_tad_file",
    "load_tad_source",
    "annotate_with_all_tad_sources",
    "annotate_svs_with_all_tad_sources",
    "mirror_and_save_all_domains",
    "get_pam50_matched_sources",
    "annotate_with_pam50_matched",
    "build_boundaries_enriched_table",
]
