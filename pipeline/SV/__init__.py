"""
SV (Structural Variant) Pipeline Module

A modular pipeline for processing structural variants with:
- VEP annotation
- Gene/element/lncRNA spatial mapping
- FIMO motif scanning
- BND neojunction motif annotation
- Evidence integration

Main entry points:
    - run_sv_pipeline(): Complete pipeline
    - run_sv_vcf_processing(): VCF processing only
    - run_sv_motif_scanning(): Motif scanning only
    - run_sv_evidence_merge(): Evidence merging only

Submodules:
    - vcf_loader: Load and parse Manta VCF files
    - sv_filtering: Apply strict/lenient filtering
    - spatial_mapping: Map SVs to genes, elements, lncRNAs
    - vep_annotation: Parse and integrate VEP annotations
    - bed_intervals: Generate BED intervals for FASTA extraction
    - motif_scanning: FIMO integration and motif hit processing
    - evidence_merger: Combine all evidence into final tables
"""

from .vcf_loader import (
    load_manta_sv_vcf,
    get_bnd_remote_coords,
    breakend_to_vcf_alt,
)

from .sv_filtering import (
    get_strict_sv_set,
    get_lenient_sv_set,
)

from .spatial_mapping import (
    map_svs_to_genes,
    map_svs_to_elements,
    map_svs_to_lncrnas,
    compute_signed_distance_with_overlap,
    classify_span_sv_gene_hit,
    classify_span_sv_elem_hit,
)

from .vep_annotation import (
    add_vep_hits_columns,
    parse_vep_csq,
)

from .bed_intervals import (
    build_sv_flanks_and_overlaps_bed,
    create_bed_from_sv_csv,
)

from .motif_scanning import (
    extract_selected_motifs,
    recombine_sv_fimo,
    parse_fimo_tsv_to_bed,
    load_all_ccre_fimo,
    annotate_bnd_neojunction_motifs,
    annotate_neojunction_motifs_from_directory,
    summarize_neojunction_motifs,
)

from .pipeline import (
    run_sv_pipeline,
    run_sv_vcf_processing,
    run_sv_motif_scanning,
    process_single_vcf,
)

__version__ = "0.1.0"

__all__ = [
    # VCF loading
    "load_manta_sv_vcf",
    "get_bnd_remote_coords",
    "breakend_to_vcf_alt",
    # Filtering
    "get_strict_sv_set",
    "get_lenient_sv_set",
    # Spatial mapping
    "map_svs_to_genes",
    "map_svs_to_elements",
    "map_svs_to_lncrnas",
    "compute_signed_distance_with_overlap",
    "classify_span_sv_gene_hit",
    "classify_span_sv_elem_hit",
    # VEP
    "add_vep_hits_columns",
    "parse_vep_csq",
    # BED intervals
    "build_sv_flanks_and_overlaps_bed",
    "create_bed_from_sv_csv",
    # Motif scanning
    "extract_selected_motifs",
    "recombine_sv_fimo",
    "parse_fimo_tsv_to_bed",
    # Neojunction motif annotation
    "load_all_ccre_fimo",
    "annotate_bnd_neojunction_motifs",
    "annotate_neojunction_motifs_from_directory",
    "summarize_neojunction_motifs",
    # Pipeline
    "run_sv_pipeline",
    "run_sv_vcf_processing",
    "run_sv_motif_scanning",
    "process_single_vcf",
]