"""
RPPA Integration Module for Immune Visibility Analysis.

This module provides comprehensive RPPA (Reverse Phase Protein Array)
data processing and analysis for the immune visibility pipeline.

Main entry points:
    - load_rppa_dataset(): Complete data loading pipeline
    - build_expression_matrix(): Construct sample × target matrix
    - compute_all_panel_scores(): Calculate pathway activation scores
    - detect_signaling_blocks(): Identify pathway disruptions
    - compute_protein_rna_discordance(): Compare protein vs RNA levels

Data structures:
    - RPPAMarkerPanels: Curated marker panel definitions
    - RPPAPathConfig: File path configuration
    - RPPAThresholdConfig: Analysis thresholds

Key concepts:
    
    1. VALIDATION STATUS: Antibodies have validation statuses:
       - "Valid": Fully validated, include by default
       - "Caution": Use with care, include by default with flag
       - "Under Evaluation": Exclude by default
    
    2. PHOSPHO PAIRS: Many proteins have total + phospho-specific
       antibodies. Activation ratio = phospho / total in log space.
    
    3. PANEL SCORES: Composite scores from biologically-related
       markers (e.g., DDR_activation, IFN_signaling, checkpoint).
    
    4. SIGNALING BLOCKS: Patterns suggesting pathway disruption
       (e.g., high DDR but low STING = sensing failure).

Example usage:
    
    >>> from rppa import load_rppa_dataset, compute_all_panel_scores
    >>> from rppa import get_default_panels, detect_signaling_blocks
    >>> 
    >>> # Load data
    >>> data = load_rppa_dataset(
    ...     sample_dir="path/to/samples",
    ...     annotation_path="path/to/annotation.csv",
    ... )
    >>> 
    >>> # Compute panel scores
    >>> panels = get_default_panels()
    >>> scores = compute_all_panel_scores(
    ...     data["expression_matrix"],
    ...     panels,
    ... )
    >>> 
    >>> # Detect signaling blocks
    >>> blocks = detect_signaling_blocks(data["expression_matrix"])
"""

__version__ = "0.1.0"

# Loader functions
from .rppa_loader import (
    load_rppa_dataset,
    load_rppa_annotation,
    load_single_sample_rppa,
    load_all_sample_rppa,
    build_expression_matrix,
    build_target_mappings,
    identify_phospho_pairs,
    integrate_sample_metadata,
    ValidationStatus,
)

# Panel definitions
from .rppa_panels import (
    RPPAMarkerPanels,
    get_default_panels,
    compute_zscore_panel_score,
    compute_pca_panel_score,
    compute_activation_ratio,
    compute_all_panel_scores,
    compute_all_activation_ratios,
    summarize_panel_coverage,
)

# Analysis functions
from .rppa_analysis import (
    detect_signaling_blocks,
    classify_ddr_ifn_quadrant,
    classify_immune_visibility_state,
    compute_protein_rna_discordance,
    identify_post_transcriptional_regulation,
    compute_immune_visibility_score,
    compute_antigen_presentation_capacity,
    test_score_association,
    correlate_with_outcomes,
)

# Configuration
from .rppa_config import (
    ImmuneVisibilityMarkers,
    RPPA_PATHS,
    RPPA_THRESHOLDS,
    IMMUNE_VISIBILITY_MARKERS,
    RPPA_TARGET_TO_GENE_OVERRIDES,
    KEY_PHOSPHO_PAIRS,
    get_rppa_config,
)

# Schemas
from .rppa_schemas import (
    empty_rppa_sample_entry,
    empty_rppa_target_entry,
    empty_rppa_gene_evidence,
    empty_rppa_panel_scores,
    empty_signaling_blocks,
    validate_rppa_expression_matrix,
    validate_rppa_annotation,
    validate_phospho_pairs,
    build_rppa_gene_evidence,
    ensure_rppa_panel_scores,
    ensure_signaling_blocks,
)

__all__ = [
    # Version
    "__version__",
    
    # Loader
    "load_rppa_dataset",
    "load_rppa_annotation",
    "load_single_sample_rppa",
    "load_all_sample_rppa",
    "build_expression_matrix",
    "build_target_mappings",
    "identify_phospho_pairs",
    "integrate_sample_metadata",
    "ValidationStatus",
    
    # Panels
    "RPPAMarkerPanels",
    "get_default_panels",
    "compute_zscore_panel_score",
    "compute_pca_panel_score",
    "compute_activation_ratio",
    "compute_all_panel_scores",
    "compute_all_activation_ratios",
    "summarize_panel_coverage",
    
    # Analysis
    "detect_signaling_blocks",
    "classify_ddr_ifn_quadrant",
    "classify_immune_visibility_state",
    "compute_protein_rna_discordance",
    "identify_post_transcriptional_regulation",
    "compute_immune_visibility_score",
    "compute_antigen_presentation_capacity",
    "test_score_association",
    "correlate_with_outcomes",
    
    # Config
    "RPPAPathConfig",
    "RPPAThresholdConfig",
    "ImmuneVisibilityMarkers",
    "RPPA_PATHS",
    "RPPA_THRESHOLDS",
    "IMMUNE_VISIBILITY_MARKERS",
    "RPPA_TARGET_TO_GENE_OVERRIDES",
    "KEY_PHOSPHO_PAIRS",
    "get_rppa_config",
    
    # Schemas
    "empty_rppa_sample_entry",
    "empty_rppa_target_entry",
    "empty_rppa_gene_evidence",
    "empty_rppa_panel_scores",
    "empty_signaling_blocks",
    "validate_rppa_expression_matrix",
    "validate_rppa_annotation",
    "validate_phospho_pairs",
    "build_rppa_gene_evidence",
    "ensure_rppa_panel_scores",
    "ensure_signaling_blocks",
]
