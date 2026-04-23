# """
# Regulatory Element Pipeline

# A modular pipeline for building regulatory element tables with
# multi-evidence gene links for immune visibility gene analysis.

# Main entry points:
#     - run_full_pipeline(): Complete pipeline
#     - run_genes_only(): Gene/lncRNA processing only
#     - run_distance_matching_only(): cCRE-gene matching only
#     - run_evidence_only(): Evidence collection only

# Modules:
#     - config: Configuration and parameters
#     - schemas: Data structure definitions
#     - utils: Shared utility functions
#     - genes: Gene and lncRNA processing
#     - regulatory_elements: cCRE loading and distance matching
#     - evidence: Evidence source processing (SCREEN, ABC, HiChIP)
# """

# from .main import (
#     run_full_pipeline,
#     run_genes_only,
#     run_distance_matching_only,
#     run_evidence_only,
# )

# from .config import (
#     PATHS,
#     BIOSAMPLES,
#     THRESHOLDS,
#     PRIMARY_GENES,
# )

# __version__ = "0.1.0"

# __all__ = [
#     # Main functions
#     "run_full_pipeline",
#     "run_genes_only",
#     "run_distance_matching_only",
#     "run_evidence_only",
#     # Config
#     "PATHS",
#     "BIOSAMPLES",
#     "THRESHOLDS",
#     "PRIMARY_GENES",
# ]




"""
Regulatory Element Pipeline
"""

__version__ = "0.1.0"
