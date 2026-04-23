"""
SNV/Indel processing module.

Functions for:
- Loading Mutect2 VCF files (VEP-annotated)
- High-confidence somatic filtering
- VEP consequence annotation parsing
- cCRE overlap matching for regulatory impact
- Optional MEME FIMO on reference windows per SNV (``annotate_snvs_with_fimo`` / ``run_fimo`` on load)
- Optional ChIP overlap on the variant base (``snv_chip.annotate_snvs_with_chip`` / ``run_chip`` on load)
"""

from .vcf_loader import (
    load_mutect_snv_vcf,
    load_mutect_snv_batch,
    save_snv_outputs,
    load_snv_csv,
)

from .somatic_filter import (
    high_conf_somatic_mask,
    apply_somatic_filter,
)

from .vep_parser import (
    add_vep_hits_columns,
    parse_csq_description,
    VEP_GENE_FIELDS,
    VEP_REGULATORY_FIELDS,
    VEP_MOTIF_FIELDS,
)

from .ccre_matching import (
    match_snvs_to_ccres,
)

from .snv_chip import annotate_snvs_with_chip
from .snv_fimo import annotate_snvs_with_fimo

__all__ = [
    # Main loaders
    "load_mutect_snv_vcf",
    "load_mutect_snv_batch",
    "save_snv_outputs",
    "load_snv_csv",
    # Filtering
    "high_conf_somatic_mask",
    "apply_somatic_filter",
    # VEP parsing
    "add_vep_hits_columns",
    "parse_csq_description",
    "VEP_GENE_FIELDS",
    "VEP_REGULATORY_FIELDS",
    "VEP_MOTIF_FIELDS",
    # cCRE matching
    "match_snvs_to_ccres",
    # Reference-window FIMO (optional)
    "annotate_snvs_with_fimo",
    # ChIP overlap on variant base (optional)
    "annotate_snvs_with_chip",
]
