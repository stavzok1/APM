"""
Schema definitions for methylation data structures.

Consistent with the pipeline's schema pattern - provides empty/default
structures and validation helpers for methylation-specific data.
"""

from typing import Dict, Any, List, Optional


# =============================================================================
# COLUMN NAME CONSTANTS
# =============================================================================

class METHYLATION_COLUMNS:
    """Standard column names used across methylation module."""
    
    # Probe reference columns
    PROBE_ID = "probeID"
    CHROM = "chrom"
    START = "start"
    END = "end"
    STRAND = "strand"
    
    # GDC reference original columns
    GENES_UNIQ = "genesUniq"
    GENE_NAMES = "geneNames"
    TRANSCRIPT_TYPES = "transcriptTypes"
    TRANSCRIPT_IDS = "transcriptIDs"
    DIST_TO_TSS = "distToTSS"
    CGI = "CGI"
    CGI_POSITION = "CGIposition"
    
    # Derived annotation columns
    IN_CGI = "in_CGI"
    CGI_CONTEXT = "CGI_context"
    GENE_LIST = "gene_list"
    
    # Promoter annotations
    IN_PROMOTER = "in_promoter"
    PROMOTER_GENES = "promoter_genes"
    PROMOTER_GENE_IDS = "promoter_gene_ids"
    
    # Gene body annotations
    IN_GENE_BODY = "in_gene_body"
    GENE_BODY_GENES = "gene_body_genes"
    
    # cCRE annotations
    OVERLAPPING_CCRES = "overlapping_ccres"
    CCRE_TYPES = "ccre_types"
    N_OVERLAPPING_CCRES = "n_overlapping_ccres"
    
    # lncRNA annotations
    IN_LNCRNA_PROMOTER = "in_lncrna_promoter"
    LNCRNA_PROMOTER_GENES = "lncrna_promoter_genes"
    
    # ATAC annotations
    OVERLAPPING_ATAC_PEAKS = "overlapping_atac_peaks"
    N_OVERLAPPING_ATAC = "n_overlapping_atac"
    
    # TAD annotations (dict column)
    TAD_DOMAINS = "TAD_domains"
    
    # Sample-level columns
    BETA = "beta"
    M_VALUE = "m_value"
    
    # Aggregation output columns - Gene level
    PROMOTER_BETA_MEAN = "promoter_beta_mean"
    PROMOTER_BETA_MEDIAN = "promoter_beta_median"
    PROMOTER_BETA_STD = "promoter_beta_std"
    PROMOTER_BETA_MIN = "promoter_beta_min"
    PROMOTER_BETA_MAX = "promoter_beta_max"
    PROMOTER_BETA_RANGE = "promoter_beta_range"
    PROMOTER_N_PROBES = "promoter_n_probes"
    PROMOTER_FRAC_HYPERMETH = "promoter_frac_hypermeth"
    PROMOTER_FRAC_HYPOMETH = "promoter_frac_hypometh"
    PROMOTER_CGI_BETA_MEAN = "promoter_CGI_beta_mean"
    PROMOTER_SHORE_BETA_MEAN = "promoter_shore_beta_mean"
    GENE_BODY_BETA_MEAN = "gene_body_beta_mean"
    GENE_BODY_N_PROBES = "gene_body_n_probes"
    
    # Aggregation output columns - cCRE level
    CCRE_BETA_MEAN = "ccre_beta_mean"
    CCRE_BETA_MEDIAN = "ccre_beta_median"
    CCRE_BETA_STD = "ccre_beta_std"
    CCRE_N_PROBES = "ccre_n_probes"
    CCRE_FRAC_HYPERMETH = "ccre_frac_hypermeth"
    CCRE_FRAC_HYPOMETH = "ccre_frac_hypometh"
    CCRE_CGI_OVERLAP = "ccre_CGI_overlap"
    CCRE_CGI_BETA_MEAN = "ccre_CGI_beta_mean"


# =============================================================================
# CGI CONTEXT VALUES
# =============================================================================

class CGI_CONTEXTS:
    """Valid CpG island context values."""
    ISLAND = "Island"
    N_SHORE = "N_Shore"
    S_SHORE = "S_Shore"
    N_SHELF = "N_Shelf"
    S_SHELF = "S_Shelf"
    OPEN_SEA = "OpenSea"
    
    SHORES = {"N_Shore", "S_Shore"}
    SHELVES = {"N_Shelf", "S_Shelf"}
    
    @classmethod
    def is_shore(cls, context: str) -> bool:
        return context in cls.SHORES
    
    @classmethod
    def is_shelf(cls, context: str) -> bool:
        return context in cls.SHELVES
    
    @classmethod
    def is_island(cls, context: str) -> bool:
        return context == cls.ISLAND


# =============================================================================
# PROBE ANNOTATION SCHEMAS
# =============================================================================

def empty_probe_annotation() -> Dict[str, Any]:
    """
    Empty/default probe annotation structure.
    
    Used when creating new probe entries or filling gaps.
    """
    return {
        # Core coordinates
        "probeID": None,
        "chrom": None,
        "start": None,
        "end": None,
        "strand": None,
        
        # CGI context
        "in_CGI": False,
        "CGI_context": "OpenSea",
        "distToTSS": None,
        
        # Gene context
        "gene_list": [],
        "in_promoter": False,
        "promoter_genes": [],
        "promoter_gene_ids": [],
        "in_gene_body": False,
        "gene_body_genes": [],
        
        # cCRE context
        "overlapping_ccres": [],
        "ccre_types": [],
        "n_overlapping_ccres": 0,
        
        # lncRNA context
        "in_lncrna_promoter": False,
        "lncrna_promoter_genes": [],
        
        # ATAC context
        "overlapping_atac_peaks": [],
        "n_overlapping_atac": 0,
        
        # TAD context (per-biosample)
        "TAD_domains": {},
    }


def empty_probe_gene_link() -> Dict[str, Any]:
    """Empty structure for a probe-to-gene link."""
    return {
        "gene_name": None,
        "gene_id": None,
        "region": None,  # "promoter" | "gene_body" | "intergenic"
        "dist_to_tss": None,
        "in_CGI": False,
        "CGI_context": "OpenSea",
    }


def empty_probe_ccre_link() -> Dict[str, Any]:
    """Empty structure for a probe-to-cCRE link."""
    return {
        "cCRE_id": None,
        "ENCODE_id": None,
        "ccre_type": None,
        "overlap_bp": 0,
    }


# =============================================================================
# SAMPLE-LEVEL SCHEMAS
# =============================================================================

def empty_sample_probe_entry() -> Dict[str, Any]:
    """
    Empty sample-level probe entry.
    
    Per-sample probe data with beta value and key annotations.
    """
    return {
        "probeID": None,
        "beta": None,
        "m_value": None,
        
        # Key annotations carried from reference
        "chrom": None,
        "start": None,
        "in_CGI": False,
        "CGI_context": "OpenSea",
        "in_promoter": False,
        "promoter_genes": [],
        "overlapping_ccres": [],
    }


# =============================================================================
# GENE-LEVEL METHYLATION SCHEMAS
# =============================================================================

def empty_gene_methylation_entry() -> Dict[str, Any]:
    """
    Empty gene-level methylation aggregation.
    
    Includes promoter and gene body metrics.
    """
    return {
        "gene_name": None,
        "gene_id": None,
        
        # Promoter metrics
        "promoter_beta_mean": None,
        "promoter_beta_median": None,
        "promoter_beta_std": None,
        "promoter_beta_min": None,
        "promoter_beta_max": None,
        "promoter_beta_range": None,
        "promoter_n_probes": 0,
        "promoter_frac_hypermeth": None,
        "promoter_frac_hypometh": None,
        "promoter_CGI_beta_mean": None,
        "promoter_shore_beta_mean": None,
        
        # Gene body metrics
        "gene_body_beta_mean": None,
        "gene_body_n_probes": 0,
        
        # Probe lists (for traceability)
        "promoter_probes": [],
        "gene_body_probes": [],
    }


def empty_lncrna_methylation_entry() -> Dict[str, Any]:
    """
    Empty lncRNA-level methylation aggregation.
    
    Same structure as gene but for lncRNAs.
    """
    return {
        "lncrna_name": None,
        "lncrna_id": None,
        
        # Promoter metrics (lncRNAs have promoters too)
        "promoter_beta_mean": None,
        "promoter_beta_median": None,
        "promoter_beta_std": None,
        "promoter_n_probes": 0,
        "promoter_frac_hypermeth": None,
        "promoter_frac_hypometh": None,
        "promoter_CGI_beta_mean": None,
        
        # Body metrics
        "body_beta_mean": None,
        "body_n_probes": 0,
        
        # Probe lists
        "promoter_probes": [],
        "body_probes": [],
    }


# =============================================================================
# CCRE-LEVEL METHYLATION SCHEMAS
# =============================================================================

def empty_ccre_methylation_entry() -> Dict[str, Any]:
    """
    Empty cCRE-level methylation aggregation.
    """
    return {
        "cCRE_id": None,
        "ENCODE_id": None,
        "ccre_type": None,
        
        # Methylation metrics
        "ccre_beta_mean": None,
        "ccre_beta_median": None,
        "ccre_beta_std": None,
        "ccre_n_probes": 0,
        "ccre_frac_hypermeth": None,
        "ccre_frac_hypometh": None,
        
        # CGI context within cCRE
        "ccre_CGI_overlap": False,
        "ccre_CGI_beta_mean": None,
        "ccre_n_CGI_probes": 0,
        
        # Probe list
        "probes": [],
    }


# =============================================================================
# COHORT-LEVEL SCHEMAS
# =============================================================================

def empty_cohort_gene_methylation() -> Dict[str, Any]:
    """
    Empty cohort-level gene methylation summary.
    
    Aggregated across samples.
    """
    return {
        "gene_name": None,
        "n_samples_with_data": 0,
        "mean_promoter_beta": None,
        "std_promoter_beta": None,
        "frac_samples_hypermeth": None,
        "frac_samples_hypometh": None,
    }


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

def validate_probe_annotation(probe: Dict[str, Any]) -> List[str]:
    """
    Validate a probe annotation dict and return list of issues.
    Returns empty list if valid.
    """
    issues = []
    
    if not isinstance(probe, dict):
        return ["probe is not a dict"]
    
    # Required fields
    required = ["probeID", "chrom", "start"]
    for field in required:
        if field not in probe or probe[field] is None:
            issues.append(f"missing required field: {field}")
    
    # Type checks
    if "promoter_genes" in probe and not isinstance(probe["promoter_genes"], list):
        issues.append("promoter_genes must be a list")
    
    if "overlapping_ccres" in probe and not isinstance(probe["overlapping_ccres"], list):
        issues.append("overlapping_ccres must be a list")
    
    if "TAD_domains" in probe and not isinstance(probe["TAD_domains"], dict):
        issues.append("TAD_domains must be a dict")
    
    return issues


def validate_gene_methylation(entry: Dict[str, Any]) -> List[str]:
    """
    Validate a gene methylation entry.
    """
    issues = []
    
    if not isinstance(entry, dict):
        return ["entry is not a dict"]
    
    if "gene_name" not in entry or entry["gene_name"] is None:
        issues.append("missing gene_name")
    
    # Check beta values in valid range
    beta_fields = [
        "promoter_beta_mean", "promoter_beta_median",
        "gene_body_beta_mean"
    ]
    for field in beta_fields:
        if field in entry and entry[field] is not None:
            if not (0 <= entry[field] <= 1):
                issues.append(f"{field} out of range [0, 1]: {entry[field]}")
    
    # Check fraction fields
    frac_fields = ["promoter_frac_hypermeth", "promoter_frac_hypometh"]
    for field in frac_fields:
        if field in entry and entry[field] is not None:
            if not (0 <= entry[field] <= 1):
                issues.append(f"{field} out of range [0, 1]: {entry[field]}")
    
    return issues


def validate_sample_beta_df(df) -> List[str]:
    """
    Validate a sample beta DataFrame.
    """
    import pandas as pd
    
    issues = []
    
    if not isinstance(df, pd.DataFrame):
        return ["not a DataFrame"]
    
    if "probeID" not in df.columns:
        issues.append("missing probeID column")
    
    if "beta" not in df.columns:
        issues.append("missing beta column")
    
    if "beta" in df.columns:
        beta = df["beta"].dropna()
        if len(beta) > 0:
            if beta.min() < 0 or beta.max() > 1:
                issues.append(f"beta values out of range: min={beta.min()}, max={beta.max()}")
    
    return issues


# =============================================================================
# SCHEMA COERCION HELPERS
# =============================================================================

def ensure_probe_annotation(val: Any) -> Dict[str, Any]:
    """
    Coerce any value into a valid probe annotation dict.
    """
    if not isinstance(val, dict):
        return empty_probe_annotation()
    
    template = empty_probe_annotation()
    for key, default in template.items():
        if key not in val:
            val[key] = default
        elif isinstance(default, list) and not isinstance(val[key], list):
            val[key] = []
        elif isinstance(default, dict) and not isinstance(val[key], dict):
            val[key] = {}
    
    return val


def ensure_gene_methylation(val: Any) -> Dict[str, Any]:
    """
    Coerce any value into a valid gene methylation entry.
    """
    if not isinstance(val, dict):
        return empty_gene_methylation_entry()
    
    template = empty_gene_methylation_entry()
    for key, default in template.items():
        if key not in val:
            val[key] = default
    
    return val


def ensure_ccre_methylation(val: Any) -> Dict[str, Any]:
    """
    Coerce any value into a valid cCRE methylation entry.
    """
    if not isinstance(val, dict):
        return empty_ccre_methylation_entry()
    
    template = empty_ccre_methylation_entry()
    for key, default in template.items():
        if key not in val:
            val[key] = default
    
    return val
