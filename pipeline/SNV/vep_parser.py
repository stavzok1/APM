"""
VEP (Variant Effect Predictor) annotation parsing.

Extracts structured consequence annotations from VEP CSQ field,
organized by feature type (Transcript, RegulatoryFeature, MotifFeature).
"""

from typing import Any, Callable, Dict, List, Mapping, Optional, Set

import numpy as np
import pandas as pd


# =============================================================================
# VEP FIELD DEFINITIONS
# =============================================================================

VEP_GENE_FIELDS: List[str] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "VARIANT_CLASS",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "GENE_PHENO",
    "SIFT",
    "PolyPhen",
    "DOMAINS",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    # Optional / plugin-dependent CSQ fields (parsed when present in the VCF CSQ header).
    "CADD_PHRED",
    "REVEL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "MUTATION_MD",
    "LoF",
    "LoF_filter",
    "LoF_info",
]

VEP_REGULATORY_FIELDS: List[str] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "Existing_variation",
    "DISTANCE",
    "VARIANT_CLASS",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
]

VEP_MOTIF_FIELDS: List[str] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "Existing_variation",
    "VARIANT_CLASS",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
]

# Consequence term sets for classification
MISSENSE_TERMS: Set[str] = {"missense_variant"}
NONSENSE_TERMS: Set[str] = {"stop_gained"}
FRAMESHIFT_TERMS: Set[str] = {"frameshift_variant"}
SPLICE_TERMS: Set[str] = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
}

# Default CSQ format string (VEP v115)
DEFAULT_CSQ_FORMAT: str = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|"
    "CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|"
    "TREMBL|UNIPARC|ARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|"
    "HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|"
    "gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_MID_AF|"
    "gnomADe_NFE_AF|gnomADe_REMAINING_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|"
    "gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|"
    "gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|MAX_AF|"
    "MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|"
    "MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
)


# =============================================================================
# PARSING HELPERS
# =============================================================================

def parse_csq_description(csq_description: str) -> List[str]:
    """
    Parse VEP CSQ field format description to extract column names.
    
    Args:
        csq_description: VEP CSQ header description string
        
    Returns:
        List of field names in order
    """
    if "Format:" in csq_description:
        fmt = csq_description.split("Format:")[1].strip().strip('"').strip()
    else:
        fmt = csq_description.strip().strip('"').strip()
    
    return fmt.split("|")


def _get_field(
    parts: List[str],
    csq_index: Dict[str, int],
    name: str,
    default: Optional[str] = None,
) -> Optional[str]:
    """
    Safely extract a field from parsed CSQ parts.
    
    Handles missing columns and truncated entries gracefully.
    """
    idx = csq_index.get(name)
    if idx is None:
        return default
    if idx >= len(parts):
        return default
    val = parts[idx]
    return default if val == "" else val


def _has_consequence(
    gene_hits: List[Dict[str, Any]],
    target_terms: Set[str],
) -> bool:
    """
    Check if ANY transcript has a consequence in target_terms.
    
    VEP uses '&' to separate multiple SO terms in one 'Consequence' field.
    """
    for hit in gene_hits:
        cons = hit.get("Consequence") or ""
        for term in cons.split("&"):
            if term.strip() in target_terms:
                return True
    return False


def _has_consequence_filtered(
    gene_hits: List[Dict[str, Any]],
    target_terms: Set[str],
    filter_field: str,
    filter_predicate: Callable[[Any], bool],
) -> bool:
    """
    Check if ANY transcript satisfying a filter has a consequence in target_terms.
    
    Args:
        gene_hits: List of transcript hit dictionaries
        target_terms: Set of SO consequence terms to check
        filter_field: Field to filter on (e.g., "CANONICAL", "MANE_SELECT")
        filter_predicate: Function returning True if the field value passes filter
        
    Examples:
        # Canonical transcripts only
        _has_consequence_filtered(hits, MISSENSE_TERMS, "CANONICAL", lambda v: v == "1")
        
        # MANE SELECT transcripts only
        _has_consequence_filtered(hits, MISSENSE_TERMS, "MANE_SELECT", lambda v: v is not None)
    """
    for hit in gene_hits:
        if not filter_predicate(hit.get(filter_field)):
            continue
        cons = hit.get("Consequence") or ""
        for term in cons.split("&"):
            if term.strip() in target_terms:
                return True
    return False


# =============================================================================
# MAIN VEP PARSING
# =============================================================================

def add_vep_hits_columns(
    df: pd.DataFrame,
    csq_description: str,
    primary_genes: List[str],
    csq_column: str = "CSQ",
    gene_symbol_mapping: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Parse VEP CSQ annotations and add structured hit columns.
    
    For each row (variant + ALT allele), adds:
    
    Hit lists:
        - gene_hits: List[Dict] - Transcript-level annotations for primary genes
        - regulatory_hits: List[Dict] - RegulatoryFeature annotations
        - motif_hits: List[Dict] - MotifFeature annotations
    
    Summary columns:
        - gene_symbols: Comma-separated unique gene symbols hit
        - hits_canonical: Whether any canonical transcript is hit
    
    Impact flags (any transcript):
        - has_missense, has_nonsense, has_frameshift, has_splice_effect
    
    Impact flags (canonical transcripts only, CANONICAL == "1"):
        - has_missense_canonical, has_nonsense_canonical, 
          has_frameshift_canonical, has_splice_effect_canonical
    
    Impact flags (MANE SELECT transcripts only):
        - has_missense_mane, has_nonsense_mane,
          has_frameshift_mane, has_splice_effect_mane
    
    Args:
        df: DataFrame with one row per (variant, ALT) and CSQ column
        csq_description: VEP CSQ format description string
        primary_genes: List of gene symbols to filter for (others ignored)
        csq_column: Name of CSQ column in df
        gene_symbol_mapping: Optional old_symbol -> canonical (HGNC/UCSC/legacy).
            When set, VEP ``SYMBOL`` values are resolved to panel symbols before filtering.
    
    Returns:
        DataFrame with additional VEP annotation columns
    """
    from ..genes.symbol_normalization import resolve_symbol_to_panel

    # Parse CSQ column names
    csq_cols = parse_csq_description(csq_description)
    csq_index = {name: i for i, name in enumerate(csq_cols)}
    
    primary_gene_set = set(primary_genes)
    
    def parse_csq_for_row(row: pd.Series) -> pd.Series:
        """Parse CSQ for a single row."""
        alt_allele = str(row.get("alt", ""))
        csq_val = row.get(csq_column)
        
        # Handle missing CSQ
        if csq_val is None or (isinstance(csq_val, float) and np.isnan(csq_val)):
            return _empty_vep_result()
        
        # CSQ might be a list or a comma-separated string
        entries = csq_val if isinstance(csq_val, list) else str(csq_val).split(",")
        
        gene_hits = []
        regulatory_hits = []
        motif_hits = []
        
        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue
            
            parts = entry.split("|")
            
            # Check allele matches this ALT
            allele = _get_field(parts, csq_index, "Allele")
            if allele is None or allele != alt_allele:
                continue
            
            feature_type = _get_field(parts, csq_index, "Feature_type")
            if feature_type is None:
                continue
            
            if feature_type == "Transcript":
                gene_name = _get_field(parts, csq_index, "SYMBOL")
                if gene_symbol_mapping is not None:
                    canonical = resolve_symbol_to_panel(
                        gene_name, primary_gene_set, gene_symbol_mapping
                    )
                else:
                    canonical = gene_name if gene_name in primary_gene_set else None
                if canonical is None:
                    continue

                hit = {fname: _get_field(parts, csq_index, fname) for fname in VEP_GENE_FIELDS}
                hit["SYMBOL"] = canonical
                gene_hits.append(hit)
                
            elif feature_type == "RegulatoryFeature":
                hit = {fname: _get_field(parts, csq_index, fname) for fname in VEP_REGULATORY_FIELDS}
                regulatory_hits.append(hit)
                
            elif feature_type == "MotifFeature":
                hit = {fname: _get_field(parts, csq_index, fname) for fname in VEP_MOTIF_FIELDS}
                motif_hits.append(hit)
        
        # Build summary fields
        gene_symbols = sorted({h.get("SYMBOL") for h in gene_hits if h.get("SYMBOL")})
        gene_symbols_str = ",".join(gene_symbols) if gene_symbols else ""
        
        hits_canonical = any(h.get("CANONICAL") == "1" for h in gene_hits)
        
        # Global consequence flags
        has_missense = _has_consequence(gene_hits, MISSENSE_TERMS)
        has_nonsense = _has_consequence(gene_hits, NONSENSE_TERMS)
        has_frameshift = _has_consequence(gene_hits, FRAMESHIFT_TERMS)
        has_splice = _has_consequence(gene_hits, SPLICE_TERMS)
        
        # Canonical-specific flags
        canonical_pred = lambda v: v == "1"
        has_missense_canonical = _has_consequence_filtered(gene_hits, MISSENSE_TERMS, "CANONICAL", canonical_pred)
        has_nonsense_canonical = _has_consequence_filtered(gene_hits, NONSENSE_TERMS, "CANONICAL", canonical_pred)
        has_frameshift_canonical = _has_consequence_filtered(gene_hits, FRAMESHIFT_TERMS, "CANONICAL", canonical_pred)
        has_splice_canonical = _has_consequence_filtered(gene_hits, SPLICE_TERMS, "CANONICAL", canonical_pred)
        
        # MANE SELECT-specific flags
        mane_pred = lambda v: v is not None
        has_missense_mane = _has_consequence_filtered(gene_hits, MISSENSE_TERMS, "MANE_SELECT", mane_pred)
        has_nonsense_mane = _has_consequence_filtered(gene_hits, NONSENSE_TERMS, "MANE_SELECT", mane_pred)
        has_frameshift_mane = _has_consequence_filtered(gene_hits, FRAMESHIFT_TERMS, "MANE_SELECT", mane_pred)
        has_splice_mane = _has_consequence_filtered(gene_hits, SPLICE_TERMS, "MANE_SELECT", mane_pred)
        
        return pd.Series({
            "gene_hits": gene_hits,
            "regulatory_hits": regulatory_hits,
            "motif_hits": motif_hits,
            "gene_symbols": gene_symbols_str,
            "hits_canonical": hits_canonical,
            "has_missense": has_missense,
            "has_nonsense": has_nonsense,
            "has_frameshift": has_frameshift,
            "has_splice_effect": has_splice,
            "has_missense_canonical": has_missense_canonical,
            "has_nonsense_canonical": has_nonsense_canonical,
            "has_frameshift_canonical": has_frameshift_canonical,
            "has_splice_effect_canonical": has_splice_canonical,
            "has_missense_mane": has_missense_mane,
            "has_nonsense_mane": has_nonsense_mane,
            "has_frameshift_mane": has_frameshift_mane,
            "has_splice_effect_mane": has_splice_mane,
        })
    
    # Apply parsing to all rows
    hits_df = df.apply(parse_csq_for_row, axis=1)
    return pd.concat([df, hits_df], axis=1)


def _empty_vep_result() -> pd.Series:
    """Return empty VEP parsing result."""
    return pd.Series({
        "gene_hits": [],
        "regulatory_hits": [],
        "motif_hits": [],
        "gene_symbols": "",
        "hits_canonical": False,
        "has_missense": False,
        "has_nonsense": False,
        "has_frameshift": False,
        "has_splice_effect": False,
        "has_missense_canonical": False,
        "has_nonsense_canonical": False,
        "has_frameshift_canonical": False,
        "has_splice_effect_canonical": False,
        "has_missense_mane": False,
        "has_nonsense_mane": False,
        "has_frameshift_mane": False,
        "has_splice_effect_mane": False,
    })


# =============================================================================
# CONSEQUENCE IMPACT CLASSIFICATION
# =============================================================================

IMPACT_HIERARCHY: Dict[str, int] = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
}


def get_max_impact(gene_hits: List[Dict[str, Any]]) -> Optional[str]:
    """Get the highest impact level from gene hits."""
    if not gene_hits:
        return None
    
    max_level = 0
    max_impact = None
    
    for hit in gene_hits:
        impact = hit.get("IMPACT")
        level = IMPACT_HIERARCHY.get(impact, 0)
        if level > max_level:
            max_level = level
            max_impact = impact
    
    return max_impact


def classify_variant_type(
    gene_hits: List[Dict[str, Any]],
) -> str:
    """
    Classify variant by most severe consequence type.
    
    Returns one of: "nonsense", "frameshift", "splice", "missense", 
                    "synonymous", "intronic", "regulatory", "other"
    """
    if _has_consequence(gene_hits, NONSENSE_TERMS):
        return "nonsense"
    if _has_consequence(gene_hits, FRAMESHIFT_TERMS):
        return "frameshift"
    if _has_consequence(gene_hits, SPLICE_TERMS):
        return "splice"
    if _has_consequence(gene_hits, MISSENSE_TERMS):
        return "missense"
    if _has_consequence(gene_hits, {"synonymous_variant"}):
        return "synonymous"
    if _has_consequence(gene_hits, {"intron_variant"}):
        return "intronic"
    
    # Check for any regulatory consequences
    regulatory_terms = {
        "regulatory_region_variant",
        "TF_binding_site_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
    }
    if _has_consequence(gene_hits, regulatory_terms):
        return "regulatory"
    
    return "other"
