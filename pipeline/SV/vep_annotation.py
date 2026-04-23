"""
VEP annotation parsing and integration.

Functions for:
- Parsing VEP CSQ format
- Extracting gene, regulatory, and motif annotations
- Adding consequence flags
"""

from typing import Any, Callable, Dict, List, Mapping, Optional, Set

import numpy as np
import pandas as pd


# =============================================================================
# VEP FIELD DEFINITIONS
# =============================================================================

VEP_GENE_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "VARIANT_CLASS", "CANONICAL", "MANE",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP",
    "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO",
    "SIFT", "PolyPhen", "DOMAINS", "AF", "gnomADe_AF", "gnomADg_AF",
    "MAX_AF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
]

VEP_REGULATORY_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "Existing_variation", "DISTANCE", "VARIANT_CLASS",
    "AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF", "CLIN_SIG", "SOMATIC",
    "PHENO", "PUBMED",
]

VEP_MOTIF_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "Existing_variation", "VARIANT_CLASS", "MOTIF_NAME",
    "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS",
    "AF", "gnomADe_AF", "gnomADg_AF", "MAX_AF",
]

# Consequence term sets
MISSENSE_TERMS = {"missense_variant"}
NONSENSE_TERMS = {"stop_gained"}
FRAMESHIFT_TERMS = {"frameshift_variant"}
SPLICE_TERMS = {"splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"}


# =============================================================================
# CSQ PARSING
# =============================================================================

def parse_vep_csq(csq_description: str) -> Dict[str, int]:
    """
    Parse VEP CSQ format description into column index mapping.
    """
    if "Format:" in csq_description:
        fmt = csq_description.split("Format:")[1].strip().strip('"').strip()
    else:
        fmt = csq_description.strip().strip('"').strip()

    csq_cols = fmt.split("|")
    return {name: i for i, name in enumerate(csq_cols)}


def _get_field(parts: List[str], name: str, csq_index: Dict[str, int], default=None):
    """Safely extract a field from CSQ parts."""
    idx = csq_index.get(name)
    if idx is None or idx >= len(parts):
        return default
    val = parts[idx]
    return default if val == "" else val


def _has_consequence(gene_hits: List[Dict], target_terms: Set[str]) -> bool:
    """Check if any transcript has a consequence in target_terms."""
    for hit in gene_hits:
        cons = hit.get("Consequence") or ""
        for term in cons.split("&"):
            if term.strip() in target_terms:
                return True
    return False


def _has_consequence_filtered(
    gene_hits: List[Dict],
    target_terms: Set[str],
    filter_field: str,
    filter_predicate: Callable,
) -> bool:
    """Check if any transcript matching filter has a consequence in target_terms."""
    for hit in gene_hits:
        if not filter_predicate(hit.get(filter_field)):
            continue
        cons = hit.get("Consequence") or ""
        for term in cons.split("&"):
            if term.strip() in target_terms:
                return True
    return False


# =============================================================================
# MAIN ANNOTATION FUNCTION
# =============================================================================

def add_vep_hits_columns(
    df: pd.DataFrame,
    csq_description: str,
    primary_genes: List[str],
    lncrna_names: Optional[List[str]] = None,
    gene_symbol_mapping: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """
    Add VEP annotation columns to SV DataFrame.
    
    Adds columns:
    - gene_hits_vep: List of transcript-level annotations
    - regulatory_hits_vep: List of regulatory feature annotations
    - motif_hits: List of motif feature annotations
    - gene_symbols: Comma-separated gene symbols
    - hits_canonical: Boolean flag for canonical transcript hit
    - has_missense, has_nonsense, has_frameshift, has_splice_effect
    - has_*_canonical: Same flags for canonical transcripts only
    - has_*_mane: Same flags for MANE transcripts only
    """
    if lncrna_names is None:
        lncrna_names = []

    from ..genes.symbol_normalization import resolve_symbol_to_panel

    csq_index = parse_vep_csq(csq_description)
    all_genes = set(primary_genes) | set(lncrna_names)

    def parse_csq_for_row(row: pd.Series) -> pd.Series:
        val = row.get("CSQ")
        
        empty_result = pd.Series({
            "gene_hits_vep": [], "regulatory_hits_vep": [], "motif_hits": [],
            "gene_symbols": "", "hits_canonical": False,
            "has_missense": False, "has_nonsense": False,
            "has_frameshift": False, "has_splice_effect": False,
            "has_missense_canonical": False, "has_nonsense_canonical": False,
            "has_frameshift_canonical": False, "has_splice_effect_canonical": False,
            "has_missense_mane": False, "has_nonsense_mane": False,
            "has_frameshift_mane": False, "has_splice_effect_mane": False,
        })
        
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return empty_result

        entries = val if isinstance(val, list) else str(val).split(",")

        gene_hits, regulatory_hits, motif_hits = [], [], []

        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue

            parts = entry.split("|")
            feature_type = _get_field(parts, "Feature_type", csq_index)
            
            if feature_type is None:
                continue

            if feature_type == "Transcript":
                gene_name = _get_field(parts, "SYMBOL", csq_index)
                if gene_symbol_mapping is not None:
                    canonical = resolve_symbol_to_panel(
                        gene_name, all_genes, gene_symbol_mapping
                    )
                else:
                    canonical = gene_name if gene_name in all_genes else None
                if canonical is None:
                    continue
                hit = {fname: _get_field(parts, fname, csq_index) for fname in VEP_GENE_FIELDS}
                hit["SYMBOL"] = canonical
                gene_hits.append(hit)
            elif feature_type == "RegulatoryFeature":
                hit = {fname: _get_field(parts, fname, csq_index) for fname in VEP_REGULATORY_FIELDS}
                regulatory_hits.append(hit)
            elif feature_type == "MotifFeature":
                hit = {fname: _get_field(parts, fname, csq_index) for fname in VEP_MOTIF_FIELDS}
                motif_hits.append(hit)

        gene_symbols = ",".join(sorted({h.get("SYMBOL") for h in gene_hits if h.get("SYMBOL")}))
        hits_canonical = any(h.get("CANONICAL") == "1" for h in gene_hits)

        canonical_pred = lambda v: v == "1"
        mane_pred = lambda v: v is not None

        return pd.Series({
            "gene_hits_vep": gene_hits,
            "regulatory_hits_vep": regulatory_hits,
            "motif_hits": motif_hits,
            "gene_symbols": gene_symbols,
            "hits_canonical": hits_canonical,
            "has_missense": _has_consequence(gene_hits, MISSENSE_TERMS),
            "has_nonsense": _has_consequence(gene_hits, NONSENSE_TERMS),
            "has_frameshift": _has_consequence(gene_hits, FRAMESHIFT_TERMS),
            "has_splice_effect": _has_consequence(gene_hits, SPLICE_TERMS),
            "has_missense_canonical": _has_consequence_filtered(gene_hits, MISSENSE_TERMS, "CANONICAL", canonical_pred),
            "has_nonsense_canonical": _has_consequence_filtered(gene_hits, NONSENSE_TERMS, "CANONICAL", canonical_pred),
            "has_frameshift_canonical": _has_consequence_filtered(gene_hits, FRAMESHIFT_TERMS, "CANONICAL", canonical_pred),
            "has_splice_effect_canonical": _has_consequence_filtered(gene_hits, SPLICE_TERMS, "CANONICAL", canonical_pred),
            "has_missense_mane": _has_consequence_filtered(gene_hits, MISSENSE_TERMS, "MANE_SELECT", mane_pred),
            "has_nonsense_mane": _has_consequence_filtered(gene_hits, NONSENSE_TERMS, "MANE_SELECT", mane_pred),
            "has_frameshift_mane": _has_consequence_filtered(gene_hits, FRAMESHIFT_TERMS, "MANE_SELECT", mane_pred),
            "has_splice_effect_mane": _has_consequence_filtered(gene_hits, SPLICE_TERMS, "MANE_SELECT", mane_pred),
        })

    hits_df = df.apply(parse_csq_for_row, axis=1)
    return pd.concat([df, hits_df], axis=1)


def summarize_vep_consequences(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize VEP consequences across all SVs."""
    summary = {
        "total_svs": len(df),
        "with_vep_gene_hits": (df["gene_hits_vep"].apply(len) > 0).sum() if "gene_hits_vep" in df.columns else 0,
        "with_regulatory_hits": (df["regulatory_hits_vep"].apply(len) > 0).sum() if "regulatory_hits_vep" in df.columns else 0,
        "with_motif_hits": (df["motif_hits"].apply(len) > 0).sum() if "motif_hits" in df.columns else 0,
        "hits_canonical": df["hits_canonical"].sum() if "hits_canonical" in df.columns else 0,
        "has_missense": df["has_missense"].sum() if "has_missense" in df.columns else 0,
        "has_nonsense": df["has_nonsense"].sum() if "has_nonsense" in df.columns else 0,
        "has_frameshift": df["has_frameshift"].sum() if "has_frameshift" in df.columns else 0,
        "has_splice_effect": df["has_splice_effect"].sum() if "has_splice_effect" in df.columns else 0,
    }
    return pd.DataFrame([summary])


def extract_high_impact_svs(df: pd.DataFrame) -> pd.DataFrame:
    """Extract SVs with high-impact VEP consequences."""
    if df.empty:
        return df.copy()
    
    mask = (
        df.get("has_missense", False) |
        df.get("has_nonsense", False) |
        df.get("has_frameshift", False) |
        df.get("has_splice_effect", False)
    )
    return df[mask].copy()
