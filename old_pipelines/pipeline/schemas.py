"""
Schema definitions and builders for nested data structures.

These helpers ensure consistent structure across all evidence types,
even when data is missing for certain biosamples or assays.
"""

from typing import Dict, Any, List, Optional


# =============================================================================
# SCREEN EVIDENCE SCHEMAS
# =============================================================================

def empty_assay_entry() -> Dict[str, Any]:
    """Empty entry for a single assay type within a biosample."""
    return {
        "score": None,
        "p_value": None,
        "strength": "none",
    }


def empty_biosample_assays(assay_types: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty dict of assays for a single biosample."""
    return {assay: empty_assay_entry() for assay in assay_types}


def empty_conservation_entry() -> Dict[str, int]:
    """Empty conservation counts for a single assay type."""
    return {
        "n_biosamples": 0,
        "n_strong": 0,
        "n_weak": 0,
    }


def empty_conservation_block(assay_types: List[str]) -> Dict[str, Dict[str, int]]:
    """Empty conservation dict across all assay types."""
    return {assay: empty_conservation_entry() for assay in assay_types}


def empty_screen_block(biosamples: List[str], assay_types: List[str]) -> Dict[str, Any]:
    """
    Empty SCREEN evidence block structure:
    {
        "per_biosample": {biosample: {assay: {score, p_value, strength}, ...}, ...},
        "conservation_global": {assay: {n_biosamples, n_strong, n_weak}, ...},
        "conservation_breast": {assay: {n_biosamples, n_strong, n_weak}, ...},
    }
    """
    return {
        "per_biosample": {bio: empty_biosample_assays(assay_types) for bio in biosamples},
        "conservation_global": empty_conservation_block(assay_types),
        "conservation_breast": empty_conservation_block(assay_types),
    }


def ensure_screen_block(
    val: Any,
    biosamples: List[str],
    assay_types: List[str],
) -> Dict[str, Any]:
    """
    Coerce any value into a valid SCREEN block, filling missing keys.
    Use this when merging data that might have gaps.
    """
    if not isinstance(val, dict):
        return empty_screen_block(biosamples, assay_types)

    val.setdefault("per_biosample", {})
    val.setdefault("conservation_global", {})
    val.setdefault("conservation_breast", {})

    # Ensure all biosamples present
    for bio in biosamples:
        if bio not in val["per_biosample"] or not isinstance(val["per_biosample"][bio], dict):
            val["per_biosample"][bio] = empty_biosample_assays(assay_types)
        else:
            for assay in assay_types:
                val["per_biosample"][bio].setdefault(assay, empty_assay_entry())

    # Ensure conservation blocks complete
    for cons_key in ["conservation_global", "conservation_breast"]:
        if not isinstance(val[cons_key], dict):
            val[cons_key] = empty_conservation_block(assay_types)
        else:
            for assay in assay_types:
                val[cons_key].setdefault(assay, empty_conservation_entry())

    return val


# =============================================================================
# ABC EVIDENCE SCHEMAS
# =============================================================================

def empty_abc_celltype_entry() -> Dict[str, Any]:
    """Empty ABC entry for a single cell type."""
    return {
        "ABC_score": None,
        "ABC_num": None,
        "activity": None,
        "distance": None,
        "element_class": None,
        "is_self_promoter": False,
        "hic_pl_scaled": None,
        "powerlaw_score": None,
        "gene_expr": None,
        "promoter_activity_q": None,
        "gene_is_expressed": False,
        "rank_within_gene": None,
        "is_present": False,
        "is_strong": False,
    }


def empty_abc_block(celltypes: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty ABC block across all cell types."""
    return {ct: empty_abc_celltype_entry() for ct in celltypes}


def empty_abc_enhancer_entry(celltypes: List[str]) -> Dict[str, Any]:
    """Empty entry for a single ABC enhancer mapped to a cCRE."""
    return {
        "start": None,
        "end": None,
        "ABC_full": empty_abc_block(celltypes),
    }


# =============================================================================
# HICHIP EVIDENCE SCHEMAS
# =============================================================================

def empty_hichip_loop_entry() -> Dict[str, Any]:
    """Empty entry for a single HiChIP loop."""
    return {
        "loop_id": None,
        "counts": None,
        "score": None,
        "n_reps": 0,
        "partner_genes": [],
        "partner_cCREs": [],
        "anchor_a": {"chr": None, "start": None, "end": None},
        "anchor_b": {"chr": None, "start": None, "end": None},
    }


def empty_hichip_celltype_entry() -> Dict[str, Any]:
    """Empty HiChIP summary for a single cell type."""
    return {
        "n_loops": 0,
        "max_counts": None,
        "max_score": None,
        "loops": [],
    }


def empty_hichip_block(celltypes: List[str]) -> Dict[str, Dict[str, Any]]:
    """Empty HiChIP block across all cell types."""
    return {ct: empty_hichip_celltype_entry() for ct in celltypes}


# =============================================================================
# GENE LINKS SCHEMA (COMBINED)
# =============================================================================

def empty_gene_link_entry(
    screen_exp_biosamples: List[str],
    screen_exp_assays: List[str],
    screen_comp_biosamples: List[str],
    screen_comp_assays: List[str],
    abc_celltypes: Optional[List[str]] = None,
    hichip_celltypes: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Complete empty gene_link entry with all evidence types.
    
    Structure:
    {
        "gene_id": str,
        "gene_type": str,
        "screen_exp": {per_biosample, conservation_global, conservation_breast},
        "screen_comp": {per_biosample, conservation_global, conservation_breast},
        "ABC_enhancers": [{start, end, ABC_full}, ...],
        "hichip": {celltype: {n_loops, max_counts, max_score, loops}, ...},
    }
    """
    entry = {
        "gene_id": None,
        "gene_type": None,
        "screen_exp": empty_screen_block(screen_exp_biosamples, screen_exp_assays),
        "screen_comp": empty_screen_block(screen_comp_biosamples, screen_comp_assays),
        "ABC_enhancers": [],
        "hichip": {},
    }
    
    if hichip_celltypes:
        entry["hichip"] = empty_hichip_block(hichip_celltypes)
    
    return entry


# =============================================================================
# CELL LINE SIGNAL SCHEMAS
# =============================================================================

def empty_cell_line_signals(signal_names: List[str]) -> Dict[str, Any]:
    """
    Empty signal dict for a single cell line.
    
    Structure:
    {
        "in_{cell_line}": False,
        "H3K27ac": None,
        "CTCF": None,
        ...
    }
    """
    return {sig: None for sig in signal_names}


def ensure_cell_line_signals(
    val: Any,
    signal_names: List[str],
    in_key: str,
) -> Dict[str, Any]:
    """Ensure cell line signal dict has all expected keys."""
    if not isinstance(val, dict):
        result = empty_cell_line_signals(signal_names)
        result[in_key] = False
        return result
    
    val.setdefault(in_key, False)
    for sig in signal_names:
        val.setdefault(sig, None)
    
    return val


# =============================================================================
# MIRNA SCHEMA
# =============================================================================

def empty_mirna_entry() -> Dict[str, Any]:
    """Empty miRNA target entry."""
    return {
        "miRNA": None,
        "weighted_context_score": None,
        "num_sites": 0,
        "UTR_starts": [],
        "UTR_ends": [],
    }


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

def validate_gene_links_structure(gene_links: Dict[str, Dict]) -> List[str]:
    """
    Validate a gene_links dict and return list of issues.
    Returns empty list if valid.
    """
    issues = []
    
    if not isinstance(gene_links, dict):
        return ["gene_links is not a dict"]
    
    for gene_name, entry in gene_links.items():
        if not isinstance(entry, dict):
            issues.append(f"{gene_name}: entry is not a dict")
            continue
        
        # Check required keys
        for key in ["screen_exp", "screen_comp", "ABC_enhancers", "hichip"]:
            if key not in entry:
                issues.append(f"{gene_name}: missing key '{key}'")
        
        # Check screen blocks
        for screen_key in ["screen_exp", "screen_comp"]:
            if screen_key in entry and isinstance(entry[screen_key], dict):
                for sub_key in ["per_biosample", "conservation_global", "conservation_breast"]:
                    if sub_key not in entry[screen_key]:
                        issues.append(f"{gene_name}.{screen_key}: missing '{sub_key}'")
    
    return issues
