"""
RPPA schema definitions for consistent data structures.

Provides empty/default structures and validation helpers
for RPPA data integration.
"""

from typing import Dict, Any, List, Optional


# =============================================================================
# RPPA SAMPLE ENTRY SCHEMAS
# =============================================================================

def empty_rppa_sample_entry() -> Dict[str, Any]:
    """
    Empty entry for RPPA data per sample.
    
    Structure mirrors protein-level evidence in gene_links.
    """
    return {
        "sample_id": None,
        "expression": {},  # target → value
        "activation_ratios": {},  # gene → ratio
        "panel_scores": {},  # panel_name → score
        "signaling_blocks": {},  # block_type → bool
        "visibility_state": None,
    }


def empty_rppa_target_entry() -> Dict[str, Any]:
    """
    Empty entry for a single RPPA target.
    """
    return {
        "target": None,
        "gene_name": None,
        "expression": None,
        "z_score": None,
        "percentile": None,
        "validation_status": None,
        "is_phospho": False,
        "phospho_site": None,
    }


def empty_rppa_gene_evidence() -> Dict[str, Any]:
    """
    Empty RPPA evidence block for a gene.
    
    Structure:
    {
        "total_protein": {target, expression, z_score},
        "phospho_forms": [{target, site, expression, z_score, activation_ratio}, ...],
        "validation_status": "Valid" | "Caution",
    }
    """
    return {
        "total_protein": {
            "target": None,
            "expression": None,
            "z_score": None,
        },
        "phospho_forms": [],
        "validation_status": None,
    }


def empty_rppa_panel_scores() -> Dict[str, Optional[float]]:
    """
    Empty panel scores dictionary.
    """
    return {
        "IFN_activated": None,
        "DDR_activation": None,
        "cGAS_STING": None,
        "DNA_repair": None,
        "checkpoint": None,
        "PI3K_AKT": None,
        "mTOR": None,
        "STAT3_suppressive": None,
        "lymphocyte_infiltration": None,
        "cytolytic": None,
        "proliferation": None,
        "apoptosis_balance": None,
    }


def empty_signaling_blocks() -> Dict[str, bool]:
    """
    Empty signaling blocks dictionary.
    """
    return {
        "ddr_sting_block": False,
        "sting_irf1_block": False,
        "stat3_override": False,
        "checkpoint_escape": False,
        "pi3k_activated": False,
        "ddr_ifn_chain_intact": False,
        "immune_desert": False,
    }


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

def validate_rppa_expression_matrix(
    expr_matrix,  # pd.DataFrame
) -> List[str]:
    """
    Validate RPPA expression matrix structure.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    # Check type
    import pandas as pd
    if not isinstance(expr_matrix, pd.DataFrame):
        return ["Expression matrix is not a DataFrame"]
    
    # Check not empty
    if expr_matrix.empty:
        issues.append("Expression matrix is empty")
    
    # Check index (samples)
    if expr_matrix.index.duplicated().any():
        n_dup = expr_matrix.index.duplicated().sum()
        issues.append(f"Expression matrix has {n_dup} duplicate sample IDs")
    
    # Check columns (targets)
    if expr_matrix.columns.duplicated().any():
        n_dup = expr_matrix.columns.duplicated().sum()
        issues.append(f"Expression matrix has {n_dup} duplicate target names")
    
    # Check values are numeric
    non_numeric = []
    for col in expr_matrix.columns:
        if not pd.api.types.is_numeric_dtype(expr_matrix[col]):
            non_numeric.append(col)
    if non_numeric:
        issues.append(f"Non-numeric columns: {non_numeric[:5]}...")
    
    return issues


def validate_rppa_annotation(
    annotation,  # pd.DataFrame
) -> List[str]:
    """
    Validate RPPA annotation DataFrame.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    import pandas as pd
    if not isinstance(annotation, pd.DataFrame):
        return ["Annotation is not a DataFrame"]
    
    # Required columns
    required = ["AGID", "peptide_target", "gene_name", "validation_status"]
    missing = [c for c in required if c not in annotation.columns]
    if missing:
        issues.append(f"Missing required columns: {missing}")
    
    # Check for duplicates
    if "peptide_target" in annotation.columns:
        if annotation["peptide_target"].duplicated().any():
            n_dup = annotation["peptide_target"].duplicated().sum()
            issues.append(f"Annotation has {n_dup} duplicate peptide_targets")
    
    # Validate validation_status values
    if "validation_status" in annotation.columns:
        valid_statuses = {"Valid", "Caution", "Under Evaluation"}
        actual = set(annotation["validation_status"].str.strip().str.title().unique())
        unknown = actual - valid_statuses
        if unknown:
            issues.append(f"Unknown validation statuses: {unknown}")
    
    return issues


def validate_phospho_pairs(
    phospho_pairs: Dict[str, Dict],
) -> List[str]:
    """
    Validate phospho pairs structure.
    
    Returns:
        List of issues (empty if valid)
    """
    issues = []
    
    if not isinstance(phospho_pairs, dict):
        return ["phospho_pairs is not a dict"]
    
    for gene, pair in phospho_pairs.items():
        if not isinstance(pair, dict):
            issues.append(f"{gene}: pair entry is not a dict")
            continue
        
        if "total" not in pair:
            issues.append(f"{gene}: missing 'total' key")
        if "phospho" not in pair:
            issues.append(f"{gene}: missing 'phospho' key")
        elif not isinstance(pair["phospho"], list):
            issues.append(f"{gene}: 'phospho' is not a list")
    
    return issues


# =============================================================================
# SCHEMA BUILDERS
# =============================================================================

def build_rppa_gene_evidence(
    gene_name: str,
    expr_matrix,  # pd.DataFrame
    target_to_gene: Dict[str, str],
    phospho_pairs: Dict[str, Dict],
    sample_id: str,
) -> Dict[str, Any]:
    """
    Build RPPA evidence structure for a specific gene and sample.
    
    Args:
        gene_name: Gene of interest
        expr_matrix: Expression matrix (samples × targets)
        target_to_gene: Mapping from target to gene
        phospho_pairs: Phospho pair information
        sample_id: Sample to extract data for
    
    Returns:
        RPPA evidence dict for this gene-sample pair
    """
    evidence = empty_rppa_gene_evidence()
    
    if sample_id not in expr_matrix.index:
        return evidence
    
    # Find targets for this gene
    gene_targets = [t for t, g in target_to_gene.items() if g == gene_name]
    
    if not gene_targets:
        return evidence
    
    # Get phospho pair info
    pair_info = phospho_pairs.get(gene_name, {})
    total_target = pair_info.get("total")
    phospho_targets = pair_info.get("phospho", [])
    phospho_sites = pair_info.get("phospho_sites", {})
    
    # Total protein
    if total_target and total_target in expr_matrix.columns:
        val = expr_matrix.loc[sample_id, total_target]
        col_vals = expr_matrix[total_target]
        z = (val - col_vals.mean()) / col_vals.std() if col_vals.std() > 0 else 0
        
        evidence["total_protein"] = {
            "target": total_target,
            "expression": float(val) if not pd.isna(val) else None,
            "z_score": float(z) if not pd.isna(z) else None,
        }
    
    # Phospho forms
    for p_target in phospho_targets:
        if p_target not in expr_matrix.columns:
            continue
        
        val = expr_matrix.loc[sample_id, p_target]
        col_vals = expr_matrix[p_target]
        z = (val - col_vals.mean()) / col_vals.std() if col_vals.std() > 0 else 0
        
        # Activation ratio if total available
        activation_ratio = None
        if total_target and total_target in expr_matrix.columns:
            total_val = expr_matrix.loc[sample_id, total_target]
            if not pd.isna(val) and not pd.isna(total_val):
                activation_ratio = float(val - total_val)  # Log-space ratio
        
        evidence["phospho_forms"].append({
            "target": p_target,
            "site": phospho_sites.get(p_target, ""),
            "expression": float(val) if not pd.isna(val) else None,
            "z_score": float(z) if not pd.isna(z) else None,
            "activation_ratio": activation_ratio,
        })
    
    return evidence


def ensure_rppa_panel_scores(
    val: Any,
) -> Dict[str, Optional[float]]:
    """
    Ensure panel scores dict has all expected keys.
    """
    if not isinstance(val, dict):
        return empty_rppa_panel_scores()
    
    template = empty_rppa_panel_scores()
    for key in template:
        val.setdefault(key, None)
    
    return val


def ensure_signaling_blocks(
    val: Any,
) -> Dict[str, bool]:
    """
    Ensure signaling blocks dict has all expected keys.
    """
    if not isinstance(val, dict):
        return empty_signaling_blocks()
    
    template = empty_signaling_blocks()
    for key in template:
        val.setdefault(key, False)
    
    return val


# Need pandas for validation
import pandas as pd
