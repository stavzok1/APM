"""
RPPA analysis functions for immune visibility assessment.

Provides:
- Signaling pathway block detection
- DDR-IFN state stratification
- Protein-RNA discordance analysis
- Immune visibility scoring
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats


# =============================================================================
# SIGNALING BLOCK DETECTION
# =============================================================================

def detect_signaling_blocks(
    expr_matrix: pd.DataFrame,
    quantile_high: float = 0.75,
    quantile_low: float = 0.25,
) -> pd.DataFrame:
    """
    Identify samples with specific signaling pathway blocks.
    
    Detects patterns suggesting disrupted signal transduction:
    
    1. DDR → STING block:
       High γH2AX (DDR active) but low STING (sensing absent)
       Suggests: cGAS-STING pathway disrupted despite DNA damage
    
    2. STING → IRF block:
       High pIRF3 (sensing active) but low IRF1 (downstream silent)
       Suggests: IFN signaling block between sensing and transcription
    
    3. STAT3 override:
       High pSTAT3 alongside IRF1
       Suggests: suppressive signal counteracting IFN response
    
    4. Checkpoint escape:
       High IRF1 (IFN active) and high PD-L1
       Suggests: tumor adapted to immune pressure via checkpoint
    
    5. PI3K escape:
       Low PTEN and high pAKT
       Suggests: PI3K pathway activation (immune-suppressive)
    
    Args:
        expr_matrix: samples × targets DataFrame
        quantile_high: Threshold for "high" expression
        quantile_low: Threshold for "low" expression
    
    Returns:
        DataFrame with boolean columns for each block type
    """
    blocks = pd.DataFrame(index=expr_matrix.index)
    
    # Helper to check if target is high/low
    def is_high(target):
        if target not in expr_matrix.columns:
            return pd.Series(False, index=expr_matrix.index)
        return expr_matrix[target] > expr_matrix[target].quantile(quantile_high)
    
    def is_low(target):
        if target not in expr_matrix.columns:
            return pd.Series(False, index=expr_matrix.index)
        return expr_matrix[target] < expr_matrix[target].quantile(quantile_low)
    
    def is_median_plus(target):
        if target not in expr_matrix.columns:
            return pd.Series(False, index=expr_matrix.index)
        return expr_matrix[target] > expr_matrix[target].median()
    
    # 1. DDR → STING block
    # High DDR activation but low STING
    ddr_active = is_high("H2AX_pS139") | is_high("ATM_pS1981")
    sting_low = is_low("STING")
    blocks["ddr_sting_block"] = ddr_active & sting_low
    
    # 2. STING → IRF1 block
    # pIRF3 active but IRF1 not induced
    irf3_active = is_high("IRF-3_pS396")
    irf1_low = is_low("IRF1")
    blocks["sting_irf1_block"] = irf3_active & irf1_low
    
    # 3. STAT3 suppressive override
    # High pSTAT3 even with IFN signaling markers present
    stat3_high = is_high("STAT3_pY705")
    irf1_present = is_median_plus("IRF1")
    blocks["stat3_override"] = stat3_high & irf1_present
    
    # 4. Checkpoint escape
    # Tumor with IFN response but high checkpoint
    irf1_high = is_high("IRF1")
    pdl1_high = is_high("PDL1")
    blocks["checkpoint_escape"] = irf1_high & pdl1_high
    
    # 5. PI3K pathway activation
    # PTEN loss or high AKT phosphorylation
    pten_low = is_low("PTEN")
    akt_high = is_high("AKT_pS473") | is_high("AKT_pT308")
    blocks["pi3k_activated"] = pten_low | akt_high
    
    # 6. Complete DDR-IFN chain intact
    # DDR → STING → IRF3 → IRF1 all active
    sting_present = is_median_plus("STING")
    irf3_present = is_median_plus("IRF-3_pS396")
    blocks["ddr_ifn_chain_intact"] = (
        ddr_active & sting_present & irf3_present & irf1_present
    )
    
    # 7. Immune desert (everything low)
    blocks["immune_desert"] = (
        is_low("IRF1") & 
        is_low("PDL1") & 
        is_low("Granzyme-B") &
        is_low("CD45")
    )
    
    return blocks


def classify_ddr_ifn_quadrant(
    ddr_score: pd.Series,
    ifn_score: pd.Series,
) -> pd.Series:
    """
    Stratify samples into DDR-IFN quadrants.
    
    Quadrants (based on median split):
    
    - hot_primed: DDR high, IFN high
      → DDR driving immune activation (your "hot" state)
    
    - sensing_failure: DDR high, IFN low
      → DDR present but not transduced to IFN (potential cGAS-STING block)
    
    - external_ifn: DDR low, IFN high
      → IFN from microenvironment, not tumor-intrinsic DDR
    
    - cold: DDR low, IFN low
      → Neither DDR nor IFN active (immune invisible)
    
    Args:
        ddr_score: DDR activation score per sample
        ifn_score: IFN signaling score per sample
    
    Returns:
        Series with quadrant labels
    """
    ddr_med = ddr_score.median()
    ifn_med = ifn_score.median()
    
    # Align indices
    common_idx = ddr_score.index.intersection(ifn_score.index)
    
    result = pd.Series("unknown", index=common_idx)
    
    for idx in common_idx:
        ddr_high = ddr_score.loc[idx] > ddr_med
        ifn_high = ifn_score.loc[idx] > ifn_med
        
        if ddr_high and ifn_high:
            result.loc[idx] = "hot_primed"
        elif ddr_high and not ifn_high:
            result.loc[idx] = "sensing_failure"
        elif not ddr_high and ifn_high:
            result.loc[idx] = "external_ifn"
        else:
            result.loc[idx] = "cold"
    
    return result


def classify_immune_visibility_state(
    expr_matrix: pd.DataFrame,
    panel_scores: pd.DataFrame,
) -> pd.DataFrame:
    """
    Comprehensive immune visibility state classification.
    
    Integrates multiple signals to classify each sample:
    
    States:
    - hot_functional: High IFN, high cytolytic, low checkpoint
    - hot_suppressed: High IFN, high cytolytic, high checkpoint
    - warm_infiltrated: Moderate IFN, some cytolytic activity
    - cold_ddr_active: Low IFN but DDR active (sensing failure)
    - cold_true: Low everything (immune invisible)
    
    Returns:
        DataFrame with classification and supporting evidence
    """
    result = pd.DataFrame(index=expr_matrix.index)
    
    # Get relevant scores
    ifn = panel_scores.get("IFN_activated", pd.Series(0, index=expr_matrix.index))
    ddr = panel_scores.get("DDR_activation", pd.Series(0, index=expr_matrix.index))
    cyto = panel_scores.get("cytolytic", pd.Series(0, index=expr_matrix.index))
    chkpt = panel_scores.get("checkpoint", pd.Series(0, index=expr_matrix.index))
    
    # Median thresholds
    ifn_med = ifn.median()
    ddr_med = ddr.median()
    cyto_med = cyto.median()
    chkpt_med = chkpt.median()
    
    # Classify
    def classify_sample(idx):
        i = ifn.get(idx, 0)
        d = ddr.get(idx, 0)
        c = cyto.get(idx, 0)
        p = chkpt.get(idx, 0)
        
        ifn_high = i > ifn_med
        ddr_high = d > ddr_med
        cyto_high = c > cyto_med
        chkpt_high = p > chkpt_med
        
        if ifn_high and cyto_high and not chkpt_high:
            return "hot_functional"
        elif ifn_high and cyto_high and chkpt_high:
            return "hot_suppressed"
        elif ifn_high or cyto_high:
            return "warm_infiltrated"
        elif ddr_high:
            return "cold_ddr_active"
        else:
            return "cold_true"
    
    result["visibility_state"] = [classify_sample(idx) for idx in expr_matrix.index]
    
    # Add supporting scores
    result["IFN_score"] = ifn
    result["DDR_score"] = ddr
    result["cytolytic_score"] = cyto
    result["checkpoint_score"] = chkpt
    
    return result


# =============================================================================
# PROTEIN-RNA DISCORDANCE
# =============================================================================

def compute_protein_rna_discordance(
    rppa_expr: pd.DataFrame,
    rna_expr: pd.DataFrame,
    target_to_gene: Dict[str, str],
    genes_of_interest: Optional[List[str]] = None,
    discordance_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Identify samples with protein/RNA discordance.
    
    High discordance suggests post-transcriptional regulation:
    - High RNA, low protein: degradation, translation block
    - Low RNA, high protein: protein stabilization
    
    This is where RPPA adds unique value over transcriptomics alone.
    
    Args:
        rppa_expr: Protein expression matrix (samples × targets)
        rna_expr: RNA expression matrix (samples × genes)
        target_to_gene: Mapping from RPPA target to gene name
        genes_of_interest: Specific genes to analyze
        discordance_threshold: Z-score difference threshold
    
    Returns:
        DataFrame with columns per gene:
        - {gene}_rna_z: Z-scored RNA expression
        - {gene}_protein_z: Z-scored protein expression
        - {gene}_discordance: Protein - RNA (z-score space)
        - {gene}_discordance_dir: "protein_high", "protein_low", "concordant"
    """
    if genes_of_interest is None:
        # Focus on key visibility genes
        genes_of_interest = [
            "IRF1", "STAT3", "PTEN", "CD274",  # Direct visibility
            "ATM", "ATR", "CHEK1",  # DDR
            "AKT1", "MTOR",  # Suppressive
        ]
    
    # Build reverse map: gene → protein targets
    gene_to_targets = {}
    for target, gene in target_to_gene.items():
        if gene not in gene_to_targets:
            gene_to_targets[gene] = []
        gene_to_targets[gene].append(target)
    
    results = pd.DataFrame(index=rppa_expr.index)
    
    for gene in genes_of_interest:
        # Get RNA
        if gene not in rna_expr.columns:
            continue
        
        rna = rna_expr[gene]
        
        # Get protein (prefer non-phospho "total" target)
        if gene not in gene_to_targets:
            # Try common name variations
            gene_variants = [gene, gene.upper(), gene.lower(), gene.capitalize()]
            found = False
            for variant in gene_variants:
                if variant in gene_to_targets:
                    gene = variant
                    found = True
                    break
            if not found:
                continue
        
        targets = gene_to_targets[gene]
        
        # Prefer non-phospho target for total protein comparison
        protein_target = None
        for t in targets:
            if "_p" not in t.lower():
                protein_target = t
                break
        if protein_target is None:
            protein_target = targets[0]
        
        if protein_target not in rppa_expr.columns:
            continue
        
        protein = rppa_expr[protein_target]
        
        # Z-score both
        rna_z = (rna - rna.mean()) / rna.std()
        protein_z = (protein - protein.mean()) / protein.std()
        
        # Align indices
        common = rna_z.index.intersection(protein_z.index)
        if len(common) == 0:
            continue
        
        rna_z = rna_z.loc[common]
        protein_z = protein_z.loc[common]
        
        # Discordance
        discordance = protein_z - rna_z
        
        # Direction classification
        def classify(d):
            if pd.isna(d):
                return "unknown"
            if d > discordance_threshold:
                return "protein_high"
            elif d < -discordance_threshold:
                return "protein_low"
            return "concordant"
        
        results.loc[common, f"{gene}_rna_z"] = rna_z
        results.loc[common, f"{gene}_protein_z"] = protein_z
        results.loc[common, f"{gene}_discordance"] = discordance
        results.loc[common, f"{gene}_discordance_dir"] = discordance.apply(classify)
    
    return results


def identify_post_transcriptional_regulation(
    discordance_df: pd.DataFrame,
    min_discordant_samples: int = 5,
) -> Dict[str, Dict]:
    """
    Identify genes showing systematic post-transcriptional regulation.
    
    Args:
        discordance_df: Output from compute_protein_rna_discordance
        min_discordant_samples: Minimum samples with discordance
    
    Returns:
        Dict[gene] = {
            "n_protein_high": count,
            "n_protein_low": count,
            "n_concordant": count,
            "dominant_pattern": "stabilization" / "degradation" / "mixed" / "concordant",
            "samples_protein_high": [sample_ids],
            "samples_protein_low": [sample_ids],
        }
    """
    # Find genes with discordance columns
    genes = set()
    for col in discordance_df.columns:
        if col.endswith("_discordance_dir"):
            gene = col.replace("_discordance_dir", "")
            genes.add(gene)
    
    results = {}
    
    for gene in genes:
        dir_col = f"{gene}_discordance_dir"
        if dir_col not in discordance_df.columns:
            continue
        
        directions = discordance_df[dir_col].dropna()
        
        n_high = (directions == "protein_high").sum()
        n_low = (directions == "protein_low").sum()
        n_conc = (directions == "concordant").sum()
        
        # Determine dominant pattern
        total_discordant = n_high + n_low
        if total_discordant < min_discordant_samples:
            pattern = "concordant"
        elif n_high > n_low * 2:
            pattern = "stabilization"  # Protein higher than expected from RNA
        elif n_low > n_high * 2:
            pattern = "degradation"  # Protein lower than expected from RNA
        else:
            pattern = "mixed"
        
        samples_high = discordance_df.index[
            discordance_df[dir_col] == "protein_high"
        ].tolist()
        samples_low = discordance_df.index[
            discordance_df[dir_col] == "protein_low"
        ].tolist()
        
        results[gene] = {
            "n_protein_high": n_high,
            "n_protein_low": n_low,
            "n_concordant": n_conc,
            "dominant_pattern": pattern,
            "samples_protein_high": samples_high,
            "samples_protein_low": samples_low,
        }
    
    return results


# =============================================================================
# IMMUNE VISIBILITY SCORING
# =============================================================================

def compute_immune_visibility_score(
    panel_scores: pd.DataFrame,
    expr_matrix: pd.DataFrame,
    weights: Optional[Dict[str, float]] = None,
) -> pd.Series:
    """
    Compute composite immune visibility score.
    
    Combines multiple indicators:
    - IFN signaling (positive)
    - cGAS-STING sensing (positive)
    - Cytolytic activity (positive)
    - Checkpoint expression (context-dependent)
    - STAT3 suppression (negative)
    - PI3K activation (negative)
    
    Args:
        panel_scores: Output from compute_all_panel_scores
        expr_matrix: Raw expression for specific markers
        weights: Optional custom weights for components
    
    Returns:
        Series of visibility scores (higher = more visible)
    """
    if weights is None:
        weights = {
            "IFN_activated": 1.0,
            "cGAS_STING": 0.8,
            "cytolytic": 0.7,
            "checkpoint": 0.3,  # Ambiguous - can indicate both visibility and evasion
            "STAT3_suppressive": -0.5,
            "PI3K_AKT": -0.4,
        }
    
    score = pd.Series(0.0, index=panel_scores.index)
    
    for component, weight in weights.items():
        if component in panel_scores.columns:
            # Z-score the component first
            vals = panel_scores[component]
            z = (vals - vals.mean()) / vals.std()
            score += weight * z.fillna(0)
    
    # Normalize to 0-100 scale
    score = (score - score.min()) / (score.max() - score.min()) * 100
    
    return score


def compute_antigen_presentation_capacity(
    expr_matrix: pd.DataFrame,
    rna_expr: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Estimate antigen presentation capacity from available markers.
    
    Note: RPPA lacks direct APM markers (TAP1, TAP2, B2M, HLA-I, PSMB8/9/10).
    This function provides indirect estimation through:
    - Trans-regulator activity (IRF1)
    - Suppressive pathways (STAT3, PI3K)
    - If RNA provided: direct APM gene expression
    
    Returns:
        DataFrame with:
        - trans_regulator_score: IRF1, IRF3 activity
        - suppressive_score: STAT3, PI3K-AKT
        - estimated_apm_capacity: Combined estimate
        - apm_rna_score: Direct APM RNA (if rna_expr provided)
    """
    result = pd.DataFrame(index=expr_matrix.index)
    
    # Trans-regulator activity (proxy for APM induction)
    trans_markers = ["IRF1", "IRF-3_pS396", "CIITA"]
    available_trans = [m for m in trans_markers if m in expr_matrix.columns]
    
    if available_trans:
        trans_z = expr_matrix[available_trans].apply(
            lambda x: (x - x.mean()) / x.std()
        )
        result["trans_regulator_score"] = trans_z.mean(axis=1)
    else:
        result["trans_regulator_score"] = np.nan
    
    # Suppressive signals (negative for APM)
    supp_markers = ["STAT3_pY705", "AKT_pS473", "AKT_pT308"]
    available_supp = [m for m in supp_markers if m in expr_matrix.columns]
    
    if available_supp:
        supp_z = expr_matrix[available_supp].apply(
            lambda x: (x - x.mean()) / x.std()
        )
        result["suppressive_score"] = supp_z.mean(axis=1)
    else:
        result["suppressive_score"] = np.nan
    
    # Combined estimate
    if "trans_regulator_score" in result and "suppressive_score" in result:
        # Higher trans, lower suppressive = better APM capacity
        result["estimated_apm_capacity"] = (
            result["trans_regulator_score"] - 0.5 * result["suppressive_score"]
        )
    
    # Direct APM RNA if available
    if rna_expr is not None:
        apm_genes = [
            "TAP1", "TAP2", "B2M", "HLA-A", "HLA-B", "HLA-C",
            "PSMB8", "PSMB9", "PSMB10", "NLRC5"
        ]
        available_apm = [g for g in apm_genes if g in rna_expr.columns]
        
        if available_apm:
            # Align indices
            common = result.index.intersection(rna_expr.index)
            apm_z = rna_expr.loc[common, available_apm].apply(
                lambda x: (x - x.mean()) / x.std()
            )
            result.loc[common, "apm_rna_score"] = apm_z.mean(axis=1)
    
    return result


# =============================================================================
# STATISTICAL TESTS
# =============================================================================

def test_score_association(
    scores: pd.Series,
    groups: pd.Series,
    test: str = "auto",
) -> Dict[str, float]:
    """
    Test association between continuous score and categorical groups.
    
    Args:
        scores: Continuous scores (e.g., visibility score)
        groups: Categorical groups (e.g., PAM50 subtype)
        test: Test type ("auto", "kruskal", "anova", "mannwhitney")
    
    Returns:
        Dict with test results
    """
    # Align
    common = scores.index.intersection(groups.index)
    scores = scores.loc[common].dropna()
    groups = groups.loc[scores.index]
    
    unique_groups = groups.unique()
    n_groups = len(unique_groups)
    
    if n_groups < 2:
        return {"error": "Need at least 2 groups"}
    
    # Auto-select test
    if test == "auto":
        test = "kruskal" if n_groups > 2 else "mannwhitney"
    
    group_values = [scores[groups == g].values for g in unique_groups]
    
    if test == "kruskal":
        stat, pval = stats.kruskal(*group_values)
        test_name = "Kruskal-Wallis"
    elif test == "anova":
        stat, pval = stats.f_oneway(*group_values)
        test_name = "One-way ANOVA"
    elif test == "mannwhitney":
        if n_groups != 2:
            return {"error": "Mann-Whitney requires exactly 2 groups"}
        stat, pval = stats.mannwhitneyu(group_values[0], group_values[1])
        test_name = "Mann-Whitney U"
    else:
        return {"error": f"Unknown test: {test}"}
    
    # Effect size (eta-squared for multi-group)
    grand_mean = scores.mean()
    ss_between = sum(
        len(g) * (g.mean() - grand_mean) ** 2 
        for g in group_values
    )
    ss_total = ((scores - grand_mean) ** 2).sum()
    eta_squared = ss_between / ss_total if ss_total > 0 else 0
    
    return {
        "test": test_name,
        "statistic": stat,
        "p_value": pval,
        "n_groups": n_groups,
        "n_samples": len(scores),
        "eta_squared": eta_squared,
        "group_means": {g: scores[groups == g].mean() for g in unique_groups},
        "group_ns": {g: (groups == g).sum() for g in unique_groups},
    }


def correlate_with_outcomes(
    scores_df: pd.DataFrame,
    outcomes_df: pd.DataFrame,
    method: str = "spearman",
) -> pd.DataFrame:
    """
    Correlate RPPA scores with clinical/molecular outcomes.
    
    Args:
        scores_df: RPPA panel scores
        outcomes_df: Outcome variables (cytolytic, infiltration, etc.)
        method: Correlation method ("spearman" or "pearson")
    
    Returns:
        DataFrame with correlations and p-values
    """
    results = []
    
    # Align indices
    common = scores_df.index.intersection(outcomes_df.index)
    scores_aligned = scores_df.loc[common]
    outcomes_aligned = outcomes_df.loc[common]
    
    for score_col in scores_aligned.columns:
        for outcome_col in outcomes_aligned.columns:
            x = scores_aligned[score_col].dropna()
            y = outcomes_aligned[outcome_col].loc[x.index].dropna()
            
            common_xy = x.index.intersection(y.index)
            if len(common_xy) < 10:
                continue
            
            x = x.loc[common_xy]
            y = y.loc[common_xy]
            
            if method == "spearman":
                r, p = stats.spearmanr(x, y)
            else:
                r, p = stats.pearsonr(x, y)
            
            results.append({
                "score": score_col,
                "outcome": outcome_col,
                "correlation": r,
                "p_value": p,
                "n": len(common_xy),
            })
    
    return pd.DataFrame(results)
