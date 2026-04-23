"""
RPPA marker panel definitions and composite score computation.

Defines biologically-informed marker panels for:
- IFN-γ signaling cascade (trans-regulatory axis)
- DDR → cGAS-STING → IFN sensing chain
- Checkpoint/evasion molecules
- Suppressive signaling pathways
- Immune context markers

Provides scoring functions for pathway activation assessment.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Union
from pathlib import Path

import pandas as pd
import numpy as np


# =============================================================================
# MARKER PANEL DEFINITIONS
# =============================================================================

@dataclass
class RPPAMarkerPanels:
    """
    Curated marker panels organized by biological function.
    
    Each panel maps to a specific component of the immune visibility model:
    
    1. IFN signaling: STAT1 → IRF1 → NLRC5 → APM transcription
       (Note: STAT1 is missing from RPPA, use IRF1 as functional proxy)
    
    2. DDR-to-IFN: DNA damage → cGAS-STING → Type I IFN → amplification
    
    3. Checkpoint/evasion: Direct immune visibility output
    
    4. Suppressive signaling: PI3K-AKT-mTOR, STAT3 (counter-regulatory)
    
    5. Immune context: For stratification and adjustment
    """
    
    # -------------------------------------------------------------------------
    # IFN-γ signaling cascade (trans-regulatory axis)
    # -------------------------------------------------------------------------
    ifn_signaling: Dict[str, List[str]] = field(default_factory=lambda: {
        # Phospho-activated forms (signaling ON)
        "activated": [
            "IRF1",              # Key readout - downstream of STAT1
            "IRF-3_pS396",       # Innate IFN activation
            "STAT3_pY705",       # Counter-regulatory (high = suppressive)
            "STAT5ALPHA",        # Alternative STAT
            "JAK2",              # Upstream kinase
            "NFKBP65_pS536",     # NF-κB activation
        ],
        # Total protein forms
        "total": [
            "IRF-3",
            "Stat3",
        ],
    })
    
    # -------------------------------------------------------------------------
    # DDR → cGAS-STING → IFN chain
    # Validates: is DDR actually triggering innate sensing?
    # -------------------------------------------------------------------------
    ddr_to_ifn_sensing: Dict[str, List[str]] = field(default_factory=lambda: {
        # DDR activation markers (phosphorylated = active)
        "ddr_activation": [
            "H2AX_pS139",        # γH2AX - active DSBs (primary marker)
            "H2AX_pS140",        # Alternative γH2AX site
            "ATM_pS1981",        # ATM kinase activated
            "ATR_pS428",         # Replication stress active
            "CHK1_pS345",        # Checkpoint active
            "CHK1_pS296",        # Alternative CHK1 site
            "CHK2_pT68",         # CHK2 activated
            "53BP1",             # DSB repair marker
        ],
        # DDR total protein
        "ddr_total": [
            "ATM",
            "ATR",
            "CHK1",
            "CHK2",
        ],
        # DNA sensing machinery
        "sensing": [
            "cGAS",              # DNA sensor
            "STING",             # Adaptor protein
            "IRF-3_pS396",       # Activation readout
            "IRF-3",             # Total IRF3
        ],
        # Repair machinery context
        "repair_machinery": [
            "MRE11",             # MRN complex
            "RAD50",             # MRN complex
            "RAD51",             # Homologous recombination
            "PARP1",             # PARP (total)
            "PARPCLEAVED",       # Cleaved PARP (apoptosis indicator)
            "KU80",              # NHEJ
            "DNA-Ligase-IV",     # NHEJ
            "XRCC1",             # BER/SSB repair
        ],
    })
    
    # -------------------------------------------------------------------------
    # Checkpoint/evasion molecules (visibility output layer)
    # -------------------------------------------------------------------------
    checkpoint_evasion: List[str] = field(default_factory=lambda: [
        "PDL1",          # CD274 - primary checkpoint
        "CTLA4",         # Co-inhibitory receptor
        "B7-H3",         # CD276 - checkpoint
        "B7-H4",         # VTCN1 - checkpoint
        "IDO",           # Metabolic immunosuppression
        "PDCD1",         # PD-1 (mostly from infiltrating T cells)
    ])
    
    # -------------------------------------------------------------------------
    # Suppressive signaling (alternative evasion mechanisms)
    # PI3K-AKT-mTOR axis associated with immune exclusion
    # -------------------------------------------------------------------------
    suppressive_signaling: Dict[str, List[str]] = field(default_factory=lambda: {
        "pi3k_akt": [
            "PTEN",              # Loss = pathway activation
            "AKT_pS473",         # AKT activated (mTORC2 site)
            "AKT_pT308",         # AKT activated (PDK1 site)
            "PDK1_pS241",        # PDK1 activated
            "PI3KP110ALPHA",     # PI3K catalytic subunit
            "PI3KP85",           # PI3K regulatory subunit
        ],
        "mtor": [
            "MTOR_pS2448",       # mTOR activated
            "P70S6K_pT389",      # S6K1 activated (mTORC1 readout)
            "S6_pS235S236",      # S6 phosphorylated (mTORC1 readout)
            "S6_pS240S244",      # Alternative S6 sites
            "4EBP1_pT37T46",     # 4E-BP1 phosphorylated
            "4EBP1_pS65",        # 4E-BP1 alternative site
            "PRAS40_pT246",      # PRAS40 (AKT target)
        ],
        "total": [
            "AKT",
            "MTOR",
            "P70S6K1",
            "S6",
            "4EBP1",
            "PRAS40",
        ],
        # STAT3 as suppressive counter-signal
        "stat3_suppressive": [
            "STAT3_pY705",       # Activated STAT3
            "Stat3",             # Total STAT3
            "IL-6",              # STAT3 activator
        ],
    })
    
    # -------------------------------------------------------------------------
    # Immune infiltration/context markers
    # For stratification and adjustment - NOT confounders per se
    # -------------------------------------------------------------------------
    immune_context: Dict[str, List[str]] = field(default_factory=lambda: {
        "lymphocyte": [
            "CD45",          # Pan-leukocyte (PTPRC)
            "CD4",           # T-helper
            "CD20",          # B-cell
            "LCK",           # T-cell signaling
            "ZAP-70",        # T/NK signaling
        ],
        "cytolytic": [
            "Granzyme-B",    # Cytolytic activity readout
            "SYK",           # Immune cell signaling
        ],
        "myeloid_related": [
            "CD86",          # Costimulatory (APC)
            "CIITA",         # MHC-II transactivator
            "HLA-DQA1",      # MHC-II (infiltration proxy)
        ],
        "other_immune": [
            "CD31",          # Endothelial/platelets
            "CD44",          # Multiple contexts
            "CD38",          # Multiple lineages
        ],
    })
    
    # -------------------------------------------------------------------------
    # Proliferation markers (for normalization/covariate)
    # -------------------------------------------------------------------------
    proliferation: List[str] = field(default_factory=lambda: [
        "PCNA",              # S-phase marker
        "CYCLINB1",          # G2/M marker
        "CDK1",              # Cell cycle kinase
        "CDK1_pY15",         # Inhibitory phosphorylation
        "Aurora-A",          # Mitotic kinase
        "PLK1",              # Mitotic kinase
    ])
    
    # -------------------------------------------------------------------------
    # Apoptosis/cell death (for ICD context)
    # -------------------------------------------------------------------------
    apoptosis: Dict[str, List[str]] = field(default_factory=lambda: {
        "pro_apoptotic": [
            "BAX",
            "BAK",
            "BID",
            "BIM",
            "Puma",              # BBC3
            "CASPASE3",
            "CASPASE7CLEAVEDD198",
            "CASPASE8",
            "CASPASE9",
            "PARPCLEAVED",
        ],
        "anti_apoptotic": [
            "BCL2",
            "BCLXL",
            "Mcl-1",
            "XIAP",
            "CIAP",
        ],
    })
    
    def get_all_targets(self) -> Set[str]:
        """Get all unique targets across all panels."""
        targets = set()
        
        # IFN signaling
        for lst in self.ifn_signaling.values():
            targets.update(lst)
        
        # DDR-IFN
        for lst in self.ddr_to_ifn_sensing.values():
            targets.update(lst)
        
        # Checkpoint
        targets.update(self.checkpoint_evasion)
        
        # Suppressive
        for lst in self.suppressive_signaling.values():
            targets.update(lst)
        
        # Immune context
        for lst in self.immune_context.values():
            targets.update(lst)
        
        # Proliferation
        targets.update(self.proliferation)
        
        # Apoptosis
        for lst in self.apoptosis.values():
            targets.update(lst)
        
        return targets
    
    def get_panel_targets(self, panel_name: str) -> List[str]:
        """Get all targets for a specific panel (flattened)."""
        panel = getattr(self, panel_name, None)
        if panel is None:
            raise ValueError(f"Unknown panel: {panel_name}")
        
        if isinstance(panel, list):
            return panel
        elif isinstance(panel, dict):
            targets = []
            for lst in panel.values():
                targets.extend(lst)
            return list(set(targets))
        else:
            raise ValueError(f"Unexpected panel type: {type(panel)}")
    
    def check_availability(
        self,
        available_targets: Set[str],
    ) -> Dict[str, Dict[str, List[str]]]:
        """
        Check which panel targets are available in the dataset.
        
        Returns:
            Dict[panel_name] = {
                "available": [targets],
                "missing": [targets],
                "coverage": fraction
            }
        """
        report = {}
        
        panels = [
            "ifn_signaling", "ddr_to_ifn_sensing", "checkpoint_evasion",
            "suppressive_signaling", "immune_context", "proliferation",
            "apoptosis"
        ]
        
        for panel_name in panels:
            all_targets = self.get_panel_targets(panel_name)
            available = [t for t in all_targets if t in available_targets]
            missing = [t for t in all_targets if t not in available_targets]
            
            report[panel_name] = {
                "available": available,
                "missing": missing,
                "coverage": len(available) / len(all_targets) if all_targets else 0,
            }
        
        return report


# =============================================================================
# SCORE COMPUTATION
# =============================================================================

def compute_zscore_panel_score(
    expr_matrix: pd.DataFrame,
    targets: List[str],
    min_available: int = 2,
) -> pd.Series:
    """
    Compute panel score as mean of z-scored targets.
    
    Args:
        expr_matrix: samples × targets DataFrame
        targets: List of target names for this panel
        min_available: Minimum targets required to compute score
    
    Returns:
        Series of panel scores (NaN for samples with insufficient data)
    """
    available = [t for t in targets if t in expr_matrix.columns]
    
    if len(available) < min_available:
        return pd.Series(np.nan, index=expr_matrix.index)
    
    subset = expr_matrix[available].copy()
    
    # Z-score each target
    zscored = (subset - subset.mean()) / subset.std()
    
    # Mean across targets
    return zscored.mean(axis=1)


def compute_pca_panel_score(
    expr_matrix: pd.DataFrame,
    targets: List[str],
    min_available: int = 3,
    n_components: int = 1,
) -> pd.Series:
    """
    Compute panel score as first principal component.
    
    Args:
        expr_matrix: samples × targets DataFrame
        targets: List of target names for this panel
        min_available: Minimum targets required
        n_components: Number of PCs to use (1 = first PC only)
    
    Returns:
        Series of panel scores
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    available = [t for t in targets if t in expr_matrix.columns]
    
    if len(available) < min_available:
        return pd.Series(np.nan, index=expr_matrix.index)
    
    subset = expr_matrix[available].copy()
    
    # Drop samples with missing values for PCA
    clean = subset.dropna()
    
    if len(clean) < 10:
        return pd.Series(np.nan, index=expr_matrix.index)
    
    # Standardize
    scaler = StandardScaler()
    scaled = scaler.fit_transform(clean)
    
    # PCA
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(scaled)
    
    result = pd.Series(np.nan, index=expr_matrix.index)
    result.loc[clean.index] = scores[:, 0]
    
    return result


def compute_activation_ratio(
    expr_matrix: pd.DataFrame,
    phospho_target: str,
    total_target: str,
) -> pd.Series:
    """
    Compute activation ratio (phospho / total) in log space.
    
    In log-transformed expression data:
        ratio = log(phospho) - log(total) = log(phospho/total)
    
    Returns:
        Series of activation ratios
    """
    if phospho_target not in expr_matrix.columns:
        return pd.Series(np.nan, index=expr_matrix.index)
    if total_target not in expr_matrix.columns:
        return pd.Series(np.nan, index=expr_matrix.index)
    
    phospho = expr_matrix[phospho_target]
    total = expr_matrix[total_target]
    
    # Ratio in log space = difference
    return phospho - total


def compute_all_panel_scores(
    expr_matrix: pd.DataFrame,
    panels: RPPAMarkerPanels,
    method: str = "zscore",  # "zscore" or "pca"
) -> pd.DataFrame:
    """
    Compute scores for all predefined panels.
    
    Args:
        expr_matrix: samples × targets DataFrame
        panels: RPPAMarkerPanels instance
        method: Scoring method ("zscore" or "pca")
    
    Returns:
        DataFrame with panel scores as columns
    """
    score_func = (
        compute_zscore_panel_score if method == "zscore" 
        else compute_pca_panel_score
    )
    
    scores = pd.DataFrame(index=expr_matrix.index)
    
    # IFN signaling
    scores["IFN_activated"] = score_func(
        expr_matrix, panels.ifn_signaling["activated"]
    )
    
    # DDR activation
    scores["DDR_activation"] = score_func(
        expr_matrix, panels.ddr_to_ifn_sensing["ddr_activation"]
    )
    
    # cGAS-STING sensing
    scores["cGAS_STING"] = score_func(
        expr_matrix, panels.ddr_to_ifn_sensing["sensing"]
    )
    
    # DNA repair
    scores["DNA_repair"] = score_func(
        expr_matrix, panels.ddr_to_ifn_sensing["repair_machinery"]
    )
    
    # Checkpoint/evasion
    scores["checkpoint"] = score_func(
        expr_matrix, panels.checkpoint_evasion
    )
    
    # PI3K-AKT
    scores["PI3K_AKT"] = score_func(
        expr_matrix, panels.suppressive_signaling["pi3k_akt"]
    )
    
    # mTOR
    scores["mTOR"] = score_func(
        expr_matrix, panels.suppressive_signaling["mtor"]
    )
    
    # STAT3 suppressive
    scores["STAT3_suppressive"] = score_func(
        expr_matrix, panels.suppressive_signaling["stat3_suppressive"]
    )
    
    # Immune infiltration
    scores["lymphocyte_infiltration"] = score_func(
        expr_matrix, panels.immune_context["lymphocyte"]
    )
    
    # Cytolytic activity
    scores["cytolytic"] = score_func(
        expr_matrix, panels.immune_context["cytolytic"]
    )
    
    # Proliferation
    scores["proliferation"] = score_func(
        expr_matrix, panels.proliferation
    )
    
    # Apoptosis balance (pro - anti)
    pro_score = score_func(expr_matrix, panels.apoptosis["pro_apoptotic"])
    anti_score = score_func(expr_matrix, panels.apoptosis["anti_apoptotic"])
    scores["apoptosis_balance"] = pro_score - anti_score
    
    return scores


def compute_all_activation_ratios(
    expr_matrix: pd.DataFrame,
    phospho_pairs: Dict[str, Dict],
    genes_of_interest: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Compute phospho/total activation ratios for key signaling nodes.
    
    Args:
        expr_matrix: samples × targets DataFrame
        phospho_pairs: From identify_phospho_pairs()
        genes_of_interest: Specific genes to compute ratios for
    
    Returns:
        DataFrame with activation ratios as columns
    """
    if genes_of_interest is None:
        # Default to key signaling nodes
        genes_of_interest = [
            "STAT3", "AKT", "MTOR", "IRF3", "ATM", "ATR",
            "CHK1", "CHK2", "S6", "4EBP1", "ERK", "MEK1",
        ]
    
    ratios = pd.DataFrame(index=expr_matrix.index)
    
    for gene in genes_of_interest:
        if gene not in phospho_pairs:
            continue
        
        pair = phospho_pairs[gene]
        total = pair.get("total")
        phospho_list = pair.get("phospho", [])
        
        if total is None or not phospho_list:
            continue
        
        for p_target in phospho_list:
            # Extract site from target name
            site = pair.get("phospho_sites", {}).get(p_target, "")
            col_name = f"{gene}_activation_{site}" if site else f"{gene}_activation"
            
            ratios[col_name] = compute_activation_ratio(
                expr_matrix, p_target, total
            )
    
    return ratios


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_default_panels() -> RPPAMarkerPanels:
    """Get default marker panel definitions."""
    return RPPAMarkerPanels()


def summarize_panel_coverage(
    expr_matrix: pd.DataFrame,
    panels: Optional[RPPAMarkerPanels] = None,
) -> pd.DataFrame:
    """
    Create summary table of panel coverage in the dataset.
    
    Returns:
        DataFrame with panel coverage statistics
    """
    if panels is None:
        panels = get_default_panels()
    
    available = set(expr_matrix.columns)
    report = panels.check_availability(available)
    
    summary = []
    for panel_name, info in report.items():
        summary.append({
            "panel": panel_name,
            "n_available": len(info["available"]),
            "n_missing": len(info["missing"]),
            "coverage": f"{info['coverage']:.1%}",
            "missing_targets": ", ".join(info["missing"][:5]) + (
                "..." if len(info["missing"]) > 5 else ""
            ),
        })
    
    return pd.DataFrame(summary)
