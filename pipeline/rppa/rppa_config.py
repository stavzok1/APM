"""
RPPA configuration settings.

Extends the main pipeline config with RPPA-specific paths,
thresholds, and marker panel definitions.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set
from pathlib import Path
from ..config  import PATHS

# =============================================================================
# RPPA PATH CONFIGURATION
# =============================================================================

@dataclass
class RPPAPathConfig:
    """RPPA-specific file paths."""
    
    # Defaults mirror pipeline.config.PathConfig (rppa_* fields)
    sample_data_dir: Path = Path(f"{PATHS.working_dir}/rppa/samples")

    annotation_file: Path = Path(f"{PATHS.annotations_dir}/rppa/antibody_annotation.csv")

    sample_metadata_tumor: Optional[Path] = Path(f"{PATHS.annotations_dir}/rppa/samples.tsv")

    output_dir: Path = Path(f"{PATHS.working_dir}/rppa/processed")
    
    def __post_init__(self):
        """Convert strings to Path objects if needed."""
        for field_name, field_value in self.__dict__.items():
            if isinstance(field_value, str):
                setattr(self, field_name, Path(field_value))


# =============================================================================
# RPPA THRESHOLD CONFIGURATION
# =============================================================================

@dataclass
class RPPAThresholdConfig:
    """RPPA-specific thresholds."""
    
    # Validation status inclusion
    include_caution: bool = True
    include_under_evaluation: bool = False
    
    # QC thresholds
    min_samples_per_target: int = 10
    min_targets_per_sample: int = 50
    
    # Discordance detection
    discordance_zscore_threshold: float = 1.0
    min_discordant_samples: int = 5
    
    # Quantile thresholds for signaling block detection
    quantile_high: float = 0.75
    quantile_low: float = 0.25
    
    # Panel score computation
    min_markers_for_panel: int = 2
    
    # Correlation significance
    correlation_pvalue_threshold: float = 0.05


# =============================================================================
# IMMUNE VISIBILITY MARKERS
# =============================================================================

@dataclass
class ImmuneVisibilityMarkers:
    """
    Markers directly relevant to antigen presentation and immune visibility.
    
    Note: Many APM genes (TAP1/2, B2M, HLA, PSMB8/9/10) are NOT in RPPA.
    This defines available proxies and related pathway members.
    """
    
    # Trans-regulators of APM (proxies for APM activity)
    trans_regulators: List[str] = field(default_factory=lambda: [
        "IRF1",          # Critical: downstream of STAT1, induces APM
        "IRF-3",         # Innate IFN pathway
        "IRF-3_pS396",   # Activated IRF3
        "CIITA",         # MHC-II master regulator (infiltration proxy)
    ])
    
    # JAK-STAT pathway (upstream of APM)
    jak_stat: List[str] = field(default_factory=lambda: [
        # Note: STAT1 is MISSING from RPPA
        "Stat3",         # Counter-regulatory
        "STAT3_pY705",   # Activated STAT3 (suppressive)
        "STAT5ALPHA",    # Alternative STAT
        "JAK2",          # Upstream kinase
    ])
    
    # Checkpoint molecules (direct visibility modulators)
    checkpoint: List[str] = field(default_factory=lambda: [
        "PDL1",          # CD274 - main checkpoint
        "CTLA4",         # Co-inhibitory
        "B7-H3",         # CD276
        "B7-H4",         # VTCN1
        "IDO",           # Metabolic suppression
    ])
    
    # DNA sensing (links DDR to immune)
    dna_sensing: List[str] = field(default_factory=lambda: [
        "cGAS",          # DNA sensor
        "STING",         # Adaptor
    ])
    
    # Cytolytic markers (functional outcome)
    cytolytic: List[str] = field(default_factory=lambda: [
        "Granzyme-B",    # Direct cytolytic readout
    ])
    
    # Missing from RPPA but critical for APM
    missing_critical: List[str] = field(default_factory=lambda: [
        "STAT1",         # Master IFN regulator - use IRF1 as proxy
        "NLRC5",         # MHC-I master regulator
        "TAP1", "TAP2",  # Peptide transport
        "B2M",           # MHC-I component
        "HLA-A", "HLA-B", "HLA-C",  # MHC-I
        "PSMB8", "PSMB9", "PSMB10",  # Immunoproteasome
    ])


# =============================================================================
# RPPA TO GENE MAPPING OVERRIDES
# =============================================================================

# Some RPPA targets have non-obvious gene name mappings
RPPA_TARGET_TO_GENE_OVERRIDES: Dict[str, str] = {
    # Phospho forms → base gene
    "PDL1": "CD274",
    "Granzyme-B": "GZMB",
    "b-Actin": "ACTB",
    "b-Catenin_pT41_S45": "CTNNB1",
    "BETACATENIN": "CTNNB1",
    "c-Abl_pY412": "ABL1",
    "CABL": "ABL1",
    "CD45": "PTPRC",
    "CD20": "MS4A1",
    "CD31": "PECAM1",
    "CD29": "ITGB1",
    "CD49B": "ITGA2",
    "CMET": "MET",
    "CMYC": "MYC",
    "CRAF": "RAF1",
    "Cyclin-D3": "CCND3",
    "ECADHERIN": "CDH1",
    "ERK2": "MAPK1",
    "Erk5": "MAPK7",
    "FAK": "PTK2",
    "GSK-3B": "GSK3B",
    "Grp75": "HSPA9",
    "HSP27": "HSPB1",
    "HSP60": "HSPD1",
    "HSP70": "HSPA1A",
    "IRF-3": "IRF3",
    "IRF-3_pS396": "IRF3",
    "IL-6": "IL6",
    "NCADHERIN": "CDH2",
    "PCADHERIN": "CDH3",
    "P53": "TP53",
    "P21": "CDKN1A",
    "P27": "CDKN1B",
    "P16INK4A": "CDKN2A",
    "P63": "TP63",
    "Stat3": "STAT3",
    "Tau": "MAPT",
    "Twist": "TWIST1",
    "VEGFR2": "KDR",
    "ZAP-70": "ZAP70",
}


# =============================================================================
# KEY PHOSPHO PAIRS FOR ACTIVATION RATIOS
# =============================================================================

# Pre-defined phospho pairs of high biological relevance
KEY_PHOSPHO_PAIRS: Dict[str, Dict[str, str]] = {
    # Gene: {total: target, phospho: target, site: description}
    "STAT3": {
        "total": "Stat3",
        "phospho": "STAT3_pY705",
        "site": "Y705",
        "interpretation": "High ratio = STAT3 pathway active (often immunosuppressive)",
    },
    "AKT": {
        "total": "AKT",
        "phospho": "AKT_pS473",  # Also AKT_pT308
        "site": "S473",
        "interpretation": "High ratio = PI3K pathway active",
    },
    "IRF3": {
        "total": "IRF-3",
        "phospho": "IRF-3_pS396",
        "site": "S396",
        "interpretation": "High ratio = innate sensing active",
    },
    "ATM": {
        "total": "ATM",
        "phospho": "ATM_pS1981",
        "site": "S1981",
        "interpretation": "High ratio = DDR active",
    },
    "ATR": {
        "total": "ATR",
        "phospho": "ATR_pS428",
        "site": "S428",
        "interpretation": "High ratio = replication stress active",
    },
    "CHK1": {
        "total": "CHK1",
        "phospho": "CHK1_pS345",
        "site": "S345",
        "interpretation": "High ratio = checkpoint engaged",
    },
    "CHK2": {
        "total": "CHK2",
        "phospho": "CHK2_pT68",
        "site": "T68",
        "interpretation": "High ratio = DDR checkpoint active",
    },
    "MTOR": {
        "total": "MTOR",
        "phospho": "MTOR_pS2448",
        "site": "S2448",
        "interpretation": "High ratio = mTOR pathway active",
    },
    "S6": {
        "total": "S6",
        "phospho": "S6_pS235S236",
        "site": "S235S236",
        "interpretation": "High ratio = mTORC1 active (translation)",
    },
}


# =============================================================================
# DEFAULT CONFIGURATION INSTANCES
# =============================================================================

RPPA_PATHS = RPPAPathConfig()
RPPA_THRESHOLDS = RPPAThresholdConfig()
IMMUNE_VISIBILITY_MARKERS = ImmuneVisibilityMarkers()


def get_rppa_config() -> Dict[str, any]:
    """Get all RPPA configuration objects."""
    return {
        "paths": RPPA_PATHS,
        "thresholds": RPPA_THRESHOLDS,
        "markers": IMMUNE_VISIBILITY_MARKERS,
        "target_overrides": RPPA_TARGET_TO_GENE_OVERRIDES,
        "key_phospho_pairs": KEY_PHOSPHO_PAIRS,
    }
