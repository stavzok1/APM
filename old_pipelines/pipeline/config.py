"""
Central configuration for the regulatory pipeline.
All paths, thresholds, gene panels, and biosamples defined here.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
from pathlib import Path


# =============================================================================
# PATH CONFIGURATION
# =============================================================================

@dataclass
class PathConfig:
    """All file paths used by the pipeline."""
    
    # Working directory
    working_dir: Path = Path("/home/stavz/masters/gdc/APM/test1")
    
    # Input files
    gencode_gtf_csv: Path = field(default_factory=lambda: Path("/home/stavz/masters/gdc/APM/test/gencode.v49.annotation.gtf.csv"))
    lncrna_csv: Path = field(default_factory=lambda: Path("/home/stavz/masters/gdc/APM/test/lncRNAs_all_features.csv"))
    ccre_csv: Path = field(default_factory=lambda: Path("/home/stavz/masters/gdc/APM/test/GRCh38-cCREs.csv"))

    # TAD processed data directory
    tads_processed: Path = Path("/home/stavz/masters/gdc/APM/data/tad_processed")
    
    # SCREEN data
    screen_exp_zip: Path = Path("/home/stavz/masters/gdc/APM/data/regulatory_elements_data/SCREEN_exp/Human-Gene-Links.zip")
    screen_exp_inner_file: str = "V4-hg38.Gene-Links.3D-Chromatin.txt"
    screen_comp_gz: Path = Path("/home/stavz/masters/gdc/APM/data/regulatory_elements_data/SCREEN_comp/V4-hg38.Gene-Links.Computational-Methods.txt.gz")
    
    # ABC model
    abc_predictions: Path = Path("/home/stavz/masters/gdc/APM/data/regulatory_elements_data/ABC_output/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt")
    
    # HiChIP
    hichip_dir: Path = Path("/home/stavz/masters/gdc/APM/data/HiCHIP")
    
    # Cell line signals
    cell_lines_dir: Path = Path("/home/stavz/masters/gdc/APM/data/regulatory_elements_data/cCRE_signals")
    
    # miRNA
    targetscan_predictions: Path = Path("/home/stavz/masters/gdc/APM/data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt")
    
    def __post_init__(self):
        """Convert strings to Path objects if needed."""
        for field_name, field_value in self.__dict__.items():
            if isinstance(field_value, str) and field_name != 'screen_exp_inner_file':
                setattr(self, field_name, Path(field_value))


# =============================================================================
# GENE PANEL
# =============================================================================

PRIMARY_GENES: List[str] = [
    # Immunoproteasome
    "PSMB8", "PSMB9", "PSMB10",
    "PSME1", "PSME2", "PSME4",
    
    # Peptide transport & loading
    "TAP1", "TAP2", "TAPBP", "TAPBPL",
    "CALR", "PDIA3", "PDIA6",
    
    # MHC-I
    "B2M",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
    
    # NK ligands
    "MICA", "MICB",
    "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L",
    "PVR", "NECTIN2", "NCR3LG1",
    
    # Negative regulators
    "SOCS1", "SOCS3",
    "NFKBIA", "TNFAIP3",
    
    # JAK-STAT pathway
    "STAT1", "STAT3",
    "JAK1", "JAK2",
    
    # Transcription factors
    "IRF1", "IRF2",
    "NLRC5",
    
    # Checkpoint
    "CD274",      # PD-L1
    "PDCD1LG2",   # PD-L2
    
    # Interferon
    "IFNG",
    "IFNGR1", "IFNGR2",
    
    # Shedding
    "ADAM10", "ADAM17",
    
    # Immunosuppressive
    "TGFB1",
    "IL6", "IL6R", "IL6ST",
    
    # STAT inhibitors
    "PIAS1", "PIAS3",
    "PTPN2", "PTPN11",
    
    # Signaling
    "EGFR",
    "SRC",
    "SMAD3",
    
    # Transcriptional coactivators & chemokines
    "EP300", "CREBBP", 
    "CCL5", "CXCL9", "CXCL10", "CXCL11",
]


# =============================================================================
# BIOSAMPLE / CELL LINE PANELS
# =============================================================================

@dataclass
class BiosampleConfig:
    """Biosample panels for different data sources."""
    
    # SCREEN experimental (3D chromatin)
    screen_exp: List[str] = field(default_factory=lambda: [
        "MCF-7",
        "MCF_10A",
        "Breast Mammary Tissue",
        "breast_epithelium_tissue_female_adult_(53_years)",
        "mammary_epithelial_cell_female_adult_(19_years)",
    ])
    
    # SCREEN computational
    screen_comp: List[str] = field(default_factory=lambda: [
        "MCF-7",
        "MCF_10A",
        "T47D",
        "mammary_epithelial_cell",
        "breast_epithelium",
    ])
    
    # ABC model cell types
    abc_celltypes: List[str] = field(default_factory=lambda: [
        "MCF-7-ENCODE",
        "MCF10A-Ji2017",
        "MCF10A_treated_with_TAM24hr-Ji2017",
        "MDA-MB-231",
        "mammary_epithelial_cell-Roadmap",
        "breast_epithelium-ENCODE",
    ])
    
    # HiChIP cell lines
    hichip_panel: List[str] = field(default_factory=lambda: [
        "B80T5", "HMEC", "K5plusK19plus", "MCF7", "MDA-MB-231", "T47D", "ZR751"
    ])
    
    # Cell lines for cCRE signal integration
    ccre_signal_cell_lines: List[str] = field(default_factory=lambda: [
        "MCF7", "HMEC1"
    ])


# =============================================================================
# THRESHOLD CONFIGURATION
# =============================================================================

@dataclass
class ThresholdConfig:
    """All numeric thresholds used in the pipeline."""
    
    # Distance windows
    lncrna_window_bp: int = 1_000_000
    ccre_window_bp: int = 1_000_000
    
    # Promoter definition (relative to TSS)
    promoter_upstream_bp: int = 2000
    promoter_downstream_bp: int = 500
    
    # Distance tiers for cCRE-gene matching
    tier_edges_bp: List[int] = field(default_factory=lambda: [0, 100_000, 250_000, 500_000, 1_000_000])
    tier_labels: List[str] = field(default_factory=lambda: ["0–100kb", "100–250kb", "250–500kb", "500–1000kb"])
    
    # SCREEN strength quantiles
    screen_weak_quantile: float = 0.90
    screen_strong_quantile: float = 0.95
    intact_hic_pvalue_threshold: float = 1e-6
    
    # ABC thresholds
    abc_present_threshold: float = 0.015
    abc_strong_threshold: float = 0.05
    
    # miRNA thresholds
    mirna_top_n: int = 200
    mirna_weight_threshold: float = -0.2
    
    # Streaming chunk sizes
    screen_chunksize: int = 300_000
    abc_chunksize: int = 500_000
    
    # Binning for fast prefilter
    coarse_bin_size: int = 100_000


# =============================================================================
# SIGNAL CONFIGURATION
# =============================================================================

# Signals to extract from cell line directories
WANTED_SIGNALS: List[str] = ["H3K27ac", "H3K4me3", "CTCF", "DNase"]

# HiChIP experiment type
HICHIP_EXPERIMENT: str = "H3K27ac"

# Column index for cCRE ID in BED files
CCRE_ID_COL_IN_BED: int = 3


# =============================================================================
# OUTPUT SUBDIRECTORIES
# =============================================================================

OUTPUT_SUBDIRS: Dict[str, str] = {
    "lncrna_matching": "lncRNA_matching",
    "regulatory_elements": "regulatory_elements_matching",
    "mirna": "miRNA",
    "evidence": "evidence",
}


# =============================================================================
# CONVENIENCE: DEFAULT CONFIGS
# =============================================================================

def get_default_config():
    """Return default configuration objects."""
    return {
        "paths": PathConfig(),
        "biosamples": BiosampleConfig(),
        "thresholds": ThresholdConfig(),
    }


# For quick imports
PATHS = PathConfig()
BIOSAMPLES = BiosampleConfig()
THRESHOLDS = ThresholdConfig()
