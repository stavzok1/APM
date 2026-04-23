"""
TAD sources configuration.

Defines available TAD biosamples, their metadata (PAM50 subtype, source study),
and provides auto-discovery of processed TAD directories.

This integrates with the main pipeline config but is kept separate for clarity
since TAD sources are numerous and have their own metadata structure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import os


# =============================================================================
# BIOSAMPLE METADATA
# =============================================================================

@dataclass
class TADBiosample:
    """Metadata for a single TAD biosample."""
    name: str                          # Folder name (e.g., "Kim_T47D")
    study: str                         # Source study (e.g., "Kim", "Rao")
    cell_line: str                     # Cell line or sample name
    pam50: Optional[str] = None        # PAM50 subtype if known
    tissue_type: str = "cell_line"     # "cell_line", "tumor_tissue", "normal_tissue", "healthy_breast"
    genome_build: str = "hg38"         # After liftover, should all be hg38
    notes: Optional[str] = None


# Complete biosample registry with metadata
TAD_BIOSAMPLE_REGISTRY: Dict[str, TADBiosample] = {
    # Kim et al. (GSE167150) - breast cancer cell lines and tissues
    "Kim_T47D": TADBiosample(
        name="Kim_T47D", study="Kim", cell_line="T47D",
        pam50="LumA", tissue_type="cell_line"
    ),
    "Kim_HMEC": TADBiosample(
        name="Kim_HMEC", study="Kim", cell_line="HMEC",
        pam50=None, tissue_type="normal_tissue",
        notes="Human mammary epithelial cells (normal)"
    ),
    "Kim_HCC70": TADBiosample(
        name="Kim_HCC70", study="Kim", cell_line="HCC70",
        pam50="Basal", tissue_type="cell_line"
    ),
    "Kim_HCC1954": TADBiosample(
        name="Kim_HCC1954", study="Kim", cell_line="HCC1954",
        pam50="HER2", tissue_type="cell_line"
    ),
    "Kim_BT549": TADBiosample(
        name="Kim_BT549", study="Kim", cell_line="BT549",
        pam50="Basal", tissue_type="cell_line"
    ),
    "Kim_ZR7530": TADBiosample(
        name="Kim_ZR7530", study="Kim", cell_line="ZR75-30",
        pam50="LumA", tissue_type="cell_line"
    ),
    "Kim_normal_tissue": TADBiosample(
        name="Kim_normal_tissue", study="Kim", cell_line="normal_tissue",
        pam50=None, tissue_type="normal_tissue",
        notes="Normal breast tissue sample"
    ),
    "Kim_TNBC_tissue1": TADBiosample(
        name="Kim_TNBC_tissue1", study="Kim", cell_line="TNBC_tissue1",
        pam50="Basal", tissue_type="tumor_tissue"
    ),
    "Kim_TNBC_tissue2": TADBiosample(
        name="Kim_TNBC_tissue2", study="Kim", cell_line="TNBC_tissue2",
        pam50="Basal", tissue_type="tumor_tissue"
    ),
    "Kim_TNBC_tissue3": TADBiosample(
        name="Kim_TNBC_tissue3", study="Kim", cell_line="TNBC_tissue3",
        pam50="Basal", tissue_type="tumor_tissue"
    ),
    
    # Rao et al. (GSE63525) - reference HMEC
    "Rao_HMEC": TADBiosample(
        name="Rao_HMEC", study="Rao", cell_line="HMEC",
        pam50=None, tissue_type="normal_tissue",
        notes="High-resolution Arrowhead domains"
    ),
    
    # Le Dily et al. (GSE109229) - breast cancer cell lines
    "LeDily_BT474": TADBiosample(
        name="LeDily_BT474", study="LeDily", cell_line="BT474",
        pam50="LumB", tissue_type="cell_line"
    ),
    "LeDily_MCF10A": TADBiosample(
        name="LeDily_MCF10A", study="LeDily", cell_line="MCF10A",
        pam50=None, tissue_type="cell_line",
        notes="Non-tumorigenic epithelial"
    ),
    "LeDily_SKBR3": TADBiosample(
        name="LeDily_SKBR3", study="LeDily", cell_line="SKBR3",
        pam50="HER2", tissue_type="cell_line"
    ),
    
    # Golloshi et al. (GSE143678) - MDA-MB-231 with different insulation thresholds
    "Golloshi_MDA231_control": TADBiosample(
        name="Golloshi_MDA231_control", study="Golloshi", cell_line="MDA-MB-231",
        pam50="Basal", tissue_type="cell_line",
        notes="Control insulation parameters"
    ),
    "Golloshi_MDA231_top10": TADBiosample(
        name="Golloshi_MDA231_top10", study="Golloshi", cell_line="MDA-MB-231",
        pam50="Basal", tissue_type="cell_line",
        notes="Top 10% insulation threshold"
    ),
    "Golloshi_MDA231_bottom10": TADBiosample(
        name="Golloshi_MDA231_bottom10", study="Golloshi", cell_line="MDA-MB-231",
        pam50="Basal", tissue_type="cell_line",
        notes="Bottom 10% insulation threshold"
    ),
    
    # van den Brand et al. (GSE273999) - healthy and primary breast
    "vandenBrand_HB1": TADBiosample(
        name="vandenBrand_HB1", study="vandenBrand", cell_line="HB1",
        pam50=None, tissue_type="healthy_breast",
        notes="Healthy breast tissue"
    ),
    "vandenBrand_HB2": TADBiosample(
        name="vandenBrand_HB2", study="vandenBrand", cell_line="HB2",
        pam50=None, tissue_type="healthy_breast"
    ),
    "vandenBrand_PB1": TADBiosample(
        name="vandenBrand_PB1", study="vandenBrand", cell_line="PB1",
        pam50=None, tissue_type="tumor_tissue",
        notes="Primary breast cancer"
    ),
    "vandenBrand_PB2": TADBiosample(
        name="vandenBrand_PB2", study="vandenBrand", cell_line="PB2",
        pam50=None, tissue_type="tumor_tissue"
    ),
    "vandenBrand_PB3": TADBiosample(
        name="vandenBrand_PB3", study="vandenBrand", cell_line="PB3",
        pam50=None, tissue_type="tumor_tissue"
    ),
    "vandenBrand_PB4": TADBiosample(
        name="vandenBrand_PB4", study="vandenBrand", cell_line="PB4",
        pam50=None, tissue_type="tumor_tissue"
    ),
    "vandenBrand_PB5": TADBiosample(
        name="vandenBrand_PB5", study="vandenBrand", cell_line="PB5",
        pam50=None, tissue_type="tumor_tissue"
    ),
}


# =============================================================================
# PAM50 SUBTYPE GROUPINGS
# =============================================================================

def get_biosamples_by_pam50(pam50: str) -> List[str]:
    """Get all biosamples matching a PAM50 subtype."""
    return [
        name for name, meta in TAD_BIOSAMPLE_REGISTRY.items()
        if meta.pam50 == pam50
    ]


def get_biosamples_by_tissue_type(tissue_type: str) -> List[str]:
    """Get all biosamples matching a tissue type."""
    return [
        name for name, meta in TAD_BIOSAMPLE_REGISTRY.items()
        if meta.tissue_type == tissue_type
    ]


def get_biosamples_by_study(study: str) -> List[str]:
    """Get all biosamples from a study."""
    return [
        name for name, meta in TAD_BIOSAMPLE_REGISTRY.items()
        if meta.study == study
    ]


# Convenience groupings
PAM50_GROUPS = {
    "Basal": get_biosamples_by_pam50("Basal"),
    "LumA": get_biosamples_by_pam50("LumA"),
    "LumB": get_biosamples_by_pam50("LumB"),
    "HER2": get_biosamples_by_pam50("HER2"),
}

TISSUE_GROUPS = {
    "cell_line": get_biosamples_by_tissue_type("cell_line"),
    "tumor_tissue": get_biosamples_by_tissue_type("tumor_tissue"),
    "normal_tissue": get_biosamples_by_tissue_type("normal_tissue"),
    "healthy_breast": get_biosamples_by_tissue_type("healthy_breast"),
}


# =============================================================================
# TAD SOURCES DISCOVERY
# =============================================================================

@dataclass
class TADSourcePaths:
    """Paths to TAD files for a single biosample."""
    biosample: str
    base_dir: Path
    domains: Path = field(init=False)
    boundaries: Path = field(init=False)
    flanks: Path = field(init=False)
    
    def __post_init__(self):
        self.domains = self.base_dir / "domains.parquet"
        self.boundaries = self.base_dir / "boundaries.parquet"
        self.flanks = self.base_dir / "flanks.parquet"
    
    def exists(self) -> bool:
        """Check if all required files exist."""
        return (
            self.domains.exists() and
            self.boundaries.exists() and
            self.flanks.exists()
        )
    
    def exists_any_format(self) -> Tuple[bool, str]:
        """Check if files exist in any format (parquet, csv, tsv)."""
        for fmt in ["parquet", "csv", "tsv"]:
            domains_path = self.base_dir / f"domains.{fmt}"
            if domains_path.exists():
                return True, fmt
        return False, ""


def discover_tad_sources(
    processed_dir: Path,
    required_files: bool = True,
) -> Dict[str, TADSourcePaths]:
    """
    Auto-discover TAD sources from processed directory.
    
    Args:
        processed_dir: Path to TADs/processed/ directory
        required_files: If True, only return sources with all files present
    
    Returns:
        Dict mapping biosample name → TADSourcePaths
    """
    sources = {}
    
    if not processed_dir.exists():
        return sources
    
    for subdir in processed_dir.iterdir():
        if not subdir.is_dir():
            continue
        
        biosample = subdir.name
        paths = TADSourcePaths(biosample=biosample, base_dir=subdir)
        
        if required_files:
            if paths.exists():
                sources[biosample] = paths
            else:
                # Check for CSV/TSV fallback
                exists, fmt = paths.exists_any_format()
                if exists:
                    # Update paths to use discovered format
                    paths.domains = subdir / f"domains.{fmt}"
                    paths.boundaries = subdir / f"boundaries.{fmt}"
                    paths.flanks = subdir / f"flanks.{fmt}"
                    sources[biosample] = paths
        else:
            sources[biosample] = paths
    
    return sources


def get_tad_sources_for_pam50(
    processed_dir: Path,
    pam50: str,
) -> Dict[str, TADSourcePaths]:
    """Get TAD sources matching a PAM50 subtype."""
    all_sources = discover_tad_sources(processed_dir)
    matching = get_biosamples_by_pam50(pam50)
    return {k: v for k, v in all_sources.items() if k in matching}


def get_tad_sources_for_study(
    processed_dir: Path,
    study: str,
) -> Dict[str, TADSourcePaths]:
    """Get TAD sources from a specific study."""
    all_sources = discover_tad_sources(processed_dir)
    matching = get_biosamples_by_study(study)
    return {k: v for k, v in all_sources.items() if k in matching}


# =============================================================================
# DEFAULT SELECTIONS
# =============================================================================

# Recommended TAD sources for different use cases
DEFAULT_CELL_LINE_SOURCES = [
    "Kim_T47D",       # LumA reference
    "Kim_HCC70",      # Basal reference
    "Kim_HCC1954",    # HER2 reference
    "LeDily_BT474",   # LumB reference
]

DEFAULT_NORMAL_SOURCES = [
    "Rao_HMEC",       # High-quality normal reference
    "Kim_HMEC",       # Kim study normal
]

# For PAM50-matched analysis
PAM50_REFERENCE_MAP = {
    "Basal": ["Kim_HCC70", "Kim_BT549", "Golloshi_MDA231_control"],
    "LumA": ["Kim_T47D", "Kim_ZR7530"],
    "LumB": ["LeDily_BT474"],
    "HER2": ["Kim_HCC1954", "LeDily_SKBR3"],
    "Normal": ["Rao_HMEC", "Kim_HMEC", "LeDily_MCF10A"],
}


def get_best_tad_source_for_pam50(
    pam50: str,
    available_sources: Dict[str, TADSourcePaths],
) -> Optional[str]:
    """
    Get the best available TAD source for a PAM50 subtype.
    
    Returns first available from the reference map, or None if none available.
    """
    if pam50 not in PAM50_REFERENCE_MAP:
        return None
    
    for candidate in PAM50_REFERENCE_MAP[pam50]:
        if candidate in available_sources:
            return candidate
    
    return None
