"""
Central configuration for the regulatory pipeline.
All paths, thresholds, gene panels, and biosamples defined here.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def _repo_root() -> Path:
    """Project root (parent of ``pipeline/``)."""
    return Path(__file__).resolve().parent.parent


def _env_path(name: str) -> Optional[Path]:
    v = os.environ.get(name, "").strip()
    return Path(v).expanduser() if v else None


def _default_working_dir() -> Path:
    """
    Large inputs and pipeline outputs under ``<repo>/data`` by default.

    Override on any machine (cluster login node, HPC scratch, etc.) with::

        export APM_WORKING_DIR=/path/to/your/APM/data
    """
    return _env_path("APM_WORKING_DIR") or (_repo_root() / "data")


def _default_annotations_dir() -> Path:
    """Annotation tables under ``<repo>/annotations`` unless ``APM_ANNOTATIONS_DIR`` is set."""
    return _env_path("APM_ANNOTATIONS_DIR") or (_repo_root() / "annotations")


def _default_vep_cache_dir() -> Path:
    return _env_path("APM_VEP_CACHE_DIR") or (Path.home() / ".vep")


def _default_sv_reference_fasta() -> Path:
    if (p := _env_path("APM_SV_REFERENCE_FASTA")) is not None:
        return p
    return Path.home() / ".vep" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"


# Captured once at import so all derived defaults stay consistent.
_APM_WDIR = _default_working_dir()
_APM_ADIR = _default_annotations_dir()


# =============================================================================
# PATH CONFIGURATION
# =============================================================================

@dataclass
class PathConfig:
    """All file paths used by the pipeline.

    **Machine-local roots** (no need to edit this file when moving hosts):

    - ``APM_WORKING_DIR``: data tree (default ``<repo>/data``).
    - ``APM_ANNOTATIONS_DIR``: annotation tables (default ``<repo>/annotations``).
    - ``APM_VEP_CACHE_DIR``: VEP cache (default ``~/.vep``).
    - ``APM_SV_REFERENCE_FASTA``: GRCh38 FASTA for SV steps (default ``~/.vep/Homo_sapiens...``).
    """

    # Working directory
    working_dir: Path = _APM_WDIR
    annotations_dir: Path = _APM_ADIR
    
    # External Input files
    gencode_gtf_pq: Path = Path(f"{_APM_WDIR}/gencode.v49.slim.parquet")
    # Full multi-feature annotation (recommended runtime source for exon_id/exon_number metadata)
    gencode_gtf_full_pq: Path = Path(f"{_APM_WDIR}/gencode.v49.annotation.parquet")
    # lncRNA–gene matching *output* (lncRNA-centric, one row per lncRNA gene + nearby-gene lists).
    # Consumed by Methylation / SV when they need matched context; **not** the upstream lncRNA
    # universe for ``run_full_pipeline`` (that is derived from ``gencode_gtf_pq`` lncRNA rows).
    lncrnas_genes_centric: Path = Path(f"{_APM_WDIR}/lncRNA_matching/lncRNAs_with_genes_1000000bp.csv")
    # Back-compat alias for ``lncrnas_genes_centric`` (SV default path; same matching artifact).
    lncrna_csv: Path = Path(f"{_APM_WDIR}/lncRNA_matching/lncRNAs_with_genes_1000000bp.csv")
    ccre_csv: Path = Path(f"{_APM_WDIR}/GRCh38-cCREs.csv")
    
    # Clinical annotations (BRCA)
    brca_clinical: Path = Path(f"{_APM_ADIR}/BRCA_clinical")

    # HGNC-style gene-symbol alias table (tab header, comma-separated rows; mis-suffixed as .xls)
    hgnc_alias_table: Path = Path(f"{_APM_ADIR}/Alias_v5.22.xls")


    # TAD processed data directory
    tads_processed: Path = Path(f"{_APM_WDIR}/TADs")
    
    # SCREEN data
    screen_exp_zip: Path = Path(f"{_APM_WDIR}/regulatory_elements_data/SCREEN_exp/Human-Gene-Links.zip")
    screen_exp_inner_file: str = "V4-hg38.Gene-Links.3D-Chromatin.txt"
    screen_comp_gz: Path = Path(f"{_APM_WDIR}/regulatory_elements_data/SCREEN_comp/V4-hg38.Gene-Links.Computational-Methods.txt.gz")
    
    # ABC model
    abc_predictions: Path = Path(f"{_APM_WDIR}/regulatory_elements_data/ABC_output/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt")
    
    # HiChIP
    hichip_dir: Path = Path(f"{_APM_WDIR}/HiCHIP")
    
    # Cell line signals
    cell_lines_dir: Path = Path(f"{_APM_WDIR}/regulatory_elements_data/cCRE_signals")
    
    # miRNA
    targetscan_predictions: Path = Path(f"{_APM_WDIR}/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt") # Computational
    mirtarbase_csv: Path = Path(f"{_APM_WDIR}/miRNA/mirtar.csv") # Experimental
    mir_family_info: Path = Path(f"{_APM_WDIR}/miRNA/miR_Family_Info.txt")    
    mirna_gff: Path = Path(f"{_APM_WDIR}/miRNA/hsa.gff")
    # miRNA loci tables for overlap (built from mirna_gff)
    mirna_precursor_loci_csv: Path = Path(f"{_APM_WDIR}/miRNA/cnv_miRNA.csv")  # pre-miRNA hairpin loci (+ metadata)
    mirna_mature_loci_csv: Path = Path(f"{_APM_WDIR}/miRNA/mirna_mature_loci.csv")  # mature arms (MIMAT) with coords
    
    # ATAC peaks
    atac_peaks_csv: Path = Path(f"{_APM_WDIR}/TCGA_ATAC/TCGA-ATAC_PanCancer_PeakSet.txt")

    
    # SNV/Mutect2 data (GDC manifest names match raw under vcfs_snv; batch loader uses VEP subset under vep_vcfs)
    snv_raw_vcf_dir: Path = Path(f"{_APM_WDIR}/SNV/vcfs_snv")
    snv_vcf_dir: Path = Path(f"{_APM_WDIR}/SNV/vep_vcfs")

    # SV pipeline paths (inputs: somatic SV VCFs under sv_vcf_dir)
    sv_vcf_dir: Path = Path(f"{_APM_WDIR}/SV/raw_somatic_sv")
    sv_output_root: Path = Path(f"{_APM_WDIR}/SV/pipeline_output")
    sv_reference_fasta: Path = field(default_factory=_default_sv_reference_fasta)
    sv_meme_file: Path = Path(f"{_APM_WDIR}/SV/motifs/core_ifn_stat_ctcf_ap1_nlrc5.plus_custom.meme")
    # GDC-style manifest (same pattern as annotations/SNV/samples.tsv)
    sv_samples_tsv: Path = Path(f"{_APM_ADIR}/SV/samples.tsv")
    all_ccre_fimo_tsv: Path = Path(f"{_APM_WDIR}/SV/all_ccre_fimo.tsv")
    vep_cache_dir: Path = field(default_factory=_default_vep_cache_dir)

    # CNV pipeline paths
    cnv_dir: Path = Path(f"{_APM_WDIR}/CNV_TCGA/CNV_extracted")
    cnv_genes: Path = Path(f"{_APM_WDIR}/cnv_genes.csv")
    # Back-compat alias used across modules; prefer mirna_mature_loci_csv when you need arm resolution by coordinates.
    mirna_path: Path = Path(f"{_APM_WDIR}/miRNA/cnv_miRNA.csv")
    cnv_annotations_path: Path = Path(f"{_APM_ADIR}/CNV/samples.tsv")
    cnv_output_dir: Path = Path(f"{_APM_WDIR}/CNV_TCGA/CNV_annotated")

    # Methylation bundle (reference / per_sample / cohort live under this tree)
    methylation_probe_reference: Path = Path(f"{_APM_ADIR}/Methylation/probe_reference.tsv")
    methylation_samples_dir: Path = Path(f"{_APM_WDIR}/Methylation/raw_samples/Tumor/")
    methylation_sample_manifest: Path = Path(f"{_APM_ADIR}/Methylation/sample_manifest.tsv")

    methylation_output_dir: Path = Path(f"{_APM_WDIR}/Methylation/")

    # RPPA (wide inputs; cohort-level outputs)
    rppa_samples_dir: Path = Path(f"{_APM_WDIR}/rppa/samples")
    rppa_processed_dir: Path = Path(f"{_APM_WDIR}/rppa/processed")
    rppa_antibody_annotation_csv: Path = Path(f"{_APM_ADIR}/rppa/TCGA_antibodies_descriptions.gencode.v36.tsv")

    # RNA expression
    rna_expression_raw : Path = Path(f"{_APM_WDIR}/RNAexp_TCGA/TCGA-BRCA.star_tpm_mapped.tsv") #TODO: figure out tpm location
    rna_expression: Path = Path(f"{_APM_WDIR}/RNAexp_TCGA/TCGA-BRCA.star_tpm_processed.tsv")
    rna_gene_col: str = "gene_symbol"
    rna_gencode_gene_probemap: Path = Path(
        f"{_APM_ADIR}/RNA/gencode.v36.annotation.gtf.gene.probemap"
    )

    # miRNA expression (wide matrix; Xena arm-specific; rows are MIMAT accessions; columns are TCGA barcodes)
    mirna_expression_tsv: Path = Path(f"{_APM_WDIR}/miRNA/XENA_mirna_arm_specific.tsv")
    # Legacy TCGA miRNA matrix previously used (kept for reference / back-compat in ad-hoc analysis)
    mirna_expression_legacy_tsv: Path = Path(f"{_APM_WDIR}/miRNA/TCGA-BRCA.mirna.tsv")
    
    # immune subtype annotations
    immune_subtype_annotations: Path = Path(f"{_APM_ADIR}/BRCA_immune_subtypes_advanced.tsv")
    thornsson_immune_table: Path = Path(f"{_APM_ADIR}/Thornsson_immune_table.tsv")
    
    # Unified clinical + immune annotations
    brca_clinical_immune_unified: Path = Path(f"{_APM_ADIR}/BRCA_clinical_immune_unified.tsv")
    
    # Normalized external annotation outputs (generated under annotations/_normalized)
    immune_subtype_annotations_normalized: Path = Path(
        f"{_APM_ADIR}/_normalized/BRCA_immune_subtypes_advanced.normalized.tsv"
    )
    thornsson_immune_table_normalized: Path = Path(
        f"{_APM_ADIR}/_normalized/Thornsson_immune_table.normalized.tsv"
    )
    hla_types_normalized: Path = Path(
        f"{_APM_ADIR}/_normalized/HLA_types.normalized.tsv"
    )
    rna_expression_sample_metadata: Path = Path(
        f"{_APM_ADIR}/_normalized/RNA_expression.sample_metadata.tsv"
    )
    atac_case_level_sample_metadata: Path = Path(
        f"{_APM_ADIR}/_normalized/ATAC_case_level.sample_metadata.tsv"
    )

    # Wide case-level matrix (comma-separated; leading columns are peak coords, then TCGA barcodes)
    atac_case_level_matrix: Path = Path(
        f"{_APM_WDIR}/TCGA_ATAC/TCGA_NORMAL_LOG_CPM_QN_BRCA_case_level.csv"
    )

    # Raw HLA types (comma-separated; includes aliquot_id)
    hla_types_tsv: Path = Path(f"{_APM_ADIR}/HLA_types.tsv")

    # GDC-style manifests derived from matrix columns (see scripts/annotations/build_annotation_sample_manifests.py)
    rna_samples_tsv: Path = Path(f"{_APM_ADIR}/RNA/samples.tsv")
    atac_samples_tsv: Path = Path(f"{_APM_ADIR}/ATAC/samples.tsv")
    hla_samples_tsv: Path = Path(f"{_APM_ADIR}/HLA/samples.tsv")
    hichip_samples_tsv: Path = Path(f"{_APM_ADIR}/HiCHIP/samples.tsv")
    mirna_samples_tsv: Path = Path(f"{_APM_ADIR}/miRNA/samples.tsv")

    # HiChIP TCGA-wide processed matrix (columns are TCGA participants, starting at 7th col)
    hichip_tcga_processed_csv: Path = Path(f"{_APM_WDIR}/HiCHIP/TCGA_BRCA_PROCESSED.csv")

    # Output Paths
  
    # gene tables
    # ``primary_*`` paths are legacy names: rows match ``PIPELINE_GENE_PANEL`` (full integration
    # when ``APM_USE_EXTENDED_GENE_PANEL`` is on), not the frozen 66 ``PRIMARY_GENES`` alone.
    genes_all_features: Path = Path(f"{_APM_WDIR}/primary_genes_all_features.csv")
    genes_only: Path = Path(f"{_APM_WDIR}/primary_genes_only.csv")

    # Tier-specific gene tables (gene-level + multifeature GENCODE rows), written by ``save_gene_tables``.
    tier2_medium_genes_only: Path = Path(f"{_APM_WDIR}/tier2_medium_genes_only.csv")
    tier2_medium_genes_all_features: Path = Path(f"{_APM_WDIR}/tier2_medium_genes_all_features.csv")
    tier3_cnv_only_genes_only: Path = Path(f"{_APM_WDIR}/tier3_cnv_only_genes_only.csv")
    tier3_cnv_only_genes_all_features: Path = Path(f"{_APM_WDIR}/tier3_cnv_only_genes_all_features.csv")
    tier4_readout_genes_only: Path = Path(f"{_APM_WDIR}/tier4_readout_genes_only.csv")
    tier4_readout_genes_all_features: Path = Path(f"{_APM_WDIR}/tier4_readout_genes_all_features.csv")

    # GENCODE-style multi-feature rows for panel lncRNAs (written by ``save_gene_tables``); not the default for mapping.
    lncrnas_all_features: Path = Path(f"{_APM_WDIR}/lncRNAs_genes_all_features.csv")
    
    # regulatory elements table (STEP 12 writes parquet primary + optional CSV).
    # Consumers should use ``pipeline.regulatory_elements.load_regulatory_element_focus``
    # which prefers ``regulatory_element_focus_with_evidence.parquet`` and falls back to CSV.
    regulatory_elements_table: Path = Path(f"{_APM_WDIR}/regulatory_elements_matching/regulatory_element_focus.csv")
    regulatory_elements_table_with_evidence_parquet: Path = Path(
        f"{_APM_WDIR}/regulatory_elements_matching/regulatory_element_focus_with_evidence.parquet"
    )
    regulatory_elements_table_with_evidence_csv: Path = Path(
        f"{_APM_WDIR}/regulatory_elements_matching/regulatory_element_focus_with_evidence.csv"
    )

    # atac peaks
    atac_peaks_mapping_csv: Path = Path(f"{_APM_WDIR}/atac_peaks/peak_id_mapping.csv")

    # ChIP-seq peaks (working_dir / ENCODE/ and CHIP_ATLAS/ subdirs)
    chip_dir: Path = Path(f"{_APM_WDIR}/CHIP")

    # Unified concatenated chip peak table (cached output of chip_loader)
    chip_unified: Path = Path(f"{_APM_WDIR}/CHIP/unified_chip_peaks.parquet")

    # Populated in __post_init__ (see OUTPUT_SUBDIRS["snv"], SV step 9)
    snv_output_dir: Path = field(init=False)
    sv_chip_enriched_dir: Path = field(init=False)
    cnv_gene_tables_dir: Path = field(init=False)

    def __post_init__(self):
        """Convert strings to Path objects if needed; attach derived module output roots."""
        # Not every ``str`` field is a filesystem path (e.g. matrix column names).
        _str_fields_keep_as_str = frozenset({"screen_exp_inner_file", "rna_gene_col"})
        for field_name, field_value in list(self.__dict__.items()):
            if isinstance(field_value, str) and field_name not in _str_fields_keep_as_str:
                setattr(self, field_name, Path(field_value))
        subs = globals().get("OUTPUT_SUBDIRS") or {}
        self.snv_output_dir = self.working_dir / subs.get("snv", "snv_somatic_annotated")
        self.sv_chip_enriched_dir = self.sv_output_root / "09_chip_enriched"
        self.cnv_gene_tables_dir = self.working_dir / "CNV_TCGA" / "CNV_gene_tables"

    def gene_table_output_paths(self) -> Tuple[Path, ...]:
        """
        Every CSV path written by ``pipeline.genes.gene_loader.save_gene_tables``.

        Order: pipeline/CNV panel multifeature, pipeline/CNV gene-only, CNV union gene-only,
        panel lncRNA multifeature, then Tier 2 / 3 / 4 (gene-only, multifeature each).
        """
        return (
            self.genes_all_features,
            self.genes_only,
            self.cnv_genes,
            self.lncrnas_all_features,
            self.tier2_medium_genes_only,
            self.tier2_medium_genes_all_features,
            self.tier3_cnv_only_genes_only,
            self.tier3_cnv_only_genes_all_features,
            self.tier4_readout_genes_only,
            self.tier4_readout_genes_all_features,
        )


# =============================================================================
# GENE PANEL
# =============================================================================
#
# Tiering follows ``research_plan/01_gene_panel_extended.md`` (canonical; ``extended_genes.md`` is a pointer).
# - ``PRIMARY_GENES``: frozen 66-gene core.
# - ``PIPELINE_GENE_PANEL``: default for ``run_full_pipeline`` (extended unless disabled).
# - ``CNV_GENES``: single deduped union (full integration + Tier 2 + Tier 3 only).
# - ``TIER4_READOUT_GENES``: labels only; included in RNA alias seeding, not CNV/pipeline panel.
# - lncRNA intervals (``main`` / ``run_genes_only``): all GENCODE ``gene_type == lncRNA`` gene rows
#   from ``gencode_gtf_pq``. ``lncrna_csv`` / ``lncrnas_genes_centric`` name the matching *output*
#   consumed by Methylation/SV. ``APM_LNCRNA_INPUT_LEGACY_CSV=1`` restores reading that CSV as input.
#
# Env: ``APM_USE_EXTENDED_GENE_PANEL=0`` forces the 66-gene pipeline panel.
# RNA HGNC scan: ``APM_HGNC_ALIAS_MAX_LINES`` optional cap for dev speed.
# Symbol normalization (VEP, miRTarBase, methylation probe gene lists, RNA):
# ``APM_USE_GENE_SYMBOL_MAPPING=0`` disables HGNC/UCSC/legacy remaps in those loaders.

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


def _dedupe_preserve_order(symbols: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for sym in symbols:
        if sym and sym not in seen:
            seen.add(sym)
            out.append(sym)
    return out


# Tier 1 — full integration (GENCODE v49 ``gene_name`` symbols; see research_plan/01_gene_panel_extended.md)
EXTENDED_PRIMARY_GENES: List[str] = [
    # IFN negative regulation (type-I brake; complements SOCS/PTPN/PIAS in core)
    "USP18",
    # Immune synapse / adhesion (visibility beyond antigen presentation)
    "ICAM1",
    "CD58",
    # MYC amplification program (full integration: immune evasion + enhancer architecture)
    "MYC",
    # MHC-I enhanceosome
    "CIITA",
    "RFX5",
    "RFXAP",
    "RFXANK",
    "NFYA",
    "NFYB",
    "NFYC",
    "POMP",
    "CANX",
    # cGAS–STING axis (GENCODE uses CGAS / STING1)
    "CGAS",
    "STING1",
    "TBK1",
    "IRF3",
    "IRF7",
    "STAT2",
    "IRF9",
    "IFNAR1",
    "IFNAR2",
    # Peptide trimming / HLA-E axis
    "ERAP1",
    "ERAP2",
    # lncRNAs (first-class regulatory subjects)
    "HCP5",
    "HCG18",
    "HLA-F-AS1",
    "MICB-DT",
    "LINC01149",
    "NEAT1",
    "MALAT1",
    "HOTAIR",
]

# Subset of ``EXTENDED_PRIMARY_GENES`` that are lncRNA symbols (same order as in that list).
TIER1_LNCRNA_GENES: List[str] = [
    "HCP5",
    "HCG18",
    "HLA-F-AS1",
    "MICB-DT",
    "LINC01149",
    "NEAT1",
    "MALAT1",
    "HOTAIR",
]

# Tier 2 — medium depth (expression + promoter methylation + CNV; not full enhancer subjects)
TIER2_MEDIUM_GENES: List[str] = [
    # Tryptophan–kynurenine axis (IFN-inducible immunosuppression + downstream sensor)
    "IDO1",
    "AHR",
    # Adenosine axis (tumor-state immunosuppression)
    "NT5E",   # CD73
    # "Don't eat me" (tumor-cell surface; macrophage receptor is readout-tier)
    "CD47",
    # Additional checkpoints / ligands (medium depth; avoid full enhancer investment)
    "LGALS9",  # TIM-3 ligand
    "CD276",   # B7-H3
    "DNMT1",
    "DNMT3A",
    "DNMT3B",
    "TET1",
    "TET2",
    "TET3",
    "EZH2",
    "SUZ12",
    "EED",
    "KDM6A",
    # H3K9/H3K36 regulators (balance Tier-2 demethylase families)
    "KDM4A",
    "KDM4B",
    "KDM4C",
    "KDM5A",
    "KDM5B",
    "KDM5C",
    "SETD2",
    "HDAC1",
    "HDAC2",
    "HDAC3",
    # Class II HDACs (family balance; immune-silencing literature often spans HDAC classes)
    "HDAC4",
    "HDAC5",
    "HDAC6",
    "HDAC7",
    "KAT2A",
    # Phosphatase family balance (Tier-0 has PTPN2/PTPN11)
    "PTPN1",
    "PTPN6",
]

# Tier 3 — CNV / dosage only (single definition; folded into CNV_GENES below)
TIER3_CNV_ONLY_GENES: List[str] = [
    "NFKB1",
    "NFKB2",
    "RELA",
    "RELB",
    "REL",
    "KPNA6",
    "XPO1",
]

# Stratifier tier — genes we want in both CNV tables and SNV-lite summaries/covariates.
# These are not “full integration” enhancer/TAD subjects; they’re tumor-state anchors.
TIER_SNV_CNV_STRATIFIER_GENES: List[str] = [
    # HRD anchors (hot/cold hypothesis conditions on HRD genotypes directly)
    "BRCA1",
    "BRCA2",
    # HRD core anchors (mechanistically central; committee-defensible)
    "PALB2",
    "ATM",
    "TP53",
    "PTEN",
    "PIK3CA",
]

# Tier 4 — readout / phenotype labels (not regulatory pipeline subjects; kept for RNA/Thorsson alignment)
TIER4_READOUT_GENES: List[str] = [
    # Immune geography / myeloid programs (expression-only readouts)
    "CXCL13",
    "CCL2",
    "CSF1",
    # Adenosine axis (often immune-cell dominated in BRCA; track as readout)
    "ENTPD1",  # CD39
    # "Don't eat me" receptor (macrophage; readout)
    "SIRPA",
    "GZMA",
    "PRF1",
    "CD8A",
    "CD8B",
    "CD4",
    "FOXP3",
    "NKG7",
    "KLRD1",
    "KLRK1",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "CTLA4",
]

FULL_INTEGRATION_GENES: List[str] = _dedupe_preserve_order(
    list(PRIMARY_GENES) + list(EXTENDED_PRIMARY_GENES)
)

# One CNV definition: full-integration + medium-depth + CNV-only tiers (research_plan/01)
CNV_GENES: List[str] = _dedupe_preserve_order(
    list(FULL_INTEGRATION_GENES)
    + list(TIER2_MEDIUM_GENES)
    + list(TIER3_CNV_ONLY_GENES)
    + list(TIER_SNV_CNV_STRATIFIER_GENES)
)

# When True, ``PIPELINE_GENE_PANEL`` matches Tier-1 extended integration (≈100 genes).
# Set env ``APM_USE_EXTENDED_GENE_PANEL=0`` to force the legacy 66-gene default.
USE_EXTENDED_PRIMARY_PANEL: bool = os.environ.get(
    "APM_USE_EXTENDED_GENE_PANEL", "1"
).strip().lower() not in ("0", "false", "no")

PIPELINE_GENE_PANEL: List[str] = (
    list(FULL_INTEGRATION_GENES)
    if USE_EXTENDED_PRIMARY_PANEL
    else list(PRIMARY_GENES)
)


def use_legacy_lncrna_intervals_csv_input() -> bool:
    """
    When True, ``main`` / ``run_genes_only`` / rebuild scripts load lncRNA intervals from
    ``PATHS.lncrna_csv`` instead of deriving them from the GENCODE genes table.

    Env: ``APM_LNCRNA_INPUT_LEGACY_CSV=1`` (old circular default: input == matching output).
    """
    return os.environ.get("APM_LNCRNA_INPUT_LEGACY_CSV", "0").strip().lower() in (
        "1",
        "true",
        "yes",
    )


RNA_EXPRESSION_ALIAS_SEED_SYMBOLS: frozenset[str] = frozenset(
    _dedupe_preserve_order(
        list(FULL_INTEGRATION_GENES)
        + list(TIER2_MEDIUM_GENES)
        + list(TIER3_CNV_ONLY_GENES)
        + list(TIER_SNV_CNV_STRATIFIER_GENES)
        + list(TIER4_READOUT_GENES)
    )
)

# Preferred name for the same seed set (RNA + all other modules use this frozenset).
PANEL_ALIAS_SEED_SYMBOLS = RNA_EXPRESSION_ALIAS_SEED_SYMBOLS

# Rare symbols in external matrices that should fold into GENCODE names above.
LEGACY_DATASET_SYMBOL_RENAMES: Dict[str, str] = {
    "TMEM173": "STING1",
    "MB21D1": "CGAS",
    "STING": "STING1",
}

# Mature-arm miRNA ids prioritized in ``research_plan/01_gene_panel_extended.md`` (§9 miRNA table).
# Expression matrices use **MIMAT** row ids; resolve those first via ``mirna_mature_loci.csv`` in
# ``pipeline.genes.panel_alias_registry``. These strings are the **arm-level** canonical labels.
PANEL_MIRNA_TIER_ARM_IDS: Tuple[str, ...] = (
    "hsa-miR-148a-3p",
    "hsa-miR-152-3p",
    "hsa-miR-9-5p",
    "hsa-miR-125a-5p",
    "hsa-miR-125b-5p",
    "hsa-miR-27a-3p",
    "hsa-miR-34a-5p",
    "hsa-miR-346-5p",
    "hsa-miR-200c-3p",
    "hsa-miR-155-5p",
    "hsa-miR-146a-5p",
)


# =============================================================================
# UCSC RNA-seq gene symbol changes
# =============================================================================

UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE: Dict[str, str] = {
    "DKFZp686O24166": "NCR3LG1",
    "PVRL2": "NECTIN2",
}


# =============================================================================
# BIOSAMPLE / CELL LINE PANELS
# =============================================================================

@dataclass
class BiosampleConfig:
    """Biosample panels for different data sources."""
    
    # SCREEN experimental (3D chromatin) — canonical keys in outputs / gene_links
    # (upstream ENCODE "Biosample" strings are mapped in biosample_names.canonicalize_screen_biosample)
    screen_exp: List[str] = field(default_factory=lambda: [
        "MCF7",
        "MCF10A",
        "breast_mammary_tissue",
        "breast_epithelium_tissue_female_adult_53_years",
        "mammary_epithelial_cell_female_adult_19_years",
    ])
    
    # SCREEN computational — same: values here are canonical output keys
    screen_comp: List[str] = field(default_factory=lambda: [
        "MCF7",
        "MCF10A",
        "T47D",
        "mammary_epithelial_cell",
        "breast_epithelium",
    ])
    
    # ABC model cell types — exact strings from the ABC file "CellType" column (filter only).
    # Nested ABC_full keys use canonical names from biosample_names.canonicalize_abc_cell_type.
    abc_celltypes: List[str] = field(default_factory=lambda: [
        "MCF-7-ENCODE",
        "MCF10A-Ji2017",
        "MCF10A_treated_with_TAM24hr-Ji2017",
        "MDA-MB-231",
        "mammary_epithelial_cell-Roadmap",
        "breast_epithelium-ENCODE",
    ])
    
    # HiChIP cell lines — subdirectory names under PATHS.hichip_dir (disk layout).
    # Output columns / gene_links.hichip keys use biosample_names.canonical_hichip_output_key.
    hichip_panel: List[str] = field(default_factory=lambda: [
        "B80T5", "HMEC", "K5plusK19plus", "MCF7", "MDA-MB-231", "T47D", "ZR751"
    ])
    
    # Cell lines for cCRE signal integration — subdirectory names under PATHS.cell_lines_dir (disk).
    # DataFrame columns use canonical_ccre_signal_column_name (e.g. MCF-7 on disk → MCF7 column).
    ccre_signal_cell_lines: List[str] = field(default_factory=lambda: [
        "MCF7", "HMEC1", "HMEC2", "MCF10A", "breast_tissue"
    ])

    # Canonical spellings; ChIP loaders normalize aliases (e.g. MCF-7 -> MCF7).
    chip_brca_celltypes: List[str] = field(default_factory=lambda: [
        "MCF7",
        "T47D",
        "HCC1954",
        "MDA-MB-231",
        "MDA-MB-468",
        "SKBR3",
        "HMEC",
        "MCF10A",
    ])

    # ChIP-Atlas unified filter: keep peaks whose cell line matches ``chip_brca_celltypes``
    # (after ``normalize_cell_line_label``) OR whose raw label mentions breast/mammary / keywords.
    chip_atlas_mammary_tissue_keywords: Tuple[str, ...] = (
        "mammary",
        "breast",
        "ductal",
        "luminal",
        "basal",
        "epithelium",
        "MCF",
        "MDA-MB",
        "HCC11",
        "HCC70",
        "HCC193",
        "T47D",
        "SKBR",
        "BT474",
        "ZR75",
        "CAMA",
        "UACC",
        "Hs578",
        "DU4475",
    )

    # Normalized cell-line label -> interpretable subtype/state (ChIP-Atlas / concordance).
    # Keys are matched after ``normalize_cell_line_label``; extend as needed.
    # Keys are ``normalize_cell_line_label`` outputs (canonical tokens).
    chip_cell_line_subtype_map: Dict[str, str] = field(default_factory=lambda: {
        "MCF7": "luminal_breast_cell_line",
        "T47D": "luminal_breast_cell_line",
        "HCC1954": "luminal_breast_cell_line",
        "MDA-MB-231": "basal_breast_cell_line",
        "MDA-MB-468": "basal_breast_cell_line",
        "SKBR3": "HER2amp_breast_cell_line",
        "HMEC": "normal_mammary_epithelium",
        "HMEC1": "normal_mammary_epithelium",
        "HMEC2": "normal_mammary_epithelium",
        "MCF10A": "immortalized_mammary_epithelium",
        "HCC1187": "basal_breast_cell_line",
        "HCC1143": "basal_breast_cell_line",
        "HCC70": "basal_breast_cell_line",
        "HCC1937": "basal_breast_cell_line",
        "HCC1569": "basal_breast_cell_line",
        "HCC1395": "basal_breast_cell_line",
        "HCC1500": "basal_breast_cell_line",
        "MDA-MB-361": "luminal_breast_cell_line",
        "MDA-MB-453": "basal_breast_cell_line",
        "MDA-MB-157": "basal_breast_cell_line",
        "BT474": "HER2amp_breast_cell_line",
        "ZR751": "luminal_breast_cell_line",
        "CAMA1": "luminal_breast_cell_line",
        "UACC812": "luminal_breast_cell_line",
        "MCF10DCIS.com": "DCIS_mammary_model",
    })


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

    # CNV -> gene summarization thresholds
    # Minimum % of gene body overlapped by a CNV segment to call loss/gain/amp at gene-level.
    cnv_gene_min_body_overlap_pct: float = 30.0
    # ASCAT3 / GDC gene-level TSV: regulatory neighborhood min-CN window (bp) + output stem.
    cnv_ascat3_regulatory_window_bp: int = 250_000
    cnv_ascat3_gene_table_stem: str = "cnv_gene_calls_ascat3"
    
    # Streaming chunk sizes
    screen_chunksize: int = 300_000
    abc_chunksize: int = 500_000
    
    # Binning for fast prefilter
    coarse_bin_size: int = 100_000
    
    # ATAC peak thresholds
    atac_gene_window_bp: int = 1_000_000
    atac_ccre_max_distance: int = 0  # 0 = overlap only
    atac_lncrna_window_bp: int = 1_000_000
    # Step 10 ATAC: TAD passes merge nested dicts on every peak × every biosample. Row-chunking
    # runs **all** biosamples on each chunk then ``pd.concat`` — same result as a single pass,
    # but peak RAM stays ~O(chunk_rows × n_biosamples) during the inner biosample loop.
    # Default is conservative for 16 GB machines; raise (e.g. 20000) if you have headroom.
    # Set ``0`` to disable chunking (legacy one-shot; very high RAM).
    # Env: ``APM_ATAC_TAD_CHUNK_ROWS`` (integer; ``0`` = disable).
    atac_tad_chunk_rows: int = 4_096
    # When ``build_atac_peak_table(..., tad_parquet_stream=...)`` streams TAD chunks to disk,
    # each annotated chunk is sliced to this many rows before ``pa.Table.from_pandas``.
    # Env: ``APM_ATAC_TAD_STREAM_WRITE_ROWS`` (integer; minimum 256).
    atac_tad_stream_write_rows: int = 2_048

    # SNV somatic filtering thresholds
    snv_min_tumor_vaf: float = 0.05          # Minimum tumor VAF (5%)
    snv_max_normal_vaf: float = 0.02         # Maximum normal VAF (2%)
    snv_min_tlod: float = 6.0                # Minimum TLOD score
    snv_min_popaf: float = 3.0               # Minimum -log10(pop AF) (i.e., AF ≤ 1e-3)
    snv_min_tumor_dp: int = 20               # Minimum tumor read depth
    snv_min_normal_dp: int = 10              # Minimum normal read depth
    snv_use_popaf: bool = True               # Whether to filter on population AF
    snv_require_tlod: bool = False           # Whether TLOD must be present
    snv_require_pass: bool = True            # Require FILTER == "PASS"
    # Optional SNV FIMO (``load_mutect_snv_vcf(..., run_fimo=True)``): window half-width and cap
    snv_fimo_flank_bp: int = 30
    snv_fimo_max_hits_per_variant: int = 50

    # SV pipeline thresholds
    sv_gene_window_bp: int = 1_000_000
    sv_element_window_bp: int = 1_000_000
    sv_flank_size_bp: int = 150
    sv_promoter_upstream_bp: int = 2000
    sv_promoter_downstream_bp: int = 500
    sv_proximal_window_bp: int = 5000
    neojunction_window_bp: int = 500_000

    
    # SV filtering thresholds (strict)
    sv_min_tumor_sr_alt: int = 2
    sv_min_tumor_alt: int = 8
    sv_max_normal_alt: int = 1
    sv_min_somatic_score: int = 25
    
    # SV filtering thresholds (lenient)
    sv_lenient_min_tumor_sr_alt: int = 1
    sv_lenient_min_tumor_pr_alt: int = 5
    sv_lenient_min_somatic_score: int = 15
    
    # FIMO thresholds
    fimo_pvalue_threshold: float = 1e-4

    # FIMO ↔ SV recombine (STEP 5 in ``pipeline/SV/pipeline.py``): ``bedtools intersect -wa -wb``
    # can emit tens of millions of rows. The merged BED is streamed in row chunks; motif hits are
    # retained per TF (short token before the first ``.`` in FIMO ``motif_id``) with best p-value
    # within each (SV, TF) or (SV, regulatory elem, TF) bucket—see ``recombine_sv_fimo``.
    sv_fimo_recombine_chunk_rows: int = 150_000
    sv_fimo_max_flank_hits_per_tf: int = 400
    sv_fimo_max_elem_hits_per_tf: int = 300
    # If the merged BED is at least this size (bytes), STEP 5 never loads the full SV CSV at once:
    # rows are read/written in ``sv_fimo_recombine_sv_csv_chunk_rows`` chunks after motif attach.
    # Set to 0 to always use the single in-memory SV dataframe (legacy behavior).
    sv_fimo_recombine_stream_if_merged_bed_bytes: int = 300_000_000
    sv_fimo_recombine_sv_csv_chunk_rows: int = 800
    # Processed SV CSV stems (filename without ``.csv``, e.g. ``..._strict_sv_set``) to skip for
    # FIMO FASTA/FIMO/intersect/recombine: avoids OOM or huge ``06`` merges. ``07_final`` gets a
    # copy of the corresponding ``02_processed_sv_csv`` rowset (no motif attachment). Empty the
    # list to process every sample again.
    # sv_fimo_recombine_skip_csv_basenames: List[str] = field(default_factory=lambda: [
    #     "TCGA-3C-AALI-01A_strict_sv_set",
    # ])

    # SV-ChIP intersection window. Mirrors sv_gene_window_bp default;
    sv_chip_window_bp: int = 1_000_000

    # Methylation thresholds
    meth_hypermeth_threshold: float = 0.7
    meth_hypometh_threshold: float = 0.3
    meth_intermediate_low: float = 0.3
    meth_intermediate_high: float = 0.7
    
    # Promoter definition for methylation (can be same as gene promoter) !!!!!!!!!!
    meth_promoter_upstream_bp: int = 2000
    meth_promoter_downstream_bp: int = 500
    
    # M-value conversion offset (to avoid log(0))
    meth_m_value_offset: float = 0.001
    
    # QC thresholds
    meth_min_valid_probe_pct: float = 50.0  # Minimum % valid probes per sample
    meth_detection_pval_threshold: float = 0.01  # If detection p-values available

   

# =============================================================================
# SV MOTIF SCANNING CONFIGURATION
# =============================================================================

# Target TF symbols for motif extraction from MEME databases (substring match on ``MOTIF`` name).
# Synced with **HOCOMOCO v14 H14CORE** (and the same names appear in H14INVIVO / H14INVITRO / H14RSNP
# for these factors). See ``pipeline/md/module_specific_processing_md/sv_→_motif_pipeline_vep_fimo.md``
# and ``scripts/motifs/hocomoco_multi_bundle_tf_report.py``.
# **Not** present as named motifs in H14 bundles (fill from JASPAR/CIS-BP/other if needed): ``RELA``,
# ``NLRC5``, ``CIITA``, ``EZH2``, ``SUZ12``, ``BRD4``.
SV_TARGET_TF_SYMBOLS: List[str] = [
    # STATs
    "STAT1", "STAT2", "STAT3",
    # IRFs (IRF9 in H14; CIITA absent — add externally if required)
    "IRF1", "IRF2", "IRF3", "IRF9",
    # NF-κB family (RELA / p65 not separate motif names in H14; REL + RELB + NFKB* are)
    "NFKB1", "NFKB2", "RELB", "REL",
    # AP-1 / bZIP
    "FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3",
    # CTCF family
    "CTCF", "CTCFL",
    "MYC",
    # NF-Y (all three subunits named in H14)
    "NFYA", "NFYB", "NFYC",
    # MHC-II / breast context (H14)
    "RFX5",
    "FOXA1", "ESR1", "GATA3",
    # TGF-β (SMAD1, SMAD2, … in HOCOMOCO names)
    "SMAD",
]

# Allowed SV types for processing
ALLOWED_SV_TYPES: List[str] = ["DEL", "DUP", "INV", "INS", "BND"]


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
# VEP ANNOTATION CONFIGURATION
# =============================================================================

# Default VEP CSQ format string (VEP v115 with --everything flag)
VEP_CSQ_FORMAT: str = (
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

# Consequence terms for classification
SNV_CONSEQUENCE_TERMS: Dict[str, List[str]] = {
    "missense": ["missense_variant"],
    "nonsense": ["stop_gained"],
    "frameshift": ["frameshift_variant"],
    "splice": ["splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"],
    "synonymous": ["synonymous_variant"],
    "regulatory": [
        "regulatory_region_variant",
        "TF_binding_site_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
    ],
}



# =============================================================================
# OUTPUT SUBDIRECTORIES
# =============================================================================

OUTPUT_SUBDIRS: Dict[str, str] = {
    "lncrna_matching": "lncRNA_matching",
    "regulatory_elements": "regulatory_elements_matching",
    "mirna": "miRNA",
    "evidence": "evidence",
    "atac_peaks": "atac_peaks",
    # Per-sample SNV tables (VEP + cCRE + miRNA overlap); see SNV.save_snv_outputs
    # Keep all SNV artifacts (inputs and outputs) under data/SNV/ for locality.
    "snv": "SNV/somatic_annotated",
    # SV pipeline subdirectories (numeric steps live under PATHS.sv_output_root)
    "sv_vep_vcfs": "sv_pipeline/01_vep_vcfs",
    "sv_processed_csv": "sv_pipeline/02_processed_sv_csv",
    "sv_bed": "sv_pipeline/03_sv_bed",
    "sv_fasta": "sv_pipeline/04_sv_fasta",
    "sv_fimo_tsv": "sv_pipeline/05_fimo_tsv",
    "sv_fimo_merged": "sv_pipeline/06_sv_fimo_merged",
    "sv_final": "sv_pipeline/07_final_sv_with_fimo",
    "sv_logs": "sv_pipeline/logs",
    # Methylation subdirectories (all under …/Methylation/)
    "methylation": "Methylation",
    "methylation_reference": "Methylation/reference",
    "methylation_per_sample": "Methylation/per_sample",
    "methylation_cohort": "Methylation/cohort",
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
