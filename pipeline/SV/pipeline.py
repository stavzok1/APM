"""
Main SV pipeline orchestrator.

Coordinates all SV processing modules for:
- VCF loading and filtering
- Spatial mapping to genes, elements, lncRNAs
- VEP annotation integration
- FIMO motif scanning
- Evidence aggregation

When paths are not provided, defaults are loaded from the main pipeline config.
"""

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple

import pandas as pd

from pipeline.genes.symbol_normalization import normalize_annotation_gene_names


# =============================================================================
# IMPORTS FROM OTHER MODULES
# =============================================================================

from .vcf_loader import load_manta_sv_vcf, breakend_to_vcf_alt
from .sv_filtering import get_strict_sv_set, get_lenient_sv_set
from .spatial_mapping import map_svs_to_genes, map_svs_to_elements, map_svs_to_lncrnas, map_svs_to_mirnas
from .vep_annotation import add_vep_hits_columns
from .bed_intervals import build_sv_flanks_and_overlaps_bed, create_beds_from_directory
from .motif_scanning import (
    recombine_all_sv_fimo,
    summarize_motif_hits,
    load_all_ccre_fimo,
    annotate_bnd_neojunction_motifs,
    summarize_neojunction_motifs,
    _sv_csv_stem_from_flanks_bed_stem,
    _sv_fimo_recombine_skip_csv_stems,
)
from ..CHIP.chip_loader import load_unified_chip
from ..CHIP.sv_chip_intersect import annotate_sv_csvs_with_chip
from ..sample_ids import add_tcga_id_columns_inplace
from ..utils import harmonize_chrom_column


def _load_sv_lncrna_feature_tables(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """
    Load lncRNAs for SV mapping.

    Accepts either a GENCODE-style multi-row table (``feature``, ``frame``) or the
    gene-centric ``lncRNAs_with_genes_1000000bp.csv`` export (one row per lncRNA gene).
    """
    raw = normalize_annotation_gene_names(
        pd.read_csv(path, low_memory=False),
        ("gene_name",),
    )
    if "feature" in raw.columns and "frame" in raw.columns:
        lncrnas = raw[raw["frame"].isin(["0", "."])].copy()
        lncrnas_only = lncrnas[lncrnas["feature"] == "gene"].copy()
    else:
        raw, _ = harmonize_chrom_column(raw)
        lncrnas_only = raw
        lncrnas = raw.copy()
    names = lncrnas_only["gene_name"].astype(str).tolist()
    return lncrnas, lncrnas_only, names


def _build_sv_vcf_feature_bundle(
    genes_path: Path,
    elements_path: Path,
    samples_tsv_path: Path,
    lncrnas_path: Optional[Path],
    mirna_path: Optional[Path],
) -> Dict[str, Any]:
    """Load gene/element/lncRNA/miRNA tables once for SV VCF batch processing."""
    from ..regulatory_elements import load_regulatory_element_focus

    # Prefer parquet when available (full GENCODE feature table is large).
    if str(genes_path).endswith((".parquet", ".pq")):
        genes_raw = pd.read_parquet(genes_path)
    else:
        genes_raw = pd.read_csv(genes_path, low_memory=False)

    genes = normalize_annotation_gene_names(
        genes_raw,
        ("gene_name",),
    )
    genes = genes[genes["frame"].isin(["0", "."])]
    genes_only = genes[genes["feature"] == "gene"]
    primary_genes = genes_only["gene_name"].tolist()
    elements = load_regulatory_element_focus(elements_path)
    lncrnas: Optional[pd.DataFrame] = None
    lncrnas_only: Optional[pd.DataFrame] = None
    lncrna_names: List[str] = []
    if lncrnas_path is not None and Path(lncrnas_path).exists():
        lncrnas, lncrnas_only, lncrna_names = _load_sv_lncrna_feature_tables(Path(lncrnas_path))
    mirnas: Optional[pd.DataFrame] = None
    if mirna_path is not None and Path(mirna_path).exists():
        mirnas = pd.read_csv(mirna_path)
    annotations = pd.read_csv(samples_tsv_path, sep="\t")
    return {
        "genes": genes,
        "genes_only": genes_only,
        "elements": elements,
        "lncrnas": lncrnas,
        "lncrnas_only": lncrnas_only,
        "mirnas": mirnas,
        "primary_genes": primary_genes,
        "lncrna_names": lncrna_names,
        "annotations": annotations,
    }


# =============================================================================
# CONFIG LOADING
# =============================================================================

def _get_config():
    """
    Try to import config from main pipeline.
    Returns (PATHS, THRESHOLDS, VEP_CSQ_FORMAT) or (None, None, None) if unavailable.
    """
    try:
        # Try importing from parent pipeline package
        import sys
        module_dir = Path(__file__).parent.parent
        if str(module_dir) not in sys.path:
            sys.path.insert(0, str(module_dir))
        
        from config import PATHS, THRESHOLDS
        
        # Try to get VEP_CSQ_FORMAT if defined
        try:
            from config import VEP_CSQ_FORMAT
        except ImportError:
            VEP_CSQ_FORMAT = None
        
        return PATHS, THRESHOLDS, VEP_CSQ_FORMAT
    except ImportError:
        return None, None, None


def _get_path_with_default(provided: Optional[Path], config_attr: str, default: Optional[Path] = None) -> Optional[Path]:
    """Get path from provided value, config, or default."""
    if provided is not None:
        return Path(provided)
    
    PATHS, _, _ = _get_config()
    if PATHS is not None and hasattr(PATHS, config_attr):
        return getattr(PATHS, config_attr)
    
    return default


def _get_threshold_with_default(provided: Optional[Any], config_attr: str, default: Any) -> Any:
    """Get threshold from provided value, config, or default."""
    if provided is not None:
        return provided
    
    _, THRESHOLDS, _ = _get_config()
    if THRESHOLDS is not None and hasattr(THRESHOLDS, config_attr):
        return getattr(THRESHOLDS, config_attr)
    
    return default


# =============================================================================
# SHELL SCRIPT PATHS (relative to module)
# =============================================================================

MODULE_DIR = Path(__file__).parent


def _get_script_path(script_name: str) -> Path:
    """Get path to shell script in module directory."""
    return MODULE_DIR / "scripts" / script_name


# =============================================================================
# SINGLE VCF PROCESSING
# =============================================================================

def process_single_vcf(
    vcf_path: Path,
    genes: pd.DataFrame,
    genes_only: pd.DataFrame,
    elements: pd.DataFrame,
    lncrnas: Optional[pd.DataFrame] = None,
    lncrnas_only: Optional[pd.DataFrame] = None,
    mirnas: Optional[pd.DataFrame] = None,
    primary_genes: Optional[List[str]] = None,
    lncrna_names: Optional[List[str]] = None,
    csq_description: Optional[str] = None,
    window: int = 1_000_000,
    filtering: str = "strict",
) -> pd.DataFrame:
    """
    Process a single VCF file through the SV pipeline.
    
    Args:
        vcf_path: Path to VCF file
        genes: Full gene features DataFrame
        genes_only: Gene-level only DataFrame (feature == "gene")
        elements: Regulatory elements DataFrame with elem_id column
        lncrnas: Optional full lncRNA features DataFrame
        lncrnas_only: Optional lncRNA-level only DataFrame
        primary_genes: List of primary gene symbols
        lncrna_names: List of lncRNA names
        csq_description: VEP CSQ format description
        window: Maximum distance for mapping
        filtering: "strict" or "lenient"
    
    Returns:
        Processed SV DataFrame
    """
    print(f"Processing: {vcf_path}")
    
    # Load VCF
    df, normal_sample, tumor_sample = load_manta_sv_vcf(str(vcf_path))
    
    # Filter
    if filtering == "strict":
        sv_df = get_strict_sv_set(df)
    else:
        sv_df = get_lenient_sv_set(df)
    
    print(f"  Filtered to {len(sv_df)} SVs ({filtering})")
    
    if sv_df.empty:
        return sv_df
    
    # Map to genes
    sv_df = map_svs_to_genes(sv_df, genes, genes_only, window=window)
    
    # Map to elements
    sv_df = map_svs_to_elements(sv_df, elements, window=window)
    
    # Map to lncRNAs if provided
    if lncrnas is not None and lncrnas_only is not None:
        sv_df = map_svs_to_lncrnas(sv_df, lncrnas, lncrnas_only, window=window)

    # Map to miRNAs if provided
    if mirnas is not None:
        sv_df = map_svs_to_mirnas(sv_df, mirnas, window=window)
    
    # Convert BND alt format
    sv_df = sv_df.rename(columns={"alt": "orig_alt"})
    sv_df["alt"] = sv_df["orig_alt"].apply(breakend_to_vcf_alt)
    sv_df["alt"] = sv_df["alt"].fillna(sv_df["orig_alt"])
    
    # Add VEP annotations if CSQ column exists
    if "CSQ" in sv_df.columns and csq_description is not None:
        if primary_genes is None:
            primary_genes = []
        if lncrna_names is None:
            lncrna_names = []
        gene_symbol_mapping = None
        if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        ):
            from ..genes.symbol_normalization import default_symbol_mapping

            gene_symbol_mapping = default_symbol_mapping()
        sv_df = add_vep_hits_columns(
            sv_df,
            csq_description,
            primary_genes,
            lncrna_names,
            gene_symbol_mapping=gene_symbol_mapping,
        )
    else:
        # VEP was skipped or the input Manta VCF lacks CSQ. Keep the atlas-aligned
        # schema by emitting the VEP columns with empty defaults so downstream code
        # can rely on them (and the non-VEP case is obvious: all False / empty).
        if not sv_df.empty:
            _SV_VEP_LIST_COLS = ("gene_hits_vep", "regulatory_hits_vep", "motif_hits")
            _SV_VEP_FLAG_COLS = (
                "hits_canonical",
                "has_missense", "has_nonsense", "has_frameshift", "has_splice_effect",
                "has_missense_canonical", "has_nonsense_canonical",
                "has_frameshift_canonical", "has_splice_effect_canonical",
                "has_missense_mane", "has_nonsense_mane",
                "has_frameshift_mane", "has_splice_effect_mane",
            )
            for col in _SV_VEP_LIST_COLS:
                if col not in sv_df.columns:
                    sv_df[col] = [[] for _ in range(len(sv_df))]
            if "gene_symbols" not in sv_df.columns:
                sv_df["gene_symbols"] = ""
            for col in _SV_VEP_FLAG_COLS:
                if col not in sv_df.columns:
                    sv_df[col] = False

    # ``flank_motif_hits`` is attached by the motif scanning step; ensure an empty
    # placeholder exists here so step 1 CSVs are already atlas-aligned.
    if "flank_motif_hits" not in sv_df.columns:
        sv_df["flank_motif_hits"] = [[] for _ in range(len(sv_df))] if not sv_df.empty else pd.Series(dtype=object)

    return sv_df


# =============================================================================
# VCF PROCESSING (BATCH)
# =============================================================================

def run_sv_vcf_processing(
    vcf_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    genes_path: Optional[Path] = None,
    elements_path: Optional[Path] = None,
    samples_tsv_path: Optional[Path] = None,
    lncrnas_path: Optional[Path] = None,
    mirna_path: Optional[Path] = None,
    csq_description: Optional[str] = None,
    window: Optional[int] = None,
    filtering: str = "strict",
) -> List[Path]:
    """
    Process all VCFs in a directory.
    
    Args:
        vcf_dir: Directory containing VCF files (default: config.PATHS.sv_vcf_dir)
        output_dir: Directory for output CSVs (default: config.PATHS.sv_output_root / "02_processed_sv_csv")
        genes_path: Path to gene features CSV (default: config.PATHS.gencode_gtf_pq)
        elements_path: Path to regulatory elements CSV (default: config.PATHS.working_dir / "regulatory_elements_matching/coding_element_focus.csv")
        samples_tsv_path: Path to sample annotations TSV (default: config.PATHS.sv_samples_tsv)
        lncrnas_path: Optional path to lncRNA features CSV (default: config.PATHS.lncrna_csv)
        csq_description: VEP CSQ format description (default: config.VEP_CSQ_FORMAT)
        window: Maximum distance for mapping (default: config.THRESHOLDS.sv_gene_window_bp)
        filtering: "strict" or "lenient"
    
    Returns:
        List of created output file paths
    """
    # Get defaults from config
    PATHS, THRESHOLDS, VEP_CSQ_FORMAT_DEFAULT = _get_config()
    
    vcf_dir = _get_path_with_default(vcf_dir, "sv_vcf_dir")
    if vcf_dir is None:
        raise ValueError("vcf_dir must be provided or configured in config.PATHS.sv_vcf_dir")
    vcf_dir = Path(vcf_dir)
    
    if output_dir is None:
        sv_output_root = _get_path_with_default(None, "sv_output_root")
        if sv_output_root is not None:
            output_dir = sv_output_root / "02_processed_sv_csv"
        else:
            raise ValueError("output_dir must be provided or configured in config.PATHS.sv_output_root")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # SV spatial mapping supports CDS/UTR/exon/start/stop flags, which require a feature-rich table
    # (not the gene-only panel export). Prefer the panel-sized features CSV when available; fall back
    # to the full GENCODE parquet for genome-wide mapping.
    genes_path = (
        _get_path_with_default(genes_path, "genes_all_features")
        or _get_path_with_default(None, "gencode_gtf_pq")
        or _get_path_with_default(None, "genes_only")
    )
    if genes_path is None:
        raise ValueError(
            "genes_path must be provided or configured in config.PATHS.gencode_gtf_pq / genes_all_features / genes_only"
        )
    
    elements_path = _get_path_with_default(elements_path, "regulatory_elements_table")
    if elements_path is None:
        # Try working_dir path
        working_dir = _get_path_with_default(None, "working_dir")
        if working_dir is not None:
            elements_path = working_dir / "regulatory_elements_matching" / "coding_element_focus.csv"
    if elements_path is None:
        raise ValueError("elements_path must be provided")
    
    samples_tsv_path = _get_path_with_default(samples_tsv_path, "sv_samples_tsv")
    if samples_tsv_path is None:
        raise ValueError("samples_tsv_path must be provided or configured in config.PATHS.sv_samples_tsv")
    
    lncrnas_path = _get_path_with_default(lncrnas_path, "lncrna_csv")
    # For arm-resolution by coordinates, prefer mature-arm loci table when present.
    mirna_path = _get_path_with_default(mirna_path, "mirna_mature_loci_csv")
    if mirna_path is None:
        mirna_path = _get_path_with_default(mirna_path, "mirna_path")
    
    if csq_description is None:
        csq_description = VEP_CSQ_FORMAT_DEFAULT
    
    window = _get_threshold_with_default(window, "sv_gene_window_bp", 1_000_000)

    bundle = _build_sv_vcf_feature_bundle(
        Path(genes_path),
        Path(elements_path),
        Path(samples_tsv_path),
        Path(lncrnas_path) if lncrnas_path is not None else None,
        Path(mirna_path) if mirna_path is not None else None,
    )
    genes = bundle["genes"]
    genes_only = bundle["genes_only"]
    elements = bundle["elements"]
    lncrnas = bundle["lncrnas"]
    lncrnas_only = bundle["lncrnas_only"]
    mirnas = bundle["mirnas"]
    primary_genes = bundle["primary_genes"]
    lncrna_names = bundle["lncrna_names"]
    annotations = bundle["annotations"]

    created_files: List[Path] = []

    for vcf_file in sorted(os.listdir(vcf_dir)):
        if not (vcf_file.endswith(".vcf") or vcf_file.endswith(".vcf.gz")):
            continue

        vcf_path = vcf_dir / vcf_file

        if vcf_file.endswith(".vcf.gz"):
            ann_row = annotations[annotations["File Name"] == vcf_file]
            if ann_row.empty:
                ann_row = annotations[annotations["File Name"] == vcf_file.replace(".vcf.gz", ".vcf")]
        else:
            ann_row = annotations[annotations["File Name"] == vcf_file + ".gz"]
            if ann_row.empty:
                ann_row = annotations[annotations["File Name"] == vcf_file]

        if ann_row.empty:
            sample_id = (
                vcf_file.replace(".vcf.gz", "").replace(".vcf", "")
                if vcf_file.endswith((".vcf", ".vcf.gz"))
                else Path(vcf_file).stem
            )
        else:
            tumor_desc = ann_row["Tumor Descriptor"].iloc[0]
            sample_ids = ann_row["Sample ID"].iloc[0]
            if str(tumor_desc).startswith("Not"):
                sample_id = sample_ids.split(",")[1].strip()
            else:
                sample_id = sample_ids.split(",")[0].strip()

        sv_df = process_single_vcf(
            vcf_path=vcf_path,
            genes=genes,
            genes_only=genes_only,
            elements=elements,
            lncrnas=lncrnas,
            lncrnas_only=lncrnas_only,
            mirnas=mirnas,
            primary_genes=primary_genes,
            lncrna_names=lncrna_names,
            csq_description=csq_description,
            window=window,
            filtering=filtering,
        )
        add_tcga_id_columns_inplace(sv_df, raw_id=sample_id)
        output_path = output_dir / f"{sample_id}_{filtering}_sv_set.csv"
        sv_df.to_csv(output_path, index=False)
        created_files.append(output_path)
        print(f"  Saved: {output_path}")

    return created_files


# =============================================================================
# MOTIF SCANNING
# =============================================================================

def run_sv_motif_scanning(
    sv_csv_dir: Optional[Path] = None,
    output_root: Optional[Path] = None,
    ref_fasta: Optional[Path] = None,
    meme_file: Optional[Path] = None,
    flank_size: Optional[int] = None,
    fimo_threshold: Optional[float] = None,
    threads: int = 4,
    all_ccre_fimo_path=None,      
    ccre_table_path=None,           
    neojunction_window=None,
    chip_dir: Optional[Path] = None,
    chip_unified_path: Optional[Path] = None,
    chip_window: Optional[int] = None,
    chip_celltype_whitelist: Optional[List[str]] = None,
    chip_tf_whitelist: Optional[List[str]] = None,
    skip_chip: bool = False,
    resume_from_fimo: bool = False,
) -> Dict[str, Path]:
    """
    Run motif scanning on processed SVs.
    
    This function orchestrates:
    1. BED interval generation
    2. FASTA extraction (via bedtools)
    3. FIMO motif scanning
    4. FIMO-SV intersection (via bedtools)
    5. Recombination of motif hits
    6. BND neojunction motif annotation (when cCRE paths are configured)
    7. ChIP-seq peak disruption annotation (unless skip_chip)
    
    Args:
        sv_csv_dir: Directory with processed SV CSVs (default: config.PATHS.sv_output_root / "02_processed_sv_csv")
        output_root: Root directory for all outputs (default: config.PATHS.sv_output_root)
        ref_fasta: Path to reference FASTA (default: config.PATHS.sv_reference_fasta)
        meme_file: Path to MEME motif file (default: config.PATHS.sv_meme_file)
        flank_size: Flank size in bp (default: config.THRESHOLDS.sv_flank_size_bp)
        fimo_threshold: FIMO p-value threshold (default: config.THRESHOLDS.fimo_pvalue_threshold)
        threads: Number of threads for parallel operations
        resume_from_fimo: If True, skip steps 1–3 and use existing ``05_fimo_tsv/*_fimo.tsv``
            plus ``03_sv_bed/*.bed``; run steps 4–7 only (merge, recombine, neojunction, ChIP).
    
    Returns:
        Dict with paths to output directories
    """
    # Get defaults from config
    PATHS, THRESHOLDS, _ = _get_config()
    
    output_root = _get_path_with_default(output_root, "sv_output_root")
    if output_root is None:
        raise ValueError("output_root must be provided or configured in config.PATHS.sv_output_root")
    output_root = Path(output_root)
    
    if sv_csv_dir is None:
        sv_csv_dir = output_root / "02_processed_sv_csv"
    sv_csv_dir = Path(sv_csv_dir)
    
    if not resume_from_fimo:
        ref_fasta = _get_path_with_default(ref_fasta, "sv_reference_fasta")
        if ref_fasta is None:
            raise ValueError("ref_fasta must be provided or configured in config.PATHS.sv_reference_fasta")
        ref_fasta = Path(ref_fasta).expanduser()

        meme_file = _get_path_with_default(meme_file, "sv_meme_file")
        if meme_file is None:
            raise ValueError("meme_file must be provided or configured in config.PATHS.sv_meme_file")
        meme_file = Path(meme_file)
    else:
        ref_fasta = None
        meme_file = None

    flank_size = _get_threshold_with_default(flank_size, "sv_flank_size_bp", 150)
    fimo_threshold = _get_threshold_with_default(fimo_threshold, "fimo_pvalue_threshold", 1e-4)
    
    # Output directories
    dirs = {
        "bed": output_root / "03_sv_bed",
        "fasta": output_root / "04_sv_fasta",
        "fimo_tsv": output_root / "05_fimo_tsv",
        "fimo_merged": output_root / "06_sv_fimo_merged",
        "final": output_root / "07_final_sv_with_fimo",
    }
    
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    fimo_skip_stems = _sv_fimo_recombine_skip_csv_stems()

    if resume_from_fimo:
        n_fimo = len(list(dirs["fimo_tsv"].glob("*_fimo.tsv")))
        print(
            "\n[SKIP steps 1–3] resume_from_fimo=True: "
            f"using existing BED/FASTA/FIMO under {output_root} "
            f"({n_fimo} *_fimo.tsv in 05_fimo_tsv)."
        )
        if n_fimo == 0:
            print(f"  [WARN] No *_fimo.tsv files in {dirs['fimo_tsv']}; steps 4–5 will produce little or no output.")
    else:
        # Step 1: Create BED intervals
        print("\n[STEP 1] Creating BED intervals...")
        bed_files = create_beds_from_directory(sv_csv_dir, dirs["bed"], flank_size)
        print(f"  Created {len(bed_files)} BED files")

        # Step 2: Extract FASTA (requires bedtools)
        print("\n[STEP 2] Extracting FASTA sequences...")
        for bed_file in bed_files:
            base = bed_file.stem
            if _sv_csv_stem_from_flanks_bed_stem(base) in fimo_skip_stems:
                print(
                    f"  [SKIP] Step 2–4 for {_sv_csv_stem_from_flanks_bed_stem(base)!r} "
                    "(THRESHOLDS.sv_fimo_recombine_skip_csv_basenames)"
                )
                continue
            out_fa = dirs["fasta"] / f"{base}.fa"

            cmd = [
                "bedtools", "getfasta",
                "-fi", str(ref_fasta),
                "-bed", str(bed_file),
                "-name",
                "-fo", str(out_fa),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                print(f"  Created: {out_fa}")
            except subprocess.CalledProcessError as e:
                print(f"  [ERROR] Failed to extract FASTA for {bed_file}: {e}")

        # Step 3: Run FIMO
        print("\n[STEP 3] Running FIMO motif scanning...")
        for fa_file in dirs["fasta"].glob("*.fa"):
            base = fa_file.stem
            if _sv_csv_stem_from_flanks_bed_stem(base) in fimo_skip_stems:
                continue
            fimo_out_dir = dirs["fimo_tsv"] / f"{base}_fimo_run"
            fimo_out_dir.mkdir(parents=True, exist_ok=True)

            cmd = [
                "fimo",
                "--oc", str(fimo_out_dir),
                "--thresh", str(fimo_threshold),
                str(meme_file),
                str(fa_file),
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                # Move fimo.tsv to expected location
                fimo_tsv = fimo_out_dir / "fimo.tsv"
                if fimo_tsv.exists():
                    out_tsv = dirs["fimo_tsv"] / f"{base}_fimo.tsv"
                    fimo_tsv.rename(out_tsv)
                    print(f"  Created: {out_tsv}")
            except subprocess.CalledProcessError as e:
                print(f"  [ERROR] FIMO failed for {fa_file}: {e}")
            except FileNotFoundError:
                print("  [ERROR] FIMO not found in PATH. Install MEME suite.")
                break

    # Step 4: Intersect FIMO with SV BEDs
    print("\n[STEP 4] Intersecting FIMO hits with SV intervals...")
    for fimo_tsv in dirs["fimo_tsv"].glob("*_fimo.tsv"):
        base = fimo_tsv.stem.replace("_fimo", "")
        if _sv_csv_stem_from_flanks_bed_stem(base) in fimo_skip_stems:
            print(
                f"  [SKIP] Step 4 for {_sv_csv_stem_from_flanks_bed_stem(base)!r} "
                "(THRESHOLDS.sv_fimo_recombine_skip_csv_basenames)"
            )
            continue

        # Find matching BED
        bed_file = dirs["bed"] / f"{base}.bed"
        if not bed_file.exists():
            continue
        
        # Convert FIMO TSV to BED
        from .motif_scanning import parse_fimo_tsv_to_bed
        fimo_bed = dirs["fimo_tsv"] / f"{base}_fimo_hits.bed"
        parse_fimo_tsv_to_bed(fimo_tsv, fimo_bed)
        
        # Intersect
        merged_bed = dirs["fimo_merged"] / f"{base}_fimo_merged.bed"
        cmd = [
            "bedtools", "intersect",
            "-a", str(bed_file),
            "-b", str(fimo_bed),
            "-wa", "-wb",
        ]
        
        try:
            with open(merged_bed, "w") as f:
                subprocess.run(cmd, check=True, stdout=f)
            print(f"  Created: {merged_bed}")
        except subprocess.CalledProcessError as e:
            print(f"  [ERROR] Intersection failed for {base}: {e}")
    
    # Step 5: Recombine motif hits
    print("\n[STEP 5] Recombining motif hits into SV tables...")
    recombine_all_sv_fimo(dirs["fimo_merged"], sv_csv_dir, dirs["final"])


    # Resolve neojunction params
    all_ccre_fimo_path = _get_path_with_default(all_ccre_fimo_path, "all_ccre_fimo_tsv")
    if all_ccre_fimo_path is None:
        raise ValueError("all_ccre_fimo_path must be provided or configured in config.PATHS.all_ccre_fimo_tsv")
    ccre_table_path = _get_path_with_default(ccre_table_path, "ccre_csv")
    if ccre_table_path is None:
        raise ValueError("ccre_table_path must be provided or configured in config.PATHS.ccre_csv")
    neojunction_window = _get_threshold_with_default(
        neojunction_window, "neojunction_window_bp", 500_000
    )
 
    # Step 6: Neojunction motif annotation (BND-specific)
    if all_ccre_fimo_path is not None and ccre_table_path is not None:
        print("\n[STEP 6] Annotating BND neojunction motifs...")
 
        neo_dir = output_root / "08_neojunction_enriched"
        neo_dir.mkdir(parents=True, exist_ok=True)
 
        # Load once (full TSV + overlap pass — large RAM / long CPU; not per-sample)
        print(
            "  Loading cCRE table and genome-wide cCRE FIMO (can take minutes and use many GB RAM)...",
            flush=True,
        )
        ccre_table = pd.read_csv(ccre_table_path)
        ccre_fimo_lookup = load_all_ccre_fimo(all_ccre_fimo_path, ccre_table)
        print("  cCRE FIMO lookup ready.", flush=True)

        final_csvs = sorted(dirs["final"].glob("*.csv"))
        for j, csv_file in enumerate(final_csvs, 1):
            print(f"  [{j}/{len(final_csvs)}] {csv_file.name} ...", flush=True)
            sv_df = pd.read_csv(csv_file)
            n_bnd = (sv_df["SVTYPE"] == "BND").sum() if "SVTYPE" in sv_df.columns else 0
 
            if n_bnd > 0:
                print(f"    {n_bnd} BNDs — annotating neojunction motifs...", flush=True)
                sv_df = annotate_bnd_neojunction_motifs(
                    sv_df, ccre_fimo_lookup, ccre_table, window=neojunction_window,
                )
                # Serialise for CSV
                sv_df["neojunction_motif_hits"] = sv_df["neojunction_motif_hits"].apply(
                    lambda x: "[]" if (x is None or x == []) else repr(x)
                )
 
            out_path = neo_dir / csv_file.name
            sv_df.to_csv(out_path, index=False)
            print(f"  Wrote: {out_path}")
 
        dirs["neojunction"] = neo_dir
    else:
        print("\n[STEP 6] Skipping neojunction annotation "
              "(all_ccre_fimo_path or ccre_table_path not provided)")
    # Step 7: ChIP-seq SV-as-disruptor annotation
    if not skip_chip:
        print("\n[STEP 7] Annotating SVs with ChIP-seq peak disruption...")
 
        # Resolve ChIP params from config defaults
        chip_dir = _get_path_with_default(chip_dir, "chip_dir")
        chip_unified_path = _get_path_with_default(chip_unified_path, "chip_unified")
        chip_window = _get_threshold_with_default(
            chip_window, "sv_chip_window_bp", 1_000_000
        )
 
        # Cell-type whitelist: caller arg > config > None (= keep all)
        if chip_celltype_whitelist is None:
            _, _, _ = _get_config()  # keep pattern; biosamples not in helper
            try:
                from config import BIOSAMPLES
                if hasattr(BIOSAMPLES, "chip_brca_celltypes"):
                    chip_celltype_whitelist = BIOSAMPLES.chip_brca_celltypes
            except ImportError:
                pass
 
        if chip_dir is None and chip_unified_path is None:
            print("  [SKIP] No chip_dir or chip_unified_path provided/configured")
        else:
            # Load: prefer cached unified parquet if present, else build from BEDs
            chip_df = None
            if chip_unified_path is not None and Path(chip_unified_path).exists():
                print(f"  Loading cached unified ChIP table: {chip_unified_path}")
                if Path(chip_unified_path).suffix == ".parquet":
                    chip_df = pd.read_parquet(chip_unified_path)
                else:
                    chip_df = pd.read_csv(chip_unified_path)
            elif chip_dir is not None:
                print(f"  Building unified ChIP table from {chip_dir}")
                chip_df = load_unified_chip(
                    working_dir=chip_dir,
                    output_path=chip_unified_path,  # cache for next run
                )
            else:
                print("  [SKIP] chip_unified_path missing and chip_dir not set")
 
            if chip_df is not None and not chip_df.empty:
                # Source dir for Step 7: prefer neojunction output, else final
                source_dir = dirs.get("neojunction", dirs["final"])
                chip_dir_out = output_root / "09_chip_enriched"
                chip_dir_out.mkdir(parents=True, exist_ok=True)
 
                annotate_sv_csvs_with_chip(
                    sv_csv_dir=source_dir,
                    chip=chip_df,
                    output_dir=chip_dir_out,
                    window=chip_window,
                    cell_type_whitelist=chip_celltype_whitelist,
                    tf_whitelist=chip_tf_whitelist,
                )
                dirs["chip"] = chip_dir_out
            else:
                print("  [SKIP] ChIP table empty or unavailable")
    else:
        print("\n[STEP 7] Skipping ChIP annotation (skip_chip=True)")
 
    return dirs
 
  

# =============================================================================
# FULL PIPELINE
# =============================================================================

def run_sv_pipeline(
    vcf_dir: Optional[Path] = None,
    output_root: Optional[Path] = None,
    genes_path: Optional[Path] = None,
    elements_path: Optional[Path] = None,
    samples_tsv_path: Optional[Path] = None,
    ref_fasta: Optional[Path] = None,
    meme_file: Optional[Path] = None,
    lncrnas_path: Optional[Path] = None,
    csq_description: Optional[str] = None,
    window: Optional[int] = None,
    flank_size: Optional[int] = None,
    fimo_threshold: Optional[float] = None,
    all_ccre_fimo_path: Optional[Path] = None,
    ccre_table_path: Optional[Path] = None,
    neojunction_window: Optional[int] = None,
    chip_dir: Optional[Path] = None,
    chip_unified_path: Optional[Path] = None,
    chip_window: Optional[int] = None,
    chip_celltype_whitelist: Optional[List[str]] = None,
    chip_tf_whitelist: Optional[List[str]] = None,
    filtering: str = "strict",
    skip_vep: bool = False,
    skip_motifs: bool = False,
    skip_chip: bool = False,
    resume_motifs_from_fimo: bool = False,
    vep_conda_env: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the complete SV pipeline.
    
    All path arguments will use defaults from config.py if not provided.
    
    Args:
        vcf_dir: Directory containing VCF files (default: config.PATHS.sv_vcf_dir)
        output_root: Root directory for all outputs (default: config.PATHS.sv_output_root)
        genes_path: Path to gene features CSV (default: config.PATHS.gencode_gtf_csv)
        elements_path: Path to regulatory elements CSV
        samples_tsv_path: Path to sample annotations TSV (default: config.PATHS.sv_samples_tsv)
        ref_fasta: Path to reference FASTA (default: config.PATHS.sv_reference_fasta)
        meme_file: Path to MEME motif file (default: config.PATHS.sv_meme_file)
        lncrnas_path: Optional path to lncRNA features CSV (default: config.PATHS.lncrna_csv)
        csq_description: VEP CSQ format description (default: config.VEP_CSQ_FORMAT)
        window: Maximum distance for mapping (default: config.THRESHOLDS.sv_gene_window_bp)
        flank_size: Flank size in bp for motif scanning (default: config.THRESHOLDS.sv_flank_size_bp)
        fimo_threshold: FIMO p-value threshold (default: config.THRESHOLDS.fimo_pvalue_threshold)
        filtering: "strict" or "lenient"
        skip_vep: Skip VEP annotation step
        skip_motifs: Skip motif scanning step
        resume_motifs_from_fimo: If True (and motifs are not skipped), skip SV VCF/VEP processing
            and reuse ``02_processed_sv_csv``; run motif steps 4–7 only using existing
            ``05_fimo_tsv`` and ``03_sv_bed`` (no new FIMO).
    
    Returns:
        Dict with pipeline results and paths
    """
    # Get defaults from config
    PATHS, THRESHOLDS, VEP_CSQ_FORMAT_DEFAULT = _get_config()
    
    # Resolve all paths with config defaults
    vcf_dir = _get_path_with_default(vcf_dir, "sv_vcf_dir")
    if vcf_dir is None:
        raise ValueError("vcf_dir must be provided or configured in config.PATHS.sv_vcf_dir")
    vcf_dir = Path(vcf_dir)
    
    output_root = _get_path_with_default(output_root, "sv_output_root")
    if output_root is None:
        # Fall back to working_dir / sv_pipeline
        working_dir = _get_path_with_default(None, "working_dir")
        if working_dir is not None:
            output_root = working_dir / "sv_pipeline"
        else:
            raise ValueError("output_root must be provided or configured in config.PATHS.sv_output_root")
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    if vep_conda_env is None:
       PATHS, _, _ = _get_config()
       if PATHS is not None and hasattr(PATHS, "vep_conda_env"):
           vep_conda_env = PATHS.vep_conda_env
       else:
           vep_conda_env = "vep_env"
    
    genes_path = _get_path_with_default(genes_path, "gencode_gtf_csv")
    elements_path = _get_path_with_default(elements_path, "regulatory_elements_table")
    samples_tsv_path = _get_path_with_default(samples_tsv_path, "sv_samples_tsv")
    ref_fasta = _get_path_with_default(ref_fasta, "sv_reference_fasta")
    meme_file = _get_path_with_default(meme_file, "sv_meme_file")
    lncrnas_path = _get_path_with_default(lncrnas_path, "lncrna_csv")
    
    if csq_description is None:
        csq_description = VEP_CSQ_FORMAT_DEFAULT
    
    window = _get_threshold_with_default(window, "sv_gene_window_bp", 1_000_000)
    flank_size = _get_threshold_with_default(flank_size, "sv_flank_size_bp", 150)
    fimo_threshold = _get_threshold_with_default(fimo_threshold, "fimo_pvalue_threshold", 1e-4)
    
    results = {
        "vcf_dir": vcf_dir,
        "output_root": output_root,
        "processed_csvs": [],
        "final_csvs": [],
    }
    
    # Output directories
    proc_csv_dir = output_root / "02_processed_sv_csv"
    proc_csv_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 60)
    print("SV PIPELINE")
    print("=" * 60)
    print(f"  VCF directory: {vcf_dir}")
    print(f"  Output root: {output_root}")
    print(f"  Filtering: {filtering}")
    print(f"  Skip VEP: {skip_vep}")
    print(f"  Skip motifs: {skip_motifs}")
    print(f"  Skip ChIP: {skip_chip}")
    print(f"  Resume motifs from FIMO (skip SV VCF→CSV; motif steps 4–7 only): {resume_motifs_from_fimo}")

    # Step 1: Process VCFs (optionally with VEP), unless resuming motif work from existing FIMO
    if resume_motifs_from_fimo:
        print("\n" + "-" * 40)
        print("STEP 1: Skipped (resume_motifs_from_fimo=True)")
        print("-" * 40)
        pat = f"*_{filtering}_sv_set.csv"
        processed_files = sorted(proc_csv_dir.glob(pat))
        print(f"  Using {len(processed_files)} existing CSV(s) matching {pat!r} under {proc_csv_dir}")
        if not processed_files:
            print(
                "  [WARN] No matching processed SV CSVs under 02_processed_sv_csv; "
                "downstream motif merge/recombine may do nothing."
            )
        results["processed_csvs"] = processed_files
    else:
        print("\n" + "-" * 40)
        print("STEP 1: Processing VCFs")
        print("-" * 40)

        # If VEP annotation is needed and not skipped, run VEP first
        active_vcf_dir = vcf_dir
        if not skip_vep:
            vep_dir = output_root / "01_vep_vcfs"
            vep_dir.mkdir(parents=True, exist_ok=True)

            # Run VEP (assumes vep is in PATH)
            print("  Running VEP annotation...")
            vcfs = sorted(vcf_dir.glob("*.vcf")) + sorted(vcf_dir.glob("*.vcf.gz"))
            for vcf_file in vcfs:
                out_vcf = vep_dir / vcf_file.name
                cmd = [
                    "vep",
                    "--input_file", str(vcf_file),
                    "--output_file", str(out_vcf),
                    "--offline", "--cache",
                    "--fasta", str(ref_fasta),
                    "--species", "homo_sapiens",
                    "--assembly", "GRCh38",
                    "--format", "vcf", "--vcf",
                    "--force_overwrite", "--everything",
                    "--symbol", "--variant_class", "--check_ref",
                    "--fork", "4",
                ]
                try:
                    subprocess.run(cmd, check=True, capture_output=True)
                    print(f"    VEP completed: {vcf_file.name}")
                except (subprocess.CalledProcessError, FileNotFoundError) as e:
                    print(f"    [WARN] VEP failed/not found for {vcf_file.name}")

            active_vcf_dir = vep_dir

        processed_files = run_sv_vcf_processing(
            vcf_dir=active_vcf_dir,
            output_dir=proc_csv_dir,
            genes_path=genes_path,
            elements_path=elements_path,
            samples_tsv_path=samples_tsv_path,
            lncrnas_path=lncrnas_path,
            csq_description=csq_description,
            window=window,
            filtering=filtering,
        )
        results["processed_csvs"] = processed_files
    
    # Step 2: Motif scanning (optional)
    if not skip_motifs:
        print("\n" + "-" * 40)
        print("STEP 2: Motif Scanning")
        print("-" * 40)
        
        motif_dirs = run_sv_motif_scanning(
            sv_csv_dir=proc_csv_dir,
            output_root=output_root,
            ref_fasta=ref_fasta,
            meme_file=meme_file,
            flank_size=flank_size,
            fimo_threshold=fimo_threshold,
            all_ccre_fimo_path=all_ccre_fimo_path,
            ccre_table_path=ccre_table_path,
            neojunction_window=neojunction_window,
            chip_dir=chip_dir,
            chip_unified_path=chip_unified_path,
            chip_window=chip_window,
            chip_celltype_whitelist=chip_celltype_whitelist,
            chip_tf_whitelist=chip_tf_whitelist,
            skip_chip=skip_chip,
            resume_from_fimo=resume_motifs_from_fimo,
        )
        # Prefer ChIP-enriched tables (09_chip_enriched) as the SV "final" deliverable when present.
        if "chip" in motif_dirs:
            results["chip_enriched_dir"] = motif_dirs["chip"]
            results["final_csvs"] = sorted(motif_dirs["chip"].glob("*.csv"))
        elif "neojunction" in motif_dirs:
            results["neojunction_dir"] = motif_dirs["neojunction"]
            results["final_csvs"] = sorted(motif_dirs["neojunction"].glob("*.csv"))
        else:
            results["final_csvs"] = sorted(motif_dirs["final"].glob("*.csv"))
        results["motif_dirs"] = motif_dirs
    else:
        results["final_csvs"] = processed_files
    
    # Summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"  Processed CSVs: {len(results['processed_csvs'])}")
    print(f"  Final CSVs: {len(results['final_csvs'])}")
    print(f"  Output directory: {output_root}")
    
    return results


# =============================================================================
# SUMMARY AND REPORTING
# =============================================================================

def summarize_sv_results(
    final_csv_dir: Path,
) -> pd.DataFrame:
    """
    Generate summary statistics across all processed SV files.
    """
    final_csv_dir = Path(final_csv_dir)
    
    summaries = []
    
    for csv_file in final_csv_dir.glob("*.csv"):
        df = pd.read_csv(csv_file)
        
        sample_id = csv_file.stem.replace("_strict_sv_set", "").replace("_lenient_sv_set", "")
        
        summary = {
            "sample_id": sample_id,
            "n_svs": len(df),
            "n_del": (df["SVTYPE"] == "DEL").sum() if "SVTYPE" in df.columns else 0,
            "n_dup": (df["SVTYPE"] == "DUP").sum() if "SVTYPE" in df.columns else 0,
            "n_ins": (df["SVTYPE"] == "INS").sum() if "SVTYPE" in df.columns else 0,
            "n_bnd": (df["SVTYPE"] == "BND").sum() if "SVTYPE" in df.columns else 0,
            "n_with_gene_hits": (df["gene_hits"].apply(lambda x: len(x) > 2 if isinstance(x, str) else False)).sum() if "gene_hits" in df.columns else 0,
            "n_with_elem_hits": (df["elem_hits"].apply(lambda x: len(x) > 2 if isinstance(x, str) else False)).sum() if "elem_hits" in df.columns else 0,
            "n_with_motif_hits": (df["flank_motif_hits"].apply(lambda x: len(x) > 2 if isinstance(x, str) else False)).sum() if "flank_motif_hits" in df.columns else 0,
            "n_with_neojunction_hits": (df["neojunction_motif_hits"].apply(lambda x: len(x) > 2 if isinstance(x, str) else False)).sum() if "neojunction_motif_hits" in df.columns else 0,
            "n_with_chip_hits": (df["chip_hits"].apply(lambda x: len(x) > 2 if isinstance(x, str) else False)).sum() if "chip_hits" in df.columns else 0,
        }
        
        summaries.append(summary)
    
    return pd.DataFrame(summaries)