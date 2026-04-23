"""
Mutect2 VCF loading and processing.

Main entry point for loading VEP-annotated Mutect2 VCF files.
Handles:
- VCF parsing with vcfpy
- Multi-allelic expansion (one row per ALT)
- Per-sample AD/DP extraction
- VAF computation
- VEP annotation parsing
- Somatic filtering
- cCRE overlap matching
"""

import io
import json
import os
import subprocess
from pathlib import Path
from typing import Tuple, List, Optional, Union
import glob

import numpy as np
import pandas as pd

try:
    import vcfpy
except ImportError:
    vcfpy = None

from .somatic_filter import high_conf_somatic_mask
from .vep_parser import add_vep_hits_columns, DEFAULT_CSQ_FORMAT, _empty_vep_result
from .ccre_matching import match_snvs_to_ccres
from .mirna_matching import match_snvs_to_mirnas


# =============================================================================
# ATLAS-ALIGNED COLUMN DEFAULTS
# =============================================================================

# Columns that ``add_vep_hits_columns`` (vep_parser) adds when CSQ is present.
# We mirror them here so the SNV table always exposes a stable atlas-aligned
# schema — even for non-VEP input VCFs (downstream code can then branch on
# ``gene_symbols == ''`` instead of missing-column errors).
SNV_VEP_ATLAS_COLUMNS: Tuple[str, ...] = (
    "gene_hits",
    "regulatory_hits",
    "motif_hits",
    "gene_symbols",
    "hits_canonical",
    "has_missense",
    "has_nonsense",
    "has_frameshift",
    "has_splice_effect",
    "has_missense_canonical",
    "has_nonsense_canonical",
    "has_frameshift_canonical",
    "has_splice_effect_canonical",
    "has_missense_mane",
    "has_nonsense_mane",
    "has_frameshift_mane",
    "has_splice_effect_mane",
)


def _ensure_snv_atlas_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Guarantee the atlas-required VEP columns exist on the SNV table.

    Non-destructive: existing values are left untouched. Missing columns are
    filled with the empty defaults defined by ``vep_parser._empty_vep_result``.
    """
    if df.empty:
        for col in SNV_VEP_ATLAS_COLUMNS:
            if col not in df.columns:
                df[col] = pd.Series(dtype=object)
        return df
    empty = _empty_vep_result()
    for col in SNV_VEP_ATLAS_COLUMNS:
        if col in df.columns:
            continue
        default = empty.get(col)
        if isinstance(default, list):
            df[col] = [[] for _ in range(len(df))]
        else:
            df[col] = default
    return df


# =============================================================================
# CONFIG HELPERS
# =============================================================================

def _get_default_elements_path() -> Optional[Path]:
    """Get default elements path from config if available."""
    try:
        from ..config import PATHS
        if hasattr(PATHS, 'regulatory_elements_table'):
            return PATHS.regulatory_elements_table
    except ImportError:
        pass
    return None


def _get_default_mirna_path() -> Optional[Path]:
    """Get default miRNA features path from config if available."""
    try:
        from ..config import PATHS
        # Prefer mature-arm loci table when present for true arm-resolution by coordinates.
        if hasattr(PATHS, "mirna_mature_loci_csv"):
            return PATHS.mirna_mature_loci_csv
        if hasattr(PATHS, "mirna_path"):
            return PATHS.mirna_path
    except ImportError:
        pass
    return None

def _get_default_output_dir() -> Optional[Path]:
    """Get default output directory from config if available."""
    try:
        from ..config import PATHS, OUTPUT_SUBDIRS
        if hasattr(PATHS, 'working_dir') and 'snv' in OUTPUT_SUBDIRS:
            return PATHS.working_dir / OUTPUT_SUBDIRS['snv']
    except ImportError:
        pass
    return None


def _get_default_primary_genes() -> List[str]:
    """Default gene panel for SNV annotation (matches ``PIPELINE_GENE_PANEL``)."""
    try:
        from ..config import PIPELINE_GENE_PANEL

        return list(PIPELINE_GENE_PANEL)
    except ImportError:
        return []


def _extract_primary_sample_id_from_samples_tsv(
    samples_tsv: Union[str, Path],
    file_name: str,
) -> Optional[str]:
    """
    Resolve a canonical TCGA sample ID for a VCF using a GDC-style samples.tsv.

    Expected columns:
      - File Name
      - Tumor Descriptor
      - Sample ID  (often "TCGA-...-01A, TCGA-...-10A")
    Returns:
      Tumor sample id (usually first entry, unless Tumor Descriptor begins with 'Not').
    """
    try:
        samples_tsv = Path(samples_tsv)
        ann = pd.read_csv(samples_tsv, sep="\t")
        if "File Name" not in ann.columns or "Sample ID" not in ann.columns:
            return None
        row = ann[ann["File Name"] == file_name]
        if row.empty:
            return None
        row = row.iloc[0]
        tumor_desc = str(row.get("Tumor Descriptor", ""))
        sample_id_field = str(row.get("Sample ID", ""))
        parts = [p.strip() for p in sample_id_field.split(",") if str(p).strip()]
        if not parts:
            return None
        if tumor_desc.startswith("Not") and len(parts) >= 2:
            return parts[1]
        return parts[0]
    except Exception:
        return None


def _clean_vcf_header(vcf_path: Union[str, Path]) -> io.StringIO:
    """
    Clean VCF header by removing non-standard comment lines.
    
    Some VEP output includes non-standard header lines (e.g., "## ENSEMBL VARIANT...")
    that vcfpy cannot parse. This function strips them.
    """
    cleaned_lines = []
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                # Only keep proper "##key=value" meta lines
                if "=" in line:
                    cleaned_lines.append(line)
                continue
            
            # First non-## line we care about is "#CHROM ..."
            if line.startswith("#CHROM"):
                cleaned_lines.append(line)
                break
        
        # Copy remaining variant records
        for line in f:
            cleaned_lines.append(line)
    
    return io.StringIO("".join(cleaned_lines))


def load_mutect_snv_vcf(
    vcf_path: Union[str, Path],
    primary_genes: Optional[List[str]] = None,
    elements_path: Optional[Union[str, Path]] = None,
    elements_df: Optional[pd.DataFrame] = None,
    mirna_df: Optional[pd.DataFrame] = None,
    csq_description: Optional[str] = None,
    clean_header: bool = True,
    apply_filter: bool = True,
    filter_kwargs: Optional[dict] = None,
    output_dir: Optional[Union[str, Path]] = None,
    save_outputs: bool = False,
    run_fimo: bool = False,
    fimo_work_dir: Optional[Union[str, Path]] = None,
    ref_fasta: Optional[Union[str, Path]] = None,
    meme_file: Optional[Union[str, Path]] = None,
    fimo_flank_bp: Optional[int] = None,
    fimo_threshold: Optional[float] = None,
    max_fimo_hits_per_variant: Optional[int] = None,
    run_chip: bool = True,
    chip_unified_path: Optional[Union[str, Path]] = None,
) -> Tuple[pd.DataFrame, Optional[str], Optional[str]]:
    """
    Load a Mutect2-style SNV/indel VCF (optionally VEP-annotated).
    
    Args:
        vcf_path: Path to VCF file
        primary_genes: List of gene symbols to filter VEP annotations.
                      If None, uses ``PIPELINE_GENE_PANEL`` from config.
        elements_path: Path to regulatory elements CSV (for cCRE matching).
                      If None, uses PATHS.regulatory_elements_table from config.
        elements_df: Pre-loaded elements DataFrame (alternative to elements_path)
        mirna_df: Pre-loaded miRNA features DataFrame (alternative to mirna_path)
        csq_description: VEP CSQ format string (auto-detected if None)
        clean_header: Whether to clean non-standard VCF header lines
        apply_filter: Whether to apply somatic filtering (uses THRESHOLDS from config)
        filter_kwargs: Override specific filter thresholds (see high_conf_somatic_mask)
        output_dir: Directory for saving outputs. If None and save_outputs=True,
                   uses PATHS.working_dir/snv from config.
        save_outputs: Whether to save DataFrame and summary to output_dir
        run_fimo: If True, run MEME FIMO on a fixed-width **reference** window around each
                  retained variant (requires ``bedtools``; ``fimo`` via ``resolve_fimo_argv()``).
                  Adds ``fimo_hits`` (list of dicts); see ``SNV.snv_fimo.annotate_snvs_with_fimo``.
        fimo_work_dir: Directory for BED/FASTA/FIMO artifacts. Default: ``output_dir/snv_fimo/<vcf_stem>``
                       or a temp-tree under ``/tmp`` when ``output_dir`` is None.
        ref_fasta: Reference FASTA (default ``PATHS.sv_reference_fasta``).
        meme_file: MEME motif file (default ``PATHS.sv_meme_file``).
        fimo_flank_bp: Half-width in bp (default ``THRESHOLDS.snv_fimo_flank_bp``).
        fimo_threshold: ``fimo --thresh`` (default ``THRESHOLDS.fimo_pvalue_threshold``).
        max_fimo_hits_per_variant: Cap per row (default ``THRESHOLDS.snv_fimo_max_hits_per_variant``).
        run_chip: If True, overlap variants with ``PATHS.chip_unified`` (strict POS overlap; no window).
        chip_unified_path: Override unified ChIP parquet path (default ``PATHS.chip_unified``).
    
    Returns:
        Tuple of (DataFrame, normal_sample_name, tumor_sample_name)
        
    DataFrame contains:
        - Basic variant info: chrom, pos, id, ref, alt, qual, filter
        - INFO fields flattened
        - Per-sample AD/DP values
        - Computed VAFs: normal_vaf, tumor_vaf
        - VEP parsed columns: gene_hits, regulatory_hits, motif_hits, etc.
        - cCRE_hits: list of overlapping regulatory elements
        - snv_chip_hits / snv_chip_aggregate: ChIP peak overlap (see ``SNV.snv_chip``)
    
    Config integration:
        - Filter thresholds: THRESHOLDS.snv_min_tumor_vaf, snv_max_normal_vaf, etc.
        - Default elements: PATHS.regulatory_elements_table
        - Default genes: ``PIPELINE_GENE_PANEL``
        - Default output: PATHS.working_dir / OUTPUT_SUBDIRS['snv']
    """
    if vcfpy is None:
        raise ImportError("vcfpy is required for VCF loading. Install with: pip install vcfpy")
    
    vcf_path = Path(vcf_path)
    
    # Apply config defaults
    if primary_genes is None:
        primary_genes = _get_default_primary_genes()
        if not primary_genes:
            raise ValueError(
                "primary_genes must be provided or PIPELINE_GENE_PANEL must be set in config"
            )
    
    # Load elements if needed for cCRE matching
    elements = None
    if elements_df is not None:
        elements = elements_df
    else:
        from ..regulatory_elements import load_regulatory_element_focus
        if elements_path is not None:
            elements = load_regulatory_element_focus(elements_path)
        else:
            default_elements = _get_default_elements_path()
            if default_elements is not None:
                try:
                    elements = load_regulatory_element_focus(default_elements)
                    print(f"Using default elements from config: {default_elements}")
                except FileNotFoundError:
                    elements = None
    
    # Apply default output_dir if saving
    if save_outputs and output_dir is None:
        output_dir = _get_default_output_dir()
    
    # Parse VCF
    if clean_header:
        buffer = _clean_vcf_header(vcf_path)
        reader = vcfpy.Reader.from_stream(buffer)
    else:
        reader = vcfpy.Reader.from_path(str(vcf_path))
    
    samples = reader.header.samples.names
    print(f"Samples in VCF: {samples}")
    
    # Identify normal/tumor samples
    normal_sample, tumor_sample = _identify_samples(samples)
    print(f"Assuming normal sample: {normal_sample}")
    print(f"Assuming tumor sample:  {tumor_sample}")
    
    # Extract CSQ description from header if not provided
    if csq_description is None:
        csq_description = _extract_csq_description(reader.header)
    
    # Parse variant records
    rows = _parse_vcf_records(reader, samples)
    df = pd.DataFrame(rows)
    
    if df.empty:
        print("Warning: No variants found in VCF")
        return df, normal_sample, tumor_sample
    
    # Convert numeric columns
    df = _convert_numeric_columns(df, samples)
    
    # Compute VAFs
    df = _compute_vafs(df, normal_sample, tumor_sample)
    
    # Parse VEP annotations
    if "CSQ" in df.columns:
        print("Parsing VEP annotations...")
        gene_symbol_mapping = None
        if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        ):
            from ..genes.symbol_normalization import default_symbol_mapping

            gene_symbol_mapping = default_symbol_mapping()
        df = add_vep_hits_columns(
            df,
            csq_description,
            primary_genes,
            gene_symbol_mapping=gene_symbol_mapping,
        )
    else:
        print(
            f"  [WARN] {vcf_path.name}: no CSQ in INFO — this VCF is NOT VEP-annotated. "
            "gene_hits / regulatory_hits / motif_hits and has_* flags will be empty. "
            "Feed a *.vep.vcf file to enable full atlas columns."
        )

    # Always ensure atlas-required columns exist (empty defaults when CSQ absent)
    df = _ensure_snv_atlas_columns(df)
    
    # Drop exact duplicates
    df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt"])
    
    # Apply somatic filtering
    if apply_filter:
        print("Applying somatic filter...")
        filter_kwargs = filter_kwargs or {}
        
        # Detect sample-specific depth columns
        tumor_dp_col = _find_dp_column(df, tumor_sample)
        normal_dp_col = _find_dp_column(df, normal_sample)
        
        if tumor_dp_col:
            filter_kwargs.setdefault("tumor_dp_col", tumor_dp_col)
        if normal_dp_col:
            filter_kwargs.setdefault("normal_dp_col", normal_dp_col)
        
        mask = high_conf_somatic_mask(df, **filter_kwargs)
        n_before = len(df)
        df = df[mask].copy()
        print(f"  Retained {len(df)}/{n_before} variants after filtering")
    
    # Match to cCREs
    if elements is not None:
        print("Matching variants to cCREs...")
        df = match_snvs_to_ccres(df, elements)
        n_with_hits = (df["cCRE_hits"].apply(len) > 0).sum()
        print(f"  {n_with_hits}/{len(df)} variants overlap cCREs")
    

    if mirna_df is not None:
        print("Matching variants to miRNAs...")
        df = match_snvs_to_mirnas(df, mirna_df)
        n_with_hits = (df["mirna_hits"].apply(len) > 0).sum()
        print(f"  {n_with_hits}/{len(df)} variants overlap miRNAs")

    if run_chip and len(df):
        from ..config import PATHS
        from .snv_chip import annotate_snvs_with_chip

        cup = Path(chip_unified_path).expanduser() if chip_unified_path else PATHS.chip_unified
        print(f"Matching variants to ChIP peaks (strict overlap, {cup})...")
        try:
            df = annotate_snvs_with_chip(df, cup)
            n_chip = (df["snv_chip_hits"].apply(len) > 0).sum()
            print(f"  {n_chip}/{len(df)} variants overlap ≥1 ChIP peak")
        except Exception as e:
            print(f"  [WARN] SNV ChIP overlap failed ({e}); writing empty snv_chip_* columns")
            df["snv_chip_hits"] = [[] for _ in range(len(df))]
            df["snv_chip_aggregate"] = [{"by_tf_source_stratum": []} for _ in range(len(df))]
    elif len(df) and "snv_chip_hits" not in df.columns:
        df["snv_chip_hits"] = [[] for _ in range(len(df))]
        df["snv_chip_aggregate"] = [{"by_tf_source_stratum": []} for _ in range(len(df))]

    if run_fimo and len(df):
        from ..config import PATHS, THRESHOLDS
        from .snv_fimo import annotate_snvs_with_fimo

        flank = int(
            fimo_flank_bp
            if fimo_flank_bp is not None
            else THRESHOLDS.snv_fimo_flank_bp
        )
        max_h = int(
            max_fimo_hits_per_variant
            if max_fimo_hits_per_variant is not None
            else THRESHOLDS.snv_fimo_max_hits_per_variant
        )
        thresh = float(
            fimo_threshold
            if fimo_threshold is not None
            else THRESHOLDS.fimo_pvalue_threshold
        )
        rf = Path(ref_fasta or PATHS.sv_reference_fasta).expanduser()
        mf = Path(meme_file or PATHS.sv_meme_file)
        if fimo_work_dir is None:
            stem = vcf_path.stem.replace(".vep", "").replace(".vcf", "")
            if output_dir is not None:
                wd = Path(output_dir) / "snv_fimo" / stem
            else:
                import tempfile

                wd = Path(tempfile.gettempdir()) / "apm_snv_fimo" / stem
        else:
            wd = Path(fimo_work_dir)
        print(f"Running SNV FIMO (work_dir={wd}, flank={flank} bp, thresh={thresh})...")
        try:
            df = annotate_snvs_with_fimo(
                df,
                ref_fasta=rf,
                meme_file=mf,
                work_dir=wd,
                flank_bp=flank,
                fimo_threshold=thresh,
                max_hits_per_variant=max_h,
            )
            nhit = (df["fimo_hits"].apply(len) > 0).sum()
            print(f"  {nhit}/{len(df)} variants have ≥1 FIMO hit")
        except FileNotFoundError as e:
            print(f"  [WARN] SNV FIMO skipped: {e}")
            df["fimo_hits"] = [[] for _ in range(len(df))]
        except subprocess.CalledProcessError as e:
            print(f"  [WARN] SNV FIMO failed ({e}); writing empty fimo_hits")
            df["fimo_hits"] = [[] for _ in range(len(df))]

    
    # Save outputs if requested
    if save_outputs and output_dir is not None:
        # Extract sample prefix from filename
        prefix = vcf_path.stem.replace(".vep", "").replace(".vcf", "")
        print(f"Saving outputs to {output_dir}...")
        save_snv_outputs(df, output_dir, prefix=prefix)
    
    return df, normal_sample, tumor_sample


def _identify_samples(samples: List[str]) -> Tuple[Optional[str], Optional[str]]:
    """Identify normal and tumor samples from sample names."""
    normal_sample = None
    tumor_sample = None
    
    for sname in samples:
        su = sname.upper()
        if "NORMAL" in su or "CONTROL" in su:
            normal_sample = sname
        if "TUMOR" in su or "TUMOUR" in su:
            tumor_sample = sname
    
    # Fallback to position
    if normal_sample is None and len(samples) >= 1:
        normal_sample = samples[0]
    if tumor_sample is None and len(samples) >= 2:
        tumor_sample = samples[-1]
    
    return normal_sample, tumor_sample


def _extract_csq_description(header) -> str:
    """Extract CSQ format description from VCF header."""
    for info_line in header.info_ids():
        if info_line == "CSQ":
            info = header.get_info_field_info("CSQ")
            if hasattr(info, "description"):
                return info.description
    
    # Return default if not found
    return DEFAULT_CSQ_FORMAT


def _parse_vcf_records(reader, samples: List[str]) -> List[dict]:
    """Parse all VCF records into list of row dicts."""
    rows = []
    
    for rec in reader:
        chrom = rec.CHROM
        pos = rec.POS
        vid = rec.ID[0] if rec.ID else None
        ref = rec.REF
        qual = rec.QUAL
        filt = ";".join(rec.FILTER) if rec.FILTER else "PASS"
        
        # Flatten INFO
        info_flat = {}
        for key, value in rec.INFO.items():
            if isinstance(value, list) and len(value) == 1:
                info_flat[key] = value[0]
            else:
                info_flat[key] = value
        
        # Expand multi-allelic: one row per ALT
        for alt_idx, alt_obj in enumerate(rec.ALT):
            alt_allele = getattr(alt_obj, "value", str(alt_obj))
            
            row = {
                "chrom": chrom,
                "pos": pos,
                "id": vid,
                "ref": ref,
                "alt": alt_allele,
                "alt_index": alt_idx,
                "qual": qual,
                "filter": filt,
            }
            row.update(info_flat)
            
            # Per-sample FORMAT fields
            for sample in samples:
                call = rec.call_for_sample[sample]
                data = call.data
                
                ad = data.get("AD")
                dp = data.get("DP")
                
                if ad is not None and len(ad) > 0:
                    ad_ref = ad[0]
                    ad_alt = ad[alt_idx + 1] if len(ad) > (alt_idx + 1) else np.nan
                else:
                    ad_ref = np.nan
                    ad_alt = np.nan
                
                row[f"{sample}_AD_ref"] = ad_ref
                row[f"{sample}_AD_alt"] = ad_alt
                row[f"{sample}_DP"] = dp if dp is not None else np.nan
            
            rows.append(row)
    
    return rows


def _convert_numeric_columns(df: pd.DataFrame, samples: List[str]) -> pd.DataFrame:
    """Convert numeric columns to appropriate dtypes."""
    for sample in samples:
        for col in [f"{sample}_AD_ref", f"{sample}_AD_alt", f"{sample}_DP"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _compute_vafs(
    df: pd.DataFrame,
    normal_sample: Optional[str],
    tumor_sample: Optional[str],
) -> pd.DataFrame:
    """Compute VAFs for normal and tumor samples."""
    if normal_sample is not None:
        ref_col = f"{normal_sample}_AD_ref"
        alt_col = f"{normal_sample}_AD_alt"
        if ref_col in df.columns and alt_col in df.columns:
            ref = df[ref_col]
            alt = df[alt_col]
            df["normal_vaf"] = alt / (ref + alt)
        else:
            df["normal_vaf"] = np.nan
    else:
        df["normal_vaf"] = np.nan
    
    if tumor_sample is not None:
        ref_col = f"{tumor_sample}_AD_ref"
        alt_col = f"{tumor_sample}_AD_alt"
        if ref_col in df.columns and alt_col in df.columns:
            ref = df[ref_col]
            alt = df[alt_col]
            df["tumor_vaf"] = alt / (ref + alt)
        else:
            df["tumor_vaf"] = np.nan
    else:
        df["tumor_vaf"] = np.nan
    
    return df


def _find_dp_column(df: pd.DataFrame, sample: Optional[str]) -> Optional[str]:
    """Find depth column for a sample."""
    if sample is None:
        return None
    
    candidate = f"{sample}_DP"
    if candidate in df.columns:
        return candidate
    
    return None


# =============================================================================
# BATCH LOADING
# =============================================================================

def load_mutect_snv_batch(
    vcf_dir: Union[str, Path],
    primary_genes: Optional[List[str]] = None,
    elements_path: Optional[Union[str, Path]] = None,
    mirna_path: Optional[Union[str, Path]] = None,
    pattern: str = "*.vcf*",
    output_dir: Optional[Union[str, Path]] = None,
    samples_tsv: Optional[Union[str, Path]] = None,
    save_per_sample: bool = False,
    save_combined: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Load multiple VEP-annotated Mutect2 VCF files from a directory.
    
    Args:
        vcf_dir: Directory containing VCF files
        primary_genes: List of gene symbols to filter VEP annotations.
                      If None, uses ``PIPELINE_GENE_PANEL`` from config.
        elements_path: Path to regulatory elements CSV.
                      If None, uses PATHS.regulatory_elements_table from config.
        mirna_path: Path to miRNA features CSV.
                      If None, uses PATHS.mirna_csv from config.
        pattern: Glob pattern for VCF files (default: "*.vcf")
        output_dir: Directory for saving outputs.
                   If None, uses PATHS.working_dir/snv from config.
        save_per_sample: Save individual sample outputs
        save_combined: Save combined DataFrame
        **kwargs: Additional arguments passed to load_mutect_snv_vcf
    
    Returns:
        Combined DataFrame with all variants and a 'source_file' column
    
    Config integration:
        Uses THRESHOLDS, PATHS, and PIPELINE_GENE_PANEL from pipeline.config
    """
    vcf_dir = Path(vcf_dir)
    vcf_files = sorted(vcf_dir.glob(pattern))
    
    if not vcf_files:
        print(f"No VCF files found matching {pattern} in {vcf_dir}")
        return pd.DataFrame()
    
    print(f"Found {len(vcf_files)} VCF files")
    
    # Apply config defaults
    if primary_genes is None:
        primary_genes = _get_default_primary_genes()
        if not primary_genes:
            raise ValueError(
                "primary_genes must be provided or PIPELINE_GENE_PANEL must be set in config"
            )
    
    if output_dir is None and (save_per_sample or save_combined):
        output_dir = _get_default_output_dir()
    
    # Setup output directory
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        per_sample_dir = output_dir / "per_sample" if save_per_sample else None
        if per_sample_dir:
            per_sample_dir.mkdir(exist_ok=True)
    else:
        per_sample_dir = None
    
    # Load elements and miRNA features once - use config default if not provided
    elements_df = None
    mirna_df = None
    from ..regulatory_elements import load_regulatory_element_focus
    if elements_path is not None:
        elements_df = load_regulatory_element_focus(elements_path)
    else:
        default_elements = _get_default_elements_path()
        if default_elements is not None:
            try:
                elements_df = load_regulatory_element_focus(default_elements)
                print(f"Using default elements from config: {default_elements}")
            except FileNotFoundError:
                elements_df = None
    
    if mirna_path is not None:
        mirna_df = pd.read_csv(mirna_path)
    else:
        default_mirna = _get_default_mirna_path()
        if default_mirna is not None and default_mirna.exists():
            print(f"Using default miRNA features from config: {default_mirna}")
            mirna_df = pd.read_csv(default_mirna)
    
    all_dfs = []
    sample_info = []

    for vcf_path in vcf_files:
        print(f"\nProcessing: {vcf_path.name}")
        try:
            df, normal, tumor = load_mutect_snv_vcf(
                vcf_path,
                primary_genes=primary_genes,
                elements_df=elements_df,
                mirna_df=mirna_df,
                **kwargs,
            )

            if not df.empty:
                df["source_file"] = vcf_path.name
                all_dfs.append(df)
                sample_info.append({
                    "file": vcf_path.name,
                    "normal_sample": normal,
                    "tumor_sample": tumor,
                    "n_variants": len(df),
                })

                if save_per_sample and per_sample_dir is not None:
                    sample_id = None
                    if samples_tsv is not None:
                        sample_id = _extract_primary_sample_id_from_samples_tsv(samples_tsv, vcf_path.name)

                    if sample_id:
                        out_dir = per_sample_dir / sample_id
                        save_snv_outputs(df, out_dir, prefix="")
                    else:
                        prefix = vcf_path.stem.replace(".vep", "").replace(".vcf", "")
                        save_snv_outputs(df, per_sample_dir, prefix=prefix)
        except Exception as e:
            print(f"  Error processing {vcf_path.name}: {e}")
            continue
    
    if not all_dfs:
        print("No variants loaded from any file")
        return pd.DataFrame()
    
    combined = pd.concat(all_dfs, ignore_index=True)
    
    # Print summary
    print(f"\n{'='*50}")
    print(f"Batch loading complete:")
    print(f"  Files processed: {len(sample_info)}")
    print(f"  Total variants: {len(combined)}")
    
    # Save combined output
    if save_combined and output_dir is not None:
        print(f"\nSaving combined outputs to {output_dir}...")
        save_snv_outputs(combined, output_dir, prefix="combined")
        
        # Also save sample manifest
        import json
        manifest_path = output_dir / "sample_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(sample_info, f, indent=2)
        print(f"  Saved sample manifest: {manifest_path}")
    
    return combined


# =============================================================================
# SAVE FUNCTIONS
# =============================================================================

def _dataframe_for_parquet_snv(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy safe for ``to_parquet`` / PyArrow.

    VCF INFO-derived columns (e.g. ``UNITIGS``) sometimes mix Python lists with
    scalars in the same column, which raises ArrowInvalid. Columns with that
    pattern are JSON-serialized per cell; other columns are unchanged.
    """
    pdf = df.copy()
    for col in pdf.columns:
        ser = pdf[col]
        if ser.dtype != object:
            continue
        non_na = ser.dropna()
        if non_na.empty:
            continue
        is_container = non_na.map(lambda x: isinstance(x, (list, tuple, dict)))
        if not (bool(is_container.any()) and bool((~is_container).any())):
            continue

        def _norm_cell(v: object) -> Optional[str]:
            if v is None:
                return None
            if isinstance(v, (list, tuple, dict)):
                return json.dumps(v, default=str)
            try:
                if pd.isna(v):
                    return None
            except (ValueError, TypeError):
                pass
            return str(v)

        pdf[col] = ser.map(_norm_cell)

    return pdf


def save_snv_outputs(
    df: pd.DataFrame,
    output_dir: Union[str, Path],
    prefix: str = "",
    save_parquet: bool = True,
    save_csv: bool = True,
    save_summary: bool = True,
) -> dict:
    """
    Save SNV DataFrame and summary statistics.
    
    Args:
        df: SNV DataFrame to save
        output_dir: Directory for output files
        prefix: Optional prefix for filenames (e.g., sample ID)
        save_parquet: Save as parquet when pyarrow can serialize the frame (skipped on engine/schema errors)
        save_csv: Save as CSV (nested dicts saved as strings, 
                  recoverable with ast.literal_eval)
        save_summary: Save summary statistics JSON
    
    Returns:
        Dict with paths to saved files
    
    Note:
        To reload CSV with nested structures:
        ```python
        import ast
        df = pd.read_csv(path, converters={
            "gene_hits": ast.literal_eval,
            "regulatory_hits": ast.literal_eval,
            "motif_hits": ast.literal_eval,
            "cCRE_hits": ast.literal_eval,
            "mirna_hits": ast.literal_eval,
            "fimo_hits": ast.literal_eval,
            "snv_chip_hits": ast.literal_eval,
            "snv_chip_aggregate": ast.literal_eval,
        })
        ```
    """
    import json
    from ..sample_ids import add_tcga_id_columns_inplace
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Ensure every saved SNV table carries normalized TCGA join keys.
    # Prefer an existing aliquot/sample_id column if present; else infer from output_dir name.
    inferred = output_dir.name if output_dir.name.startswith("TCGA-") else None
    if "aliquot_id" in df.columns:
        add_tcga_id_columns_inplace(df, raw_id_col="aliquot_id")
    elif "sample_id" in df.columns:
        add_tcga_id_columns_inplace(df, raw_id_col="sample_id")
    elif inferred is not None:
        add_tcga_id_columns_inplace(df, raw_id=inferred)
    
    # Build filename prefix
    file_prefix = f"{prefix}_" if prefix else ""
    saved_files = {}
    
    # Save parquet (preserves list/dict columns natively); optional if pyarrow/fastparquet installed
    if save_parquet:
        parquet_path = output_dir / f"{file_prefix}snv_variants.parquet"
        try:
            _dataframe_for_parquet_snv(df).to_parquet(parquet_path, index=False)
            saved_files["parquet"] = parquet_path
            print(f"  Saved parquet: {parquet_path}")
        except ImportError as e:
            print(f"  [WARN] Skipping parquet (no parquet engine): {e}")
        except Exception as e:
            # Last resort if Arrow still rejects the schema (e.g. new odd INFO types).
            print(f"  [WARN] Skipping parquet; CSV still written ({type(e).__name__}): {e}")
    
    # Save CSV (nested dicts become string repr, use ast.literal_eval to recover)
    if save_csv:
        csv_path = output_dir / f"{file_prefix}snv_variants.csv"
        df.to_csv(csv_path, index=False)
        saved_files["csv"] = csv_path
        print(f"  Saved CSV: {csv_path}")
    
    # Save summary statistics
    if save_summary:
        summary = _compute_snv_summary(df)
        summary_path = output_dir / f"{file_prefix}snv_summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        saved_files["summary"] = summary_path
        print(f"  Saved summary: {summary_path}")
    
    return saved_files


def load_snv_csv(
    csv_path: Union[str, Path],
    nested_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Load SNV CSV with nested dict/list columns properly parsed.
    
    Args:
        csv_path: Path to CSV file
        nested_columns: Columns containing nested structures.
                       Defaults to standard SNV nested columns.
    
    Returns:
        DataFrame with nested structures restored
    """
    import ast
    
    if nested_columns is None:
        nested_columns = [
            "gene_hits",
            "regulatory_hits",
            "motif_hits",
            "cCRE_hits",
            "mirna_hits",
            "fimo_hits",
            "snv_chip_hits",
            "snv_chip_aggregate",
        ]
    
    def safe_literal_eval(x):
        if pd.isna(x) or x == "":
            return []
        if isinstance(x, (list, dict)):
            return x
        try:
            return ast.literal_eval(str(x))
        except (ValueError, SyntaxError):
            return []

    header_cols = set(pd.read_csv(csv_path, nrows=0).columns)
    nested_present = [c for c in nested_columns if c in header_cols]
    converters = {col: safe_literal_eval for col in nested_present}

    return pd.read_csv(csv_path, converters=converters)


def _compute_snv_summary(df: pd.DataFrame) -> dict:
    """Compute summary statistics for SNV DataFrame."""
    summary = {
        "n_total_variants": len(df),
        "n_unique_positions": df[["chrom", "pos"]].drop_duplicates().shape[0] if len(df) > 0 else 0,
    }
    
    # Consequence counts
    for col in ["has_missense", "has_nonsense", "has_frameshift", "has_splice_effect"]:
        if col in df.columns:
            summary[f"n_{col.replace('has_', '')}"] = int(df[col].sum())
    
    # Canonical/MANE counts
    for col in ["hits_canonical", "has_missense_canonical", "has_nonsense_canonical"]:
        if col in df.columns:
            summary[f"n_{col}"] = int(df[col].sum())
    
    # Gene symbols affected
    if "gene_symbols" in df.columns:
        all_genes = set()
        for genes_str in df["gene_symbols"].dropna():
            if genes_str:
                all_genes.update(g.strip() for g in str(genes_str).split(",") if g.strip())
        summary["genes_affected"] = sorted(all_genes)
        summary["n_genes_affected"] = len(all_genes)
    
    # cCRE hit counts
    if "cCRE_hits" in df.columns:
        n_with_ccre = (df["cCRE_hits"].apply(lambda x: len(x) if isinstance(x, list) else 0) > 0).sum()
        summary["n_variants_in_ccres"] = int(n_with_ccre)
        
        # cCRE type breakdown
        ccre_types = {}
        for hits in df["cCRE_hits"]:
            if isinstance(hits, list):
                for hit in hits:
                    etype = hit.get("elem_type", "unknown")
                    if etype:
                        primary_type = etype.split(",")[0]
                        ccre_types[primary_type] = ccre_types.get(primary_type, 0) + 1
        summary["ccre_type_counts"] = ccre_types

    if "fimo_hits" in df.columns:
        n_fimo = (
            df["fimo_hits"].apply(lambda x: len(x) if isinstance(x, list) else 0) > 0
        ).sum()
        summary["n_variants_with_fimo_hits"] = int(n_fimo)

    if "snv_chip_hits" in df.columns:
        n_chip = (
            df["snv_chip_hits"].apply(lambda x: len(x) if isinstance(x, list) else 0) > 0
        ).sum()
        summary["n_variants_with_chip_hits"] = int(n_chip)

    # VAF statistics
    if "tumor_vaf" in df.columns:
        tumor_vaf = df["tumor_vaf"].dropna()
        if len(tumor_vaf) > 0:
            summary["tumor_vaf_median"] = float(tumor_vaf.median())
            summary["tumor_vaf_mean"] = float(tumor_vaf.mean())
    
    # Chromosome distribution
    if "chrom" in df.columns:
        summary["variants_per_chrom"] = df["chrom"].value_counts().to_dict()
    
    return summary


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_coding_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to variants with coding consequences (missense, nonsense, frameshift)."""
    mask = (
        df.get("has_missense", False) |
        df.get("has_nonsense", False) |
        df.get("has_frameshift", False)
    )
    return df[mask].copy()


def get_splice_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to variants with splice consequences."""
    return df[df.get("has_splice_effect", False)].copy()


def get_canonical_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to variants affecting canonical transcripts."""
    return df[df.get("hits_canonical", False)].copy()


def get_regulatory_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to variants overlapping regulatory elements."""
    has_ccre = df.get("cCRE_hits", pd.Series([[]]*len(df))).apply(len) > 0
    has_reg_hit = df.get("regulatory_hits", pd.Series([[]]*len(df))).apply(len) > 0
    has_motif_hit = df.get("motif_hits", pd.Series([[]]*len(df))).apply(len) > 0
    has_fimo = df.get("fimo_hits", pd.Series([[]]*len(df))).apply(len) > 0
    has_chip = df.get("snv_chip_hits", pd.Series([[]]*len(df))).apply(len) > 0

    return df[has_ccre | has_reg_hit | has_motif_hit | has_fimo | has_chip].copy()