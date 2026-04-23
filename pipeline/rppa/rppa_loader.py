"""
RPPA data loading and expression matrix construction.

Handles:
- Loading per-sample RPPA files
- Constructing sample × target expression matrix
- Loading and validating antibody annotations
- Handling validation status filtering
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Set
from dataclasses import dataclass, field
import warnings

import pandas as pd
import numpy as np


# =============================================================================
# VALIDATION STATUS HANDLING
# =============================================================================

@dataclass
class ValidationStatus:
    """Validation status categories and handling rules."""
    
    VALID: str = "Valid"
    CAUTION: str = "Caution"
    UNDER_EVALUATION: str = "Under Evaluation"
    
    # Default inclusion rules
    include_valid: bool = True
    include_caution: bool = True  # Include but flag
    include_under_evaluation: bool = False  # Exclude by default
    
    @classmethod
    def get_included_statuses(
        cls,
        include_caution: bool = True,
        include_under_evaluation: bool = False,
    ) -> Set[str]:
        """Get set of validation statuses to include."""
        statuses = {cls.VALID}
        if include_caution:
            statuses.add(cls.CAUTION)
        if include_under_evaluation:
            statuses.add(cls.UNDER_EVALUATION)
        return statuses


# =============================================================================
# ANNOTATION LOADING
# =============================================================================

def load_rppa_annotation(
    annotation_path: Path,
    include_caution: bool = True,
    include_under_evaluation: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load and validate RPPA antibody annotation.
    
    Args:
        annotation_path: Path to annotation CSV
        include_caution: Include "Caution" status antibodies
        include_under_evaluation: Include "Under Evaluation" antibodies
    
    Returns:
        Tuple of (filtered_annotation, excluded_annotation)
        
    Expected columns:
        AGID, peptide_target, gene_name, antibody_origin, source,
        catalog_number, validation_status, gene_id
    """
    annot = pd.read_csv(annotation_path, sep="\t")
    
    # Validate required columns
    required_cols = ["AGID", "peptide_target", "gene_name", "validation_status"]
    missing = set(required_cols) - set(annot.columns)
    if missing:
        raise ValueError(f"Annotation missing required columns: {missing}")
    
    # Normalize validation status (handle case variations)
    annot["validation_status"] = annot["validation_status"].str.strip().str.title()
    
    # Get valid statuses
    valid_statuses = ValidationStatus.get_included_statuses(
        include_caution=include_caution,
        include_under_evaluation=include_under_evaluation,
    )
    
    # Split by validation status
    mask_included = annot["validation_status"].isin(valid_statuses)
    
    included = annot[mask_included].copy()
    excluded = annot[~mask_included].copy()

    if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    ):
        from pipeline.genes.symbol_normalization import (
            apply_symbol_mapping_series,
            default_symbol_mapping,
        )

        included["gene_name"] = apply_symbol_mapping_series(
            included["gene_name"], default_symbol_mapping()
        )
    
    # Report
    print(f"RPPA annotation loaded:")
    print(f"  Total antibodies: {len(annot)}")
    for status in annot["validation_status"].unique():
        count = (annot["validation_status"] == status).sum()
        included_str = "✓" if status in valid_statuses else "✗"
        print(f"    {status}: {count} [{included_str}]")
    print(f"  Included: {len(included)}, Excluded: {len(excluded)}")
    
    return included, excluded


def build_target_mappings(annotation: pd.DataFrame) -> Dict[str, Dict]:
    """
    Build comprehensive mappings from annotation.
    
    Returns:
        Dict with keys:
        - target_to_gene: peptide_target → gene_name
        - target_to_agid: peptide_target → AGID
        - gene_to_targets: gene_name → [peptide_targets]
        - agid_to_target: AGID → peptide_target
        - caution_targets: set of targets with "Caution" status
    """
    mappings = {
        "target_to_gene": {},
        "target_to_agid": {},
        "gene_to_targets": {},
        "agid_to_target": {},
        "caution_targets": set(),
    }
    
    for _, row in annotation.iterrows():
        target = row["peptide_target"]
        gene = row["gene_name"]
        agid = row["AGID"]
        status = row.get("validation_status", "Valid")
        
        mappings["target_to_gene"][target] = gene
        mappings["target_to_agid"][target] = agid
        mappings["agid_to_target"][agid] = target
        
        if gene not in mappings["gene_to_targets"]:
            mappings["gene_to_targets"][gene] = []
        mappings["gene_to_targets"][gene].append(target)
        
        if status == "Caution":
            mappings["caution_targets"].add(target)
    
    return mappings


# =============================================================================
# PER-SAMPLE FILE LOADING
# =============================================================================

def sample_id_from_rppa_filename(file_path: Path) -> Optional[str]:
    """
    Derive a specimen-level id from a TCGA-style RPPA filename.

    GDC per-sample files often name like ``TCGA-...-01A-..._RPPA_data.tsv`` while
    the ``lab_id`` column holds many internal replicate ids per row — it must
    not be used as ``sample_id`` for pivoting.
    """
    from ..sample_ids import normalize_tcga_id

    stem = Path(file_path).stem
    for suf in ("_RPPA_data", "_rppa_data", "_RPPA", "_rppa", "_normalized"):
        if len(stem) > len(suf) and stem.lower().endswith(suf.lower()):
            stem = stem[: -len(suf)]
            break
    tid = normalize_tcga_id(stem.strip())
    if tid.sample_vial:
        return str(tid.sample_vial)
    if tid.sample:
        return str(tid.sample)
    if tid.aliquot:
        return str(tid.aliquot)
    if tid.participant:
        return str(tid.participant)
    return None


def load_single_sample_rppa(
    file_path: Path,
    valid_targets: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """
    Load RPPA data from a single sample file.
    
    Expected format:
        AGID, lab_id, catalog_number, set_id, peptide_target, protein_expression
    
    Args:
        file_path: Path to sample RPPA file
        valid_targets: If provided, filter to only these targets
    
    Returns:
        DataFrame with columns: [peptide_target, protein_expression, lab_id]
    """
    df = pd.read_csv(file_path, sep="\t")
    
    # Validate required columns
    required = ["peptide_target", "protein_expression", "lab_id"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Sample file {file_path} missing columns: {missing}")
    
    # Filter to valid targets if provided
    if valid_targets is not None:
        df = df[df["peptide_target"].isin(valid_targets)].copy()

    return df[["peptide_target", "protein_expression", "lab_id"]]


def load_all_sample_rppa(
    sample_dir: Path,
    annotation: pd.DataFrame,
    metadata: Optional[pd.DataFrame] = None,
    file_pattern: str = "*.tsv",
    sample_id_source: str = "metadata",
    # also: "filename" (TCGA stem / sample_vial), "lab_id" (legacy; first row only)
    metadata_file_col: str = "File Name",
    metadata_sample_col: str = "Sample ID",
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Load RPPA data from all sample files in a directory.
    
    Args:
        sample_dir: Directory containing per-sample RPPA files
        annotation: Filtered annotation DataFrame
        file_pattern: Glob pattern for sample files
        sample_id_source: How to derive sample ID
            - "metadata": Use manifest columns (``metadata_file_col`` → ``metadata_sample_col``)
            - "filename": Parse TCGA id from file path stem (recommended without manifest)
            - "lab_id": First ``lab_id`` in file (legacy; wrong for many TCGA exports)
    
    Returns:
        Tuple of:
        - Long-format DataFrame: [sample_id, peptide_target, protein_expression]
        - Dict mapping filename → sample_id
    """
    sample_dir = Path(sample_dir)
    valid_targets = set(annotation["peptide_target"])
    
    all_data = []
    file_to_sample = {}
    
    files = list(sample_dir.glob(file_pattern))
    if not files:
        raise ValueError(f"No files matching '{file_pattern}' in {sample_dir}")
    
    use_manifest = (
        sample_id_source == "metadata"
        and metadata is not None
        and not metadata.empty
        and metadata_file_col in metadata.columns
        and metadata_sample_col in metadata.columns
    )
    if sample_id_source == "metadata" and not use_manifest:
        print(
            "  Note: metadata manifest missing or incomplete; "
            "using lab_id from each RPPA file as sample_id."
        )

    print(f"Loading RPPA data from {len(files)} sample files...")
    
    for i, fpath in enumerate(files):
        try:
            df = load_single_sample_rppa(fpath, valid_targets)
            
            # Determine sample ID
            if use_manifest:
                rows = metadata.loc[metadata[metadata_file_col] == fpath.name]
                if rows.empty:
                    warnings.warn(
                        f"No metadata row where {metadata_file_col!r} == {fpath.name!r}; "
                        f"falling back to filename-derived sample_id."
                    )
                    sid = sample_id_from_rppa_filename(fpath)
                    sample_id = sid if sid else str(df["lab_id"].iloc[0])
                else:
                    sample_id = str(rows[metadata_sample_col].iloc[0])
            elif sample_id_source == "lab_id":
                sample_id = str(df["lab_id"].iloc[0])
            else:
                # "filename" or default when manifest unused
                sid = sample_id_from_rppa_filename(fpath)
                sample_id = sid if sid else str(df["lab_id"].iloc[0])
            
            df["sample_id"] = sample_id
            file_to_sample[fpath.name] = sample_id
            
            all_data.append(df[["sample_id", "peptide_target", "protein_expression"]])
            
            if (i + 1) % 50 == 0:
                print(f"  Loaded {i + 1}/{len(files)} files...")
                
        except Exception as e:
            warnings.warn(f"Error loading {fpath}: {e}")
            continue
    
    if not all_data:
        raise ValueError("No sample data could be loaded")
    
    combined = pd.concat(all_data, ignore_index=True)
    print(f"  Total: {len(combined)} measurements from {combined['sample_id'].nunique()} samples")
    
    return combined, file_to_sample


# =============================================================================
# EXPRESSION MATRIX CONSTRUCTION
# =============================================================================

def build_expression_matrix(
    long_data: pd.DataFrame,
    annotation: pd.DataFrame,
    handle_duplicates: str = "mean",  # "mean", "first", "error"
    min_samples_per_target: int = 10,
    min_targets_per_sample: int = 50,
) -> Tuple[pd.DataFrame, Dict[str, any]]:
    """
    Construct sample × target expression matrix from long-format data.
    
    Args:
        long_data: Long-format DataFrame [sample_id, peptide_target, protein_expression]
        annotation: Antibody annotation for validation
        handle_duplicates: How to handle duplicate measurements
        min_samples_per_target: Drop targets with fewer samples
        min_targets_per_sample: Drop samples with fewer targets
    
    Returns:
        Tuple of:
        - Expression matrix (samples × targets)
        - QC stats dict
    """
    # Check for duplicates
    dup_check = long_data.groupby(["sample_id", "peptide_target"]).size()
    has_dups = (dup_check > 1).any()
    
    if has_dups:
        n_dups = (dup_check > 1).sum()
        if handle_duplicates == "error":
            raise ValueError(f"Found {n_dups} duplicate sample-target pairs")
        elif handle_duplicates == "mean":
            print(f"  Averaging {n_dups} duplicate measurements")
            long_data = (
                long_data
                .groupby(["sample_id", "peptide_target"], as_index=False)
                ["protein_expression"]
                .mean()
            )
        else:  # first
            print(f"  Taking first of {n_dups} duplicate measurements")
            long_data = long_data.drop_duplicates(
                subset=["sample_id", "peptide_target"],
                keep="first"
            )
    
    # Pivot to wide format
    expr_matrix = long_data.pivot(
        index="sample_id",
        columns="peptide_target",
        values="protein_expression"
    )
    
    # QC: target coverage
    samples_per_target = expr_matrix.notna().sum(axis=0)
    targets_to_drop = samples_per_target[samples_per_target < min_samples_per_target].index
    
    if len(targets_to_drop) > 0:
        print(f"  Dropping {len(targets_to_drop)} targets with <{min_samples_per_target} samples")
        expr_matrix = expr_matrix.drop(columns=targets_to_drop)
    
    # QC: sample coverage
    targets_per_sample = expr_matrix.notna().sum(axis=1)
    samples_to_drop = targets_per_sample[targets_per_sample < min_targets_per_sample].index
    
    if len(samples_to_drop) > 0:
        print(f"  Dropping {len(samples_to_drop)} samples with <{min_targets_per_sample} targets")
        expr_matrix = expr_matrix.drop(index=samples_to_drop)
    
    # Compute QC stats
    denom = int(expr_matrix.size)
    missing_rate = float(expr_matrix.isna().sum().sum() / denom) if denom > 0 else 0.0
    qc_stats = {
        "n_samples": len(expr_matrix),
        "n_targets": len(expr_matrix.columns),
        "missing_rate": missing_rate,
        "targets_dropped": list(targets_to_drop),
        "samples_dropped": list(samples_to_drop),
        "samples_per_target": samples_per_target.to_dict(),
        "targets_per_sample": targets_per_sample.to_dict(),
    }
    
    print(f"  Final matrix: {qc_stats['n_samples']} samples × {qc_stats['n_targets']} targets")
    print(f"  Missing rate: {qc_stats['missing_rate']:.2%}")
    
    return expr_matrix, qc_stats


# =============================================================================
# PHOSPHO-FORM IDENTIFICATION
# =============================================================================

def identify_phospho_pairs(
    annotation: pd.DataFrame,
) -> Dict[str, Dict[str, Union[str, List[str]]]]:
    """
    Identify total/phospho antibody pairs for ratio computation.
    
    Phospho-specific antibodies typically contain patterns like:
    - _pS123 (phospho-serine at position 123)
    - _pT456 (phospho-threonine)
    - _pY789 (phospho-tyrosine)
    - _pT123Y456 (dual phosphorylation)
    
    Returns:
        Dict[gene_name] = {
            "total": target_name or None,
            "phospho": [phospho_target_names],
            "phospho_sites": {target: site_string}
        }
    """
    import re
    
    # Pattern to detect phospho-specific antibodies
    phospho_pattern = re.compile(r"_p([STY]\d+)", re.IGNORECASE)
    
    pairs = {}
    
    for gene, group in annotation.groupby("gene_name"):
        targets = group["peptide_target"].tolist()
        
        total_candidates = []
        phospho_targets = []
        phospho_sites = {}
        
        for t in targets:
            match = phospho_pattern.search(t)
            if match:
                phospho_targets.append(t)
                # Extract all phospho sites
                sites = phospho_pattern.findall(t)
                phospho_sites[t] = "_".join(sites)
            else:
                total_candidates.append(t)
        
        # Select best total (prefer exact gene name match, shorter names)
        total = None
        if total_candidates:
            # Sort by length (shorter = more likely to be total)
            total_candidates.sort(key=len)
            total = total_candidates[0]
        
        if total or phospho_targets:
            pairs[gene] = {
                "total": total,
                "phospho": phospho_targets,
                "phospho_sites": phospho_sites,
            }
    
    # Report
    n_with_pairs = sum(
        1 for p in pairs.values() 
        if p["total"] is not None and len(p["phospho"]) > 0
    )
    print(f"Identified {n_with_pairs} genes with total + phospho pairs")
    
    return pairs


# =============================================================================
# SAMPLE METADATA INTEGRATION
# =============================================================================

def integrate_sample_metadata(
    expr_matrix: pd.DataFrame,
    metadata_path: Path,
    sample_id_col: str = "sample_id",
    merge_on: str = "index",  # or column name
) -> pd.DataFrame:
    """
    Integrate sample metadata with expression matrix.
    
    Args:
        expr_matrix: Expression matrix (samples × targets)
        metadata_path: Path to sample metadata TSV
        sample_id_col: Column name for sample ID in metadata
        merge_on: How to match samples
            - "index": Match metadata sample_id_col to expr_matrix index
            - column name: Match on specific column
    
    Returns:
        Metadata DataFrame with same index as expr_matrix
    """
    metadata = pd.read_csv(metadata_path, sep="\t")
    
    if sample_id_col not in metadata.columns:
        raise ValueError(f"Metadata missing column: {sample_id_col}")
    
    # Set sample ID as index for matching
    metadata = metadata.set_index(sample_id_col)
    
    # Align with expression matrix
    common_samples = expr_matrix.index.intersection(metadata.index)
    
    if len(common_samples) == 0:
        warnings.warn("No matching samples between expression matrix and metadata")
        return pd.DataFrame(index=expr_matrix.index)
    
    if len(common_samples) < len(expr_matrix):
        missing = set(expr_matrix.index) - set(common_samples)
        warnings.warn(
            f"{len(missing)} samples in expression matrix not in metadata"
        )
    
    return metadata.loc[metadata.index.intersection(expr_matrix.index)]


# =============================================================================
# MAIN LOADING FUNCTION
# =============================================================================

def load_rppa_dataset(
    sample_dir: Path,
    annotation_path: Path,
    metadata_path: Optional[Path] = None,
    include_caution: bool = True,
    include_under_evaluation: bool = False,
    file_pattern: str = "*.tsv",
    min_samples_per_target: int = 10,
    min_targets_per_sample: int = 50,
    metadata_file_col: str = "File Name",
    metadata_sample_id_col: str = "Sample ID",
) -> Dict[str, any]:
    """
    Complete RPPA dataset loading pipeline.
    
    Args:
        sample_dir: Directory with per-sample RPPA files
        annotation_path: Path to antibody annotation
        metadata_path: Optional path to sample metadata (TSV with file manifest columns)
        metadata_file_col: Column in metadata matching each per-sample RPPA filename
        metadata_sample_id_col: Column in metadata with the TCGA-style sample id
        include_caution: Include "Caution" validation status
        include_under_evaluation: Include "Under Evaluation" status
        file_pattern: Glob pattern for sample files
        min_samples_per_target: QC threshold
        min_targets_per_sample: QC threshold
    
    Returns:
        Dict with:
        - expression_matrix: samples × targets DataFrame
        - annotation: filtered annotation DataFrame
        - excluded_annotation: excluded antibodies
        - mappings: target/gene mapping dicts
        - phospho_pairs: total/phospho pair information
        - qc_stats: quality control statistics
        - metadata: sample metadata (if provided)
        - caution_flag: DataFrame flagging caution targets
    """
    print("=" * 60)
    print("RPPA DATASET LOADING")
    print("=" * 60)
    
    # 1. Load annotation
    print("\n[1/5] Loading annotation...")
    annotation, excluded = load_rppa_annotation(
        annotation_path,
        include_caution=include_caution,
        include_under_evaluation=include_under_evaluation,
    )
    
    mappings = build_target_mappings(annotation)
    phospho_pairs = identify_phospho_pairs(annotation)

    manifest_df: Optional[pd.DataFrame] = None
    meta_path = Path(metadata_path) if metadata_path is not None else None
    if meta_path is not None and meta_path.is_file():
        print("\n[2/5] Loading sample metadata manifest...")
        manifest_df = pd.read_csv(meta_path, sep="\t", low_memory=False)
        if metadata_file_col not in manifest_df.columns or metadata_sample_id_col not in manifest_df.columns:
            warnings.warn(
                f"Metadata at {meta_path} missing {metadata_file_col!r} or "
                f"{metadata_sample_id_col!r}; will use lab_id from each RPPA file."
            )
            manifest_df = None
        else:
            print(f"  Manifest rows: {len(manifest_df)}")
    elif meta_path is not None:
        warnings.warn(f"Metadata path does not exist or is not a file: {meta_path}")

    print("\n[3/5] Loading sample files...")
    long_data, file_map = load_all_sample_rppa(
        sample_dir,
        annotation,
        manifest_df,
        file_pattern,
        sample_id_source="metadata" if manifest_df is not None else "filename",
        metadata_file_col=metadata_file_col,
        metadata_sample_col=metadata_sample_id_col,
    )

    print("\n[4/5] Building expression matrix...")
    expr_matrix, qc_stats = build_expression_matrix(
        long_data,
        annotation,
        min_samples_per_target=min_samples_per_target,
        min_targets_per_sample=min_targets_per_sample,
    )

    metadata_aligned: Optional[pd.DataFrame] = None
    if meta_path is not None and meta_path.is_file():
        try:
            print("\n[5/5] Aligning metadata to expression matrix...")
            metadata_aligned = integrate_sample_metadata(
                expr_matrix, meta_path, sample_id_col=metadata_sample_id_col
            )
            print(f"  Metadata rows aligned to matrix: {len(metadata_aligned)}")
        except Exception as e:
            warnings.warn(f"Could not align metadata to expression matrix: {e}")
            metadata_aligned = None

    # Create caution flag matrix
    caution_flags = pd.DataFrame(
        False,
        index=expr_matrix.index,
        columns=expr_matrix.columns
    )
    for target in mappings["caution_targets"]:
        if target in caution_flags.columns:
            caution_flags[target] = True
    
    print("\n" + "=" * 60)
    print("LOADING COMPLETE")
    print("=" * 60)
    
    return {
        "expression_matrix": expr_matrix,
        "annotation": annotation,
        "excluded_annotation": excluded,
        "mappings": mappings,
        "phospho_pairs": phospho_pairs,
        "qc_stats": qc_stats,
        "metadata": metadata_aligned,
        "caution_flags": caution_flags,
        "file_map": file_map,
    }
