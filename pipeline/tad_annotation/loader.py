"""
TAD data loading and multi-biosample annotation orchestration.

Provides functions to:
    - Load TAD domains, boundaries, flanks from processed files
    - Run annotation across multiple TAD biosamples
    - Handle different file formats (parquet, csv, tsv)

This is the main entry point for TAD annotation in the pipeline.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

from .tad_config import (
    TADSourcePaths,
    TAD_BIOSAMPLE_REGISTRY,
    discover_tad_sources,
    get_best_tad_source_for_pam50,
)
from .annotator import (
    annotate_genes_with_tads,
    annotate_ccres_with_tads,
    annotate_lncrnas_with_tads,
    annotate_svs_with_tads,
)
from .mirroring import mirror_all_features_into_domains, mirror_ccres_into_boundaries


def _step4_parallel_feature_tables() -> bool:
    """
    When True, annotate genes / cCREs / lncRNAs for one biosample on up to three threads.

    Each task gets its own shallow copy of the (small) TAD tables so ``annotate_df_with_tads``
    can normalize ``chrom`` columns without races. Disable with ``APM_TAD_STEP4_PARALLEL=0``.
    """
    return os.environ.get("APM_TAD_STEP4_PARALLEL", "1").strip().lower() not in (
        "0",
        "false",
        "no",
    )


# =============================================================================
# FILE LOADING
# =============================================================================

def load_tad_file(path: Path) -> pd.DataFrame:
    """
    Load a TAD file (domains, boundaries, or flanks) with format auto-detection.
    
    Supports: .parquet, .csv, .tsv
    """
    suffix = path.suffix.lower()
    
    if suffix == ".parquet":
        return pd.read_parquet(path)
    elif suffix == ".csv":
        return pd.read_csv(path)
    elif suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    else:
        # Try parquet first, then csv
        try:
            return pd.read_parquet(path)
        except Exception:
            return pd.read_csv(path)


def load_tad_source(paths: TADSourcePaths) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load all TAD data for a single biosample.
    
    Args:
        paths: TADSourcePaths with domains, boundaries, flanks paths
    
    Returns:
        Tuple of (tad_domains, boundaries, flanks) DataFrames
    """
    domains = load_tad_file(paths.domains)
    boundaries = load_tad_file(paths.boundaries)
    flanks = load_tad_file(paths.flanks)
    
    return domains, boundaries, flanks


# =============================================================================
# MULTI-BIOSAMPLE ANNOTATION
# =============================================================================

def annotate_with_all_tad_sources(
    genes_df: pd.DataFrame,
    ccre_df: Optional[pd.DataFrame] = None,
    lncrnas_df: Optional[pd.DataFrame] = None,
    *,
    processed_dir: Path,
    biosamples: Optional[List[str]] = None,
    skip_missing: bool = True,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Annotate features with TAD context from multiple biosamples.
    
    Each biosample's annotation is stored under TAD_domains[biosample_name].
    
    Args:
        genes_df: Genes DataFrame to annotate
        ccre_df: Optional cCREs DataFrame to annotate
        lncrnas_df: Optional lncRNAs DataFrame to annotate
        processed_dir: Path to TADs/processed/ directory
        biosamples: List of biosample names to use (None = all available)
        skip_missing: Skip biosamples with missing files instead of erroring
        verbose: Print progress
    
    Returns:
        Tuple of annotated (genes_df, ccre_df, lncrnas_df)
    
    Example:
        >>> genes_df, ccre_df, lncrnas_df = annotate_with_all_tad_sources(
        ...     genes_df,
        ...     ccre_df,
        ...     lncrnas_df,
        ...     processed_dir=Path("/home/stavz/masters/gdc/TADs/processed"),
        ...     biosamples=["Kim_T47D", "Kim_HCC70", "Rao_HMEC"],
        ... )
    """
    # Discover available sources
    available = discover_tad_sources(processed_dir, required_files=True)
    
    if verbose:
        print(f"Found {len(available)} TAD sources in {processed_dir}")
    
    # Determine which biosamples to process
    if biosamples is None:
        biosamples = list(available.keys())
    else:
        # Filter to requested biosamples that exist
        if skip_missing:
            biosamples = [b for b in biosamples if b in available]
        else:
            missing = [b for b in biosamples if b not in available]
            if missing:
                raise FileNotFoundError(f"TAD sources not found: {missing}")
    
    if verbose:
        print(f"Processing {len(biosamples)} biosamples: {biosamples[:5]}{'...' if len(biosamples) > 5 else ''}")
    
    # Annotate with each biosample
    for i, biosample in enumerate(biosamples):
        if verbose:
            print(f"  [{i+1}/{len(biosamples)}] {biosample}...", end=" ")
        
        paths = available[biosample]
        domains, boundaries, flanks = load_tad_source(paths)

        use_parallel = _step4_parallel_feature_tables() and (
            ccre_df is not None or lncrnas_df is not None
        )
        if use_parallel:
            max_workers = 1 + int(ccre_df is not None) + int(lncrnas_df is not None)
            with ThreadPoolExecutor(max_workers=max_workers) as ex:
                fut_genes = ex.submit(
                    annotate_genes_with_tads,
                    genes_df,
                    domains.copy(),
                    flanks.copy(),
                    boundaries.copy(),
                    biosample,
                )
                fut_ccre = (
                    ex.submit(
                        annotate_ccres_with_tads,
                        ccre_df,
                        domains.copy(),
                        flanks.copy(),
                        boundaries.copy(),
                        biosample,
                    )
                    if ccre_df is not None
                    else None
                )
                fut_lnc = (
                    ex.submit(
                        annotate_lncrnas_with_tads,
                        lncrnas_df,
                        domains.copy(),
                        flanks.copy(),
                        boundaries.copy(),
                        biosample,
                    )
                    if lncrnas_df is not None
                    else None
                )
                genes_df = fut_genes.result()
                if fut_ccre is not None:
                    ccre_df = fut_ccre.result()
                if fut_lnc is not None:
                    lncrnas_df = fut_lnc.result()
        else:
            # Annotate genes (always)
            genes_df = annotate_genes_with_tads(
                genes_df, domains, flanks, boundaries, biosample
            )

            # Annotate cCREs if provided
            if ccre_df is not None:
                ccre_df = annotate_ccres_with_tads(
                    ccre_df, domains, flanks, boundaries, biosample
                )

            # Annotate lncRNAs if provided
            if lncrnas_df is not None:
                lncrnas_df = annotate_lncrnas_with_tads(
                    lncrnas_df, domains, flanks, boundaries, biosample
                )
        
        if verbose:
            print("done")
        del domains, boundaries, flanks
        import gc
        gc.collect()
    
    return genes_df, ccre_df, lncrnas_df


def annotate_svs_with_all_tad_sources(
    sv_df: pd.DataFrame,
    *,
    processed_dir: Path,
    biosamples: Optional[List[str]] = None,
    skip_missing: bool = True,
    ins_mode: str = "point",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Annotate SVs with TAD context from multiple biosamples.
    
    Args:
        sv_df: SVs DataFrame to annotate
        processed_dir: Path to TADs/processed/ directory
        biosamples: List of biosample names (None = all available)
        skip_missing: Skip missing biosamples
        ins_mode: "point" or "span" for insertions
        verbose: Print progress
    
    Returns:
        Annotated sv_df
    """
    available = discover_tad_sources(processed_dir, required_files=True)
    
    if biosamples is None:
        biosamples = list(available.keys())
    elif skip_missing:
        biosamples = [b for b in biosamples if b in available]
    
    if verbose:
        print(f"Annotating SVs with {len(biosamples)} TAD sources")
    
    for biosample in biosamples:
        paths = available[biosample]
        domains, boundaries, flanks = load_tad_source(paths)
        
        sv_df = annotate_svs_with_tads(
            sv_df, domains, flanks, boundaries, biosample,
            ins_mode=ins_mode,
        )
    
    return sv_df


# =============================================================================
# MIRRORING WITH ALL SOURCES
# =============================================================================

def mirror_and_save_all_domains(
    genes_df: pd.DataFrame,
    ccre_df: Optional[pd.DataFrame] = None,
    lncrnas_df: Optional[pd.DataFrame] = None,
    *,
    processed_dir: Path,
    biosamples: Optional[List[str]] = None,
    output_suffix: str = "_with_hits",
    verbose: bool = True,
) -> None:
    """
    Mirror feature annotations back into each TAD domain table and save.
    
    For each biosample:
        1. Load its domain table
        2. Mirror gene/cCRE/lncRNA hits from the annotated feature tables
        3. Save enriched domain table
    
    Args:
        genes_df: Annotated genes DataFrame (with TAD_domains column)
        ccre_df: Optional annotated cCREs DataFrame
        lncrnas_df: Optional annotated lncRNAs DataFrame
        processed_dir: Path to TADs/processed/ directory
        biosamples: List of biosample names (None = all available)
        output_suffix: Suffix for output files (e.g., "_with_hits")
        verbose: Print progress
    """
    available = discover_tad_sources(processed_dir, required_files=True)
    
    if biosamples is None:
        biosamples = list(available.keys())
    
    if verbose:
        print(f"Mirroring to {len(biosamples)} domain tables")
    
    for biosample in biosamples:
        paths = available[biosample]
        domains, boundaries, flanks = load_tad_source(paths)

        # Some TAD sources name ids as `name` rather than `*_id`.
        if "domain_id" not in domains.columns and "name" in domains.columns:
            domains = domains.copy()
            domains["domain_id"] = domains["name"]
        if "boundary_id" not in boundaries.columns and "name" in boundaries.columns:
            boundaries = boundaries.copy()
            boundaries["boundary_id"] = boundaries["name"]
        
        # Mirror all features
        domains = mirror_all_features_into_domains(
            domains,
            genes_df=genes_df,
            lncrnas_df=lncrnas_df,
            ccre_df=ccre_df,
            biosample=biosample,
            mode="primary",
        )
        ctcf_df = ccre_df
        if ccre_df is not None and "type" in ccre_df.columns:
            ctcf_df = ccre_df[ccre_df["type"].astype(str).str.contains("CTCF", na=False)]
        boundaries = mirror_ccres_into_boundaries(
            boundaries,
            ctcf_df,
            biosample=biosample,
            mode="primary",
        )

        def _nullify_empty_dict_cells(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
            """
            Writing a dict-typed column with *only empty dicts* can cause pyarrow to infer a
            struct with zero child fields and raise ArrowNotImplementedError. To keep outputs
            writable, store empty dict cells as nulls.
            """
            df = df.copy()
            for c in cols:
                if c not in df.columns:
                    continue
                df[c] = df[c].apply(lambda x: None if isinstance(x, dict) and len(x) == 0 else x)
            return df
        # Save enriched domain table
        output_path = paths.base_dir / f"domains{output_suffix}.parquet"
        domains = _nullify_empty_dict_cells(domains, ["gene_hits", "lncRNA_hits", "cCRE_hits"])
        domains.to_parquet(output_path)

        output_path = paths.base_dir / f"boundaries_with_ctcf_hits.parquet"
        boundaries = _nullify_empty_dict_cells(boundaries, ["cCRE_hits", "gene_hits", "lncRNA_hits"])
        boundaries.to_parquet(output_path)
        
        if verbose:
            print(f"  Saved {output_path.name} for {biosample}")


# =============================================================================
# PAM50-MATCHED ANNOTATION
# =============================================================================

def get_pam50_matched_sources(
    sample_pam50: str,
    processed_dir: Path,
    fallback_to_normal: bool = True,
) -> List[str]:
    """
    Get TAD sources matching a sample's PAM50 subtype.
    
    Args:
        sample_pam50: PAM50 subtype of the sample
        processed_dir: Path to TADs/processed/
        fallback_to_normal: If no match, fall back to normal tissue sources
    
    Returns:
        List of matching biosample names
    """
    from .tad_config import PAM50_REFERENCE_MAP
    
    available = discover_tad_sources(processed_dir, required_files=True)
    
    # Get candidates for this PAM50
    candidates = PAM50_REFERENCE_MAP.get(sample_pam50, [])
    matches = [c for c in candidates if c in available]
    
    if not matches and fallback_to_normal:
        normal_candidates = PAM50_REFERENCE_MAP.get("Normal", [])
        matches = [c for c in normal_candidates if c in available]
    
    return matches


def annotate_with_pam50_matched(
    genes_df: pd.DataFrame,
    sample_pam50: str,
    *,
    ccre_df: Optional[pd.DataFrame] = None,
    lncrnas_df: Optional[pd.DataFrame] = None,
    processed_dir: Path,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Annotate features using only PAM50-matched TAD sources.
    
    Args:
        genes_df: Genes DataFrame
        sample_pam50: PAM50 subtype to match
        ccre_df: Optional cCREs DataFrame
        lncrnas_df: Optional lncRNAs DataFrame
        processed_dir: Path to TADs/processed/
        verbose: Print progress
    
    Returns:
        Tuple of annotated DataFrames
    """
    matched_sources = get_pam50_matched_sources(sample_pam50, processed_dir)
    
    if verbose:
        print(f"PAM50 '{sample_pam50}' matched to: {matched_sources}")
    
    if not matched_sources:
        if verbose:
            print("  Warning: No matching TAD sources found")
        return genes_df, ccre_df, lncrnas_df
    
    return annotate_with_all_tad_sources(
        genes_df,
        ccre_df,
        lncrnas_df,
        processed_dir=processed_dir,
        biosamples=matched_sources,
        verbose=verbose,
    )
