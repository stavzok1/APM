"""
Main RPPA pipeline orchestrator.

Coordinates all RPPA processing steps and integrates with
the main regulatory element pipeline.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings

import pandas as pd
import numpy as np

from .rppa_loader import (
    load_rppa_dataset,
    build_target_mappings,
    identify_phospho_pairs,
)
from .rppa_panels import (
    RPPAMarkerPanels,
    get_default_panels,
    compute_all_panel_scores,
    compute_all_activation_ratios,
    summarize_panel_coverage,
)
from .rppa_analysis import (
    detect_signaling_blocks,
    classify_ddr_ifn_quadrant,
    classify_immune_visibility_state,
    compute_protein_rna_discordance,
    identify_post_transcriptional_regulation,
    compute_immune_visibility_score,
    compute_antigen_presentation_capacity,
    correlate_with_outcomes,
)
from ..sample_ids import add_tcga_id_columns_from_index_inplace
from .rppa_config import (
    RPPAPathConfig,
    RPPAThresholdConfig,
    RPPA_TARGET_TO_GENE_OVERRIDES,
)


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_rppa_pipeline(
    sample_dir: Path = RPPAPathConfig.sample_data_dir,
    annotation_path: Path = RPPAPathConfig.annotation_file,
    output_dir: Path = RPPAPathConfig.output_dir,
    metadata_path: Path = RPPAPathConfig.sample_metadata_tumor,
    rna_expr: Optional[pd.DataFrame] = None,
    include_caution: bool = True,
    include_under_evaluation: bool = False,
    scoring_method: str = "zscore",
    min_samples_per_target: int = 10,
    min_targets_per_sample: int = 50,
    save_outputs: bool = True,
) -> Dict[str, Any]:
    """
    Run the complete RPPA analysis pipeline.
    
    Args:
        sample_dir: Directory containing per-sample RPPA files
        annotation_path: Path to antibody annotation CSV
        output_dir: Directory for output files
        rna_expr: Optional RNA expression matrix for discordance analysis
        metadata_path: Optional sample metadata CSV
        include_caution: Include "Caution" validation status
        include_under_evaluation: Include "Under Evaluation" status
        scoring_method: Panel scoring method ("zscore" or "pca")
        min_samples_per_target: QC threshold
        min_targets_per_sample: QC threshold
        save_outputs: Whether to save outputs to files
    
    Returns:
        Dict containing all processed data and results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 60)
    print("RPPA ANALYSIS PIPELINE")
    print("=" * 60)
    
    # =========================================================================
    # STEP 1: Load Data
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 1: Loading RPPA data")
    print("-" * 40)
    
    data = load_rppa_dataset(
        sample_dir=sample_dir,
        annotation_path=annotation_path,
        metadata_path=metadata_path,
        include_caution=include_caution,
        include_under_evaluation=include_under_evaluation,
        min_samples_per_target=min_samples_per_target,
        min_targets_per_sample=min_targets_per_sample,
    )
    
    expr_matrix = data["expression_matrix"]
    annotation = data["annotation"]
    mappings = data["mappings"]
    phospho_pairs = data["phospho_pairs"]
    metadata = data["metadata"]
    
    # =========================================================================
    # STEP 2: Panel Coverage Check
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 2: Checking panel coverage")
    print("-" * 40)
    
    panels = get_default_panels()
    coverage = summarize_panel_coverage(expr_matrix, panels)
    print(coverage.to_string(index=False))
    
    if save_outputs:
        coverage.to_csv(output_dir / "panel_coverage.csv", index=False)
    
    # =========================================================================
    # STEP 3: Compute Panel Scores
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 3: Computing panel scores")
    print("-" * 40)
    
    panel_scores = compute_all_panel_scores(
        expr_matrix, panels, method=scoring_method
    )
    
    print(f"  Computed {len(panel_scores.columns)} panel scores")
    print(f"  Non-null rates:")
    for col in panel_scores.columns:
        rate = panel_scores[col].notna().mean()
        print(f"    {col}: {rate:.1%}")
    
    if save_outputs:
        panel_scores.to_csv(output_dir / "panel_scores.csv")
    
    # =========================================================================
    # STEP 4: Compute Activation Ratios
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 4: Computing activation ratios")
    print("-" * 40)
    
    activation_ratios = compute_all_activation_ratios(
        expr_matrix, phospho_pairs
    )
    
    print(f"  Computed {len(activation_ratios.columns)} activation ratios")
    
    if save_outputs:
        activation_ratios.to_csv(output_dir / "activation_ratios.csv")
    
    # =========================================================================
    # STEP 5: Detect Signaling Blocks
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 5: Detecting signaling blocks")
    print("-" * 40)
    
    blocks = detect_signaling_blocks(expr_matrix)
    
    print("  Block frequencies:")
    for col in blocks.columns:
        freq = blocks[col].mean()
        print(f"    {col}: {freq:.1%}")
    
    if save_outputs:
        blocks.to_csv(output_dir / "signaling_blocks.csv")
    
    # =========================================================================
    # STEP 6: DDR-IFN Stratification
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 6: Stratifying by DDR-IFN state")
    print("-" * 40)
    
    ddr_score = panel_scores.get("DDR_activation", pd.Series(dtype=float))
    ifn_score = panel_scores.get("IFN_activated", pd.Series(dtype=float))
    
    if not ddr_score.isna().all() and not ifn_score.isna().all():
        quadrants = classify_ddr_ifn_quadrant(ddr_score, ifn_score)
        print("  DDR-IFN quadrant distribution:")
        print(quadrants.value_counts().to_string())
    else:
        quadrants = pd.Series("unknown", index=expr_matrix.index)
        print("  Could not compute DDR-IFN quadrants (missing scores)")
    
    # =========================================================================
    # STEP 7: Immune Visibility Classification
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 7: Classifying immune visibility state")
    print("-" * 40)
    
    visibility = classify_immune_visibility_state(expr_matrix, panel_scores)
    
    print("  Visibility state distribution:")
    print(visibility["visibility_state"].value_counts().to_string())
    
    if save_outputs:
        visibility.to_csv(output_dir / "visibility_classification.csv")
    
    # =========================================================================
    # STEP 8: Compute Visibility Score
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 8: Computing visibility score")
    print("-" * 40)
    
    visibility_score = compute_immune_visibility_score(panel_scores, expr_matrix)
    
    print(f"  Score range: {visibility_score.min():.1f} - {visibility_score.max():.1f}")
    print(f"  Mean: {visibility_score.mean():.1f}, Median: {visibility_score.median():.1f}")
    
    # =========================================================================
    # STEP 9: Protein-RNA Discordance (if RNA provided)
    # =========================================================================
    discordance = None
    post_transcriptional = None
    
    if rna_expr is not None:
        print("\n" + "-" * 40)
        print("STEP 9: Computing protein-RNA discordance")
        print("-" * 40)
        
        discordance = compute_protein_rna_discordance(
            expr_matrix,
            rna_expr,
            mappings["target_to_gene"],
        )
        
        post_transcriptional = identify_post_transcriptional_regulation(
            discordance
        )
        
        # Report genes with systematic regulation
        regulated_genes = [
            g for g, info in post_transcriptional.items()
            if info["dominant_pattern"] in ["stabilization", "degradation"]
        ]
        print(f"  Found {len(regulated_genes)} genes with systematic post-transcriptional regulation")
        
        if save_outputs:
            discordance.to_csv(output_dir / "protein_rna_discordance.csv")
            
            # Save post-transcriptional summary
            pt_summary = pd.DataFrame([
                {"gene": g, **info}
                for g, info in post_transcriptional.items()
            ])
            pt_summary.to_csv(output_dir / "post_transcriptional_regulation.csv", index=False)
    else:
        print("\n  Skipping protein-RNA discordance (no RNA data provided)")
    
    # =========================================================================
    # STEP 10: Antigen Presentation Capacity
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 10: Estimating antigen presentation capacity")
    print("-" * 40)
    
    apm_capacity = compute_antigen_presentation_capacity(
        expr_matrix, rna_expr
    )
    
    if "estimated_apm_capacity" in apm_capacity.columns:
        cap = apm_capacity["estimated_apm_capacity"]
        print(f"  APM capacity range: {cap.min():.2f} - {cap.max():.2f}")
    
    if save_outputs:
        apm_capacity.to_csv(output_dir / "apm_capacity.csv")
    
    # =========================================================================
    # STEP 11: Save Combined Results
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 11: Saving combined results")
    print("-" * 40)
    
    # Combine key results into single DataFrame
    combined = pd.DataFrame(index=expr_matrix.index)
    combined["visibility_score"] = visibility_score
    combined["visibility_state"] = visibility["visibility_state"]
    combined["ddr_ifn_quadrant"] = quadrants
    
    # Add panel scores
    for col in panel_scores.columns:
        combined[f"panel_{col}"] = panel_scores[col]
    
    # Add key activation ratios
    for col in activation_ratios.columns:
        combined[f"ratio_{col}"] = activation_ratios[col]
    
    # Add signaling blocks
    for col in blocks.columns:
        combined[f"block_{col}"] = blocks[col]
    
    # Add APM capacity
    if "estimated_apm_capacity" in apm_capacity.columns:
        combined["estimated_apm_capacity"] = apm_capacity["estimated_apm_capacity"]
    
    if save_outputs:
        # Ensure outputs carry normalized TCGA join keys (sample_vial/sample/participant)
        add_tcga_id_columns_from_index_inplace(combined)
        combined.to_csv(output_dir / "rppa_analysis_combined.csv")
        
        # Also save as parquet for faster loading with nested data
        combined.to_parquet(output_dir / "rppa_analysis_combined.parquet")
    
    print(f"\n  Saved combined results: {len(combined)} samples × {len(combined.columns)} features")
    
    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("RPPA PIPELINE COMPLETE")
    print("=" * 60)
    print(f"\nOutputs saved to: {output_dir}")
    print(f"  - panel_coverage.csv")
    print(f"  - panel_scores.csv")
    print(f"  - activation_ratios.csv")
    print(f"  - signaling_blocks.csv")
    print(f"  - visibility_classification.csv")
    print(f"  - apm_capacity.csv")
    print(f"  - rppa_analysis_combined.csv/parquet")
    if discordance is not None:
        print(f"  - protein_rna_discordance.csv")
        print(f"  - post_transcriptional_regulation.csv")
    
    return {
        "expression_matrix": expr_matrix,
        "annotation": annotation,
        "mappings": mappings,
        "phospho_pairs": phospho_pairs,
        "panel_scores": panel_scores,
        "activation_ratios": activation_ratios,
        "signaling_blocks": blocks,
        "ddr_ifn_quadrants": quadrants,
        "visibility_classification": visibility,
        "visibility_score": visibility_score,
        "apm_capacity": apm_capacity,
        "protein_rna_discordance": discordance,
        "post_transcriptional_regulation": post_transcriptional,
        "combined": combined,
        "qc_stats": data["qc_stats"],
        "metadata": data.get("metadata"),
        "caution_flags": data["caution_flags"],
    }


# =============================================================================
# INTEGRATION WITH MAIN PIPELINE
# =============================================================================

def integrate_rppa_with_samples(
    rppa_results: Dict[str, Any],
    sample_ids: List[str],
) -> pd.DataFrame:
    """
    Extract RPPA features for a specific set of samples.
    
    Args:
        rppa_results: Output from run_rppa_pipeline
        sample_ids: Sample IDs to extract
    
    Returns:
        DataFrame with RPPA features for specified samples
    """
    combined = rppa_results["combined"]
    
    # Find matching samples
    common = set(sample_ids).intersection(combined.index)
    
    if not common:
        warnings.warn("No matching samples between RPPA and requested list")
        return pd.DataFrame(index=sample_ids)
    
    if len(common) < len(sample_ids):
        missing = set(sample_ids) - common
        warnings.warn(f"{len(missing)} samples not found in RPPA data")
    
    return combined.loc[list(common)].copy()


def add_rppa_to_gene_table(
    gene_table: pd.DataFrame,
    rppa_results: Dict[str, Any],
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """
    Add RPPA-derived features to a per-gene table.
    
    Adds columns for:
    - Visibility score
    - Panel scores
    - DDR-IFN quadrant
    - Signaling block flags
    
    Args:
        gene_table: DataFrame with sample_col column
        rppa_results: Output from run_rppa_pipeline
        sample_col: Column containing sample identifiers
    
    Returns:
        Gene table with RPPA columns added
    """
    if sample_col not in gene_table.columns:
        raise ValueError(f"Gene table missing column: {sample_col}")
    
    combined = rppa_results["combined"]
    
    # Create mapping for samples in gene_table
    result = gene_table.copy()
    
    # Select key columns to add
    rppa_cols = [
        "visibility_score",
        "visibility_state",
        "ddr_ifn_quadrant",
        "panel_IFN_activated",
        "panel_DDR_activation",
        "panel_cGAS_STING",
        "panel_checkpoint",
        "panel_PI3K_AKT",
        "panel_STAT3_suppressive",
        "panel_cytolytic",
        "block_ddr_sting_block",
        "block_sting_irf1_block",
        "block_stat3_override",
        "block_checkpoint_escape",
        "block_pi3k_activated",
        "estimated_apm_capacity",
    ]
    
    available_cols = [c for c in rppa_cols if c in combined.columns]
    
    for col in available_cols:
        result[f"rppa_{col}"] = result[sample_col].map(
            combined[col].to_dict()
        )
    
    return result


# =============================================================================
# QUICK ANALYSIS FUNCTIONS
# =============================================================================

def quick_rppa_summary(
    sample_dir: Path,
    annotation_path: Path,
) -> pd.DataFrame:
    """
    Quick summary of RPPA data without full pipeline.
    
    Returns basic statistics about the dataset.
    """
    data = load_rppa_dataset(
        sample_dir=sample_dir,
        annotation_path=annotation_path,
    )
    
    expr = data["expression_matrix"]
    
    summary = {
        "n_samples": len(expr),
        "n_targets": len(expr.columns),
        "missing_rate": f"{expr.isna().mean().mean():.1%}",
        "n_valid_antibodies": len(data["annotation"]),
        "n_excluded_antibodies": len(data["excluded_annotation"]),
        "n_phospho_pairs": sum(
            1 for p in data["phospho_pairs"].values()
            if p["total"] and p["phospho"]
        ),
    }
    
    return pd.DataFrame([summary])


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run RPPA analysis pipeline")
    parser.add_argument("--sample-dir", type=str, required=True,
                        help="Directory with per-sample RPPA files")
    parser.add_argument("--annotation", type=str, required=True,
                        help="Path to antibody annotation CSV")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Output directory")
    parser.add_argument("--rna-expr", type=str, default=None,
                        help="Optional RNA expression matrix")
    parser.add_argument("--metadata", type=str, default=None,
                        help="Optional sample metadata")
    parser.add_argument("--include-caution", action="store_true", default=True,
                        help="Include Caution status antibodies")
    parser.add_argument("--exclude-caution", action="store_true",
                        help="Exclude Caution status antibodies")
    parser.add_argument("--scoring-method", choices=["zscore", "pca"],
                        default="zscore", help="Panel scoring method")
    
    args = parser.parse_args()
    
    # Load RNA if provided
    rna = None
    if args.rna_expr:
        rna = pd.read_csv(args.rna_expr, index_col=0)
    
    run_rppa_pipeline(
        sample_dir=Path(args.sample_dir),
        annotation_path=Path(args.annotation),
        output_dir=Path(args.output_dir),
        rna_expr=rna,
        metadata_path=Path(args.metadata) if args.metadata else None,
        include_caution=not args.exclude_caution,
        scoring_method=args.scoring_method,
    )
