"""
Main pipeline orchestrator.

Coordinates all modules to build the complete regulatory element table
with multi-evidence gene links.
"""

import os
from pathlib import Path
from typing import List, Optional

import pandas as pd

# Configuration
from .config import (
    PATHS, BIOSAMPLES, THRESHOLDS,
    PRIMARY_GENES, WANTED_SIGNALS, HICHIP_EXPERIMENT,
    OUTPUT_SUBDIRS,
)

# Gene modules
from .genes import (
    load_genes,
    load_lncrnas,
    filter_genes_by_names,
    filter_lncrnas,
    add_promoter_columns,
    create_genes_bed,
    match_lncrnas_to_genes,
    save_lncrna_matching,
    get_mirna_targets,
)

from .genes.lncrna_matching import get_matched_lncrna_names
from .genes.gene_loader import harmonize_multiple_dfs, save_gene_tables

# Regulatory element modules
from .regulatory_elements import (
    load_ccres,
    add_multiple_cell_line_signals,
    match_ccres_to_genes,
    save_all_matching_outputs,
    build_element_focus_table,
    build_gene_summary_table,
)
from .regulatory_elements.element_table import initialize_gene_links_column

# Evidence modules
from .evidence import (
    build_screen_exp_links,
    build_screen_comp_links,
    collapse_screen_to_nested,
    build_abc_links,
    map_abc_to_ccres,
    build_hichip_links,
    integrate_hichip_to_element_table,
    merge_all_evidence,
)
from .evidence.abc_links import collapse_abc_per_link
from .evidence.evidence_merger import (
    attach_gene_links_to_elements,
    print_evidence_summary,
)

# TAD annotation
from .tad_annotation import (
    annotate_with_all_tad_sources,
    mirror_and_save_all_domains,
)


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_full_pipeline(
    working_dir: Optional[Path] = None,
    gene_panel: Optional[List[str]] = None,
    skip_hichip: bool = False,
    skip_abc: bool = False,
    skip_mirna: bool = False,
    skip_tads: bool = False,
    tad_biosamples: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Run the complete regulatory element pipeline.
    
    Args:
        working_dir: Output directory (defaults to config)
        gene_panel: List of genes to process (defaults to PRIMARY_GENES)
        skip_hichip: Skip HiChIP processing (faster for testing)
        skip_abc: Skip ABC processing
        skip_mirna: Skip miRNA processing
        skip_tads: Skip TAD annotation
        tad_biosamples: Specific TAD biosamples to use (None = all available)
    
    Returns:
        Final element focus table with all evidence
    """
    # Setup
    working_dir = Path(working_dir or PATHS.working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    gene_panel = gene_panel or PRIMARY_GENES
    
    print("\n" + "=" * 60)
    print("REGULATORY ELEMENT PIPELINE")
    print("=" * 60)
    print(f"Working directory: {working_dir}")
    print(f"Gene panel size: {len(gene_panel)}")
    
    # =========================================================================
    # STEP 1: Load and filter genes
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 1: Loading genes and lncRNAs")
    print("-" * 40)
    
    genes_all = load_genes(PATHS.gencode_gtf_csv)
    lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
    ccres = load_ccres(PATHS.ccre_csv)
    
    # Harmonize chromosomes
    [genes_all, lncrnas_all, ccres], _ = harmonize_multiple_dfs(
        [genes_all, lncrnas_all, ccres]
    )
    
    # Filter to gene panel
    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(genes, THRESHOLDS.promoter_upstream_bp, THRESHOLDS.promoter_downstream_bp)
    
    lncrnas = filter_lncrnas(lncrnas_all)
    lncrnas = add_promoter_columns(lncrnas, THRESHOLDS.promoter_upstream_bp, THRESHOLDS.promoter_downstream_bp)
    
    print(f"  Filtered genes: {len(genes)}")
    print(f"  Total lncRNAs: {len(lncrnas)}")
    
    # =========================================================================
    # STEP 2: lncRNA matching
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 2: Matching lncRNAs to genes")
    print("-" * 40)
    
    pairs_df, genes_with_lnc, lncrnas_with_genes = match_lncrnas_to_genes(
        genes, lncrnas, THRESHOLDS.lncrna_window_bp
    )
    
    lncrna_output_dir = working_dir / OUTPUT_SUBDIRS["lncrna_matching"]
    save_lncrna_matching(
        lncrna_output_dir, pairs_df, genes_with_lnc, lncrnas_with_genes,
        THRESHOLDS.lncrna_window_bp
    )
    
    matched_lncrna_names = get_matched_lncrna_names(lncrnas_with_genes)
    lncrnas_matched = lncrnas[lncrnas["gene_name"].isin(matched_lncrna_names)]
    print(f"  Matched lncRNAs: {len(matched_lncrna_names)}")
    
    # =========================================================================
    # STEP 3: cCRE-gene distance matching
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 3: Matching cCREs to genes by distance")
    print("-" * 40)
    
    # Add cell-line signals to cCREs
    ccres = add_multiple_cell_line_signals(
        ccres,
        PATHS.cell_lines_dir,
        BIOSAMPLES.ccre_signal_cell_lines,
        WANTED_SIGNALS,
    )
    
    # Match cCREs to protein-coding genes
    pair_df_coding = match_ccres_to_genes(
        genes, ccres,
        window_bp=THRESHOLDS.ccre_window_bp,
        tier_edges=THRESHOLDS.tier_edges_bp,
        tier_labels=THRESHOLDS.tier_labels,
    )
    
    # Build element focus table
    elem_focus = build_element_focus_table(
        ccres, pair_df_coding, THRESHOLDS.tier_labels,
        cell_line_cols=BIOSAMPLES.ccre_signal_cell_lines,
    )
    
    gene_summary = build_gene_summary_table(pair_df_coding, THRESHOLDS.tier_labels)
    
    # Save outputs
    reg_output_dir = working_dir / OUTPUT_SUBDIRS["regulatory_elements"]
    save_all_matching_outputs(
        reg_output_dir, pair_df_coding, elem_focus, gene_summary, "coding"
    )
    
    # Also match lncRNAs (optional output)
    if len(lncrnas_matched) > 0:
        pair_df_lncrna = match_ccres_to_genes(
            lncrnas_matched, ccres,
            window_bp=THRESHOLDS.ccre_window_bp,
        )
        elem_focus_lnc = build_element_focus_table(
            ccres, pair_df_lncrna, THRESHOLDS.tier_labels
        )
        gene_summary_lnc = build_gene_summary_table(pair_df_lncrna, THRESHOLDS.tier_labels)
        save_all_matching_outputs(
            reg_output_dir, pair_df_lncrna, elem_focus_lnc, gene_summary_lnc, "lncRNA"
        )
    
    # =========================================================================
    # STEP 4: TAD annotation
    # =========================================================================
    if not skip_tads:
        print("\n" + "-" * 40)
        print("STEP 4: Annotating features with TAD context")
        print("-" * 40)
        
        genes, ccres, lncrnas_matched = annotate_with_all_tad_sources(
            genes,
            ccres,
            lncrnas_matched if len(lncrnas_matched) > 0 else None,
            processed_dir=PATHS.tads_processed,
            biosamples=tad_biosamples,
            verbose=True,
        )
        
        # Optionally mirror back and save enriched domain tables
        # mirror_and_save_all_domains(
        #     genes, ccres, lncrnas_matched,
        #     processed_dir=PATHS.tads_processed,
        #     biosamples=tad_biosamples,
        # )
    else:
        print("\n  Skipping TAD annotation (skip_tads=True)")
    
    # =========================================================================
    # STEP 5: SCREEN evidence
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 5: Building SCREEN evidence")
    print("-" * 40)
    
    # Experimental links
    links_exp, assays_exp = build_screen_exp_links(
        zip_path=PATHS.screen_exp_zip,
        inner_file=PATHS.screen_exp_inner_file,
        gene_list=gene_panel,
        panel_biosamples=BIOSAMPLES.screen_exp,
        breast_biosamples=BIOSAMPLES.screen_exp,
        intact_pvalue_threshold=THRESHOLDS.intact_hic_pvalue_threshold,
        weak_quantile=THRESHOLDS.screen_weak_quantile,
        strong_quantile=THRESHOLDS.screen_strong_quantile,
        chunksize=THRESHOLDS.screen_chunksize,
    )
    links_exp = collapse_screen_to_nested(
        links_exp, BIOSAMPLES.screen_exp, assays_exp, "screen_exp"
    )
    
    # Computational links
    links_comp, assays_comp = build_screen_comp_links(
        gz_path=PATHS.screen_comp_gz,
        gene_list=gene_panel,
        panel_biosamples=BIOSAMPLES.screen_comp,
        breast_biosamples=BIOSAMPLES.screen_comp,
        weak_quantile=THRESHOLDS.screen_weak_quantile,
        strong_quantile=THRESHOLDS.screen_strong_quantile,
        chunksize=THRESHOLDS.screen_chunksize,
    )
    links_comp = collapse_screen_to_nested(
        links_comp, BIOSAMPLES.screen_comp, assays_comp, "screen_comp"
    )
    
    # =========================================================================
    # STEP 6: ABC evidence
    # =========================================================================
    if not skip_abc:
        print("\n" + "-" * 40)
        print("STEP 6: Building ABC evidence")
        print("-" * 40)
        
        abc_raw = build_abc_links(
            abc_path=PATHS.abc_predictions,
            gene_list=gene_panel,
            celltypes=BIOSAMPLES.abc_celltypes,
            present_threshold=THRESHOLDS.abc_present_threshold,
            strong_threshold=THRESHOLDS.abc_strong_threshold,
            chunksize=THRESHOLDS.abc_chunksize,
        )
        abc_mapped = map_abc_to_ccres(abc_raw, ccres)
        abc_collapsed = collapse_abc_per_link(abc_mapped)
    else:
        print("\n  Skipping ABC evidence (skip_abc=True)")
        abc_collapsed = pd.DataFrame(columns=["ENCODE_id", "gene_name", "ABC_enhancers"])
    
    # =========================================================================
    # STEP 7: Merge evidence into gene_links
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 7: Merging all evidence")
    print("-" * 40)
    
    gene_links_df = merge_all_evidence(
        links_exp=links_exp,
        links_comp=links_comp,
        abc_collapsed=abc_collapsed,
        screen_exp_biosamples=BIOSAMPLES.screen_exp,
        screen_exp_assays=assays_exp,
        screen_comp_biosamples=BIOSAMPLES.screen_comp,
        screen_comp_assays=assays_comp,
    )
    
    elem_focus = attach_gene_links_to_elements(elem_focus, gene_links_df)
    
    # =========================================================================
    # STEP 8: HiChIP evidence
    # =========================================================================
    if not skip_hichip:
        print("\n" + "-" * 40)
        print("STEP 8: Building HiChIP evidence")
        print("-" * 40)
        
        hichip_data = build_hichip_links(
            hichip_dir=PATHS.hichip_dir,
            cell_types=BIOSAMPLES.hichip_panel,
            ccres=ccres,
            genes=genes,
            experiment=HICHIP_EXPERIMENT,
        )
        
        elem_focus = integrate_hichip_to_element_table(elem_focus, hichip_data)
    else:
        print("\n  Skipping HiChIP evidence (skip_hichip=True)")
    
    # =========================================================================
    # STEP 9: miRNA targets
    # =========================================================================
    if not skip_mirna:
        print("\n" + "-" * 40)
        print("STEP 9: Processing miRNA targets")
        print("-" * 40)
        
        mirna_output_dir = working_dir / OUTPUT_SUBDIRS["mirna"]
        get_mirna_targets(
            genes=genes,
            targetscan_path=PATHS.targetscan_predictions,
            top_n=THRESHOLDS.mirna_top_n,
            weight_threshold=THRESHOLDS.mirna_weight_threshold,
            output_dir=mirna_output_dir,
        )
    else:
        print("\n  Skipping miRNA processing (skip_mirna=True)")
    
    # =========================================================================
    # STEP 10: Save final outputs
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 10: Saving final outputs")
    print("-" * 40)
    
    # Save gene tables
    save_gene_tables(genes_all, gene_panel, matched_lncrna_names, working_dir)
    
    # Save BED file
    create_genes_bed(genes, working_dir / "apm_genes.bed")
    
    # Save final element table
    final_output = reg_output_dir / "regulatory_element_focus_with_evidence.csv"
    elem_focus.to_csv(final_output, index=False)
    print(f"  Saved final element table: {final_output}")
    
    # Print summary
    print_evidence_summary(elem_focus)
    
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    
    return elem_focus


# =============================================================================
# MODULAR PIPELINE FUNCTIONS
# =============================================================================

def run_genes_only(
    working_dir: Optional[Path] = None,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Run only the gene/lncRNA processing steps.
    
    Useful for quick iterations on gene filtering logic.
    """
    working_dir = Path(working_dir or PATHS.working_dir)
    gene_panel = gene_panel or PRIMARY_GENES
    
    print("Running gene processing only...")
    
    genes_all = load_genes(PATHS.gencode_gtf_csv)
    lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
    
    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(genes)
    
    lncrnas = filter_lncrnas(lncrnas_all)
    lncrnas = add_promoter_columns(lncrnas)
    
    pairs_df, genes_with_lnc, lncrnas_with_genes = match_lncrnas_to_genes(
        genes, lncrnas, THRESHOLDS.lncrna_window_bp
    )
    
    lncrna_output_dir = working_dir / OUTPUT_SUBDIRS["lncrna_matching"]
    save_lncrna_matching(
        lncrna_output_dir, pairs_df, genes_with_lnc, lncrnas_with_genes,
        THRESHOLDS.lncrna_window_bp
    )
    
    return genes


def run_distance_matching_only(
    genes: pd.DataFrame,
    working_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Run only the cCRE-gene distance matching.
    
    Useful for iterating on distance tier logic.
    """
    working_dir = Path(working_dir or PATHS.working_dir)
    
    print("Running distance matching only...")
    
    ccres = load_ccres(PATHS.ccre_csv)
    
    pair_df = match_ccres_to_genes(
        genes, ccres,
        window_bp=THRESHOLDS.ccre_window_bp,
    )
    
    elem_focus = build_element_focus_table(
        ccres, pair_df, THRESHOLDS.tier_labels
    )
    
    return elem_focus


def run_evidence_only(
    gene_panel: Optional[List[str]] = None,
    include_hichip: bool = True,
    include_abc: bool = True,
) -> pd.DataFrame:
    """
    Run only the evidence collection steps.
    
    Returns merged evidence without integrating into element table.
    """
    gene_panel = gene_panel or PRIMARY_GENES
    
    print("Running evidence collection only...")
    
    # SCREEN experimental
    links_exp, assays_exp = build_screen_exp_links(
        zip_path=PATHS.screen_exp_zip,
        inner_file=PATHS.screen_exp_inner_file,
        gene_list=gene_panel,
        panel_biosamples=BIOSAMPLES.screen_exp,
    )
    links_exp = collapse_screen_to_nested(
        links_exp, BIOSAMPLES.screen_exp, assays_exp, "screen_exp"
    )
    
    # SCREEN computational
    links_comp, assays_comp = build_screen_comp_links(
        gz_path=PATHS.screen_comp_gz,
        gene_list=gene_panel,
        panel_biosamples=BIOSAMPLES.screen_comp,
    )
    links_comp = collapse_screen_to_nested(
        links_comp, BIOSAMPLES.screen_comp, assays_comp, "screen_comp"
    )
    
    # ABC
    if include_abc:
        abc_raw = build_abc_links(
            abc_path=PATHS.abc_predictions,
            gene_list=gene_panel,
            celltypes=BIOSAMPLES.abc_celltypes,
        )
        abc_mapped = map_abc_to_ccres(abc_raw, load_ccres(PATHS.ccre_csv))
        abc_collapsed = collapse_abc_per_link(abc_mapped)
    else:
        abc_collapsed = pd.DataFrame(columns=["ENCODE_id", "gene_name", "ABC_enhancers"])
    
    # Merge
    gene_links_df = merge_all_evidence(
        links_exp=links_exp,
        links_comp=links_comp,
        abc_collapsed=abc_collapsed,
        screen_exp_biosamples=BIOSAMPLES.screen_exp,
        screen_exp_assays=assays_exp,
        screen_comp_biosamples=BIOSAMPLES.screen_comp,
        screen_comp_assays=assays_comp,
    )
    
    return gene_links_df


def run_tad_annotation_only(
    genes: pd.DataFrame,
    ccres: Optional[pd.DataFrame] = None,
    lncrnas: Optional[pd.DataFrame] = None,
    biosamples: Optional[List[str]] = None,
) -> tuple:
    """
    Run only the TAD annotation step.
    
    Useful for testing TAD integration separately.
    
    Args:
        genes: Genes DataFrame
        ccres: Optional cCREs DataFrame
        lncrnas: Optional lncRNAs DataFrame
        biosamples: Specific biosamples to use (None = all)
    
    Returns:
        Tuple of annotated (genes, ccres, lncrnas)
    """
    print("Running TAD annotation only...")
    
    return annotate_with_all_tad_sources(
        genes,
        ccres,
        lncrnas,
        processed_dir=PATHS.tads_processed,
        biosamples=biosamples,
        verbose=True,
    )


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run regulatory element pipeline")
    parser.add_argument("--working-dir", type=str, help="Output directory")
    parser.add_argument("--skip-hichip", action="store_true", help="Skip HiChIP processing")
    parser.add_argument("--skip-abc", action="store_true", help="Skip ABC processing")
    parser.add_argument("--skip-mirna", action="store_true", help="Skip miRNA processing")
    parser.add_argument("--skip-tads", action="store_true", help="Skip TAD annotation")
    parser.add_argument("--genes-only", action="store_true", help="Run only gene processing")
    parser.add_argument("--evidence-only", action="store_true", help="Run only evidence collection")
    
    args = parser.parse_args()
    
    if args.genes_only:
        run_genes_only(working_dir=args.working_dir)
    elif args.evidence_only:
        run_evidence_only(include_hichip=not args.skip_hichip, include_abc=not args.skip_abc)
    else:
        run_full_pipeline(
            working_dir=args.working_dir,
            skip_hichip=args.skip_hichip,
            skip_abc=args.skip_abc,
            skip_mirna=args.skip_mirna,
            skip_tads=args.skip_tads,
        )