"""
Main pipeline orchestrator.

Coordinates all modules to build the complete regulatory element table
with multi-evidence gene links.

Run from the **repository root** (``APM/``), either:

- ``python -m pipeline.main``  (preferred; keeps package context), or
- ``python pipeline/main.py``  (supported: repo root is prepended to ``sys.path`` below).
"""

import gc
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# ``python pipeline/main.py`` puts ``.../pipeline`` on ``sys.path``, not the repo root, so
# ``from .config`` fails (no package parent). Absolute ``pipeline.*`` imports need ``APM/``.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pandas as pd

# Configuration
from pipeline.config import (
    PATHS,
    BIOSAMPLES,
    THRESHOLDS,
    PIPELINE_GENE_PANEL,
    CNV_GENES,
    WANTED_SIGNALS,
    HICHIP_EXPERIMENT,
    OUTPUT_SUBDIRS,
    UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE,
    use_legacy_lncrna_intervals_csv_input,
)

# Gene modules
from pipeline.genes import (
    load_genes,
    load_lncrnas,
    filter_genes_by_names,
    filter_lncrnas,
    add_promoter_columns,
    create_genes_bed,
    match_lncrnas_to_genes,
    save_lncrna_matching,
    get_mirna_targets,
    get_mirtarbase_targets,
)

from pipeline.genes.lncrna_matching import get_matched_lncrna_names
from pipeline.genes.gene_loader import (
    harmonize_multiple_dfs,
    lncrna_gene_intervals_from_annotation,
    save_gene_tables,
)
from pipeline.biosample_names import ccre_signal_output_columns

# Regulatory element modules
from pipeline.regulatory_elements import (
    load_ccres,
    add_multiple_cell_line_signals,
    match_ccres_to_genes,
    save_all_matching_outputs,
    build_element_focus_table,
    build_gene_summary_table,
    save_regulatory_element_focus_evidence_csv,
    save_regulatory_element_focus_evidence_parquet,
)


def _rss_mb() -> float:
    """Process RSS (MB) via ``/proc/self/status``; returns ``-1.0`` if unavailable."""
    try:
        with open("/proc/self/status", "r") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    parts = line.split()
                    return float(parts[1]) / 1024.0
    except Exception:
        pass
    return -1.0


def _log_rss(label: str) -> None:
    rss = _rss_mb()
    if rss >= 0:
        print(f"  [RSS] {label}: {rss:,.0f} MB")
# Evidence modules
from pipeline.evidence import (
    build_screen_exp_links,
    build_screen_comp_links,
    collapse_screen_to_nested,
    build_abc_links,
    map_abc_to_ccres,
    build_hichip_links,
    integrate_hichip_to_element_table,
    merge_all_evidence,
)
from pipeline.evidence.abc_links import collapse_abc_per_link
from pipeline.evidence.evidence_merger import (
    attach_gene_links_to_elements,
    print_evidence_summary,
)

# TAD annotation
from pipeline.tad_annotation import (
    annotate_with_all_tad_sources,
    mirror_and_save_all_domains,
    build_boundaries_enriched_table,
)

from pipeline.atac_peaks import build_atac_peak_table, save_atac_outputs, annotate_df_with_peaks

# ChIP-seq peak integration
from pipeline.CHIP.chip_hits import integrate_chip_hits_to_element_table

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
    skip_atac: bool = False,
    tad_biosamples: Optional[List[str]] = None,
    atac_peaks_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Run the complete regulatory element pipeline.
    
    Args:
        working_dir: Output directory (defaults to config)
        gene_panel: List of genes to process (defaults to ``PIPELINE_GENE_PANEL``)
        skip_hichip: Skip HiChIP processing (faster for testing)
        skip_abc: Skip ABC processing
        skip_mirna: Skip miRNA processing
        skip_tads: Skip TAD annotation
        skip_atac: Skip ATAC peaks processing
        tad_biosamples: Specific TAD biosamples to use (None = all available)
        atac_peaks_path: Path to ATAC peaks file (defaults to config)
    
    Returns:
        Final element focus table with all evidence
    """
    # Setup
    working_dir = Path(working_dir or PATHS.working_dir)
    working_dir.mkdir(parents=True, exist_ok=True)
    gene_panel = gene_panel or PIPELINE_GENE_PANEL
    
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
    
    genes_all = load_genes(PATHS.gencode_gtf_pq)
    if use_legacy_lncrna_intervals_csv_input():
        print("  lncRNA intervals: legacy CSV (APM_LNCRNA_INPUT_LEGACY_CSV=1)")
        lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
    else:
        print("  lncRNA intervals: GENCODE gene rows (gene_type=lncRNA) from gencode_gtf_pq")
        lncrnas_all = lncrna_gene_intervals_from_annotation(genes_all)
    ccres = load_ccres(PATHS.ccre_csv)
    
    # Harmonize chromosomes
    [genes_all, lncrnas_all, ccres], _ = harmonize_multiple_dfs(
        [genes_all, lncrnas_all, ccres]
    )
    
    # Filter to gene panel
    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(genes, THRESHOLDS.promoter_upstream_bp, THRESHOLDS.promoter_downstream_bp)

    cnv_genes = filter_genes_by_names(genes_all, CNV_GENES)
    cnv_genes = add_promoter_columns(cnv_genes, THRESHOLDS.promoter_upstream_bp, THRESHOLDS.promoter_downstream_bp)
    
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
        cell_line_cols=ccre_signal_output_columns(BIOSAMPLES.ccre_signal_cell_lines),
    )
    
    gene_summary = build_gene_summary_table(pair_df_coding, THRESHOLDS.tier_labels)
    
    # Save outputs
    reg_output_dir = working_dir / OUTPUT_SUBDIRS["regulatory_elements"]
    save_all_matching_outputs(
        reg_output_dir, pair_df_coding, elem_focus, gene_summary, "coding"
    )

    # =========================================================================
    # STEP 3b: Add ChIP-seq peak hits to cCREs
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 3b: Annotating cCREs with ChIP-seq peak hits")
    print("-" * 40)
    if PATHS.chip_unified.exists():
        elem_focus = integrate_chip_hits_to_element_table(
            elem_focus=elem_focus,
            ccres_full=ccres,
            chip_unified=PATHS.chip_unified,
        )
        n_with_chip = elem_focus["chip_hits"].apply(lambda x: len(x) > 0).sum()
        print(f"  Added chip_hits to elem_focus: {n_with_chip}/{len(elem_focus)} elements have ChIP overlaps")
    else:
        print(f"  chip_unified not found at {PATHS.chip_unified} — skipping chip_hits")
    
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
    
    # genes.to_csv(working_dir / "genes.csv", index=False)
    # elem_focus.to_csv(working_dir / "elem_focus.csv", index=False)
    # lncrnas_matched.to_csv(working_dir / "lncrnas_matched.csv", index=False)


    # =========================================================================
    # STEP 4: TAD annotation
    # =========================================================================
    if not skip_tads:
        print("\n" + "-" * 40)
        print("STEP 4: Annotating features with TAD context")
        print("-" * 40)
        
        genes, ccres, lncrnas_matched = annotate_with_all_tad_sources(
            genes,
            elem_focus,
            lncrnas_matched if len(lncrnas_matched) > 0 else None,
            processed_dir=PATHS.tads_processed,
            biosamples=tad_biosamples,
            verbose=True,
        )

        # Optional: mirror feature hits back into TAD domain tables + enrich boundaries.
        #
        # This writes enriched per-biosample domain tables under each TAD source folder,
        # plus a single combined boundary table with overlaps (cCRE/ATAC/probes + adjacent-domain genes).
        do_mirror = os.environ.get("APM_TAD_MIRROR", "1").strip().lower() not in ("0", "false", "no")
        do_boundary_enrich = os.environ.get("APM_TAD_BOUNDARY_ENRICH", "1").strip().lower() not in (
            "0",
            "false",
            "no",
        )

        if do_mirror:
            print("\n" + "-" * 40)
            print("STEP 4b: Mirroring feature hits into TAD domains/boundaries")
            print("-" * 40)
            mirror_and_save_all_domains(
                genes_df=genes,
                ccre_df=ccres,
                lncrnas_df=lncrnas_matched if (lncrnas_matched is not None and len(lncrnas_matched) > 0) else None,
                processed_dir=PATHS.tads_processed,
                biosamples=tad_biosamples,
                verbose=True,
            )
        else:
            print("\n  Skipping TAD mirroring (APM_TAD_MIRROR=0)")

        if do_boundary_enrich:
            print("\n" + "-" * 40)
            print("STEP 4c: Building unified enriched boundary table")
            print("-" * 40)
            try:
                bdf = build_boundaries_enriched_table(
                    processed_dir=PATHS.tads_processed,
                    biosamples=tad_biosamples,
                    gene_panel=gene_panel,
                )
                out = Path(PATHS.tads_processed) / "boundaries_enriched.parquet"
                out.parent.mkdir(parents=True, exist_ok=True)
                bdf.to_parquet(out, index=False)
                print(f"  Saved: {out} (rows={len(bdf)}, cols={bdf.shape[1]})")
            except Exception as e:
                print(f"  [WARN] boundary enrichment failed: {e}")
        else:
            print("\n  Skipping boundary enrichment table (APM_TAD_BOUNDARY_ENRICH=0)")
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
        quantiles_out_path=(working_dir / OUTPUT_SUBDIRS["evidence"] / "screen_exp_score_quantiles.csv"),
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
        quantiles_out_path=(working_dir / OUTPUT_SUBDIRS["evidence"] / "screen_comp_score_quantiles.csv"),
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
        print("STEP 9a: Processing miRNA targets")
        print("-" * 40)
        
        mirna_output_dir = working_dir / OUTPUT_SUBDIRS["mirna"]
        get_mirna_targets(
            genes=genes,
            gtf=genes_all,
            targetscan_path=PATHS.targetscan_predictions,
            top_n=THRESHOLDS.mirna_top_n,
            weight_threshold=THRESHOLDS.mirna_weight_threshold,
            output_dir=mirna_output_dir,
        )

        print("\n" + "-" * 40)
        print("STEP 9b: Processing miRTarBase validated targets")
        print("-" * 40)
 
        mirtarbase_output_dir = working_dir / OUTPUT_SUBDIRS["mirna"] / "mirtarbase"
        mirtarbase_tables = get_mirtarbase_targets(
            mirtarbase_csv=PATHS.mirtarbase_csv,
            family_info_tsv=PATHS.mir_family_info,
            gene_panel=gene_panel,
            output_dir=mirtarbase_output_dir,
        )
    else:
        print("\n  Skipping miRNA processing (skip_mirna=True)")
    
    # =========================================================================
    # STEP 10: ATAC peaks processing
    # =========================================================================
    atac_table = None
    if not skip_atac:
        print("\n" + "-" * 40)
        print("STEP 10: Processing ATAC peaks")
        print("-" * 40)

        atac_path = atac_peaks_path or PATHS.atac_peaks_csv
        if atac_path.exists():
            atac_output_dir = working_dir / OUTPUT_SUBDIRS["atac_peaks"]
            atac_output_dir.mkdir(parents=True, exist_ok=True)
            atac_parquet_path = atac_output_dir / "atac_peaks_annotated.parquet"
            atac_table = build_atac_peak_table(
                peaks=atac_path,
                genes=genes,
                ccres=elem_focus,
                lncrnas=lncrnas_matched if len(lncrnas_matched) > 0 else None,
                gene_panel=gene_panel,
                gene_window_bp=THRESHOLDS.atac_gene_window_bp,
                lncrna_window_bp=THRESHOLDS.atac_lncrna_window_bp,
                ccre_max_distance=THRESHOLDS.atac_ccre_max_distance,
                tier_edges=THRESHOLDS.tier_edges_bp,
                tier_labels=THRESHOLDS.tier_labels,
                tad_processed_dir=PATHS.tads_processed if not skip_tads else None,
                tad_biosamples=tad_biosamples,
                tad_parquet_stream=atac_parquet_path if not skip_tads else None,
                verbose=True,
            )

            # Saving ATAC table as CSV requires JSON-encoding large nested columns
            # and can easily double memory. Prefer parquet by default.
            gc.collect()
            save_atac_outputs(
                atac_table,
                atac_output_dir,
                save_parquet=True,
                save_csv=False,
                reuse_main_parquet_if=atac_parquet_path if not skip_tads else None,
            )
            # Peak table nested columns + parquet serialization can peak RAM; drop columns
            # not needed for downstream distance annotation (still on disk in parquet).
            _atac_keep = {
                "peak_id", "chrom", "start", "end", "center", "score", "percentGC", "annotation",
                "n_ccres_total", "n_genes_total", "n_lncrnas_total",
            }
            _drop = [c for c in atac_table.columns if c not in _atac_keep]
            if _drop:
                atac_table.drop(columns=_drop, inplace=True)
            gc.collect()
        else:
            print(f"  ATAC peaks file not found: {atac_path}")
            print("  Skipping ATAC processing")

        if atac_table is not None:
            atac_table = atac_table[
                (atac_table["n_ccres_total"] > 0)
                | (atac_table["n_genes_total"] > 0)
                | (atac_table["n_lncrnas_total"] > 0)
            ]
            atac_reduced = atac_table[
                ["peak_id", "chrom", "start", "end", "center", "score", "percentGC", "annotation"]
            ]

            genes = annotate_df_with_peaks(
                genes, atac_reduced, kind="gene", window_bp=THRESHOLDS.atac_gene_window_bp
            )
            lncrnas = annotate_df_with_peaks(
                lncrnas, atac_reduced, kind="lncrna", window_bp=THRESHOLDS.atac_lncrna_window_bp
            )
            elem_focus = annotate_df_with_peaks(
                elem_focus,
                atac_reduced,
                kind="ccre",
                window_bp=THRESHOLDS.atac_ccre_max_distance,
                copy_df=False,
            )
    else:
        print("\n  Skipping ATAC processing (skip_atac=True)")

    # =========================================================================
    # STEP 11: RNA expression normalization
    # =========================================================================

    # print("\n" + "-" * 40)
    # print("STEP 11: Normalizing RNA expression gene symbols")
    # change_str = "\n".join([f"  {k} -> {v}" for k, v in UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE.items()])
    # print(f"  Normalizing using UCSC RNA-seq gene symbol changes:\n{change_str}")
    # print("-" * 40)
    # rna = pd.read_csv(PATHS.rna_expression_raw, sep="\t", index_col=0)
    # rna = normalize_expression_mat(rna)
    # rna.to_csv(PATHS.rna_expression, sep="\t", index=False)
    
    
    # =========================================================================
    # STEP 12: Save final outputs
    # =========================================================================
    print("\n" + "-" * 40)
    print("STEP 12: Saving final outputs")
    print("-" * 40)
    _log_rss("STEP 12 enter")

    # ---------------------------------------------------------------------
    # FREE memory before heavy writes. ``elem_focus`` is the only object we
    # must keep for STEP 12; large nested companions are already on disk
    # (parquet / CSV outputs from earlier steps) and can be dropped here.
    # Assigning ``locals()[name] = None`` does **not** release Python locals
    # (CPython returns a snapshot), so we issue explicit ``del`` statements.
    # ---------------------------------------------------------------------
    try: del atac_table  # noqa: F821
    except Exception: pass
    try: del atac_reduced  # noqa: F821
    except Exception: pass
    try: del lncrnas  # noqa: F821
    except Exception: pass
    try: del lncrnas_matched  # noqa: F821
    except Exception: pass
    try: del lncrnas_all  # noqa: F821
    except Exception: pass
    try: del hichip_data  # noqa: F821
    except Exception: pass
    try: del ccres  # noqa: F821
    except Exception: pass
    try: del gene_links_df  # noqa: F821
    except Exception: pass
    try: del links_exp  # noqa: F821
    except Exception: pass
    try: del links_comp  # noqa: F821
    except Exception: pass
    try: del abc_collapsed  # noqa: F821
    except Exception: pass
    try: del pair_df_coding  # noqa: F821
    except Exception: pass
    try: del pair_df_lncrna  # noqa: F821
    except Exception: pass
    gc.collect()
    _log_rss("after pre-STEP12 memory release")

    # ---------------------------------------------------------------------
    # Save final ``elem_focus`` as parquet FIRST (the most expensive artifact),
    # so even if later steps fail we keep the main output. Parquet is always
    # written; CSV is opt-in (APM_WRITE_ELEM_FOCUS_CSV=1) — the 70k-row CSV
    # with huge nested ``gene_links`` / ``TAD_domains`` / ``atac_peaks`` is
    # what previously triggered the OOM kill at STEP 12.
    # ---------------------------------------------------------------------
    parquet_final = reg_output_dir / "regulatory_element_focus_with_evidence.parquet"
    save_regulatory_element_focus_evidence_parquet(elem_focus, parquet_final, verbose=True)
    gc.collect()
    _log_rss("after elem_focus parquet write")

    # Save gene tables (needs ``genes_all``; drop it right after).
    save_gene_tables(
        genes,
        cnv_genes,
        gene_panel,
        matched_lncrna_names,
        working_dir,
        genes_all_harmonized=genes_all,
    )
    try:
        del genes_all  # noqa: F821 — released now that Tier 2-4 exports are on disk
    except Exception:
        pass
    gc.collect()
    _log_rss("after save_gene_tables + drop genes_all")

    # Save BED file
    create_genes_bed(genes, working_dir / "apm_genes.bed")

    write_csv_env = os.environ.get("APM_WRITE_ELEM_FOCUS_CSV", "0").strip().lower()
    if write_csv_env in ("1", "true", "yes", "on"):
        csv_final = reg_output_dir / "regulatory_element_focus_with_evidence.csv"
        gc.collect()
        save_regulatory_element_focus_evidence_csv(elem_focus, csv_final, verbose=True)
        _log_rss("after elem_focus CSV write")
    else:
        print(
            "  Skipping regulatory_element_focus_with_evidence.csv "
            "(set APM_WRITE_ELEM_FOCUS_CSV=1 to also write CSV; parquet is the default)."
        )

    # Print summary
    print_evidence_summary(elem_focus)

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    _log_rss("PIPELINE COMPLETE")

    return elem_focus


# =============================================================================
# MODULAR PIPELINE FUNCTIONS
# =============================================================================

def run_genes_only(
    working_dir: Optional[Path] = None,
    gene_panel: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Run only the gene/lncRNA processing steps."""
    working_dir = Path(working_dir or PATHS.working_dir)
    gene_panel = gene_panel or PIPELINE_GENE_PANEL
    
    print("Running gene processing only...")
    
    genes_all = load_genes(PATHS.gencode_gtf_pq)
    if use_legacy_lncrna_intervals_csv_input():
        lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
    else:
        lncrnas_all = lncrna_gene_intervals_from_annotation(genes_all)

    genes = filter_genes_by_names(genes_all, gene_panel)
    genes = add_promoter_columns(genes)

    cnv_genes = filter_genes_by_names(genes_all, CNV_GENES)
    cnv_genes = add_promoter_columns(cnv_genes)
    
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

def run_rna_expression_only(
    *,
    rna_input: Optional[Path] = None,
    rna_output: Optional[Path] = None,
    gene_col: str = PATHS.rna_gene_col,
    verbose: bool = True,
) -> Path:
    """
    Run only STEP 11 (RNA expression gene symbol normalization).

    Reads `PATHS.rna_expression_raw` by default and writes `PATHS.rna_expression`.
    This is intentionally standalone so you can re-run RNA normalization without
    running the full main pipeline.
    """
    import pandas as pd

    in_path = Path(rna_input or PATHS.rna_expression_raw)
    out_path = Path(rna_output or PATHS.rna_expression)
    print(in_path)


    max_alias_lines = os.environ.get("APM_HGNC_ALIAS_MAX_LINES", "").strip()
    max_alias_n = int(max_alias_lines) if max_alias_lines.isdigit() else None

    from pipeline.RNA_exp.normalize_expression_mat import normalize_expression_mat

    if os.environ.get("APM_USE_GENE_SYMBOL_MAPPING", "1").strip().lower() in (
        "0",
        "false",
        "no",
    ):
        symbol_mapping = dict(UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE)
        map_note = "UCSC manual map only (APM_USE_GENE_SYMBOL_MAPPING=0)"
    else:
        from pipeline.genes.symbol_normalization import get_panel_symbol_mapping

        symbol_mapping = get_panel_symbol_mapping(max_alias_n)
        map_note = f"full panel registry (HGNC + manual + UCSC + legacy); {PATHS.hgnc_alias_table}"

    if verbose:
        print("\n" + "-" * 40)
        print("RNA ONLY: STEP 11 (Normalize RNA expression gene symbols)")
        print("-" * 40)
        sample_keys = list(symbol_mapping.items())[:12]
        change_str = "\n".join([f"  {k} -> {v}" for k, v in sample_keys])
        more = "" if len(symbol_mapping) <= 12 else f"\n  ... ({len(symbol_mapping)} total entries)"
        print(f"  Input : {in_path}")
        print(f"  Output: {out_path}")
        print(f"  Map mode: {map_note}")
        print(f"  Using gene_col={gene_col!r}; symbol mapping sample:\n{change_str}{more}")

    rna = pd.read_csv(in_path, sep="\t")
    print(' '.join(rna.columns[:2].tolist()))
    rna = normalize_expression_mat(rna, gene_col=gene_col, mapping=symbol_mapping)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    rna.to_csv(out_path, sep="\t", index=False)

    if verbose:
        print(f"  Wrote: {out_path}")
        print(f"  Rows: {len(rna)}  Cols: {len(rna.columns)}")

    return out_path

def run_cnv_genes_only(
    working_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run only the CNV gene processing."""
    working_dir = Path(working_dir or PATHS.working_dir)
    
    print("Running CNV gene processing only...")
    genes_all = load_genes(PATHS.gencode_gtf_pq)
    cnv_genes = filter_genes_by_names(genes_all, CNV_GENES)
    cnv_genes = add_promoter_columns(cnv_genes)
    
    return cnv_genes


def run_distance_matching_only(
    genes: pd.DataFrame,
    working_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Run only the cCRE-gene distance matching."""
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
    """Run only the evidence collection steps."""
    gene_panel = gene_panel or PIPELINE_GENE_PANEL
    
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
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Run only the TAD annotation step."""
    print("Running TAD annotation only...")
    
    return annotate_with_all_tad_sources(
        genes,
        ccres,
        lncrnas,
        processed_dir=PATHS.tads_processed,
        biosamples=biosamples,
        verbose=True,
    )


def run_atac_only(
    genes: Optional[pd.DataFrame] = None,
    ccres: Optional[pd.DataFrame] = None,
    gene_panel: Optional[List[str]] = None,
    atac_peaks_path: Optional[Path] = None,
    working_dir: Optional[Path] = None,
    include_tads: bool = True,
    tad_biosamples: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Run only the ATAC peaks processing.
    
    Args:
        genes: Pre-loaded genes DataFrame (loads from config if None)
        ccres: Pre-loaded cCREs DataFrame (loads from config if None)
        gene_panel: List of genes to filter to
        atac_peaks_path: Path to ATAC peaks file
        working_dir: Output directory
        include_tads: Whether to include TAD annotations
        tad_biosamples: Specific TAD biosamples to use
    
    Returns:
        ATAC peak table. With ``include_tads=True``, the main nested columns live in
        ``<working_dir>/atac_peaks/atac_peaks_annotated.parquet``; the returned frame is slim.
    """
    print("Running ATAC peaks processing only...")
    
    working_dir = Path(working_dir or PATHS.working_dir)
    gene_panel = gene_panel or PIPELINE_GENE_PANEL
    atac_path = atac_peaks_path or PATHS.atac_peaks_csv
    atac_output_dir = working_dir / OUTPUT_SUBDIRS["atac_peaks"]
    atac_output_dir.mkdir(parents=True, exist_ok=True)
    atac_parquet_path = atac_output_dir / "atac_peaks_annotated.parquet"

    # Load data if not provided
    if genes is None:
        genes_all = load_genes(PATHS.gencode_gtf_pq)
        genes = filter_genes_by_names(genes_all, gene_panel)
        genes = add_promoter_columns(genes, THRESHOLDS.promoter_upstream_bp, THRESHOLDS.promoter_downstream_bp)
    
    if ccres is None:
        ccres = load_ccres(PATHS.ccre_csv)
    
    # Build ATAC table
    atac_table = build_atac_peak_table(
        peaks=atac_path,
        genes=genes,
        ccres=ccres,
        gene_panel=gene_panel,
        gene_window_bp=THRESHOLDS.atac_gene_window_bp,
        ccre_max_distance=THRESHOLDS.atac_ccre_max_distance,
        tier_edges=THRESHOLDS.tier_edges_bp,
        tier_labels=THRESHOLDS.tier_labels,
        tad_processed_dir=PATHS.tads_processed if include_tads else None,
        tad_biosamples=tad_biosamples,
        tad_parquet_stream=atac_parquet_path if include_tads else None,
        verbose=True,
    )

    # Save outputs
    save_atac_outputs(
        atac_table,
        atac_output_dir,
        save_parquet=True,
        save_csv=False,
        reuse_main_parquet_if=atac_parquet_path if include_tads else None,
    )
    print(f"  Saved ATAC peak table: {atac_output_dir}")
    
    return atac_table


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # run_full_pipeline()
    import argparse
    
    parser = argparse.ArgumentParser(description="Run regulatory element pipeline")
    parser.add_argument("--working-dir", type=str, help="Output directory")
    parser.add_argument("--skip-hichip", action="store_true", help="Skip HiChIP processing")
    parser.add_argument("--skip-abc", action="store_true", help="Skip ABC processing")
    parser.add_argument("--skip-mirna", action="store_true", help="Skip miRNA processing")
    parser.add_argument("--skip-tads", action="store_true", help="Skip TAD annotation")
    parser.add_argument("--skip-atac", action="store_true", help="Skip ATAC peaks processing")
    parser.add_argument("--genes-only", action="store_true", help="Run only gene processing")
    parser.add_argument("--evidence-only", action="store_true", help="Run only evidence collection")
    parser.add_argument("--atac-only", action="store_true", help="Run only ATAC processing")
    parser.add_argument("--atac-peaks", type=str, help="Path to ATAC peaks file")
    parser.add_argument("--rna-only", action="store_true", help="Run only RNA expression processing")
    parser.add_argument(
        "--log-file",
        type=str,
        default=None,
        help="Append stdout/stderr to this file; on exit/SIGINT/SIGTERM append dmesg kill/oom hints (Linux)",
    )

    args = parser.parse_args()

    if args.log_file:
        from pipeline.run_logging import install_run_tee

        _lf = Path(args.log_file).expanduser()
        install_run_tee(_lf)
        print(f"[pipeline] tee logging to {_lf}", flush=True)

    if args.genes_only:
        run_genes_only(working_dir=args.working_dir)
    elif args.evidence_only:
        run_evidence_only(include_hichip=not args.skip_hichip, include_abc=not args.skip_abc)
    elif args.atac_only:
        atac_path = Path(args.atac_peaks) if args.atac_peaks else None
        run_atac_only(
            atac_peaks_path=atac_path,
            working_dir=args.working_dir,
            include_tads=not args.skip_tads,
        )
    elif args.rna_only:
        run_rna_expression_only(
            rna_input=PATHS.rna_expression_raw,
            rna_output=PATHS.rna_expression,
            gene_col=PATHS.rna_gene_col,
            verbose=True,
        )
    else:
        atac_path = Path(args.atac_peaks) if args.atac_peaks else None
        run_full_pipeline(
            working_dir=args.working_dir,
            skip_hichip=args.skip_hichip,
            skip_abc=args.skip_abc,
            skip_mirna=args.skip_mirna,
            skip_tads=args.skip_tads,
            skip_atac=args.skip_atac,
            atac_peaks_path=atac_path,
        )
    
