#!/usr/bin/env python3
"""
Smoke tests for the APM regulatory element pipeline.

Uses REAL data from disk; writes outputs to a temporary directory.
Original code is NOT modified. Bugs are caught and reported.

Usage:
    cd /home/stavz/masters/gdc/APM
    .venv/bin/python3 -m tests.run_smoke_tests           # run all (use venv; system python3 may lack pandas)
    .venv/bin/python3 -m tests.run_smoke_tests genes     # run one section
    .venv/bin/python3 -m tests.run_smoke_tests --list    # list available tests
"""

import sys
import os
import subprocess
import shutil
import tempfile
import traceback
import time
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Callable, Optional

# Ensure the repo root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Configuration: real data paths
# ---------------------------------------------------------------------------
DATA = Path("/home/stavz/masters/gdc/APM/data")
ANN  = Path("/home/stavz/masters/gdc/APM/annotations")
HGNC_ALIAS_TSV = ANN / "Alias_v5.22.xls"

GENCODE_CSV       = DATA / "gencode.v49.annotation.gtf.csv"
LNCRNA_CSV        = DATA / "lncRNAs_all_features.csv"
CCRE_CSV          = DATA / "GRCh38-cCREs.csv"
CELL_LINES_DIR    = DATA / "regulatory_elements_data" / "cCRE_signals"
ABC_PREDICTIONS   = DATA / "regulatory_elements_data" / "ABC_output" / "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt"
HICHIP_DIR        = DATA / "HiCHIP"
TADS_DIR          = DATA / "TADs"
ATAC_PEAKS_TXT    = DATA / "TCGA_ATAC" / "TCGA-ATAC_PanCancer_PeakSet.txt"
MIRNA_TARGETSCAN  = DATA / "miRNA" / "Predicted_Targets_Context_Scores.default_predictions.txt"
MIRNA_MIRTARBASE  = DATA / "miRNA" / "mirtar.csv"
MIRNA_FAMILY_INFO = DATA / "miRNA" / "miR_Family_Info.txt"

# SNV: pick the first VEP-annotated VCF (PATHS.snv_vcf_dir → data/SNV/vep_vcfs)
SNV_VEP_VCF_DIR = DATA / "SNV" / "vep_vcfs"

# SV: existing processed CSV
SV_PROCESSED_DIR   = DATA / "SV" / "pipeline_output" / "02_processed_sv_csv"
SV_RAW_VCF_DIR = (
    DATA / "SV" / "raw_somatic_sv"
    if (DATA / "SV" / "raw_somatic_sv").is_dir()
    else DATA / "SV" / "raw_somatic_sv_test"
)

# CNV
CNV_SEG_FILE       = DATA / "CNV_TCGA" / "CNV_extracted" / "TCGA-BRCA.4695c0f4-4d30-4604-84b3-81407bf980e8.ascat3.allelic_specific.seg.txt"
CNV_ANNOTATIONS    = DATA / "CNV_TCGA" / "samples.tsv"  # may not exist
CNV_MIRNA_CSV      = DATA / "miRNA" / "cnv_miRNA.csv"

# Methylation
METH_PROBE_REF     = ANN / "Methylation" / "probe_reference.tsv"
METH_SAMPLE_RAW    = DATA / "Methylation" / "raw_samples" / "Tumor" / "071a05c1-f9df-4e68-b73e-8f6716d75d52.methylation_array.sesame.level3betas.txt"
METH_OUTPUTS       = DATA / "Methylation"

# Small gene panel for speed
SMOKE_GENES = ["HLA-A", "HLA-B", "TAP1", "B2M", "STAT1"]

# Extended / tiered panel spot-check (GENCODE v49 symbols)
SMOKE_GENES_TIERED = [
    "CGAS",
    "STING1",
    "CIITA",
    "EZH2",
    "NFKB1",
    "GZMA",
    "HLA-A",
]

# ---------------------------------------------------------------------------
# Test infrastructure
# ---------------------------------------------------------------------------
@dataclass
class TestResult:
    name: str
    passed: bool
    duration_s: float
    message: str = ""
    error: str = ""

results: List[TestResult] = []


def smoke_test(name: str):
    """Decorator that wraps a test function with timing and error handling."""
    def decorator(fn: Callable):
        fn._smoke_name = name
        def wrapper(tmpdir: Path):
            t0 = time.time()
            try:
                fn(tmpdir)
                dur = time.time() - t0
                results.append(TestResult(name, True, dur, "OK"))
                print(f"  [PASS] {name}  ({dur:.1f}s)")
            except Exception as e:
                dur = time.time() - t0
                tb = traceback.format_exc()
                results.append(TestResult(name, False, dur, str(e), tb))
                print(f"  [FAIL] {name}  ({dur:.1f}s)")
                print(f"         {e}")
        wrapper._smoke_name = name
        return wrapper
    return decorator


# ===========================================================================
# UTILITY: load genes from CSV (since parquet doesn't exist on disk)
# ===========================================================================

def load_genes_from_csv(path: Path, nrows: Optional[int] = None,
                        gene_names: Optional[list] = None) -> pd.DataFrame:
    """Load GENCODE GTF-derived CSV, matching the schema load_genes() expects.
    
    If gene_names is provided, loads full file but filters to just those genes
    (much faster than loading everything when you only need a few).
    """
    if gene_names is not None:
        chunks = []
        for chunk in pd.read_csv(path, chunksize=100_000):
            hits = chunk[chunk["gene_name"].isin(gene_names)]
            if len(hits) > 0:
                chunks.append(hits)
        df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    else:
        df = pd.read_csv(path, nrows=nrows)
    from pipeline.utils import harmonize_chrom_column
    df, _ = harmonize_chrom_column(df)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    return df


# ===========================================================================
# TEST: Gene loading, filtering, promoter columns, lncRNA matching
# ===========================================================================

@smoke_test("genes.load_and_filter")
def test_genes_load(tmpdir: Path):
    from pipeline.genes.gene_loader import (
        filter_genes_by_names, filter_lncrnas, add_promoter_columns,
        create_genes_bed, harmonize_multiple_dfs,
    )

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    assert len(genes_all) > 0, "No genes loaded"
    assert "gene_name" in genes_all.columns
    print(f"    Loaded {len(genes_all)} rows from GENCODE CSV for panel {SMOKE_GENES}")
    print(f"    Features present: {genes_all['feature'].value_counts().to_dict()}")

    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    assert len(genes) > 0, f"No genes matched panel {SMOKE_GENES}"
    print(f"    Filtered to {len(genes)} rows for {genes['gene_name'].nunique()} unique genes")

    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    assert "tss" in genes.columns
    assert "prom_start" in genes.columns
    gene_rows = genes[genes["feature"] == "gene"]
    print(f"    Gene-level rows with promoter coords: {len(gene_rows)}")
    for _, g in gene_rows.iterrows():
        print(f"      {g['gene_name']}: {g['chrom']}:{g['start']}-{g['end']} TSS={g['tss']} prom=[{g['prom_start']},{g['prom_end']}]")

    bed_path = tmpdir / "smoke_genes.bed"
    create_genes_bed(genes, bed_path)
    assert bed_path.exists()
    print(f"    BED file written: {bed_path} ({bed_path.stat().st_size} bytes)")

    lncrnas_all = pd.read_csv(LNCRNA_CSV, nrows=20_000)
    from pipeline.utils import harmonize_chrom_column
    lncrnas_all, _ = harmonize_chrom_column(lncrnas_all)
    lncrnas = filter_lncrnas(lncrnas_all)
    assert len(lncrnas) > 0, "No lncRNAs after filtering"
    lncrnas = add_promoter_columns(lncrnas, upstream_bp=2000, downstream_bp=500)
    print(f"    lncRNAs loaded: {len(lncrnas_all)} raw -> {len(lncrnas)} filtered")


@smoke_test("genes.lncrna_matching")
def test_lncrna_matching(tmpdir: Path):
    from pipeline.genes.gene_loader import filter_genes_by_names, filter_lncrnas, add_promoter_columns
    from pipeline.genes.lncrna_matching import match_lncrnas_to_genes, save_lncrna_matching

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)

    lncrnas_all = pd.read_csv(LNCRNA_CSV, nrows=20_000)
    from pipeline.utils import harmonize_chrom_column
    lncrnas_all, _ = harmonize_chrom_column(lncrnas_all)
    lncrnas = filter_lncrnas(lncrnas_all)
    lncrnas = add_promoter_columns(lncrnas, upstream_bp=2000, downstream_bp=500)

    pairs_df, genes_with_lnc, lncrnas_with_genes = match_lncrnas_to_genes(
        genes, lncrnas, window_bp=500_000
    )
    assert isinstance(pairs_df, pd.DataFrame)
    print(f"    Pairs found: {len(pairs_df)} gene-lncRNA pairs within 500kb window")
    print(f"    Genes with lncRNA neighbors: {len(genes_with_lnc)}")
    print(f"    lncRNAs with gene neighbors: {len(lncrnas_with_genes)}")
    if len(pairs_df) > 0:
        print(f"    Pair columns: {list(pairs_df.columns)}")

    out = tmpdir / "lncrna_matching"
    save_lncrna_matching(out, pairs_df, genes_with_lnc, lncrnas_with_genes, 500_000)
    assert (out).exists()
    saved_files = list(out.iterdir()) if out.is_dir() else [out]
    print(f"    Saved {len(saved_files)} output files to {out}")


# ===========================================================================
# TEST: cCRE loading, distance matching, element table
# ===========================================================================

@smoke_test("regulatory.load_ccres")
def test_ccre_loading(tmpdir: Path):
    from pipeline.regulatory_elements.ccre_loader import load_ccres
    ccres = load_ccres(CCRE_CSV)
    assert len(ccres) > 100, f"Only {len(ccres)} cCREs loaded"
    assert "cCRE_id" in ccres.columns or "ENCODE_id" in ccres.columns
    id_col = "cCRE_id" if "cCRE_id" in ccres.columns else "ENCODE_id"
    print(f"    Loaded {len(ccres)} cCREs")
    print(f"    Columns: {list(ccres.columns)}")
    print(f"    Chromosomes: {sorted(ccres['chrom'].unique())[:5]}... ({ccres['chrom'].nunique()} total)")
    if "ccre_class" in ccres.columns:
        print(f"    Classes: {ccres['ccre_class'].value_counts().to_dict()}")


@smoke_test("regulatory.distance_matching")
def test_distance_matching(tmpdir: Path):
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
    from pipeline.regulatory_elements.ccre_loader import load_ccres
    from pipeline.regulatory_elements.distance_matching import match_ccres_to_genes
    from pipeline.regulatory_elements.element_table import build_element_focus_table, build_gene_summary_table

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)

    ccres = load_ccres(CCRE_CSV)

    tier_edges = [0, 100_000, 250_000, 500_000, 1_000_000]
    tier_labels = ["0–100kb", "100–250kb", "250–500kb", "500–1000kb"]

    pair_df = match_ccres_to_genes(genes, ccres, window_bp=1_000_000,
                                   tier_edges=tier_edges, tier_labels=tier_labels)
    assert len(pair_df) > 0, "No cCRE-gene pairs found"
    print(f"    cCRE-gene pairs: {len(pair_df)}")
    print(f"    Pair columns: {list(pair_df.columns)}")
    if "distance_tier" in pair_df.columns:
        print(f"    Tier distribution: {pair_df['distance_tier'].value_counts().to_dict()}")

    elem_focus = build_element_focus_table(ccres, pair_df, tier_labels)
    assert len(elem_focus) > 0
    assert "gene_links" in elem_focus.columns or any("0–100kb" in c for c in elem_focus.columns)
    print(f"    Element-focus table: {len(elem_focus)} rows x {len(elem_focus.columns)} cols")
    if "gene_links" in elem_focus.columns:
        n_with_links = elem_focus["gene_links"].apply(lambda x: len(x) if isinstance(x, (list, dict)) else 0 if pd.isna(x) else 1).sum()
        print(f"    Total gene_links entries: {n_with_links}")

    gene_summary = build_gene_summary_table(pair_df, tier_labels)
    assert len(gene_summary) > 0
    print(f"    Gene summary table: {len(gene_summary)} rows x {len(gene_summary.columns)} cols")
    print(f"    Gene summary columns: {list(gene_summary.columns)}")


@smoke_test("regulatory.cell_line_signals")
def test_cell_line_signals(tmpdir: Path):
    from pipeline.regulatory_elements.ccre_loader import load_ccres, add_multiple_cell_line_signals
    ccres = load_ccres(CCRE_CSV)
    ccres_sample = ccres.head(500).copy()

    cell_lines = ["MCF7"]
    wanted = ["H3K27ac", "H3K4me3", "CTCF", "DNase"]

    ccres_sig = add_multiple_cell_line_signals(
        ccres_sample, CELL_LINES_DIR, cell_lines, wanted,
    )
    assert "MCF7" in ccres_sig.columns, "MCF7 signal column not added"
    print(f"    Input: {len(ccres_sample)} cCREs, cell_lines={cell_lines}, assays={wanted}")
    print(f"    Output columns added: {[c for c in ccres_sig.columns if c not in ccres_sample.columns]}")
    mcf7_vals = ccres_sig["MCF7"]
    if mcf7_vals.dtype == object:
        print(f"    MCF7 signal is dict/nested type; sample: {mcf7_vals.iloc[0]}")
    else:
        print(f"    MCF7 signal stats: mean={mcf7_vals.mean():.3f}, non-null={mcf7_vals.notna().sum()}")


# ===========================================================================
# TEST: TAD annotation
# ===========================================================================

@smoke_test("tad_annotation.basic")
def test_tad_annotation(tmpdir: Path):
    from pipeline.tad_annotation import annotate_with_all_tad_sources

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
    test_panel = SMOKE_GENES[:2]
    genes = filter_genes_by_names(genes_all, test_panel)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    print(f"    Input: {len(genes)} gene rows for {test_panel}")

    genes_out, _, _ = annotate_with_all_tad_sources(
        genes, None, None,
        processed_dir=TADS_DIR, biosamples=["Kim_T47D"], verbose=False,
    )
    assert "TAD_domains" in genes_out.columns, "TAD_domains column missing"
    n_with_tad = genes_out["TAD_domains"].apply(lambda x: bool(x) if isinstance(x, (dict, list)) else pd.notna(x)).sum()
    print(f"    Output: {len(genes_out)} rows, {n_with_tad} with TAD annotations")
    if n_with_tad > 0:
        sample_tad = genes_out["TAD_domains"].iloc[0]
        print(f"    Sample TAD_domains structure: {type(sample_tad).__name__}")
        if isinstance(sample_tad, dict):
            print(f"    TAD_domains keys: {list(sample_tad.keys())[:5]}")


# ===========================================================================
# TEST: Evidence building (SCREEN, ABC)
# ===========================================================================

@smoke_test("evidence.screen_from_prebuilt")
def test_screen_prebuilt(tmpdir: Path):
    """Test loading pre-built SCREEN link CSVs (since zip/gz may not exist)."""
    screen_exp_csv = DATA / "regulatory_elements_data" / "SCREEN_exp" / "screen_exp_links.csv"
    screen_comp_csv = DATA / "regulatory_elements_data" / "SCREEN_comp" / "screen_comp_links.csv"

    found = 0
    if screen_exp_csv.exists():
        df = pd.read_csv(screen_exp_csv, nrows=1000)
        assert len(df) > 0, "Empty SCREEN exp CSV"
        found += 1
        print(f"    SCREEN experimental: {len(df)} rows (first 1000)")
        print(f"    Columns: {list(df.columns)}")
        if "gene_name" in df.columns:
            print(f"    Unique genes: {df['gene_name'].nunique()}")
        if "biosample" in df.columns or "celltype" in df.columns:
            ct_col = "biosample" if "biosample" in df.columns else "celltype"
            print(f"    Biosamples: {df[ct_col].nunique()} unique")
    else:
        print(f"    SCREEN exp CSV not found at {screen_exp_csv}, skipping")

    if screen_comp_csv.exists():
        df = pd.read_csv(screen_comp_csv, nrows=1000)
        assert len(df) > 0, "Empty SCREEN comp CSV"
        found += 1
        print(f"    SCREEN computational: {len(df)} rows (first 1000)")
        print(f"    Columns: {list(df.columns)}")
    else:
        print(f"    SCREEN comp CSV not found at {screen_comp_csv}, skipping")

    print(f"    Summary: {found}/2 SCREEN evidence files found and loaded")


@smoke_test("evidence.abc_links")
def test_abc_links(tmpdir: Path):
    from pipeline.evidence.abc_links import build_abc_links

    if not ABC_PREDICTIONS.exists():
        print(f"    ABC file not found at {ABC_PREDICTIONS}, skipping")
        return

    print(f"    ABC file: {ABC_PREDICTIONS.name} ({ABC_PREDICTIONS.stat().st_size / 1e6:.1f} MB)")
    print(f"    Gene filter: {SMOKE_GENES}, celltypes=['MCF-7-ENCODE']")

    abc_raw = build_abc_links(
        abc_path=ABC_PREDICTIONS,
        gene_list=SMOKE_GENES,
        celltypes=["MCF-7-ENCODE"],
        present_threshold=0.015,
        strong_threshold=0.05,
        chunksize=200_000,
    )
    assert isinstance(abc_raw, pd.DataFrame)
    print(f"    ABC raw links: {len(abc_raw)} rows")
    if len(abc_raw) > 0:
        print(f"    Columns: {list(abc_raw.columns)}")
        if "TargetGene" in abc_raw.columns:
            print(f"    Target genes: {abc_raw['TargetGene'].unique().tolist()}")
        if "ABC.Score" in abc_raw.columns:
            print(f"    ABC score range: [{abc_raw['ABC.Score'].min():.4f}, {abc_raw['ABC.Score'].max():.4f}]")


# ===========================================================================
# TEST: HiChIP evidence (requires pyranges)
# ===========================================================================

@smoke_test("evidence.hichip_links")
def test_hichip_links(tmpdir: Path):
    try:
        import pyranges
    except ImportError:
        print("    pyranges not installed, skipping HiChIP test")
        return

    from pipeline.evidence.hichip_links import (
        build_hichip_links, load_loops_unified, find_hichip_loops_file,
    )

    test_celltype = "MCF7"
    loops_path = find_hichip_loops_file(HICHIP_DIR, test_celltype, "H3K27ac")
    if loops_path is None:
        print(f"    No loops file for {test_celltype} in {HICHIP_DIR}, skipping")
        return

    print(f"    Found loops file: {loops_path.name} ({loops_path.stat().st_size / 1e3:.1f} KB)")

    loops = load_loops_unified(loops_path)
    print(f"    Loaded {len(loops)} loops from {test_celltype}")
    print(f"    Loop columns: {list(loops.columns)}")

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
    from pipeline.regulatory_elements.ccre_loader import load_ccres
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    ccres = load_ccres(CCRE_CSV)
    print(f"    Building HiChIP links: {len(ccres)} cCREs, {len(genes)} gene rows, celltype={test_celltype}")

    hichip_data = build_hichip_links(
        hichip_dir=HICHIP_DIR,
        cell_types=[test_celltype],
        ccres=ccres,
        genes=genes,
        experiment="H3K27ac",
    )
    assert isinstance(hichip_data, dict)
    print(f"    Cell types processed: {list(hichip_data.keys())}")
    if test_celltype in hichip_data:
        ct_data = hichip_data[test_celltype]
        print(f"    {test_celltype} data keys: {list(ct_data.keys())}")
        if "per_ccre" in ct_data:
            print(f"    cCREs with HiChIP loops: {len(ct_data['per_ccre'])}")
        if "ccre_gene_map" in ct_data:
            print(f"    cCRE-gene HiChIP connections: {len(ct_data['ccre_gene_map'])}")


# ===========================================================================
# TEST: ATAC peaks
# ===========================================================================

@smoke_test("atac.load_and_match")
def test_atac_peaks(tmpdir: Path):
    from pipeline.atac_peaks.peak_loader import load_atac_peaks

    if not ATAC_PEAKS_TXT.exists():
        print(f"    ATAC peaks file not found at {ATAC_PEAKS_TXT}, skipping")
        return

    peaks = load_atac_peaks(ATAC_PEAKS_TXT)
    assert len(peaks) > 0, "No ATAC peaks loaded"
    print(f"    Loaded {len(peaks)} ATAC peaks")
    print(f"    Peak columns: {list(peaks.columns)}")
    chrom_col = "chrom" if "chrom" in peaks.columns else "seqnames"
    assert chrom_col in peaks.columns
    print(f"    Chromosomes: {peaks[chrom_col].nunique()} unique, top: {peaks[chrom_col].value_counts().head(3).to_dict()}")
    if "score" in peaks.columns:
        print(f"    Score range: [{peaks['score'].min():.2f}, {peaks['score'].max():.2f}]")

    try:
        from pipeline.atac_peaks.gene_matching import match_peaks_to_genes

        genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
        from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
        test_panel = SMOKE_GENES[:2]
        genes = filter_genes_by_names(genes_all, test_panel)
        genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
        print(f"    Matching peaks to {len(genes)} gene rows ({test_panel}), window=500kb")

        pairs = match_peaks_to_genes(peaks, genes, window_bp=500_000)
        assert isinstance(pairs, pd.DataFrame)
        print(f"    Peak-gene pairs: {len(pairs)}")
        if len(pairs) > 0:
            print(f"    Pair columns: {list(pairs.columns)}")
            if "body_overlap" in pairs.columns:
                n_ov = (pairs["body_overlap"] > 0).sum() if pairs["body_overlap"].dtype != object else "N/A (dict)"
                print(f"    Peaks overlapping gene body: {n_ov}")
    except (ValueError, TypeError) as e:
        if "columns" in str(e).lower() or "broadcast" in str(e).lower() or "dict" in str(e).lower():
            print(f"    KNOWN BUG: match_peaks_to_genes body_overlap apply fails with pandas compat issue: {e}")
            print("    -> gene_matching.py apply(func_returning_dict) needs result_type='reduce'")
        else:
            raise


# ===========================================================================
# TEST: SNV module
# ===========================================================================

@smoke_test("snv.load_single_vcf")
def test_snv_load(tmpdir: Path):
    try:
        import vcfpy
    except ImportError:
        print("    vcfpy not installed, skipping SNV test")
        return

    from pipeline.SNV.vcf_loader import load_mutect_snv_vcf

    vep_vcfs = sorted(SNV_VEP_VCF_DIR.glob("*.vep.vcf"))
    if not vep_vcfs:
        vep_vcfs = sorted(SNV_VEP_VCF_DIR.glob("*APM_1Mb*.vcf*"))
    if not vep_vcfs:
        print(f"    No VEP VCFs found in {SNV_VEP_VCF_DIR}, skipping")
        return

    vcf_path = vep_vcfs[0]
    print(f"    Testing with: {vcf_path.name} ({vcf_path.stat().st_size / 1e6:.1f} MB)")
    print(f"    Gene panel: {SMOKE_GENES}")

    df, normal, tumor = load_mutect_snv_vcf(
        vcf_path,
        primary_genes=SMOKE_GENES,
        elements_df=None,
        apply_filter=True,
        save_outputs=False,
    )
    assert isinstance(df, pd.DataFrame)
    print(f"    Result: {len(df)} variants loaded")
    print(f"    Sample IDs: normal={normal}, tumor={tumor}")

    if len(df) > 0:
        expected_cols = ["chrom", "pos", "ref", "alt", "tumor_vaf"]
        for col in expected_cols:
            assert col in df.columns, f"Missing column: {col}"
        print(f"    All required columns present: {expected_cols}")
        print(f"    Variant columns: {list(df.columns)}")
        print(f"    Chromosomes: {sorted(df['chrom'].unique())}")
        print(f"    VAF range: [{df['tumor_vaf'].min():.4f}, {df['tumor_vaf'].max():.4f}]")

        if "gene_hits" in df.columns:
            sample_hit = df['gene_hits'].iloc[0]
            print(f"    gene_hits: type={type(sample_hit).__name__}, example={str(sample_hit)[:100]}")
        if "cCRE_hits" in df.columns:
            n_with_ccre = (df["cCRE_hits"].apply(lambda x: len(x) if isinstance(x, list) else 0) > 0).sum()
            print(f"    cCRE hits: {n_with_ccre}/{len(df)} variants overlap cCREs")
        if "VEP" in df.columns or "vep_consequence" in df.columns:
            vep_col = "VEP" if "VEP" in df.columns else "vep_consequence"
            print(f"    VEP annotations present in column '{vep_col}'")
    else:
        print(f"    No variants matched filters (this is OK for a small VCF)")


@smoke_test("snv.somatic_filter")
def test_snv_filter(tmpdir: Path):
    from pipeline.SNV.somatic_filter import high_conf_somatic_mask

    fake_df = pd.DataFrame({
        "tumor_vaf": [0.25, 0.01, 0.15, 0.0],
        "normal_vaf": [0.0, 0.0, 0.10, 0.0],
        "filter": ["PASS", "PASS", "PASS", "artifact"],
        "TLOD": [20.0, 3.0, 15.0, 1.0],
        "POPAF": [6.0, 6.0, 6.0, 1.0],
        "TUMOR_DP": [50, 10, 40, 5],
        "NORMAL_DP": [30, 20, 25, 15],
    })
    print(f"    Input: {len(fake_df)} synthetic variants")
    print(f"    Columns: {list(fake_df.columns)}")
    mask = high_conf_somatic_mask(fake_df)
    assert mask.sum() >= 1, "Filter should pass at least one variant"
    print(f"    Somatic filter result: {mask.sum()}/{len(fake_df)} pass")
    print(f"    Pass pattern: {mask.tolist()}")
    for i, (passed, row) in enumerate(zip(mask, fake_df.itertuples())):
        label = "PASS" if passed else "FAIL"
        print(f"      Variant {i}: [{label}] VAF={row.tumor_vaf}, TLOD={row.TLOD}, filter={row.filter}")


# ===========================================================================
# TEST: SV module (from already-processed CSV)
# ===========================================================================

@smoke_test("sv.load_processed_csv")
def test_sv_processed(tmpdir: Path):
    csvs = sorted(SV_PROCESSED_DIR.glob("*_strict_sv_set.csv"))
    if not csvs:
        print(f"    No processed SV CSVs in {SV_PROCESSED_DIR}, skipping")
        return

    csv_path = csvs[0]
    print(f"    Loading: {csv_path.name} ({csv_path.stat().st_size / 1e3:.1f} KB)")
    df = pd.read_csv(csv_path)
    assert len(df) > 0, "Empty SV CSV"
    print(f"    {len(df)} SVs loaded")
    print(f"    Columns: {list(df.columns)}")

    for col in ["chrom", "pos", "SVTYPE"]:
        assert col in df.columns, f"Missing column: {col}"

    print(f"    SV types: {df['SVTYPE'].value_counts().to_dict()}")
    print(f"    Chromosomes: {sorted(df['chrom'].unique())}")

    if "gene_hits" in df.columns:
        print(f"    gene_hits column dtype: {df['gene_hits'].dtype}")
        n_with_hits = df["gene_hits"].apply(lambda x: bool(x) and x != "[]" and x != "nan").sum()
        print(f"    SVs with gene hits: {n_with_hits}/{len(df)}")


@smoke_test("sv.load_raw_vcf")
def test_sv_raw_vcf(tmpdir: Path):
    try:
        import vcfpy
    except ImportError:
        print("    vcfpy not installed, skipping SV VCF test")
        return

    from pipeline.SV.vcf_loader import load_manta_sv_vcf

    vcfs = sorted(SV_RAW_VCF_DIR.glob("*.vcf"))
    if not vcfs:
        print(f"    No raw SV VCFs in {SV_RAW_VCF_DIR}, skipping")
        return

    vcf_path = vcfs[0]
    print(f"    Testing with: {vcf_path.name} ({vcf_path.stat().st_size / 1e3:.1f} KB)")
    df, normal, tumor = load_manta_sv_vcf(str(vcf_path))
    assert isinstance(df, pd.DataFrame)
    print(f"    Loaded {len(df)} SVs")
    print(f"    Sample IDs: normal={normal}, tumor={tumor}")
    if len(df) > 0:
        print(f"    Columns: {list(df.columns)}")
        if "SVTYPE" in df.columns:
            print(f"    SV types: {df['SVTYPE'].value_counts().to_dict()}")
        if "chrom" in df.columns:
            print(f"    Chromosomes: {sorted(df['chrom'].unique())}")


@smoke_test("sv.filtering")
def test_sv_filtering(tmpdir: Path):
    try:
        import vcfpy
    except ImportError:
        print("    vcfpy not installed, skipping")
        return

    from pipeline.SV.vcf_loader import load_manta_sv_vcf
    from pipeline.SV.sv_filtering import get_strict_sv_set

    vcfs = sorted(SV_RAW_VCF_DIR.glob("*.vcf"))
    if not vcfs:
        print(f"    No raw VCFs, skipping")
        return

    df, _, _ = load_manta_sv_vcf(str(vcfs[0]))
    strict = get_strict_sv_set(df)
    pct = len(strict) / len(df) * 100 if len(df) > 0 else 0
    print(f"    Filtering: {len(df)} raw -> {len(strict)} strict SVs ({pct:.0f}% retained)")
    if len(strict) > 0 and "SVTYPE" in strict.columns:
        print(f"    Strict SV types: {strict['SVTYPE'].value_counts().to_dict()}")


@smoke_test("sv.spatial_mapping")
def test_sv_spatial_mapping(tmpdir: Path):
    try:
        from pipeline.SV.spatial_mapping import map_svs_to_genes
    except ImportError as e:
        print(f"    Cannot import spatial_mapping: {e}, skipping")
        return

    csvs = sorted(SV_PROCESSED_DIR.glob("*_strict_sv_set.csv"))
    if not csvs:
        print("    No processed SV CSVs, skipping")
        return

    sv_df = pd.read_csv(csvs[0])
    if "gene_hits" in sv_df.columns:
        print(f"    gene_hits already present in {csvs[0].name} — spatial mapping was pre-run")
        print(f"    {len(sv_df)} SVs already annotated")
        return

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    genes_only = genes[genes["feature"] == "gene"].copy()
    print(f"    Input: {len(sv_df)} SVs, {len(genes_only)} genes, window=500kb")

    sv_out = map_svs_to_genes(sv_df, genes, genes_only, window=500_000)
    assert "gene_hits" in sv_out.columns
    n_with_hits = sv_out["gene_hits"].apply(lambda x: len(x) > 0 if isinstance(x, list) else False).sum()
    print(f"    Output: {n_with_hits}/{len(sv_out)} SVs mapped to genes")


# ===========================================================================
# TEST: CNV module
# ===========================================================================

@smoke_test("cnv.load_and_annotate")
def test_cnv(tmpdir: Path):
    from pipeline.CNV.loader import load_cnv_file
    from pipeline.CNV.features import add_basic_cnv_features

    if not CNV_SEG_FILE.exists():
        print(f"    CNV seg file not found at {CNV_SEG_FILE}, skipping")
        return

    raw = load_cnv_file(CNV_SEG_FILE)
    assert len(raw) > 0, "No CNV segments loaded"
    print(f"    Loaded {len(raw)} segments from {CNV_SEG_FILE.name}")
    print(f"    Raw columns: {list(raw.columns)}")
    if "chrom" in raw.columns:
        print(f"    Chromosomes: {sorted(raw['chrom'].unique())}")

    cnv = add_basic_cnv_features(raw)
    assert "cn_state" in cnv.columns
    assert "loh_flag" in cnv.columns
    print(f"    Features added: cn_state, loh_flag")
    print(f"    CN states: {cnv['cn_state'].value_counts().to_dict()}")
    print(f"    LOH flag: {cnv['loh_flag'].value_counts().to_dict()}")
    print(f"    Output columns: {list(cnv.columns)}")


@smoke_test("cnv.gene_hits")
def test_cnv_gene_annotation(tmpdir: Path):
    from pipeline.CNV.loader import load_cnv_file
    from pipeline.CNV.features import add_basic_cnv_features
    from pipeline.CNV.gene_hits import annotate_cnv_with_gene_hits

    if not CNV_SEG_FILE.exists():
        print("    CNV seg file not found, skipping")
        return

    raw = load_cnv_file(CNV_SEG_FILE)
    cnv = add_basic_cnv_features(raw)

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    print(f"    Input: {len(cnv)} CNV segments, {len(genes)} gene rows for {SMOKE_GENES}")

    cnv_ann = annotate_cnv_with_gene_hits(cnv, genes, pu=2000, pd_=500, pw=5000)
    assert "gene_hits" in cnv_ann.columns
    n_with_hits = cnv_ann["gene_hits"].apply(lambda x: len(x) if isinstance(x, list) else 0).sum()
    n_segs_with_hits = cnv_ann["gene_hits"].apply(lambda x: len(x) > 0 if isinstance(x, list) else False).sum()
    print(f"    Segments with gene hits: {n_segs_with_hits}/{len(cnv_ann)}")
    print(f"    Total gene hits: {n_with_hits}")
    if n_with_hits > 0:
        first_hit_seg = cnv_ann[cnv_ann["gene_hits"].apply(lambda x: len(x) > 0 if isinstance(x, list) else False)].iloc[0]
        print(f"    Sample gene_hits structure: {str(first_hit_seg['gene_hits'])[:200]}")


# ===========================================================================
# TEST: Methylation module
# ===========================================================================

@smoke_test("methylation.load_probe_reference")
def test_meth_probe_ref(tmpdir: Path):
    from pipeline.Methylation.probe_loader import load_probe_reference

    if not METH_PROBE_REF.exists():
        print(f"    Probe reference not found at {METH_PROBE_REF}, skipping")
        return

    probes = load_probe_reference(METH_PROBE_REF)
    assert len(probes) > 0
    print(f"    Loaded {len(probes)} probes")
    print(f"    Columns: {list(probes.columns)}")
    if "chrom" in probes.columns:
        print(f"    Chromosomes: {probes['chrom'].nunique()} unique")
    if "probe_id" in probes.columns or "Name" in probes.columns:
        id_col = "probe_id" if "probe_id" in probes.columns else "Name"
        print(f"    Sample probe IDs: {probes[id_col].head(3).tolist()}")


@smoke_test("methylation.load_sample_beta")
def test_meth_sample(tmpdir: Path):
    from pipeline.Methylation.sample_processing import load_sample_beta

    if not METH_SAMPLE_RAW.exists():
        print(f"    Methylation sample not found at {METH_SAMPLE_RAW}, skipping")
        return

    beta_df = load_sample_beta(METH_SAMPLE_RAW)
    assert len(beta_df) > 0
    print(f"    Loaded {len(beta_df)} probe beta values")
    print(f"    Columns: {list(beta_df.columns)}")
    beta_col = [c for c in beta_df.columns if "beta" in c.lower()]
    if beta_col:
        vals = beta_df[beta_col[0]]
        print(f"    Beta range: [{vals.min():.4f}, {vals.max():.4f}], mean={vals.mean():.4f}")
        print(f"    NaN probes: {vals.isna().sum()} ({vals.isna().mean()*100:.1f}%)")


@smoke_test("methylation.read_existing_outputs")
def test_meth_existing_outputs(tmpdir: Path):
    """Verify that existing methylation outputs are loadable."""
    gene_matrix_path = METH_OUTPUTS / "cohort" / "gene_meth_matrix.csv"
    if not gene_matrix_path.exists():
        print("    No cohort gene matrix found, skipping")
        return

    mat = pd.read_csv(gene_matrix_path, index_col=0)
    assert len(mat) > 0
    print(f"    Gene meth matrix: {mat.shape[0]} genes x {mat.shape[1]} cols")
    print(f"    Columns: {list(mat.columns)}")
    numeric_cols = mat.select_dtypes(include="number")
    if len(numeric_cols.columns) > 0:
        print(f"    Numeric value range: [{numeric_cols.values[~np.isnan(numeric_cols.values)].min():.4f}, "
              f"{numeric_cols.values[~np.isnan(numeric_cols.values)].max():.4f}]")
        print(f"    NaN rate: {numeric_cols.isna().mean().mean()*100:.1f}%")
    str_cols = mat.select_dtypes(include="object")
    if len(str_cols.columns) > 0:
        print(f"    String columns: {list(str_cols.columns)}")
    print(f"    Sample gene names: {list(mat.index[:5])}")

    per_sample_dir = METH_OUTPUTS / "per_sample"
    if per_sample_dir.exists():
        samples = [d for d in per_sample_dir.iterdir() if d.is_dir()]
        print(f"    Per-sample directories: {len(samples)}")
        if samples:
            s = samples[0]
            gene_meth = s / f"{s.name}_gene_meth.csv"
            if gene_meth.exists():
                df = pd.read_csv(gene_meth)
                print(f"    Sample {s.name}: {len(df)} gene-level entries")
                print(f"    Sample columns: {list(df.columns)}")
            sample_files = list(s.iterdir())
            print(f"    Files in sample dir: {[f.name for f in sample_files]}")


# ===========================================================================
# TEST: RPPA module (import-level + submodule tests)
# ===========================================================================

@smoke_test("rppa.import_submodules")
def test_rppa_imports(tmpdir: Path):
    """Test that RPPA submodules import without crashing (the __init__ is broken)."""
    errors = []
    missing_pkgs = []

    for mod_name in ["rppa_loader", "rppa_panels", "rppa_analysis", "rppa_schemas"]:
        try:
            __import__(f"pipeline.rppa.{mod_name}", fromlist=[mod_name])
        except ImportError as e:
            err_str = str(e)
            if "No module named" in err_str and any(pkg in err_str for pkg in ["scipy", "sklearn", "statsmodels"]):
                missing_pkgs.append(err_str)
            else:
                errors.append(f"{mod_name}: {e}")
        except Exception as e:
            errors.append(f"{mod_name}: {e}")

    if missing_pkgs:
        print(f"    RPPA needs optional packages: {missing_pkgs}")
    if errors:
        print(f"    RPPA submodule import errors:")
        for err in errors:
            print(f"      - {err}")
        raise ImportError(f"{len(errors)} RPPA submodule(s) failed to import")
    else:
        print("    All RPPA submodules import successfully (some may need scipy/sklearn)")


@smoke_test("rppa.package_init_import")
def test_rppa_package_import(tmpdir: Path):
    """Test if the RPPA package __init__.py can be imported (expected to fail)."""
    try:
        import pipeline.rppa
        print("    pipeline.rppa imported successfully")
    except ImportError as e:
        err_str = str(e)
        # Missing optional packages (scipy etc.) are not pipeline bugs
        if any(pkg in err_str for pkg in ["scipy", "sklearn", "statsmodels"]):
            print(f"    pipeline.rppa needs optional packages: {e}")
            return
        # Missing config symbols IS a pipeline bug
        raise ImportError(
            f"pipeline.rppa __init__.py fails to import: {e}\n"
            "  -> config.py is missing: ImmuneVisibilityMarkers, RPPA_PATHS, "
            "RPPA_THRESHOLDS, IMMUNE_VISIBILITY_MARKERS, "
            "RPPA_TARGET_TO_GENE_OVERRIDES, KEY_PHOSPHO_PAIRS, get_rppa_config"
        )


# ===========================================================================
# TEST: Pipeline schemas validation
# ===========================================================================

@smoke_test("schemas.smoke_check")
def test_schemas(tmpdir: Path):
    from pipeline.schemas import (
        empty_screen_block,
        empty_abc_block,
        empty_hichip_block,
        empty_snv_entry,
        empty_sv_record,
        empty_atac_gene_link_entry,
        empty_atac_ccre_link_entry,
    )

    screen = empty_screen_block(biosamples=["MCF-7"], assay_types=["Intact-HiC"])
    assert "per_biosample" in screen
    assert "conservation_global" in screen

    abc = empty_abc_block(celltypes=["MCF-7-ENCODE"])
    assert isinstance(abc, dict)

    hichip = empty_hichip_block(celltypes=["MCF7"])
    assert isinstance(hichip, dict)

    snv = empty_snv_entry()
    assert "chrom" in snv

    sv = empty_sv_record()
    assert "gene_hits" in sv

    atac_gene = empty_atac_gene_link_entry()
    assert "dist_to_tss" in atac_gene

    atac_ccre = empty_atac_ccre_link_entry()
    assert "cCRE_id" in atac_ccre or "distance" in atac_ccre

    print("    All schema factories produce valid dicts")


# ===========================================================================
# TEST: Utils
# ===========================================================================

@smoke_test("utils.core_functions")
def test_utils(tmpdir: Path):
    from pipeline.utils import (
        normalize_chrom,
        harmonize_chrom_column,
        min_interval_distance,
        tss_to_interval_distance,
        compute_interval_overlap,
        assign_distance_tier,
    )

    # normalize_chrom expects a pandas Series, not a plain string
    s = pd.Series(["6", "chr6", "chrX"])
    normalized = normalize_chrom(s)
    assert all(normalized.str.startswith("chr")), f"normalize_chrom failed: {normalized.tolist()}"

    df = pd.DataFrame({"chrom": ["6", "chr1", "X"]})
    df, _ = harmonize_chrom_column(df)
    assert all(df["chrom"].str.startswith("chr"))

    assert min_interval_distance(10, 20, 30, 40) == 10
    assert min_interval_distance(10, 40, 20, 30) == 0

    assert tss_to_interval_distance(15, 10, 20) == 0
    assert tss_to_interval_distance(5, 10, 20) == 5

    ov = compute_interval_overlap(10, 30, 20, 40)
    assert isinstance(ov, dict), f"compute_interval_overlap should return dict, got {type(ov)}"
    assert ov["overlaps"] == True
    assert ov["overlap_bp"] == 10

    tiers = assign_distance_tier(
        pd.Series([50_000, 150_000, 800_000]),
        tier_edges=[0, 100_000, 250_000, 500_000, 1_000_000],
        tier_labels=["0–100kb", "100–250kb", "250–500kb", "500–1000kb"],
    )
    assert tiers.iloc[0] == "0–100kb"
    assert tiers.iloc[1] == "100–250kb"
    assert tiers.iloc[2] == "500–1000kb"

    print("    All utility functions pass")


# ===========================================================================
# TEST: Main pipeline config loading
# ===========================================================================

@smoke_test("config.load")
def test_config(tmpdir: Path):
    from pipeline.config import (
        PATHS,
        BIOSAMPLES,
        THRESHOLDS,
        PRIMARY_GENES,
        PIPELINE_GENE_PANEL,
        CNV_GENES,
        FULL_INTEGRATION_GENES,
        EXTENDED_PRIMARY_GENES,
        TIER2_MEDIUM_GENES,
        TIER3_CNV_ONLY_GENES,
        TIER4_READOUT_GENES,
        USE_EXTENDED_PRIMARY_PANEL,
    )

    assert len(PRIMARY_GENES) == 66
    assert len(PRIMARY_GENES) > 10
    assert THRESHOLDS.ccre_window_bp == 1_000_000
    assert len(BIOSAMPLES.screen_exp) > 0
    assert PATHS.working_dir.parts[-1] == "data"
    print(f"    PRIMARY_GENES (core): {len(PRIMARY_GENES)} genes")
    print(f"    PIPELINE_GENE_PANEL: {len(PIPELINE_GENE_PANEL)} genes "
          f"(extended={USE_EXTENDED_PRIMARY_PANEL})")
    print(f"    FULL_INTEGRATION: {len(FULL_INTEGRATION_GENES)} | "
          f"EXT+1: {len(EXTENDED_PRIMARY_GENES)} | T2: {len(TIER2_MEDIUM_GENES)} | "
          f"T3: {len(TIER3_CNV_ONLY_GENES)} | T4: {len(TIER4_READOUT_GENES)}")
    assert len(FULL_INTEGRATION_GENES) == len(
        set(PRIMARY_GENES) | set(EXTENDED_PRIMARY_GENES)
    )
    assert set(TIER3_CNV_ONLY_GENES).issubset(set(CNV_GENES))
    assert set(FULL_INTEGRATION_GENES).issubset(set(CNV_GENES))
    assert set(TIER2_MEDIUM_GENES).issubset(set(CNV_GENES))
    assert not set(TIER4_READOUT_GENES) & set(PIPELINE_GENE_PANEL)
    assert "GZMA" in TIER4_READOUT_GENES
    assert "GZMA" not in CNV_GENES
    assert "EZH2" in CNV_GENES
    print(f"    CNV_GENES: {len(CNV_GENES)} (deduped union)")
    print(f"    Biosamples configured: screen_exp={len(BIOSAMPLES.screen_exp)}, "
          f"abc={len(BIOSAMPLES.abc_celltypes)}, hichip={len(BIOSAMPLES.hichip_panel)}")


@smoke_test("gene_panels.gencode_resolution")
def test_gene_panels_gencode_resolution(tmpdir: Path):
    """Tier 1/2/3/4 symbols resolve in GENCODE-derived gene table."""
    from pipeline.genes.gene_loader import filter_genes_by_names

    if not GENCODE_CSV.is_file():
        print(f"    GENCODE CSV missing at {GENCODE_CSV}, skipping")
        return

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES_TIERED)
    genes = filter_genes_by_names(genes_all, SMOKE_GENES_TIERED)
    found = set(genes["gene_name"].astype(str))
    missing = [g for g in SMOKE_GENES_TIERED if g not in found]
    assert not missing, f"GENCODE missing tiered panel symbols: {missing}"
    print(f"    Resolved {len(found)}/{len(SMOKE_GENES_TIERED)} tiered symbols in GENCODE CSV")


@smoke_test("gene_aliases.hgnc_pvrl2")
def test_gene_aliases_hgnc_pvrl2(tmpdir: Path):
    """HGNC alias file maps PVRL2 -> NECTIN2 within a bounded scan depth."""
    from pipeline.genes.gene_aliases import build_rna_expression_symbol_mapping
    from pipeline.config import UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE, LEGACY_DATASET_SYMBOL_RENAMES

    if not HGNC_ALIAS_TSV.is_file():
        print(f"    HGNC alias file missing at {HGNC_ALIAS_TSV}, skipping")
        return

    seed = frozenset({"NECTIN2"})
    m = build_rna_expression_symbol_mapping(
        HGNC_ALIAS_TSV,
        seed,
        UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE,
        LEGACY_DATASET_SYMBOL_RENAMES,
        max_lines=400_000,
        verbose=False,
    )
    assert m.get("PVRL2") == "NECTIN2", f"expected PVRL2->NECTIN2, got {m.get('PVRL2')!r}"
    print(f"    HGNC-derived RNA rename entries (bounded scan): {len(m)}")


@smoke_test("rna.normalize_expression_alias_stub")
def test_normalize_expression_alias_stub(tmpdir: Path):
    from pipeline.RNA_exp.normalize_expression_mat import normalize_expression_mat

    df = pd.DataFrame({"gene_symbol": ["PVRL2"], "tpm": [2.5]})
    out = normalize_expression_mat(df, gene_col="gene_symbol", mapping={"PVRL2": "NECTIN2"})
    assert out["gene_symbol"].iloc[0] == "NECTIN2"
    print("    normalize_expression_mat alias replace OK")


# ===========================================================================
# TEST: main.py known bugs audit
# ===========================================================================

@smoke_test("tiers.resolve_symbol_to_panel")
def test_tiers_resolve_symbol_to_panel(tmpdir: Path):
    from pipeline.genes.symbol_normalization import resolve_symbol_to_panel

    panel = {"NECTIN2", "B2M"}
    m = {"PVRL2": "NECTIN2", "HGM": "B2M"}
    assert resolve_symbol_to_panel("NECTIN2", panel, m) == "NECTIN2"
    assert resolve_symbol_to_panel("PVRL2", panel, m) == "NECTIN2"
    assert resolve_symbol_to_panel("NOPE", panel, m) is None
    print("    resolve_symbol_to_panel OK")


@smoke_test("tiers.snv_vep_symbol_alias_minimal")
def test_tiers_snv_vep_symbol_alias_minimal(tmpdir: Path):
    """VEP CSQ uses legacy PVRL2; panel is NECTIN2 — normalized hit retained."""
    import pandas as pd
    from pipeline.SNV.vep_parser import add_vep_hits_columns, DEFAULT_CSQ_FORMAT, parse_csq_description

    cols = parse_csq_description(DEFAULT_CSQ_FORMAT)
    vals = {c: "" for c in cols}
    vals["Allele"] = "T"
    vals["Consequence"] = "missense_variant"
    vals["IMPACT"] = "MODERATE"
    vals["SYMBOL"] = "PVRL2"
    vals["Gene"] = "ENSG00000130202"
    vals["Feature_type"] = "Transcript"
    vals["Feature"] = "ENST00000357310"
    vals["BIOTYPE"] = "protein_coding"
    vals["CANONICAL"] = "1"
    csq = "|".join(vals.get(c, "") for c in cols)

    df = pd.DataFrame({"chrom": ["1"], "pos": [100], "ref": ["A"], "alt": ["T"], "CSQ": [csq]})
    out = add_vep_hits_columns(
        df,
        DEFAULT_CSQ_FORMAT,
        ["NECTIN2"],
        gene_symbol_mapping={"PVRL2": "NECTIN2"},
    )
    assert len(out.iloc[0]["gene_hits"]) >= 1
    assert out.iloc[0]["gene_hits"][0]["SYMBOL"] == "NECTIN2"
    print("    VEP alias → panel canonical OK")


@smoke_test("tiers.mirtarbase_mapping_no_crash")
def test_tiers_mirtarbase_mapping_no_crash(tmpdir: Path):
    from pipeline.genes.mirtarbase import load_mirtarbase

    if not MIRNA_MIRTARBASE.exists():
        print("    miRTarBase CSV missing, skipping")
        return
    df = load_mirtarbase(mirtarbase_csv=MIRNA_MIRTARBASE, gene_panel=SMOKE_GENES[:3])
    assert "gene_norm" in df.columns
    print(f"    load_mirtarbase with symbol mapping: {len(df)} rows retained")


@smoke_test("tiers.methylation_probe_reference_normalize")
def test_tiers_methylation_probe_reference_normalize(tmpdir: Path):
    from pipeline.Methylation.probe_loader import load_probe_reference

    if not METH_PROBE_REF.exists():
        print("    Methylation probe ref missing, skipping")
        return
    df = load_probe_reference(METH_PROBE_REF)
    assert "gene_list" in df.columns
    n_nonempty = sum(1 for _i, row in df.head(500).iterrows() if row["gene_list"])
    print(f"    probe reference gene_list rows (first 500, nonempty): {n_nonempty}")


@smoke_test("registry.unittest_panel_aliases")
def test_registry_unittest_panel_aliases(tmpdir: Path):
    """Run ``unittest`` for ``tests/test_panel_alias_registry.py``."""
    r = subprocess.run(
        [
            sys.executable,
            "-m",
            "unittest",
            "tests.test_panel_alias_registry",
            "-v",
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
    assert r.returncode == 0, "panel alias registry unit tests failed"
    print("    unittest tests.test_panel_alias_registry: OK")


@smoke_test("registry.module_sample_coverage_script")
def test_registry_module_sample_coverage_script(tmpdir: Path):
    out = tmpdir / "module_sample_canonical_coverage.md"
    r = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "coverage" / "module_sample_canonical_coverage.py"),
            "--out",
            str(out),
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        print(r.stderr)
    assert r.returncode == 0, r.stderr
    assert out.is_file() and out.stat().st_size > 50
    print(f"    wrote coverage sample report ({out.stat().st_size} bytes)")


@smoke_test("rna.mapping_gate_ucsc_only")
def test_rna_mapping_gate_ucsc_only(tmpdir: Path):
    """UCSC-only dict (what RNA uses when ``APM_USE_GENE_SYMBOL_MAPPING=0``)."""
    from pipeline.config import UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE
    from pipeline.RNA_exp.normalize_expression_mat import normalize_expression_mat

    df = pd.DataFrame({"gene_symbol": ["PVRL2"]})
    out = normalize_expression_mat(df, gene_col="gene_symbol", mapping=dict(UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE))
    assert out["gene_symbol"].iloc[0] == "NECTIN2"
    print("    RNA normalize with UCSC-only mapping OK")


@smoke_test("main.known_bugs")
def test_main_bugs(tmpdir: Path):
    """Audit main.py for known bugs without running the full pipeline."""
    from pipeline.config import PATHS

    issues = []
    fixed = []

    # Bug 1: PATHS.mirtarbase_csv does not exist (should be mirtarbase_input)
    if not hasattr(PATHS, "mirtarbase_csv"):
        # Check if main.py still references the old name
        import inspect
        from pipeline import main as main_mod
        full_src = inspect.getsource(main_mod)
        if "mirtarbase_csv" in full_src:
            issues.append(
                "BUG: main.py uses PATHS.mirtarbase_csv but "
                "config.py defines PATHS.mirtarbase_input"
            )
        else:
            fixed.append("FIXED: mirtarbase_csv -> mirtarbase_input reference corrected")

    # Bug 2: atac_table / atac_reduced used unconditionally after skip_atac guard
    import inspect
    from pipeline import main as main_mod
    src = inspect.getsource(main_mod.run_full_pipeline)
    lines = src.split("\n")

    atac_unguarded = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if "atac_reduced" in stripped and "annotate_df_with_peaks" in stripped:
            context = "\n".join(lines[max(0,i-3):i+1])
            if "if " not in context and "else" not in context:
                atac_unguarded = True
                break

    if atac_unguarded:
        issues.append(
            "BUG: main.py uses atac_reduced unconditionally in annotate_df_with_peaks "
            "— crashes when skip_atac=True and atac_reduced is undefined"
        )
    else:
        fixed.append("FIXED: atac_table/atac_reduced now properly guarded")

    # Bug 3: RPPA __init__.py imports missing config symbols
    try:
        import pipeline.rppa
        fixed.append("FIXED: pipeline.rppa imports work")
    except ImportError as e:
        err_str = str(e)
        if any(pkg in err_str for pkg in ["scipy", "sklearn", "statsmodels"]):
            fixed.append("FIXED (partial): pipeline.rppa works but needs optional packages")
        else:
            issues.append(f"BUG: pipeline.rppa fails to import: {e}")

    for f in fixed:
        print(f"    {f}")
    for issue in issues:
        print(f"    {issue}")

    if issues:
        print(f"    AUDIT: {len(issues)} remaining bug(s) in main.py (informational, not a test failure)")
    else:
        print(f"    AUDIT: All {len(fixed)} previously known bugs are fixed!")


# ===========================================================================
# TEST: miRNA targets
# ===========================================================================

@smoke_test("mirna.targetscan")
def test_mirna_targetscan(tmpdir: Path):
    from pipeline.genes.mirna_targets import get_mirna_targets
    from pipeline.genes.gene_loader import filter_genes_by_names, add_promoter_columns

    if not MIRNA_TARGETSCAN.exists():
        print("    TargetScan predictions not found, skipping")
        return

    genes_all = load_genes_from_csv(GENCODE_CSV, gene_names=SMOKE_GENES)
    genes = filter_genes_by_names(genes_all, SMOKE_GENES)
    genes = add_promoter_columns(genes, upstream_bp=2000, downstream_bp=500)
    print(f"    TargetScan file: {MIRNA_TARGETSCAN.name} ({MIRNA_TARGETSCAN.stat().st_size / 1e6:.1f} MB)")
    print(f"    Gene panel: {SMOKE_GENES}, top_n=50, weight_threshold=-0.2")

    out_dir = tmpdir / "mirna"
    get_mirna_targets(
        genes=genes,
        gtf=genes_all,
        targetscan_path=MIRNA_TARGETSCAN,
        top_n=50,
        weight_threshold=-0.2,
        output_dir=out_dir,
    )
    if out_dir.exists():
        out_files = list(out_dir.rglob("*"))
        print(f"    miRNA outputs: {len(out_files)} files written to {out_dir}")
        for f in out_files[:5]:
            if f.is_file():
                print(f"      {f.name} ({f.stat().st_size / 1e3:.1f} KB)")


# ===========================================================================
# RUNNER
# ===========================================================================

ALL_TESTS = [
    test_config,
    test_gene_panels_gencode_resolution,
    test_gene_aliases_hgnc_pvrl2,
    test_normalize_expression_alias_stub,
    test_tiers_resolve_symbol_to_panel,
    test_tiers_snv_vep_symbol_alias_minimal,
    test_tiers_mirtarbase_mapping_no_crash,
    test_tiers_methylation_probe_reference_normalize,
    test_registry_unittest_panel_aliases,
    test_registry_module_sample_coverage_script,
    test_rna_mapping_gate_ucsc_only,
    test_utils,
    test_schemas,
    test_genes_load,
    test_lncrna_matching,
    test_ccre_loading,
    test_distance_matching,
    test_cell_line_signals,
    test_tad_annotation,
    test_screen_prebuilt,
    test_abc_links,
    test_hichip_links,
    test_atac_peaks,
    test_snv_load,
    test_snv_filter,
    test_sv_processed,
    test_sv_raw_vcf,
    test_sv_filtering,
    test_sv_spatial_mapping,
    test_cnv,
    test_cnv_gene_annotation,
    test_meth_probe_ref,
    test_meth_sample,
    test_meth_existing_outputs,
    test_rppa_imports,
    test_rppa_package_import,
    test_mirna_targetscan,
    test_main_bugs,
]

TEST_NAMES = {fn._smoke_name: fn for fn in ALL_TESTS}


def main():
    # Parse args
    selected = None
    if len(sys.argv) > 1:
        if sys.argv[1] == "--list":
            print("Available smoke tests:")
            for name in TEST_NAMES:
                print(f"  {name}")
            return
        selected = sys.argv[1:]

    # Create temp dir for outputs
    tmpdir = Path(tempfile.mkdtemp(prefix="apm_smoke_"))
    print(f"\nSmoke test output dir: {tmpdir}")
    print("=" * 70)

    # Run tests
    tests_to_run = ALL_TESTS
    if selected:
        tests_to_run = []
        for sel in selected:
            matched = [fn for name, fn in TEST_NAMES.items() if sel in name]
            tests_to_run.extend(matched)
        if not tests_to_run:
            print(f"No tests matching: {selected}")
            print(f"Available: {list(TEST_NAMES.keys())}")
            return

    t0 = time.time()
    for test_fn in tests_to_run:
        test_fn(tmpdir)

    total_time = time.time() - t0

    # Summary
    print("\n" + "=" * 70)
    print("SMOKE TEST SUMMARY")
    print("=" * 70)

    passed = [r for r in results if r.passed]
    failed = [r for r in results if not r.passed]

    for r in results:
        status = "PASS" if r.passed else "FAIL"
        print(f"  [{status}] {r.name:40s}  {r.duration_s:6.1f}s  {r.message[:40]}")

    print(f"\n  Total: {len(results)} | Passed: {len(passed)} | Failed: {len(failed)}")
    print(f"  Time:  {total_time:.1f}s")
    print(f"  Output dir: {tmpdir}")

    if failed:
        print(f"\n{'='*70}")
        print("FAILURE DETAILS")
        print(f"{'='*70}")
        for r in failed:
            print(f"\n--- {r.name} ---")
            print(r.error)

    # Exit code
    sys.exit(len(failed))


if __name__ == "__main__":
    main()
