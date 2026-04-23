# CNV Pipeline Module

## Overview

CNV annotation sub-package for the APM regulatory pipeline.
Drop `cnv/` into `pipeline/` alongside existing modules.

## Structure

```
pipeline/cnv/
  __init__.py    Public API exports
  loader.py      File I/O + sample ID resolution
  features.py    CN state, LOH flag, segment length
  geometry.py    Signed distance + region classification
  gene_hits.py   Gene / lncRNA / miRNA annotation
  elem_hits.py   cCRE annotation
  runner.py          Orchestrator (single + batch)
  gene_level_ascat.py  GDC ASCAT3 gene-level TSV → per-vial gene tables (LOH, window min CN, promoters)
  CNV_MODULE.md      This file
```

## ASCAT3 gene-level copy number (GDC)

Files whose names contain ``ascat3.gene_level_copy_number`` or ``gene_level_copy_number`` (for example ``TCGA-BRCA.<uuid>.ascat3.gene_level_copy_number.v36.tsv``) in ``PATHS.cnv_dir`` are **not** segment files. When present, ``process_cnv_directory`` reads them with ``load_ascat_gene_level_tsv``, resolves the tumor vial via the same ``samples.tsv`` manifest, and writes enriched tables under ``gene_tables_root`` (typically ``PATHS.cnv_gene_tables_dir``) as ``{sample_vial}_{THRESHOLDS.cnv_ascat3_gene_table_stem}.csv`` (default stem ``cnv_gene_calls_ascat3``), plus ``cnv_gene_calls_all_samples_ascat3.csv``. Segment-derived gene calls remain ``{sid}_cnv_gene_calls.csv`` when ``*.seg`` / compatible segment TSVs are processed.

Columns added beyond the GDC file include ``cn_state``, ``cn_minor`` / ``cn_major`` / ``allele_delta``, ``loh_flag`` (minor allele 0 with major ≥ 1), ``loh_only`` (diploid total CN with LOH), ``regulatory_window_min_cn`` / ``regulatory_window_bp``, and optional ``tss`` / ``prom_start`` / ``prom_end`` / ``strand`` when ``genes_path`` carries those fields (e.g. ``cnv_genes.csv`` after promoter annotation).

## Config: add to PathConfig in config.py

```python
cnv_dir: Path = Path(".../CNV/raw")
cnv_genes: Path = Path(".../gencode.v49.annotation.gtf.csv")
cnv_annotations_path: Path = Path(".../CNV/samples.tsv")
cnv_output_dir: Path = Path(".../CNV/annotated")
```
## How to Run

### A) Batch (replaces notebook)
```python
from pipeline.CNV import process_cnv_directory
from pipeline.config import PATHS
process_cnv_directory(
    cnv_dir=str(PATHS.cnv_dir),
    genes_path=str(PATHS.cnv_genes),
    lncrnas_path=str(PATHS.lncrna_csv),
    mirna_path=str(PATHS.mirna_path),
    elements_path=str(PATHS.regulatory_elements_table),
    ann_path=str(PATHS.cnv_annotations_path),
    out_dir=str(PATHS.cnv_output_dir))
```

### B) Single-sample
```python
from pipeline.cnv import load_cnv_file, annotate_single_sample
raw = load_cnv_file("path/to/sample.seg.txt")
result = annotate_single_sample(raw, genes, lncrnas, mirna, elements)
```

## Output: gene_hits, lncRNAs_hits, mir_hits, elem_hits (list[dict])
