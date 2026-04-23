# Methylation Module Documentation

## Overview

The methylation module integrates DNA methylation data (450K/EPIC array) with the regulatory pipeline. It provides:

1. **Probe Reference Table** - Master annotation of all probes with genomic context
2. **Per-Sample Processing** - Beta value loading, validation, and enrichment
3. **Feature-Level Aggregation** - Gene, lncRNA, and cCRE methylation summaries
4. **Cohort Matrices** - Features × Samples matrices for downstream analysis

---

## Module Structure

```
methylation/
├── __init__.py              # Clean exports
├── meth_schemas.py          # Schema definitions & validation
├── probe_loader.py          # Probe reference loading & annotation
├── sample_processing.py     # Per-sample beta handling
├── aggregation.py           # Gene/lncRNA/cCRE aggregation
├── methylation_table.py     # Main orchestrator
└── config_additions.py      # Config additions reference
```

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        METHYLATION PIPELINE                              │
└─────────────────────────────────────────────────────────────────────────┘

STEP 1: Build Probe Reference (once)
    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
    │ GDC Probe   │     │   Genes     │     │   cCREs     │
    │ Reference   │     │ (GENCODE)   │     │  (ENCODE)   │
    └──────┬──────┘     └──────┬──────┘     └──────┬──────┘
           │                   │                   │
           ▼                   ▼                   ▼
    ┌─────────────────────────────────────────────────────┐
    │              Probe Annotation Pipeline               │
    │  • Promoter overlap (gene panel)                    │
    │  • Gene body overlap                                │
    │  • cCRE overlap                                     │
    │  • lncRNA promoter overlap                          │
    │  • ATAC peak overlap                                │
    │  • TAD domain context (per biosample)               │
    └─────────────────────────────────────────────────────┘
           │
           ▼
    ┌─────────────────────────────────────────────────────┐
    │         probe_annotations.parquet                    │
    │  ~450K-850K probes with full genomic context        │
    └─────────────────────────────────────────────────────┘

STEP 2: Per-Sample Processing
    ┌─────────────┐     ┌─────────────────────┐
    │ Sample Beta │     │  Probe Reference    │
    │   File      │     │   (from Step 1)     │
    └──────┬──────┘     └──────────┬──────────┘
           │                       │
           ▼                       ▼
    ┌─────────────────────────────────────────────────────┐
    │           Sample Processing Pipeline                 │
    │  • Load & validate beta values                      │
    │  • Compute M-values                                 │
    │  • Enrich with reference annotations                │
    │  • Aggregate to genes/lncRNAs/cCREs                 │
    └─────────────────────────────────────────────────────┘
           │
           ├──► {sample}_probes.parquet
           ├──► {sample}_gene_meth.csv
           ├──► {sample}_lncrna_meth.csv
           └──► {sample}_ccre_meth.csv

STEP 3: Cohort Matrices
    ┌─────────────────────────────────────────────────────┐
    │         All Sample Aggregations                      │
    └─────────────────────────────────────────────────────┘
           │
           ▼
    ┌─────────────────────────────────────────────────────┐
    │  • gene_meth_matrix.parquet  (Genes × Samples)      │
    │  • lncrna_meth_matrix.csv    (lncRNAs × Samples)    │
    │  • ccre_meth_matrix.parquet  (cCREs × Samples)      │
    │  • sample_qc_summary.csv                            │
    └─────────────────────────────────────────────────────┘
```

---

## Probe Reference Annotations

Each probe gets these annotations:

### Core Coordinates
| Column | Type | Description |
|--------|------|-------------|
| `probeID` | str | Illumina probe ID (e.g., cg00000029) |
| `chrom` | str | Chromosome (chr1-chr22, chrX, chrY) |
| `start` | int | CpG start position |
| `end` | int | CpG end position |
| `strand` | str | Strand (+/-) |

### CGI Context
| Column | Type | Description |
|--------|------|-------------|
| `in_CGI` | bool | Probe overlaps CpG island |
| `CGI_context` | str | Island/N_Shore/S_Shore/N_Shelf/S_Shelf/OpenSea |
| `distToTSS` | int | Distance to nearest TSS (from GDC reference) |

### Gene Context
| Column | Type | Description |
|--------|------|-------------|
| `in_promoter` | bool | Overlaps any panel gene promoter |
| `promoter_genes` | list[str] | Panel genes whose promoters overlap |
| `promoter_gene_ids` | list[str] | Corresponding ENSEMBL IDs |
| `in_gene_body` | bool | Overlaps any panel gene body (excl. promoter) |
| `gene_body_genes` | list[str] | Panel genes whose bodies overlap |

### cCRE Context
| Column | Type | Description |
|--------|------|-------------|
| `overlapping_ccres` | list[str] | cCRE IDs that overlap this probe |
| `ccre_types` | list[str] | Types of overlapping cCREs (pELS, dELS, etc.) |
| `n_overlapping_ccres` | int | Count of overlapping cCREs |

### lncRNA Context
| Column | Type | Description |
|--------|------|-------------|
| `in_lncrna_promoter` | bool | Overlaps any matched lncRNA promoter |
| `lncrna_promoter_genes` | list[str] | lncRNA names whose promoters overlap |

### ATAC Context
| Column | Type | Description |
|--------|------|-------------|
| `overlapping_atac_peaks` | list[str] | ATAC peak IDs that overlap |
| `n_overlapping_atac` | int | Count of overlapping peaks |

### TAD Context
| Column | Type | Description |
|--------|------|-------------|
| `TAD_domains` | dict | Per-biosample TAD context. **Build default:** *slim* payload (`APM_METH_TAD_SLIM` defaults on; set `APM_METH_TAD_SLIM=0` for full gene-style `domains` + boundary blobs). |

When **reusing** `reference/probe_annotations.parquet`, the pipeline loads only columns needed for enrichment/aggregation (see `METHYLATION_PROBE_REFERENCE_LOAD_COLUMNS` in `methylation_table.py`), so an old file with huge full `TAD_domains` does not OOM on read. Set `APM_METH_LOAD_FULL_PROBE_REFERENCE=1` if you truly need every column in memory.

Parquet is written in slices (`APM_METH_PROBE_PARQUET_CHUNK`, default 50000 rows). **Free disk space** before building the reference.

---

## Gene-Level Methylation Metrics

For each gene in the panel:

### Promoter Metrics
| Column | Type | Description |
|--------|------|-------------|
| `promoter_beta_mean` | float | Mean beta across promoter probes |
| `promoter_beta_median` | float | Median beta |
| `promoter_beta_std` | float | Standard deviation |
| `promoter_beta_min` | float | Minimum beta |
| `promoter_beta_max` | float | Maximum beta |
| `promoter_beta_range` | float | Max - Min |
| `promoter_n_probes` | int | Number of probes in promoter |
| `promoter_frac_hypermeth` | float | Fraction of probes with beta > 0.7 |
| `promoter_frac_hypometh` | float | Fraction of probes with beta < 0.3 |
| `promoter_CGI_beta_mean` | float | Mean beta of CGI probes only |
| `promoter_CGI_n_probes` | int | Number of CGI probes |
| `promoter_shore_beta_mean` | float | Mean beta of shore probes only |

### Gene Body Metrics
| Column | Type | Description |
|--------|------|-------------|
| `gene_body_beta_mean` | float | Mean beta across gene body probes |
| `gene_body_n_probes` | int | Number of probes in gene body |

---

## cCRE-Level Methylation Metrics

For each cCRE with overlapping probes:

| Column | Type | Description |
|--------|------|-------------|
| `cCRE_id` | str | cCRE identifier |
| `ccre_beta_mean` | float | Mean beta across cCRE probes |
| `ccre_beta_median` | float | Median beta |
| `ccre_beta_std` | float | Standard deviation |
| `ccre_n_probes` | int | Number of probes in cCRE |
| `ccre_frac_hypermeth` | float | Fraction hypermethylated |
| `ccre_frac_hypometh` | float | Fraction hypomethylated |
| `ccre_CGI_overlap` | bool | cCRE contains CGI probes |
| `ccre_CGI_beta_mean` | float | Mean beta of CGI probes only |
| `ccre_CGI_n_probes` | int | Number of CGI probes |

---

## Usage Examples

### 1. Build Probe Reference (One-Time)

```python
from methylation import build_probe_reference_table
from pipeline.config import PATHS, THRESHOLDS, PRIMARY_GENES

# Load your genes and cCREs (already loaded in main pipeline)
genes = pd.read_csv(PATHS.gencode_gtf_csv)
ccres = pd.read_csv(PATHS.ccre_csv)
lncrnas = pd.read_csv(PATHS.lncrna_csv)

# Build annotated probe reference
probe_reference = build_probe_reference_table(
    probe_reference_path=PATHS.methylation_probe_reference,
    genes=genes,
    ccres=ccres,
    lncrnas=lncrnas,
    gene_panel=PRIMARY_GENES,
    upstream_bp=THRESHOLDS.promoter_upstream_bp,
    downstream_bp=THRESHOLDS.promoter_downstream_bp,
    output_path=PATHS.working_dir / "methylation/reference/probe_annotations.parquet",
)
```

### 2. Process Single Sample

```python
from methylation import build_sample_methylation_tables

results = build_sample_methylation_tables(
    sample_path=Path("/path/to/sample_beta.tsv"),
    probe_reference=probe_reference,
    sample_id="TCGA-A8-A09K",
    gene_panel=PRIMARY_GENES,
    lncrna_panel=matched_lncrna_names,
    output_dir=Path("/output/methylation/per_sample/TCGA-A8-A09K"),
)

# Access results
gene_meth = results["genes"]      # Gene-level aggregation
ccre_meth = results["ccres"]      # cCRE-level aggregation
probes = results["probes"]        # Per-probe data with annotations
```

### 3. Build Cohort Matrices

```python
from methylation import build_cohort_matrices

# Sample manifest: DataFrame with columns [sample_id, beta_path]
sample_manifest = pd.read_csv("/path/to/manifest.tsv", sep="\t")

cohort = build_cohort_matrices(
    sample_paths=[Path(p) for p in sample_manifest["beta_path"]],
    sample_ids=sample_manifest["sample_id"].tolist(),
    probe_reference=probe_reference,
    gene_panel=PRIMARY_GENES,
    output_dir=Path("/output/methylation/cohort"),
)

# Access matrices
gene_matrix = cohort["genes"]     # Genes × Samples
ccre_matrix = cohort["ccres"]     # cCREs × Samples
qc_summary = cohort["qc_summary"] # Sample QC metrics
```

### 4. Run Full Pipeline

```python
from methylation import run_methylation_pipeline

# Prepare sample manifest
manifest = pd.DataFrame({
    "sample_id": ["TCGA-A8-A09K", "TCGA-BH-A0B7", ...],
    "beta_path": ["/path/to/sample1.tsv", "/path/to/sample2.tsv", ...],
})

results = run_methylation_pipeline(
    probe_reference_path=PATHS.methylation_probe_reference,
    sample_manifest=manifest,
    genes=genes,
    ccres=ccres,
    gene_panel=PRIMARY_GENES,
    working_dir=PATHS.working_dir / "methylation",
    lncrnas=lncrnas,
    lncrna_panel=matched_lncrna_names,
    build_reference=True,
    build_per_sample=True,
    build_cohort=True,
)
```

### 5. Integrate with Element Table

```python
from methylation import integrate_methylation_with_element_table

# After processing a sample
element_table_with_meth = integrate_methylation_with_element_table(
    element_table=elem_focus,
    sample_ccre_meth=results["ccres"],
)

# Now element_table has columns:
# meth_ccre_beta_mean, meth_ccre_n_probes, etc.
```

---

## Output Directory Structure

```
methylation/
├── reference/
│   ├── probe_annotations.parquet    # Full annotated reference (use this)
│   └── probe_annotations.csv        # Simplified CSV for inspection
│
├── per_sample/
│   └── {sample_id}/
│       ├── {sample_id}_probes.parquet      # Per-probe data
│       ├── {sample_id}_gene_meth.csv       # Gene aggregation
│       ├── {sample_id}_lncrna_meth.csv     # lncRNA aggregation
│       ├── {sample_id}_ccre_meth.csv       # cCRE aggregation
│       └── {sample_id}_qc.json             # QC metrics
│
├── cohort/
│   ├── gene_meth_matrix.parquet     # Genes × Samples (promoter beta mean)
│   ├── gene_meth_matrix.csv
│   ├── lncrna_meth_matrix.csv       # lncRNAs × Samples
│   ├── ccre_meth_matrix.parquet     # cCREs × Samples
│   ├── ccre_meth_matrix.csv
│   └── sample_qc_summary.csv        # Per-sample QC
│
└── pipeline_summary.json            # Pipeline run metadata
```

---

## Integration with Other Data Types

### Comparison Table

| Data Type | Gene Level | cCRE Level | Per-Sample | Cohort Matrix |
|-----------|------------|------------|------------|---------------|
| **Methylation** | Promoter + body mean | cCRE mean | ✓ | ✓ |
| **Expression (RNA-seq)** | TPM/counts | - | ✓ | ✓ |
| **ATAC-seq** | Peak accessibility near TSS | Peak overlap | ✓ | ✓ |
| **Copy Number** | Gene dosage | - | ✓ | ✓ |
| **Structural Variants** | Gene/promoter hits | Enhancer hits | ✓ | - |
| **TADs** | Domain context | Domain context | Reference | Reference |
| **SCREEN** | - | Evidence scores | Reference | Reference |
| **ABC** | - | Enhancer scores | Reference | Reference |
| **HiChIP** | - | Loop connections | Reference | Reference |

### Typical Analysis Workflow

```python
# 1. Build all reference tables (once)
#    - Gene annotations
#    - cCRE annotations with SCREEN/ABC/HiChIP
#    - Probe reference with all overlaps

# 2. Per sample, for each data type:
#    - RNA-seq: gene expression
#    - Methylation: gene/cCRE methylation
#    - ATAC: peak accessibility
#    - CNV: gene dosage
#    - SV: gene/enhancer hits

# 3. Combine into gene-centric table:
gene_table = genes[["gene_name", "gene_id", "chrom", "start", "end"]]
gene_table = integrate_methylation_with_gene_table(gene_table, gene_meth)
gene_table = integrate_expression_with_gene_table(gene_table, expr_data)
# ... etc

# 4. Combine into cCRE-centric table:
ccre_table = elem_focus  # Already has SCREEN, ABC, HiChIP
ccre_table = integrate_methylation_with_element_table(ccre_table, ccre_meth)
# ... etc
```

---

## Biological Interpretation Notes

### Promoter vs Gene Body Methylation

- **Promoter hypermethylation** → Gene silencing
- **Gene body methylation** → Often correlates with *active* transcription
- These often have *inverse* relationships with expression

### CGI Context Matters

- **CGI probes**: More stable, often reflect long-term silencing
- **Shore probes**: More dynamic, often reflect recent regulatory changes
- Consider analyzing CGI and non-CGI probes separately

### Integration with APM Hypothesis

For the immune visibility (APM) analysis:

1. **Promoter methylation of APM genes** → Direct silencing mechanism
2. **Enhancer (cCRE) methylation** → May block TF binding (STAT1, IRF1)
3. **CGI context** → May indicate different silencing mechanisms
4. **Combined with SVs** → SV disrupts enhancer + methylation silences promoter = redundant silencing
