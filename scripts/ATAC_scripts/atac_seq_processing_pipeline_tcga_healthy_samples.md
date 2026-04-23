# ATAC-seq Processing Pipeline (TCGA + Healthy Samples)

This directory contains a modular pipeline for processing ATAC-seq data by harmonizing TCGA Pan-Cancer ATAC profiles with healthy (ENCODE) samples, followed by normalization and replicate collapsing.

The pipeline is designed to be:

- **Modular** – individual steps can be re-run independently  
- **Reproducible** – all paths and outputs are centrally configured  
- **Scalable** – supports multiple healthy samples and large cohorts

---

## Pipeline Overview

The pipeline consists of four main steps:

1. **Recount healthy ATAC BAMs** on the TCGA Pan-Cancer peak set  
2. **Harmonize** TCGA raw ATAC counts with healthy sample vectors  
3. **Normalize** ATAC counts using logCPM and quantile normalization  
4. **Collapse** TCGA technical replicates to case-level profiles

A master shell script orchestrates the workflow and supports checkpointing to avoid unnecessary recomputation.

---

## Scripts

### 1. `recount_atac.sh`

**Purpose**  
Recount a single healthy / ENCODE ATAC-seq BAM file using the TCGA Pan-Cancer ATAC peak set, producing a per-peak raw count vector compatible with TCGA ATAC matrices.

**Arguments**

- `TCGA_ATAC_PeakSet.bed`
- `ENCODE_sample.bam`
- `out_dir`

**Usage**

    recount_atac.sh <TCGA_ATAC_PeakSet.bed> <ENCODE_sample.bam> <out_dir>

**Inputs**

- Pan-Cancer ATAC peak set (BED)
- Healthy / ENCODE ATAC BAM file

**Outputs** (written to `out_dir`)

- `<bam_basename>.vector` — per-peak raw count vector
- `<bam_basename>.counts.txt` — BEDTools coverage output
- Intermediate files (sorted BAM, BAM index, genome file, sorted peak BED)

---

### 2. `harmonize_atac_counts.py`

**Purpose**  
Create a unified raw ATAC counts matrix by filtering TCGA samples to a selected cancer type and appending healthy ATAC vectors.

**Key operations**

- Filter TCGA samples using `_primary_disease`
- Map ATAC libraries via `bam_prefix` to `Case_ID`
- Append healthy `.vector` files

**Arguments (CLI)**

    python harmonize_atac_counts.py \
      --paths-json atac_paths.json \
      --first-sample-col 7 \
      --primary-disease "breast invasive carcinoma" \
      --tcga-cancer-type BRCA

**Inputs (from `atac_paths.json`)**

- `TCGA_RAW_ATAC_COUNTS_PATH`
- `TCGA_SAMPLES_TABLE_PATH`
- `ATAC_META_PATH`
- `HEALTHY_COUNTS_DIR`

**Outputs**

- `OUT_PATH` — combined raw ATAC counts matrix (TSV)

---

### 3. `normalize_atac_edger.R`

**Purpose**  
Normalize ATAC counts using edgeR logCPM and quantile normalization.

**Arguments**

- input matrix
- output matrix
- first sample column (0-based)

**Usage**

    normalize_atac_edger.R <input> <output> <first_sample_col>

**Outputs**

- LogCPM + quantile-normalized ATAC matrix

---

### 4. `collapse_atac_replicates.py`

**Purpose**  
Collapse TCGA technical replicates into case-level ATAC profiles.

**Rules**

- TCGA replicates averaged by `Case_ID`
- Non-TCGA samples unchanged

**Arguments (CLI)**

    python collapse_atac_replicates.py \
      --paths-json atac_paths.json \
      --first-sample-col 7

**Inputs**

- `ATAC_LOGCPM_QN_PATH`
- `ATAC_META_PATH`

**Outputs**

- `ATAC_CASE_LEVEL_PATH`

---

## Master Pipeline Script

### `run_atac_pipeline.sh`

**Purpose**  
Run the complete ATAC-seq processing pipeline with stage-level checkpointing.

**Stages**

- recount
- harmonize
- normalize
- collapse

**Usage**

    run_atac_pipeline.sh --peakset <TCGA_ATAC_PeakSet.bed> [options]

**Outputs**

- Healthy `.vector` files
- Raw combined ATAC matrix
- Normalized ATAC matrix
- Case-level ATAC matrix

---

## Configuration: `atac_paths.json`

Example keys:

- `TCGA_RAW_ATAC_COUNTS_PATH`
- `HEALTHY_COUNTS_DIR`
- `TCGA_SAMPLES_TABLE_PATH`
- `ATAC_META_PATH`
- `OUT_PATH`
- `ATAC_LOGCPM_QN_PATH`
- `ATAC_CASE_LEVEL_PATH`

---

## Notes

- File existence is used as a checkpointing mechanism.
- The pipeline is designed for easy extension to additional cancer types.
