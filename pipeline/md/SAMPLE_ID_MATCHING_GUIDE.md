# Sample ID Matching Guide (TCGA BRCA APM pipeline)

This file explains **how sample IDs appear in each module/input**, how they map onto the pipeline’s
**normalized join keys**, and the recommended join strategy across modules.

---

## Normalized TCGA join keys (canonical)

The pipeline normalizes any TCGA-like identifier into these keys (see `pipeline/sample_ids.py`):

- **`participant`**: `TCGA-XX-YYYY`
- **`sample`**: `TCGA-XX-YYYY-SS` where \(SS\) is the 2-digit sample type (e.g. `01`, `10`, `11`)
- **`sample_vial`**: `TCGA-XX-YYYY-SSA` (adds vial letter, e.g. `01A`, `10A`, `11A`)
- **`aliquot`**: full aliquot barcode when present (extra hyphen groups)

### Recommended canonical keys

- **Cross-module “omics” joins (SNV/SV/CNV/Methylation/RPPA/HLA/ATAC)**: join on **`sample_vial`**
  - Rationale: these are assay-specific outputs and naturally track `...-01A`/`...-10A`/`...-11A`.
- **Case-level joins (RNA expression + immune tables)**: join on **`sample`**
  - Rationale: most immune tables and RNA matrices are indexed at `...-01` granularity.
- **Participant-only data (HiChIP columns in your processed file)**: join on **`participant`**

---

## Where the ID comes from (per module / input)

### SNV (Mutect2 VCFs → VEP VCFs)

- **Input VCF “sample labels” inside the VCF**: typically **`TUMOR` / `NORMAL`**
  - These are **not TCGA IDs**.
- **Canonical TCGA mapping source**: `annotations/SNV/samples.tsv`
  - Column: **`Sample ID`** contains a tumor+normal pair like: `TCGA-...-01A, TCGA-...-10A`
  - The pipeline picks the tumor entry based on `Tumor Descriptor`.
- **Output identity**:
  - Per-sample output directories are named by **`sample_vial`** (e.g. `TCGA-XX-YYYY-01A`)
- **Normalization behavior**:
  - `save_snv_outputs()` adds normalized columns to saved SNV tables:
    - Prefer `aliquot_id` column if present, else `sample_id` column, else infer from the per-sample output folder name.

### SV (Manta SV VCFs)

- **Input VCF internal labels**: often `TUMOR` / `NORMAL` (not joinable)
- **Canonical TCGA mapping source**: `annotations/SV/samples.tsv`
  - `Sample ID` provides **`sample_vial`**.
- **Output identity**:
  - Saved files named like `{sample_vial}_strict_sv_set.csv` (or lenient)
- **Normalization behavior**:
  - `pipeline/SV/pipeline.py` injects normalized columns using `raw_id=sample_id` before saving.

### CNV (ASCAT segments)

- **Canonical TCGA mapping source**: `annotations/CNV/samples.tsv`
  - `pipeline/CNV/loader.py::extract_sample_id_from_annotations()` returns a tumor **`sample_vial`**
    using `Tumor Descriptor` logic.
- **Output identity**:
  - Saved files named like `{sample_vial}_cnv_annotated.csv`
- **Normalization behavior**:
  - `pipeline/CNV/runner.py` injects normalized columns using `raw_id=sid` before saving.

### Methylation

- **Canonical sample id**:
  - Your methylation per-sample processing already uses **`sample_vial`** (`TCGA-...-01A`) as `sample_id`
    and/or output folder naming.
- **Output identity**:
  - Per-sample tables: `{sample_vial}_probes.csv`, `{sample_vial}_gene_meth.csv`, `{sample_vial}_ccre_meth.csv`, etc.
- **Normalization behavior**:
  - `pipeline/Methylation/methylation_table.py` injects normalized columns using `raw_id=sample_id`
    before writing each per-sample table.

### RPPA

- **Metadata source**: `annotations/rppa/samples.tsv` (TCGA `Sample ID` column is `sample_vial`)
- **Matrix identity**:
  - RPPA expression matrices typically use `sample_vial` as the index (and sometimes also have a numeric lab id upstream).
- **Normalization behavior**:
  - `pipeline/rppa/rppa_main.py` injects normalized columns derived from the **index** before saving
    `rppa_analysis_combined.csv/parquet`.

---

## External inputs (RNA / immune / ATAC / HLA / HiCHIP / miRNA)

These are *not* produced by a pipeline module, but they are key join sources.

### RNA expression (wide matrix)

- **Paths (config)**: `PATHS.rna_expression` / `PATHS.rna_expression_raw` (wide TSV; TCGA tokens in **column headers**).
- **Where the ID is**: **column headers** (not row labels), optionally with the first column used as gene / feature index.
- **Observed shape**: often `sample_vial`-level (e.g. `TCGA-...-01A`) in this file
- **Column-derived manifest (recommended for tooling)**:
  - `PATHS.rna_samples_tsv` → `annotations/RNA/samples.tsv`
  - Regenerate with `python scripts/annotations/build_annotation_sample_manifests.py` after the matrix changes.
  - Each row mirrors GDC-style fields; `Sample ID` is the parsed tumor barcode (`sample_vial` when available).
- **Normalization strategy**:
  - Do **not** reshape the full matrix.
  - Optional companion table: `annotations/_normalized/RNA_expression.sample_metadata.tsv`
    - One row per sample column header
    - Contains normalized keys (`participant`, `sample`, `sample_vial`, `aliquot`)

### miRNA expression (wide matrix)

- **Path (config)**: `PATHS.mirna_expression_tsv`
  - Default: `data/miRNA/XENA_mirna_arm_specific.tsv`
- **Format**: tab-separated TSV
- **Units**: log2(RPM+1) in your current file
- **Where the ID is**: **column headers** (samples)
- **Observed ID shape**: `sample_vial`-level (e.g. `TCGA-...-01A`) in this file
- **Column-derived manifest**:
  - `PATHS.mirna_samples_tsv` → `annotations/miRNA/samples.tsv`
  - Regenerate with `python scripts/annotations/build_annotation_sample_manifests.py`

### Immune subtype annotations (tidy tables)

1) **Advanced immune subtypes**
- **Path (config)**: `PATHS.immune_subtype_annotations`
  - `annotations/BRCA_immune_subtypes_advanced.tsv`
- **ID column**: `TCGA Participant Barcode` (typically `sample` like `TCGA-...-01`)
- **Normalized output written**:
  - `annotations/_normalized/BRCA_immune_subtypes_advanced.normalized.tsv`

2) **Thornsson/Thorsson immune table**
- **Path (config)**: `PATHS.thornsson_immune_table`
  - `annotations/Thornsson_immune_table.tsv`
- **ID column**: `sample_id` (typically `sample` like `TCGA-...-01`)
- **Normalized output written**:
  - `annotations/_normalized/Thornsson_immune_table.normalized.tsv`

### ATAC case-level matrix (wide)

- **Path (config)**: `PATHS.atac_case_level_matrix`
  - Default: `data/TCGA_ATAC/TCGA_NORMAL_LOG_CPM_QN_BRCA_case_level.csv`
- **Format**: comma-separated CSV; TCGA sample columns typically start after the leading peak / coordinate columns (e.g. column 11+).
- **Where the ID is**: **column headers**
- **Observed ID shape**: aliquot-like headers (e.g. `TCGA-...-01A-31-A615-42`)
  - Normalization derives:
    - `sample_vial` = `TCGA-...-01A`
    - `sample` = `TCGA-...-01`
    - `participant` = `TCGA-...`
    - `aliquot` = the full header string
- **Column-derived manifest**:
  - `PATHS.atac_samples_tsv` → `annotations/ATAC/samples.tsv` (from matrix headers; same generator script as RNA).
- **Normalized companion table written**:
  - `annotations/_normalized/ATAC_case_level.sample_metadata.tsv`

### HLA types (tidy)

- **Path (config)**: `PATHS.hla_types_tsv` → `annotations/HLA_types.tsv`
- **ID column**: `aliquot_id` (full aliquot barcode)
- **Delimiter note**: this file is comma-separated in this repo; the normalizer auto-detects.
- **Aliquot-derived manifest**:
  - `PATHS.hla_samples_tsv` → `annotations/HLA/samples.tsv` (one row per **unique** `aliquot_id`; `Sample ID` = parsed `sample_vial` when possible).
  - Regenerate with `python scripts/annotations/build_annotation_sample_manifests.py`.
- **Normalized output written**:
  - `annotations/_normalized/HLA_types.normalized.tsv`

### HiCHIP processed matrix (wide, TCGA participants)

- **Path (config)**: `PATHS.hichip_tcga_processed_csv`
  - Default: `data/HiCHIP/TCGA_BRCA_PROCESSED.csv`
- **Format**: comma-separated CSV
- **Where the ID is**: **column headers**, starting at the **7th column** (first 6 columns are coordinates)
- **Observed ID shape**: `participant` (e.g. `TCGA-A7-A0CH`)
- **Column-derived manifest**:
  - `PATHS.hichip_samples_tsv` → `annotations/HiCHIP/samples.tsv`
  - Regenerate with `python scripts/annotations/build_annotation_sample_manifests.py`
  - Note: `Sample ID` in this manifest is **participant** (no sample type/vial in the source headers).

---

## Practical join recipes (what to do in analysis code)

### Join RNA (matrix) to module outputs

1) Load module output (SNV/SV/CNV/Methylation/RPPA/HLA/ATAC-derived tables) → has `sample_vial`.
2) Load `RNA_expression.sample_metadata.tsv`.
3) Join on `sample_vial` (if RNA headers are `...-01A`) **or** join on `sample` (if RNA headers are `...-01`).

### Join immune tables to module outputs

Immune tables are mostly `sample`-level (`...-01`), while module outputs are `sample_vial` (`...-01A`).

Recommended:
- Reduce module tables to case-level by using `sample` (drop vial letter), then join to immune on `sample`.

### Join HLA to anything else

HLA provides `aliquot_id` → derive `sample_vial` and join on `sample_vial`.

### Join ATAC to anything else

ATAC headers look aliquot-like → derive `sample_vial` from the header and join on `sample_vial`.

---

## Cell line / biosample label normalization

ENCODE and internal file layouts often disagree on spelling (`MCF-7` vs `MCF7`, tissue names with
spaces vs snake_case). The pipeline keeps **two layers**:

1. **On-disk / upstream identity** (unchanged strings you must match when reading inputs):
   - `BiosampleConfig.ccre_signal_cell_lines`: subdirectory names under `PATHS.cell_lines_dir`.
   - `BiosampleConfig.hichip_panel`: subdirectory names under `PATHS.hichip_dir`.
   - `BiosampleConfig.abc_celltypes`: exact `CellType` values in the ABC predictions TSV (filter only).

2. **Canonical keys in outputs** (nested dicts, `elem_focus` columns, gene_links):
   - ChIP: `normalize_cell_line_label` on `cell_type` (whitelist `chip_brca_celltypes` is canonical).
   - cCRE nested signal columns: `canonical_ccre_signal_column_name(disk_subdir)` so a folder
     named `MCF-7` would still load from disk but store under column `MCF7` (see `ccre_loader.add_multiple_cell_line_signals`).
   - SCREEN: raw `Biosample` from the file → `canonicalize_screen_biosample`; `BiosampleConfig.screen_exp` /
     `screen_comp` list those **canonical** keys (order = column order in collapsed tables).
   - ABC: after filtering on upstream `CellType`, rows are relabeled with `canonicalize_abc_cell_type`
     before aggregation (e.g. `MCF-7-ENCODE` → `MCF7_ENCODE`).
   - HiChIP: loops are read from `hichip_panel` disk names; `elem_focus` / `gene_links["hichip"]` keys use
     `canonical_hichip_output_key` (usually identical to the folder name for lines like `MCF7`).

All helpers live in `pipeline/biosample_names.py`.

### PAM50 / subtype labels (tumor vs healthy)

For **biology-aware grouping** (LumA, LumB, HER2, Basal, plus non-tumor buckets), use
`pipeline/cell_line_subtype_map.py`:

- ``subtype_for_key(canonical_biosample_or_line)`` → ``SubtypeRecord`` with ``pam50_group``
  in ``{LumA, LumB, HER2, Basal, Normal_like, Healthy_tissue, Unknown}``.
- Tumor cell-line annotations are aligned with ``../TADs/config.py`` ``CELL_LINE_METADATA``
  (TADs uses ``Her2`` spelling; this module exposes ``HER2`` for the four-class label).
- ``Healthy_tissue`` covers primary / organoid-style ENCODE biosample keys (e.g. SCREEN
  ``breast_mammary_tissue``); ``Normal_like`` covers benign cultures (HMEC, MCF10A, …).
- ``PAM50_GROUP_TO_CELL_LINES`` mirrors TADs ``PAM50_TO_CELL_LINES`` with the same spelling
  normalization and a few APM-only lines (e.g. ``MDA-MB-468``, ``HMEC1`` / ``HMEC2``).

---

## Generated normalized external files

Run:

```bash
python scripts/annotations/normalize_external_annotations.py
```

Outputs:
- `annotations/_normalized/BRCA_immune_subtypes_advanced.normalized.tsv`
- `annotations/_normalized/Thornsson_immune_table.normalized.tsv`
- `annotations/_normalized/HLA_types.normalized.tsv`
- `annotations/_normalized/RNA_expression.sample_metadata.tsv`
- `annotations/_normalized/ATAC_case_level.sample_metadata.tsv`

