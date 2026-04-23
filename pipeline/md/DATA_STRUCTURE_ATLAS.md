# Pipeline Data Structure Atlas

Visual reference for every output table and its nested column structures.

**Related**:
- [Thresholds and weak/strong signal rules](THRESHOLDS_AND_SIGNAL_STRENGTH.md) — defaults in `ThresholdConfig` (`pipeline/config.py`) and how SCREEN, ABC, methylation, and filters use them.
- [APM session log (science + paths + commands)](APM_SESSION_LOG.md) — human-readable changelog for long working sessions.
- [Gene symbol normalization](GENE_SYMBOL_NORMALIZATION.md) — how HGNC/UCSC/legacy/Ensembl aliases fold into canonical panel symbols.

---

## Table of Contents

1. [Element Focus Table (Main Pipeline Output)](#1-element-focus-table)
2. [ATAC Peak Table](#2-atac-peak-table)
3. [SNV Table](#3-snv-table)
4. [SV Table](#4-sv-table)
5. [CNV Segment Table](#5-cnv-segment-table)
6. [Methylation Tables](#6-methylation-tables)
7. [RPPA Tables](#7-rppa-tables)
8. [Gene / lncRNA Tables](#8-gene--lncrna-tables)
9. [RNA Expression Matrix (log2(TPM+1))](#9-rna-expression-matrix-log2tpm1)
10. [ATAC Case-Level Sample Matrix (`TCGA_NORMAL_LOG_CPM_QN_BRCA_case_level`)](#10-atac-case-level-sample-matrix-tcga_normal_log_cpm_qn_brca_case_level)
11. [Unified Clinical + Immune Table (BRCA)](#11-unified-clinical--immune-table-brca)
12. [miRTarBase Output Tables](#12-mirtarbase-output-tables)
13. [HiChIP TCGA BRCA processed matrix (external)](#hichip-tcga-brca-processed)
14. [miRNA Xena arm-specific expression (external)](#mirna-xena-arm-specific)
15. [Cohort covariates table (post-processing)](#15-cohort-covariates-table-post-processing)
16. [External annotation inputs (from `config.PATHS`)](#16-external-annotation-inputs-from-configpaths)

---

## 1. Element Focus Table

**Long-format scanning tables (deferred)**: The pipeline emits **flat** `scan_*` columns on element-focus and ATAC peak outputs for scanner/query performance. Long tables at the full \((\text{element} \times \text{gene} \times \text{biosample} \times \text{assay})\) grain are **not** written by default; introduce them only when workloads need joins that flat columns cannot serve.

**Primary file**: `regulatory_element_focus_with_evidence.parquet`
(STEP 12 default; written in row chunks via `pyarrow.ParquetWriter`, `zstd`-compressed.
Optional legacy CSV is written **only** when `APM_WRITE_ELEM_FOCUS_CSV=1`.)

**Load with**:

```python
from pipeline.regulatory_elements import load_regulatory_element_focus
elem = load_regulatory_element_focus()  # decodes JSON-encoded nested columns
```

Nested object columns (`gene_links`, `TAD_domains`, `TAD_boundary_overlaps`,
`atac_peaks`, `chip_hits`, `ABC_enhancers`, `hichip`, `screen_exp`,
`screen_comp`) are **JSON-encoded strings on disk** to keep PyArrow /
write RAM bounded; the loader decodes them back into `dict` / `list`.

### Scanning columns (derived, scalar)

For the combinatorial scanning engine, the pipeline also materializes **pushdown-friendly scalar columns** derived from the nested JSON bundles (env `APM_WRITE_SCANNING_COLUMNS=0` disables).

These columns are prefixed `scan_` and are intended for fast filtering without decoding nested JSON:

- **Gene links**
  - `scan_gene_links_n_genes`
  - `scan_gene_links_any_panel_gene`
  - `scan_gene_links_any_tier1_lncrna`
- **SCREEN (row-level)**
  - `scan_screen_exp_n_biosamples_any`, `scan_screen_exp_n_biosamples_strong`
  - `scan_screen_exp_MCF7_max_score`, `scan_screen_exp_MCF10A_max_score`
  - `scan_screen_exp_breast_n_strong`
- **ABC**
  - `scan_ABC_max_score_any_celltype`, `scan_ABC_max_score_MCF7_ENCODE`
  - `scan_ABC_n_celltypes_strong`
  - `scan_ABC_self_promoter_any`
- **HiChIP**
  - `scan_hichip_n_celltypes_with_loops`
  - `scan_hichip_max_score_any`
  - `scan_hichip_n_loops_MCF7`
- **ChIP**
  - `scan_chip_hits_n_TFs`
  - `scan_chip_hits_has_STAT1`, `scan_chip_hits_has_CTCF`
  - `scan_chip_hits_MCF7_n_TFs`

### `scan_*` column contract (derivation + threshold source)

All definitions are implemented in `pipeline/scanning_columns.py` (`derive_elem_focus_scanning_columns` / `derive_atac_peak_scanning_columns`). Runtime contract tests: `tests/test_scanning_columns_spec.py`.

| Column | Derivation rule | Threshold / symbol source |
|--------|-----------------|---------------------------|
| `scan_gene_links_n_genes` | `len(gene_links)` when `gene_links` is a `dict` of gene → evidence bundle; else `0` | Structural (no threshold) |
| `scan_gene_links_any_panel_gene` | `any(gene in PIPELINE_GENE_PANEL for gene in gene_links.keys())` | `pipeline.config.PIPELINE_GENE_PANEL` (env `APM_USE_EXTENDED_GENE_PANEL` switches legacy 66 vs extended) |
| `scan_gene_links_any_tier1_lncrna` | `any(gene in TIER1_LNCRNA_GENES for gene in gene_links.keys())` | `pipeline.config.TIER1_LNCRNA_GENES` |
| `scan_screen_exp_n_biosamples_any` | Count biosamples under `screen_exp["per_biosample"]` with any assay `strength != "none"` | Strength strings from SCREEN linker (`pipeline/evidence/screen_links.py`); quantile thresholds from `ThresholdConfig` / per-run CSV under `working_dir/evidence/` |
| `scan_screen_exp_n_biosamples_strong` | Count biosamples with any assay `strength == "strong"` | Same as above |
| `scan_screen_exp_{MCF7,MCF10A}_max_score` | Max `score` across assays for that biosample in `screen_exp` | Raw score only; strength is separate |
| `scan_screen_exp_breast_n_strong` | Sum of `n_strong` across assays in `screen_exp["conservation_breast"]` | Nested counter from SCREEN export |
| `scan_ABC_max_score_any_celltype` | Max `ABC_score` across `ABC_enhancers[]` entries / cell types | `is_strong` on each cell type uses ABC thresholds in `ThresholdConfig` (see `THRESHOLDS_AND_SIGNAL_STRENGTH.md`) |
| `scan_ABC_n_celltypes_strong` | Number of cell types with `is_strong == True` | ABC strong rule (config) |
| `scan_ABC_self_promoter_any` | Any enhancer with `is_self_promoter == True` | Boolean from ABC record |
| `scan_ABC_max_score_MCF7_ENCODE` | Max score restricted to `MCF7_ENCODE` cell type key | Same as ABC score |
| `scan_hichip_n_celltypes_with_loops` | Count `hichip` dict entries with `n_loops > 0` | HiChIP loop presence is count-based (no universal score threshold in this column) |
| `scan_hichip_max_score_any` | Max `max_score` across cell types in `hichip` | Per-cell-type upstream processing |
| `scan_hichip_n_loops_MCF7` | `hichip["MCF7"]["n_loops"]` (or `0`) | Named biosample key |
| `scan_chip_hits_n_TFs` | `len(chip_hits)` top-level TF keys | Structural |
| `scan_chip_hits_has_STAT1` / `has_CTCF` | TF key present in `chip_hits` | Structural |
| `scan_chip_hits_MCF7_n_TFs` | Count TF keys that have a `MCF7` cell-type subtree | Named cell line key |

**Nested JSON “same shape” caveat**: In this atlas, “same shape” applies to **closely related spatial hit families** (e.g. SV `gene_hits` vs `lncRNA_hits` built by the same mapper). It does **not** mean SNV VEP `gene_hits` (CSQ transcript dicts) share keys with SV spatial `gene_hits`. When in doubt, treat schemas as **table- and modality-specific**; SV span-level key stability is guarded by `tests/test_sv_span_gene_hit_schema_keys.py`.

**Legacy file**: `regulatory_element_focus_with_evidence.csv` (same schema;
nested columns are JSON strings that the loader also decodes).

**Granularity**: One row per cCRE (candidate cis-Regulatory Element)

**Units**:
- Genomic coordinates (`start`, `end`): **bp**
- Distance windows / distances: **bp**
- Fractions (overlap fractions, normalized positions): **0–1**

### Flat columns

```
┌─────────────────────────────────────────────────────────────────────┐
│  FLAT COLUMNS (one value per row)                                   │
├────────────────┬──────────┬─────────────────────────────────────────┤
│ Column         │ Type     │ Example                                 │
├────────────────┼──────────┼─────────────────────────────────────────┤
│ cCRE_id        │ str      │ "EH38E1594103.1"                        │
│ ENCODE_id      │ str      │ "EH38E1594103"                          │
│ type           │ str      │ "pELS"                                  │
│ raw_type       │ str      │ "pELS,CTCF-bound"                       │
│ chrom          │ str      │ "chr6"                                  │
│ start          │ int      │ 29942000                                │
│ end            │ int      │ 29942800                                │
│ 0–100kb        │ str      │ "HLA-A,HLA-G"  (comma-sep gene names)  │
│ 100–250kb      │ str      │ "HLA-B"                                 │
│ 250–500kb      │ str      │ ""                                      │
│ 500–1000kb     │ str      │ ""                                      │
│ min_dist_to_   │ int      │ 5230                                    │
│   any_gene     │          │                                         │
│ genes_by_      │ str      │ "HLA-A:5230,HLA-G:89420"               │
│   exact_dist   │          │                                         │
└────────────────┴──────────┴─────────────────────────────────────────┘
```

### Cell-line signal columns (one per cell line)

```
Column name = cell line (e.g. "MCF7", "HMEC1", "breast_tissue")
Each cell is a dict:

MCF7 ─────────────────────────────────────────────
│
├── "in_MCF7"     : bool        ← is this cCRE active in MCF7?
├── "H3K27ac"     : float|None  ← signal value
├── "H3K4me3"     : float|None
├── "CTCF"        : float|None
└── "DNase"       : float|None
```

### `chip_hits` column (ChIP-seq overlaps per cCRE)

**Source**: `PATHS.chip_unified` (cached unified ChIP peak table), built by `pipeline/CHIP/chip_loader.py`.
There is **no** cell-line whitelist here: every peak in the parquet that overlaps the cCRE is counted.
(Optional **SV-only** ChIP mapping in `pipeline/CHIP/sv_chip_intersect.map_svs_to_chip` *does* support
`cell_type_whitelist` / `tf_whitelist` when building SV `chip_hits`.)
The unified table includes a `sample_id` (derived from each BED filename stem) so multiple
samples/replicates for the same TF×cell line stay distinguishable (e.g. `CTCF_MCF7_1`, `CTCF_MCF7_2`).

Each cell is a nested dict (TF hierarchy) aggregated over overlapping ChIP peaks:
`tf → cell_type → source → summary`:

**Cell line keys (`cell_type`)** use canonical labels (e.g. `MCF7`, `T47D`); raw ChIP tables may say `MCF-7` — see `pipeline/biosample_names.normalize_cell_line_label`. SCREEN `per_biosample` keys follow `canonicalize_screen_biosample` (see `BiosampleConfig.screen_*`); ABC `ABC_full` keys follow `canonicalize_abc_cell_type`.

```
chip_hits ─────────────────────────────────────────────────────────────
│
├── "CTCF" ────────────────────────────────────────────────────────────
│   ├── "MCF7" ────────────────────────────────────────────────────────
│   │   ├── "ENCODE" ──────────────────────────────────────────────────
│   │   │   ├── "n_peaks"        : int
│   │   │   ├── "n_samples"      : int              ← distinct `sample_id`s with ≥1 hit
│   │   │   ├── "sample_ids"     : list[str]        ← e.g. ["CTCF_MCF7_1","CTCF_MCF7_2"]
│   │   │   ├── "max_overlap_bp" : int
│   │   │   ├── "max_score_norm" : float|None  (ENCODE signalValue; CHIP_ATLAS is None)
│   │   │   └── "mean_score_norm": float|None
│   │   │
│   │   │   └── "per_sample" ──────────────────────────────────────────
│   │   │       │ dict keyed by `sample_id` with the same summary keys:
│   │   │       ├── "CTCF_MCF7_1": {n_peaks, max_overlap_bp, max_score_norm, mean_score_norm}
│   │   │       └── "CTCF_MCF7_2": { ... }
│   │   └── "CHIP_ATLAS" ── (same summary keys)
│   └── "T47D" ── (same shape)
│
├── "STAT1" ── (same shape)
└── ... per TF observed overlapping this cCRE
```

**Where `chip_hits` is integrated**

| Surface | Role |
|--------|------|
| **Element focus table** (`regulatory_element_focus_with_evidence.*`) | Column `chip_hits`: TF → cell line → source → overlap summaries per cCRE. Written in `pipeline/main.py` via `integrate_chip_hits_to_element_table` (`pipeline/CHIP/chip_hits.py`). |
| **SV processed tables** | Column `chip_hits`: list of overlap dicts per SV / breakpoint (see §4). Built in `pipeline/CHIP/sv_chip_intersect.py` during the SV pipeline. |

**Where `chip_hits` is not used**

CNV segment outputs, methylation probe/gene tables, RPPA tables, RNA expression matrices, and HiChIP / miRNA standalone exports do **not** attach ChIP peak bundles as `chip_hits`.

**SNV:** uses separate columns ``snv_chip_hits`` / ``snv_chip_aggregate`` (strict overlap with ``PATHS.chip_unified``; not the cCRE nested ``chip_hits`` shape).

### `TAD_domains` column (nested dict, per biosample)

```
TAD_domains ──────────────────────────────────────────────────────────
│
├── "Kim_T47D" ───────────────────────────────────────────────────────
│   │
│   ├── "domains" ────────────────────────────────────────────────────
│   │   │   dict of { domain_id → spatial_relation }
│   │   │
│   │   ├── "Kim_T47D_domain_001"  : "contains"
│   │   ├── "Kim_T47D_domain_002"  : "overlap_right"
│   │   └── ...
│   │
│   └── "primary" ────────────────────────────────────────────────────
│       │   the single best-matching domain
│       │
│       ├── "domain_id"   : "Kim_T47D_domain_001"
│       ├── "rel"         : "contains"
│       │
│       ├── "domain" ─────────────────────────────────────────────────
│       │   ├── "chrom"   : "chr6"
│       │   ├── "start"   : 1000000
│       │   ├── "end"     : 2000000
│       │   └── "len"     : 1000000
│       │
│       ├── "feature" ────────────────────────────────────────────────
│       │   ├── "chrom"   : "chr6"
│       │   ├── "start"   : 29942000
│       │   ├── "end"     : 29942800
│       │   └── "mid"     : 29942400
│       │
│       ├── "boundaries" ─────────────────────────────────────────────
│       │   ├── "left" ───────────────────────────────────────────────
│       │   │   ├── "boundary_id" : "Kim_T47D_b_001"
│       │   │   ├── "dist_bp"     : 500000
│       │   │   └── "overlap"     : false
│       │   │
│       │   ├── "right" ──────────────────────────────────────────────
│       │   │   ├── "boundary_id" : "Kim_T47D_b_002"
│       │   │   ├── "dist_bp"     : 498000
│       │   │   └── "overlap"     : false
│       │   │
│       │   └── "nearest" ────────────────────────────────────────────
│       │       ├── "side"    : "right"
│       │       └── "dist_bp" : 498000
│       │
│       └── "normalized" ─────────────────────────────────────────────
│           ├── "frac_from_left"             : 0.501
│           ├── "frac_from_right"            : 0.499
│           └── "frac_to_nearest_boundary"   : 0.499
│
├── "Kim_HCC70"    ─── (same structure as above)
├── "Rao_HMEC"     ─── (same structure)
└── ... up to 24 biosamples
```

### `gene_links` column (the main multi-evidence bundle)

```
gene_links ───────────────────────────────────────────────────────────
│
│   dict keyed by gene_name
│
├── "HLA-A" ──────────────────────────────────────────────────────────
│   │
│   ├── "gene_id"    : "ENSG00000206503"
│   ├── "gene_type"  : "protein_coding"
│   │
│   ├── "screen_exp" ─── [SCREEN experimental evidence] ──────────────
│   │   │
│   │   ├── "per_biosample" ──────────────────────────────────────────
│   │   │   │
│   │   │   ├── "MCF7" ──────────────────────────────────────────────
│   │   │   │   │
│   │   │   │   ├── "Intact-HiC" ────────────────────────────────────
│   │   │   │   │   ├── "score"    : 12.45
│   │   │   │   │   ├── "p_value"  : 1.2e-8
│   │   │   │   │   └── "strength" : "strong" | "weak" | "none"
│   │   │   │   │
│   │   │   │   ├── "RNAPII-ChIAPET" ── (same shape) ────────────────
│   │   │   │   └── ...per assay type found
│   │   │   │
│   │   │   ├── "MCF10A" ── (same shape per assay)
│   │   │   └── ...per biosample in panel
│   │   │
│   │   ├── "conservation_global" ────────────────────────────────────
│   │   │   │   per assay type:
│   │   │   ├── "Intact-HiC" ────────────────────────────────────────
│   │   │   │   ├── "n_biosamples" : 3
│   │   │   │   ├── "n_strong"     : 2
│   │   │   │   └── "n_weak"       : 1
│   │   │   └── ...
│   │   │
│   │   └── "conservation_breast" ── (same shape as conservation_global)
│   │
│   ├── "screen_comp" ─── [SCREEN computational evidence] ────────────
│   │   │
│   │   │   (IDENTICAL SHAPE to screen_exp, with different
│   │   │    biosamples and assay types)
│   │   │
│   │   ├── "per_biosample"
│   │   │   ├── "MCF7"
│   │   │   │   ├── "ABC" ── {score, p_value, strength}
│   │   │   │   ├── "EPI"
│   │   │   │   └── ...
│   │   │   └── ...
│   │   ├── "conservation_global"
│   │   └── "conservation_breast"
│   │
│   │   (Provenance) The per-run SCREEN score quantiles used to assign `strength`
│   │   are persisted by the main pipeline under `working_dir/evidence/`:
│   │     - `screen_exp_score_quantiles.csv` (+ `.meta.json`)
│   │     - `screen_comp_score_quantiles.csv` (+ `.meta.json`)
│   │
│   ├── "ABC_enhancers" ─── [Activity-by-Contact predictions] ────────
│   │   │
│   │   │   list of enhancer entries mapped to this cCRE
│   │   │
│   │   ├── [0] ─────────────────────────────────────────────────────
│   │   │   ├── "start"    : 29940000
│   │   │   ├── "end"      : 29942500
│   │   │   │
│   │   │   └── "ABC_full" ─── dict keyed by celltype ───────────────
│   │   │       │
│   │   │       ├── "MCF7_ENCODE" ────────────────────────────────────
│   │   │       │   ├── "ABC_score"             : 0.035
│   │   │       │   ├── "ABC_num"               : 0.042
│   │   │       │   ├── "activity"              : 15.6
│   │   │       │   ├── "distance"              : 5230
│   │   │       │   ├── "element_class"         : "intergenic"
│   │   │       │   ├── "is_self_promoter"      : false
│   │   │       │   ├── "hic_pl_scaled"         : 0.8
│   │   │       │   ├── "powerlaw_score"        : 0.02
│   │   │       │   ├── "gene_expr"             : 8.2
│   │   │       │   ├── "promoter_activity_q"   : 0.75
│   │   │       │   ├── "gene_is_expressed"     : true
│   │   │       │   ├── "rank_within_gene"      : 0.85
│   │   │       │   ├── "is_present"            : true
│   │   │       │   └── "is_strong"             : false
│   │   │       │
│   │   │       ├── "MCF10A-Ji2017" ── (same keys)
│   │   │       └── ...
│   │   │
│   │   └── [1] ── (another enhancer entry, same shape)
│   │
│   └── "hichip" ─── [H3K27ac HiChIP loop evidence] ────────────────
│       │
│       │   dict keyed by cell type
│       │
│       ├── "MCF7" ───────────────────────────────────────────────────
│       │   ├── "n_loops"    : 3
│       │   ├── "max_counts" : 42.0
│       │   ├── "max_score"  : 8.5
│       │   │
│       │   └── "loops" ── list ──────────────────────────────────────
│       │       ├── [0] ──────────────────────────────────────────────
│       │       │   ├── "loop_id"      : 142
│       │       │   ├── "counts"       : 42.0
│       │       │   ├── "score"        : 8.5
│       │       │   ├── "n_reps"       : 2
│       │       │   ├── "anchor_ccre"  : "A"
│       │       │   ├── "anchor_prom"  : "B"
│       │       │   └── ...
│       │       └── [1] ── ...
│       │
│       ├── "HMEC" ── (same shape)
│       └── ...
│
├── "HLA-G" ── (same full structure as "HLA-A" above)
└── ...per gene linked to this cCRE
```

---

## External sample-ID normalization (recommended join keys)

Across modules and annotation tables, append these columns whenever a table has *one row per sample*:

- `participant`: `TCGA-XX-YYYY`
- `sample`: `TCGA-XX-YYYY-SS` (e.g. `...-01`)
- `sample_vial`: `TCGA-XX-YYYY-SSA` (e.g. `...-01A`)
- `aliquot`: full aliquot barcode when available (`TCGA-...-01A-..-..-..`)

For wide matrices where samples are **column headers** (RNA / ATAC), do **not** reshape the matrix.
Instead, generate a companion table `*.sample_metadata.tsv` with one row per sample column and the same
normalized keys (script: `scripts/annotations/normalize_external_annotations.py`).

## 2. ATAC Peak Table

**File**: `atac_peaks_annotated.csv` / `.parquet`
**Granularity**: One row per ATAC peak

**Units**:
- Genomic coordinates (`start`, `end`, `center`, `length`): **bp**
- Distances to TSS / windows: **bp**
- `percentGC`: fraction **0–1**

### Flat columns

```
┌─────────────────────────┬──────────┬──────────────────────────────┐
│ Column                  │ Type     │ Example                      │
├─────────────────────────┼──────────┼──────────────────────────────┤
│ peak_id                 │ str      │ "chr6:29940000-29941200"     │
│ chrom                   │ str      │ "chr6"                       │
│ start                   │ int      │ 29940000                     │
│ end                     │ int      │ 29941200                     │
│ center                  │ int      │ 29940600                     │
│ length                  │ int      │ 1200                         │
│ original_name           │ str      │ "Peak_42381"                 │
│ score                   │ float    │ 7.82                         │
│ annotation              │ str      │ "Intergenic"                 │
│ percentGC               │ float    │ 0.52                         │
│ n_genes_total           │ int      │ 4                            │
│ n_genes_overlapping     │ int      │ 1                            │
│ n_lncrnas_total         │ int      │ 2                            │
│ n_lncrnas_overlapping   │ int      │ 0                            │
│ n_ccres_total           │ int      │ 3                            │
│ n_ccres_overlapping     │ int      │ 2                            │
│ linked_genes            │ list     │ ["HLA-A","HLA-G","MICB"]     │
│ linked_lncrnas          │ list     │ ["LINC01149","HCG18"]        │
└─────────────────────────┴──────────┴──────────────────────────────┘
```

### Scanning columns (derived, scalar)

ATAC outputs also include a small set of `scan_` columns derived from nested objects to support fast filtering:

- `scan_gene_links_n_genes`, `scan_gene_links_any_panel_gene`
- `scan_lncrna_links_n_lncrnas`, `scan_lncrna_links_any_tier1_lncrna`
- `scan_ccre_links_n_ccres`
- `scan_TAD_domains_n_biosamples`
- `scan_TAD_boundary_overlaps_any`

### `genes_by_tier` / `lncrnas_by_tier`

```
genes_by_tier ────────────────────────────────────
│
├── "0–100kb"     : ["HLA-A", "HLA-G"]
├── "100–250kb"   : ["MICB"]
├── "250–500kb"   : []
└── "500–1000kb"  : ["MICA"]
```

### `gene_links` column (ATAC-specific)

```
gene_links ───────────────────────────────────────────────────────────
│
├── "HLA-A" ──────────────────────────────────────────────────────────
│   ├── "gene_id"       : "ENSG00000206503"
│   ├── "gene_type"     : "protein_coding"
│   ├── "dist_to_tss"   : 5000
│   ├── "tier"          : "0–100kb"
│   ├── "tss_position"  : 29942470
│   ├── "strand"        : "+"
│   │
│   └── "body_overlap" ───────────────────────────────────────────────
│       ├── "overlaps"               : true
│       ├── "overlap_bp"             : 800
│       ├── "overlap_interval"       : [29940000, 29940800]
│       ├── "overlap_frac_of_peak"   : 0.67
│       ├── "overlap_frac_of_gene"   : 0.001
│       └── "overlap_type"           : "promoter" | "gene_body" | null
│
├── "HLA-G" ── (same shape)
└── ...
```

### `lncrna_links` column

```
(Identical shape to gene_links above,
 keyed by lncRNA name instead of gene name)
```

### `ccre_links` column

```
ccre_links ── list ───────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "cCRE_id"    : "EH38E1234567.1"
│   ├── "ENCODE_id"  : "EH38E1234567"
│   ├── "raw_type"   : "pELS,CTCF-bound"
│   ├── "distance"   : 0
│   │
│   └── "overlap" ────────────────────────────────────────────────────
│       ├── "overlaps"               : true
│       ├── "overlap_bp"             : 500
│       ├── "overlap_interval"       : [29940200, 29940700]
│       ├── "overlap_frac_of_peak"   : 0.42
│       └── "overlap_frac_of_ccre"   : 0.25
│
├── [1] ── (same shape)
└── ...
```

### `ccre_types` column

```
ccre_types ───────────────────────────────────────
│
├── "pELS"      : 3
├── "dELS"      : 1
├── "CTCF-only" : 1
└── "PLS"       : 0
```

### `TAD_domains` (same nested layout as in element-focus `TAD_domains`)

### `TAD_boundary_overlaps` column (ATAC-specific)

```
TAD_boundary_overlaps ────────────────────────────────────────────────
│
├── "Kim_T47D" ───────────────────────────────────────────────────────
│   ├── "overlaps_boundary" : true
│   ├── "n_boundaries"      : 1
│   │
│   └── "boundaries" ── list ─────────────────────────────────────────
│       └── [0] ──────────────────────────────────────────────────────
│           ├── "boundary_id"       : "Kim_T47D_boundary_042"
│           ├── "overlap_bp"        : 200
│           ├── "overlap_interval"  : [29940000, 29940200]
│           └── "boundary_strength" : "strong"
│
├── "Kim_HCC70" ── (same shape)
└── ...per biosample
```

---

## 16. External Annotation Inputs (from `config.PATHS`)

These are not pipeline outputs; they are external metadata tables referenced by
`pipeline/config.py` and typically joined by TCGA sample barcode.

### `PATHS.immune_subtype_annotations`

**File**: `annotations/BRCA_immune_subtypes_advanced.tsv`  
**Join key**: `sample_id` (e.g. `TCGA-3C-AAAU-01`)

**Observed header (first columns)**:

```
sample_id
PAM50_recomputed
PAM50Call_RNAseq
PAM50_final
thornsson_subtype
CPE
```

Notes:
- `PAM50Call_RNAseq` may be empty for some rows.
- `thornsson_subtype` appears present but is often empty in this table.

### `PATHS.thornsson_immune_table`

**File**: `annotations/Thornsson_immune_table.tsv`  
**Join key**: `TCGA Participant Barcode` (e.g. `TCGA-A2-A0CO-01`)

This is a wide table combining immune subtype (C1..C6) with many immune and clinical
features. The observed header includes (non-exhaustive):

```
TCGA Participant Barcode
TCGA Study
Immune Subtype
TCGA Subtype
pathologic_stage
pathologic_T
pathologic_N
pathologic_M
PAM50_final
... many immune infiltration / neoantigen / survival fields ...
Wound Healing
```

Notes:
- Column names include spaces and punctuation (keep as-is when loading with pandas).
- Many fields are numeric fractions/scores; missing values are common.

### `PATHS.brca_clinical`

**File**: `annotations/BRCA_clinical` (tab-separated; despite name, this is a file)  
**Join key**: `sampleID` (e.g. `TCGA-3C-AAAU-01`)

This is a very wide clinical table with TCGA clinical fields and multiple “nature2012”
cluster assignments. The observed header begins with:

```
sampleID
AJCC_Stage_nature2012
Age_at_Initial_Pathologic_Diagnosis_nature2012
CN_Clusters_nature2012
... (many TCGA clinical form fields) ...
vital_status
year_of_initial_pathologic_diagnosis
... multiple _GENOMIC_ID_* columns ...
```

Notes:
- There are hundreds of columns; treat it as “sample-level metadata”.
- Many columns are sparsely populated.

<a id="hichip-tcga-brca-processed"></a>

### HiChIP TCGA BRCA processed matrix (external)

**Path (config)**: `PATHS.hichip_tcga_processed_csv`  
**On-disk name (typical)**: `data/HiCHIP/TCGA_BRCA_PROCESSED.csv` (BRCA cohort–level HiChIP / FitHiChIP–style processed interactions; TCGA-wide columns)

**Granularity**: one row per interaction (anchor pair); sample signal in wide columns.

**Format**: comma-separated CSV.

**Leading columns (first six; coordinate + anchor pair)**:

```
chr1   : str   (e.g. "chr1")
s1     : int   (bp, anchor 1 start)
e1     : int   (bp, anchor 1 end)
chr2   : str
s2     : int
e2     : int
```

**Sample columns**: from the **7th column** onward, headers are **TCGA case / participant** barcodes (e.g. `TCGA-A7-A0CH`) — not `...-01` sample-type ids. Values are numeric scores / normalized counts for that interaction in that case (units follow your upstream FitHiChIP / processing).

**Join key**: **`participant`** (derive from each `TCGA-XX-YYYY-...` column if you need `sample` / `sample_vial`; see `pipeline/md/SAMPLE_ID_MATCHING_GUIDE.md`).

**Column-derived manifest**: `PATHS.hichip_samples_tsv` → `annotations/HiCHIP/samples.tsv` (from `scripts/annotations/build_annotation_sample_manifests.py`).

**Relationship to cell-line HiChIP loops**: per-line HiChIP loop BED/BEDPE files under `PATHS.hichip_dir` (e.g. `MCF7/…`) feed `gene_links["hichip"]` in the element table; this **TCGA matrix** is a separate cohort resource for tumor-level chromatin interaction scores.

<a id="mirna-xena-arm-specific"></a>

### miRNA Xena arm-specific expression (external)

**Path (config)**: `PATHS.mirna_expression_tsv`  
**On-disk name (typical)**: `data/miRNA/XENA_mirna_arm_specific.tsv`

**Granularity**: wide matrix — **rows** = mature-arm features (often **MIMAT** accessions or arm labels, depending on how the matrix was built); **columns** = TCGA samples.

**Format**: tab-separated TSV.

**Header pattern**: first column is commonly `sample` (row id / feature id column); remaining columns are TCGA barcodes (often **`sample`**-level ids like `TCGA-…-01`, sometimes vial-level — treat headers as opaque until normalized).

**Units**: typically **log-transformed** read counts (e.g. log2(RPM+1)) per Xena / upstream pipeline; confirm against your file README if you need exact units.

**Join key**: use normalized keys from **`annotations/_normalized/`** companion tables when present, or column headers as **`raw_sample_id`** and map via `pipeline/sample_ids.normalize_tcga_id` (see `SAMPLE_ID_MATCHING_GUIDE.md`).

**Column-derived manifest**: `PATHS.mirna_samples_tsv` → `annotations/miRNA/samples.tsv` (from `scripts/annotations/build_annotation_sample_manifests.py`).

**Legacy note**: `PATHS.mirna_expression_legacy_tsv` points at an older BRCA miRNA matrix for reference only; the pipeline default expression path for arm-level work is **`mirna_expression_tsv`** above.

---

## 3. SNV Table

**File**: `<prefix>_snv_variants.csv` / `.parquet`
**Granularity**: One row per ALT allele per variant site (after somatic filtering)

**Units**:
- Genomic positions (`pos`): **bp**
- Depth / allele depths (`*_DP`, `*_AD_*`): **reads**
- VAF (`normal_vaf`, `tumor_vaf`): fraction **0–1**
- `qual`, `TLOD`: unitless caller scores

### Flat columns

```
┌──────────────────────────┬──────────┬──────────────────────────────┐
│ Column                   │ Type     │ Example                      │
├──────────────────────────┼──────────┼──────────────────────────────┤
│ chrom                    │ str      │ "chr6"                       │
│ pos                      │ int      │ 29942470                     │
│ id                       │ str|None │ "rs12345"                    │
│ ref                      │ str      │ "C"                          │
│ alt                      │ str      │ "T"                          │
│ qual                     │ float    │ 30.0                         │
│ filter                   │ str      │ "PASS"                       │
│ TLOD                     │ float    │ 12.5                         │
│ POPAF                    │ float    │ 6.0                          │
│ <NORMAL>_AD_ref          │ int      │ 45                           │
│ <NORMAL>_AD_alt          │ int      │ 0                            │
│ <NORMAL>_DP              │ int      │ 45                           │
│ <TUMOR>_AD_ref           │ int      │ 38                           │
│ <TUMOR>_AD_alt           │ int      │ 12                           │
│ <TUMOR>_DP               │ int      │ 50                           │
│ normal_vaf               │ float    │ 0.0                          │
│ tumor_vaf                │ float    │ 0.24                         │
│ gene_symbols             │ str      │ "HLA-A"                      │
│ hits_canonical           │ bool     │ true                         │
│ has_missense             │ bool     │ true                         │
│ has_nonsense             │ bool     │ false                        │
│ has_frameshift           │ bool     │ false                        │
│ has_splice_effect        │ bool     │ false                        │
│ has_missense_canonical   │ bool     │ true                         │
│ has_nonsense_canonical   │ bool     │ false                        │
│ has_frameshift_canonical │ bool     │ false                        │
│ has_splice_effect_canon. │ bool     │ false                        │
│ source_file              │ str      │ "sample1.vep.vcf"  (batch)   │
└──────────────────────────┴──────────┴──────────────────────────────┘
```

**Coding / CDS hit note (SNV):**
- There is **no single** `cds_flag` column on SNVs. The SNV module represents “coding hit” via **VEP consequence classes**:
  - `has_missense`, `has_nonsense`, `has_frameshift` (and `has_splice_effect` as a separate splice proxy).
- Convenience filter: `pipeline/SNV/vcf_loader.py:get_coding_variants()` returns variants where any of {missense, nonsense, frameshift} is True.
- Transcript-level `gene_hits` entries do include `CDS_position` / `Protein_position` when VEP provides them (see below), but that is **annotation metadata**, not a boolean “CDS overlap” flag.

### `gene_hits` column

```
gene_hits ── list ────────────────────────────────────────────────────
│
├── [0] ── (one entry per affected transcript) ───────────────────────
│   ├── "Allele"           : "T"
│   ├── "Consequence"      : "missense_variant"
│   ├── "IMPACT"           : "MODERATE"
│   ├── "SYMBOL"           : "HLA-A"
│   ├── "Gene"             : "ENSG00000206503"
│   ├── "Feature_type"     : "Transcript"
│   ├── "Feature"          : "ENST00000376802"
│   ├── "BIOTYPE"          : "protein_coding"
│   ├── "EXON"             : "3/8"
│   ├── "INTRON"           : null
│   ├── "HGVSc"            : "c.123C>T"
│   ├── "HGVSp"            : "p.Ala41Val"
│   ├── "cDNA_position"    : "456"
│   ├── "CDS_position"     : "123"
│   ├── "Protein_position" : "41"
│   ├── "Amino_acids"      : "A/V"
│   ├── "Codons"           : "gCc/gTc"
│   ├── "CANONICAL"        : "YES"
│   ├── "MANE_SELECT"      : "NM_002116.8"
│   ├── "SIFT"             : "deleterious(0.01)"
│   ├── "PolyPhen"         : "probably_damaging(0.95)"
│   ├── "gnomADe_AF"       : 0.0001
│   ├── "gnomADg_AF"       : 0.00012
│   ├── "MAX_AF"           : 0.0002
│   └── "CLIN_SIG"         : null
│
├── [1] ── (another transcript of same or different gene)
└── ...
```

### `regulatory_hits` column

```
regulatory_hits ── list ──────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "Allele"        : "T"
│   ├── "Consequence"   : "regulatory_region_variant"
│   ├── "IMPACT"        : "MODIFIER"
│   ├── "Feature_type"  : "RegulatoryFeature"
│   ├── "Feature"       : "ENSR00000123456"
│   ├── "BIOTYPE"       : "promoter"
│   ├── "DISTANCE"      : null
│   └── "VARIANT_CLASS" : "SNV"
│
└── [1] ── ...
```

### `motif_hits` column

```
motif_hits ── list ───────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "Allele"                 : "T"
│   ├── "Consequence"            : "TF_binding_site_variant"
│   ├── "Feature_type"           : "MotifFeature"
│   ├── "Feature"                : "ENSM00000012345"
│   ├── "MOTIF_NAME"            : "STAT1::STAT2"
│   ├── "MOTIF_POS"             : 8
│   ├── "HIGH_INF_POS"          : "Y"
│   ├── "MOTIF_SCORE_CHANGE"    : -0.15
│   └── "TRANSCRIPTION_FACTORS" : "STAT1,STAT2"
│
└── [1] ── ...
```

### `fimo_hits` column (reference-window MEME FIMO)

**Source**: `pipeline/SNV/snv_fimo.py` → `annotate_snvs_with_fimo`, invoked from `load_mutect_snv_vcf` / `load_mutect_snv_batch` when `run_fimo=True` (cohort runner: `--snv-fimo`, default on).

**Orthogonal to `motif_hits`**: `motif_hits` comes from VEP **MotifFeature** CSQ lines (alternate allele / TF binding consequence). `fimo_hits` scans a **fixed-width hg38 reference** window around each retained variant with MEME **FIMO** and the same motif database as SV FIMO (`PATHS.sv_meme_file` by default). PWMs are selected using **`SV_TARGET_TF_SYMBOLS`** in `pipeline/config.py` when (re)building that MEME file (see `pipeline/md/module_specific_processing_md/sv_→_motif_pipeline_vep_fimo.md`). Bin resolution: `resolve_fimo_argv()` (`APM_FIMO_BIN`, `PATH`, `VEP_ENV` / `vep_env`, or common `~/miniforge3/envs/vep_env/bin/fimo`, etc.).

Each list element is one motif occurrence overlapping the extracted window (genomic coordinates are **0-based start**, **half-open end**). There is **no** `distance_to_pos` field; use `variant_pos` (1-based VCF `POS`) vs `start`/`end` if you need an offset.

```
fimo_hits ── list ─────────────────────────────────────────────────────
│
├── [0] ───────────────────────────────────────────────────────────────
│   ├── "TF"                 : str   (prefix of MEME motif_id, e.g. "CTCF")
│   ├── "motif_id"           : str   (e.g. "MA0139.1")
│   ├── "rel_start"          : int   (1-based within FIMO window sequence)
│   ├── "rel_stop"           : int   (1-based inclusive in FIMO TSV; mapped to half-open hg38 end)
│   ├── "start"              : int   (hg38 0-based inclusive)
│   ├── "end"                : int   (hg38 0-based exclusive)
│   ├── "strand"             : str|None
│   ├── "p_value"            : float|None
│   ├── "q_value"            : float|None
│   ├── "score"              : float|None
│   ├── "variant_pos"        : int   (1-based VCF POS for this row)
│   └── "matched_sequence"   : str|None
│
└── [1] ── ...
```

### `snv_chip_hits` / `snv_chip_aggregate` (unified ChIP overlap)

**Source**: ``PATHS.chip_unified`` (same parquet as cCRE ``chip_hits`` / SV ``chip_hits``), attached in ``load_mutect_snv_vcf`` when ``run_chip=True`` (default in cohort runner: ``--snv-chip``). Implementation: ``pipeline/SNV/snv_chip.py``.

**Overlap rule:** strict **single-base** overlap with VCF ``POS`` (1-based): a peak ``[start, end)`` in 0-based BED coordinates overlaps iff ``start <= POS - 1`` and ``end >= POS``. **No** distance window.

``snv_chip_hits`` — list (one dict per overlapping peak), analogous to SV ``chip_hits`` entries:

```
snv_chip_hits ── list ────────────────────────────────────────────────
│
├── [0]
│   ├── "tf"                              : str
│   ├── "cell_type"                       : str   (canonical label)
│   ├── "cell_subtype"                    : str   (ChIP-Atlas subtype map; may be "")
│   ├── "source"                          : "ENCODE" | "CHIP_ATLAS"
│   ├── "sample_id"                       : str|None   (BED stem in unified table)
│   ├── "score_norm"                      : float|None
│   ├── "chrom"                           : str
│   ├── "peak_start" / "peak_end"        : int   (0-based half-open peak)
│   ├── "variant_pos"                     : int   (1-based VCF POS)
│   ├── "overlap_bp"                      : 1
│   ├── "overlap_start_0based"            : int
│   ├── "overlap_end_0based_exclusive"    : int
│   └── "stratum"                         : str   (subtype:* / context:* — see code)
```

``snv_chip_aggregate`` — per-row dict summarizing the same hits by **TF × source × stratum** (``mean_score_norm``, ``max_score_norm``, ``n_peaks``, ``n_chip_samples``, ``chip_sample_ids``).

### `cCRE_hits` column

```
cCRE_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "cCRE_id"            : "EH38E1594103.1"
│   ├── "elem_type"          : "pELS,CTCF-bound"
│   ├── "chrom"              : "chr6"
│   ├── "elem_start"         : 29942000
│   ├── "elem_end"           : 29942800
│   └── "genes_by_exact_dist": "HLA-A:5230,HLA-G:89420"
│
└── [1] ── ...
```

### `mirna_hits` column (if miRNA matching enabled)

```
mirna_hits ── list ───────────────────────────────────────────────────
│
├── [0] ── (one entry per overlapping miRNA feature)
│   ├── fields vary by miRNA reference structure
│   └── ...
└── ...
```

### Batch summary file (`<prefix>_snv_summary.json`)

```json
{
  "n_total_variants": 342,
  "n_unique_positions": 340,
  "n_missense": 12,
  "n_nonsense": 2,
  "n_frameshift": 1,
  "n_splice": 0,
  "n_hits_canonical": 8,
  "genes_affected": ["HLA-A", "HLA-B", "TAP1"],
  "n_genes_affected": 3,
  "n_variants_in_ccres": 45,
  "n_variants_with_fimo_hits": 0,
  "n_variants_with_chip_hits": 0,
  "ccre_type_counts": {"pELS": 20, "dELS": 15, "CTCF-only": 10},
  "tumor_vaf_median": 0.18,
  "tumor_vaf_mean": 0.22,
  "variants_per_chrom": {"chr6": 150, "chr1": 50, ...}
}
```

---

## 4. SV Table

**File**: `<sample>_strict_sv_set.csv` (in `07_final_sv_with_fimo/`)
**Granularity**: One row per structural variant call

**Units**:
- Genomic positions (`pos`, `END`, breakpoint positions): **bp**
- `SVLEN`: **bp**
- Read support counts (`*_alt`, `*_sr_alt`, `*_pr_alt`): **reads**
- Signed distances / overlap sizes: **bp**

### Flat columns

```
┌──────────────────────────────┬──────────┬──────────────────────────┐
│ Column                       │ Type     │ Example                  │
├──────────────────────────────┼──────────┼──────────────────────────┤
│ id                           │ str      │ "MantaDEL:42:0:1:0:0:0"  │
│ chrom                        │ str      │ "chr6"                   │
│ pos                          │ int      │ 29940000                 │
│ END                          │ int      │ 29960000                 │
│ SVTYPE                       │ str      │ "DEL" | "DUP" | "INV"   │
│                              │          │   | "INS" | "BND"       │
│ SVLEN                        │ int|None │ -20000                   │
│ ref                          │ str      │ "A"                      │
│ alt                          │ str      │ "<DEL>"                  │
│ qual                         │ float    │ 999.0                    │
│ filter                       │ str      │ "PASS"                   │
│ SOMATICSCORE                 │ int      │ 45                       │
│ normal_alt                   │ int      │ 0                        │
│ tumor_alt                    │ int      │ 12                       │
│ tumor_sr_alt                 │ int      │ 5                        │
│ tumor_pr_alt                 │ int      │ 7                        │
│ bnd_remote_chrom             │ str|None │ null (or "chr9" for BND) │
│ bnd_remote_pos               │ int|None │ null                     │
│ gene_symbols                 │ str      │ "HLA-A,HLA-G"           │
│ hits_canonical               │ bool     │ true                     │
│ has_missense                 │ bool     │ false                    │
│ has_nonsense                 │ bool     │ false                    │
│ has_frameshift               │ bool     │ true                     │
│ has_splice_effect            │ bool     │ true                     │
│ has_*_canonical              │ bool     │ ...                      │
│ has_*_mane                   │ bool     │ ...                      │
└──────────────────────────────┴──────────┴──────────────────────────┘
```

### `gene_hits` column (SV-specific spatial mapping)

```
gene_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "gene_name"          : "HLA-A"
│   ├── "gene_id"            : "ENSG00000206503"
│   ├── "strand"             : "+"
│   ├── "signed_dist"        : 0         ← 0 means overlap
│   ├── "overlap_start"      : 29942000
│   ├── "overlap_end"        : 29945000
│   ├── "overlap_bp"         : 3000
│   ├── "overlap_percent"    : 0.15
│   │
│   │   region classification flags (0 or 1):
│   ├── "promoter_flag"      : 1
│   ├── "gene_body_flag"     : 1
│   ├── "cds_flag"           : 1         ← **only** when a CDS feature row is overlapped
│   ├── "utr_flag"           : 0
│   ├── "exon_flag"          : 1
│   ├── "mane_cds_flag"      : 0|1       ← best-effort MANE CDS overlap (if MANE is annotated in features)
│   ├── "intron_only_flag"   : 0
│   ├── "upstream_5kb_flag"  : 0
│   ├── "downstream_5kb_flag": 0
│   ├── "start_codon_flag"   : 0
│   ├── "stop_codon_flag"    : 0
│   │
│   ├── "region_hit"         : "promoter+exon"
│   ├── "hit_side"           : "span" | "point" | "bp1" | "bp2"
│   ├── "transcript_id"      : "ENST00000376802"
│   └── "transcript_type"    : "protein_coding"
│
│   feature-interval attribution (lists; may be empty):
│   ├── "exon_interval_ids"        : list[str]  (e.g. ["ENST...:exon:40691570-40691867", ...])
│   ├── "exon_ids"                 : list[str]  (stable GTF exon IDs when available; else empty)
│   ├── "cds_interval_ids"         : list[str]
│   ├── "utr_interval_ids"         : list[str]
│   ├── "start_codon_interval_ids" : list[str]
│   └── "stop_codon_interval_ids"  : list[str]
│
├── [1] ── ...
└── ...
```

### `elem_hits` column (SV-specific)

```
elem_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "elem_id"           : "EH38E1594103.1"
│   ├── "elem_type"         : "pELS"
│   ├── "chrom"             : "chr6"
│   ├── "elem_start"        : 29942000
│   ├── "elem_end"          : 29942800
│   ├── "signed_dist"       : 0
│   ├── "overlap_start"     : 29942000
│   ├── "overlap_end"       : 29942800
│   ├── "overlap_bp"        : 800
│   ├── "overlap_percent"   : 1.0
│   │
│   │   region flags:
│   ├── "overlaps_flag"     : 1
│   ├── "proximal_flag"     : 0
│   ├── "distal_flag"       : 0
│   ├── "region_hit"        : "overlaps"
│   ├── "hit_side"          : "span"
│   │
│   │   FIMO motif hits (added after motif scanning step):
│   └── "motif_hits" ── list ─────────────────────────────────────────
│       ├── [0] ──────────────────────────────────────────────────────
│       │   ├── "start"           : 29942100
│       │   ├── "end"             : 29942118
│       │   ├── "TF"             : "CTCF"
│       │   ├── "motif_id"       : "MA0139.1"
│       │   ├── "score"          : 22.5
│       │   ├── "p_value"        : 3.2e-6
│       │   ├── "q_value"        : 1.5e-4
│       │   ├── "strand"         : "+"
│       │   └── "distance_to_pos": 100
│       └── [1] ── ...
│
├── [1] ── ...
└── ...
```

### `lncRNA_hits` column

```
(Same shape as gene_hits, keyed for lncRNA features)
```

**Feature identity / transcript tracking note (SV):**
- SV `gene_hits` currently stores **flags** (`cds_flag`, `exon_flag`, etc.) and a **single** `transcript_id`/`transcript_type` chosen from overlapped `transcript` rows (prefers MANE transcript if annotated).
- It does **not** store a list of overlapped **exon IDs** (stable exon identifiers) or the exact exon numbers. If you want exon-level attribution, we’d need to extend SV mapping to emit the specific overlapped feature row(s) (e.g. transcript_id + feature start/end, or exon number if available in the GENCODE features table).

---

## 15. Cohort covariates table (post-processing)

**Doc**: `pipeline/md/COHORT_COVARIATES.md`  
**Build**: `pipeline/covariates/build_covariates.py` (`build_cohort_covariates`)  
**Primary output**: `data/covariates/<run_id>/cohort_covariates.parquet` (and optional CSV)

**Optional / extended providers** (columns appear when enabled at build time; see `pipeline/md/COHORT_COVARIATES.md`):

- **`--ddr-scores`**: BRCA-filtered **`DDRscores.txt`‑style** table merged as additional numeric/categorical columns (HRD / purity / signature‑3 / TP53 score, etc.).
- **HiChIP TCGA processed matrix** (`PATHS.hichip_tcga_processed_csv`): participant-level summaries **lifted to `sample_vial`** when HiChIP provider is enabled.
- **Pipeline-native SNV** (`PATHS.snv_output_dir` per-sample `*_snv_summary.json`): best-effort driver-oriented columns (joined using `annotations/SNV/samples.tsv` to recover the true tumor vial from UUID-ish run folders).
- **Sample coverage**: prefers `analysis/sample_coverage/output/current/` when present (stable mirror), else the latest `run_*` directory.

**Replication context (base index)**: the covariates base table can include **`n_vials_per_sample`**, **`n_vials_per_participant`**, and related flags so downstream models see TCGA vial duplication explicitly.

### SV→gene “mirroring” outputs (gene-centric disruption summaries)

**General-purpose mirroring code**: `pipeline/SV/gene_mirroring.py`  
**Covariates integration**: `pipeline/covariates/providers/sv_disruptions.py` (wraps the general mirroring and joins by `sample_vial`)  
**Input**: SV strict tables under `data/SV/pipeline_output/07_final_sv_with_fimo/*_strict_sv_set.csv`

**Global SV disruption summary columns** (always present; prefixed when merged into covariates table):
- `SV_any_gene__n_unique_genes_hit`
- `SV_any_gene__n_hits`, `SV_any_gene__n_bp_hits`, `SV_any_gene__n_span_hits`
- `SV_any_gene__n_cds_hits` (if `cds_flag` is available/populated), `SV_any_gene__n_exon_hits`
- `SV_any_gene__n_cds_bp_hits`, `SV_any_gene__n_exon_bp_hits`
- `SV_any_gene__n_promoter_hits`
- `SV_any_gene__n_start_stop_hits`

**Per-gene SV disruption columns** (only for configured genes-of-interest; e.g. TP53/PTEN/PIK3CA/BRCA1/BRCA2):
- `SV_<GENE>__any_hit`, `SV_<GENE>__any_promoter_hit`
- `SV_<GENE>__any_exon_hit`, `SV_<GENE>__any_exon_bp_hit`
- `SV_<GENE>__any_cds_hit`, `SV_<GENE>__any_cds_bp_hit`
- counts: `SV_<GENE>__n_hits`, `__n_bp_hits`, `__n_promoter_hits`, `__n_exon_hits`, `__n_exon_bp_hits`, `__n_cds_hits`, `__n_cds_bp_hits`

### `mir_hits` column (SV ↔ miRNA loci; optional)

When miRNA loci are provided (defaults to `PATHS.mirna_path`, built from `data/miRNA/hsa.gff`),
SVs are also mapped to nearby/overlapping **pre-miRNA loci**. Each hit may include mature-arm
identifiers for joining to arm-specific expression (Xena `MIMAT...`) and miRTarBase.

```
mir_hits ── list ──────────────────────────────────────────────────────
│
├── [0]
│   ├── "gene_name"         : "Hsa-Mir-8-P1a_pre"
│   ├── "gene_id"           : "MI0000342"
│   ├── "mature_names"      : "Hsa-Mir-8-P1a_3p,Hsa-Mir-8-P1a_5p*" | ""   (optional)
│   ├── "mature_accessions" : "MIMAT0000318,MIMAT0004571" | ""             (optional; comma-separated)
│   ├── "strand"            : "+" | "-"
│   ├── "signed_dist"       : int
│   ├── "overlap_start"     : int|None
│   ├── "overlap_end"       : int|None
│   ├── "overlap_bp"        : int
│   ├── "overlap_percent"   : float
│   ├── "region_hit"        : "overlaps" | "proximal" | "distal"
│   └── "hit_side"          : "span" | "point" | "bp1" | "bp2"
└── ...
```

### `gene_hits_vep` / `regulatory_hits_vep` / `motif_hits` columns

```
(Same shapes as the SNV gene_hits / regulatory_hits / motif_hits
 described in Section 3 above)
```

### `flank_motif_hits` column

```
flank_motif_hits ── list ─────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "start"           : 29939850
│   ├── "end"             : 29939868
│   ├── "TF"             : "STAT1"
│   ├── "motif_id"       : "MA0137.3"
│   ├── "score"          : 18.3
│   ├── "p_value"        : 5.1e-5
│   ├── "q_value"        : 2.0e-3
│   ├── "strand"         : "-"
│   ├── "distance_to_pos": -150
│   └── "flank_side"     : "left" | "right"
│
└── [1] ── ...
```

### `chip_hits` column (SV ↔ ChIP-seq disruption evidence)

ChIP hits are added by `pipeline/CHIP/sv_chip_intersect.py` Step 7 in the SV pipeline.

**Type**:
- In memory / parquet: `list[dict]`
- When written to CSV by the batch annotator: serialized with `repr(...)` (so it will be a string like `"[]"` or `"[{...}, {...}]"`).

```
chip_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "tf"            : str   (e.g. "CTCF")
│   ├── "cell_type"     : str   (e.g. "MCF7")
│   ├── "source"        : str   (e.g. "ENCODE" | "CHIP_ATLAS")
│   ├── "score_norm"    : float|None
│   ├── "chrom"         : str
│   ├── "peak_start"    : int
│   ├── "peak_end"      : int
│   ├── "signed_dist"   : int
│   ├── "overlap_start" : int|None
│   ├── "overlap_end"   : int|None
│   ├── "overlap_bp"    : int
│   ├── "overlaps_flag" : int   (0/1)
│   └── "hit_side"      : "span" | "bp1" | "bp2"
│
└── ... (for BND SVs, breakpoint-specific keys may also appear)
    ├── "bp_index"  : 1|2
    ├── "bp_chrom"  : str
    ├── "bp_pos"    : int
    ├── "mate_chrom": str|None
    └── "mate_pos"  : int|None
```

---

## 5. CNV Segment Table

**File**: per-sample annotated CSV in `PATHS.cnv_output_dir`
**Granularity**: One row per CNV segment

**Units**:
- Segment coordinates (`start`, `end`): **bp**
- `segment_length`: **bp**
- `num_probes`: **probes**
- `segment_mean`: typically **log2(copy ratio)** (unitless; depends on upstream ASCAT/GDC segment format)
- Signed distances / overlap sizes: **bp**

### Flat columns

```
┌──────────────────────────┬──────────┬──────────────────────────────┐
│ Column                   │ Type     │ Example                      │
├──────────────────────────┼──────────┼──────────────────────────────┤
│ chrom                    │ str      │ "chr6"                       │
│ start                    │ int      │ 29000000                     │
│ end                      │ int      │ 31000000                     │
│ num_probes               │ int      │ 450                          │
│ segment_mean             │ float    │ -0.85                        │
│ sample_id                │ str      │ "TCGA-A2-A0CM"              │
│                          │          │                              │
│ (derived features from features.py):                               │
│ cn_state                 │ str      │ "loss" | "gain" | "neutral"  │
│ loh_flag                 │ bool     │ false                        │
│ segment_length           │ int      │ 2000000                      │
└──────────────────────────┴──────────┴──────────────────────────────┘
```

### `gene_hits` column

```
gene_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "gene_name"        : "HLA-A"
│   ├── "gene_id"          : "ENSG00000206503"
│   ├── "strand"           : "+"
│   ├── "signed_dist"      : 0
│   ├── "overlap_start"    : 29942000
│   ├── "overlap_end"      : 29945000
│   ├── "overlap_bp"       : 3000
│   ├── "overlap_percent"  : 1.0
│   │
│   │   region flags:
│   ├── "promoter_flag"    : 1
│   ├── "gene_body_flag"   : 1
│   ├── "exon_flag"        : 1
│   ├── "intron_only_flag" : 0
│   ├── "upstream_5kb_flag": 0
│   ├── "downstream_5kb_flag": 0
│   └── "region_hit"       : "gene_body+promoter"
│
└── [1] ── ...
```

### `lncRNA_hits`, `mir_hits` columns

```
(Same spatial hit structure as `gene_hits` above, but different feature identity.)

lncRNAs_hits ── list ─────────────────────────────────────────────────
│
├── [0]
│   ├── "gene_name"        : "<lncRNA name>"
│   ├── "gene_id"          : "<ENSG...>" | null
│   ├── "strand"           : "+" | "-"
│   ├── "signed_dist"      : int
│   ├── "overlap_start"    : int
│   ├── "overlap_end"      : int
│   ├── "overlap_bp"       : int
│   ├── "overlap_percent"  : float
│   ├── region flags (same keys as gene_hits; derived from GTF subfeatures when available)
│   └── "region_hit"       : str  (e.g. "gene_body+promoter", "upstream", "downstream", ...)
└── ...

mir_hits ── list ──────────────────────────────────────────────────────
│
├── [0] ── (one entry per miRNA locus within the CNV window)
│   ├── "gene_name"        : "Hsa-Mir-551-P1_pre"
│   ├── "gene_id"          : "MI0003556"
│   ├── "mature_names"     : "Hsa-Mir-551-P1_3p,Hsa-Mir-551-P1_5p*" | ""   (optional)
│   ├── "mature_accessions": "MIMAT0003214" | ""                               (optional; comma-separated)
│   ├── "strand"           : "+" | "-"
│   ├── "signed_dist"      : int
│   ├── "overlap_start"    : int
│   ├── "overlap_end"      : int
│   ├── "overlap_bp"       : int
│   ├── "overlap_percent"  : float
│   └── "region_hit"       : str
└── ...

Notes:
- miRNA features come from `PATHS.mirna_path` (CNV-ready locus table, `data/miRNA/cnv_miRNA.csv`).
- Because miRNA loci do not have exon/intron/promoter subfeatures in this reference,
  the subfeature flags used for genes/lncRNAs are typically not populated.
```

### `elem_hits` column

```
elem_hits ── list ────────────────────────────────────────────────────
│
├── [0] ──────────────────────────────────────────────────────────────
│   ├── "elem_id"          : "EH38E1594103.1"
│   ├── "elem_type"        : "pELS"
│   ├── "chrom"            : "chr6"
│   ├── "elem_start"       : 29942000
│   ├── "elem_end"         : 29942800
│   ├── "signed_dist"      : 0
│   ├── "overlap_start"    : 29942000
│   ├── "overlap_end"      : 29942800
│   ├── "overlap_bp"       : 800
│   ├── "overlap_percent"  : 1.0
│   ├── "overlaps_flag"    : 1
│   ├── "proximal_flag"    : 0
│   ├── "distal_flag"      : 0
│   └── "region_hit"       : "overlaps"
│
└── [1] ── ...
```

---

## 6. Methylation Tables

**Units**:
- `beta_value`: fraction **0–1**
- `m_value`: logit(beta) (unitless)
- `detection_pval`: p-value **0–1**
- Aggregated methylation features (`*_beta_mean/median/std`): beta **0–1**

### 6a. Probe Reference Table

**File**: `reference/probe_annotations.csv`
**Granularity**: One row per methylation probe (450K or EPIC array)

```
┌───────────────────────────┬──────────┬──────────────────────────────┐
│ Column                    │ Type     │ Example                      │
├───────────────────────────┼──────────┼──────────────────────────────┤
│ probeID                   │ str      │ "cg00000029"                 │
│ chrom                     │ str      │ "chr16"                      │
│ start                     │ int      │ 53468112                     │
│ end                       │ int      │ 53468113                     │
│ strand                    │ str      │ "+"                          │
│ in_CGI                    │ bool     │ true                         │
│ CGI_context               │ str      │ "Island" | "Shore" |         │
│                           │          │   "Shelf" | "OpenSea"        │
│ distToTSS                 │ int      │ 1200                         │
│                           │          │                              │
│ in_promoter               │ bool     │ true  (within panel gene     │
│                           │          │         promoter region)      │
│ promoter_gene             │ str      │ "TAP1"                       │
│                           │          │                              │
│ in_gene_body              │ bool     │ true                         │
│ gene_body_gene            │ str      │ "TAP1"                       │
│                           │          │                              │
│ n_overlapping_ccres       │ int      │ 2                            │
│ overlapping_ccre_ids      │ str      │ "EH38E123,EH38E456"         │
│                           │          │                              │
│ in_lncrna_promoter        │ bool     │ false                        │
│ lncrna_promoter_gene      │ str      │ ""                           │
│                           │          │                              │
│ n_overlapping_atac        │ int      │ 1                            │
│ overlapping_atac_ids      │ str      │ "chr16:53468000-53469000"    │
│                           │          │                              │
│ TAD annotations (optional, one set of cols per biosample)           │
└───────────────────────────┴──────────┴──────────────────────────────┘
```

### 6b. Per-Sample Probe Table

**File**: `per_sample/<sample_id>/<sample_id>_probes.csv`
**Granularity**: One row per probe (joined with reference annotations)

```
┌─────────────────────┬──────────┬──────────────────────────────┐
│ Column              │ Type     │ Example                      │
├─────────────────────┼──────────┼──────────────────────────────┤
│ probeID             │ str      │ "cg00000029"                 │
│ beta_value          │ float    │ 0.82                         │
│ m_value             │ float    │ 2.18  (log2 transform)       │
│ detection_pval      │ float    │ 0.001                        │
│ is_valid            │ bool     │ true                         │
│                     │          │                              │
│ (+ all probe reference columns from 6a above)                │
└─────────────────────┴──────────┴──────────────────────────────┘
```

### 6c. Per-Sample Gene Aggregation

**File**: `per_sample/<sample_id>/<sample_id>_gene_meth.csv`
**Granularity**: One row per panel gene

```
┌────────────────────────────┬──────────┬──────────────────────────────┐
│ Column                     │ Type     │ Example                      │
├────────────────────────────┼──────────┼──────────────────────────────┤
│ gene_name                  │ str      │ "HLA-A"                      │
│ sample_id                  │ str      │ "TCGA-A2-A0CM"              │
│ promoter_beta_mean         │ float    │ 0.12                         │
│ promoter_beta_median       │ float    │ 0.10                         │
│ promoter_beta_std          │ float    │ 0.05                         │
│ promoter_n_probes          │ int      │ 8                            │
│ promoter_frac_hypermeth    │ float    │ 0.0  (fraction > 0.7)        │
│ promoter_frac_hypometh     │ float    │ 1.0  (fraction < 0.3)        │
│ promoter_CGI_beta_mean     │ float    │ 0.08                         │
│ promoter_shore_beta_mean   │ float    │ 0.15                         │
│ gene_body_beta_mean        │ float    │ 0.65                         │
│ gene_body_n_probes         │ int      │ 12                           │
└────────────────────────────┴──────────┴──────────────────────────────┘
```

### 6d. Per-Sample cCRE Aggregation

**File**: `per_sample/<sample_id>/<sample_id>_ccre_meth.csv`

```
┌────────────────────────────┬──────────┬──────────────────────────────┐
│ Column                     │ Type     │ Example                      │
├────────────────────────────┼──────────┼──────────────────────────────┤
│ cCRE_id                    │ str      │ "EH38E1594103.1"             │
│ sample_id                  │ str      │ "TCGA-A2-A0CM"              │
│ ccre_beta_mean             │ float    │ 0.45                         │
│ ccre_beta_median           │ float    │ 0.42                         │
│ ccre_n_probes              │ int      │ 3                            │
│ ccre_frac_hypermeth        │ float    │ 0.0                          │
│ ccre_frac_hypometh         │ float    │ 0.33                         │
│ ccre_CGI_overlap           │ bool     │ true                         │
│ ccre_CGI_beta_mean         │ float    │ 0.38                         │
└────────────────────────────┴──────────┴──────────────────────────────┘
```

### 6e. Per-Sample ATAC / lncRNA Aggregation

```
(Similar shape to cCRE aggregation:
 atac_peak_id + atac_beta_mean + atac_n_probes + ...
 lncrna_name + promoter_beta_mean + promoter_n_probes + ...)
```
### 6f. Cohort Matrix (Gene example)

**File**: `cohort/gene_meth_matrix.csv`
**Layout**: Genes (rows) × Samples (columns) + meta column

```
┌───────────┬──────────┬──────────┬──────────┬───────────────────────┐
│ gene_name │ TCGA-A2  │ TCGA-B6  │ TCGA-E2  │ meta                  │
│           │ -A0CM    │ -A0RH    │ -A1L9    │                       │
├───────────┼──────────┼──────────┼──────────┼───────────────────────┤
│ HLA-A     │ 0.12     │ 0.08     │ 0.65     │ {see below}           │
│ HLA-B     │ 0.15     │ NaN      │ 0.55     │                       │
│ TAP1      │ 0.05     │ 0.04     │ 0.72     │                       │
└───────────┴──────────┴──────────┴──────────┴───────────────────────┘

"meta" column structure (per row):

meta ─────────────────────────────────────────────
│
├── "TCGA-A2-A0CM" ──────────────────────────────
│   └── "n_probes" : 8
│
├── "TCGA-B6-A0RH"
│   └── "n_probes" : null    ← NaN beta → null
│
├── "TCGA-E2-A1L9"
│   └── "n_probes" : 6
│
└── ...per sample in matrix
```

### 6g. QC Summary

**File**: `cohort/sample_qc_summary.csv`

```
┌─────────────────┬──────────┬──────────────────────────────┐
│ Column          │ Type     │ Example                      │
├─────────────────┼──────────┼──────────────────────────────┤
│ sample_id       │ str      │ "TCGA-A2-A0CM"              │
│ n_total         │ int      │ 485577                       │
│ n_valid         │ int      │ 480123                       │
│ pct_valid       │ float    │ 98.9                         │
│ mean_beta       │ float    │ 0.48                         │
│ median_beta     │ float    │ 0.52                         │
└─────────────────┴──────────┴──────────────────────────────┘
```

---

## 7. RPPA Tables

**Units**:
- RPPA “expression” values and downstream analysis tables are typically **z-scores** (unitless; relative)
- Ratio features are unitless (often log-ratios depending on upstream preprocessing)

### 7a. Expression Matrix

**File**: internal (loaded via `load_rppa_dataset`)
**Layout**: Samples (rows) × Antibody Targets (columns)

```
┌──────────────┬──────────┬──────────┬──────────┬──────────┐
│ sample_id    │ STAT1    │ p-STAT1  │ B2M      │ PD-L1    │
│ (index)      │ _Y701    │          │          │          │
├──────────────┼──────────┼──────────┼──────────┼──────────┤
│ TCGA-A2-A0CM │ 0.35     │ -0.22    │ 1.12     │ -0.85    │
│ TCGA-B6-A0RH │ -0.15    │ 0.45     │ 0.88     │ 0.12     │
└──────────────┴──────────┴──────────┴──────────┴──────────┘
```

### 7b. Panel Scores

**File**: `panel_scores.csv`
**Granularity**: One row per sample

```
┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│ sample_id    │ IFN_     │ DDR_     │ cGAS_    │ check-   │ PI3K_    │
│ (index)      │ activated│ activ.   │ STING    │ point    │ AKT      │
├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│ TCGA-A2-A0CM │ 1.23     │ 0.85     │ 0.42     │ -0.31    │ -0.15    │
│ TCGA-B6-A0RH │ -0.65    │ 1.22     │ -0.20    │ 0.55     │ 0.88     │
└──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
```

Full panel list (must match `pipeline/rppa/rppa_schemas.empty_rppa_panel_scores` key order):

<!-- GENERATED:rppa_panel_score_columns BEGIN -->
IFN_activated, DDR_activation, cGAS_STING, DNA_repair, checkpoint, PI3K_AKT, mTOR, STAT3_suppressive, lymphocyte_infiltration, cytolytic, proliferation, apoptosis_balance
<!-- GENERATED:rppa_panel_score_columns END -->

### 7c. Signaling Blocks

**File**: `signaling_blocks.csv`

```
┌──────────────┬───────────┬───────────┬──────────┬───────────┬──────────┐
│ sample_id    │ ddr_sting │ sting_    │ stat3_   │ checkpoint│ pi3k_    │
│              │ _block    │ irf1_block│ override │ _escape   │ activated│
├──────────────┼───────────┼───────────┼──────────┼───────────┼──────────┤
│ TCGA-A2-A0CM │ false     │ false     │ false    │ false     │ false    │
│ TCGA-B6-A0RH │ true      │ true      │ false    │ true      │ true     │
└──────────────┴───────────┴───────────┴──────────┴───────────┴──────────┘

Additional blocks: ddr_ifn_chain_intact, immune_desert
```

### 7d. Combined Analysis Table

**File**: `rppa_analysis_combined.csv` / `.parquet`
**Granularity**: One row per sample (the final "wide" RPPA output)

```
┌───────────────────────────────────────────────────────────────────────────┐
│  COMBINED TABLE — all RPPA results per sample                             │
├───────────────────────────────┬──────────┬────────────────────────────────┤
│ Column                        │ Type     │ Description                    │
├───────────────────────────────┼──────────┼────────────────────────────────┤
│ visibility_score              │ float    │ composite immune vis. score    │
│ visibility_state              │ str      │ "visible" | "partially_visible"│
│                               │          │   | "invisible" | "unknown"   │
│ ddr_ifn_quadrant              │ str      │ "DDR+IFN+" | "DDR+IFN-" |     │
│                               │          │   "DDR-IFN+" | "DDR-IFN-"    │
│ panel_IFN_activated           │ float    │ z-score panel score            │
│ panel_DDR_activation          │ float    │                                │
│ panel_cGAS_STING              │ float    │                                │
│ panel_checkpoint              │ float    │                                │
│ panel_PI3K_AKT                │ float    │                                │
│ panel_STAT3_suppressive       │ float    │                                │
│ panel_cytolytic               │ float    │                                │
│ ... (all panel_* columns)     │          │                                │
│ ratio_<gene>_activation       │ float    │ phospho / total log ratio      │
│ ... (all ratio_* columns)     │          │                                │
│ block_ddr_sting_block         │ bool     │ pathway disruption flag        │
│ block_sting_irf1_block        │ bool     │                                │
│ block_stat3_override          │ bool     │                                │
│ block_checkpoint_escape       │ bool     │                                │
│ block_pi3k_activated          │ bool     │                                │
│ estimated_apm_capacity        │ float    │ antigen presentation estimate  │
└───────────────────────────────┴──────────┴────────────────────────────────┘
```

### 7e. Protein-RNA Discordance (optional)

**File**: `protein_rna_discordance.csv`

```
┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┐
│ sample_id    │ gene     │ protein  │ rna      │ discord- │ direction│
│              │          │ z_score  │ z_score  │ ance     │          │
├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│ TCGA-A2-A0CM │ HLA-A    │ 1.2      │ -0.5     │ 1.7      │ stabili- │
│              │          │          │          │          │ zation   │
└──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┘
```

---

## 8. Gene / lncRNA Tables

**Units**:
- Genomic coordinates (`start`, `end`, `tss`, promoter bounds): **bp**
- Distances / windows (e.g. lncRNA matching windows): **bp**

### 8a. Gene Tables (from `save_gene_tables`)

**Config**: `pipeline/config.py` → `PathConfig` fields (iterate via `PATHS.gene_table_output_paths()`).

**Gene symbol sets** (same module): `PIPELINE_GENE_PANEL` (Tier 1 full integration when extended mode is on; legacy 66 otherwise), `CNV_GENES`, `TIER2_MEDIUM_GENES`, `TIER3_CNV_ONLY_GENES`, `TIER4_READOUT_GENES`. The on-disk filenames below are fixed; the **`primary_*`** names are historical — rows for `primary_genes_*.csv` follow **`PIPELINE_GENE_PANEL`**, not the frozen `PRIMARY_GENES` list alone.

**Files** under `PATHS.working_dir` (default `data/`):

| Role | `PathConfig` attribute | On-disk basename |
|------|------------------------|------------------|
| Tier 1 pipeline panel — multifeature GENCODE rows | `genes_all_features` | `primary_genes_all_features.csv` |
| Tier 1 pipeline panel — `feature == "gene"` + promoters | `genes_only` | `primary_genes_only.csv` |
| CNV module union (`CNV_GENES`) — gene rows + promoters | `cnv_genes` | `cnv_genes.csv` |
| Matched panel lncRNAs — multifeature rows | `lncrnas_all_features` | `lncRNAs_genes_all_features.csv` |
| Tier 2 medium-depth — gene rows + promoters | `tier2_medium_genes_only` | `tier2_medium_genes_only.csv` |
| Tier 2 — multifeature | `tier2_medium_genes_all_features` | `tier2_medium_genes_all_features.csv` |
| Tier 3 CNV-only — gene rows + promoters | `tier3_cnv_only_genes_only` | `tier3_cnv_only_genes_only.csv` |
| Tier 3 — multifeature | `tier3_cnv_only_genes_all_features` | `tier3_cnv_only_genes_all_features.csv` |
| Tier 4 readout — gene rows + promoters | `tier4_readout_genes_only` | `tier4_readout_genes_only.csv` |
| Tier 4 — multifeature | `tier4_readout_genes_all_features` | `tier4_readout_genes_all_features.csv` |

Tier 2–4 CSVs are written only when `save_gene_tables` receives a non-empty **`genes_all_harmonized`** frame (full `load_genes` output after chromosome harmonization). Otherwise the pipeline logs a warning and skips those six files.

```
┌───────────────────────┬──────────┬──────────────────────────────┐
│ Column                │ Type     │ Example                      │
├───────────────────────┼──────────┼──────────────────────────────┤
│ chrom                 │ str      │ "chr6"                       │
│ start                 │ int      │ 29942470                     │
│ end                   │ int      │ 29945884                     │
│ strand                │ str      │ "+"                          │
│ gene_name             │ str      │ "HLA-A"                      │
│ gene_id               │ str      │ "ENSG00000206503"            │
│ gene_type             │ str      │ "protein_coding"             │
│ feature               │ str      │ "gene" | "transcript" | ...  │
│ transcript_id          │ str|None │ "ENST00000376802.8"          │
│ transcript_type        │ str|None │ "protein_coding"             │
│ transcript_name        │ str|None │ "HLA-A-201"                  │
│ exon_number            │ int|None │ 3   (only on exon rows)      │
│ exon_id                │ str|None │ "ENSE0000…" (exon rows)      │
│ tss                   │ int      │ 29942470                     │
│ prom_start            │ int      │ 29940470                     │
│ prom_end              │ int      │ 29942970                     │
│                       │          │                              │
│ (if TAD annotation ran):                                        │
│ TAD_domains           │ dict     │ (same structure as §1)       │
│                       │          │                              │
│ (if ATAC annotation ran):                                       │
│ peak_links            │ dict/list│ nearby peaks + distances      │
│                       │          │                              │
│ (lncRNA matching columns, if present):                          │
│ lncRNAs_within_1000kb │ list     │ ["LINC01149", "HCG18"]       │
└───────────────────────┴──────────┴──────────────────────────────┘
```

### 8b. lncRNA Matching Outputs

**Upstream lncRNA interval universe** (Steps 1–2 in `pipeline/main.py`): By default, lncRNA loci are **not** loaded from `lncRNA_matching/lncRNAs_with_genes_1000000bp.csv`. They are derived from **`PATHS.gencode_gtf_pq`** as every annotation row with `feature == "gene"` and `gene_type == "lncRNA"` (`lncrna_gene_intervals_from_annotation` in `pipeline/genes/gene_loader.py`). That file is still written **after** matching and is referenced as **`PATHS.lncrnas_genes_centric`** / **`PATHS.lncrna_csv`** for downstream modules (Methylation, SV) that expect the lncRNA-centric matched table. To restore the previous behaviour (use that CSV as the input interval list as well), set **`APM_LNCRNA_INPUT_LEGACY_CSV=1`**.

**Files** in `lncRNA_matching/`:

- `genes_lncRNAs_1000000bp_distances.csv` — all pair rows
- `genes_with_lncRNAs_1000000bp.csv` — gene-centric
- `lncRNAs_with_genes_1000000bp.csv` — lncRNA-centric

```
Pairs table:
┌───────────────┬───────────────┬─────────────┬─────────────────────┐
│ gene_name     │ lncRNA_name   │ chrom       │ min_distance_bp     │
├───────────────┼───────────────┼─────────────┼─────────────────────┤
│ HLA-A         │ HCG18         │ chr6        │ 45230               │
│ HLA-A         │ LINC01149     │ chr6        │ 230500              │
│ TAP1          │ HCG18         │ chr6        │ 12300               │
└───────────────┴───────────────┴─────────────┴─────────────────────┘
```

### 8c. lncRNA Interaction Tables (miRNAs + RBPs; ENCORI + POSTAR3)

This module builds compact interaction artifacts for a **deterministic lncRNA set**:
- Tier‑1 panel lncRNAs (`TIER1_LNCRNA_GENES`)
- plus N extra lncRNAs closest to any of the frozen 66 `PRIMARY_GENES` (via §8b pairs table; default N=20)

Canonical build entrypoint: `python -m pipeline.lncRNA_interactions.build_all`
(module doc: `pipeline/md/module_specific_processing_md/LNCRNA_INTERACTIONS_MODULE.md`).

**Output folder**: `data/lncRNA_interactions/`

| Role | On-disk basename | Notes |
|------|------------------|-------|
| Selected lncRNA list | `selected_lncrnas.txt` | one gene symbol per line (reproducibility) |
| ENCORI miRNA→lncRNA interactions | `encori_mirna_targets.parquet` | miRNAs are **arm-level** (e.g. `hsa-miR-339-5p`) + MiRBase `miRNAid` |
| ENCORI fetch diagnostics | `encori_mirna_targets_diagnostics.csv` | per selected lncRNA: rows returned for **symbol** vs **GENCODE gene_id** fallback |
| ENCORI RBP→lncRNA targets (curated subset) | `encori_rbp_targets.parquet` | RBPTarget rows for **top-N POSTAR3 RBPs** × `selected_lncrnas` (not a global `RBP=all,target=all` dump) |
| ENCORI RBP fetch diagnostics | `encori_rbp_targets_diagnostics.csv` | per (RBP, lncRNA): rows returned for **symbol** vs **GENCODE gene_id** fallback |
| POSTAR3 peaks overlapping selected lncRNA loci | `postar3_overlaps.parquet` | peak-level rows (chrom/start/end/rbp/assay/cell_tissue/score) + provenance columns when non-empty |
| POSTAR3 per-RBP summary (across selected lncRNAs) | `postar3_rbp_summary.parquet` | `rbp`, `n_peaks`, `n_assays`, `n_cell_tissue`, `assays`, `cell_tissue` |
| POSTAR3 per-lncRNA locus summary | `postar3_lncrna_summary.parquet` | `gene_name`, `chrom`, `start`, `end`, `n_overlapping_peaks` |
| Convenience shortlist of RBPs | `recommended_rbps_from_postar3.csv` | same as `postar3_rbp_summary` but minimal columns |

**POSTAR3 master peak table (large; upstream of §8c overlaps)**  
- **File**: `data/RBP-RNA/POSTAR3.parquet` (column schema aligned with `scripts/RBP-RNA/build_postar3_parquet.py`: chrom/start/end/peak_id/strand/rbp/assay/cell_tissue/source_accession/score).  
- **Build**: stream `data/RBP-RNA/POSTAR3.txt` → parquet (do not load full TXT in RAM).  
- **Overlap geometry**: `build_all` / `postar3_summary.py` support `--postar3-region-mode` **`gene` | `exons` | `introns` | `promoter`** (promoter window defaults are configurable on the CLI). When overlaps are written, extra columns include **`__selected_set`** (e.g. `tier1_plus_20close`) and **`__region_mode`** (the mode used for that build).

**Predicted miRNA→lncRNA targets (RNAhybrid; local tool path)**  
All paths are under `data/lncRNA_interactions/predicted_targets/rnahybrid/` (default `PATHS.working_dir`):

| Role | Basename | Notes |
|------|----------|--------|
| Exon BED (selected lncRNAs) | `selected_lncrnas.exons.bed` | BED6; chrom names are **normalized** to match the genome FASTA (e.g. `chr1` → `1` for Ensembl primary assembly FASTA). |
| Exon FASTA | `selected_lncrnas.exons.fa` | `bedtools getfasta`; companion stderr log `selected_lncrnas.exons.fa.bedtools.stderr.txt` on failure |
| Spliced transcript FASTA | `selected_lncrnas.spliced_exons.fa` | one sequence per `gene_name` (exons concatenated in transcript order) |
| Truncated targets (optional) | `selected_lncrnas.spliced_exons.3prime_trunc.fa` | last **N** bp per lncRNA (default **N=8000** via `--max-target-bp`; `0` disables). Avoids RNAhybrid “target too long” / aborts on megabase lncRNAs. |
| miRNA subset FASTA | `mirnas.subset.fa` | written when using `--mirna-subset-from-encori` (miRs present in `encori_mirna_targets.parquet`) |
| miRNA subset diagnostics | `mirna_subset_from_encori.csv` | which miRNA names were found in the mature FASTA |
| RNAhybrid stdout | `rnahybrid.raw.txt` | full tool output (includes structure lines under `-c`) |
| Parsed hits table | `rnahybrid.hits.parquet` | one row per **compact** hit line: `target`, `target_hit_start`, `miRNA`, `seed_len`, `mfe`, `raw_line`, plus `stderr`, `exit_code`, `cmd`, **`max_target_bp_truncation`** (string or null) |

CLI entrypoint: `python -m pipeline.lncRNA_interactions.predicted_targets_cli` (see `analysis/COMMANDS.md` and `LNCRNA_INTERACTIONS_MODULE.md`).

### 8d. Unified ChIP peak table (cached input to `chip_hits`)

**File**: `PATHS.chip_unified` → e.g. `data/CHIP/unified_chip_peaks.parquet`  
**Build**: `python scripts/build_unified_chip_peaks.py` (reads all BEDs under `PATHS.chip_dir/ENCODE` and `.../CHIP_ATLAS`).  
**Grain**: one row per peak; includes **`sample_id`** (BED filename stem) so replicate TF×cell-line tracks stay distinguishable (e.g. `CTCF_MCF7_1`). Consumed by regulatory `chip_hits` integration when the parquet exists (`pipeline/main.py` does not auto-rebuild this file).

### 8e. BRCA isoform expression (filtered wide matrix; post-processing)

**Inputs**: pan-cancer isoform TPM TSV (default `data/RNAexp_TCGA/tcga_rsem_isoform_tpm`), BRCA clinical / participant whitelist, GENCODE transcript probemap under `annotations/RNA/`.  
**Scripts**: `pipeline/RNA_exp/filter_isoform_expression_brca.py`, `pipeline/RNA_exp/report_isoform_brca_genes.py`  
**Typical outputs** (paths configurable via CLI):  
- `data/RNAexp_TCGA/isoform_brca_panel.parquet` (or `.tsv`) — BRCA columns only, transcript rows filtered to the analysis gene panel (and optional `--genes-only` mode)  
- `data/RNAexp_TCGA/isoform_brca_panel_qc.json` — chunking / match-rate QC  

---

## 9. RNA Expression Matrix (log2(TPM+1))

**Files**:
- Raw: `PATHS.rna_expression_raw` (e.g. `data/RNAexp_TCGA/TCGA-BRCA.star_tpm_mapped.tsv`)
- Processed: `PATHS.rna_expression` (e.g. `data/RNAexp_TCGA/TCGA-BRCA.star_tpm_processed.tsv`)

**Layout**: Genes (rows) × Samples (columns)

### Header / columns

```
gene_symbol       : str    (e.g. "TSPAN6")
Ensembl_ID        : str    (e.g. "ENSG00000000003.15")
<sample columns>  : float  expression value in **log2(TPM+1)** (unitless)
```

Sample column headers are TCGA sample barcodes (commonly `sample_vial`, e.g. `TCGA-XX-YYYY-01A`).

Notes:
- This pipeline’s RNA “STEP 11” (in `pipeline/main.py`) normalizes gene symbols; it does **not**
  change the expression units. In the provided BRCA matrix, values look like log2(TPM+1).

### Normalized sample metadata companion table

Because the sample IDs live in the **column headers**, we keep the matrix wide and generate a
separate, tidy table with normalized join keys:

- `annotations/_normalized/RNA_expression.sample_metadata.tsv`
  - `raw_sample_id`, `participant`, `sample`, `sample_vial`, `aliquot`, `source`

---

## 10. ATAC Case-Level Sample Matrix (`TCGA_NORMAL_LOG_CPM_QN_BRCA_case_level`)

**File**: `PATHS.atac_case_level_matrix`  
Example: `data/TCGA_ATAC/TCGA_NORMAL_LOG_CPM_QN_BRCA_case_level.csv`  
Column-derived sample list: `PATHS.atac_samples_tsv` (`annotations/ATAC/samples.tsv`, from `scripts/annotations/build_annotation_sample_manifests.py`).

**Layout**: Peaks (rows) × Samples (columns), with peak coordinates in leading columns.

### Leading columns (per peak)

```
seqnames   : str   (chrom, e.g. "chr1")
start      : int
end        : int
name       : str   (peak label)
score      : float
annotation : str   (e.g. "Promoter")
GC         : float
```

### Sample columns (units)

After the leading peak columns, each remaining column is a sample.
In this file, sample column headers often look like aliquot-like barcodes:
`TCGA-XX-YYYY-01A-..-..-..`

Values are **log CPM** (counts-per-million), **quantile normalized** (per the file naming and ATAC module docs).
This is a per-peak per-sample accessibility measure (unitless after log transform).

### Normalized sample metadata companion table

- `annotations/_normalized/ATAC_case_level.sample_metadata.tsv`
  - `raw_sample_id`, `participant`, `sample`, `sample_vial`, `aliquot`, `source`

---

## 11. Unified Clinical + Immune Table (BRCA)

**File**: `annotations/BRCA_clinical_immune_unified.tsv`  
**Granularity**: One row per TCGA case-level sample (`sample` like `TCGA-...-01`)

**Units**:
- Clinical stage/T/N/M: categorical strings
- `CPE`: unitless score (0–1 in observed data)

### Columns

```
sample_id         : str  (e.g. "TCGA-3C-AALJ-01")
pathologic_stage  : str
pathologic_T      : str
pathologic_N      : str
pathologic_M      : str
PAM50_final       : str
CPE              : float

# normalized join keys (added by pipeline)
participant       : str
sample            : str
sample_vial       : str | null
aliquot           : str | null
```

---

## 11a. Thorsson Immune Table (PanCancer; C1–C6)  *(input annotation)*

**Files**:
- `annotations/Thornsson_immune_table.tsv`  *(raw)*
- `annotations/_normalized/Thornsson_immune_table.normalized.tsv` *(same content with normalized join keys when available)*
- Related subsets: `annotations/TCGA_immune_subtypes.tsv`, `annotations/BRCA_immune_subtypes.tsv`,
  `annotations/BRCA_immune_subtypes_advanced.tsv`

**Granularity**: One row per TCGA participant/sample row from Thorsson et al. immune subtypes.

**Key concept**: the column **`Immune Subtype`** is the canonical **C1–C6** immune subtype label.

**Recommended usage (Tier 4 deferral)**: treat this table as the canonical source for **immune composition / stromal context** (e.g. CIBERSORT-like cell fractions, `Leukocyte Fraction`, `Stromal Fraction`, `TIL Regional Fraction`). The gene panel’s `TIER4_READOUT_GENES` is intentionally compact (lymphoid/cytolytic/checkpoint readouts) and does **not** attempt to enumerate myeloid/stromal markers exhaustively.

### Minimal join keys / identifiers

```
TCGA Participant Barcode  : str   (participant / case; e.g. "TCGA-XX-YYYY")
TCGA Study                : str   (e.g. "BRCA")
Immune Subtype            : str   (C1..C6)
TCGA Subtype              : str   (study-specific; e.g. PAM50-ish / TCGA subtype when present)
```

### Commonly used score blocks (high-signal for QC / sanity)

These columns are used frequently for quick biological orientation and should be **non-empty** in typical TCGA cohorts:

- **Immune / stromal fractions**
  - `Leukocyte Fraction`, `Stromal Fraction`, `Intratumor Heterogeneity`
  - `Lymphocyte Infiltration Signature Score`, `TIL Regional Fraction` (when present)
- **Immune program scores**
  - `IFN-gamma Response`, `TGF-beta Response`, `Wound Healing`, `Proliferation`, `Macrophage Regulation`
- **Mutation / neoantigen burden**
  - `SNV Neoantigens`, `Indel Neoantigens`, `Silent Mutation Rate`, `Nonsilent Mutation Rate`
- **CNV / genome instability**
  - `Number of Segments`, `Fraction Altered`, `Aneuploidy Score`, `Homologous Recombination Defects`
- **Receptor diversity**
  - `BCR Evenness`, `BCR Richness`, `BCR Shannon` (and `TCR *` analogues when present)
- **Cell-type fractions (CIBERSORT-like)**
  - A large family of columns with names like:
    `B Cells Memory`, `B Cells Naive`, `T Cells CD8`, `T Cells CD4 Memory Activated`, `NK Cells Activated`,
    `Macrophages M0/M1/M2`, `Dendritic Cells Activated/Resting`, `Mast Cells Activated/Resting`, `Neutrophils`, `Monocytes`, …

### Notes on width (+700 columns)

This table is intentionally **wide** (hundreds of columns). For analysis/QC, treat it as a set of *named blocks*:
- Clinical / subtype labels (stage, PAM50, ER/PR/HER2 fields)
- Global burden scores (mutation/CNV/aneuploidy/HRD)
- Immune program scores (IFNG, TGFb, wound healing, proliferation, etc.)
- Immune cell fractions + receptor diversity metrics

If you need a programmatic column inventory (grouped by block), generate it from the header rather than
hand-maintaining 700+ names. (The sanity suite under `analysis/sanity/` is the intended place for that.)

---

## 12. miRTarBase Output Tables

**Source**: `PATHS.mirtarbase_csv` (`data/miRNA/mirtar.csv`) + optional `PATHS.mir_family_info`  
**Module**: `pipeline/genes/mirtarbase.py`  
**Output directory**: `{PATHS.working_dir}/miRNA/mirtarbase/` (see `mirtarbase_module.md`)

**Units / types**:
- Binary flags: **0/1**
- Counts (`n_*`): **integer study or row counts**
- `evidence_score`: **integer** weighted sum (see below)
- JSON columns: **JSON string** in CSV (use `json.loads`)

### Support type slugs (column name suffix)

```
functional_mti          ← "Functional MTI"
functional_mti_weak     ← "Functional MTI (Weak)"
nonfunctional_mti       ← "Non-Functional MTI"
nonfunctional_mti_weak  ← "Non-Functional MTI (Weak)"
```

### Experiment classes (six)

```
reporter | binding | protein | rna | perturbation | other
```

Cross columns use `has_<experiment>__<support_slug>` and `n_<experiment>__<support_slug>_studies`.

---

### 12a. `mirtar_interaction_study_collapsed.csv`

**Granularity**: one row per **(miRNA, gene, study_id)**

```
┌────────────────────┬──────────┬─────────────────────────────────────────┐
│ Column             │ Type     │ Notes                                   │
├────────────────────┼──────────┼─────────────────────────────────────────┤
│ miRNA              │ str      │ e.g. hsa-miR-155-5p                     │
│ miRNA_family       │ str      │ TargetScan family or name proxy       │
│ gene               │ str      │ gene symbol                             │
│ study_id           │ str      │ PMID or REF:: / MIRTAR:: pseudo-id      │
│ has_functional_mti │ int 0/1  │                                         │
│ has_functional_mti_weak │ int 0/1 │                                  │
│ has_nonfunctional_mti   │ int 0/1 │                                  │
│ has_nonfunctional_mti_weak │ int 0/1 │                               │
│ has_reporter … has_other │ int 0/1 │ experiment class flags         │
│ has_<exp>__<support_slug> │ int 0/1 │ 24 cross flags (6×4)          │
└────────────────────┴──────────┴─────────────────────────────────────────┘
```

---

### 12b. `mirtar_interaction_summary.csv`

**Granularity**: one row per **(miRNA, gene)** — aggregated over studies

```
┌──────────────────────────────┬──────────┬──────────────────────────────┐
│ Column                       │ Type     │ Notes                        │
├──────────────────────────────┼──────────┼──────────────────────────────┤
│ miRNA, miRNA_family, gene    │ str      │                              │
│ n_studies                    │ int      │ unique studies               │
│ n_functional_mti_studies     │ int      │ … + weak / nonfunctional …   │
│ n_<experiment>_studies       │ int      │ six experiment totals       │
│ n_<exp>__<support>_studies    │ int      │ 24 cross-count columns       │
│ experiment_support_counts_json │ str │ nested dict as JSON string    │
│ evidence_score               │ int      │ weighted evidence (below)    │
└──────────────────────────────┴──────────┴──────────────────────────────┘
```

**`evidence_score`** (integer, per miRNA×gene):

```
evidence_score = 3×Σ(reporter studies)
               + 3×Σ(binding studies)
               + 2×Σ(protein studies)
               + 1×Σ(rna studies)
               + 1×Σ(perturbation studies)
```

**`experiment_support_counts_json`** (conceptual tree):

```
experiment_support_counts_json ── (JSON string)
│
└── "<experiment>" ─────────────────────────────────────────────────────
    ├── "Functional MTI"            : int
    ├── "Functional MTI (Weak)"    : int
    ├── "Non-Functional MTI"       : int
    └── "Non-Functional MTI (Weak)" : int
```

---

### 12c. `mirtar_gene_summary.csv`

**Granularity**: one row per **gene** (panel)

**Scalar columns** (representative):

```
gene
n_unique_miRNAs
n_unique_families
n_total_studies
n_<support_slug>_studies          (4 support types)
n_<experiment>_studies            (6 experiment classes)
n_<experiment>__<support_slug>_studies  (24 columns)
```

**JSON columns** (each cell is a **JSON string**):

```
experiment_support_counts_json ── same nested shape as §12b

mirna_study_counts_json ──────────
│
└── "<miRNA>" : int   (studies targeting this gene via that miRNA)

family_study_counts_json ─────────
│
└── "<family>" : int

support_to_mirnas_json ───────────
│
└── "<Support Type>" : list[str]   (miRNA names)

experiment_to_mirnas_json ────────
│
└── "<experiment>" : list[str]

experiment_support_to_mirnas_json
│
└── "<experiment>"
    └── "<Support Type>" : list[str]

support_to_families_json / experiment_to_families_json / experiment_support_to_families_json
    (same pattern at family level)

mirna_detail_json ─────────────────────────────────────────────────────
│
└── "<miRNA>"
    ├── "n_studies" : int
    ├── "support_counts" : { "<Support Type>": int, ... }
    └── "experiment_support_counts" : { "<experiment>": { "<Support Type>": int, ... } }
```

---

### 12d. `mirtar_mirna_summary.csv`

**Granularity**: one row per **miRNA**

Mirror of gene summary: scalar counts + JSON dicts (`target_gene_study_counts_json`, `target_gene_detail_json`, etc.) — same nesting ideas as §12c but gene-centric keys.

---

### 12e. `mirtar_family_summary.csv`

**Granularity**: one row per **miRNA family**

```
miRNA_family
n_unique_miRNAs
n_unique_targets
n_total_studies
n_<support_slug>_studies
n_<experiment>_studies
```

---

## Cross-Reference: Shared Structure Patterns

### Gene symbol normalization & alias registry (pipeline-wide)

**Why this exists**: external tables mix symbol conventions (HGNC aliases, legacy names, UCSC quirks, Ensembl IDs).
The pipeline normalizes these into the **canonical symbols used in `pipeline/config.py`** so joins are stable across
RNA/SNV/SV/CNV/methylation/miRNA modules.

**Canonical implementation**:
- **Registry**: `pipeline/genes/panel_alias_registry.py`
  - `get_gene_flat_alias_map(max_lines=None)` → **`alias → canonical`**
  - `get_gene_canonical_to_aliases(max_lines=None)` → **`canonical → {aliases}`**
- **Entry-point helpers**: `pipeline/genes/symbol_normalization.py`
  - `default_symbol_mapping()` (full registry; slow once, cached)
  - `get_panel_symbol_mapping(max_lines)` (optionally capped)

**Alias sources (merged)**:
- **HGNC alias table**: `annotations/Alias_v5.22.xls` (streamed; header is tab-separated, body is CSV)
- **UCSC fixes**: `UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE` (`pipeline/config.py`)
- **Legacy remaps**: `LEGACY_DATASET_SYMBOL_RENAMES` (`pipeline/config.py`; e.g. `TMEM173→STING1`, `MB21D1→CGAS`)
- **Hand-curated extras**: `MANUAL_CANONICAL_TO_ALIASES` (`panel_alias_registry.py`)
- **Ensembl IDs**: `ENSG…` (+ versionless) → canonical `gene_name` from `PATHS.gencode_gtf_pq`

**Performance knob (important)**:
- `APM_HGNC_ALIAS_MAX_LINES=<N>` limits how much of the HGNC alias table is scanned (dev speed).
  This can **miss aliases** that appear later in the file (e.g. tokens like `NF-YA` or `p65` may not map if truncated).
  For correctness-sensitive runs, prefer leaving it unset so the registry performs a full scan (cached).

**Disable mapping**:
- `APM_USE_GENE_SYMBOL_MAPPING=0` disables HGNC/UCSC/legacy remaps in the major loaders that support it.

### Overlap info dict (used in ATAC, SV, CNV, utils)

```
overlap ──────────────────────────────────────────
├── "overlaps"            : bool
├── "overlap_bp"          : int
├── "overlap_interval"    : [start, end] | null
├── "overlap_frac_of_a"   : float   (fraction of first interval)
└── "overlap_frac_of_b"   : float   (fraction of second interval)
```

### TAD_domains dict (used in elem_focus, genes, ATAC peaks)

```
(See Section 1 TAD_domains for full tree — identical structure
 everywhere it appears)
```

### Spatial hit dict pattern (used in SV gene_hits, SV elem_hits, CNV)

```
signed_dist + overlap_start + overlap_end + overlap_bp
  + overlap_percent + region classification flags
  + region_hit string + hit_side
```

### VEP consequence hit pattern (used in SNV, SV)

```
Allele + Consequence + IMPACT + SYMBOL + Gene + Feature_type
  + Feature + BIOTYPE + CANONICAL + MANE_SELECT + SIFT + PolyPhen
  (gene_hits variant)

or

Allele + Consequence + Feature_type + Feature + MOTIF_NAME
  + MOTIF_SCORE_CHANGE + TRANSCRIPTION_FACTORS
  (motif_hits variant)
```

