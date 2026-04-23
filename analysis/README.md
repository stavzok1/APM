# Analysis utilities

Scripts here are **downstream of the pipeline**: they read `pipeline.config.PATHS` and annotation tables under `annotations/`, then write organized outputs under `analysis/sample_coverage/output/`.

## Sample × module coverage (see `analysis/sample_coverage/`)

**Script:** `sample_module_coverage.py`

**Prerequisite (RNA / ATAC / HLA / HiCHIP / miRNA):** generate column-derived manifests once (or after matrix updates):

```bash
python scripts/annotations/build_annotation_sample_manifests.py
```

This writes `annotations/RNA/samples.tsv`, `annotations/ATAC/samples.tsv`, `annotations/HLA/samples.tsv`, `annotations/HiCHIP/samples.tsv`, and `annotations/miRNA/samples.tsv` (paths on `PATHS.*_samples_tsv`). The coverage script prefers these files, then falls back to `annotations/_normalized/*.tsv`, then scans wide matrix headers directly.

**What it does**

1. Loads sample identifiers from each omics / annotation source (SNV, SV, CNV, methylation, RPPA, **miRNA**, RNA, ATAC, HLA, **HiCHIP**; plus immune tables and unified clinical metadata).
2. Normalizes IDs using `pipeline.sample_ids.normalize_tcga_id` (same rules as `SAMPLE_ID_MATCHING_GUIDE.md`).
3. For GDC-style manifests (`samples.tsv`), picks the **tumor** aliquot using `Tissue Type` when comma-aligned with `Sample ID`; otherwise falls back to the same rule as `pipeline/SNV/vcf_loader.py` (first vs second token when `Tumor Descriptor` starts with `Not`).
4. Writes:
   - **`sample_vial_presence.tsv`**: one row per tumor `sample_vial`, boolean columns per assay keyed at vial level.
   - **`participant_presence.tsv`**: full participant matrix (omics + metadata columns).
   - **`participant_presence_omics.tsv`** / **`participant_presence_metadata.tsv`**: same rows, split columns.
   - **`module_counts_omics.tsv`** / **`module_counts_metadata.tsv`**: counts split by omics vs metadata; **`module_counts.tsv`** keeps the combined legacy layout.
   - **`intersections_sample_vial.tsv`**: omics-only multi-way overlaps at vial level.
   - **`intersections_participant_omics.tsv`** / **`intersections_participant_metadata.tsv`**: participant intersections split.
   - **`pairwise_overlap_sample_vial.tsv`**, **`pairwise_overlap_participant_omics.tsv`**, **`pairwise_overlap_participant_metadata.tsv`**
   - **`module_notes.tsv`**: non–patient-local reference modules (if any).
5. Saves **figures** (requires `matplotlib`), in two bundles with filenames `*_omics_<run_key>.png` and `*_metadata_<run_key>.png` (vial-level plots only in the omics bundle).
6. **Clinical stratification** (after a coverage run): `clinical_omics_stratification.py --presence <run>/participant_presence_omics.tsv` merges with `BRCA_clinical_immune_unified.tsv` (needs `pathologic_stage_collapsed`; run `scripts/annotations/update_brca_clinical_unified_stage.py` once) and writes PAM50 / stage tables plus heatmaps under a `stratification_*` folder.

**HLA cohort filter:** `scripts/annotations/filter_hla_to_brca_clinical_cohort.py` archives the full comma table to `annotations/HLA/HLA_types.unfiltered.comma.csv` (once) and rewrites `annotations/HLA_types.tsv` to participants present in `PATHS.brca_clinical`. Then rerun `scripts/annotations/normalize_external_annotations.py` and `scripts/annotations/build_annotation_sample_manifests.py` so `_normalized` and `HLA/samples.tsv` stay in sync.

**One-shot (WSL):** from repo root, `bash scripts/run_full_annotation_coverage_pipeline.sh` runs filter → unified stage column → normalize → manifests → coverage → stratification under the latest `analysis/sample_coverage/output/run_*`.

**Run** (from repository root):

```bash
python analysis/sample_module_coverage.py
```

Optional custom output directory:

```bash
python analysis/sample_module_coverage.py --out analysis/output/my_run
```

**Interpretation**

- Use **`sample_vial`** tables to compare SNV/SV/CNV/methylation/RPPA/RNA/ATAC/HLA (omics-style IDs).
- Use **`participant`** tables when mixing in **immune** annotations that are indexed at `TCGA-…-01` style or participant barcodes.
- HiChIP in the regulatory pipeline is **not** a TCGA per-patient assay here; see `module_notes.tsv`.
