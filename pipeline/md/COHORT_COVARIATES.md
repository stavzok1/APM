# Cohort covariates (post-processing module)

This repo intentionally separates:

- **Element-/gene-centric artifacts** produced by `pipeline/main.py` (cCRE table + `gene_links` evidence), from
- **Sample-centric covariates** computed after a cohort is assembled.

The cohort covariates builder lives under `pipeline/covariates/` and is **not** called from `pipeline/main.py`.

## Output

By default the builder writes under:

- `data/covariates/<run_id>/cohort_covariates.parquet` (canonical)
- `data/covariates/<run_id>/cohort_covariates.csv` (optional)
- `data/covariates/<run_id>/metadata.json` (provenance + provider diagnostics)

## Join key policy

The covariates table index is **`sample_vial`** (`TCGA-XX-YYYY-01A` when available).
Normalized columns are included to support joining to sources that only have:

- `sample` (`TCGA-XX-YYYY-01`)
- `participant` (`TCGA-XX-YYYY`)
- `aliquot` (full barcode)

Some sources in this repo are naturally at **sample-level** (clinical/immune unified, RNA-wide TPM). When those are added, they are **lifted** onto `sample_vial` by matching on the derived `sample` field, so the same sample-level value may appear on multiple vials if present in the cohort.

## What is included (default providers)

- **Coverage/presence**: reads `analysis/sample_coverage/output/current/tables/omics/sample_vial_presence.tsv` when present (written with `python analysis/sample_module_coverage.py --write-current`); otherwise falls back to the latest `analysis/sample_coverage/output/run_*/tables/omics/sample_vial_presence.tsv` (already `sample_vial`-keyed).
- **Clinical/immune**: reads `annotations/BRCA_clinical_immune_unified.tsv` (sample-keyed) and lifts to `sample_vial`.
- **RPPA**: reads `data/rppa/processed/rppa_analysis_combined.parquet` (panel scores, blocks).
- **RNA signatures**: reads `data/RNAexp_TCGA/...processed.tsv` and computes gene-set scores (see `pipeline/RNA_exp/signatures.py`).

Optional add-ons:

- **DDR/HRD score table** (`--ddr-scores`): `DDRscores.txt`-style table with HRD metrics, purity/ploidy, mutSig3, TP53 score, etc. (BRCA-only filtered by acronym when present).
- **Mutation MAF** (`--maf`): derives TP53 and PIK3CA hotspot covariates from a MAF-like table.
- **External covariate tables** (`--external name=path[:sep][:id_col]`): e.g., HRD scores exported elsewhere.

## Signatures (RNA)

The post-processing module uses `pipeline/RNA_exp/signatures.default_gene_sets()` which includes:

- Basic immune visibility readouts (APM class I, CD8, NK, IFNG axis, CYT)
- **Hypoxia**: Buffa hypoxia metagene (`HYPOXIA_BUFFA_mean`) scored as mean
- **Autophagy (proxy)**: `AUTOPHAGY_core_mean` (transcriptional state proxy, not flux)
- **MHC-II readout**: `MHCII_readout_mean` over `HLA-DRA`, `CD74`, `CIITA` (when present)

## CLI

Run from repo root:

```bash
.venv/bin/python3 -m pipeline.covariates.build_covariates --run-id run_cov_001
```

Add a MAF and an external HRD table:

```bash
.venv/bin/python3 -m pipeline.covariates.build_covariates \\
  --run-id run_cov_002 \\
  --maf /path/to/somatic.maf \\
  --external hrd=/path/to/hrd_scores.tsv:\\t:sample_id
```

Use the canonical DDRscores table in `data/covariates`:

```bash
.venv/bin/python3 -m pipeline.covariates.build_covariates \\
  --run-id run_cov_003 \\
  --ddr-scores data/covariates/DDRscores.txt
```

