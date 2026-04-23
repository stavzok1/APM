## Sanity checks / QC (quick orientation)

This folder contains **repeatable sanity tests** for:
- **Structural QC**: suspicious NaN/empty columns, nested payloads that look undecoded, per-sample file presence.
- **Biological QC**: "does the world look sane" checks across modalities (RNA, CNV, methylation, ATAC-derived methylation).

### Inputs

These scripts are designed to work on a cohort subset processed via
`analysis/trial_max_coverage_cohort/run_cohort_module_processing.py`.

They accept either:
- a **scratch** folder containing `cohort_processing_outputs.json`, or
- a methylation subset manifest (TSV) listing sample IDs.

### Run (WSL)

Structural QC:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m analysis.sanity.structural_qc \
  --scratch-json analysis/trial_max_coverage_cohort/processing_scratch_20260420_123616/cohort_processing_outputs.json
```

Produces (under `analysis/sanity/output/structural_qc_<timestamp>/`):
- `sample_ids.txt`: sample IDs used (from the methylation subset manifest in the scratch json)
- `probe_reference_subset__suspicious_columns.csv`, `probe_reference_subset__nested_cells.csv`
- `methylation_*__suspicious_columns.csv` + `methylation_*__nested_cells.csv` for cohort matrices
- `methylation_per_sample_presence_and_qc.csv`
- `cnv_gene_calls_presence.csv`
- Per-module scans:
  - `snv__per_file_qc.csv`, `snv__aggregated_suspicious_columns.csv`
  - `sv_final__per_file_qc.csv`, `sv_final__aggregated_suspicious_columns.csv` (plus `sv_neojunction`, `sv_chip_enriched`, `sv_processed` when present)
  - `cnv_annotated__per_file_qc.csv`, `cnv_annotated__aggregated_suspicious_columns.csv` (annotated segments)
  - `rppa__per_file_qc.csv`, `rppa__aggregated_suspicious_columns.csv`

Biological orientation:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m analysis.sanity.biological_qc \
  --scratch-json analysis/trial_max_coverage_cohort/processing_scratch_20260420_123616/cohort_processing_outputs.json
```

Produces (under `analysis/sanity/output/biological_qc_<timestamp>/`):
- `rna_signatures.csv`: CYT, CD8/NK/IFNG means, APM class I mean (computed from TPM)
- `rna_apm_key_genes_tpm.csv`: TPM for HLA-A/B/C, B2M, TAP1/2, PSMB8/9, NLRC5
- `methylation_promoter_beta_key_genes.csv`: promoter beta means for those genes (from methylation cohort matrix)
- `cnv_log2_key_genes.csv`: per-gene CNV log2 values (from per-sample CNV gene tables when present)
- `corr_cnv_vs_expression_key_genes.csv`: per-gene CNV↔expression correlation across the subset
- `corr_promoter_beta_vs_expression_key_genes.csv`: per-gene promoter beta↔expression correlation across the subset

### Thorsson immune table column inventory

This is an **annotation inventory** tool (not a pipeline-output sanity check), so it lives under `scripts/`:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 scripts/coverage/profile_thorsson_immune_table.py \
  --input annotations/Thornsson_immune_table.tsv
```

By default it writes:

`pipeline/md/contract_tables/Thorsson_immune_table.profile.md`

### Outputs

Each script writes a timestamped folder under:

`analysis/sanity/output/<run_id>/`

with CSV summaries you can inspect or plot in notebooks.

