# Tests index (what exists and where)

This repo mixes **unit tests** (fast, deterministic) and **smoke / integration tests** (use real data under `data/`).

If you are running commands from Windows-hosted shells, prefer the WSL bridge pattern in `.cursor/rules/wsl-shell-bridge.mdc`.

---

## Fast unit tests (no big data required)

### QC tests and thresholds (definition doc)

See `analysis/qc/QC_TESTS_AND_THRESHOLDS.md`.

### Gene panel doc sync

- **Test**: `tests/test_docs_sync.py`  
  Verifies `research_plan/01_gene_panel_extended.md` generated blocks and gene-symbol mentions stay in sync with `pipeline/config.py`.
- **Regenerator**: `scripts/panel/regen_gene_panel_doc_fragments.py`

### Atlas + scanning contracts

- `tests/test_atlas_generated_blocks.py` — RPPA panel name list in `DATA_STRUCTURE_ATLAS.md` matches `empty_rppa_panel_scores()`.
- `tests/test_scanning_columns_spec.py` — `scan_gene_links_n_genes` / `any_panel_gene` match `gene_links` dict semantics.
- `tests/test_sv_span_gene_hit_schema_keys.py` — stable key set for SV span `classify_span_sv_gene_hit` output.
- `tests/test_hla_loh_ascat_segments.py` — segment-based HLA LOH (minor copy = 0) on synthetic overlaps.
- `tests/test_hla_snv_vep.py` — HLA-focused SNV summaries from synthetic ``gene_hits`` (LoF vs high missense).

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m unittest -q tests.test_docs_sync
.venv/bin/python3 scripts/panel/regen_gene_panel_doc_fragments.py --check
```

### In-code schema snapshot (builders in `pipeline/schemas.py`)

- **Generator**: `scripts/schema/dump_pipeline_schema_snapshot.py`
- **Fixture**: `tests/fixtures/pipeline_schema_snapshot.json`
- **Test**: `tests/test_schema_snapshot.py`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 scripts/schema/dump_pipeline_schema_snapshot.py --check
.venv/bin/python3 -m unittest -q tests.test_schema_snapshot
```

---

## Output-schema snapshots (use actual on-disk outputs)

These tests read a small set of columns from existing parquet outputs and snapshot the **nested key structure** of JSON-encoded columns.

- **Generator**: `scripts/schema/dump_output_schema_snapshot.py`
- **Fixture**: `tests/fixtures/output_schema_snapshot.json`
- **Test**: `tests/test_output_schema_snapshot.py`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 scripts/schema/dump_output_schema_snapshot.py --check
.venv/bin/python3 -m unittest -q tests.test_output_schema_snapshot
```

Notes:
- Requires the element-focus parquet at `pipeline.config.PATHS.regulatory_elements_table_with_evidence_parquet`.
- Optionally snapshots ATAC peak outputs if `data/atac_peaks/atac_peaks_annotated.parquet` exists.

---

## Smoke tests (real data; slower)

### Pipeline smoke suite (multi-module)

- **Runner**: `tests/run_smoke_tests.py`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m tests.run_smoke_tests --list
.venv/bin/python3 -m tests.run_smoke_tests
```

### Selected unit tests

There are additional unit tests for specific modules (examples):

- `tests/test_panel_alias_registry.py`
- `tests/test_snv_fimo.py`
- `tests/test_snv_chip.py`
- `tests/test_chip_atlas_loader.py`
- `tests/test_covariates_module.py`

Run all unit tests (excluding smoke runner):

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m unittest -q
```

---

## Cross-modal QC (real data; cohort run)

- **Runner**: `analysis/qc/run_cross_modal_qc.py`
- **Inputs**:
  - `--scratch-json` from `analysis/trial_max_coverage_cohort/.../cohort_processing_outputs.json`
  - `annotations/Thornsson_immune_table.tsv`
  - `data/Methylation/cohort/gene_meth_matrix.csv` (via scratch context)
  - `data/CNV_TCGA/CNV_gene_tables/*_cnv_gene_calls_*.csv` (via scratch context)
  - element focus parquet at `pipeline.config.PATHS.regulatory_elements_table_with_evidence_parquet`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/qc/run_cross_modal_qc.py \
  --scratch-json analysis/trial_max_coverage_cohort/processing_scratch_20260421_192526/cohort_processing_outputs.json
```

---

## Full QC (pipeline/qc modules; pass/warn/fail findings)

- **Runner**: `analysis/qc/run_full_qc.py`
- **Design**: compute()+assert_cohort() per module under `pipeline/qc/`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/qc/run_full_qc.py \
  --scratch-json analysis/trial_max_coverage_cohort/processing_scratch_20260421_192526/cohort_processing_outputs.json
```

---

## QC suite (one command)

- **Runner**: `analysis/qc/run_qc_suite.py`
- Builds cohort RPPA outputs (so `panel_scores.csv` exists under `ctx.rppa_output_dir`), then runs:
  - `analysis/qc/run_cross_modal_qc.py`
  - `analysis/qc/run_full_qc.py`

Run:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/qc/run_qc_suite.py \
  --scratch-json analysis/trial_max_coverage_cohort/processing_scratch_20260421_192526/cohort_processing_outputs.json
```

