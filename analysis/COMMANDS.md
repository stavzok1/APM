## APM analysis commands

This file is a “copy/paste” index of common commands you’ll run while iterating on cohort
processing + QC from the repo root.

### Build a max-coverage cohort (50 tumor vials)

```bash
cd /home/stavz/masters/gdc/APM
OUT=analysis/trial_max_coverage_cohort/run_$(date -u +%Y%m%d_%H%M%S)
.venv/bin/python3 analysis/trial_max_coverage_cohort/build_trial_cohort.py --n 50 --out "$OUT"
```

This writes (among others): `"$OUT/cohort_summary.tsv"`.

### Run cohort module processing (all modules; canonical outputs)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv"
```

### Run a single module

```bash
cd /home/stavz/masters/gdc/APM

# SNV only
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --snv-only

# SV only
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --sv-only

# CNV only
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --cnv-only

# Methylation only
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --methylation-only

# RPPA only
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --rppa-only
```

Only one `--*-only` flag can be set per invocation.

### Methylation: reuse probe reference

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --methylation-only --methylation-reuse-probe-reference
```

### SV: control VEP & motifs

```bash
# Run SV with VEP enabled (slower)
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --sv-only --no-sv-skip-vep

# Resume motif work from existing FIMO (reuse 02/03/05; rebuild 06/07/08/09)
.venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \
  --cohort-summary "$OUT/cohort_summary.tsv" \
  --sv-only --sv-resume-motifs-from-fimo
```

### SNV: recreate VEP VCFs (`vcfs_snv` → `vep_vcfs`)

Raw Mutect2 VCFs live under:
- `data/SNV/vcfs_snv` (or legacy `data/SNV/vcfs_SNV`)

VEP outputs are written to:
- `data/SNV/vep_vcfs`

The canonical script in this repo:
- `scripts/SNV_scripts/run_vep_apm_region.sh`

```bash
cd /home/stavz/masters/gdc/APM

# This script:
# 1) builds `data/SNV/apm_genes_1Mb.bed` if missing (from pipeline gene panel)
# 2) bedtools-intersects each VCF to that ±1Mb bed
# 3) runs VEP to produce: *.APM_1Mb.vep.vcf under data/SNV/vep_vcfs/
bash scripts/SNV_scripts/run_vep_apm_region.sh
```

Notes:
- Requires `bedtools` and `vep` available on PATH, and a VEP cache + FASTA at `~/.vep/`.
- The intermediate `*.APM_1Mb.vcf` is **pre-VEP** (no `CSQ`). The pipeline expects `*.APM_1Mb.vep.vcf` when present.

### Structural QC (schema/empties/nested-cells)

This reads the canonical outputs recorded in `cohort_processing_outputs.json` and produces a new
run folder under `analysis/sanity/output/`.

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/sanity/structural_qc.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"
```

### Biological QC (orientation sanity checks)

This uses the same `--scratch-json` pointer and writes a new run folder under `analysis/sanity/output/`.

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/sanity/biological_qc.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"
```

### Cohort QC report (analysis module; markdown output)

Writes `analysis/qc/output/<run_id>/report.md` plus supporting CSVs.

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/qc/run_cohort_qc_report.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"
```

### QC suite (recommended one-command runner)

Builds cohort RPPA outputs (if needed) and runs both QC reports:
- `analysis/qc/run_cross_modal_qc.py`
- `analysis/qc/run_full_qc.py`

Writes new run folders under `analysis/qc/output/`.

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 analysis/qc/run_qc_suite.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"
```

### QC components (run individually)

```bash
cd /home/stavz/masters/gdc/APM

# Build cohort RPPA outputs (panel_scores, blocks, etc.)
.venv/bin/python3 analysis/qc/build_rppa_outputs_for_cohort.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"

# Cross-modal QC report (Tier 1 structural + cross-modal correlations)
.venv/bin/python3 analysis/qc/run_cross_modal_qc.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"

# Full QC report (pipeline/qc modules; pass/warn/fail findings)
.venv/bin/python3 analysis/qc/run_full_qc.py \
  --scratch-json "analysis/trial_max_coverage_cohort/processing_scratch_<UTC>/cohort_processing_outputs.json"
```

### TAD boundaries enriched table (beyond CTCF-only boundary hits)

Builds a single parquet across biosamples with:
- overlapping cCRE IDs
- overlapping ATAC peak IDs
- overlapping methylation probe IDs
- adjacent-domain gene lists (left/right)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m pipeline.tad_annotation.boundary_enrich_cli \
  --processed-dir "$(python3 -c 'from pipeline.config import PATHS; print(PATHS.tads_processed)')"
```

### Debugging input linking (why a sample didn’t run)

Every cohort run writes a scratch folder containing:
- `cohort_processing_prepare.log`: per-sample MISS diagnostics + a JSON dict of what linked
- `cohort_processing_outputs.json`: canonical output paths (used by sanity scripts)

Open the prepare log and search for `MISS cnv` / `MISS sv` / `MISS snv` / `MISS rppa`.

### ChIP (ChIP-Atlas + ENCODE → unified peaks)

From repo root (`PATHS.chip_unified` defaults to `data/CHIP/unified_chip_peaks.parquet`).

```bash
cd /home/stavz/masters/gdc/APM

# ChIP-Atlas UCSC / gffTags multi-column BEDs → narrow 6-column (in-place; *.pre_narrow.bak backups)
.venv/bin/python3 scripts/chip/preprocess_chip_atlas_beds.py --dry-run
.venv/bin/python3 scripts/chip/preprocess_chip_atlas_beds.py

# Merge ENCODE + ChIP-Atlas into one parquet (run after preprocess if Atlas was UCSC-style)
.venv/bin/python3 scripts/chip/build_unified_chip_peaks.py

# TF × source: how many experiments (distinct sample_id / BED stem) per TF
.venv/bin/python3 scripts/chip/summarize_chip_tf_experiments.py
.venv/bin/python3 scripts/chip/summarize_chip_tf_experiments.py --list-samples
```

### Isoform expression (BRCA column whitelist + panel transcripts)

```bash
cd /home/stavz/masters/gdc/APM

# Full filter (default paths under data/RNAexp_TCGA/ …)
.venv/bin/python3 pipeline/RNA_exp/filter_isoform_expression_brca.py \
  --out-table data/RNAexp_TCGA/isoform_tpm_brca_panel.parquet \
  --qc-json data/RNAexp_TCGA/isoform_tpm_brca_panel_qc.json

# Protein-coding / primary-gene universe only (no lncRNA distance or lnc GENCODE features)
.venv/bin/python3 pipeline/RNA_exp/filter_isoform_expression_brca.py --genes-only \
  --out-table data/RNAexp_TCGA/isoform_tpm_brca_genes_only.parquet \
  --qc-json data/RNAexp_TCGA/isoform_tpm_brca_genes_only_qc.json

# Retention / genes / MANE / tags report (genes-only rules; add --cohort-summary to match a trial cohort)
.venv/bin/python3 pipeline/RNA_exp/report_isoform_brca_genes.py --out-json /tmp/isoform_brca_genes_report.json

# Quick scan (first N chunks only)
.venv/bin/python3 pipeline/RNA_exp/report_isoform_brca_genes.py --max-chunks 5 --top-genes 20
```

### lncRNA interactions (ENCORI / starBase API)

Writes cached raw TSVs under `data/external_cache/encori/…` and the aggregated parquets under
`data/lncRNA_interactions/`.

Notes:
- ENCORI **RBP-Target** supports only RBPs that have CLIP datasets in ENCORI (classic RBPs like `FUS`,
  `HNRNPC`, `ELAVL1`, `AGO2`, etc.). Many chromatin modifiers (e.g. `EZH2`, `BRD4`, `EP300`) are **not**
  available as “RBP” in ENCORI.
- miRNA targets are restricted to: the 8 Tier-1 panel lncRNAs + 20 additional closest-to-primary-gene
  lncRNAs from `data/lncRNA_matching/genes_lncRNAs_1000000bp_distances.csv` (deterministic).

```bash
cd /home/stavz/masters/gdc/APM

# RBP→lncRNA targets (cellType=all, geneType=lncRNA)
.venv/bin/python3 -m pipeline.lncRNA_interactions.encori_cli \
  --rbp-list EZH2,SUZ12,EED,BRD4,EP300,HDAC1,HDAC2 \
  --clip-exp-num 2

# miRNA→lncRNA targets for the pipeline’s selected lncRNAs
.venv/bin/python3 -m pipeline.lncRNA_interactions.encori_cli \
  --do-mirna --clip-exp-num 2 --n-extra-close-lncrnas 20
```

### lncRNA interactions (unified outputs: ENCORI miRNAs + POSTAR3 RBPs)

This produces one folder of canonical interaction tables under `data/lncRNA_interactions/`
(selected lncRNA list + ENCORI miRNA-target rows + POSTAR3-derived RBP summaries).

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m pipeline.lncRNA_interactions.build_all \
  --n-extra-close-lncrnas 20 \
  --encori-clip-exp-num 2 \
  --encori-rbp-clip-exp-num 2 \
  --encori-rbp-top-n 12 \
  --postar3-region-mode exons \
  --postar3-chunk-rows 500000
```

### lncRNA predicted miRNA targets (tool-based; RNAhybrid)

Requires:
- `bedtools`
- `RNAhybrid` (Ubuntu package: `rnahybrid`)
- A mature miRNA FASTA (arm-level names), e.g. `data/miRNA/mature.fa`
- A genome FASTA (defaults to `PATHS.sv_reference_fasta`)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m pipeline.lncRNA_interactions.predicted_targets_cli \
  --mirna-fasta data/miRNA/mature.fa \
  --n-extra-close-lncrnas 20
```

### SNV cohort run (FIMO + ChIP overlap)

`run_cohort_module_processing.py` passes **`--snv-fimo`** by default (use `--no-snv-fimo` to skip). Requires `bedtools`, `fimo` (see `pipeline/SNV/snv_fimo.resolve_fimo_argv`), and `PATHS.sv_meme_file`.

It also passes **`--snv-chip`** by default (use **`--no-snv-chip`** to skip): strict **POS** overlap with **`PATHS.chip_unified`** → columns `snv_chip_hits` / `snv_chip_aggregate` (see `pipeline/md/DATA_STRUCTURE_ATLAS.md` §3).

SV + motif pipeline (stages, thresholds, streaming, neojunction, ChIP): `pipeline/md/module_specific_processing_md/sv_→_motif_pipeline_vep_fimo.md` (canonical). `pipeline/SV/SV_PIPELINE_README.md` is a pointer to that file.

### Curated SV/SNV MEME library (rebuild from a database dump)

PWM databases are **downloaded manually** from HOCOMOCO / JASPAR / CIS-BP (or your mirror); the repo only **filters and merges** them in MEME text form via `extract_selected_motifs()` (`pipeline/SV/motif_scanning.py`).

The checked-in library is `data/SV/motifs/core_ifn_stat_ctcf_ap1_nlrc5.meme` (see `PATHS.sv_meme_file`). The PWM **name filter** is `SV_TARGET_TF_SYMBOLS` in `pipeline/config.py` (documented in `pipeline/md/module_specific_processing_md/sv_→_motif_pipeline_vep_fimo.md`). To rebuild from a larger HOCOMOCO/JASPAR MEME file:

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -c "
from pathlib import Path
from pipeline.SV.motif_scanning import extract_selected_motifs
from pipeline.config import SV_TARGET_TF_SYMBOLS, PATHS
src = Path('path/to/hocomoco_or_jaspar_full.meme')   # <-- your download
extract_selected_motifs(src, PATHS.sv_meme_file, SV_TARGET_TF_SYMBOLS)
print('wrote', PATHS.sv_meme_file)
"
```

