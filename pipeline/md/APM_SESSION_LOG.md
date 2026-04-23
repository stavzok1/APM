# APM session log (reconstructed)

**Purpose**: Single place for **(a)** what changed scientifically, **(b)** where artifacts live on disk, and **(c)** how to rerun key commands. Written for your own continuity (not a handoff doc).

**Scope**: This file summarizes the **long multi-topic working session** whose detailed tool trace lives in the Cursor transcript [APM multi-topic session](281c0a87-4090-4ea4-93ba-f344add95be0), plus follow-on fixes in the lncRNA / RNAhybrid lane.

**Companion**: Canonical output shapes remain in `pipeline/md/DATA_STRUCTURE_ATLAS.md` (updated in the same pass as this log to include late **lncRNA interactions** and a few other tables that were missing).

---

## (a) Science and integration outcomes (what actually moved)

### Data registry and documentation

- Confirmed **HiChIP TCGA BRCA processed** matrix and **Xena arm-specific miRNA** matrix are first-class `PATHS` inputs and feed annotation manifest builders; fixed stale docs that still referenced older miRNA paths.
- Expanded **`DATA_STRUCTURE_ATLAS.md`** so late outputs (lncRNA module, RNAhybrid intermediates, POSTAR3 master parquet, unified ChIP parquet, BRCA isoform filter outputs, richer covariates notes) are discoverable.

### ChIP-seq (large refactor)

- Unified ENCODE + ChIP-Atlas peaks now carry a per-file **`sample_id`** so **replicate TF×cell-line** tracks stay distinct (`CTCF_MCF7_1`, `CTCF_MCF7_2`, …) while summaries still expose aggregate **`n_samples`** / **`sample_ids`** under `chip_hits` in the element table.
- **`chip_loader.py`** now tolerates **ChIP-Atlas UCSC BED9+** exports (track line + gffTags), not only legacy six-column tables.
- **`scripts/build_unified_chip_peaks.py`** is the supported way to (re)build `PATHS.chip_unified`; regulatory `main` **consumes** the cached parquet if present but does **not** rebuild it automatically.

### SNV: FIMO on reference windows

- **`load_mutect_snv_vcf`** (and batch loaders) support **`run_fimo=True`**, adding **`fimo_hits`** per retained variant by scanning a fixed-width **hg38 reference** window with MEME **FIMO** (orthogonal to VEP `motif_hits`). Cohort wiring uses the existing analysis flags documented in SNV module docs.

### Isoform expression at TCGA scale

- Added a **chunked, disk-first** workflow to subset a **pan-cancer isoform TPM** table to **BRCA participants** and to transcripts relevant to the **analysis gene panel**, with transcript→gene tagging via the GENCODE probemap (`pipeline/RNA_exp/filter_isoform_expression_brca.py`, reporting helper alongside).

### Cohort covariates (post-`main`)

- Extended **`pipeline/covariates/build_covariates.py`** with additional providers and join logic:
  - **DDRscores** table ingestion (`--ddr-scores`).
  - **HiChIP participant** matrix summarized and **lifted to `sample_vial`**.
  - **Pipeline-native SNV** summaries from per-sample JSON, with a **filename→tumor vial** mapping fix via `annotations/SNV/samples.tsv` (UUID-ish intermediate directory names are not sample IDs).
  - Explicit **SNV+CNV stratifier anchors** (`TP53`, `PTEN`, `PIK3CA`) are used as genes-of-interest in SV mirroring and are supported via MAF-style mutation covariates and DDR-style tables where available.
- **Sample coverage canonical path**: `analysis/sample_coverage/output/current/` (opt-in mirror via `--write-current` on the coverage script); covariates coverage provider prefers this when present.
- **Vial plurality**: base covariate index can carry **`n_vials_per_sample`**, **`n_vials_per_participant`**, and related replication flags.

### Research panel governance (tiers)

- Encoded an “iterated panel verdict” into **`pipeline/config.py`** tier lists and synchronized narrative in **`research_plan/01_gene_panel_extended.md`** (Tier 1–4 additions; explicit separation of **panel subjects** vs **covariates / stratifiers**).

### lncRNA interactions (ENCORI + POSTAR3 + predicted miRNA targets)

- **ENCORI**: more robust TSV parsing (comment lines, one-line errors), stable empty schemas when an RBP is not in ENCORI’s CLIP universe, **ENSG-style `gene_id` fallback** for RBPTarget with diagnostics, and practical documentation in `analysis/COMMANDS.md`.
- **POSTAR3**: streaming **POSTAR3.txt → parquet**; overlap summaries support **feature-aware region modes** (`gene` / `exons` / `introns` / `promoter`) with provenance columns on overlap rows.
- **RNAhybrid**: end-to-end local prediction for **miRNA vs spliced lncRNA exon sequences**, including **genome FASTA vs GENCODE `chr` prefix harmonization**, **3′ truncation** for very long transcripts, **auto `-m/-n`**, **ENCORI-restricted miRNA FASTA** option, and a **parser aligned to RNAhybrid `-c` compact lines**.

### Repository hygiene

- Fixed stale **`scripts/...` paths** after script reorganization (notably module sample coverage under `scripts/coverage/…` and smoke-test registry expectations).

---

## (b) Artifacts and paths (disk truth)

Paths below use **`PATHS.working_dir`** defaults (`data/` …) where applicable. If you changed `working_dir`, prepend your configured root.

### ChIP

| Artifact | Typical path | Notes |
|-----------|--------------|--------|
| Unified peaks parquet | `data/CHIP/unified_chip_peaks.parquet` | `sample_id` per BED stem; rebuild via `scripts/build_unified_chip_peaks.py` |
| ENCODE / ChIP-Atlas raw BEDs | `data/CHIP/ENCODE/*.bed`, `data/CHIP/CHIP_ATLAS/*.bed` | filenames define `sample_id` |

### POSTAR3 (RBP peaks)

| Artifact | Typical path | Notes |
|-----------|--------------|--------|
| Master peak table | `data/RBP-RNA/POSTAR3.parquet` | Built from `data/RBP-RNA/POSTAR3.txt` via `scripts/RBP-RNA/build_postar3_parquet.py` |

### lncRNA interactions

| Artifact | Path under `data/lncRNA_interactions/` | Notes |
|-----------|----------------------------------------|--------|
| Selected lncRNAs | `selected_lncrnas.txt` | deterministic 8 + N closest |
| ENCORI miRNA targets | `encori_mirna_targets.parquet` | + `encori_mirna_targets_diagnostics.csv` |
| ENCORI RBP targets | `encori_rbp_targets.parquet` | + `encori_rbp_targets_diagnostics.csv` |
| POSTAR3 overlaps / summaries | `postar3_overlaps.parquet`, `postar3_rbp_summary.parquet`, `postar3_lncrna_summary.parquet` | overlaps include `__selected_set`, `__region_mode` when written |
| RBP shortlist | `recommended_rbps_from_postar3.csv` | convenience |
| RNAhybrid tree | `predicted_targets/rnahybrid/*` | see **§8c** in `DATA_STRUCTURE_ATLAS.md` for the full file list |

### Isoform expression (BRCA subset)

| Artifact | Typical path | Notes |
|-----------|--------------|--------|
| Filtered isoform table | `data/RNAexp_TCGA/isoform_brca_panel.parquet` (or `.tsv`) | CLI chooses |
| QC JSON | `data/RNAexp_TCGA/isoform_brca_panel_qc.json` | match rates, column coverage |

### Covariates

| Artifact | Path | Notes |
|-----------|------|--------|
| Covariates table | `data/covariates/<run_id>/cohort_covariates.parquet` | optional `.csv` |
| Run metadata | `data/covariates/<run_id>/metadata.json` | provider diagnostics |
| DDRscores input (example) | `data/covariates/DDRscores.txt` | pass with `--ddr-scores` |

### Sample coverage (omics presence)

| Artifact | Path | Notes |
|-----------|------|--------|
| Stable mirror | `analysis/sample_coverage/output/current/` | populate via `python analysis/sample_module_coverage.py --write-current` |
| Timestamped runs | `analysis/sample_coverage/output/run_*/` | fallback when `current/` absent |

### ENCORI caches

| Artifact | Path | Notes |
|-----------|------|--------|
| Raw ENCORI TSVs | `data/external_cache/encori/...` | per-query cache |

---

## (c) Commands (rerun / reproduce)

Run from repo root on WSL with **`PYTHONPATH`** set to the repo (or use **`./.venv/bin/python3`** after `cd`).

### Unified ChIP peaks

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 scripts/build_unified_chip_peaks.py
.venv/bin/python3 scripts/build_unified_chip_peaks.py --inspect-only
```

### lncRNA interactions (ENCORI + POSTAR3)

```bash
cd /home/stavz/masters/gdc/APM
PYTHONPATH=/home/stavz/masters/gdc/APM .venv/bin/python3 -m pipeline.lncRNA_interactions.build_all \
  --postar3-region-mode exons
```

ENCORI-only CLIs remain available as documented in **`analysis/COMMANDS.md`** (`encori_cli`).

### RNAhybrid predicted miRNA targets

```bash
cd /home/stavz/masters/gdc/APM
PYTHONPATH=/home/stavz/masters/gdc/APM .venv/bin/python3 -m pipeline.lncRNA_interactions.predicted_targets_cli \
  --mirna-fasta data/miRNA/mature.fa \
  --mirna-subset-from-encori \
  --n-extra-close-lncrnas 20 \
  --max-target-bp 8000
```

Use **`--max-target-bp 0`** to disable 3′ truncation (expect failures or extreme runtime on very long lncRNAs).

### POSTAR3.txt → parquet (streaming)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 scripts/RBP-RNA/build_postar3_parquet.py \
  --in data/RBP-RNA/POSTAR3.txt \
  --out data/RBP-RNA/POSTAR3.parquet \
  --chunksize 200000
```

### BRCA isoform filter (chunked)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 pipeline/RNA_exp/filter_isoform_expression_brca.py \
  --out-table data/RNAexp_TCGA/isoform_brca_panel.parquet \
  --qc-json data/RNAexp_TCGA/isoform_brca_panel_qc.json
```

### Cohort covariates

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m pipeline.covariates.build_covariates --run-id run_cov_001
.venv/bin/python3 -m pipeline.covariates.build_covariates --run-id run_cov_ddr --ddr-scores data/covariates/DDRscores.txt
```

### Smoke tests (use venv interpreter)

```bash
cd /home/stavz/masters/gdc/APM
.venv/bin/python3 -m tests.run_smoke_tests
```

---

## Timeline: this long session (workstreams, chronological-ish)

1. **Data structure + atlas** for HiChIP BRCA processed + Xena miRNA; fix doc/TOC drift; clarify venv vs system `python3`.
2. **ChIP multi-sample `sample_id`**, aggregation under `chip_hits`, rebuild script, loader fixes for new ChIP-Atlas BED format.
3. **Operational clarification**: regulatory `main` does not rebuild `chip_unified`; SV path may.
4. **ChIP consumer audit** (where cell-type whitelists apply).
5. **SNV FIMO** integration + smoke execution on a real VEP VCF.
6. **Isoform pan-cancer → BRCA** chunked filter + QC/report scripts.
7. **Covariates** extensions (DDR, HiChIP participant lift, native SNV join fixes, coverage `current/`, vial plurality).
8. **Gene tier panel** alignment in config + research plan markdown; covariate/signature **inventory discussion** (no code change when requested).
9. **Repo path sweep** after `scripts/` reorg (smoke tests + analysis references).
10. **ENCORI** hardening + diagnostics + RBPTarget availability reality check.
11. **POSTAR3** parquet builder + lncRNA overlap **region modes**.
12. **RNAhybrid** pipeline + WSL/tooling + reference contig naming + truncation + hit parsing.
13. **Atlas update** (this log’s sibling deliverable) to register the above outputs explicitly.

---

## Timeline: APM overall (since project start; high level, reconstructed)

This is a **narrative reconstruction** of how the repo reads today—not a git-blame audit. Treat dates as **relative phases**.

| Phase | What landed (conceptually) |
|-------|----------------------------|
| **A. Regulatory backbone** | cCRE-centric **`regulatory_element_focus_with_evidence`** pipeline; nested JSON columns; distance rings; thresholds doc; GENCODE-driven gene and lncRNA interval logic. |
| **B. Multi-omics evidence layers** | RNA TPM matrices + normalization helpers; ATAC matrices + SCREEN/ABC hooks where configured; methylation probe→gene logic; RPPA composites and “block” heuristics. |
| **C. Germline/somatic structure variants** | SV pipeline with element/gene/lncRNA/miRNA hits, TAD boundaries, optional motif scanning; strict SV tables for downstream mirroring. |
| **D. SNV + CNV** | VEP-driven SNV tables with rich CSQ parsing; optional **FIMO** on reference windows; CNV segments with gene/lncRNA/miRNA/element hits. |
| **E. Gene-centric “mirroring”** | SV (and related) disruption summaries lifted into **sample-centric** covariates and QC surfaces. |
| **F. Immune / clinical annotation** | Unified BRCA clinical + immune tables; Thorsson and immune subtype satellites; miRTarBase processing outputs. |
| **G. Larger “world model” inputs** | HiChIP TCGA processed matrix registration; Xena miRNA matrix; ChIP peak corpus growth → unified parquet + replicate-aware `chip_hits`. |
| **H. lncRNA program** | lncRNA matching tables; **ENCORI** + **POSTAR3** interaction bundle; **RNAhybrid** predicted miRNA targeting; atlas wiring. |
| **I. Expression at transcript resolution** | Isoform TPM subsetting for BRCA + panel genes (streaming). |
| **J. Operational maturity** | `DATA_STRUCTURE_ATLAS` as the column bible; `analysis/COMMANDS.md`; smoke tests; WSL bridge rules; script directory hygiene. |

If you later want a **machine-generated** “timeline since first commit,” we can add a short script once the repo is initialized as a git worktree in this environment.
