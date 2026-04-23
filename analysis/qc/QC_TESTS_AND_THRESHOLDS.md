## QC tests and thresholds (current suite)

This document is the **single place** that enumerates what QC checks are currently implemented, **what files they read**, and the **definitions + thresholds** used to emit pass/warn/fail (or to hard-stop Tier 1 structural checks).

Code entry points:
- `analysis/qc/run_qc_suite.py`: one command runner (RPPA build + QC reports)
- `analysis/qc/run_cross_modal_qc.py`: cross-modal QC report (Tier 1 + selected Tier 3)
- `analysis/qc/run_full_qc.py`: modular QC report (Tier 2 within-modality + Tier 3 SV/LOH)

Output root:
- `analysis/qc/output/<run_id>/`

---

### Tier 1 — structural integrity (hard-stop)

#### 1) Scanning-column integrity (element-focus)
- **Where**: `analysis/qc/run_cross_modal_qc.py::_qc_scanning_columns_integrity`
- **Inputs**: `PATHS.regulatory_elements_table_with_evidence_parquet`
- **Definition**: Load a small slice of the element-focus parquet and recompute a subset of `scan_*` columns from nested JSON (`gene_links`, `screen_exp`, `ABC_enhancers`, `hichip`, `chip_hits`) using `pipeline.scanning_columns.derive_elem_focus_scanning_columns()`. Compare stored vs recomputed per row.
- **Threshold**:
  - **Hard fail** if `n_scan_cols_checked > 0` **and** `min_ok_rate < 0.999`
  - If scan columns are missing entirely, the check is treated as vacuously OK (no hard fail).

---

### Tier 2 — within-modality sanity (warning-oriented)

#### 2) Methylation promoter-beta distribution sanity
- **Where**: `pipeline/qc/within_modality.py::compute_within_modality_sanity`
- **Inputs**: `ctx.methylation_working_dir/cohort/gene_meth_matrix.csv`
- **Definition** (per sample): compute `mean`, `median`, `p05`, `p95` across genes.
- **Thresholds**:
  - If the matrix is missing: emit `WARN` (`within_modality.methylation_gene_matrix_present`)
  - No numeric “fail” thresholds yet (this is descriptive QC).

#### 3) RNA sanity (subset)
- **Where**: `pipeline/qc/within_modality.py`
- **Inputs**: `PATHS.rna_expression` via `read_tpm_wide_subset(...)`
- **Gene set**: APM class-I set from `pipeline.RNA_exp.signatures.GeneSets().apm_class_i`
- **Definition** (per sample): compute `mean_log2tpm1`, `std_log2tpm1` across that gene set.
- **Thresholds**:
  - If stats can’t be computed at all: emit `WARN` (`within_modality.rna_subset_present`)
  - If median per-sample std dev is extremely low: emit `WARN` (`within_modality.rna_subset_variance`)
    - Current cutoff: `median(std_log2tpm1) < 0.05`

#### 4) CNV gene-call sanity (presence/scale)
- **Where**: `pipeline/qc/within_modality.py`
- **Inputs**: per-sample gene-call files under `ctx.cnv_gene_tables_dir`:
  - `<sample>_cnv_gene_calls_ascat3.csv` (preferred) else `<sample>_cnv_gene_calls.csv`
- **Definition** (per sample): locate a CNV numeric column in (`cnv_log2`, `segment_mean`, `log2_copy_ratio`, `log2_ratio`) or compute `log2(copy_number/2)` if `copy_number` exists; compute mean/std across genes.
- **Thresholds**:
  - If no gene-call files found for any sample: emit `WARN` (`within_modality.cnv_gene_calls_present`)
  - No numeric “fail” thresholds yet (this is descriptive QC).

#### 5) RPPA panel-score sanity (gross degeneracy only)
- **Where**: `pipeline/qc/within_modality.py::assert_within_modality_sanity`
- **Inputs**: `ctx.rppa_output_dir/panel_scores.csv`
- **Important definition note**:
  - RPPA “panel scores” are computed as **mean of per-target z-scores** (see `pipeline/rppa/rppa_panels.py::compute_zscore_panel_score`).
  - Therefore, a panel score’s std dev is **not expected** to be ~1.0 in general.
- **Thresholds**:
  - If `panel_scores.csv` missing: emit `WARN` (`within_modality.rppa_panel_scores_present`)
  - For a small set of score columns (IFN/DDR/cGAS if present, else first ~8 cols):
    - Emit `WARN` if `std` is outside **[0.15, 2.5]** (`within_modality.rppa_panel_score_scale`)
    - Emit `WARN` if `abs(mean) > 0.3` (`within_modality.rppa_panel_score_centering`)

---

### Tier 3 — cross-modal biological consistency (warning-oriented unless noted)

#### 6) Promoter methylation ↔ RNA expression
- **Where**: `analysis/qc/run_cross_modal_qc.py::_qc_methylation_vs_rna`
- **Inputs**:
  - Methylation: `ctx.methylation_working_dir/cohort/gene_meth_matrix.csv`
  - RNA: `PATHS.rna_expression` via `read_tpm_wide_subset(...)`
- **Definition** (per gene): Spearman + Pearson correlation across samples between promoter beta and RNA log2(TPM+1).
- **Discordance counters** (per sample; used for “swap-like” signals):
  - Count genes where `beta > 0.7` and expression is above gene’s 75th percentile
  - Count genes where `beta < 0.2` and expression is below gene’s 25th percentile
- **Thresholds**:
  - Currently no automated warn/fail thresholds; reported for interpretation.

#### 7) CNV (gene-level log2) ↔ RNA dosage
- **Where**: `analysis/qc/run_cross_modal_qc.py::_qc_cnv_vs_rna`
- **Inputs**:
  - CNV: per-sample gene calls (`copy_number` → `cnv_log2 = log2(cn/2)` if needed)
  - RNA: `PATHS.rna_expression`
- **Definition**:
  - Per gene: Spearman + Pearson correlation across samples between `cnv_log2` and expression
  - Per sample: `cnv_expr_concordant_frac` heuristic (gain→expr above median; loss→below)
- **Thresholds**:
  - Currently no automated warn/fail thresholds; reported for interpretation.

#### 8) IFN axis ↔ APM response
- **Where**: `analysis/qc/run_cross_modal_qc.py::_qc_ifn_axis`
- **Inputs**:
  - Thorsson: `annotations/Thornsson_immune_table.tsv` (`IFN-gamma Response`)
  - RNA: `PATHS.rna_expression` gene-set scores
- **Definition**:
  - `apm_response_to_ifn_score = APM_classI_mean - median(APM_classI_mean)`
  - Flag `ifn_high_apm_low_flag` if IFN score ≥ 75th percentile AND APM response ≤ 25th percentile
- **Thresholds**:
  - The flag is produced; no additional warn/fail thresholds yet.

#### 9) ATAC accessibility ↔ RNA in cis (linked peaks)
- **Where**: `analysis/qc/run_cross_modal_qc.py::_qc_atac_vs_rna_cis`
- **Inputs**:
  - Peak↔gene pairs: `PATHS.regulatory_elements_table_with_evidence_parquet` columns `gene_links` + `atac_peak_links`
  - Peak accessibility: `PATHS.atac_case_level_matrix` (streamed; TCGA aliquot columns are prefix-matched and then normalized to TCGA sample ids)
  - RNA: `PATHS.rna_expression`
- **Definitions**:
  - Per (gene, peak_id): Spearman/Pearson correlation across samples between peak accessibility and RNA(gene)
  - Evidence tags per pair (from element-focus `gene_links[gene]` bundle):
    - `ccre_has_screen_strong` (any SCREEN “strong” in `screen_exp.per_biosample.*.*`)
    - `ccre_has_hichip` (non-empty `hichip`)
    - `ccre_has_abc` (non-empty `ABC_enhancers`) **tag only**
    - `ccre_evidence_supported = has_screen_strong OR has_hichip`
    - `ccre_gene_tier`, `ccre_gene_dist_to_tss` (when present)
  - Gene-level aggregation:
    - For each gene, per sample: `gene_atac_max = max(ATAC across linked peaks)`
    - Correlate `gene_atac_max` vs RNA(gene)
    - Report `ccre_evidence_supported_rate` (fraction of linked cCREs with supported evidence)
- **Thresholds**:
  - No automated warn/fail thresholds yet; outputs are intended for cohort-specific calibration.

#### 10) SV disruption ↔ RNA expression loss (per gene; adjusted)
- **Where**: `pipeline/qc/sv_vs_rna.py::compute_sv_disruption_vs_rna`
- **Inputs**:
  - SV per-sample strict tables: `ctx.sv_output_root/07_final_sv_with_fimo/<sample>_strict_sv_set.csv`
  - RNA: `PATHS.rna_expression` subset (APM class-I by default)
  - CNV: gene calls `copy_number` → `cnv_log2 = log2(cn/2)`
  - IFN covariate: Thorsson `IFN-gamma Response` (preferred), else RNA signature `IFNG_mean`
- **Definitions** (per gene):
  - `SV_hit`: sample has any SV gene hit for gene (promoter or CDS flags)
  - `delta_median_unadjusted = median(expr | SV_hit) - median(expr | !SV_hit)`
  - CNV-neutral subset: `abs(cnv_log2) <= 0.2`
  - OLS regression (when enough samples): `expr ~ SV_hit + cnv_log2 + IFN`
- **Thresholds**:
  - Minimum sample sizes to compute adjusted comparisons:
    - `min_hit = 5`, `min_no = 10`
  - CNV-neutral band: `abs(cnv_log2) <= 0.2`
  - No cohort-wide warn/fail threshold is enforced yet; outputs are per-gene for interpretation.

#### 11) HLA LOH ↔ neoantigen / infiltration (join sanity)
- **Where**: `pipeline/qc/loh_vs_neoantigen.py` (+ `pipeline/qc/hla_loh_ascat_segments.py`)
- **Inputs**:
  - **Preferred**: ASCAT3 **allelic** segment files (`PATHS.cnv_dir` / `*.ascat3.allelic_specific.seg.txt`), indexed from `PATHS.cnv_annotations_path` (`annotations/CNV/samples.tsv`) by tumor vial id.
  - **Fallback**: CNV gene calls `loh_flag` for `HLA-A`, `HLA-B`, `HLA-C`, `B2M` in `*_cnv_gene_calls_ascat3.csv` (gene-level ASCAT export can miss allelic imbalance when `min_copy_number == max_copy_number`).
  - **Optional**: HLA OptiType-style table `PATHS.hla_types_tsv` (`annotations/HLA_types.tsv`), joined by **prefix match** (`aliquot_id.startswith(<tumor_vial>)`).
  - **Optional**: per-sample SNV tables under ``<snv_output>/per_sample/*_snv_variants.csv`` (indexed from ``annotations/SNV/samples.tsv``), using parsed VEP **transcript** rows in ``gene_hits`` for `HLA-A` / `HLA-B` / `HLA-C` / `B2M`.
  - Thorsson immune table: neoantigens + infiltration proxies (if present)
- **Definitions**:
  - Gene body intervals for `HLA-A/B/C` and `B2M` are read from the per-sample gene table (same file as fallback), then overlapped with allelic segments. The segment with **largest overlap bp** with each gene body is chosen; **segment LOH** if that segment has `Minor_Copy_Number == 0` and overlap covers **≥ 30%** of the gene interval length.
  - `hla_loh_gene_table_any = any(loh_flag > 0)` on HLA genes (fallback path).
  - `hla_loh_seg_any = any(HLA-A/B/C/B2M segment LOH)` when a segment file is resolved for the vial.
  - `hla_loh_any`: **segment-based** `hla_loh_seg_any` when a segment file is available; otherwise falls back to `hla_loh_gene_table_any`.
  - **SNV / VEP (HLA loci, somatic table)** — `pipeline/qc/hla_snv_vep.py`:
    - **Audit (unfiltered)**: `hla_snv_lof_any`, `hla_snv_high_missense_any`, `hla_snv_clin_path_any` (same definitions as before: LoF-style consequences or LoFtee-style **LoF** flags; `IMPACT == HIGH` + `missense_variant`; ClinVar-style `CLIN_SIG`).
    - **Damage gates** (used for `hla_damage_cnv_or_snv_any`): max population AF over `gnomADe_AF` / `gnomADg_AF` / `gnomAD_AF` / `MAX_AF` / `AF` must be **≤ 0.01** when any such field is numeric; **no** such fields parsed ⇒ treat as eligible (somatic exports often omit MAF).
    - `hla_snv_lof_rare_any`: LoF-like hit passes the MAF gate above.
    - `hla_snv_high_missense_damage_any`: HIGH missense passes the MAF gate **and** either (i) no CADD/REVEL/SpliceAI scores are present on the hit, or (ii) **CADD_PHRED ≥ 20**, or **REVEL ≥ 0.5**, or max(**SpliceAI_pred_DS_AG**, **SpliceAI_pred_DS_AL**) **≥ 0.5**, or pathogenic **CLIN_SIG** (same strings as `hla_snv_clin_path_any`).
    - `hla_snv_clin_path_damage_any`: pathogenic **CLIN_SIG** and the same MAF gate (excludes common “pathogenic” tagging with high AF when gnomAD is populated).
    - `hla_snv_min_pop_af_seen` / `hla_snv_max_tumor_vaf`: min of gnomAD/MAX/AF fields and max `tumor_vaf` over rows that carry an HLA transcript hit (for manual QC).
    - `hla_damage_cnv_or_snv_any`: `hla_loh_any` **or** `hla_snv_lof_rare_any` **or** `hla_snv_high_missense_damage_any` **or** `hla_snv_clin_path_damage_any` (presentation-restriction **proxy**, not phasing-aware “LOH”).
  - `hla_type_*` columns: `A1`, `A2`, `B1`, `B2`, `C1`, `C2`, `Reads`, `Objective`, plus `hla_typing_aliquot_id` when a typing row matches.
  - Join to Thorsson columns (neoantigens, leukocyte fraction, etc.)
- **Thresholds**:
  - If Thorsson table missing: emit `WARN` (`hla_loh_vs_neoantigen.thornsson_present`)
  - If no allelic segment files are indexed from the manifest: emit `WARN` (`hla_loh.ascat_seg_index`)
  - If no per-sample SNV files are indexed from the SNV manifest: emit `WARN` (`hla_snv.per_sample_index`)
  - No controlled regression is implemented yet (this is currently a join/coverage sanity check).

---

### How SCREEN “strength” is assigned (what it means)

SCREEN strength labels are assigned in `pipeline/evidence/screen_links.py::_apply_strength_classification`:

- Scores are grouped by **assay_type only** and compared to empirical quantiles:
  - `weak` if `score >= q_weak` where `q_weak = quantile(score, screen_weak_quantile)` (default 0.90)
  - `strong` if `score >= q_strong` where `q_strong = quantile(score, screen_strong_quantile)` (default 0.95)
  - `none` otherwise
- For **Intact-HiC** only (experimental links), `strong` additionally requires:
  - `p_value <= intact_hic_pvalue_threshold` (default 1e-6)

Important interpretation caveats:
- “weak/strong” is a **relative rank within the assay’s score distribution** in the *extracted subset* (after filtering to the pipeline’s gene list).
- It is **not** recomputed per biosample, so biosamples with systematically different score scales are not normalized separately.
- Because quantiles are computed after gene filtering, the absolute meaning of “top 5%” can shift if you change the gene panel.

Reasonableness assessment (current defaults):
- **Pros**:
  - Assay-specific quantiles avoid mixing incomparable assays.
  - The Intact-HiC p-value gate makes “strong” more defensible for that assay.
  - For evidence integration/QC, a rank-based “strong” tag is often enough to prioritize or stratify links.
- **Potential pitfalls**:
  - Quantiles computed on a *panel-filtered* subset can make “strong” less stable across runs/panels.
  - Per-assay-only quantiles can under/over-call strong for particular biosamples if score distributions vary by biosample.
  - Computational links have dummy p-values (NaN), so the p-value gate does not apply there.

If we want more “absolute” or stable behavior, the next steps would be:
- Precompute assay quantiles on a broader SCREEN baseline (or cache per-run quantiles with provenance).
- Optionally compute quantiles **per (assay_type, biosample)** and store both a global and biosample-normalized strength tag.
- Prefer using the raw `score` + `p_value` fields for downstream modeling, treating `"strong"/"weak"` as a coarse convenience label.

Implementation note (now implemented):
- The main pipeline persists the **per-run per-assay quantiles** it used to:
  - `working_dir/evidence/screen_exp_score_quantiles.csv` (+ `.meta.json`)
  - `working_dir/evidence/screen_comp_score_quantiles.csv` (+ `.meta.json`)
  These capture the quantiles used for that run and basic provenance (source paths + gene set size).

