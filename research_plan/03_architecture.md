# 03 — Analysis architecture

Layer contract for the analysis layer built on top of `pipeline/`. **No code in this document** — only the interface.

---

## L0 — Atlas tables (provided by `pipeline/`)

Canonical definitions live in `pipeline/md/DATA_STRUCTURE_ATLAS.md`. The L0 tables consumed by the analysis layer:

| Table | Granularity | Key nested fields |
|---|---|---|
| Element Focus | per cCRE | `gene_links` (SCREEN/ABC/HiChIP/ChIP per biosample + conservation), `chip_hits`, `TAD_domains` × 24 |
| ATAC Peak | per peak | `gene_links` with `body_overlap`, `TAD_boundary_overlaps` |
| ATAC case-level | peak × sample | log-CPM QN (~74 BRCA) |
| SNV | per allele | `motif_hits` with Δ, `cCRE_hits`, `gene_hits` (VEP) |
| SV | per call | `elem_hits.motif_hits` (FIMO), `flank_motif_hits`, `chip_hits`, `gene_hits` |
| CNV | per segment | `gene_hits`, `elem_hits`, `cn_state`, `loh_flag` |
| Methylation | per sample | 6c gene-level (`promoter_beta_*`), 6d cCRE-level (`ccre_beta_*`) |
| RPPA | per sample | `panel_scores`, `signaling_blocks`, `protein_rna_discordance`, `estimated_apm_capacity` |
| RNA | gene × sample | log2(TPM+1) |
| 13 per-tumor HiChIP | loop × sample (13 samples) | Calibration anchor for conservation prior |
| miRTarBase | per (miR, gene) | `evidence_score`, exp × support cross counts |
| Thorsson | per sample | C1–C6, activated/resting cell fractions |
| Clinical + immune unified | per sample | stage, PAM50, vital status, CPE |

---

## L1 — Per-sample × per-gene feature builders

Each builder returns a tidy table `(sample, gene, feature_value)` or a dense matrix with sample index and gene column. One function, one feature. Names below are the contract — future code must match.

### Raw layer features

| Builder | Returns |
|---|---|
| `build_expression_matrix()` | `samples × genes`, log2(TPM+1) |
| `build_cnv_dosage(gene_set)` | `samples × genes`, gene-averaged log2 ratio + `cn_state` + `loh_flag` |
| `build_promoter_methylation(gene_set)` | `promoter_beta_mean`, `promoter_beta_std`, `promoter_frac_hypermeth`, `promoter_CGI_beta_mean` |
| `build_genebody_methylation(gene_set)` | `gene_body_beta_mean` |
| `build_enhancer_methylation(gene_set)` | Enhancer-level β via `gene_links` → `cCRE_id` → `ccre_beta_mean` |
| `build_enhancer_atac(gene_set, case_level=True)` | Mean / sum accessibility at gene's enhancers |
| `build_rppa_panel()` | RPPA `panel_scores`, `signaling_blocks`, `estimated_apm_capacity` |
| `build_protein_rna_discordance(gene_set)` | `PRD_g` signed |
| `build_ifn_signature()` | Hallmark IFNG ssGSEA + IFNα |
| `build_thorsson_states()` | `C1..C6`, activated/resting fractions, `CYT`, CD8a |
| `build_clinical()` | stage, PAM50, OS/PFS, treatment flags, purity |

### Element-resolved novel features

| Builder | Returns | Used by |
|---|---|---|
| `build_crml_snv(gene_set, tf_set, chip_required=True)` | `CRML_snv_g` | H1, H2, H3 |
| `build_crml_sv(gene_set, tf_set, chip_required=True)` | `CRML_sv_g` | H1, H3 |
| `build_enhancer_redundancy(gene_set, sources=[ABC,HiChIP,SCREEN])` | `R_g` | H7, H8, H54 |
| `build_boundary_conservation_disruption(gene_set)` | `BCD_g` conservation-weighted | H9, H55 |
| `build_mirna_pressure(gene_set, min_evidence=4)` | `MIR_g` | H26–H31 |
| `build_lncrna_cis_panel(gene_set, window=1_000_000)` | lncRNA cis-matrix | H32–H37 |

### Within-panel composite builders

| Builder | Returns | Used by |
|---|---|---|
| `build_jak_stat_intactness()` | `JSI` from CNV, SV, damaging SNV, β of {JAK1, JAK2, STAT1} | H12 |
| `build_stat1_stat3_ratio(source="rna"\|"rppa")` | `R_STAT` | H40 |
| `build_irf1_irf2_ratio()` | `R_IRF` | H41 |
| `build_feedback_tension()` | SOCS1 residual given IFN response | H42 |
| `build_pathway_coherence()` | PC1 eigenratio across APM | H44 |
| `build_peptide_pipeline_coherence()` | Cross-station coherence | H47 |
| `build_dual_axis_visibility()` | MHC-I idx, NKL idx, shedding; `QUAD` | H13, H14 |
| `build_cgas_sting_intactness()` | From Tier-1 extended panel | H48, H56 |
| `build_ifn_receptor_ratio()` | `log(IFNGR·/IFNAR·)` | H58 |

---

## L1.5 — Conservation-normalized counterparts

For every feature with a multi-biosample or subtype-consensus version, a residual column is added:

| Feature | Conservation-normalized counterpart |
|---|---|
| `enhancer_atac[g]` | `enhancer_atac[g] − atac_consensus_subtype(g)` |
| `hichip_edge[e→g]` | Calibrated via 13 tumor HiChIP; residual vs subtype consensus |
| `abc_score[e→g]` | Residual vs breast-biosample consensus |
| `boundary_strength[b]` | Residual vs 24-biosample consensus |
| `chip_binding[tf, e]` | Residual vs same-subtype ChIP average |

Details in `04_conservation_prior.md`.

---

## L2 — Composite scores

| Score | Inputs | Used by |
|---|---|---|
| `APM_intrinsic` | RNA residualized on IFN-γ + purity + PAM50 | H2, H16, H17 |
| `IVI` (Intrinsic Visibility Index) | sparse PLS of {`CRML`, `R_g`, `BCD_g`, `MIR`, `PRD`, `APC`} on `APM_intrinsic` | H18 |
| `CHEMCHAIN_score` | Serial-mediation score | H16 |
| `COMB_HLA` | PC1 of locus-β in HLA region | H19 |
| `ARCH` | RPPA block-vector cluster label | H10, H11 |
| `COHER` | Pathway coherence | H44 |
| `BOTTLE` | `min` rule output | H38 |

---

## L3 — Hypothesis engine

A single mediation / interaction engine routes all hypotheses. Signatures:

```
fit_regression(data, y, x, covariates, interactions=None) → coef + CI + diagnostics
fit_mediation(data, x, m, y, covariates, moderator=None, n_boot=2000) → ACME, ADE, prop_mediated, CI
fit_serial_mediation(data, chain=[x, m1, m2, ..., y], covariates) → per-link indirect, total indirect, CI
fit_survival(data, score, covariates, endpoint, treatment=None) → HR + CI + Imai sensitivity
```

Each hypothesis file (`hypotheses/hXX_*.py`) is a **column-mapping script** plus an interpretation hook. No bespoke stats per hypothesis.

---

## L4 — Reports

One generator per figure. All figures in `research_plan/05_timeline.md` are mapped to a generator name. Figure panels regenerate from frozen feature tables; a change upstream invalidates a cache flag.

---

## Feature contract rules (to freeze first)

1. **Sample key**: use `sample` = `TCGA-XX-YYYY-SS` (per Atlas §External sample-ID normalization). Every L1 builder emits this column.
2. **Gene key**: use `gene_name` (HUGO symbol as in `PRIMARY_GENES`).
3. **Missingness policy**: per feature, encode as `NaN` and provide a companion `*_is_imputed` flag when imputation was applied (especially for consensus priors).
4. **Units**: declare units per column in the builder docstring; conform to Atlas units (β 0–1, log2(TPM+1) for expression, log-CPM-QN for ATAC, z-score for RPPA).
5. **Caching**: feature tables write to `analysis/output/features/<feature>.parquet`; every builder reads the latest cache unless `--refresh` is passed.

Until these rules are codified, downstream hypothesis code is fragile. **Freezing this contract is the first item on the timeline.**
