# Canonical architecture — axes × models × analyses × results

> **GENERATED** by `analyses/ops/gen_architecture.py` — do NOT hand-edit. A materialized VIEW over the
> one-home docs (registry, catalog, STATE_OF_PLAY): every finding is an `MH-##` (→ registry), every
> module a path. Regenerate after tagging changes. Join key: `docs/derived/axis_assignment.tsv`.
> Axes: **12** · findings mapped: **177** · modules tagged: **223**.

## The learned model  `model`
*what the estimator IS*

- **model:** ONE dense learned-τ² non-negative Gibbs posterior/gene; two readouts (π≡1 dense→coupling; evidence-π→discovery). Lasso RETIRED.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **12** findings · **26** modules
- **results (12):** [MH-102f](DISCOVERY_REGISTRY.md) · [MH-102e](DISCOVERY_REGISTRY.md) · [MH-102c](DISCOVERY_REGISTRY.md) · [MH-93](DISCOVERY_REGISTRY.md) · [MH-92](DISCOVERY_REGISTRY.md) · [MH-30](DISCOVERY_REGISTRY.md) · [MH-32](DISCOVERY_REGISTRY.md) · [MH-124b](DISCOVERY_REGISTRY.md) · [MH-126c](DISCOVERY_REGISTRY.md) · [MH-141](DISCOVERY_REGISTRY.md) · [MH-145](DISCOVERY_REGISTRY.md) · [MH-151](DISCOVERY_REGISTRY.md)
- **analyses (26):** `eval/_niter_floor_test.py`, `eval/_share_drift_check.py`, `eval/_niter_bias_test.py`, `eval/grid_oof.py`, `eval/within_family_showdown.py`, `learned/channel_cn.py`, `learned/attribution_eb.py`, `learned/kd.py`, `learned/gauge.py`, `learned/spike_slab.py`, `learned/occupancy.py`, `learned/readouts.py`, `learned/hierarchical.py`, `learned/families.py` …

## Edge existence & coupling  `edge-existence`
*do the curated/predicted edges exist AND act (coupling|C)?*

- **model:** partial-Spearman(arm, mRNA | C) via the Gibbs posterior; C = CPE+target_cn+mal_prolif (+8 deconv fractions).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **29** findings · **50** modules
- **results (29):** [MH-101](DISCOVERY_REGISTRY.md) · [MH-100](DISCOVERY_REGISTRY.md) · [MH-80](DISCOVERY_REGISTRY.md) · [MH-78](DISCOVERY_REGISTRY.md) · [MH-71](DISCOVERY_REGISTRY.md) · [MH-1](DISCOVERY_REGISTRY.md) · [MH-2](DISCOVERY_REGISTRY.md) · [MH-3](DISCOVERY_REGISTRY.md) · [MH-6](DISCOVERY_REGISTRY.md) · [MH-7](DISCOVERY_REGISTRY.md) · [MH-17](DISCOVERY_REGISTRY.md) · [MH-22](DISCOVERY_REGISTRY.md) · [MH-44](DISCOVERY_REGISTRY.md) · [MH-45](DISCOVERY_REGISTRY.md) · [MH-46](DISCOVERY_REGISTRY.md) · [MH-53](DISCOVERY_REGISTRY.md) · [MH-59](DISCOVERY_REGISTRY.md) · [MH-57](DISCOVERY_REGISTRY.md) · [MH-47](DISCOVERY_REGISTRY.md) · [MH-H1](DISCOVERY_REGISTRY.md) · [MH-H2](DISCOVERY_REGISTRY.md) · [MH-H3](DISCOVERY_REGISTRY.md) · [MH-111](DISCOVERY_REGISTRY.md) · [MH-112](DISCOVERY_REGISTRY.md) · [MH-122](DISCOVERY_REGISTRY.md) · [MH-130e](DISCOVERY_REGISTRY.md) · [MH-132b](DISCOVERY_REGISTRY.md) · [MH-136](DISCOVERY_REGISTRY.md) · [MH-154](DISCOVERY_REGISTRY.md)
- **analyses (50):** `coupling_predictor_comparison.py`, `gene_roles.py`, `brca_deconvolution.py`, `seed_family_justification.py`, `arm_expression.py`, `tcga_batch.py`, `ago_gate.py`, `hallmark_interaction.py`, `coupling_inference.py`, `stats.py`, `eval/held_out_tuning.py`, `eval/core_coupling_deconv_retest.py`, `eval/site_free_null.py`, `eval/core_coupling_composition_retest.py` …

## CN-dose causal  `cn-causal`
*CN-dose causal identification (instrument, exclusion)*

- **model:** CN-locus 2SLS instrument + Hansen-J over-ID + T1 screen. ⛔ instruments RETRACTED; DESIGN live (highest-value open item).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **15** findings · **4** modules
- **results (15):** [MH-99](DISCOVERY_REGISTRY.md) · [MH-91](DISCOVERY_REGISTRY.md) · [MH-89](DISCOVERY_REGISTRY.md) · [MH-87](DISCOVERY_REGISTRY.md) · [MH-4](DISCOVERY_REGISTRY.md) · [MH-5](DISCOVERY_REGISTRY.md) · [MH-28](DISCOVERY_REGISTRY.md) · [MH-29](DISCOVERY_REGISTRY.md) · [MH-43](DISCOVERY_REGISTRY.md) · [MH-119/120](DISCOVERY_REGISTRY.md) · [MH-126](DISCOVERY_REGISTRY.md) · [MH-124r](DISCOVERY_REGISTRY.md) · [MH-128c](DISCOVERY_REGISTRY.md) · [MH-129](DISCOVERY_REGISTRY.md) · [MH-133](DISCOVERY_REGISTRY.md)
- **analyses (4):** `mirna_locus_cnv.py`, `eval/within_family_cn_iv.py`, `analyses/cnv_locus/cn_dosage_attribution.py`, `analyses/cnv_locus/mirna_cnv_subtype_depth.py`

## Attribution / identity  `attribution`
*WHO owns a gene's regulation — identity vs magnitude?*

- **model:** Shapley/LMG on R² (identity, collinearity-fair) vs Gibbs β (magnitude). `share`=β_f/Σβ is NOT identity.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **17** findings · **8** modules
- **results (17):** [MH-97](DISCOVERY_REGISTRY.md) · [MH-96](DISCOVERY_REGISTRY.md) · [MH-95](DISCOVERY_REGISTRY.md) · [MH-94](DISCOVERY_REGISTRY.md) · [MH-88](DISCOVERY_REGISTRY.md) · [MH-86](DISCOVERY_REGISTRY.md) · [MH-41](DISCOVERY_REGISTRY.md) · [MH-42](DISCOVERY_REGISTRY.md) · [MH-110](DISCOVERY_REGISTRY.md) · [MH-118](DISCOVERY_REGISTRY.md) · [MH-121](DISCOVERY_REGISTRY.md) · [MH-127b](DISCOVERY_REGISTRY.md) · [MH-138](DISCOVERY_REGISTRY.md) · [MH-140](DISCOVERY_REGISTRY.md) · [MH-146](DISCOVERY_REGISTRY.md) · [MH-150](DISCOVERY_REGISTRY.md) · [MH-152](DISCOVERY_REGISTRY.md)
- **analyses (8):** `eval/attribution_showdown.py`, `learned/attribution.py`, `learned/analyses/attribution_identity.py`, `output/learned/mh127/code/build_strata.py`, `analyses/nmf/within_pathway_nmf.py`, `analyses/pressure_dev/denominator_attribution_sweep.py`, `analyses/pressure_dev/pressure_attribution_validation.py`, `analyses/misc/genomewide_promiscuity.py`

## Decoy / specificity  `decoy`
*does curation beat an abundance/variance-matched fake?*

- **model:** site-free / abundance+variance-matched decoy; Hungarian on signed loadings+dose+variance. Gap ≈ −0.012 (two designs agree).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **22** findings · **12** modules
- **results (22):** [MH-65](DISCOVERY_REGISTRY.md) · [MH-40](DISCOVERY_REGISTRY.md) · [MH-124](DISCOVERY_REGISTRY.md) · [MH-125](DISCOVERY_REGISTRY.md) · [MH-127](DISCOVERY_REGISTRY.md) · [MH-127c](DISCOVERY_REGISTRY.md) · [MH-128](DISCOVERY_REGISTRY.md) · [MH-128a](DISCOVERY_REGISTRY.md) · [MH-128b](DISCOVERY_REGISTRY.md) · [MH-130](DISCOVERY_REGISTRY.md) · [MH-130a](DISCOVERY_REGISTRY.md) · [MH-130b](DISCOVERY_REGISTRY.md) · [MH-131](DISCOVERY_REGISTRY.md) · [MH-131b](DISCOVERY_REGISTRY.md) · [MH-134](DISCOVERY_REGISTRY.md) · [MH-135](DISCOVERY_REGISTRY.md) · [MH-136b](DISCOVERY_REGISTRY.md) · [MH-139](DISCOVERY_REGISTRY.md) · [MH-144](DISCOVERY_REGISTRY.md) · [MH-147](DISCOVERY_REGISTRY.md) · [MH-148](DISCOVERY_REGISTRY.md) · [MH-149](DISCOVERY_REGISTRY.md)
- **analyses (12):** `eval/decoy_bench.py`, `learned/seed_rarity.py`, `learned/seq_specificity.py`, `learned/analyses/gene_atlas.py`, `output/learned/mh127/code/universe.py`, `output/learned/mh127/code/atk_leak.py`, `output/learned/mh127/code/ladder_oof.py`, `output/learned/mh127/code/build_fakes.py`, `output/learned/mh127/code/stratlens_final.py`, `output/learned/mh127/code/permnull.py`, `output/learned/mh127/code/mh127_synth_cptac.py`, `analyses/misc/mir301_family_network.py`

## Discovery  `discovery`
*novel edges beyond curation (its own lane)*

- **model:** site-free empirical null (heavy-tailed) → Simes-within-family → BH; per-edge EMPTY, signal is set-level + convergent-evidence.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **11** findings · **16** modules
- **results (11):** [MH-23](DISCOVERY_REGISTRY.md) · [MH-24](DISCOVERY_REGISTRY.md) · [MH-26](DISCOVERY_REGISTRY.md) · [MH-38](DISCOVERY_REGISTRY.md) · [MH-39](DISCOVERY_REGISTRY.md) · [MH-55](DISCOVERY_REGISTRY.md) · [MH-123](DISCOVERY_REGISTRY.md) · [MH-142](DISCOVERY_REGISTRY.md) · [MH-143](DISCOVERY_REGISTRY.md) · [MH-153](DISCOVERY_REGISTRY.md) · [MH-156](DISCOVERY_REGISTRY.md)
- **analyses (16):** `eval/_e6_production_check.py`, `eval/targetscan_orphan_coupling.py`, `method_dev/site_ladder/discovery_site_evidence.py`, `learned/discovery.py`, `learned/eval/retrofit.py`, `learned/analyses/enrichment.py`, `analyses/ops/axis_bootstrap.py`, `analyses/ops/rerun_normal_excluded_batch.py`, `analyses/dcis_ev/dcis_orphan_rigor.py`, `analyses/dcis_ev/ev_mirna_gse255660_screen.py`, `analyses/dcis_ev/ev_mirna_deep.py`, `analyses/mir301/mir301_family_depth.py`, `analyses/mir301/mir301_focus_genes.py`, `analyses/evidence_dev/orphan_noncircular_test.py` …

## CPTAC / protein  `protein`
*does it hold at the protein layer?*

- **model:** CPTAC prospective (n=101); β_source = a·β_target gauge; composition confound reframe (epithelial survive, stromal collapse).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **20** findings · **26** modules
- **results (20):** [MH-77](DISCOVERY_REGISTRY.md) · [MH-31](DISCOVERY_REGISTRY.md) · [MH-33](DISCOVERY_REGISTRY.md) · [MH-34](DISCOVERY_REGISTRY.md) · [MH-35](DISCOVERY_REGISTRY.md) · [MH-36](DISCOVERY_REGISTRY.md) · [MH-103](DISCOVERY_REGISTRY.md) · [MH-104](DISCOVERY_REGISTRY.md) · [MH-105](DISCOVERY_REGISTRY.md) · [MH-106](DISCOVERY_REGISTRY.md) · [MH-107](DISCOVERY_REGISTRY.md) · [MH-108](DISCOVERY_REGISTRY.md) · [MH-109](DISCOVERY_REGISTRY.md) · [MH-113](DISCOVERY_REGISTRY.md) · [MH-114](DISCOVERY_REGISTRY.md) · [MH-115](DISCOVERY_REGISTRY.md) · [MH-116/117](DISCOVERY_REGISTRY.md) · [MH-126b](DISCOVERY_REGISTRY.md) · [MH-130c](DISCOVERY_REGISTRY.md) · [MH-130d](DISCOVERY_REGISTRY.md)
- **analyses (26):** `cptac_batch.py`, `evidence_scoring.py`, `eval/buffa_validation.py`, `eval/protein_compartment_null.py`, `eval/grid_full.py`, `eval/coupling_grid.py`, `eval/cptac_validation.py`, `learned/confounders.py`, `learned/cptac_data.py`, `learned/within_family.py`, `learned/genomic_context.py`, `learned/eval/ood_protein.py`, `learned/eval/bayes_parity.py`, `learned/eval/dossier.py` …

## Progression / state  `progression`
*GTEx→NAT→tumor trajectory / state*

- **model:** nested state model M_t = a·M_H + Δ; GTEx-healthy anchor (NAT is field-effect-contaminated).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **16** findings · **22** modules
- **results (16):** [MH-102](DISCOVERY_REGISTRY.md) · [MH-102b](DISCOVERY_REGISTRY.md) · [MH-102g](DISCOVERY_REGISTRY.md) · [MH-102d](DISCOVERY_REGISTRY.md) · [MH-79](DISCOVERY_REGISTRY.md) · [MH-72](DISCOVERY_REGISTRY.md) · [MH-70](DISCOVERY_REGISTRY.md) · [MH-66](DISCOVERY_REGISTRY.md) · [MH-64](DISCOVERY_REGISTRY.md) · [MH-9](DISCOVERY_REGISTRY.md) · [MH-27](DISCOVERY_REGISTRY.md) · [MH-64](DISCOVERY_REGISTRY.md) · [MH-58](DISCOVERY_REGISTRY.md) · [MH-50](DISCOVERY_REGISTRY.md) · [MH-132](DISCOVERY_REGISTRY.md) · [MH-137](DISCOVERY_REGISTRY.md)
- **analyses (22):** `healthy_anchor.py`, `gtex_mirna_matrix.py`, `family_normal_reference.py`, `pressure_engine.py`, `mirna_state_class.py`, `estimate_scores.py`, `method_dev/landscape/cross_state_top50.py`, `learned/states.py`, `learned/card.py`, `learned/ago_loading.py`, `learned/analyses/parallel_view.py`, `analyses/ops/rerun_parallel.py`, `analyses/dcis_ev/dcis_ev_progression.py`, `analyses/edge_panels/edge_acquired_pressure_panels.py` …

## Subtype  `subtype`
*PAM50-stratified coupling / who-is-pressured*

- **model:** per-PAM50 coupling test (learned β). ⬜ full subtype contrasts + subtype-stratified discovery = board §Y.2/#1.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **15** findings · **26** modules
- **results (15):** [MH-68](DISCOVERY_REGISTRY.md) · [MH-67](DISCOVERY_REGISTRY.md) · [MH-8](DISCOVERY_REGISTRY.md) · [MH-10](DISCOVERY_REGISTRY.md) · [MH-11](DISCOVERY_REGISTRY.md) · [MH-14](DISCOVERY_REGISTRY.md) · [MH-15](DISCOVERY_REGISTRY.md) · [MH-12](DISCOVERY_REGISTRY.md) · [MH-13](DISCOVERY_REGISTRY.md) · [MH-18](DISCOVERY_REGISTRY.md) · [MH-19](DISCOVERY_REGISTRY.md) · [MH-16](DISCOVERY_REGISTRY.md) · [MH-20](DISCOVERY_REGISTRY.md) · [MH-21](DISCOVERY_REGISTRY.md) · [MH-25](DISCOVERY_REGISTRY.md)
- **analyses (26):** `hybrid_pressure.py`, `stratum_characterization.py`, `pam50_gene_resolution.py`, `run_all.py`, `data_loaders.py`, `subtype_contrasts.py`, `robustness_checks.py`, `visibility_archetype_contrasts.py`, `MIRACLE/make_miracle_figures.py`, `learned/data.py`, `learned/subtype.py`, `learned/analyses/lineage_verdict.py`, `analyses/figures/make_grant_figure.py`, `analyses/figures/make_seed_grant_figures.py` …

## Outcome / prognostic  `outcome`
*does it predict survival?*

- **model:** survival on pressure/β. Pressure magnitude non-prognostic (data-bound: TCGA=OS only). ⬜ learned-β prognostic untested.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **2** findings · **8** modules
- **results (2):** [MH-61](DISCOVERY_REGISTRY.md) · [MH-60](DISCOVERY_REGISTRY.md)
- **analyses (8):** `analyses/cnv_locus/mirna_locus_sv_overlap.py`, `analyses/pressure_dev/pressure_prognostic_gene_centric.py`, `analyses/outcome/outcome_survival.py`, `analyses/outcome/outcome_program_role_topology.py`, `analyses/outcome/outcome_extra_axes.py`, `analyses/outcome/outcome_pressure_survival.py`, `analyses/outcome/outcome_subtype.py`, `analyses/outcome/outcome_pressure_variants.py`

## External cohorts  `external`
*independent-cohort replication*

- **model:** Buffa / METABRIC / SCAN-B / CCLE. Cross-cohort loss ~80% (cohort boundary), ~0% (mRNA→protein layer).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **12** findings · **9** modules
- **results (12):** [MH-76](DISCOVERY_REGISTRY.md) · [MH-75](DISCOVERY_REGISTRY.md) · [MH-74](DISCOVERY_REGISTRY.md) · [MH-73](DISCOVERY_REGISTRY.md) · [MH-69](DISCOVERY_REGISTRY.md) · [MH-68](DISCOVERY_REGISTRY.md) · [MH-67](DISCOVERY_REGISTRY.md) · [MH-66](DISCOVERY_REGISTRY.md) · [MH-65](DISCOVERY_REGISTRY.md) · [MH-37](DISCOVERY_REGISTRY.md) · [MH-63](DISCOVERY_REGISTRY.md) · [MH-62](DISCOVERY_REGISTRY.md)
- **analyses (9):** `pressure_signature_curated_info.py`, `presentation/make_figures.py`, `analyses/cnv_locus/ccle_cn_expr_subtype_depth.py`, `analyses/cnv_locus/ccle_mirna_cn_concordance.py`, `analyses/builders/_build_depmap_dependency.py`, `analyses/pressure_dev/pressure_prognostic_improve.py`, `analyses/outcome/outcome_escape_signature.py`, `analyses/outcome/outcome_metabric.py`, `analyses/outcome/outcome_buffa_pressure.py`

## DCIS / EV / pre-malignant  `dcis-ev`
*the DCIS / extracellular-vesicle lane*

- **model:** in-situ→invasive coupling; selective EV export; MH-55 corroborates MH-114 16 days early.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **6** findings · **16** modules
- **results (6):** [MH-56](DISCOVERY_REGISTRY.md) · [MH-54](DISCOVERY_REGISTRY.md) · [MH-52](DISCOVERY_REGISTRY.md) · [MH-51](DISCOVERY_REGISTRY.md) · [MH-49](DISCOVERY_REGISTRY.md) · [MH-48](DISCOVERY_REGISTRY.md)
- **analyses (16):** `analyses/dcis_ev/dcis_risom_outcome_niche.py`, `analyses/dcis_ev/ev_mirna_replication.py`, `analyses/dcis_ev/ev_mirna_screening.py`, `analyses/dcis_ev/ev_mirna_gse270497_ts.py`, `analyses/dcis_ev/dcis_caf_tgfb_gse196354.py`, `analyses/dcis_ev/dcis_mrna_coupling.py`, `analyses/dcis_ev/dcis_geo_loader.py`, `analyses/dcis_ev/dcis_mir29c_gse93740.py`, `analyses/dcis_ev/dcis_lcm_compartment_gse162670.py`, `analyses/dcis_ev/ev_mirna_ts_pressure.py`, `analyses/dcis_ev/ev_mirna_multi_cohort.py`, `analyses/dcis_ev/ev_mirna_followup.py`, `analyses/dcis_ev/dcis_timing.py`, `analyses/dcis_ev/dcis_htan_compartment.py` …
