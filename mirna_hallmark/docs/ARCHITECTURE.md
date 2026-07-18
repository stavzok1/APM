# Canonical architecture — axes × models × analyses × results

> **GENERATED** by `analyses/ops/gen_architecture.py` — do NOT hand-edit. A materialized VIEW over the
> one-home docs (registry, catalog, STATE_OF_PLAY): every finding is an `MH-##` (→ registry), every
> module a path. Regenerate after tagging changes. Join key: `docs/derived/axis_assignment.tsv`.
> Axes: **12** · findings mapped: **177** · modules tagged: **254**.

## The learned model  `model`
*what the estimator IS*

- **model:** ONE dense learned-τ² non-negative Gibbs posterior/gene; two readouts (π≡1 dense→coupling; evidence-π→discovery). Lasso RETIRED.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **12** findings · **30** modules
- **results (12):** [MH-102f](DISCOVERY_REGISTRY.md) · [MH-102e](DISCOVERY_REGISTRY.md) · [MH-102c](DISCOVERY_REGISTRY.md) · [MH-93](DISCOVERY_REGISTRY.md) · [MH-92](DISCOVERY_REGISTRY.md) · [MH-30](DISCOVERY_REGISTRY.md) · [MH-32](DISCOVERY_REGISTRY.md) · [MH-124b](DISCOVERY_REGISTRY.md) · [MH-126c](DISCOVERY_REGISTRY.md) · [MH-141](DISCOVERY_REGISTRY.md) · [MH-145](DISCOVERY_REGISTRY.md) · [MH-151](DISCOVERY_REGISTRY.md)
- **analyses (30):** `analyses/cptac/composition_pinpoint_lineage.py`, `analyses/outcome/outcome_famous_compare.py`, `analyses/outcome/outcome_mir200_deepdive.py`, `analyses/spatial/spatial_xenium.py`, `eval/_niter_bias_test.py`, `eval/_niter_floor_test.py`, `eval/_share_drift_check.py`, `eval/grid_oof.py`, `eval/within_family_showdown.py`, `learned/analyses/eb_shrink.py`, `learned/analyses/shrinkage_compare.py`, `learned/attribution_eb.py`, `learned/calibration.py`, `learned/card.py` …

## Edge existence & coupling  `edge-existence`
*do the curated/predicted edges exist AND act (coupling|C)?*

- **model:** partial-Spearman(arm, mRNA | C) via the Gibbs posterior; C = CPE+target_cn+mal_prolif (+8 deconv fractions).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **29** findings · **43** modules
- **results (29):** [MH-101](DISCOVERY_REGISTRY.md) · [MH-100](DISCOVERY_REGISTRY.md) · [MH-80](DISCOVERY_REGISTRY.md) · [MH-78](DISCOVERY_REGISTRY.md) · [MH-71](DISCOVERY_REGISTRY.md) · [MH-1](DISCOVERY_REGISTRY.md) · [MH-2](DISCOVERY_REGISTRY.md) · [MH-3](DISCOVERY_REGISTRY.md) · [MH-6](DISCOVERY_REGISTRY.md) · [MH-7](DISCOVERY_REGISTRY.md) · [MH-17](DISCOVERY_REGISTRY.md) · [MH-22](DISCOVERY_REGISTRY.md) · [MH-44](DISCOVERY_REGISTRY.md) · [MH-45](DISCOVERY_REGISTRY.md) · [MH-46](DISCOVERY_REGISTRY.md) · [MH-53](DISCOVERY_REGISTRY.md) · [MH-59](DISCOVERY_REGISTRY.md) · [MH-57](DISCOVERY_REGISTRY.md) · [MH-47](DISCOVERY_REGISTRY.md) · [MH-H1](DISCOVERY_REGISTRY.md) · [MH-H2](DISCOVERY_REGISTRY.md) · [MH-H3](DISCOVERY_REGISTRY.md) · [MH-111](DISCOVERY_REGISTRY.md) · [MH-112](DISCOVERY_REGISTRY.md) · [MH-122](DISCOVERY_REGISTRY.md) · [MH-130e](DISCOVERY_REGISTRY.md) · [MH-132b](DISCOVERY_REGISTRY.md) · [MH-136](DISCOVERY_REGISTRY.md) · [MH-154](DISCOVERY_REGISTRY.md)
- **analyses (43):** `ago_gate.py`, `analyses/architecture/geneset_architecture.py`, `analyses/cnv_locus/sc_mirna_target_k562.py`, `analyses/edge_panels/edge_breast_context.py`, `analyses/edge_panels/edge_partial_corr_panels.py`, `analyses/edge_panels/edge_pressure_vs_corr_scatter.py`, `analyses/evidence_dev/evidence_scoring_sensitivity.py`, `analyses/misc/_oof_arm_vs_family_eiv.py`, `analyses/misc/compartment_annotate_network.py`, `analyses/misc/decoupling_validation.py`, `analyses/misc/escape_mechanism.py`, `analyses/misc/gene_force_concentration.py`, `analyses/misc/target_combined_anticorr.py`, `analyses/nmf/nmf_advanced.py` …

## CN-dose causal  `cn-causal`
*CN-dose causal identification (instrument, exclusion)*

- **model:** CN-locus 2SLS instrument + Hansen-J over-ID + T1 screen. ⛔ instruments RETRACTED; DESIGN live (highest-value open item).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **15** findings · **4** modules
- **results (15):** [MH-99](DISCOVERY_REGISTRY.md) · [MH-91](DISCOVERY_REGISTRY.md) · [MH-89](DISCOVERY_REGISTRY.md) · [MH-87](DISCOVERY_REGISTRY.md) · [MH-4](DISCOVERY_REGISTRY.md) · [MH-5](DISCOVERY_REGISTRY.md) · [MH-28](DISCOVERY_REGISTRY.md) · [MH-29](DISCOVERY_REGISTRY.md) · [MH-43](DISCOVERY_REGISTRY.md) · [MH-119/120](DISCOVERY_REGISTRY.md) · [MH-126](DISCOVERY_REGISTRY.md) · [MH-124r](DISCOVERY_REGISTRY.md) · [MH-128c](DISCOVERY_REGISTRY.md) · [MH-129](DISCOVERY_REGISTRY.md) · [MH-133](DISCOVERY_REGISTRY.md)
- **analyses (4):** `analyses/cnv_locus/cn_dosage_attribution.py`, `analyses/cnv_locus/mirna_cnv_subtype_depth.py`, `eval/within_family_cn_iv.py`, `mirna_locus_cnv.py`

## Attribution / identity  `attribution`
*WHO owns a gene's regulation — identity vs magnitude?*

- **model:** Shapley/LMG on R² (identity, collinearity-fair) vs Gibbs β (magnitude). `share`=β_f/Σβ is NOT identity.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **17** findings · **7** modules
- **results (17):** [MH-97](DISCOVERY_REGISTRY.md) · [MH-96](DISCOVERY_REGISTRY.md) · [MH-95](DISCOVERY_REGISTRY.md) · [MH-94](DISCOVERY_REGISTRY.md) · [MH-88](DISCOVERY_REGISTRY.md) · [MH-86](DISCOVERY_REGISTRY.md) · [MH-41](DISCOVERY_REGISTRY.md) · [MH-42](DISCOVERY_REGISTRY.md) · [MH-110](DISCOVERY_REGISTRY.md) · [MH-118](DISCOVERY_REGISTRY.md) · [MH-121](DISCOVERY_REGISTRY.md) · [MH-127b](DISCOVERY_REGISTRY.md) · [MH-138](DISCOVERY_REGISTRY.md) · [MH-140](DISCOVERY_REGISTRY.md) · [MH-146](DISCOVERY_REGISTRY.md) · [MH-150](DISCOVERY_REGISTRY.md) · [MH-152](DISCOVERY_REGISTRY.md)
- **analyses (7):** `analyses/misc/genomewide_promiscuity.py`, `analyses/nmf/within_pathway_nmf.py`, `analyses/pressure_dev/denominator_attribution_sweep.py`, `analyses/pressure_dev/pressure_attribution_validation.py`, `eval/attribution_showdown.py`, `learned/analyses/attribution_identity.py`, `learned/attribution.py`

## Decoy / specificity  `decoy`
*does curation beat an abundance/variance-matched fake?*

- **model:** site-free / abundance+variance-matched decoy; Hungarian on signed loadings+dose+variance. Gap ≈ −0.012 (two designs agree).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **22** findings · **5** modules
- **results (22):** [MH-65](DISCOVERY_REGISTRY.md) · [MH-40](DISCOVERY_REGISTRY.md) · [MH-124](DISCOVERY_REGISTRY.md) · [MH-125](DISCOVERY_REGISTRY.md) · [MH-127](DISCOVERY_REGISTRY.md) · [MH-127c](DISCOVERY_REGISTRY.md) · [MH-128](DISCOVERY_REGISTRY.md) · [MH-128a](DISCOVERY_REGISTRY.md) · [MH-128b](DISCOVERY_REGISTRY.md) · [MH-130](DISCOVERY_REGISTRY.md) · [MH-130a](DISCOVERY_REGISTRY.md) · [MH-130b](DISCOVERY_REGISTRY.md) · [MH-131](DISCOVERY_REGISTRY.md) · [MH-131b](DISCOVERY_REGISTRY.md) · [MH-134](DISCOVERY_REGISTRY.md) · [MH-135](DISCOVERY_REGISTRY.md) · [MH-136b](DISCOVERY_REGISTRY.md) · [MH-139](DISCOVERY_REGISTRY.md) · [MH-144](DISCOVERY_REGISTRY.md) · [MH-147](DISCOVERY_REGISTRY.md) · [MH-148](DISCOVERY_REGISTRY.md) · [MH-149](DISCOVERY_REGISTRY.md)
- **analyses (5):** `analyses/misc/mir301_family_network.py`, `eval/decoy_bench.py`, `learned/analyses/gene_atlas.py`, `learned/seed_rarity.py`, `learned/seq_specificity.py`

## Discovery  `discovery`
*novel edges beyond curation (its own lane)*

- **model:** site-free empirical null (heavy-tailed) → Simes-within-family → BH; per-edge EMPTY, signal is set-level + convergent-evidence.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **11** findings · **16** modules
- **results (11):** [MH-23](DISCOVERY_REGISTRY.md) · [MH-24](DISCOVERY_REGISTRY.md) · [MH-26](DISCOVERY_REGISTRY.md) · [MH-38](DISCOVERY_REGISTRY.md) · [MH-39](DISCOVERY_REGISTRY.md) · [MH-55](DISCOVERY_REGISTRY.md) · [MH-123](DISCOVERY_REGISTRY.md) · [MH-142](DISCOVERY_REGISTRY.md) · [MH-143](DISCOVERY_REGISTRY.md) · [MH-153](DISCOVERY_REGISTRY.md) · [MH-156](DISCOVERY_REGISTRY.md)
- **analyses (16):** `analyses/dcis_ev/dcis_orphan_rigor.py`, `analyses/dcis_ev/ev_mirna_deep.py`, `analyses/dcis_ev/ev_mirna_gse255660_screen.py`, `analyses/evidence_dev/orphan_noncircular_report.py`, `analyses/evidence_dev/orphan_noncircular_test.py`, `analyses/evidence_dev/orphan_noncircular_universe.py`, `analyses/mir301/mir301_family_depth.py`, `analyses/mir301/mir301_focus_genes.py`, `analyses/ops/axis_bootstrap.py`, `analyses/ops/rerun_normal_excluded_batch.py`, `eval/_e6_production_check.py`, `eval/targetscan_orphan_coupling.py`, `learned/analyses/enrichment.py`, `learned/discovery.py` …

## CPTAC / protein  `protein`
*does it hold at the protein layer?*

- **model:** CPTAC prospective (n=101); β_source = a·β_target gauge; composition confound reframe (epithelial survive, stromal collapse).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **20** findings · **22** modules
- **results (20):** [MH-77](DISCOVERY_REGISTRY.md) · [MH-31](DISCOVERY_REGISTRY.md) · [MH-33](DISCOVERY_REGISTRY.md) · [MH-34](DISCOVERY_REGISTRY.md) · [MH-35](DISCOVERY_REGISTRY.md) · [MH-36](DISCOVERY_REGISTRY.md) · [MH-103](DISCOVERY_REGISTRY.md) · [MH-104](DISCOVERY_REGISTRY.md) · [MH-105](DISCOVERY_REGISTRY.md) · [MH-106](DISCOVERY_REGISTRY.md) · [MH-107](DISCOVERY_REGISTRY.md) · [MH-108](DISCOVERY_REGISTRY.md) · [MH-109](DISCOVERY_REGISTRY.md) · [MH-113](DISCOVERY_REGISTRY.md) · [MH-114](DISCOVERY_REGISTRY.md) · [MH-115](DISCOVERY_REGISTRY.md) · [MH-116/117](DISCOVERY_REGISTRY.md) · [MH-126b](DISCOVERY_REGISTRY.md) · [MH-130c](DISCOVERY_REGISTRY.md) · [MH-130d](DISCOVERY_REGISTRY.md)
- **analyses (22):** `analyses/cptac/composition_pinpoint_report.py`, `analyses/cptac/cptac_acquired_pressure.py`, `analyses/cptac/cptac_orphan_confound_pilot.py`, `analyses/cptac/cptac_orphan_discovery.py`, `analyses/cptac/cptac_orphan_evidence_table.py`, `analyses/cptac/cptac_resid_survivors.py`, `analyses/cptac/cptac_target_specificity.py`, `analyses/dcis_ev/dcis_caf_nf_gse37527.py`, `analyses/evidence_dev/fetch_encori_mrna.py`, `analyses/spatial/spatial_mibi_anchor.py`, `eval/buffa_validation.py`, `eval/coupling_grid.py`, `eval/cptac_validation.py`, `eval/grid_full.py` …

## Progression / state  `progression`
*GTEx→NAT→tumor trajectory / state*

- **model:** nested state model M_t = a·M_H + Δ; GTEx-healthy anchor (NAT is field-effect-contaminated).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **16** findings · **17** modules
- **results (16):** [MH-102](DISCOVERY_REGISTRY.md) · [MH-102b](DISCOVERY_REGISTRY.md) · [MH-102g](DISCOVERY_REGISTRY.md) · [MH-102d](DISCOVERY_REGISTRY.md) · [MH-79](DISCOVERY_REGISTRY.md) · [MH-72](DISCOVERY_REGISTRY.md) · [MH-70](DISCOVERY_REGISTRY.md) · [MH-66](DISCOVERY_REGISTRY.md) · [MH-64](DISCOVERY_REGISTRY.md) · [MH-9](DISCOVERY_REGISTRY.md) · [MH-27](DISCOVERY_REGISTRY.md) · [MH-64](DISCOVERY_REGISTRY.md) · [MH-58](DISCOVERY_REGISTRY.md) · [MH-50](DISCOVERY_REGISTRY.md) · [MH-132](DISCOVERY_REGISTRY.md) · [MH-137](DISCOVERY_REGISTRY.md)
- **analyses (17):** `analyses/cross_state/cross_state_coupling.py`, `analyses/cross_state/cross_state_deep_dive.py`, `analyses/cross_state/cross_state_landscape.py`, `analyses/dcis_ev/dcis_ev_progression.py`, `analyses/edge_panels/edge_acquired_pressure_panels.py`, `analyses/edge_panels/edge_transition_pressure_panels.py`, `analyses/misc/mirna_comovement.py`, `analyses/misc/nat_tumor_umap.py`, `analyses/ops/rerun_parallel.py`, `analyses/pressure_dev/pressure_prognostic_signature.py`, `analyses/pressure_dev/program_pressure_by_role.py`, `family_normal_reference.py`, `learned/ago_loading.py`, `learned/analyses/parallel_view.py` …

## Subtype  `subtype`
*PAM50-stratified coupling / who-is-pressured*

- **model:** per-PAM50 coupling test (learned β). ⬜ full subtype contrasts + subtype-stratified discovery = board §Y.2/#1.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **15** findings · **22** modules
- **results (15):** [MH-68](DISCOVERY_REGISTRY.md) · [MH-67](DISCOVERY_REGISTRY.md) · [MH-8](DISCOVERY_REGISTRY.md) · [MH-10](DISCOVERY_REGISTRY.md) · [MH-11](DISCOVERY_REGISTRY.md) · [MH-14](DISCOVERY_REGISTRY.md) · [MH-15](DISCOVERY_REGISTRY.md) · [MH-12](DISCOVERY_REGISTRY.md) · [MH-13](DISCOVERY_REGISTRY.md) · [MH-18](DISCOVERY_REGISTRY.md) · [MH-19](DISCOVERY_REGISTRY.md) · [MH-16](DISCOVERY_REGISTRY.md) · [MH-20](DISCOVERY_REGISTRY.md) · [MH-21](DISCOVERY_REGISTRY.md) · [MH-25](DISCOVERY_REGISTRY.md)
- **analyses (22):** `MIRACLE/make_miracle_figures.py`, `analyses/cnv_locus/ccle_breast_target_anticorr.py`, `analyses/cnv_locus/cell_line_subtype_fidelity.py`, `analyses/cnv_locus/mirna_cnv_genome_maps.py`, `analyses/cnv_locus/mirna_locus_genome_dispersion.py`, `analyses/cross_state/cross_state_expression_panels.py`, `analyses/dcis_ev/ev_mirna_301a_target_depth.py`, `analyses/edge_panels/edge_pressure_panels.py`, `analyses/figures/make_grant_figure.py`, `analyses/figures/make_seed_grant_figures.py`, `analyses/misc/gene_expression_explainability.py`, `analyses/misc/spine_claim_audit.py`, `analyses/nmf/nmf_sample_signatures.py`, `analyses/nmf/within_subtype_signatures.py` …

## Outcome / prognostic  `outcome`
*does it predict survival?*

- **model:** survival on pressure/β. Pressure magnitude non-prognostic (data-bound: TCGA=OS only). ⬜ learned-β prognostic untested.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **2** findings · **8** modules
- **results (2):** [MH-61](DISCOVERY_REGISTRY.md) · [MH-60](DISCOVERY_REGISTRY.md)
- **analyses (8):** `analyses/cnv_locus/mirna_locus_sv_overlap.py`, `analyses/outcome/outcome_extra_axes.py`, `analyses/outcome/outcome_pressure_survival.py`, `analyses/outcome/outcome_pressure_variants.py`, `analyses/outcome/outcome_program_role_topology.py`, `analyses/outcome/outcome_subtype.py`, `analyses/outcome/outcome_survival.py`, `analyses/pressure_dev/pressure_prognostic_gene_centric.py`

## External cohorts  `external`
*independent-cohort replication*

- **model:** Buffa / METABRIC / SCAN-B / CCLE. Cross-cohort loss ~80% (cohort boundary), ~0% (mRNA→protein layer).
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **12** findings · **9** modules
- **results (12):** [MH-76](DISCOVERY_REGISTRY.md) · [MH-75](DISCOVERY_REGISTRY.md) · [MH-74](DISCOVERY_REGISTRY.md) · [MH-73](DISCOVERY_REGISTRY.md) · [MH-69](DISCOVERY_REGISTRY.md) · [MH-68](DISCOVERY_REGISTRY.md) · [MH-67](DISCOVERY_REGISTRY.md) · [MH-66](DISCOVERY_REGISTRY.md) · [MH-65](DISCOVERY_REGISTRY.md) · [MH-37](DISCOVERY_REGISTRY.md) · [MH-63](DISCOVERY_REGISTRY.md) · [MH-62](DISCOVERY_REGISTRY.md)
- **analyses (9):** `analyses/builders/_build_depmap_dependency.py`, `analyses/cnv_locus/ccle_cn_expr_subtype_depth.py`, `analyses/cnv_locus/ccle_mirna_cn_concordance.py`, `analyses/outcome/outcome_buffa_pressure.py`, `analyses/outcome/outcome_escape_signature.py`, `analyses/outcome/outcome_metabric.py`, `analyses/pressure_dev/pressure_prognostic_improve.py`, `presentation/make_figures.py`, `pressure_signature_curated_info.py`

## DCIS / EV / pre-malignant  `dcis-ev`
*the DCIS / extracellular-vesicle lane*

- **model:** in-situ→invasive coupling; selective EV export; MH-55 corroborates MH-114 16 days early.
- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **6** findings · **16** modules
- **results (6):** [MH-56](DISCOVERY_REGISTRY.md) · [MH-54](DISCOVERY_REGISTRY.md) · [MH-52](DISCOVERY_REGISTRY.md) · [MH-51](DISCOVERY_REGISTRY.md) · [MH-49](DISCOVERY_REGISTRY.md) · [MH-48](DISCOVERY_REGISTRY.md)
- **analyses (16):** `analyses/dcis_ev/dcis_caf_tgfb_gse196354.py`, `analyses/dcis_ev/dcis_geo_loader.py`, `analyses/dcis_ev/dcis_htan_compartment.py`, `analyses/dcis_ev/dcis_lcm_compartment_gse162670.py`, `analyses/dcis_ev/dcis_mir29c_gse93740.py`, `analyses/dcis_ev/dcis_mrna_coupling.py`, `analyses/dcis_ev/dcis_risom_outcome_niche.py`, `analyses/dcis_ev/dcis_timing.py`, `analyses/dcis_ev/ev_mirna_followup.py`, `analyses/dcis_ev/ev_mirna_gse270497_ts.py`, `analyses/dcis_ev/ev_mirna_multi_cohort.py`, `analyses/dcis_ev/ev_mirna_replication.py`, `analyses/dcis_ev/ev_mirna_screening.py`, `analyses/dcis_ev/ev_mirna_ts_pressure.py` …

## Shared / infrastructure  `shared`
*axis-agnostic loaders, builders, batch/CN helpers, and the retired-pressure builders — serve every axis*

- **modules (55):** `analyses/builders/_build_meth_chunked.py`, `analyses/builders/_build_tf_census.py`, `analyses/cnv_locus/mirna_locus_chr19_megacluster_gallery.py`, `analyses/cnv_locus/mirna_locus_sv_case_plots.py`, `analyses/edge_panels/edge_prior_refinement.py`, `analyses/evidence_dev/encori_evidence_sensitivity.py`, `analyses/evidence_dev/encori_mirtar_comparison.py`, `analyses/evidence_dev/encori_share_sensitivity.py`, `analyses/misc/_oof_arm_vs_family_eiv2.py`, `analyses/misc/_oof_arm_vs_family_eiv3.py`, `analyses/misc/concordance_target_join.py`, `analyses/misc/genome_wide_promiscuity.py`, `analyses/misc/hub_gene_methylation.py`, `analyses/ops/gen_architecture.py`, `analyses/pressure_dev/pressure_evidence_sensitivity.py`, `analyses/spatial/spatial_coexpression_niche.py`, `analyses/spatial/spatial_deconv_c2l.py`, `analyses/spatial/spatial_visium_public.py`, `arm_expression.py`, `brca_deconvolution.py` …
