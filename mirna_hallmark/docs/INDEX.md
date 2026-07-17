# Docs index — `mirna_hallmark`

Map of every document under `mirna_hallmark/docs/` (and the notable markdown that
lives outside it). This is the human-facing entry point; [AGENTS.md](../AGENTS.md)
is the agent-onboarding entry point. Links are relative markdown so they render on
GitHub **and** resolve in Obsidian (vault set to relative-link mode).

**Status vocabulary** used below:

| status | meaning |
|--------|---------|
| `active` | current canonical reference or living result doc — trust it |
| `living` | continuously updated (ledger / registry / what's-next) |
| `parked` | design/plan intentionally not built yet — do **not** auto-run |
| `derived` | generated rendition of a canonical doc (paper / slides / PDF) |
| `historical` | one-shot session or handoff, work now **consumed** — read for provenance, not current state |
| `reference` | external literature / data catalog |

`refs` = how many other markdown files in the subproject link to it (in-degree).
Orphans (`refs 0`) are flagged.

---

## 1. Canonical references — read these first

The load-bearing spine. Every number, formula, and component traces back here.

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [MODELING_FRAMEWORK.md](MODELING_FRAMEWORK.md) | active | 7 | Conceptual **and** implemented framework — the connective layer between concepts, formulas, and code. Start here for "what are we modeling and why". |
| [FORMULAS.md](FORMULAS.md) | active | 8 | Code-faithful exact formula for every quantity (pressure / share / specificity / coupling / family / validation-priority). The *how*. |
| [METHODS.md](METHODS.md) | active | 8 | Methods narrative — data → pressure → gate → coupling. |
| [DATA_SOURCES.md](DATA_SOURCES.md) | active | 4 | Every dataset the subproject reads, what it is, and its preprocessing basis. |
| [MIRNA_GENOMIC_CONTEXT_AXIS.md](MIRNA_GENOMIC_CONTEXT_AXIS.md) | active | 4 | miRNA locus genomic-context axis (MH-101): strand-aware host classification (47% coding-host), arm→locus GENCODE canonicalization, the `.N` universe-redefinition, host-target relatedness, + the open downstream unification of MH-99/MH-100 under a per-edge host-gene lens. |
| [ANALYSES_CATALOG.md](ANALYSES_CATALOG.md) | living | 5 | One row per analysis component (module, purpose, inputs, outputs). Orient here before grepping. |
| [MODULE_MAP.md](MODULE_MAP.md) | active | 1 | Post-reorg navigation (top-level 160 → 36 `.py`). Which module does what. |
| [EDGE_QUESTION_TAXONOMY.md](EDGE_QUESTION_TAXONOMY.md) | active | 7 | What a single miRNA→target edge can be asked. Backs the `apm-edge-question` skill. |
| [GENE_QUESTION_TAXONOMY.md](GENE_QUESTION_TAXONOMY.md) | active | 9 | A gene's total incoming miRNA regulation. Backs the `apm-gene-question` skill. |

## 2. Living operational docs

Updated continuously — the current state of the project lives here.

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [DISCOVERY_REGISTRY.md](DISCOVERY_REGISTRY.md) | living | 12 | Every finding (MH-*) with evidence tier + source table. The record of what's established. |
| [ANALYSIS_RUN_LEDGER.md](ANALYSIS_RUN_LEDGER.md) | living | 5 | One row per run: when, status, runtime, catalog alignment. Maintained via `apm-analysis-ledger`. |
| [WHATS_NEXT.md](WHATS_NEXT.md) | living | 5 | Prioritized forward TODO. |
| [PROGRAM_FORWARD_BOARD.md](PROGRAM_FORWARD_BOARD.md) | living | 1 | **Consolidated** open-work board across the WHOLE subproject (fusion/Bayes · CN · progression · discovery · validation · infra) — one place to see everything; points into WHATS_NEXT / LEARNED_MODEL_WHATS_NEXT / CHANNEL_FUSION for depth. |
| [ISOMIR_AWARE_MODELING.md](ISOMIR_AWARE_MODELING.md) | active | 1 | Design: per-sample seed-resolved family predictors (isomiR-corrected X_fam) + propagation to CN/coupling; Phase 0 built (MH-96 per-sample QC), 1-4 planned. |
| [handoffs/HANDOFF_PROTEIN_AND_COMPOSITION.md](handoffs/HANDOFF_PROTEIN_AND_COMPOSITION.md) | **active** | 1 | **⭐ START HERE for the protein axis + composition machinery (session 2026-07-12, MH-103…113).** Verified state, the 5 approved follow-ups, live outputs — and **8 SELF-RETRACTIONS listed up front so they are not re-derived** (βᵗ was a LEAK; 'discordance forces HMC' false; 'protein 3× more compartment-exposed' a unit+cohort artifact; 'deconv=True everywhere' would IMPOSE an over-control; …). Key yield: `learned/readouts.py` (the FOUR readouts genome-wide + the composition RETENTION tag, validated at p=9e-12) and the **confounder-ARCHITECTURE principle** (proliferation/host localise to an ARM = LOCUS properties; composition localises to the FAMILY = an EXPRESSION property ⇒ the unit predicts which tool can fix it). |
| [CPTAC_PROTEIN_CHANNEL_PLAN.md](CPTAC_PROTEIN_CHANNEL_PLAN.md) | living | 1 | **⛔ CENTREPIECE FALSIFIED (MH-103, see its §0a) — `βᵗ` is NOT supported at n=101 and the prior evidence was a LEAK; four leak-free framings agree, incl. the max-power aggregate (BH q<0.10 in 1/17 genes; PTEN dead). Coherent with MH-34: the translational residual is ~1% of genes, below a 101-patient cohort's resolution. SURVIVING: `locus_cn_cptac` (BUILT, r=0.997) · δ-transportability · Bar-5 (validates β) · the unbiased `a`.** Also established (and unaffected): protein carries only **≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier ****≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**)** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**) of the mRNA channel's Fisher information about β (≤7.6% at any `a_g`) ⇒ never a coupling lever either; the miRNA CONFOUNDS the mRNA→protein slope (marginal `a` biased up to 38%); **ENGINE VERDICT: Gibbs, no HMC** — corrects CHANNEL_FUSION §7/§L/§J's "discordance forces HMC". |

## 3. Result reports — established findings

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [REPORT.md](REPORT.md) | active | 10 | Detailed results report (TCGA-BRCA), pressure spine refreshed 2026-06-17. |
| [BIOLOGICAL_SYNTHESIS.md](BIOLOGICAL_SYNTHESIS.md) | active | 9 | Biology-first distilled surface of what's established. |
| [LANDSCAPE_REPORT.md](LANDSCAPE_REPORT.md) | active | 5 | Cross-state pressure landscape: GTEx-healthy → TCGA-NAT → tumor. |
| [MIRNA_CNV_DOSAGE_REPORT.md](MIRNA_CNV_DOSAGE_REPORT.md) | active | 6 | miRNA CNV dosage (universe copy-number → expression concordance). |
| [MH124_ANTICOUPLING_VALIDITY.md](MH124_ANTICOUPLING_VALIDITY.md) | active | 2 | **MH-124: are anti-couplings meaningful?** Edge EXISTENCE (decoy-controlled + CN-instrumented) vs WEIGHT IMPORTANCE (attribution at chance) vs the S+L refutation. ⚠ **§5's "exogenous CN instrument" existence leg is DOWNGRADED by MH-126** — read that first. |
| [MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md](MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md) | active | 1 | **MH-126: aneuploidy in the CN instrument · what the learned model is FOR · the GRADED attribution question.** Ploidy is a **calibration** fix (impossible-edge negative reduced form 86.4% → 50.0%), **not** a power fix; **⛔ `pi_causal` is NOT exogenous**; β is de-confounded + transferable but **cannot name** the canonical family (rank = chance, abundance beats it) though it does carry a weak, family-FE-borderline **evidence gradient**. Carries the **w-circularity gate** (β is bit-identical under shuffled w). |
| [MH127_DECOY_MODEL_GENE_BUDGET.md](MH127_DECOY_MODEL_GENE_BUDGET.md) | active | 1 | **MH-127: the DECOY MODEL — is the learning doing anything beyond FITTING?** Learn β on a MATCHED FAKE (abundance-matched, site-free, non-family, non-cluster, non-collinear) regulator set and compare the gene budget. **In-cohort: NO** — the fitted decoy reaches **79–91%** of the real model's OOF field (real-vs-fake **p=0.21/0.58, n.s.**) and **beats real unfitted abundance** ⇒ **MH-115's headline is RESTATED (a FITTING claim, not curation)**. **What survives: OUT-OF-COHORT TRANSFER** (frozen β predict CPTAC mRNA, the decoy's do not; paired p=0.030/0.027) and a **fit-free high-evidence-weight tertile** (q=0.006, FAKE1 only). The stratified rescue (regulator count · Shapley class) **FAILS**. ⚠ **Its stratified numbers and its CPTAC-protein null are AMENDED/CORRECTED by [MH130_GENE_LENS.md](MH130_GENE_LENS.md).** |
| [MH130_GENE_LENS.md](MH130_GENE_LENS.md) | active | 1 | **MH-130/131/132: THE GENE LENS — the per-gene competence map.** Answers, per gene and **a priori**: can this gene's miRNA regulation be **measured at all** (**17.6% cannot**; 48.1% have ONE seed family, where β ≡ uniform exactly)? does the **curated edge set** beat a matched fake there? does it **transfer**? **The 78.9% real-vs-fake "tie" is NOISE-FLOOR DILUTION** (at the floor the decoy *beats* the real model). **THE MODEL'S HONEST DOMAIN OF VALIDITY — fit-free, 27% of the universe: `n_fam ≥ 3` AND `w_max > median`** (gap −0.0376/−0.0282 in both fake sets; **complement +0.0002, p=0.72**). Also: the **SNR strata RETRACTED** (double-log bug), the **CPTAC-protein discrepancy SETTLED** (it was MH-127's decoy trim) — **but a site-free fitted fake reaches 99% of the protein transfer**, so that cell is *not* curation evidence. Artifacts: `output/learned/gene_atlas.tsv`, `output/learned/gene_competence_map.tsv`. |

## 4. Framework renditions (derived from §1)

Generated presentations of [MODELING_FRAMEWORK.md](MODELING_FRAMEWORK.md) — keep in
sync with the canonical, don't edit independently.

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [MODELING_FRAMEWORK_EXTERNAL.md](MODELING_FRAMEWORK_EXTERNAL.md) | derived | 2 | External / paper version (also `.pdf`, `.docx`) — abstract-led, NotebookLM source. |
| [MODELING_FRAMEWORK_SLIDES_COMPANION.md](MODELING_FRAMEWORK_SLIDES_COMPANION.md) | derived | **0** | Per-slide build spec + worked numbers for the NotebookLM deck from _EXTERNAL. Orphan — reachable only from here. |

## 5. Parked design docs — not built, do not auto-run

Intentional forward designs. Each self-declares parked status.

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md](ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md) | parked | 2 | Identity (who) vs magnitude (force) attribution frame + implementation TODOs. |
| [ATTRIBUTION_CONTEXT_AXIS.md](ATTRIBUTION_CONTEXT_AXIS.md) | parked | 1 | Bridge from edge-discovery/coupling to a context axis. |
| [PRESSURE_FUTURE_OPTIONS.md](PRESSURE_FUTURE_OPTIONS.md) | parked | 2 | Deferred pressure-formula options A–E (evidence-scoring companion). |
| [PRESSURE_PROGNOSTIC_DESIGN.md](PRESSURE_PROGNOSTIC_DESIGN.md) | parked | 1 | Parsimonious functional miRNA-pressure prognostic design (needs METABRIC-full, EGA-gated). |
| [PRIMARY_TUMOR_ROADMAP.md](PRIMARY_TUMOR_ROADMAP.md) | parked | 1 | Applying the DCIS/EV rigor lens to the primary-tumor arc. |
| [DUAL_SPINE_COMPARISON_PLAN.md](DUAL_SPINE_COMPARISON_PLAN.md) | parked→executed | **0** | miRTar M0 vs ENCORI M0′ plan. **Executed** — outcome is `output/dual_spine_comparison/CANONICAL_DECISION.md`. Orphan; supersede or link. |

## 6. Learned-model cluster (BUILT — production, 2026-07-06)

> **⚠ CURRENT MODEL STATE — read this ONE line before reasoning about "the learned model", and trust nothing that
> contradicts it (2026-07 convergence, MH-85):** the model is **ONE dense learned-τ² non-negative Bayesian posterior
> per gene, two readouts** — π≡1 dense → coupling/attribution/identifiability; evidence-π → discovery — and **the
> adaptive lasso is RETIRED to a baseline.** Canonical statement: **[`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md §6b`](LEARNED_MODEL_DISCOVERY_SYNTHESIS.md)** (+ §0 thesis).
> **Any text describing lasso-as-primary, "Bayes = uncertainty layer only", or the Decision A–J / Bars / "no-payoff"
> gap ledger is PRE-CONVERGENCE (≤2026-07-05) and SUPERSEDED** — that includes `LEARNED_MODEL_DESIGN_RESPONSE.md`
> (historical) and `LEARNED_MODEL_WHATS_NEXT.md §"DESIGN vs AS-BUILT"`. Never ground the *current* model on a
> `historical`-tagged doc or a pre-convergence section; if a section's claim about "what the model is" is undated or
> older than §6b, it does not describe the current model. (Guard added 2026-07-10 after a grounding error.)

The learned miRNA→gene regulatory matrix is **implemented and in production** in `learned/` (not parked). The
canonical trio is **ESTIMATOR_MAP** (what model for what), **LEARNED_MODEL_METHODS** (formulas), **VALIDATION**
(does it work); **DISCOVERY_SYNTHESIS** is the capstone (§6b = the converged model in one line); **WHATS_NEXT** is
the living TODO. The three design-phase docs are now **historical** (their work is built — read for provenance).
Inclusion migrated to POOLED-HE 2026-07-06. *Consolidation candidate: fold LEARNED_MODEL_METHODS formulas into
FORMULAS.md and promote the trio into §1 canonical, so the learned model isn't a separate "cluster" now that it's core.*

| doc | status | refs | role in cluster |
|-----|--------|-----:|------|
| [LEARNED_MODEL_ESTIMATOR_MAP.md](LEARNED_MODEL_ESTIMATOR_MAP.md) | active | 2 | **Canonical map** — what estimator for what job, and why over alternatives. Current (pooled-HE, loss lens). |
| [LEARNED_MODEL_METHODS.md](LEARNED_MODEL_METHODS.md) | active | 1 | Formula-level spec of every estimator (§0–§17). Current. |
| [LEARNED_MODEL_RATIONALE.md](LEARNED_MODEL_RATIONALE.md) | active | 1 | **Conceptual companion** — the *why* behind every §of METHODS (why this construction, what breaks otherwise). §-aligned to METHODS. Current. |
| [LEARNED_MODEL_DISCOVERY_SYNTHESIS.md](LEARNED_MODEL_DISCOVERY_SYNTHESIS.md) | active | 1 | **Capstone synthesis** — the whole learned-τ² + discovery + retrofit + dossier architecture in one doc (one fit → 4 jobs; deconv handling; candidacy TS∪K_D; 3-layer validation → 268 triple-validated edges). MH-82/83/84. |
| [LEARNED_MODEL_VALIDATION.md](LEARNED_MODEL_VALIDATION.md) | active | 1 | Validation dossier — numbers pending pooled-HE regen. |
| [LEARNED_MODEL_WHATS_NEXT.md](LEARNED_MODEL_WHATS_NEXT.md) | living | 1 | Cluster-local living TODO. |
| [LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md](LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md) | parked | 1 | **Design (no code yet)** — the broader object §6b is a slice of: one latent M fused from many observation *channels*; CN channel = Gap-2B Bayes (between-family exclusion + within-family δ as one posterior, a conjugate Gaussian term on β). Grounded in `_gibbs_posterior`. |
| [LEARNED_MODEL_STATE_CHANNEL_PLAN.md](LEARNED_MODEL_STATE_CHANNEL_PLAN.md) | active | 1 | **PLAN (awaiting review) — the STATE/progression channel `M_T = a·M_H + Δ`.** Carries the SETTLED deconvolution/confounder-panel decision (β is panel-invariant, ρ=0.94 across atlases; the CIBERSORTx "reconciliation blocker" is RETIRED — zero new runs) and the real cross-state problem (the matched-C policy; C's tumour-only terms). Binds the CPTAC session's C too. |
| [LEARNED_REGULATORY_MATRIX_DESIGN.md](LEARNED_REGULATORY_MATRIX_DESIGN.md) | active | 2 | Design note on learning M (companion to MODELING_FRAMEWORK). |
| [LEARNED_MODEL_DESIGN_RESPONSE.md](LEARNED_MODEL_DESIGN_RESPONSE.md) | historical | 3 | Design report (decision log). Work now built — provenance. |
| [LEARNED_MODEL_DISCUSSION_PROMPT.md](LEARNED_MODEL_DISCUSSION_PROMPT.md) | historical | 1 | The framing prompt the design answered. Consumed. |
| [LEARNED_MODEL_BUILD_PLAN.md](LEARNED_MODEL_BUILD_PLAN.md) | historical | 3 | Reuse-vs-rebuild + phased plan. Built. |

## 7. Historical — consumed session artifacts & handoffs

One-shot docs whose work is **done and folded into the registry/reports**. Read for
provenance; they do **not** reflect current state. Candidates to move to a
`docs/archive/` folder (see index footer).

| doc | status | refs | what it captured (now consumed) |
|-----|--------|-----:|------|
| [HANDOFF_NORMAL_EXCLUDED_BATCH_RERUN.md](HANDOFF_NORMAL_EXCLUDED_BATCH_RERUN.md) | historical | 1 | Normal-excluded + batch-corrected rerun — **done** (registry MH-53). |
| [ORPHAN_DISCOVERY_HANDOFF.md](ORPHAN_DISCOVERY_HANDOFF.md) | historical | **0** | CPTAC orphan-edge discovery → next-analysis handoff — consumed (orphan eCLIP validation done). |
| [CPTAC_VALIDATION_SESSION_SUMMARY.md](CPTAC_VALIDATION_SESSION_SUMMARY.md) | historical | 1 | CPTAC protein validation session summary — consumed. |
| [DCIS_EV_SYNTHESIS.md](DCIS_EV_SYNTHESIS.md) | historical | 1 | DCIS / pre-malignant / EV synthesis — done (MH-48..54). |

## 8. Reference scans

| doc | status | refs | what it is |
|-----|--------|-----:|------------|
| [CRC_PORT_LITERATURE_SCAN.md](CRC_PORT_LITERATURE_SCAN.md) | reference | **0** | CRC-port literature scan (also `.pdf`). Orphan — link from WHATS_NEXT or a port-planning doc. |

---

## Markdown outside `docs/`

Not part of this index's link graph, but part of the subproject knowledge base:

- **Subproject roots:** [AGENTS.md](../AGENTS.md) (agent onboarding), [README.md](../README.md), [CLAUDE.md](../CLAUDE.md)
- **method_dev/** (method-development notes):
  [aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md](../method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md),
  [aggregate_pressure/AGGREGATE_PRESSURE_FINDINGS.md](../method_dev/aggregate_pressure/AGGREGATE_PRESSURE_FINDINGS.md),
  [aggregate_pressure/FOLLOWUPS.md](../method_dev/aggregate_pressure/FOLLOWUPS.md),
  [arm_expression/ARM_EXPRESSION_FLOOR.md](../method_dev/arm_expression/ARM_EXPRESSION_FLOOR.md),
  [arm_expression/SILENT_ARM_REMOVAL.md](../method_dev/arm_expression/SILENT_ARM_REMOVAL.md),
  [site_ladder/SITE_FILTER_LADDER_PLAN.md](../method_dev/site_ladder/SITE_FILTER_LADDER_PLAN.md),
  [method_dev/README.md](../method_dev/README.md)
- **output/** (per-run decisions/validation):
  [dual_spine_comparison/CANONICAL_DECISION.md](../output/dual_spine_comparison/CANONICAL_DECISION.md),
  [encori_share_sensitivity/VALIDATION.md](../output/encori_share_sensitivity/VALIDATION.md)
- **Grant/**: [DRAFT_proposal_mirna_hallmark.md](../Grant/DRAFT_proposal_mirna_hallmark.md), [DRAFT_seed_grant_mirna_subtype.md](../Grant/DRAFT_seed_grant_mirna_subtype.md)
- **presentation/**: [presentation/README.md](../presentation/README.md)

---

## Housekeeping flags (as of this index build, 2026-07-05)

- **Broken link in [AGENTS.md](../AGENTS.md):** the parent-repo reference
  `analysis/COHORT_AND_JOIN_CONVENTIONS.md` was missing its `docs/` segment (file is at
  `analysis/docs/COHORT_AND_JOIN_CONVENTIONS.md`). *(Fixed in the pass that created this index.)*
- **Orphans** (`refs 0`, reachable only from this index): `MODELING_FRAMEWORK_SLIDES_COMPANION`,
  `DUAL_SPINE_COMPARISON_PLAN` (executed), `ORPHAN_DISCOVERY_HANDOFF` (consumed),
  `CRC_PORT_LITERATURE_SCAN`.
- **Proposed archive:** the §7 historical docs could move to `docs/archive/` to shrink
  the active surface. Not done automatically — pending decision.
