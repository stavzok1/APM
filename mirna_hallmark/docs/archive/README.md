# Archive — consumed session artifacts & superseded design docs

> **Goal:** hold documents whose work is **done and folded into the canonical record**, so they stay
> readable for provenance without sitting on the active surface.
> **What belongs here:** one-shot handoffs/session summaries whose task completed; design docs whose
> design is now built; derived reading-guides for bundles that no longer exist. **NOT** anything a
> reader could mistake for current state.
> **Update trigger:** when a doc's work is consumed, or a design it describes is built/abandoned.
> **Sync-partner:** `../INDEX.md` (which must stop listing an archived doc as active).

**Nothing here describes the current model or the current findings.** For that, read
[`../STATE_OF_PLAY.md`](../STATE_OF_PLAY.md). For the finding-by-finding record, read
[`../DISCOVERY_REGISTRY.md`](../DISCOVERY_REGISTRY.md).

Archived 2026-07-17. All files retain full git history (`git log --follow`).

---

## The `MH12x_*.md` series — the per-finding docs

**The pattern is retired and now banned** (`../../CLAUDE.md` → documentation protocol §3). Six docs / 184 KB,
abandoned at MH-130 while the science ran on to MH-137 — and precisely where the retractions rotted (MH-127
still read as current while MH-137 had superseded its central number four times).
**The registry row is the record of a finding.** Every finding below has one, verified before archiving.

| doc | registry rows that now carry it | what the doc still holds |
|-----|-------------------------------|--------------------------|
| `MH124_ANTICOUPLING_VALIDITY.md` | MH-124 · 124b · **124r** | ⭐ §4b's **estimand argument** (the literature's claim is a *level* claim; β is a *covariance* functional — asking β to recover it is an estimand mismatch, not a model failure). **§5/§5b are VOID** — do not cite. |
| `MH125_WHAT_SURVIVES.md` | **MH-125** *(backfilled 2026-07-17)* | The adversarial re-audit of 19 claims: 1 survives / 6 untouched / 7 retest / 5 die. |
| `MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md` | MH-126 · 126b · 126c | The attribution **rank** result (cite this, not §4b) + the **w-circularity gate**. Ploidy = calibration fix, not power fix. |
| `MH127_DECOY_MODEL_GENE_BUDGET.md` | MH-127 · 127b · 127c | ⛔ Magnitudes **superseded by MH-137**; its "in-cohort tie" verdict is **reversed**. |
| `MH128_DECOY_RESOLUTION_AND_CN_GOLD.md` | **MH-128** *(parent, backfilled)* · 128a · 128b · 128c | Its own phase-3 headline was refuted by its own three adversaries. |
| `MH130_GENE_LENS.md` | **MH-130** *(parent, backfilled)* · 130a–e | The **atlas** half stands (17.6% unmeasurable · 48.1% one-family · the 27% domain). The **gap** half sits on a decoy since discredited. |

⚠ **`MH-129` is an ID COLLISION** — the same id names two unrelated pieces of work (the variance-match gate +
FAKE3 build; and the within-arm two-way-FE CN instrument). Both legs are independently retracted. See registry
row MH-129.

---

## Merged into a survivor (2026-07-17) — content moved, not lost

Each of these had its **live** content folded into the doc that owns that question, then was archived.
The pre-convergence / retired material was deliberately **left behind here** rather than carried forward.

| doc | live content went to | what stayed behind |
|-----|----------------------|--------------------|
| `LEARNED_MODEL_ESTIMATOR_MAP.md` | **`LEARNED_MODEL_METHODS{,_FORMAL}.md` §18** (job → estimator → chosen-over) + the pooled-HE inclusion policy | The whole **pre-convergence** frame: lasso-as-primary, "Hierarchical Bayes = uncertainty layer", Bar-1/Bar-2 comparators, "spike-and-slab loses, lasso ships". ⛔ Also its **row 7 "CN-locus instrument = causality"** — there is no live causality job (MH-124r/126, MH-133). |
| `ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md` | **`LEARNED_MODEL_RATIONALE.md` §15 + §15a** — the identity-vs-magnitude frame, the loud-passenger/quiet-owner cases, and ⭐ **why `share` (β_f/Σβ) is NOT identity** | **§1–§12, the entire softmax era** — retired estimator. |
| `ATTRIBUTION_PRIMITIVE.md` | *(nothing merged)* | The `attribute()` primitive was **never built**; its §5/§8 "live bugs" are themselves stale. |
| `ATTRIBUTION_CONTEXT_AXIS.md` | *(nothing merged)* | Parked in its entirety. |
| `COLLINEARITY_AND_IDENTIFIABILITY.md` | **`LEARNED_MODEL_RATIONALE.md` §4a** — Kind A vs Kind B collinearity, the resolution ladder, the identifiability ceiling as a **limit not a fitting failure** | Its §5 code map (belongs to `ANALYSES_CATALOG.md` under the one-home rule); hierarchical δ-pooling stays **parked**. |
| `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` | **`LEARNED_MODEL_RATIONALE.md`** — ~3 paragraphs: M as a latent parameter, every source a noisy **channel**, and the honesty clause (fusion pays **only** where a channel is non-redundant/exogenous) | The **~60-axis map** (largely unbuilt), the CN channel algebra, ⛔ and its **§7/§L/§J engine claims — "the discordance link forces HMC" is WRONG**; the verdict is Gibbs. |
| `MODULE_MAP.md` | **`ANALYSES_CATALOG.md`** — the orientation section: code tree, the dotted-path run rule, the 35 top-level modules classified spine/baseline/pipeline | Its `learned/` build narrative (**pre-convergence**: "adaptive-lasso" estimator, Bar-3/Bar-5, "Bayes' value = uncertainty, not coupling") and the one-time wave-1/wave-2 migration mechanics. ⚠ Its folder counts were **stale**; the catalog carries freshly measured ones. |

**Not archived, deliberately:** `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` stays active — only its §0/§6b was folded
into `MODELING_FRAMEWORK.md`; its **discovery-lane architecture** (one fit → 4 jobs, candidacy TS∪K_D, the
3-layer validation) is live and has no other home.

---

## Consumed session artifacts & superseded design docs

| doc | why archived |
|-----|--------------|
| `LEARNED_MODEL_DESIGN_RESPONSE.md` | Design-phase report (2026-07-04). **Self-declares HISTORICAL/SUPERSEDED** in its own banner. The model converged after it: one dense learned-τ² posterior, two readouts, lasso retired (`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md §6b`). Its Decision A–J / Bars framing is pre-convergence. |
| `LEARNED_MODEL_DISCUSSION_PROMPT.md` | The framing prompt the design answered. Consumed. |
| `LEARNED_MODEL_BUILD_PLAN.md` | Reuse-vs-rebuild + phased plan (2026-07-04). Built. ⚠ Its one still-live line is recorded in `STATE_OF_PLAY.md`: Bar-5 independent-cohort (Buffa) was **never run on the learned M**. |
| `LEARNED_MODEL_PDF_INDEX.md` | Reading guide for a 4-PDF bundle. The PDFs now live in `../derived/`; the guide indexed a distribution format, not a doc. |
| `CPTAC_VALIDATION_SESSION_SUMMARY.md` | Session summary, 2026-06-23/24 (MH-33→36). Consumed — and its conclusions are **superseded**: the protein validation ran with no composition block (MH-107/114). |
| `ORPHAN_DISCOVERY_HANDOFF.md` | Orphan-edge discovery handoff. Consumed. ⚠ The orphan nomination list it fed collapsed **594 → 23** under composition adjustment (MH-114). |
| `HANDOFF_NORMAL_EXCLUDED_BATCH_RERUN.md` | ⚠ **The doc says "this task is mid-flight … the full rerun did not finish." That is stale.** The ledger records the rerun as **settled 2026-06-26** (MH-53), with refreshed counts throughout. Archived as done; kept because it documents the code changes. |
| `MODELING_FRAMEWORK_SLIDES_COMPANION.md` | Per-slide build spec for the NotebookLM deck derived from `MODELING_FRAMEWORK_EXTERNAL.md`. Derived-of-a-derived; was orphaned (in-degree 0). |
| `DUAL_SPINE_COMPARISON_PLAN.md` | miRTar M0 vs ENCORI M0′ plan. **Executed** — the outcome memo is now tracked at [`../decisions/DUAL_SPINE_CANONICAL_DECISION.md`](../decisions/DUAL_SPINE_CANONICAL_DECISION.md) (it had been stranded inside the gitignored `output/`). The plan's own header still says "Status: planned"; that is stale. |

---

## Parked design docs whose PREMISE died (archived 2026-07-17)

Each was carried as "parked / revisit later". Checked before archiving: in every case the thing they
plan for has since been retired, answered, or blocked on something that was never started. Their
surviving threads are absorbed into `../PROGRAM_FORWARD_BOARD.md` §Z.

| doc | why archived |
|-----|--------------|
| `PRESSURE_FUTURE_OPTIONS.md` | Deferred options A–E for the **pressure formula**. ⛔ **Dead by supersession:** MH-115 **RETIRED the pressure heuristic** — on CPTAC protein it is p=0.11, on CPTAC mRNA p=0.83, while the learned β gives −0.036 (p=5.3e-13). *"The evidence weighting adds nothing over raw abundance, ever."* Tuning a retired estimator is a category error. |
| `LEARNED_REGULATORY_MATRIX_DESIGN.md` | The design note for **learning M**, self-declared *"design / parked. Not in any spine."* **The learned model is production** (`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b). Historical. |
| `PRESSURE_PROGNOSTIC_DESIGN.md` | The prognostic construct. **Blocked on "METABRIC-full (EGA-pending)" — a DAR that was never filed** (no accession, no submission date, no owner anywhere in the repo). And **MH-76 already answered the generalization question**: the frozen panel gives DFI +0.002 / OS +0.006, panel-alone C 0.48–0.53 ≈ random. The last live lead (coherence→OS) is on the board. |
| `PRIMARY_TUMOR_ROADMAP.md` | Take the DCIS/EV rigor lens to primary tumours. **Overtaken by events** — its central ask (the composition retest) *happened*: MH-107/114 ran it, found ~57% of the protein validation is compartment arithmetic, and proved it with an orientation-stratified shuffled null. Its one live thread — the **definitive orphan retest at n=1,041** — is on the board §Z. |
| `CRC_PORT_LITERATURE_SCAN.md` | Literature landscape for porting the machine to colorectal cancer (2026-06-28). ⛔ **Blocked on the obvious:** the breast machine's own existence claim rests on **ONE observational line** (both CN instruments retracted, MH-133). Porting an unvalidated machine to a second tissue multiplies the claim, not the evidence. |

---

## Dead CHANNELS whose measured survivors were extracted first (archived 2026-07-17)

Each doc describes a channel that is **dead**. Each also held **live, measured results with no other home**.
Those were extracted to their proper homes *before* archiving — the docs go, the measurements stay.

| doc | the channel's death | survivors, and where they now live |
|-----|---------------------|-----------------------------------|
| `CN_INSTRUMENT.md` | ⛔ **BOTH CN instruments RETRACTED.** `pi_causal = γ·b` with **b an OBSERVATIONAL OLS slope** ⇒ it re-derived the correlation it was built to validate (γ is site-blind: HE +0.19931 vs decoy +0.19844, p=0.20; clean reduced form arm-clustered **p=0.115**). The within-arm two-way-FE replacement's control class was **71% real binders**; against a genuinely site-free control **τ=−0.0007, p=0.84**, failing the site-type efficacy ladder. **Refutations at n=216k–235k, not power failures.** | **`LEARNED_MODEL_METHODS{,_FORMAL}.md` §8a** — the **F>10 ∧ T1-clean admission gate** (exclusion FAILS for ~73% of well-instrumented edges; "F>10" alone overstates the causal-usable set ~4×; T1 is CONSERVATIVE so 27%-clean is a LOWER bound) and the **soft F-weight** (`_ADMIT_SOFT_C=10`: same admitted set, 4.4× less null leakage). ⚠ *These were never uniquely in this doc — the finding's home is registry **MH-89/90**, the formula was in the (already archived) fusion design. What was missing was the estimator SPEC; §8a is that.* |
| `LEARNED_MODEL_STATE_CHANNEL_PLAN.md` | ⛔ **MEASURED AND CANCELLED** (MH-102d). τ ≈ 0; info **0.6–0.7%**; `learned/channel_state.py` never built. Structural: precision ∝ a² and a≈0.11, so **even unlimited GTEx donors add ~1%**. | **`METHODS §1a`** — the **panel decision** (β is ρ=0.94 across *entirely different* atlases vs ρ=0.81 with **no C at all** ⇒ *which* control you use matters ~10× less than *whether* you control; and the tumour C **cannot** move to HBCA in principle, since `mal_prolif` is defined on the Cancer-Epithelial fraction and a normal-breast atlas has no malignant class) + the **matched-C policy** (CPE/target_cn/mal_prolif are tumour-only; NAT/GTEx must have NO purity term — a participant-keyed join silently hands 107/113 NAT samples their own patient's **tumour** CPE). **`METHODS §18a`** — the **β gauge convention** (raw-`r` for the MODEL, z-Y GAUGE-ONLY). ⚠ *Its §10 still cites the **VOID** 1.7%/≤8.8% Fisher figures; that void number is archived with it.* |
| `CPTAC_PROTEIN_CHANNEL_PLAN.md` | ⛔ **Centrepiece `βᵗ` FALSIFIED at n=101** (MH-103) — the prior positive evidence was a **mediator LEAK** (`a` fit on all samples, then OOF-ed). | **`LEARNED_MODEL_VALIDATION.md` §6.1–§6.5** — `locus_cn_cptac` (r=0.997, but fidelity ≠ power: F>10 for only **59/685** arms) · **δ-transportability** (same member dominates **84.5%** vs 43.7% chance; member-share r=0.991 — scope: δ's *abundance* input only) · **Bar-5 PASSES** (581 genes; and it **locates the limit**: the COHORT jump costs ~80%, the mRNA→PROTEIN jump ~0%) · **`a_g`** (complex subunits buffered, 0.117 vs 0.437, p=1.4e-29 — ⚠ with MH-106's refinement that this is a **large-obligate-complex** effect and the RPL/RPS proxy **over-states** the general claim) · the **Fisher bound** (protein carries ~4–6%, ≤7.6% at any `a_g` ⇒ **never a coupling lever**). Its §6 VOID numbers were corrected before archiving. |
