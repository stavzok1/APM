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
