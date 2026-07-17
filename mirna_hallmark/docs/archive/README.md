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
