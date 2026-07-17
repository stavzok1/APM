# mirna_hallmark — Claude Code

## Before any task here

Read **`mirna_hallmark/AGENTS.md`** first. This subproject has its **own** analysis catalog,
discovery registry, and run ledger under `mirna_hallmark/docs/` — do **not** mix them with the
main `analysis/` docs, `analysis/DISCOVERY_REGISTRY.md`, or `analysis/docs/ANALYSIS_RUN_LEDGER.md`.

## What this subproject is

A focused study of **miRNA regulation of MSigDB Hallmark programs** in TCGA-BRCA. Reuses parent
APM data and loaders; writes all outputs under `mirna_hallmark/output/`.

## Essential commands

- **Python:** always `.venv/bin/python3` from **repo root** (not from this subdir)
- **Run a module:** `.venv/bin/python3 -m mirna_hallmark.<module>`
- **Full pipeline:** `.venv/bin/python3 -m mirna_hallmark.run_all` (add `--skip-mirna-cnv` once CNV cache exists)

## Module spine (run order)

```
config.py                        paths, strata, AGO gate params, cutoffs
hallmark_sets.py                 parse annotations/GSEA/h.all.v2026.1.Hs.txt
build_edges.py                   vectorized Hallmark-scoped miRTarBase edges
data_loaders.py                  clinical, miRNA arms, RNA, AGO/RISC, target CNV (cached)
stats.py                         bh_fdr, kruskal_across_strata, hypergeom, zscore_rows
ago_gate.py                      per-sample RISC capacity + saturating gate
mirna_locus_cnv.py               miRNA-universe CNV by stratum
hallmark_interaction.py          AGO-gated pressure + anti-correlation (CORE)
run_all.py                       orchestrator → output/run_manifest.json
```

## DOCUMENTATION PROTOCOL — MANDATORY, do not skip

**Documentation is a BLOCKING step at the end of every arc, not an afterthought.** An "arc" = a set of related
analyses/builds that reaches a conclusion, verdict, or deliverable.

**At the end of ANY arc you MUST, before treating it as done:**
1. **ASK the user**: *"This concludes the &lt;arc&gt; — shall I document it, and confirm the scope?"* — unless they
   already told you to document. **Never silently finish an arc.** (This was added because arcs kept ending undocumented.)
2. **Update EVERY applicable canonical home** (usually several, not one):
   - `docs/ANALYSIS_RUN_LEDGER.md` — an `MH-##` row for any run/build that changed a finding or output.
   - `docs/DISCOVERY_REGISTRY.md` — a claim row (with an S/A/P/R strength tag) for any new finding/verdict.
   - Learned-arc docs when the learned model changed: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` (narrative),
     **`LEARNED_MODEL_METHODS.md` + `LEARNED_MODEL_METHODS_FORMAL.md` (TWIN docs — same spec, plain-text vs LaTeX;
     ALWAYS edit BOTH)**, `LEARNED_MODEL_RATIONALE.md` (why), `LEARNED_MODEL_VALIDATION.md` (numbers),
     `LEARNED_MODEL_ESTIMATOR_MAP.md`.
   - `docs/DATA_SOURCES.md` — any new dataset/source (+ where in code it's used).
   - Memory — durable cross-session facts.
3. **Route by conclusion type:** a RESULT → VALIDATION + narrative; a METHOD → the two METHODS twins; a WHY →
   RATIONALE. A *modeling* conclusion does **NOT** go in `BIOLOGICAL_SYNTHESIS.md` (biology-findings only).
4. **Retract/annotate** any doc claim this arc superseded — never leave stale/contradicted text.
5. Then verify: grep the trackers for the new terms; confirm each deliverable appears in its home(s).

**Creating a NEW doc — it MUST open with this header block** (defines the doc's contract):
```
> **Goal:** <one line — what this doc is for>
> **What belongs here:** <which findings/cases go into this doc — and which do NOT>
> **Update trigger:** <when this doc must be updated>
> **Sync-partner:** <twin doc(s) that must be edited together, or "none">
```
Then register it in `docs/INDEX.md` with its role tag (`active`/`derived`/…).

Follow the ledger skill for the mechanics: `.claude/skills/apm-analysis-ledger/SKILL.md`.

## Working axioms — MANDATORY (not optional polish)

These are completion gates, on the same footing as the documentation protocol:

1. **VERIFY BEFORE ASSERT.** Never state a fact about build-status, identity/equivalence, safety, or a bottleneck
   ("X is unbuilt/done", "bit-identical", "no-op", "the bottleneck is Y") until you have RUN the check: grep the
   **ledger** for build-status, run a **diff/control** for equivalence, `cProfile` for bottlenecks. Tag claims
   VERIFIED vs INFERRED; if you can't cite the check, don't state it as fact. (memory `verify-before-assert`)

   **1a. THE MEASURED-ONLY GATE (the teeth — added 2026-07-12 after this failed FOUR times in one session).**
   A **number** or a **direction** (bigger/smaller/stronger/understated/overstated/helps/hurts) may enter a
   PERSISTED artifact — registry, ledger, board, memory, docstring, or a summary stated to the user as fact —
   **only if you personally produced it by running something.** Otherwise it stays in chat, tagged `[INFERRED]`.
   - ⚠ **Arithmetic on someone else's number is NOT a measurement.** Rescaling, extrapolating, or
     unit-converting a reported statistic to manufacture a new one is the single most common way a fabricated
     number enters a doc. **Never transform a reported statistic — RE-DERIVE it from the data.** (The session's
     worst case: multiplying a published Shapley width `0.41 × 1.73 → ±0.71` — invalid, it exceeds the [0,1]
     bound, AND the true correction went the *opposite* direction. Re-deriving gave 0.243.)
   - **Every quantitative claim in a registry/ledger cell must be traceable to that row's `evidence (output)`
     column.** If you cannot name the file or run that produced the number, it does not go in.
   - **PREDICT BEFORE YOU TEST.** Write the prediction down first. A wrong prediction is a signal that reasoning
     is unreliable *in this domain* — from then on, every downstream claim in that domain needs a check, not an
     argument. (Every one of this session's four errors was a confident prediction that the test refuted.)
2. **TRACE THE DOWNSTREAM RIPPLE.** When a change touches a shared output/function/schema, `grep -rln` every
   consumer, re-run the stale persisted ones, check for shared-mutation, and hunt residual old-pattern uses —
   *before* declaring done. (memory `downstream-ripple-canonical`)

   **2a. THE CONFOUNDER BLOCK IS NOT A "MORE IS BETTER" KNOB — KNOW WHICH `C` EACH READOUT USES, AND REPORT
   RETENTION INSTEAD OF SILENTLY CONDITIONING.** (Added 2026-07-12 after MH-107; **corrected same day — the first
   version wrongly said "always use the same C everywhere", which would IMPOSE an over-control.**)
   - **The design is deliberate, not sloppy** (`LEARNED_MODEL_METHODS.md §1`, `DESIGN_RESPONSE` Decision B):
     *"tumour attribution = fractions; **gate/coupling = core only unless deconv requested**"*, and
     *"**Cancer-Epithelial is deliberately EXCLUDED from the `deconv` block** — conditioning on the compartment the
     target is expressed in **over-controls the signal**"*. **A miRNA acting through a cell-STATE shift PRODUCES a
     composition change — the composition is partly the MECHANISM, not a nuisance.** Same trap as MH-100
     (proliferation) and MH-101 (host gene): a global control **over-controls the majority**.
   - **THE RULE: never silently condition a plausible-mechanism covariate away. Run BOTH blocks and report the
     RETENTION** (`retention = β_deconv/β_core`; `card.py` / `learned/readouts.py`:
     **cell_intrinsic ≥0.7 · partial · composition_explained <0.4**). **FLAG, don't delete** — the reader decides.
     Measured on the HE universe: **~20% of family couplings are `composition_explained`.**
   - **The MH-107 error, stated correctly:** `cptac_validation` / `ood_protein` / `dossier.tier3_protein` had **NO
     composition variant at all**, so they **never computed retention** — and a compartment-driven result was
     presented as cell-intrinsic. (EMT proteins ↔ CAF fraction **r=+0.509**; ZEB1 protein **+0.768**; miR-200, an
     EPITHELIAL miRNA, ↔ CAF **−0.35** ⇒ the flagship "miR-200 represses EMT protein" is compartment arithmetic.
     Cost: gene-level FDR hits **66→3**; **27%** of the "268 triple-validated" edges sign-flip.)
   - **Adjudicate confound-vs-mechanism, don't assume:** test whether the effect survives the adjustment in the
     **well-powered** cohort (miR-200→ZEB1 mRNA in TCGA: ρ=−0.209, **p=8.7e-12** after composition ⇒ CONFOUND for
     that edge — the edge is real, only its bulk-protein magnitude was inflated). Per-gene verdicts: MH-100's
     OOF 2×2. (memory `cptac-protein-composition-confound`)
3. **RUNTIME IS A DELIVERABLE.** Proactively `cProfile` + optimize any batch/pipeline run — don't wait for a
   timeout or a user request. The win is usually redundant data loading, not the numerics; cache DERIVED objects,
   not just raw reads. Verify optimizations are output-identical (a control run). (memory `profile-before-batch-loops`)

4. **A NULL MUST BE ABLE TO CARRY THE ARTIFACT — AND "IT REPLICATES" IS NOT EVIDENCE OF BIOLOGY.**
   (Added 2026-07-12 after MH-114, where BOTH of these gave a confidently wrong verdict before the right test
   overturned it. They are the default way a confounded result survives scrutiny.)
   - **Before trusting a null or a replication, state what the ARTIFACT would do to it.** If the artifact and
     the biology predict the SAME outcome, **the test is not an arbiter** — no matter how rigorous it looks.
   - **A shuffled/permutation null compared on the MEAN is blind to a SIGN-SYMMETRIC artifact.** A random
     miRNA–target pair is equally likely to be same-compartment (artifact **+**) or opposite (artifact **−**)
     ⇒ the artifact **cancels in the mean** (measured: shuffled mean ρ = +0.0001) and lives in the **spread**.
     ⇒ **STRATIFY THE NULL BY THE CONFOUND'S OWN AXIS.** Restratified by compartment orientation, **shuffled
     non-edges reproduced the real edges' gradient EXACTLY (−0.1290 vs −0.1272)** — pairs that *cannot* repress
     produced the same bulk anti-correlation ⇒ the signal WAS arithmetic.
   - **Cross-cohort reproducibility cannot arbitrate confounding when the cohorts SHARE the confound** — the
     artifact replicates. (Composition adjustment *lowered* prospective↔TCGA-105 concordance +0.263 → +0.123,
     because the shared CAF-rich stroma was the reproducible part.)
   - **The identification strategy is an ASYMMETRY the confound is blind to.** Here: compartment arithmetic
     doesn't know about 3'UTR sites; repression doesn't know about cell type. Find that asymmetry.
     (memory `null-design-and-shared-confounds`)

## Key docs (this folder only)

| Doc | Purpose |
|-----|---------|
| `docs/INDEX.md` | **Doc map** — every doc grouped by role (canonical / living / reports / parked / historical) with status + link graph. Find any doc here. |
| `README.md` | Quickstart + output map |
| `docs/MODELING_FRAMEWORK.md` | Conceptual + implemented framework (resolution hierarchy, edge universe, pressure, architecture, coupling/confounding, protein, validation) — "how is this modeled and why" |
| `docs/BIOLOGICAL_SYNTHESIS.md` | Biology-first findings surface (start here for "what did we learn") |
| `docs/DCIS_EV_SYNTHESIS.md` | DCIS→invasive + EV arc synthesis (MH-48..56): two-compartment/two-phase model, stroma reframe, caveats, what-to-pursue |
| `docs/PRIMARY_TUMOR_ROADMAP.md` | Phased roadmap to take the rigor/compartment lens to primary tumours + the core framework (composition-retest, spatial, orphan retrofit) |
| `docs/MIRNA_CNV_DOSAGE_REPORT.md` | miRNA CNV dosage report |
| `docs/WHATS_NEXT.md` | Forward-looking extensions and caveats |
| `docs/REPORT.md` | Detailed results |
| `docs/METHODS.md` | Formal methods: pressure, gate, enrichment, stats |
| `docs/DATA_SOURCES.md` | Every dataset read (core spine, priors, validation, DCIS/EV corpus) + preprocessing provenance |
| `docs/ANALYSES_CATALOG.md` | One row per analysis component (inputs/outputs/CLI) |
| `docs/DISCOVERY_REGISTRY.md` | Claims + strength tags |
| `docs/ANALYSIS_RUN_LEDGER.md` | Per-component last-run timestamp + status |

## Environment

Running on Linux at repo root `/sci/labs/michall/stavzok/APM`. No WSL bridge needed here.
