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
2. **⭐ THE ONE-HOME RULE — write each fact to EXACTLY ONE doc. Everything else LINKS.**
   *(Rewritten 2026-07-17. The previous rule said "update EVERY applicable canonical home (usually several,
   not one)" — a **fan-out rule with no consolidation rule**. It is the measured cause of this subproject's
   doc rot: a finding lands in N docs, and when it is later refuted the retraction lands in 1 of N.
   **This is not a prediction — it is the history.** MH-38's "108 triple-validated orphans" was retracted
   inside MH-114's row and left live with an S/R tag in its own row and the ledger for 5 days, while the
   talk still cited it. The same protein Fisher-information number exists in three docs at three values,
   one of which is VOID. **Do not restate a number in a second doc. Link to its home.**)*

   | the fact you produced | its ONE home |
   |---|---|
   | a finding / verdict / retraction | **`docs/DISCOVERY_REGISTRY.md`** — one `MH-##` row, S/A/P/R tag, `evidence (output)` filled |
   | a run (when, status, runtime) | **`docs/ANALYSIS_RUN_LEDGER.md`** — one row |
   | the CURRENT state of an axis | **`docs/STATE_OF_PLAY.md`** — the verdict only, ≤1 block, no derivation |
   | a formula / estimator spec | **`FORMULAS.md`**, or the **`LEARNED_MODEL_METHODS{,_FORMAL}.md` TWINS** (same spec, plain vs LaTeX — **ALWAYS edit BOTH**) |
   | *why* a design is the way it is | **`LEARNED_MODEL_RATIONALE.md`** |
   | a dataset / source | **`docs/DATA_SOURCES.md`** |
   | a component (inputs/outputs/CLI) | **`docs/ANALYSES_CATALOG.md`** |
   | what to do next | **`docs/PROGRAM_FORWARD_BOARD.md`** |
   | a biology finding | **`BIOLOGICAL_SYNTHESIS.md`** — biology only; a *modeling* conclusion does **NOT** go here |
   | a durable cross-session fact | memory |

3. **A NEW DOC IS THE EXCEPTION, NOT THE DEFAULT — and a per-finding `MH###_*.md` is now BANNED.**
   That pattern produced six docs / 184 KB (MH-124…130), was abandoned at MH-130, and is precisely where the
   retractions rotted: MH-127 still reads as current while MH-137 has superseded its central number four times.
   **The registry row IS the record of a finding.** Write a new doc only for a *durable reference* that has a
   home in the table above and no owner yet — and if you do, it MUST carry the header block below and be
   registered in `docs/INDEX.md`. **If you cannot name the doc's non-overlapping purpose in one line, do not
   create it.**

4. **⭐ A RETRACTION GOES IN THE RETRACTED CLAIM'S OWN HOME — not in the retractor's.**
   Recording "MH-38 is dead" inside MH-114's row is **not** a retraction; it is a footnote nobody reading MH-38
   will ever see. **Edit the original row**: prepend `⛔ RETRACTED (MH-###, YYYY-MM-DD)` or
   `⚠ SUPERSEDED BY MH-###`, downgrade its strength tag, and *then* cross-link. Same for the ledger row and for
   `STATE_OF_PLAY.md`. **Grep for the claim's headline number across `docs/` and fix every hit, or you have not
   retracted it.**

5. **Then verify, don't assume:** grep the trackers for the new terms; confirm the deliverable appears in its
   ONE home and that no *second* doc restates it; confirm every superseded claim carries its banner.

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
   - **The design is deliberate, not sloppy** (`LEARNED_MODEL_METHODS.md §1`; the `DESIGN_RESPONSE` Decision B
     that this used to cite is now `docs/archive/`, historical):
     *"tumour attribution = fractions; **gate/coupling = core only unless deconv requested**"*.
     **A miRNA acting through a cell-STATE shift PRODUCES a composition change — the composition is partly the
     MECHANISM, not a nuisance.** Same trap as MH-100 (proliferation) and MH-101 (host gene): a global control
     **over-controls the majority**. ⇒ **the RULE below stands.**
   - ⛔ **BUT ITS STATED RATIONALE IS FALSE — do not repeat it (MH-111/114, 2026-07-12).** `METHODS §1` justifies
     the block by *"Cancer-Epithelial is deliberately EXCLUDED — conditioning on the compartment the target is
     expressed in over-controls the signal"*. **That is void:** the 9 Wu-major fractions sum to **exactly
     1.000000** and **R²(Cancer-Epithelial ~ the 8 "held-out") = 1.000000 in BOTH cohorts** (TCGA n=1059, CPTAC
     n=133) ⇒ conditioning on the 8 **algebraically determines** tumour content. **The hold-out is a simplex
     illusion.** The practice is nevertheless safe — but for a *different, measured* reason: purity contributes
     only **0.9%** of the retention loss (CPE in core C already absorbs it) while the **stromal MIX contributes
     33.1%**, and a block with tumour content removed reproduces the result. **Fixing `METHODS §1`'s rationale is
     open correctness debt** (`STATE_OF_PLAY.md` → CPTAC axis).
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

**Read in this order. Do not skip 1 — it exists so you do not have to reconstruct the state from 48 docs.**

| # | Doc | Purpose |
|---|-----|---------|
| **1** | **`docs/STATE_OF_PLAY.md`** | **⭐ WHERE WE ARE NOW.** Per axis (model · edge existence · attribution · decoy · CPTAC · progression · Buffa): what STANDS, what is DEAD, what is OPEN — with MH ids and dates. Carries the **doc-traps table**: which docs read as current but are not. |
| **2** | `docs/DISCOVERY_REGISTRY.md` | **THE SOURCE OF TRUTH.** One `MH-##` row per finding + strength tag. **Where this and any other doc disagree, THIS WINS.** ⚠ It runs AHEAD of the docs — MH-131…137 exist only here. |
| **3** | `docs/ANALYSIS_RUN_LEDGER.md` | One row per run: when, status, runtime. |
| 4 | `docs/INDEX.md` | The doc map — role + currency per doc. Find anything here. |
| 5 | `README.md` | Quickstart + output map |
| 6 | `docs/MODELING_FRAMEWORK.md` | Conceptual + implemented framework — "how is this modeled and why" |
| 7 | `docs/FORMULAS.md` · `docs/METHODS.md` | The exact formula for every quantity; the methods narrative |
| 8 | `docs/BIOLOGICAL_SYNTHESIS.md` | Biology-first findings surface (biology only — **not** modeling conclusions) |
| 9 | `docs/DATA_SOURCES.md` | Every dataset read + preprocessing provenance |
| 10 | `docs/ANALYSES_CATALOG.md` | One row per analysis component (inputs/outputs/CLI) |
| 11 | `docs/PROGRAM_FORWARD_BOARD.md` | What to do next |

Reports (`REPORT.md`, `LANDSCAPE_REPORT.md`, `MIRNA_CNV_DOSAGE_REPORT.md`, `DCIS_EV_SYNTHESIS.md`) and the
per-axis docs are in `docs/INDEX.md`. **Consumed / superseded docs live in `docs/archive/` — nothing there
describes current state.**

## Environment

Running on Linux at repo root `/sci/labs/michall/stavzok/APM`. No WSL bridge needed here.
