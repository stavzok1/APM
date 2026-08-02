# mirna_hallmark — Claude Code

> **⭐ THIS FILE IS THE CANONICAL onboarding for `mirna_hallmark/`** — orientation, the reuse contract, the
> module spine, the documentation protocol, the working axioms, and the doc map. Claude Code auto-loads it; other
> tools (Cursor) are sent here by the thin `AGENTS.md`. There is no second onboarding doc to read.
> **First, get oriented:** `docs/ARCHITECTURE.md` (the axes×models×analyses×results map) → `docs/STATE_OF_PLAY.md`
> (per-axis verdicts).

## Before any task here

This subproject is **SEPARATE** — it has its **own** analysis catalog, discovery registry, and run ledger under
`mirna_hallmark/docs/`. Do **not** mix them with the main `analysis/` docs, `analysis/DISCOVERY_REGISTRY.md`, or
`analysis/docs/ANALYSIS_RUN_LEDGER.md`.

## What this subproject is

A focused study of **miRNA regulation of MSigDB Hallmark programs** in TCGA-BRCA. Reuses parent APM data and
loaders; writes all outputs under `mirna_hallmark/output/`. Central questions: which Hallmark programs are most
targeted by high-evidence miRNAs; does evidence-weighted miRNA pressure anti-correlate with target expression
(AGO/RISC-gated); how do the miRNA loci and target genes vary by copy number across strata.

## Reuse, don't duplicate (import from the parent, don't re-implement)

| Need | Reused from |
|------|-------------|
| miRTarBase normalize | `pipeline.genes.mirtarbase.load_mirtarbase` (+ constants) — **NOT** the slow `get_mirtarbase_targets`; use `build_edges.compute_interaction_summary_fast` (vectorized) |
| miRNA pressure | `mirna_hallmark.pressure_build` / `pressure_engine` (spine `softmax_z_logrpm + evidence_mass`) — ⚠ the pressure heuristic is **§6b-RETIRED** (MH-115); the canonical estimator is the learned Gibbs (`learned/`) |
| miRNA-universe CNV | `analysis.cohort_landscapes.cnv.dosage_landscape_cnv` |
| clinical / RNA / partial Spearman | `analysis.utils.common.loaders` |
| ⭐ **ANY per-gene question** ("where does X help / which genes / what predicts it") | **`learned/gene_axes.py`** — build the axes, `scan` them under FDR, `contrast` the extremes, `sign_analysis` the gain. **Consult it BEFORE hand-rolling a stratification; it encodes two traps that have already cost errors** (see axiom 8) |

*(Confounder-block `C` policy, module spine, and design decisions are NOT restated here — they have homes:
axiom 2a below, `docs/ARCHITECTURE.md` + `docs/ANALYSES_CATALOG.md`, and `docs/MODELING_FRAMEWORK.md`.)*

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

   **3a. PARALLELIZING A BATCH ON THIS NFS BOX — READ THIS BEFORE YOU WRITE A `Pool` (learned the hard way,
   2026-07-18; the genome-wide card took a full session of dead ends).** In order:
   1. **Memoize EVERY static-file read at module level FIRST.** A `pd.read_csv` that is negligible for one item
      is *fatal* at N_items × N_workers on NFS (thrash). Audit the whole per-item call path — `cProfile` on one
      item HIDES this (each read looks small; only N×workers makes it dominate). This is the load-bearing fix.
   2. **Single-thread BLAS BEFORE `import numpy`** (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`) —
      else a fork after a BLAS-heavy step inherits a locked OpenMP threadpool → **futex DEADLOCK** at ~2% CPU,
      and 7 workers × unbounded threads oversubscribe 8 cores.
   3. **Prefer INDEPENDENT SUBPROCESS SHARDS over in-process `multiprocessing.Pool`.** The Pool deadlocked/hung
      **three different ways** here (fork-after-threadpool; spawn NFS thundering-herd → dead worker → `map`
      hang; COW-growth stall). Independent shards + concat have nothing shared to deadlock. Isolate a native
      crasher (`double free` on a degenerate design, e.g. GAPDH) by running the shard's items one-per-subprocess.
   4. **A small-scale test gives FALSE confidence** — NFS contention is *nonlinear* (40 genes passed, 1,436
      thrashed on identical code). **Verify the MECHANISM at FULL scale:** workers must be COMPUTE-BOUND
      (`ps STAT R`, ~90% CPU) — not wall-time on a subset. Diagnose stalls via `/proc/<pid>/wchan`
      (`futex_wait`=lock/threadpool; `pipe_read`=blocked on a worker), never guess. (memory `batch-nfs-per-unit-reads`)

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

5. **⭐ NEVER HEADLINE A RATIO OR A THRESHOLD COUNT WITHOUT LOOKING AT WHERE THE MASS SITS.**
   *(Added 2026-07-17 after this fired **FOUR times in one day**, each time producing a confident number that a
   trivial perturbation moved. It is now the single most common way a non-finding enters this registry.)*
   - **The tell:** the quantity underneath is *fine* — it is the **summary statistic** that is unstable.
     MH-144 is the cleanest case: two implementations whose per-gene ceilings agree at **corr 0.9956** report
     **17.6% vs 25.5%** "not measurable", because `ceiling ≤ 0` is a threshold and **39.8% of genes sit within
     ±0.01 of it**. A shift of **0.06% of the scale** reclassified 8% of the universe.
   - **The four instances.** ① **MH-138** — a per-gene median on a **BIMODAL** distribution (21% near 0, 43% at
     exactly 1.0) moved **55.5% → 70.9%** under a change licensed as *noise*; the 50th percentile lives in the
     empty valley between the modes. ② **MH-144** — the threshold count above. ③ **MH-119** — `share = β_f/Σβ`
     read **999%** because a gene's βs cancelled to Σβ ≈ 0. ④ **MH-146** — a ratio of two SDs whose
     **denominator is 0 for 48% of rows**; the tool's printed headline was **112.86×**, the honest gated value
     **~1.4×**.
   - **THE RULE.** Before a ratio or a threshold count enters a doc: **(a)** plot/quantify where the mass sits
     relative to the boundary or the denominator; **(b)** if there is a pile-up, **gate the denominator**
     (`identity > 0.05`, `|z| > 2`, `net_pressure ≥ 0.5` — the `readouts.add_reliability` pattern) and report
     the **gated** statistic; **(c)** if a threshold is unavoidable, **sweep it** and report a value from a
     region where the two arms agree (MH-144: `≤ 0.02` reproduces at ~52% where `≤ 0` does not); **(d)** state
     the robust alternative next to the fragile one, or drop the fragile one.
   - **A statistic evaluated where its denominator vanishes is not a finding; it is a coin-flip with decimals.**
     ⇒ registry rows MH-138, MH-144, MH-146; the gates in `learned/readouts.py::add_reliability`.

6. **⭐ BEFORE FIXING A DISCREPANCY, VERIFY THE TWO THINGS ARE THE SAME OBJECT — MOST "BUGS" HERE ARE NAMING COLLISIONS.**
   *(Added 2026-07-17 after this fired **THREE times in one day**. In all three a "fix" would have changed
   CORRECT code, and in two of them the fix would have silently destroyed a load-bearing design feature.)*
   - **The tell:** two artifacts disagree, someone labels it a bug, and the label propagates faster than the
     evidence. **Read the evidence column, not the title.**
   - **The three instances.** ① **MH-138** — "§5 (bagged NNLS) vs §6b (dense posterior)" looked like a
     contradiction about *canonical attribution*. It was **ONE WORD DOING TWO JOBS**: MAGNITUDE (the Gibbs
     mean) vs IDENTITY (Shapley fed NNLS weights). Both docs were right. ② **MH-142** — E6's *"SELF-PARTIALLING
     BUG **(live)**"*: production scores **orphans**, disjoint from `he_agg` by construction (**0 of 6,744**
     dossier pairs are HE edges). Its own evidence column said *"class-matched estimator comparison"* — the
     effect was an artifact of **MH-124's own experiment**. Removing `he_agg` would have converted the
     discovery lane from *"couples BEYOND curation"* to *"couples at all"*. ③ **MH-145** — the deferred
     *"`coupling_grid` field bug"*: `grid_full`/`grid_oof` call `LR.fit_gene` = **the adaptive lasso §6b
     RETIRED**, while `coupling_grid` reads the Gibbs. **Two estimators, one label.** Nothing to reconcile;
     MH-116/117 needed re-**labelling**, not re-running.
   - **THE RULE.** Given "A and B disagree, it's a bug": **(a)** name the ESTIMAND each answers — different
     estimands are not a bug (axiom 5's cousin: *ask what question the number answers*); **(b)** name the
     ESTIMATOR each runs — grep the call, do not trust the docstring or the arm's label; **(c)** name the
     UNIVERSE/UNIT each scores — MH-142 died on `0 of 6,744`, a one-line check; **(d)** read the claim's
     **evidence column**, which records what was actually measured, over its **status tag**, which records what
     someone believed. **A `(live)` tag is not evidence.** *(Same failure as MH-38's "108 triple-validated
     orphans", which propagated for a session on its `S/R` tag while its input had collapsed 594→23.)*
   - **If the check is cheap, run it before writing the fix.** All three collapsed to a single grep or a
     one-line set intersection. ⇒ registry rows MH-138, MH-142, MH-145; guard `eval/_e6_production_check.py`.

7. **⭐ CLASSIFY THE QUANTITY, NOT ITS CONTAINER — a unit/estimator belongs to the CALL that produced it.**
   *(Added 2026-08-01 after this fired **THREE times in one day**, each time producing a confident, wrong
   statement about what a number IS.)*
   - **The tell:** you are about to say what a number means, and your evidence is the **doc** that
     describes it, the **module** it lives in, the **artifact** it sits on, or an **aggregate statistic**
     computed over it. All four are containers. None carries the unit reliably.
   - **The three instances.** ① **MH-191** — `beta`'s rung inferred from the docs (*"estimated per family,
     broadcast to arms"*) and from an aggregate (*"one distinct β in 90.6% of cells"*). Both wrong:
     `readouts.run(level="arm")` fits it **PER ARM**, and the 90.6% was simply the ~90% of cells that are
     SINGLETONS. `readouts.py:185` settles it in ten seconds. ② **MH-192** — "these labels are wrong"
     inferred from **my own verifier**, whose `nunique(dropna=False)` counted NaN as a value, so every
     partially-covered column looked like a violation: **4 false positives, and the check was never
     checked.** ③ **MH-194** — `static_owner_family` classified by the **MODULE** (`realization.py`,
     "mostly per-arm") instead of the **COLUMN**; it is an attribution claim at the wrong rung, and the
     fix moves the owner for **11.7%** of the genes where the rung can matter.
   - **THE RULE.** Before asserting a quantity's unit, estimator or rung: **(a)** grep the CALL that
     produces it — not the docstring, not the module name, not the column name (axiom 6's cousin);
     **(b)** ask what the UNIVERSE'S SHAPE does to any statistic you use as evidence — a skewed universe
     (90% singletons, 30% NaN) is happily consistent with the opposite fact (axiom 5's cousin);
     **(c)** when a TOOL flags something, validate the tool on a case whose answer you already know.
   - **A property is real; the container is just the wrong place to read it from.**
     ⇒ registry rows MH-191, MH-192, MH-194; guard `learned/card_rungs.py` (`--check`).

8. **A PER-GENE QUESTION IS A MULTI-AXIS QUESTION — AND AN AGGREGATE WIN IS NOT PER-GENE ACTIONABLE UNTIL
   YOU SHOW THE HARM IS GATEABLE.** (MH-201, 2026-08-02.) ⇒ **use `learned/gene_axes.py`; do not hand-roll
   a stratification.** A single number for "does X help?" is close to uninterpretable: in MH-201 the headline
   margin (β −0.0186 vs abundance −0.0022) only became a finding once asked along many gene axes, and **three
   of the four things that came out were on axes I had not thought to build.**
   - **BUILD ALL FOUR AXIS FAMILIES, not just the gene's own.** ⭐ **The REGULATOR-ENSEMBLE axes are usually
     the strongest and are the ones most often forgotten** — `reg_dose_hhi` was the largest effect in the
     whole test (q=2e-05). The mechanism generalises: **a weighted estimator cannot beat an unweighted sum
     when one member dominates the ensemble's abundance, because the sum already IS that member.**
   - **DISPERSION, NOT LEVEL.** A gene's / regulator's **dynamic range** predicts (`self_sd`, `buffa_sd`
     q≈0.001); its **MEAN expression is nothing** (q=0.45). A flat feature cannot correlate however high it
     sits. Always carry a dispersion term or you will conclude expression doesn't matter.
   - ⛔ **THE DEGENERACY TRAP — always ask where your candidate is MATHEMATICALLY INERT and SPLIT on it.**
     With ONE unit, `Σw·Z = w·Z` and the unweighted reference is `Z`; Spearman is scale-invariant, so they
     are **IDENTICAL BY CONSTRUCTION** (verified max|Δ| = 0.0). That was **43%** of MH-201's universe, and
     pooling it produced an incoherent *"median Δ = +0.0000 with Wilcoxon p = 1.3e-27"* and made an
     ABUNDANCE effect read as a β effect. `gene_axes.mask_degenerate()`.
   - ⛔ **THE MOVING-SUPPORT TRAP (axiom 5's cousin).** Raw HHI is bounded below by `1/k`, so `corr(k, HHI)`
     is strongly negative **whatever the biology does** — I reported **−0.667** as evidence that big designs
     are diffuse; floor-corrected it is **−0.075**. `gene_axes.hhi()` normalises by default.
   - ⭐ **ASK WHETHER THE GAIN IS SIGN CORRECTION, NOT SHARPENING.** They mean different things and a mean
     hides it. MH-201: the unweighted reference had the right sign in **51.1% — chance** — vs β's 58.5%;
     net **+49** rescued; and the correction was **ASYMMETRIC** — it moved wrong answers (p=1.5e-16) and
     moved right ones by **exactly 0.0000**. `gene_axes.sign_analysis()`.
   - ⭐⭐ **THEN ACCOUNT FOR THE HARM SYMMETRICALLY AND TEST WHETHER IT IS GATEABLE — AND TEST THE GATE ON
     THE AXIS FAMILY THAT PREDICTS THE *SIGN*, NOT THE ONE THAT PREDICTS THE *MAGNITUDE*. THEY ARE
     DIFFERENT FAMILIES, AND CONFLATING THEM COST ME A WRONG CONCLUSION.** MH-201: β helped 55.7% of genes
     but **HURT 34.7%** (gains beat losses only **1.55×**). I gated on the ENSEMBLE axis (`reg_dose_hhi`),
     found the hurt-rate flat at 34–36%, and concluded *"not gateable, β is an aggregate bet only"*.
     **WRONG.** The **COMPOSITION** axes are the ones that predict the SIGN — the flagship failure FN1 sits
     at the **89th percentile** of `comp_tcga_mrna_driver_share` and was fully predictable, and that axis
     separates the tails at **p=0.00999** where the ensemble axis gives p=0.10. **Both gates together:
     hurt 34.7% → 27.6%, sign-correct 58.5% → 69.3%, mean margin 2.3×, while retaining EVERY showcase gene.**
     ⇒ **ensemble axes = how much it gains where it gains; confounder/retention axes = whether it gains at
     all.** Ask which question your candidate gate answers before concluding anything from it. Report a
     showcase gene (MYC: +0.099 → −0.340) only alongside its counterexample (FN1: −0.216 → +0.121) **and
     the gate that separates them.**
   ⇒ registry row MH-201; module `learned/gene_axes.py` (self-checked against every number above).

## Key docs (this folder only)

**Read in this order. Do not skip 1 — it exists so you do not have to reconstruct the state from 48 docs.**

| # | Doc | Purpose |
|---|-----|---------|
| **0** | **`docs/ARCHITECTURE.md`** | **⭐ HOW IT ALL FITS** — the generated axes×models×analyses×results map (12 axes + a `shared` bucket). Start here for orientation; it links into the docs below. Interactive view: `docs/derived/architecture.html`. Regenerate: `analyses/ops/gen_architecture.py`. |
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

## Guardrail agents (`.claude/agents/`) — use them

- **`rigor-auditor`** (opus) — run it on a finding BEFORE it enters `DISCOVERY_REGISTRY`. It wraps
  `apm-rigor-protocol` + the six axioms and returns PASS / CONCERNS / FAIL. The gate against a confounded or
  unmeasured claim getting a row.
- **`consistency-sweeper`** (sonnet, read-only) — run it before a wrap-up. It audits axis-tag coverage,
  `ARCHITECTURE.md` freshness, catalog↔modules, ledger↔outputs, `/tmp` reads, and import orphans, and returns a
  drift report. (AST imports, not regex; scope-aware catalog — see its prompt.)

The architecture map, its join keys (`docs/derived/{axis_assignment,module_axis}.tsv`), and the map generator
(`analyses/ops/gen_architecture.py`) are maintained BY these + the regenerate-not-hand-edit discipline.

## Environment

Running on Linux at repo root `/sci/labs/michall/stavzok/APM`. No WSL bridge needed here.
