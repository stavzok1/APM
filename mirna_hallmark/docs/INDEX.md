# Docs index — `mirna_hallmark`

> **Goal:** find the right doc fast, and know whether to trust it.
> **What belongs here:** one row per doc — role, currency, one line on what it is. **NOT** verdicts
> (→ [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md)), **NOT** findings (→ [`DISCOVERY_REGISTRY.md`](DISCOVERY_REGISTRY.md)).
> **Update trigger:** a doc is added, archived, or changes role.
> **Sync-partner:** [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) (its "Doc traps" table must agree with the ⚠ rows here).

**Index rebuilt 2026-07-17** (registry current through **MH-137**).

---

## 🚩 Start here

| # | doc | why |
|---|-----|-----|
| **1** | **[STATE_OF_PLAY.md](STATE_OF_PLAY.md)** | **The one doc that answers "where are we".** Per axis: what stands, what's dead, what's open — with MH ids and dates. Read this before anything else. |
| **2** | [DISCOVERY_REGISTRY.md](DISCOVERY_REGISTRY.md) | **The source of truth.** Append-only, current, one row per finding with a strength tag. Where this and any other doc disagree, **this wins**. |
| **3** | [ANALYSIS_RUN_LEDGER.md](ANALYSIS_RUN_LEDGER.md) | One row per run: when, status, runtime, what changed. |

> ⚠ **The registry runs AHEAD of the docs.** MH-131…137 — the most consequential fortnight in the
> project (both CN instruments retracted; the decoy control rebuilt; the first real positive control)
> — exist **only as registry rows**. There is no `MH13x_*.md`. A doc's silence is not evidence.

---

## 1. Canonical references — the load-bearing spine

| doc | currency | what it is |
|-----|----------|------------|
| [MODELING_FRAMEWORK.md](MODELING_FRAMEWORK.md) | active | **"What are we modeling and why" + ⭐ THE CANONICAL MODEL STATEMENT** (folded from SYNTHESIS §6b, 2026-07-17). ⚠ §9.3/§9.4's MH-38 "108 triple-validated orphans" is **RETRACTED** in place. |
| [FORMULAS.md](FORMULAS.md) | active | Code-faithful exact formula for every quantity. The *how*. |
| [METHODS.md](METHODS.md) | active | Methods narrative — data → pressure → gate → coupling. |
| [DATA_SOURCES.md](DATA_SOURCES.md) | active | Every dataset read, and its preprocessing basis. |
| [ANALYSES_CATALOG.md](ANALYSES_CATALOG.md) | living | **What code exists + how to run it.** Opens with an orientation section (the code tree, the dotted-path run rule, the 35 top-level modules classified spine/baseline/pipeline), then one row per component. ✅ Absorbed `MODULE_MAP.md` 2026-07-17. |
| [DECOY_PREMISE_REGISTER.md](DECOY_PREMISE_REGISTER.md) | active | ⭐ **Is the decoy control sound?** One table of every ASSUMPTION the site-free negative control rests on (fake validity · matching · scoring · what the result may be used for), each with a verdict (HOLDS / INERT / BOUND / OPEN / CLOSED-BY-DATA) and a pointer to the `MH-##` row that tested it. **A map, not a source — it deliberately carries no numbers**; construction lives in the `eval/decoy_bench.py` docstring, the headline in `STATE_OF_PLAY` Axis 4. |
| [EDGE_QUESTION_TAXONOMY.md](EDGE_QUESTION_TAXONOMY.md) | active | What a single miRNA→target edge can be asked. Backs `apm-edge-question`. |
| [GENE_QUESTION_TAXONOMY.md](GENE_QUESTION_TAXONOMY.md) | active | A gene's total incoming regulation. Backs `apm-gene-question`. |
| [MIRNA_GENOMIC_CONTEXT_AXIS.md](MIRNA_GENOMIC_CONTEXT_AXIS.md) | active | Strand-aware host classification; the per-edge host-gene lens. |

## 2. The learned model

> **The converged statement in one line:** ONE dense learned-τ² non-negative Bayesian posterior per
> gene, TWO readouts (π≡1 dense → coupling/attribution; evidence-π → discovery). **Lasso RETIRED.**
> Canonical: [`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`](LEARNED_MODEL_DISCOVERY_SYNTHESIS.md) §6b.
> Anything describing lasso-as-primary, "Bayes = uncertainty layer only", or Decision A–J / Bars is
> **pre-convergence (≤2026-07-05) and superseded.**

| doc | currency | role |
|-----|----------|------|
| [LEARNED_MODEL_DISCOVERY_SYNTHESIS.md](LEARNED_MODEL_DISCOVERY_SYNTHESIS.md) | active | **Capstone** — §6b = the converged model in one line. |
| [LEARNED_MODEL_METHODS.md](LEARNED_MODEL_METHODS.md) | active | Formula-level spec + **§18 estimator selection** (job → estimator → chosen-over). **Twin of** `_FORMAL` — always edit BOTH. ✅ Absorbed `ESTIMATOR_MAP` 2026-07-17. ⛔ §8 (CN causality) and 🔬 §5 (canonical attribution) now carry banners — read them. |
| [LEARNED_MODEL_METHODS_FORMAL.md](LEARNED_MODEL_METHODS_FORMAL.md) | active | The LaTeX twin. |
| [LEARNED_MODEL_RATIONALE.md](LEARNED_MODEL_RATIONALE.md) | active | The *why* behind every § of METHODS — incl. **§15/§15a identity-vs-magnitude** (⭐ why `share` is NOT identity), **§4a collinearity/identifiability**, and the **channel frame**. ✅ Absorbed the attribution trio + COLLINEARITY + CHANNEL_FUSION's live frame 2026-07-17. |
| [LEARNED_MODEL_VALIDATION.md](LEARNED_MODEL_VALIDATION.md) | active | Validation dossier. ⚠ §1 still carries the MH-127 restatement, not MH-134/135/137. |

## 3. Findings — the record is the REGISTRY

> ✅ **The six `MH12x_*.md` docs were archived 2026-07-17** — the per-finding-doc pattern is retired and now
> **banned** (`../CLAUDE.md` §3). Every finding was verified to have a registry row first (MH-125/128/130 were
> backfilled). Read [DISCOVERY_REGISTRY.md](DISCOVERY_REGISTRY.md); the docs below are deep-dive appendices in
> [`archive/`](archive/README.md), kept for provenance.

| doc | currency | what it holds |
|-----|----------|---------------|
| [MH124_ANTICOUPLING_VALIDITY.md](archive/MH124_ANTICOUPLING_VALIDITY.md) | ⚠ **mixed** | §4b (the attribution estimand argument) **STANDS**. **§5/§5b are VOID — do not cite** (`pi_causal` not exogenous). |
| [MH125_WHAT_SURVIVES.md](archive/MH125_WHAT_SURVIVES.md) | ⚠ mixed | The adversarial re-audit. Its exogenous leg is retracted (MH-133). |
| [MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md](archive/MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md) | active | **Cite this for attribution rank.** Ploidy = a calibration fix, not a power fix. Carries the w-circularity gate. |
| [MH127_DECOY_MODEL_GENE_BUDGET.md](archive/MH127_DECOY_MODEL_GENE_BUDGET.md) | ⛔ **superseded** | Magnitudes superseded by **MH-137**. Its "in-cohort tie" verdict is **reversed**. |
| [MH128_DECOY_RESOLUTION_AND_CN_GOLD.md](archive/MH128_DECOY_RESOLUTION_AND_CN_GOLD.md) | ⛔ superseded | §3's CN leg retracted by MH-133. |
| [MH130_GENE_LENS.md](archive/MH130_GENE_LENS.md) | ⚠ **half** | The **atlas** half (ceilings, 17.6% unmeasurable, 48.1% one-family) **STANDS**. The **gap** half sits on a discredited decoy. |

## 4. Findings surface + detail reports

[BIOLOGICAL_SYNTHESIS.md](BIOLOGICAL_SYNTHESIS.md) is the **surface** — read it for "what did we learn"
(biology only; modeling conclusions live elsewhere). The three detail reports were demoted to
[`reports/`](reports/) on 2026-07-17: they hold real results but overlap the surface.

| doc | currency | what it is |
|-----|----------|------------|
| [BIOLOGICAL_SYNTHESIS.md](BIOLOGICAL_SYNTHESIS.md) | active | **Biology-first surface.** Start here for findings. |
| [reports/REPORT.md](reports/REPORT.md) | detail | Detailed TCGA-BRCA results. |
| [reports/LANDSCAPE_REPORT.md](reports/LANDSCAPE_REPORT.md) | ⚠ detail | Cross-state landscape. **Its per-edge FDR counts rest on the uncalibrated null** (MH-123/124) — banner pending. |
| [reports/MIRNA_CNV_DOSAGE_REPORT.md](reports/MIRNA_CNV_DOSAGE_REPORT.md) | detail | miRNA CNV dosage. |
| [reports/PROGRESSION_LANDSCAPE_REPORT.md](reports/PROGRESSION_LANDSCAPE_REPORT.md) | detail | **⭐ Consolidated progression-landscape characterization** (17-thread map: structure · mechanism · players · validation). Descriptive; sync-partner `BIOLOGICAL_SYNTHESIS.md §12`; tested spine = MH-158/160/162/163. |
| [DCIS_EV_SYNTHESIS.md](DCIS_EV_SYNTHESIS.md) | active | DCIS/EV arc (MH-48..56). ⭐ Its MH-55 result **independently corroborates MH-114** 16 days early. |

## 5. Channels & axes

> ⛔ **2026-07-17: two dead-channel docs were ARCHIVED** — `STATE_CHANNEL_PLAN` (measured and cancelled, τ≈0)
> and `CPTAC_PROTEIN_CHANNEL_PLAN` (`βᵗ` falsified at n=101). **Their measured survivors were extracted first**:
> METHODS §1a/§18a and VALIDATION §6.1–§6.5. See [`archive/README.md`](archive/README.md) for each address.
> ⚠ **`CN_INSTRUMENT.md` STAYS ACTIVE.** Its instruments are retracted, but its **design** — the exclusion
> restriction, multi-IV + Hansen-J over-ID, the CN/expression asymmetry — is what a **revival** builds on, and
> that revival is the programme's highest-value open item (board §A). *A dead result is archivable; a live
> design is not.*

| doc | currency | what it is |
|-----|----------|------------|
| [CN_INSTRUMENT.md](CN_INSTRUMENT.md) | ⚠ **design LIVE, instruments DEAD** | The CN-locus instrument's design + identifiability argument. ⛔ Both instruments retracted (MH-124r/126, MH-133) — **do not cite it as causal evidence**. ✅ **Kept for the REVIVAL**: §2 exclusion restriction · §5 multi-IV + Hansen-J over-ID · §6 the CN/expression asymmetry · §7 the corrected architecture · §9 built-vs-gaps. The estimator spec (F>10 ∧ T1-clean, soft F-weight) moved to `LEARNED_MODEL_METHODS §8a`. |
| [SEED_SITE_INSTRUMENT_DESIGN.md](SEED_SITE_INSTRUMENT_DESIGN.md) | 🔧 **design ONLY, nothing built** | The **target-side** edge-existence instrument (board §A↔§E): a seed-disrupting 3′UTR SNV/editing/APA perturbs the *site*, not arm dose, so it is blind to the confounds that killed CN. Carries the causal DAG, the `site×dose` interaction identification, and the exclusion-violation ledger. P0a run (128 overlap events → SNV-alone is set-level-only). No registry row. |

| doc | currency | what it is |
|-----|----------|------------|
| [ISOMIR_AWARE_MODELING.md](ISOMIR_AWARE_MODELING.md) | active | isomiR-corrected X_fam. Phase 2–4 built, **default-OFF** (coupling wash). |

## 6. Forward planning — ONE doc

> ✅ **2026-07-17: the five parked design docs were archived** — each planned for something since retired,
> answered, or blocked on an unstarted step. Their live threads are absorbed into the board's **§Z**.

| doc | currency |
|-----|----------|
| **[PROGRAM_FORWARD_BOARD.md](PROGRAM_FORWARD_BOARD.md)** | **active, current through MH-137.** The single forward doc. Opens with a **🔨 DO THESE FIRST** block. ✅ **Consolidated 2026-07-17**: absorbed `WHATS_NEXT.md` + `LEARNED_MODEL_WHATS_NEXT.md` (86 KB across three docs → 32 KB in one); both are now in [`archive/`](archive/README.md). |

## 7. Subfolders

| folder | what |
|--------|------|
| [archive/](archive/README.md) | **Consumed / superseded.** Nothing here describes current state. 9 docs + a README explaining each. |
| [derived/](derived/) | Generated renditions (PDF/DOCX) of the canonical docs. Don't edit — regenerate. |
| [decisions/](decisions/) | Locked decision memos. ⭐ **Rescued from the gitignored `output/`** — they were untracked. |
| [ARCHITECTURE.md](ARCHITECTURE.md) | **⭐ THE MAP** (generated) — axes × models × analyses × results, joined by reference. Start here for "how does it all fit". Regenerate: `analyses/ops/gen_architecture.py`. |
| [ARCHITECTURE_PLAN.md](ARCHITECTURE_PLAN.md) | **The plan for the canonical axes×models×analyses×results map** (2026-07-18). A materialized VIEW over the one-home docs; generated, non-redundant. Draft axis taxonomy inside. |
| [handoffs/](handoffs/) | Session handoffs. **`HANDOFF_DISCOVERY_NULL.md` (2026-07-17) is the current one** — discovery lane / site-free null / FDR (MH-154/155). ⚠ `HANDOFF_PROTEIN_AND_COMPOSITION.md` is superseded on 3 counts — see STATE_OF_PLAY. |

## 8. Reference & renditions

| doc | what |
|-----|------|
| [MODELING_FRAMEWORK_EXTERNAL.md](derived/MODELING_FRAMEWORK_EXTERNAL.md) | derived — external/paper version, NotebookLM source. |

---

## Markdown outside `docs/`

- **Roots:** [AGENTS.md](../AGENTS.md) · [CLAUDE.md](../CLAUDE.md) · [README.md](../README.md)
- **method_dev/**: see [method_dev/README.md](../method_dev/README.md) — 4 arcs (aggregate_pressure, arm_expression, site_ladder, landscape).
- **Grant/**, **presentation/** — ⚠ `presentation/` F25 + deck are **stale** (figure from 2026-06-24; prose fixed 2026-07-12).

## Known structural debt

1. **MH-131…137 have no doc and (for 133/134/136) no ledger row.** The registry is ahead of everything.
2. **Three forward docs** need consolidating into one.
3. **`DISCOVERY_REGISTRY.md` (338KB) and `ANALYSIS_RUN_LEDGER.md` (262KB)** are append-only monoliths;
   single cells run 3,000+ chars and retractions live inside the rows they retract.
4. **The doc protocol fans out with no consolidation rule** — "update EVERY applicable canonical home"
   means a finding lands in N places and its later retraction lands in 1. That is the mechanism behind
   MH-38 being dead in MH-114's row and alive in its own.
