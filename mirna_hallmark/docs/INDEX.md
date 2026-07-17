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
| [MODELING_FRAMEWORK.md](MODELING_FRAMEWORK.md) | active | Conceptual + implemented framework. "What are we modeling and why." |
| [FORMULAS.md](FORMULAS.md) | active | Code-faithful exact formula for every quantity. The *how*. |
| [METHODS.md](METHODS.md) | active | Methods narrative — data → pressure → gate → coupling. |
| [DATA_SOURCES.md](DATA_SOURCES.md) | active | Every dataset read, and its preprocessing basis. |
| [ANALYSES_CATALOG.md](ANALYSES_CATALOG.md) | living | One row per analysis component. Orient here before grepping. |
| [MODULE_MAP.md](MODULE_MAP.md) | active | Which module does what. |
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
| [LEARNED_MODEL_METHODS.md](LEARNED_MODEL_METHODS.md) | active | Formula-level spec. **Twin of** `_FORMAL` — always edit BOTH. |
| [LEARNED_MODEL_METHODS_FORMAL.md](LEARNED_MODEL_METHODS_FORMAL.md) | active | The LaTeX twin. |
| [LEARNED_MODEL_RATIONALE.md](LEARNED_MODEL_RATIONALE.md) | active | The *why* behind every § of METHODS. |
| [LEARNED_MODEL_VALIDATION.md](LEARNED_MODEL_VALIDATION.md) | active | Validation dossier. ⚠ §1 still carries the MH-127 restatement, not MH-134/135/137. |
| [LEARNED_MODEL_ESTIMATOR_MAP.md](LEARNED_MODEL_ESTIMATOR_MAP.md) | ⚠ **part-stale** | **TRAP** — billed canonical, but self-dates 2026-07-06 and its **main table is pre-convergence**. Only its banner UPDATEs + the pooled-HE policy are current. |
| [LEARNED_MODEL_WHATS_NEXT.md](LEARNED_MODEL_WHATS_NEXT.md) | ⚠ stale | Cluster-local TODO. Predates MH-133…137. |
| [LEARNED_REGULATORY_MATRIX_DESIGN.md](LEARNED_REGULATORY_MATRIX_DESIGN.md) | active | Design note on learning M. |
| [COLLINEARITY_AND_IDENTIFIABILITY.md](COLLINEARITY_AND_IDENTIFIABILITY.md) | active | Why identity is hard. Hierarchical δ-pooling = parked. |

## 3. Findings — deep dives (the registry is the record; these are appendices)

| doc | currency | what it holds |
|-----|----------|---------------|
| [MH124_ANTICOUPLING_VALIDITY.md](MH124_ANTICOUPLING_VALIDITY.md) | ⚠ **mixed** | §4b (the attribution estimand argument) **STANDS**. **§5/§5b are VOID — do not cite** (`pi_causal` not exogenous). |
| [MH125_WHAT_SURVIVES.md](MH125_WHAT_SURVIVES.md) | ⚠ mixed | The adversarial re-audit. Its exogenous leg is retracted (MH-133). |
| [MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md](MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md) | active | **Cite this for attribution rank.** Ploidy = a calibration fix, not a power fix. Carries the w-circularity gate. |
| [MH127_DECOY_MODEL_GENE_BUDGET.md](MH127_DECOY_MODEL_GENE_BUDGET.md) | ⛔ **superseded** | Magnitudes superseded by **MH-137**. Its "in-cohort tie" verdict is **reversed**. |
| [MH128_DECOY_RESOLUTION_AND_CN_GOLD.md](MH128_DECOY_RESOLUTION_AND_CN_GOLD.md) | ⛔ superseded | §3's CN leg retracted by MH-133. |
| [MH130_GENE_LENS.md](MH130_GENE_LENS.md) | ⚠ **half** | The **atlas** half (ceilings, 17.6% unmeasurable, 48.1% one-family) **STANDS**. The **gap** half sits on a discredited decoy. |
| [CN_INSTRUMENT.md](CN_INSTRUMENT.md) | ⚠ mixed | The instrument's design + the 73% exclusion-failure verdict. ⛔ Both instruments now retracted (MH-133). |

## 4. Result reports

| doc | currency | what it is |
|-----|----------|------------|
| [BIOLOGICAL_SYNTHESIS.md](BIOLOGICAL_SYNTHESIS.md) | active | Biology-first surface. "What did we learn." |
| [REPORT.md](REPORT.md) | active | Detailed TCGA-BRCA results. |
| [LANDSCAPE_REPORT.md](LANDSCAPE_REPORT.md) | ⚠ | Cross-state landscape. **Its per-edge FDR counts rest on the uncalibrated null** (MH-123/124). |
| [MIRNA_CNV_DOSAGE_REPORT.md](MIRNA_CNV_DOSAGE_REPORT.md) | active | miRNA CNV dosage. |
| [DCIS_EV_SYNTHESIS.md](DCIS_EV_SYNTHESIS.md) | active | DCIS/EV arc (MH-48..56). ⭐ Its MH-55 result **independently corroborates MH-114** 16 days early. |

## 5. Channels & axes

| doc | currency | what it is |
|-----|----------|------------|
| [LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md](LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md) | parked | The ~60-axis fusion map. ⚠ Its §7/§L/§J "discordance forces HMC" is **wrong** — corrected to Gibbs. |
| [LEARNED_MODEL_STATE_CHANNEL_PLAN.md](LEARNED_MODEL_STATE_CHANNEL_PLAN.md) | ⛔ **CLOSED** | The state channel was **measured and cancelled** (τ≈0, info 0.6%). §10 = the post-mortem. **Do not rebuild.** Carries the settled panel decision. |
| [CPTAC_PROTEIN_CHANNEL_PLAN.md](CPTAC_PROTEIN_CHANNEL_PLAN.md) | ⛔ **centrepiece dead** | `βᵗ` falsified at n=101. ⚠ §6 still carries VOID PAM50 numbers contradicting §1. |
| [ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md](ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md) | active §0 | §0 = the identity-vs-magnitude frame. §1–§12 (softmax era) **RETIRED**. |
| [ATTRIBUTION_PRIMITIVE.md](ATTRIBUTION_PRIMITIVE.md) | ⚠ part-stale | The `attribute()` primitive — **still unbuilt**. Its §5/§8 "live bugs" are themselves stale. |
| [ATTRIBUTION_CONTEXT_AXIS.md](ATTRIBUTION_CONTEXT_AXIS.md) | parked | Bridge to a context axis. |
| [ISOMIR_AWARE_MODELING.md](ISOMIR_AWARE_MODELING.md) | active | isomiR-corrected X_fam. Phase 2–4 built, **default-OFF** (coupling wash). |

## 6. Forward planning

> ⚠ **Three competing forward docs, all predating MH-133…137. Consolidation pending.**

| doc | currency |
|-----|----------|
| [PROGRAM_FORWARD_BOARD.md](PROGRAM_FORWARD_BOARD.md) | ⚠ stale (2026-07-12) — the consolidated board; still the best of the three. |
| [WHATS_NEXT.md](WHATS_NEXT.md) | ⚠ stale — subproject-wide, pressure-arc era. |
| [PRESSURE_FUTURE_OPTIONS.md](PRESSURE_FUTURE_OPTIONS.md) · [PRESSURE_PROGNOSTIC_DESIGN.md](PRESSURE_PROGNOSTIC_DESIGN.md) · [PRIMARY_TUMOR_ROADMAP.md](PRIMARY_TUMOR_ROADMAP.md) | parked. ⚠ PROGNOSTIC_DESIGN is blocked on "METABRIC-full (EGA-pending)" — **no DAR was ever filed**. |

## 7. Subfolders

| folder | what |
|--------|------|
| [archive/](archive/README.md) | **Consumed / superseded.** Nothing here describes current state. 9 docs + a README explaining each. |
| [derived/](derived/) | Generated renditions (PDF/DOCX) of the canonical docs. Don't edit — regenerate. |
| [decisions/](decisions/) | Locked decision memos. ⭐ **Rescued from the gitignored `output/`** — they were untracked. |
| [handoffs/](handoffs/) | Session handoffs. ⚠ `HANDOFF_PROTEIN_AND_COMPOSITION.md` is superseded on 3 counts — see STATE_OF_PLAY. |

## 8. Reference & renditions

| doc | what |
|-----|------|
| [MODELING_FRAMEWORK_EXTERNAL.md](MODELING_FRAMEWORK_EXTERNAL.md) | derived — external/paper version, NotebookLM source. |
| [CRC_PORT_LITERATURE_SCAN.md](CRC_PORT_LITERATURE_SCAN.md) | reference — CRC-port literature scan. |

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
