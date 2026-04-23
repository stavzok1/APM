# 05 — PhD timeline

Three scenarios, same scope. The only real lever is **how fast hypotheses can be batch-tested via the automated engine in L3**. Everything else is planning.

Assumes:
- Pipeline + Atlas are **done** (they are).
- Analysis layer (`analysis/`) exists partially; feature-builder contract not yet frozen.
- Agent-assisted hypothesis generation and testing is available (see `06_automation_plan.md`).
- Committee expects three first-author papers or one + thesis chapters of equivalent weight.

---

## Decision gates

Work is gated by five critical dependencies, in order:

| Gate | Meaning | Blocks |
|---|---|---|
| **G1 — Contract freeze** | L1 feature builder contract in `03_architecture.md` frozen, sample-key normalization in place | Everything downstream |
| **G2 — Consensus prior calibrated** | 13 HiChIP calibrated against subtype-consensus prior (`04_conservation_prior.md`) | All enhancer-gene hypotheses (most of the catalog) |
| **G3 — Extended panel integrated** | Tier-1 additions in `01_gene_panel_extended.md` (cGAS–STING, CIITA cluster, ERAP, lncRNAs) flowing through the pipeline | H48, H56, H50, H32–H34 |
| **G4 — Per-allele layer** | Allele-aware RNA + LOHHLA + haplotype-phased methylation built per `08_hla_allele_layer.md` | Family M (H46, H63–H65) |
| **G5 — External cohort access** | METABRIC + at least one ICB cohort (IMvigor210 / IMpassion130 / Liu / Van Allen) accessible and harmonized | H6, H18 replication, H59 / H64 validation, Aim 2 paper acceptance |

---

## Fast path (~12 months to defense-ready results)

| Month | Deliverable | Gated by |
|---|---|---|
| 0 | Research plan written (this folder), hypotheses catalogued | — |
| 1 | G1 met: L1 builders specified + sample normalization + parquet caching | — |
| 1–2 | G2 met: 13 HiChIP vs subtype-consensus calibration report, feature `*_resid` variants live | G1 |
| 1–2 | G4 started: allele-aware RNA builder + LOHHLA stand-up; Family M features begin flowing | G1 |
| 2 | G3 met: Tier-1 panel integrated into pipeline config + rerun outputs | G1 |
| 2–4 | Aim 1 atlas: variance decomposition for all 66 + extended; H1, H19, H38, H41, H44 results; Figure 1–3 | G1, G2 |
| 2–4 | **Wet-lab W2/W6 pilot** (allele-specific reporters): quick methodological check | G1 |
| 3–4 | G4 met: per-allele outcome variables live; H46 and H63 first pass | G1 |
| 4–5 | Aim 1 paper draft; preprint | — |
| 5–7 | Aim 2 scores: CRML, BCD_g, R_g, MIR, PRD, IVI, COMB_HLA constructed; Cox PH TCGA; four-class hot/cold (H59); H4, H9, H12, H14, H46, H18, H63 | G2, G4 |
| 5–9 | **Wet-lab W1 (CRISPRi) + W3 (demeth/EZH2i reversal)** run in parallel | — |
| 7–8 | G5 met: METABRIC + ICB cohort ingested | — |
| 8–9 | Aim 2 external replication; TNBC hot/cold tested in ICB cohort (H59, H64); H5, H6 | G5 |
| 9–10 | Aim 2 paper draft with wet-lab figures; preprint | W1 + W3 readouts |
| 10–11 | Aim 3 platform: DDR panel + EMT panel fingerprint; H17 | — |
| 10–11 | **HLA regulatory-loss methods paper** (standalone: H46 + H63 + 13-HiChIP calibration) | G2, G4 |
| 11–12 | Aim 3 paper draft; thesis narrative assembled | — |
| 12 | Defense-ready: 3–4 preprints, 1–2 submitted, wet-lab figures in place | — |

Risk: assumes engine fires cleanly through ~40 Tier-1 hypotheses with minimal re-scoping, G5 does not slip, and the wet-lab bench has dCas9-KRAB + NK92 or equivalent available. The HLA per-allele layer (Family M) fits comfortably into the fast path because most of the ingredients are already present.

---

## Realistic path (~18 months to defense-ready)

Slippage budgets added to the fast path:

- G1 + G2 take 3 months total (contract iterations + calibration sensitivity analyses).
- Aim 1 takes 4 months (reviewer rounds + supplementary).
- Aim 2 external replication adds 2 months (cohort harmonization always takes more than expected).
- Aim 3 takes 4 months (second panel requires its own literature dive and a light rerun).
- Thesis writing: 2 months overlap with aim 3.

| Quarter | Deliverables |
|---|---|
| Q1 | G1, G2 met; 20+ Tier-1 hypotheses run; Aim 1 figures drafted |
| Q2 | Aim 1 preprint; G3 met; CRML/BCD/IVI scores defined |
| Q3 | Aim 2 cohort validation; Aim 2 preprint; G4 met |
| Q4 | Aim 3 second panel run; Aim 3 preprint |
| Q5 | Revisions from reviewers + defense preparation |
| Q6 | Defense |

---

## Conservative path (~24 months, "enough" results)

"Enough" = committee-passable thesis with one published paper + two strong preprints.

- Aim 1 paper accepted mid-cycle.
- Aim 2 + Aim 3 as preprints.
- External ICB validation replaced by METABRIC-only survival validation.
- HLA per-allele layer (Family M) remains in-scope because typing is already available; H46 + H63 should land even on this path. H64 + H65 may be deferred.
- Wet-lab: **W2 or W6** (allele-specific reporter) as minimum backstop if CRISPRi / demeth experiments slip.
- Tier-2 / Tier-3 hypothesis families partially executed; the 15 lead-set hypotheses guaranteed.

Key change: if an entire aim stalls, the lead hypothesis set alone (15 hypotheses, all Tier-1 or Tier-2) is sufficient for a defensible thesis.

---

## Hypothesis throughput model

Let:

- \( N \) = number of hypotheses in the lead + priority-1/2 set ≈ 35
- \( H \) = human-in-the-loop hours per hypothesis = 4 (interpretation + figure + writing)
- \( E \) = engine compute cost per hypothesis = 2 hours amortized

If `fit_mediation` and `fit_regression` are productionized with stable feature caches (`analysis/output/features/*.parquet`), 35 hypotheses × 6 hours = 210 hours ≈ **6 weeks full-time** of actual analytical work, excluding cohort ingestion and writing.

Meaning: the rate-limiting step is not statistical testing. It is **writing and external cohort ingestion**. Budget accordingly.

---

## Writing strategy

Three parallel documents from month 3 onward:

| Document | Owner format | Updates |
|---|---|---|
| Aim 1 manuscript | LaTeX in `writing/aim1/` | Every completed hypothesis adds a paragraph + figure panel |
| Aim 2 manuscript | LaTeX in `writing/aim2/` | Every score builder commit triggers a method-section update |
| Thesis outline | Markdown in `writing/thesis/` | Section per aim, grows by copy-edit from papers |

Avoid writing thesis separately from papers. Write papers first, thesis absorbs them.

---

## Committee checkpoint targets

| Month (realistic) | Checkpoint |
|---|---|
| 3 | Hypothesis catalog signed off; feature contract frozen |
| 6 | Aim 1 preprint to committee + draft of Fig 1–3 |
| 12 | Aim 2 preprint + external cohort results |
| 18 | Aim 3 preprint + full thesis draft |
| 20 | Committee revisions |
| 24 | Defense |

---

## What makes the timeline fast

1. Pipeline is already the hard part; element-resolved tables exist.
2. Hypotheses are rank-ordered in `06_automation_plan.md`; Tier-1 can batch-execute.
3. Scores compose from builders; no bespoke per-hypothesis code.
4. Writing is incremental, not end-loaded.
5. Conservation-prior methodology upgrades small-sample modules without waiting for new data.

## What makes the timeline slip

1. External cohort access (G4) is the single most common slippage source.
2. Tier-3 hypotheses (H6, H16, H24, H46) are honest unknowns; budget for two of four to succeed.
3. Reviewer rounds: preprint promptly; don't treat revisions as serial.
4. HLA-typing pipeline (for H46) is its own project — either include as a lemma or defer.
