# Handoff — CPTAC/protein axis + the COMPOSITION machinery (session of 2026-07-12)

> **Goal:** hand a fresh session the VERIFIED state of the protein axis (MH-103…109) and the composition/retention
> machinery (MH-110…113), the five approved follow-ups, and — critically — the list of claims I made and then had
> to RETRACT, so they are not re-derived.
> **What belongs here:** verified state, live outputs, approved next steps, retractions. NOT method detail
> (→ the registry/ledger rows named below), NOT the original plan (→ `CPTAC_PROTEIN_CHANNEL_PLAN.md`).
> **Update trigger:** when a follow-up lands or a claim here is superseded.
> **Sync-partner:** `CPTAC_PROTEIN_CHANNEL_PLAN.md`, `DISCOVERY_REGISTRY.md` (MH-103…113), `PROGRAM_FORWARD_BOARD.md`.

---

## ⭐ THE BOTTOM LINE (MH-114 — read this before anything else)

**The CPTAC bulk-PROTEIN validation of miRNA repression is ~57% COMPARTMENT ARITHMETIC. This is PROVEN, not
argued:** degree-preserving **SHUFFLED non-edges** — pairs that *cannot* repress, by construction —
reproduce the real curated edges' compartment gradient **EXACTLY** (**−0.1290 vs −0.1272**). Bulk protein and
bulk miRNA are both compartment-weighted averages, so a miRNA merely *living in a different cell type* than a
target anti-correlates with it with **zero** cell-autonomous repression.

- **Composition adjustment removes ARTIFACT, not biology** (it kills 90% of that gradient in real AND fake
  edges alike). ⛔ **I first concluded the opposite** — see trap #1 below.
- **The real effect SURVIVES**, as the compartment-BLIND real-minus-shuffled contrast: **≈−0.011** unadjusted,
  **≈−0.0065** adjusted (MW **p=2.0e-03**), *the same in both orientations* — as true repression must be.
- **The FDR collapse is a POWER WALL, not evidence of absence.** Genes **34 → 3**, orphan candidates
  **594 → 23** — but both sets stay significantly negative *as a set* (genes MW **p=1.8e-08**; orphans MW
  **p=4.2e-141**, 87% keep sign). At the true \|ρ\|≈0.14, n=101 has **~25% power**.
- ⇒ **The protein layer validates miRNA repression as a SET, at ~40% of the advertised magnitude. It does NOT
  validate individual edges at n=101.**

**TWO METHODOLOGICAL TRAPS, both of which fooled me (now axiom 4 in `mirna_hallmark/CLAUDE.md`):**
1. **A MEAN-based shuffled null is blind to a SIGN-SYMMETRIC artifact.** Random pairs are equally likely
   same-compartment (artifact **+**) or opposite (artifact **−**) ⇒ **it cancels in the mean** (shuffled mean
   ρ = +0.0001) ⇒ I wrongly concluded "adjustment removes biology". **STRATIFY THE NULL BY THE CONFOUND'S OWN
   AXIS.**
2. **"It replicates in an independent cohort" is NOT evidence of biology when both cohorts SHARE the
   confound.** Adjustment *lowered* prospective↔TCGA-105 concordance (**+0.263 → +0.123**) — because the
   shared CAF-rich stroma **was** the reproducible part. **The artifact replicates.**

**The identification strategy is an ASYMMETRY the confound is blind to:** compartment arithmetic doesn't know
about 3'UTR sites; repression doesn't know about cell type.

---

## 0. READ THIS FIRST — I was wrong a lot. Do not re-derive these.

This session produced **9 self-retractions**. Every one was caught by *running a check*, never by reasoning.
Treat my intuitions in this domain as unreliable; **test, don't argue.**

| # | I claimed | Truth | Where |
|---|---|---|---|
| 1 | βᵗ (translational latent) is real — OOF 9/9, PTEN z=3.1 | **LEAK** — the mediator `a` was fitted on ALL samples, then OOF-ed. Fit inside the fold ⇒ signal gone. **βᵗ FALSIFIED.** | MH-103 |
| 2 | "The discordance link forces HMC" | **False.** Discordance is ADDITIVE ⇒ linear ⇒ conjugate ⇒ **stays on Gibbs.** | plan §3.1 |
| 3 | STAR is the right CPTAC mRNA | **Wrong** — LinkedOmics is better on two non-circular measures. Default flipped. | MH-104 |
| 4 | "Unbiased `a` is up to 38% off" | **Immaterial genome-wide** (paired median bias −0.001; 52.4% coin-flip). From a 9-gene panel. | MH-105 |
| 5 | "miRNA effect reaches protein only where protein tracks mRNA" (the "mechanistic clincher") | **Does not hold.** Continuous test: a_IV ρ=+0.03 (p=0.65); a pure nuisance (instrument-strength F) predicts equally well. | MH-109 |
| 6 | "Buffering is 2.4× stronger under IV" | **Selection artifact** (OLS-on-all-genes vs IV-on-F>10). Same set: OLS −0.157 vs IV −0.136 — IV slightly WEAKER (2SLS is less efficient). | MH-108-C |
| 7 | "Bulk protein is ~3× more compartment-exposed than mRNA" | **UNIT+COHORT artifact.** Same cohort/unit/C: mRNA 38.5% vs protein 33.7%, p=0.99. **No layer effect.** | MH-113 |
| 8 | "C-parity ⇒ deconv=True everywhere" | **Would IMPOSE an over-control.** `deconv` is deliberately NOT a default — composition may BE the signal. | MH-107-C |
| 9 | "Composition adjustment removes BIOLOGY (over-control)" — from a shuffled null whose shuffled mean was ≈0, and from adjustment *lowering* cross-cohort concordance | **BACKWARDS.** The mean CANCELS a sign-symmetric artifact, and the artifact REPLICATES across two cohorts that share the stroma. Stratified by orientation, **fake edges reproduce the real gradient exactly** ⇒ adjustment removes **ARTIFACT**. | MH-114 |

**⚠ Also: `retention` rates are n-DEPENDENT** (TCGA n=1041 → 14%; CPTAC n=101 → 38.5%, same layer). It is an
unstable RATIO when β_core is small. **Never compare retention rates across cohorts of different n.**

---

## 1. VERIFIED STATE — the protein axis

**The plan's centrepiece is DEAD.** βᵗ (a per-family translational-repression latent) is **not supported at n=101**,
in four leak-free framings incl. the max-power aggregate (BH q<0.10 in **1/17** genes; **PTEN dead**, p=0.82).
Coherent with MH-34: the translational residual is a **~1%-of-genes** phenomenon; n=101 resolves ~one gene's worth.
**There was never a βᵗ field. Only n fixes this.** (MH-103)

**What DID land:**
- **Bar-5 PASSES** (581 genes): the learned M carries to an independent cohort AND layer (p=1.3e-4). ⭐ **The loss is
  at the COHORT boundary (~80%), not the LAYER boundary (~0%).** Stratified: retention rises with regulator count and
  in-cohort coupling strength (10–20 families → 0.287; ceiling ρ<−0.25 → 0.217, p=3.5e-5) — **a MIXTURE, not a
  uniform loss.** (MH-104, MH-106)
- **Dose-delivery TRANSPORTS** (unlike β): same dominant member 84.5% vs 43.7% chance; member-share r=0.991.
  The "inherited δ" assumption is **verified**. (MH-104)
- **`a_g` (mRNA→protein slope) causal via CN-IV**: `a_OLS` 0.561 → **`a_IV` 0.893** ⇒ mRNA reliability λ=0.63.
  Protein tracks mRNA ~1:1 causally. **Complex buffering is CAUSAL** (monomer 0.943 → ≥31 subunits 0.526, p=8.8e-20).
  ⚠ **NOT NOVEL** — this is Gonçalves 2017 (*Cell Systems*) + Csárdi 2015 (*PLoS Genet*) re-derived. (MH-105/106/108)
- **`locus_cn_cptac` BUILT** (r=0.997 vs ASCAT) but **WEAK as an instrument** (F>10 for only 8.6% of arms;
  5 usable families). Fidelity ≠ power.
- **MH-107: the CPTAC protein validation never computed the composition-RETENTION tag** ⇒ a compartment-driven result
  (EMT/ZEB1) was presented as cell-intrinsic. **Honest verdict (MH-107-C): ~1/3 of hits CELL-INTRINSIC, ~38%
  COMPARTMENT ARITHMETIC, ~30% partial.** NOT "the result is dead". The **edge is real** (miR-200→ZEB1, TCGA mRNA,
  composition-adjusted, n=1041, **p=8.7e-12**) — composition INFLATES it, does not invent it.

---

## 2. VERIFIED STATE — the composition machinery (the session's real yield)

**`learned/readouts.py` — THE FOUR READOUTS, genome-wide** (closes the board's "canonical M rerun"):
1,549 genes × 5,117 (gene,family) rows, 8.3 min on 8 workers. ONE posterior, read four ways:
**COUPLING** (β, β_sd) · **ATTRIBUTION-magnitude** (β) · **ATTRIBUTION-identity** (Bayesian Shapley `share ± share_sd`;
for a linear aggregate Shapley(f)=β_f exactly ⇒ per-draw share β_f/Σβ) · **IDENTIFIABILITY** (|z|>2, 29.9% of
families) · **DISCOVERY** (evidence-π PIP, 20.6%). Median top-share **0.973** ⇒ a gene's repression is highly
**concentrated in ONE seed family**. (MH-110)

**⚠ `deconv` is NOT a default — composition may BE the signal** (`METHODS §1`: *"gate/coupling = core only unless
deconv requested"*). `readouts.py` runs **BOTH** blocks and reports **`retention = β_deconv/β_core`**
(cell_intrinsic ≥0.7 / partial / composition_explained <0.4). Edge-level: **64.4% / 21.6% / 14.0%**.

**RETENTION IS VALIDATED** (not asserted): per-family compartment loading predicts it **monotonically** (least-locked
→ retention 0.980, 5.4% explained; most-locked → 0.677, **28.4%**). Canonical case passes: **miR-200 (the textbook
EPITHELIAL miRNA) median retention 0.594 vs 0.874 for others, p=9.1e-12.** (MH-110)

**PINPOINTED (MH-111):**
- **⚠ The Wu-major fractions are a SIMPLEX** (sum = exactly 1.0) ⇒ **R²(Cancer-Epithelial ~ the 8 "held-out"
  lineages) = 1.0000** ⇒ the design's "hold the malignant fraction out" is **mathematically an ILLUSION**.
  **BUT BENIGN:** purity contributes **0.9%** of the loss (CPE in core C already absorbs it); the **stromal MIX**
  contributes **33.1%**.
- **WHICH CELLS:** **CAFs 10.7% · Endothelial 6.3% · Myeloid 5.6% · PVL 4.1%** · Normal-Epi 1.5% · B-cells 0.6% ·
  **T-cells 0.0% · Plasmablasts 0.0%** ⇒ a **MESENCHYMAL/VASCULAR/MYELOID axis, NOT lymphocytes.**
- **GENE ALLOCATION (Shapley is additive ⇒ exact):** `gene_retention = Σβ_deconv/Σβ`. 1,502 genes: median 0.800;
  cell_intrinsic 60.5% / partial 29.5% / composition_explained 10.0%. **The worst single family carries 95% (median)
  of a gene's composition loss.**
- **ARM → FAMILY: an INFORMATIVE NEGATIVE.** Within-family compartment-lock spread = **0.0121** ⇒ seed-family members
  are compartment-HOMOGENEOUS; all arm summaries predict retention identically ⇒ **arm-first buys NOTHING for
  composition; the FAMILY is the right unit.**

**⭐ THE ARCHITECTURE PRINCIPLE (the session's most generalisable finding):**

| confound | localises to | why | the handle |
|---|---|---|---|
| proliferation (MH-100) | **ARM** (miR-93-5p, conc. 0.85) | a **LOCUS** property (intronic in MCM7) | host mRNA |
| host gene (MH-101) | **ARM** (dose-dominant) | a **LOCUS** property | host conditioning |
| **composition (MH-111)** | **FAMILY** (members homogeneous) | an **EXPRESSION** property | retention tag |

*Same-seed members sit at different LOCI (⇒ diverge on host/prolif) but are CO-EXPRESSED (⇒ agree on cell type).*
**The unit at which a confound acts predicts which tool can fix it.**

**RETENTION IS A CLEAN ARTIFACT FILTER (MH-112):** of 537 `composition_explained` edges, **486 (90.5%) are
OPPOSITE-compartment** (miRNA one side, target the other ⇒ bulk-mixing arithmetic) vs a **52.8% coin-flip baseline**
for cell-intrinsic edges (**Fisher OR=0.09, p=1.6e-87**). ⇒ **The "hidden stromal network" hypothesis is REFUTED —
retention flags arithmetic, not biology.** **Residue = 47 both-STROMAL candidate compartment-intrinsic edges; top hit
`RBPJ ← miR-21-5p/590-5p`** (miR-21 = canonical CAF miRNA; RBPJ = Notch effector).

---

## 3. THE FIVE APPROVED FOLLOW-UPS (user-approved 2026-07-12, in priority order)

- **A. Wire `retention` + compartment ORIENTATION into the discovery lane / dossier.** The evidence-π PIP is blind to
  composition. An edge that is *discovered* **AND** `composition_explained` **AND** opposite-compartment is almost
  certainly an artifact (p=1.6e-87 backing). **Cleans the 568-novel-candidate list and the gold/dossier sets.**
- **B. Emit a per-(gene, family) COMPOSITION DOWN-WEIGHT** — exactly parallel to MH-99's `coding_down_weight`
  (→ inflate s²). Since **95%** of a gene's loss sits in ONE family, this is surgical, not a blanket filter.
- **C. The 47 compartment-intrinsic candidates** (led by miR-21→RBPJ). Small, specific, biologically motivated —
  the one place stromal regulation may be real. A targeted look, not a program.
- **D. Reconcile with the DCIS/EV "stroma reframe"** (`DCIS_EV_SYNTHESIS.md`) — that arc already concluded stroma
  mattered; we now have quantitative per-edge compartment orientation. May explain or sharpen it.
- **E. Document the composition machinery in the METHODS twins** (`LEARNED_MODEL_METHODS.md` +
  `_FORMAL.md` — **ALWAYS edit BOTH**), and **correct §1's Cancer-Epithelial rationale** (the hold-out is a simplex
  illusion; the practice is safe because CPE absorbs purity).

**✅ MH-107 CORRECTNESS DEBT — CLOSED (MH-114).** `cptac_validation` re-run for **BOTH** cohorts into the
canonical paths (pre-fix preserved at `output/cptac_validation_PRE_MH107/`); `ood_protein` re-run (hub survives
7/7). **A SECOND composition-blind lane was found and fixed:** `analyses/cptac/cptac_orphan_discovery` built
`purity+cin` straight off the clinical frame and never called `_covariates` ⇒ the whole **orphan nomination**
deliverable was screened with no composition block (**594 → 23** candidates). Both lanes now route through the
shared `_covariates`. A **manifest provenance bug** was also fixed (`conf_label` was a literal written *before*
the composition join ⇒ the manifest said "purity + CIN" while 11 covariates were fitted).

**⚠ STILL OPEN:**
- **`dossier.tier3_protein`** (MH-84's "268 triple-validated" — **27% sign-flip**) — code is fixed, **not yet
  re-run**.
- **`buffa_validation`** consumes `cptac_validation/orphan_discovery/` ⇒ the **"108 orphans triple-validated"**
  claim (talk, Validation 3) rests on the OLD orphan list. **Re-run it.**
- **`presentation/`** — F25 + `talk.qmd` were UPDATED to the honest story (the old slide claimed "53/1143 genes
  FDR-neg; ZEB1 ρ=−0.55" and annotated **CALD1/MMP2, which are CAF markers**). **Figures need regenerating**
  (`make_figures.py`) and the deck re-rendering.
- The talk's **"cell-composition confound pilot: 24 → 17 survive"** (71%) is contradicted by the full run
  (594 → 23 = 4%). Worth reconciling — `analyses/cptac/cptac_orphan_confound_pilot.py`.

---

## 4. CODE / OUTPUTS

**New modules:** `learned/readouts.py` (the four readouts + retention) · `learned/cptac_data.py` (CPTAC Model-3
assembly; **default RNA source = `linkedomics`**) · `scripts/analysis/build_locus_cn_from_gene_cn.py` ·
`scripts/cptac/build_cptac_cibersortx_mixture.py`.

**Fixed (shared spine):** `data._arm_name_map` global-cache **POISONING** (a CPTAC-first call silently handed TCGA the
CPTAC map; fingerprinted on index+shape; **verified strict no-op**) · `gauge.beta_table` dropped `.N` arms
**ASYMMETRICALLY** (TCGA 19 / CPTAC 58 ⇒ the cross-cohort gauge compared DIFFERENT arm sets) ·
`cptac_validation._covariates` / `ood_protein` / `dossier.tier3_protein` now composition-aware (`*_RAW` kept).
**⚠ The state session's gauge numbers (0.6%, a=1.125) predate the arm fix — worth re-running.**

**Key outputs** (`output/learned/`): `readouts_{edges,genes}.tsv` · `family_compartment_loading.tsv` ·
`arm_compartment_loading.tsv` · `family_arm_dose_compartment.tsv` · `composition_pinpoint.tsv` ·
`gene_composition_allocation.tsv` · `edge_compartment_orientation.tsv` · `layer_retention_cptac.tsv` ·
`bar5_ladder.tsv` · `cptac_a_causal_iv.tsv` · `cptac_a_by_cn_direction.tsv` · `cptac_propagation_slope_a.tsv` ·
`cptac_protein_composition_rerun.tsv` (carries `retention` + `composition_class`).
**Data:** `data/CORUM/humanComplexes.txt` (CORUM 5.3; ⚠ host TLS chain broken → fetched with verification relaxed,
payload validated — see manifest).

**Registry/ledger:** MH-103 · 104 · 105 · 106 · 107 · **107-C** · 108 · **108-C** · 109 · 110 · 111 · 112 · **113**.
