# MH-127 — THE DECOY MODEL: does the learned gene budget beat a MATCHED FAKE regulator set?

> **Goal:** answer one question with a control that was missing from every previous benchmark — **when you
> learn weights β on a MATCHED, IMPOSSIBLE regulator set, does the gene budget `Σ_f β_f·X_f` predict the gene
> as well as the real one?** I.e. **is the learning doing anything beyond FITTING?** MH-115 compared a
> *fitted* β against *unfitted* abundance; that comparison rewards fitting, not curation. This document
> supplies the missing arm: **fitted-REAL vs fitted-FAKE.**
> **What belongs here:** the decoy design + its match QC, the four-rung ladder (real/fake × unfitted/fitted)
> in TCGA-OOF and in CPTAC (mRNA + protein, both C blocks), the stratified answer (regulator count · Shapley
> class · curated evidence weight), the corrected framing of MH-115, the adversarial verification log, and the
> strict separation of *underpowered* from *negative*.
> **What does NOT belong here:** any number not produced in this arc; any statistic transformed from another
> document rather than re-derived; biology claims (this is a MODELING conclusion — see
> `LEARNED_MODEL_VALIDATION.md` / `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`, **not** `BIOLOGICAL_SYNTHESIS.md`).
> **Update trigger:** any new decoy-controlled benchmark of the learned model; any change to the gene-field
> construction (`eval/coupling_grid._learned_M`, `learned/families.collapse_by_family`,
> `learned/attribution_eb._prep`); any re-run of the CPTAC transfer.
> **Sync-partner:** `ANALYSIS_RUN_LEDGER.md` (MH-127 row), `DISCOVERY_REGISTRY.md` (MH-127 / MH-127b /
> MH-127c, and the **MH-115 annotation**), `LEARNED_MODEL_VALIDATION.md` (§1 coupling gate + §4 negatives),
> `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` (narrative), `MH124_ANTICOUPLING_VALIDITY.md` (the per-edge decoy,
> of which this is the gene-level counterpart), `INDEX.md`.

**STATUS: COMPLETE.** Phase 1 (build the matched decoys + the strata), Phase 2 (the four-rung ladder, TCGA
OOF + CPTAC transfer), and **six adversarial lenses** (2× leakage, 2× match-quality, 2× stratum-artifact).
**Where a verifier refuted a result, the refutation is what is recorded below.**

---

## ⭐ THE VERDICT, FIRST SENTENCE

**IN-COHORT, NO: fitted-REAL does not beat fitted-FAKE.** On the abundance-matched, site-free, non-family,
non-cluster, non-collinear decoy design, the fitted decoy reaches **78.9%** of the real model's out-of-fold
gene-level anti-correlation over all 899 testable genes (**91.3%** under the structure-matched decoy on
multi-regulator genes), the paired gap is **Δρ = −0.0098 [95% CI −0.0213, +0.0018], real wins on 52.5% of
genes, Wilcoxon p = 0.211** (FAKE2: Δ = −0.0060, p = 0.577), and **the direction is stable across all 5 decoy
draws while never approaching significance.** The gap sits *below* the paired minimum detectable effect
(MDE = 0.0161) ⇒ **not significant, NOT proven zero.**

**BUT the curated set is not inert — it separates from the decoy in exactly two places, and both are
fit-free or out-of-cohort:**

1. **OUT-OF-COHORT TRANSFER (the real result).** Frozen TCGA β applied to CPTAC (n=101, never refit):
   the **REAL** field predicts CPTAC mRNA (**ρ = −0.0189, p = 8.9e-04**, 685 genes) and the **matched FAKE
   field does not** (**−0.0056, p = 0.19**); paired **REAL-vs-FAKE p = 0.030** (FAKE1) / **0.027** (FAKE2),
   and under FAKE2 it **survives composition adjustment** (Δ = −0.0126, **p = 0.018**). A fitted decoy buys
   in-cohort fit; it does **not** buy transfer.
2. **THE HIGH-EVIDENCE TERTILE (fit-free stratifier).** Splitting genes by *curated evidence weight* `w`
   (which no model touched), the top tertile shows real > fake **both fitted (Δ = −0.058, p = 3.0e-04,
   BH q = 0.0059) and unfitted (Δ = −0.058, p = 3.3e-04, q = 0.0059)** — the only stratum in the whole arc
   that separates. ⚠ It does **not** replicate under the structure-matched decoy FAKE2 (q = 0.20) ⇒ **a lead,
   not a result.**

**And the fitted DECOY beats REAL UNFITTED ABUNDANCE** (Δ = −0.0100, p = 0.054 over all genes;
**Δ = −0.0357, p = 3.6e-04 on multi-regulator genes**) — which is why **MH-115's headline must be restated
(§4).**

---

## §1 THE DESIGN

**THE FLAW BEING FIXED** (user, 2026-07-13): *"β is FITTED to the gene and abundance is NOT."* In bulk tumour
data every miRNA and every mRNA loads on a few global axes (proliferation, purity/composition), so **any**
predictor set regressed on Y will predict Y out-of-sample. The control is a **decoy model**: learn β over a
matched fake regulator set and compare **fitted-REAL vs fitted-FAKE**.

**UNIVERSE.** 1,421 pooled-HE genes (miRTarBase-HE ∪ TarBase-LT-functional, arm in `X`, arm above the
canonical `arm_expression` floor), n = 1,041 primary tumours, 555 candidate decoy arms.

**DECOY ELIGIBILITY** (the MH-124 §2 "strict decoy v3" cascade, reused not reinvented). A fake arm `a` may
replace a real regulator of gene `g` only if: `a` is expressed **and** covered by the genome-wide scanMiR scan
**and** has miRNA coords; **not** a real regulator of `g`; **SITE-FREE** on `g` (no TargetScan context++ site
**and** no scanMiR duplex, any row); **not** a seed-family mate of any real regulator; **not** a genomic-cluster
mate (50 kb single-linkage); and **|r| < 0.30 against every real regulator** (so it is not a collinear proxy).

**MATCHING = OPTIMAL ASSIGNMENT** (`scipy.optimize.linear_sum_assignment`, cost = |Δ median log2RPM| inside a
0.25 caliper, else 1e3+|Δ| — the LP provably maximises the number of in-caliper matches first). Two decoy
sets, 5 seeds each: **FAKE1** (arm→arm, abundance-matched) and **FAKE2** (FAKE1 + a family-level assignment
that reproduces `n_families` and the family *size profile* exactly, i.e. the collinear block structure).

**⭐ THE PHASE-1 FINDING THAT SHAPED EVERYTHING: a fully abundance-matched decoy is STRUCTURALLY IMPOSSIBLE at
the top of the abundance range.** On the full replacement set the abundance match fails hard (real median 7.07
vs fake 5.02 log2RPM, **MWU p = 6.4e-148**), and the deficit is monotone in abundance (real arms >10 log2RPM:
only 12% matched, mean Δ = −3.31). Abundant miRNAs are themselves curated regulators of most genes and carry a
scanMiR duplex on almost everything (78% of all arm×gene pairs), so an *abundant but impossible* regulator
barely exists in the transcriptome. **That bias points IN FAVOUR OF REAL** — i.e. the naive full-replacement
decoy would rig the test toward the exact artifact the task exists to detect.

**⇒ THE DEPLOYED DESIGN IS THE `in_caliper` TRIMMED SET:** keep only the real arms that have an in-caliper
twin, **symmetrically in both arms**. Abundance is then unbiased (**MWU p = 0.975**, signed Δ = −0.003,
median |Δ| = 0.071) and multi-regulator genes survive (899 genes ≥1 arm, 344 ≥2, 80 ≥4). **This changes the
estimand** — see §7.

**THE FOUR RUNGS** (all on the identical family-collapsed design, only the WEIGHTS and the ARMS differ):

| rung | arms | weights | question |
|---|---|---|---|
| **R1** | REAL | unfitted (uniform) | mere abundance on the curated regulators (the MH-115 baseline) |
| **R2** | FAKE | unfitted (uniform) | mere abundance on impossible regulators (the "global axes" floor) |
| **R3** | REAL | **fitted β** | the learned model |
| **R4** | FAKE | **fitted β** | ⭐ **the decoy model — the missing control** |

**THE GENE FIELD WAS BUILT FROM SCRATCH.** `eval/coupling_grid._learned_M` was **not** reused (it was under
audit at the time; it has since been fixed under MH-125). Construction, exactly:
`data.assemble_gene(gene, w_prior_source="ledger", deconv=False)` → `families.collapse_by_family` (the **TRUE
RPM pool** `log2(1+Σ(2^x−1))`, **not** a `groupby.mean`) → `attribution_eb._prep` (yr = −resid(Y|C), Xz
z-scored, variance-floored) → `spike_slab._gibbs_posterior(pi=ones, 2000/700)` ⇒ **dense β = coupling**.
Field = `Xz @ β` (**on the scale β was fitted on**, never z-scored β into raw doses). Core C = CPE +
`target_cn` + `mal_prolif`.

**OOF LEAKAGE CONTROL.** 5-fold `KFold(shuffle, random_state=0)` over patients, identical folds for every
gene/rung/seed/decoy so every comparison is paired on the same split. Per fold, **fitted on TRAIN only**: the
C coefficients, the z-scoring mean/sd, the sd<0.1 variance floor, and β. Held-out Y is residualised with the
TRAIN C-coefs. Scoring = Spearman(field, resid_test), averaged over folds. **CPTAC = a true transfer** (β
frozen from TCGA, standardised design rebuilt *within* CPTAC).

---

## §2 ⭐ THE FOUR-RUNG TABLE (THE DELIVERABLE)

All numbers below re-derived by me from the persisted per-gene rows
(`output/learned/mh127/{ladder_oof.tsv.gz, p2_transfer_results.tsv.gz}` → `mh127_tcga_rungs.tsv`,
`mh127_cptac_rungs.tsv`). FAKE = mean over its 5 seeds, per gene. Negative ρ = repressive = better.

### 2a. TCGA — 5-fold OOF (the in-cohort ceiling; n ≈ 1,041 patients)

| decoy | gene set | **R1** real·unfitted | **R2** fake·unfitted | **R3** real·FITTED | **R4** fake·FITTED | fake **reach** (R4/R3) |
|---|---|---|---|---|---|---|
| FAKE1 | all (899) | −0.0268 | −0.0157 | **−0.0465** | **−0.0367** | **78.9%** |
| FAKE1 | n_fam ≥ 2 (325) | −0.0368 | −0.0152 | **−0.0917** | **−0.0725** | **79.1%** |
| FAKE2 | all (834) | −0.0225 | −0.0161 | **−0.0396** | **−0.0336** | **84.9%** |
| FAKE2 | n_fam ≥ 2 (278) | −0.0306 | −0.0221 | **−0.0822** | **−0.0750** | **91.3%** |

**The paired tests (Wilcoxon; ties excluded from win-rates; bootstrap CI on the paired Δ):**

| test | what it asks | FAKE1 (all, 899) | FAKE1 (n_fam≥2, 325) | FAKE2 (all, 834) | FAKE2 (n_fam≥2, 278) |
|---|---|---|---|---|---|
| **R3 vs R4** | ⭐ **is the learning more than fitting?** | Δ **−0.0098** [−0.0213,+0.0018], 52.5% win, **p = 0.211** (MDE 0.0161) | Δ −0.0192 [−0.0379,−0.0010], 55.1%, **p = 0.105** (MDE 0.0266) | Δ −0.0060, 51.0%, **p = 0.577** | Δ −0.0072, 51.8%, **p = 0.587** |
| **R1 vs R2** | curation ALONE (unfitted) | Δ −0.0110, **p = 0.112** | Δ −0.0216, **p = 0.050** | Δ −0.0065, p = 0.405 | Δ −0.0085, p = 0.372 |
| **R3 vs R1** | the **fitting gain** on REAL | Δ −0.0198, **p = 4.9e-18** (571 exact ties) | Δ −0.0549, p = 1.5e-18 | Δ −0.0171, p = 9.7e-13 | Δ −0.0516, p = 3.3e-13 |
| **R4 vs R2** | ⭐ the **fitting gain on the DECOY** | Δ **−0.0210**, **p = 3.2e-30** (72.8% of non-tied) | Δ −0.0573, p = 1.4e-30 | Δ −0.0176, p = 6.9e-21 | Δ −0.0529, p = 1.4e-21 |
| **R4 vs R1** | ⭐ **fitted FAKE vs real ABUNDANCE** | Δ −0.0100, **p = 0.054** | Δ **−0.0357**, **p = 3.6e-04** | Δ −0.0111, **p = 0.038** | Δ **−0.0444**, **p = 5.4e-05** |

**Read the last two rows together.** The decoy gains **as much from fitting (−0.0210)** as the curated set
does (**−0.0198**) — **fitting is a GENERIC gain, not a curation gain** — and the fitted decoy then **beats
the real curated abundance sum**, decisively on multi-regulator genes.

**A structural fact that governs the aggregate** (predicted before the run, then confirmed): for a gene whose
trimmed design has **one family**, `Xz@β` is a monotone rescaling of the abundance pool ⇒ **R3 ≡ R1 exactly**
(measured: 5,620/5,650 one-family units are exact ties, 99.5%). 574/899 (64%) of FAKE1 genes are one-family
after trimming and contribute **zero** fitting signal. The aggregate near-tie is therefore **not** a dilution
by "diffuse" genes — it is 64% of genes where no fitting happens at all, plus a real-vs-fake gap that never
separates in the rest.

### 2b. CPTAC — the TRANSFER (frozen TCGA β, never refit; n = 101; both C blocks)

Per-gene Spearman on CPTAC samples; one-sample Wilcoxon over genes (is the rung predictive **at all**?), and
the paired **R3 vs R4** decisive test. "core" = CPE + target_cn + mal_prolif; "comp" = + the Wu-major
composition block (axiom 2a: report **both**, never silently condition).

| decoy | layer | C | **R1** real·unf | **R2** fake·unf | **R3** real·FIT | **R4** fake·FIT | **R3 vs R4** (paired) |
|---|---|---|---|---|---|---|---|
| FAKE1 (685 g) | mRNA | core | −0.0052 (p=0.21) | +0.0044 (p=0.37) | ⭐ **−0.0189 (p = 8.9e-04)** | −0.0056 (p = 0.19) | **Δ −0.0132, p = 0.030** (MDE 0.020) |
| FAKE1 | mRNA | comp | +0.0016 (p=0.92) | +0.0072 (p=0.14) | −0.0052 (p = 0.18) | +0.0019 (p = 0.68) | Δ −0.0071, p = 0.130 |
| FAKE1 | **protein** | core | +0.0042 (p=0.42) | +0.0036 (p=0.74) | −0.0063 (**p = 0.38, NULL**) | −0.0042 (p = 0.30) | Δ −0.0021, **p = 0.86** (MDE 0.020) |
| FAKE1 | **protein** | comp | +0.0071 (p=0.14) | +0.0029 (p=0.59) | +0.0019 (**p = 0.75, NULL**) | −0.0005 (p = 0.93) | Δ +0.0023, **p = 0.86** (MDE 0.018) |
| FAKE2 (632 g) | mRNA | core | −0.0073 (p=0.12) | +0.0042 (p=0.43) | ⭐ **−0.0193 (p = 5.9e-04)** | −0.0057 (p = 0.26) | **Δ −0.0137, p = 0.027** |
| FAKE2 | mRNA | comp | −0.0026 (p=0.42) | +0.0080 (p=0.089) | **−0.0092 (p = 0.043)** | +0.0031 (p = 0.48) | ⭐ **Δ −0.0126, p = 0.018** |
| FAKE2 | **protein** | core | +0.0053 (p=0.31) | +0.0035 (p=0.50) | −0.0033 (**p = 0.71, NULL**) | −0.0032 (p = 0.62) | Δ −0.0003, **p = 0.96** |
| FAKE2 | **protein** | comp | +0.0073 (p=0.14) | +0.0043 (p=0.26) | +0.0034 (**p = 0.55, NULL**) | +0.0010 (p = 0.70) | Δ +0.0022, **p = 0.95** |

**Three things to take from this table.**
1. **mRNA transfer is where curation earns its keep.** The real fitted field is the **only** rung that
   transports (p = 8.9e-04 / 5.9e-04); real unfitted abundance does not (p = 0.21), the fitted decoy does not
   (p = 0.19), and the paired real-vs-fake test is significant under **both** decoys (p = 0.030 / 0.027) —
   and under FAKE2 it **survives composition** (p = 0.018).
2. **The protein arm cannot arbitrate.** Under this from-scratch, correctly-pooled field the **real** model is
   itself null on CPTAC protein (p = 0.38 core, 0.75 comp) ⇒ you cannot compare two models on a readout where
   the reference model is null. It is **also** underpowered for a ceiling-sized gap (MDE 0.018–0.020 vs an
   in-cohort gap of ~0.017). **Do not report this as "real = fake on protein".** See §5 and §7·4.
3. The composition block halves-to-kills the mRNA transfer for FAKE1 (−0.0189 → −0.0052) but leaves a
   significant real-vs-fake gap under FAKE2. Composition remains a live confound (MH-107/MH-114); it is not
   the whole story here.

---

## §3 THE HEAD-TO-HEAD PARTIAL (max power; does each field carry signal the other lacks?)

Partial Spearman on held-out samples, TCGA OOF, 503 genes (n_fam ≥ 2):

| partial | ρ | p | genes negative |
|---|---|---|---|
| **REAL·fit \| FAKE·fit** | **−0.0564** | **4.0e-16** | 66.6% |
| **FAKE·fit \| REAL·fit** | **−0.0378** | **1.5e-12** | 63.0% |
| FAKE·fit \| real abundance | −0.0431 | 1.6e-14 | — |

**Both are true and they answer different questions — quoting only one would misrepresent the arc.**
*Marginal*: cannot reject real = fake (reach 79–91%, p = 0.11). *Partial*: the real field carries signal the
decoy does not (p = 4e-16) — **and the decoy carries an independent chunk of comparable size (67% of it)**.
The honest reading: **curation contributes something real; so does a matched impossible design.**

In CPTAC the partial is **asymmetric**, which is the cleaner statement: REAL·fit | FAKE·fit = **−0.0198,
p = 3.6e-04**; FAKE·fit | REAL·fit = −0.0041, **p = 0.29**.

---

## §4 ⭐ THE CORRECTED FRAMING OF MH-115 — IT MUST BE RESTATED

MH-115 concluded: *"the LEARNED β is the only estimator that predicts CPTAC protein, and abundance is
REDUNDANT given β."* The comparison it ran was **fitted β vs unfitted abundance**. This arc shows that
comparison **rewards fitting, not curation**:

- **A fitted DECOY beats real unfitted abundance** (TCGA OOF: Δ = −0.0100, p = 0.054 all genes;
  **Δ = −0.0357, p = 3.6e-04** on multi-regulator genes; FAKE2: p = 0.038 / 5.4e-05).
- **The decoy gains from fitting exactly as much as the real set does** (−0.0210 vs −0.0198).
- **MH-115's own "β is additional to abundance" partial survives** under my field (REAL·fit | real abundance =
  −0.0441, p = 5.2e-11; abundance | β = **+0.0187**, i.e. no independent repressive signal) — **but the same
  test against the FAKE fitted model dies on protein** (p = 0.93) and only survives on TCGA (p = 4e-16) and
  CPTAC mRNA (p = 3.6e-04).

**⇒ MH-115's headline as stated is WRONG, and the registry row is annotated accordingly.** Proposed corrected
wording (now carried in `DISCOVERY_REGISTRY.md`):

> *"The learned β out-predicts the evidence-weighted pressure heuristic and raw abundance — **but so does a
> matched, impossible, abundance-matched decoy regulator set once it is FITTED** (MH-127). The advantage over
> abundance is therefore a **FITTING** effect, not a curation effect. What survives the decoy control is
> narrower and more interesting: **the curated design TRANSPORTS out-of-cohort and the decoy does not**
> (CPTAC mRNA: real −0.019 p=8.9e-04 vs fake −0.006 p=0.19; paired p=0.030/0.027)."*

**⚠ ONE DISCREPANCY TO RECONCILE, NOT A RETRACTION.** MH-115's CPTAC **protein** claim was produced by
`eval/coupling_grid` (which has since been fixed under MH-125 and re-run: its persisted
`coupling_grid_summary.tsv`, 2026-07-13 18:22, still reports learned·protein·core mean ρ = −0.035,
p = 7.8e-12 over 1,162 genes). **My from-scratch field does not reproduce a protein signal (p = 0.38 / 0.75)
on the 685 decoy-testable genes.** The two constructions differ in (a) gene universe — mine is the
decoy-testable **trimmed** subset, which drops the most abundant arms (§7·1) — and (b) the design the β
multiply. **Someone must diff the two field constructions.** [INFERRED, not measured: the trim is the likely
cause.] Until then, the protein leg of MH-115 is **flagged, not withdrawn**.

---

## §5 THE STRATIFIED ANSWER — the rescue FAILS (and I will not soften it)

**THE PREDICTION** (task, and mine): *for SINGLE_DRIVER / few-regulator genes REAL should beat FAKE
decisively; for DIFFUSE / many-regulator genes FAKE should match REAL.* **MEASURED: the opposite ordering, and
no stratum is significant.**

Shapley-concentration classes were built from the canonical `attribution.shapley_identity` (LMG on R², the
collinearity-fair readout), applied only where non-degenerate (n_fam ≥ 3; **top-1 share ≥ 0.5 by construction
at p = 2**, so 100% of 2-family genes are vacuously SINGLE_DRIVER, and 48% of the universe has ONE family and
no identity at all).

**R3 vs R4 by stratum (FAKE1, TCGA OOF, paired Wilcoxon):**

| stratum | n | Δ (real − fake) | p |
|---|---|---|---|
| regulator count 1–3 | 710 | −0.0069 | 0.599 |
| 4–8 | 155 | −0.0195 | 0.173 |
| 9–15 | 29 | −0.0209 | 0.508 |
| 16–25 | 5 | — | too few |
| **26+ (incl. PTEN, 82 curated regulators)** | **0** | — | **NOT DECOY-CONTROLLABLE AT ALL** |
| SINGLE_DRIVER | 440 | −0.0145 | 0.273 |
| CO_DRIVERS | 36 | −0.0362 | 0.176 |
| DIFFUSE | 27 | −0.0320 | 0.248 |

**The gap is SMALLEST where credit is concentrated and LARGEST where it is diffuse — the reverse of the
prediction — and nothing is significant** (FAKE2 agrees: SINGLE_DRIVER p = 0.85, CO_DRIVERS p = 0.061,
DIFFUSE p = 0.62). **The "real but scoped" rescue does not hold: there is no subgroup in which the real
fitted model demonstrably beats the matched fitted decoy.**

**⛔ AND SECTION-C-STYLE "reach % by class" ORDERINGS ARE WITHDRAWN.** An adversarial stratum lens measured
that the reach-% ordering across classes is **indistinguishable from randomly shuffled class labels
(permutation p = 0.92)**, that **no** stratification (real-fit-derived, fit-free, DoF-matched, or
FAKE-fit-derived) produces heterogeneity (Kruskal–Wallis n.s.), and that per-stratum p-values run 0.067–0.99.
**The class ordering is NOISE. Do not quote it.** (It also re-fit the strata on a disjoint half and found the
classes are not a leak — the strata are *innocent*, just uninformative.)

**⭐ THE ONE STRATIFIER THAT DOES SEPARATE — and it is FIT-FREE.** Splitting genes into tertiles of **curated
evidence weight `w`** (a property of the literature, not of any fit;
`output/learned/mh127/stratlens_all_tests.tsv`):

| decoy | tertile | R3 vs R4 (fitted) | R1 vs R2 (**unfitted**) |
|---|---|---|---|
| FAKE1 | w-lo (114) | +0.0001, p = 0.71 | −0.0114, p = 0.70 |
| FAKE1 | w-mid (114) | −0.0064, p = 0.85 | −0.0023, p = 0.92 |
| **FAKE1** | **w-hi (115)** | **−0.0580, p = 3.0e-04, BH q = 0.0059** | **−0.0575, p = 3.3e-04, q = 0.0059** |
| FAKE2 | w-hi (101) | −0.0419, p = 0.022, **q = 0.20** | −0.0290, p = 0.067, q = 0.34 |

**This is the single most pro-curation number in the arc, and note WHERE it lives: in the UNFITTED rung as
much as the fitted one (−0.0575 vs −0.0580).** In the high-evidence tertile the *edge set itself* — not the
learning — carries the signal. **⚠ But it does not replicate under the structure-matched decoy FAKE2
(q = 0.20), and the arc's own rule is "report under both decoys; if the gap differs, it is not robust."
⇒ SUGGESTIVE LEAD, NOT AN ESTABLISHED RESULT.**

---

## §6 WHAT THE LEARNING ACTUALLY BUYS — the adjudicated list

MH-126b claimed three things for the learned model. The decoy control adjudicates them:

| claim | verdict under the decoy control | evidence |
|---|---|---|
| **De-confounding** (β is estimated *given* C; a raw correlation is not) | ⚠ **TRUE BUT NOT CURATION-SPECIFIC.** The decoy is fitted given the same C and gets the same benefit. De-confounding is a property of the ESTIMATOR, and it is *why* a fitted decoy beats real abundance — it cannot be used as evidence that the curated edges are real. | R4 vs R1: p = 0.054 (all) / 3.6e-04 (multi-reg) |
| **Joint apportionment under collinearity** | ⛔ **NOT CURATION-SPECIFIC.** The decoy gains from fitting as much as the real set (−0.0210 vs −0.0198), and the structure-matched decoy FAKE2 (identical family-size profile) reaches **91.3%** of the real model on multi-regulator genes. | R4 vs R2 p = 3.2e-30; reach table §2a |
| **Out-of-sample / out-of-layer TRANSFER** | ✅ **SURVIVES — this is the one that is curation-specific.** Frozen TCGA β predict **CPTAC mRNA** (real p = 8.9e-04; matched decoy p = 0.19; paired p = 0.030 / 0.027, and p = 0.018 composition-adjusted under FAKE2). ⛔ **DIES on CPTAC protein — for BOTH arms** (the real model is itself null there in this construction, and the test is underpowered). | §2b, §3 |
| **In-cohort gene-budget prediction (the MH-115 claim)** | ⛔ **DOES NOT SURVIVE.** A matched impossible design reaches 79–91% of it. | §2a |
| **The curated EDGE SET carries signal a matched impossible set does not** | 🟡 **PARTIAL.** Yes in the CPTAC partial (asymmetric: p = 3.6e-04 vs 0.29) and in the high-evidence tertile (fit-free, q = 0.006) — but the same TCGA partial shows the **decoy also carries independent signal** (p = 1.5e-12), and the marginal in-cohort test cannot separate them. | §3, §5 |

**The one-line summary of what the learning buys:** *de-confounded, jointly-apportioned fitting — which a fake
regulator set gets for free — plus a transferable component that the fake set does not get.* **Edge existence
does not rest on this model** (MH-124's site/abundance-matched enrichment and MH-126's instrument work are the
existence evidence; note MH-124r retracted the `pi_causal` leg).

---

## §7 UNDERPOWERED vs NEGATIVE (kept strictly separate)

1. **R3 vs R4 (the decisive test) is UNDERPOWERED for a small true effect, NOT proven zero.** Bootstrap CI
   **[−0.0213, +0.0018]**, paired **MDE = 0.0161** (80% power, α = 0.05). **I can exclude a real-model
   advantage larger than ~0.021 ρ; I cannot exclude a smaller one.** The point estimate favours real in
   **all 5 seeds and both decoy sets** — consistent direction, indistinguishable magnitude. A leakage lens
   further measured that the **per-gene SE is understated ~1.8×** by the shared-fold structure ⇒ the true
   interval is **wider**, so the data exclude even less. *Honest headline: whatever real-vs-fake advantage
   exists is smaller than the generic gain from fitting (−0.021), which the decoy gets for free.*
2. **CPTAC PROTEIN is UNDERPOWERED AND REFERENCE-NULL — it is not a null of the decoy question.**
   Observed real-vs-fake Δ = +0.0023, **MDE = 0.018** (comp) / 0.020 (core); power to detect the
   ceiling-sized gap (−0.017) ≈ **0.56**. And the reference arm (real) is itself null (p = 0.75).
3. **The 26+ regulator bin (12 genes, incl. PTEN) is NOT decoy-controllable at all** (0–3 eligible fake arms:
   hub genes have target-rich 3′UTRs, so essentially no arm is site-free on them). The graded analysis **tops
   out at 9–15 regulators (n = 29)**; 16–25 (n = 5) is too thin. This is a **structural limit of the control**,
   not a null.
4. **No multiple-testing correction was applied across the ~30 stratum tests in §5.** If it were, **nothing**
   would survive except the fit-free w-hi tertile (which carries its own BH q = 0.0059 within its 36-test
   family, and fails under FAKE2).

---

## §8 ADVERSARIAL VERIFICATION LOG (six lenses; refutations recorded)

| lens | verdict on the attacked result | what it measured |
|---|---|---|
| **Leakage ×2** (independent re-implementations) | **SURVIVES** | Both reproduced the ladder from scratch (one bit-exact, max\|Δ\| = 0.0e+00 on all four rungs). A **200-permutation global patient-relabel null run through the pipeline centres both arms at ≈0** (real arm mean ρ = −0.00001, SD 0.0036) ⇒ fitting gains **nothing** under a null Y. A permuted-Y probe returns a **POSITIVE** (conservative) ρ (+0.011) and is **arm-symmetric** (+0.01127 real vs +0.01121 fake). The maximal conceivable β/C/z leak is worth ~0.007 ρ and **helps the decoy as much as the real model**. |
| **Match-quality ×2** | **SURVIVES** | The decoys are **not collinear proxies** (subspace R² median **0.028**; 0% of designs > 0.5). Abundance matched (MWU p = 0.96–0.97); **detection rate matched** (p = 0.30–0.59); zero columns hit the variance floor. The axes that are **mismatched all HANDICAP THE FAKE**: SD (real 1.03 vs 0.98, p = 3.5e-17), HE degree (median 20 vs 3 targets), TargetScan promiscuity (504 vs 0 genes), \|r\| with miRNA PC1. A **size-matched** decoy at the deployed design still reaches **76.1%** and beats real abundance at **p = 7e-9**. Also measured the unmeasured caveat: **the trim costs the REAL model 1.8×** while the decoy's reach stays ≈0.79. |
| **Stratum lens (a)** | **⛔ REFUTED the ladder's "n.s. in EVERY stratum"** | Built a **fit-free** stratifier (curated evidence weight) and found real > fake at **BH q = 0.0059** in the top tertile, **fitted AND unfitted** (§5). Its own DoF-artifact and collinearity hypotheses were **refuted by its own measurements** (FAKE2 is 100% DoF-matched; the DoF gradient runs the opposite way). It **confirms** the aggregate near-tie (its full-universe R3vR4 = −0.0119, p = 0.146) and the "fitting is a generic gain" core. |
| **Stratum lens (b)** | **CONFIRMS the headline; ⛔ KILLS the class ordering** | The reach-% ordering across Shapley classes is **indistinguishable from shuffled labels (perm p = 0.92)**; no stratifier produces heterogeneity; the classes are **not** a leak (disjoint-half refit agrees) — they are simply **uninformative**. |

---

## §9 CAVEATS — what this test is and is not

1. **THE TRIMMED ESTIMAND.** This asks *"does the curated design beat a matched impossible design **on the
   subset of regulators for which a matched control exists**"* — not on the full curated set. Arm retention is
   77% in the 1–3 bin but 38–48% for 4–15-regulator genes, and **the dropped arms are the most abundant ones**
   (32% of curated arms are unmatchable at any caliper; median 7.80 log2RPM). ⇒ the trimmed test may
   **understate** the real model (measured: the trim costs the real model 1.8×). The full-set test is
   **abundance-biased toward real** and was deliberately not used as a headline.
2. **IRREDUCIBLE AMBIGUITY OF THE WHOLE DESIGN** (the deepest thing the match-quality lens found, and it cuts
   both ways): **abundance is simultaneously the NUISANCE (shared global axes) and the CAUSAL DOSE.** Matching
   on it matches away part of the causal channel; not matching it confounds. For the most abundant curated
   arms **no valid decoy exists at all**.
3. **The 5 seeds are NOT 5 independent draws** (mean Jaccard 0.69 / 0.76; median eligible pool 40 arms). Seed
   spread is a **stability** measure, not a confidence interval.
4. **The Shapley classes are derived FROM the real fitted model** — legitimate as a stratifier (both arms are
   scored inside the same stratum), **not** evidence that the classes are biology. Per MH-124, β names the
   literature regulator **at chance**, so "SINGLE_DRIVER" means *the model concentrates its credit*, not
   *one miRNA drives this gene*.
5. **Effect sizes are TINY everywhere** (|ρ| 0.02–0.09 at the ceiling; median fitted R² ≈ 0.028). This is a
   comparison of two weak predictors.
6. **Fit under CORE C.** β were fitted under core C in TCGA (the canonical coupling block); the composition
   block was applied at **scoring** time in CPTAC. A fully composition-**fitted** decoy ladder was not run.
7. **NOT MEASURED:** the full (untrimmed) ladder; a per-edge (rather than per-gene) decoy score under this
   design; any CN-instrument cross-check of the genes where real does beat fake; a "no *meaningful* site"
   relaxation of the site-free rule (which would enlarge the eligible pool and could rescue hub genes — at the
   cost of weakening the control).

---

## §10 REPRODUCE INDEX

All artifacts persisted to **`mirna_hallmark/output/learned/mh127/`**:

| file | what |
|---|---|
| `fake_regulator_sets.tsv.gz` | 29,920 rows — (gene, seed, fake_set, real_arm → fake_arm, abundances, `in_caliper`). **⭐ Filter `in_caliper=True` to get the deployed TRIMMED design.** |
| `fake_set_match_quality.tsv.gz`, `trimmed_design_quality.tsv.gz` | per-(gene,seed,set) match QC (abundance, n_families, within-design mean \|r\|, max \|r\| to any real regulator) |
| `gene_strata.tsv.gz` | 1,421 genes × 28 — dense-Gibbs β summary, `shapley_identity` vector, top1/entropy/gini, Shapley class, regulator-count bins, decoy-testability flags |
| `ladder_oof.tsv.gz` | per (gene, seed, fake_set) TCGA 5-fold OOF ρ for all four rungs |
| `p2_transfer_results.tsv.gz` | per (gene, seed, fake_set, REAL/FAKE) CPTAC transfer ρ — mRNA + protein × {raw, core, comp} |
| `mh127_tcga_rungs.tsv`, `mh127_cptac_rungs.tsv` | **the §2 tables**, re-derived from the two files above |
| `stratlens_all_tests.tsv` | the §5 stratified tests incl. the fit-free evidence-weight tertiles |
| `code/` | `universe.py` · `build_fakes.py` · `build_strata.py` · `ladder_oof.py` · `p2_transfer.py` · `stratlens_final.py` · `permnull.py` · `atk_leak.py` · `mh127_synth_{tcga,cptac}.py` |

No repo module was changed by this arc.
