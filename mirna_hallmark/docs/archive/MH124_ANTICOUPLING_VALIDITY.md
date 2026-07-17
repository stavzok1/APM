# MH-124 — ARE ANTI-COUPLINGS MEANINGFUL? Edge existence, true importance, and the arm-vs-family weight question

> **Goal:** settle, with orthogonal evidence and no spurious conclusions, whether miRNA→target anti-couplings
> in TCGA bulk carry real repression signal — for HE (curated) arms vs abundance-matched site-free decoys —
> and whether the LEARNED WEIGHTS (per-arm and per-family) recover what biology/literature says they should.
> **What belongs here:** the design of each test, the CONTROL it carries, the measured result, and the verdict.
> Negative and refuted results are first-class. Nothing enters without a control.
> **What does NOT belong here:** anything not measured in this task; anything from the retracted MH-118..121 arc.
> **Update trigger:** any new measurement on these questions.
> **Sync-partner:** `ANALYSIS_RUN_LEDGER.md` (MH-124 row), `DISCOVERY_REGISTRY.md` (MH-124 claim).


> # ⛔⛔ CRITICAL RETRACTION (2026-07-13, MH-126/128) — **§5 AND §5b ARE VOID AS STATED.**
> **`instrument.exclusion`'s `pi_causal` IS NOT EXOGENOUS.** Read the code (`instrument.py:1118-1136`):
> ```python
> a        = (zs @ xr) / zz          # first stage  CN_seg -> X_fam        EXOGENOUS ✓
> b, delta = lstsq([xr, zs], yr)     # b = X_fam -> Y | Z_s                ⚠ OBSERVATIONAL OLS SLOPE
> pi_causal = a * b                  # <- PRODUCT-OF-COEFFICIENTS MEDIATION estimator, NOT an IV
> ```
> `b` is the **observational regression of the target on the miRNA** and carries every confound. We measured
> that HE edges anti-correlate 5× more than decoys (−0.050 vs −0.011) ⇒ **`pi_causal = a·b` separates HE from
> decoy FOR THAT REASON ALONE, with zero causal content.** A verifier reproduced the entire "3× validation"
> from a **zero-causal-effect simulation**.
>
> **THE CLEAN TEST (`scratchpad/clean_reduced_form.py`; the REDUCED FORM `corr(resid(CN|C), resid(Y|C))` is the
> only quantity that touches Y solely through copy number; F-gated on the FIRST STAGE only, since the exclusion
> gate keys on δ_s which is itself derived from the observational regression):**
>
> | C block | DECOY reduced form (must be ≈0 for validity) | HE vs DECOY |
> |---|---|---|
> | C (as used in §5) | **−0.105, 86.7% neg — NULL BROKEN, TEST INVALID** | −0.111 vs −0.106, p=0.32 |
> | **C + ploidy** (valid) | **−0.0021, 49.9% neg** ✅ | **−0.0096 vs −0.0025, p=0.030** |
>
> **…but under the CORRECT UNIT OF ANALYSIS (1,793 edges are only ~330 ARMS): clustered by arm, p = 0.115 —
> NOT SIGNIFICANT.** Seed-stability: p<0.05 in only **2/10** matching draws (median p=0.068).
> For contrast, the **observational** coupling on the SAME edges separates at **p=5.4e-16**, and the first stage
> is IDENTICAL for real and fake (p=0.76) — exactly as it must be.
>
> ⇒ **THE EXOGENOUS INSTRUMENT DOES NOT DISTINGUISH REAL EDGES FROM IMPOSSIBLE ONES.**
> ⇒ **THE "THREE INDEPENDENT VALIDATIONS" HEADLINE IS FALSE.** The other two legs (set-level OR 1.41; the
> decoy-controlled PIP/prediction test) are **BOTH OBSERVATIONAL** and both rest on the same anti-correlation —
> they were never independent of each other. **Edge existence now rests on ONE line of evidence: curated edges
> anti-correlate more than abundance-matched site-free ones (p=1.1e-11).** Real, but one line, not three.
>
> ⚠ **UNDERPOWERED, NOT REFUTED.** The reduced-form point estimate points the right way and is *quantitatively
> consistent* with a real effect: first stage 0.184 × an implied causal effect −0.052 = −0.0096, exactly what is
> observed. With only ~157 HE arms carrying a usable instrument, the reduced form is attenuated below what this
> cohort can resolve. **"The data cannot say" — NOT "the answer is no."**
>
> **DO NOT CITE §5/§5b. DO NOT CALL `pi_causal` EXOGENOUS.**


> ### ⛔ SECOND RETRACTION (2026-07-16, MH-129/133) — **THE WITHIN-ARM REPLACEMENT ALSO FAILS.**
> After `pi_causal` was retracted (banner above), a **within-arm two-way fixed-effects** design was built as the
> honest replacement: `reduced(m,g) = α_m + γ_g + τ·site(m,g)`, where α_m absorbs ALL locus pleiotropy (a FE, so
> no observational slope is ever touched) and τ is identified by **site presence = sequence = exogenous to CN**.
> It gave τ = −0.0039 (p=0.029, 2-way clustered), a placebo at exactly 0, and a Wald of −0.019 [−0.034, −0.003].
> **IT DOES NOT SURVIVE ITS OWN FALSIFICATION PROGRAMME:**
> * **F1 — the control class was 71% REAL BINDERS.** My "non-site" class was merely *not-curated*; the
>   genome-wide scanMiR scan (10.4M pairs) shows **71.1%** of the panel carries a duplex. Against a **genuinely**
>   site-free control (no TargetScan site **AND** no scanMiR duplex): **τ = −0.0007, p = 0.84** on the very
>   F>10 subset that produced the headline. **Nothing.**
> * **F2 — the SITE-TYPE EFFICACY LADDER FAILS.** Copy number cannot see a site type, so genuine
>   site-mediated repression MUST follow 8mer > 7mer-m8 > 7mer-A1 (Bartel). Measured vs the genuine site-free
>   control: **−0.0017 / −0.0001 / −0.0006** — none significant, **NOT monotone**, ladder contrast p=0.26.
> * **F3 — guide-vs-passenger confirms in only 1 of 3 specifications** (unweighted p=0.010; F-weighted p=0.17;
>   confident-guide p=0.29). The passenger estimate swings −0.0025 → −0.0034 → −0.0167 as n falls 1,011 → 356 →
>   109 — the signature of noise.
> ⇒ **BOTH CN INSTRUMENTS ARE RETRACTED. EDGE EXISTENCE HAS NO SURVIVING EXOGENOUS VALIDATION.**
> It rests on the OBSERVATIONAL contrast alone (**MH-131: gap −0.0215, p=6.5e-04, an upper bound**).
> ⚠ **UNDERPOWERED vs REFUTED:** F1/F2 are genuine failures against a proper control, not power problems
> (n = 216k–235k pairs). F3 is underpowered. Do not cite MH-129's τ, Wald, or "predicted sites beat curation".


**STATUS: COMPLETE — but see BOTH RETRACTIONS above.** SETTLED: Q1, Q1b (§2.4) · Q2 (§5) · Q3a (§5b.1) · Q4, Q5 (§4) · Q6 (E8) · **Q7 (§6)**. UNDERPOWERED (not negative): **Q3b (§5b.2)** — the literature set (n=16), not the estimator, is the constraint. **Plus a correctness BUG in the canonical Gibbs found, fixed and verified (§6b).**

> ## ⛔ AMENDED BY MH-126 (2026-07-13) — READ BEFORE USING §5
> **The third existence leg — "the EXOGENOUS CN instrument separates real from decoy" (`pi_causal` −0.0235 vs
> −0.0075) — IS NOT EXOGENOUS AND IS DOWNGRADED.** `instrument.exclusion` returns `pi_causal = γ_s · b_fam`,
> where `b_fam` is the **OLS partial slope of Y on X_fam given [C, Z_s]** — the instrument-**orthogonal**
> (endogenous) variation; it reproduces the plain no-copy-number partial correlation at ρ≈0.93, and the whole
> "3× separation" was **reproduced from a simulation with ZERO causal effect**. The only exogenous factor,
> **`γ_s`, is SITE-BLIND (HE +0.1993 vs DECOY +0.1984, ratio 1.004×, p=0.20)**, and the genuinely exogenous
> **reduced form** shows **no HE-vs-decoy effect** at the **locus** unit (τ=−0.0015, p=0.67).
> **§5's other finding — that the RAW instrument fails because of aneuploidy — is CONFIRMED and EXPLAINED**
> (a ploidy sign-flip, not co-amplification; the fix takes the impossible-edge negative-reduced-form rate from
> 86.4% to 50.0%, the null). ⇒ **EDGE EXISTENCE now rests on the two OBSERVATIONAL legs (§2.4 enrichment,
> §5's decoy-controlled model test) and has NO surviving exogenous validation.**
> **Q3b (attribution) is UNCHANGED and still UNDERPOWERED** (aneuploidy-controlled: 4/16, p=0.18).
> See **`MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md`** §1 and `CN_INSTRUMENT.md` §9x.
> MH-126 also answers the **graded** version of Q3b: β does **not** rank the canonical family above chance
> either (0.518, p=0.66) — but it *does* carry a weak evidence gradient (+0.068, permutation p=0.005).

### ⭐ THE HEADLINE SO FAR
**⛔ RETRACTED: "EDGE EXISTENCE is validated THREE independent ways"** — the CN leg was NOT exogenous and the other two are both observational. **ONE line of evidence stands.** Original text: — (1) set-level enrichment vs abundance-matched
site-free (**OR 1.41, p=2.5e-13**, survives composition); (2) decoy-controlled model test (**PIP 23.0% vs
15.5%, p=1.1e-04**; prediction **0.146 vs 0.104, p=2.1e-07**); (3) the **exogenous CN instrument**,
exclusion-gated and abundance-matched (**π_causal −0.0235 vs −0.0075, p<1e-4**).
**IMPORTANCE / ATTRIBUTION — ⚠ CORRECTED (§4b, user-caught).** The original "at chance" verdict used **β** and
**argmax**, both of which the doctrine forbids for this question. With the doctrine's own estimator, the model
**DOES rank the canonical family above chance** (Shapley identity rank 0.319, **p=0.0076**; DOSE 0.258,
**p=0.0003**; the two are INDEPENDENT) — while **β is at chance (p=0.48)**, because the canonical regulator is
**high-DOSE, not high-variance**. ⇒ **the model RANKS the right family high but does not NAME it** (argmax n.s.).
Separately: decoys still absorb **36%** of the total |β| and reach **71%** of the predictive field.
⇒ **THE EDGES ARE REAL. THE WEIGHTS ARE NOT TRUSTWORTHY AS ATTRIBUTION. Existence ≠ importance, and the model
currently conflates them.**

**⭐ AND THE TASK HAS ONE UNIFYING RESULT (§6):** the anti-coupling signal and its background **occupy the SAME
subspace**. A rank decomposition (`X = S + L`) does NOT separate them — removing the low-rank part kills the
REAL signal *faster* than the decoys, monotonically, until separation → 0. The top miRNA factors ARE
proliferation (r=0.33) and purity (r=0.41) — but that is also **where the biology acts**. ⇒ **No endogenous
re-basis of the miRNA×sample matrix can identify these edges.** Identification must come from **outside** it —
**an asymmetry the confound is blind to** (axiom 4). **Locus copy number is exactly that**, and it is precisely
the lens that works (§5, §5b).

This file is the single source of truth for the task; it is written to survive a
context compaction. Each section states DESIGN → CONTROL → RESULT → VERDICT.

---

## 0. The questions (user-set, 2026-07-13)

1. **Are anti-couplings meaningful for HE arms vs no-site pairs — or is it just ABUNDANCE?** Test with
   **abundance-matched** decoys, decoys chosen carefully (no seed-family mate, no genomic-cluster mate,
   near-uncorrelated with every real regulator).
2. **Orthogonally validate EDGE EXISTENCE and TRUE IMPORTANCE via the CN instrument** — applied carefully,
   with the cluster/co-location distinctions the infra already encodes.
3. **The same for FAMILY→target edges** — the CN instrument aggregates several loci, by the arms' loci.
   **Infra exists (`instrument.family_multi_iv`, `cn_copy`, `exclusion`, `between_family_exclusion`). DO NOT
   REINVENT.**
4. **Do the per-ARM weights WITHIN a family match literature/expectation?**
5. **Is using n_arms predictors (instead of n_families) itself the problem?** Does the FAMILY-weighted model
   at lower resolution give the coupling weight to the **right literature-wise families**?
6. **Attribution within family to arms** — by the mechanisms this project devised (dose / chimeric / CN).
7. **P.S. Does the arm-resolution model need a SPARSE + LOW-RANK (S+L) decomposition?**

---

## 1. WHAT IS ALREADY ESTABLISHED (carried in, with its control)

| # | finding | control that established it |
|---|---|---|
| E1 | **The per-edge theoretical t-null is wrong by 3–4× in scale.** Site-free pairs (no TargetScan site, no scanMiR duplex ⇒ cannot repress) are "FDR-significant" at **35.3% (core C) / 25.0% (deconv)** against a nominal 5%. Fitted null: **σ₀ = 0.083–0.132 vs theoretical 0.031**; **μ₀ NEGATIVE**. | 60,000 site-free pairs, Efron empirical null (fit location+scale; ⚠ exceedance-counting is resolution-limited, min p = 1/N) |
| E2 | **Curated HE edges DO beat site-free pairs** — abundance-matched (log2RPM>4): **OR 1.41, Fisher p=2.5e-13**, and it **SURVIVES composition adjustment** (OR 1.43 core → 1.41 deconv, essentially unchanged). ⇒ **the curated signal was never compartment arithmetic.** | abundance-matched site-free class |
| E3 | **Sequence candidacy is mostly compartment arithmetic.** Predicted-site vs site-free, abundance-matched: **core OR 1.36 (p=7e-24) → deconv OR 1.11 (p=8e-04)**. HE vs candidate goes from n.s. (p=0.12) to **OR 1.27 (p=5e-08)**. ⇒ composition adjustment **SEPARATES**: curation survives it, sequence prediction does not. | same, under both C blocks |
| E4 | ⇒ **The HE signal is a DISTRIBUTIONAL SHIFT, not per-edge significance.** A typical HE edge has ρ≈−0.15 against σ₀=0.09 ⇒ z≈−1.4 ⇒ real, but nowhere near surviving FDR across 4,000 tests. Recalibrated: HE **31.3% → 0.0%**. | Efron null |
| E5 | **The orphan universe is 100% expression-selected** (0/6,744 have a pure-sequence route) ⇒ its 84% TCGA hit-rate is **circular**. An in-fold-selection OOF *is* computable and does not collapse (92.5% split-half) — **but it is not an arbiter: site-free FAKE edges pass it at 77.5%, and abundance-matched they are identical (90.7% vs 90.7%).** Honest selection-free orphan rate: **19.2%**. | in-fold selection OOF + site-free fake class + abundance matching |
| E6 | **⚠ SELF-PARTIALLING BUG (live):** `discovery`/`dossier` partial every edge on `C + he_agg`, where `he_agg` is built FROM the HE arms ⇒ an HE arm is partialled on a covariate **containing itself** (+ its seed-family mates; 22.6% of HE arms share a family with another HE arm of the same gene). HE mean \|ρ\| shrinks **−40.1%** (vs −6.9% cand, −5.3% fake); HE rate **0.9%** (with he_agg) → 23.2% (LOO) → **32.6%** (C-only). **Options: (2) leave-one-FAMILY-out, (3) two labelled estimators (marginal + incremental), (5) the joint Gibbs — the principled fix, already built.** | class-matched estimator comparison |
| E7 | **Design §F stands (within-family):** same seed ⇒ **site counts 100% identical** across members; within-family scanMiR K_D ratio **median 1.06×**. ⚠ **BUT scanMiR K_D is a function of nt 1–10, NOT the seed** — same-seed members get *different* K_D tables in **83.9%** of families, so "K_D is constant" is FALSE as stated; the differences are simply **small** (median within-(family,gene) range 0.021 on a mean \|repression\| 0.44 ≈ 5%). **And scanMiR sees NOTHING of AGO loading or 5′-isomiR seed shift** ⇒ it cannot test the premise `w_m = θ_m·a_m` in full. | all 1,203 same-seed KdModel pairs; genome-wide scan |
| E8 | **Within-family ATTRIBUTION is a DOSE question, and every coupling lens fails it.** Against **three independent ground truths**: measured chimeric occupancy (n=207, permutation null) — dose **63.3%** vs arm β **44.9% (chance)**; **literature-known repressing member (n=14)** — abundance **92.9% (p=1e-4)**, chimeric 84.6%, CN dose ~75%, **arm β 42.9% (p=0.38, CHANCE)**. | permutation nulls throughout |
| E9 | **Shapley must be NESTED** (family Shapley × within-family share), not one joint game over all arms: the joint game **UP-credits** multi-arm families by **+4.9 pp (p=1.3e-80)**, driven by family SIZE (a degrees-of-freedom artifact). | joint vs nested comparison |

---

## 2. THE DECOY TEST — is β repression, or transportable background?

### 2.1 Design
For each gene, augment its real (HE) regulator design with the **same number of DECOY arms** and fit ONE
canonical Gibbs posterior (dense π≡1) on the joint design in TCGA.

**⚠ DECOY DEFINITION MATTERS (user-flagged, and my v1 was WRONG).**
* **v1 (INVALID):** excluded only arms already in the design, arms with a scanMiR site on that gene, and
  curated edges. **It did NOT exclude seed-family mates or genomic-cluster mates of a real regulator** —
  which are **co-transcribed and collinear**, so they absorb β by collinearity and "predict" CPTAC by
  **proxying the real regulator**. That alone could produce the whole result.
* **v2 (STRICT, the valid one):** additionally requires the decoy to be **NOT a seed-family mate**, **NOT a
  genomic-cluster mate** (`instrument._genomic_clusters`, 50 kb single-linkage), and **|r| < 0.30 against
  EVERY real regulator** of that gene.
* **v3 (REQUIRED, pending):** additionally **ABUNDANCE-MATCHED** to the real regulators — this is the user's
  core question and the dominant confound.

### 2.2 v1 result (⚠ INVALID — collinear decoys not excluded; recorded for provenance only)
| | REAL | DECOY (v1) |
|---|---|---|
| median β | +0.00400 | +0.00369 |
| \|z\|>2 "identified" | 20.7% | 18.4% |
| decoy share of gene's total \|β\| | — | **43.4%** |
| Spearman(β, CPTAC-protein edge ρ), core | −0.2510 | **−0.2528** |
| " , composition-adjusted | −0.1203 | −0.1044 |

**DO NOT CITE.** Superseded by v2/v3.

### 2.3 v2 (STRICT) — decoys not family/cluster mates, |r|<0.30 vs every real arm
816 genes, 4,728 real / 4,728 strict decoy arms.
* **PREDICTION LANE (the key panel):** REAL-only gene field mean OOF partial ρ **+0.1592** vs DECOY-only
  **+0.1164** (real better in **62%** of genes, Wilcoxon **p=3.67e-13**). ⇒ **decoys reach 73% of the real
  field's predictive power.**
* PIP: with the dense π≡1 prior, PIP is **100% for everything** (degenerate — π≡1 means "always include";
  PIP carries no information there). ⚠ **My first PIP comparison was CONFOUNDED** — real arms got their
  individual evidence priors while decoys got the median. Fixed in v3 with a matched uniform prior.

### 2.4 ⭐ v3 (STRICT **+ ABUNDANCE-MATCHED** + MATCHED PRIOR) — **the answer to Q1b**
482 genes, 1,663 real / 1,663 decoy arms. Decoy differs from its real twin in **EXACTLY ONE THING: no binding
site.** Abundance matched (**real median 7.07 vs decoy 6.98, MWU p=0.80** ✓); not a seed-family mate; not a
genomic-cluster mate; |r| < 0.30 against every real regulator; **matched uniform prior π=0.3 for BOTH**, so any
PIP gap is **the DATA**, not the prior.

| | REAL (HE) | **DECOY (no site)** | test |
|---|---|---|---|
| **median β** | +0.00002 | **+0.00002** | identical |
| mean β | +0.01142 | +0.00757 | MWU p=1.89e-04 |
| **PIP>0.5** (matched prior ⇒ **the DATA**) | **23.0%** | **15.5%** | **MWU p=1.10e-04** |
| \|z\|>2 | 19.4% | 12.7% | |
| **decoy share of the gene's total \|β\|** | — | **36.3%** | |
| **OOF prediction** (real-only vs decoy-only field) | **+0.1460** | **+0.1035** | **71% of real**; 59% of genes; Wilcoxon **p=2.05e-07** |

**VERDICT (Q1 + Q1b): NOT just abundance — but SMALL.** With abundance matched and every collinearity route
closed, **the data DOES separate real regulators from impossible ones** (PIP 23.0% vs 15.5%, p=1.1e-04;
prediction 0.146 vs 0.104, p=2.1e-07). **But decoys reach 71% of the predictive power, absorb 36% of the
repression weight, and have an IDENTICAL median β.**
⇒ **The anti-coupling signal is real and site-driven — ~29% of the predictive signal riding on ~71%
background.** The model detects that curated regulators **as a SET** carry genuine information beyond
abundance; it **cannot** say **WHICH** one matters (§4).
Consistent with E2 (set-level OR 1.41, p=2.5e-13) and E4 (nothing survives per-edge). Three framings, one
conclusion.
Outputs: `output/learned/decoy_v3_abundance_matched.tsv`, `decoy_v3_prediction.tsv`.

---

## 3. ⭐ THE INCLUSION-CRITERION INSIGHT (settles "how did we get true regulators?")

`assemble_gene(he_only=True)` builds the design **from the gene's CURATED HE regulators only**. **Every column
the model can weight is already a literature-validated regulator of that gene.** So METHODS §5's *"recovers
curated drivers (PTEN miR-21/103a/181/182/96)"* is **TRUE BY CONSTRUCTION — there is nothing else in the
design to pick.**

⇒ **The biology enters through the INCLUSION CRITERION (curation), not through the learning.** The data's job
is only to *distribute weight among* an already-correct set — and that distribution fails every ground-truth
test (E8). **This is not a contradiction; it is the resolution.**

**⇒ THE OPEN QUESTION THIS FORCES (§4 below): does the FAMILY model — the identified estimand (Design §F) —
give the coupling weight to the RIGHT literature-wise FAMILY? That has NEVER been tested.**

---

## 4. ⛔ BETWEEN-FAMILY ATTRIBUTION vs LITERATURE — **THE MODEL DOES NOT RECOVER THE CANONICAL FAMILY**

### 4.1 Hand-curated classic MTIs (n=19 genes with the literature family IN the design)
| lens | names the literature family | chance (1/n_fam) | binomial p |
|---|---|---|---|
| **FAMILY β** (Design §F's IDENTIFIED estimand) | **1/19 = 5.3%** | 8.2% | **0.80 — BELOW chance** |
| DOSE (abundance) | 3/19 = 15.8% | 8.2% | 0.20 |

Only **COL1A1 ← miR-29** is recovered. The misses are textbook: **ZEB1 ← miR-200** (β says miR-590-3p) ·
**VEGFA ← miR-200** (β says **miR-718**, an arm detected in 1.8% of patients) · **MYC ← let-7** (β says
miR-184) · **CDK6 ← miR-34** (β says miR-21) · **E2F1 ← miR-17** (β says miR-1205).

### 4.2 Well-powered, objective literature proxy (top **evidence_score** regulator per gene; n=173)
**NON-CIRCULAR:** the dense posterior is called with `np.ones(p)` — it never sees the evidence prior.
| lens | hit | chance | p |
|---|---|---|---|
| ARM β | 25.4% | 20.1% | 0.052 |
| ARM dose | 24.3% | 20.1% | 0.10 |
| **FAMILY β** | 24.3% | 21.9% | **0.25** |
| FAMILY dose | 25.4% | 21.9% | 0.15 |
**All at chance.** (The borderline arm-β 0.052 is uncorrected across 4 tests, on genes with a **median of only
5 regulators** where chance is already 20%.)

### 4.3 DECOY **FAMILIES** (strict: no family/cluster overlap, |r|<0.30 vs every real family dose)
516 genes, 3,558 real / 3,558 decoy families. median β REAL +0.00440 vs **DECOY +0.00401**; |z|>2 **24.0% vs
19.1%**; **decoy share of the gene's total |β| = 36.2%.**

### 4.4 VERDICT (Q5)
**⛔ NO. The FAMILY model — the identified estimand — does NOT give the coupling weight to the
literature-canonical family. It is at or below chance, at BOTH resolutions.** Lowering the resolution from
arm to family does **not** fix the attribution problem. ⚠ n=19 excludes a *strong* effect but cannot prove
"exactly chance"; the n=173 evidence-ranked test agrees.
Outputs: `output/learned/between_family_literature.tsv`, `decoy_family_betas.tsv`, `literature_power_test.tsv`.

Literature ground-truth set (classic MTIs, gene → canonical repressing seed family):
ZEB1/ZEB2/CDH1/VEGFA ← miR-200 · CCND1/BCL2 ← miR-15/16 · CDKN1A ← miR-106 · E2F1/BCL2L11 ← miR-17 ·
PTEN ← miR-19 · CDKN1B/ESR1 ← miR-221/222 · KRAS/MYC/HMGA2 ← let-7 · PDCD4/RECK/THBS1/TPM1 ← miR-21 ·
DNMT3A/DNMT3B/MCL1/COL1A1 ← miR-29 · SIRT1/CDK6/MET/NOTCH1/CD44 ← miR-34 · SOCS1/CEBPB ← miR-155.

Tests: (a) does the FAMILY β name the literature family, vs a 1/n_families permutation null; (b) how does it
compare to DOSE and to Shapley-identity; (c) **decoy FAMILIES** (strictly non-overlapping) — do they get β as
large as real families?

---

## 4b. ⛔⛔ CORRECTION TO §4 / Q4 / Q5 — **THE ATTRIBUTION VERDICT USED THE WRONG READOUT** (user-caught, 2026-07-13)

> **User:** *"if you are looking for a canonical arm that doesn't have a variance — it won't be there?"*

**§4's "attribution is at chance" was measured with `beta` and scored by `argmax`. BOTH ARE WRONG**, and the
subproject's own doctrine says so:
* `ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md`: identity = **Shapley/LMG on R²** (`attribution.shapley_identity`) —
  **NOT β**. (Additive Shapley on β splits NOTHING under collinearity.)
* `CN_INSTRUMENT.md §7`: *"a high-abundance, **LOW-VARIANCE** member delivers most of the **DOSE** yet earns ≈0
  variation-share ... **level/measurement is the right functional and variation-Shapley is the WRONG one here**."*
* The literature's *"miR-X represses gene Y"* comes from **knockdown/overexpression — a perturbation of LEVEL.**
  Asking a **covariance** estimator (β) to recover a **level** claim is an **ESTIMAND MISMATCH**, not a model failure.

**MEASURED (`scratchpad/variance_attribution.py` → `output/learned/variance_attribution.tsv`; 21 literature genes;
canonical guide arm → seed family via `families.family_of`; rank 0 = top of the gene's families, 0.5 = chance):**

| readout for the CANONICAL family | mean rank | Wilcoxon vs chance |
|---|---|---|
| **DOSE / level** (family mean) | **0.258** | **p = 0.0003** ✅ |
| **Shapley IDENTITY** (LMG on R²) | **0.319** | **p = 0.0076** ✅ |
| β (the coefficient) | 0.486 | p = 0.48 — **CHANCE** |
| *the canonical family's own VARIANCE* | 0.578 | p = 0.95 — **NOT a variance outlier** |

**⭐ The canonical regulator is HIGH-DOSE, NOT HIGH-VARIANCE** ⇒ a covariance functional cannot see it.
**And the two winning signals are INDEPENDENT:** Spearman(identity rank, dose rank) = **−0.023 (p=0.92)**;
dose-adjusted identity rank at a chance-level dose = **0.290, 95% CI [0.050, 0.455], bootstrap p=0.0068.**

**⚠ The variance-FLOOR hypothesis is REFUTED as the mechanism:** **0/21** canonical families are zeroed by the
`sd < 0.1` floor (their sd runs 0.67–1.21). The canonical family is in the design and gets a β — the β is just
**at chance**, because β is the wrong functional.

### ⇒ THE CORRECTED VERDICT

**The model DOES concentrate credit on the literature-canonical family, well above chance** — read with the
doctrine's own estimator (`shapley_identity`) and scored by **RANK**. It places it in the **top ~32%** of a
gene's 20–60 candidate families. **But it does NOT reliably make it #1** (argmax: identity 2/21, dose 3/21;
neither significant). ⇒ **RANKS IT HIGH; DOES NOT NAME IT.** That is materially better than "at chance" and
materially weaker than "recovers the literature". **Q4/Q5's "at/below chance" verdicts are RETRACTED as stated.**

⚠ This is the **THIRD** time in this project β was used where the doctrine specifies `shapley_identity`
(cf. MH-122, memory `read-the-doctrine-before-inventing`).

---

## 5. ⭐ CN-INSTRUMENT ORTHOGONAL VALIDATION — **EDGE EXISTENCE IS CONFIRMED** (Q2)

### 5.1 The raw instrument is ALSO confounded — it CANNOT validate on its own
`edge_instrument` over the v3 set (HE 818 / decoy 794 usable, F>10):
| | HE | abundance-matched DECOY | |
|---|---|---|---|
| usable (F>10) | 50.0% | 53.8% | decoys **more** instrumentable |
| ρ_reduced < 0 | 89.2% | **86.5%** | OR 1.29, p=0.055 |
| **mean ρ_reduced** | −0.1136 | −0.1087 | **p=0.21 — NO DIFFERENCE** |

**86.5% of SITE-FREE edges have a negative CN→target reduced form.** That is **ANEUPLOIDY** — a miRNA locus's
CN tracks global CN/CIN, which moves the target through paths C does not close. ⇒ **the RAW reduced form is
NOT an orthogonal validator.** (`CN_INSTRUMENT.md` says so: ~73% exclusion failure.)
⚠ Also: the abundance match **BREAKS** on the instrumentable subset (instrumentability itself tracks
abundance) — it must be **re-matched** after gating.

### 5.2 ⚠ DO NOT REINVENT — my hand-rolled δ_s FAILED verification
A hand-rolled single-arm δ_s differed from `instrument.exclusion()`'s by **|diff| up to 0.185** (CCND1). The
infra's `exclusion()` works on the **family aggregate X_fam** and groups instruments **by CN segment**.
**All hand-rolled exclusion numbers were discarded.**
✅ **The valid fast path (VERIFIED bit-identical, max|diff| = 0.00e+00, 9.1 s → 0.01 s):** call
`INS.exclusion(gene, [arm], assembled=(Y, X_custom, C, w), coding=True, host=True)` where Y and C come from
`assemble_gene(he_only=True)` (identical) and `X_custom` holds the arm's expression from the global X matrix
(so DECOY arms, absent from the HE design, can be instrumented).

### 5.3 ⭐ THE RESULT — exclusion-gated, abundance re-matched
Infra exclusion verdict (`coding=True, host=True`): **HE 34.5% clean / 65.5% pleiotropic**; decoy 31.0% / 69.0%
(median |δ_s| 0.078 / 0.086). *Canonical edges are not spared:* **CCND1←miR-16 is flagged `pleiotropic`,
δ_s = −0.23** (the RB1@13q14 co-amplification the docs name).

**On EXCLUSION-CLEAN edges, abundance RE-MATCHED (HE 234 vs DECOY 234; abundance MWU p=0.88 ✓):**
| | HE (curated) | **DECOY (site-free)** | test |
|---|---|---|---|
| **π_causal** (CN→target **THROUGH the miRNA**) | **−0.02349** (76% neg) | **−0.00754** (62% neg) | **MWU p < 1e-4** |
| π_raw (reduced form) | −0.03168 (78% neg) | −0.01779 (72% neg) | p = 0.0004 |

**VERDICT (Q2): ✅ EDGE EXISTENCE IS CONFIRMED BY THE EXOGENOUS INSTRUMENT.** The causal channel is **3×
stronger** for curated edges than for abundance-matched impossible ones (p<1e-4). **The signal ONLY appears
after the exclusion gate** — the raw instrument gives a FALSE NEGATIVE (p=0.21). **The exclusion stack is not
optional.**
Outputs: `output/learned/{cn_edge_validation,cn_infra_validation}.tsv`.

---

## 5b. FAMILY-LEVEL CN — **✅ EXISTENCE CONFIRMED · ❌ ATTRIBUTION UNDERPOWERED** (Q3)

**USE THE EXISTING INFRA. DO NOT REINVENT** (a prior attempt invented `argmax|wald_m|`, which was DEGENERATE:
co-located members share ONE CN column ⇒ identical γ/π/Wald ⇒ argmax decided by ROW ORDER; retracted).

`instrument.exclusion(gene, members, assembled=…, coding=True, host=True)` **groups instruments by CN SEGMENT**
— which IS the "aggregate several loci by the arms' loci" family instrument. Per gene-family, π_causal is the
**exclusion-weight–weighted mean over the family's segments** (F>10). Scripts: `cn_family.py`, `cn_fam_exist.py`,
`cn_fam_match.py`. Outputs: `cn_family_effects.tsv` (1,522 gene-families / 456 genes; 33.3% exclusion-clean),
`cn_family_existence.tsv`, `cn_family_literature.tsv`, `literature_headtohead.tsv`.

### 5b.1 ✅ Q3a — EXISTENCE: the family instrument separates HE from site-free decoy FAMILIES

Decoy family = strict (no member has a TargetScan site in the gene · not a cluster mate · |r|<0.30 vs every real
regulator), **abundance re-matched** on the instrumentable subset (the raw match BREAKS there, MWU p=0.000 —
exactly as at the arm level; abundance and repression predict the SAME sign, so it is not an arbiter until
matched — **axiom 4**).

| subset | HE π_causal | decoy π_causal | MWU p | seed-stability |
|---|---|---|---|---|
| all instrumented (684 vs 684, match p=0.78) | **−0.01639** (68% neg) | −0.00766 (63% neg) | **2.8e-05** | 10/10 draws p<0.05 |
| **exclusion-clean** (219 vs 219, match p=0.87) | **−0.02272** (73% neg) | −0.00831 (63% neg) | **7.7e-05** | 10/10 draws p<0.05 |

⇒ **Family→target edge EXISTENCE is confirmed by the exogenous instrument, at the SAME effect size as the arm
level** (arm: −0.02349 / 76% neg vs −0.00754 / 62% neg). Both resolutions carry real, site-dependent,
CN-driven repression. As at the arm level, **the signal only appears after the exclusion gate.**

### 5b.2 ❌ Q3b — TRUE IMPORTANCE: does the CN-CAUSAL estimate name the LITERATURE family?

The sharp question: the observational β failed (§4). Does the **exogenous** estimate succeed where it could not?
Ground truth = the canonical guide arm of each classic MTI, mapped to its seed family via `families.family_of`
(exact family match — **no substring heuristics**). Pruned: CDH1 (miR-200 acts *indirectly*, via ZEB1) and VEGFA
(not a canonical miR-200 target).

| estimator | names the literature family | chance | binom p |
|---|---|---|---|
| **learned β** (observational, full design) | **1/16 = 6%** | 8.1% | 0.74 |
| dose (abundance, full design) | 3/16 = 19% | 8.1% | 0.13 |
| **CN-causal π** (exogenous, instrumented only) | **4/16 = 25%** | **21.6%** | 0.46 |
| CN-causal π, exclusion-clean | 4/16 = 25% | 21.6% | 0.46 |

Paired: only-β 1 · **only-CN 4** · neither 11 (McNemar p=0.19). CN-causal's four hits are ALL canonical:
**ZEB1←miR-200bc-3p/429**, **BCL2L11←miR-17~92**, **MET←miR-34**, **CD44←miR-34**.

**⛔ VERDICT: NOT SETTLED — UNDERPOWERED, NOT NEGATIVE.** No estimator beats chance at n=16. The CN-causal π
numerically leads and wins the paired comparison 4-to-1, but **its chance baseline is 21.6%, not 8.1%** —
precisely *because* only ~7 of ~20 families per gene are instrumentable, which shrinks the candidate set.
**The binding constraint is the literature set (16 usable genes), not the estimator.**

> **⚠ RETRACTED (self-caught, same session):** a first pass reported **CN-causal 5/16 = 31.2%, p=0.09** using a
> substring matcher. That matcher (a) counted **PDCD4 ← miR-21-*3p*/3591-3p** as a hit for literature "miR-21",
> though the guide is miR-21-**5p** and the 3p family has a DIFFERENT SEED and different targets; and (b) a
> "strict" fix that blanket-skipped `-3p` families was *also* wrong — **for miR-200b/c the guide IS the -3p arm**,
> so it threw away the real ZEB1 hit. Only the **family_of** mapping is correct. Both earlier numbers are void.

* **EDGE EXISTENCE (arm level):** `instrument.edge_instrument(arm, gene)` → `first_stage_F` (usable = F>10),
  `rho_reduced` (<0 ⇒ genetic dose represses the target), `causal_concordant` (reduced<0 ∧ observational<0 ∧
  F>10). Compare the **causal-concordant rate for HE vs abundance-matched SITE-FREE decoys** — this is the
  orthogonal existence test, and it is immune to the observational background entirely.
* **PLEIOTROPY (all four kinds), via the infra:** `locus_origin` (silent-source loci are NOT valid
  instruments) · `pleiotropy()` (co-located miRNAs that also target the gene) · `between_family_exclusion()`
  (co-located OTHER-SEED families) · `exclusion(coding=True, host=True)` (co-amplified coding gene / intronic
  host, joint down-weight; **ACT ONLY ON REDUCTIONS**).
* **FAMILY level:** `family_multi_iv` (each member's focal-locus CN as a separate instrument; γ = dose
  delivery, π = reduced form) + `cn_copy` (the co-location-correct dose share) + `exclusion` (T1 δ_s gate).
* **Doc limits to respect:** CN's finest unit is the **CLUSTER**; it CANNOT separate -5p/-3p of one hairpin,
  nor hairpins in one polycistron (miR-18↔miR-20a CN ρ=1.00). It CAN separate across clusters/loci.

---

## 6. S + L (SPARSE + LOW-RANK) DECOMPOSITION — **⛔ REFUTED. The signal LIVES IN the low-rank part.** (Q7)

**Hypothesis (user):** the arm-resolution design may need `X = S + L` — sparse targeting signal + low-rank
shared structure. The low-rank part would be the background this task measured (compartment, proliferation,
global programs — the thing that makes σ₀ 3–4× too wide, and that lets site-free decoys reach 71% of the real
predictive field). If so, projecting `L` out should shrink the DECOYS far more than the REAL arms.

### 6.1 Design (and the trap it avoids)

**⚠ A naive S+L reproduces the E6 self-partialling bug**: factors computed from ALL arms CONTAIN the real
regulators, so projecting them out partials each HE arm on a covariate containing itself. The clean route is
**RUV** (Gagnon-Bartsch & Speed 2012): estimate the unwanted factors from **NEGATIVE-CONTROL FEATURES only**.

* **L** = top-k left singular vectors of the **negative-control arms** — arms that *cannot* regulate this gene
  (not in the design · not a seed-family mate · not a genomic-cluster mate of any design arm), computed
  **per gene, leave-the-design-out**.
* **S** = `X − proj_L(X)` — the sparse part. `L` is also appended to `C` for the outcome.
* Fit through the **canonical** path: `AE._prep(Y, S, [C|L])` → `spike_slab._gibbs_posterior(learn_pi0=True)`
  (prior-free base rate ⇒ real and decoy get the SAME prior — the decoy design's "matched prior").
* **BASELINE REPRODUCES §2.4**: REAL PIP>0.5 **21.1%** vs DECOY **14.3%** (§2.4: 23.0% / 15.5%).

**⚠ First attempt was DEGENERATE and caught:** run with `π≡1` (the *dense/coupling* readout) every PIP = 100%
by construction — the inclusion prior IS 1. PIP is only meaningful under `learn_pi0`. Those numbers are void.

### 6.2 Result — MONOTONE, and the OPPOSITE of the prediction

| k (RUV factors removed) | REAL PIP>0.5 | DECOY PIP>0.5 | **separation (REAL−DECOY)** | decoy \|β\| share |
|---|---|---|---|---|
| **0 (baseline)** | **21.1%** | 14.3% | **+6.8%** (p=8.2e-05) | 49.0% |
| 2 | 15.5% | 13.0% | +2.5% | 46.7% |
| 5 | 8.8% | 7.5% | +1.3% | 49.0% |
| 10 | 4.5% | 4.2% | **+0.3%** | 47.6% |

**Removing the low-rank component destroys the REAL signal FASTER than the decoys**, monotonically in k, and
drives the real-vs-decoy separation to **zero**. The decoy |β| share **never moves** (49% throughout) — S+L
does not reclaim a single point of the β the decoys absorb.

**⚠ MY PREDICTION (P2: "decoys will collapse, real will hold") WAS WRONG.** Per the measured-only gate, that
downgrades my reasoning in this domain: rank-structure intuitions here must be tested, not argued.
(P3 — "S+L will not rescue attribution" — is moot: it doesn't rescue existence either.)

### 6.3 ⭐ Why — and what it means for the whole task

The top miRNA-matrix factors **ARE the background axes**: mPC1 (14.2% of variance) ↔ **proliferation** r=0.33;
mPC2 (6.4%) ↔ **purity/CPE** r=0.41. But **that is also where the biology is.** miRNA repression *acts through*
the global programs — proliferation and EMT are not merely nuisance, they are the substrate the Hallmark miRNAs
regulate. This is **axiom 2a at the design-matrix level**: a plausible-mechanism covariate, removed globally,
**over-controls**. Rank cannot tell "the program that confounds" from "the program that is repressed" — **the
signal and the artifact occupy the SAME subspace**.

⇒ **The identification CANNOT come from inside the miRNA×sample matrix.** No endogenous re-basis of X — S+L or
otherwise — separates them, because there is no direction in that matrix that only one of them occupies.

⇒ **It must come from an EXOGENOUS source, and it does: the CN instrument (§5, §5b).** That is exactly
**axiom 4's** prescription — *"the identification strategy is an ASYMMETRY the confound is blind to"*. Germline/
somatic locus copy number is blind to compartment and proliferation; the decoy's 3′UTR has no site. **This is
why the instrument separates HE from decoys (π_causal −0.0227 vs −0.0083) where the rank decomposition cannot.**

---

## 6b. ⛔ A CORRECTNESS BUG IN THE CANONICAL GIBBS — FOUND, FIXED, VERIFIED (and it overturns MH-119's *cause*)

**`spike_slab._rtnorm_pos` silently violated its own support.** The half-normal slab asserts **β ≥ 0**
(repression ⇒ positive weight, METHODS §5). The sampler drew from `TN_{[0,∞)}(mu, s²)` by inverse-CDF:

```python
a = -mu / s;  lo = ndtr(a);  u = lo + (1-lo)*rng.random()
u = min(max(u, 1e-12), 1 - 1e-12);  return mu + s*ndtri(u)      # ← BUG
```
For `mu/s < -7.0345`, `ndtr(-mu/s)` **saturates to exactly 1.0** in float64 ⇒ `lo = 1.0` ⇒ `u` clips to
`1-1e-12` ⇒ the "draw" becomes the **deterministic constant** `mu + s·ndtri(1-1e-12) = mu + 7.0345·s`, which is
**NEGATIVE**. **MEASURED** (`ndtri(1-1e-12) = 7.0344869100`):

| mu/s | frac of draws NEGATIVE | unique values in 200 draws |
|---|---|---|
| −7.034 | 0% | **1** (already degenerate) |
| **−7.05** | **100%** | **1** |
| −9.0 | 100% | 1 |

**Persisted contamination:** **161/5117 = 3.15%** of `readouts_edges.tsv` β (129 genes) and **179/5802 = 3.09%**
of `readouts_arm_edges.tsv` (138 genes) are **negative — impossible under the model**.

**FIX:** Robert (1995) exact one-sided truncated-normal sampler (rejection from N(0,1) when the mode is inside
the support; **exponential proposal** with optimal rate `α = (a+√(a²+4))/2` in the far-left tail).
**VERIFIED two ways:** (1) **support** — 0 negative draws at every `mu/s` from −20 to +5, and 4000 unique values
(the determinism collapse is gone); (2) **correctness** — KS test against the exact `scipy.stats.truncnorm`
**PASSES at every mu** (D ≤ 0.016, all p > 0.01; sample means match the analytic TN means). As `mu/s → −∞` the
draw now → 0⁺, which is what the model means by *"this edge is off"*.

### ⭐ It overturns the CAUSE recorded for MH-119 (not the fix)

MH-119 found `share` (β_f/Σβ) exploding to **43.7** and recorded the cause as **"SIGN CANCELLATION — 100% of bad
rows have an anti-repressive β"**. The observation is right; **the causal claim is WRONG.** Those anti-repressive
β **are this bug**. The half-normal model *cannot* produce a negative coefficient — **the impossible value WAS
the alarm, and it was read as biology.** (The `share_sd` = 316 alarm was noted and under-read for the same
reason.) The MH-119 **deliverable stands** (a bounded `share_abs` + reliability flags are correct regardless);
only the mechanism is re-attributed. ⚠ `ratio-readouts-need-a-denominator-gate` (memory) corrected accordingly.

### It did NOT manufacture any MH-124 conclusion — RE-VERIFIED

| MH-124 result | buggy sampler | **FIXED sampler** |
|---|---|---|
| decoy test (§2.4): REAL vs DECOY PIP>0.5 | 21.1% / 14.3% (gap +6.8%) | **21.0% / 14.4% (gap +6.6%, p=8.9e-05)** |
| S+L (§6, k=5) separation | +1.3% | **+1.1%** |
| decoy \|β\| share | 49.0% | **48.7%** |
| negative β in the refit | (3%) | **0 / 3326** |

**Q1b and Q7 hold unchanged.** §5/§5b (the CN instrument) are **2SLS and never touch the Gibbs** ⇒ unaffected by
construction.

**DOWNSTREAM RIPPLE (axiom 2):** `readouts_{edges,genes}.tsv` + `readouts_arm_{edges,genes}.tsv` **REGENERATED**
with the fixed sampler. Consumers to re-check: `eval/coupling_grid.py`, `eval/attribution_showdown.py`,
`eval/within_family_cn_iv.py`, `analyses/cptac/composition_pinpoint_lineage.py`.

---

## 7. VERDICTS (filled as each test lands)

| question | verdict | evidence |
|---|---|---|
| Q1 anti-couplings meaningful vs no-site? | **YES at SET level; NO per edge** (E2, E4, §2.4) | OR 1.41 abundance-matched, survives composition; per-edge z≈−1.4 |
| Q1b **is it just abundance?** | **⭐ NO — but the margin is SMALL** (§2.4) | abundance-MATCHED strict decoys: PIP 23.0% vs 15.5% (p=1.1e-04), prediction 0.146 vs 0.104 (p=2.1e-07) — **but decoys reach 71% of the real field and take 36% of the |β|** |
| Q2 **CN existence validation (arm)** | **⛔ RETRACTED — `pi_causal` is NOT exogenous (a·b, b = observational OLS). Clean reduced form: arm-clustered p=0.115, n.s.** (see banner) | exclusion-gated + abundance-matched: **π_causal −0.0235 (HE) vs −0.0075 (decoy), p<1e-4**. ⚠ the RAW instrument FAILS (p=0.21) — aneuploidy makes CN→target negative for 86.5% of impossible edges. **The exclusion gate is essential.** |
| Q3a **CN existence validation (FAMILY)** | **⛔ RETRACTED — same `pi_causal` defect** (see banner) | exclusion-gated + abundance-matched: **π_causal −0.0227 (HE) vs −0.0083 (decoy), p=7.7e-05**, seed-stable 10/10. Same effect size as the arm level. |
| Q3b **CN → right family? (true importance)** | **❓ UNDERPOWERED — not negative** (§5b.2) | exact seed-family match, n=16: β 6% · dose 19% · **CN-causal 25% vs 21.6% chance (p=0.46)**; paired only-CN 4 vs only-β 1 (McNemar p=0.19). Hits all canonical (ZEB1←miR-200, BCL2L11←miR-17, MET/CD44←miR-34). **Literature set is the constraint, not the estimator.** |
| Q4 per-arm weights vs literature | **⚠ RETRACTED as stated (§4b)** — the β/argmax test was the WRONG READOUT | β at chance (42.9%, p=0.38) **but that is a covariance functional on a LEVEL claim** |
| Q5 **family model → right family?** | **⭐ CORRECTED (§4b): YES BY RANK, with the RIGHT readout** | **Shapley IDENTITY rank 0.319 (p=0.0076)** and **DOSE rank 0.258 (p=0.0003)**, independent of each other; **β is at chance (0.486, p=0.48)** — β is a COVARIANCE functional and the canonical regulator is **high-DOSE, not high-variance** (variance rank 0.578, p=0.95). **Ranks it high; does not name it** (argmax still n.s.). |
| Q6 within-family attribution | **DOSE, not coupling** (E8) | 3 independent ground truths |
| Q7 **S+L (sparse + low-rank)?** | **⛔ REFUTED — the signal LIVES IN the low-rank part** (§6) | RUV factors from negative-control arms; removing them destroys REAL faster than DECOY, **monotonically**: separation +6.8% (k=0) → +2.5% (k=2) → +1.1% (k=5) → **+0.3% (k=10)**; decoy \|β\| share never moves (49%). **My prediction was WRONG.** ⇒ signal and background occupy the SAME subspace ⇒ **identification must be EXOGENOUS** (which is why the CN instrument works). |
| **BUG: `_rtnorm_pos` broke the β≥0 support** | **⛔ FOUND + FIXED + VERIFIED** (§6b) | `ndtr` saturates at `mu/s < −7.0345` ⇒ deterministic NEGATIVE draw; **3.15%** of persisted β impossible. Fixed (Robert 1995); KS-verified. **Overturns MH-119's *cause*** (it was the bug, not "sign-cancellation biology"). **MH-124's conclusions re-verified unchanged.** |
