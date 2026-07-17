# MH-130 — THE GENE LENS: a per-gene competence map for the learned miRNA model

> **Goal:** Answer, **per gene and a priori**, three questions the aggregate never could — (a) can this gene's
> miRNA regulation be **measured at all**? (b) if yes, does the **curated edge set** beat a matched fake there,
> and does the **learning** add anything? (c) does it **transfer**? Genes that fail (a) are reported as
> **NOT MEASURABLE** and are never averaged in as evidence.
> **What belongs here:** the gene-resolution lens on the learned model — the noise floor, the competence
> classes, the model's honest domain of validity, and the adjudication of the variance-match / SNR / dose-
> variance / protein-discrepancy questions that this arc settled. **Does NOT belong here:** biology claims
> (→ `BIOLOGICAL_SYNTHESIS.md`), estimator specs (→ `LEARNED_MODEL_METHODS*.md` twins), or the decoy
> construction itself (→ `MH127_DECOY_MODEL_GENE_BUDGET.md`, `MH128_DECOY_RESOLUTION_AND_CN_GOLD.md`).
> **Update trigger:** any change to the gene universe, the decoy construction, the C blocks, or the OOF
> ladder; any new cohort/layer added to the transfer grid; any re-scoring of FAKE3.
> **Sync-partner:** `MH127_DECOY_MODEL_GENE_BUDGET.md` (its stratified numbers are amended here) and
> `docs/DISCOVERY_REGISTRY.md` rows MH-127/MH-127b/MH-127c/MH-130/MH-131/MH-132.

**Every number below carries its five coordinates: RESOLUTION (arm/family) × WEIGHTING (uniform/dose/beta/
shapley_identity) × COHORT (TCGA-OOF / CPTAC) × LAYER (mRNA / protein) × C BLOCK (core / composition).**
A number without them is not usable and does not appear here.

C-block vocabulary: **core** = CPE/purity + proliferation + target_cn (MH-127's block, `mh127_core`) or
purity + CIN (`coupling_grid`'s block, `grid_core`) — both are reported wherever they differ.
**composition** = core + the 8 deconvolution fractions.

---

## 0. The one-paragraph verdict

The learned model is **not** a model of the HE universe. It is a model of **~27% of it**, and that 27% is
identifiable **a priori, without touching Y or any fit**: genes with **≥3 curated seed-family regulators AND
above-median curated edge evidence (`w_max`)**. Inside that class the curated edge set beats an
abundance-matched, site-free decoy (TCGA-OOF / mRNA / family / β / core: gap **−0.0376, p=0.0038** FAKE1;
**−0.0282, p=0.043** FAKE2) and the frozen β transports to CPTAC. Outside it the gap is **+0.0002, p=0.72 —
exactly zero**, and in the 17.6% of genes with **no explainable ceiling at all** the *decoy out-predicts the
real model*. The famous **78.9% real-vs-fake "tie" is dilution**: it is the average of a real effect over a
quarter of the universe and a null (sometimes reversed) over the rest. Two things this arc did **not**
establish: that the decoy's SD mismatch is harmless (that rebuttal fails arc-wide BH), and that the CPTAC
**protein** transfer is evidence for the *curated* edge set (a site-free fitted fake reaches 99% of it).

---

## 1. ⭐ LEAD — WAS THE SNR STRATIFICATION CONFOUNDED BY A VARIANCE MISMATCH? **NO — BUT IT IS DEAD ANYWAY.**

The question this arc was convened to settle: MH-127's fakes are matched on **abundance** but never on **SD**.
If the fakes were systematically lower-variance, then "real beats fake in the low-regulator-variance stratum"
would be an artifact of the decoy being even lower-variance *there*.

**Measured (MH-129, TCGA, n=1041, X=log2(RPM+1), arm and family resolution, trimmed `in_caliper` design):**

| question | answer |
|---|---|
| Are the fakes SD-matched? | **NO.** `build_fakes._cost` = \|Δabundance\| only; SD never enters. Trimmed in-caliper pairs: FAKE1 SD **−4.3%**, p=3.35e-18 (n=7905); FAKE2 **−3.8%**, p=4.68e-19. Abundance (the matched axis) Δ=−0.003. |
| Did that mismatch manufacture the SNR strata? | **NO — it runs the WRONG WAY.** Inside the exact `reg_var` Q1 stratum I had claimed, the **fake's SD is +16.2% HIGHER** than the real's (arm, p=1.4e-09; FAKE2 +17.3%). Inside the `reg_abund` Q4 stratum: SD gap **−1.9%, p=0.78 = no mismatch**. Gene-level spearman(gap_fit, signed d_sd) = **+0.001, p=0.97**. |
| So is the SNR result safe? | **NO. It is REFUTED by a bug in my own annotation code.** `snr_strat.py` ran `np.log2(X+1)` on a matrix **already** in log2(RPM+1). Consequences, measured: spearman(true VAR, double-logged VAR) = **+0.154** (the axis was not variance); spearman(abundance, VAR) **+0.557 → −0.587 = SIGN FLIP**; in `mh127_snr_strata.tsv` spearman(reg_abund, reg_var) = **−0.834** and varQ1 ∩ abundQ4 Jaccard = 0.585. **The "two independent SNR findings" were ONE finding, and "low variance" was "high abundance" in disguise.** |

**RETRACTED — the corrected numbers** (TCGA-OOF / mRNA / **family** / **β** / **core**, 5-seed mean, trimmed
design; gap = ρ_real − ρ_fake, negative = real wins):

| claim (as published) | corrected |
|---|---|
| "real beats fake at HIGH regulator abundance", gap −0.036, **p=0.008** | **gap −0.0200, p=0.086** (FAKE1) / **−0.0213, p=0.156** (FAKE2) — **NOT SIGNIFICANT** |
| "real beats fake at LOW regulator variance", **p=0.028–0.049** | **gap −0.0045, p=0.539** (FAKE1) / **−0.0035, p=0.668** (FAKE2) — **GONE** |
| "fake beats real at LOW regulator abundance", +0.020, p=0.040 | not reproduced on the correct scale |

The only surviving stratum on the corrected axis (true variance **Q4 high**, gap −0.0341, p=0.023) appears in
FAKE1 only (FAKE2 p=0.155) **and** is exactly where the decoy is most SD-handicapped (−44.3%, p=1.8e-29) —
so it is variance-confounded and must not be used. **Nothing in the SNR arc survives. Do not carry it forward.**
Independently re-measured in MH-131 (quartiles on the corrected axis, both fake sets, BH-corrected): abundance
quartiles monotone but not significant; variance quartiles non-monotone and null.

---

## 2. ⭐ THE NOISE FLOOR — HOW MANY GENES CANNOT BE MEASURED AT ALL

**Definition.** `ceiling_R2_oof` = 5-fold OOF R² (KFold 5, shuffle, seed 0 — the *same* folds as the ladder) of
an **unrestricted OLS** fit of the gene's **full real regulator design** on the C-residualised target, with C,
z-scoring and the sd<0.1 floor all decided on TRAIN only. Unrestricted OLS is the upper bound for **any linear
weighting** (real β, fake β, dose, uniform) — so a gene with ceiling ≤ 0 cannot be explained by *any* of them.

**TCGA / mRNA / n_eff=1041 / family resolution:**

| | core C | composition C |
|---|---|---|
| ceiling ≤ 0 (**NOT MEASURABLE**) | **250 / 1421 = 17.6%** | **385 / 1421 = 27.1%** |
| ceiling ≤ 0.02 | **734 / 1421 = 51.7%** | — |
| median ceiling | **+0.0185** | **+0.0058** |

Arm resolution gives 17.6% / 27.0% — **resolution is immaterial; the C block is what moves the floor**
(composition adjustment eats ~69% of the median ceiling and adds 135 noise-floor genes).

**The real headline is the tail, not the zero count.** Only 17.6% of genes are strictly unexplainable, but
**51.7% have a ceiling ≤ 2%** — half the HE universe can have at most 2% of its post-C residual variance
explained by any linear weighting of its curated regulators.

**And half the universe has no *shape* at all:** median n_arms = 2, median n_fam = 2; **684/1421 (48.1%) genes
have exactly ONE seed family**. For those, `Xz@β` is a monotone rescaling of the single column ⇒ **β ≡ uniform
identically** (verified here: max\|ρ_fam_beta_real − ρ_fam_uniform_real\| over all n_fam=1 genes = **0.0
exactly**). The question "does the learning add anything?" is **undefined** for half the universe — it is not
a null, it is a non-question.

⚠ The ceiling is **substantially design capacity, not biology**: spearman(ceiling, n_fam) = **+0.564**. Use it
as a **dilution filter** (where the model works), never as evidence for *why*. My declared target-detection
floor (median log2(TPM+1) ≥ 1 AND dropout ≤ 0.20) fails only **112/1421 = 7.9%** of genes — **target detection
is not the binding constraint; regulator-design width is.**

---

## 3. ⭐⭐ THE 78.9% TIE IS DILUTION

MH-127's headline (TCGA-OOF / mRNA / family / β / core, 899 genes, 5-seed mean): ρ_real −0.0465 vs ρ_fake
−0.0367 ⇒ **the fake reaches 78.9% of the real model, paired p=0.211 (n.s.)**. That number is a **mixture of a
strong effect and a null, dominated by the null.**

**Measured here on the competence map (same coordinates, FAKE1):**

| stratum | n | ρ_real | ρ_fake | fake as % of real | gap | p (Wilcoxon) |
|---|---:|---:|---:|---:|---:|---:|
| **ALL ladder genes** (MH-127's headline) | 899 | −0.0465 | −0.0367 | **78.9%** | −0.0098 | 0.211 |
| MEASURABLE only (ceiling > 0) | 724 | −0.0568 | −0.0387 | 68.3% | −0.0180 | **0.020** |
| MEASURABLE **and** n_fam ≥ 2 | 452 | −0.0755 | −0.0498 | 65.9% | −0.0258 | **0.0091** |
| **A_COMPETENT** (n_fam ≥ 3 **and** w_max > median, measurable) | 230 | −0.0963 | −0.0570 | **59.2%** | **−0.0393** | **0.0034** |
| **NOT MEASURABLE** (ceiling ≤ 0) | 175 | −0.0043 | −0.0284 | **661%** | **+0.0241** | **0.019** |

FAKE2: all-genes 84.9% (gap p=0.58) → A_COMPETENT 66.5% (gap **−0.0286, p=0.052**).

**The monotone ladder is the finding.** And at the noise floor the decoy does not merely tie — **it beats the
real model** (gap **+0.0241, p=0.019** FAKE1; **+0.0232, p=0.021** FAKE2), which is what a fitted-on-noise
predictor does when the real design has no signal to protect it.

MH-131 reached the same conclusion on an **exchangeable, decoy-symmetric** stratifier
(`ceil_sym = (ceil_real + ceil_fake)/2`, for which E[gap \| ceil_sym] = 0 under the null **by construction**,
so it is not selection on the outcome): noise floor 33.1% of genes, fake/real **255%**; top tertile 22.4% of
genes, fake/real **51.1%**, gap −0.0712, p=1e-4, cluster-bootstrap p=0.0005. That stratifier **survived** a
Y-permutation calibration of its own selection operator (artifact share ≤5%) and a fully non-circular
patient sample-split (gap −0.055, boot p=0.0018), and the 18-cell headline grid survives **one global BH over
all 771 cells of the arc** (q = 0.0014–0.042).

---

## 4. ⭐⭐ THE DELIVERABLE — THE PER-GENE COMPETENCE MAP

**Artifacts (both persisted):**
- `mirna_hallmark/output/learned/gene_atlas.tsv` — **1421 genes × 98 columns**. The a-priori atlas: ceilings
  (in-sample / adjusted / OOF R² / OOF ρ) × {arm, family} × {core, composition}; regulator dose & variance
  distributions (Gini/entropy on linear-RPM shares and on log2); dose–variance decoupling; collinearity
  (mean/max \|r\|, condition number); curated evidence weights (`w_max`, `w_sum`); Shapley concentration.
- `mirna_hallmark/output/learned/gene_competence_map.tsv` — **1421 genes × 130 columns**. The atlas **joined to
  the outcome**: TCGA-OOF ρ (real / fake, both fake sets, 5-seed mean, arm & family, β & uniform), the
  real-vs-fake gap, the fitting gain, and the CPTAC transfer (mRNA & protein × β & uniform × core & composition)
  — plus the class labels below.

**The classes.** (b) and (c) are decided by **fit-free** axes only — `n_fam` (design width, arm-count-matched
in the decoys) and `w_max` (external curated edge evidence). (a) uses Y but is blind to the real-vs-fake
comparison.

| class | rule | n | what it means |
|---|---|---:|---|
| **NOT_MEASURABLE** | ceiling_R2_oof (fam/core) ≤ 0 | **250 (17.6%)** | No linear weighting of the curated regulators explains anything OOF. **Never average these in as evidence** — the decoy beats the real model here. |
| **D_UNDEFINED** | n_fam ≤ 1 (and measurable) | **500 (35.2%)** | One seed family ⇒ β ≡ uniform **exactly**. The "does learning add anything" question cannot be asked. |
| **C_WEAK** | n_fam = 2 and w_max ≤ median | **122 (8.6%)** | Measurable, but the decoy ties (fake = 99.2% of real, gap −0.0005, p=0.68). |
| **B_PARTIAL** | n_fam ≥ 3 **or** w_max > median (not both) | **176 (12.4%)** | gap −0.0199, **p=0.23** — not significant. |
| **A_COMPETENT** | n_fam ≥ 3 **and** w_max > median | **373 (26.2%)** | The domain of validity — see §5. |

(Of the 384 genes meeting the fit-free A rule, only **11** fail the measurability gate ⇒ **the fit-free axes
nearly imply measurability**, which is why the class can be declared a priori.)

---

## 5. ⭐⭐ THE MODEL'S HONEST DOMAIN OF VALIDITY

**It exists, it is a-priori, and it is 27% of the universe.**

> **CLASS A — `n_fam ≥ 3` AND `w_max > median(w_max)`. 384/1421 genes = 27.0%.**
> Both axes are **fit-free** (design width + curated evidence weight); neither touches Y, β, or the decoy.

**Measured inside it** (all five coordinates given per row):

| readout | resolution | weighting | cohort | layer | C | value | p | n |
|---|---|---|---|---|---|---|---|---:|
| real-vs-fake gap, FAKE1 | family | β | TCGA-OOF | mRNA | core | **−0.0376** | **0.0038** | 238 |
| real-vs-fake gap, FAKE2 | family | β | TCGA-OOF | mRNA | core | **−0.0282** | **0.043** | 193 |
| ρ_real / ρ_fake (FAKE1) | family | β | TCGA-OOF | mRNA | core | −0.0944 / −0.0568 (fake = **60.2%**) | — | 238 |
| transfer, β | family | β | CPTAC | mRNA | core (grid) | **−0.0779** | **1.0e-11** | 300 |
| transfer, β | family | β | CPTAC | **protein** | core (grid) | **−0.0727** | **9.8e-11** | 300 |
| transfer, uniform | family | uniform | CPTAC | protein | core (grid) | −0.0098 | 0.30 (**null**) | 300 |
| paired β − uniform | family | — | CPTAC | protein | core (grid) | **−0.0629** | **4.3e-16** | 300 |
| transfer, β | family | β | CPTAC | protein | **composition** | **−0.0172** | **0.034** | 300 |

**The complement (1037 genes) is empty:** TCGA-OOF / mRNA / family / β / core gap = **+0.0002, p=0.72**.
CPTAC / protein / family / β / core = −0.0153, p=0.0088 — a weak residual, ~5× smaller than in class A.

**Read this precisely.**
1. **In-cohort (TCGA-OOF, mRNA), class A is where the CURATED EDGE SET wins** — the decoy is beaten in both
   fake sets, and this replicates. But **the learning does not add anything in-cohort**: the paired
   (gap_β − gap_uniform) is −0.0104 (arm, p=0.225) / −0.0063 (fam, p=0.62) in the strongest stratum. ⚠ That is
   an **underpowered null, not a zero** — its 95% CI admits up to ~⅓ of the separation (MH-132).
2. **Out-of-cohort, the LEARNING is what carries it** — uniform is flat on CPTAC protein (−0.0098, p=0.30)
   while β is at −0.0727. This is the same asymmetry MH-127 found (transfer needs **both** the curated edge set
   and the learned weights) and it is sharper inside class A.
3. **Composition survives in class A** (protein −0.0172, p=0.034) where it is null on the whole universe
   (−0.0073, p=0.064) — but at ~24% of the core-C magnitude. The honest reading remains: **most of the core-C
   protein transfer is compartment arithmetic**, and I have not separated the surviving fraction per gene.
4. ⛔ **The protein cell is NOT decoy-controlled — see §7.** Class A's mRNA transfer *is* (MH-131 §5: REAL
   −0.0281 p=0.016 vs FAKE +0.0035, gap p=0.012 FAKE1 / 0.0038 FAKE2, replicated; survived MH-132's
   match-artifact lens). Its **protein** transfer is not.

---

## 6. THE VARIANCE-MATCH CONFOUND — WHERE IT ACTUALLY STANDS

Three routes were run. **They do not fully agree, and the refuter wins where it won.**

| route | result |
|---|---|
| **(i) FAKE3 — a variance-matched decoy, built and scored** (2-D optimal assignment: abundance caliper 0.25 **and** log2-SD caliper 0.25; match achieved: abundance Δ=+0.1%, SD Δ=−0.3%) | On a **common-support** design where the REAL design is byte-identical across decoys (551 genes, 5 seeds): paired FAKE3 − FAKE1 delta = **−0.0003, p=0.63** (family/β); all 18 cells \|Δ\| ≤ 0.0036, all q_BH = 0.97. **Variance-matching the decoy changes the gap by nothing.** ⚠ **But this test is weak**: on the common-support pairs FAKE1's SD mismatch is only +1.6% (p=0.076), and the subset has no real-win to begin with (gap +0.007, n.s.) ⇒ close to a null-on-null. |
| **(ii) \|d_sd\| as a covariate** (full power, 899 genes — a caliper is **infeasible**: only 2.1% of top-abundance real arms have *any* legal 2-D-matched twin) | The mismatch **is** a genuine artifact leg: partial spearman(\|d_sd\|, gap \| real_sd) = **−0.0878, p=0.0084** (FAKE1) / **−0.0970, p=0.0051** (FAKE2), surviving conditioning on the real design's own SD; whereas partial(real_sd, gap \| \|d_sd\|) = **nothing** (p=0.66 / 0.54). ⚠ This **corrects MH-129**, which measured only the *signed* d_sd (null, p=0.90) and over-concluded "the SD mismatch does not predict the gap AT ALL". The **signed** mismatch does not; the **magnitude** does. |
| **(iii) robust median split** inside the unbiased ceil_sym top tertile: does the win survive where the decoy is well SD-matched? | Well-matched half gap **−0.0650, p=0.0103** (FAKE1, n=100) / −0.0547, **p=0.066** (FAKE2). ⛔ **THIS REBUTTAL FAILS.** Under one global BH over the arc's 771 cells it lands at **q=0.082** (FAKE1) and FAKE2 does not clear 0.05 raw ⇒ it **fails the arc's own both-fake-sets rule** (MH-132). |

**Status, stated without spin:** the decoy is **not** SD-matched (−4.1 to −4.3%, p≈1e-17, re-measured twice
independently). The mismatch is a **real but minor artifact leg**. The rebuttal that the win survives where the
decoy is well-matched is **underpowered, not established**. The artifact share of the real-win is bounded
somewhere in **16–68%** (the median split says ~16%, a linear extrapolation to \|d_sd\|=0 says 53–68%; the two
routes disagree and were **not** adjudicated). **What is NOT in doubt:** the SD mismatch cannot explain the
**stratum structure** (§3) — inside the ceiling strata the mismatch is flat, and inside the retracted SNR strata
it ran the wrong way.

⛔ **STRUCTURAL LIMIT (escalate).** A caliper-matched decoy test of the **high-abundance** stratum is
**infeasible with this pool**: the abundance caliper binds (SD-caliper feasible 75.8% of real arms, abundance
caliper 36.3%, 2-D 20.3%; in the top abundance quintile 2-D feasibility is **1.0%**). The abundant curated
regulators with **no legal twin** include miR-21-5p, miR-22-3p, miR-10b-5p, miR-30a-5p, miR-148a-3p, miR-182-5p.
Any future variance-controlled test must **adjust** (SD as a covariate), not **match**.

---

## 7. THE PROTEIN DISCREPANCY — RESOLVED, THEN THE CURATION READING REFUTED

**Two runners disagreed:** `eval/coupling_grid.py` said the TCGA-frozen β predicts CPTAC protein;
MH-127's `p2_transfer.py` said protein is null (p=0.38).

**Cause, isolated (MH-131, same 685 genes, same 101 participants, same rung):** it is **MH-127's decoy TRIM**,
not the gene set and not the C block. Restricting the grid's protein cell to MH-127's own 685 genes still gives
ρ=**−0.0377, p=1.2e-08**. Same genes / same samples / same weighting, protein, family, β:
**TRIM −0.0058 (p=0.30) → FULL −0.0228 (p=1.3e-04)** under MH-127's own C; **TRIM −0.0127 → FULL −0.0326
(p=1.2e-06)** under grid-core C. The trim drops **56.2% of arms**.

⇒ **MH-127's "the protein channel dies" is an ARTIFACT of the decoy-matching trim. Corrected: on the FULL
design the TCGA-frozen dense β transports to CPTAC protein at ρ = −0.0310 (p=4.0e-09, family/β/grid_core,
n=1098), and protein is the STRONGER layer (mRNA −0.0186).** The MH-127c registry row is corrected accordingly.
The *mechanism* the user proposed (the trim deletes the abundant arms, which is where the signal is) is **half
right**: it does delete them (+5.55 log2RPM, p=1.7e-221) — but a **random** trim of the same size costs a
statistically indistinguishable amount (MH127TRIM vs RANDTRIM p=0.129) ⇒ the loss is **design WIDTH**.
⚠ MH-132 downgrades even this: the width-not-abundance contrast is an **underpowered null**, not a demonstrated
equivalence. Carry "the trim cripples the design", not a mechanism.

⛔ **AND THEN THE CURATION READING IS REFUTED.** MH-132 built the protein-cell decoy MH-131 declared impossible
(a fitted **site-free** fake on the **FULL** design) and it reaches **99% of the real model's −0.031**
(the fake's own protein transfer p=9e-8). **The −0.031 is real arithmetic, but it is NOT evidence that the
curated edge set is right.** The protein transfer establishes only that a **fitted** miRNA field transports —
not that **these** edges do. **The one decoy-controlled transfer result in the arc is the CPTAC mRNA cell inside
class A** (§5.4).

---

## 8. THE DOSE–VARIANCE AXIS (MH-124 §4b) — STRUCTURE CONFIRMED, CONSEQUENCE REFUTED

**CONFIRMED at universe scale.** `var_rank_of_top_dose` (0 = most variable, 1 = least) has mean **0.5555**
(arm) / **0.5549** (family) against a **within-gene SD-shuffle permutation null** (2000 perms, null mean 0.500):
**p=0.0005 (arm) / 0.0015 (family)**, and **the effect grows with regulator count** (family, ≥3 regs 0.562,
n=457 → ≥6 regs **0.596**, n=171, p=0.0005). **The top-dose regulator is systematically among the LEAST
variable** — PTEN's miR-21-5p var_rank 0.877; ESR1's miR-22-3p **1.000** (the least variable of its 13
families); CDKN1A's miR-22-3p 0.852; ZEB1's miR-200c-3p 0.733. A covariance functional **is** structurally
blind to these regulators.

**Also measured:** within a gene's own regulator set, dose and variance are **decoupled** — median
spearman(dose, SD) = **+0.039** (arm, n=494) / **+0.091** (family, n=457), 46%/44% of genes negative. The global
across-arm relation (+0.557) does **not** transfer into the within-gene design.

**REFUTED — the blindness costs nothing measurable, and the proposed fix makes things worse.**
- gap vs `fam_var_rank_of_top_dose`: spearman **−0.056, p=0.211** (FAKE1) / **−0.047, p=0.32** (FAKE2); after
  partialling on n_fam (mandatory — ceiling ~ n_fam is +0.564): p=0.33 / 0.31. **Null, sign inconsistent.**
- **The operative test:** does a DOSE- or IDENTITY-weighted field succeed where β is blind? Paired
  (ρ_dose_real − ρ_beta_real), TCGA-OOF / mRNA / family / core, in the **BLIND** stratum (var_rank ≥ 0.75):
  **+0.0349, p<1e-4 — POSITIVE means DOSE IS WORSE.** ρ_β = −0.0806 vs ρ_dose = −0.0458: **β is 76% better than
  dose exactly where β is supposed to be blind.** shapley_identity − β = +0.0219, p<1e-4 (also worse).
  Replicated in FAKE2 (dose − β = +0.0297, p<1e-4).

⚠ **Scope.** This refutes the extension of §4b to the **coupling** estimand (TCGA-OOF, mRNA, Spearman of the
field vs the C-residualised target). MH-124 §4b's original claim was about **attribution** (ranking the
canonical regulator) — a **different estimand**, and it is **NOT overturned**. Do not carry "dose weighting is
dead" into the attribution arc.

---

## 9. UNDERPOWERED ≠ NEGATIVE — the explicit separation

| claim | status |
|---|---|
| Class A real-vs-fake win (TCGA-OOF/mRNA/fam/β/core) | **POSITIVE.** Replicates in both fake sets (p=0.0038 / 0.043); survives arc-wide BH on the headline grid. |
| The 78.9% tie is dilution | **POSITIVE.** Survives an exchangeable stratifier, a Y-permutation calibration of the selection operator, a non-circular patient split, and one global BH over 771 cells. |
| "The learning adds nothing **in-cohort**" | **UNDERPOWERED NULL.** p=0.19–0.62; 95% CI admits up to ~⅓ of the separation. Not a demonstrated zero. |
| "The SD mismatch does not manufacture the win" | **UNDERPOWERED.** Route (iii) fails global BH (q=0.082) and fails the both-fake-sets rule. FAKE3 changes nothing, but that test is near-null-on-null. |
| "The trim costs WIDTH, not abundance" | **UNDERPOWERED NULL** (RANDTRIM contrast p=0.129). The *premise* — the trim is not neutral — is **POSITIVE**. |
| The SNR strata (abundance / variance) | **REFUTED** (scale bug + corrected re-run). Not underpowered — **wrong**. |
| MH-124 §4b's *consequence* ("β fails on blind genes") | **REFUTED** (β beats dose *in* the blind stratum, p<1e-4). Its *structure* is **CONFIRMED**. |
| CPTAC **protein** transfer as evidence for the curated edges | **REFUTED** (a site-free fitted fake reaches 99%). The transfer *number* is positive; its *curation* reading is not. |

---

## 10. OPEN HOLES

1. **The protein cell has never faced a decoy that could lose.** Fakes exist only on the TRIM design; the trim
   is exactly what kills protein. The one decoy that reached the FULL design **tied**. This is the arc's biggest
   hole.
2. **The high-abundance stratum is not decoy-testable** with this pool (2.1% 2-D feasibility). Needs a different
   null: permuted edges within abundance bins, or SD entered as a covariate.
3. **The composition-adjusted universe** is where everything weakens (protein −0.0073, p=0.064 globally;
   mRNA sign-flips). Class A survives at −0.0172 (p=0.034) but at a quarter of the magnitude. Per-gene
   composition-class stratification **not run**.
4. **`ceil_sym` / `ceiling_R2_oof` is design capacity, not biology** (ρ=+0.564 with n_fam). It tells you *where*
   the model works, never *why*.
5. **Gene dependence:** the across-gene t-tests treat genes as independent. Measured ICC ≈ 0 and patient-bootstrap
   deff = 1.00 (MH-132) ⇒ the dependence attack failed, but the absolute p-values (1e-11 etc.) are still not
   calibrated for shared designs.

---

## 11. PROVENANCE

**Persisted artifacts.**
`mirna_hallmark/output/learned/gene_atlas.tsv` (1421 × 98) ·
`mirna_hallmark/output/learned/gene_competence_map.tsv` (1421 × 130) ·
`mirna_hallmark/output/learned/mh127/` (`ladder_oof.tsv.gz`, `p2_transfer_results.tsv.gz`,
`fake_regulator_sets.tsv.gz`, `gene_strata.tsv.gz`).

**Design (identical across the whole arc).** Gene universe = the 1421 HE genes in MH-127's `universe.pkl`.
Design = `learned.data.assemble_gene(w_prior_source="ledger", deconv=False, edges=<pooled-HE ∩ RNA-matrix ∩
arm-expression floor>)`; family resolution = `learned.families.collapse_by_family` (TRUE RPM pool
log2(1+Σ(2^x−1)), **not** a groupby mean); z-scoring and the sd<0.1 floor = `learned.attribution_eb._prep`
verbatim; β = `learned.spike_slab._gibbs_posterior` dense (π≡1, the coupling readout);
identity = `learned.attribution.shapley_identity`. 5-fold OOF = `KFold(5, shuffle, random_state=0)` over
participants, with C, z-params and the sd-floor **all decided on TRAIN only** — the same folds in every phase.
X = log2(RPM+1) throughout (the MH-129 double-log bug is not repeated anywhere downstream).

**Verification gates passed.** `gene_strata`'s Shapley columns re-derived from scratch on 30 random genes:
**bit-exact** (max\|diff\| = 0.0000). The generalized ladder reproduces MH-127's `ladder_oof` **bit-exactly**
(max\|diff\| = 0.00e+00, pearson = 1.000000 on all 4 legacy columns). MH-132's independent harness reproduces
MH-131 bit-exactly on 7 columns. The CPTAC engine reproduces `mh131_transfer_grid.tsv.gz` on 30 checked cells
to 5 decimals.

**Phases.** MH-129 (variance-match gate + FAKE3 build) · MH-130 (gene atlas + this competence map) ·
MH-131 (stratified ladder; protein-discrepancy settlement) · MH-132 (four adversarial refutation lenses:
match-artifact, selection-and-ceiling, multiplicity-and-power, protein-cell decoy).
Scratchpad code: `mh129_audit.py`, `mh129_true_variance_axis.py`, `mh129_build_fake3.py`, `mh130_gene_atlas.py`,
`mh130_competence_map.py`, `mh131_ladder.py`, `mh131_ceilfake.py`, `mh131_settle.py`, `mh131_whytrim.py`,
`mh131_transfer_grid.py`, `mh131_declass.py`, `mh132_*.py`.
