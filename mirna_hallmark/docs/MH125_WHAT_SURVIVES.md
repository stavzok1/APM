# MH-125 — WHAT SURVIVES: the honest ledger after the decoy/null audit

> **Goal:** answer, in one place and without spin, the question *"if fake edges anti-correlate almost as well as
> real ones, do any of our results survive?"* — with a mechanism, a claim-by-claim verdict, and a re-run list.
> **What belongs here:** the direct answer to the mechanism question; the measured real-vs-fake separation on
> every readout we have; the SURVIVES / NEEDS_RETEST / DIES table for the audited claims; the existence-vs-
> attribution distinction; what must be re-run; and what is genuinely unknown (as distinct from negative).
> **What does NOT belong here:** new biology (→ `BIOLOGICAL_SYNTHESIS.md`); model design (→ the METHODS twins);
> the design of the decoy tests themselves (→ `MH124_ANTICOUPLING_VALIDITY.md`, `MH127_…`, `MH128_…`);
> any number not produced in this session or read from a named persisted output.
> **Update trigger:** any new decoy-controlled, null-calibrated, or out-of-cohort test that changes a verdict
> in the table below; any re-run of a claim currently marked NEEDS_RETEST.
> **Sync-partner:** `DISCOVERY_REGISTRY.md` (the audited rows), `ANALYSIS_RUN_LEDGER.md` (MH-125 row),
> `MH124_ANTICOUPLING_VALIDITY.md`, `MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md`,
> `MH127_DECOY_MODEL_GENE_BUDGET.md`, `MH128_DECOY_RESOLUTION_AND_CN_GOLD.md`.

---

## ⭐ THE ONE THING TO KNOW

**The edge SET is real. The WEIGHTS are not attribution. And most of what died is a class of READOUTS, not the
biology.**

A fake (site-free) miRNA arm anti-correlates with a gene almost as easily as a real regulator does — **and that
is a fact about bulk tumour data, not evidence that our edges are fake.** It was always true; we just never
measured it. What it destroys is the *inference we drew from a bare anti-correlation* — FDR counts of
"significantly coupled" programs, the evidence-weighted pressure heuristic, and any statement of the form
*"miR-X is the driver of gene Y"*. What it does **not** destroy is the finding that curated regulators, as a
set, carry genuine site-driven repression signal beyond abundance and beyond composition — that survives every
control we have been able to build.

Of the **19 registry claims that were adversarially re-audited** in this arc: **1 survives as stated, 6 are
untouched, 7 need a retest, 5 die.** The ones that die are all built on the retired pressure heuristic + an
uncalibrated null. The learned model's core out-of-sample claim (β predicts CPTAC protein; abundance and the
heuristic do not) **is confirmed in the regenerated grid** — with one important qualifier in §2.4.

---

## 1. "How is that even possible?" — the mechanism, in plain language

### 1.1 Bulk tumour expression is dominated by a few global axes

Every sample in TCGA-BRCA differs from every other in **how much it is proliferating** and **how much of it is
tumour vs stroma vs immune** (purity / composition). Those two axes move *thousands* of genes and *hundreds* of
miRNAs at once. Measured, not assumed: the top factors of the miRNA × sample matrix **are** those axes —
mPC1 correlates with proliferation at **r = 0.33** and mPC2 with purity at **r = 0.41**
[`MH124_ANTICOUPLING_VALIDITY.md` §6].

So take any miRNA and any mRNA at random. If the miRNA happens to load positively on proliferation and the gene
loads negatively, they will anti-correlate across 1,041 patients — **with no physical interaction whatsoever.**

### 1.2 The size of the effect, measured

| quantity | value | source |
|---|---|---|
| SD of the per-edge coupling statistic **expected under independence** | **0.031** | theoretical t-null |
| SD **actually observed** among 60,000 site-free pairs (no TargetScan site, no scanMiR duplex ⇒ *cannot* repress) | **0.083 – 0.132** | MH-124 E1 |
| ⇒ the null we had been testing against was | **3–4× too narrow** | |
| fraction of **impossible** pairs that pass "FDR q<0.05" against that null | **25.0% (composition-adj) – 35.3% (core)** | MH-124 E1 |
| fraction of curated HE edges that survive an **honest** empirical null, per-edge | **0.0%** (from a nominal 31.3%) | MH-124 E4 |

**Read that last row carefully.** It does not say the edges are fake. It says **no single edge is individually
significant in bulk** — a typical HE edge sits at ρ ≈ −0.15 against a null SD of 0.09, i.e. z ≈ −1.4. The
signal is a **distributional shift of the whole curated set**, not a per-edge fact. Every claim we ever made of
the form "N of these edges are FDR-significant" was counting noise at the 25–35% rate.

### 1.3 Why the decoys are *supposed* to be good

The decoys are not random junk. By construction they are **abundance-matched**, **not seed-family mates**, **not
genomic-cluster mates**, **carry no binding site of any kind** — and, in the corrected build, they are matched
on the **global factor loadings** too. In other words a decoy is *an arm that looks exactly like a real
regulator in every respect except that it physically cannot bind the target.* That is the point. If such an arm
still anti-correlates, the anti-correlation was never about binding.

**And it does.** [`MH128_DECOY_RESOLUTION_AND_CN_GOLD.md` §1, un-handicapped fake, 3,104 head-to-head edges]:

| set | realized partial ρ (arm vs target \| C) |
|---|---|
| REAL curated HE edges | **−0.0353** |
| FAKE, axis-matched (impossible edges) | **−0.0227** |
| FAKE, MH-124's older strict rule (handicapped) | −0.0069 |

⇒ **Anti-correlation in bulk is cheap. The honest site-driven effect is ≈ −0.02, not the ≈ −0.039 MH-124
implied — every real-vs-decoy gap MH-124 published is ~2.25× too large** (its decoy was accidentally *stripped*
of global-axis loading by the `|r|<0.30` cap, which flattered the real set).

**This is the direct answer to the question.** Yes — an impossible edge reaches roughly two-thirds of a real
edge's anti-correlation. It is possible because both arms ride the same two global axes. It is a property of
bulk data. It would be true in any bulk cohort, for any pair of molecules.

---

## 2. What *does* separate real from fake

### 2.1 The readouts that separate (all observational)

| readout | REAL | FAKE | test | source |
|---|---|---|---|---|
| set-level enrichment, abundance-matched | — | — | **OR 1.41, p=2.5e-13**; survives composition (1.43 → 1.41) | MH-124 E2 |
| realized coupling, axis-matched fake | −0.0353 | −0.0227 | **MWU p=3.0e-04**; adjusted for abundance+loadings+degree, b(is_real) = **−0.0173, p=5.8e-07** | MH-128 §1 |
| posterior inclusion (matched prior ⇒ the gap is the DATA) | PIP>0.5 **23.0%** | **15.5%** | MWU p=1.1e-04 | MH-124 §2.4 |
| held-out gene-field prediction | 0.146 | 0.104 | Wilcoxon p=2.1e-07 — but see below | MH-124 §2.4 |
| decoy share of a gene's total \|β\| | — | **36.3%** | — | MH-124 §2.4 |

**Two corrections to how those last two are usually quoted:**
* the "**decoy reaches 71% of the real field**" figure is a *handicapped-decoy* number. Denominated correctly
  against the full learned model, the honest figure is **~0.75 (β) [0.69, 0.81]** [MH-128 §2a].
* the separation is **modest and set-level in every one of these lanes.** None of them licenses a per-edge or
  per-arm statement.

### 2.2 ⛔ The exogenous CN instrument does NOT separate them — do not cite the "3×"

This must be said plainly because it was, until very recently, our headline. **`instrument.exclusion`'s
`pi_causal` is not exogenous.** It is `γ_s · b_fam` — a first stage times an **observational OLS slope** — i.e.
a product-of-coefficients mediation estimator, not an IV. The "π_causal −0.0235 vs −0.0075, 3× separation" was
**reproduced from a simulation with zero causal effect** [`MH126_…` §2; retraction banner at the top of
`MH124_…`]. The only genuinely exogenous factor, `γ_s`, is **site-blind** (HE +0.1993 vs decoy +0.1984,
p=0.20), and the genuinely exogenous **reduced form** shows **no HE-vs-decoy effect** at the correct
(arm/locus) unit of analysis (**p = 0.115**, significant in only 2/10 matching draws).

⇒ **There is currently NO exogenous validation of edge existence.** Existence rests on **one** observational
line (§2.1), not three. *(That reduced form is **underpowered, not negative** — see §6.)*

### 2.3 The three decisive tests run this session

**T1 — the CPTAC decoy transfer test.** *(Does a decoy's β carry out-of-cohort, out-of-modality signal? The
program's prediction was NO.)* **REFUTED.** In CPTAC-2 prospective (n=101 independent patients, TMT-11;
2,585 pairs, 1,293 real / 1,292 decoy):
* decoy β predicts CPTAC **protein** at Spearman(β, realized ρ) = **−0.243, p=8.6e-19** (core C) and
  **−0.128, p=3.7e-06** (composition-adjusted). Real reaches −0.304 / −0.198.
* real beats decoy on β-transport **only at the protein layer** (Δ = −0.061 core / −0.070 comp; gene-cluster
  bootstrap one-sided **p = 0.048 / 0.040, both 95% CIs cross zero**) and **not at all at mRNA** (Δ = −0.011,
  p = 0.40).
* on **realized per-edge coupling**, after composition adjustment, real curated edges have **zero (slightly
  positive) mean coupling in CPTAC at both layers** and are **statistically indistinguishable from impossible
  edges** (protein p=0.32, mRNA p=0.14). Re-matching on *CPTAC-measured* abundance erodes even the core-C
  separation (p 0.003 → 0.057).

⇒ **"It replicates out-of-cohort and out-of-modality" is NOT an arbiter of edge specificity.** TCGA and CPTAC
are both bulk breast tumour: they *share the confound*, so **the artifact replicates too**. This is exactly the
trap recorded in memory as *null-design-and-shared-confounds* / MH-114.

**T2 — the `coupling_grid._learned_M` code audit + MH-115 re-derivation.** Four real defects were found and
**have since been fixed and the grid regenerated** [`mirna_hallmark/eval/coupling_grid.py`, carries the
"⚠ MH-125 FIX" notes; `output/learned/coupling_grid_summary.tsv`, regenerated 2026-07-13 18:22, **read
directly this session**]. The claim that the bug understated the protein effect by 2.8× was **refuted** (the
true move is ~1.05×). **MH-115's headline holds in the corrected, regenerated grid:**

| CPTAC, out-of-sample | learned β | abundance | heuristic |
|---|---|---|---|
| protein, core C | **−0.0347, p=7.8e-12** | −0.0110, p=0.023 | −0.0079, p=0.11 |
| protein, composition C | **−0.0131, p=7.4e-04** | −0.0009, p=0.83 | +0.0031, p=0.44 |
| mRNA, core C | **−0.0197, p=1.8e-05** | +0.0085, p=0.053 | −0.0009, p=0.83 |
| mRNA, composition C | −0.0054, p=0.12 | +0.0068, p=0.052 | +0.0023, p=0.50 |

⇒ **only the learned β predicts held-out CPTAC protein, and it survives composition adjustment** (p=7.4e-04).
The pressure heuristic is indistinguishable from mere abundance, and both are null — MH-115's retirement of the
heuristic stands. *(One MH-115 sub-claim flips: in TCGA the heuristic (−0.0459) now beats abundance (−0.0417),
not the reverse. And note the in-session scratchpad reconstruction of the corrected field gave the
protein/composition cell as p=0.055; the repo's own regenerated pipeline gives p=7.4e-04. **The persisted file
is the number of record; the disagreement is unresolved and worth one re-check.**)*

**T3 — the random-regulator null.** Substituting *every* real HE arm with an abundance-matched, non-family,
non-cluster, near-uncorrelated decoy and refitting the canonical Gibbs:
* the **gene-level FDR count** is **91–93% reproduced** by fake regulators (REAL 513/1,206 sig vs NULL 477).
  **The count carries essentially no information.** There *is* a real specificity component (paired real−null
  Δ = −0.012 / −0.035, Wilcoxon p = 9.7e-03 / 1.3e-03) but it is ~10–27% of the effect.
* the CPTAC-protein readout **does not rescue it**: the abundance-matched random-regulator null reaches
  **72%** of the real edge-level β↔protein coefficient, difference not distinguishable (gene-clustered
  bootstrap **p = 0.246**).

### 2.4 ⚠ THE THREE TESTS DISAGREE — and the disagreement is informative

* **T2 says** the learned β beats abundance and the heuristic on held-out protein, decisively.
* **T1 and T3 say** a *fitted fake* regulator set predicts held-out protein about as well as the real one.

Both are true, and the reconciliation is **measured, not argued** [MH-128 §2b]: **a fitted FAKE beats real
UNFITTED abundance** (gap −0.0617, Wilcoxon p=1.7e-19; cluster-block permutation p=5e-05), reproducing **61% of
the learned model's entire margin over abundance**. T2 compares *fitted vs unfitted*; T1/T3 compare
*fitted-real vs fitted-fake*.

⇒ **MH-115's "β beats mere abundance" must be read as a claim about FITTING, not about CURATION.** Roughly
**two-fifths** of β's advantage over abundance is attributable to the curated edge set; the rest is what any
flexible fit on the same global covariance structure would buy.

**One genuine, unresolved contradiction between our own arcs:** MH-127 measured that the *gene-field* transfer
to CPTAC **does** separate real from fake (real predicts CPTAC mRNA ρ=−0.0189, p=8.9e-04; matched fake does
not, p=0.19; paired p=0.030, and p=0.018 under FAKE2 with composition) — whereas T1 finds **no** real-vs-fake
separation at CPTAC at the **edge** level after composition. Different unit (gene field vs edge), different
decoy build. **We do not currently know which resolution is right. Do not quote either as settled.**

---

## 3. ⭐ THE SURVIVES / DIES TABLE

**Scope.** This table covers the **19 claims re-audited adversarially in this arc (MH-1 … MH-24, the Aim-1
Hallmark / pressure-heuristic block)**. Every one was verified by re-running or re-reading its own output; the
verifier's status wins where it disagreed with the first-pass auditor. **MH-4, MH-5, MH-9, MH-10 and MH-14 were
not in the audited set** and are unassessed here. The ~120 later registry rows (MH-25 …) were **not** re-audited
in this arc — but note that any of them resting on `compute_gene_pressure`, on a theoretical per-edge FDR count,
or on a β-argmax is exposed to exactly the failure modes below and should be treated as suspect until checked.

### DIES — cannot be restated faithfully and truly (5)

| claim | one-line reason |
|---|---|
| **MH-2** — "strongest anti-correlations are bile-acid / estrogen / apical / adipogenesis / myogenesis" | The **superlative is the claim**, and the ranking is not robust: under a composition-adjusted C the top-5 becomes apoptosis/p53/IFN-γ/ROS/OXPHOS (none of the five named); the claim also never matched its own file (EMT was #1 by partial ρ); estimator is the retired heuristic, which mere abundance reproduces at rank ρ=0.87. *(The anti-correlations themselves are real and composition-robust — 40/50 programs qualify — which is why the superlative is the content and the content is gone.)* |
| **MH-12** — "Basal concentrates pressure on a TSG/apoptosis node set" | It is a **curation-degree map, not a basal map**: Spearman(in-degree, basal pressure) = +0.24 (p=3.7e-19); the same hubs top the **luminal** list (ρ(luminal, non-luminal) = 0.46 heuristic / **0.9964 learned**); degree-matched, the 9-gene set as written is n.s. (p=0.15). |
| **MH-16** — "Basal enriched for EMT/apoptosis/TSG pressure; per-subtype driver arms named" | Basal pressure is positive for **82.4% of ALL genes** (a global shift); against that background EMT shows **no** enrichment (p=0.913). Drivers are **not subtype-specific**: 85.9% of genes get an identical top-5 driver set in all four subtypes. Apoptosis alone clears background (p<1e-3) and is a legitimate new question. |
| **MH-17** — "basal coupling survives and strengthens under proliferation adjustment, 8/8 programs" | A **zero-edge** predictor — mean z(log2 RPM) over all arms, no targeting information at all — also gets **8/8** (ρ −0.20 to −0.43); an arm-label-permuted decoy gets **7.67/8**. The statistic cannot see the edge set, so it cannot close a specificity objection. |
| **MH-24** — "orphan coupling clusters on seed families at TSG nodes" | Verdict delivered as DIES by the verifier; **the reason text was truncated in transmission to this document.** Do not re-assert the claim; re-read the audit output before rewriting the row. |

### NEEDS_RETEST — the number is not defensible; the phenomenon may be (7)

| claim | one-line reason |
|---|---|
| **MH-1** — "24/50 Hallmark programs show significant negative pressure↔expression coupling" | Against a **degree-preserving edge shuffle**, decoys get **27.6/50** — *above* the real 25/50 (8/8 draws ≥ real). And a program's own pressure explains its own expression no better than a foreign program's (cognate ρ −0.065 vs non-cognate −0.063; cognate rank 22/49 = chance). The count carries zero wiring information. Retest with the learned β + a shuffle null + a cognate-vs-foreign specificity test. |
| **MH-7** — "proliferation programs persist, immune programs collapse" | The split is a **covariate-set artifact**: purity is in the model, proliferation is not. With proliferation actually adjusted (MKI67), **E2F persists (+0.125, p=1e-4) but G2M does not (p=0.10)** and MITOTIC_SPINDLE sign-flips. The immune half is *stronger* than claimed (it reverses to significant negative). Also the registry text says "CPE+PAM50"; the file is CPE+HRD+batch. |
| **MH-8** — "estrogen/metabolic anti-correlations are driven by Normal-like tumours" | **Contradicted by its own regenerated file**: with Normal-like *deleted*, ER-early is still ρ=−0.286 (q=1.2e-20) and bile-acid −0.357 (q=4.4e-32). And "within LumA/LumB ρ≈0" is half false — cholesterol/bile-acid survive in LumA at q=3.4e-06. The Normal-like stratum no longer exists in any current output ⇒ literally unmeasurable without a rerun. |
| **MH-13** — "the most impactful miRNA per Hallmark is subtype-dependent; intra-family arms oppose" | Only the word **"impactful"** fails: the argmax is curation-depth × abundance (curation flips the winner in 30/50 luminal, 44/50 non-luminal; bootstrap top-1 reproducibility 72–85%). The **subtype-dependence** (50/50 flips vs a null of 25.1, p<0.005) and the **intra-family opposition** (miR-29a +0.42 / miR-29c −0.90, p=3.4e-23) are **verified** — reword as a DOSE claim, do not delete. |
| **MH-21** — "miR-17~92 coupling is basal-maximal; only top-load family coupled (ρ −0.45, 44/50)" | The **comparative core survives** a matched-decoy null (real 44/50, ρ −0.456 vs decoy 16/50, −0.097; p=0.032 — so it is not cheap background). But **"basal-maximal" and the −0.45 magnitude are manufactured by the single proliferation covariate**: drop it and Basal is −0.080, *no deeper than LumA*, with Her2 deepest. `prolif_verdict.tsv` classifies 3/5 of this claim's own nodes as **over_control**. Retest = adjudicate proliferation as confounder vs mediator (miR-17~92 is a direct E2F/MYC target *and* represses cell-cycle brakes ⇒ collider-adjacent). |
| **MH-22** — "CDKN1A/TGFBR2/VIM survive target-CN + dynamic-range controls (p<0.01)" | Per-gene p-values against the null that is **3–4× too narrow**; at σ₀≈0.09–0.13, VIM (ρ=−0.17) does not survive an honest null and CDKN1A is borderline. The negative verdicts (BIM, IRF1) are safe. No composition block. |
| **MH-23** — "TargetScan orphans add no Basal-specific layer (both reach 42/50)" | The **site-free decoy floor of this readout is 45–46/50** — *above* both numbers being contrasted. The test has no demonstrated power, so it can neither detect nor exclude additivity. (The claim's *direction* is not refuted; the downstream advice — don't flat-union HE ∪ TS-orphan — is unharmed.) |

### SURVIVES — re-tested against a harder null and still standing (1)

| claim | one-line reason |
|---|---|
| **MH-19** — "per-gene robustness splits the hub: a proliferation-robust TSG core survives; IRF1/BCL2L11 do not" | **Re-scored against an empirical abundance-matched site-free null** (32,215 pairs; measured σ₀=0.127 vs theoretical 0.076): **TGFBR2←miR-17-5p/miR-20a-5p and PTEN/CDKN1A←miR-106b-5p clear it** (BH q<0.05, E2F/MYC-independent proliferation proxy, CPE+HRD+batch). **IRF1 and BCL2L11 sit exactly at the decoy null under every proxy (0/5)** — the negative half is fully confirmed. ⚠ **Amend the numbers: 9/14 → 4/14 (2/14 under the E2F-independent proxy); VIM DIES (q=0.084); CDKN1A/PTEN survive via miR-106b-5p ONLY, not "all three arms".** |

### UNAFFECTED — nothing in this arc touches them (6)

| claim | one-line reason |
|---|---|
| **MH-3** — canonical programs enriched for HE target genes | Pure hypergeometric overlap of the **curated annotation** with Hallmark sets — no expression, no coupling, no null in dispute. |
| **MH-6** — the AGO gate rescales but rarely reorders | A relative/ordering statement about a monotone rescaling; asserts no significance. |
| **MH-11** — the D76/IRF1 route reproduces in bulk basal contrasts | Already downgraded to pattern-strength (**P**) by MH-19; the artifact uses no learned quantity. No *new* death — but do not re-inflate it either. |
| **MH-15** — IRF1 route visible in the archetype panel | Untouched by this arc — but **independently rotten**: IRF1 is one of the six genes that *define* the archetype label (`visibility_archetypes.py:53`), so the gene-level half is **circular by construction**. Keep only the descriptive Basal-vs-Luminal clause, flagged as unadjusted abundance. |
| **MH-18** — the basal hub is not a curation artifact (miR-17-5p #1 under three weightings) | A statement about **regulatory load** (a degree/weight property of the curated universe), not about coupling or the learned β. |
| **MH-20** — miR-17~92/106b~25/106a~363 carry ~22% of basal repressive load | Computation is current and reproducible — but **the registry's numbers are stale**: re-derived from `family_load_share.tsv` the shares are **21.5% / 28.1% (TargetScan) / 8.9% / 3.9%**, over **702** ranked arms (not "~515"). Fix the prose. Also: this is a *definitional share of a defined metric*; load ≠ coupling. |

---

## 4. The load-bearing distinction: EXISTENCE vs IMPORTANCE

### 4.1 Where the biology actually enters the model

`assemble_gene(he_only=True)` builds each gene's design **from that gene's curated high-evidence regulators**.
**The biology enters through the INCLUSION CRITERION — which edges are in the design — not through the
learning** [MH-124 §3]. Decades of luciferase assays, westerns and proteomics chose the columns. The Gibbs
sampler then apportions covariance among columns it was handed.

That is why the two halves of the model have completely different standing:

| question | status | evidence |
|---|---|---|
| **EXISTENCE** — "is this a real regulatory edge?" | **VALIDATED** (set-level, observational, one line) | curated edges beat abundance-matched impossible ones: OR 1.41 (p=2.5e-13, survives composition); realized ρ −0.0353 vs −0.0227 (p=3.0e-04); PIP 23.0% vs 15.5% (p=1.1e-04) |
| **IMPORTANCE / ATTRIBUTION** — "which regulator matters for this gene?" | **NOT VALIDATED** | β names the literature-canonical regulator **at chance** (arm-level 42.9%, p=0.38; graded rank 0.518 = chance, p=0.66). **Mere abundance names it far better (0.240, p<1e-4).** Decoys absorb **36% of |β|** and reproduce **~75%** of the field. |

**Nuance that must not be lost:** the model **ranks** the canonical family above chance with the doctrine's own
estimators — Shapley identity rank 0.319 (**p=0.0076**) and dose 0.258 (**p=0.0003**) [MH-124 §4b, user-caught
correction to the original "at chance" verdict]. **It ranks the right family high; it cannot NAME it.** β
specifically is the *worst* functional for this — it is also the functional most willing to hand credit to a
decoy [MH-128 §2a].

### 4.2 The one-sentence version

**The edge set is curated-real; the weights are untrustworthy as attribution; and the model's out-of-sample
advantage over abundance is ~40% curation and ~60% fitting.**

Use β as a **de-confounded, jointly-apportioned, transferable coupling estimator** — it is the only aggregator
in this subproject that predicts CPTAC protein out-of-sample and survives composition. Do **not** use it to say
who drives a gene.

---

## 5. What must be re-run, in dependency order

**Already done in this arc (verified by file mtime, this session):**
0. `readouts_{edges,genes}` + `readouts_arm_{edges,genes}` regenerated 2026-07-13 15:56/15:59 after the
   `spike_slab._rtnorm_pos` support bug (3.15% of persisted β were negative — impossible for a repression-only
   prior). Current β has 0.0% negatives.
1. `eval/coupling_grid.py` fixed (true RPM pool; `fillna(0)`; per-gene member sets; z-scored design + variance
   floor) and `coupling_grid_{summary,genes,by_composition}.tsv` regenerated 2026-07-13 18:22.
   ⚠ **`eval/coupling_grid.py` is currently UNTRACKED in git — commit it.**

**Still owed, in order:**
2. **Null recalibration for every FDR count in the Aim-1 block.** No count of "N/50 programs" or "N/14 routes"
   may be re-asserted against the theoretical t-null. The replacement null is a **degree-preserving edge
   shuffle** (set level) or an **abundance-matched site-free decoy** (edge level). Affects
   `hallmark_interaction/`, `robustness/aim1_proliferation/`, `robustness/targetscan_orphan/`.
3. **Re-derive, or retire, everything built on `pressure_build.compute_gene_pressure`.** MH-115 retired it (no
   better than abundance, zero out-of-sample signal) and T2 re-confirms it (CPTAC protein p=0.11, mRNA p=0.83).
   Consumers: `hallmark_interaction`, `subtype_contrasts`, `pam50_gene_resolution`,
   `visibility_archetype_contrasts`, `robustness/*`, `targetscan_orphan`. Replacement estimator = the learned β
   **with a composition block**.
4. **Add a composition block wherever there is none.** `config.CONFOUNDER_NUMERIC` is only `(CPE, HRD)` — no
   deconvolution fractions, no proliferation. That is the covariate set behind MH-1/2/7/8/12/13/16/21/22.
5. **Propagate the 2026-07-11 arm-name `.N` universe fix** (~18 HE arms / 300 edges recovered, incl. canonical
   miR-101/124/126-3p). **Every output dated 2026-06-26 predates it** — that is the entire Aim-1 block above.
6. **Registry text fixes** that need no compute: MH-20's shares (21.5/28.1/8.9/3.9, 702 arms); MH-7's covariate
   set (CPE+HRD+batch, not CPE+PAM50); MH-19's counts (4/14, VIM out).

---

## 6. What is genuinely unknown — *underpowered ≠ negative*

These are **open questions, not refutations.** Presenting them as refutations would be as wrong as presenting
them as results.

1. **Within-family attribution against the literature (n = 16).** The canonical-MTI list has **sixteen** usable
   entries. Aneuploidy-controlled, the model hits **4/16 (p=0.18)**. That is a test with almost no power. It
   licenses "we cannot demonstrate attribution", **not** "attribution is wrong". The constraint is **the
   literature set, not the estimator** [MH-124 §5b.2, MH-126 §1].
2. **The exogenous CN reduced form.** Only **~157 HE arms** carry a usable instrument. Clustered at the correct
   unit, p = 0.115 (significant in 2/10 matching draws). But the point estimate points the right way and is
   *quantitatively consistent* with a real effect (first stage 0.184 × implied causal −0.052 = −0.0096 =
   exactly what is observed). **"The data cannot say" — not "the answer is no."** A better-powered exogenous
   instrument remains the single highest-value thing we could build, because it is the only lens the shared
   bulk confound is blind to.
3. **Gene-field vs edge-level CPTAC transfer** — MH-127 says the real field separates from a fake one
   (p=0.030); T1 says it does not at edge resolution after composition. **Unreconciled.**
4. **The protein-only real-vs-decoy β-transport gap** (Δ = −0.061/−0.070, p ≈ 0.04–0.048, **both CIs cross
   zero**). Mechanistically coherent — protein carries translational repression, which a site-free decoy cannot
   produce, and the gap is absent at mRNA. **A lead, not a result. Do not headline it.**
5. **Whether Hallmark-program specificity exists at all.** Cognate-vs-foreign specificity is currently at
   *exact chance* with the heuristic — but it has **never been tested with the learned β + composition**. The
   50 pressure rows are near-collinear (median pairwise ρ = 0.807), so today we are testing one global miRNA
   axis fifty times. A properly-estimated version could still find program structure. It could also confirm
   there is none. **We do not know.**
6. **Whether the surviving MH-19 routes (TGFBR2←miR-17/20a, PTEN/CDKN1A←miR-106b) are causal.** They clear the
   hardest null we have built, which is an *observational* null. No exogenous evidence exists for them either.

---

## 7. Provenance

| tag | what it means |
|---|---|
| T1 / T2 / T3 | run in this session (MH-125); scripts in `scratchpad/{decoy_cptac*,step1_repro,step2_corrected,step3_mh115,vz_*}.py`, run with `PYTHONPATH=/sci/labs/michall/stavzok/APM .venv/bin/python3`. **No repo module was modified by those runs.** |
| MH-124 §x / MH-126 / MH-127 / MH-128 | measured and recorded in the named doc; each carries its own control and adversarial log |
| audit table (§3) | each row re-verified by re-running or re-reading that claim's own persisted output; the verifier's status wins where it disagreed with the first-pass auditor |
| grid table (§2.3, T2) | read directly this session from `mirna_hallmark/output/learned/coupling_grid_summary.tsv` (regenerated 2026-07-13 18:22, post-fix) |

Numbers not produced by a run in this session or read from a named persisted output **are not in this
document.**
