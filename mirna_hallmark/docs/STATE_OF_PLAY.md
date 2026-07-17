# State of play — `mirna_hallmark`

> **Goal:** the ONE doc you read to know where the project stands. Per axis: what **stands**, what is
> **dead**, what is **open** — each with its MH id and date, so you can go deeper without guessing.
> **What belongs here:** one short block per axis, verdicts only. **NOT** methods (→ `FORMULAS.md`,
> `LEARNED_MODEL_METHODS.md`), **NOT** the finding-by-finding record (→ `DISCOVERY_REGISTRY.md`),
> **NOT** run history (→ `ANALYSIS_RUN_LEDGER.md`), **NOT** forward plans (→ `PROGRAM_FORWARD_BOARD.md`).
> **Update trigger:** any arc that changes a verdict. If you write an `MH-` registry row that supersedes
> or retracts something, update the matching block here **in the same pass**.
> **Sync-partner:** `DISCOVERY_REGISTRY.md` (the record of record — this doc is its executive summary).

**Last updated: 2026-07-17** · covering the registry through **MH-143**.

---

## How to read this

**The `DISCOVERY_REGISTRY.md` is the source of truth.** It is append-only and current. Every other
doc — including the ones labelled "canonical" — may lag it. Where this doc and any other doc
disagree, the registry wins, then this doc, then everything else.

⚠ **The registry runs ahead of the docs.** MH-131…143 (the most consequential fortnight in the
project) exist **only as registry rows**. There is no `MH13x_*.md` for them, and that is now policy —
the per-finding-doc pattern is **banned** (`../CLAUDE.md` §3). Do not conclude from a doc's silence
that nothing happened.

---

## The one-paragraph version

The last two weeks were a **self-falsification campaign**, and it worked. Most of the June headlines
are dead. What survives is narrower, better-founded, and honestly bounded. The model is a good
**coupling estimator** and is **not** an attribution engine. Edge existence rests on **one
observational line** — both copy-number instruments are retracted. The curated edge set beats a
properly-built fake by **≈0.012–0.013 ρ** — a value two independent designs now agree on (MH-137/139),
and which *shrank* at every control improvement before converging there. The
single strongest result in the program is a **negative control that correctly reports zero**
(MH-136). Two whole channels — state and protein — were measured and cancelled, for the same
structural reason: information, not effort.

---

## Axis 1 — The learned model (what it *is*)

**STANDS.** One **dense learned-τ² non-negative Bayesian posterior per gene, two readouts** —
π≡1 dense → coupling/attribution/identifiability; evidence-π → discovery. Half-normal slab
(`θ_m ~ N⁺(0,ν²)`) encodes "miRNAs repress". **ν² is learned, not guessed** — that, not L2, is the
whole edge over the lasso (`nnridge_cv ≈ lasso`; OOF −0.168 vs −0.152, p=9e-16). **Adaptive lasso is
RETIRED to a baseline.** Canonical: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b (2026-07-13).

**Engine: Gibbs, not HMC** — measured, not aesthetic. Every live channel is Gaussian-conjugate;
Student-t is a scale mixture; the protein discordance link is *additive* ⇒ conjugate; HMC was only
warranted for the A3 field, which the A4 probe gated out (pooling is a wash, Δρ −0.0006). Gibbs is
also the better estimator (split-half ρ=0.822 vs bagged-NNLS 0.729). **Carve-out:** fit the *gauge*
with NNLS — Gibbs's heavy-tailed posterior SDs break the errors-in-variables correction.

**What it is FOR (MH-126b):** a **de-confounded, jointly-apportioned, transferable coupling
estimator**. **What it is NOT:** an attribution engine. *It must not be sold as an answer to "which
miRNA represses this gene".*

**OPEN / known limits:**
- **Posterior widths are ~25% too narrow** — bagged NNLS **1.37×**, Gibbs **1.29×**, Student-t ν=7
  **1.13×** (`CPTAC_PROTEIN_CHANNEL_PLAN.md` §"CORRECTED 2026-07-12").
  ⚠ **`LEARNED_MODEL_STATE_CHANNEL_PLAN.md` §11(a) still says "~40% (NNLS 28%, Gibbs 39%)" and "SBC
  promoted" — SUPERSEDED.** The earlier figure was measured on a biased gene subset (the
  scale-dependent sd-floor bug); **SBC is RETIRED as the wrong tool.**
  ⭐ **MH-92's `nu=7` substantially FIXES the calibration (0.77 → 0.89)** — so Student-t is
  "robustify-not-lever" *for coupling*, but it **is** the calibration fix. Turn it on for any power calc.
  ⚠ **Two things do NOT inherit the inflation:** **Shapley shares are RATIOS, so a common β inflation
  cancels entirely** (re-derived: MH-94's width is **1.8× too WIDE**, not too narrow) — and note
  `readouts.py` does not pass `nu`, so the shipped table is Gaussian.
- **Identifiability ceiling:** median posterior SD/|β| = **0.799**; only **27.8% of 5,117 units
  identified** (|z|>2). Not a fitting failure — same-seed arms share the binding site.
- **The n≈1000 ceiling:** every internal lever (CN channel, δ-pooling, isomiR refit, t-likelihood,
  cross-gene pooling) lands "immaterial at n≈1000".
- ✅ **E6 "self-partialling bug (live)" — NOT LIVE; the tag is RETRACTED (MH-142, 2026-07-17).**
  Verified: every call site scores **orphans**, which are disjoint from `he_agg` by construction —
  **0 of dossier's 6,744 candidate pairs are HE edges** (vs 5,938 HE edges). Nothing is partialled on
  a covariate containing itself in production. E6 **stands as a trap for CLASS-MATCHED comparisons**
  (its evidence column: *"class-matched estimator comparison"*) — push HE arms through this estimator,
  as MH-124's own experiment did, and you get a ~36× depressed HE rate **as an artifact of your design**.
  ⛔ **Do not "fix" it by removing `he_agg`** — that term is what makes an orphan's coupling *partial on
  what curation already captures*, i.e. the point of the lane. Guard: `eval/_e6_production_check.py`.
- The **shipped readouts do not include the CN channel** (`readouts.py` calls `_gibbs_posterior`
  without `channels=`) — the production table is mRNA-likelihood + evidence-prior only.

**FIXED, don't re-litigate:** the `_rtnorm_pos` **sampler bug** (inverse-CDF saturated for
`mu/s < −7.0345` ⇒ deterministic *negative* draws, impossible under a half-normal slab; contaminated
3.15% of β). Replaced with Robert (1995) rejection sampling, KS-verified. ⚠ It **overturns MH-119's
recorded cause** — the `share`=43.7 explosion was logged as "anti-repressive biology"; it was this bug.

---

## Axis 2 — Edge existence

**⛔ It rests on ONE observational line:** *curated edges anti-correlate more than abundance-matched
site-free ones.* The "three independent validations" headline is **FALSE** — the other two legs were
both observational and never independent of each other.

**Both CN instruments are RETRACTED:**
- `pi_causal = γ·b` where **b is an observational OLS slope** ⇒ a product-of-coefficients mediation
  estimator, not an IV. It re-derived the anti-correlation it was meant to validate. γ — the only
  exogenous factor — is **site-blind** (HE +0.19931 vs decoy +0.19844, p=0.20). The clean reduced
  form, arm-clustered: **p=0.115, n.s.** (MH-124r/126, 2026-07-13). The retraction is in the
  production code at `learned/instrument.py:1128–1148`.
- The within-arm two-way-FE replacement: its control class was **71% real binders**. Against a
  genuinely site-free control **τ=−0.0007, p=0.84**, and it fails the site-type efficacy ladder
  (8mer/7mer-m8/7mer-A1 not monotone, p=0.26) — copy number cannot see a site type, so genuine
  site-mediated repression *must* follow it (MH-133, 2026-07-16). **Refutations at n=216k–235k pairs,
  not power failures.**

**⚠ Also bounding this axis:** the per-edge null is **3–4× too narrow**. Site-free (cannot-repress)
pairs pass "FDR q<0.05" at **25–35%**; under an honest empirical null **0.0%** of curated HE edges
survive per-edge (was a nominal 31.3%). **This does not say the edges are fake** — a typical HE edge
is ρ≈−0.15 against σ₀≈0.09. It says the signal is a **set-level distributional shift, not per-edge
significance** (MH-123/124). Any FDR count in any doc predating this rests on the uncalibrated null.

**OPEN:** restoring an exogenous existence validation is, by both docs' own assessment, **the single
highest-value open item in the program.**

---

## Axis 3 — Attribution / Shapley

**STANDS: the model RANKS the canonical family high; it does NOT NAME it.**

Cite **MH-126 Test 1** (n=32, 2026-07-13), not MH-124 §4b (n=21) — same result, better powered,
correctly clustered. Mean normalised rank (0 = top, 0.5 = chance):

| readout | rank | verdict |
|---|---|---|
| **abundance** | **0.240** | p<1e-4 ✅ |
| **`shapley_identity`** | **0.370** | p=0.018 ✅ (cluster p=0.023–0.047) |
| **β** | **0.518** | p=0.66 — **chance** |

The 32 genes hold only **13 distinct families**, and abundance is family-constant ⇒ the independent
unit is the **family**. Cluster-bootstrapped, shapley survives. **Argmax is at chance at every k.**

**MH-124 §4b and MH-126 do NOT conflict** — both say β is at chance; §4b's positive claim was about
`shapley_identity` and dose, which MH-126 reproduces. §4b survives as the **estimand argument**: the
literature's "miR-X represses gene Y" is a **level** claim (knockdown/overexpression); β is a
**covariance** functional. Asking β to recover it is an estimand mismatch, not a model failure.

**WITHDRAWN:** shapley's **evidence gradient** — 85% is family fame (family-preserving null p=0.297).
Only β's gradient survives, barely (b≈+0.031, p=0.074–0.097 under family FE). Shapley wins the
**rank** test and loses the **gradient** test; β vice versa. Dissociation, not inconsistency.

**The w-circularity gate PASSES, bit-identically:** β is unchanged (max|Δβ| = **0.0**, 0/48 genes)
under shuffled *and* constant `w`, with the positive control firing (`w_fam` moved 42/48). Mechanism,
not luck: `X_fam` is pure abundance and the dense Gibbs is called with `π = ones(p)`.
⚠ **But `pip_discovery` and `prior_pi` ARE w-contaminated by construction** — any "do canonical
regulators score higher?" test on those columns is **circular**.

✅ **FIXED 2026-07-17 (MH-140) — the table now has a TRUE identity column, and it disagrees with
magnitude in 1 of every 4 multi-family genes.** `readouts.py` emits **`identity`** = `shapley_identity`
(Shapley/LMG on **R²(Xz·m, yr)** — a NON-additive value function, so collinear overlap is genuinely split),
with **m = bagged NNLS**. The old `share` (= β_f/Σβ) is renamed **`beta_frac`** and demoted to MAGNITUDE:
for the ADDITIVE aggregate Shapley(f) ≡ β_f, so it split nothing (MH-122), and `top_family_identity` was
`argmax β` in disguise — **identical to `top_family_magnitude` in 99.35% of genes**.
**Now: WHO ≠ HOW-MUCH.** Agreement **85.3%** genome-wide, **75.2% among the 819 multi-family genes**.
**The mechanism, at genome scale:** `m_nnls` is exactly 0 in **1,623/5,117** families (31.7%); `beta` is
exactly 0 in **0/5,117** — the half-normal slab cannot zero an un-informed family, NNLS can.
⭐ **Identity can also say UNDEFINED:** 213/1,549 genes (13.8%) have NaN identity, and in **100%** of them
NNLS zeroed every family ⇒ no coupling to distribute. `beta` structurally cannot make that statement.
**Control: purely additive** — `beta`/`beta_sd`/`z`/`pip_dense`/`pip_discovery` bit-identical (max|Δ|=0.0).
⚡ Three independent wins, all controlled: `shapley_identity` **8.2–9.7×** (output-identical); **1 of 4 Gibbs
chains per gene was DEAD** (the evidence-π chain ran on both C blocks; `dec["pip_disc"]` was never read) —
gated, output-identical; and **N_ITER/BURN 2000/700 → 1000/350** (MH-141), licensed against the sampler's own
reseed jitter and verified NOISE-not-BIAS. **Run 5.9 → 1.6 min (3.7×) while emitting MORE.** ⭐ `identity` is
**bit-identical** under the chain cut — it rides the NNLS weights, not the Gibbs mean.

**✅ SETTLED (MH-138): "§5 bagged NNLS vs §6b dense posterior" was never a conflict — one word, two
jobs.** MAGNITUDE = the Gibbs posterior mean (`beta`) ⇒ §6b. IDENTITY = `shapley_identity`, which
needs **fixed weights** from `canonical_M` = **bagged NNLS** ⇒ §5, for that job. And NNLS is
*required* there, not vestigial: the half-normal slab has a **strictly positive mean**, so an
un-informed family **cannot be zeroed** — measured, `beta == 0` in **0/5,117**, 100% positive —
whereas NNLS returns exact zeros. The built `bayes_shapley_identity` already does the right thing:
**NNLS fixes support+signs, Gibbs draws supply the width.** Rule of thumb elsewhere: ***bagged NNLS
for the GAUGE, Gibbs for the MODEL.***

**OPEN:** **scale the literature set** — this, not the estimator, gates attribution. n went 16→21→32;
the rank test became decisive while argmax stayed at chance. Needs a versioned, auditable ground-truth
pull. Until then: *"β does not attribute" is MEASURED; "the CN instrument cannot attribute" is
UNDERPOWERED, not refuted.*

---

## Axis 4 — Decoy validation (does curation beat a matched fake?)

**MH-137 is the row of record, ✅ CORROBORATED by MH-139 (both 2026-07-17).**

**Core gap −0.0119 (p=3.4e-05); deconv −0.0090; retention 0.76; 1,322 genes** (dose-matched by a hard
|Δdose|<1 caliper). **MH-139 re-measured it the other way** — caliper OFF, all 1,349 genes kept, dose
residual removed post hoc by regressing gap on Δdose with **b re-derived on this decoy (+0.01122,
95% CI [+0.0018,+0.0198])**: **gap at Δdose=0 = −0.0129 [−0.0205, −0.0052]**. **−0.0119 lies inside
that CI** ⇒ the caliper's 14.4% universe restriction (which drops precisely the HIGH-DOSE regulators,
1,164 vs 31 RPM median, p=2e-179) **cost nothing material.** Two designs, one answer: **≈−0.012 to −0.013.**

⚠ **The trajectory runs DOWN as controls improve** — this is the single most important thing to
internalise on this axis:

> **−0.045** (v3, loading-handicapped) → **−0.0306** (6mer-contaminated) → **−0.0147** (dose-weighted)
> → **−0.0119** (evidence-excluded + calipered)

Three defects fixed in MH-137: a **126× evidence hole** (only 2% of miRTarBase was excluded, so
CLIP/CLASH-supported *real* edges served as "decoys"; now 1,790,439 pairs excluded); the **6mer gate**
(Poisson-gated, not presence — a single 6mer is 24% chance noise at the median UTR); the **dose
caliper** (Δdose −3.05 → −0.105).

**The honest statement is a BOUND, not a win:** *real beats fake by ≤0.012 ρ, and every control
improvement so far has shrunk it.* Whether the next fix (Manakov chimeric eCLIP + POSTAR3, still
uncovered) takes it to zero is **genuinely open**.

**⭐ STANDS — MH-136, the arc's first real positive control, and its strongest result:** using the
genome-wide 6mer map, **187 genes whose curated edges are ALL seedless show gap +0.0006 — exactly
zero**; **810 all-seeded genes show −0.0325**. A priori, sequence-only, non-circular, and the model
behaves correctly on **both** sides. *No site → no repression → no coupling → no gap.*
It also closes the **21% positive-coupling mystery** (MH-132): controlled for regulator count,
`frac_seedless` predicts positive coupling (+0.042, p=2.2e-04).

**The internal null is back:** the **1-family stratum** (where β ≡ uniform exactly) gives
**−0.0007, p=0.15** on the clean control. MH-135's retraction of it was an artifact of the broken decoy.

**STANDS: the advantage over *unfitted abundance* is a FITTING effect, not curation** (MH-115/127).
A fitted matched fake also beats real unfitted abundance. **Rule: an abundance baseline is not a
control. Benchmark any aggregator against a fitted matched decoy.**

**Domain of validity — half-confirmed:**
- **Width axis (`n_fam ≥ 3`): CONFIRMED** and sharpened — the gap concentrates at 3–4 families
  (−0.0314, p=0.0006). ⚠ Not monotone under MH-137: it **peaks at 3–4 and falls at 5+**.
- **Evidence axis (`w_max > median`): NEVER tested on a clean control.** Its only support fails the
  arc's own both-fake-sets rule (q=0.006 FAKE1, q=0.20 FAKE2).
- ⇒ Carry MH-130's "27% domain" as a **live hypothesis**, not an established partition.

**The gene lens (`gene_atlas.tsv` STANDS; the gap columns need a re-run):**
- **17.6% of genes are NOT MEASURABLE at all** (OOF R² ceiling of an unrestricted OLS ≤ 0); **51.7%
  have a ceiling ≤ 2%**; 27.1% not measurable under composition C. The ceiling is **design capacity,
  not biology** (spearman(ceiling, n_fam) = +0.564).
- **48.1% of genes have ONE seed family, where β ≡ uniform exactly** ⇒ *"does the learning add
  anything?"* is **UNDEFINED for half the universe — a non-question, not a null.**

---

## Axis 5 — CPTAC / protein

**⛔ DEAD — do not rebuild:** `βᵗ` (the translational-repression latent, the plan's centrepiece) is
**not supported at n=101**, and all prior positive evidence was a **LEAK** (mediator slope `a` fitted
on all samples, then OOF-ed; fit `a` inside the fold and it vanishes). Four leak-free framings agree;
decisively the aggregate (1 df, not 60): **BH q<0.10 in 1/17 genes — BCL2 alone; PTEN d=−0.006,
p=0.82** (MH-103).

**⚠ MH-34 is NOT overturned** — 13/17 stay repression-directed (binomial p≈0.025). The falsification
is of the **modelling object**, not the biology. The translational residual is ~1% of genes; 17 genes
predicts ~0.2 hits; we got 1. **There was never a βᵗ field. Only n fixes this.**

**STANDS — protein can never be a coupling lever:** it carries **4–6%** of the mRNA channel's Fisher
information about β, and **≤7.6% at any `a_g`** (pre-registered ceiling; a 5× correction to the
number moved the conclusion zero). Converges with the state channel's independent 0.6%.

**STANDS — MH-114, the best-designed test in the program:** ~**57%** of the protein validation is
**compartment arithmetic**. Degree-preserving shuffled null, stratified by compartment orientation:
**shuffled non-edges reproduce the real gradient exactly (−0.1290 vs −0.1272).** Pairs that *cannot*
repress produce the same bulk anti-correlation. The true compartment-blind effect is **≈−0.011**,
surviving adjustment at **≈−0.0065** (p=2.0e-03). Damage: gene FDR-neg **34→3**, hallmark **9→0**,
EMT **−0.330 → +0.162 (sign flip)**, orphan candidates **594→23**.

⚠ **The protein transfer is NOT decoy-controlled** — a site-free fitted fake reaches **99%** of it
(MH-130d). It is real arithmetic, but **not evidence the curated edge set is right.** The only
decoy-controlled transfer in the arc is the **CPTAC mRNA** cell inside class A.

⚠ **MH-114's "34→3" was measured on the retired heuristic** (`compute_gene_pressure`), not the
learned posterior (MH-115) — and its FDR counts rest on the uncalibrated null (MH-123). The
*conclusion* is unharmed; **do not reprint the counts.**

**The confounder-ARCHITECTURE principle (MH-111) — the most generalisable thing here:** proliferation
and host localise to an **ARM** (LOCUS properties); composition localises to the **FAMILY** (an
EXPRESSION property, because same-seed members sit at different loci but are co-expressed).
**The unit at which a confound acts predicts which tool can fix it.**

---

## Axis 6 — Progression / state

**Two different objects. Do not conflate them.**

**⛔ The state CHANNEL is CLOSED — measured and cancelled** (MH-102d, 2026-07-12).
`learned/channel_state.py` **was never built** (verified absent). τ — the acquired-vs-constitutive
payload, the axis's entire deliverable — is **indistinguishable from zero** (GTEx τ=0.0009, info
**0.6%**; NAT 0.7%; CPTAC 0.7%). Structural: channel precision ∝ **a²**, and a≈0.11 ⇒ **even
unlimited GTEx donors add ~1%**. *Composition attenuation does not merely bias the healthy β — it
destroys its information content.* The additive `M_T = M_H + Δ` form **would have faked the
headline** (Δ>0 on nearly every edge, manufactured out of attenuation). **Do not rebuild it.**

**✅ The abundance-level cross-state LANDSCAPE is built and run** — GTEx-healthy → TCGA-NAT → tumor,
721 arms / 5,108 edges. Verified against the persisted 2026-06-27 artifact: **acquired_realized 640
(129 arms) · lost 209 (89 arms) · stable 2,141 · acquired_unrealized 1,775**. ⚠ Its **per-edge FDR
class labels** rest on the null MH-124 measured as 3–4× too narrow; `LANDSCAPE_REPORT.md` carries no
banner about this. ⚠ NAT is structurally underpowered (n≈104): `nat_decoupled` = **1 edge**.

**STANDS — the panel question is SETTLED, zero new CIBERSORTx runs needed:** β is **ρ=0.94** across
entirely different atlases (Wu-9 vs HBCA-11) and ρ=0.99 across reference resamples, but only ρ=0.81
vs **no C at all** (mean |β| halves). ⇒ ***which* composition control you use moves β ~10× less than
*whether* you control composition.** The "reference reconciliation blocker" is **RETIRED**.

**What survives the cancellation:** `learned/gauge.py` and `learned/confounders.py` — *they are what
killed the axis, which is what falsifiable infrastructure is for.* **Standing rule:** any future
cross-cohort channel runs `gauge.calibrate()` → read `info_ratio` + τ **BEFORE** it is built.

---

## Axis 7 — Buffa / independent cohort / outcome

**⚠ "METABRIC" in MH-73/74/75/76 means Buffa, n=207 — verified in code**
(`pressure_prognostic_signature.py:166`, `pressure_prognostic_gene_centric.py:97` both import
`load_buffa`). Buffa's miRNA arm **is** the METABRIC invasive miRNA subset (GSE22216 miRNA +
GSE22219 mRNA, 207 paired tumours, **DRFS only, 77 events**, no TCGA overlap). ⇒ the
+0.056 / +0.060 / +0.028 / +0.019 are **four views of ONE 207-patient result**, not four validations.

**STANDS (the well-earned negatives):**
- **Pressure magnitude is exhaustively non-prognostic** — 0 recurrence FDR across 3,368 features.
- **Functional > magnitude** — realized +0.060 vs received-pressure −0.016 (a *within-cohort*
  contrast that cohort-specificity cannot explain).
- **miR-210 → DRFS reproduces Buffa 2011** — a working **positive control**, which is what makes the
  pressure nulls credible rather than merely underpowered.

**⛔ DEAD — the "TCGA is sparse so it cannot arbitrate" defense.** MH-76 ran the **frozen-panel**
test (Buffa-trained, frozen, scored on TCGA), which removes overfitting as an explanation entirely:
**DFI +0.002, OS +0.006; panel-alone C 0.48–0.53 ≈ random.** MH-73/74 still assert the superseded
framing and need annotating.

**⛔ SETTLED BY MEASUREMENT 2026-07-17 — "108 orphans triple-validated" IS 3.** Re-ran
`eval/buffa_validation` against the composition-adjusted screen (594 → 23 candidates since 2026-07-12,
MH-114): **`triple_validated` 108 → 3; ECM/collagen among them 30 → 0** (all rows 64 → 0). The three
survivors are all miRGeneDB **guide** arms — **miR-195-5p→PSMD7 · miR-181b-5p→IRS2 · miR-22-3p→STK39**
— and only **1** is neg-sig in Buffa. **The miR-29→collagen/ECM headline is gone.** The pre-registered
expectation ("single digits, possibly zero ECM") was correct.
⚠ **The module had been DEAD, not merely stale.** `GEO_DIR` was a `__file__`-relative hop
(`parent.parent/"data"`), correct when the module sat at the subproject top level — but it was moved
into `eval/`, which added a directory level and pointed it at a non-existent `mirna_hallmark/data/`.
Every input raised `FileNotFoundError`; the GEO cache was at the **repo root** all along. Fixed to
`C.REPO_ROOT`. **That is why this lane was never re-run after MH-114: it could not run.** Ripple
checked — every other `__file__`-relative path in the subproject is correct (`learned/readouts.py`,
`learned/cptac_data.py`, `presentation/make_figures.py` use `parents[2]`; `config.py` correctly uses
`parent.parent` at the top level).
✅ **CONTROL — the HE-edge legs reproduce EXACTLY** (concordance **+0.319** · sign-concordance
**0.593** · TCGA neg-sig same-sign **0.673** · neg-sig-in-Buffa **0.128**), so only the orphan lane
moved. Pre-fix outputs preserved at `output/buffa_validation_PRE_MH114/`.

**⚠ Buffa is NOT in the learned arc at all** — one docstring line in `learned/eval/__init__.py`
(marked `STUB`); there is no `ood_cohort.py`. **β has never touched Buffa.** This is the most
interesting gap in the program: MH-104 measured the model's loss as **~80% at the COHORT boundary
and ~0% at the LAYER boundary**, and Buffa is the only independent-patient cohort in the repo with
both miRNA and mRNA.

**⚠ "METABRIC-full (EGA-pending)" is doing unearned work across five docs.** There is **no accession,
no DAR id, no submission date, no owner** anywhere in the repo. The evidence supports *"never applied
for"*, not *"pending"*. It would also be a **power** upgrade, not an **independence** one — same
platform, same consortium, Buffa is a subset of it.

---

## The four things worth doing next

1. **`β(TCGA) → Buffa mRNA`** — the model's only untested boundary is the cohort boundary; this is the
   only cohort that can test it; the inputs are already on disk. ⚠ **Pre-register MH-114's
   compartment-orientation stratification** — both cohorts are bulk breast and share the CAF confound,
   so a clean replication proves nothing on its own.
2. **Re-run `buffa_validation` (~90 s) and annotate MH-38 / MH-55 / MH-73 / MH-74.** Correctness debt:
   a dead claim is live in the talk.
3. **Re-run the competence map's gap columns through `eval/decoy_bench.py`** (2:22 full universe).
   Every real-vs-fake number in `gene_competence_map.tsv` sits on a discredited control; the atlas
   half is fine.
4. **Restore an exogenous existence validation** — both `CN_INSTRUMENT.md` and MH-126 name this as the
   highest-value open item, and nothing else changes the program's standing as much.

---

## Doc traps — read this before trusting a doc

| doc | trap |
|---|---|
| `LEARNED_MODEL_ESTIMATOR_MAP.md` | Billed as the "canonical map", but self-dates to **2026-07-06** and its main table is **pre-convergence** (lasso-as-primary, "Bayes = uncertainty layer"). Only its two banner UPDATEs and the pooled-HE policy are current. |
| `MH127_*` / `MH128_*` / `MH130_*` | Read as current; their central magnitudes are **superseded by MH-137**. MH-130's *atlas* half stands; its *gap* half does not. |
| `MH124_*` §5 / §5b | **VOID** (`pi_causal` not exogenous). **Do not cite.** §4b is untouched and stands. |
| `CPTAC_PROTEIN_CHANNEL_PLAN.md` | §6 still carries **VOID** PAM50-contaminated numbers (1.7% / ≤8.8%), contradicting §1 on the same page. §1's correction strings are self-nested and unreadable. |
| `handoffs/HANDOFF_PROTEIN_AND_COMPOSITION.md` | Dated 2026-07-12; **12+ rows have landed since**. Superseded on the estimator (MH-115), the decoy control (MH-127/130/132), and the null calibration (MH-123/124). |
| `PROGRAM_FORWARD_BOARD.md`, `WHATS_NEXT.md`, `LEARNED_MODEL_WHATS_NEXT.md` | **Three competing forward docs, all predating MH-133…137.** All stale on the CN instrument. Consolidation pending. |
| `LANDSCAPE_REPORT.md` | Its per-edge FDR class counts rest on the **uncalibrated null** (MH-123/124). No banner yet. |
| `presentation/` F25 + `talk.pptx` | `talk.qmd` prose was fixed 2026-07-12 and is honest; **the F25 figure beside it is from 2026-06-24** and shows the old numbers. Deck needs re-rendering. |
