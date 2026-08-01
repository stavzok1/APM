# State of play — `mirna_hallmark`

> **Goal:** the ONE doc you read to know where the project stands. Per axis: what **stands**, what is
> **dead**, what is **open** — each with its MH id and date, so you can go deeper without guessing.
> **What belongs here:** one short block per axis, verdicts only. **NOT** methods (→ `FORMULAS.md`,
> `LEARNED_MODEL_METHODS.md`), **NOT** the finding-by-finding record (→ `DISCOVERY_REGISTRY.md`),
> **NOT** run history (→ `ANALYSIS_RUN_LEDGER.md`), **NOT** forward plans (→ `PROGRAM_FORWARD_BOARD.md`).
> **Update trigger:** any arc that changes a verdict. If you write an `MH-` registry row that supersedes
> or retracts something, update the matching block here **in the same pass**.
> **Sync-partner:** `DISCOVERY_REGISTRY.md` (the record of record — this doc is its executive summary).

**Last updated: 2026-07-18** · covering the registry through **MH-159**.

---

## How to read this

**The `DISCOVERY_REGISTRY.md` is the source of truth.** It is append-only and current. Every other
doc — including the ones labelled "canonical" — may lag it. Where this doc and any other doc
disagree, the registry wins, then this doc, then everything else.

⚠ **The registry runs ahead of the docs.** MH-131…153 (the most consequential fortnight in the
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
  site-mediated repression *must* follow it (MH-133, 2026-07-16). ⚠ **LANGUAGE CORRECTED (MH-157,
  2026-07-18): these are UNDERPOWERED / can't-tell, NOT "refutations."** The ladder's MDE@80% (0.00442)
  exceeds the whole 8mer effect (0.00174); no bin is distinguishable from 0. The "τ=−0.0007, p=0.84" was
  the F>10-**gated** spec — F-weighted on the same control gives p=0.057. **The two named instruments stay
  dead regardless** (their retraction is unaffected); only the framing is downgraded.

**⛔ AND A THIRD DESIGN — guide/passenger within-hairpin — IS ALSO FORECLOSED (MH-157, 2026-07-18).** The
one confound-immune-*looking* CN design (two arms of a hairpin share the identical locus CN ⇒ ploidy &
co-amplification difference out; contrast identified only by which seed targets which gene) does **not**
revive the instrument. It is (a) **conceptually foreclosed** — CN cannot separate hairpin arms (`CN_INSTRUMENT
§3`), so the contrast is a *gene-set* comparison, not edge identification; and (b) **empirically an abundance
confound** — it vanishes when arms are comparable-abundance (|ΔabundE|<1 gap +0.004, p=0.65) and is ~half
co-cluster leak (guide targets co-targeted by a polycistron neighbour 43% vs passenger 16%); the clean
guide-own residual (−0.006, p=0.12) fails the site-efficacy ladder (context++ slope p=0.53). **Edge existence
still has NO surviving exogenous validation.**

**⚠ Also bounding this axis:** the per-edge null is **3–4× too narrow**. Site-free (cannot-repress)
pairs pass "FDR q<0.05" at **25–35%**; under an honest empirical null **0.0%** of curated HE edges
survive per-edge (was a nominal 31.3%). **This does not say the edges are fake** — a typical HE edge
is ρ≈−0.15 against σ₀≈0.09. It says the signal is a **set-level distributional shift, not per-edge
significance** (MH-123/124). Any FDR count in any doc predating this rests on the uncalibrated null.
**⛔ And a permutation null does NOT rescue it (MH-154):** `coupling_permutation`'s 2000 Freedman–Lane
draws were measured to reproduce the **theoretical** null (σ₀=0.0309 vs 0.0311, R²=0.9997) — **permuting
destroys the very structure that inflates the real null.** Nor does **BY**, which corrects *dependence*,
never *scale*: MH-45's BY count (**33.3%**) sits inside the **35.3%** rate for *impossible* edges.
**The calibration standard is `eval/site_free_null.py`** (rescued from `/tmp`, reproduces MH-123 exactly).
It is **conservative** — the site-free class holds some real non-canonical targets.

**⭐ THE DISCOVERY LANE, RESOLVED (MH-155/156).** Per-edge and per-family discovery in bulk are **EMPTY**
under the honest FDR (empirical heavy-tailed null → Simes within family → BH across families; `q_family_emp<0.05
= 0` everywhere). BY was demoted to a worst-case sensitivity bound — it corrected the wrong defect (dependence,
not the heavy tail; it even kept 2 tail false-positives the empirical FDR caught). **The defensible claim is
CONVERGENT-EVIDENCE, not per-edge:** coupling co-varies with independent biochemistry, abundance-controlled, on
multiple axes — **site-level chimeric overlap partial ρ=−0.096 (p=2.7e−15), scanMiR K_D −0.091 (p=7e−14),
site-confidence rung −0.050 (p=3.9e−5); site COUNT n.s.** (quality not quantity) — and the site–duplex
coincidence climbs the site-confidence ladder monotonically (9%→27%). Deliverable = **A1∩site-level-chimeric,
157 edges** (miR-18a/miR-17~92 cluster). Effects are SMALL by construction (bulk is weak per-edge); the strength
is the convergence. **⚠ NOT a per-edge licence.** Home: `learned/discovery.py`, `method_dev/site_ladder/`.

**OPEN:** restoring an exogenous existence validation is, by both docs' own assessment, **the single
highest-value open item in the program** — and it is now **TESTED, not aspirational (MH-159, 2026-07-18).**
**⛔ GERMLINE cis-eQTL MR built + run + refuted.** Germline genotype is genuinely exogenous (randomized at
meiosis; no reverse causation / tumour-state confound), so a germline cis-eQTL instruments a miRNA's dose.
Built at n=1,075 (near-complete overlap): SNP6 birdseed → hg38 cis-SNP matrix (1,075×128k). But the
instruments are **honestly WEAK** (cross-fitted OOF: only **10 arms F>10**; the single-SNP "38" was
winner's-curse inflation) and **NON-SPECIFIC** — the reduced form is repression-directed for the miRNA's
curated targets **no more than for random non-target genes** (curated +0.0027 vs control +0.0049, paired
p=0.66; still null with the strong multi-SNP instruments, p=0.42). ⇒ **germline fixes EXOGENEITY but not
EXCLUSION** (cis-SNPs still hit host genes/clusters/paralogs) and the miRNA cis-heritability is small; SNP
aggregation doesn't rescue it. **Edge existence still rests on ONE observational line.** ⬜ Only untested
extension: full genotype imputation (unlikely to change it — low OOF-F + specificity null). Home:
`method_dev/cn_guide_passenger/GERMLINE_FINDINGS.md`.
For the discovery lane specifically: correlation-matched/scrambled-seed
family null · unify the two lanes by family size · subtype-stratified (PAM50) run.

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

**Domain of validity — ⛔ THE "27% DOMAIN" IS RETIRED AS STATED (MH-147, 2026-07-17).** Re-tested on the
canonical decoy (1,349 genes):
- ✅ **Width axis (`n_fam ≥ 3`) HOLDS** — n=519, gap **−0.0264, p=4.2e-06**.
- ⛔ **Evidence axis (`w_max > median`) ADDS NOTHING GIVEN WIDTH** — within `n_fam ≥ 3`, high (−0.0298)
  vs low (−0.0130): **MWU p=0.293**. It separates alone only because it correlates with width.
  **The rule reduces to `n_fam ≥ 3`.**
- ⛔ **The complement is NOT zero** — MH-130's load-bearing "outside: +0.0002, p=0.72" becomes
  **−0.0114, p=6.3e-03 (n=936)**. The model works outside its claimed domain, just weaker.
- ⇒ **There is no partition. Competence is a GRADIENT in design width** — consistent with
  `spearman(ceiling, n_fam)=+0.551`. **The defensible claim:** the gap SCALES WITH REGULATOR-DESIGN
  WIDTH, is ~2.3× larger at `n_fam ≥ 3` than in the complement, and is **nowhere exactly zero**.

**The gene lens — ✅ REPRODUCIBLE 2026-07-17 (MH-144), and one headline RETIRED:**
- ⛔ **"17.6% NOT MEASURABLE" is FRAGILE — do not cite it.** It is a threshold count evaluated exactly
  where the distribution **piles up**: 39.8% of genes sit within ±0.01 of zero ceiling. A faithful
  re-implementation gives **25.5% on MH-130's own genes** while the per-gene ceilings agree at
  **corr 0.9956** — a systematic shift of **−0.00136** (0.06% of the R² scale) reclassifies ~8% of the
  universe. **The ROBUST statement: `ceiling ≤ 0.02` → ~52–53%**, which reproduces across implementations.
  ⚠ Same failure mode as MH-138's per-gene median (see Axis 3) — **twice in one day**.
- ✅ **What reproduces** (current 1,549-gene universe): ONE seed family **47.1%** · `A_COMPETENT`
  **27.6%** · spearman(ceiling, n_fam) **+0.551** ⇒ the lens's STRUCTURE is sound.
- ⚠ **MH-130's universe was stale + undocumented** — 1,421 genes, a strict subset of today's 1,549; the
  150 it never saw are 96.9% single-family. It almost certainly built on the decoy-eligible set (cf.
  MH-127's 1,436), which its doc never records.
- The ceiling is **design capacity, not biology** — target detection fails only 7.9% of genes;
  **regulator-design width is the binding constraint.**
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

**✅ NEW — the WITHIN-PATIENT PAIRED axis is a THIRD object (MH-158, 2026-07-18), and it works where the
cross-cohort channel couldn't:** differencing tumour − **own-NAT** (n=103 matched pairs) removes the patient
baseline, bypassing the composition attenuation that killed the channel. Engine `learned/realization.py`
(M = `M_complement`, fit on the 975 non-matched tumours — leak-free). **DOSE/potential:** the acquired dose is
~46% patient-specific (`own_specific_frac` median **0.46**) ⇒ own-NAT ≠ cohort/subtype anchor; within-gene
share shifts recover canonical grip (miR-21→PTEN). **COUPLING/realized:** a WEAK set-level signal — real edges
beat matched site-free decoys by **~0.017 ρ in the OPPOSITE-compartment stratum after composition-Δ**
(gene-clustered p=2.7e-4; SAME n.s.). Own-NAT gives a stronger raw within-patient ρ than cohort (p=2.8e-6) **but a decoy control DEMOTES that advantage to SHARED-baseline** (decoys show the same own−cohort gap, MWU p=0.92) — the repression-specific signal is the real-vs-decoy gap alone, reference-independent. ⚠ Tiny
(~1/6 single-gene SE), **NOT per-edge**, single cohort; **NAT≠healthy** ⇒ final malignant step only.
rigor-auditor **A** (CONCERNS→re-scoped: null-referenced gap not −0.050/64%; retention median 0.61 not 38%).
✅ **Phase 2/3 DONE (2026-07-18, rigor-auditor PASS, all 9 predictions pre-registered+confirmed).** Res-3 family:
nonlinear pooling adds no coupling (`family_agg` −0.046 ≈ arm-level −0.050). **Owner-convergence** (Res-4): the
realization-Shapley owner tracks the static-identity owner above a within-gene random-owner null (z=6.1) **but the
coincidence is partly MECHANICAL (shared M vector), realization~dose is n.s. (z=1.6), and realization~static is not
distinguishable from static~dose ⇒ the realization-SPECIFIC ownership content is NULL.** Retention×realization
near-null. Genomic-context (P): host-coupled retain+realize more than intergenic (antisense exception noted).
Hallmark (P, descriptive): metabolic/differentiation programs most realized. Res-5 within-family at chance (prereg
null). Own-retaining 7% not clinically distinct (fragile purity hit, n=5). Details in MH-158's registry row.
✅ **CLASS → REALIZATION (MH-160, rigor-auditor A):** the cross-state `shift_class` PREDICTS within-patient
realization, ordered by repression-LIVENESS not dose — live-in-NAT edges {constitutive, field_established}
realize more than acquired-only ones (Δ −0.039, gene-clustered) despite LOWER dose; anti-dose (dose+coupling_tum
controlled −0.032), decoy-arbitrated (direct paired diff −0.024, p=0.01); partial circularity (ρ=0.234 with the
defining coupling_nat) quantified + defused. Enabled by the integrated `progression_edge_card`.
✅ **`shift_class` REBUILT → CALIBRATED TWO-AXIS (MH-166, 2026-07-19):** the old coupling-only −0.1 label is
SUPERSEDED by a two-axis annotation = per-state **calibrated** coupling (`site_free_null` per state, replaces −0.1;
tumour null 2.67× reproduces MH-123) × **same-platform** `arm_lfc_NAT_TUM` dose (replaces the soft QN healthy leg).
New quantified cols (`coupling_p_*`, `coupling_z_tum`, `realization_score`, `dose_class`, `wiring_frac`) + NEW class
`dose_acquired_uncoupled` (571 edges: dose up, no calibrated coupling). ⭐ **At EQUAL dose the calibrated coupling
axis alone separates realization** — `dose_acquired_uncoupled` +0.002 vs `acquired_realized` −0.101 (dose 1.13 vs
1.19). MH-160 RE-VALIDATED: core ordering holds/strengthens (−0.046) but the decoy-arbitrated increment is now
MARGINAL (−0.026, p=0.098) as `constitutive` collapsed 127→6 under calibration. Home: `card.py::_shift_class`.
⚠ **Healthy-leg fallacy (MH-166 follow-up, 2026-07-19):** the GTEx coupling leg is BLIND where the uniquely-mappable
pipeline collapsed multi-mapping arms. Added `healthy_leg`/`healthy_potential`/graded `healthy_uninformative`
(miTED-aware, 305 flagged; contamination 24%→13%; MH-160 robust). Flagged arms = dose-constitutive, coupling-blind;
recoverable via same-seed co-transcription surrogate (miR-17 via miR-20a ⇒ coupling-ACQUIRED). Genome-wide surrogate
node + per-target GTEx→NAT→TUM identity trajectory = NEXT ARC. See FORMULAS §11g + DISCOVERY_REGISTRY MH-166.
⭐⭐ **EXTERNAL VALIDATION ARC (MH-166, 2026-07-19) — the big current-state update:**
• **PRECISION validated externally, TWO orthogonal axes:** McGeary&Bartel functional (19 arms, Fisher **p=7.96e-33**,
  HE **5.6×** enriched vs null) + breast chimeric binding GSE263552 (**27×**). The edges we use are real.
• **Family-collapse validated FUNCTIONALLY** (miR-29a≈29b, per-target ρ 0.87–0.95) — but only the seed-HOMOGENEOUS
  case; within-family TCGA arm-corr is HETEROGENEOUS (med 0.47, 6 tight/6 loose; miR-29 raw 0.44 ≈ partial 0.43,
  so co-expression ≠ the 0.87 seed-driven interchangeability). ≈25% seed-heterogeneous families still gated per-family.
• **⭐ COUPLING IS REALIZATION, NOT EXISTENCE (supersedes the "weak set-level signal" framing above):** breast TCGA
  coupling is INDEPENDENT of edge reality (Fisher OR=0.99 vs McGeary-functional) and ~13% SENSITIVE — it measures
  realized DOMINANT regulation across patients (observational, variance-driven), a harder/different question than
  existence. High-specificity (2.9–5× null) / low-sensitivity. ⇒ judge EXISTENCE on functional/binding; use coupling
  as REALIZATION (presence = dominant realized regulator; absence ≠ not-real). Do NOT use coupling as an existence arbiter.
• **RECALL → DISCOVERY (expansion pending):** HE captures ~6% of functional targets (precision-by-design). Built the
  tiered candidacy + `model_expansion_list.tsv` (functional × seed-type/Kd × context × chimeric × family-mate) tied to
  the project's coupling-discovery lane → **21 gold∩external crown jewels** (miR-17~92→TMBIM6/AHNAK/XIAP; miR-182→MET
  is literature-VALIDATED in breast, the rest novel-but-assay-backed). ⚠ **abundance is NOT a "tilt"/artifact — abundant
  arms being dominant realized regulators IS the biology; a "real edge" = structure(seed) + physical(binding) +
  impact(coupling), all three. The abundance association only limits RECALL (the low-abundance tail), not precision —
  do NOT down-weight coupling.** ⬜ **OPEN: MODEL EXPANSION + Gibbs REFIT (board item 0) — fold candidates in, re-fit,
  read new β + gene pressure.** Files in `output/learned/realization/{model_expansion_list,missed_edges_*,breast_chimeric_*}.tsv`.

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
| `CPTAC_PROTEIN_CHANNEL_PLAN.md` | ✅ **FIXED + ARCHIVED 2026-07-17.** §6's VOID PAM50 numbers (1.7% / ≤8.8%) corrected to **4–6% / ≤7.6%**, the self-nested correction strings de-nested, and a third live citation of the void figure found at §0a. The doc is now in [`archive/`](archive/README.md) — its channel is dead (`βᵗ` falsified at n=101) and its **survivors were extracted to `LEARNED_MODEL_VALIDATION.md` §6.1–§6.5** (`locus_cn_cptac` · δ-transportability · Bar-5 · `a_g` · the Fisher bound). |
| `handoffs/HANDOFF_PROTEIN_AND_COMPOSITION.md` | Dated 2026-07-12; **12+ rows have landed since**. Superseded on the estimator (MH-115), the decoy control (MH-127/130/132), and the null calibration (MH-123/124). |
| `PROGRAM_FORWARD_BOARD.md`, `WHATS_NEXT.md`, `LEARNED_MODEL_WHATS_NEXT.md` | **Three competing forward docs, all predating MH-133…137.** All stale on the CN instrument. Consolidation pending. |
| `LANDSCAPE_REPORT.md` | Its per-edge FDR class counts rest on the **uncalibrated null** (MH-123/124). No banner yet. |
| `presentation/` F25 + `talk.pptx` | `talk.qmd` prose was fixed 2026-07-12 and is honest; **the F25 figure beside it is from 2026-06-24** and shows the old numbers. Deck needs re-rendering. |
