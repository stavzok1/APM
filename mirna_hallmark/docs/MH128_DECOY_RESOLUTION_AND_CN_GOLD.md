# MH-128 — THE DECOY, RESOLVED · and the CN instrument is NOT the attribution gold standard

> **Goal:** settle what the "decoy" result of MH-124/MH-127 actually means — was the decoy HANDICAPPED, does a
> FITTED FAKE regulator set beat the real one, and can the CN instrument replace the literature MTI list as the
> arbiter of attribution.
> **What belongs here:** the three questions, each with a verdict sentence first; the corrected magnitudes that
> supersede MH-124's real-vs-decoy gaps; every adversarial refutation that beat a phase result (the refutation
> wins); and a strict separation of NEGATIVE from UNDERPOWERED.
> **What does NOT belong here:** biology (→ `BIOLOGICAL_SYNTHESIS.md`); the model's design (→ the METHODS twins);
> anything not measured in this arc.
> **Update trigger:** any new measurement on decoy construction, the real-vs-fake ladder, or the CN instrument's
> role as an attribution arbiter.
> **Sync-partner:** `MH124_ANTICOUPLING_VALIDITY.md` (§2.4/§4b/§5 — annotated by this doc),
> `MH127_DECOY_MODEL_GENE_BUDGET.md`, `CN_INSTRUMENT.md` §9x/§9y, `ATTRIBUTION_PRIMITIVE.md` §6,
> `ANALYSIS_RUN_LEDGER.md` (MH-128), `DISCOVERY_REGISTRY.md` (MH-128a/b/c).

---

## ⭐ THE ONE THING TO KNOW

**The decoys were never "impossible arms the model stupidly loves". They are arms that carry the same global
covariance structure as the real regulators — and MH-124's strict decoy was HANDICAPPED on exactly that
structure, so every real-vs-decoy gap MH-124 published is ~2.2× too large.**

Remove the handicap and three things happen at once:

1. **The real edge set still separates** from an impossible one (realized anti-correlation −0.0353 vs −0.0227,
   MWU p=3.0e-04) — **edge existence survives**, but the honest site-driven effect is **≈ −0.02**, not the
   ≈ −0.039 implied by MH-124's −0.050-minus-−0.011.
2. **A FITTED FAKE reproduces ~75% of the learned model's held-out prediction** (β, vs the full real model;
   95% CI [0.69, 0.81]) **and decisively beats REAL UNFITTED ABUNDANCE** (p=1.7e-19). So MH-115's "the learned β
   beats mere abundance" is, to ~61% of its margin, a claim about **FITTING**, not about **CURATION** — which is
   exactly what MH-127 concluded on a different decoy construction. Two independent decoy rules, same verdict.
3. **The CN instrument cannot rescue attribution.** The quantity that was scored — `pi_causal = γ_s·b_fam` —
   contains **β's own observational slope**. β agrees with `b_fam` (the pure observational slope, *no copy number
   anywhere in it*) at **+0.873**, *higher* than with `pi_causal` (+0.736), and a Y-permutation null reproduces
   **+0.618** of that +0.736. Against the correct null, **IDENTITY (Shapley) agrees with the CN allocation at
   least as well as β does** ⇒ **MH-124 §4b is NOT refuted**, and the CN instrument does **not** supersede the
   literature list.

**Read the model accordingly: it is a decent DETECTOR of the curated edge SET and a poor ALLOCATOR of credit
inside it — and its in-cohort predictive advantage is mostly the fit, not the biology.**

---

## 0. Provenance of every number below

| block | produced by | artifacts |
|---|---|---|
| Phase 1 (fake-set rebuild, match quality, realized ρ) | this arc, phase-1 agent | `scratchpad/fakeset_{NEW_axis_matched,OLD_strict_rawcorr,REAL_regulators}.tsv`, `fakeaxis_pairs_ALL.tsv`, `fakeaxis_{build,analyse,finish}.py`, `PREDICTIONS_fakeaxis.md` |
| Phase 2 (the 4-rung OOF ladder) | this arc, phase-2 agent | `scratchpad/l2_ladder.tsv`, `l2_pergene.tsv`, `l2_h2h_pairs.tsv`, `l2_gene_clusters.tsv`, `l2_ladder.py` |
| Phase 3 (CN gold standard) | this arc, phase-3 agent | `scratchpad/mh127/{core.py, gate.py, analyse.py, rows_pass1.tsv, REPORT.txt, PREDICTION.md}` |
| Adversarial R1 (handicap direction, ladder) | this arc, refuter | `scratchpad/adv_*`, `adv_bg.tsv` |
| Adversarial R2 (leakage / circularity / capacity) | this arc, refuter; **re-run and reproduced by the synthesist** | `scratchpad/kill_{perm,capacity,pool,select,rand,fullreal,final}.py`, `kill_{perm,rand,fullreal}.{log,tsv}` |
| Adversarial R3 (power vs negative) | this arc, refuter; **re-run and reproduced by the synthesist** | `scratchpad/pwr{1,2,3,4}.py`, `pwr_{mde,concentration,q5,ratio,nneeded}.tsv` |
| Adversarial R4/R5/R6 (CN circularity + null exchangeability) | this arc, refuters; **`analyse_null.py` re-run and reproduced by the synthesist** | `scratchpad/mh127/{attack1,attack2,attack3,analyse_null}.py`, `null_ylane.tsv`, `bfam.log`, `attack2.log` |

**No repo module was modified by this arc.** Scale throughout: 763 genes / 4,249 real (gene,arm) HE edges;
head-to-head subset **753 genes / 3,104 edges**; TCGA n=1041 (Normal-like excluded).

> ⚠ **Do not mix scales.** This arc's REAL realized anti-correlation is **−0.0353** (head-to-head) / **−0.0408**
> (all 4,249 edges), not MH-124's **−0.050**, because the site rule here is *stricter* (TargetScan context++ site
> **AND** scanMiR duplex **AND** the curated edge) and the head-to-head subset restricts REAL to arms for which a
> legal fake exists. **Every ratio quoted below is measured on the SAME rows in the SAME run** and is internally
> valid; the absolute levels are not interchangeable with MH-124's.

---

## 1. Q1 — WAS THE STRICT DECOY HANDICAPPED?

### VERDICT: **YES. CONFIRMED.** The old decoy under-loads the global axes, and it does so *because of the
### `|r|<0.30` raw-correlation cap*. Every MH-124 real-vs-decoy gap on a coupling readout is **≈2.2× too large**
### — and the real set **still separates** from an un-handicapped fake.

**The handicap, measured.** Both fake sets were rebuilt through **one** code path so the rules differ in exactly
one thing (site-free ∧ non-seed-family ∧ non-cluster ∧ abundance caliper 0.25 is common to both; the OLD rule
adds a `|r_raw|<0.30` cap and picks at random, the NEW rule caps on the *sparse* residual correlation and
minimises the distance to the real arm's factor **loadings**). Loadings are read out on a global miRNA-PC basis
(descriptive only — never a covariate); the matching itself used a **leave-the-design-out** per-gene factor basis.

| head-to-head, 3,104 edges / 753 genes | REAL | OLD (strict) | NEW (axis-matched) | OLD p | NEW p |
|---|---|---|---|---|---|
| abundance (log2 RPM) | 5.329 | 5.320 | 5.315 | 0.99 ✓ | 0.89 ✓ |
| **\|loading\| mPC1 (proliferation axis)** | **0.323** | **0.254** | **0.329** | **6.4e-68** ⛔ | **0.22** ✓ |
| \|loading\| mPC2 (purity axis) | 0.242 | 0.195 | 0.200 | 1.9e-46 ⛔ | 5.8e-32 ⛔ **(still unfixed)** |
| detection rate | 0.885 | 0.883 | 0.883 | 0.16 ✓ | 0.95 ✓ |
| HE degree (promiscuity) | 26.27 | 9.56 | 10.60 | ~0 ⛔ | ~0 ⛔ (both) |
| yield (real arms that got a legal fake) | — | 70.7% | 86.0% | | |

**⭐ The consequence — realized anti-correlation, partial Spearman(arm, target \| C):**

| set | ρ | gap vs REAL | test |
|---|---|---|---|
| **REAL** | **−0.0353** (n=3,104) | — | |
| OLD strict decoy | −0.0069 (n=14,993) | **−0.0285** | MWU p=3.8e-20 |
| **NEW axis-matched fake** | **−0.0227** (n=15,457) | **−0.0126** | MWU **p=3.0e-04**; paired-Wilcoxon p=3.1e-17 |
| NEW vs OLD (the fake is 3.3× more anti-correlated once un-handicapped) | | | MWU p=4.5e-26 |

⇒ **OVERSTATEMENT FACTOR = 2.25×** (all-pairs panel: 2.12× — same answer).
⇒ **REAL STILL SEPARATES** from the un-handicapped fake (p=3.0e-04).

**The residual gap is BINDING, not degree.** OLS with gene-clustered SEs: b(is_real) OLD −0.0285 raw →
−0.0218 (+abundance, loadings) → −0.0217 (+degree), p=5.4e-06; NEW −0.0126 raw → −0.0198 → **−0.0173,
p=5.8e-07**. Two independent routes — matching *or* adjustment — **converge on a site-driven effect of ≈ −0.02**.

**The arms actually shipped in MH-124.** Measured on `output/learned/decoy_v3_abundance_matched.tsv`
(482 genes / 1,663 pairs): |mPC1| REAL 0.342 vs DECOY **0.275** (ratio 0.80, p=1.8e-21); |mPC2| 0.270 vs **0.221**
(0.82, p=9.9e-18); abundance matched (p=0.80 ✓); detection matched (p=0.32 ✓). ⇒ **the decoys MH-124 actually
used carry a ~20% relative loading deficit on BOTH global axes.**

### ⚠ ADVERSARIAL R1 — the direction audit (CONFIRMED, and it *located* the handicap)

Which fake-set rule biases which way? Measured over the site-free abundance-eligible base pool (725 genes;
un-handicapped pool background ρ = **−0.01842**, REAL arms −0.03554):

1. **The OLD `|r_raw|<0.30` cap FAVOURS REAL** — the pass-set sits at −0.00963, i.e. the cap **strips
   +0.00879** of background (p=1.2e-10, CI [+0.0042, +0.0136]) ≈ **48%** of the pool's background; the arms it
   *rejects* carry **−0.03150**, nearly the real level. **This is the handicap, and it is the raw-correlation
   cap, not the abundance match.**
2. **The NEW sparse cap favours the FAKE**, microscopically (pass-set −0.01985, Δ **−0.00143**, p=1.4e-04).
3. **The loading match itself is NEUTRAL** (‖λ‖ REAL 0.5375 vs NEW 0.5394).

⇒ Both residual asymmetries in the NEW design **flatter the fake**, so "real still beats fake" is
**CONSERVATIVE, not inflated**.

### ⚠ ADVERSARIAL R2 — a SECOND selection, pointing the other way (and it re-denominates the recovery ratio)

Requiring every real arm to *have* a legal fake **deletes 34.3% of real arms**, and the deleted ones are **not
random**: within the 753 ladder genes, RETAINED real arms have ρ **−0.0353** and DROPPED ones **−0.0532**
(t-p=2.5e-03, MWU p=3.0e-03), and DROPPED arms are far more abundant (median log2 RPM 10.29 vs 5.33) and more
promiscuous (HE degree 42.1 vs 26.3). ⇒ **the head-to-head REAL set is a WEAKENED real set** ⇒ any
`|fake|/|real|` recovery ratio computed on it is **INFLATED**. (This is the mirror image of MH-127's
"abundance-matching is structurally impossible at the top of the range" finding — the same arms, from the other
side.) **Fix: quote decoy recovery against the FULL real model.** See §2.

⇒ **NET:** the 2.25× overstatement of the *anti-correlation gap* stands (identical REAL rows on both sides);
recovery *ratios* must be denominated on the full model.

---

## 2. Q2 — DOES THE FITTED-REAL GENE BUDGET BEAT THE FITTED-FAKE ONE?

### VERDICT: **Barely, and not at every resolution × weighting once the honest null is used. What IS decisive is
### that a FITTED FAKE beats REAL UNFITTED ABUNDANCE (p=1.7e-19).** And the predicted "wrong functional" headline
### is **REFUTED in the opposite direction**: β is the functional most willing to hand credit to a decoy — the
### decoy recovers *more* of the real model under β than under identity or dose, not less.

**The 4-rung ladder** (753 genes, 3,104 head-to-head edges, 3 fake draws, strict 5-fold OOF over *patients*;
C fit on train, z-scoring + the `sd<0.1` floor decided on train; `_gibbs_posterior(pi=ones)` = β, exact-LMG
`shapley_identity` = identity, train-mean level = dose; family unit = `collapse_by_family` TRUE RPM pool).
Mean per-gene OOF Spearman; negative = the repressive field works.

| rung | | value |
|---|---|---|
| **R1** REAL, UNFITTED (abundance pool) | | **−0.0581** |
| **R2** FAKE, UNFITTED | OLD / NEW | −0.0201 / **−0.0378** |
| **R3** REAL, FITTED (arm) | β / shapley / dose | **−0.1326** / −0.0955 / −0.0638 |
| **R4** FAKE, FITTED (arm, NEW) | β / shapley / dose | **−0.1152** / −0.0756 / −0.0450 |

### 2a. ⭐ DECOY RECOVERY — correctly denominated (R2's re-run, 753 genes, arm resolution, identical folds)

| weighting | vs the head-to-head (capacity-matched, 3.98-arm) real model | **vs THE LEARNED MODEL (6.05-arm, all HE regulators)** |
|---|---|---|
| **β** | 0.881 [0.814, 0.958] | **0.747 [0.691, 0.807]** |
| shapley (identity) | 0.791 [0.666, 0.933] | **0.623 [0.526, 0.729]** |
| dose | 0.713 [0.565, 0.876] | **0.673 [0.533, 0.823]** |

*(gene bootstrap, 4,000 reps.)* The OLD (handicapped) decoy reproduces MH-124's published "decoy = 71% of the
real field" almost exactly (0.697–0.720 at β on the head-to-head) — **so MH-124's 71% is a HANDICAPPED-decoy
number, and the honest figure against the full learned model is ~0.75.**

**The functional ordering, honestly:** the robust part is that **β is where the decoy does BEST** (all
denominators agree). The *dose < shapley < β* ordering claimed in phase 2 is **NOT resolved**: a paired
cluster-bootstrap gives **P(full order) = 0.79 (family) / 0.83 (arm)**, and against the full-model denominator
the middle two **swap** (shapley 0.623 < dose 0.673 < β 0.747). ⇒ **Report only: β hands the most credit to a
decoy.** ⛔ The pre-registered story ("the field was built with the wrong functional; identity or dose would
separate") is **REFUTED in direction** — β does *not* under-separate relative to identity; it over-shares.

### 2b. ⭐ THE DECISIVE RESULT — a fitted FAKE beats real UNFITTED abundance

| contrast (753 genes) | value | p |
|---|---|---|
| fitted FAKE (arm/β) vs the **FULL** real abundance pool | −0.1171 vs −0.0554, **gap −0.0617** | **1.7e-19** (Wilcoxon); cluster-block permutation **p=5e-05** |
| fitted FAKE (arm/shapley) vs FULL abundance pool | −0.0758 vs −0.0554, gap −0.0204 | 2.9e-03 |
| fitted FAKE (arm/dose) vs FULL abundance pool | −0.0457 vs −0.0554, gap **+0.0097** | 0.085 — **abundance wins** |
| **the real model's own margin over the same pool** | −0.1568 vs −0.0554, margin −0.1014 | — |
| ⇒ **the decoy reproduces 61% of the learned model's margin over mere abundance** | | |

Independently reproduced by R1 on a **fully unhandicapped** decoy (no correlation cap at all): gap −0.0528,
p<1e-4. ⇒ **MH-115's "β beats mere abundance" must be read as a FITTING claim.** This is the *third* independent
construction (MH-127's trimmed decoy; this arc's axis-matched decoy; R1's uncapped decoy) to reach it.

### 2c. Real-beats-fake — where it actually holds, under the honest null

Genes are **not independent** (regulator-set Jaccard → Louvain, 43 blocks). Measured **ICC = 0.026–0.090**,
**DEFF = 1.23–2.48** ⇒ n_eff = 305–611 of 753 genes. The naive Wilcoxon is **10–40× anti-conservative**.
Cluster-block **sign-flip permutation** (20k reps) is the honest test:

| Q2 cell (REAL-fitted − FAKE-fitted, NEW axis-matched fake) | gap | Wilcoxon p | **block-perm p** |
|---|---|---|---|
| **arm / β** | −0.0174 | 3.3e-05 | **0.017** ✅ |
| arm / shapley | −0.0199 | 8.2e-04 | 0.054 |
| arm / dose | −0.0187 | 8.8e-04 | 0.076 |
| fam / β | −0.0130 | 1.4e-03 | 0.055 |
| fam / dose | −0.0167 | 3.1e-03 | 0.094 |
| fam / shapley | −0.0149 | 7.7e-03 | 0.132 |
| **Q4** REAL-unfitted − FAKE-unfitted (NEW) | −0.0203 | 2.2e-03 | **0.063** |

⇒ **Only ONE of six cells (arm/β) survives at α=0.05**, and **"curation carries signal with no fitting"
(Q4-NEW) FAILS** the honest null. Real beats fake on only **53–59%** of genes — a coin flip plus a nudge.
**⚠ UNDERPOWERED, NOT NEGATIVE:** genes needed for 80% power at the observed effects = **697–1,619** (we have
753). The direction is stable in all 3 draws (per-draw sd of the gap 0.0005–0.0025) and 87% of the variance is
**true gene heterogeneity**, so more patients/iterations/draws cannot help — **only more genes can.**
(Against the OLD handicapped decoy, every cell "passes" — that is the 2.2× overstatement showing up again:
Q2 arm/β OLD block-perm p=0.017 at gap −0.0402 vs NEW −0.0174.)

### 2d. ⛔ ARM > FAMILY — **WITHDRAWN. It is a DESIGN ARTIFACT.**

Phase 2 reported "the ARM design separates real-from-fake ~1.75× better than FAMILY" (and it *does* survive the
cluster permutation, block-perm p=5e-05). **But the cause is the fake's construction, not the resolution.**
Decoy rule (e) forbids a fake from being a **seed-family mate** of any real regulator ⇒ **the fake design is
structurally incapable of collapsing under family pooling**:

* On the 199 "designs-differ" genes: **REAL loses 1.524 columns** to family pooling; the **NEW fake loses 0.077**
  (OLD 0.051). On that subset the fake family model gets **+1.447 extra columns** on average (95.7% of gene-draws).
* Exact decomposition (β): `arm_gap − fam_gap = −0.0174` = **pooling cost to REAL (−0.0179, p=7.3e-08) minus
  pooling cost to FAKE (−0.0004, p=0.47)**. The identity closes to 4 decimals for β, shapley and dose.
* On the **19** gene-draws where the two designs are genuinely capacity-matched (`nfam_fake == nfam_real`), the
  arm-minus-family contrast is **n.s.** (p=0.62 / 0.49 / 0.57). *(Tiny n — but the decomposition identity above
  does not need power.)*

⇒ **"FAMILY resolution" pools the real model and leaves the fake model alone.** Q5 measured the **cost of
pooling to the real design**, not a resolution's ability to separate. **The ARM-vs-FAMILY clause is retracted.**

---

## 3. Q3 — DOES THE MODEL'S CREDIT ALLOCATION AGREE WITH THE CN-CAUSAL ALLOCATION?

### VERDICT: **The CN instrument CANNOT be the model-blind gold standard — the quantity that was scored is β's
### own estimator.** What survives a leak-free null: **the model's covariance allocations DO agree with the
### exogenous CN reduced form — and IDENTITY agrees at least as well as β.** ⇒ **"β > identity > dose" is
### REFUTED; MH-124 §4b stands.** And the reframe's own premise dies with it: **the CN instrument is not a LEVEL
### perturbation in the estimand sense** — its reduced form is a cross-patient covariance functional, so it can
### only arbitrate covariance allocations. **It does not supersede the literature MTI list; they answer different
### questions.**

### 3a. What the phase-3 run measured (466 genes = a randomised 30% subsample; exclusion-clean, F>10, C+ploidy)

| GT = `pi_causal` (137 genes / 444 families) | per-gene ρ | Wilcoxon | phase-3's within-gene shuffle null |
|---|---|---|---|
| **β** | **+0.736** | 2.5e-17 | perm p=0.001 |
| **identity** (Shapley/LMG) | **+0.453** | 2.4e-07 | perm p=0.001 |
| **dose** | **+0.093** | 0.227 | perm p=0.099 |

Phase 3 concluded "β > identity > dose ⇒ MH-124 §4b refuted". **That conclusion is wrong, and three independent
refuters killed it the same way.**

### 3b. ⛔ THE LEAK — `pi_causal` IS β's OWN ESTIMATOR

`instrument.exclusion` (`instrument.py:1116-1136`) returns **`pi_causal = γ_s · b_fam`**, where `b_fam` is the
**observational** OLS slope of Y on X_fam given [C, Z_s] — fit on the **same Y, same X, same samples, same C**
as β. (`attribution_eb._prep` fits β on `yr = −resid(Y|C)` ⇒ β ≈ −b.) **This is exactly the defect MH-126 already
established** — and it re-surfaces here from inside the attribution question.

**Measured (`scratchpad/mh127/bfam.log`, the clean set, 137 genes / 444 pairs):**

| model | vs **`b_fam`** — the purely observational slope, **NO copy number anywhere in it** | vs `pi_causal` | vs `pi_raw` (the reduced form) |
|---|---|---|---|
| β | **+0.873** (p=4.0e-24) | +0.736 | +0.465 |
| identity | +0.523 (p=8.1e-09) | +0.453 | +0.326 |
| dose | +0.031 (p=0.69) | +0.093 | +0.091 |

**β agrees BETTER with a quantity containing zero copy number than with `pi_causal`.** ⇒ the "CN-causal
allocation" ranking is largely **an identity of the regression**, not a causal test.

### 3c. ⭐ THE CORRECT NULL — and it REVERSES the ordering

The three allocations do **not** share a chance baseline. A **Y-residual permutation null** (design, X, CN matrix
and C held 100% fixed; only `resid(Y|C)` permuted across patients — so the *only* thing destroyed is the true
signal, while the shared-Y-noise route between β and `pi_causal` is preserved) gives, over 25 permutations
(`mh127/analyse_null.py` → `null_ylane.tsv`; re-run and reproduced by the synthesist):

| GT | model | **observed** | **null mean** | null sd | p | **excess over null** |
|---|---|---|---|---|---|---|
| `pi_causal` | **β** | +0.736 | **+0.618** | 0.053 | 0.038 | **+0.118** |
| `pi_causal` | **identity** | +0.453 | **+0.104** | 0.089 | 0.038 | **+0.349** |
| `pi_causal` | dose | +0.093 | +0.006 | 0.083 | 0.231 | +0.087 (n.s.) |
| **`pi_raw`** (exogenous reduced form) | **β** | +0.458 | +0.149 | 0.062 | 0.038 | **+0.309** |
| **`pi_raw`** | **identity** | +0.334 | +0.010 | 0.085 | 0.038 | **+0.324** |
| **`pi_raw`** | dose | +0.099 | +0.008 | 0.065 | 0.115 | +0.091 (n.s.) |
| `b_fam` (no CN at all) | β | +0.873 | +0.685 | 0.053 | 0.038 | +0.188 |
| `b_fam` | identity | +0.523 | +0.121 | 0.097 | 0.038 | +0.402 |

*(Re-gating the exclusion filter inside each permutation — "noise + selection" — changes nothing: β/`pi_causal`
null +0.618, identity null +0.114.)*

**Read the LAST column, not the first.** β's chance baseline against `pi_causal` is **+0.618**; identity's is
**+0.104**. On the *shared* scale of excess-over-null, **identity (+0.349) beats β (+0.118) against `pi_causal`**,
and on the **fully exogenous `pi_raw`** they are **effectively tied (+0.324 vs +0.309)**. A refuter's independent
**sample-split** (fit β on one half, score the CN allocation on the other) agrees: β excess **+0.391 (p=0.004)**
on `pi_causal`; on `pi_raw`, β **+0.382 (p=1.7e-05)** and identity **+0.321 (p=1.9e-04)**.

⇒ **THE ANSWERS:**
* **Does the model's credit allocation agree with the CN-causal allocation? YES — beyond any leak.** The
  exogenous reduced form `pi_raw` agrees with **both** covariance allocations (β and identity), and a
  sample-split confirms it. This is a **real, leak-free agreement** and it is the arc's positive CN result.
* **Does IDENTITY agree better than BETA? At least as well — and clearly better against `pi_causal`.** The
  phase-3 "β > identity" ordering was **the leak talking**. ⇒ **MH-124 §4b (identity/dose rank the canonical
  family above chance; β is at chance) is NOT refuted by the CN instrument.**
* **DOSE shows no detectable agreement** with the CN allocation on any arbiter — **but this is UNDERPOWERED, NOT
  NEGATIVE**: the bootstrap CI on dose's `pi_causal` agreement is **[−0.048, +0.235]** with **MDE 0.201**. The
  data cannot exclude a true dose agreement up to ≈+0.24. It can only say **dose < β**.

### 3d. ⛔ AND THE PREMISE ITSELF — "the CN instrument is nature's knockdown, the same estimand as the literature"

**REFUTED.** Copy number perturbs dose **across patients**, so its reduced form is a **covariance / slope**
functional — and empirically it ranks families the way a covariance functional does (β +0.31, identity +0.32
excess) and **not** the way DOSE does (+0.09, n.s.). The literature MTI list ranks the *opposite* way (MH-124
§4b, n=21 genes: dose rank 0.258 p=3e-4; identity 0.319 p=0.0076; **β at chance**, 0.486). **The two gold
standards are not the same estimand.** The right conclusion is **not** "β wins" but:

> **Each arbiter measures the functional it is built from. The CN instrument vindicates the model's covariance
> allocation; it cannot adjudicate a LEVEL claim any more than β can.**

The user's reframe — *the literature list is a biased, estimand-mismatched gold standard* — **stands**. But the
CN instrument is **not** the unbiased replacement for it. There is currently **no** exogenous, level-perturbing,
well-powered arbiter of attribution in this project.

### 3e. ✅ ONE CLEAN, USABLE POSITIVE — put ploidy in C for every CN readout

Adding a per-sample **mean ASCAT3 CN over the 506 miRNA hairpin loci** (`instrument.locus_cn`) to C, measured on
the same 466-gene subsample:

| | C = core | **C + ploidy** |
|---|---|---|
| well-instrumented families (F>10) | 1,051 (47.8%) | 1,162 (52.8%) |
| **exclusion-clean** (F>10, \|δ\|≤0.05) | 355 (16.1%) | **598 (27.2%)  (+68%)** |
| genes with ≥2 clean families | 80 | **137** |
| **frac `pi_raw` < 0 (the aneuploidy sign artifact)** | **0.878** | **0.566** |

⇒ **USE IT.** This **independently reproduces MH-126's** aneuploidy finding from a completely different entry
point, and it is a *coverage* win (+68% clean families), not only a calibration one. Every conclusion in §3
holds under both C blocks (core is weaker).

---

## 4. WHAT THIS ARC CHANGES IN THE STANDING RECORD

| standing claim | status after MH-128 |
|---|---|
| MH-124 §2.4: "decoys reach **71%** of the real predictive field; absorb **36%** of \|β\|; PIP 23.0% vs 15.5%" | **UPPER BOUNDS on separation, measured against a HANDICAPPED decoy.** The 71% is reproduced by the OLD decoy here (0.70–0.72) and becomes **~0.75 [0.69, 0.81]** against the full learned model with an axis-matched fake. ⚠ **NOT MEASURED:** the PIP and \|β\|-share panels were **not** re-run against the axis-matched fake. |
| MH-124: "real edges anti-correlate **5×** more than decoys (−0.050 vs −0.011)" | **The GAP is ~2.2× overstated.** Un-handicapped: −0.0353 vs −0.0227 (p=3.0e-04). **The honest site-driven effect is ≈ −0.02, not ≈ −0.039.** **Edge existence SURVIVES** (still separates, and survives conditioning on degree: b=−0.0173, p=5.8e-07). |
| MH-124 §4b: "identity/dose rank the canonical family above chance; β is at chance" | **STANDS.** The CN instrument does **not** refute it (§3c). |
| MH-124 §5 / Q2 / Q3a: the CN-instrument existence leg | **already retracted by MH-126**; this arc **independently re-derives the same defect** from the attribution side (β↔`b_fam` +0.873, higher than β↔`pi_causal`). |
| MH-115: "the learned β beats mere abundance" | **A FITTING claim, ~61% reproducible by a fake regulator set** (decoy margin −0.0617 vs the real model's −0.1014 over the same abundance pool; p=1.7e-19). Third independent decoy construction to say so. |
| MH-127: "the fitted decoy reaches 79–91% of the real model in-cohort" | **REPLICATED with a different decoy rule**: 0.747 [0.691, 0.807] at β against the full learned model. |
| **NEW (this arc):** "ARM resolution separates real-from-fake better than FAMILY" | **NEVER TRUE — withdrawn as a design artifact** before it was published anywhere but here (§2d). |

---

## 5. CAVEATS — read before quoting anything above

1. **NOT MEASURED:** the MH-124 **PIP**, **\|β\|-share** and **CPTAC-protein-transport** panels were **not**
   re-run against the axis-matched fake. The 2.25× applies to the **realized-anti-correlation** diagnostic and
   (via the ladder) to the **held-out-prediction** readout. Whether the PIP/transport gaps shrink by the same
   factor is **untested**.
2. **mPC2 is STILL handicapped in the NEW fake** (0.200 vs REAL 0.242, p=5.8e-32 — barely better than OLD's
   0.195). ⇒ **2.25× is a LOWER bound** on the overstatement.
3. **Promiscuity/degree is unmatched in BOTH fake sets** (HE degree 26.3 vs ~10). Partly unavoidable (an arm
   site-free on gene *g* necessarily has fewer sites) — and the real-vs-fake gap **survives conditioning on it**.
4. **The cluster bootstrap has only 43 blocks**; its CIs are coarse. That is the honest price of gene dependence,
   and it is what turns naive-Wilcoxon "significance" into 5-of-6 non-significant cells.
5. **Phase 3 ran on 30% of the gene universe** (466/1,571; gene order randomised, rows written incrementally ⇒
   unbiased subsample, but a subsample). Samplers were trimmed for compute (Gibbs 800/300 vs canonical 2000/700;
   Shapley 120 perms) — this adds noise to β and identity; **it does not favour either.**
6. **The Y-permutation null in §3c is a GLOBAL null** (no miRNA effect anywhere), and it has a **p-floor of
   1/26 = 0.038** at 25 permutations. Its *means* — the quantities the argument rests on — are estimated to
   ~0.01–0.02, far tighter than the gaps involved (+0.618 vs +0.006).
7. **The `coding`/`host` pleiotropy gate was NOT run** in phase 3 (node saturation). From reading
   `instrument.exclusion` (**INFERRED**), `coding`/`host` only set an aggregation *weight*, and 96.0% of
   well-instrumented families have exactly one strong segment ⇒ the weight is inert for them; its only real role
   would be as an extra filter. **This should be run to close the loop.**
8. **"Underpowered" is not "negative"** anywhere in this doc. The two live underpowered results are:
   real-fitted > fake-fitted at family resolution (needs ~1,100–1,600 genes), and dose's agreement with the CN
   allocation (CI up to +0.24).

---

## 6. WHAT TO DO NEXT (ordered by what would actually move a verdict)

1. **Re-run MH-124's PIP / |β|-share / CPTAC-transport panels against the axis-matched fake.** They are the last
   MH-124 numbers still quoted against a handicapped decoy.
2. **Scale the gene universe** (753 → ~1,600 with a fake). It is the *only* lever on the real-vs-fake gap:
   87% of the variance is true gene heterogeneity, so more patients / draws / Gibbs iterations cannot help.
3. **Build an mPC2-matched fake** — the current one still under-loads the purity axis, so the reported 2.25× is
   a lower bound.
4. **Stop scoring attribution against `pi_causal`.** If the CN instrument is used as an arbiter at all, the only
   admissible quantity is the **reduced form `pi_raw`** with **ploidy in C** — and even then it arbitrates a
   *covariance* claim, not a *level* claim.
5. **The level/dose estimand still has no well-powered arbiter.** The literature MTI list (n=21) is biased and
   tiny; the CN instrument is a covariance functional. A genuine level perturbation (a miRNA
   overexpression/knockdown panel with matched transcriptomes) is the missing data, not a missing estimator.
