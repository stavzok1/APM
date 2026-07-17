# MH-126 — ANEUPLOIDY IN THE CN INSTRUMENT, WHAT THE LEARNED MODEL IS *FOR*, AND THE GRADED ATTRIBUTION QUESTION

> **Goal:** answer, with measured evidence and adversarial verification, the three questions raised after
> MH-124: (1) is **aneuploidy** the missing covariate that makes the raw CN instrument fail, and does
> *modelling* it (rather than gating it) fix the instrument **and buy power**? (2) if the edge set comes from
> curation and the weights do not attribute, **what is the learned model actually for**? (3) when the model is
> fit over the curated set, **do the βs give the canonical families heavier weights** — the *graded* question
> (rank / top-k / evidence-gradient), not the argmax one.
> **What belongs here:** the design of each test, the CONTROL it carries, the measured result, the adversarial
> refutation where one landed, and the verdict. Refuted and underpowered results are first-class, and the two
> are kept strictly distinct.
> **What does NOT belong here:** anything not measured in this task; any number transformed from another
> document rather than re-derived; the retracted MH-118..121 arc.
> **Update trigger:** any new measurement on the CN instrument's covariate block, on the model's out-of-sample
> value, or on the evidence-gradient / attribution question.
> **Sync-partner:** `ANALYSIS_RUN_LEDGER.md` (MH-126 row), `DISCOVERY_REGISTRY.md` (MH-126 / MH-126b /
> MH-126c claims), `CN_INSTRUMENT.md` §9 (the covariate-block gap), `MH124_ANTICOUPLING_VALIDITY.md`
> (whose §5 CN-existence leg this document **downgrades**), `INDEX.md`.

**STATUS: COMPLETE.** Three questions asked; three answered. **Two of the three answers are partly negative,
and one of them retracts a leg of MH-124.** Every headline below was put through an adversarial verifier
(circularity lens · confound lens · power/unit-of-analysis lens); **where the verifier refuted a result, the
refutation is what is recorded.**

---

## ⭐ THE HEADLINE

1. **ANEUPLOIDY IS REAL AND IT IS A CALIBRATION FIX, NOT A POWER FIX.** A per-sample **ploidy** term (mean
   ASCAT3 CN over the 506 miRNA hairpin loci) removes a systematic sign bias that made **86.4% of *impossible*
   (site-free decoy) edges show a NEGATIVE CN→target reduced form**; with it, the decoy rate falls to **50.0%
   — exactly the null.** But the mechanism is **not** the co-amplification story that motivated it (that has
   the wrong sign), it does **not** improve HE-vs-decoy **discrimination** (ΔAUC = −0.002, n.s.), and it does
   **not** make attribution answerable (4/16 → 4/16, p=0.18).
2. **⛔ AND THE SAME WORK REFUTES MH-124's THIRD EXISTENCE LEG.** The statistic MH-124 called "the exogenous
   CN instrument" — `pi_causal` — **is not exogenous**: it is `γ_s · b_fam`, first stage × the **observational**
   OLS slope (ρ≈0.93 with the plain no-CN partial correlation; the entire "3× separation" was reproduced from a
   simulation with **zero causal effect**). The **only** genuinely exogenous factor, `γ_s`, is **site-blind**
   (HE +0.199 vs DECOY +0.198, ratio 1.004×, p=0.20), and the one genuinely exogenous statistic (the CN→target
   reduced form) shows **no HE-vs-decoy effect** once the unit of analysis is the **locus** (τ=−0.0015, p=0.67).
3. **THE LEARNED MODEL IS NOT AN ATTRIBUTION ENGINE — IT IS A DE-CONFOUNDED, JOINTLY-APPORTIONED,
   TRANSFERABLE COUPLING ESTIMATOR.** That is not rhetoric: β is the **only** aggregator in this subproject
   that predicts CPTAC protein out-of-sample, it **survives composition adjustment**, and it is **entirely
   additional to mere abundance** while abundance is **redundant given β** (MH-115, verified against the
   registry). What it cannot do is **name** the causal regulator — and that is an **identifiability** limit,
   not a fitting failure.
4. **THE GRADED QUESTION: YES, BUT WEAKLY, AND NOT ROBUSTLY.** Canonical (heavily-curated) families do get
   somewhat heavier β within a gene — but the effect is **small** (per-gene Spearman(β, w) = **+0.068**,
   permutation p=0.005), it **only just survives** a full family fixed effect (b≈+0.031, **p≈0.07–0.10 —
   borderline, not significant**), and **β still cannot NAME the canonical family** (mean normalized rank
   **0.518** = exactly chance, p=0.66). **Mere abundance names it far better (0.240, p<1e-4).**

---

## 0. What was run, and the circularity gate that licenses Q3

All scripts are in this session's scratchpad
(`/tmp/claude-207054/-sci-labs-michall-stavzok-APM/c6b5af90-f9bb-4c16-a3b3-2c02f8fc9b32/scratchpad/`):
`common.py`, `characterise.py`, `run_main.py`, `analyse.py`, `lit.py`, `sens.py` (Q1);
`w_circularity.py` (the gate); `mh126_build.py`, `mh126_lit.py`, `mh126_tests.py` (Q3);
`adv_*`, `circ_*`, `pow*` (the adversarial lenses). **No repo module was modified by this task.**

**THE GATE (CONFIRMED).** The graded test in §3 is only meaningful if the learned β never saw the evidence
weight `w`. It did not:

* **Fit path.** Over 48 genes (n=1041 participants; p = 2–40 families), refitting the canonical dense/coupling
  posterior (`readouts._posteriors` → `collapse_by_family` → `attribution_eb._prep` →
  `spike_slab._gibbs_posterior(π≡1)`) under **real vs within-gene-SHUFFLED vs CONSTANT** `w`:
  **max|Δβ| = 0.000000e+00** (0/48 genes moved at all), same for `beta_sd`, `pip_dense`, `X_fam`.
  Positive control fired: `w_fam` moved in 42/48 genes (max|Δ| = 11.77).
* **Source path.** Monkeypatching `learned.evidence.ledger.edge_weights` itself (so `assemble_gene`'s `w`,
  `collapse_by_family`'s `w_fam` **and** `priors.edge_priors` all see a permuted ledger): **max|Δβ| = 0.0**,
  while `prior_pi` moved in 42/48 (max 0.926) and `pip_discovery` in 41/48 (max 0.866).
* The adversarial verifier repeated this end-to-end through `readouts.gene_readouts` on PTEN/ZEB1/CDKN1A/BCL2
  (175 family units, 66.4% of 469,409 ledger rows changed): **max|Δbeta| = max|Δbeta_sd| = max|Δshare| =
  max|Δz| = max|Δpip_dense| = 0.0 exactly.**
* **Why:** `collapse_by_family` builds `X_fam = log2(1 + Σ_members RPM)` — pure abundance — and the `wf` it
  returns **is a dead variable in `readouts._posteriors`**; the dense Gibbs is called with `π = np.ones(p)`.

⚠ **TWO COLUMNS IN THE SAME FILE ARE w-CONTAMINATED BY CONSTRUCTION** — `pip_discovery` and `prior_pi` (the
evidence-π prior *is* `w`). Any "do canonical regulators score higher?" test run on those is **circular**.
**Use `beta`.** Likewise `spike_slab.fit_gene_ss` and `eb_shrink.fit_gene_eb` (reached via `mvp.oof_gate`)
**do** eat `w` and must not be substituted.
⚠ **Residual circularity that the gate does NOT close:** the candidate SET is evidence-defined (pooled-HE).
β is free of the graded weight; the universe is not free of the evidence **threshold**.

**Persisted-artifact provenance (verified in this task):** `output/learned/readouts_edges.tsv` — 5,117
gene×family rows, 1,549 genes, **0 negative β**, regenerated **2026-07-13 15:56** (i.e. post the MH-124b
`_rtnorm_pos` fix). All §3 numbers are from that file.

---

## 1. Q1 — "WE NEED ANEUPLOIDY"

### VERDICT

> **Aneuploidy is the missing covariate — and modelling it FIXES the instrument's CALIBRATION, but it does
> NOT buy power, and it does NOT make attribution answerable. The diagnosis behind the request was also
> wrong: the mechanism is a ploidy-driven SIGN FLIP, not co-amplification. And chasing it exposed that the
> statistic MH-124 validated existence with (`pi_causal`) was never exogenous in the first place.**

### 1.1 The baseline, reproduced

n = 3,111 edges / 482 genes (`decoy_v3_abundance_matched.tsv`); the re-implemented reduced form was verified
**bit-identical** to `instrument.edge_instrument` (max|diff| = 0.00e+00).
Decoy %neg reduced form **86.4%** (MH-124 reported 86.5%); HE 89.3%. HE exclusion-clean **34.6%** (MH-124:
34.5%). Clean π_causal HE −0.0221 vs decoy −0.0087, p=1.55e-05. **The MH-124 baseline reproduces.**

### 1.2 The mechanism — measured, and it is NOT co-amplification

* corr(CN_miRNA-locus, CN_target-gene) across participants: **HE median +0.563, DECOY +0.569, 100% positive in
  both, HE-vs-DECOY p=0.31.** The chromosomal tide is universal and completely **site-blind**.
* Co-amplification alone predicts a **POSITIVE** reduced form. The observed decoy reduced form is **NEGATIVE**.
  The measured mechanism is a **sign flip**:
  * corr(resid(CN_locus | C), ploidy) = **+0.437**, 100% positive (n=3,111)
  * corr(resid(Y | C), ploidy) = **−0.135**, negative in **94.6% of 482 genes**
  * corr(CN_target, ploidy) raw = +0.727, 99.4% positive
  i.e. once the gene's **own cis-dose** (`target_cn`, already in C) is removed, a higher-ploidy tumour has
  **lower residual expression of essentially every gene** (compositional/TPM dilution + over-adjustment on
  `target_cn`), while **any** miRNA locus's CN goes **up** with ploidy. Product = a negative reduced form for
  **any** pair, site or no site.
* **One per-gene scalar explains most of it:** R²(rho_reduced ~ corr(resid(Y|C), ploidy)) = **0.554** (n=1,627
  usable edges); Spearman(rz_ploidy × ry_ploidy, rho_reduced) = **+0.728, p=2e-268**.

### 1.3 The fix, and what it buys

`C_ANEU = C + [aneuploidy_score, fraction_genome_altered, ploidy]`, added to **both** the first stage and the
reduced form. (`C` from `assemble_gene` = `[CPE, target_cn, mal_prolif]` — **verified to contain no aneuploidy
term**. `ploidy` is DERIVED: per-participant mean ASCAT3 CN over the 506 miRNA hairpin loci in
`instrument.locus_cn()`.)

| readout | before | after |
|---|---|---|
| DECOY %neg reduced form | 86.4% | **50.0%** (= the null) |
| HE %neg reduced form | 89.3% | 53.6% |
| HE exclusion-clean | 34.6% | **55.6%** |
| DECOY exclusion-clean | 30.6% | **57.7%** |
| median \|δ_s\| (HE) | 0.078 | 0.044 |
| HE instrumentable (F>10) | 50.6% | 54.9% |

**Which term does the work:** the *official* aneuploidy measures barely help (Taylor arm-level aneuploidy score
alone: decoy %neg 81.5%; FGA alone 82.9%; both 82.2%). **The derived per-sample ploidy does all of it alone
(49.9%)** — and it *strengthens* the first stage (mean first-stage ρ +0.096 → +0.118), so it is removing
nuisance from the instrument, not signal.

**CONTROL (retention, paired on edges instrumented under both C variants):** π_causal **grows** for both
classes (HE 1.32×, DECOY 1.98×) ⇒ **this is not an over-control.** But it also does not preferentially kill
decoys.

### 1.4 ⛔ WHAT THE ADVERSARIAL VERIFIER REFUTED (this is what stands)

**(a) `pi_causal` IS NOT EXOGENOUS — and neither, therefore, was MH-124's "exogenous CN validation".**
`instrument.py` computes `pi_causal = γ_s · b_fam`, where `γ_s` is the first stage but **`b_fam` is the OLS
partial slope of Y on X_fam given [C, Z_s]** — i.e. estimated from exactly the variation in X_fam that is
**orthogonal to the instrument**, the **endogenous** variation. Verified by re-running `INS.exclusion`:
`pi_causal == γ_s · b_fam` exactly, and `b_fam` reproduces the plain **no-copy-number-anywhere** partial
correlation to ~2 decimals on every edge checked.
* Decomposing the headline (exclusion-clean, abundance-matched, n=422/side, C_ANEU): **`γ_s` — the only
  exogenous factor — is site-blind: HE +0.19931 vs DECOY +0.19844, ratio 1.004×, p=0.20.** All of the
  separation lives in `b_fam`, the observational slope.
* An independent verifier **reproduced the entire "3× exogenous CN validation" from a simulation with ZERO
  causal effect.**

**(b) THE ONE GENUINELY EXOGENOUS STATISTIC (the CN→target reduced form) DOES NOT SEPARATE HE FROM DECOY.**
Under a proper unit of analysis it dies: **locus as the unit, τ = −0.0015, p=0.67**; and matching arms on
their own site-free background turns the "rescued" raw instrument (p=0.038) into **p=0.71**.

**(c) EVERY HE-vs-DECOY p-VALUE IN §1.3 IS EDGE-LEVEL PSEUDO-REPLICATION.** 3,111 edges but only **482 genes /
409 arms / 277 loci**; the abundance-matching variable is a **per-arm constant** (409/409 arms have exactly one
value), so it can never remove a locus-level confound. Variance of π_causal explained by cluster: gene 0.372,
arm 0.352, locus 0.300. Under **two-way (gene + arm) fixed effects**, **59–72% of the claimed HE-vs-decoy gap
is gene+arm composition**, and **P5 ("the raw instrument is rescued") dies: p 0.038 → 0.70.**

**(d) THE FIX DOES NOT IMPROVE DISCRIMINATION AT ALL: ΔAUC = −0.002 (n.s.).**

**(e) THE ATTRIBUTION PAYOFF IS REFUTED (and it is a NEGATIVE, not a "cannot say", only in the narrow sense
that nothing MOVED).** Family-level exclusion on the 22 classic-MTI genes: instrumented 222 → 258,
exclusion-clean 80 → 124; but **answerable genes 15 → 16**, and the CN-causal π names the literature family
**4/15 = 27% (chance 15.8%, p=0.20) → 4/16 = 25% (chance 14.2%, p=0.18) — UNCHANGED.**
⚠ **This is UNDERPOWERED, not evidence of absence.** The binding constraint is the **literature set (n=16)**,
exactly as MH-124 said. Aneuploidy did not relax it.

### 1.5 WHAT SURVIVES (all three verifiers agree)

* The **calibration fix** (decoy %neg 86.4% → 50.0%) — cluster-robust, reproduced independently, and untouched
  by every attack.
* The **ploidy sign-flip mechanism** (+0.437 / −0.135 / R²=0.55; 94.6% of 482 genes) — a per-gene statistic,
  not power-vulnerable.
* The **exclusion-clean expansion** (34.6% → 55.6%) and the **first-stage strengthening**.
* The **attribution refutation** (4/16, p=0.18).

### 1.6 ⚠ Behavioural side-effect

**After the fix the exclusion GATE stops discriminating** (57.7% of decoys "clean" vs 55.6% of HE). The gate is
no longer a filter; whatever discrimination remains lives in π_causal — which §1.4(a) just showed is
observational. **Any analysis that uses the exclusion gate as a validity filter changes behaviour under
`C_ANEU`.** This is a contract change, not a tuning knob.

---

## 2. Q2 — "THEN HOW DOES LEARNING HELP US?"

### VERDICT

> **The learned model is not, and should not be sold as, an attribution engine. Its measured value is that
> it is the ONLY estimator in this subproject that (a) is DE-CONFOUNDED by construction, (b) APPORTIONS
> jointly under collinearity, and (c) TRANSFERS OUT OF SAMPLE — to a different cohort AND a different
> molecular layer — and it is entirely ADDITIONAL to mere abundance, which is redundant given it. What it
> cannot do is NAME the causal regulator among same-seed candidates, and that is an IDENTIFIABILITY limit
> (same-seed arms share the same binding site), not a fitting failure.**

The case is made from **registered, verified** measurements. Every number below was read out of
`DISCOVERY_REGISTRY.md` (MH-115 / MH-116-117) and `ANALYSIS_RUN_LEDGER.md` in this task, **not** taken from a
prompt or transformed by arithmetic.
⚠ **CAVEAT ON PROVENANCE:** a `coupling_grid` **field bug is under separate investigation** in a parallel
session. The MH-115/116/117 numbers below are quoted **as registered**; they were **not** re-derived here. If
that audit moves them, this section must be re-checked.

### (a) DE-CONFOUNDING — the reason β survives where the heuristic collapses

β is estimated **given C** (purity/CPE · `target_cn` · proliferation; + the composition block when requested).
A raw correlation, and the evidence-weighted **pressure heuristic** built on it, are not.
**MH-115 (registered):** on the *same* family-dose design, only the WEIGHT differing —

| aggregator | CPTAC protein, core C | CPTAC protein, **composition-adjusted** |
|---|---|---|
| abundance | −0.0060 (**p=0.20**) | +0.0040 (p=0.30) |
| heuristic pressure | −0.0079 (**p=0.11**) | +0.0031 (p=0.44) |
| **learned β** | **−0.0364 (p=5.3e-13)** | **−0.0126 (p=9.6e-04)** |

⇒ **only the learned β transports, and only the learned β survives composition.** The heuristic — the
aggregator behind MH-34, F25 and the talk — is p=0.83 on CPTAC mRNA (exactly zero) and is **beaten by raw
abundance even in-sample on TCGA** (−0.056 vs −0.046).

### (b) JOINT APPORTIONMENT under collinearity

The dense posterior fits **all** of a gene's family regulators jointly; the heuristic and the abundance
baseline score each regulator marginally. Measured consequence (MH-115, edge level, 2,773 HE edges vs a
measured CPTAC-protein partial-ρ):

* learned β Spearman **−0.278 (p=2.1e-50)** core / **−0.139 (p=1.8e-13)** composition-adjusted
* heuristic evidence −0.032 (p=0.096) · abundance −0.010 (p=0.61)
* **PARTIALS: β given abundance −0.140 (p=1.5e-13); abundance given β −0.012 (p=0.52)**
  ⇒ **β is entirely ADDITIONAL to abundance; abundance is REDUNDANT given β.**
* Only β shows a monotone dose-response across its quartiles; evidence and abundance quartiles are flat noise.

### (c) OUT-OF-SAMPLE TRANSFER — including to a universe curation cannot reach

* **OOF over patients** (MH-116/117, all 1,041 TCGA patients, 5-fold): gene HE **−0.126** OOF vs −0.177
  in-sample; edge HE −0.609 vs −0.712 ⇒ **~85% of the TCGA signal is retained out-of-fold — it is not
  circularity.**
* **ORPHAN universe** (miRTarBase-absent edges — **no curation exists, so the heuristic cannot even be
  formed**): the TCGA-mRNA-fitted β still predicts realized coupling in **CPTAC protein** (different patients,
  different layer, uncurated): **−0.189 core / −0.050 composition-adjusted.** **Abundance FLIPS SIGN there
  (+0.071 at TCGA edge level — the wrong direction)** ⇒ abundance is not merely weak off-curation, it is
  actively misleading.

### (d) The honest case AGAINST

* **It does not ATTRIBUTE.** MH-124: β names the literature-canonical regulator at chance (per-arm 42.9%,
  p=0.38; per-family 1/16 = 6% vs 8.1% chance). **§3 of this document extends that: β does not even RANK the
  canonical family above chance (0.518, p=0.66) — and mere abundance ranks it far better (0.240).**
* **Why — and this matters:** MH-124 §6 showed that signal and confound occupy the **same subspace** (removing
  low-rank structure kills the real signal *faster* than decoys, monotonically). Same-seed arms share the same
  binding site, so **no endogenous re-weighting can separate them.** This is an **identifiability** ceiling,
  and §3.4 measures it: median posterior SD/|β| = **0.799**; only **27.8% of 5,117 units are identified**
  (|z|>2).
* **And the exogenous escape hatch is now weaker than we thought:** §1.4 refutes the exogeneity of the CN
  statistic that MH-124 leaned on.

**⇒ The model earns its keep as a COUPLING and RANKING-OF-EDGES estimator with de-confounding and transfer.
It does not earn its keep as an answer to "which miRNA represses this gene". Do not sell it as one.**

---

## 3. Q3 — "DO THE βs GIVE THE CANONICAL FAMILIES HEAVIER WEIGHTS?" (the GRADED question)

### VERDICT

> **YES — but only as a WEAK GRADIENT, and it is borderline once family identity is absorbed. β does track
> curation weight within a gene (per-gene Spearman +0.068, permutation p=0.005; it beats a
> family-identity-preserving null at p=0.010) — but under a full FAMILY fixed effect the coefficient is
> b≈+0.031, p≈0.07–0.10, i.e. NOT significant. And β still CANNOT NAME the canonical family: its mean
> normalized rank is 0.518, exactly chance. Mere abundance names it much better (0.240, p<1e-4). "Heavier on
> average" and "the heaviest one" are different claims — only the first is even weakly supported.**

### 3.0 Setup

5,117 gene×family rows from the current post-bugfix `readouts_edges.tsv`; `w` = `ledger.edge_weights`
(max member weight = exactly `collapse_by_family`'s `w_fam`); within-gene CV of `w` median 0.373.
**Literature ground truth: 32 genes** (double MH-124's 16) — MH-124's 22-gene classic-MTI set plus 24 textbook
first-report MTIs, mapped gene → canonical **guide arm** → seed family via `families.family_of` (no substring
matching). **Sanity: `w`'s argmax family IS the literature family in 23/32 = 72%** (mean normalized w-rank
0.031) ⇒ the evidence gradient (Test 3, n=524 genes) is the **well-powered version** of the literature test.

### 3.1 TEST 1 — does β rank the canonical family? **NO.**

Normalized rank of the literature family (0 = heaviest, chance = 0.5), n=32 genes:

| readout | mean rank | median | Wilcoxon |
|---|---|---|---|
| **β** | **0.518** | 0.500 | **p=0.66 — EXACTLY CHANCE** |
| `shapley_identity` | 0.370 | 0.274 | p=0.0180 |
| **abundance** | **0.240** | 0.243 | **p<1e-4** |
| `w` (sanity) | 0.031 | — | — |

TOP-k (hits/32 vs Poisson-binomial chance): β **1/32 top-1** (chance 9.4%, p=0.96) · 9/32 top-3 (p=0.58) ·
16/32 top-5 (p=0.21) — **β is at or below chance at every k.** Paired: abundance − β = **−0.278** (better in
23/32, p=5e-04); shapley − β = −0.148 (18/32, p=0.0225).

**Robustness (adversarial power lens).** The 32 genes contain only **13 distinct canonical families**
(miR-34 ×5, miR-21 ×5, miR-15/16 ×4 …) and `abundance` is family-constant ⇒ the independent unit is the
**family**, not the gene. Re-tested by **cluster bootstrap over the 13 families (20k)**: abundance rank 0.240,
CI [0.126, 0.356], **p=0.0000** (design effect 2.73 — the SE really is 1.65× the iid one, and it still
survives); shapley 0.370, CI [0.237, 0.498], p=0.023; β's ranks are **not** clustered (DE 0.86) and remain at
chance. **Cluster Wilcoxon (unit = family, n=13): abundance p=0.0002 · shapley p=0.047 · β p=0.675.**
⇒ **Test 1's verdict is leverage-proof: β = chance; abundance beats it.**

### 3.2 TEST 2 — does the CN-causal π do better on the literature set? **CANNOT SAY.**

On the *same* instrumented candidate set: all-instrumented (n=20 genes) π mean rank 0.469 (p=0.35) vs β 0.478
(p=0.30), paired diff −0.009 (p=0.38); exclusion-clean (n=7) π 0.292 vs β 0.379 (p=0.41).
**No significant advantage — and this is UNDERPOWERED (n=20 / n=7), NOT a negative.** It also sits oddly
against Test 3, where the CN π has the strongest *raw* gradient. **Resolve by scaling the literature set, not
by trusting the n=7** — and note §1.4(a): π_causal is not exogenous anyway, so a win here would not have been
the exogenous validation it looks like.

### 3.3 TEST 3 — the evidence GRADIENT (3,797 gene×family rows; 524 genes with ≥3 curated families)

**(a) Per-gene Spearman(readout, w)**, mean [median], Wilcoxon:

| readout | mean ρ | p |
|---|---|---|
| **β** | **+0.0683** [+0.0982] | 4.2e-03 |
| shapley | +0.0715 [+0.1000] | 3.3e-03 |
| **abundance** | +0.1710 [+0.2428] | 5.2e-12 |
| CN π_causal (n=209) | +0.1345 | 6.2e-04 |
| `pip_discovery` | +0.3899 | **CIRCULAR — positive control only** |

**(b) Pooled within-gene rank regression, gene FE, cluster-robust SE (coefficient on rank(w)):**

| readout | raw | \| abundance | \| abundance + promiscuity + family-size | **\| FULL FAMILY FE** |
|---|---|---|---|---|
| **β** | +0.0516 (p=5.5e-03) | +0.0731 (p=8.5e-05) | +0.0506 (p=9.8e-03) | **b≈+0.031, p=0.074–0.097 ⚠ BORDERLINE** |
| shapley | +0.0695 (p=7.5e-05) | +0.0364 (p=0.042) | +0.0112 (**p=0.55**) | dead |
| CN π | +0.1148 (p=2.9e-04) | +0.1130 (p=4.7e-04) | +0.0787 (p=0.024) | **p=0.121 — dead** |

**(c) NULLS.** Permuting `w` **within gene** (1000×): β +0.0683 vs null 0.0006±0.0240, **p=0.005**; abundance
p=0.001; CN π p=0.001. **But the sharper null is the adversary's:** permuting `w` **within family, across its
genes** (2000×; preserves each family's `w` multiset exactly and destroys only the *edge-specific* part) —
**β obs +0.0676 vs null +0.0186±0.0204, p=0.0100 → SURVIVES** (family "fame" alone accounts for 27% of the
observed effect); **shapley obs +0.0708 vs null +0.0603, p=0.297 → DEAD (85% is family fame)**; CN π p=0.028.
Variance of the within-gene-demeaned ranks explained by a family main effect: β 0.326 · shapley 0.283 ·
abundance 0.731 · CN π 0.269 · w 0.208.

**(d) The study-bias channel is real but is a SUPPRESSOR, not a confound.** Spearman(abundance, w) within gene
= **+0.1710** (well-studied miRNAs *are* more abundant) — but Spearman(β, abundance) within gene = **−0.0801
(p=9.1e-04)**, so conditioning on abundance makes β's gradient **go up**, not down. **"Canonical = abundant" is
FALSE for β.**

**⇒ TEST 3 VERDICT:** β's evidence gradient is **real but small**, and its **edge-specific** component
(the part that is not just "famous family") is **borderline** once family identity is fully absorbed. The
author's original "✅ SURVIVES (b=+0.0506, p=9.8e-03)" must be read as **b≈+0.031, p≈0.07–0.10.**
`shapley`'s gradient and the CN π's "strongest gradient" claim are both **withdrawn** — they are family fame.

### 3.4 TEST 4 — the ceiling ("the data cannot say" is a live partial explanation)

Median posterior SD/|β| = **0.799** (IQR 0.450–0.935); only **27.8% of 5,117 units are identified** (|z|>2).
Typical designs are **not** near-singular (median condition number 1.5; max pairwise |r| between family
predictors median 0.287, >0.8 in only 4.2% of genes) — **but big designs are**: ≥10 families (n=98 genes) →
cond 8.8, max|r| 0.687, only **17% identified** — and that is exactly where the gradient decays (per-gene
ρ(β, w) = **+0.045, p=0.11** on big designs vs **+0.122, p=0.0086** on well-identified genes).
⚠ The moderation test itself is **n.s. (MWU p=0.119)** ⇒ I am **NOT** claiming identifiability significantly
moderates the gradient. **"The data cannot say" and "canonical families genuinely do not get the heaviest
weight" are NOT distinguishable here.**

---

## 4. WHAT SHOULD NOW BE RE-RUN — and what should NOT

**Aneuploidy did NOT buy power. Do not schedule a power-motivated re-run.** What it bought is *calibration*,
and the honest to-do list follows from that:

1. **PERSIST THE PLOIDY CONTROL — for calibration, not for power.** Add a `ploidy` column to the CN
   instrument's confounder block (`confounders.build_C` or an instrument-local `C`). **Prerequisites before it
   lands:** (a) **leave-one-out ploidy** (the focal locus currently contributes 1/506 to its own control — the
   LOO check was NOT run); (b) a **VIF/conditioning check** (`ploidy` is r=0.828 with `target_cn`, already in
   C); (c) decide whether the *official* aneuploidy scores stay at all (they contribute almost nothing:
   81.5% / 82.9% alone vs ploidy's 49.9%).
2. **RIPPLE (canonical completion gate).** If it lands in the shared `C`, **every persisted CN-instrument
   output goes stale**: `exclusion.tsv`, `cn_edge_validation.tsv`, `cn_family_effects.tsv`,
   `cn_family_existence.tsv`, `cn_family_literature.tsv`, `channel_cn`. **The ripple has NOT been traced.**
   ⚠ And note §1.6: the **exclusion gate stops discriminating** under `C_ANEU` — consumers that use it as a
   validity filter change behaviour.
3. **⛔ STOP CALLING `pi_causal` EXOGENOUS.** It is `γ_s · b_fam` = first stage × **observational** slope. Any
   claim that needs exogeneity must use the **reduced form** or a proper 2SLS/LATE ratio — and when it does,
   the HE-vs-decoy separation **does not survive** a locus-level unit of analysis (τ=−0.0015, p=0.67).
   **`MH124_ANTICOUPLING_VALIDITY.md` §5's third existence leg is hereby DOWNGRADED.** MH-124's other two
   existence legs (set-level abundance-matched enrichment OR 1.41; the decoy-controlled model test) are
   **untouched** — but they are both *observational*, so **edge existence currently has no surviving exogenous
   validation.** Restoring one is now the single highest-value open item on the CN arc.
4. **FIX THE UNIT OF ANALYSIS in every HE-vs-decoy test.** 3,111 edges = 482 genes / 409 arms / 277 loci;
   ICCs 0.30–0.37. Use **gene + arm two-way FE** (or cluster at the locus). Per-arm abundance matching is a
   per-arm constant and cannot remove a locus-level confound.
5. **SCALE THE LITERATURE SET — this, not the instrument, is what gates attribution.** n went 16 → 32 here and
   the argmax test stayed at chance while the *rank* test became decisive for abundance/shapley. A versioned,
   auditable ground-truth pull (rather than a hand-curated list) is the prerequisite for any further
   attribution claim. Until then, **"β does not attribute" is measured; "the CN instrument cannot attribute"
   is UNDERPOWERED, not refuted.**
6. **ADD A FAMILY FIXED-EFFECT ROW to any within-gene gradient test.** Gene FE does **not** remove the family
   main effect, and it is worth ~half of the raw gradient.

---

## 5. NOT MEASURED (state these as gaps, not as negatives)

* LOO-ploidy; VIF/conditioning of `ploidy` against `target_cn`; a per-chromosome/arm mean-CN control (which
  would absorb the tide more completely but carries a genuine over-control risk — the instrument's own
  variation *is* arm-level CN).
* Whether the ploidy sign-flip exists in **CPTAC** (everything here is TCGA-only).
* An **arm-resolution** version of Test 1/3 (`readouts_arm_edges.tsv` exists, but §F/§8 doctrine says the arm
  β is a *nomination*, not the estimand).
* Whether the β↔evidence gradient holds **out of sample** (CPTAC).
* Re-derivation of the MH-115/116/117 numbers used in §2 (quoted as registered; a `coupling_grid` field bug is
  under separate investigation).
