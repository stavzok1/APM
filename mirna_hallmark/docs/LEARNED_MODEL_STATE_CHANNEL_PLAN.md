# Plan — the STATE / progression channel (`M_T = a·M_H + Δ`)

> **Goal:** the plan for folding the healthy/progression axis into the §6b posterior as a channel — the
> decisions, the evidence that settled them, the phasing, and the seams with the CPTAC session.
> **What belongs here:** the state channel's design decisions + rationale + phase plan, and the Phase-0
> probe results that justify them (incl. the settled deconvolution/confounder-panel decision, which also
> binds the CPTAC session's C). NOT results of the built channel (→ `LEARNED_MODEL_VALIDATION.md`), NOT the
> abundance-level cross-state arc (already built + documented: `parallel_view.py`, `states.py`,
> `analyses/cross_state/`).
> **Update trigger:** when a phase lands, or a decision here is overturned by a probe.
> **Sync-partner:** `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` (§Q state-as-dimension, §I2 state channel, §7
> catalog row, §E gauge, §J engine, §N confounder block), `PROGRAM_FORWARD_BOARD.md` §C/§F.

> # ⛔ VERDICT (2026-07-12, MH-102d): **THE STATE CHANNEL WAS MEASURED AND CANCELLED. `channel_state.py` was NOT built.**
> On a **calibrated** gauge (split-half control `a=1.125 [1.00,1.28]`, cal✓) under the **E3 standardized-β** convention:
> **τ — the per-edge acquired-vs-constitutive readout, i.e. this axis's ENTIRE deliverable — is indistinguishable from ZERO**,
> and the coupling contribution is **2.4%** of the mRNA likelihood's own precision. §10 has the full post-mortem.
> Everything below §0–§9 is the *plan as written before the measurement*; it is kept as the record of what was tested and why.

**Status: CLOSED (was: PLAN — awaiting review).** All numbers below are read-only probes run
2026-07-12 (scratchpad: `probe_panels.py`, `panel_precision.py`, `settle.py`, `beta_invariance.py`).

---

## 0. What the axis delivers (expectation-setting, up front)

The state channel is **not** a coupling lever. If wiring were conserved and the gauge perfect, adding GTEx's
n=327 to the tumour's n≈1000 shrinks the β SE by ~13% at best — and the healthy β is attenuated and noisier
than that ideal (§2). Every internal lever in this program has landed "immaterial-at-n≈1000"; a weak, attenuated
second sample of the *same* regression will not be the exception.

**The deliverable is the Δ readout:** a per-edge posterior on `Δ = M_T − a·M_H` → {**constitutive** (Δ≈0, wired
the same in healthy breast) / **acquired** (Δ>0, repression built in the tumour) / **lost** (Δ<0, healthy
repression released)}, with honest width, on the same latent the rest of the model uses. Secondarily a
**π/inclusion lift** (healthy repression as *existence* evidence) — where MH-94 found the larger effect lives.

**The axis's structural advantage:** it is the only channel whose exogeneity is **not an assumption**. CN needs
an exclusion restriction (and fails it for ~73% of well-instrumented edges); the healthy data is a **disjoint
cohort** — nothing to exclude. Its threat is cross-cohort commensurability (§2), which is tractable.

---

## 1. Phase-0 evidence — what was probed, and what it settled

### 1a. THE DECIDING TEST — does the confounder panel move β? (n=1444 coefficients)
β_H fit on the 327 paired GTEx donors (canonical estimator: `_prep` + bagged NNLS), 147 genes, under four
confounder blocks. Spearman ρ of the resulting β vectors (on both-nonzero edges in parentheses):

| vs `wu300` | ρ | what it perturbs |
|---|---|---|
| `wu150` — same atlas, 150-cell signature | **0.986 (0.971)** | reference cell-subsample (pure technical replicate) |
| `hbca` — a **completely different atlas** | **0.941 (0.851)** | the reference itself |
| **no C at all** | **0.812 (0.579)** | composition control entirely |

Mean \|β\| **halves** once any composition control is added: `none` 0.0205 → `wu300` 0.0099.

⇒ **The panel choice is immaterial. Having a C at all is everything.** *Which* reference you use moves β ~10×
less than *whether* you control composition. **The "reference reconciliation" that the handoff framed as this
axis's blocking first task is not a blocker.** Any of the three panels works. (§3)

⇒ And the corollary is the real defect on this axis: **`state.cross_state_transfer` fits M_H with no confounder
block at all** — `state._fit_M(Y, Xz)` takes only Y and X (verified by reading the code), while the tumour target
it is scored against *is* C-residualised. So ~half of today's "healthy wiring" signal is composition artifact,
and it also still uses the retired non-negative lasso. **This is Phase 1, and it is a wiring fix, not a run.**

### 1b. The reliability decomposition (what the 150-cell reruns bought)
The 150-cell reruns are a **pure technical replicate** — same reference, same samples, same method, same cell-type
definitions; only the signature's cell subsample differs. That is the **reliability ceiling**: nothing can agree
across references better than a type agrees with *itself* under reference resampling.

| Wu-9 type | 300 vs 150 (GTEx) | 300 vs 150 (NAT) | NAT-alone vs NAT-in-pool |
|---|---|---|---|
| PVL | **0.998** | **0.999** | 0.995 |
| Endothelial / CAFs / Myeloid / Normal Epithelial | 0.99 | 0.99 | 0.90–0.98 |
| Plasmablasts / B-cells / Cancer Epithelial | 0.96–0.98 | 0.89–0.97 | 0.79–0.90 |
| **T-cells** | **0.685** | **0.731** | 0.742 |

**Only `T-cells` is technically unreliable.** Everything else reproduces at ≥0.89, most ≥0.97.

### 1c. Wu-9 vs HBCA-11 disagree on PVL and lymphoid — but that is ANNOTATION, not noise
Same samples, two atlases, coarse lineages: epithelial +0.94/+0.97, endothelial +0.96/+0.96, fibroblast
+0.89/+0.88, myeloid +0.84/+0.87 — but **perivascular −0.46/−0.54** and **lymphoid −0.56/−0.29** (GTEx n=514 /
NAT n=113). Since PVL is *technically* reproducible at 0.998 (§1b), the anti-correlation is the two atlases
**partitioning the mesenchymal and immune compartments differently** — each panel is internally self-consistent.
*(An earlier reading of these as "unmeasured compartments / reference artifacts" was wrong and is retracted, §9.)*

### 1d. Conditioning is fine — kappa was never the problem
kappa(HBCA-11) = **8.1**; kappa(Wu-9) = **12.0**. No signature-column pair exceeds r=0.90 (max: LASP~DDC2 = 0.80).
Nothing in the panels in use is degenerate. **DDC1/DDC2 exist only in HBCA-11 and appear in no Wu-9 run** — and
the whole HBCA-vs-Wu swap moves β by only ρ=0.94, so they are irrelevant to the model. **Nothing to rerun.**

### 1e. C *completeness* for immune is a separate, still-open question
Wu-9's lymphoid tracks Thorsson's lymphocyte-infiltration score at only **r=0.31**, where **LM22 reaches 0.70**
(n≈1050; BayesPrism 0.34; Wu-9-vs-LM22 within-immune share 0.39). §1a proves C's *span is stable across panels*;
it does **not** prove C is *sufficient*. Same shape as MH-100, where C's `mal_prolif` proved materially incomplete.
⇒ a **completeness test**, not a panel question (§4.3). *(Caveat: Thorsson's LIS is itself expression-derived, so it
is a consensus comparator, not ground truth — the fair reading is LM22 vs Wu-9 against the same yardstick.)*

### 1f. Coverage, and the adipose hole
GTEx: 346 miRNA × 514 mRNA donors → **327 paired**, and **all 327 are in `gtex_wu_major`**. 812 arms, of which
**411 carry usable dose variation** (sd>0.3) — so ~half the arm universe is speakable-about at all; the rest get
no channel term (axis P1). GTEx breast is adipose-dominant and **no reference has an adipocyte class** (not Wu-9,
not HBCA-11, not LM22): on HBCA the adipose smears into Fibroblast 0.202 / Vascular-endo 0.166 / PVL 0.126, and
the fit degrades (median R 0.77 vs 0.82 on NAT). GTEx-on-Wu is worse still — 12.6% *plasmablasts* in normal
breast, 31% total immune vs 12% on HBCA — **yet its median fit R is higher (0.843)**, so goodness-of-fit does not
certify fractions. **No deconvolution run can close this; only a marker metagene can.** (§4.2)

### 1g. Dead ends (recorded so no one re-walks them)
- **Local re-solve of CIBERSORTx** from its emitted adjusted signature + mixture (which would have made fine panels
  free and docker-independent): **fails** — overall r=0.68, mean\|diff\| 0.048 vs the tool's own output. Abandoned.
- **A fine 42–49-way panel:** cannot help. β is 0.97–0.99 invariant to reference *detail* and 0.94 invariant to the
  *whole atlas*; a finer C cannot move what a different atlas doesn't move. And the one column a fine immune panel
  would subdivide (`T-cells`) is the only technically-unstable one. **Do not run it.**

---

## 2. The model — why `M_T = a·M_H + Δ`, not `M_T = M_H + Δ`

**The composition-attenuation trap.** GTEx breast is ~32% epithelial (coarse); the tumour is epithelial-dominant.
A repression edge acting in epithelium, measured in bulk, is attenuated roughly in proportion to the epithelial
share — *even under perfectly conserved wiring*. Add a platform change (GTEx TPM vs TCGA RPM — relative scale only,
memory `parallel-view-abundance-vs-wiring`) and non-identical C blocks (§4.1), and β_H sits on a systematically
compressed scale. Under the additive form `Δ = M_T − M_H`, that compression appears as **Δ>0 for essentially every
edge** — it would manufacture the axis's own headline ("acquired wiring") out of a composition artifact. Direction-
known ⇒ design it out, don't check for it.

```
M_T = a · M_H + Δ            a = ONE global scalar, shared by all edges
```
`a` absorbs everything systematic and uninteresting — composition attenuation, the TPM↔RPM platform scale, and the
residual difference between the C blocks. The channel's information is then the **per-edge deviation Δ from the
global healthy→tumour map**, which is the biological quantity we want and is invariant to the panel decision.
The global scale is a nuisance; we never claim it.

- `a` is estimated **once, across edges** (errors-in-variables ⇒ Deming/TLS through the origin, on the frozen
  bagged-NNLS estimates), **out-of-fold w.r.t. the genes it is later applied to**, then frozen — so the channel
  cannot inform the global mean it was calibrated on. `a = 0` (healthy carries no information) is a testable null.
- `τ_Δ` gets an inverse-gamma hyperprior and is sampled (EB method-of-moments as fallback). **`τ_Δ → ∞` ⇒ the
  channel is inert ⇒ bit-identical to §6b** — validation-ladder item (i), for free.
- **`a` cannot absorb *per-edge* differential confounding.** That is why the C blocks must be built the *same way*
  in each state (§4.1) even though `a` handles the global offset.

---

## 3. DECISION — the panel question is CLOSED

**TCGA stays Wu-major. GTEx and NAT use Wu-9 too. No new CIBERSORTx run — zero.**

Three independent reasons, in order of force:
1. **β doesn't care** (§1a): ρ=0.94 across entirely different atlases, 0.99 across reference resamples.
2. **The model's tumour C cannot move to HBCA even in principle:** `data._malignant_prolif` is *defined* as the
   proliferation metagene residualised on the **Cancer Epithelial** fraction — and HBCA, a *normal* breast atlas,
   has no malignant class. A swap would delete the covariate `mal_prolif` is built from.
3. **Ripple:** `_DECONV_COLS` feeds C for every module; a swap invalidates every persisted output.

The handoff's premise ("GTEx/NAT are on HBCA, TCGA is on Wu, there is no `tcga_hbca`") was overtaken by runs that
landed afterwards: **`gtex_wu_major` exists (514 samples, covering 327/327 paired donors).** Reconciliation went the
other direction. **The §C/§F "CIBERSORTx blocker" is retired**, and the CPTAC session's C is untouched — the
"one panel, shared" coordination fight does not exist.

---

## 4. DECISION — the confounder block (this is the real problem, not the panel)

### 4.1 The matched-C policy — the actual "reconciliation"
`CPE`, `target_cn` and `mal_prolif` are **all tumour-only** — they do not exist in healthy tissue. So the C blocks
**cannot be made identical across states no matter which panel is chosen.** This, not the panel, is the genuine
cross-state confounder problem, and it needs an explicit policy:

- **Shared core (all three states):** the 8 non-malignant Wu-9 fractions (`_DECONV_COLS`, Cancer Epithelial held
  out) + a proliferation metagene + the adipose metagene (§4.2).
- **Tumour-only extras** (`CPE`, `target_cn`, `mal_prolif`): kept in the tumour fit, because dropping validated
  confounders to buy symmetry is the wrong trade.
- **The asymmetry is then TESTED, not assumed:** refit β_T under the *reduced, healthy-matchable* C and measure the
  shift. If β_T is stable, the asymmetry is harmless and `a` absorbs the rest. If not, the channel's `pihat` is
  built from the reduced-C β_T and the difference is reported.

### 4.2 The three C blocks — all exist, none need running
**BUILT: `learned/confounders.py` — `build_C(cohort, participants)` + `availability()`.** One builder, four cohorts,
verified **drop-in identical** to both existing consumers (`data._deconv()` max|diff|=0.0 over 1059 participants;
`states._cibersortx_state_cov()` max|diff|=0.0 over 113) ⇒ migrating TCGA and NAT onto it is a **no-op**; only GTEx
and CPTAC actually change — the two that were broken.

| cohort | source | n | note |
|---|---|---|---|
| **tumour** | `tcga_tumor_nat_pooled` (tumour rows) | 1059 | ≡ the wired `tcga_cibersortx_fractions.tsv`. **Unchanged, zero ripple** |
| **NAT** | the **`-NAT` rows of the SAME pooled run** | 113 | ⚠ **NOT** the standalone `tcga_nat_alone`. β_NAT is ρ=0.957 between them (second-order), but the standalone run's decomposition is worse: with no cancer cells to anchor the malignant signature, it misassigns ~10% of each NAT sample into `Cancer Epithelial` (0.033→0.101), eaten out of `Normal Epithelial` (0.245→0.175) — a real C column. *(An earlier recommendation to use `tcga_nat_alone` is RETRACTED, §9.)* The pooled rows are what `states.py` already reads ⇒ no-op |
| **GTEx** | **`gtex_wu_major`** | 514 (327/327 paired donors) | replaces **no C at all** in `state.cross_state_transfer` and the metagene *fallback* in `states.assemble_state` (whose docstring flagged CIBERSORTx-on-GTEx as "the documented future upgrade" — it has arrived) |
| **CPTAC** | **`cptac_wu_major`** | 133 (covers **104/104** of the prospective miRNA arm matrix) | replaces **no composition at all** in `cptac_validation` |

**Purity stays cohort-native and is NOT harmonised** — TCGA `CPE`, CPTAC `ESTIMATE`, NAT/GTEx **none**. Justified:
purity is nearly redundant given composition (β ρ=0.984 with vs without CPE), and CPE is *not* the malignant fraction
(r=0.554, means 0.732 vs 0.235 — not substitutes). ⚠ NAT/GTEx must have **no** purity term: the clinical table is
*participant*-keyed, so a naive join silently hands **107/113** NAT samples their own patient's **tumour** CPE — a
covariate describing a different specimen. `build_C` returns `None` to block it.

**The adipose metagene is DEFAULT-OFF, gated on the over-control test — not mandatory** *(revised; an earlier
"mandatory" is retracted, §9)*. Measured on the 514 GTEx samples: the 8 Wu fractions already explain **69.6%** of the
adipose axis, and the metagene is **r=−0.71 with `Normal Epithelial`** / −0.54 with `Cancer Epithelial` — i.e. largely
an **inverse-epithelial axis**. Residualising on it risks removing the compartment the target is expressed in, the exact
failure mode `data._latent`'s global PCs hit. Only 30.4% of it is new information. Gate it (P0b), don't assume it.

### 4.3 The immune-completeness test (the one open deconv-adjacent item)
Does C actually control immune composition, or merely appear to? Score candidate immune covariates in **tumour**,
where two independent yardsticks exist (**LM22** `hires_out`, 1106 samples; **Thorsson** immune scores):
Wu-9's four immune columns vs an immune marker metagene vs LM22-derived immune. If C under-controls, **the fix must
be the metagene** — `states._state_metagene_cov` already contains one — because it is the only immune covariate
computable identically in GTEx/NAT/tumour. **LM22 is a validator, never a covariate** (tumour-only; and its
fractions are immune-*relative*, not bulk). Shared item with the CPTAC session.

---

## 5. DECISION — PRIOR form (a plug-in channel term), not a 2nd likelihood

The prior form **is already a term in the existing `channels=` API** — no sampler changes:
```python
{"members": [(j, 1.0)],                     # loading 1 on family j's own β  (cf. channel_cn's γ)
 "pihat":   beta_H[j] / a,                  # healthy β, mapped onto the tumour gauge
 "s2":      se_H[j]**2 / a**2 + tau_d**2}   # its SE + the acquired-wiring spread
```
`_gibbs_posterior` adds `N(π̂; β_j, s²)` — Gaussian-conjugate ⇒ **stays on Gibbs (J1)**.

The joint 2nd-likelihood form co-estimates β_H and β_T from both cohorts' raw data. Because the cohorts are
**disjoint samples**, the two-stage plug-in is nearly statistically equivalent; the joint form only recovers the
healthy posterior's non-Gaussianity (the half-normal spike at 0), which is second-order — for real sampler surgery.
**Recommend the prior form**, with the joint form as a gated fallback if most β_H pile on the spike.

**Engine:** this channel never forces HMC. If CPTAC moves the posterior to J2, this term ports as a one-line
Gaussian prior. They should not wait for us; we will not force J2 on them.

---

## 6. Phases

**Phase 0 — the shared infrastructure (MH-102). ✅ LARGELY DONE (2026-07-12) — no CIBERSORTx runs.**
- ✅ **P0a — the canonical C builder.** `learned/confounders.py` — `build_C(cohort, participants)` + `availability()`
  for TCGA · NAT · GTEx · CPTAC. Verified **drop-in** (max|diff|=0.0 vs *both* `data._deconv()` and
  `states._cibersortx_state_cov()`) ⇒ TCGA/NAT migration is a **no-op**; only GTEx and CPTAC change.
- ✅ **P0b — the adipose metagene: GATED OFF** (§4.2). 69.6% explained by the existing Wu fractions, r=−0.71 with
  `Normal Epithelial` ⇒ over-control risk. Available via `adipose=True`; does not fire by default.
- ✅ **P0e — the gauge (axis E).** `learned/gauge.py` — `fit_gauge` / `a_for` (out-of-fold) / `to_channel_terms` /
  `beta_table` / `shuffled_gauge`. **The ladder (all three PASS their shuffled null):** CPTAC **a=2.369** ·
  NAT **a=0.528** · GTEx **a=0.231**. **It factorizes** — state ≈2× × platform ≈2.3× ≈ 4.3 ≈ 1/0.231, which *is*
  the NAT platform-bridge decomposition. Two estimator bugs were caught **only** by the null (registry MH-102b).
- ✅ **P0d — the immune-completeness test.** `learned/immune_completeness.py`. C is immune-**incomplete** (LM22-lymphoid
  0.505 / IFNG-response 0.763 of variance survives C) but the fix is **INERT** (β ρ=0.967, |β| −5%, same on immune and
  non-immune targets) ⇒ **C UNCHANGED**. Closes the last open C question for *both* channels.
- ✅ **The cross-cohort shuffle-null primitive** — `gauge.shuffled_gauge()`, shared (analogue of `exclusion(shuffle_cn=)`).
- ✅ **P0c — the matched-C policy test.** `CPE`/`mal_prolif` are tumour-only, so the C blocks cannot be column-identical.
  Refitting β_T under the reduced, healthy-matchable C (composition only) gives **ρ=0.900 (0.837 both-nonzero)**,
  mean|β| 0.0072→0.0070. So the asymmetry is **real but second-order** — the same magnitude as an atlas swap (0.94) or
  a refsample swap (0.92). **It is already ABSORBED**, not ignored: the gauge is fit *source → β_T-under-the-FULL-C*
  (the model's actual latent), so the global part lands in `a` and the per-edge part in `τ`, which inflates `s²` and
  makes the channel self-down-weight. **No reduced-C variant of β_T is needed.**
- ✅ **P0f — coverage / identifiability.** Each source can speak about ~**52% of target genes** and **54–59% of target
  edges** (NAT 138 genes/1282 edges · GTEx 137/1178 · CPTAC 137/1251). The rest get **no term** and revert to
  mRNA-only (axis P1).
- ✅ **The last reference mismatch is CLOSED** (user added `cptac_wu_alt`, 2026-07-12): **all four cohorts now sit on
  REF-A** (300-cell). Measured immaterial, as predicted — CPTAC fractions agree r=0.88–0.996 between REF-A and REF-B
  with identical fit R (0.866/0.864), and the gauge moves only **a: 2.369 → 2.314**. `build_C` now sources CPTAC from
  `cptac_wu_alt`; `cptac_wu_major` (REF-B) is retained as a sensitivity run.

**Phase 0 is CLOSED.** Remaining before the channel:
- ⬜ **Migrate the three call sites onto `build_C`.** TCGA/NAT are verified no-ops; **GTEx and CPTAC are not** —
  `cptac_validation` currently has *no composition term*, so its prospective results (the 34 FDR-negative edges)
  **will move**. That is a correction, not a regression, but it is a re-run with a ripple: trace it, don't fold it in.

**Phase 1 — retire the broken healthy fit.** Rebuild `state.cross_state_transfer` on the canonical estimator with a
real C_GTEx (kills the lasso *and* the no-C fit, §1a). Its current Δ=0/retention verdict is composition-confounded
and must be re-derived before it is cited anywhere.

**Phase 2 — the channel** (`learned/channel_state.py`, mirroring `channel_cn.py`): `state_terms(gene, cols, members)`
→ the `channels=` list. **Validation ladder, all three:** (i) **nests §6b** — `τ_Δ→∞` ⇒ bit-identical; (ii) **gauge**
— single strong clean edge, channel-only β == mRNA-only β; (iii) **shuffled null** — permute GTEx donor labels
between the miRNA and mRNA matrices (`state_terms(shuffle_state=True)`, mirroring `exclusion(shuffle_cn=)`) ⇒ β_H→0
⇒ the channel's weight must collapse.

**Phase 3 — the readout (the deliverable).** Per-edge Δ posterior → {constitutive / acquired / lost} with honest
width; joined into the per-edge attribution card (board §D). Then the **π/inclusion half** (healthy repression as
existence evidence, mirroring MH-94's `between_family_bayes`) — no double-count with the curated evidence prior
(axis M): that prior is *curation*, this is *data*, from a disjoint cohort.

**Phase 4 — fold in the trajectory infra (conditional, low priority).** **Keep** `parallel_view.change_trajectory` /
`state_pressure_attribution` / `states.budget_shift` as readouts — they are mature; only their M source changes. Use
`states.realization` (111 within-patient NAT→tumour pairs, **baseline-free** — the cleanest causal test in the
subproject) as the channel's **out-of-sample validation**, not as a channel.

**Gated OUT:** **Q2** (state as a *latent dimension* with a smoothness prior) — per §Q's own diagnostic, only if
Phase 3 shows ΔM is real. **Fine subtype panels** — §1g, cannot help C.

---

## 7. Coordination with the CPTAC / protein session
1. **Composition: resolved, and it costs them nothing.** TCGA C stays Wu-major (§3). No shared-panel decision is
   pending; no run is needed. The one shared item is the **immune-completeness test** (§4.3/P0d) — if C
   under-controls immune, it affects protein discordance too (often composition-driven).
2. **Engine:** this channel is Gaussian-conjugate (J1). Theirs may force J2 — **their call alone**; our term ports
   trivially either way.
3. **(protein × state)** — the natural joint object if both channels land (is protein discordance state-dependent?).
   **Flagged, not built**, by either session solo.
4. **Gauge:** each session runs its own single-edge check; ours is the `a`-constant of §2.

## 8. Runtime (the deliverable-not-afterthought axiom)
The healthy fit is one extra 327×~10 regression per gene — trivial. The cost is **loads**: the GTEx mRNA GCT and
miRNA matrix hoisted/cached once (`state._CACHE` already does), and the frozen β_H table persisted
(`output/learned/state_beta_healthy.tsv`) so the channel is a lookup, not a refit. Batch over ~1500 genes → profile
first; parallelize as `prolif_verdict` does (8 workers).

## 9. Retractions (claims made during planning and withdrawn on measurement)
Recorded per the `verify-before-assert` axiom — these were asserted before being checked, and the checks refuted them:
- ❌ *"The cross-solver disagreement is concentrated in near-degenerate compartments."* → kappa(HBCA-11)=8.1,
  kappa(Wu-9)=12.0, no signature pair >0.90. **Nothing is degenerate.**
- ❌ *"Perivascular and lymphoid are unmeasured / reference artifacts."* → PVL's technical reproducibility is
  **0.998/0.999**. The Wu-vs-HBCA anti-correlation is an **annotation/partition difference between atlases**, not
  noise. Only `T-cells` is genuinely unstable (0.685/0.731).
- ❌ *"Fine immune subtypes are unreproducible, based on GTEx+NAT."* → both are normal breast (one regime tested
  twice, not two). Tested properly in tumour instead; the fine-panel verdict now rests on §1a/§1g, not on this.
- ❌ *"Reference reconciliation is the axis's first task."* (the handoff's premise) → **retired**: `gtex_wu_major`
  exists, and β is panel-invariant anyway (§1a/§3).
- ❌ *"Use `tcga_nat_alone` for the NAT C, not the pooled rows."* → **reversed on measurement.** That advice rested on
  fraction shifts and a subspace rotation; β_NAT was never tested. It is ρ=0.957 between the two sources, while the
  standalone run's malignant negative control degrades 3× (0.033→0.101) at the direct expense of `Normal Epithelial`,
  a real C column. **Use the pooled `-NAT` rows** — which `states.py` already did.
- ❌ *"An adipocyte metagene is mandatory for GTEx."* → **default-OFF**: it is 69.6% explained by the existing 8 Wu
  fractions and is r=−0.71 with `Normal Epithelial` (an inverse-epithelial axis) ⇒ genuine over-control risk. Gated.
- ❌ *"Harmonise the refsamples with 3 CIBERSORTx re-runs."* → **no runs.** The filenames lie; there are only TWO Wu
  references (md5 of the inferred refsample), the production mismatch is CPTAC-alone, and it costs β ρ=0.986 —
  absorbed by the gauge constant + s² floor. Also: S-mode re-adjusts per mixture regardless, so identical refsample
  files would *not* have produced identical signatures anyway.
- ❌ *"CPTAC needs a marker-metagene C (no deconvolution exists)."* → superseded: `cptac_wu_major` (133) now exists and
  covers **104/104** of the prospective miRNA cohort. It also makes `mal_prolif` computable in CPTAC (it needs the
  Cancer-Epithelial fraction). Metagenes survive only for adipose (gated) and possibly immune-completeness.


---

## 10. POST-MORTEM — why the axis is closed (MH-102d)

> ⚠ **NUMBERS BELOW SUPERSEDED — **FINAL (both sd-floors scale-free): GTEx a=0.137 ρ=0.188 τ=0.0009 info=0.6% · NAT a=0.326 ρ=0.277 τ=0.0221 info=0.7% · CPTAC a=0.333 ρ=0.280 τ=0.0023 info=0.7%.** Conclusion UNCHANGED (τ≈0, info ~0.6–0.7%) ⇒ state channel stays CANCELLED. The two sd-floors (`min_sd=0.2` on Y, `sd<0.1` on X) were ABSOLUTE and the cohorts are on different scales; found by the CPTAC session. The verdict is unchanged.**


**The verdict.** Under a gauge that passes its own split-half control (`a=1.125 [1.00,1.28]`) and the **E3
standardized-β** convention:

| source | a | ρ | **τ (the payload)** | **info added to β_T** | shuffle-null |
|---|---|---|---|---|---|
| *split-half TCGA (control)* | *1.125* | *0.72* | — | *124%* | — |
| **GTEx** | 0.112 | 0.169 | **0.0008** | **2.4%** | dead ✓ |
| **CPTAC** (mRNA) | 0.336 | 0.313 | **0.0024** | **0.6%** | dead ✓ |
| **NAT** | 0.327 | 0.287 | 0.0268 | **0.4%** | ⚠ **NOT dead** (a_null=0.076) |

**1. No payload.** τ is the spread of per-edge deviations — the acquired-vs-constitutive wiring readout, and the
*only* thing this axis was ever going to deliver. It is ~0 for GTEx and CPTAC; a simulated τ_true=0 world produces
*larger* τ̂ by chance. NAT's apparent τ=0.027 is retired: it survives a matched-C refit (so it is not the C
asymmetry) but it **vanishes under the Gibbs estimator** (bagged NNLS's boundary zeros give `se=0`, understating
`mean(se²)` and hence overstating τ) **and its shuffle-null fails**.

**2. No coupling gain, for a structural reason.** The channel's precision scales as **a²**, and `a≈0.11` ⇒
`a²≈0.012`. Even with unlimited GTEx donors it would add ~1%. **Composition attenuation does not merely bias the
healthy β — it destroys its information content.**

**3. The handoff's form would have faked the headline.** With β_H ≈ 0.11·β_T, the additive `Δ = M_T − M_H` gives
Δ>0 on nearly every edge ⇒ "acquired wiring" manufactured out of attenuation. The `a·M_H` correction was necessary;
it just revealed there was nothing underneath.

**4. Two estimator pathologies, caught by the CONTROLS and not by reasoning.**
- Fed the **Gibbs posterior's** β/SD, the split-half control returns **a=4.1** (3.5e6 in one configuration) where the
  truth is 1.0. ⚠ **CORRECTED FRAMING — this is NOT a defect in Gibbs.** Gibbs is the **BETTER estimator** (split-half
  reproducibility ρ=**0.822** vs bagged NNLS's **0.729**) and **the model keeps it**. The fault is *this gauge's*
  errors-in-variables correction: it divides by `Var(β̂)−mean(se²)`, and Gibbs's `mean(se²)` is dominated by a **heavy
  tail of a few enormous posterior SDs** (`sqrt(mean(se²))=0.055` vs a typical se of 0.015) ⇒ reliability **0.17** vs
  NNLS's **0.72**. `_MIN_RELIABILITY` + `Gauge.usable` now **refuse** the gauge instead of returning 4.1.
  **Bagged NNLS for the GAUGE; Gibbs for the MODEL.**
- The **raw-`r` convention silently absorbed the Y-scale**: median sd(resid Y) is **TCGA 0.237 vs NAT 0.600 /
  GTEx 0.612** — normal tissue is ~2.5× more variable — and `a_raw/a_zY = 2.25 ≈ the sd ratio 2.58`. So the earlier
  "healthy β is 23% of tumour β = composition attenuation" was **mostly Y-scale**. Fixed by **E3** (`zscore_y=True`,
  CHANNEL_FUSION §E — which the design doc had anticipated and this plan skipped).

**5. The one estimator-robust quantity** is the rank concordance ρ: split-half TCGA **0.72** (the reproducibility
ceiling) · NAT 0.27–0.39 · GTEx 0.10–0.18. Healthy β shares ~a fifth of the reproducible per-edge structure with
tumour β. That is a *relative* statement and cannot found a channel — a channel needs a **scale**.

**Independent convergence.** The parallel **CPTAC/protein session** reached the same structural conclusion by a
different route: protein carries **1.7%** of the mRNA channel's Fisher information about β (≤8.8% even at a perfect
slope) ⇒ "protein cannot move β." Both sessions then pivoted the same way: **the value of an exogenous source is a
NEW LATENT (a different question), never a coupling gain on the same β.** (→ `CPTAC_PROTEIN_CHANNEL_PLAN.md`)

**What survives.** `learned/confounders.py` (build_C) and `learned/gauge.py` (calibrate / fit_gauge / shuffled_gauge
/ info_ratio) are the durable output — **they are what killed the axis**, which is what falsifiable infrastructure is
for. Any future cross-cohort channel should be run through `calibrate()` → `info_ratio` → `τ` **before** it is built.


---

## 11. TWO MODEL-WIDE SPIN-OFFS (from the split-half control, 2026-07-12)

**(a) CALIBRATION — the model's posterior widths are ~40% too narrow.** Measured against independent half-cohorts, the
reported uncertainties **understate the true sampling variability**: bagged NNLS by **28%** (se 0.0055 vs true sampling
sd 0.0079), the dense Gibbs posterior by **39%** (0.0145 vs 0.0237). And the se distribution is **heavy-tailed**
(`sqrt(mean(se²))=0.055` vs a typical 0.015), so a minority of edges are *under*-confident while the typical edge is
*over*-confident — a two-sided failure. This feeds **PIP**, the **Bayesian-Shapley widths** (MH-94), **δ-pooling
confidence** (MH-98), and **every channel's `s²`**. ⇒ the board's **§E BN simulation-based calibration** item is no
longer hypothetical and has been **promoted**.

**(b) THE β GAUGE CONVENTION — DECIDED: keep raw-`r` for the model; z-Y is GAUGE-ONLY.**

| | global ρ(A,B) | top-25 | top-50 | **top-100** |
|---|---|---|---|---|
| **raw-`r`** (current) | 0.728 | **12**/25 | 25/50 | **63**/100 |
| z-Y | 0.728 | 11/25 | 25/50 | 55/100 |

Within-gene family rankings are **exactly invariant** (ρ=1.0000) — z-scoring Y rescales all of a gene's β by one
constant. Cross-gene is near-invariant (ρ=0.992). The decision rests on **reproducibility**: raw-`r`'s discovery
ranking is more stable across independent halves, because z-Y divides β by `sd(resid Y)` (a **7.2×** spread across
genes) which is smallest exactly where the residual is most noise-dominated.
⚠ **Honest caveat:** raw-`r`'s top-50 **is** tilted toward high-variance genes (median sd(residY) **0.369** vs 0.243
overall). z-Y removes the tilt but costs stability. If a variance-neutral view is ever wanted, z-Y is an alternative
*view*, not a replacement.
⚠ **The cross-cohort GAUGE must still use z-Y** (`gauge.beta_table(zscore_y=True)`): there the Y-scales differ
*systematically* (TCGA 0.237 vs GTEx 0.612) and the whole factor lands in `a`.
⚠ **Discovery-lane reality check:** top-100 edge reproducibility across half-cohorts is only **~63%** under EITHER
convention.
