# CPTAC / PROTEIN channel — the plan (βᵗ, the translational-repression latent)

> **Goal:** the agreed, evidence-checked build plan for the CPTAC/protein axis — what gets built, what it
> yields, which models it uses, what it inherits vs. adds, and the prerequisites it is blocked on.
> **What belongs here:** the protein-axis plan, its verified feasibility numbers, its phase gates, and the
> retractions/corrections made while forming it. NOT results (→ VALIDATION/registry once run), NOT the general
> fusion frame (→ CHANNEL_FUSION_DESIGN), NOT the state axis (→ HANDOFF_PROGRESSION_STATE).
> **Update trigger:** when a phase lands, a gate fires, or a prerequisite is unblocked.
> **Sync-partner:** `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` (§7 protein row, §I1/§L/§J/§E — **this doc
> CORRECTS §7's conjugacy and engine claims**), `handoffs/HANDOFF_CPTAC_PROTEIN.md` (the orientation this answers),
> `PROGRAM_FORWARD_BOARD.md` §A③/§E.

**Status:** PLAN, user-approved 2026-07-12. Not started. ~~Blocked on P1 (CIBERSORTx on CPTAC).~~

> ✅ **P1 IS UNBLOCKED (state session, 2026-07-12).** CPTAC is deconvolved: **`data/external/cibersortx/cptac_wu_alt/`**
> (133 samples, Wu-major **REF-A** — the *same* reference as TCGA/NAT/GTEx; covers **104/104** of the prospective miRNA
> cohort). Build C with **`learned.confounders.build_C("cptac", parts)`** — one builder, four cohorts, verified drop-in
> (max|diff|=0.0 vs `data._deconv()`). It also makes **`mal_prolif` computable in CPTAC** for the first time (it needs
> the Cancer-Epithelial fraction). ⚠ Wiring composition into `cptac_validation` **will move** its prospective results
> (the 34 FDR-negative edges) — its C currently has **no composition term at all**. That is a correction, not a
> regression, but it is a re-run with a ripple.
>
> 🔁 **INDEPENDENT CONFIRMATION of this plan's §1 headline (state session, MH-102d).** Measured on a *calibrated*
> cross-cohort gauge (`learned/gauge.py`; split-half control a=1.125 ✓): **CPTAC-mRNA→TCGA carries 0.6%** of the mRNA
> likelihood's precision on β (GTEx 2.4%, NAT 0.4%) — i.e. the *easiest possible* cross-cohort transfer (same modality,
> same disease state) already caps out at ~1%, converging with this plan's Fisher-information figure (⛔ **the 1.7%
> once cited here is VOID — PAM50-contaminated C; §1's re-measured value is ≈4–6%, ceiling ≤7.6%**) by a different
> route. **The state channel was cancelled on this evidence** (registry MH-102d). It corroborates: an exogenous source's
> value is a **new latent** (`βᵗ`), never a coupling gain on the same β.
>
> ⚠ **Two GAUGE pathologies to avoid** if you estimate any cross-cohort slope — note both are about the *gauge*, not the model:
> (a) **use bagged NNLS to fit the gauge, NOT the Gibbs posterior.** Fed Gibbs's β/SD the split-half control returns a=4.1
> where the truth is 1.0 — not because Gibbs is wrong (**it is the BETTER estimator**: split-half ρ=0.822 vs NNLS 0.729, and
> the model KEEPS it) but because the errors-in-variables correction divides by `Var(β̂)−mean(se²)`, and Gibbs's `mean(se²)`
> is dominated by a heavy tail of huge posterior SDs (reliability 0.17 vs 0.72). Run `gauge.calibrate()` and honour
> `Gauge.usable` — it now refuses rather than returning a mis-scaled `a`.
> (b) fit the gauge on the **standardized-β convention** (`zscore_y=True`, the default in `gauge.beta_table`), or the slope
> silently absorbs the cohort's Y-scale (median sd(resid Y): TCGA 0.237 vs GTEx 0.612). ⚠ This is a **gauge-only** convention
> — the MODEL stays on raw-`r` (tested: within-gene rankings are exactly invariant, and raw-`r`'s cross-gene ranking is
> *more* reproducible, 63/100 vs 55/100 top-100 split-half stability).
>
> 📌 **A model-wide calibration finding you should factor into your `βᵗ` power calc:** the reported uncertainties
> UNDERSTATE the true sampling variability — bagged NNLS by **28%** (0.72×), Gibbs by **39%** (0.61×), measured against
> independent half-cohorts. So posterior widths (and any power calculation resting on them) are ~40% too narrow. At n≈101
> this matters for βᵗ. (→ board §E, SBC.)

---

## ⛔ 0a. PHASE-0 VALIDATION FALSIFIED THE PLAN'S CENTREPIECE (2026-07-12) — READ FIRST

**`βᵗ` as a per-family latent is NOT SUPPORTED at n=101. The evidence I built this plan on was a LEAK.**

**The leak.** Every earlier βᵗ result in this doc (OOF 9/9 positive; PTEN z=3.1/2.7/2.3; PEBP1 z=3.6/2.5) built the
target as `protein_resid = P − a·M` with **`a` fit on ALL samples**, then OOF-ed `X → protein_resid`. The mediator
was fit using the test folds. **Fit `a` INSIDE the fold and the signal disappears.**

**The leak-free test** (does `X` improve out-of-fold prediction of `P` **beyond the gene's own mRNA `M`**?
null = permute `X` only), run with the **model's own estimator** (NNLS, non-negative — correctly regularized for
p≈63 > n=101):

| estimator / predictor set | `d_incr` > 0 | median `d_incr` | **z > 2** |
|---|---|---|---|
| ridge, all families, best-λ **chosen on the test metric** (optimistic) | 5/10 | **+0.001** | — |
| **NNLS (the model's own), all ~60 families** | 3/10 | **−0.0088** | **0/10 — not even PTEN** |
| NNLS + **top-3 families by TCGA β** (non-circular prior; K post-hoc) | 6/10 | +0.0062 | **3/10: BCL2 3.62 · CDKN1A 2.30 · PDCD4 2.23** |

**What this establishes:** at CPTAC's n, the miRNA family doses do **not** improve out-of-sample prediction of
protein beyond the transcript. It is a **`p` problem** — a strong non-circular prior (top-3 TCGA families) recovers
a weak signal (d = 0.01–0.05) in **3/10** genes — but the full per-family βᵗ posterior the plan proposed (width,
PIP, Shapley, genome-wide discovery) **would be prior-dominated noise.**

**⚠ PTEN — cited as the headline throughout this doc — does NOT survive (z=0.53).** BCL2 / CDKN1A / PDCD4 do.

**What this does NOT overturn:** MH-34/35's marginal **association** results (aggregate pressure ↔ `protein_resid`,
11–29/1132 genes FDR; MH-40's decoy control). Association ≠ out-of-sample prediction; a weak real association can
fail to yield predictive gain at n=101 with p≈60. The biology of translational repression is not in question — our
**power to resolve it per-family at n=101** is.

### The verdict is SETTLED — four independent leak-free framings, including the maximally-powered one

The per-family failure could have been a `p>n` artifact. It is not. The **AGGREGATE** frame spends **1 degree of
freedom, not 60** — the predictor is `X_cptac · β_TCGA`, i.e. **weights FIXED from another cohort AND another layer,
zero parameters fitted in CPTAC**. It is the most powerful test available, and it is what MH-34 actually did:

| framing (all leak-free) | result |
|---|---|
| ridge, all ~60 families, λ swept (best-λ picked **on the test metric**) | median gain **+0.001** |
| **NNLS (model's own), all ~60 families** | **z>2 in 0/10** |
| NNLS + top-3 families by TCGA β (K post-hoc ⇒ optimistic) | z>2 in 3/10 (d = 0.01–0.05) |
| **AGGREGATE — TCGA-fixed weights, 1 df, BH-FDR over 17 genes** | **q<0.10 in 1/17 — only BCL2** (d=+0.028, q≈0.085, marginal). **PTEN: d=−0.006, p=0.82.** |

**⇒ NOT a leak, NOT an estimator choice, NOT a p>n artifact. It is the TRUE EFFECT SIZE meeting the cohort size.**

**And this is exactly what our own biology already said.** MH-34 established repression is **predominantly
mRNA-mediated**, with a translational residual in **11/1132 genes ≈ 1%**. On a 17-gene panel that predicts ~0.2
hits; we found 1. **There was never a βᵗ *field* to fit.** βᵗ is real (Bartel/Selbach; MH-35's survivors) but ~1%
prevalent, and **n=101 resolves approximately one gene's worth of it.** The plan's premise — "a new dimension of M"
— mistook a 1%-prevalence phenomenon for a latent field. **No method fixes this; only n does.**

---

## ✅ 0b. PHASE-1 RESULTS — the surviving axis, RUN (MH-104, 2026-07-12)

### (A) BAR-5 SYSTEMATIC — **PASSES**, with a proper ceiling. And the loss is at the COHORT boundary, not the LAYER.
`β` fit on TCGA mRNA only (never sees CPTAC, never sees protein), scored on 581 genes. The essential control I
initially omitted: **fit β on TCGA half-A, score on TCGA half-B** — the same-cohort out-of-sample *ceiling*, which
tells us whether a flat cross-cohort result means "β doesn't transport" or "my aggregate is broken."

| rung | median ρ | % repression-directed | binomial p | **retention vs ceiling** |
|---|---|---|---|---|
| **(i) TCGA held-out mRNA** — same cohort, same layer (the CEILING) | **−0.117** | **88.8%** | **2e-88** | 1.000 |
| (ii-a) CPTAC mRNA — **STAR** | −0.017 | 54.6% | 1.5e-2 | 0.148 |
| (ii-b) CPTAC mRNA — **LinkedOmics** | −0.023 | 56.1% | 1.8e-3 | **0.193** |
| **(iii) CPTAC protein** — cross-cohort **+** cross-layer | **−0.023** | **57.7%** | **1.3e-4** | **0.201** |

**Three findings:**
1. **The aggregate is NOT broken** (88.8% negative in TCGA held-out, p=2e-88) ⇒ the flat cross-cohort result is REAL.
2. **BAR-5 PASSES**: the learned M carries to an **independent cohort AND an independent layer** (p=1.3e-4). This is
   the honest, FDR-controlled version of MH-83's hand-picked 7-gene "7/7".
3. ⭐ **THE LOSS IS ENTIRELY AT THE COHORT BOUNDARY, NOT THE LAYER BOUNDARY.** Crossing cohorts costs **~80%**
   (0.117 → 0.023); crossing mRNA→protein costs **~0%** (0.193 → 0.201, i.e. protein ≈ mRNA). This **independently
   re-derives the βᵗ falsification from a different direction** (protein adds nothing beyond the transcript) **and**
   it locates the real barrier of the whole cross-cohort programme. Corroborates the state session's 0.6% (MH-102d).

### (A2) ⭐ STRATIFIED — the ~80% loss is a MEDIAN hiding a MIXTURE. Transfer is HETEROGENEOUS and PREDICTABLE. (MH-106)
The flat retention figure in (A) is a median over all genes and it **understates what actually happens**. Stratifying
(587 genes; `n_fam` distribution: median **3**, IQR 2–6, range 2–63):

| stratum | n | retention vs ceiling | % repression-directed | binom p |
|---|---|---|---|---|
| **2–5 regulator families** | 386 | **0.082** | 52.6% | 0.17 **(ns)** |
| 5–10 families | 124 | 0.159 | 59.7% | 1.9e-2 |
| **10–20 families** | 48 | **0.287** | **68.8%** | 6.6e-3 |
| **in-cohort ceiling ρ < −0.25** (strongest) | 80 | **0.217** | **72.5%** | **3.5e-5** |
| in-cohort ceiling ρ > −0.08 (weak) | 150 | ~0 | 51.3% | ns |
| **`a_g` < 0.15** (protein DECOUPLED from mRNA) | 50 | **−0.139** | **42%** | ns — **NO transfer** |
| **`a_g` 0.55–0.80** (protein TRACKS mRNA) | 118 | **0.392** | **66.1%** | **3.0e-4** |

**MIXTURE, not a shift:** observed sd(ρ_protein) **0.136** vs **0.100** under a pure null at n=101 ⇒ **1.36× excess
dispersion**; and corr(ceiling, ρ_CPTAC_protein) = **+0.176 (p=1.8e-5)** ⇒ **a real SUBSET of genes transfers.**

⛔ **RETRACTED (MH-109) — "the miRNA effect reaches protein only where protein tracks mRNA" DOES NOT HOLD.**
I called this "the mechanistic clincher". It was an artifact of coarse binning + an unstable ratio-of-medians
"retention" statistic. **The continuous test says nothing:** on the 204 genes with a usable CN instrument, Spearman
with protein-repression strength is **a_OLS +0.010 (p=0.88) · a_IV +0.032 (p=0.65) · λ −0.003 (p=0.97)**; partial
`a_IV | λ` = +0.031 (p=0.66). On all 490 genes `a_OLS` gives ρ=+0.103 (p=0.022) — **but instrument strength F (pure
CN variance, a total nuisance) predicts equally well (ρ=+0.114, p=0.011)** ⇒ not interpretable as biology.
**UNSUPPORTED, not disproven:** the miRNA→protein signal is AT THE NOISE FLOOR (median ρ=−0.02 vs a per-gene SE of
~0.10 at n=101) — you cannot detect a MODULATION of a signal indistinguishable from zero.

**⇒ CORRECTED CLAIM:** it is **NOT** that "β doesn't transport". **β transports for the genes where β is
WELL-DETERMINED** (many regulators · strong in-cohort coupling · unbuffered protein) — **retention 0.22–0.39** —
**and not at all for the poorly-determined majority.** The cohort barrier is real but **heterogeneous and
predictable**, not uniform.

### (B) ⚠ RNA-SOURCE DEFAULT CORRECTED — `linkedomics`, not `star`. My earlier call was WRONG.
Two independent non-circular measurements both favour **LinkedOmics**: (1) agreement with the *independent protein
assay* (0.3713 vs 0.3678); (2) **Bar-5 rung (ii)** — retention of the TCGA ceiling **0.193 vs 0.148** (p 1.8e-3 vs
1.5e-2). The `star` default rested on a **gauge-commensurability** argument that is now **moot** (`zscore_y=True`
standardizes Y anyway, and the cross-cohort β pooling it served is cancelled). `star` remains the ONLY source usable
by CIBERSORTx (needs absolute linear TPM) and is retained for the both-ways robustness check.
`cptac_data.mrna()` default flipped; **this also explains the "protein beats mRNA" anomaly** — STAR was simply the
noisier mRNA; with LinkedOmics, protein (0.201) ≈ mRNA (0.193), as biology requires.

### (C) V1 — DOSE-DELIVERY TRANSPORTS. The plan's "inherited δ" assumption is now VERIFIED, not assumed.
Member dose-share (which member carries the family's dose = δ's primary input), TCGA n=1041 vs CPTAC n=101, on the
**58 multi-member families expressed in both** (canonical floor RPM≥10):
- **Same dominant member: 84.5%** vs **43.7% chance** ⇒ **1.93× enrichment**
- **median member-share correlation 0.991**; median total-variation distance **0.145**
- (all 260 families incl. noise-dominated: 73.5% / TV 0.230 — the degradation is entirely in undetectable arms)
⇒ **δ's abundance component transports.** (Scope: this validates the *abundance* input to `delta_pooling`, not its
full mRNA-width + CN + chimeric fusion.) `output/learned/v1_dose_delivery_transport.tsv`.

### (D) The CPTAC CN instrument is USABLE BUT WEAK — V1's full-δ leg is 5-family-limited.
First stage (locus CN → arm dose), F: **TCGA (n=961) median 3.16, 234 arms F>10 (32%)** vs **CPTAC (n=101) median
1.10, 59 arms F>10 (8.6%)**; **49 arms F>10 in BOTH**. First-stage γ replicates across cohorts at Spearman **0.349**.
Families resolvable (≥2 distinct loci at F>10): **TCGA 31 · CPTAC 5 · BOTH 5** (let-7 · miR-15/16 · miR-17~92 ·
miR-181 · miR-25/92). ⇒ the full CN-fused δ comparison is feasible but **low-powered**; the abundance-based test (C)
is the well-powered one and is what we report.

---

### PIVOT — what the protein axis IS (all four survive, none depend on βᵗ)
1. **`locus_cn_cptac`** — BUILT + VALIDATED (median **r=0.997** vs ASCAT; coordinate hole patched, 482/482 loci;
   CPTAC 122×414, median CN 1.99). An independent second CN-instrument cohort.
2. **V1 — δ transportability.** The first direct test of the cross-cohort assumption the whole program rests on.
   Needs (1), not βᵗ.
3. **Bar-5 systematic.** Validates **β** at the protein layer — the `a·β` path, a DIFFERENT and much easier claim
   than βᵗ, and MH-83 already shows it works (7/7). This is the "highest-value validation left", and it is intact.
4. **`a`, the propagation slope — VALIDATED (MH-105).** One coefficient per gene ⇒ well-powered; median **0.432**.
   ⭐ **Falsification test PASSED:** obligate-complex subunits (RPL/RPS/PSM/SNRP, n=139) have median `a_g` **0.117** vs
   **0.437** for all other genes (n=9057), **Mann-Whitney p=1.4e-29** — protein-complex **BUFFERING**, recovered
   genome-wide, exactly as predicted ⇒ `a_g` IS the propagation slope.
   ⚠ **RETRACTED:** "the unbiased `a` is up to 38% off the marginal one" — that came from a 9-gene hand-picked panel.
   Genome-wide the miRNA-confound correction is **IMMATERIAL** (paired median bias **−0.001**; adjusted<marginal in
   **52.4%** = coin-flip, 597 genes). It matters only on strongly miRNA-coupled genes (ZEB1 −38%), not globally.

**The honest headline this axis now delivers:** *the protein layer CONFIRMS the mRNA-mediated model and BOUNDS the
translational residual below the per-gene resolution of a 101-patient cohort.* A real, publishable, negative result.

### ⚠ Precision on the claim — the translational component EXISTS; it is not per-gene resolvable
An earlier wording here ("miRNA adds **nothing** to protein beyond the transcript") **over-claimed** and is corrected.
That statement is about *out-of-sample prediction*. The **association is real but sub-threshold**:
- **MH-34's own frame re-run on this panel** (partial Spearman(aggregate pressure, protein | mRNA, C), proper C):
  **13/17 genes NEGATIVE** = repression-directed (binomial **p≈0.025**) — a consistent directional signal —
  but **p<0.05 in 1/17 (KRT8) and BH q<0.10 in 0/17.**
- ⇒ **The partial-correlation frame agrees with the OOF frame.** It does NOT contradict it, and **MH-34 is NOT
  overturned**: MH-34 ran on **1,132 genes and found 11–29 hits ≈ 1–2%**; on a 17-gene panel that predicts **0.2–0.4**
  hits and we observe 0–1. **Quantitatively consistent.**

**Correct claim:** the translational component is **real, weak, and widespread-in-direction**, but **NOT resolvable
per gene at n=101**, and therefore **cannot support a βᵗ latent** (per-family, PIP, Shapley, discovery). Association
(MH-34/35, valid, FDR on 1132 genes) ≠ per-gene resolution ≠ a fittable latent field. **The falsification is of the
MODELING OBJECT, not of the biology or of MH-34.**

### ⚠ NAME COLLISION — two different `a`s (do not conflate across the two doc sets)
| | `gauge.Gauge.a` (state session) | **`a_g` (this doc)** |
|---|---|---|
| model | `β_source ≈ a·β_target + c` | `protein_g ≈ a_g · mRNA_g` |
| meaning | cross-cohort **scale NUISANCE** (composition/platform/C), explicitly "never claimed" | **mRNA→protein propagation slope** |
| resolution | **ONE GLOBAL scalar** (applied out-of-fold over genes, still global) | **PER-GENE** — median **0.397**, IQR 0.23–0.62; 11% of genes have a_g<0.1 |
Its per-gene *variation is the biology* (low `a_g` = a protein-buffered gene ⇒ the CORUM complex-buffering axis BM).
⚠ On the "unbiased `a`": the miRNA-confound correction to `a` is **IMMATERIAL GENOME-WIDE** (paired median bias **−0.001**, adjusted<marginal in **52.4%** = coin-flip, 597 genes) — it is material only on the few STRONGLY miRNA-coupled genes (ZEB1 −38%), NOT a global bias (ZEB1 0.551→0.340, BCL2 −25%, VIM −19% are real but unrepresentative).

### ✅ `gauge`'s `min_sd` scale-dependence — FIXED by the state session (`min_sd_rel`, relative to each cohort's median sd).

---

## 0. The one-paragraph version  *(⚠ SUPERSEDED by §0a — retained for provenance)*

Protein is **not** the lever that beats the mRNA floor — **measured, it carries **≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**) of the mRNA channel's
information about β, and ≤7.6% even if the mRNA→protein slope were perfect.** What it uniquely gives is a
**new latent: βᵗ, translational repression that leaves the transcript intact and is therefore invisible to
RNA-seq at any sample size.** βᵗ is verified real **and survives real composition adjustment** (out-of-fold vs a
30× shuffled null, under `build_C("cptac")`: PTEN z=2.7, PEBP1 z=2.5 — independently re-deriving MH-35's Class-B
survivors by a fully multivariate route) but low-power (n=101). The build is a **cohort-internal mediation model
in CPTAC** (`P ~ X + M + C`), fused to TCGA only through a **hierarchical β with a LEARNED cohort-heterogeneity
variance τ²**. Everything stays **Gaussian-conjugate ⇒ Gibbs. The engine does NOT move to HMC.**

---

## 1. Verified feasibility (every number checked, not asserted)

> **⚠ ALL NUMBERS BELOW RE-DERIVED under `learned.confounders.build_C("cptac")` (2026-07-12).** An earlier
> version of this table used a hand-built "rich C" that **included PAM50 dummies** — exactly what the state
> session's `lineage_verdict` prohibits (PAM50 is computed *from* the mRNA; 27/50 of its classifier genes are our
> targets; it costs −36% of |β|). Those numbers (`a_g` 0.437/0.390, ratio 1.7%, σ² ratio 0.91) are **VOID** and
> were **re-measured, not rescaled** (measured-only gate). The verdict did not change — **it strengthened.**

| quantity | measured (correct C) | consequence |
|---|---|---|
| `a_g` (protein~mRNA slope) | **0.397** (median, 9,196 genes) | the link coefficient |
| σ²_mRNA / σ²_protein | **0.79** (0.643 / 0.819) | |
| **information ratio protein→β** | **≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**) — the formula shown is the superseded `a_g`=0.397 derivation: = 0.397²·(101/1041)·0.79; **≤7.6% even at a_g=1.0** | **protein cannot move β. Robust to every plausible a_g.** Converges with the state session's independently-derived **0.6%** for CPTAC-mRNA→TCGA. |
| βᵗ out-of-fold vs 30× shuffled null | **2/9 z>2 = PTEN (2.7), PEBP1 (2.5)**; median OOF +0.202 | βᵗ **real, survives composition, low-power** |
| βᵗ under composition (the load-bearing re-check) | strong hits **HOLD**; weak ones **die** (ZEB1 0.084→−0.023, PDCD4 0.193→−0.017) | composition removes false-positive discordance, not the real signal |
| **protein OOF > mRNA OOF in CPTAC** | 0.348 (8/9 z>2) vs 0.206 (3/9) | **predicted by the model** (protein sees a·β+βᵗ) |
| CPTAC n | **101** (miRNA∩protein∩clinical) | the binding constraint |
| gene coverage | 9,556 protein / 9,196 resid | **not** a constraint |
| p/n | **0.28** median, **0.63** PTEN (TCGA 0.027) | 10× worse conditioning ⇒ wide βᵗ |
| family coverage (with `.N` fix) | **314/324 complete (96.9%)**, 7 genuinely-absent arms | clean |
| **locus-CN proxy** (nearest-5 genes) | **r = 0.998** vs ASCAT truth; **no gene deserts** (median 28 kb, zero >1 Mb) | CPTAC CN instrument is REAL |
| β transportability TCGA↔CPTAC | median corr **0.175** — **UNDERPOWERED, uninterpretable** | ⇒ **learn τ², don't assume** |
| paired normals | **18** (PDC000120), **PROTEOME-ONLY** (verified: 0 GDC files on their normal samples) | Phase 3 = paired protein-baseline test |

---

## 1b. PREREQUISITES — what is needed from the user

| # | need | why | who |
|---|---|---|---|
| **P1** | **CIBERSORTx fractions for the 122 CPTAC samples (Wu-major panel)** | CPTAC's C today is `purity + CIN` (2 covariates) vs TCGA's 11. Composition is **load-bearing for βᵗ** (mediator-conditioning ⇒ any unmeasured mRNA↔protein confounder biases it). **BLOCKS Phase 0/2.** | **USER / state session.** CIBERSORTx runs via **docker**, which is **not available in this shell**. The state session ran `gtex_wu_major` today (2026-07-12) ⇒ they have a working docker path + the Stanford token (in `gtex_wu_major/run.log`). **I can prepare the CPTAC mixture-TPM file + the exact command; someone with docker executes one command.** |
| **P2** | **Panel decision: Wu-major** (do NOT re-deconvolve TCGA onto HBCA) | `data._DECONV_COLS` is Wu-major and baked into C ⇒ switching invalidates every persisted learned output (downstream-ripple rule). **Effectively already decided** — `gtex_wu_major/` exists as of 2026-07-12. | confirm with state session |
| P3 | MS **technical** load proxy (TIC / total protein / plex norm factor) | the only *non*-expression-derived translational-load term (see §3.3 — PC1-of-proteome is REJECTED) | I'll look in PDC; **may not exist**. Non-blocking — falls back to the plex term. |

### 1b.1 RNA SOURCE — a real fork, decided by role (checked 2026-07-12)
Two CPTAC RNA sources exist and they are **only r=0.81 correlated per-gene across samples** (median; 10% of
expressed genes <0.5), with a **1.8× SD difference**. Neither is a better measurement — the non-circular test
(which RNA agrees better with the *independent* protein assay) is a **wash**: LinkedOmics median corr 0.3713 vs
STAR 0.3678 (Δ=0.008, p=7e-18 but 2% relative; STAR wins on 45% of genes). So the choice is decided by role:

| use | source | why |
|---|---|---|
| **C / CIBERSORTx** | **GDC STAR TPM** (mandatory) | LinkedOmics `.cct` is **gene-median-centered log2 with negative values** ⇒ absolute between-gene scale destroyed ⇒ **unusable** by CIBERSORTx. Built: `scripts/cptac/build_cptac_cibersortx_mixture.py` → `output/cibersortx_transfer/mixture_cptac_brca.txt` (53,929 × 133; **101/101 of the analysis cohort**). |
| **model M / β_CPTAC** | **GDC STAR TPM** (`log2(TPM+1)`) | **GAUGE.** The hierarchical pooling `β_CPTAC ~ N(β_TCGA, τ²)` requires both βs on ONE scale. TCGA's Y is `log2(TPM+1)` from GDC STAR; LinkedOmics has **1.8× the SD** ⇒ τ² would absorb a pure SCALE artifact instead of real cohort heterogeneity, corrupting the quantity that decides how much we borrow from TCGA. |
| **legacy validation lane** (`load_cptac_layers`; MH-33/34/35/36/40/77/83/84) | **LinkedOmics (unchanged)** | marginally better there; within-CPTAC correlations need no TCGA gauge; switching would invalidate published numbers for no gain. **Do NOT switch it.** |

**⚠ MANDATORY ROBUSTNESS CHECK (not optional):** **βᵗ is the coefficient on X _conditional on M_.** Conditioning
on a *different* M (r=0.81!) gives a *different* βᵗ. Since neither source is better, silently picking one is a
hidden researcher degree of freedom **on the headline result**. ⇒ **βᵗ is computed under BOTH RNA sources.** If it
does not survive the swap, that is a **finding and a caveat**, not something to bury.

### 1b.2 CIBERSORTx refsample — use `scref_wu_major.txt`
The two existing Wu-major runs used **different reference subsamples**: `tcga_wu_major/` (→ the model's C) used
**`scref_wu_major.txt`**; `gtex_wu_major/` (state session, 2026-07-12) used **`wu_major_300cells.txt`**. Checked:
**same 9 Wu-major lineages** (NOT a panel mismatch — an earlier framing of mine, corrected), but the derived
signatures differ (3317 vs 3425 genes, Jaccard 0.69; per-lineage profile r **0.92–0.97**) ⇒ similar-but-not-identical
fractions. **CPTAC is pinned to `scref_wu_major.txt`** (matches the model's C). For strict comparability the cheap
fix is for the **state session to re-run GTEx on `scref_wu_major.txt`** — NOT to re-derive TCGA's, which would
ripple through every persisted learned output.

**NOT needed from the user (all local or self-fetchable):** CPTAC CN (local, gene-level, proxy r=0.998) ·
the 18 NAT proteomes (PDC GraphQL, client exists) · `target_cn` (local) · TMT plex (already fetched) ·
segment-level CN — **does not exist, checked three ways** (2026-07-12): (i) **repo/superfolder filesystem search**:
the only `.seg` files are TCGA ASCAT3 (`data/CNV_TCGA/CNV_extracted/`); CPTAC has gene-level only
(`HS_CPTAC_BRCA_2018_CNA.cct`); (ii) **GDC**: CPTAC-2 Breast has **no Copy Number Variation data category at all**;
(iii) **cBioPortal `brca_cptac_2020`**: gene-level GISTIC/log2 only, segment endpoint empty. Recoverable only by
re-calling CN from controlled-access WGS/WXS BAMs — and **unnecessary**: the nearest-5-genes proxy hits **r=0.998**
vs ASCAT truth with **no gene deserts**.

---

## 2. The model space (the frame — memory `protein-mediation-model-space`)

Per gene: **X** = seed-family doses · **M** = the gene's mRNA · **P** = its protein · **C** = confounders.

```
        ┌────── β ──────▶  M ────── a ──────┐
   X ───┤                                    ├──▶  P
        └─────────────── βᵗ ────────────────┘

M = −X·β  + e_M                β  = mRNA destabilization  (what §6b fits)
P = a·M − X·βᵗ + e_P           βᵗ = TRANSLATIONAL repression (mRNA-invisible)
⇒ P = −X·(a·β + βᵗ) + …        a  = mRNA→protein propagation
```

**The three core regressions — not redundant:**

| # | regression | gives | status |
|---|---|---|---|
| 1 | `M ~ X + C` | **β** | ✅ the current §6b model (TCGA, n=1041) |
| 2 | `P ~ X + C` | **τ = a·β + βᵗ** (total) | ✅ but only MARGINAL (`protein_anticorr`, `ood_protein`) |
| **3** | **`P ~ X + M + C`** | **βᵗ** (coef on X) **and `a`** (coef on M) | ❌ **NOT BUILT — the core build** |

**The miRNA is a CONFOUNDER of the mRNA→protein slope** (a common cause of M and P): it pushes M down *and*
P down, so the marginal `P~M` slope absorbs that shared push and is **biased up** on strongly-coupled genes — ZEB1
0.551→0.340 (**−38%**), BCL2 −25%, VIM −19% — **but ⚠ this is IMMATERIAL GENOME-WIDE (MH-105: paired median bias
−0.001; 52.4% = coin-flip over 597 genes). The mechanism is real; the global correction is not.** Old text follows:
BCL2 −25%, VIM −19%. Model 3 fixes it and yields βᵗ and `a` from ONE fit.

Existing `protein_resid` is a 2-step approximation of Model 3. **Checked: it is ≈unbiased for βᵗ** (ratio to
joint = 0.983 — an earlier `(1−R²)` attenuation prediction was WRONG and is retracted). But it does not give
`a`, gives no joint uncertainty, and is marginal/redundancy-blind.

**Extensions (named; 6 is the important one):** 4 bivariate/SUR · **5 mediation decomposition** (proportion
mediated `a·β/(a·β+βᵗ)` per-edge WITH width — upgrades MH-34's descriptive claim) · **6 target-gene CN as an
INSTRUMENT for the mediator M** · 7 `locus_cn_cptac` as instrument for X · 8 state-conditional · 9 saturating
link (**GATED OUT** — failed held-out gate) · 10 CORUM complex-buffering · 11 reverse P→M (unidentifiable).

---

## 3. Architecture — hierarchical β with a LEARNED τ² (supersedes "don't fuse")

βᵗ is estimated from `P ~ X + M`, which conditions on CPTAC's **own** M ⇒ **cohort-transportability-free**.
(Subtracting TCGA's β instead would misattribute any cohort difference `β_TCGA ≠ β_CPTAC` to βᵗ — manufacturing
translational repression out of cohort heterogeneity.)

**But conditioning on M strips X's variance** (R²(M~X)=0.63 for PTEN) ⇒ βᵗ is noisy at n=101. And the
alternative (βᵗ = τ − a·β, borrowing TCGA's n=1041) is *more precise but biased* if β doesn't transport.
**A genuine bias–variance tradeoff — and the transportability test is UNDERPOWERED (median corr 0.175, but
β_CPTAC is itself mostly noise at p/n=0.63, which attenuates the correlation even if β transports perfectly).**

⇒ **Do not decide it by hand. Let the model learn it:**

```
β_TCGA   ←  TCGA mRNA likelihood (n=1041)
β_CPTAC  ~  N(β_TCGA, τ²_cohort)          ← τ²_cohort LEARNED
βᵗ       ←  CPTAC protein equation
```

Nests both extremes and picks automatically: cohorts agree ⇒ τ²→0 ⇒ borrow TCGA's n=1041; cohorts differ ⇒
τ²→∞ ⇒ fall back to CPTAC-only. **Gaussian-conjugate ⇒ stays on Gibbs.** Transportability becomes an
*estimated quantity*, not an assumption.

### 3.1 ENGINE VERDICT (for the state session): **Gibbs. No HMC.**
`CHANNEL_FUSION §7/§L/§J` claims the protein "discordance link" is non-conjugate and forces J2/HMC. **That is
wrong and is corrected here:** discordance is an **additive** term (`P = a·M − X·βᵗ`), hence **linear, hence
Gaussian-conjugate**. Only *saturation* (L2, `a/(a+K)`) breaks conjugacy — and that link already failed its
held-out gate on mRNA. **The state session can build against Gibbs.**

### 3.2 The hazard, and the fix
Conditioning on M is **conditioning on a mediator** ⇒ any unmeasured **U → M and U → P** (composition, purity,
plex, translational capacity, RNA degradation) biases `a` and hence βᵗ. **Fix = instrument the mediator with
the target gene's own CN** (`target_cn` drives M strongly; no plausible path to P except through M). Kills
*all* unmeasured mRNA↔protein confounding. **CPTAC gene-level CN available 122/122, on disk.**

### 3.3 ⚠ The global protein-load factor is REJECTED as a C term (three precedents)
PC1-of-the-proteome is an **expression-derived factor sharing noise with the outcome P** — structurally the
same object rejected three times: **Decision B** (expression-derived purity → over-control, used CPE),
**MH-100** (fuller proliferation control → over_control 2.3:1 over confound genome-wide → global control
REJECTED, per-gene flag), **MH-101** (host lens → mean dC −0.032, net over-control → gene flag).
⇒ Prefer a **technical** load proxy (MS TIC / total protein / plex normalization factor — not proteome-derived).
Falsification test before adopting anything: **does it kill PTEN (z=3.1) and PEBP1 (z=3.6)?** If yes, it's
over-controlling. If an expression-derived term is unavoidable → MH-100 OOF-2×2 per-gene flag, never a global term.

### 3.4 Confounder block C — channel-specific (axis N2)
TCGA's learned C = `CPE, 8× Wu-major deconv, target_cn, mal_prolif` — **no batch** (batch lives in the pressure
spine, not the learned regression). CPTAC C mirrors it: `purity, 8× deconv (P1), target_cn, mal_prolif, CIN`.
**TMT plex enters the PROTEIN equation ONLY** — it is an *assay* artifact (not expression-derived ⇒ no
over-control hazard) and empirically SHARPENS (MH-34: `protein_resid` survivors 11→29). Run with/without as
a sensitivity, as `batch_plex/` already does.

---

## 4. Inherit / validate / add

**INHERITS unchanged (no new machinery):** `spike_slab._gibbs_posterior` (incl. `channels=`, `nu=`,
`return_samples`) · `families.collapse_by_family` (so βᵗ is commensurate with β) · `attribution.bayes_shapley_identity`
(the MH-94 per-draw Shapley is linear in βᵗ ⇒ applies untouched) · `within_family.delta_pooling` · the evidence
prior · `eval/cptac_validation.load_cptac_layers` · `learned/eval/ood_protein.py` (built + hub-validated, MH-83).

**βᵗ is FAMILY-level** (§8 applies verbatim: translational repression acts through the same seed-matched 3′UTR
sites ⇒ same-seed members share it). **Member attribution is INHERITED, not re-estimated:**
`M_member^protein = δ_share_m × βᵗ_f` — dose-delivery is a property of the locus, not the outcome layer.

**VALIDATES:**
- **V1 — δ transportability (the prize).** `locus_cn_cptac` gives CPTAC its OWN exogenous instrument ⇒ compute
  δ_m independently in CPTAC and compare to TCGA's. **A direct, non-circular test of the assumption the whole
  cross-cohort program rests on.** Did not exist before.
- **V2 — the multivariate collapse test.** MH-84's "59% of gold edges discordance-coupled" is **marginal and
  redundancy-blind**. Model 3 is multivariate. **Predicted: a substantial fraction collapses** (cf. MH-83 Q2,
  miR-520/C19MC −0.16→~0). Auditing our own headline.
- **V3 — Bar-5 systematic.** `ood_protein` 7 genes → universe + shuffled-protein null + FDR.

**GENUINELY NEW:**
1. **βᵗ** — a per-family translational-repression latent with width/PIP/Shapley. The mRNA channel gives
   **zero** information about it ⇒ a **new dimension of M**, not a better estimate of the old one.
2. **`locus_cn_cptac`** — validated CN instrument (r=0.998, 482/482 loci once the coordinate hole is patched).
3. **Protein > mRNA coupling in CPTAC** — new; **predicted by the two-latent model before measurement**.
4. **`a_g` VALIDATED as the propagation slope (MH-105, refined with OFFICIAL CORUM 5.3 — MH-106).**
   ⚠ **REFINED:** buffering is a **LARGE-OBLIGATE-COMPLEX** phenomenon, **not generic complex membership**.
   CORUM subunits (n=3519) median `a_g` **0.410** vs non-CORUM (n=5677) **0.443** — real but **tiny** (p=4.5e-4).
   The evidence is the **DOSE-RESPONSE with complex size**: largest-complex **2–3 members → a_g 0.407** vs
   **≥31 members → 0.261**; **Spearman(complex size, a_g) = −0.079, p=2.6e-6**. My earlier RPL/RPS/PSM/SNRP proxy
   (0.117 vs 0.437, p=1.4e-29) over-stated the *general* claim precisely because it IS the large-complex set.
   ⚠ The *"unbiased a is 38% off"* deliverable is **RETRACTED** (immaterial genome-wide: paired bias −0.001, 52.4%).
5. **Bug found:** the `.N` arm-name bug bites the CPTAC path (17 families in a 10-gene probe).

---

## 5. Phases and gates

**Phase 0 — prerequisites, C block, gauge (BLOCKING).**
Matched CPTAC C (needs **P1**) · `target_cn` (free) · plex (protein eq. only) · patch the 137-locus coordinate
hole (`mirna_mature_loci.csv` via `pre_gene_id` → 482/482) · **per-cohort cache key for `_arm_name_map`**
(global `_CACHE` would poison the TCGA map) · CPTAC CN log2-ratio → absolute (`2·2^log2r`) · §E single-edge
gauge **on tcga105** (same patients ⇒ isolates scale from transportability).
**Gate:** gauge equality, else nothing downstream is valid.

**Phase 1 — Bar-5 systematic + the βᵗ scan (the deliverable).**
Scale `ood_protein` to the universe (shuffled null, FDR) · fit Model 3 genome-wide · V2 collapse test.
**Gate:** does βᵗ beat its shuffled null genome-wide? **Predicted yield: tens–low-hundreds of genes, NOT
thousands.** Posteriors will be wide; the width is the honest output.

**Phase 2 — `locus_cn_cptac` + hierarchical τ² + V1.**
Build the CN instrument; add the learned-τ² cross-cohort β; run the δ-transportability test.

**Phase 3 — the 18 paired normals: a PAIRED-BASELINE test (scope settled, not blocked).**
**VERIFIED (GDC, 3 independent queries):** the 18 NAT cases have **ZERO** GDC files on their *normal* sample;
CPTAC-2 Breast has 134 cases and exactly **1** Solid-Tissue-Normal sample registered; the LinkedOmics RNA and
protein matrices are both exactly 122 columns with **identical, tumor-only** sample sets. Krug 2020 profiled
**"122 treatment-naive primary breast cancers"** — the **18 NATs are PROTEOME-ONLY comparators** (PDC000120).
⇒ **No NAT miRNA and no NAT mRNA exist.** So: no model can be *fit* in normal tissue, and **`a_g`-in-normal is
NOT estimable** (an earlier promise of "`a_g` normal vs tumor" is **retracted**).

**What the 18 CAN do — and it is worth having.** Their matched *tumors* are inside the n=101 (full miRNA + mRNA +
protein). So use the NAT purely as a **per-patient protein BASELINE**:
```
ΔP_g = P_tumor,g − P_NAT,g        (18 paired patients — removes each patient's own protein baseline)
test:  do genes with high βᵗ lose protein relative to THEIR OWN normal?
```
A paired design kills inter-patient protein variance, so it is orthogonal to everything else here and it tests
βᵗ against a *normal-tissue reference* rather than a cohort mean. n=18 ⇒ powered only for strong effects
(the PTEN/PEBP1 tier), but it is a real, independent corroboration axis. Fetch: PDC GraphQL (client exists,
`scripts/cptac/fetch_cptac_brca_tmt_plex.py` already talks to PDC000120) — **no user action needed.**

**GATED OUT:** the saturating link (L2 → HMC). Predicted never triggered.

---

## 6. Ceilings — what this will NOT deliver

> **⛔ THE "1.7% / ≤8.8%" NUMBERS IN THIS SECTION WERE VOID AND CONTRADICTED §1 ON THE SAME PAGE (fixed 2026-07-17).**
> They were computed with **PAM50 in the confounder block** — exactly what `lineage_verdict` **prohibits** (§6b.1:
> PAM50 is computed *from* the mRNA; **27/50** of its classifier genes are our targets; it costs **−36% of |β|**).
> §1 and §7.8 had already **re-measured** them (measured-only gate: re-measured, **not rescaled**) under
> `build_C("cptac")`, but **this section was never updated** — so the doc has been carrying two different
> information ratios in two places since 2026-07-12. Corrected below to §1's values. ⭐ **THE VERDICT IS UNCHANGED
> EITHER WAY — at 1.7% or at 4–6%, protein CANNOT move β**; the correction makes the ceiling *tighter*, not looser
> (≤7.6% vs the void ≤8.8%). Void wording retained verbatim at the end of this section for provenance.

- **It will not improve coupling.** **≈4–6%** of the mRNA channel's information about β (⚠ **CORRECTED, MH-108**:
  the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives
  **≈6%**, and the direct Bar-5 retention² check gives **≈4%**), with a **pre-registered ceiling of ≤7.6% at
  `a_g`=1.0**. Anyone claiming protein "breaks past the mRNA floor" is wrong. *(⛔ VOID wording, retained for
  provenance: "**It will not improve coupling.** 1.7% (≤8.8% at any `a_g`)." — PAM50-contaminated C, see banner.)*
- **It will not give per-gene βᵗ for most genes.** Only 2/9 of a *hit-enriched* panel cleared the null at n=101.
- **It will not resolve within-family identity at protein.** βᵗ is family-level; members come from inherited δ.

---

## 6b. INHERITED from the state session (2026-07-12) — binding on this build

1. **PAM50 must NEVER enter C** (`lineage_verdict`). It is computed *from* the mRNA; 27/50 of its classifier genes
   are our targets; it costs **−36% of |β|** on exactly those (ESR1/miR-18 −79%). If lineage adjustment is needed,
   **IHC ER/HER2 is the non-circular instrument — and it is a per-edge FLAG, not a global control.**
   ⚠ *This caught a real error here:* an earlier `a_g` figure in this plan was computed with PAM50 in C. Retracted
   and re-measured (§1).
2. **Posterior widths are ~25% too narrow** — CORRECTED 2026-07-12 (bagged NNLS 1.37×, Gibbs **1.29×**, Student-t
> ν=7 **1.13×**). ⚠ My earlier note to you said "35–70%" and "Student-t does not fix it" — **both RETRACTED**: they were
> measured on a biased gene subset (the same scale-dependent sd-floor bug you found). **MH-92's `nu=7` substantially
> FIXES the calibration (0.77 → 0.89) — turn it on for your βᵗ power calc.** Two things do NOT propagate: **Shapley
> shares are RATIOS so a common β inflation cancels entirely** (re-derived: MH-94's width is 1.8× too *wide*), and all
> OOF/permutation results are untouched. SBC remains the wrong tool.

---

## 7. Retractions made while forming this plan (all mine)

1. "The discordance link forces HMC" — **wrong**, it's additive ⇒ conjugate ⇒ Gibbs.
2. "Protein enters as a `channels=` pseudo-observation" — **wrong** for protein (no instrument; p/n=0.63 makes
   the reduced-form vector unestimable). It is a **second likelihood block**.
3. "No CPTAC CN counterpart exists" — **wrong**; gene-level CN → locus CN at **r=0.998**.
4. "`protein_resid` is attenuated by (1−R²)" — **wrong**; measured ratio 0.983.
5. "Gene deserts limit CN coverage" — **wrong**; no deserts (median 28 kb), it was a coordinate hole in our map.
6. "Don't fuse into TCGA at all" — **superseded** by the hierarchical learned-τ² form (user-driven).
7. "25 families absent from CPTAC" — **wrong**; the `.N` bug. Real answer: 7 arms.
8. **`a_g` = 0.437 / 0.390-under-rich-C, information ratio 1.7%, σ² ratio 0.91** — **VOID**: the "rich C" contained
   **PAM50**, which `lineage_verdict` prohibits. **Re-measured (not rescaled) under `build_C("cptac")`: a_g 0.397,
   σ² ratio 0.79, ratio **≈4–6%**, ≤7.6% at a_g=1.0.** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the
   ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5
   retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**)
   Verdict unchanged and *strengthened*. ⚠ The **1.7% / ≤8.8%** figures in **§6** were part of this same VOID set
   and went uncorrected until 2026-07-17 — see the banner there.
9. **PREDICTION FAILED (logged per the measured-only gate).** I predicted real composition would broadly weaken βᵗ
   (median OOF 0.193 → ~0.12–0.15). **It did not** — median went *up* (0.202), the two strong hits held (PTEN
   z=2.7, PEBP1 z=2.5), and only the *weak* ones collapsed (ZEB1, PDCD4 → ~0). My intuition about
   composition-driven discordance is therefore unreliable in this domain ⇒ **check, don't argue, on anything
   downstream of it.**
