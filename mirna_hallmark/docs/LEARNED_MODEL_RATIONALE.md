# Learned model — RATIONALE (the *why* behind every estimator)

The conceptual companion to `LEARNED_MODEL_METHODS.md`. METHODS is the **spec** — notation, objective,
formula, one line of "what it computes" — deliberately kept clean. This doc holds the **reasoning** that
does not belong in a formula sheet: *why this construction and not the naive one, what breaks otherwise,
what the estimator is really identifying, and how to read its output honestly.* Section numbers mirror
METHODS exactly (§0–§17), so a reader can hold the two side by side.

- **What is the formula?** → `LEARNED_MODEL_METHODS.md`
- **Which estimator for which job, and chosen over what?** → `LEARNED_MODEL_ESTIMATOR_MAP.md`
- **Does it actually work (numbers, positive controls, honest negatives)?** → `LEARNED_MODEL_VALIDATION.md`
- **Why this design at all (decision log)?** → `LEARNED_MODEL_DESIGN_RESPONSE.md`

---

## The one idea that organizes everything

Every estimator here answers a question about the map *miRNA-arm abundance → target mRNA*, and the whole
design turns on **matching the estimator to the estimand's robustness**. Two questions that look identical
have opposite stability under the central nuisance — **seed-family collinearity** (same-seed arms move
together):

- **Aggregate / selection / prediction** questions (does coupling exist? does the model beat abundance?
  does a new edge enter?) read the *bundle* `X·M` or an in/out decision. Collinearity inflates *coefficient*
  variance but **not** fitted-value variance — the aggregate is stable even when the individual weights are
  not.
- **Coefficient-read / attribution** questions (how much does *this* arm contribute? how does its weight
  shift across states?) read the individual `M(m,g)` — precisely the object collinearity destabilizes
  (single-fit bootstrap correlation ≈ 0.03).

> ⚠ **ESTIMATOR ASSIGNMENTS UPDATED 2026-07-17 — the distinction above is durable; the estimators it used to
> name are not.** This block previously routed the aggregate to *"the **model** (adaptive lasso)"* and the
> coefficient read to *"**bagged NNLS** on a fixed family support"*. Both are **pre-convergence**. Under §6b
> the model is **ONE dense learned-τ² non-negative Bayesian posterior per gene, with two readouts of it**
> (π≡1 dense → coupling/attribution/identifiability; evidence-π → discovery), and **the adaptive lasso is
> RETIRED to a baseline** — it is beaten on OOF coupling genome-wide (−0.168 vs −0.152, wins 58%, Wilcoxon
> p=9e-16), and since `nnridge_cv ≈ lasso`, the edge is the **learned τ²**, not L2. So **both** bullets are
> now served by **the same posterior**; the robustness distinction survives as *why the two readouts are read
> differently*, not as two different fits. **Bagged NNLS is not retired** — it remains correct for the
> **cross-cohort gauge**, because Gibbs's heavy-tailed posterior SDs break the errors-in-variables correction
> (it returns a=4.1 where the split-half truth is 1.0). *Rule of thumb: **bagged NNLS for the GAUGE, Gibbs for
> the MODEL**;* on the model itself Gibbs is simply better (split-half reproducibility ρ=0.822 vs 0.729).
> 🔬 **Open, not decided here:** §5 still names unpenalized bagged NNLS as *canonical attribution* and
> `attribution.shapley_identity` still accepts `canonical_M` (bagged NNLS) **or** a Gibbs β draw, while §6b
> assigns attribution to the dense posterior. **Which is canonical for attribution is unresolved — do not
> infer it from this block.** Current state: [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) Axis 1/3.

If you ever feel two sections are "doing the same thing twice," this is why they are not: they read
different functionals of the same fit, with different robustness, so they need different machinery.

---

## The frame behind the frame — M is a latent, every data source is a CHANNEL

*(Absorbed from `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` §0 — the durable conceptual frame only. That
doc's ~60-axis design map is parked and largely unbuilt; its engine claims are superseded — the engine
verdict is **Gibbs**, see `STATE_OF_PLAY.md` Axis 1.)*

The generative reading that organizes what counts as "adding information" to this model: **`M` is a
latent regulatory parameter, and every data source is a noisy *channel* observing it through its own
link.** The shipped model is one slice of that object — it fuses exactly **two** channels, a
curated-evidence *prior* and the tumour-mRNA *likelihood*, and everything else (CN, CPTAC protein,
chimeric/K_D) sits *outside* as validation. The frame says those could instead **feed** the posterior,
with the model's several "jobs" staying **readouts** of one latent rather than separate estimators.
The β latent is **per seed family** — same-seed members share it (§4) — so "which arm" is a
dose-delivery readout, not a coordinate of the latent.

**The honesty clause, which is the whole reason this is a map and not a mandate.** Re-estimating the
*same* observational data with a Bayesian prior earns nothing (verified Δρ = −0.002, pre-convergence).
**Fusion pays only where a channel is non-redundant / exogenous**, or where instrument-sharing adds
identification. Where a channel is redundant, the fused model *reduces to the current one* — no gain,
and it takes on link-misspecification risk for nothing. So every candidate channel is gated against the
frequentist interim plus a shuffled null before it is built, and the burden of proof is on the channel
to show it is carrying information the mRNA likelihood does not already have.

---

## §0 Notation & identification — why residual space, why POOLED-HE

**Why everything lives in the C-residual space `(I − P_C)`.** We never want to know "does a high-miRNA
patient have a low target" — that correlation is polluted by purity, copy number, and proliferation, all of
which move miRNA and mRNA jointly. We want the *within-stratum* question: holding those fixed, does more of
this arm mean less target. Residualizing **both** predictor and target on `C` (Frisch–Waugh–Lovell) makes
every downstream coupling a **partial** association by construction, so no later step has to "remember" to
adjust — the confounders are already gone from the coordinates.

**Why POOLED-HE is the candidate set, not all predicted edges.** The universe of *predicted* miRNA→target
pairs is enormous and mostly false; feeding it to the model would drown real edges in selection noise. We
restrict candidacy to edges with **experimental** support (miRTarBase high-evidence ∪ TarBase-v9
low-throughput functional). This is a *prior on existence*, imposed before any TCGA data is seen — which is
what keeps the model non-circular. The migration from miRTarBase-HE-only to POOLED-HE was checked to be
gate-neutral on anchors (Δρ ≤ 0.003): it adds coverage without moving the headline.

## §1 Confounder block C — why each term, and the composition trap

The subtle confounder is **cell composition**, and it is a trap because it hides in *both* modalities. A
tumour that is 80% epithelial vs 40% epithelial differs in bulk miRNA and bulk mRNA for reasons that have
nothing to do with regulation. So `C` is built to strip composition **without** stripping the biology:

- `CPE` (purity) and `target_cn` (the gene's own copy number) remove the two dosage confounders that would
  otherwise fake a miRNA→mRNA link.
- `mal_prolif` is the delicate one: proliferation drives both, but "proliferation" as measured in bulk is
  partly *just* the epithelial fraction. So we residualize the proliferation metagene on the
  Cancer-Epithelial fraction first — leaving **per-malignant-cell** proliferation, the biological part, and
  discarding the composition part.
- **Why Cancer-Epithelial is excluded from the `deconv` block** (the deliberate asymmetry): the target is
  *expressed in* the epithelial compartment, so conditioning on that compartment's fraction over-controls —
  it regresses away the very signal we want. We strip the *non-malignant* fractions (they are pure
  nuisance) but keep the compartment the biology lives in.

**Why the deconv default flips by task.** Gate/coupling runs with **core C only**: we are asking "is there
coupling at all," and a leaner `C` has more power. The canonical *attribution* runs with **deconv ON**: when
we name a driver we want a **cell-intrinsic** weight, not one an analyst can dismiss as stroma — so we pay
the power cost to buy the stronger claim. Same model, different burden of proof, different `C`.

**Why batch is post-fit only.** An in-lasso ComBat-style batch correction removes biological heterogeneity
along with plate effects (it cannot tell them apart inside the fit). So batch conditions the **significance
test** (`oof_gate`'s ρ computation) but never the **model** — we correct the ruler, not the object.

## §2 The prior w — inverse penalty, not a location (and where spike-and-slab would go)

The evidence prior is the project's most misread object, so state it plainly: **`w_m` is an inverse penalty,
not a coefficient.** A strong-evidence edge is not pushed toward a large weight; it is made **harder to zero
out**. Formally `λ_m = α / w_m`, so the prior sets *ordering and selection* (who survives shrinkage) while
the **data set magnitude** (§3). This is why the source of `w` is empirically non-critical (ledger ≈
evidence ≈ fused within 0.005 on the gate): the prior only has to get the *ranking* roughly right, and any
reasonable evidence fusion does.

**Why fused/deduped magnitude and not binary inclusion.** Binary "is it curated" throws away that some edges
have one qPCR and some have ten reporter assays. `w_m = Σ_class w_class·log1p(#distinct PMIDs)` keeps that
gradient, deduplicated across databases so a paper both DBs cite counts once, and `log1p` so the tenth PMID
matters less than the second. `w_class` grades by assay *directness* (reporter > western > qPCR > CLIP)
because a luciferase construct is closer to "this arm represses this gene" than a binding footprint.

**Where the spike-and-slab fits — and why the lasso is what ships.** The parked Bayesian design
(`LEARNED_MODEL_DESIGN_RESPONSE §D/E`) replaces the lasso's *implicit* L1 selection with an **explicit,
evidence-graded inclusion indicator** `z_{g,m} ~ Bernoulli(π_{g,m})` — a true "does this edge exist" random
variable, separate from magnitude. It is now **built and benchmarked** (`learned/spike_slab.py`; §2a below).
Verdict: it is competitive and recovers the curated drivers, but it does **not** beat the adaptive lasso's
soft-L1 selection on out-of-fold coupling — which is exactly why the lasso ships and the spike-and-slab
stays the documented alternative rather than the production model.

### §2a Spike-and-slab benchmark — the honest head-to-head
The model: `r = Xβ + ε`, `β_m = z_m·θ_m`, `z_m ~ Bernoulli(π_m)` with `π_m` monotone in the evidence prior,
`θ_m ~ N⁺(0,ν²)` (half-normal slab → keeps the `M ≥ 0` repression box). Evidence enters **only** through
inclusion `π`, so the benchmark isolates the design's actual claim. Estimated by componentwise SSVS Gibbs
(`spike_slab._gibbs_ss`); the gate reads the posterior-mean coefficient `M_m = E[z_m·θ_m]` and PIP `E[z_m]`,
scored on the **same folds / C / family predictors** as the lasso.

- **Hub panel:** beats abundance **4/5**, recovers curated drivers at PIP≈1.0 (ZEB1 miR-200/429, ESR1
  miR-18, PTEN miR-17~92), but matches/beats the lasso only **2/5** (mean ρ −0.313 vs lasso −0.351). It
  *wins* where the true support is a tight cooperative family (ZEB1) and *loses* where many weak arms add
  up (PTEN, CDKN1A) — the lasso's graded shrinkage keeps more of that diffuse signal than a hard in/out
  indicator.
- **Genome-wide sample (107 genes):** the gap holds and is decisive — spike-and-slab is FDR-significant on
  **35%** of genes vs the lasso's **45%**, mean ρ **−0.078 vs −0.125**, and matches/beats the lasso on only
  **21%** of genes. (Run: `spike_slab.gate_fdr(limit=120)` →
  `output/learned/gate_fdr_spike_slab.tsv`; lasso ρ on this sample = −0.125, matching the published
  genome-wide −0.123, so the sample is representative.)
- **Takeaway:** an explicit inclusion prior buys **interpretability** (a real PIP per edge) but not
  **coupling** over the lasso. It is the right tool when the question is "does this edge exist, with a
  posterior probability" and the wrong tool when the question is "predict the target as well as possible."

**PIP = Posterior Inclusion Probability** = `E[z_m]`, the fraction of post-burn Gibbs samples in which
edge `m` was included. It is the model's "probability this edge is a real regulator, given data + prior";
PIP = 1.00 means every sample kept it (e.g. ESR1 miR-18, PTEN miR-17~92). This is the object the lasso
*cannot* give — a lasso coefficient is either zero or not, with no calibrated probability attached.

**Where it wins, and by how much.** The spike-and-slab beats the lasso exactly where the true support is a
tight cooperative seed-family and hard inclusion is *correct*: **ZEB1 ρ −0.465 vs lasso −0.438 (~0.03
deeper, ~6%)**, GATA3 −0.204 vs −0.200. It loses on the diffuse genes where many weak arms each contribute
a little (PTEN, CDKN1A) — a Bernoulli in/out discards signal the lasso's graded shrinkage keeps.

### §2b Tuning tried (and why none of it beats the lasso globally)
Two principled modifications were run head-to-head on the hub panel:
1. **Full evidence-graded slab** (`priors.edge_priors` → π from the ledger + a slab *scale* `1+gain·μ·τ`
   widened by scanMiR biochemical affinity μ and evidence depth τ; `spike_slab.fit_gene_ss(use_priors=True)`):
   **a wash** — mean ρ −0.309 vs the uniform-slab −0.313. This is the same lesson as §2/§3 for the lasso:
   *the prior sets ordering, not magnitude* — a fancier magnitude prior does not move coupling. (The scale
   prior is deliberately location-at-0, so it never biases the coefficient up; it only loosens shrinkage.)
2. **Base inclusion rate p0** (the real lever): the diffuse genes want *more* inclusion — PTEN improves
   monotonically (p0 0.3→0.5→0.7 gives ρ −0.273 → −0.311 → **−0.393**, finally *beating* lasso −0.334),
   but CDKN1A peaks at p0=0.5 (−0.193) then falls. So the optimum is **per-gene**, and pushing p0 up
   globally turns the spike-and-slab into an approximate ridge — surrendering the sparse, interpretable PIP
   that was its only advantage over the lasso. There is no single global p0 that beats the lasso across
   genes. The one modification that *could* help — an empirical-Bayes / per-gene adaptive p0 — is essentially
   re-deriving the lasso's continuous shrinkage in Bayesian clothing, which is why it is not pursued.

### §2c Learned shrinkage TOWARD the prior — empirical Bayes (the τ² question)
A distinct and principled idea: instead of using the curated prior as an inverse *penalty* (ordering) and
*guessing* the shrinkage strength α, make the prior a *location* the coefficients shrink toward and let the
**data learn how hard** via a between-effect variance τ² — small τ² where the prior agrees with the data
(shrink hard, trust the prior), large τ² where it disagrees (shrink little, trust the data). Built as
`eb_shrink.py`: per gene `β_m ~ N(c·ŵ_m, τ²)`, ŵ = unit-normalised curated weight, `c` and `τ²` learned by
conjugate Gibbs; `c` = how much the prior's *direction* is used at all. **This is genuinely different from
§10** (which learns τ² but anchors to a *data-estimated* cross-gene mean μ_m, not the curated prior).
Result on the hub panel — the data's own verdict on the prior:

| gene | ρ_eb | ρ_lasso | learned τ² | c | τ²/c² |
|---|---|---|---|---|---|
| PTEN | −0.204 | −0.334 | 0.033 | 0.011 | 273 |
| GATA3 | −0.185 | −0.200 | 2.365 | −0.006 | 56870 |
| ESR1 | −0.473 | −0.536 | 0.175 | **0.270** | **2.4** |
| ZEB1 | **−0.596** | −0.438 | 0.110 | 0.039 | 71 |
| CDKN1A | −0.260 | −0.246 | 0.062 | 0.012 | 431 |

- **The learned τ² CONFIRMS the design empirically.** `τ²/c²` (data dispersion ÷ prior-anchored part) is
  71–57000 in 4/5 genes — the data-estimated between-effect variance *swamps* the prior-anchored component,
  so the model chooses to shrink little toward the prior and let the data set magnitude. "Prior sets
  ordering, not magnitude" is thus **derived from the data** (τ² learned large), not assumed. ESR1 (ratio
  2.4, c=0.27) is the one gene where the prior magnitude genuinely pulls.
- **Coupling is a wash** (mean ρ −0.344 vs lasso −0.351): learning the anchor strength doesn't beat guessing
  α, because coupling is insensitive to shrinkage strength (the flat α-sweep, §3). Learning τ² *justifies*
  the guess rather than improving on it — real methodological value (removes the "α is arbitrary" objection),
  zero headline movement.
- **Side-finding:** where EB *wins* (ZEB1 −0.596 vs −0.438, large) it is **not** the prior (c≈0) — it is that
  a **Gaussian ridge beats L1 on tight cooperative families** (keeps all correlated members instead of
  selecting one). Same direction as the spike-slab ZEB1 win, larger. This lead was chased at full scale — §2d.

### §2d Learned-τ² ridge vs the lasso — GENOME-WIDE (the one thing that moved the needle)
Chasing §2c's ZEB1 lead: is a **learned-τ² non-negative ridge** (dense, shrinkage strength learned by the
between-effect variance ν², repress-box intact) a better *coupling* estimator than the shipped lasso? Built
as `shrinkage_compare.py` (all estimators in ONE leak-free harness — train-only residualisation, adaptive
w-scaling, raw-X aggregate — differing only in penalty + λ-source). Disentangling the two axes the ZEB1
teaser conflated (penalty type vs the M≥0 box):

- **Penalty type is inert; λ-SOURCE is everything.** `nnridge_cv` (L2, λ by CV) ≈ `nnlasso` (−0.1525 vs
  −0.1518 genome-wide) — swapping L1→L2 does nothing. But `ebridge_nn` (L2, λ **learned** by ν²) beats the
  lasso: **mean −0.1683 vs −0.1518, wins 379/655 multi-regulator genes (58%), Wilcoxon p=9×10⁻¹⁶.** So it is
  the **learned τ²**, not L2, and CV was the wrong stand-in (it optimises prediction MSE, which over-shrinks
  for coupling). Robust to the ν² hyperprior (InvGamma(a0,b0) swept 0.01→20: PTEN −0.388↔−0.395, always
  beating lasso −0.219) — the win is NOT a tuned knob (§2c's "n≫p pins ν²"). MacKay type-II ML (`ebridge_unc`)
  is the CV-free version and agrees.
- **Two honest caveats full scale exposes (n=5 hid both):** (1) **tail-driven** — mean Δ −0.0165 but
  **median Δ −0.003**; the typical gene barely moves, the lift is concentrated in genes with MANY regulators
  (CCL5, MET n=24, E2F2, cytokines/MMPs) where lasso's sparse selection discards signal. (2) The
  **cooperative-family/collinearity hypothesis FAILS at scale** — Spearman(collin, Δ)=−0.074, non-monotonic
  across tertiles. The mechanism is **regulator-count, not collinearity**: at n≫p, dense learned-τ² shrinkage
  beats sparse selection *when there are many arms to select among*. The ZEB1 collinearity story was small-n.
- **Verdict:** learned-τ² NN-ridge is a genuinely better *coupling* estimator (significant, biology-intact),
  but its edge is concentrated in many-regulator genes, and it is **dense** (no selection) — so it is a better
  fit for the *prediction/coupling* job, not a drop-in for the lasso's selection/discovery role (match
  estimator to estimand). `output/learned/shrinkage_compare.tsv`.

### §2e Learned-τ² for ATTRIBUTION — EB posterior vs bagged NNLS  (`attribution_eb.py`)
Can the same learned-τ² model also serve the *attribution* estimand (the per-edge weight bagged-NNLS reads,
§5), and does it add calibrated uncertainty? Per gene, the non-negative learned-τ² **posterior** (half-normal
Gibbs, π≡1 → `_gibbs_posterior` mean+sd) vs bagged NNLS, on the SAME z-scored family predictors:
- **Stabilises as well as bagging:** reproducibility EB **0.996** ≥ NNLS 0.983 (single-lasso ≈0.03);
  agreement ρ 0.844.
- **Recovers curated drivers** via the identified set |mean/sd|>2: PTEN miR-21/103/181/182 (4/4), GATA3
  miR-27, ESR1 miR-18 — and **correctly declines** the §9-unidentified ones (ESR1 miR-221/222; ZEB1's
  miR-200 split, z≈1.7). The EB |z|>2 line **reproduces §9's conditioned-partial identifiability verdict for
  free** — one fit gives weight + identifiability, consolidating two estimators.
- **Caveat:** the raw EB *mean* is biased UP for un-informed edges (the half-normal prior has a positive
  mean → zero-signal families drift to it, not to 0 like NNLS). Must read through the |z|>2 filter, or switch
  on the evidence-graded inclusion π (spike-and-slab) to zero them at source. So "EB mean" is not a naked
  drop-in for "NNLS weight"; the deployable form is **NNLS-style point estimate + EB posterior sd**.
- **Stabilisation is the SHRINKAGE, not the prior** (works prior-free); the prior's only distinct attribution
  role — tie-breaking collinear members — is deliberately NOT used (it would undo §9's honest "unidentified").

### §2f Combining τ² with the LASSO — the Bayesian lasso (redundant, and why)
Natural idea: learned-τ² is not ridge-by-definition (§1's L1/L2 point) — a **Laplace** prior gives the
**Bayesian lasso** (Park & Casella), which *should* keep sparse selection AND learn the shrinkage AND give a
posterior. Built (`spike_slab._gibbs_blasso`, non-negative: per-edge τ_j² ~ Exp(λ²/2), global λ² learned by
Gamma) and wired into both harnesses. Result — **it collapses to the ridge on all three jobs:**
- **Coupling:** hub mean −0.363 ≈ ebridge_nn −0.362 (no gain).
- **Attribution:** reproducibility 0.996, agreement-with-NNLS 0.87, driver recovery via |z|>2 — all identical
  to the ridge. Does **not** fix the un-informed-edge mean bias — marginally *worse* (PTEN spurious miR-718
  bl_mean 0.170 > ridge 0.149: the Laplace's heavy tail lets a signal-free edge draw a large local τ_j² →
  larger truncated positive mean). |z|>2 kills it either way (bl_z 0.99).
- **Why:** we read the **posterior MEAN**, which for a Laplace prior is **not sparse** — only the MODE (MAP)
  is. Bayesian model-averaging erases the L1 sparsity, so blasso-mean ≈ ridge-mean on both aggregate and
  per-edge. **Genuine Bayesian selection needs an explicit inclusion indicator (the spike-and-slab, §2a) —
  not a Laplace prior — and that estimator LOSES on coupling** (sparse hurts the aggregate at n≫p).

**The fundamental split this exposes, and its reconciliation:**
- **coupling wants DENSE** (ridge / blasso-mean / learned-τ²) — selection hurts;
- **attribution wants a CLEAN/sparse point estimate** — but the sparse spike-and-slab hurts coupling;
- **the |z|>2 identified-set reconciles them** — read the *dense* learned-τ² model through the identifiability
  filter → clean attribution (drivers recovered, §9 reproduced) *from* the same dense model that wins coupling.

**So the one unified estimator is the learned-τ² NON-NEGATIVE model** (ridge; blasso is equivalent, adds
nothing), with three readouts: **aggregate → coupling (§2d), posterior mean±sd → attribution (§5), |z|>2 →
identifiability (§9)**. Combining τ² with the lasso is *redundant*, not additive. See §2g for the one genuine
tension this leaves.

### §2g Adding the INCLUSION indicator (spike-slab-lasso) — the tension that has no single-fit fix
The last combination: inclusion `z~Bernoulli(π_evidence)` + Laplace slab + learned λ² — genuine selection
AND Bayesian shrinkage (`spike_slab._gibbs_ss_blasso`), tested for BOTH coupling and attribution.
- **Coupling:** the inclusion indicator **destroys the ridge's win** — ssblasso mean −0.312, *worse* than the
  dense ridge −0.362 and ≈ the plain lasso. Selection discards the diffuse many-arm signal the dense ridge
  keeps. **Confirmed: you cannot have selection AND best coupling** (the spike-and-slab's coupling loss, §2a,
  in the unified harness).
- **Attribution/identifiability:** the inclusion **PIP is the BEST signal on the hard co-driver genes** —
  recovers ZEB1 miR-200/429 (2/2) and CDKN1A miR-15/16+miR-17~92 (2/2), which the dense ridge `|z|>2` **misses
  entirely (0/2, 0/2)**: the evidence-π rescues evidence-supported drivers that pure data-driven
  identifiability can't resolve under collinearity. But it *under*-recovers PTEN's many drivers (2/4 vs 4/4)
  and the raw mean stays partly contaminated (un-informed edges leak on their included sweeps).
- **The irreducible tension (now fully mapped):** coupling wants **dense** (no inclusion); identifiability on
  co-driver genes wants **inclusion** (evidence-informed). They pull opposite ways ⇒ **no single fit is optimal
  for all three.** But it is the SAME Bayesian family — flip ONE knob:
  - **π ≡ 1 (dense)** → coupling (§2d winner) + attribution point estimate (|z|>2, §2e);
  - **evidence-π (inclusion)** → identifiability, esp. the hard co-driver genes (PIP>0.5).
  **One framework, two π-settings — not one fit.** That is the honest, defensible consolidation of the arc:
  a single learned-τ² non-negative Bayesian model, run dense for coupling/attribution and inclusion-on for
  co-driver identifiability. Full estimator grid in `output/learned/shrinkage_compare.tsv` +
  `attribution_eb.run()`.

**The full 2×2 (slab shape × selection), coupling on the hub — collapses to ONE axis:**

| | ridge slab (Gaussian) | lasso slab (Laplace) |
|---|---|---|
| **dense** (π≡1) | ebridge_nn −0.362 | blasso −0.363 |
| **inclusion** (spike) | ssridge −0.304 | ssblasso −0.312 |

Slab shape (ridge vs lasso) is **irrelevant** — rows are flat (−0.362≈−0.363; −0.304≈−0.312). Only **dense vs
inclusion** matters, and inclusion *hurts* coupling (both spike variants drop below even the plain lasso
−0.320). So the entire estimator space reduces to a single knob: **dense → coupling/attribution, inclusion →
co-driver identifiability.** Gaussian-vs-Laplace and CV-vs-EB are second-order; selection is the axis.

### §2h Where seed-rarity enters — inclusion not slab (holds); but the gain is NOT validated (term OFF)

The specificity/rarity prior (`seed_rarity.py`) had two candidate loci and a free gain; `rarity_bench.py` +
an edge-coupling audit (2026-07-09, soups {PTEN, CCND1, VEGFA} vs clean {ESR1, CDH1}) resolved the locus but
**overturned** the gain.

- **Why not the slab (holds).** Letting a rare-seed edge get a *looser slab* so it shows up in the **dense**
  coupling/attribution fit does nothing: dense OOF coupling **ΔρOOF = +0.000** every gain/gene, attribution
  magnitudes identical to 3 decimals. The §2b/§2c mechanism at its sharpest — **at n≫p the likelihood dominates
  the slab prior** (doubling a rare edge's slab 1.36→2.96 moves its mean by 3e-5). A scale prior has no purchase
  where there is data ⇒ specificity cannot act on model A. Rarity's only possible lever is the discrete inclusion
  decision in π.
- **Why the inclusion gain FAILS (retracted).** It is tempting that rarity in π lifts the rarest edge's PIP past
  0.5 in every soup — but that PIP-crossing is a **proxy the prior trivially satisfies**, not biological truth.
  Audited against independent signals, the lifted edges are wrong: PTEN miR-105 **unexpressed** (RPM 0.48), CCND1
  miR-99/100 **non-coupling** (ρ −0.02, q_BY 1.0), VEGFA miR-718 non-coupling (ρ −0.07, q_BY 0.45). Across all
  covered edges **Spearman(rarity, edge-ρ) = +0.10/+0.11/−0.01** (wrong sign; on PTEN the commonest seeds couple
  2× harder than the rarest). **Seed rarity ≠ realized coupling.** A global rarity channel selects for obscurity
  (tiny/artifactual targetomes of unexpressed orphans), not specificity. The one real case (miR-18a in miR-17~92)
  is *moderate*-rarity + expressed + inside a coupling family — rarity as a *within-family* tiebreak, not a global
  promoter.
- **State:** `rarity_gain` defaults to **0 (off)**; plumbing retained. Re-enabling requires an **expression gate**
  + **within-cooperative-family scope** + **coupling validation** — a redesign, not a re-tune (SYNTHESIS §6a).

**Consequence for "one fit, N jobs":** the dense/inclusion split is intact; rarity, if it ever ships, lives on the
inclusion side only and is provably absent from coupling/attribution — but it does not ship as-is.

## §3 The model — why non-negative, why adaptive, why α barely matters

**Non-negativity** is not a convenience; it encodes the biology "miRNAs repress." Without the `M ≥ 0` box,
OLS/lasso will happily assign a positive coefficient to an anti-correlated-by-confounding arm and call it
activation — noise wearing a sign. The box turns "explain the target" into "explain the *repression* of the
target," which is the only claim we can defend.

**Why adaptive (Zou reparametrization) and not plain lasso.** Plain lasso shrinks all edges equally, so it
drops a well-curated edge and a random one on the same footing. The adaptive penalty lets curated edges
resist shrinkage — the prior earns its keep exactly here. Implementation is the standard trick: scale
column `m` by `w_m`, run a *uniform*-α non-negative lasso, unscale `M = b·w`. Nothing exotic; the "adaptive"
behavior is entirely in the column scaling.

**Why α is a sparsity knob, not a coupling knob.** Across α ∈ [0.001, 0.02] the mean gate ρ is flat
(−0.340) and only the *count* of nonzero edges moves (9→7); only α=0.05 over-shrinks. So α is not a
hyperparameter we agonize over — it tunes how many edges survive, not whether the aggregate couples. The
canonical *attribution* (§5) sidesteps α entirely by being unpenalized bagged NNLS.

## §4 Seed-family collapse — the estimand is family→gene

This is a statement about **what is identifiable**, not a modeling convenience. Same-seed arms share the
targeting seed and are near-collinear in expression; no estimator can split "was it miR-200b or miR-200c"
from bulk data — the design matrix simply lacks the variation. Rather than report a fabricated split, we
**collapse the family into one predictor** and state the honest estimand: *this seed family → this gene*.
The individual arm is then a **nomination**, not a claim. Pooling in true RPM space (`log2(1+Σ(2^x−1))`)
and not in log space matters because abundances add in linear, not log, units. The family prior is the
**max** member weight (the family is at least as evidenced as its best-studied arm).

### §4a Two kinds of collinearity — why they need different resolutions
*(Absorbed from `COLLINEARITY_AND_IDENTIFIABILITY.md`.)*

The starting fact is the one the whole design turns on: when predictors co-vary, the **aggregate** `X·M`
is stable but the **individual member weight** is not. `ρ(X·M_lasso, X·M_baggedNNLS)` is within 0.03 and
the aggregate bootstrap SD is ~0.01, while a single fit's *per-member* coefficient has bootstrap
correlation **0.03** across resamples. So the gate / discovery / aggregate-coupling read the *stable*
object and are fine; everything that asks **"who owns g"** — attribution, budget, identity, driver
nomination — is a collinearity problem. There are **two biologically distinct kinds**, and conflating
them is the trap.

**Kind A — SEED-FAMILY collinearity (shared *targetome*).** Arms with the same seed (miR-200b/200c/429)
have *identical predicted targets*, so they repress `g` **through the same sites**. This is the hard
kind: confounded not just in abundance but in **targeting** — even with independent abundance you could
not tell which arm's binding cut `g`, because it is the same binding site. So the **estimand changes**:
we do not fight the collinearity, we make the family the unit (§4 collapse), apportion realized pressure
to arms by linear-RPM share as an explicit *abundance nomination* rather than a resolved coefficient
(§5a), and resolve drivers only by conditioned coupling (§9). A truly collinear member has ~no residual
variance ⇒ it lands in `family-only` automatically, no threshold needed. **Honest floor:** when members
are genuinely collinear, the abundance apportionment is the only split available, and it is flagged as a
nomination, not a coefficient.

**Kind B — GENOMIC-CLUSTER collinearity (shared *transcription*, distinct targetomes).** miRNAs
co-transcribed from one locus (miR-183/96/182; miR-17~92) or co-regulated across loci (the two miR-200
clusters) have **different seeds ⇒ different targetomes**, but **correlated abundance**. This is the
*nuisance* kind — the abundance collinearity is real, the biology is arm-specific — and it is
**fundamentally easier than A because they do not share targets**: typically only one cluster-mate has a
seed match in `g`'s 3′UTR, so **sequence disambiguates what abundance cannot**. The co-transcribed
non-targeter is riding along, not regulating `g`. Two sub-cases matter: **B1** (same polycistron — one
transcript, one locus, one CN ⇒ only the target side and conditioned expression separate them) and
**B2** (separate loci, co-expressed ⇒ each locus's copy number independently perturbs *its* cluster, so
CN can act as the separator). Hence the resolution ladder, cheapest first: **sequence-first prune**
(among abundance-collinear pairs `|ρ|≥0.7`, drop the one with ~zero biochemical potential on `g`), then
**conditioned coupling** across cluster-mates, then **CN disambiguation for B2** — which is deliberately
*not* a forced single owner: CN speaks only where its locus instrument is strong, and a same-polycistron
pair is guarded back to `CN-blind`. Genome-wide, Kind B is a minority of the curated attribution set —
only **60/1571 pooled-HE genes (4%)** carry a cross-seed collinear cluster — and the sequence-prune
barely fires there (3 riders); it earns its keep in the **discovery/orphan** set, where site-less
predicted edges appear.

**One line:** Kind A you resolve by *collapsing* (the family **is** the answer); Kind B you resolve from
the *target side* (sequence tells you which co-transcribed mate is real), with CN as the exogenous
separator when the mates sit on different loci.

**Why the "resolve to members" endpoint stays parked.** The single principled method covering both kinds
is **hierarchical δ-pooling** — instead of collapse (which is the infinite-shrinkage limit), fit a
family/cluster weight plus a **shrunk per-member δ with posterior width**, turning "unidentifiable" into
a quantified credible interval. It is **parked**. NB do not confuse it with the *cross-gene* pooling that
IS built (§10, `β(m,g) ~ N(μ_m, τ²)`), which was tested and did not help (Δρ −0.002; effects are
gene-specific) — that is an **orthogonal pooling axis**, and its null does not predict that within-family
δ-pooling would fail.

**The identifiability ceiling is a LIMIT, not a fitting failure.** This is the honest reading of the
model's own uncertainty: median posterior SD/|β| = **0.799**, and only **27.8% of 5,117 units** are
identified at |z|>2 (`STATE_OF_PLAY.md` Axis 1). Same-seed arms share the binding site, so the design
matrix simply lacks the variation to split them — no estimator, and no amount of fitting effort, creates
information the data does not contain. Reporting a split anyway would be fabrication (§9).

## §5 Canonical attribution — why bagged NNLS instead of the lasso's coefficients

The single-lasso coefficient `M(m,g)` is the **unstable functional** (bootstrap corr 0.03 under
collinearity), so reading it directly to attribute a driver would report noise. Two fixes stack: **NNLS**
(unpenalized, so no α-dependent shrinkage distorting the *relative* weights) and **bagging** (average over
40 bootstraps → variance down ~1/B → reproducible at corr ≈ 0.99). The **variance floor** (`sd < 0.1` → zero
weight) enforces a rule the biology demands: an arm undetectable-in-this-state cannot be a driver in this
state, regardless of what the fit wants to do with its tiny residual variance. z-scoring is not cosmetic —
level-scale NNLS is ill-conditioned and recovers high-variance noise arms; z-scoring puts every family on
the same gauge so the weight reflects coupling, not dynamic range.

**Why the same fixed support across states (§5's ΔM).** To compare a weight in tumour vs NAT vs GTEx, the
*support* must be identical — otherwise "the weight changed" is confounded with "a different edge was
selected." Fixing the family support and re-fitting per state makes `dM = M_state1 − M_state2` a comparable
estimand.

### §5a Budget contribution — why family-size correction
A naive "M·X" would let a three-arm family count its abundance three times. Apportioning by each arm's share
of the family's pooled abundance guarantees `Σ_arm contribution = M_fam·X_fam` — the family contributes once,
split among its arms by realized abundance. A singleton is unchanged (its share is 1). This keeps
"realized force" additive and honest at the family level.

## §6 Validation gate — why out-of-fold, and what genome-wide FDR means

**Why OOF and not in-sample fit.** An in-sample `X·M` correlates with `Y` partly because it was *fit* to
`Y` — that is circular and always "passes." Holding out folds and scoring `X·M_train` on unseen patients
makes the gate a genuine *predictive* test: the weights must generalize, not memorize.

**Why three bars, not one.** A significant coupling is not enough — it must beat the two baselines it could
be a trivial restatement of. **Bar-1 (abundance):** if equal-weighting every arm couples just as well, the
*learning* added nothing. **Bar-2 (curated fixed-M):** if freezing the prior weights couples just as well,
the *fitting* added nothing (this is the heuristic the whole program nests and must improve on). Passing
requires `rho_model < rho_abund` **and** `rho_model ≤ rho_curated`.

**What "genome-wide FDR" actually means here.** Each gene yields one pass/fail; testing ~1,100–1,500 genes
at p<0.05 would produce dozens of chance "passes." The genome-wide FDR (`gate_fdr`) computes, per gene, a
one-sided partial-*t* p-value for `rho_model < 0` on `df ≈ n − 8` (the −8 = the confounder columns removed),
then controls the **expected fraction of genes we call coupled that are actually null** across the whole
scan. **Benjamini–Yekutieli, not just BH**, because the gene tests are **dependent through shared
regulators** (one hub arm drives many genes); BY is valid under arbitrary dependence, BH would be
anti-conservative. Caveat baked in: the test is **only meaningful for ≥2-regulator genes** — for a
singleton `X·M ≡ abundance` up to scale, so there is nothing for "learning beats abundance" to detect.

## §7 Per-edge coupling — what it adds over the gate, and why permute

The gate (§6) is an *aggregate* verdict on the whole regulator bundle. §7 zooms to **one (arm → gene) edge**
and asks whether *that* arm's abundance anti-correlates with the target after partialling out `C` and going
to ranks (partial Spearman). This is the edge-level analog of the gene-level gate.

**Why a permutation null instead of a textbook p-value — two independent reasons:**
1. **Family dependence.** Same-seed arms are not independent tests. A parametric p-value assumes an
   independent-sample null and is miscalibrated. The permutation builds the null empirically and — crucially
   — applies the **same permutation π to every arm row**, so the null inherits the *same* cross-arm
   correlation structure as the observed statistic. Apples to apples.
2. **No clean analytic null for a partial-rank statistic.** After residualizing ranks on `C`, the sampling
   distribution of the partial Spearman is not a standard Spearman/t distribution. The permutation gives the
   exact finite-sample null for *this* statistic.

The two reported p-values do different jobs: `p_perm` is the raw permutation tail; `p_z = Φ((obs−μ)/σ)` is a
smooth Gaussian-tail approximation from the null moments, more stable in the extreme tail where few
permutations land. Multiplicity (BH / BY / family-Simes) then runs across edges, with the **seed family as
the hypothesis unit** (§14) so collinear arms are not counted as independent evidence.

## §8 Causality (CN-locus instrument) — from correlation toward cause

§6–§7 establish *association*; they cannot rule out a third cause driving both arm and target. §8 reaches
for causation with an **instrument**: the copy number of the miRNA's own locus, `CN(m)`. The logic of IV —
if a variable that moves the arm *only through the arm* (genetic dosage) also moves the target, that is much
harder to explain by confounding than a plain correlation.

- **First stage / F > 10:** the instrument must actually move the arm; a weak instrument (F ≤ 10) gives
  unreliable IV and is discarded.
- **Reduced form (`rho_reduced`):** does residualized `CN(m)` anti-correlate with residualized `Y`? This is
  the *genetic-dose* channel — expected < 0 (more locus copies → more arm → less target).
- **`rho_obs` (the comparison baseline):** the plain observational partial-Spearman of the arm's **measured
  abundance** vs the target, same `C`. It is the correlational channel. `causal_concordant` demands
  **both** `rho_reduced < 0` **and** `rho_obs < 0` (plus F > 10): the correlation you see is *corroborated*
  by the dosage channel. If `rho_obs < 0` but `rho_reduced ≥ 0`, the observational anti-correlation is not
  backed by genetics → likely confounded, not causal.

**Why cluster-exclusion is mandatory, not optional.** CN moves the *whole genomic cluster*, so a co-located
miRNA that also targets `g` violates the IV exclusion restriction — the instrument is no longer arm-specific.
The fix conditions each candidate on its co-targeters: `arm_unique` requires the arm's own conditioned
partial to survive (< −0.1) while **no co-targeter** does. `strong_causal = concordant AND (cluster_clean OR
arm_unique)` — a causal claim that survives the "was it really this arm or its neighbour" test.

## §9 Driver resolution inside a family — honesty about identifiability

Given a family (§4) that the model selected, §9 asks *which member(s) carry the repression* by measuring
each member's anti-correlation **net of its family-mates** (+C). The design is built to **refuse to force a
winner** when the data cannot support one — the three outcomes *are* the finding:

- **0 survivors → "family-only":** the members are so collinear that conditioning removes each one's residual
  variance (`part → 0`); the split is genuinely unidentified (e.g. ESR1 miR-221/222). Reporting "miR-221
  drives ESR1" here would be fabrication.
- **1 survivor → "single-driver"** (e.g. PTEN miR-181b).
- **≥2 survivors → "co-drivers":** the repression is carried by several members, reported as a *set* (PTEN
  miR-17~92 = miR-106b + miR-20a; ZEB1 = miR-429 + miR-200c).

This is the guardrail that stops family collapse (§4) from silently becoming an over-claim about a single
arm. The empirical ~13:7:11 single:co:family-only split says a large minority of "drivers" are honestly
unresolvable — and the model says so.

## §10 Uncertainty (hierarchical Bayes) — error bars, deliberately not point estimates

The lasso/NNLS give points, not posteriors, and fit each gene alone. §10 adds **calibrated uncertainty** by
fitting a program-wise hierarchical Bayesian regression that pools each miRNA's coefficient across all the
genes it targets (`β(m,g) ~ N(μ_m, τ²)`), yielding a posterior mean ± sd per edge and an "identified edge"
rule `|β/sd| > 2`. Partial pooling shrinks noisy single-gene estimates toward the miRNA's cross-gene
behavior, so an edge reads as identified only if the data support it beyond its own uncertainty.

**Two deliberate scope limits.** (1) Cross-gene magnitude pooling does **not** improve OOF coupling here
(the honest negative in VALIDATION §1) — so §10 is used purely as the **uncertainty object**, not as a
better predictor; the coupling verdict stays per-gene (§6/§7). (2) It answers "how sure," not "who" (§9) or
"is it causal" (§8). Reading it as a driver-selector would overreach — its job is error bars.

## §11 Subtype tests — three different questions, not one

"Does regulation differ by subtype" hides three separable questions, and conflating them is the classic
error: **LEVEL** (does pressure magnitude `X·M` differ across PAM50 — Kruskal), **SLOPE** (does the
*coupling* differ — a subtype×pressure interaction F-test), and **WIRING** (do the *weights* `M` themselves
differ). A subtype can have higher pressure with identical wiring (level differs, wiring conserved) — so you
must test the right one for the claim you want to make.

**Why the pooled wiring estimator corrects the independent one.** Fitting `M` *independently* per subtype
differences two small-n noisy fits and **over-calls remodeling** — it flagged PTEN as rewired. The pooled
estimator (Decision H) makes the whole-cohort `M_all` the **prior mean** and estimates only a *regularized
deviation* `δ_s` per subtype (ridge toward `M_all`), so small-n noise shrinks to zero instead of masquerading
as rewiring. Result: PTEN is borderline (~1.35× across λ, not robust), only RB1-Her2 modestly remodels ⇒
**subtype wiring is mostly conserved**; the robust subtype signals are LEVEL and common-M SLOPE, not
rewiring. This is the design's `M_T = M_H + Δ` nesting in its subtype flavour — the same shrink-the-delta
logic that guards the state comparison.

## §12 Discovery — why a new edge needs three gates, not one

An orphan (m,g) beyond the curated set must clear three independent hurdles, each blocking a different false
positive: (1) **conditioned partial-Spearman** below threshold, conditioning on `C` *and the curated-HE
aggregate* — so a "discovery" cannot be a shadow of the known regulators; (2) **scanMiR biochemical
support** — a real duplex must be plausible, not just a correlation; (3) **composition retention** under
deconv — the edge must survive stripping stromal/immune composition, i.e. be cell-intrinsic. Edge-level FDR
comes from a permutation null-scan (real vs permuted hit counts), and hypotheses collapse to seed-family
→gene (§14). The three gates are AND-ed because each alone has a characteristic failure mode (statistical
shadow / no mechanism / composition artefact).

## §13 Enrichment — a hypergeometric over the right universe

The only subtlety is the **universe**: enrichment of a target set in a Hallmark program is computed over the
**HE universe** `U`, not all genes. Testing against the whole genome would inflate every fold-change (the HE
set is already biased toward studied genes); conditioning on `U` asks the honest question — *given the edges
we could have found, is this program over-represented.* BH across programs.

## §14 Multiple-testing machinery — why the hypothesis unit is the family

The recurring principle: **the seed family is one hypothesis, not N arms.** Counting family-mate arms as
independent tests would both inflate the test count (over-conservative multiplicity) and double-count
correlated evidence (anti-conservative discovery). The toolkit matches the dependence structure: **BH** when
tests are near-independent; **Benjamini–Yekutieli** (the `1/H_n` factor) under arbitrary dependence, e.g.
the gene-level gate where hubs induce dependence; **Simes-per-family** to roll family-mate p-values into one
family q; and a **family-wise min-statistic permutation null** for the permutation tests. Which one applies
is dictated by *where the dependence is*, not by taste.

## §15 Identity vs magnitude — two attribution objects that must never be conflated

Over a gene's active families there are two legitimate but *different* "who matters" answers, and reporting
one as the other is a category error:

- **MAGNITUDE (budget)** = `M_fam·X_fam` normalized — *realized force*, abundance-included. Rewards an arm
  for being abundant. This is *the budget*: who delivers the dose.
- **IDENTITY** = **Shapley** credit for the fixed-weight aggregate's explained variance (LMG for the linear
  model) — *who explains the target*, abundance-removed. Rewards on-target coupling, and is
  **collinearity-fair** (`attribution.shapley_identity`).

The `gap = identity − magnitude` is the interesting object: **+** = a quiet on-target owner (explains more
than its budget — ESR1 miR-18: identity 0.54 / budget 0.20), **−** = a loud passenger (high budget, explains
little — ESR1 miR-22: budget 0.26 / identity 0.01). The gap is *primarily* abundance-vs-coupling and
*secondarily* collinearity (co-varying families split shared R²).

**The two measured reasons they diverge.** (i) **Coupling strength** — an abundant family whose pressure
barely tracks `Y` is a *passenger*. ESR1 miR-22 is the **most abundant** family (abund 16.0) → budget 0.26,
but its pressure barely tracks ESR1 (individual ρ −0.12) → identity 0.01. The converse is the **quiet
owner**: ESR1 miR-18 is scarce (abund 3.0) yet its pressure tracks ESR1 strongly (individual ρ −0.43) → it
owns 0.54 of the coupling on a 0.20 budget. Note the ESR1 families are **not** collinear (abundance
|ρ| ≤ 0.29) — here the gap is *purely* abundance-vs-coupling. (ii) **Collinearity** — co-varying families
split the shared credit; this is the smaller effect, and it appears where families genuinely co-vary (e.g.
ZEB1's miR-200bc/429 vs miR-141/200a). The `indiv_rho` column exposes the mechanism directly: **identity
tracks indiv_rho, budget tracks abundance.**

**The identity↔magnitude rank-correlation is a budget-faithfulness diagnostic** — PTEN **+0.91** means
budget ≈ ownership; ESR1 **+0.45** means budget over-credits loud-but-off-target arms.

**Budget ≠ identity is NOT a weight-calibration failure.** Budget = `M·X̄` = mean *force* (depends on
abundance); identity = Shapley R² = *variance explained* (depends on how the family's variation tracks `Y`).
Abundance ≠ variance-explained, so **even a perfect `M`** keeps a loud-steady/loud-weak family high-budget
and low-identity. `M` is already calibrated to coupling — the fit maximises the aggregate's fit — and the
budget merely re-multiplies by `X̄`. The gap is the *signal* (loud passenger vs quiet owner), not an error,
which is exactly why both are kept.

### §15a ⭐ Why `share` (= β_f / Σβ) is NOT an identity estimator
The single most important thing to get right in this section, and a **live mislabel** in the shipped table
(`learned/readouts.py` labels `share` as *"ATTRIB / IDENTITY — the Bayesian Shapley"* in genome-wide
`readouts_edges.tsv`; see `STATE_OF_PLAY.md` Axis 3). It is not identity — anyone reading
`readouts_edges.share` as identity is reading **β**.

The reason is a property of Shapley, not an implementation slip. **Shapley's whole power lives in the value
function.** For an **ADDITIVE** value function, the Shapley value is trivially `β_f` — each player's marginal
contribution is the same in every coalition, so the ordering-average has nothing to average over. It
therefore **splits nothing under collinearity**, and normalising it just gives `β_f / Σβ`. That is `share`:
a renormalised coefficient wearing a Shapley label.

The doctrine's identity uses the **NON-ADDITIVE** value function `v(S) = R²(X_S · M_S, Y)`. Here a collinear
pair's *joint* R² is far less than the sum of its members' marginal R²s — the coalition is worth less than
its parts because the parts overlap. That sub-additivity is precisely what gives the ordering-average
something to do: it genuinely **divides the overlap** instead of double-counting it. That is
`attribution.shapley_identity` = **LMG on R²**, and it is why identity is collinearity-fair while `share` is
not.

**The discipline that keeps this from being abused:** identity is **⊥ coupling** (§6/§7). A high identity
share is the **distribution of the coupling that already exists**; it is **never** evidence that a driver
exists. **Identity says *who*, magnitude says *how much*, NEITHER says *whether*.** You establish existence
with §6/§7/§8 — a driver call still needs FDR-significant coupling — and only then use identity to apportion
credit among established drivers.

**Why the precursor softmax identity was retired.** The precursor (abundance-removed, softmax-share)
identity was shown to be **⊥ coupling** in the wrong way: it surfaced *silenced specialists*
(miR-137→NCOA3, ρ = −0.02 NS) while the abundant arm was the one that actually coupled (miR-17, ρ = −0.25).
That "identity" was a literature/structure hypothesis, not a driver. The learned Shapley identity is
**coupling-based** — it can only credit families that contribute to the *real* coupling, so it cannot
surface a non-coupling specialist. The structural "who is designed to repress this" question did not
disappear; it moved to an explicit expression-free object (§16) rather than being smuggled through a
softmax.

## §16 Structural identity — the expression-free "who is designed to repress this"

Every object above needs the arm to be *expressed* to score. A specialist that is **silenced** in tumour has
zero abundance, zero budget, zero Shapley — it vanishes from the realized picture even though it is the arm
*designed* to own the gene. §16 recovers it with an **expression-free** score: biochemical `potential`
(scanMiR/TargetScan strength × evidence) × `specificity` (the family's pooled evidence-mass fraction aimed at
`g`) × `confidence` (a depth multiplier — no evidence ⇒ no ownership). Because there is **no biochemical
specialist** (every arm hits ~700 targets roughly uniformly), specificity is a *pure evidence construct* —
we are honest that "who specifically owns this gene" is a curation fact, not a biochemical one. Pooling the
**whole seed family's** evidence (not one representative arm) stops a thinly-studied member from faking
concentration. `lost_specialist` overlays this with baseline-active-but-tumour-silenced to surface the arms
a purely expression-based readout would miss.

## §17 Loss/gain mechanism — why promoter methylation, and why bidirectional

If §16 says an arm is silenced, §17 asks *by what mechanism* — and promoter methylation is the testable one.
The **no-cCRE-hop** detail matters: the gene-centric probe annotation misses miRNA loci (which is why the
older pipeline fell back to cCREs), so §17 maps CpG probes **directly** over the pri-miRNA promoter window.
**Bidirectional** because silencing is only half the story: `hyper` (Δβ up, arm silenced → repression *lost*
→ target up) and `hypo` (Δβ down, arm de-repressed → repression *gained* → target down) are symmetric
mechanistic calls. The TSG-vs-oncomiR validation (Δβ +0.061 vs −0.070; miR-124 +0.193, miR-21 −0.176)
confirms the direction convention points the right way.

---

## Reading guide — which doc answers your question

| You are asking… | Go to |
|---|---|
| What is the exact objective/formula? | `LEARNED_MODEL_METHODS.md` (§ same number) |
| *Why* is it built this way / what breaks otherwise? | **this doc** (§ same number) |
| Which estimator do I use for job X, and vs what? | `LEARNED_MODEL_ESTIMATOR_MAP.md` |
| Does it work — numbers, controls, honest negatives? | `LEARNED_MODEL_VALIDATION.md` |
| Why was the whole thing designed this way? | `LEARNED_MODEL_DESIGN_RESPONSE.md` |

*Every §here mirrors the same § in METHODS. If they ever drift, METHODS (the spec) is authoritative for
formulas and this doc is authoritative for reasoning — reconcile, don't let them diverge.*
