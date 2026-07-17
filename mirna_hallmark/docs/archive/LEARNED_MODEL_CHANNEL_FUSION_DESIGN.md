# Channel-fused posterior — design (the one model that absorbs the causal/exclusion machinery)

> **Goal:** specify the *broader* Bayesian object the converged model is a slice of — a **latent regulatory
> field fused from multiple observation channels** — and the concrete, code-grounded first build: the **CN
> exogenous channel** (Gap-2B Bayes: the between-family exclusion as **family-β covariance** in one posterior; the
> arm-level "which arm" question is dose-delivery in Ring-1 L2 — a *separate* latent that can also go Bayesian via
> δ-pooling, but is NOT a repression-β in this channel).
> **What belongs here:** the generative frame, the math of adding a channel to the current Gibbs, the phased build,
> and the falsification. NOT the *implementations* of the **frequentist interims** this Bayes absorbs — those are
> plural (CN between-family = MH-88 `CN_INSTRUMENT §9`; Shapley **identity** = `attribution.py`; **identifiability**
> partial-corr = `parallel_view`) — the doc *compares against* the relevant one per readout, but each lives in its
> own home. NOT results (→ VALIDATION/registry).
> **Update trigger:** when a channel is built/benched, or the frame changes. Move rows to `built` as they land.
> **Sync-partner:** `CN_INSTRUMENT.md §9` (the CN instrument + the MH-88 CN interim this Bayes-ifies),
> `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md §6b` (the current model = the per-gene two-channel slice this generalizes),
> `LEARNED_MODEL_ESTIMATOR_MAP.md` (add the channel-fusion row when built).

> **Status: CN CHANNEL BUILT (2026-07-10 … 2026-07-12); the rest of the frame is still design.** The first exogenous
> channel — CN — is fully built and run:
> - **between-family covariance** form (`channel_cn.multi_family_terms`, **MH-93**) + **continuous-inclusion PIP** (**MH-94**);
> - **scope-3 exclusion gate** (`instrument.exclusion`/`run_exclusion`, **MH-90**) — 45,187 seg-rows; each strong segment
>   now carries `seg_dose_share`;
> - **JOINT coding+HOST pleiotropy** down-weight (`coding_pleiotropy`, **MH-99** + **MH-101**): host folded in guard-exempt
>   (retired MAX); **PER-LOCUS host** (`locus_host_map`→`locus_context.tsv`, **MH-101 per-locus** — each CN locus conditioned
>   on its OWN host) + the **FULL-universe coding scan** (2026-07-12, `do_scan=coding` — **3059 pure-coding + 36 host**
>   down-weights) + `grade_host_downweights` (the ambiguous `repression_direction` set resolved via the OOF lens);
> - **within-family δ-pooling** (`within_family.delta_pooling`, **MH-98**) — dose-delivery, orthogonal to β.
>
> Hot path is 6.4× / coding-scan 17× faster (pure caching+vectorization, outputs verified identical — memory
> `profile-before-batch-loops`). **There is NO unfinished "fold δ-pooling into the channel" remainder** (I earlier
> mischaracterized it as one): δ-pooling is a DELIBERATELY separate, orthogonal latent (F2, §6-F — member-level *dose-delivery*
> δ_m; same-seed members share β_f per §8, so there is no per-member repression β to add as a β-term), AND the fully-joint
> Gibbs form was **tested and found WORSE in MH-98** (mRNA collinearity at n≈1000 buries the CN-resolved member). Folding
> member-level latents into ONE joint posterior only becomes meaningful at the **A3 de-collapsed (arm×gene) field** — i.e. it
> is subsumed by the A-axis continuation below, not a standalone CN build.
>
> **FUNDAMENTAL continuations (§6):**
> 1. **The A-axis / Bayes field (fit UNIT, §6.A) — ✅ RESOLVED (2026-07-12): keep A1 (per-gene); the field does NOT pay.**
>    Probe #1 (`exclusion.tsv`): A1 leaves a large weak-instrument territory on the table (512 weak segments hit ≥3 genes).
>    Probe #2 (`hierarchical.oof`, the A4 μ_m-pooling diagnostic): a WASH even on weak genes (mean Δρ −0.0006; weak-gene
>    5/10 = coin-flip) — μ_m is a poor shrinkage target (β_g heterogeneous). Since A3 (HMC field) is warranted only if A4
>    pays, the cheap diagnostic gates OUT the expensive field build. "One posterior over the whole field" doesn't beat
>    per-gene at n≈1000. So there is **no field to build** — A1 is confirmed as the fit unit. (§6.A ★ has the numbers.)
> 2. **The exogenous PROTEIN channel (§6.I1, CPTAC) — now the SOLE remaining fundamental lever.** With the A-axis resolved
>    to A1 and every internal refinement (CN channel, δ-pooling, coding/host gate, isomiR refit, t-likelihood, cross-gene
>    pooling) landing "immaterial-at-n≈1000", a NON-redundant second likelihood (protein) is the only move left that can
>    reach past the tumour-mRNA floor — validate against CPTAC protein-vs-mRNA discordance (§I1).
>
> **NB — isomiR refit is NOT the A-axis.** It is **axis AS2** (§10 AS): a 5′-isomiR shifts the seed → a *different targetome*
> → it re-DEFINES the edge (same locus, different β), orthogonal to A's fit-unit question. It is necessary (17 seed-heterogeneous
> families, `within-family-structure`) but gates a **headline revalidation** because refitting `X_fam` changes coupling/CN/
> Shapley/discovery — see `ISOMIR_AWARE_MODELING.md`. Design-first per user; grounded in the *current* model (§6b), not `DESIGN_RESPONSE`.

---

## 0. Thesis

The converged model (`DISCOVERY_SYNTHESIS §6b`) is **one dense learned-τ² non-negative Bayesian posterior over the
edge weights M, per gene**, read four ways (coupling / attribution / identifiability / discovery). Read precisely,
it has **two implicit limits**:

1. **It fuses only TWO channels** — a curated-evidence *prior* + the tumour-mRNA *likelihood*. Everything else is
   *external*: the CN instrument is external 2SLS validation, CPTAC protein is external OOD (Bar-5), chimeric/K_D is
   identity-only.
2. **It is fit PER GENE.** But a locus CN event is **one exogenous shock** perturbing many arms *and* testable
   against **every gene those arms hit** — so causal identification is intrinsically **cross-gene**.

**The broader object relaxes both:** *M is a latent causal repression parameter; every data source is a noisy
**channel** observing it through its own link, fused in ONE posterior over the (**family** × gene [× state]) field.*
(The β latent is **per seed family** — same-seed members share it, §8; "which arm" is a dose-delivery readout, not a
latent coordinate — see §1.) The
four "jobs" stay readouts; the new move is that **more sources FEED the posterior** instead of validating it from
outside. Every current tool is then a marginal/conditional readout of this object; **Gap-2B's "CN-anchor on the
cross-family covariance" is the trivial consequence of the CN channel and the mRNA channel sharing the latent.**

**Honesty / scope (the A/F lesson, correctly scoped):** re-estimating the *same* observational data with a
Bayesian prior earns nothing (verified Δρ=−0.002, `WHATS_NEXT §DESIGN-vs-AS-BUILT`, pre-convergence). Fusion pays
**only** where a channel is **non-redundant / exogenous**, or where instrument-sharing adds identification — the
collinear-flat directions and weak-per-gene-instrument regimes. Where redundant it reduces to the current model
(no gain, and link-misspecification risk). **So this frame is a MAP, not a mandate to build the whole field.**
Each channel/unit-widening is gated against the frequentist interim + a shuffled null, exactly as before.

---

## 1. The channel is a Gaussian term on the same β (why this is a MINIMAL extension)

The current dense sampler `spike_slab._gibbs_posterior` fits, per gene, on C-residualised z-scored predictors that
are **already collapsed to seed families** (`FAM.collapse_by_family`) — so its columns are **families** and the
coefficient it samples is **per-family** `β_f = z_f·θ_f`, `θ_f ~ N⁺(0,ν²)`. **There is no per-arm β:** same-seed
members repress identically (§8), so they pool into the family dose `X_fam,f` and share `β_f`; "which arm" is a
dose-delivery *readout* (§2.2), not a coefficient. The per-coordinate update solves a Gaussian conditional with
precision `A = X_fam,f·X_fam,f/σ²` and cross-term `B = (X_fam,f·r/σ²)` (lines 149–151).

**A channel = one more Gaussian observation of `β_f`, so it just adds `(A_c, B_c)` to that same conditional.** For
the CN instrument the algebra is exact and already ours — all per **family**, not per arm:

- **first stage** (gene-independent): `X_fam,f = γ_{f,s}·Z_s + Cδ` — `γ_{f,s}` = how segment-*s* copy number moves
  family *f*'s **pooled** dose `X_fam,f`. Estimated ONCE per (family, segment) from the arm/CN data
  (`family_multi_iv`, the `a` in `exclusion`) and **reused across every gene** — where cross-gene strength pools.
- **structural** = the latent: `Y = Σ_f β_f·X_fam,f` (`β_f` = M, the thing we're sampling; sum over **families**).
- ⇒ **reduced form** `Y = π_s·Z_s` with **`π_s = Σ_{f on s} γ_{f,s}·β_f`** — the segment's reduced-form coefficient
  is a *linear combination of the β's of the **families** on that segment that hit g* (our family AND co-located
  other-seed families — that sum IS the between-family confound, first-class). *Within* a family the arms have
  already pooled into `X_fam,f` with one `β_f`, so the sum is over families, never arms.

So the instrument supplies **`π̂_s ~ N(Σ_{f on s} γ_{f,s}·β_f, s²_π,s)`**, `s²_π,s` = the reduced-form variance
(large when weak ⇒ the channel self-downweights = the graded Bayesian F>10 gate). Its contribution to **family f's**
conditional, holding the other families' β:

```
e_s   = π̂_s − Σ_{f'≠f} γ_{f',s} β_{f'}      # segment reduced form minus co-located FAMILIES' explained part
A_cn  = γ_{f,s}² / s²_π,s                     # add to A
B_cn  = γ_{f,s} · e_s / s²_π,s                # add to B
```

Drop `A_cn, B_cn` into lines 149–151 and the rest of the sampler is unchanged. **Everything Gap-2B needs falls
out of this one term:**

- **Between-family exclusion** = automatic. Co-located other-seed families **and co-amplified coding regulators of Y**
  are *other terms in the same segment sum* `π_s = Σ_{f'} γ_{f',s} β_{f'} + Σ_c γ_{c,s} β_c` (coding `c` via
  OmniPath/DepMap, H2/BP); the shared `π̂_s` is explained jointly, and the posterior distributes it by each term's γ
  + its own mRNA evidence. The induced **posterior COVARIANCE across the co-located β's** (they trade off to match the
  shared reduced form) is *literally* "the CN reduced form as a causal anchor on the cross-family posterior
  covariance." No binary exonerate, no NaN — graded, with width.
- **Within-family: `β_f` (repression) is family-level, but the DOSE-DELIVERY split can be Bayesian too.** The CN
  channel constrains only the family `β_f` (via `X_fam,f`); same-seed members share `β_f` (§8), so there is no
  per-member *repression* β. But "which arm delivered the dose" (Ring-1 L2) has a Bayesian form — **δ-pooling**
  `M_member = M_fam + δ_m`, a shrunk per-member δ **with posterior width** (wide/collapsing for the ~8% `|ρ|≥0.95`
  co-transcribed pairs = the honest resolution limit). This does **NOT** contradict §8 (it grades *dose-delivery
  resolution*, not repression identity — no exoneration, every same-seed member is a real regulator). It is the
  Bayesian upgrade of the frequentist L2 point (chimeric/abundance, co-primary), **orthogonal to `β_f`**, and it is
  the parked "full family pooling" (Decision F / hierarchical δ-pooling). So the CN channel is family-`β_f`; δ-pooling
  is a *separate* L2 latent that can share the posterior.
- **Weak instrument** = `s²_π,s` large ⇒ `A_cn,B_cn → 0` ⇒ the edge reverts to the mRNA-only posterior. Reduces to
  §6b exactly when CN is uninformative — the "no loss where redundant" guarantee.

**Why a Gaussian term (and when it isn't the right one).** `π̂_s` is a **regression-coefficient estimate** (the
mediated reduced form `pi_causal`), so at n~1000 it is **asymptotically normal** (CLT) — a Gaussian `N(π̂_s; Σγβ, s²_π)`
is the correct asymptotic likelihood for it, *and* it is conjugate with the Gaussian Gibbs (the minimal extension).
The honest upgrade is a **Student-t**: `s²_π` is *estimated*, and a t (heavier tails) both accounts for that and
guards against collider-biased outlier edges — worth it if the `δ_s²`-inflated `s²_π` is itself uncertain. Non-Gaussian
is needed for **other** channels, not this one: count data (binding → Poisson/NB), the protein-discordance link, and
ceRNA/occupancy **saturation** all break Gaussianity/conjugacy → HMC (axes **AB / L / J**). So: Gaussian for the CN
channel (justified + conjugate), t as the robust option, non-Gaussian where the *link* or *support* demands it.

---

## 2. The CN channel (Gap-2B Bayes) — the shared mechanism (unit/options in §6)

> This section gives the *mechanism* that every §6 option shares; it uses the per-gene framing (§6 **A1**) only for
> concreteness — the unit (A1 per-gene / A2 neighborhood / A3 field / A4 hier-prior), the π̂ form, γ, gauge, and
> hierarchy are all **open axes mapped in §6, not decided here.**

The channel term is identical whatever the unit; inputs all exist in `instrument.py`:

| quantity | source |
|---|---|
| segment membership (independent-CN units on g's regulators) | `_genomic_clusters` / focal-locus + `family_causal_attribution`'s segment grouping |
| `γ_{f,s}` (first stage CN_s→**X_fam,f**), gene-independent | `exclusion`'s `gamma_s` (`a`) / `family_multi_iv` (per **family**, per segment) |
| `π̂_s = pi_causal`, `s²_π,s` (pleiotropy-corrected reduced form) | **`exclusion`** (`pi_causal`, `s2_pi`) — scope-3, built |
| co-located other-seed **families** on the segment | `pleiotropy()` → `family_of` (as in `between_family_exclusion`) |

Wire: a `channel_cn(gene)` builder returns per-segment `(families, γ_{f,s}, π̂_s=pi_causal, s²_π)` **from `exclusion`**;
`_gibbs_posterior` gains an optional `channels=[...]` arg that augments the per-**family** `(A,B)` per the §1 formulas.
Readouts: posterior mean/sd of each **family's** `β_f`; the **between-family share = posterior mass on our family's β /
segment total**, replacing MH-88's Shapley point estimate with a distribution. (Arm-level "which arm" is *not* here —
it is the L2 dose-delivery readout, §2.2.)

**Coding-gene pleiotropy is a first-class term in the segment sum, not just miRNA families.** A CN event moves
*everything* on the co-amplified segment — including **protein-coding genes** — so if a co-amplified coding gene
regulates Y directly, its path `γ_{c,s}·β_c` is in the reduced form right alongside the miRNA families'
`Σ_f γ_{f,s}·β_f`. This is **not a minor correction**: MH-89 found it is often the *dominant* confound — for the
canonical miR-16→CCND1 edge, **RB1 at 13q14** (a cell-cycle regulator co-deleted with the miR-15/16-1 locus) absorbs
~23% of the CN→CCND1 reduced form, a bigger piece than the miRNA between-family confound. We do **not** estimate `β_c`
as a channel latent — we **condition the CN channel on the coding gene's expression `X_c`** (observed RNA-seq; the
co-amplified regulators of Y enumerated via **OmniPath/DepMap**), which blocks the `γ_{c,s}·β_c` path exactly like the
frequentist H2/BP escalation (condition on RB1). It is already *inside* the T1 direct-effect `δ_s` (which captures all
non-family pleiotropy); the coding-gene step *attributes* that `δ_s` and confirms it. So between-family exclusion in
this channel means **against co-located other-seed families AND co-amplified coding regulators** — both are terms in
the same covariance object.

**On cross-gene pooling (RETRACTING the earlier "joint neighborhood" overstatement).** I'd claimed fitting the
genes a segment hits *jointly* is a distinct win — it is not. The instrument's strength pools across genes only two
ways, both already accounted for: (a) the **gene-independent first stage `γ_{f,s}`** (estimated once from arm/CN
data, reused across every gene — no joint fit needed); (b) **`μ_f` β-pooling** (a family's typical effect across its
targets — the existing `hierarchical.py` μ, axis AF — *not* locus-specific). Genes that merely share a locus have
**no special joint structure beyond these**: `π_{s,g}` is per-gene (needs `Y_g`), each `β_{f,g}` is its own. So
"neighborhood joint fit" is subsumed by shared-γ + μ_f, not a separate broader win.

### 2.1 Feasibility-probe RESULTS (2026-07-10 — E-gauge + T1 exclusion screen, scratchpad)

Ran the axis-E gauge check + a per-edge exclusion screen before any sampler code. **Three findings that shape the build:**

1. **GAUGE = PASS (on a valid instrument).** On the exclusion-clean canonical edge **miR-92a-3p → PTEN** (F=150),
   the IV-implied standardized β (**+0.319** = π/γ) ≈ the observational β (**+0.273**), gap **0.045**. So the CN
   channel term `N(π̂; γ·β)` and the mRNA likelihood **land on the same β** — the fusion premise is sound, the units
   are coherent (axis **E1** confirmed).
2. **The Wald π/γ is unstable → use the precision-weighted term (never π/γ).** Across one family's loci β_IV=π/γ
   ran +1.6 → −12 (small-γ blow-up). The channel term `A_cn=γ²/s², B_cn=γ·e/s²` **divides by nothing** (small γ →
   weak pull) — confirmed the right form empirically (axis **B1**).
3. **EXCLUSION is the binding constraint, and it fails for the MAJORITY.** A per-edge **T1 test** — does the CN→target
   reduced form SURVIVE conditioning on the arm's own abundance (the mediator)? — screened all 403 testable well-
   instrumented edges: **only 108 (27%) are exclusion-clean** (`pleiotropy_frac = π_cond/π_raw < 0.4`); **median frac
   = 0.75**; F>10 overstates the causal-usable set ~4×. The **canonical miR-16→CCND1 FAILS** (frac 0.91): its locus is
   the chr13q14 RB1/DLEU cell-cycle CN hotspot (T1 0.242→0.220; RB1↔CCND1 ρ=+0.57 absorbs ~23%; MH-88 exonerates the
   miR-15/16 family entirely on CCND1). **miR-92a→PTEN passes** (frac 0.16) because 13q31 amplification hits PTEN
   *through* the miRNA. **Validity is locus-specific, exactly as theory predicts.**

**Design consequences (folded into the axes):**
- **New REQUIRED admission gate for the CN channel:** `F>10 **AND** T1-clean` (mediator-conditioning), not F>10 alone.
  T1 is cheap and per-edge — it defines the channel's honest support (~108 edges). *This is a first-class element, not
  a refinement.* (Sits under axis **C/gating** as the channel's inclusion rule + a new instrument-validity screen.)
- **Coding-gene pleiotropy (H/BP) PROMOTED from peripheral to prerequisite** — RB1 at 13q14 is a live confounder on a
  canonical edge; MH-88's miRNA-only between-family exclusion is **insufficient**, coding-gene conditioning is required.
- **Gauge/MVP anchor = miR-92a→PTEN** (clean+strong+canonical), *not* miR-16→CCND1 (a lesson: "canonical + strong F"
  ≠ valid instrument).
- **Honest headline:** the CN instrument's exclusion restriction fails for ~54–73% of its own well-instrumented edges
  → the channel is real but its *support is small and must be earned per-edge* (T1). This is the arc's key empirical fact.

### 2.2 The exclusion gate — precise spec (for scope-3)

**Two sides of the IV — do not conflate them** (my earlier drafts did):
- **Instrument side (which CN is exogenous):** per-segment / per-locus CN — focal locus (`_arm_focal_locus`), the segment =
  what a CN event moves (`_genomic_clusters` / the |ρ|>0.999 independent-CN grouping); a multi-locus arm → **multi-IV** (§5,
  its loci jointly, each net of the others). **All segment / per-locus structure lives here.**
- **Exposure/mediator side (what it instruments):** the **whole-family dose `X_fam = log2 Σ_m 2^{x_m}`** (all members, all
  loci) — because **coupling is family-level** (`Y = β_family·X_fam + …`, §8). This is exactly `family_multi_iv`'s exposure.
  The mediator is **never segment-scoped**; off-segment same-family arms are *constituents of `X_fam`*, not controls to
  partial (the whole "on-segment pool + partial off-segment" muddle is retired — it summed to `X_fam` anyway).

**Two complementary exclusion tests on that object:**
- **Over-ID (Hansen-J) — PRIMARY where a family has ≥2 independent-CN segments.** Fit one 2SLS β, test `E[Z_s(Y−βX_fam)]=0`
  across segments. **Collider-free** (never conditions on the mediator) and — because it uses the GMM moment, not per-instrument
  `π_s/γ_s` — **immune to the Wald-ratio small-γ blow-up**. Already built in `family_multi_iv`. Catches **heterogeneous**
  pleiotropy (segments disagree). *Blind spot:* if all the family's segments share the *same* confounder (a co-amplification
  hotspot), they agree → over-ID passes falsely.
- **T1 (mediator conditioning) — runs on EVERY edge (not a fallback).** It is the *only* test that catches **homogeneous**
  pleiotropy (all segments confounded the same way — e.g. proliferation-driven aneuploidy correlating CN genome-wide with a
  cell-cycle target), which over-ID structurally cannot see; and it is the sole test when a family has one instrument. Condition `Y~Z_s` on `X_fam`(+C);
  the residual `Z_s` coefficient is the **direct pleiotropic effect `δ_s`**. Catches **homogeneous** pleiotropy over-ID misses.
  **Read `δ_s` against its SE (is the direct effect ≠ 0?), NOT the ratio `frac = π_cond/π_raw`** — `frac` inherits the *same*
  small-denominator instability as the Wald ratio (small `π_raw`, when `βγ_s` and `δ_s` cancel → the "frac 4.99" cases, which
  are ratio artifacts, **not** a separate collider category). Keep `frac` only as a readable normalisation ("fraction of the
  reduced form that is direct") *when `π_raw` is well-determined*.
  - **T1 is CONSERVATIVE.** Conditioning on `X_fam` (a descendant of `Z_s`) opens a collider path if a shared program
    `U → X_fam, U → Y` is not fully in `C`, inflating `δ_s` → biasing toward "pleiotropic". Safe for an admission gate (clean
    verdicts stay trustworthy) but it means **the "27% clean" is a LOWER BOUND** (some flags are residual-collider, not real
    pleiotropy). `X_fam` (pooled) also *reduces* mediator measurement-error vs a single arm; target-CN pleiotropy is
    pre-controlled (`target-CN ∈ C`).

**T1 always; over-ID additionally wherever ≥2 independent instruments exist — run BOTH, they are complementary not ordered**
(over-ID = collider-free heterogeneous; T1 = homogeneous, the blind spot). A family passes only if **neither** flags it:
admission = `F>10 ∧ T1 δ_s≈0 ∧ over-ID-clean (where ≥2 indep instruments)`.

**Escalation — attribute a flagged `δ_s ≠ 0`.** There is **no "rescue" step**: with `X_fam` the mediator we are at the
family unit *by construction* (the earlier "within-family sibling rescue" was an artifact of the retired on-segment-pool
mediator). So the ladder is pure attribution of the residual direct effect:
1. **Between-family other-seed co-targeters.** Condition additionally on the co-located other-seed families' `X_fam`;
   if `δ_s` collapses ⇒ their (different-mechanism) repression carried it → **down-weight ours**. **STATUS (this session):**
   frequentist interim = MH-88; the **Bayes forms are now built** — the cross-family β-COVARIANCE (§2.4, MH-93, weak) and
   the continuous-**inclusion PIP** (§2.5, MH-94, less-weak: existence is coarser than magnitude). So "between-family" is no
   longer just the MH-88 point Shapley — it's a graded posterior (covariance + PIP-entry).
2. **Co-located coding regulators (H2/BP; OmniPath/DepMap; RB1@13q14).** Condition on the coding gene's expression; if `δ_s`
   collapses ⇒ horizontal coding pleiotropy → **down-weight**. **STATUS: BUILDING (2026-07-10, MH-99).** Mechanism validated
   on the anchor — conditioning miR-16→CCND1's `δ_s` on **RB1@13q14** (co-deleted, CCND1-co-varying +0.57) absorbs **35%** of
   it (> MH-89's ~23%). The general gate auto-enumerates co-located coding genes (GENCODE window) whose expression confounds Y.
3. **Residual after 1–2 = unknown pleiotropy** → exclude, or lean on over-ID with an added instrument.

**Arm-level attribution ("which arm within the family") is ORTHOGONAL to exclusion** — it is the §8 question, answered by
**per-segment CN** (same-seed members on *different* segments) or **Ring-1 L2** (chimeric/abundance, for same-segment members
CN cannot separate), and is **never** mixed into the exclusion test. **STATUS: the Bayes form is now built** — **δ-pooling**
(MH-98, `within_family.delta_pooling`): fuses mRNA + per-segment CN + chimeric into a per-member share with honest width (the
graded upgrade of the discrete `resolution_report`). Grades *dose-delivery*, not repression (§8) — so still orthogonal to exclusion.

In the channel-fusion posterior the gate is not a preprocessing filter but the **admission weight** on the CN term (→ `s²_π`):
pass → full `s²_π`; down-weighted → inflated `s²_π` (weaker pull); unknown/excluded → drop the CN term (mRNA-only). (Ties axis
**C** channel-noise + the `F>10 ∧ over-ID-clean ∧ T1-δ_s≈0` admission rule.)

**AS-BUILT (scope-3, `instrument.exclusion` / `run_exclusion`, 2026-07-10) — the gate is CONTINUOUS, not a filter (verified on
miR-92a→PTEN clean vs miR-16→CCND1 pleiotropic; reproduces the scratchpad `π_raw` 0.241/0.116).** Four things the build+verification
forced, each replacing a cruder first cut:**
- **Standardise (C-residualise + z-score) everything** — raw slopes are on incomparable CN↔log-TPM scales; z-scoring makes `δ_s` a
  standardised partial coefficient and the decomposition `pi_causal = π_raw − δ_s` (= `a·b`, verified to 1e-4) hold.
- **The channel gets `(pi_causal, s²_π, γ_s)`, per segment** — `pi_causal = a·b` (`a=γ_s` = Z_s→X_fam, `b` = X_fam→Y|Z_s) is the
  pleiotropy-CORRECTED mediated reduced form. A "pleiotropic" edge is **not dropped** — its corrected signal (miR-16→CCND1: −0.024)
  still gives `β = pi_causal/γ_s`, just down-weighted.
- **PLEIOTROPY handling is CONTINUOUS (no binary filter)** — this was the over-conservatism: the binary `seg_admit` at a hand-set
  `δ_thr` discarded valid-but-modest instruments; now the correction rides in `pi_causal`/`s²_π` continuously. **BUT the weak-instrument
  cut is NOT optional** (I'd wrongly said "no F>10 needed"): the admission **weight = `A_cn = γ_s²/s²_π`** (the channel's β-precision,
  §1) **F>10-gated**, NOT `1/s²_π`. The shuffled-CN null showed `1/s²_π` admits a **noise floor** — under the null `γ²` is matched by
  `se_a²` inside `s²_π` and they cancel (dead instrument keeps weight ~`1/b²`). The `F>10` weak-IV cut is the null-separator
  (**null-validated: shuffled-CN → total weight collapses 128×, strong-count 86→2**). So: *continuous within the strong set
  (pleiotropy), F>10 gate at the noise floor (weak instrument)* — two different jobs, not one threshold.
- **`s²_π = Sobel (b²·SE_a² + a²·SE_b²) + δ_s²`.** The **`δ_s²` term is essential and parameter-free**: the pleiotropy correction is
  precisely estimated but COLLIDER-biased, so without it a pleiotropic edge gets *higher* precision than a clean one (verified
  backwards); `δ_s²` is the honest bias-variance floor (if all of `δ_s` were collider artifact, the error in `pi_causal` is `δ_s`).
  Result: clean miR-92a→PTEN precision ≈ **689** vs pleiotropic miR-16→CCND1 ≈ **21** (~33× down-weight, continuous).
- **Over-ID (Hansen-J) is a family CONTEXT flag only** — over-powered at n~1000 (rejects trivial heterogeneity, p=1e-4) → heterogeneity
  heads-up, not a per-edge verdict. Per-segment `(pi_causal, s²_π)` is what the channel consumes.

### 2.3 MVP RESULT — the channel is real but WEAK; the CN's value is DIAGNOSTIC (MH-91, 2026-07-10)

Wired `exclusion`'s `(pi_causal, s²_π, γ_s)` into `spike_slab._gibbs_posterior(channels=)` (the §1 `A_cn/B_cn` term;
`channels=None` = bit-identical to §6b) and benched fused-vs-mRNA-only β on 6 cluster-genes (234 channel-families):
- **Directionally correct:** the channel moves β **10× more on flat/collinear families** (mean |Δβ| 0.0030) than
  already-identified ones (0.0003) — it acts exactly where the mRNA can't identify, as designed.
- **But quantitatively weak:** even on flat directions |Δβ| ≈ 0.003 vs an mRNA-sd of 0.017 (**~0.2σ**), with **~zero
  posterior tightening**; material on only a handful of edges (max |Δβ| = 0.049, TGFBR2).
- **Why (fundamental, not a bad function).** The IV weight is `A_cn = γ²/s²_π ≈ γ²·n/σ² = γ² × A_mRNA`, and `γ` (first
  stage CN→abundance) is small (~0.2) ⇒ only ~`γ²≈4%` of the exposure variance is CN-driven ⇒ the IV carries only
  ~`γ²·n ≈ 40` *effective* observations against OLS's 995 — **~25× weaker by construction.** The `F=150` vs "weak"
  paradox: **F tests `γ≠0` (significance); `γ²` is the *variance explained* (strength)** — a locus can be
  significantly-but-weakly instrumenting. `β·γ` (the reduced form) is intrinsically small (a product of two moderate
  effects) yet carries the full `1/√n` noise ⇒ ~2σ. **`β·γ` IS the right object** — it equals 2SLS (same precision),
  so it extracts the *full* exogenous information; there is no higher-SNR exogenous handle (γ alone has no target, the
  observational β is confounded). The weakness is the small exogenous slice, not the function.
- **Verdict:** the channel-fusion is the **exogenous-data completion of the §6b "priors don't move coupling"** result
  — real, non-redundant, but too weak at n=995 to materially move the fused β. So the **CN's value is DIAGNOSTIC**: the
  exclusion GATE (which edges are causally-clean) + the `pi_causal` pleiotropy correction, **not** a broad β-mover.
  Fixed-n ceiling: `s²_π ∝ 1/n` ⇒ a larger cohort (METABRIC-full, EGA-gated) would strengthen it; multi-IV pools
  segments but k is small (marginal). This concludes the channel-fusion **build** arc.

### 2.4 CROSS-FAMILY posterior-covariance form — BUILT + verdict (MH-93, 2026-07-10)

The MVP (§2.3) wired `exclusion`'s **per-family** `(pi_causal, s²_π, γ_s)` one-family-per-term → it never induced the
cross-family covariance. **Built `learned/channel_cn.py`** (`multi_family_terms(gene, cols, members)`): groups all of
g's families' member arms into **genomic clusters** (`_genomic_clusters`, single-linkage `window`) — each cluster = one
shared CN segment `Z_c` — and emits ONE channel term per cluster with `members = [(col_f, γ_{f,c})…]`, `pihat = Z_c→r`,
`s²_π = se(π̂_c)² + δ_c²`. A ≥2-member term = the cross-family case; the `_gibbs_posterior` `members` path (`e_s = π̂ −
Σ_{other} γβ`) then induces the **cross-family posterior covariance** = between-family exclusion with honest width. Added
`return_samples=` to read that covariance from the joint β draws. Unifies the MVP (1-member) + cross-family (≥2-member).

- **GAUGE FIX (axis E, make-or-break) — `π̂` on the RAW `r` scale, not z-scored.** The dense `_prep` leaves `r =
  −resid(Y|C)` **un-z-scored** (sd≈0.15) and z-scores only X; so the channel must compute `π̂ = Z→r` on the **raw r**
  (z-scoring `Z` and each family dose to match `Xz`). A first cut z-scored `r` too → `π̂` ~6× off, silently killing the
  channel. **Verified:** on clean single-family repressors the IV-β (`π̂/γ`) ≈ the standardized single-edge obs-β **both
  on the `r=−resid(Y)` convention** (positive = repression), gaps 0.02–0.08 (miR-103 0.150/0.226, miR-181 0.169/0.208);
  sign-diffs (miR-205 IV +0.15 / obs −0.10) are pleiotropic edges the gate flags, not sign bugs. **The r-scale is a pure
  rescaling: standardized |Δβ|/σ is gauge-INVARIANT** (raw 0.0144 vs z-scored-r 0.0140, corr 0.99) ⇒ the §2.3 weak
  verdict is scale-independent; the fix only makes `channel_cn` consistent with production `_prep`.
- **Falsification ladder:** ① nests §6b (channel off = bit-identical) ✅; gauge (E1) ✅; ④ shuffled-CN null → term
  collapses to mRNA-only ✅. **③ resolves the collinear cluster** = **directionally yes, materially no**: on the canonical
  clusters (miR-17~92/106b-25 → PTEN/CCND1) the channel *does* add cross-family anti-correlation in the right direction
  (PTEN miR-106b↔miR-17 −0.068→−0.069; CCND1 −0.040→−0.043) and shifts the between-family share (~0.49→0.48), but only by
  **~0.7% channel weight (~0.1σ)**.
- **VERDICT (MH-93):** the cross-family Bayes form **inherits the §2.3/MH-91 weakness**. The between-family split is
  **dominated by the mRNA likelihood** (which already induces the collinear anti-correlation in the dense fit); the CN
  anchor adds only its ~0.7% exogenous slice. So the design's "the Bayes covariance *replaces* MH-88's point Shapley" is
  **half right**: the Bayes form supplies the honest *width*, but the CN does **not** materially improve the *resolution*
  — MH-88's frequentist Shapley (conditioning on co-located families' **abundance** = full mRNA-scale info) is the
  stronger between-family resolver; the CN-anchor covariance is its exogenous-but-weak complement. Same "real-but-
  diagnostic" shape as everything CN at n≈1000. **⚠ SCOPE:** this closes only the **magnitude-share (covariance) half**.
  The `_gibbs_posterior` ran **DENSE (π≡1, every family forced in)**, so it grades the *share/magnitude* among all-in
  families but does **NOT** touch the *entry* — MH-88's binary "does the family enter the coalition" gate. The continuous
  replacement for that gate is the **inclusion PIP** (axis G2, §2.5), which is a *distinct and less-weak* result. Remaining
  Gap-2B Bayes items: **G2 continuous-inclusion** (§2.5, done MH-94) + **within-family δ-pooling** (§6-F, the general
  same-seed-arm resolver). Code: `learned/channel_cn.py`, `spike_slab._gibbs_posterior(return_samples=)`.

### 2.5 CONTINUOUS-INCLUSION between-family exclusion (G2) — the binary-entry → PIP replacement (MH-94, 2026-07-10)

MH-93's covariance form grades the **magnitude** but leaves the **entry** binary (dense π≡1). The Bayesian replacement for
MH-88's binary "does family f enter the exclusion coalition" is the **inclusion indicator** `z_f ~ Bernoulli(π_f)` → its
posterior **PIP** = continuous P(f is a genuine co-targeter). Run in **INCLUSION mode** (evidence-π) the CN channel already
lifts the PIP through `log_L` (`B_cn` enters the marginal likelihood → the `expit(logit_π + log_L)` inclusion draw) — so
reading the PIP from an inclusion-mode fit **with** the CN channel IS G2, with **no prior/likelihood double-count** (axis M:
CN stays in θ, flows to inclusion via the likelihood, not a separate prior-π term).

- **RESULT — G2 is where the CN is LESS weak (existence is coarser than magnitude).** On the AMBIGUOUS co-located families
  (0.15<PIP<0.85, where between-family exclusion actually bites), the CN moves the **entry** by **ΔPIP 0.03–0.09** — an
  order of magnitude more than it moved the magnitude (~0.1σ / 0.7%). Canonical case: **CCND1 miR-17~92 PIP 0.485→0.575
  (+0.091)** — the CN existence-evidence lifts the confound family's entry from a coin-flip toward "in," directionally
  right; TGFBR2 reshuffles miR-18↓/miR-19↑. Modest but real, and materially bigger than the magnitude effect — because
  inclusion asks the coarser "is the edge real?" that a weak instrument can still answer.
- **GENERAL (not CN-specific) — Bayesian Shapley over the inclusion posterior.** The binary-entry→PIP + point-Shapley→
  posterior-Shapley reframe applies to **every** Shapley in the program (identity_vs_magnitude, between-family, driver),
  not just CN: (1) continuous entry = `z_f` (PIP-weighted; a usually-out family gets ~0); (2) attribution per posterior
  draw → a *distribution* of shares (vs MH-88's ±0.01 MC point) — for the linear reduced form `π=Σγ_fβ_f`, Shapley is exact
  and cheap (`φ_f=γ_fβ_f` → per-draw share `γ_fβ_f/Σγβ`); for the R²-variance case (identity) it's the same principle with
  per-draw LMG; (3) any channel (CN) is just an informer of the posterior. Code: `channel_cn.between_family_bayes`.
- **PROPAGATED + VALIDATED on the identity Shapley (MH-95, `attribution.bayes_shapley_identity`).** Applied the reframe to
  `identity_vs_magnitude`: its binary `_exonerate_between` rider gate → continuous **PIP** entry, its point Shapley on fixed M
  → **per-posterior-draw Shapley ± sd**. The honest width **validates AND discriminates**: ESR1 miR-18 ownership robust
  (0.598±0.022, PIP 1.0, confirms the "quiet owner" finding) but PTEN miR-18 **sd 0.118 > mean 0.117** at PIP 0.495 = genuinely
  can't tell it owns any PTEN (same arm, opposite reliability — only width shows it); and continuous PIP is *more discriminating*
  than the binary gate (ESR1 kept all 7 vs PIP separating 3 solid / 4 marginal). = the "driver call carries its credible
  interval" upgrade (axis AQ). Confirms the reframe is **channel-agnostic and general** — CN between-family was one instance.

---

## 3. Falsification (must clear all)

1. **Nests §6b:** with the CN channel off (or `s²_π=∞`) the posterior = the current dense posterior, bit-for-bit.
2. **Reproduces the well-posed frequentist interim:** where MH-88 is well-posed (single family, strong F), the
   channel-fused posterior mean β must agree with the 2SLS/portion point estimate.
3. **Resolves what the interim CANNOT:** on the collinear co-located clusters (17~92 → PTEN/CCND1), MH-88 gives a
   binary exonerate / a Monte-Carlo Shapley (±0.01 noise) / NaN partials; the posterior must give a **stable graded
   ownership with honest width** — the payoff that justifies the Bayes over the interim.
4. **Shuffled-CN null:** permuting CN labels collapses `π̂_s → 0` ⇒ the channel exerts no pull (β returns to mRNA-
   only). If a "gain" survives the shuffle it was leakage, not the instrument.
5. **Not the A/F null:** the gain must live on the **flat/collinear directions** an observational prior can't move
   (else it is the −0.002 re-estimation with extra steps).
6. **Admits only exclusion-clean edges (§2.2):** the CN term must be down-weighted/dropped on edges failing over-ID
   (heterogeneous pleiotropy) or T1 (`δ_s ≠ 0`, homogeneous) — a channel that pulls β on the ~73% pleiotropic edges is
   injecting confound, not causal signal (the MH-89 finding). Cross-check over-ID and T1 agree where both are powered.

---

## 4. Roadmap — one frame, channels added one at a time (each gated on §3-style wins)

- **CN channel** (this doc, §2) — Gap-2B Bayes; between-family exclusion as **family-β covariance** in one posterior (arm-level = L2 dose-delivery, not here).
- **Protein channel** (CPTAC) — the abundance→protein link as a second *likelihood* (the translational residual
  mRNA can't see), not external Bar-5. Pays on discordant genes.
- **Cross-gene pooling** — via the gene-independent `γ_{f,s}` (already) + `μ_f` β-pooling (axis AF); *not* a locus-specific "neighborhood" fit (retracted, §2).
- **State channel** — healthy GTEx fit as the prior mean for tumour (`M_T = M_H + Δ`), i.e. the healthy channel;
  Decision H's principled home (subtype flavour already applied).
- **Binding channel** (chimeric/K_D) — already a prior (`priors.μ`, identity); the frame just names it a channel.

**What stays OUT unless a result demands it:** occupancy/ceRNA saturation substrate (failed the held-out gate),
and any channel whose link is not identifiable from data (mis-specification > information).

---

## 5. Code touchpoints (for when it leaves `parked`)

`spike_slab._gibbs_posterior` (+`channels=`, +`return_samples=` args) · `instrument.{family_multi_iv (γ),
family_causal_attribution (π̂ + add SE), _genomic_clusters, pleiotropy, locus_cn, _arm_focal_locus}` · **`learned/channel_cn.py`
BUILT (MH-93, `multi_family_terms`)** · bench vs `between_family_exclusion.tsv` (MH-88) + a shuffled-CN null. Register the
estimator in `ESTIMATOR_MAP` on build.

**✅ SCOPE-3 DONE (2026-07-10, MH-90): `instrument.exclusion` / `run_exclusion` → `output/learned/instrument/exclusion.tsv`**
(45k segment-rows; anchors verified; universe 25% clean/75% pleiotropic). Built CONTINUOUS per §2.2 + the AS-BUILT refinements
in §2.1 (`pi_causal, s²_π=Sobel+δ_s², γ_s`; over-ID context; no filter). The spec below is the as-built contract. (a) **mediator/exposure = the whole-family
`X_fam`** (per `family_multi_iv`, coupling is family-level) — NOT segment-scoped; segment/per-locus structure is instrument-side
only (focal locus + `_genomic_clusters`; multi-IV §5 for multi-locus arms). (b) **Two complementary exclusion tests:** over-ID
(Hansen-J, already in `family_multi_iv`) **primary where ≥2 independent-CN segments** (collider-free, Wald-stable, catches
heterogeneous pleiotropy) **additionally where ≥2 indep instruments**, + **T1 on EVERY edge** (the only test for homogeneous
pleiotropy — over-ID's blind spot — and the sole test at 1 instrument); run BOTH, complementary not ordered; **read
the absolute `δ_s = π_cond` against its SE, NOT the ratio `frac`** (frac has the Wald small-`π_raw` instability — the "4.99"
cases are ratio artifacts, not colliders); T1 is CONSERVATIVE (collider over-flags ⇒ clean is a lower bound). (c) **Escalation**
= attribute a flagged `δ_s` (no rescue step): between-family (MH-88 frequentist; **now also Bayes covariance MH-93 +
inclusion-PIP MH-94**, §2.4/§2.5) → **coding-gene H2/BP (✅ DONE MH-99: `coding_pleiotropy`, RB1→CCND1 absorbs 42%; full-universe scan run 2026-07-12 → 3059 pure-coding down-weights — §H)** → residual/over-ID;
emit an admission **weight** (→ `s²_π`), not binary. (d) **Arm attribution is ORTHOGONAL** (per-segment CN or Ring-1 L2, §8;
**Bayes form = δ-pooling MH-98**), not part of exclusion. Admission = `F>10 ∧ over-ID-clean ∧ T1-δ_s≈0`. §2.1 = scratchpad probes.

---

## 6. Design-space map — options kept in play (NOTHING pruned)

Convention (from `DESIGN_RESPONSE`): every axis lists its options with **buys / costs / the DIAGNOSTIC that would
later separate them**. We are *mapping*, not deciding — each option stays live until a named cheap test kills it.
Axes are roughly orthogonal (a build picks one per axis); dependencies noted.

### A. Fit UNIT / support — "per-gene vs joint neighborhood" (the Limit-2 axis)

> **★ THIS IS THE FUNDAMENTAL CONTINUATION (2026-07-12).** Every exogenous piece the CN channel needs is now built
> (between-family covariance, δ-pooling, coding+host gate, full scan); what remains at the *frame* level is the unit —
> **"one posterior over the whole field" is still A1 (per-gene) in practice.** The next fundamental move is not more CN
> plumbing but resolving A: **does widening the fit unit recover identification A1 leaves on the table?** Resolve it by
> DIAGNOSTIC, not assertion — the concrete first test is cheap and falsifiable (no field build needed to run it):
> 1. **A1-leaves-on-the-table probe — ✅ RUN (2026-07-12, from `exclusion.tsv`): the trigger FIRES.** Of 845 (family,segment)
>    instruments, **796 (94%) hit ≥3 genes**, and **512 (64% of multi-gene segments) are WEAK per-gene (median F<10) yet hit
>    ≥3 genes** — 27,402 segment-gene tests in the weak-but-multiply-testable regime (miR-3613-3p hits 392 genes at medF 6;
>    477 weak segments hit ≥5 genes). So by the stated criterion A1 IS leaving cross-gene identification on the table — the
>    A question is LIVE, not the "expected no" the A2 bullet guessed. **⚠ NECESSARY not sufficient:** the first stage `γ_{f,s}`
>    is already gene-independent/shared in A1, so this large territory only *pays* if joint fitting sharpens the per-gene
>    **β_g** — which is exactly probe #2. Don't build A3 off probe #1 alone.
> 2. **A4 vs A1 (μ_m magnitude-pooling) — ✅ RUN (2026-07-12, `hierarchical.oof`): WASH → A4 does NOT pay.** On a broad
>    representative sample (40 multi-regulator genes, full n): mean Δρ(pooled−solo) **−0.0006**, `pool_helps` 19/40 (coin-flip);
>    **on the WEAK genes where pooling should help most (|ρ_solo|<0.15): mean Δρ −0.0002, 5/10 helped** = also a coin-flip. On
>    the strong canonical panel (PTEN/CDKN1A/BCL2…) pooling actively HURTS (−0.034, and WORSE at small n: −0.18 at n=40) — the
>    OPPOSITE of the strength-borrow hypothesis, i.e. μ_m is a POOR shrinkage target because the per-gene β_g are genuinely
>    heterogeneous. So the weak-instrument territory (probe #1) does NOT convert to coupling gain via cross-gene pooling.
> **★ A-AXIS VERDICT (2026-07-12): A1 (per-gene) is CONFIRMED; the field (A2/A3/A4) does not pay for coupling.** Probe #1's
> territory is real but probe #2 shows the cheap way to exploit it is a wash — and since A3 (HMC field) is warranted ONLY IF
> A4 pays, **the cheap diagnostic gates OUT the expensive A3 build.** "One posterior over the whole field" does not beat
> per-gene at n≈1000 (β_g heterogeneity + the program's recurring immaterial-at-n≈1000 floor). This RESOLVES the last frame-
> level DOF: keep A1. (A2 joint-reduced-form pooling is even less likely — γ is already shared; left un-run as dominated.)
> ~~Only if A2/A4 pay does the A3 field get built.~~
> **isomiR is a DIFFERENT axis (AS2, §10):** it changes the *edge/seed*, not the fit unit — run it as its own headline-gated
> refit (`ISOMIR_AWARE_MODELING.md`), and note it INTERACTS with A (a refit `X_fam` re-instruments every unit).

- **A1 per-gene, shared γ (MVP).** One posterior per gene; the CN channel added; γ_m precomputed once, reused.
  *Buys:* smallest change to §6b; cross-gene pooling already happens through γ's precision. *Costs:* the
  reduced-form direction π̂_s is estimated from one gene → weak-per-gene instruments stay weak. *Diagnostic:* count
  genes where the segment F<10 but the segment hits ≥3 genes — if many, A1 is leaving identification on the table.
- **A2 joint shared-locus NEIGHBORHOOD.** Fit all genes a segment hits jointly. *Buys:* **likely nothing beyond A1+AF**
  (retracted — see §2): the instrument's cross-gene strength already pools through the **gene-independent `γ_{f,s}`**
  (A1 reuses it) and through **`μ_f` β-pooling** (axis AF, not locus-specific); each gene's `β_{f,g}` and `π_{s,g}`
  are its own. *Diag:* does a joint fit sharpen β sd on weak-instrument genes vs A1-with-shared-γ-and-μ_f? Expected
  **no** — kept only to be *tested*, not asserted as a win.
- **A3 full (arm×gene) field.** One posterior over the whole matrix + all segments + μ_m miRNA pooling. *Buys:*
  maximal strength-borrow; every readout coherent. *Costs:* heaviest inference; mis-spec propagates globally. *Diag:*
  only worth it if A2 shows neighborhood-pooling helps AND μ_m cross-gene pooling helps *together* (else A2 suffices).
- **A4 hierarchical-prior compromise.** Keep per-gene fits but link β_{m,·} across a miRNA's genes via the existing
  `hierarchical.py` μ_m level, and let the CN channel enter each. *Buys:* pooling without a joint sampler. *Costs:*
  the CN cross-gene identification isn't jointly used (only the magnitude prior is shared). *Diag:* A4 vs A2 on the
  weak-instrument genes — does magnitude-pooling alone recover what joint-reduced-form pooling does?

### B. The instrument observation — what π̂ actually is
- **B1 plug-in reduced-form coefficient + SE** (the MH-88 quantity, but **OLS-linear not Spearman** — the β·γ
  identity is linear). *Buys:* trivial to produce from `family_causal_attribution`. *Costs:* plug-in ignores joint
  (γ,β,C) uncertainty. *Diag:* compare to B2 on a few segments; if the interval barely moves, B1 is fine.
- **B2 full 2SLS / GMM moment** `E[Z'(Y−Xβ)|C]=0` as the channel (proper IV, C handled jointly). *Buys:* rigorous,
  no plug-in scale error. *Costs:* non-conjugate-ish; more code. *Diag:* does B2 change any between-family share vs B1?
- **B3 per-segment vs per-locus (multi-IV) granularity.** A multi-locus arm = several channel terms (one per active
  source, `multi_iv`); over-ID becomes automatic posterior *tension* (two independent-CN directions that disagree →
  wide/flagged β, replacing the Hansen-J reject→recurse with a graded fit). *Diag:* do multi-locus arms' posteriors
  widen exactly where MH-88's over-ID rejected?
- **B4 conditioned vs raw reduced form.** MH-88 conditions π on the other-seed pool (mediation adj.); in the Bayes
  the covariance handles co-location, so raw π̂ is the natural input. *Keep both* to check they agree (they should).

### C. Channel noise s²_π — how hard the instrument pulls
- **C1 plug-in SE** (homoskedastic) — simplest. **C2 robust/Hansen SE** (matches `multi_iv`). **C3 hyper-prior on
  s²_π** (InvGamma, data sets channel weight). *Buys (C3):* honest auto-weighting; *Costs:* can down-weight the
  channel to nothing or over-trust it. **C4 soft (all segments, SE-weighted) vs C5 hard F>10 gate** (matches MH-88).
  *Diag:* the soft/hard and plug-in/hyper choices should only matter on borderline-F segments — sweep F∈{5,10,20} and
  see where the between-family share flips.
- **✅ C4-vs-C5 RESOLVED (2026-07-12, `instrument.soft_fweight_sweep` + `exclusion(shuffle_cn=)`).** The soft null-corrected
  weight `(γ²−c·se_a²)₊/s²_π = (γ²/s²_π)·(1−c/F)₊` (since `se_a²=γ²/F`) swept over c∈{5,10,20} on the real table vs a
  full-universe SHUFFLED-CN null (which collapses strong segments 17,053→141, ~0.8% false-strong-rate). **C4 at c=10 is a
  strictly cleaner gate than the C5 hard cliff** — near-identical real ranking (ρ=0.995 with `beta_weight`) but **4.4× less
  null leakage** (14,827→3,383 weight; SNR 43→124), by down-weighting the ~11% of real weight in the borderline F∈[10,20]
  band the cliff over-credits. c=5 too permissive (SNR 22), c=20 too strict (drops 40% of real signal). **ADOPTED as the
  DEFAULT (2026-07-12, `_ADMIT_SOFT_C=10`): `beta_weight` = `(γ²/s²_π)·(1−10/F)₊`** — same admitted set as F>10, re-weighted;
  the choice is immaterial at the channel's current ~0.7% weight but it is the strictly-cleaner gate for when the CN channel
  is up-weighted (e.g. the protein-channel fusion). `output/learned/instrument/soft_fweight_sweep.tsv`.

### D. First-stage γ (CN→arm) treatment
- **D1 plug-in point γ** (precompute once; γ is gene-independent + high-n → precise). **D2 γ as a latent** (full
  errors-in-variables Bayesian IV: sample γ,β jointly — correct propagation). **D3 shared across states too?** (γ may
  differ tumour vs NAT). *Diag:* D1 vs D2 — with n~1000 patients the γ posterior is tight, so D2's extra width is
  likely negligible; confirm on one segment before paying for it.

### E. Scale / gauge commensurability (the crux — do both channels see the SAME β?)
- The mRNA channel runs on C-residualised **z-scored** predictors (`r=−resid(Y|C)`); π̂ is on CN→Y units. β is shared
  only if both are on one scale. **E1 residualise+z-score CN consistently** so `γ_{f,s}·β_f` lands on the mRNA β scale
  (estimate γ on the same z-scored X). **E2 explicit rescaling gauge constant** (risky — identifiability). **E3
  standardise everything to a repression-per-SD-dose gauge** (Decision-C flavour). *Diag (must-run first):* on a
  clean single-family strong-F edge (miR-16→CCND1), does the CN-channel-only posterior mean β equal the mRNA-only
  posterior mean β? If not, the gauge is wrong — this is the make-or-break check before any multi-family claim.

### F. Latent hierarchy (which levels β has)
- **F1 family-only β** (collapse as now) + CN channel on families — the DEFAULT, and correct for repression: same-seed
  members repress identically (§8) ⇒ **there is no per-arm repression β**, so no repression-δ pooling. **F2 arm-level
  *dose-delivery* δ-pooling** `M_member = M_fam + δ_m` (NOT repression) — the **parked Decision-F item** / hierarchical
  δ-pooling: a shrunk per-member δ with **posterior width** (wide/collapsing for the ~8% `|ρ|≥0.95` co-transcribed
  pairs), the Bayesian upgrade of the frequentist L2 point (chimeric/abundance, co-primary). A **separate** latent from
  `β_f`, orthogonal to the CN channel — it grades *resolution*, not repression, so it does NOT violate §8. **F3 +
  program level above gene** (pooling `β_f` across a program). *Diag:* F2 is worth building where a segment has ≥2
  same-seed members with divergent chimeric/dose; it never splits `β_f`. (Only "arm-level *repression* β" is retired —
  that would violate §8; dose-delivery δ-pooling is legitimate and parked.)

### G. Inclusion-mode coupling (does CN touch z, or only θ?)
- **G1 CN affects magnitude θ only** (cleanest; coupling readout). **G2 CN also lifts inclusion π** — a causally-
  confirmed edge → higher PIP, a *new discovery booster* (the instrument as evidence-of-existence, the "existence
  gate" framing in the arc memory). *Diag:* in the evidence-π (discovery) readout, does adding G2 promote any FDR-
  borderline edge that CPTAC/Manakov then supports? (ties to the discovery lane, not just coupling.)

### H. Coding-gene pleiotropy (the sibling Gap-2B item) — SAFE GATE BUILT + PRODUCTIONIZED (`instrument.coding_pleiotropy`, MH-99, 2026-07-11)
- **STATUS: validated + safe.** Conditioning miR-16→CCND1's exclusion `δ_s` on **RB1@13q14** (co-deleted, co-varying
  +0.57) absorbs **42%** (35% single-gene); miR-16→**BCL2** 25%. Correct "co-located" criterion = **CN CO-AMPLIFICATION**
  (coding CN co-varies with the miRNA-locus CN), NOT genomic distance — RB1 is 1.7 Mb away, missed by a window, caught by CN-coamp.
- **THE SAFETY BACKSTOP = "trust only reductions" (NOT the OmniPath filter).** OmniPath-direct is too restrictive — RB1 is a
  2-hop CCND1 regulator (RB→E2F→CCND1), absent from OmniPath's curated direct edges, so it would DROP the validated confounder.
  Instead: a valid backdoor-blocking conditioning **REDUCES** `δ_s` (→ a coding `down_weight`), whereas an **INCREASE** is a
  collider artifact (X_c a collider between Z_s and Y) → **IGNORE** (keep `δ_s`). Plus a proliferation-survival candidate filter.
  Result: 13q14→CCND1 down_weight **0.44** (RB1/CAB39L), the POLR2A/AURKB locus (over-conditioning) correctly **collider_ignore→0**.
  Emits per-segment `coding_down_weight`∈[0,1] = the CN-channel down-weight (inflate `s²_π`), **NOT a subtraction** of the miRNA effect.
- **⚠ Probing this exposed a MODEL-WIDE finding** (memory `proliferation-confound-vs-mechanism`): the over-conditioning came from
  proliferation markers (POLR2A/AURKB) still Y-correlated after C's `mal_prolif` ⇒ `mal_prolif` is INCOMPLETE — but the fix is NOT
  "control proliferation more completely," because proliferation is BOTH confounder AND mechanism (miR-15/16→CCND1 acts through it),
  so over-controlling proliferation-effector targets removes real signal. Right question = per-gene appropriateness. (Design §B / axis AL.)
- **✅ (a) genome-wide gene CN DONE** — `build_genomewide_gene_cn.py` (ASCAT allele-specific-seg × GENCODE intersection →
  `output/matrices/gene_cn_genomewide.tsv.gz`, **20006 coding genes × 1060 participants**); `_gene_cn()` prefers it. **✅ (b)
  WIRED into `exclusion(coding=True)`** — per segment emits `coding_down_weight` + `s²_π_coding = s²_π/(1−dw)` + `beta_weight_coding`
  (13q14→CCND1: weight 0.5→0.3). **⬜ (c) the proliferation-signal revisit** — the model-wide finding, still open (memory).
- **✅ HOST fold-in + PER-LOCUS + FULL SCAN (MH-101, 2026-07-11/12).** The intronic **HOST** gene is folded into `coding_pleiotropy`
  as a mandatory **guard-exempt** candidate (retiring the earlier separate-gates-via-MAX heuristic → ONE joint down-weight that
  de-double-counts collinear confounders). Host is now **PER-LOCUS** (`genomic_context.locus_host_map`→`locus_context.tsv`): each
  CN locus is conditioned on the host it ACTUALLY sits in (not the arm's coding-first representative), which dropped 2 silent-locus
  spurious down-weights (miR-30c→NFYC, miR-194→BPNT1). The former host-family `do_scan` restriction (a speed workaround) is **LIFTED**
  (`do_scan = coding`) now the scan is 17× faster → the generic coding scan runs on EVERY family: **3059 pure-coding + 36 host**
  down-weights (~18% of strong segments carry a co-amplified coding confounder — the coding subset of MH-89's 73% exclusion-failing rate).
- **✅ GRADE the ambiguous host down-weights (`grade_host_downweights`→`host_downweight_grades.tsv`).** A `repression_direction`
  down-weight (δ_s opposite-sign to γ) is genuine-host-confound-vs-residual-miR-leak ambiguous; resolved by joining `host_confound_lens`'s
  OOF `dC` (does conditioning on the host help/hurt held-out coupling), graded on |dC|≤0.02=immaterial: of 26 → **6 confound_confirmed
  · 6 over_control (relax candidates) · 11 host_immaterial · 3 unresolved**.
- The segment sum `Σ_{f on s} γ_{f,s} β_f` assumes only miRNA families on the segment move Y. **H1 flag-only** (as MH-88 remaining).
  **H2 add the co-amplified coding gene as another term** in the segment sum (a non-miRNA β_coding with its own γ),
  so horizontal coding-pleiotropy becomes just another channel term; OmniPath/DepMap supplies which coding genes on
  the segment regulate Y. *Buys:* the exclusion restriction's biggest threat handled inside the same object. *Diag:*
  does adding the coding term shrink the miRNA β's on segments carrying a known coding regulator (e.g. an oncogene
  co-amplified with the cluster)? **H3 combine the WHOLE co-amplified coding block on the segment as ONE nuisance factor**
  (its expression PCs) instead of cherry-picking known regulators — catches the confound *agnostically* (regulators
  OmniPath misses), **potentially stronger** than H2; *cost:* over-conditioning risk (must not soak up the family's own
  signal). **H4 (target-side) pool the SHARED pleiotropy across the target genes a segment hits** as a segment-level
  confounder latent — valid *if* the pleiotropy is low-rank (e.g. a proliferation/aneuploidy factor driving many
  cell-cycle targets). *Note — distinct from the RETRACTED "neighborhood" idea (§2):* that pooled the **signal** (β)
  and was wrong; H3/H4 pool the **nuisance**, which is legitimate.

### I. Other channels (map the pattern now; build later, each gated)
- **I1 protein (CPTAC)** — abundance→protein link as a second *likelihood* (translational residual); **non-Gaussian
  link ⇒ may break Gibbs conjugacy** (see J). **I2 state (GTEx M_H as prior mean, M_T=M_H+Δ)** — the healthy channel
  (Decision H). **I3 binding (chimeric/K_D)** — already `priors.μ`; the frame just *names* it a channel. *Diag per
  channel:* same as CN — does it move a flat/non-redundant direction, and die under its own shuffled null?

### J. Inference engine
- **J1 extend the existing Gibbs** (valid while every channel is Gaussian-conjugate — CN is). **J2 HMC/NUTS
  (numpyro/Stan)** once a non-conjugate link (protein occupancy, saturating) or the full field (A3) enters. **J3
  EB/variational** for scale. *Diag:* stay on J1 until a channel's link is non-conjugate OR the joint field is too
  big for Gibbs mixing — then J2. (Don't pre-pay for HMC.)

### K. Parallel frequentist + prerequisites ("want both")
- **K1 keep MH-88 as the fast point-estimate readout alongside the Bayes** (interim = fast/point; Bayes = full/width;
  run in parallel, as elsewhere in the program). **Prereqs to enable the map at all:** (i) add `π̂ + s²_π (SE)` and the
  `γ` solo first-stage to `family_causal_attribution`'s output; (ii) a **segment → genes-it-hits** index (for A2/A3);
  (iii) a **coding-gene-on-segment → target** map (OmniPath/DepMap, for H2). *These are cheap and unlock several axes
  — the natural first code, independent of which options win.*

**Reading the map:** the make-or-break is **E (gauge)** — run that single-edge equality check *first*; everything
else is a live menu. The cheapest unlocking work is **K's prerequisites**. The biggest conceptual fork is **A
(unit)** — per-gene vs neighborhood vs field — and it's decided by *diagnostics* (A1/A2/A4 comparisons), not upfront.

---

## 7. Channel catalog (axis I to full depth — nothing pruned)

Every source that could FEED the latent M, mapped on the same template: **latent it observes · link · Gibbs-
conjugate? (stays J1 vs forces J2 HMC) · which direction it identifies (is it non-redundant/exogenous?) · coverage
· its own shuffled null**. Conjugacy is the load-bearing property — it decides whether a channel is a cheap
Gaussian term or forces the engine to change (axis J).

| channel | latent observed | link | Gibbs-conjugate? | identifies | coverage | null |
|---|---|---|---|---|---|---|
| **mRNA coupling** (current) | `r = X·β` | linear | ✅ (is the sampler) | abundance-driven coupling | all HE edges (~1000 genes) | shuffle Y |
| **evidence prior** (current) | π/μ/τ (a PRIOR, not a likelihood) | — | ✅ prior term | a-priori existence + magnitude scale | all edges | shuffle evidence |
| **CN dose** (this doc) | `π̂_s = Σ_f γ_{f,s} β_f` | linear (reduced form) | ✅ conjugate Gaussian term | **exogenous causal + the flat collinear-family direction** | well-instrumented families (F>10; subset of genes) | shuffle CN labels |
| **protein** (CPTAC) | β on protein: `prot = a·β·(mRNA-repr) + transl. residual` | translational (mRNA→protein discordance; ~linear-ish but noisy, possibly saturating) | ⚠️ **borderline** — linear-Gaussian approx keeps J1; a discordance/saturating link forces **J2** | the **realized-at-protein** repression the mRNA channel can't see | CPTAC cohort (n-limited, ~120) | shuffle protein |
| **state / healthy** (GTEx, NAT) | `M_T = M_H + Δ` | prior-mean nesting (Decision H) — a **prior channel**, or a 2nd mRNA-likelihood on healthy Y | ✅ (Gaussian prior on Δ) | **acquired vs constitutive** wiring | GTEx-breast ~300 / NAT ~104 | shuffle state labels |
| **binding** (chimeric / K_D) | L2 arm split + μ magnitude | affinity → prior μ (already) + L2 within-family | ✅ prior term | **arm identity** + affinity ordering | K_D genome-wide 746 arms; chimeric sparse (Manakov HEK / TarBase) | shuffle arm labels |
| **methylation rewiring** | edge on/off: π depends on promoter β | state-conditional **inclusion** channel (lifts/kills z) | ✅ (modifies π, Beta-Bernoulli) | **lost/gained specialists** (epigenetic switch) | HM450 tumour/normal | shuffle β |
| **target-side SNV / APA / editing** | edge support: site created/destroyed | reference-blind **existence** channel on z | ✅ (modifies π) | **reference-blind real** edges (seed gain/loss, 3′UTR exposure) | SNV/VEP; APAatlas | shuffle SNV/site |
| **AGO loading** | arm functional-abundance reweight | within-hairpin passenger-exclusion only (cross-arm INERT) | ✅ (a prior weight) | passenger QC only — **found coupling-inert** ([[ago-loading-arm-axis]]) | GSE censuses | shuffle arm |
| **AGO capacity / ceRNA** (`ago_gate`) | per-sample RISC saturation | **nonlinear** `a/(1+Σa)` on X→repression | ❌ **non-conjugate → J2**; the occupancy substrate that **FAILED the held-out gate** | global competition / ceRNA | per-sample | — |

**Per-channel notes kept in play (not verdicts):**
- **protein** — ⚠ **THE ROW ABOVE AND THIS NOTE ARE SUPERSEDED (2026-07-12) → `CPTAC_PROTEIN_CHANNEL_PLAN.md`.**
  Two corrections, both verified: **(1) "a discordance link is non-conjugate → J2" is WRONG.** Discordance is an
  **ADDITIVE** term (`P = a·M − X·βᵗ`), hence **linear ⇒ Gaussian-conjugate ⇒ STAYS ON GIBBS.** Only *saturation*
  (L2 `a/(a+K)`) breaks conjugacy — and that link already failed its held-out gate on mRNA. **The protein channel
  does NOT move the engine.** (2) **Protein is NOT a coupling lever:** its Fisher information about β is
  `a_g²·(n_p/n_m)·(σ²_m/σ²_p)` = **≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier ****≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**)** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**) (measured under `build_C("cptac")`: `a_g`=0.397, σ² ratio 0.79, n=101 vs
  1041) and **≤7.6% even at `a_g`=1.0** — converging with the state session's independent **0.6%** (MH-102d).
  Its value is the **new latent βᵗ** (translational repression — invisible to RNA-seq at any n), estimated by
  **Model 3** (`P ~ X + M + C`) *cohort-internally in CPTAC*, not as a `channels=` term on the TCGA β (protein has
  no instrument, and at p/n≈0.63 the reduced-form vector is unestimable — it must be a **second likelihood block**).
  - *Historical (the superseded framing):* "the fork is link form — plain-Gaussian stays on Gibbs; the honest
    discordance model is non-conjugate → J2; fit linear first and graduate." The premise (discordance ⇒ non-conjugate)
    does not hold.
- **state** has two homes (axis Q): a **prior** (M_H → prior mean, Gaussian, Gibbs-cheap) or a **second likelihood**
  (healthy Y observed, so healthy and tumour β co-estimated). Map both; the prior form is the Decision-H subtype
  flavour already validated.
- **methylation / SNV / APA / editing** are all **existence/inclusion** channels (they change *whether* an edge is
  on, per state/patient) — they enter through π (the evidence-π discovery readout), NOT θ. They're the
  `WHATS_NEXT §EPIGENETIC-REWIRING` + `§PARALLEL-VIEW SNV leg` items, unified as "support-changing channels."
- **AGO capacity / ceRNA** is the one channel that is **intrinsically non-conjugate** (saturating) AND already
  failed its gate — kept in the map for completeness, flagged do-not-build-first.

---

## 8. Further axes (the dimensions §6 didn't cover — still mapping, nothing pruned)

### L. Link function form (per channel)
- **L1 linear `a·w`** (current, all channels) · **L2 saturating occupancy `a/(a+K)`** (retires D(m); FAILED the
  held-out gate for mRNA — kept for protein/capacity where saturation is real) · **L3 logit/threshold gauge** ·
  **L4 channel-specific links** (mRNA linear, protein saturating, CN linear). *Diag:* per channel, does the
  nonlinear link beat linear out-of-fold? (mRNA already says no; protein/capacity untested.)

### M. Prior-vs-likelihood placement of each source (the double-counting guard)
- The same physical fact can enter as a **prior** or a **likelihood** — never both. Evidence, binding = priors;
  mRNA, CN, protein = likelihoods; the tricky ones: **CN "existence"** (a causal edge is more likely real) could
  lift π (prior-side, axis G2) *or* stay purely in θ (likelihood-side, G1) — **not both** (double-counts the
  instrument). Methylation/SNV are likelihood-on-π. *Diag:* map each source to exactly one slot; check no source's
  PMID/observation is counted in both a prior and a likelihood term (the §6b non-circularity discipline, extended).

### N. Confounder block C — shared vs channel-specific
- **N1 one shared C** (CPE+target-CN+malignant-prolif, the §6b core) across all channels · **N2 channel-specific C**
  (CN reduced form residualises CN on C; protein needs a protein-composition C; healthy-state C differs) · **N3
  minimal-C per channel** (over-adjustment/collider risk differs by channel). *Diag:* does a channel's coupling
  survive its *own* correct C? (e.g. the CN channel must residualise CN and Y on the same C the mRNA channel uses,
  or β is on two scales — ties to E).
  - **NB (MH-100):** the shared C's proliferation term (`mal_prolif`) is materially INCOMPLETE (a residual-proliferation
    axis carries |ρ|≤0.83 with the target after C), yet a fuller *global* control over-controls the majority of genes
    (50% over_control vs 18% confound on a 28-gene panel; mean Δ_Cspace −0.032) → **C not changed**; the per-gene
    confound-vs-mechanism verdict is axis **AL** (`prolif_verdict.py`, a flag not a subtraction).

### O. Explicit prior covariance Σ over the latent (beyond what CN induces)
- The CN channel induces cross-family covariance *from data*. Independently, we can put a **prior** Σ on the family
  β's: **O1 independent** (current — diagonal) · **O2 block by co-location** (co-amplified families correlated a
  priori) · **O3 block by seed similarity / shared program** · **O4 learned Σ** (Wishart). *Buys:* structure where
  data is thin; *Costs:* a wrong Σ biases. *Diag:* does an informative Σ change any between-family share vs the
  data-only (CN-induced) covariance? If not, keep O1.

### P. Missing-data / channel coverage (the (edge × channel) matrix is sparse)
- Not every edge has every channel (CN only well-instrumented; protein only CPTAC; chimeric sparse). **P1 drop the
  absent term** (the posterior naturally does this — an edge with no CN term = mRNA-only) · **P2 partial-pool across
  coverage** (borrow the channel's typical effect for uncovered edges via a hyper-prior) · **P3 impute** (risky).
  *Diag:* map the coverage matrix first; P2 only matters if a channel covers <½ the edges it's relevant to.

### Q. State: a DIMENSION of the latent vs a CHANNEL
- **Q1 state as a channel** (healthy = prior on tumour, §7). **Q2 state as a latent dimension** `M_{m,g,state}` with
  a smoothness/nesting prior across GTEx→NAT→tumour (the field is (arm×gene×**state**)). *Buys (Q2):* the
  abundance-vs-wiring decomposition (Δx·M vs x·ΔM) becomes a posterior readout, not a separate analysis. *Diag:*
  Q2 only if wiring (ΔM) is real — the subtype flavour said mostly-conserved, so gate on the state ΔM being nonzero.

### R. Direction / sign box
- **R1 half-normal θ≥0** (repression only, current) · **R2 signed θ** for edges a **direction channel** (TarBase-v9
  Negative/Positive) flags as validated-activating (the 6 edges found, `evidence/direction.py`) · **R3 per-edge sign
  prior** from the direction channel. *Diag:* do the 6 activating edges' posteriors want positive β? (cheap check.)

### S. Cooperativity / interaction structure (Decision G)
- **S1 additive β only** (current) · **S2 primed product terms** for arm-pairs with distinct sites 8–40 nt apart ·
  **S3 shared-AGO pool** as the global-saturation mechanism. *Costs:* S2/S3 sit on the occupancy substrate (failed).
  *Diag:* gate hard on held-out improvement; do not build before mRNA-linear is beaten.

### T. ceRNA / global competition
- **T1 none** (total miRNA is RPM-constant, CV 0.05 — competition enters via capacity not load) · **T2 per-sample
  AGO-capacity latent** (`ago_gate`, the saturating gate) · **T3 explicit free-pool term**. *Diag:* the gate HURT
  coupling at every rung (MH-79) → keep T1 unless a protein/capacity result demands otherwise.

### U. Cross-channel DISAGREEMENT as a first-class readout
- With channels, **agreement/disagreement is a finding**: **U1 per-edge channel-provenance vector** (mRNA-only /
  CN-confirmed / protein-confirmed / discordant) · **U2 posterior-predictive check per channel** (does the fused
  posterior predict each channel's held-out data?) · **U3 discordance genes** (mRNA-coupled but protein-flat =
  translational buffering; CN-null but mRNA-coupled = confounded/shared-program). *Buys:* the discordances ARE the
  biology (the §6b "one posterior" gains a provenance dimension). *Diag:* none needed — this is pure readout, free.

### V. Readout / attribution space (what the fused posterior exposes)
- **V1 the four §6b jobs** (coupling/attribution/identity/discovery) · **V2 + channel provenance** (U1) · **V3 +
  per-channel contribution decomposition** (how much each channel moved each β — a Shapley-over-channels) · **V4 +
  causal vs observational split per edge** (the CN channel's share = `cn_driven_fraction`, now with width). Map: each
  is a marginal of the same posterior.

### W. Compute / tractability budget (per A×J combination)
- Gibbs per-gene (A1×J1) ≈ current cost · joint neighborhood (A2×J1) = block sampler over ~5-gene groups · full
  field (A3) or non-conjugate (J2) = HMC, hours. *Diag:* profile A1 first (memory: `profile-before-batch-loops`);
  the field is only worth its compute if A2 diagnostics say pooling pays.

### X. Evaluation reference per channel (which bar scores which)
- mRNA → OOF held-out patients (Bar 1) · CN → the MH-88 interim + shuffled-CN (Bar 4) · protein → CPTAC discordance
  (Bar 5) · state → independent cohort (Buffa/METABRIC) · binding → Manakov/CLIP. Map: **each channel is scored
  against the reference it is not derived from** (non-circularity, per channel).

**Coverage note:** this map now spans *sources* (§7 catalog), *where each enters* (M prior/likelihood), *its link*
(L), *its confounders* (N), *the latent's shape* (F hierarchy · O covariance · Q state-dim), *missingness* (P),
*sign/interaction/competition* (R/S/T), *readouts* (U/V), and *engine/compute/eval* (J/W/X). Open a new axis here
rather than collapsing one when a new possibility appears — the map grows, it doesn't prune.

---

## 9. Further axes, round 2 — sample structure, measurement, units, causal graph, meta/inference

Grouped for navigation; still mapping, nothing pruned. New axes lettered Y → onward.

### Sample / population structure
- **Y. Patient heterogeneity.** Y1 population-level β (current) · Y2 **patient/subtype random effects** `β_{m,g,i}=β_{m,g}+b_i`
  (some tumours wire an edge, others don't) · Y3 **mixture over patient classes** (responder/non-responder latent
  class per edge — ties to the "spiker" arms, e.g. miR-135b active in a subset). *Diag:* does an edge's coupling
  concentrate in a patient subset (bimodal residual)? If yes, Y3 earns its place.
- **Z. External cohorts as CHANNELS (not just validation).** Z1 external held-out only (current — Buffa/METABRIC =
  Bar) · Z2 **each cohort a second mRNA likelihood** with a cohort random effect (fold replication INTO the
  posterior; borrow strength) · Z3 cohort-specific β (test transportability). *Diag:* does folding METABRIC in as a
  channel sharpen β vs using it only as a held-out bar?
- **AA. Batch / technical nuisance.** AA1 post-fit check (current, ≤0.03) · AA2 in-model nuisance latent (plate/site
  as a random effect) · AA3 covariate-protected ComBat pre-fit. *Diag:* the batch effect was negligible (MH-46) →
  AA1 unless a channel (protein plex, cohort) reintroduces it.

### Measurement models (the likelihood's distributional form) — **most-suitable distribution PER channel**
Each channel's *data type* dictates its distribution; this is not cosmetic and **it hits the dense model itself**
(the mRNA likelihood is one of these rows), so decide it per channel, not once:

| channel | data | most-suitable | why / when it matters |
|---|---|---|---|
| **mRNA (the DENSE §6b likelihood)** | log-CPM residual | Gaussian (approx) · **NB** (count-honest) · **t** (robust) | Gaussian-on-log-CPM is the *voom-style approximation*; NB is the count-honest generative model (matters for **low-expression** genes); **t** is the cheapest real upgrade (heavy-tailed residuals from outlier tumours). **So the distribution question is the dense model's own, not a channel add-on.** |
| **CN** | a *coefficient estimate* `π̂_s` | **Gaussian → t** | asymptotically normal (CLT, n~1000); **t** is honest that `s²_π` is *estimated* + guards collider outliers |
| **protein (CPTAC)** | log-ratio abundance | log-normal + **non-linear discordance link** | the *link* (translational residual) is the issue, not the noise family |
| **binding (chimeric/CLIP)** | read **counts** | **Poisson / NB** (or binomial bound/not) | counts ⇒ not Gaussian |
| **methylation** | β ∈ [0,1] | **Beta** (bounded) or Gaussian on M-values | bounded support |

- **AB. Verdict — RESOLVED (MH-92, 2026-07-10): Student-t, NB rejected.** The diagnostic ran (85 genes × 1041 pts, exact
  dense fit): the dense mRNA residuals are **heavy-tailed** (excess-kurt 1.2 median / 3.5 pooled, t-MLE **df≈7**, P(|z|>4)=33×
  Gaussian) from amplified-subset (+skew, ERBB2/3) & near-floor (−skew, 3 genes) tumours — **NOT** count-like at low expression
  (Spearman(expr,kurt)=+0.19, wrong sign for NB). **Student-t** is now BUILT in `spike_slab._gibbs_posterior(nu=)` (scale-mixture,
  stays Gibbs; `nu=None` bit-identical). It is a **correctness/robustness upgrade, not a lever** at n≈1000 (nests Gaussian;
  robust 1.76× under contamination; OOF Δρ +0.001 n.s., attribution corr 0.999) — the same "principled-but-immaterial-at-this-
  depth" shape as the CN channel (§2.3) and the gate SEs (bootstrap confirms Sobel adequate, s²_π δ²-floor-dominated). t also
  applies to the CN channel term (`t['nu']`, guards collider-outlier segments) and NB stays reserved for the future count
  channels (binding). Non-Gaussian (counts/Beta/saturation) still forces HMC (axis J). See METHODS §3b (twins), WHATS_NEXT ①,
  VALIDATION, registry/ledger MH-92.
- **AC. Predictor X measurement / compositionality.** miRNA RPM is **compositional** (sums to ~constant, CV 0.05).
  AC1 log2(RPM+1) treated as free (current) · AC2 **CLR/ILR compositional transform** (honest simplex geometry) ·
  AC3 spike-in / absolute-scale if available. *Diag:* does the compositional transform change any coupling sign?
  (the near-constant sum says small, but the ceRNA/competition story lives exactly here — ties to T.)
- **AD. CN measurement granularity.** AD1 continuous copy number (current) · AD2 **ordinal / categorical** (amp /
  neutral / del) · AD3 **allele-specific / parent-of-origin** (ASCAT gives it; imprinted loci, DLK1-DIO3) · AD4
  focal-vs-arm-vs-broad CN level (the event scale). *Diag:* does allele-specific CN change the instrument on
  imprinted clusters? (the methylation lane already flags DLK1-DIO3.)

### Units & pooling
- **AE. Regulator unit / pooling granularity.** AE1 seed FAMILY (current) · AE2 seed hexamer / 8mer · AE3 arm ·
  AE4 hairpin/locus · AE5 cross-species-conserved-family as a distinct pooling attribute. (Interacts with F
  hierarchy — F is the *levels*, AE is the *unit definition*.) *Diag:* the family-key was already a bug-source
  (MH-86 miR-number vs seed) — pin the unit before pooling.
- **AF. Program / Hallmark pooling.** AF1 gene-independent (current) · AF2 **program as a pooling level** (μ_program
  — a miRNA's typical weight within a Hallmark) · AF3 multi-membership handling (a gene in N sets — weight-share vs
  duplicate-row). *Diag:* do genes in a shared program have correlated M for a shared regulator beyond chance?
- **AG. Cross-gene M gauge.** Is `β_{m,gA}` comparable to `β_{m,gB}` (a miRNA's weight ACROSS its targets)? AG1
  per-gene z-scored (not comparable across genes) · AG2 a shared dose gauge (comparable — enables μ_m pooling to
  mean something) · AG3 target-normalised. (Distinct from E, which is per-edge cross-CHANNEL; this is cross-GENE.)
  *Diag:* the μ_m pooling in `hierarchical.py` assumes cross-gene comparability — is it actually gauged?
- **AH. Compartment-specific M.** AH1 bulk (current + deconv-in-C) · AH2 **compartment-specific latent**
  `M_{m,g,compartment}` (cell-intrinsic vs stromal — the miR-29→ECM stromal class, MH-51) · AH3 single-cell/spatial
  channel to identify it. *Diag:* deconv-in-C already separates cell-intrinsic; AH2 only if a compartment carries a
  distinct M sign.

### Causal graph / identification
- **AI. Instrument set (CN is one).** AI1 somatic locus CN (this doc) · AI2 **miRNA-locus SNV / germline eQTL**
  (MR-style, cross-individual) · AI3 trans-eQTL / methylation-as-instrument for the miRNA · AI4 **multiple
  instruments jointly** (CN + eQTL over-identify the same β — a stronger Sargan/Hansen test). *Diag:* is a second
  instrument available and non-redundant with CN? (eQTL breaks reverse-causation differently than somatic CN.)
- **AJ. Full DAG / structural model.** AJ1 single-instrument reduced form (this doc) · AJ2 **explicit SCM**
  (CN→miRNA→mRNA→protein, with target-CN→mRNA and TF→{miRNA,mRNA} confounders as nodes) · AJ3 mediation split
  (CN→miRNA→mRNA vs CN→co-amplified-coding→mRNA — the H2 coding-pleiotropy term, as a formal mediator). *Diag:* does
  writing the DAG expose a backdoor the current C doesn't block?
- **AK. Feedback / reverse causation.** AK1 assume acyclic (current) · AK2 model miRNA↔target feedback (target is a
  TF/host of the miRNA — the instrument is exactly what licenses ignoring this, so map where it's NOT ignorable) ·
  AK3 flag feedback edges (host-gene/intragenic miRNAs). *Diag:* which HE edges are intragenic (miRNA in the
  target's intron)? those need the instrument most.
- **AL. Target's own transcription-rate latent (Decision B, "most load-bearing").** AL1 proliferation proxy only
  (current) · AL2 **intronic/pre-mRNA nascent-transcription latent** (orthogonal to mature mRNA — removes TF-co-drive
  positive-weight artifacts) · AL3 target-CN + target-methylation as the target's genomic set-point. *Diag:* the
  TF-regulon proxy over-controlled; does a nascent-transcription proxy remove positive-β artifacts WITHOUT collapsing
  real signal? (HIGH value, needs intronic-read data — check availability.)
  - **VERDICT (MH-100, `learned/prolif_verdict.py`) — AL1 is INCOMPLETE, but AL2's premise ("control more completely")
    is REJECTED, and the resolution is PER-GENE.** The residual axis `R = resid(broad cell-cycle metagene | C)`
    (target-excluded) still carries **|ρ|≤0.83** of the target after C (mean 0.42; AURKA +0.83 / TGFBR2 −0.64 / PTEN,ZEB1
    −0.54) ⇒ `mal_prolif` (E2F+G2M) is materially incomplete. **But this is the one place a modeling choice is a real
    coupling LEVER** (adding R moves β by median **−20 to −30%**, vs the sub-σ of the CN channel §2.3 / Student-t AB /
    gate-SEs) — and it acts **three ways at once**: de-confounds (miR-130/301/454→ZEB1 −97%, canonical miR-200bc/429
    spared −11%), over-controls proliferation effectors (MYB 60–89%), and **un-masks** real edges (let-7→CDK6 +350%,
    miR-25/92→PTEN +63% — proliferation was a *negative* confounder hiding them). **Non-circular arbiter = OOF 2×2**
    (train {C, C+R} × eval {C-space, C+R-space}): control helps both → confound; hurts both → mechanism/over-control;
    flips → fragile (OOF can't settle — the disagreement is about the estimand, not predictive fit). **28-gene panel:
    over_control 50% · confound 18% (PTEN/CDK4/VEGFA/THBS1/CDKN1B) · fragile 7% · neutral 25%; mean Δ_Cspace −0.032
    (net harm)** ⇒ a global fuller control over-controls the majority → **C left as-is (row-10 unchanged).** **AL2 is
    DATA-BLOCKED** — no intronic/pre-mRNA reads in the repo. Deliverable = a per-gene **FLAG, not a subtraction**
    (sibling of the coding gate's "trust only reductions", §H) → `output/learned/prolif_verdict.tsv`. Proliferation is
    BOTH confound and mechanism; the right object is per-gene appropriateness, not "capture better." (Ties axis **N**.)
  - **GENOME-WIDE + PER-EDGE + ARM-SOURCE (extension, 8-worker parallel — whole 1507-gene / 1191-coupled-edge universe in
    4.7 min).** (a) **Most genes are proliferation-NEUTRAL (980/65%)** — the panel's 50%-over_control was proliferation-
    ENRICHED; among the 527 affected, over_control beats confound **2.8:1** (360 vs 129; net Δ_Cspace −0.017) → the verdict
    HOLDS genome-wide. (b) **Per-EDGE** (`prolif_verdict_edges_genomewide.tsv`): mechanism 430 (36%) · ambiguous 352 (30%) ·
    **prolif_robust 236 (20%)** · **deconfound 88 (7%)** · **unmasked_real 85 (7%)** — controlling proliferation over-strips
    mediated signal (36%) far more than removes genuine confound (7%); 7% of edges are actively HIDDEN by proliferation
    (discovery leads); 20% robust. (c) **Arm-SOURCE (axis F / §8 dose-delivery, `arm_prolif_source`):** a family's
    proliferation confound is NOT diffuse — it localizes to ONE arm (**median concentration 0.85**), systematically the
    arm **miR-93-5p** (106b~25) across dozens of miR-17~92 targets. **MECHANISM TESTED (`confound_mechanism`) — HOST-GENE
    transcriptional, NOT CN (corrects an earlier "CN-localizable" assertion):** miR-93/106b/25 is INTRONIC in **MCM7** (a
    replication gene), so it is co-transcribed from the host promoter; conditioning on host MCM7 removes **65–97%** of the
    confound vs the locus CN only **8–10%** (miR-93 host-frac 0.72 / CN-frac 0.08; first stage ρ_abund,CN≈0.25 is too weak).
    ⇒ this is **intragenic-host co-transcription (axis BP)**, and the handle is the **host mRNA, not the CN instrument** — the
    CN channel/T1 exclusion barely moves it. Sign-disagreeing arms (miR-15a/424 + vs miR-195/497 −) EXPLAIN the fragile
    verdicts. This grades the confound SOURCE (which arm) + its MECHANISM (host transcription), NOT the verdict (family-β, §8).

### Meta / inference strategy (how we navigate the map itself)
- **AM. Model averaging vs selection over the map.** AM1 pick one option per axis (a single model) · AM2 **Bayesian
  model averaging** over live options (integrate over the axis — e.g. average over unit A1/A2, over link L1/L2) ·
  AM3 stacking / ensemble by held-out weight. *Buys (AM2/3):* no premature pruning — the map's uncertainty enters
  the posterior. *Costs:* compute × #options. *Diag:* do the top options disagree on any conclusion? if not, select.
- **AN. Sequential vs simultaneous fusion.** AN1 joint fit (all channels at once) · AN2 **sequential Bayesian
  updating** (mRNA posterior → CN update → protein update; each channel refines the prior) · AN3 modular/cut
  (prevent a mis-specified channel from contaminating others — "cut" feedback). *Buys (AN2):* interpretable "what
  each channel added"; *Costs:* order-dependence unless done right. *Diag:* does channel order change the posterior?
  (it shouldn't if the joint is well-specified — a consistency check.)
- **AO. Prior/channel sensitivity as a deliverable.** AO1 none · AO2 **leave-one-channel-out** (how much does each
  conclusion depend on each channel — the fusion's robustness) · AO3 prior-sensitivity sweep (π/μ/τ, s²_π). *Buys:*
  the honesty layer — a conclusion that flips when one channel is removed is fragile. *Diag:* free once the fusion
  runs; report per headline edge.
- **AP. Identifiability as a first-class READOUT.** AP1 implicit (posterior width) · AP2 **explicit flat-direction
  report** (which linear combinations of β the data can't separate — condition number, the collinear-family
  directions) · AP3 an identifiability score per edge (how much of β is data vs prior). *Buys:* names exactly where
  a channel (CN) is *needed* — the flat directions ARE the CN channel's job. *Diag:* free; it's the diagnostic that
  tells you which edges to instrument.
- **AQ. Decision-theoretic thresholds / uncertainty propagation.** AQ1 point readouts (current cards) · AQ2 **propagate
  posterior width** into the cards / discovery-FDR / attribution (a driver call carries its credible interval) · AQ3
  utility-weighted decisions (cost of false-positive vs false-negative per lane — discovery vs validation differ).
  *Diag:* which downstream claims would change sign under their posterior width? (the ones to re-report with intervals.)
- **AR. Prior elicitation structure.** AR1 flat / weakly-informative · AR2 the graded evidence ledger (current
  π/μ/τ, PMID-deduped, assay-class-weighted) · AR3 hierarchical prior on the ledger's own weights (learn the assay-
  class weights). *Diag:* does learning the assay-class weights beat the hand-set transfection-calibrated ones?

**Round-2 coverage:** adds sample/population structure (Y/Z/AA), measurement distributional form (AB/AC/AD), unit &
pooling granularity (AE/AF/AG/AH), the fuller causal graph & instrument set (AI/AJ/AK/AL), and the meta-inference
strategy for navigating the map without pruning (AM–AR). The map is intentionally over-complete: most axes will
resolve to their cheap default via a one-line diagnostic, but each is written so the default is a *tested* choice,
not an unexamined one.

---

## 10. Further axes, round 3 — sequence biology, universe/scope, literature, validation design, engineering

### miRNA & target sequence biology (the predictor's real complexity)
- **AS. isomiR heterogeneity.** AS1 one canonical mature per arm (current) · AS2 **5′-isomiRs as distinct edges**
  (a 5′ shift moves the seed → a *different* targetome — same locus, different β) · AS3 3′-isomiRs (stability/loading,
  same seed). *Diag:* what fraction of an arm's reads are 5′-shifted (isomiR-seq)? if material, AS2 splits the edge.
- **AT. miRNA-seed A-to-I editing (retargeting).** AT1 ignore · AT2 **edited-seed as a separate regulator** (edited
  miR-376/miR-381 retarget) — a miRNA-side support channel (sibling to the target-side SNV channel, §7). *Diag:*
  editing sites in expressed arms (REDIportal) → new seed → new edge set.
- **AU. Target 3′UTR isoform / site availability.** AU1 reference 3′UTR (current — all annotated sites present) ·
  AU2 **APA/isoform-weighted site availability** (which UTR is expressed → which sites exist per sample; short-UTR
  escapes distal sites) · AU3 isoform as a per-sample multiplier on the edge. (Target-side sibling of AS; ties to the
  APA support channel.) *Diag:* APAatlas short/long-UTR ratio on the target — does it gate the site set?
- **AV. Kinetics / steady-state assumption.** AV1 steady-state cross-sectional (current) · AV2 mRNA-decay-rate-aware
  (repression acts on degradation; fast- vs slow-turnover targets respond differently) · AV3 pseudotime/trajectory
  if a dynamic axis exists. *Diag:* mostly not identifiable from cross-sectional data — map, default AV1, flag.

### Universe / scope of the field
- **AW. Regulator universe.** AW1 HE seed-families (current coupling) · AW2 HE ∪ discovered orphans (current
  discovery support) · AW3 genome-wide TargetScan/K_D nomination · AW4 CLIP-physical-site universe (ENCORI/POSTAR3).
  (The support the field is defined over — interacts with G/§7 support-changing channels.)
- **AX2. Target-gene universe.** AX2a the 50-Hallmark union (~4,384 genes, current) · AX2b genome-wide · AX2c a
  focused causal panel (well-instrumented genes only, for the CN-channel MVP). *Diag:* the CN channel only has power
  where an instrument exists → AX2c is the natural first scope for Gap-2B Bayes.
- **AY. Discovery floor / detectability.** AY1 RPM≥10 (current discovery universe) · AY2 sweep {50/10/1} (sensitivity,
  a re-filter not a re-fit) · AY3 per-state floors. (The `WHATS_NEXT §K_D` floor-sweep item, folded in.)

### Evidence / prior sourcing
- **AZ. Literature-mining channel.** AZ1 curated DBs only (miRTarBase/TarBase, current) · AZ2 **LLM/text-mined
  edge evidence** (a graded literature-support prior beyond the curated set) · AZ3 breast-context-weighted literature.
  *Diag:* non-circularity (mined evidence must not leak TCGA) + calibration vs curated tiers.
- **BA. Evidence ledger structure.** BA1 scalar weight (current PMID-deduped Σ) · BA2 **structured per-assay-class
  vector** (reporter/western/proteomics/CLIP as separate prior dimensions) · BA3 direction-typed (TarBase Neg/Pos).
  (Refines M/AR — the ledger's internal shape as it enters the prior.)

### Validation / falsification design
- **BB. Null construction.** BB1 shuffle-CN (edge-label permutation, current plan) · BB2 parametric null (γ=0) ·
  BB3 negative-control instruments (CN of unlinked loci) · BB4 negative-control genes (non-targets). *Diag:* do the
  nulls agree on the FDR? (each catches a different leakage.)
- **BC. Cross-validation scheme.** BC1 patient-level K-fold (current OOF) · BC2 gene-level held-out (generalise to new
  genes) · BC3 cohort-level (transportability, ties to Z) · BC4 locus-level (held-out instrument). *Diag:* which
  generalisation is claimed? score against that fold.
- **BD. Multiplicity across the field.** BD1 per-gene FDR (current) · BD2 field-wide FDR (all edges jointly) · BD3
  hierarchical FDR (program → gene → edge). *Diag:* the joint field (A3) changes the multiplicity object — map it.

### Engineering / reproducibility (the axes that bit us already)
- **BE. Determinism.** BE1 seeded + order-sorted (the MH-88 fix — sampled Shapley / set-iteration must be sorted) ·
  BE2 exact where feasible (small p). *Diag:* two runs bit-identical (the check we ran for `between_family_exclusion`).
- **BF. Compute strategy.** BF1 per-gene serial · BF2 parallel over genes (BiocParallel-style, the K_D-scan pattern) ·
  BF3 background for long fits (memory: `profile-before-batch-loops`, long runs → background). *Diag:* profile a
  small sample + hoist loads out of loops before any batch (the standing rule).
- **BG. Caching / incrementality.** BG1 recompute · BG2 cache γ / π̂ / K_D as reusable `(arm,gene)` resources ·
  BG3 incremental channel add (don't refit mRNA when adding CN). *Diag:* the content-hash cache pattern already used.
- **BH. Code tracking / provenance.** BH1 the `learned/` tree is currently **git-untracked** (`??`) — a provenance
  hole flagged repeatedly; BH2 `git add` the tree before the arc grows. (Not modeling, but real; map it so it's not
  forgotten.)

**Round-3 coverage:** sequence-level biology (AS–AV: isomiR / seed-editing / target-isoform / kinetics), universe &
scope (AW/AX2/AY), evidence sourcing shape (AZ/BA), validation design (BB–BD), and engineering/reproducibility
(BE–BH). Together with §6–§8 the map now covers, end to end: **what feeds the latent, how, on what scale, over what
support & universe, with what causal identification, what sample/measurement/sequence structure, read out how, fused
by what engine, validated against what, and built reproducibly.** New possibilities open a new axis; the map never prunes.

---

## 11. Full axis index (navigation)

`§6` **A** unit · **B** instrument-obs · **C** channel-noise · **D** γ · **E** per-edge gauge · **F** hierarchy · **G**
inclusion-coupling · **H** coding-pleiotropy · **I** other-channels · **J** engine · **K** parallel+prereqs.
`§7` **channel catalog** (mRNA · evidence · CN · protein · state · binding · methylation · SNV/APA · loading · capacity).
`§8` **L** link-form · **M** prior-vs-likelihood · **N** per-channel-C · **O** Σ-prior · **P** coverage · **Q** state-dim ·
**R** sign · **S** cooperativity · **T** ceRNA · **U** disagreement-readout · **V** readout-space · **W** compute · **X** eval-ref.
`§9` **Y** patient-RE · **Z** cohorts-as-channels · **AA** batch · **AB** Y-likelihood · **AC** X-compositionality · **AD**
CN-granularity · **AE** regulator-unit · **AF** program-pooling · **AG** cross-gene-gauge · **AH** compartment-M · **AI**
instrument-set · **AJ** DAG/SCM · **AK** feedback · **AL** target-transcription-latent · **AM** model-averaging · **AN**
sequential-vs-joint · **AO** channel-sensitivity · **AP** identifiability-readout · **AQ** uncertainty-propagation · **AR**
prior-elicitation. `§10` **AS** isomiR · **AT** seed-editing · **AU** target-isoform · **AV** kinetics · **AW** regulator-
universe · **AX2** target-universe · **AY** discovery-floor · **AZ** literature-channel · **BA** ledger-structure · **BB**
null · **BC** CV-scheme · **BD** multiplicity · **BE** determinism · **BF** compute · **BG** caching · **BH** code-tracking.
`§12` **BI** network/indirect · **BJ** biogenesis-capacity · **BP** host-gene-pleiotropy · **BK** sponge-nodes · **BL**
subclonal-CN · **BM** protein-complex-buffering · **BN** SBC/calibration · **BO** foundation-model-priors · **BQ** ancestry.

**~60 axes + a 10-channel catalog.** The three make-or-break / cheapest-first anchors from across the whole map:
**E (per-edge gauge)** must pass before any multi-family claim; **K prerequisites** (π̂+SE+γ, segment→genes index,
coding→target map) are option-independent unlocks; **AX2c + A1 (well-instrumented panel, per-gene)** is the smallest
honest first slice. Everything else is a live menu resolved by named diagnostics, not upfront choice.

---

## 12. Further axes, round 4 — network, biogenesis, competition, clonal, protein-structure, calibration, priors

### Network / indirect structure
- **BI. Indirect / network path effects.** The edge miRNA→gene may be a **PATH** (miRNA→TF→gene), not a direct site.
  BI1 direct-site only (current) · BI2 model indirect cascades (network propagation over a miRNA→TF→target graph) ·
  BI3 **separate direct (site+CLIP) vs indirect (no-site, network-mediated)** coupling as distinct readouts. *Diag:*
  does the coupled edge carry a physical site? no-site-but-coupled = candidate indirect (and a CN-channel test
  distinguishes: an indirect edge's CN reduced form flows through the TF, not the miRNA — a pleiotropy signature).

### miRNA biology / biogenesis
- **BJ. Biogenesis capacity as a global X confounder.** Drosha/Dicer/DGCR8/XPO5/AGO2 status scales mature miRNA
  levels globally (and is cancer-dysregulated). BJ1 X taken as given (current) · BJ2 **biogenesis-machinery
  expression as a per-sample supply-capacity covariate/latent** (a global scaling of X, the supply-side analog of
  the `ago_gate` demand-side capacity) · BJ3 as an instrument-validity flag (a biogenesis-limited locus makes less
  mature product than its CN implies → weak first stage). *Diag:* does Drosha/Dicer in C move coupling?
- **BP. Host-gene / intragenic instrument validity.** An **intragenic** miRNA shares its locus CN *and* co-
  transcription with its **host gene** — so locus CN also moves the host, and if the host regulates the target that
  is horizontal pleiotropy through the host (the coding-pleiotropy H2 case, in its sharpest form). BP1 ignore · BP2
  **flag intragenic arms + add the host as a segment term** (like H2) · BP3 use host mRNA as the first-stage
  readout (host expression = the co-transcription proxy). *Diag:* which HE arms are intragenic (miRBase host
  annotation)? their CN instrument is host-confounded and must carry BP2/BP3.

### Competition / ceRNA nodes (extends T)
- **BK. Explicit sponge / competitor nodes.** BK1 none (current) · BK2 **lncRNA / circRNA / pseudogene sponges as
  explicit competing nodes** (the target's coupling depends on sponge abundance drawing the miRNA away) · BK3 the
  target as a competitor for its co-targets (ceRNA crosstalk within the field). *Diag:* is a known sponge expressed,
  and does conditioning on it change the edge? (mostly small — total-miRNA is RPM-constant — but node-specific for
  high-copy sponges.)

### CN / clonal resolution (extends AD)
- **BL. Subclonal / clonal CN.** BL1 bulk clonal CN (current, ASCAT tumour-average) · BL2 **subclonal / CCF-weighted
  CN** (a CN event in a *minor* clone perturbs dose in only part of the tumour → a diluted instrument) · BL3 purity-
  adjusted integer CN. *Diag:* do multi-clonal tumours attenuate the first stage? (ASCAT gives clonal/subclonal;
  ties to purity in C and to BL as an instrument-strength modifier.)

### Protein / realization structure (extends the protein channel §7)
- **BM. Protein-complex buffering / stoichiometry.** Why mRNA↔protein discordance: complex subunits are **buffered**
  (excess subunit degraded) so protein doesn't track mRNA-repression. BM1 independent-gene protein link (channel
  default) · BM2 **complex-aware** (CORUM stoichiometry — a subunit's protein is buffered by its complex) · BM3 flag
  buffered genes and attenuate the protein channel's weight there. *Diag:* are the protein-flat-but-mRNA-coupled
  genes (readout U3) enriched for complex subunits? if yes, BM2 explains the discordance rather than discarding it.

### Validation sub-axes (extends BB–BD)
- **BN. Simulation-based calibration / prior-predictive.** BN1 none · BN2 **SBC** (simulate β from the prior, refit,
  check posterior-rank uniformity — does the sampler recover known truth?) · BN3 prior-predictive checks (does the
  prior generate plausible data?). *Diag:* the honesty layer for the Bayesian machinery itself — run once per model
  *form* (e.g. once the CN channel term is added), independent of the biology.

### Prior sourcing (extends AZ/AR)
- **BO. Pretrained / foundation-model priors.** BO1 mechanistic priors only (TargetScan context++ / scanMiR K_D,
  current) · BO2 **pretrained miRNA-target deep models** as a μ/π prior channel · BO3 sequence-language-model
  embeddings of the miRNA/3′UTR. *Diag:* non-circularity (the pretrained model must NOT be trained on TCGA or on the
  same HE evidence) + does it beat scanMiR/context++ on held-out breast binding?

### Population structure (extends AI/Y — only when a germline instrument is used)
- **BQ. Ancestry / population stratification.** BQ1 ignore (somatic CN needs no ancestry control) · BQ2 **ancestry
  PCs in C** (required for the germline-eQTL instrument AI2 — population structure confounds germline dose) · BQ3
  ancestry-stratified fits. *Diag:* only live if AI2 (eQTL instrument) is built; somatic-CN (this doc) is exempt.

**Round-4 coverage:** network/indirect paths (BI), miRNA-supply biology (BJ/BP), explicit competition nodes (BK),
clonal CN resolution (BL), protein-complex structure (BM), Bayesian self-calibration (BN), foundation-model priors
(BO), and population structure for the germline-instrument branch (BQ). These are the **peripheral-but-distinct**
frontier: each is a real modeling choice, most default off for the somatic-CN MVP, but each now has a named home and
a diagnostic rather than being an unexamined omission. **BP (intragenic host pleiotropy)** and **BL (subclonal CN)**
are the two here that directly touch the CN instrument's *validity* — promote them from peripheral if the MVP hits
an intragenic or highly-subclonal arm.
