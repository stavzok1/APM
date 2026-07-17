# Attribution — Identity vs Magnitude: concepts, definitions, interpretations

> **Status: §0 is the LEARNED-MODEL realization (built 2026-07-06, `learned/attribution.py`, Design
> §Decision I). §1–§12 are the PRECURSOR (softmax-share) conceptual origin — the softmax identity is now
> RETIRED (§0).** Cross-refs: `LEARNED_MODEL_METHODS.md` §16, `GENE_QUESTION_TAXONOMY.md`,
> `MODELING_FRAMEWORK.md` §5a (share ⊥ coupling), engine `pressure_engine.py`.

## 0. LEARNED-MODEL REALIZATION (2026-07-06) — the identity we actually use
The learned model **retires the precursor softmax identity** (§1, abundance-removed specialist) and exposes
two coupling-grounded objects (`learned/attribution.py`):
- **MAGNITUDE / force** = the budget share `M(f)·X_fam` (`states.budget_shift`, METHODS §5a) — who realizes
  the most repression (abundance-included).
- **IDENTITY / who** = the **Shapley/LMG** credit for the FIXED-weight aggregate's explained variance
  R²(X·M, Y), decomposed across seed-families (`attribution.identity_vs_magnitude`) — who **fairly OWNS the
  coupling**, splitting shared credit under collinearity (Shapley = LMG for the linear aggregate).

**Why the softmax identity was retired:** §12's edge reality-check showed the precursor abundance-removed
identity ⊥ coupling — it surfaced *silenced specialists* (miR-137→NCOA3, ρ=−0.02 NS) while the abundant arm
actually coupled (miR-17, ρ=−0.25). So that "identity" was a literature/structure hypothesis, not a driver.
The learned Shapley identity is **coupling-based** — it can only credit families that contribute to the *real*
coupling, so it cannot surface a non-coupling specialist.

**Why identity ≠ budget:** budget rewards **abundance** (M·X̄ = realized force); identity rewards **on-target
coupling** (does the family's pressure explain Y). Two reasons they diverge — (i) *coupling strength* (an
abundant family whose pressure barely tracks Y is a passenger), and (ii) *collinearity* (co-varying families
split the shared credit — the smaller effect).

**What it delivers** (the identity≠magnitude gap, gene-specific):
- separates **quiet owners** (high identity, modest budget — ESR1 **miR-18**: owns 0.54 of the coupling on a
  0.20 budget; scarce [abund 3.0] but its pressure tracks ESR1 strongly, individual ρ −0.43) from **loud
  passengers** (high budget, ~zero ownership — ESR1 **miR-22**: the *most abundant* [16.0] → budget 0.26, but
  its pressure barely tracks ESR1 [individual ρ −0.12] → identity 0.01). NB the ESR1 families are **not**
  collinear (abundance |ρ| ≤ 0.29) — here the gap is purely abundance-vs-coupling, not redundancy;
- the collinearity effect appears where families co-vary (e.g. ZEB1's miR-200bc/429 vs miR-141/200a);
- the identity↔magnitude rank-corr is a **budget-faithfulness diagnostic** (PTEN +0.91 = budget≈ownership;
  ESR1 +0.45 = budget over-credits loud-but-off-target arms). CLI: `python -m mirna_hallmark.learned.attribution PTEN ESR1`.

**Discipline (unchanged, §5a):** identity share ⊥ coupling — a high identity share is the *distribution* of the
coupling that exists, never evidence a driver exists; a driver call still needs FDR-significant coupling
(METHODS §6/§7). Identity says *who*, magnitude says *how much force*, neither says *whether*.

**Budget ≠ identity is NOT a weight-calibration failure.** Budget = M·X̄ = mean *force* (depends on abundance);
identity = Shapley R² = *variance explained* (depends on how the family's variation tracks Y). Abundance ≠
variance-explained, so even a perfect M keeps a loud-steady/loud-weak family high-budget & low-identity. M is
already calibrated to coupling (the NNLS maximises the aggregate's fit); the budget merely re-multiplies by X̄.
The gap is the *signal* (loud passenger vs quiet owner), not an error — which is why we keep both. (The
`indiv_rho` column exposes the mechanism: identity tracks indiv_rho, budget tracks abundance.)

**When the PRECURSOR (softmax) identity is still meaningful:** for **silenced / lost** regulators. A specialist
miRNA methylation-silenced in tumour has no variance (→ low Shapley) and low abundance (→ low budget), so BOTH
learned metrics miss it — but it is the *structural* owner (who is designed to regulate g). The abundance-removed
softmax identity surfaces it (the TSG-miRNA-loss lens, miR-137→NCOA3). It is a structural/literature HYPOTHESIS
(⊥ current coupling, §12) — use it to *find* silenced specialists, then confirm the loss separately. So:
**precursor identity = who-by-design (incl. silenced) · budget = current force · Shapley identity = current explained-variance.**

**The silenced-specialist lens is now BUILT (2026-07-06), not just softmax:** `learned/structural_identity.py`
replaces the softmax proxy with an *expression-free structural identity* = **potential · specificity · confidence**
(candidate arms = pooled-HE ∩ scanMiR, with a TargetScan fallback for arms scanMiR lacks — recovers miR-137).
potential = scanMiR/TargetScan biochem × (1+log1p evidence), direction-filtered via TarBase-v9. specificity = the
family's POOLED evidence-mass fraction on g (pool the whole seed family; there is NO biochemical specialist to
floor to — every arm hits ~700 targets ~uniformly, so specificity is a pure evidence construct). confidence =
min(1, pooled ev_depth/10) enters as a **MULTIPLIER** (no evidence ⇒ no ownership; the earlier biochem-floor blend
is retired). `lost_specialist(gene)` then adds the loss layer (baseline-active GTEx/NAT → tumour-silenced) and the
**bidirectional methylation mechanism leg** (§8.4, `learned/methylation.py`): promoter hyper-methylation → repression
LOST, hypo → repression GAINED — "confirm the loss separately" is now a call, not a TODO.

---
*Below: the precursor (softmax) origin of the concepts — retained for the conceptual ladder + the MH-80
evidence, but the softmax share itself is superseded by §0.*

## 1. The attribution formula (current, as implemented)

Per-(gene g, arm m, sample s), the structural-attribution **contribution** is

```
c(m,g,s) = softmax_{m'∈R(g)}[ δ(m,s) ]_m  ×  w(m,g)

  δ(m,s) = log2RPM(m,s) − median_cohort( log2RPM(m) )            # the softmax LOGIT
  w(m,g) = log_ev(m,g) / log1p( Σ_g' log_ev(m,g') )             # specificity / promiscuity discount
```

- The softmax is **gene-local** (only across R(g), the gene's regulators) and **per sample**.
- The figure/lens quantity is the within-gene normalization
  `share(m,g) = mean_s|c(m,g,s)| / Σ_{m'} mean_s|c(m',g,s)|` (≡ `global_abs_share`, expr_mode="softmax",
  target_norm="evidence_mass").
- Engine refs: `expression_multiplier(·,"softmax")` returns δ; `_softmax_rows` is the gene-local softmax;
  `_edge_weights`/`mirna_mass_denominators` build w; `cohort_median_expression` is the per-arm median.

## 2. What question this asks

**The logits are NOT abundance — they are each arm's log-RPM deviation from its OWN cohort median.**
Consequence: per-arm centering **removes absolute abundance**. A constitutively-loud arm (miR-21) sits at
its own median → δ≈0; a near-silent specialist also sits at its own median → δ≈0; they get *similar*
softmax shares. The competition is over *who is transiently elevated relative to their own baseline*, not
who is loud. `w` then discounts promiscuous hubs (a gene is a small fraction of a hub's targeting).

So the attribution asks: **"Among gene g's regulators, which arm is (i) a specialist for g and (ii)
comparatively elevated vs its own baseline — i.e., whose IDENTITY is most attributable to this gene?"**
It is an **identity/specificity** question, **not a force question**.

## 3. Interpretation (i) vs (ii) — and why neither is "absolute pressure"

- **(i) arm relative to its OWN budget** = the weight `w(m,g)` (fraction of the arm's targeting on g) — a
  multiplier *inside* the contribution.
- **(ii) share of the GENE's budget** = the final normalized `share` — BUT it is a **relative** share with
  **absolute abundance removed** (median-centering). It is **not** "how much absolute pressure out of the
  budget."
- True (ii)-absolute requires restoring magnitude (×logrpm / ×z) — see §6.
- **Evidence of the identity≠magnitude gap:** miR-137→NCOA3 has abundance share **0.00** but attribution
  share **0.34**; miR-370→MAP3K8 abundance **0.00**, attribution **0.38**. A near-silent miRNA cannot be
  exerting the most molecular force — "#1 attributed" is an IDENTITY statement.

## 4. Promiscuity placement — WEIGHT (current) vs LOGIT (proposed) vs BOTH

The current weight is a **division** by the mass denominator `w = log_ev/D(m)`. The faithful logit form of
"divide by D" is therefore **subtract LOG D** (not raw D — subtracting raw D = `÷exp(λD)`, exponential and
far too harsh):

```
δ'(m,s) = δ(m,s) + log_ev(m,g) − λ·log D(m)      ⟺      share ∝ exp(δ)·ev / D(m)^λ
```

- **In the weight (current):** multiplicative, post-softmax — does not reshape the competition, only scales
  slices.
- **In the logit (proposed):** the softmax normalizer responds to promiscuity → **competitive,
  redistributive** (moves share hub→specialist). More faithful (a promiscuous miRNA divides its molecules →
  lower per-target occupancy = a competition effect). Hook: `_softmax_gene_logits` via the `log_ev_arms`
  slot (`softmax_z_logrpm_spec`/`_enc`).

**KEY NUANCE (resolves the "weight vs logit" question precisely):**
- For the **renormalized attribution SHARE** (`global_abs_share`, the figure lens): the per-sample softmax
  normalizer `Z=Σexp(δ)` is arm-constant and **cancels** on renormalization, so `softmax(δ)·(ev/D)`
  renormalized **≡** `softmax(δ+log ev−log D)`. → **weight-division and log-D-in-logit COINCIDE; no-op for
  the share.**
- For the **aggregate pressure** `P=Σ share·magnitude` (what couples): they **DIFFER** — weight form's
  normalizer sees only δ (promiscuity inert in the competition); logit form's normalizer `Z'=Σexp(δ)·ev/D`
  **responds to promiscuity**, redistributing the per-sample share before magnitude. → the difference is
  real only for the coupled aggregate.
- → **grid axis: promiscuity placement ∈ {none, weight, logit, both}** × source {degree, evidence-mass} × λ
  (grade by aggregate coupling, since the share is invariant).

## 4b. Reference is a TWO-SLOT axis — logit-reference AND magnitude-reference (separable)

Reference (cohort / NAT / healthy) can enter in two distinct places, meaning different things:

- **Logit reference** ∈ {cohort, NAT, healthy}: redefines *"elevated vs what"* = the **identity/competition**
  baseline. `δ_ref(m,s) = log2RPM(m,s) − ref-median(m)` (current = cohort_median; swap for NAT/GTEx median).
- **Magnitude reference** ∈ {cohort-z, NAT-dev, healthy-dev, raw-logrpm}: the **force** baseline ("pressure
  above normal"). Already implemented: `_dev_health_matrix`, `softmax_devhealthy*` modes.

They are **independent → a 2-slot grid.**

**Biological tension (interacts with §7 identity/magnitude):**
- **NAT/healthy in the LOGIT** surfaces tumor-**induced (gain)** arms (oncomiR high vs normal → big +δ →
  wins) and **down-weights silenced arms** (methylation-silenced miR-137 is *below* normal → δ_healthy<0 →
  loses share). So NAT-in-logit answers a **gain** question and **conflates identity with up-direction** —
  it would HIDE the silenced-specialist quadrant (miR-137/miR-370) we care about.
- **Cohort-median logit + specificity weight** keeps a **direction-neutral identity** (miR-137 at its own
  median, δ≈0, surfaces via specificity) → CAN surface silenced specialists.
- Empirically **MH-79's NAT-referenced share won at the aggregate (+12.6%)** ⇒ for **coupling/force**,
  pathological-engagement reference helps; for **faithful identity** (loss-direction included), cohort +
  specificity is right. **Hypothesis:** reference belongs in the **magnitude** slot (engagement vs normal)
  while the **identity** slot stays reference-neutral — to be gridded + interpreted, not settled.
- Open sub-question: signed vs |δ| in the logit (a down-vs-normal arm exerts *less* pressure, so should it
  lose share? identity vs engagement-direction).

## 5. Per-sample → aggregate → coupling (the implication chain)

- Bare softmax shares **sum to 1** ⇒ `Σ_m share(m,g,s) ≡ 1` per sample ⇒ the bare aggregate is
  **degenerate** (confirmed: share-only gene predictor ≈ dead, ~143 neg-sig).
- A meaningful, sample-varying aggregate needs a **non-constant magnitude**: `P(g,s)=Σ_m share·z·logrpm`.
- **Aggregation-lemma flip (MH-79 / §5a):** at the edge, per-arm-monotone transforms can't change
  rank-coupling; at the **Σ**, monotonicity breaks, so *how the share distributes the magnitude across
  arms* genuinely changes the aggregate's sample-trajectory → changes coupling. The share = the
  *combination rule*, not a force.
- Coupling test = "does the gene fall in samples where its regulators' share-weighted pressure is high?"
- **Caveat:** aggregating per-sample fractions **blends identity (which) and magnitude (how much)** →
  motivates separating the two axes (§7).

## 6. The construction ladder (interpretations)

| Construction | = | Interpretation |
|---|---|---|
| `softmax` | identity only | *who* is the specific regulator; sums to 1 → no aggregate force |
| `softmax_logrpm` | identity × absolute level | attributed **absolute** pressure; re-weights toward abundant arms → drifts to `abundance_sum` |
| `softmax_z_logrpm` (spine) | identity × level × responsiveness | **realized, responsive** attributed pressure; high when abundant AND elevated vs own norm |
| `abundance_sum` | pure magnitude | no identity weighting |

**Spectrum: pure identity (softmax) ↔ pure magnitude (abundance_sum); constructions interpolate.**
Prediction to test: adding magnitude should **drift the surfaced arm from specialist → abundant** (already
seen: CASP8 realized top reverts to the abundant miR-106b).

## 7. Identity vs Magnitude — the synthesizing frame

- **Identity axis** = softmax share (who; ⊥ coupling; abundance removed).
- **Magnitude axis** = abundance / logrpm / z (how much; tracks coupling).
- **Identity–magnitude plane** per (gene, arm): x = identity share, y = realized-pressure (or abundance) share.
  - **high identity, low magnitude** = *silenced specialist* → miR-137→NCOA3, miR-370→MAP3K8 (the
    methylation-loss signature; the specific TSG-miRNA that's been switched off).
  - **low identity, high magnitude** = *loud promiscuous hub* → miR-21, miR-17 (buffered / non-specific).
  - **on-diagonal** = loud specialist (identity = magnitude).
- The decoupling cases live off-diagonal; this plane is the conceptual centerpiece for "does pressure
  matter" and for reading the MH-80 reorderings.

## 8. PARKED implementation TODOs (handle later — do NOT auto-run)

1. **4th "realized-pressure" lens** (`softmax_z_logrpm` share) in the budget figures, side-by-side with the
   identity (`softmax`) lens, so the identity↔magnitude decoupling is visible per gene.
2. **Identity–magnitude plane figure** (x = identity share, y = realized/abundance share); label the four
   quadrants; mark miR-137/miR-370 (high-id/low-mag) and miR-21/miR-17 (low-id/high-mag).
3. **Edge-level partial-ρ** for miR-137→NCOA3 and miR-370→MAP3K8 (confirm directional coupling at the
   *edge*, not just the gene aggregate).
4. ✅ **DONE (2026-07-06) — miRNA-promoter-methylation loss mechanism** (`learned/methylation.py`). The
   silencing-mechanism test is now a DIRECT probe-over-promoter gate (the cCRE hop is retired: the pipeline's
   gene-centric probe annotation has no coverage below chr1:38 Mb → misses miRNA loci, which is *why* it fell
   back to cCRE). Uses the full sesame HM450 manifest (`PATHS.methylation_probe_reference`, 485,577 probes
   hg38); promoter window = strand-aware [TSS−2000,TSS+500] around the pri-miRNA hairpin TSS; gate =
   Δβ(tumour−normal) ≥ 0.15 ∧ tumour β ≥ 0.20 over the overlapping CpG probes (802 tumour vs 99 normal raw
   betas, cached probe×sample matrix). **Validated**: canonical CpG-silenced TSG-miRNAs mean Δβ **+0.061** vs
   oncomiRs **−0.070** (miR-124 +0.193, miR-9 +0.180, miR-129 +0.176; miR-21 −0.176, miR-155 −0.191). Wired
   into `structural_identity.lost_specialist` as `delta_beta`/`meth_confirmed` (priced only for LOST calls).
   It correctly *discriminates mechanism*: RB1←miR-335 (Δβ +0.011) & TGFBR2←miR-204/211 (Δβ −0.156) are real
   losses but NOT promoter-methylation-driven. NOT YET: within-tumour β↔arm-expression ρ (needs UUID→barcode
   map); genome-wide lost-specialist discovery (now unblocked by the cache).
5. **Promiscuity-placement grid axis** {none / weight / logit / both} × {degree, evidence-mass} × λ (penalty
   = subtract **log D**, ≡ ÷D^λ; grade by **aggregate coupling** since the share is placement-invariant);
   reuse `_softmax_gene_logits` / the spec/enc path.
5b. **Reference 2-slot grid** — logit-reference {cohort/NAT/healthy} × magnitude-reference {cohort-z/NAT-dev/
   healthy-dev/raw}, independently. Test the hypothesis that reference belongs in magnitude (engagement) and
   identity stays reference-neutral; include the signed-vs-|δ| logit question. Reuse `_dev_health_matrix` /
   `softmax_devhealthy*`; logit-reference = swap `cohort_median_expression` for NAT/GTEx median.
6. **Quantify the identity→magnitude drift**: compare surfaced arm under softmax vs softmax_logrpm vs
   softmax_z_logrpm on the canonical cases (does adding magnitude revert to the abundant arm?).
7. **Curate the 75 robust canonical cases** — "canonical" is a noisy automated proxy (PMID/evidence don't
   capture textbook-ness); separate curated textbook edges from hypothesis-grade.
9. **Per-gene whitelist saving** for the cartesian (~15–20 key configs, ≈28k rows) — recover the granularity
   the per-config summary dropped (the 91M-row OOM tradeoff).
10. **Promiscuity metric overhaul** — ✅ DONE (`genomewide_promiscuity.py`, HE-only genome-wide + breast-expr;
   rank-faithful r=0.92, 15% top-arm flips, headline cases robust). Remaining: add `promisc_he_expr` to the
   cartesian pdefs grid + an evidence-mass (not just degree) genome-wide variant.
11. **Unify the two share builders** (`competition_share_map` = no evidence weight, used by the cartesian; vs
   `compute_gene_pressure_contributions` evidence-weighted, used by the attribution figures) — or always state
   which. Add the evidence-WEIGHT as an explicit cartesian axis + λ-sweep for the logit penalty.
8. Decide whether any of the above is a NEW axis vs already covered by the running `apm_cartesian_overnight`
   grid (cartesian sweeps reference/function/promiscuity-def/magnitude/aggregate/contrib — but promiscuity
   **placement (logit vs weight)** and the identity↔magnitude *read-out* are likely new).

## 10. The grid map — three questions, 3-layer pipeline, archetype meanings

### 10a. There are THREE attributive questions (not two) — reference enters at a different layer for each

| Question | asks | per-miRNA output | reference enters at |
|---|---|---|---|
| **Q-IDENTITY** | *whose* gene is this (direction/state-neutral) | specificity/identity share | **none** (cohort or symmetric) |
| **Q-LEVEL** | who is pressing *now* (tumor) | realized pressure share | **none** (absolute tumor level) |
| **Q-GAIN** | who drove the *change* vs normal | gain attribution (signed) | **magnitude-difference** layer |

The silenced specialist (miR-137) reads correctly under all three: **high** Q-IDENTITY (the specialist),
**small** Q-LEVEL (barely expressed now), **negative** Q-GAIN (lost pressure vs normal = the methylation
event). They conflict only if a Q-GAIN construction is mislabeled as Q-IDENTITY. **Clean Q-GAIN attribution
= build realized `c(m,g,s)=share×magnitude` FIRST, then difference across reference `Δc=c_tumor−c_normal`**
— reference at the magnitude layer, identity stays reference-neutral. Putting reference *inside the logit*
conflates "who is induced" (within-competition reweight) with "how much more pressure" (a magnitude Δ).

### 10b. The construction = a 3-layer pipeline

```
LAYER 1 IDENTITY (share, within-sample, Σ=1): comp fn {softmax/temp/sparsemax/linear/mass-action}
   × logit content {δ(ref) + [−λ·log D promisc?] + [specificity bonus?]} × logit ref {cohort/NAT/healthy}
LAYER 2 MAGNITUDE (per-arm force): term {none/logrpm/z/z·logrpm/dev} × mag ref {cohort-z/NAT-dev/healthy-dev/raw}
   → contribution c(m,g,s) = share × magnitude
LAYER 3 QUESTION (read-off): reference-diff {level(tumor)/gain(tumor−normal)} × aggregate {sum/mean}
   × sign {signed/pos/abs}   → P(g,s) (couple) OR per-arm attribution vector
```
Cross-cutting: **within-sample vs cohort** (L1 softmax = within-patient competition; z/δ + coupling =
across-cohort) and **per-miRNA-attributive vs aggregate-coupling** (read the per-arm vector vs Σ-then-correlate).

### 10c. Archetype map — hypothesized meaning per construction

| Construction | Layers | Per-miRNA meaning | Aggregate / coupling | sample↔cohort |
|---|---|---|---|---|
| bare softmax (cohort) | id only | identity: specific+locally-elevated regulator | **degenerate** (Σ=1), a label not a force | within-sample, cohort logit |
| softmax, NAT-logit | id only | induction-identity: up vs normal | degenerate; gain-biased, **hides silenced** | within-sample, normal logit |
| abundance_sum | mag only | raw abundance share | total abundance → couples (**baseline**) | pure cohort magnitude |
| softmax·logrpm | id×mag | absolute attributed pressure | attributed absolute load → couples; **drifts to abundant** | share × cohort level |
| softmax·z | id×mag | dynamic engagement: who *moves* | total responsiveness | share × cohort deviation |
| softmax·z·logrpm (spine) | id×mag | realized responsive pressure | realized total → **best coupling** | share × cohort z·level |
| softmax·dev_NAT | id×mag, NAT-mag | pressure *above normal* per arm | total **acquired** pressure → tumor-specific repression | reference at magnitude (clean Q-GAIN-level) |
| Δ(tumor−NAT) realized | id×mag, L3 gain | **who drove the change** (loss=neg) | pressure **gain** | two-state differential |
| +promiscuity-in-logit (any) | L1 mod | share **identical** (renorm) | **redistributed** → diff coupling | no-op per-miRNA, real aggregate |
| comp-fn {mass-action/sparsemax/temp} | L1 mod | how winner-take-all the share is | concentration of the aggregate | within-sample competition shape |

### 10d. Hypotheses to test on the grid (mapped to layers)
1. Magnitude drifts surfaced arm **specialist→abundant** (softmax→softmax_logrpm→softmax_z_logrpm→abundance_sum
   is a continuum) — quantify.
2. Promiscuity placement = **no-op per-miRNA, real lever for coupling** (§4) → grade by aggregate.
3. Reference belongs in **magnitude** (engagement/gain), identity **reference-neutral** → NAT helped MH-79
   coupling as Q-GAIN-level, not by reshaping identity.
4. Off-diagonal of the **identity–magnitude plane** = the interesting genes (silenced specialists vs loud hubs).
5. Open (point 3): identity reference-neutral vs **symmetric |δ|** (engagement-neutral); promiscuity = **degree**
   (spreads attention) vs **evidence-mass** (spreads *confident* attention).
6. mass-action keeps a free-capacity term (≈ saturating abundance at single-regulator genes) → expected closest
   non-identity contender; sparsemax = hard arm selection; temperature tunes concentration.

## 11. Corrections & refinements from the cartesian results (2026-06-28)

- **CORRECTION — promiscuity-in-logit WAS tested.** The cartesian promiscuity axis {degree, evidence_mass}
  enters as `base_logit − disc` (disc = log1p(#targets) / evidence-mass = the log-D form), division in
  linear space for linear/massaction. It **HURT coupling** (none 451 > degree/evidence 435), at **λ=1 only**.
  So my §4 "proposed (untested)" framing was wrong for the logit; what's still untested = λ-sweep + the
  evidence-**weight** inside the shares + a clean weight-vs-logit head-to-head.
- **The cartesian shares carry NO evidence-specificity weight** (`competition_share_map` ignores evidence) —
  pure abundance-deviation softmax. So the **coupling-optimal** construction (`nat|z`, no promiscuity, no
  specificity) is NOT the **identity-faithful** attribution (the MH-80 figures use evidence-weighted
  `compute_gene_pressure_contributions`). **Coupling-optimal ≠ identity-faithful** — unify the two share
  builders or at least always state which is in use.
- **Promiscuity metric overhaul (2b genome-wide + 2c breast-expression) — DONE 2026-06-28**
  (`genomewide_promiscuity.py`, cached `output/genomewide_promiscuity/he_targetome_promiscuity.tsv`;
  wired as `promisc_he_expr` in `_ref_baselines`):
  - **2b genome-wide** + **2c breast-expr filter (TCGA-tumor median log2>1, NOT GTEx-only)**, **HE edges
    only** — `Support Type == 'Functional MTI'` (the strong set; the 3.98M `Functional MTI (Weak)` records
    are excluded — including them gave a spurious 84× inflation, rank-uncorrelated, that falsely made
    miR-137 look like a 930-target hub).
  - **Result: HE genome-wide degree is rank-faithful to the hallmark-scoped one (Spearman 0.92), ~1.5× scale**;
    expression filter trims ~10%. Swapping the denominator flips the surfaced arm for **only 14/94 (15%)** of
    canonical genes, and **all headline cases are stable** (ABCA1→miR-33b, NCOA3→miR-137, CASP8→miR-376a,
    ESR1→miR-4728, SPRY2→miR-24). So MH-80's attribution is **robust** to the corrected metric.
  - Sharpens the **identity/specificity** axis (the attribution), **not** coupling (cartesian: promiscuity
    hurts coupling). The available `promisc_he_expr` baseline can be added to the cartesian pdefs grid.
- **Per-gene grid results were not saved** (per-config summary only — the 91M-row OOM). Fix: persist per-gene
  rows for a **whitelist** of ~15–20 key configs (winners + baselines + archetypes), ≈28k rows.

## 12. Identity–magnitude diagnostic + edge-level reality check (2026-06-28, executed)

- **Identity–magnitude plane + drift** (`pressure_attribution_validation.py --identity-magnitude`;
  `identity_magnitude_plane.tsv`, `magnitude_drift.tsv`, `fig4_identity_magnitude.png`): the two axes
  **decouple** — **identity-top = the surfaced specialist arm in 94%** of canonical genes, **magnitude-top =
  the abundant (masked) arm in 77%**; the surfaced arm drifts specialist→abundant along the ladder, and
  **94% of that drift happens at the identity→level (+logrpm) step** (level→magnitude +z moves only 22%) —
  i.e. the **absolute level** drags identity back to abundance (matches the cartesian: z helps, logrpm hurts).
  Hypotheses 1 & 4 confirmed.
- **EDGE-LEVEL REALITY CHECK (claim-changing) — the surfaced specialists do NOT couple; the abundant arms do.**
  Edge partial-ρ (arm vs target | CPE+HRD+target-CN [+composition]): **NCOA3** surfaced miR-137 = −0.019
  (p=0.73, NS) vs masked **miR-17-5p = −0.25 (p=1.6e-15)**; **MAP3K8** surfaced miR-370 = −0.05 (NS) vs masked
  **miR-375 = −0.26 (p<1e-15)**. miR-137 is usable in only 342/957 samples (largely silenced). **So the
  attribution surfaces the structural specialist + literature TSG-miRNA, but TCGA coupling does NOT confirm it
  as the ACTIVE regulator — the gene's net-repression is carried by the ABUNDANT arm, and the specialist is
  largely silenced (no variance to couple).** This is §5a empirically (identity ⊥ coupling) and **tempers
  MH-80**: the surfaced "canonical" arms are *identity/literature hypotheses*, not coupling-confirmed drivers;
  ABCA1→miR-33b remains the strongest because miR-33b is both specific AND expressed. The "silenced specialist"
  read (high identity, low magnitude, no coupling because switched off) is consistent — and is exactly what the
  methylation check (§8.4) tests.

## 9. Biological throughline (the directionality of the reorderings)

Every MH-80 flip runs **loud-promiscuous-hub (abundance) → low-abundance specialist (attribution)**, and
the specialist is repeatedly a **methylation-silenced tumor-suppressor miRNA** (miR-137, miR-370; cf.
miR-33→ABCA1). So the attribution axis behaves as a detector of *"which specific, often epigenetically-
silenced, regulator owns this gene's identity"*, while abundance detects *"which oncomiR is loudest."* For a
TSG-miRNA-**loss** mechanism, the attribution axis is the biologically relevant one and abundance is
actively misleading — but each case must still be closed with the **coupling** (magnitude) and the
**literature**, never the share alone (§5a).
