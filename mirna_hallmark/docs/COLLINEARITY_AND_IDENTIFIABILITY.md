# Collinearity & identifiability — the two kinds, how each is resolved

How the learned model handles the fact that miRNA abundances co-vary, so "which arm/family repressed gene g"
is often **non-identifiable** from the data alone. Companion to `LEARNED_MODEL_ESTIMATOR_MAP.md` (which cites
"family = the identifiability unit" throughout without laying out *why* or *what's left*). Status: **active**.

---

## 0. The problem

When predictors co-vary, the **aggregate** `X·M` is stable but the **individual member weight** is not:

- collinearity inflates *coefficient* variance, not *fitted-value* variance → `ρ(X·M_lasso, X·M_baggedNNLS)` within 0.03, aggregate bootstrap SD ~0.01;
- but a single lasso fit's *per-member* coefficient has bootstrap corr **0.03** across resamples.

So the gate / discovery / aggregate-coupling read the **stable** object and are fine; everything that asks
**"who owns g"** (attribution, budget, identity, driver nomination) is a collinearity problem. There are **two
biologically distinct kinds**, with different origins, different identified estimands, and different resolutions.
Conflating them is the trap.

---

## 1. Kind A — SEED-FAMILY collinearity (shared *targetome*)

**Origin.** Arms with the **same seed** (miR-200b/200c/429; the members of a TargetScan family). Same seed ⇒
*identical predicted targets* ⇒ they repress g **through the same sites**. This is the *hard* kind: confounded
not just in abundance but in **targeting** — even with independent abundance you couldn't tell which arm's
binding cut g, because it is the same binding site.

**The estimand changes.** Because they share targets, the identified object is **family→gene**, not arm→gene.
We don't fight the collinearity — we make the family the unit. Three steps:

1. **Collapse** — `families.collapse_by_family`:
   `X_fam = log2(1 + Σ_members (2^x − 1))` (pool in *linear* RPM = biologically correct total family dose,
   re-logged); `w_fam = max` member prior. One predictor; the identifiability unit **everywhere** (gate,
   discovery, canonical M, Shapley all run on families).

2. **Apportion to arms** — `states.budget_shift` (§5a):
   `contribution(arm) = M_fam · X_fam · (2^a / Σ 2^a)` — split the family's realized pressure by each arm's
   **linear-RPM share**, so `Σ_arms = M_fam·X_fam` (no family-size inflation). An **abundance nomination**
   ("which arm carries the dose"), NOT a resolved coefficient.

3. **Resolve the driver** — `families.resolution_report` (target-driven):
   residualise the member arm **and the target Y** on (other members + C); a member is a driver only if its
   **own** conditioned partial-Spearman with g survives (`ρ < −0.1`). Output `single-driver / co-drivers /
   family-only`. A truly collinear member has ~no residual variance ⇒ auto `family-only` (no threshold needed).

**Status: fully built.** Honest floor: when members are genuinely collinear, resolution returns `family-only`
and the abundance apportionment is the only split — flagged as a nomination, not a coefficient.

---

## 2. Kind B — GENOMIC-CLUSTER collinearity (shared *transcription*, distinct targetomes)

**Origin.** miRNAs co-transcribed from one locus (miR-183/96/182; miR-17~92) or co-regulated across loci (the
two miR-200 clusters). **Different seeds ⇒ different targetomes**, but **correlated abundance** (shared
transcription). The *nuisance* kind: the abundance collinearity is real, the biology is arm-specific.

**Why it is fundamentally easier than A.** They **do not share targets**. For gene g, typically only *one*
cluster-mate has a seed match in g's 3′UTR — so **sequence disambiguates what abundance cannot**. The
co-transcribed non-targeter is riding along, not regulating g.

**Two sub-cases:**
- **B1 — same polycistron** (183/96/182 = one primary transcript, one locus, one CN): collinear in abundance
  *and* CN ⇒ only the target side + conditioned expression separate them.
- **B2 — separate loci, co-expressed** (two miR-200 clusters, different CN): the **CN instrument** separates
  them — each locus's copy number independently perturbs *its* cluster.

**Current handling:**
- **General identity/budget path** — emergent only: adaptive lasso *selects* one (unstable); Shapley *splits*
  the credit (a "shared" identity, not a single owner); `parallel_view.within_gene_ranks` now *names* the
  collinear pair (`|ρ_abund|≥0.7`). No sequence disambiguation.
- **CN-instrument path** — `instrument.cluster_attribution` **already resolves Kind B**: finds co-located
  co-targeters (`pleiotropy`), then either **cluster-clean** (no co-located co-targeter ⇒ arm-specific) or
  **attributes within the cluster** by conditioned partial anti-corr (residualise arm on co-located
  co-targeters + C; unique survivor = carrier), CN supplying exogeneity. Same §4b logic, keyed on genomic
  co-location rather than seed.

**The problem, precisely:** the cluster resolution lives **only inside the CN-instrument causal test**, scoped
to genomically co-located co-targeters. The **general** attribution/identity/discovery path for a gene never
applies it — it splits or selects, and never says "of these two co-expressed clusters, only miR-141 has a site
in g, so it is the owner; miR-200a is co-transcribed noise."

---

## 3. How to resolve Kind B (in order of cost)

1. **Sequence-first prune** (cheapest, do first): among **abundance-collinear** regulator pairs (`|ρ|≥0.7`),
   if one has biochemical potential on g (scanMiR/TargetScan site) and the other ~0, the site-less one is the
   co-transcribed rider ⇒ drop/down-weight it *before* the identity/budget split. Turns Kind-B collinearity
   into either a single owner or a genuine multi-regulator set. **→ built as `families.cluster_resolution`.**
2. **Generalized conditioned resolution**: run `resolution_report`'s conditioned partial anti-corr across
   **genomic-cluster-mates** (different seed families), not just seed-family-mates — `cluster_attribution`'s
   logic lifted out of the CN path into observational attribution. Survivor = driver; collinear rider drops.
   **→ built into `families.cluster_resolution` (the `coupling` column).**
3. **CN-instrument disambiguation for B2** — **BUILT** (2026-07-07) into `cluster_resolution`: per-family locus
   CN as independent instruments (different loci ⇒ CN separable even where expression is collinear). Joint
   reduced-form `g ~ Σ CN_locus + C` → per-family CAUSAL verdict `cn_causal` ∈ {**CN-repressor** (first-stage
   fs≥0.15 ∧ exogenous reduced-form ≤ −0.1), **CN-null** (adjudicable, no effect), **CN-blind** (weak/shared
   instrument)}. It is NOT a forced single owner — CN speaks only where its locus instrument is strong. Guarded
   against B1 (same-polycistron: if CN_f loses >60% of its variance to a cluster-mate's CN, same locus ⇒
   CN-blind). **Validated:** ZEB1 miR-200bc/429 adjudicable (mild/null); **PTEN → CN-causal repressor miR-19**
   (the canonical miR-17~92→PTEN edge, Olive 2009) where the *observational* coupling picked miR-18 — the
   exogenous layer can **correct a confounded observational call**. **Caveats:** sparse (miRNA-locus CN
   instruments are mostly weak → most clusters CN-blind); paralog-spanning seed families (chr13/chrX
   17~92/106a~363) give ambiguous locus CN. Remaining gap: FDR over the cluster reduced-forms.
4. **Hierarchical partial-pooling** (`families.py:11`, Phase 3b): instead of collapse (= infinite-shrinkage
   limit), fit a family/cluster weight + a **shrunk per-member δ with posterior width** — every member a
   shrunk-but-nonzero weight *and* an uncertainty. Turns "unidentifiable" into a quantified credible interval.
   Subsumes A and B; the honest endpoint. **Parked.** NOTE — do not confuse with the *cross-gene* pooling that
   IS built (`hierarchical.py`: `β(m,g) ~ N(μ_m, τ²)`, sharing arm m's effect ACROSS the genes it targets),
   which was tested and **did not help** (Δρ −0.002; effects are gene-specific). That is an **orthogonal pooling
   axis** — the cross-gene null does *not* predict the within-family/cluster δ-pooling (a per-gene, member-level
   shrinkage) would fail. The collinearity-relevant one is the within-family δ-pooling, still unbuilt.

**How often B fires (genome-wide, 2026-07-07):** only **60/1571 pooled-HE genes (4%)** carry a cross-seed
collinear cluster; of 154 families-in-clusters, **119 designed-not-coupling, 30 resolved owner, 3 rider,
2 couples-no-site** (`output/learned/cluster_resolution_scan.tsv`). So in the curated *attribution* set Kind-B is
a minority, resolved by conditioned coupling (step 2); the sequence-prune (step 1) barely fires (3 riders) — it
earns its keep in the **discovery/orphan** set, where site-less predicted edges appear. The 2 `couples-no-site`
are the reference-blind (SNV/editing/APA) candidates.

---

## 4. Unifying view

| | Kind A (seed family) | Kind B (genomic cluster) |
|---|---|---|
| shared | **targets** (same seed) | **transcription** (same locus/regulation) |
| identified estimand | family→gene (collapse) | arm→gene per-mate (distinct targets) |
| separator | conditioned coupling only (targets identical) | **sequence** (which mate targets g) + conditioned coupling + **CN** (B2) |
| built | collapse + apportion + resolve ✅ | `cluster_resolution`: sequence-prune + conditioned coupling + **CN-overlay (B2 causal)** ✅ |
| parked | hierarchical δ-pooling | CN-overlay FDR; hierarchical δ-pooling |

**One line:** Kind A you resolve by *collapsing* (the family **is** the answer); Kind B you resolve by the
*target side* (sequence tells you which co-transcribed mate is real), with CN as the exogenous separator when
the mates sit on different loci. The **hierarchical pool** is the single principled method covering both — the
real "resolve to members" endpoint, currently parked.

---

## 5. Code map

| Object | Module.function |
|--------|-----------------|
| seed-family map | `families.family_of` |
| collapse to family predictor | `families.collapse_by_family` |
| arm apportionment (abundance) | `states.budget_shift` |
| Kind-A driver resolution (conditioned coupling) | `families.resolution_report` |
| **Kind-B cluster resolution (sequence prune + conditioned coupling)** | `families.cluster_resolution` |
| Kind-B CN causal disambiguation (co-located co-targeters) | `instrument.cluster_attribution`, `instrument.pleiotropy` |
| collinear-pair naming (diagnostic) | `parallel_view.within_gene_ranks` |
| collinearity-fair identity split | `attribution.shapley_identity` |
| hierarchical δ-pooling (parked) | `hierarchical.py` (Phase 3b) |
