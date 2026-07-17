# Learned-τ² model + discovery — SYNTHESIS (the whole architecture, one doc)

Capstone summary of the learned-τ² Bayesian regulatory model and the discovery/validation apparatus built
on it (MH-82 → MH-84, 2026-07). Companion to the estimator-level docs: `LEARNED_MODEL_RATIONALE.md`
(the *why* per §), `LEARNED_MODEL_METHODS.md` (formulas), `LEARNED_MODEL_ESTIMATOR_MAP.md` (which
estimator), `LEARNED_MODEL_VALIDATION.md` (numbers). This doc is the **"how the pieces fit into one
system, and what it produces."**

---

## 0. Thesis in one paragraph

A **single dense learned-τ² non-negative Bayesian regression**, fit per gene over an augmented regulator
support `[HE ∪ discovered orphans]`, serves **four jobs** from one posterior — coupling, attribution,
identifiability, and (via an inclusion-mode variant) discovery. It is **non-circular** (priors from
curation + sequence + biochemistry only, never TCGA X/Y), **composition-adjusted** (CIBERSORTx
deconvolution the canonical way), **credibility-gated** (validated on independent CPTAC protein — ⚠ **see §6d:
under the MH-127 decoy control the out-of-cohort validation that survives is CPTAC *mRNA*, not protein; and the
in-cohort "beats abundance" gate is a FITTING result, not a curation result**), and
**FDR-controlled**. The discovery half nominates edges from two candidacy sources (TargetScan + scanMiR
K_D), fuses them into a consensus, and validates each at three independent layers (mRNA, protein,
literature) — yielding **268 triple-validated novel edges**. ⚠ **QUALIFIED 2026-07-17 (MH-107/114): `dossier.tier3_protein` ran with NO composition block. ~27% of the 268 SIGN-FLIP under adjustment and mean |ρ| drops −60% (−0.190 → −0.075). The code is fixed but `output/learned/dossier.log` is still dated 2026-07-09 = PRE-FIX ⇒ the count is stale. Distinct from the retracted MH-38 '108 orphans' claim, but the same confound. Re-run before citing.**

---

## 1. The estimator

Per gene `g`, on confounder-residualised, z-scored seed-family predictors:

    r = X·β + ε ,  ε ~ N(0, σ²)                       # r = −resid(Y|C): repression → positive
    β_m = z_m · θ_m                                    # coefficient = inclusion × magnitude
    θ_m ~ N⁺(0, ν²)                                    # half-normal slab — β ≥ 0 ("miRNAs repress")
    z_m ~ Bernoulli(π_m)                               # inclusion indicator
    ν², σ² ~ InvGamma                                  # ν² (=τ²) LEARNED, not guessed

One knob, `π`, switches the two modes:
- **π ≡ 1 (DENSE)** → learned-τ² ridge. Every candidate in; ν² sets shrinkage. **The core-jobs fit.**
- **evidence/uniform/learned-base-rate π (INCLUSION)** → spike-and-slab. Genuine selection (PIP). **Discovery.**

**Why it wins (genome-wide):** the dense learned-τ² ridge beats the shipped adaptive lasso on OOF coupling —
mean ρ **−0.168 vs −0.152**, wins 58%, **Wilcoxon p = 9×10⁻¹⁶** — because at n≫p the lasso's sparse
selection over-shrinks the aggregate, and dense learned-τ² keeps the diffuse many-arm signal (concentrated
in many-regulator genes). CV-ridge ≈ lasso, so it is the *learned* τ², not L2 per se. Every alternative was
tried (spike-slab, Bayesian lasso, EB-shrink-toward-prior, the full 2×2 slab×selection grid); they collapse
to one axis — **dense wins coupling, inclusion wins discovery; slab shape and CV-vs-EB are second-order.**

---

## 2. Four jobs, one fit

| job | readout of the *same* dense posterior | replaces | verdict (genome-wide) |
|---|---|---|---|
| **Coupling** | aggregate `X·E[β]` | adaptive lasso | **better** (p=9e-16); CPTAC-protein-validated |
| **Attribution** | posterior **mean ± sd** | bagged NNLS | **equal + unified** (agree 0.80, repro 0.97>0.98·NNLS) |
| **Identifiability** | **`\|z\| = \|mean/sd\| > 2`** | §9 conditioned-partial | **equal** (prec 0.86, rec 0.89; = full-conditional, shrunk) |
| **Discovery** | **PIP** (inclusion mode) | lasso + permutation scan | complementary; see §4 |

Key facts that make this *one* fit, not four:
- `|z|>2` **is** the full-conditional partial-Spearman, shrinkage-stabilised (Spearman 0.88), read free from
  the model — so identifiability needs no separate estimator.
- Attribution reads the same posterior mean; the calibrated sd comes for free (bagging only approximates it).
- **Discovery is the exception** — it needs the *inclusion* mode, because it must nominate a *representative*
  of a collinear orphan group (model-averaging), which the dense `|z|` over-conditions away.

---

## 3. Composition / deconvolution — the canonical way

**Yes, cell-composition is adjusted the standard CIBERSORTx way, throughout.** The confounder block `C` can
carry the **8 non-malignant CIBERSORTx fractions** — `CAFs, T-cells, Myeloid, B-cells, Endothelial, PVL,
NormalEpi, Plasmablasts` — plus malignant-compartment **proliferation** (E2F/G2M metagene residualised on
the Cancer-Epithelial fraction). Cancer-Epithelial itself is deliberately excluded (the target lives there —
conditioning on it over-controls).

Where it is applied:
- **Attribution** forces `deconv=True` → cell-intrinsic weights.
- **Discovery** applies it as the **robustness gate**: `robust = the edge's coupling survives adding the
  CIBERSORTx fractions to C`. **Every edge in the consensus tables is composition-robust by construction.**
- **Retrofit** uses `deconv=True`.
- **Gold-standard confirmation** (MH-70/71): 96% of edges stay negative under **BayesPrism** deconvolution —
  the result is not a CIBERSORTx artifact (BayesPrism 96% / CIBERSORTx-S 94% / metagene 93% converge, same
  miR-29→ECM stromal exception).
- **Dossier** (§5): the per-edge `realized_coupling` is reported both raw and **deconv-adjusted**
  (`realized_coupling_deconv`, gold sets). Result: of the 676 gold edges that couple raw (ρ<−0.1),
  **658 (97%) still couple after CIBERSORTx composition adjustment**; raw↔deconv coupling correlate **0.96**
  (mean attenuates −0.228 → −0.175, ~23%, i.e. a small composition-driven component removed, the signal
  intact). So the gold discoveries are **composition-robust at the per-edge level**, not just via the gate.

---

## 4. Discovery pipeline

**Candidacy (two sources).** Orphans = predicted (miRNA, gene) edges beyond the curated HE set:
- **TargetScan** (default) — context++ conserved sites. Candidacy ⊥ the scanMiR gate → independent.
- **scanMiR K_D** — top-80 strongest predicted duplexes/gene. Biophysical, reaches non-canonical seeds;
  but candidacy ≈ the scanMiR gate (same source) → **less independent** ⇒ complementary, not primary.

**The grid (per candidacy source).** Spike-and-slab, gated (scanMiR duplex + deconv robustness + permutation
FDR), across the inclusion-prior modes:

| prior | role | TargetScan (robust / FDR) |
|---|---|---|
| **learned base-rate** (Beta–Bernoulli, no hand-set value) | **canonical shortlist** | 885 / **0.001** |
| **uniform** (π=0.5) | high-recall PIP | 1225 / 0.014 |
| **densez** (`\|z\|`-over-orphans) | high-recall, complementary to PIP | 1344 / 0.044 |
| evidence-π (skeptical) | **retired** (dominated) | 1115 / 0.070 |

`|z|` and PIP are **complementary** for discovery (densez adds ~430–620 beyond PIP): `|z|` = full-conditional
(drops redundant groups), PIP = model-averaged (keeps a representative). The prior in PIP is a
**recall/precision dial** — learned = clean shortlist, uniform = high-recall; not a requirement (learned lets
the data set the inclusion rate and lands conservative on its own).

**Consensus (persisted).** `discovery_consensus_{targetscan,kd,combined}.tsv` — one row/edge, boolean per
method, so **any intersection is a filter**. Combined **TS ∪ K_D = 6,744 edges**. Key tiers:
- **all-5 ultragold = 414** (marginal partial + 3 Bayes priors + full-conditional all agree)
- **triple-gold = 431** (partial ∩ learned-PIP ∩ densez-`|z|`)
- **cross-candidacy = 503** (both TargetScan-learned AND K_D-learned) — real regardless of nomination source
- Cross-validation: **187/431 (43%)** of TargetScan gold is independently re-discovered by biophysical K_D.

**Retrofit (closes the loop).** Fold discovered orphans into the support, re-fit the dense ridge over
`[HE ∪ discovered]`: **100% of learned discoveries earn `|z|>2` in the joint fit** (cross-estimator check —
selected by PIP, confirmed by `|z|`), and **~9% of curated-regulator identifiability calls flip** (the
omitted-variable bias real discoveries were causing). One augmented-support fit then serves all four jobs.

---

## 5. Per-edge dossier + 3-layer validation

`discovery_dossier.tsv` — every one of the 6,744 edges characterised: realized coupling
(partial-Spearman arm↔target | C + HE-aggregate), detectability, biochem, literature, and a `dossier_pass`
flag (couples ρ<−0.1 **and** expressed **and** duplex).

- **Overall:** 72% couple, 37% pass. TargetScan 85% couple > K_D 65% (K_D broader = noisier per-edge);
  cross-candidacy 74% pass = highest quality.
- **The consensus is validated, non-circularly:** the more independent methods agree, the higher the
  *independently-measured* coupling — **triple-gold/ultragold couple 100%, pass ~81%** (the ~19% non-pass
  fail on expression/duplex, *not* coupling).

**Tier 3 — independent CPTAC protein** (`discovery_dossier_gold.tsv`): on the gold sets (608 testable),
**66% protein-coupled** (ρ<0 with target protein in independent prospective-CPTAC patients), 59%
discordance-coupled (beyond-mRNA). **★ 268 edges are TRIPLE-VALIDATED** — couple at **mRNA** (TCGA) + **protein**
(independent CPTAC) + **literature** (PMIDs). Led by the miR-17~92 cluster on a coherent target set:
AHNAK←miR-106b (mRNA −0.54 / protein −0.64), MAPK1/TIMP2/TNS1/PDLIM5←miR-106b, PTPRM←miR-130b,
AKAP12←miR-15b, IRS2←miR-181b, LPAR1←miR-429.

---

## 6. The priors (non-circular)

`priors.py` — `π` (inclusion), `μ` (biochemical magnitude, scanMiR), `τ` (evidence-depth confidence), and a
**seed-rarity** channel (`seed_rarity.py`). All from curation + sequence + biochemistry, **never TCGA**:
- **π** — evidence-graded inclusion, monotone in the PMID-deduped ledger weight.
- **μ / τ** — slab-*scale* channels (loosen shrinkage for strong/deeply-evidenced edges); location stays at 0,
  so the prior sets ordering, the data sets magnitude.
- **seed-rarity** — the **collinearity-immune resolver**: same-seed arms are expression-collinear, so every
  expression axis is trapped in the soup; sequence seed-rarity adjudicates the collinear direction the
  likelihood cannot (e.g. within miR-17~92, miR-18a's rare seed nominates it as PTEN's specialist).

The shuffled-evidence null confirms the *dense coupling* is prior-**independent** (evidence-shuffle doesn't
change it) — so coupling can't be circular through a prior it ignores; the prior's work is in
inclusion/discovery.

### 6a. Seed-rarity — locus decided (π, not slab); gain NOT validated → term OFF by default (`rarity_bench.py` + edge audit, 2026-07-09)

Two things were tested on collinear soups {PTEN, CCND1, VEGFA} vs clean {ESR1, CDH1}, both readouts of the one
posterior, gain sweep {0.4…1.6} × loci {π, slab-gated·τ, slab-standalone}. **One conclusion holds, one is retracted.**

- **✅ Locus = π, not the slab (holds).** Rarity in the slab-scale (the channel the *dense* model A uses) is
  **inert**: dense OOF coupling ΔρOOF = **+0.000** at every gain/gene, attribution magnitudes identical to 3
  decimals — at n≫p the likelihood dominates the slab prior (a rare edge's slab doubling 1.36→2.96 moves M by
  3e-5). *Closes "should specificity also enter model A?" — it can be wired in, but the data won't let it matter.*
  Rarity's only possible leverage is the discrete **inclusion** decision.
- **❌ Gain = 1.2 (RETRACTED).** It was tuned to a **proxy** — "does the rarest-seed edge cross PIP 0.5" — which
  rarity trivially satisfies and which is **not** biological correctness. An independent audit (edge-level
  partial-Spearman coupling + arm expression) shows the lifted "specialists" are the **wrong** edges: PTEN miR-105
  **unexpressed** (RPM 0.48), CCND1 miR-99/100 **non-coupling** (ρ −0.02, q_BY 1.0), VEGFA miR-718 non-coupling
  (ρ −0.07, q_BY 0.45). Worse, across all covered edges **Spearman(rarity, edge-ρ) = +0.10 / +0.11 / −0.01** (wrong
  sign) — on PTEN the *commonest*-seed edges couple twice as hard (−0.098) as the rarest (−0.044). Seed rarity
  ≠ realized coupling; a global high-gain rarity channel promotes globally-rarest seeds, which are
  disproportionately unexpressed/low-confidence orphans with tiny (often artifactual) targetomes.

**Why the original miR-18a case worked and this didn't:** miR-18a-in-miR-17~92 is **moderate**-rarity (0.26),
expressed, and a member of a *coupling cooperative family* — there rarity legitimately tiebreaks collinear co-seed
arms the marginal likelihood can't separate. The failure mode is treating rarity as a *global* inclusion promoter,
which selects for obscurity, not specificity.

**Current state:** `rarity_gain` **defaults to 0 (rarity OFF)** in `priors.inclusion_prior`/`edge_priors`; the
plumbing is retained. Output: `output/learned/rarity_bench.tsv`.

**Concept validation of the redesign — also negative (2026-07-09).** Before building the expression-gated,
within-cooperative-family redesign, its premise was concept-tested on the real clusters (miR-17~92, miR-15/16) it
was meant to serve, against BOTH candidate targets:
- **vs edge coupling (magnitude):** within-cluster Spearman(rarity, edge-ρ) = −0.07 (PTEN) / +0.30 / +0.50 (CCND1)
  — the *commonest* seed is the strongest coupler in every cluster (miR-17 family ρ −0.47 on PTEN, rarity 0.05).
  *Caveat:* coupling is the magnitude axis and rarity is an identity axis (ATTRIBUTION_IDENTITY_VS_MAGNITUDE), so
  this target is arguably category-confused — hence the second test.
- **vs evidence-mass specificity (identity, non-circular; `structural_identity.specificity`):** Spearman(rarity,
  ev_spec) = **−0.53 (PTEN cluster)** / +0.50 (CCND1 cluster) / +0.30 / +0.11 (global) — **inconsistent**; in the
  PTEN cluster the rarest members (miR-18, miR-193) are the *least* evidence-specific.

**Verdict on rarity:** seed-targetome-*size* rarity is not a sound specificity operationalization — it tracks
obscurity (small/artifactual targetomes of low-confidence, often unexpressed arms), not dedication, and agrees with
neither coupling nor the curated evidence-mass specificity. It is retired (OFF).

**The right statistic — per-edge affinity CONCENTRATION (validated, 2026-07-09).** [⚠ SUPERSEDED — see §6c:
concentration-*share* was itself found flawed (HE-share inflates small targetomes → wrongly crowned miR-96 over
canonical miR-200 on ZEB1); the FINAL identity statistic is affinity **PERCENTILE** "is g among the family's
strongest targets", seed-family-keyed. The convergence-with-evidence-mass result below still stands.] The fix was
a shape change: specificity must be a **per-(arm,gene)** property, not a per-arm one. Defining
`seq_spec(arm,gene) = affinity(arm→gene) / Σ_{g'∈HE(arm)} affinity(arm→g')` — the share of the arm's total
predicted repression aimed at this gene, from TargetScan **context++** (which already folds in site-type +
local-context + conservation; `Predicted relative KD` is the alternative biophysical weight), restricted to the
arm's **HE target universe** (cleaner denominator than the weak genome-wide tail, aligned to the model's gene
space) — recovers the textbook specialists rarity missed: **miR-19b/19a & miR-26a → PTEN** (Olive 2009, Huse 2009),
**miR-15/16/195 → CCND1** (Cimmino 2005), **miR-200/429 → ZEB1**, with the abundant generalist miR-17 correctly
*last*. Systematically it **converges with the independent curated evidence-mass specificity** at Spearman
**+0.49…+0.69 across 6/6 testable genes** (PTEN .56, CCND1 .61, ESR1 .49, VEGFA .50, ZEB1 .69, BCL2 .61), vs rarity
+0.04…+0.30 (ZEB1 −0.10). Two non-circular specificity signals — sequence affinity-concentration and curation
evidence-concentration — agreeing at ρ~0.6 is convergent validation of a real identity axis. *Caveat:* the two
share the **curated HE universe** (denominator), so the agreement is in the numerators (sequence vs
evidence-depth), which are genuinely independent; a genome-denominator version (needed for uncurated orphans) is
purer-sequence but noisier.

**Built + inclusion-π role tested (`seq_specificity.py`, `seq_spec_bench.py`, 2026-07-09).** The statistic is a
module (`seq_spec(gene, candidates, universe='HE'|'genome')`) and is wired into `priors.inclusion_prior`/
`edge_priors` as an **expression-gated** π channel (`spec_gain`, default 0 = off; gate = soft detectability floor
log2(11)). Retrying the goal rarity failed: gated seq-spec-π **works safely** — it lifts genuine expressed+coupling
specialists (PTEN miR-25/92 PIP 0.52→0.86, edge-ρ −0.35) and the **expression gate suppresses the unexpressed
high-seq-spec arms that rarity promoted** (miR-494 gate 0.17, miR-216a 0.10, CCND1 miR-302 0.09). But the lift is
**modest and uneven** (strong on PTEN, negligible on CCND1 where the couplers are already selected) and it can
still promote an expressed-but-non-coupling edge (PTEN miR-148, edge-ρ +0.17); coupling stays neutral (Δρ ≈ +0.01).
**Verdict:** the π role is now *defensible* (unlike rarity) but marginal → kept **default-off** pending the
discovery-precision validation. seq-spec's demonstrated strength is as an **identity/attribution readout** (role 1)
and an orphan **discovery-candidacy** signal (role 2, genome universe) — both still to be wired. AGO **loading(arm)**
remains a separate functional-abundance leg.

A separate, better-founded *functional-abundance* axis (AGO **loading(arm)** — RISC-loaded fraction, measured
non-circularly) is noted for its own leg; it discounts passenger arms on the abundance side, orthogonal to
specificity.

### 6b. The converged final model, in one line

**One dense learned-τ² non-negative Bayesian posterior per gene, two readouts of it:** run **π ≡ 1 (dense)** for
**coupling + attribution + identifiability** (the §1–2 winner, prior-inert by design), and **evidence-π
(inclusion)** for **discovery/selection**. Slab-scale carries μ·τ only. The **seed-rarity** channel is **off by
default** (§6a: locus settled = inclusion-π, but a global rarity gain promotes obscure/unexpressed seeds, not
coupling specialists — parked pending an expression-gated, coupling-validated redesign). Nothing else in the
estimator space moved the needle (§2g): dense-vs-inclusion is the sole axis, slab shape and CV-vs-EB are
second-order. This is the model to cite as *the* learned model; the adaptive lasso is retired to a baseline.

### 6c. Attribution-IDENTITY — the final model, validation, and verdict (2026-07-09)

The specificity work converged into a standalone **attribution-identity** model (`attribution_identity.py`) —
the abundance-removed "who OWNS gene g" readout — distinct from the regression's attribution-*magnitude* (M).

**Statistic (superseding concentration-share):** `identity(m,g)` uses the **affinity PERCENTILE** of g within the
family's genome targetome ("is g among the family's strongest targets"), NOT concentration-share (share is
confounded by targetome size — it crowned miR-96 over miR-200 on ZEB1; percentile is scale-free and captures a
hub-family-with-a-famous-target). context++ beat K_D (checked: K_D scores weak/6mer sites, worse; don't combine). ⚠ **SUPERSEDED (MH-87):** that check was HE-restricted — on the GENOME-WIDE universe all-site scanMiR K_D BEATS context++ at per-arm specialist recovery (0.89 vs 0.79 fair, `kd_fair_bench.tsv`), and strong-site thresholding HURTS. K_D is now the per-arm specificity (`kd.genome_affinity`, `seq_specificity.affinity_percentile_kd`) + discovery biochemical support; context++ stays L1 (family-level, where K_D's per-arm edge is muted).

**Two levels (match the model's estimand):**
- **L1 FAMILY ownership** = softmax over TargetScan **seed families** of z(affinity-percentile). Co-seed matures
  target identically, so ownership is a seed property. **Seed-key FIXED 2026-07-09:** affinity is now keyed by the
  canonical seed family via `seq_specificity.seed_family_map()` (from `data/miRNA/miR_Family_Info.txt`, = FAM
  grouping) — previously miR-number-keyed with a MAX-over-members hack that inflated multi-member families. The
  core model (`FAM.family_of`) was always correct; the bug was isolated to the identity affinity table (now fixed).
- **L2 ARM distribution** = who within the family, from **pooled chimeric provenance** (`chimeric_evidence.py`:
  Manakov + TarBase qCLASH/CLASH, per-source or overlap-resolved `mean_norm`), else expression×loading fallback.
  Combined arm ranking = P(family)·P(arm|family).

**Validation (as far as available data allows — honest):**
- **Face-valid**: recovers canonical specialists (miR-19/26→PTEN, miR-15/16→CCND1, miR-200/141→ZEB1), abundance-
  removed, and names the physically-bound mature.
- **Physical binding**: Manakov chimeric p=7e-8 **did NOT independently replicate** (Manakov-exclusive TarBase
  p=0.53) → treat as inconclusive. **Breast-matched** (Farazi MCF7 PAR-CLIP p=7.5e-7; MCF7+MDA HITS-CLIP p=5.9e-5,
  modest) — **non-circular through HE** (CLIP ∉ HE, verified <2% overlap; HE = functional evidence only) but
  **seed-circular** (percentile & CLIP both seed-driven; breast CHIMERIC, which would break it, does not exist).
- **Pan-context L1**: chimeric p=2.9e-5, AGO-CLIP p=8.4e-3 (too dense). **Arm-level L2: 91% cross-source agreement**
  (Manakov vs TarBase-chimeric on the guide arm) — the robust part. ASYMMETRY: arm-resolution (L2) reproducible;
  family-affinity (L1) modest+seed-circular.

**Model-value verdict (the honest bottom line):** identity earns **NO value in the TCGA regression** — inert on
coupling (p0×spec_gain sweep: spec_gain null, p0 owns coupling → dense; `genome_p0_spec_sweep.py`), on
attribution-magnitude (data-driven), on discovery-precision (identity⊥coupling; `seq_spec_bench.py`,
`discovery_dossier`), and on **identifiability** (`identity_payoff_genome.py`: 0% top-driver flip genome-wide;
evidence-π already saturates). It is a **descriptive standalone model** (who's *designed* to own g), a
hypothesis-generator — NOT a model-improving prior. AGO **loading** (`ago_loading.py`) is validated (correct
guide biology; passenger edges couple 2× weaker p=1e-3) but likewise coupling-inert; its value is arm-QC +
discovery passenger-filter. The seq-spec inclusion-π channel is retained but **default-off**.

Modules: `seq_specificity.py` (+ `seed_family_map`, `affinity_percentile`), `attribution_identity.py`,
`ago_loading.py`, `chimeric_evidence.py`; benches `rarity_bench`, `seq_spec_bench`, `genome_seq_fdr`,
`genome_p0_spec_sweep`, `identity_payoff_genome`. Data: `data/external/breast_parclip_farazi2014/`. See
`he-definition-and-evidence-weight` memory + `DATA_SOURCES.md` §interaction-evidence.

### 6d. ⭐⭐ THE DECOY CONTROL — what the fit actually buys (MH-127, 2026-07-13)

**The model's in-cohort advantage over abundance is a FITTING effect, not a curation effect.** Learn the same
dense β on a **matched fake** regulator design — abundance-matched, **site-free** (no TargetScan site *and* no
scanMiR duplex), not a seed-family mate, not a genomic-cluster mate, |r|<0.30 against every real regulator —
and the **fitted decoy reaches 79–91%** of the real model's 5-fold-OOF gene field, with the paired real-vs-fake
test **not significant** (Δρ = −0.0098, 52.5% of genes, **p = 0.211**; structure-matched decoy p = 0.577;
MDE 0.0161 ⇒ *underpowered, not proven zero*). The decoy **gains from fitting as much as the real set does**
(−0.0210 vs −0.0198) and then **beats real unfitted abundance** (p = 0.054 overall; **p = 3.6e-04** on
multi-regulator genes). ⇒ **MH-115's "the learned β beats mere abundance" is RESTATED: an abundance baseline
is not a control, and de-confounding + joint apportionment are properties of *any* fitted design, including an
impossible one.**

**What the curated design does buy, measured: TRANSFER.** Frozen TCGA β predict **CPTAC mRNA** (ρ = −0.0189,
**p = 8.9e-04**) while the **matched decoy's frozen β do not** (−0.0056, p = 0.19); paired **p = 0.030 / 0.027**
under the two decoys, and **p = 0.018 composition-adjusted** under the structure-matched one. The CPTAC
**protein** layer cannot arbitrate — the *real* model is itself null there under a from-scratch, correctly
family-pooled field, and the test is underpowered. **The stratified rescue fails outright** (no regulator-count
or Shapley-class stratum separates real from fake; the class ordering is indistinguishable from shuffled
labels). The **only** stratifier that separates is **fit-free** (the top tertile of curated evidence weight,
where the *unfitted* real design already beats the unfitted decoy) — and it does not replicate under the second
decoy, so it is a lead, not a result. Full account: `MH127_DECOY_MODEL_GENE_BUDGET.md`.

---

## 7. Artifacts & reproduce

**Modules** (`mirna_hallmark/learned/`): `spike_slab.py` (`_gibbs_posterior`, `_gibbs_ss`, `_gibbs_blasso`,
`_gibbs_ss_blasso`; `learn_pi0` Beta–Bernoulli), `priors.py` + `seed_rarity.py`, `regression.py` (lasso
baseline), `shrinkage_compare.py` (coupling grid), `attribution_eb.py` (attribution/identifiability/selection
genome-wide), `eval/bayes_parity.py` (CPTAC parity, shuffle null, discovery grid, K_D), `eval/retrofit.py`,
`eval/dossier.py`.

**Outputs** (`mirna_hallmark/output/learned/`): `shrinkage_compare.tsv`, `attribution_eb_full.tsv`,
`identifiability_vs_partial.tsv`, `selection_lasso_vs_bayes.tsv`, `gate_fdr_spike_slab.tsv`,
`discovery_bayes_*{,_kd}.tsv`, `discovery_consensus_{targetscan,kd,combined}.tsv`,
`discovery_dossier{,_gold}.tsv`, `retrofit_grid.tsv`.

**Provenance:** MH-82 (coupling/attribution/identifiability + estimator sweep), MH-83 (gap-closers: CPTAC
parity, shuffle null, genome-wide identifiability, discovery-vs-partial), MH-84 (discovery policy, retrofit,
K_D candidacy, consensus, dossier, tier-3 protein). Full decision log in `LEARNED_MODEL_RATIONALE.md`
§2a–§2g; numbers in `LEARNED_MODEL_VALIDATION.md` §1; run ledger `ANALYSIS_RUN_LEDGER.md`.

---

## 8. What is NOT built (honest edges)

- **Pipeline wiring** — discovery → retrofit → single fit is proven but run as separate steps, not one
  orchestrated command.
- **K_D independence caveat** — K_D candidacy shares the scanMiR source with the gate; TargetScan stays
  primary. `densez`-K_D is high-FDR (0.149) and dropped.
- **Protein-layer power** — CPTAC protein ρ is modest (mean −0.07; n~100, translational buffering); it's a
  *direction* check (66% negative ≫ chance), not a magnitude claim.
- **Dossier expensive columns** — CN-instrument causality, methylation mechanism, subtype-specificity are
  defined but run only on demand, not pre-computed over all 6,744.
