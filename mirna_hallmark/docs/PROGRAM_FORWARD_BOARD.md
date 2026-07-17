# Program forward board — the whole open-work map (`mirna_hallmark`)

> **Goal:** ONE consolidated view of every open forward thread across the whole subproject — edge existence,
> the decoy/claim-1 arc, attribution, cross-cohort validation, the discovery lane, the DCIS/EV lane, and infra —
> so nothing is lost between the detail docs. This is the *index of what's next*; the *depth* lives in the
> registry row named on each line.
> **What belongs here:** a one-to-three-line entry per open thread, with a status marker and a pointer to its
> `MH-##` row or home doc. NOT methods/rationale (→ their home docs), NOT results (→ `DISCOVERY_REGISTRY.md`),
> NOT done-and-closed items except as compact ✅/⛔ anchors that stop someone rebuilding a dead thing.
> **Update trigger:** when a thread is started/finished, or a new thread is opened anywhere.
> **Sync-partner:** `docs/STATE_OF_PLAY.md` (the current-state verdicts this board plans against; if the two
> disagree, STATE_OF_PLAY wins, and the registry wins over both).

**Last updated: 2026-07-17** · planning against the registry through **MH-148**.
*Supersedes and replaces `WHATS_NEXT.md` and `LEARNED_MODEL_WHATS_NEXT.md` (both archived; both predated MH-133…137).*

Status: ✅ done (anchor only) · 🔨 immediate/flagged · ⬜ open · 🔬 investigate · ⛔ blocked/dead · 🚫 deliberate skip.

---

## 🔨 DO THESE FIRST

Highest value ÷ cost. ⚠ Two of the original four turned out to be claims that did not survive checking
(MH-38's 108 orphans were 3; E6's "live bug" was never live) — verify the claim before building the fix.

1. ✅ **DONE 2026-07-17 — `eval/buffa_validation` re-run: the "108 triple-validated orphans" ARE 3.**
   `triple_validated` **108 → 3**; **ECM/collagen among them 30 → 0**. Survivors (all guide arms):
   **miR-195-5p→PSMD7 · miR-181b-5p→IRS2 · miR-22-3p→STK39**; 1 neg-sig in Buffa. ⚠ **It required a code fix —
   the module had been DEAD since the `eval/` reorg** (`GEO_DIR` was a `__file__`-relative hop pointing at a
   non-existent `mirna_hallmark/data/`; the cache was at the repo root). **That is why it was never re-run.**
   HE legs reproduce exactly (+0.319 / 0.593 / 0.673) = the control. `talk.qmd` updated. **(→ MH-38)**
2. 🔨 **Annotate MH-38 / MH-55 / MH-73 / MH-74** — all four still carry **S/R** strength tags with **no stale
   marker**. MH-38/55: input collapsed 594→23. MH-73/74: superseded by MH-76's frozen-panel test. Per the
   one-home rule, the retraction banner goes in *their* rows. **(→ CLAUDE.md doc protocol §4)**
3. ✅ **NOT A BUG — E6's `(live)` tag RETRACTED 2026-07-17 (MH-142).** Verified in code and data: every call
   site scores **orphans**, disjoint from `he_agg` by construction — **0 of dossier's 6,744 candidate pairs are
   HE edges** (vs 5,938 HE). Nothing self-partials in production. E6 **stands as a trap for CLASS-MATCHED
   comparisons** (its evidence: *"class-matched estimator comparison"*) — MH-124's own experiment pushed HE arms
   through this estimator and got a ~36× depressed HE rate **from its own design**. ⛔ **Do NOT remove `he_agg`**:
   it is what makes an orphan's coupling *partial on what curation captures*. Guard: `eval/_e6_production_check.py`.
   **(→ MH-142)**
4. ✅ **DONE 2026-07-17 (MH-140) — `readouts.share` relabelled AND the table gained a TRUE identity column.**
   `identity` = `shapley_identity` (LMG on R², m = bagged NNLS); `share`→`beta_frac` (demoted to MAGNITUDE);
   + `identity_deconv`, `identity_eq_magnitude`. **WHO ≠ HOW-MUCH in ~25% of multi-family genes** (agree 75.2%,
   n=819; was 99.35% when "identity" was argmax β). `m_nnls`=0 in 31.7% of families vs `beta`=0 in 0/5,117.
   213 genes (13.8%) now honestly report identity = UNDEFINED. Change is purely additive (5 columns
   bit-identical); `shapley_identity` optimized 8.2–9.7× output-identical. **(→ MH-140, MH-138)**
5. ⬜ **`β(TCGA) → Buffa mRNA`** — the learned model's **only untested boundary is the cohort boundary**
   (MH-104: **~80% loss at the cohort boundary, ~0% at the layer boundary**), and Buffa (n=207, miRNA+mRNA) is
   the only independent-patient cohort in the repo that can test it. **β has never touched Buffa** — one `STUB`
   docstring in `learned/eval/__init__.py`; no `ood_cohort.py`. Inputs are on disk.
   ⚠ **Pre-register MH-114's compartment-orientation stratification** — both cohorts are bulk breast and share
   the CAF confound, so a clean replication proves nothing on its own. **(→ axiom 4; MH-104, MH-114)**
6. ⬜ **CURATED-EDGE ADMISSIBILITY — can this HE edge repress IN BREAST AT ALL? (user-proposed 2026-07-17)**
   An HE edge earns its place by a luciferase assay somewhere — not necessarily in **this tissue**. Two a-priori,
   Y-blind, model-blind gates ask whether it *could* operate here:
   * **SITE:** does the arm have a seed site in the target's 3'UTR? (TargetScan context++ ∪ scanMiR duplex ∪
     Poisson-gated 6mer — `decoy_bench._site_maps`, already built.)
   * **EXPRESSION:** is the arm expressed in breast — **tumour (TCGA) OR healthy (GTEx)**? (`_FLOOR = log2(11)`,
     RPM≥10, `card.py`; GTEx breast is ingested via `gtex_mirna_matrix.py`. **Healthy matters**: an arm silent in
     tumour but present in healthy breast is a *lost* regulator, not an inadmissible one.)
   **An edge failing BOTH cannot repress here by any modelled mechanism** ⇒ it is curation TRANSFERRED from
   another tissue/cell line, and it is sitting in the design diluting the real arm.
   **⭐ HALF OF THIS IS ALREADY MEASURED, and it lands:** MH-136 — **187 genes whose curated edges are ALL
   seedless give gap +0.0006, EXACTLY ZERO**, vs −0.0325 for 810 all-seeded genes. *No site → no repression →
   no coupling → no gap.* MH-132 independently found the 21% positive-coupling genes are **~9/13 positive in
   GTEx NORMAL BREAST too** ⇒ *"a **curation-transferability** problem, not a model defect"* — the same claim
   from the other side. **The expression half has never been tested.**
   **THE TEST (cheap — both artifacts exist):** tag every HE edge `admissible = has_site ∧ expressed_breast`;
   re-run `decoy_bench` on the admissible subset. **PREDICTION (register before running): the gap GROWS** —
   inadmissible edges are noise in the REAL arm, so removing them should sharpen the contrast beyond MH-137's
   −0.0119 / MH-139's −0.0129.
   ⚠ **Guardrails.** (a) This is **DIAGNOSTIC first, a filter second** — axiom 2a: flag, don't delete. A gate
   that removes 30% of the universe needs its own null before it changes a headline. (b) It is **NOT** the
   decoy's site-free rule (that defines the FAKE; this classifies the REAL). (c) Report on the gated set AND
   the complement — if the complement is not ~0, the gate is not what we think. (d) `n_fam` shrinks under the
   gate, and MH-147 showed the gap **scales with design width** — so re-check the width axis inside the gated
   set, or the gain is just width. **(→ MH-136, MH-132, MH-147; `ARM_EXPRESSION_FLOOR.md`, `SILENT_ARM_REMOVAL.md`)`**

7. ⬜ **Restore an exogenous existence validation** — see §A. Not cheap, but both `CN_INSTRUMENT.md` and MH-126
   name it **the single highest-value open item in the program**, and nothing else changes the program's
   standing as much.

---

## A. Edge existence — the one open foundation

**⛔ Both CN instruments are RETRACTED. Edge existence rests on ONE observational line** (curated edges
anti-correlate more than abundance-matched site-free ones). The "three independent validations" headline is
**FALSE** — the other two legs were both observational and never independent of each other.
**(→ STATE_OF_PLAY Axis 2; MH-124r/126, MH-133)**

- ⛔ **`pi_causal` — DEAD, do not rebuild.** `pi_causal = γ·b` with **b an observational OLS slope** ⇒ a
  product-of-coefficients mediation estimator, not an IV; it re-derived the anti-correlation it was meant to
  validate. Clean reduced form, arm-clustered: **p=0.115, n.s.** Retraction is in production code at
  `learned/instrument.py:1128–1148`. `MH124_*` §5/§5b are **VOID — do not cite** (§4b is untouched).
- ⛔ **Within-arm two-way-FE replacement — DEAD.** Its control class was **71% real binders**; against a
  genuinely site-free control **τ=−0.0007, p=0.84**, and it fails the site-type efficacy ladder
  (8mer/7mer-m8/7mer-A1 not monotone, p=0.26). **Refutations at n=216k–235k pairs, not power failures.** (MH-133)
- ⬜ **THE OPEN ITEM: a new exogenous handle.** The identification strategy must be an **asymmetry the confound
  is blind to** (axiom 4). Copy number cannot see a site type — that is what killed both CN designs; whatever
  replaces them must pass the same ladder.
- 🔬 **AI2 — second instrument: germline-eQTL / methylation-as-instrument.** Mapped, unbuilt. Over-identifies
  the same β and breaks reverse-causation differently than somatic CN (a stronger Sargan/Hansen jointly).
  Needs an ancestry-PC C block (axis BQ). **The leading candidate for the item above.** **(→ CHANNEL_FUSION §9 AI)**
- ⚠ **Bounding this axis regardless:** the per-edge null is **3–4× too narrow**; site-free pairs pass "FDR
  q<0.05" at **25–35%**; under an honest empirical null **0.0%** of curated HE edges survive per-edge. The
  signal is a **set-level distributional shift, not per-edge significance** (MH-123/124). **Any FDR count in any
  doc predating this rests on the uncalibrated null** — see §G for the two docs that still need banners.
- 🚫 **Do not carry forward** the mapped-but-unbuilt CN-instrument refinements (BL subclonal/CCF-weighted CN,
  AD3 allele-specific CN on imprinted loci). They sharpen a design whose exclusion restriction is dead; they are
  not a route back to identification.
- ⚠ **Note on the CN *channel* (distinct from the CN *instrument*):** the fusion channel was only ever
  "real but weak" (~0.7% of β's information), and **the shipped readouts do not include it at all**
  (`readouts.py` calls `_gibbs_posterior` without `channels=`). The soft-F-weight gate (`_ADMIT_SOFT_C=10`) and
  the coding/host pleiotropy gate are built and correct but currently **do not touch the production table**.
  Decide deliberately whether to wire them or retire them — do not leave this ambiguous. **(→ MH-91, MH-99/101)**

## B. Decoy / claim-1 — the surviving positive result

**MH-137 is the row of record and supersedes every earlier magnitude: gap −0.0119 (p=3.4e-05), deconv −0.0090,
retention 0.76, 1,322 genes.** ⚠ **The magnitude has SHRUNK at every control fix:**
**−0.045 → −0.0306 → −0.0147 → −0.0119.** The honest statement is a **BOUND, not a win**.

- ⬜ **Re-run the competence map's gap columns through `eval/decoy_bench.py`** (2:22, full universe). Every
  real-vs-fake number in `gene_competence_map.tsv` sits on a **discredited decoy**; the `gene_atlas.tsv` half is
  fine. **(→ MH-130a, MH-137)**
- ⬜ **The next control fix: Manakov chimeric eCLIP + POSTAR3 are still uncovered by the evidence exclusion.**
  Whether closing that hole takes the gap to zero is **genuinely open** — it is the arc's decisive remaining test.
- ⬜ **The evidence axis (`w_max > median`) has NEVER been tested on a clean control.** Its only support fails
  the arc's own both-fake-sets rule (q=0.006 FAKE1, q=0.20 FAKE2). ⇒ carry MH-130's **"27% domain" as a live
  hypothesis, not an established partition**. The **width axis (`n_fam ≥ 3`) is CONFIRMED** and sharpened
  (3–4 families −0.0314, p=0.0006) — ⚠ but **not monotone under MH-137: it peaks at 3–4 and falls at 5+**.
- ✅ **Anchors — do not re-litigate:** MH-136's **seedless positive control** (187 all-seedless genes, gap
  **+0.0006 = exactly zero**; 810 all-seeded **−0.0325**) is the arc's first real positive control and stands ·
  the **1-family internal null is back** (−0.0007, p=0.15 on the clean control; MH-135's retraction of it was an
  artifact of the broken decoy) · **an abundance baseline is NOT a control — benchmark any aggregator against a
  fitted matched decoy** (MH-115/127) · MH-130b's whole **SNR-strata arc is dead** (double-log bug);
  `mh127_snr_strata.tsv` is poisoned — delete or re-derive.
- 🔬 **Aggregate "raw force" vs abundance-sum, tumour↔GTEx Δρ.** Gene-rung design in
  `method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md`. **Two gates before any run:** it sits
  on the **retired pressure heuristic**, and its comparator is an **abundance baseline, which the rule above
  says is not a control**. Re-scope against a fitted matched decoy or drop it.

## C. Attribution / Shapley

**STANDS: the model RANKS the canonical family high; it does NOT NAME it.** Cite **MH-126 Test 1** (n=32), not
MH-124 §4b (n=21). β is at chance (rank 0.518, p=0.66); `shapley_identity` 0.370 (p=0.018); abundance 0.240.
**Argmax is at chance at every k.** **(→ STATE_OF_PLAY Axis 3)**

- ⬜ **Scale the literature ground-truth set — this, not the estimator, gates attribution.** n went 16→21→32 and
  the rank test became decisive while argmax stayed at chance. Needs a **versioned, auditable ground-truth pull**.
  Until then: *"β does not attribute" is MEASURED; "the CN instrument cannot attribute" is UNDERPOWERED, not refuted.*
- ⬜ **One canonical rerun of the persisted attribution outputs** (`output/learned/programs/*`, cards) now the
  estimator is settled. Optional: canonicalize the `shift_vs_weight` per-edge weight diagnostic.
- 🔨 **See DO-FIRST #4** — the `share` mislabel is an attribution-doctrine violation shipped in the genome-wide table.
- ⚠ **Do not build a "do canonical regulators score higher?" test on `pip_discovery` or `prior_pi`** — both are
  **w-contaminated by construction**; the test is circular. (β itself passes the w-circularity gate
  bit-identically: max|Δβ| = 0.0 under shuffled *and* constant w.)
- 🚫 **Do not carry "dose weighting is dead" from MH-130e into the attribution arc** — it refuted the extension
  of §4b to the **coupling** estimand only. §4b's attribution claim is a different estimand and is not overturned.

## D. Cross-cohort & external validation

- ⬜ **`β(TCGA) → Buffa mRNA`** — see DO-FIRST #5. The program's only untested boundary.
- ⚠ **"METABRIC" in MH-73/74/75/76 means Buffa, n=207** — verified in code. The +0.056 / +0.060 / +0.028 /
  +0.019 are **four views of ONE 207-patient result, not four validations.** Do not recount them.
- ⬜ **METABRIC-full: file the DAR or stop citing it.** "EGA-pending" appears in **five docs** with **no
  accession, no DAR id, no submission date, no owner** anywhere in the repo. The evidence supports *"never
  applied for"*, not *"pending"*. ⚠ It would also be a **power** upgrade, not an **independence** one — same
  platform, same consortium; Buffa is a subset of it.
- ⬜ **CPTAC composition re-runs still pending (MH-107):** **MH-34/35/36/39/40 · MH-83 · MH-84 tier-3.** The
  code is fixed (`cptac_validation._covariates` 2→11 covariates, now RAISES rather than running confounded);
  the persisted outputs are not. ⚠ MH-114 supersedes MH-107's magnitudes — re-run against MH-114's framing
  (~57% of the protein validation is compartment arithmetic; true compartment-blind effect ≈−0.011, surviving
  adjustment at ≈−0.0065, p=2.0e-03).
- ⬜ **The CPTAC protein transfer is NOT decoy-controlled** — a site-free fitted fake reaches **99%** of it
  (MH-130d). It is real arithmetic but **not evidence the curated edge set is right**. *"The protein channel has
  never faced a decoy that could lose"* — the arc's own words. The only decoy-controlled transfer is the
  **CPTAC/mRNA** cell inside the competent class.
- ⬜ **V1 δ-transportability, remaining scope.** The *abundance* input transports (same member dominates 84.5%
  vs 43.7% chance; member-share correlation 0.991, MH-104). The **full CN/chimeric-fused δ** is untested and
  low-powered — the CPTAC CN instrument is F>10 for only **59/685 arms (8.6%)** vs TCGA's 32%; **49 arms / 5
  families** usable in both.
- 🔬 **Manakov chimeric-duplex validation of the discovery shortlist** — also the input to §B's evidence-hole fix.
- ⚠ **Standing rule for any future cross-cohort channel:** run `gauge.calibrate()` → read `info_ratio` + τ
  **BEFORE** building it. This is what killed the state channel; it is what falsifiable infrastructure is for.

## E. Discovery lane

- 🔨 **Per-edge (+ per-gene) attribution card** — one row/edge joining what is already computed: regime context
  (arm median RPM, % samples above the RPM≥10 floor, IQR, **spiker flag**; target median + IQR), budget share
  (E7/G4 `states.budget_shift`, GTEx→NAT→tumour Δ), composition tag (deconv retention: cell-intrinsic ≥0.7 /
  partial / composition-explained <0.4), and shift-class (`mirna_state_class.joint_edge_class` extended with
  paired-E5 realization + a low-variance/undetectable split). Per-gene card = the G-series roll-up (G1 net
  repression, G8 coherence+role, G4 budget concentration, composition fraction, G10 shift-class). Pure join.
  **(→ `ATTRIBUTION_CONTEXT_AXIS.md`)**
- 🔬 **Top-discovery deep-dives.** (a) the **miR-17~92 / miR-106b-25 target cluster** (miR-106b/93/17/19a →
  RABEP1, AHNAK, TGFBRAP1, IL6ST, AFF1, WWP1, M6PR — ubiquitous arms, cell-intrinsic, strong, weakly-curated):
  lit scan each, then CPTAC protein + Manakov binding + within-patient realization; is it a coordinated
  cluster→program hit? ⚠ note miR-93-5p is also the systematic **proliferation-confound arm** (MH-100) — check
  the `prolif_verdict` flag per edge before reading these as coupling. (b) **miR-135b-5p → GATA3** (fully novel,
  a **spiker**, 28%>floor ⇒ coupling lives in the miR-135b-high subset): which tumours, is the subset a clinical
  stratum? (c) the **568 fully-novel candidates** → triage (cell-intrinsic + expressed/spiker + biochemically
  supported), lit-novelty check, shortlist for protein + Manakov.
  ⚠ Every FDR-based candidate count here inherits the **uncalibrated per-edge null** (§A).
- ⬜ **Epigenetic-rewiring discovery — the bidirectional methylation lane.** Detector built
  (`structural_identity` + `methylation`; positive controls pass: miR-124/9/129/137 hyper→LOST, miR-21/155
  hypo→GAINED); needs genome-wide scale-up. **Axis A** = LOST specialist (baseline-active, tumour-silenced,
  promoter HYPER-methylated Δβ≥+0.15 → target de-repressed). **Axis B** = GAINED specialist (baseline-silent,
  tumour-active, HYPO-methylated Δβ≤−0.15 → target acquires repression). Nominate on `structural_identity`
  (`structural_share ≥ 0.10` ∧ `confidence ≥ 0.5` — expression-free, so it survives silencing); the ONE
  inferential claim is the **target consequence** (group Δ + within-tumour g-mRNA vs promoter β; permutation
  null, BY-FDR). **Must condition:** purity/composition (global tumour hypomethylation is the known Axis-B trap
  ⇒ demand CGI-promoter-local, not open-sea) · CN-vs-methylation (cross `mirna_locus_cnv`) · imprinted loci
  (DLK1-DIO3 — flag & exclude; miR-127 baseline β 0.92 is exactly this). **Prereq: the UUID→TCGA-barcode map.**
  Power-bounded (normal n=97); the tumour arm is deep.
- ⬜ **Target-side mechanism legs — the reference-blind quadrant.** `identity ≫ structural` is not only
  "opportunistic": it can be a coupling via a site the **reference genome does not encode**. Build the
  **seed-SNV leg** (APM has SNV/VEP): scan each gene's 3′UTR tumour SNVs for seed-match **gain** (acquired
  site) or **loss** (the target-side analog of the silenced specialist). Same family: **A-to-I editing** and
  **APA**. ⇒ the identity/structural quadrant is a **mechanism router** — loss side {methylation, seed-SNV-loss,
  CN-loss}; gain side {hypomethylation, seed-SNV-gain, editing/APA-exposed, uncurated-real, residual-confound}.
  Subsumes the old "cross EMT/NF-κB genes with APM SV/3′UTR disruption" item.
- 🔬 **context++ × K_D fusion.** The two are **nearly orthogonal — Spearman +0.065** over 78k (family,gene)
  pairs (K_D = scanMiR RBNS *binding*; context++ = TargetScan *repression-outcome* + conservation + site
  context) ⇒ exactly where a fusion can beat either. MH-86's "don't combine" was HE-conditioned and used the
  worse (HE) K_D. Test for specificity + discovery candidacy. **(→ MH-87)**
- ⬜ **K_D scan scope + discovery floor sweep.** Current genome-wide K_D = **746 detectable arms** (RPM≥10).
  Optional: (1) a complete **~2,600-arm reference scan** (one `scanmir_genomewide_par.R` run) as a reusable
  `(arm,gene)` resource; (2) a floor sensitivity sweep (RPM ≥ 50/10/1) — a re-*filter*, not a re-scan; the
  floor is the sensitivity/precision lever, not an FDR lever.
- ⬜ **CLIP-source candidate expansion** (ENCORI / POSTAR3 physical sites) beyond TargetScan nomination.
- ⬜ **Regulation-pattern archetypes + seed-rarity.** (a) Tabulate the per-gene pairwise-ranking concordance
  matrix over the four rankings (abundance · budget · identity · structural) → cluster genes into archetypes
  (abundance-controlled / weight-controlled / designed-but-silenced / redundant-collinear) — the systematic
  version of the ad-hoc PTEN (structural anti-abundance ρ=−0.65) / ESR1 reads. (b) Ground **structural
  specificity in targetome RARITY**: how many genes carry the arm's 7–8mer genome-wide, tied to scanMiR K_D +
  context++ (rare strong site = genuine specialist; ubiquitous seed = pseudo-specialist) — makes "designed to
  *specifically* repress g" a sequence claim. (c) **GLOBAL** rank concordance via `program_network` engagement
  (Σ_g M) — loud-but-unwired vs scarce-but-hub arms at genome scale.
- ⬜ **State-stratified discovery** — run the discovery lane **per state** (own permutation-FDR null + scanMiR +
  deconv tag), then compare: found-in-tumour-not-NAT = candidate **acquired** regulator · found-in-all =
  constitutive · found-only-NAT = field/lost. **Distinct from the wiring axis** (which uses a *fixed* support to
  *compare* known edges). Power-bounded: NAT n=104 / GTEx n=327 ≪ tumour ~1000; and *undetectability* is worst
  in healthy ⇒ "not found in healthy" ≠ "not a healthy edge". NAT is the more feasible end.

## F. Progression / state — the landscape, not the channel

**Two different objects. Do not conflate them.** **(→ STATE_OF_PLAY Axis 6)**

- ⛔ **The state CHANNEL is CLOSED — measured and cancelled. DO NOT REBUILD.** `learned/channel_state.py` was
  never built (verified absent). τ — the acquired-vs-constitutive payload, the axis's **entire deliverable** —
  is indistinguishable from zero (GTEx τ=0.0009, info **0.6%**; NAT 0.7%; CPTAC 0.7%). Structural: channel
  precision ∝ **a²**, and a≈0.11 ⇒ **even unlimited GTEx donors add ~1%**. The additive `M_T = M_H + Δ` form
  **would have faked the headline** (Δ>0 on nearly every edge, manufactured out of composition attenuation).
  **(→ MH-102d, `LEARNED_MODEL_STATE_CHANNEL_PLAN.md §10`)**
- ✅ **The abundance-level cross-state LANDSCAPE is built and run** (GTEx-healthy → TCGA-NAT → tumour, 721 arms
  / 5,108 edges; acquired_realized 640 · lost 209 · stable 2,141 · acquired_unrealized 1,775).
- 🔨 **`LANDSCAPE_REPORT.md` needs a banner** — its per-edge FDR class labels rest on the null MH-124 measured
  as **3–4× too narrow**, and it carries no warning. ⚠ NAT is structurally underpowered (n≈104):
  `nat_decoupled` = **1 edge**.
- ⬜ **Abundance-vs-wiring decomposition** (Δx·M vs x·ΔM) as a posterior readout + **subtype interaction tests
  (Gelman–Stern)** — a readout over the landscape, not a channel.
- ⬜ **GTEx logFC gauge check** — `parallel_view.change_trajectory` reports both the naive RPM−TPM (`rawFC_GN`)
  and the QN-into-TCGA (`relFC_GN`) GTEx→NAT change; validate the miRNA RPM≈TPM assumption against the QN
  version to find where they diverge (housekeeping/spike arms).

## G. Correctness debt, docs, figures, infra

- 🔨 **Fix `LEARNED_MODEL_METHODS §1`'s rationale** (and its `_FORMAL` twin — always edit both). It justifies
  the C block by *"Cancer-Epithelial is deliberately EXCLUDED — conditioning on the compartment the target is
  expressed in over-controls the signal"*. **That is VOID:** the 9 Wu-major fractions sum to exactly 1.000000
  and **R²(Cancer-Epithelial ~ the 8 held-out) = 1.0000 in both cohorts** ⇒ the hold-out is a **simplex
  illusion**. The practice is nevertheless **safe for a different, measured reason**: purity contributes only
  **0.9%** of the retention loss while the stromal MIX contributes **33.1%**. **(→ MH-111/114)**
- ⬜ **Regenerate presentation F25 + re-render the deck.** `talk.qmd` prose was fixed 2026-07-12 and is honest;
  **the F25 figure beside it is from 2026-06-24** and shows the old numbers.
- ⬜ **Posterior calibration — ✅ DOC CONFLICT RESOLVED 2026-07-17; the fix itself is still open.**
  **SETTLED VALUE: widths are ~25% too narrow** — bagged NNLS **0.73×** · Gibbs Gaussian **0.77×** ·
  **Student-t ν=7 0.89×** (MH-102e, `learned/calibration.py`). The two docs never actually disagreed on the
  measurement — this board reports **reliability** factors and `CPTAC_PROTEIN_CHANNEL_PLAN.md` reports their
  **reciprocals** as **inflation** factors (1.37× / 1.29× / 1.13×). The lone outlier was
  `LEARNED_MODEL_STATE_CHANNEL_PLAN.md` §11(a)'s **"~40% (NNLS 28%, Gibbs 39%) / SBC promoted"**, measured on a
  **biased gene subset** (the scale-dependent sd-floor bug) — now **bannered SUPERSEDED at source**.
  **SBC stays RETIRED as the wrong tool** (the fault is the MODEL — unmodeled between-participant heterogeneity —
  not the sampler, so SBC would PASS and mislead).
  **To apply:** multiply reported widths by the inflation factor; **turning `nu=7` on is a live recommendation**
  (⚠ `readouts.py` does not pass `nu` — the shipped table is Gaussian, i.e. the 0.77× rung).
  **Triage:** inverse-variance **fusion weights are SAFE** (SE *ratios* cancel); **all OOF/permutation results
  are UNTOUCHED**.
  ✅ **MH-94 WIDTH SUB-CONFLICT — DIRECTION SETTLED 2026-07-17 (MH-146); the flagship is NOT.** MEASURED
  (`calibration.shapley_resampled_width`, 64 family-rows): the posterior width **UNDERSTATES by ~1.4×** on
  families that carry identity (median 1.37 at identity>0.05) ⇒ **the board's DIRECTION was right; the CPTAC
  plan's "1.8× too WIDE" is RETIRED.** ⛔ But the tool's own printed headline ("understates by 112.86×") is a
  **zero-denominator artifact** — 48% of rows have `posterior_sd ≈ 0` (zero-identity families) where the ratio
  explodes to 3881–7655×. ⬜ **Open:** gate that summary on `identity_mean > 0.05`, and re-run MH-94's
  flagship (PTEN miR-141/200a 0.77±0.41 ⛔ **[RETRACTED — MH-150, 2026-07-17: re-run at MH-94's OWN config gives **0.003 ± 0.012, PIP 0.042**. The SHARE collapsed 250×, not just its width; the exemplar is void. PTEN's identity is now miR-17~92 (0.501) + miR-21/590 (0.223) — the canonical regulators. Prime suspect: the `_rtnorm_pos` sampler bug, which MH-94 predates. DO NOT CITE.]** vs MH-102e's re-derived 0.243) at ITS OWN config — the configs here
  are not comparable, and 0.41→0.243 points the OTHER way. **(→ MH-146)**
- ✅ **CLOSED (MH-151/152) — the isomiR cache needs NOTHING.** Coverage of the model cohort is **99.9%**
  (1,078/1,079 participants; ZERO missing). The manifest's 18 extra Primary Tumor files are **duplicate
  vials**, not missing people — a rebuild and a download both gain nothing.
- ✅ **ANSWERED (MH-153) — YES, it survives, and the risk is elsewhere.** miR-17~92's dose is **corr 0.95**
  under `isomir=True` and **PTEN's identity argmax is UNCHANGED** (miR-18-5p; identity corr 0.92). **Donor
  families are BUFFERED BY POOLING** — losing ⅔ of one member's reads is absorbed by the rest.
- ⬜ **THE REAL EXPOSURE, newly surfaced (MH-153) — audit RECEIVERS and the ORPHAN MASS, not donors:**
  (a) small families have no pool to dilute a donation — `miR-543` gains **+351%** from jumped reads;
  (b) the `orphan` bucket **DELETES 4.1% of total RPM** (654.3M → 627.4M) — reads whose shifted seed
  matches no family are dropped from EVERY family. **`X_fam` does not conserve mass and nothing flags it.**
- 🔨 **BEFORE any OOF isomiR run (MH-153): fit the shrinkage target on TRAIN ONLY.** `seed_composition`
  Dirichlet-shrinks low-read cells toward the **cohort** family-distribution — a cross-sample X-leak (no
  Y-leak). Confined to low-read cells and the path is default-OFF, so nothing shipped is affected.
  **(→ MH-153, MH-152)**
- ⬜ **Infra:** `git add` the untracked `learned/` tree (a provenance hole, flagged repeatedly) · fold
  ENCORI/POSTAR3/Manakov into the ledger (**union, not summed**) · `baselines/` re-export shims · vectorize the
  per-sample segment×locus overlap in `mirna_locus_cnv` (interval join per chromosome — the main slow step in
  `run_all`) · covariate-protected ComBat batch **only if** a channel reintroduces batch (never naive dummies).
- ⬜ **Doc consolidation, continued.** This board absorbed `WHATS_NEXT.md` + `LEARNED_MODEL_WHATS_NEXT.md`.
  Still queued: merge `LEARNED_MODEL_METHODS` → `FORMULAS.md`; merge `MODULE_MAP` → `ANALYSES_CATALOG`; archive
  `LEARNED_MODEL_{DESIGN_RESPONSE,DISCUSSION_PROMPT,BUILD_PLAN}`. See `STATE_OF_PLAY.md`'s **doc-traps table**
  for which docs read as current but are not.
- ⬜ **Low-priority coordinate/coverage debt** (deferred, cheap, no claim depends on it): extend MirGeneDB coords
  to the high-abundance expr-only MIMAT tail (~1,418 GDC MIMATs lack loci; ~24 have mean log2RPM>2; lead gap
  **miR-625-3p**) · the MI*→GDC hairpin bridge (348/506 mapped, 158 zero-bp) + per-locus `pct_mapped` reporting
  · multi-set governance (set-size-weighted / Jaccard-deduplicated Hallmark scoring).

## H. DCIS / pre-malignant / EV lane

- ⛔ **HTAN Duke (Strand) DCIS atlas — the RNA + outcome instrument, ACCESS-BLOCKED.** Two uses: gene-level
  miR-29 targets (COL1A1/3A1/4A1/5A1, LOX/LOXL2, SPARC, FBN1) across normal→DCIS→IBC **and by recurrence
  outcome**; and the MH-48 abundance-asymmetry signature + **miR-185-5p** as a DCIS→recurrence predictor — the
  one question the public cohorts (3 events) cannot answer. Patient bulk sequencing is dbGaP/EGA-controlled
  (DAR), but **derived L3/L4 are often open and are sufficient** (L1/L2 not needed) — confirm the specific open
  files on the portal before committing. **Check periodically.**
- ⬜ **EV recipient-cell delivery — GSE297447** (the EV **mRNA** companion to GSE297448): are the targets of the
  *exported* miRNAs (miR-140-3p, miR-205) de-repressed in recipient cells / present as EV cargo? Closes the
  export→delivery loop the MCF10A EVs only half-answer.
- ⬜ **Rab27A / secretion-machinery axis** — does export-selectivity (cellular-down→EV-up) track
  Rab27A/secretory-pathway expression across progression? Separates selective packaging from bulk shedding.
- ⬜ **miR-185-5p** — the one arm surviving composition adjustment as an epithelial grade-decliner (MH-51,
  adj ρ=−0.75, q_adj=0.026). Cheap in existing data: its targets, Hallmark links, and whether it is a bona fide
  epithelial progression suppressor (vs the stroma-driven miR-29c).
- ⬜ **Stroma-aware compartment extension.** The framework (MH-1..47) is epithelial-centric, but miR-29c→ECM
  showed the invasion biology lives in non-epithelial compartments; the **myoepithelial-loss program
  (miR-145/-143/-126)** is the LATE invasion leader (MH-48) and deserves its own module rather than being
  adjusted away.
- ⬜ **The desmoplastic driver, if not miR-29** (MH-54's open end). ⚠ **The remaining clean test is IN-SITU**
  (spatial / single-cell CAF), **not another cultured-CAF set** — GSE196354 was inconclusive due to **culture
  drift** (FAP/ACTA2 only +0.15 NS, collagens flat), so cultured CAFs under-represent the in-vivo myCAF state.
- ⬜ **A DCIS-STAGED sorted/LCM small-RNA-seq accession** — the one gap the field does not fill. Scan persists.
- 🚫 **Do NOT** run a broad small-RNA-seq sweep, or acquire another plasma-EV cohort (the tissue-export
  signature does **not** translate to a robust circulating marker — miR-205-5p replicated in only 1 of 6
  cohorts, miR-140-3p in 1; plasma EV is multi-source). 🚫 Do **not** use the miRNA-inference tools (SiCmiR,
  STmiR) as evidence for the miR-29c→ECM axis — they infer miRNA activity *from target expression* and are
  **circular** for exactly this question; landscape/hypothesis only.

## I. ⛔ Dead anchors — do not rebuild

One line each. The detail is in the registry row named.

- ⛔ **`βᵗ`, the CPTAC translational-repression latent** — the plan's centrepiece, **not supported at n=101**;
  all prior positive evidence was a **LEAK**. Four leak-free framings agree; the aggregate (1 df) gives BH
  q<0.10 in **1/17 genes — BCL2 alone**; **PTEN d=−0.006, p=0.82**. ⚠ **MH-34 is NOT overturned** (13/17 stay
  repression-directed, binomial p≈0.025) — the falsification is of the **modelling object**, not the biology.
  **There was never a βᵗ field. Only n fixes this.** **(→ MH-103)**
- ⛔ **Protein as a coupling lever** — it carries **4–6%** of the mRNA channel's Fisher information about β, and
  **≤7.6% at any `a_g`**. It must justify itself as a **different question**, never as a coupling gain. **(→ MH-108)**
- ⛔ **The state channel** — see §F. **(→ MH-102d)**
- ⛔ **"TCGA is sparse so it cannot arbitrate"** — MH-76's **frozen-panel** test removes overfitting as an
  explanation entirely: **DFI +0.002, OS +0.006; panel-alone C 0.48–0.53 ≈ random.**
- ✅ **Pressure magnitude is exhaustively non-prognostic** — 0 recurrence FDR across 3,368 features. A
  **well-earned negative**, credible because **miR-210 → DRFS reproduces Buffa 2011** as a working positive
  control. Functional > magnitude (realized +0.060 vs received-pressure −0.016, a *within-cohort* contrast).
- ✅ **CIBERSORTx / the panel question is SETTLED — ZERO new runs.** β is **ρ=0.94** across entirely different
  atlases (Wu-9 vs HBCA-11), ρ=0.99 across reference resamples, but only ρ≈0.58–0.81 vs **no C at all** (mean
  |β| halves) ⇒ ***which* composition control you use moves β ~10× less than *whether* you control it.** Wu-9
  stays the panel for all four cohorts. **Do NOT run** a fine 42–49-way panel, TCGA×HBCA, or refsample
  "harmonisation". Two gaps are **not** deconvolution problems: **adipose** (no reference has an adipocyte class
  and GTEx breast is adipose-dominant ⇒ needs a marker metagene) and **immune C-completeness** (Wu-9 lymphoid
  tracks Thorsson r=0.31 vs LM22's 0.70 — ⚠ LM22 + Thorsson are **validators, never covariates**).
- ✅ **The n≈1000 ceiling** — every internal lever (CN channel, δ-pooling, isomiR refit, Student-t, cross-gene
  pooling, the A3/A4 field) landed **"immaterial at n≈1000"** ⇒ **the frontier is exogenous, not more internal
  refinement.** Don't re-open them.
- ✅ **Engine: Gibbs, not HMC** — measured, not aesthetic. ⚠ **Carve-out:** fit the *gauge* with NNLS (Gibbs's
  heavy-tailed posterior SDs break the errors-in-variables correction). **Only a future non-conjugate channel
  reopens this** — binding = Poisson/NB, methylation = Beta.
- ✅ **`_rtnorm_pos` sampler bug — FIXED, don't re-litigate** (contaminated 3.15% of β; Robert (1995) rejection
  sampling, KS-verified). ⚠ It **overturns MH-119's recorded cause**: `share`=43.7 was this bug, not
  "anti-repressive biology".
- 🚫 **Occupancy / ceRNA saturation substrate** (Decision C/G) — failed its held-out gate; Phase-4 cooperativity
  sits on it (gate HARD on held-out if ever revisited). Likewise the `priors.py` spike-and-slab π/μ/τ object
  (permutation-FDR discovery already *is* the inclusion inference), the A/D/F pooling partials (verified no
  gain), and any channel whose link isn't identifiable.
- 🚫 **Fusion honesty-layer readouts** (axes U/AO/AP) and the **binding channel** (already a prior, `priors.μ` —
  the fusion frame only renames it): free readouts, but they **materialize only if a 2nd non-redundant channel
  lands**. Per §A and this section, none is on the board.
- ⬜ **The one design partial still worth having: B, the transcription-rate proxy** (the design's "most
  load-bearing"). An **intronic / pre-mRNA-read** nascent-transcription proxy as a C column, so a positive M
  isn't a TF co-driving miRNA and target. **HIGH value but RISKY** (over-control — the TF-regulon proxy failed
  exactly this way) and **data-blocked** on intronic reads. Check availability first.

---

## Open questions worth a dedicated study (not scheduled)

- 🔬 Is the basal **immune-pressure tone** a *cause* (miRNA-driven IFN silencing) or a *consequence* (low
  infiltration → low immune gene/miRNA co-program)? ⚠ post-MH-107/114, any version of this must run under the
  composition C and report retention, or it is compartment arithmetic by construction.
- 🔬 Does **8q24 amplification** (miR-151a, AGO2) create a coordinated "amplicon-driven RISC + guide" state that
  amplifies pressure non-linearly?
- 🔬 Are the **differentiation-program** repressions (estrogen / adipogenesis / myogenesis) a generic
  loss-of-identity signal shared with other carcinomas (pan-cancer test)?
- ⬜ **Per-sample hub-route OLS explainability** (`gene_expression_explainability.py`, in progress): nested OLS
  on `log2(TPM+1) ~ CN + CPE + PAM50 + promoter meth [+ route miRNA block]`. Early Basal read: CDKN1A ΔR²≈0.16
  (all) / 0.05 (rna_low); PTEN ≈0.11 / 0.05; **IRF1 curated route ≈0.006 but the TargetScan alt (miR-130b/301)
  ≈0.086 / 0.30** ⇒ the IRF1 grant route should use `curated_plus_targetscan` (combined Basal ΔR²≈0.10 / 0.32).
  ⚠ The IRF1 route **failed proliferation adjustment within Basal** (MH-19) and is exploratory until this
  resolves it.
