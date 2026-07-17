# Program forward board — the whole open-work map (`mirna_hallmark`)

> **Goal:** ONE consolidated view of every open forward thread across the whole subproject — the learned-model
> fusion arc, the CN causal machinery, the progression/state axis, the discovery lane, validation, and infra —
> so nothing is lost between the several detail docs. This is the *index of what's next*; the *depth* lives in
> the sync-partner docs.
> **What belongs here:** a one-to-three-line entry per open thread, with a status marker and a pointer to the
> doc that carries its full spec/rationale. NOT the detailed method/rationale (those live in their home docs),
> NOT results (→ VALIDATION/registry), NOT done-and-closed items except as ✅ anchors for context.
> **Update trigger:** when a thread is started/finished, or a new thread is opened anywhere; keep it in sync with
> the two source TODO docs below whenever they change.
> **Sync-partner:** `LEARNED_MODEL_WHATS_NEXT.md` (learned-model detail + design-vs-as-built), `WHATS_NEXT.md`
> (subproject-wide detail), `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` (the ~60-axis fusion map + channel catalog).

Status: ✅ done · 🔨 immediate/flagged · ⬜ open · 🔬 investigate · ⛔ blocked on external · 🚫 deliberate skip.
Where an item has a canonical home, it is named in **(→ doc)**; open that for the full spec.

---

## A. Fusion / Bayes integration — the channel-fused posterior
The converged §6b model generalized to a latent regulatory field fused from many channels. **(→ CHANNEL_FUSION_DESIGN;
its §6–§12 hold the ~60-axis map + channel catalog — this board only lists the live *builds*, not every axis.)**

- ✅ **① Observation distribution → Student-t** (MH-92) — dense mRNA residuals heavy-tailed (df≈7, NB rejected); built
  `_gibbs_posterior(nu=)`; robustify-not-lever. **(→ METHODS §3b, WHATS_NEXT ①)**
- ✅ **Between-family exclusion — Bayes DONE (both halves).** (a) **Magnitude/covariance** (MH-93, `channel_cn.py`): correct
  but WEAK (~0.7%/~0.1σ, mRNA-dominated); gauge fixed (raw-`r`) + gauge-invariant. (b) **Continuous-inclusion / G2** (MH-94,
  `between_family_bayes`): the binary-entry→PIP replacement — CN existence-evidence lifts entry ΔPIP 0.03–0.09 (~10× the
  magnitude effect; existence is coarser) + **Bayesian Shapley with honest width** (PTEN miR-141/200a 0.77±0.41 = genuine
  non-identifiability MH-88's point hid). **General beyond CN** = Bayesian Shapley over the inclusion posterior (applies to
  identity/driver Shapley too). **(→ CHANNEL_FUSION §2.4/§2.5)**
- ✅ **within-family STRUCTURE taxonomy + isomiR QC** (MH-96, `learned/within_family.py`) — looks inside the §8 collapse:
  homogeneous-collapsed / dose-partitioned / seed-heterogeneous (mis-collapse) / all-silent, from expression + 5′-isomiR
  seed-shift (RAW TCGA isomiR, full cohort fetched) + locus; `family_card(gene)` folds in the MH-94 PIP-entry. Collapse
  valid for bulk (median seed-shift 0.038) but miR-17~92/miR-25/92 seed-heterogeneous. ⚠ full-cohort cache rebuild pending.
- ✅ **within-family δ-pooling — DONE (MH-98, `within_family.delta_pooling`); the last Gap-2B Bayes item, CLOSED.** The full
  `M_member = M_fam + δ_m` resolver: **fuses mRNA (member posterior width) + per-segment CN (Ring-1, resolves between-segment)
  + chimeric (within-segment L2)** into one per-member share. PTEN miR-17: fusion surfaces miR-519d (CN-only, mRNA-collinear). =
  `family_card.delta_conf`. Method = **inverse-variance fusion of two INDEPENDENT share estimators** (the correct model); a
  fully-joint hierarchical Gibbs was tested + is WORSE for collinear predictors (mRNA dominates, τ² couples) → no refinement pending.
  (per-segment CN separating different-loci members + chimeric/loading/abundance, NOT the weak family-β slice; wide δ =
  honest non-identifiability for the ~8% |ρ|≥0.95 co-transcribed pairs). Distinct object, may not be as weak. **(→ §6-F)**
- ✅ **② Coding-gene pleiotropy gate — SAFE + PRODUCTIONIZED + WIRED (MH-99, `instrument.coding_pleiotropy` + `exclusion(coding=True)`).**
  RB1@13q14 absorbs 42% of miR-16→CCND1's δ (25% BCL2). Safety = **"trust only reductions"** (not OmniPath — too restrictive, RB1
  is 2-hop). Genome-wide gene CN built (20006×1060, ASCAT×GENCODE); `coding_down_weight`→`s²_π_coding`→`beta_weight_coding` (13q14
  weight 0.5→0.3). Only remaining piece = the proliferation-signal revisit (below). **(→ CHANNEL_FUSION §H/§2.2)**
- ✅ **PROLIFERATION signal — confound vs mechanism CHARACTERIZED + RESOLVED per-gene (MODEL-WIDE, MH-100).** C's
  `mal_prolif` is materially INCOMPLETE (residual-prolif `R=resid(cell-cycle metagene|C)` carries |ρ|≤0.83 after C)
  AND — unlike CN/t/gate-SEs — a real coupling **LEVER** (β moves −20–30%); but it acts three ways at once
  (de-confounds miR-130→ZEB1 −97% / over-controls effectors / un-masks let-7→CDK6 +350%), so a GLOBAL fuller control is
  **REJECTED** (over-controls 50% of a 28-gene panel; the clean fix AL2 nascent-transcription latent is data-blocked).
  A non-circular **OOF 2×2** (train{C,C+R}×eval{C,C+R}) classifies per gene {confound / over_control / fragile}.
  **GENOME-WIDE DONE (8-worker parallel, 1507 genes / 1191 edges, 4.7 min):** 65% neutral (panel was prolif-enriched),
  over_control:confound **2.8:1** among the affected → verdict holds; per-EDGE {mechanism 36% · prolif_robust 20% ·
  deconfound 7% · **unmasked_real 7%** = real edges hidden by proliferation}. **ARM-SOURCE DONE:** the confound localizes
  to ONE arm (median concentration 0.85), systematically **miR-93-5p**. **MECHANISM TESTED (`confound_mechanism`): HOST-GENE
  transcriptional, NOT CN** — miR-93/106b/25 is intronic in **MCM7** (replication gene); host mediates **65–97%** vs locus CN
  **8–10%** ⇒ intragenic-host co-transcription (**axis BP**), handle = host mRNA not the CN instrument (an earlier
  "CN-localizable" guess, refuted). Deliverable = per-gene/edge **FLAG** (`prolif_verdict.py` → `prolif_verdict{,_genomewide,_edges_genomewide}.tsv`,
  `prolif_arm_source.tsv`, `prolif_confound_mechanism.tsv`), C unchanged. Residual open = **host-mRNA (BP) conditioning of the
  confound-classified intronic-oncomiR edges** + independent-cohort check. **(→ CHANNEL_FUSION §9 AL / §6 N / §12 BP; memory `proliferation-confound-vs-mechanism`)**
- ✅ **A-AXIS (fit UNIT) RESOLVED — keep A1 (per-gene); the field does NOT pay (2026-07-12, MH-A-unit).** Probe #1 (`exclusion.tsv`):
  A1 leaves a large weak-instrument territory on the table (512/796 multi-gene segments weak-per-gene). Probe #2 (`hierarchical.oof`,
  the A4 μ_m cross-gene pooling): a WASH even on weak genes (Δρ −0.0006, coin-flip). Since A3 (HMC field) is warranted only if A4
  pays, the cheap probe GATES OUT the field build ⇒ the ⬜ cross-gene / μ_f β-pooling (AF) item below is CLOSED (keep A1). **(→ CHANNEL_FUSION §6.A ★)**
- ⛔ **③ Protein channel (CPTAC) — `βᵗ` FALSIFIED at n=101 (MH-103, Phase-0 validation). (→ `CPTAC_PROTEIN_CHANNEL_PLAN.md §0a`)**
  **The plan's centrepiece is DEAD, and the earlier positive evidence was a LEAK** (the mediator `a` was fitted on ALL
  samples, then `X → protein_resid` OOF-ed; fit `a` INSIDE the fold and the signal vanishes). Leak-free, **four
  independent framings agree**: ridge/all-families λ-swept (median gain **+0.001**); **NNLS all-families z>2 in 0/10**;
  NNLS + top-3-TCGA prior 3/10 (K post-hoc); and — decisively — the **AGGREGATE frame** (predictor = `X_cptac·β_TCGA`,
  weights FIXED from another cohort AND layer ⇒ **1 df, not 60**, the max-power test): **BH q<0.10 in 1/17 genes, only
  BCL2.** **PTEN — the headline throughout — is dead in every framing (aggregate d=−0.006, p=0.82).**
  ⇒ **NOT a leak / estimator / p>n artifact — it is the TRUE EFFECT SIZE meeting the cohort size.** And it matches our
  own biology: **MH-34 already said repression is predominantly mRNA-mediated, translational residual in ~1% of genes
  (11/1132)** — on 17 genes that predicts ~0.2 hits; we got 1. **There was never a βᵗ field to fit. Only n fixes this.**
  **SURVIVING deliverables (none need βᵗ):** ✅ **`locus_cn_cptac`** (BUILT+VALIDATED, median **r=0.997** vs ASCAT;
  coord hole patched 482/482; CPTAC 122×414) · ⬜ **V1 δ-transportability** (first direct test of the cross-cohort
  assumption the program rests on) · ⬜ **Bar-5 systematic** (validates **β** at protein — the `a·β` path, a different
  and EASIER claim; MH-83 shows it works 7/7) · ✅ **`a` propagation slope VALIDATED (MH-105)** — complex subunits buffered (median a_g **0.117** vs **0.437**, p=1.4e-29); ⚠ the miRNA-confound correction is **IMMATERIAL genome-wide** (median bias −0.001). **Engine verdict stands: Gibbs, no HMC.**
  **Headline the axis now delivers:** *the protein layer CONFIRMS the mRNA-mediated model and BOUNDS the translational
  residual below the PER-GENE resolution of a 101-patient cohort.* ⚠ **Precision (verified):** the component is REAL but
  sub-threshold — re-running **MH-34's own partial-correlation frame** on the panel gives **13/17 repression-directed
  (binomial p≈0.025)** yet **BH q<0.10 in 0/17**; MH-34's 11–29/1132 (≈1–2%) predicts 0.2–0.4 hits on 17 genes, so the
  two frames **AGREE**. **MH-34 is NOT overturned — the falsification is of the MODELING OBJECT, not the biology.**
- ✅ **③b PROTEIN AXIS — the SURVIVING deliverables, RUN (MH-104). (→ `CPTAC_PROTEIN_CHANNEL_PLAN.md §0b`)**
  **(A) BAR-5 SYSTEMATIC PASSES (581 genes) — and the loss is at the COHORT boundary, NOT the LAYER boundary.**
  β fit on TCGA mRNA alone (blind to CPTAC and to protein). The control that makes it readable — **fit on TCGA
  half-A → score on TCGA half-B** = the same-cohort out-of-sample CEILING: (i) **−0.117, 88.8% repression-directed,
  p=2e-88** ⇒ the aggregate is sound. Then (ii) CPTAC mRNA **−0.023, p=1.8e-3, retention 0.193**; (iii) **CPTAC
  protein −0.023, 57.7%, p=1.3e-4, retention 0.201.** ⇒ Bar-5 **PASSES** (the honest FDR-controlled successor to
  MH-83's hand-picked 7/7), but ⭐ **the COHORT jump costs ~80% while the mRNA→PROTEIN jump costs ~0%** — which
  **independently re-derives MH-103's βᵗ falsification by a different route** and locates the programme's real
  limit at the cohort boundary (converging with MH-102d's 0.6%).
  **(B) ⚠ RNA-source default CORRECTED → `linkedomics`** (my earlier `star` was WRONG): protein-agreement 0.3713 vs
  0.3678 AND Bar-5 retention 0.193 vs 0.148. Also **explains the "protein beats mRNA" anomaly** — STAR was the
  noisier mRNA; with LinkedOmics protein≈mRNA, as biology requires. (`star` stays: the ONLY CIBERSORTx-usable source.)
  **(C) V1 — DOSE-DELIVERY TRANSPORTS** (the "inherited δ" assumption **VERIFIED, not assumed**): on the 58
  multi-member families expressed in both cohorts, the **same member dominates 84.5% vs 43.7% chance (1.93×)**,
  **member-share correlation 0.991**. Scope = δ's *abundance* input, not its full CN/chimeric fusion.
  **(D) CPTAC CN instrument REAL but WEAK:** F>10 for **59/685 arms (8.6%)** vs TCGA's 32%; **49 arms / only 5
  families** usable in both ⇒ the full CN-fused δ test is low-powered.
  **BUGS FIXED (shared spine):** `data._arm_name_map` global-cache **POISONING** (CPTAC-first silently handed TCGA the
  CPTAC map, 6272→5963; fixed, **verified strict no-op**) · `gauge.beta_table` dropped `.N` arms **ASYMMETRICALLY**
  (TCGA 19 / CPTAC 58 ⇒ the cross-cohort gauge compared DIFFERENT arm sets; fixed → families TCGA 382→402, CPTAC
  376→394). ⚠ **The state session's gauge numbers (0.6%, a=1.125) predate this arm fix and are worth re-running.**
- 🗄 *(historical)* **③ Protein channel — plan as approved 2026-07-12 (superseded by the rows above).**
  **⚠ TWO PRIOR CLAIMS ON THIS ROW ARE NOW CORRECTED (verified, see the plan §1/§3.1):**
  (a) **NOT a "fundamental lever."** Measured under `build_C("cptac")`, protein carries only **≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier ****≈4–6%** (⚠ **CORRECTED, MH-108**: the earlier **1.2%** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**)** used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct Bar-5 retention² check gives ≈4%. **Verdict UNCHANGED — the pre-registered ceiling was ≤7.6% at a_g=1.0.**) of the mRNA channel's Fisher
  information about β (`a_g²·(n_p/n_m)·(σ²_m/σ²_p)` = 0.397²·(101/1041)·0.79) — and **≤7.6% even if `a_g`=1.0** ⇒ it **cannot move β**,
  robust to every plausible `a_g`. It washes like CN (0.7%), and **converges with the state session's independently-derived 0.6%**
  for CPTAC-mRNA→TCGA (MH-102d). **Its real value = a NEW LATENT `βᵗ`** (translational repression, which leaves the transcript
  intact and is therefore invisible to RNA-seq at ANY n). βᵗ verified real **and survives real composition adjustment** (OOF vs
  shuffled null under the composition-C: PTEN z=2.7, PEBP1 z=2.5 — re-deriving MH-35's Class-B survivors multivariately; weak
  hits ZEB1/PDCD4 correctly collapse) but low-power (n=101).
  ⚠ *An earlier version of this row cited `a_g`=0.437 / 1.7% / ≤8.8% — those were computed with **PAM50 in C** (which
  `lineage_verdict` prohibits) and are **VOID**; the figures above were **re-measured, not rescaled**.*
  (b) **The discordance link does NOT force HMC.** Discordance is **additive** (`P = a·M − X·βᵗ`) ⇒ **linear ⇒ Gaussian-conjugate
  ⇒ STAYS ON GIBBS (J1).** Only *saturation* (L2 `a/(a+K)`) breaks conjugacy, and that link already failed its held-out gate.
  **ENGINE VERDICT for the state session: the engine does NOT move. Build against Gibbs.**
  Build = **Model 3** (`P ~ X + M + C`, a cohort-internal mediation model in CPTAC) + **hierarchical β with a LEARNED τ²_cohort**
  (transportability becomes an *estimated* quantity, not an assumption) + **target-gene CN as an instrument for the mediator**.
  **Blocked on:** CIBERSORTx fractions for CPTAC (Wu-major panel). **(→ CHANNEL_FUSION §7, §I1, §L — those sections' conjugacy/engine
  claims are superseded by the plan)**
- ✅ **④ Soft/null-corrected F-weight — DONE (2026-07-12, `soft_fweight_sweep` + `exclusion(shuffle_cn=)`).** `w_soft(c)=(γ²−c·se_a²)₊/s²_π
  = (γ²/s²_π)·(1−c/F)₊`. Swept c∈{5,10,20} on the real table vs a full-universe SHUFFLED-CN null (null collapses strong 17053→141).
  **soft_c=10 is a strictly cleaner gate than the hard F>10 cliff** — same real ranking (ρ=0.995) but 4.4× less null leakage (SNR 43→124),
  by down-weighting the ~11% of weight in the borderline F∈[10,20] band. c=5 too loose, c=20 too strict. **ADOPTED as the DEFAULT
  (`_ADMIT_SOFT_C=10`): `beta_weight = (γ²/s²_π)·(1−10/F)₊`; batch re-run + downstream refreshed.** Immaterial at the channel's ~0.7%
  weight now, but the cleaner gate for when CN is up-weighted (protein fusion). **(→ CHANNEL_FUSION §C)**
- ✅ **Cross-gene pooling — RESOLVED (see A-axis above): μ_m/μ_f β-pooling is a WASH, keep A1.** (The retracted neighborhood fit stays retracted.) **(→ CHANNEL_FUSION §6.A/§AF)**
- ⬜ **Binding channel** (chimeric/K_D) — already a prior (`priors.μ`); the frame just names it a channel. **(→ §7, §I3)**
- ⬜ **State channel** = the progression axis (see C).
- ⬜ **Fusion honesty-layer readouts** (axes U / AO / AP) — cross-channel disagreement provenance vector (U1),
  leave-one-channel-out sensitivity (AO2), explicit flat-direction / identifiability report (AP2 — *names which edges
  the CN channel is even needed for*). Free readouts; materialize once a 2nd non-redundant channel lands. **(→ CHANNEL_FUSION §8 U, §9 AO/AP)**
- ✅ **isomiR Phase 2–4 refit — DONE + GATED OFF (2026-07-12, AS2).** Phase 2 built (`collapse_by_family(isomir=True)`, opt-in) +
  general Phase-4 gate → systematic COUPLING WASH (bidirectional per-edge, real |Δρ| up to 0.14 on jump families but signed mean ≈0)
  ⇒ default-OFF, Phase 3 NOT triggered; deliverable = per-edge re-attribution (attribution/discovery), exactly the expected
  "attribution-yes / coupling-no at n≈1000". **(→ ISOMIR_AWARE_MODELING.md §Phase 2/4, CHANNEL_FUSION §10 AS)**
- 🚫 occupancy/ceRNA saturation substrate (failed held-out gate); any channel whose link isn't identifiable.

## A2. ⚠⚠ CPTAC PROTEIN LAYER — COMPOSITION CONFOUND (MH-107). RE-RUNS PENDING.
- 🔨 **The CPTAC protein validation ran with NO cell-composition term** while **the TCGA model it validates already
  conditions on the 8 Wu-major lineages** — a validation weaker than the model under validation. **Confounding path
  (measured):** EMT proteins ↔ CAF fraction **r=+0.509** (vs +0.033 for all other proteins); **ZEB1 protein ↔ CAF
  r=+0.768**; **miR-200 (an EPITHELIAL miRNA) ↔ CAF r=−0.35** ⇒ "miR-200 represses EMT protein" is **compartment
  arithmetic**. **Graded severity:** MH-34 gene/hallmark **CATASTROPHIC (FDR-neg 66→3, −95%; EMT −0.123→−0.009)** ·
  MH-84 per-edge **MODERATE (73% keep sign, 27% SIGN-FLIP; ~1/4 of the "268 triple-validated" are artifacts)** ·
  MH-83 hub **SURVIVES 7/7** (but ZEB1 −75%; BCL2 strengthens).
  **⭐ REFRAME, NOT RETRACTION:** stroma-loaded programs collapse (EMT/COAGULATION/ALLOGRAFT) while **tumour-INTRINSIC
  programs survive or strengthen** (ESTROGEN_RESPONSE, PI3K_AKT_MTOR, BCL2, CDKN1B, ESR1) ⇒ the protein layer validates
  cell-autonomous miRNA regulation **in the EPITHELIAL compartment**.
  **⭐ CONFOUND not MEDIATOR (adjudicated):** in TCGA (whose C already has the lineages) miR-200→ZEB1 mRNA **survives
  composition at n=1041, ρ=−0.209, p=8.7e-12** ⇒ conditioning does NOT kill the edge. **The edge is REAL** (true effect
  bounded 0.209 ≤ |ρ| ≤ 0.446); only the bulk-protein magnitude was confounded.
  **✅ CODE FIXED** (`cptac_validation._covariates` 2→11 covariates, RAISES rather than running confounded;
  `ood_protein`, `dossier.tier3_protein` — all were RAW Spearman). **⬜ RE-RUNS PENDING: MH-34/35/36/39/40 · MH-83 · MH-84 tier-3.**
  **(→ registry MH-107, memory `cptac-protein-composition-confound`)**

## B. CN instrument — within / between family
- ✅ Gap 1/2A/Ring-1 (MH-87), between-family frequentist interim (MH-88), T1 exclusion gate (MH-90), channel MVP +
  weak-but-diagnostic verdict (MH-91). **(→ CN_INSTRUMENT §7/§9/§10, memory `cn-instrument-gap2-arc`)**
- ✅ within-family δ-pooling + cross-family covariance = the A Bayes-unification item **DONE** (MH-93/94/98; δ-pooling is a
  separate/orthogonal dose-delivery latent, inverse-variance not joint-Gibbs — the joint form tested-worse).
- ✅ **coding-gene pleiotropy = A② — DONE + extended (MH-99/101):** per-locus HOST fold-in (`locus_host_map`→`locus_context.tsv`),
  the host-family `do_scan` restriction LIFTED → full-universe scan (3059 pure-coding + 36 host down-weights), and the 26 ambiguous
  `repression_direction` down-weights GRADED (`grade_host_downweights`). Coding scan 17× faster. **(→ CHANNEL_FUSION §H)**
- ✅ **within-family arm attribution** (Ring-1 L2 dose-delivery) — the frequentist point the δ-pooling Bayes-ifies (MH-98).
- ✅ **BP intragenic-host pleiotropy — BUILT (MH-101 genomic-context axis, per-locus host instrument validity).** Remaining mapped-but-unbuilt
  instrument threats: 🔬 BL subclonal/CCF-weighted CN, AD3 allele-specific CN on imprinted loci. **(→ CHANNEL_FUSION §12 BP/BL, §AD)**
- 🔬 **AI2 second instrument — germline-eQTL / methylation-as-instrument** (over-identifies the same β, breaks
  reverse-causation differently than somatic CN; a stronger Sargan/Hansen jointly with CN). Mapped, unbuilt; needs ancestry-PC C (BQ). **(→ CHANNEL_FUSION §9 AI)**

## C. Progression / state axis — ⛔ **CLOSED (MH-102d, 2026-07-12). The channel was MEASURED and CANCELLED.**
> ⚠ **numbers corrected (scale-dependent sd-floor bug, found by the CPTAC session):** **FINAL (both sd-floors scale-free): GTEx a=0.137 ρ=0.188 τ=0.0009 info=0.6% · NAT a=0.326 ρ=0.277 τ=0.0221 info=0.7% · CPTAC a=0.333 ρ=0.280 τ=0.0023 info=0.7%.** Conclusion UNCHANGED (τ≈0, info ~0.6–0.7%) ⇒ state channel stays CANCELLED. 
- ⛔ **STATE CHANNEL: NOT BUILT.** On a **calibrated** gauge (split-half control a=1.125, cal✓) + the E3 standardized-β
  convention: **τ (the per-edge acquired/lost-wiring payload — the axis's ENTIRE deliverable) is indistinguishable from
  ZERO**, and the coupling contribution is **2.4%** of the mRNA likelihood's own precision (a same-cohort half gives 124%).
  Structural: precision ∝ a², and a≈0.11 ⇒ even infinite GTEx donors add ~1%. **Composition attenuation doesn't just bias
  the healthy β — it destroys its information.** ⚠ **Applies to the PROTEIN channel too:** CPTAC-mRNA→TCGA is the *easiest*
  cross-cohort transfer (same modality, same state) and caps at **0.6%** ⇒ a protein channel must justify itself as a
  DIFFERENT QUESTION, never as a coupling gain. **(→ registry MH-102d, STATE_CHANNEL_PLAN §10)**
- ✅ **SHARED CROSS-COHORT INFRASTRUCTURE — BUILT (MH-102, 2026-07-12). Both channels (state AND protein) consume this.**
  (a) **`learned/confounders.py`** — `build_C(cohort, participants)` for **four** cohorts (TCGA·NAT·GTEx·CPTAC), verified
  drop-in (max|diff|=0.0 vs `data._deconv()` *and* `states._cibersortx_state_cov()`). (b) **`learned/gauge.py`** — the
  cross-cohort gauge `β_source = a·β_target + c + δ`, out-of-fold `a`, `to_channel_terms()` (Gaussian-conjugate ⇒ **J1**),
  the shared β producer `beta_table()`, and **`shuffled_gauge()` = the cross-cohort SHUFFLE-NULL primitive**. **Ladder (all
  PASS their null): CPTAC a=2.369 · NAT a=0.528 · GTEx a=0.231**, factorizing as state ≈2× × platform ≈2.3×. ⚠ **a≠1 even
  for CPTAC (same modality, same disease state) ⇒ the protein channel CANNOT feed β_protein in raw.** (c)
  **`learned/immune_completeness.py`** — C is immune-*incomplete* but the fix is *inert* (β ρ=0.967) ⇒ **C unchanged**.
  **(→ registry MH-102/102b/102c, STATE_CHANNEL_PLAN §4/§6)**
- ⬜ **State channel `M_T = a·M_H + Δ`** (Decision-H *state* flavour, axis Q) — GTEx-breast M_H as prior mean for tumour;
  the "acquired vs constitutive wiring" readout. **PLANNED; the panel question is SETTLED (2026-07-12). (→ `LEARNED_MODEL_STATE_CHANNEL_PLAN.md`)**
  ⚠ **The additive `M_T = M_H + Δ` form is WRONG** — GTEx is ~32% epithelial vs an epithelial-dominant tumour ⇒ β_H is
  composition-**attenuated** ⇒ Δ>0 on nearly every edge = a *manufactured* "acquired wiring" artifact. A **global gauge
  constant `a`** absorbs it (composition + TPM↔RPM platform + C-block differences); the **per-edge deviation Δ** is the
  deliverable. ⚠ **"Reference reconciliation" is NOT the open piece** (see F) — the real cross-state problem is the
  **matched-C policy**: `CPE` / `target_cn` / `mal_prolif` are all **TUMOUR-ONLY**, so the C blocks cannot be identical
  across states whatever panel is chosen. **(→ CHANNEL_FUSION §Q, STATE_CHANNEL_PLAN §2/§4)**
- ⬜ **Abundance-vs-wiring decomposition** (Δx·M vs x·ΔM) as a posterior readout; **subtype interaction tests
  (Gelman–Stern)**. **(→ WHATS_NEXT REMAINING PHASES)**
- ⬜ **State-stratified discovery** (per-state permutation-FDR: found-in-tumour-not-NAT = acquired, all-states =
  constitutive, only-NAT = field/lost); power-bounded (NAT n=104 / GTEx n=327). **(→ WHATS_NEXT EXPAND)**
- ⬜ **GTEx→NAT→tumour trajectory** (`parallel_view.change_trajectory`) + the **GTEx logFC gauge check** (rawFC_GN vs
  relFC_GN — validate the miRNA RPM≈TPM assumption). **(→ WHATS_NEXT PARALLEL-VIEW)**

## D. Discovery lane
- 🔨 **Per-edge (+ per-gene) attribution card** — one row/edge joining regime context (range stats, spiker flag),
  budget share (E7/G4), composition tag (deconv retention), shift-class; per-gene = the G-series roll-up. **(→ WHATS_NEXT IMMEDIATE)**
- 🔬 **Top-discovery deep-dives:** the **miR-17~92 / 106b-25 cluster** (RABEP1/AHNAK/WWP1/…), **miR-135b→GATA3** (spiker,
  28%>floor), the **568 fully-novel candidates** (triage + lit-novelty + Bar-5/Manakov shortlist). **(→ WHATS_NEXT INVESTIGATE)**
- ⬜ **Epigenetic-rewiring discovery** (methylation lost/gained specialists, Axis A de-repression / Axis B acquired-
  repression) — detector built (`structural_identity`+`methylation`), needs genome-wide scale-up; confounders =
  purity/composition, CN-vs-methylation, imprinted loci (DLK1-DIO3); prereq = UUID→TCGA-barcode map. **(→ WHATS_NEXT EPIGENETIC-REWIRING)**
- ⬜ **Target-side SNV / APA / A-to-I editing** mechanism legs (reference-blind edges: seed gain/loss, 3′UTR exposure,
  seed editing) — the identity/structural quadrant as a mechanism router. **(→ WHATS_NEXT PARALLEL-VIEW SNV leg)**
- 🔬 **context++ × K_D fusion** (orthogonal, ρ 0.065 → a fusion may beat either) + complete ~2,600-arm K_D reference
  scan + discovery floor sweep (RPM 50/10/1). **(→ WHATS_NEXT K_D)**
- ⬜ **CLIP-source candidate expansion** (ENCORI/POSTAR3 physical sites) + **regulation-pattern archetypes** (parallel-
  view pairwise-ranking fingerprint → abundance-/weight-/designed-but-silenced/redundant-collinear classes) +
  **seed→targetome-rarity** sequence-specificity ground. **(→ WHATS_NEXT PARALLEL-VIEW)**

## E. Validation gaps
- 🔬 **CPTAC protein — SYSTEMATIC (Bar 5)** — docs flag this **the highest-value item left**; only miR-135b→GATA3 done.
  Infra: `eval/cptac_validation.py`. **(→ WHATS_NEXT VALIDATION GAPS)**
- 🔬 **External-cohort replication (METABRIC / Buffa)** of the learned M / top edges (Buffa infra exists; METABRIC-miRNA
  EGA-gated, partial subset) + Manakov chimeric validation of the shortlist.
- ✅ **POSTERIOR CALIBRATION — MEASURED + FIXED-BY-REPORTING (MH-102e, `learned/calibration.py`). ⚠ SBC was the WRONG
  TOOL and that board item is RETIRED.** The model's reported uncertainties **understate the true sampling variability**
  (measured against independent half-cohorts, 3 splits): **bagged NNLS 0.73× (inflate 1.37×) · Gibbs Gaussian 0.77× (1.29×) · Gibbs Student-t ν=7 **0.89×** (1.13×)** ⇒ **posterior widths are ~25% too narrow** *(an earlier "35–70%" was measured on a biased gene subset — RETRACTED)*.
  **MH-92's Student-t SUBSTANTIALLY FIXES it (0.77→0.89, nearly calibrated) — turning `nu=7` on is a LIVE recommendation** *(the earlier "does not fix it" is RETRACTED)* ⇒ heavy residual tails are part of the story,
  not all of it; the residue is **unmodeled between-participant heterogeneity** (the likelihood treats patients as iid
  given X and C — they are not). **⇒ SBC would PASS and mislead**: it asks "is the sampler right *given the model*?", and
  the fault is the MODEL, not the sampler. The honest instrument is resampling (`posterior_calibration()`), and the fix is
  to **multiply reported widths by `inflation_factor`**. **TRIAGE — who is affected:** width claims (**MH-94** Bayesian-
  Shapley "honest width", **MH-98** δ-pooling confidence) are **understated** — note the direction: wider intervals make
  MH-94's non-identifiability claim **STRONGER**; inverse-variance **fusion weights are SAFE** (they use SE *ratios*, which
  cancel); channel `s²` immaterial (~0.7% weight); and **all OOF/permutation results are UNTOUCHED** (`prolif_verdict`'s
  2×2, discovery FDR, coupling ρ) — which covers most of the program's headline findings. **(→ registry MH-102e)**
- ✅ **β GAUGE CONVENTION — DECIDED: keep raw-`r` for the model; z-Y is GAUGE-ONLY (2026-07-12).** Tested: within-gene
  family rankings are **exactly invariant** (ρ=1.0000) to z-scoring Y, and the cross-gene ranking is near-invariant
  (ρ=0.992) — but raw-`r`'s discovery ranking is **more reproducible** across independent halves (**63/100 vs 55/100**
  top-100), because z-Y divides β by `sd(resid Y)` (a **7.2×** spread) which is smallest exactly where the residual is
  most noise-dominated. ⚠ Honest caveat: raw-`r`'s top-50 **is** tilted toward high-variance genes (median sd(residY)
  0.369 vs 0.243 overall); z-Y removes that tilt but costs stability. ⚠ **The cross-cohort GAUGE must still use z-Y**
  (`gauge.beta_table(zscore_y=True)`) — there the cohorts' Y-scales differ *systematically* (TCGA 0.237 vs GTEx 0.612)
  and the factor lands entirely in `a`. ⚠ **Discovery-lane reality check:** top-100 edge reproducibility across
  half-cohorts is only **~63%** under EITHER convention. **(→ registry MH-102d)**

## F. Design-vs-as-built partials + composition/infra
- ⬜ **B transcription-rate proxy** (intronic/pre-mRNA nascent latent) — the design's "most load-bearing"; HIGH value but
  RISKY (over-control, the TF-proxy failure mode) + needs intronic-read data. **(→ WHATS_NEXT DESIGN-vs-AS-BUILT B)**
- ✅ **E direction channel** (TarBase-v9, 6 activating edges flagged); ✅ **Decision I Shapley** identity-vs-magnitude
  (`attribution.identity_vs_magnitude`); ✅ **H subtype** flavour. A/D/F-pooling **subsumed** (verified no gain).
- ⬜ **Phase-4 cooperativity / ceRNA (Decision G)** — sits on the failed occupancy substrate; gate HARD on held-out.
  `priors.py` spike-and-slab π/μ/τ object (LOW payoff, lasso-FDR discovery already works).
- ✅ **Canonical M** unified + speed pass; ⬜ **remaining:** one canonical rerun of the persisted attribution outputs
  (`output/learned/programs/*`, cards) now the estimator is settled; optional canonicalize the `shift_vs_weight`
  per-edge weight diagnostic. **(→ WHATS_NEXT CANONICAL M)**
- ✅ **CIBERSORTx — SETTLED (2026-07-12). ZERO new runs needed; the "reference reconciliation" blocker is RETIRED.**
  The panel choice **does not move β**: on 1444 (gene,family) coefficients fit on the 327 paired GTEx donors, β_H is
  ρ=**0.94** concordant across *entirely different atlases* (Wu-9 vs HBCA-11) and ρ=**0.99** across reference
  cell-subsamples (300 vs 150 cells) — while **no-C vs any-C is ρ=0.58 and *halves* |β|**. ⇒ *having* a C is ~10× more
  important than *which*. Wu-9 stays the panel for all three states (HBCA cannot supply Cancer-Epithelial, which
  `mal_prolif` is built from; and swapping `_DECONV_COLS` ripples through every persisted output).
  ✅ **BUILT: `learned/confounders.py`** — `build_C(cohort, participants)` + `availability()`, ONE builder for **four**
  cohorts (tcga · nat · gtex · **cptac**), verified **drop-in identical** to both existing consumers (`data._deconv()` and
  `states._cibersortx_state_cov()`, max|diff|=0.0) ⇒ TCGA/NAT migration is a **no-op**; only GTEx (had **NO C at all** in
  `state.cross_state_transfer`) and CPTAC (had **no composition at all**) change. Sources: tumour+NAT = the pooled run's
  tumour/`-NAT` rows · GTEx = `gtex_wu_major` (514, **327/327** paired donors) · CPTAC = `cptac_wu_major` (133, **104/104**
  of the prospective miRNA cohort — also makes `mal_prolif` computable there). Purity stays **cohort-native** (CPE /
  ESTIMATE / none) — it is nearly redundant given composition (β ρ=0.984) and NAT/GTEx must have **none** (a
  participant-keyed join would hand 107/113 NAT samples their own *tumour's* CPE).
  **Genuine gaps (not deconvolution):** (i) **adipose** — *no* reference has an adipocyte class (Wu-9, HBCA-11, LM22) and
  GTEx breast is adipose-dominant ⇒ needs a marker metagene; (ii) **immune C-completeness** — Wu-9 lymphoid tracks Thorsson's
  infiltration score at r=0.31 vs **LM22's 0.70** ⇒ test in tumour (LM22 + Thorsson are the *validators*, never covariates —
  they exist only for tumour); shared with the CPTAC session.
  **Do NOT run:** a fine 42–49-way panel (β is 0.94–0.99 invariant to reference detail ⇒ a finer C cannot help), TCGA×HBCA,
  or refsample "harmonisation" (the run FILENAMES lie — by md5 there are only TWO Wu references; the sole production
  mismatch is CPTAC-on-REF-B, worth β ρ=0.986, absorbed by the gauge + s² floor; and S-mode re-adjusts per mixture anyway).
  🚫 Retracted: the "degenerate compartments / unmeasured PVL+lymphoid" reading (kappa 8.1/12.0, PVL's technical
  reproducibility 0.998 — the Wu-vs-HBCA anti-correlation is an *annotation* difference); and "use `tcga_nat_alone`"
  (β ρ=0.957, but its malignant negative control degrades 3× at the expense of `Normal Epithelial` — use the pooled rows).
  **(→ STATE_CHANNEL_PLAN §1/§3/§4/§9)**
- ⬜ **Infra:** `git add` the untracked `learned/` tree (provenance hole, flagged repeatedly); doc consolidation
  (~13 learned-model docs → ~6); fold ENCORI/POSTAR3/Manakov into the ledger (union, not summed); `baselines/`
  re-export shims; proper covariate-protected ComBat batch (only if a channel reintroduces batch). **(→ WHATS_NEXT COMPOSITION/INFRA)**

---

## Top of the board (what the priority signals say — updated 2026-07-12)
The fusion arc's causal core is now COMPLETE (Bayes unification MH-93/94/98, coding+host gate MH-99/101, A-axis resolved to A1,
isomiR gated). Every INTERNAL lever has landed "immaterial-at-n≈1000" ⇒ the frontier is now exogenous, not more internal refinement.
1. **CPTAC protein channel (A③ / E, Bar-5)** — THE sole remaining fundamental lever; docs' "highest-value item left" AND the first
   channel to force HMC (axis J). A non-redundant 2nd likelihood the mRNA floor can't absorb.
2. **State/progression channel (C)** — the other exogenous axis. **No longer gated on CIBERSORTx** (F: settled, zero runs);
   planned in `LEARNED_MODEL_STATE_CHANNEL_PLAN.md`. First work = wire the three C blocks (GTEx has *no* C today) + the
   matched-C policy, not a deconvolution.
3. **Per-edge/per-gene attribution card (D)** — the flagged IMMEDIATE deliverable, pure join over what's computed.
4. **Discovery deep-dives + external replication (D/E)** — miR-17~92 cluster, 568 novel candidates, METABRIC/Manakov.
5. **Infra debt (F)** — `git add` the untracked `learned/` tree (provenance hole); doc consolidation.
   *(Done & retired from priority: Bayes unification, coding/host pleiotropy gate + per-locus + full scan + grading, A-axis / cross-gene pooling, isomiR refit.)*
