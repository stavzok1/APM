# Learned miRNA→gene regulatory model — design report

> **⚠ HISTORICAL — DESIGN-PHASE (2026-07-04), SUPERSEDED. Do NOT cite this as the current model.** The model has
> since CONVERGED: **one dense learned-τ² Bayesian posterior, two readouts, adaptive lasso RETIRED** — see
> **[`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md §6b`](../LEARNED_MODEL_DISCOVERY_SYNTHESIS.md)**. In particular this doc's
> "lasso as one baseline among several", the Decision A–J log, and the 5-Bar framing predate that convergence;
> where they say "the model is X" or "X has no payoff", read §6b instead. Kept for provenance/reasoning only.

**Responding to** `LEARNED_MODEL_DISCUSSION_PROMPT.md` (§9 deliverable).
**Scope.** How to *learn* the miRNA-arm → Hallmark-gene regulatory weights from TCGA-BRCA
expression under multi-source evidence priors, how they aggregate/attribute, and how to know they
capture real targeting. The current `softmax`/`D(m)`/`edge_w` machinery is treated as **one baseline
among several**, per §0.

**Status.** Design / not yet built. This is the reasoning-first spec to audit and act on. Companion
to `FORMULAS.md` (coded heuristic), `LEARNED_REGULATORY_MATRIX_DESIGN.md` (earlier RRR sketch, which
this supersedes and extends), `DATA_SOURCES.md`, `MODELING_FRAMEWORK.md`,
`EDGE_QUESTION_TAXONOMY.md`, `GENE_QUESTION_TAXONOMY.md`.

**Legibility convention.** I mark **[assert]** where a claim is established by literature/data cited
here, **[infer]** where I am reasoning beyond direct evidence, and **[choice]** where it is a design
preference that could go the other way. Every major decision carries an *Alternatives kept in play*
block — I prune only what is genuinely dominated.

---

## 1. Executive summary

**The model I would build.** A **program-wise (per-Hallmark) hierarchical Bayesian
occupancy-regression**. For each of the 50 Hallmark sets, learn a non-negative repression weight for
each candidate (arm → gene) edge from tumour mRNA under a confounder block, with:

- an **occupancy (saturating) link** for how arm abundance turns into repression — not a linear
  `a·w` — so the a·w gauge, the RISC-loading threshold, and ceRNA/pool competition are represented by
  the *functional form* rather than by a hand-set denominator;
- a **spike-and-slab, evidence-graded inclusion prior** `π_{g,m}` plus an informative **magnitude
  prior mean** `μ_{g,m}` and **prior strength** `τ_{g,m}`, fused from a re-graded evidence stack
  (curated MTI tiers **+ TarBase-v9 method/confidence + a single PMID-deduped AGO-CLIP axis
  (ENCORI∪POSTAR3∪TarBase-CLIP; primary) + chimeric duplex coordinates (confirmatory) + TargetScan
  context++ / biochemical K_D + breast literature**), with **site-level** evidence setting magnitude
  and interaction support and **edge-level** evidence setting inclusion and sign;
- **seed-family partial pooling** so the identified estimand is **family→gene**, with the per-arm
  split reported *with its posterior width* rather than faked;
- **primed cooperative product terms** only for arm-pairs with distinct sites 8–40 nt apart, and a
  **shared free-AGO pool** as the mechanism for global saturation / ceRNA;
- a **nested state model** `M_tumour = M_healthy + Δ` (healthy GTEx fit as the prior for the tumour
  fit) that separates *acquired abundance* from *acquired per-site weight*.

**Observation target `Y`.** **Primary = tumour mRNA level with an explicit transcription-control
block in `C`** — a **minimal, composition-robust** set: purity = **CPE/ESTIMATE**, target-CN, **one
malignant-compartment (Cancer-Epithelial CIBERSORTx) proliferation**, a transcription proxy, plus
**deconvolution compartment fractions** for stroma/immune programs; HRD and batch off by default
(Decision B) — because mRNA carries the majority of miRNA repression at steady state **[assert: Guo 2010;
Eichhorn 2014]** and is where the curated evidence and the power are. **Validation `Y` = CPTAC protein and
protein-vs-mRNA discordance** (the translational residual the mRNA fit cannot see) plus an independent
mRNA cohort (SCAN-B / METABRIC / Buffa).

**Prior scheme.** Non-circular by construction: `π, μ, τ` and the interaction support come from
**sequence + curation + external-tissue assays only**, never from TCGA `(X,Y)`; expression estimates
*magnitudes*; the learned-M and curated-M are kept side by side and disagreements are reported as
findings.

**Evaluation.** A falsification ladder with pre-named criteria, scored against the *right* reference
per claim: (1) beat **raw abundance** out-of-fold on held-out patients — the primary bar; (2) match or
extend **curated fixed-M** while nesting it; (3) beat a **shuffled-evidence null**; (4) a **miRNA-locus
CN instrument** (2SLS) for causality; (5) an **orthogonal protein layer** (CPTAC) and **independent
cohort** for OOD. Abandon the learned model for raw abundance if it cannot clear Bar 1.

**The 3–5 decisions that matter most.**
1. **Observation model + a minimal, composition-robust `C`** — regress mRNA with a *transcription-control
   block*, not naked level; validate on protein discordance. Keep `C` lean: purity, target-CN, **one
   CIBERSORTx-deconvolved (malignant-compartment) proliferation** term, transcription proxy; drop HRD and
   batch unless a diagnostic justifies them (over-adjustment/collider risk). (Everything downstream is
   biased without the transcription control; it is the biggest hole in the current heuristic.)
2. **Occupancy link + abundance-threshold gauge** — fix abundance on the `log2(RPM+1)` scale, put a
   saturating gate with a functional-copy-number threshold, let `M` carry per-site affinity. This
   resolves the a·w non-identifiability and *retires `D(m)`* (promiscuity becomes a consequence of a
   shared finite RISC pool, not an input).
3. **Evidence fusion — de-duplicate by PMID *globally* (functional AND physical), which is not done today.**
   The current score counts miRTarBase studies + *adds* ENCORI CLIP on top with **no cross-database dedup**;
   the same paper double-counts both as CLIP (ENCORI/POSTAR3/TarBase — TarBase is ~93% CLIP from just 998
   PMIDs) **and** as functional (a luciferase/western paper curated by both miRTarBase and TarBase). Fix =
   a unified **`(m,g,PMID,assay-class)` ledger built at edge-build time**, counting *distinct PMIDs* per
   class. TarBase's `regulation` is ~all-Negative (301 Positive in 4.72M) → **not** a direction prior.
   Sequence/K_D is a separate orthogonal (PMID-free) magnitude channel; K_D gives the **orphan** magnitude.
4. **Family partial pooling** — the identified unit is the seed family; per-arm resolution is a
   nomination with uncertainty, not a coefficient.
5. **Nested state model** — test whether `M` itself moves H→tumour (APA/AGO capacity) versus only `x`
   moving; this is falsifiable and it decides §7.

---

## 2. Literature synthesis (refreshed from the web, July 2026)

What the current state of the art establishes across the §1j facets, with what I **adopt** vs
**improve on**, and contradictions flagged. Citations are listed in §8.

### 2.1 Silencing mechanism & magnitude
- **mRNA destabilisation dominates at steady state.** By the time substantial repression ensues,
  mRNA decay explains ~66–90% of miRNA-mediated repression; translational repression is rapid but
  weak (~10–25% residual) **[assert: Guo et al. 2010; Eichhorn et al. 2014]**. Single-cell work shows a
  hierarchical target response to miRNA induction (Kim/Bartel line) **[assert]**.
  → **Adopt:** mRNA is a legitimate **primary** readout. → **Improve on:** the translational-only
  residual is invisible in mRNA and *only* appears in protein / protein-vs-mRNA discordance — so
  protein is not a nice-to-have, it is the layer that separates a real post-transcriptional edge from a
  transcriptional confound (§Decision B, §Decision J).
- **Effector biology.** AGO→TNRC6/GW182→CCR4-NOT/PAN2-3 deadenylation→decapping→decay, plus initiation
  block **[assert]**. → I do **not** model the effector chain explicitly; its net is the saturating
  decay-rate multiplier in the occupancy link (McGeary's `b≈1.8`, "bound AGO nearly triples decay").
- **Fine-tuning / mild effects.** Most targets < 2-fold; buffering, threshold-sharpening; one miRNA
  hits hundreds of proteins, each mildly **[assert: Selbach 2008; Baek 2008]**. → **Adopt:** low SNR
  per edge ⇒ borrowing strength (priors, groups, low rank) is mandatory, not optional.

### 2.2 Quantitative site efficacy & 3′UTR biology
- **Site-type hierarchy** 8mer > 7mer-m8 > 7mer-A1 > 6mer, plus 3′-supplementary pairing **[assert]**.
- **TargetScan context++** sums ~14 features (site type, supplementary pairing, local AU, minimum
  distance, 3′UTR length, SPS, TA, PCT, ORF length/8mer, offset-6mer) into a continuous per-site score
  **[assert: Agarwal et al. 2015]**. → **Adopt** as a per-edge continuous magnitude prior available for
  **every** predicted edge, orphans included.
- **The biochemical / K_D model is the current frontier and it changes my recommendation.**
  AGO-RBNS measures relative K_D for all ≤12-nt sequences; a Langmuir **occupancy** model
  `N ∝ a/(a+K_D)` predicts repression with r²≈0.30–0.37 (vs ≤0.14 for 30 prior algorithms), and a CNN
  generalises K_D to arbitrary miRNAs from the 5′ sequence (~+50% over TargetScan on held-out miRNAs)
  **[assert: McGeary, Bartel et al. 2019; scanMiR / Soutschek 2022 packages it]**. Two features I lean
  on hard: (i) **flanking-dinucleotide AU context spans a ~100-fold affinity range** — often larger
  than the site-type difference — so "the site exists" is a weak prior; the *graded* K_D is the strong
  one; (ii) repression scales with **occupancy**, a saturating function of **free** AGO-miRNA
  concentration and K_D. → **Adopt** the occupancy link as the functional form (§Decision C) and K_D as
  the mechanistic magnitude prior (§Decision E). → **Improve on:** McGeary fits the free-AGO scalar
  `a_g` *per transfection* and holds it fixed across cell types — i.e. the model **openly does not
  resolve the relative-vs-absolute abundance gauge**; that unresolved scalar *is* our §3f gauge problem,
  and I anchor it with the functional threshold rather than leaving it free.
- **Positional / accessibility / APA.** Functional sites enrich near UTR ends, deplete ~15 nt after the
  stop and in long-UTR interiors; open (accessible) structure predicts efficacy; APA shortens 3′UTRs in
  proliferating/tumour cells and deletes distal sites **[assert]**. → **Adopt:** APA makes **support
  itself sample/state-dependent** — a direct argument for the state-dependent `M` (§Decision H) and for
  conditioning site presence on 3′UTR isoform where PDUI is available.
- **Deep-learning site predictors** (miRAW, DeepMirTar, miTAR, TargetNet, TEC-miTarget; and 2026
  multi-instance models e.g. PAIR-Former) reach high *classification* AUC on curated positives/negatives
  **[assert]**. → I **do not** adopt these as the estimator: they predict *binary functionality* on
  curated sets, tend to overfit the negatives, and do not give a calibrated *magnitude* on the natural
  edge distribution. I adopt only the **biochemical/K_D** line, which is quantitatively grounded and
  regression-calibrated.

### 2.3 Seed families, cooperativity & competition
- **Cooperativity window 8–40 nt** between *distinct* sites gives super-additive repression, strongest
  ~15–26 nt, lost beyond ~40 nt **[assert: Saetrom 2007; Grimson 2007; Broderick 2011]**; the human
  transcriptome is enriched for cooperativity-permitting spacings **[assert: Rinck/Sætrom-line 2013]**.
  → **Adopt** as the trigger for *primed* product terms (§Decision G).
- **Seed-family redundancy** — same-seed members are near-substitutable at a shared site (sub-additive),
  the opposite sign to cooperativity **[assert]**. → **Adopt** as the biological mandate for the family
  group structure (§Decision F).
- **ceRNA is real but narrow.** Physiological competition requires near-equimolarity / low target
  abundance; Denzler showed target relief only above ~1.5×10⁵ added sites/cell, insensitive to miRNA
  level; Bosson (AGO iCLIP) found most active miRNAs are *not* ceRNA-susceptible, with low-abundance
  families (miR-25) the exception **[assert: Denzler 2014/2016; Bosson 2014]**. → **Adopt:** model
  ceRNA/global saturation as a **cohort/sample-level shared-pool** effect (rare, global), **not** a
  per-edge ceRNA term — the literature says a coarse resolution is the honest one.
- **Functional dose threshold.** Below a per-cell copy floor a miRNA is inert regardless of site
  quality; a few most-abundant miRNAs dominate AGO loading **[assert: Mullokandov 2012]**. → **Adopt**
  as the anchor that *breaks the a·w gauge degeneracy* (§Decision C).

### 2.4 Expression-based regulatory modelling
- **GenMiR++** (Bayesian, sequence-primed, marginal associations) and **TaLasso** (penalised
  regression, *conditional/partial* associations) are the canonical expression-network methods; sparse
  network-regularised NMF (**SNMNMF**) and **PIMiM** find miRNA-gene **comodules**; Bayesian **sparse
  reduced-rank regression** simultaneously does rank reduction + predictor/response selection in `p≫n`
  omics **[assert]**. → **Adopt** the conditional-not-marginal stance (TaLasso), the comodule/low-rank
  idea (as an *overlay*, §Decision A), and Bayesian sparse-RRR machinery for the universe layer.
  → **Improve on (the differentiators):** none of these push a **quality-weighted, assay-directness-
  graded, breast-audited magnitude** prior into penalty+mean (they use sequence only for *support*);
  none carry **purity/proliferation/target-CN/transcription** confounder control (composition
  confounding is the #1 failure mode); none carry an **internal instrument**; none model
  **cooperativity as primed products**. These four gaps are exactly what our data lets us close.
- **MoPC** (partial correlation of mRNA–protein *conditioned on miRNA*) infers miRNA–target
  interactions from TCGA/CPTAC multi-omics **[assert, 2024]**. → **Adopt its logic as a validation
  estimator** on the protein layer (it is precisely the protein-discordance test).

### 2.5 Causal / identification
- **Invariant Causal Prediction using PAM50 subtypes as "environments"** finds miRNA→mRNA edges whose
  relationship is *invariant across subtypes*, arguing this is robust to TF confounding; validated
  against transfection (|log2FC|>0.3); it explicitly **misses subtype-specific edges** **[assert: Meng
  et al. 2019]**. → **Adopt** the environments idea for the state/subtype axis (§Decision H) and the
  invariance-vs-specificity tradeoff it names. → **Improve on:** ICP tests invariance but doesn't
  estimate a *magnitude* under an evidence prior or carry a dose instrument.
- **miRNA Mendelian randomisation** studies use germline *eQTLs* of circulating miRNAs as instruments
  for cancer risk **[assert, many 2024–2025]**. → I flag this as a **different design** from ours: those
  are germline-eQTL, cross-individual, risk-outcome MRs. Our instrument is **somatic miRNA-locus CN
  perturbing arm dose within tumours** — a within-cohort dose instrument for the *target*, not a
  risk MR. The MR literature validates the *logic* (genetic dose instrument breaks reverse causation);
  the *estimand* is ours to define.
- **Composition/purity confounding** in tumour omics is well documented; correlation ≠ targeting
  without adjustment **[assert]**. → **Adopt** into the design `C` block.

### 2.6 Protein-level & discordance
- pSILAC established mild, widespread protein effects **[assert: Selbach/Baek 2008]**; CPTAC breast
  proteogenomics shows real mRNA–protein discordance (e.g. TP53 r≈0.36; post-transcriptional control of
  specific genes) **[assert]**. → **Adopt:** protein and discordance are the orthogonal OOD layer; CPTAC
  phospho/acetyl (now inventoried) give *effector* readouts for select programs.

### 2.7 State / field-effect biology
- Histologically normal tumour-adjacent breast tissue carries field-cancerisation defects (DNA-methylation
  drift, CNAs in ~40% of adjacent tissue, altered expression) **[assert]**. → **Adopt the brief's stance
  as literature-backed:** **NAT is its own state, never a healthy stand-in**; GTEx is the only bona-fide
  healthy anchor. Using NAT as normal understates tumour change and can invert directions.

### 2.8 Breast-specific & Hallmark miRNA biology
- Subtype-specific miRNA programs: miR-200 family high in luminal/epithelial, lost in TNBC (EMT);
  let-7 high in LumA; miR-21/210/221 up and miR-10b/145/205 down in TNBC; **AGO2-PAR-CLIP shows the
  *target regulation* itself differs by subtype**, not just miRNA levels **[assert]**. → **Adopt:** the
  model is inherently cohort/subtype-specific; subtype is a legitimate "environment" (ties to §Decision
  H and ICP). The §1h arm/target table is consistent with current literature and I use it as the
  **known-driver recovery** panel (face validity), not as fitted priors on top of the sequence evidence.

### 2.9 Evaluation & overfitting
- Sequence target predictors over-predict by ≥7-fold; CLIP shows binding, not function, and carries its
  own false positives; robust evaluation uses held-out CV, orthogonal-layer and perturbation validation
  **[assert: Pinzón 2017; Kern 2021]**. → **Adopt:** never trust in-sample fit; every promotion gate is a
  pre-named held-out / orthogonal win (§Decision J, §6).

**Net contradiction check.** Nothing in the refreshed literature contradicts the §1 biology or §2 data
assumptions. The one *sharpening*: "site exists" is a much weaker prior than the brief's site-type
hierarchy alone implies — the **graded K_D (with flanking context)** is the strong prior, and the
biochemical model's *unresolved free-AGO scalar* is a candid admission that the abundance gauge (§3f) is
genuinely unidentified from expression alone and must be anchored by external biology (the threshold).

---

## 3. Decision log

Each entry: **Decision / How / Why / Alternatives kept in play (and when each wins) / Assumptions &
risks (and how to detect failure).**

### Decision A — Model object (§3a)

**Decision.** Primary estimator = **program-wise (per-Hallmark) hierarchical Bayesian regression** with
a spike-and-slab + informative-normal prior on non-negative weights, an **occupancy link** for the
abundance→repression map, seed-family partial pooling, an unpenalised confounder block, and primed
cooperative terms. Ship a **frequentist surrogate** (non-negative adaptive group-elastic-net) as the
cheap MVP. Keep **low-rank RRR/NMF as a structural overlay** and a **feature/FM model as the orphan
prior-mean generator**, not as the primary edge estimator.

**How.** For a program `P` with genes `G_P` and candidate regulators `A_P`, fit the generative model of
§4 jointly over `G_P` (borrowing strength within the program), with `S≈1,060` tumour samples (primary
tumour, PAM50 **Normal-like excluded**; §4.1). Within a
program `|G_P|≈50–200`, `|A_P|` a few hundred candidate arms → still `p>n`, made identifiable by
sparse+group+informative priors. Estimation by mean-field VI (fast) or Gibbs with a spike-and-slab
(gold), each gene's confounder block unpenalised.

**Why.**
- *Biology (§1c, §1h):* the Hallmark program is the native unit; a program-scoped `M` reads as "the
  regulatory matrix of EMT / P53", ties to the architecture layer, and keeps `p≫n` tractable.
- *Data (§2a):* `S≈1,060` (primary tumour, Normal-like excluded), `A≈2,236` (detectable far fewer),
  `G≈4,384`. A single universe fit is
  `p≫n` at its worst and maximally circular; per-program scoping cuts `p` an order of magnitude while
  staying interpretable.
- *Identifiability:* the informative prior + non-negativity + family groups convert the ill-posed
  `p≫n` problem into a well-posed MAP/posterior; the Bayesian form is the natural home for the
  multi-source prior (Decision E) and, crucially, **reports uncertainty** so the family
  non-identifiability (Decision F) is surfaced, not hidden.

**Alternatives kept in play.**
- **Per-gene adaptive-lasso (gene-focused).** *Wins* as the MVP and for deep single-driver dives (PTEN
  ~91 arms, GATA3, ESR1, ZEB1) where you want one clean penalised regression, maximum interpretability,
  and direct comparison to the existing per-edge coupling. *Loses* cross-gene strength; if orphan
  support is opened it becomes `p≫n` per gene. Use it first; it is Phase 1.
- **Low-rank reduced-rank regression / NMF (comodules).** *Wins* when regulation is genuinely
  low-dimensional (a few latent programs) and you want a structural roll-up / borrow maximal strength;
  it is the SNMNMF/PIMiM/Bayesian-sparse-RRR family. *Loses* per-edge interpretability (a low-rank `M`
  is dense and fights the sparse curated prior). Keep it as an **overlay** that explains residual
  covariance, not as the edge estimator.
- **Feature / factorization-machine model** `M_{g,m}=φ(g)ᵀΘψ(m)`. *Wins* for **reaching orphan edges** —
  it learns a sequence/evidence→weight map that generalises beyond curated support. *Loses*
  identification (indirect) and is the most circular (it re-learns the sequence prior). Use it to
  **generate `μ` for orphans**, feeding Decision E, not as the primary.
- **Mechanistic occupancy GLM (fully nonlinear).** *Wins* for maximal mechanism fidelity. *Loses*
  scalability and needs K_D for canonical sites only. Adopt its **link** inside the regression rather
  than fitting the full nonlinear system everywhere.

**What I would *not* do.** A single universe-wide **dense** RRR as the primary estimator
(identifiability bites hardest, circularity highest, interpretability lowest), or a **deep-learning
target predictor retrained on expression** (overfits, uninterpretable, and just relearns the sequence
prior we already have, per §2.2/§2.9).

**Assumptions & risks.** Assumes within-program regulation is approximately additive-in-occupancy with
sparse support. *Fails* if a few dense latent factors dominate (detect: large residual covariance after
the sparse fit → promote the RRR overlay) or if programs are too small for stable fits (detect:
posterior inclusion probabilities never separate from the prior → fall back to gene-focused or pool
programs).

---

### Decision B — Observation target `Y` (§6) **[most load-bearing]**

**Decision.** **Primary `Y` = tumour mRNA level, regressed with an explicit transcription-control block
in `C`** (not naked level). **Validation `Y` = CPTAC protein *and* protein-vs-mRNA discordance**, plus an
independent mRNA cohort (SCAN-B / METABRIC / Buffa). Use **mRNA Δ vs a subtype centroid** as the primary
only in the state-differenced variant (Decision H), never Δ-vs-NAT.

**How.** Model `Y_{g,s} = α_g + C_sβ_g − ρ_{g,s} + ε`, where `ρ` is the repression aggregate (§4). **`C`
should be the *minimal* set of common causes, not a kitchen sink** (feedback, 2026-07-04): over-adjustment
inflates variance and, worse, risks **collider / over-control bias** if a covariate is a *consequence* of
the miRNA's effect rather than a cause of both arm and target. Reasoning through each candidate against
"is it a plausible common cause of *this* (arm→gene) pair?":
- **Purity — keep, but THINK about *which* purity** (feedback, 2026-07-04). CPE (Consensus Purity
  Estimate) and ABSOLUTE are **CN/methylation/multi-modal** — largely *orthogonal to expression*. The
  user's own **breast single-cell-atlas CIBERSORTx** run yields a **Cancer-Epithelial fraction**
  (`output/brca_deconvolution/tcga_cibersortx_fractions.tsv`) = an *expression-derived* malignant purity.
  For a confounder in an **mRNA** regression this matters: an expression-derived purity **shares noise
  with the outcome `Y`** (both are linear combinations of the same expression), so using it as the purity
  term risks **partial collinearity / over-control**. → **Decision (2026-07-04, user):** **use CPE /
  ESTIMATE for now** — the established purity scalar (CPE is multi-modal and preferred; ESTIMATE is
  expression-based but on immune/stromal signature genes largely *disjoint* from the target panel, so its
  circularity is mild). The deconvolution Cancer-Epithelial fraction is kept as a *composition* covariate
  and a robustness variant, **not** the purity term.
- **Target-gene CN — keep (target-side).** CN drives target mRNA directly and independently of the
  miRNA; omitting it inflates apparent repression on deleted genes.
- **Proliferation — keep exactly ONE, and make it *malignant-compartment*** (feedback). Proliferation
  co-drives oncomiRs and their targets via the E2F/MYC axis, but the raw E2F/G2M bulk metagene is itself
  confounded with tumour-cell fraction (it rises with purity, not only with per-cell cycling). *Correction
  to my first draft:* the clinical table's `cibersortx_*` columns are **IFNG-only (LM22, immune)** — a
  non-immune *residual*, **not** a purified malignant signal. The right tool is the user's **breast
  single-cell-atlas CIBERSORTx** (cell types incl. **Cancer Epithelial**): compute proliferation on the
  **Cancer-Epithelial HiRes GEP** (imputed malignant-cell expression — the README's "phase 2" custom
  epithelial signature), or, as the immediately-available approximation, **residualise bulk proliferation
  on the Cancer-Epithelial fraction**. Keep `thornsson_proliferation` / `rppa_panel_proliferation` and the
  three E2F/G2M/MKI67 metagenes only as *sensitivity* variants (collinear — do **not** stack). *Caveat:*
  proliferation is partly *downstream* of the E2F/G2M/MYC Hallmarks, so for those programs conditioning on
  it risks over-control — report with and without it there.
- **Composition — use tissue-matched *deconvolution compartment fractions*, and they are large-effect for
  stromal programs** (feedback + repo evidence). Replace the generic ESTIMATE Immune+Stromal / LM22 with
  the **breast single-cell-atlas deconvolution fractions** the user built — CIBERSORTx (Cancer Epithelial,
  Normal Epithelial, CAFs, PVL, Endothelial, B/T/Myeloid/Plasmablasts) or BayesPrism (MH-70). The repo's
  own **coupling deconvolution retest** (`core_coupling_deconv_retest`, 2026-06-26) shows this is the
  **gold standard and that the metagene *under-corrects***: with 8 cell fractions + proliferation,
  732/1,013 headline edges stay FDR-negative (median ρ −0.148 → −0.109), but the **most-attenuated edges
  are *all* miR-29b→ECM/collagen** (COL5A1/2/3, COL3A1, FBN1, SPARC, LOXL2, ADAM12, PDGFRB: −0.3/−0.4 →
  ≈0) — the composition confound made concrete (miR-29 and its ECM targets both track CAF/stromal
  fraction). So compartment fractions are **essential and program-conditional**: mandatory for
  stroma/immune-expressed programs (EMT/ECM, IFN, INFLAMMATORY, ANGIOGENESIS); for cell-intrinsic programs
  (E2F, G2M, MYC, DNA-repair) purity + one malignant-compartment proliferation term usually suffice.
- **HRD — drop from the default (reconsidered).** HRD is a genomic-instability score that is largely a
  *correlate* of subtype/purity/proliferation, not an independent common cause of a specific (arm→gene)
  pair; once purity + proliferation + PAM50 are in, it is near-redundant. It sits in the current spine
  (`CONFOUNDER_NUMERIC = CPE + thornsson_hrd_score`) mostly by inertia. Keep it only where an ablation
  shows it moves the estimate.
- **Batch — drop unless diagnosed.** TCGA-BRCA miRNA-seq is single-platform (BCGSC); include the plate
  factor (`plate_both`) only if a batch diagnostic (guided-PCA / PVCA on the arms) shows a real effect.
- **Transcription proxy — keep (the new, essential piece).** Not currently in the repo covariates; in
  decreasing cleanliness: (i) **intronic / pre-mRNA read signal** for the target (direct transcription-
  rate proxy if derivable from STAR), (ii) **TF-regulon activity** for the gene's known TFs, (iii)
  **latent expression factors** (PEER/PCA) soaking up global transcriptional programs.

So the **minimal default `C` = {purity = CPE/ESTIMATE (decided), target-CN, one
malignant-compartment proliferation (Cancer-Epithelial CIBERSORTx GEP), transcription proxy}**, plus
**tissue-matched deconvolution compartment fractions** for stroma/immune-expressed programs; HRD and batch
off by default.
Validation: refit-free scoring on CPTAC protein + a **MoPC-style** partial correlation of mRNA–protein
*conditioned on the arm* to isolate the translational residual.

**Why.**
- *Mechanism (§1b):* mRNA carries 66–90% of steady-state repression **[assert: Guo 2010; Eichhorn
  2014]** → mRNA is the highest-powered readout and where the curated evidence lives. But miRNA acts
  post-transcriptionally: `target = transcription − repression − other`. Regressing on **level** makes
  transcription a confounder — a TF co-driving a miRNA and its target manufactures a spurious *positive*
  weight (or masks a real negative one). This is the single largest defect in the current heuristic,
  which regresses/couples on level with only purity+HRD+prolif in `C` and **no transcription proxy**.
- *The matrix view's gift (§2 of the earlier note):* the confounder that drives both `x` and `y`
  creates omitted-variable bias in `M` unless it is in the design. Putting transcription in `C` is the
  fix, and it is *cheaper and more honest* than differencing.
- *The translational residual is only visible in protein* **[assert: §1b, Selbach]** → protein and
  especially **discordance** (mRNA present, protein down) is the layer that certifies an edge is truly
  post-transcriptional and not a surviving transcriptional confound. This makes protein the decisive
  non-circular OOD check, not a secondary nicety.

**Alternatives kept in play.**
- **Naked mRNA level (no transcription proxy).** *Wins* only on simplicity/coverage. *Loses* to
  confounding — dominated for causal claims. Keep only as an ablation to *quantify* the transcription
  confound (compare weights with/without the proxy; the delta is the confound size).
- **mRNA residualised on the transcription proxy (fit on residuals).** *Wins* when you trust the proxy
  and want the post-transcriptional part isolated *before* fitting. *Loses* if the proxy is imperfect
  (over-residualising removes real signal). Near-equivalent to putting the proxy in `C`; I prefer
  in-`C` because it estimates the adjustment jointly rather than committing to a fixed residualisation.
- **mRNA Δ vs reference (GTEx / subtype centroid).** *Wins* paired with the state-differenced `M`
  (Decision H): gene-specific baseline transcription cancels, improving identification. *Loses* power
  (differencing inflates noise) and, if referenced to **NAT**, injects field-effect contamination
  **[assert: §2.7]** — so reference to GTEx (cross-platform noise cost) or subtype centroid (platform
  constant, preferred) only.
- **Protein level as *primary*.** *Wins* if the question is function-first and CPTAC coverage suffices.
  *Loses* on `n` (smaller), MS batch, and that it *understates* the ~75% of repression realised at mRNA.
  Keep as the co-primary for the subset of genes with strong CPTAC coverage and as the validation layer.
- **Protein-vs-mRNA discordance as *primary*.** *Wins* as the cleanest post-transcriptional signal.
  *Loses* on `n` and needs paired CPTAC. Strongest as the **validation** estimand (Decision J, Bar 5).

**Assumptions & risks.** Assumes a usable transcription proxy exists. *Fails* if no proxy captures the
relevant TF programs (detect: learned positive weights on canonical repressive edges → proxy inadequate;
escalate to differencing or protein-primary). Assumes the mRNA↔protein map is informative; *fails* for
translational-only targets (detect: edge couples on protein/discordance but not mRNA → a real edge the
mRNA readout structurally misses; flag as translational-class).

---

### Decision C — The a·w scale gauge (§3f) **[must-resolve]**

**Decision.** **Fix arm abundance on the `log2(RPM+1)` scale; use a saturating occupancy link with a
functional-copy-number threshold; let `M` carry per-site affinity (`1/K`).** Report **two** quantities,
never one: the **potential weight** `M_{g,m}` (per-unit affinity, gauge-anchored) and the **realised
contribution** `w_eff(m,g,s)=gate(x_{m,s})·M_{g,m}` (abundance-gated). The functional threshold — a
physical, non-rescalable anchor — is what breaks the degeneracy.

**How.** The contribution is multiplicatively non-identifiable in the *linear* form `a·w` (rescale
`x_m→c·x_m`, `M→M/c`, identical `Y`). Note this is **not** cured merely by going saturating:
`N=a/(a+K)` is itself invariant to `(a,K)→(ca,cK)` — the *ratio* `a/K` is identified, not `a` and `K`
separately. What pins the gauge is **external biology in absolute units**:
1. **Abundance scale fixed** at `log2(RPM+1)` (no per-arm variance standardisation), so a rare arm and
   an abundant arm are *not* equalised — the abundant arm loads more RISC, as biology demands (§1e).
2. **A saturating gate with an absolute threshold** `gate(x)= (x−x0)_+ / (κ + (x−x0)_+)`-style (or the
   biochemical Langmuir occupancy), with the half-max / floor `x0` anchored to the **functional
   copy-number threshold** (~hundreds of copies; **[assert: Mullokandov 2012]**) mapped onto the RPM
   scale. Below `x0`, contribution ≈ 0 *regardless of site quality* — the biology that makes the gauge
   physical rather than conventional.
So "a strong edge" = **high per-site affinity `1/K`** *among arms that clear the loading threshold*; a
rare high-affinity arm has large `M` but ~0 realised repression.

**Why.**
- *Mechanism (§1e; §2.2):* occupancy is the biochemically correct map from abundance to repression
  **[assert: McGeary 2019]**, and the free-AGO scalar that McGeary leaves unfit *is* our gauge freedom —
  anchoring it with the Mullokandov threshold is the principled resolution.
- *It retires `D(m)`:* the promiscuity denominator was a hand-set proxy for "finite RISC spread over
  many targets." With a **shared free-AGO pool** (Decision G), a promiscuous arm mechanically dilutes
  its own per-target occupancy — promiscuity becomes a *consequence*, not an input (see Decision G).
- *It aligns the gauge with the question (coupling-invariance lemma, repo `EDGE_QUESTION_TAXONOMY`):*
  for a **single edge**, Spearman coupling is invariant to any monotone edge weight and to the gauge — so
  the gauge only matters for **attribution/identity** (who owns the gene) and **aggregate magnitude**
  (how much total). Fixing it once and exposing potential-vs-realised serves both without contaminating
  the coupling test.

**Alternatives kept in play.**
- **Per-arm variance standardisation (z-scoring `x`).** *Wins* when you only care about *dynamics*
  ("who moves across tumours") and want scale-free coupling; it is the current `z`-spine. *Loses* by
  erasing cross-arm abundance (a rare and an abundant arm both get SD=1) — biologically wrong for
  loading, and rank-invisible to the threshold. Keep for the *coupling* readout, not for attribution or
  magnitude.
- **Row/column-norm constraint on `M`** (`‖M_{·,m}‖=1`). *Wins* mathematically (removes the degeneracy
  cleanly). *Loses* biological meaning (an arbitrary norm). Dominated by the threshold gauge, which is
  physical.
- **Fix abundance, let `M` carry magnitude, no gate (linear).** *Wins* on simplicity and is what the
  RRR sketch assumed. *Loses* the threshold and saturation biology; a low-abundance arm with a great
  site gets over-credited. Keep as the unsaturated-regime approximation (`a≪K` ⇒ occ≈a/K) where it is
  valid.

**Assumptions & risks.** Assumes RPM is a usable monotone proxy for per-cell copy number *up to a
per-sample constant*. *Fails* under strong compositional distortion (a global miRNA-content shift
rescales all RPMs) — detect via the AGO/RISC-capacity covariate and via CN-instrument concordance;
mitigate by the shared-pool normalisation (Decision G) which absorbs per-sample loading. The threshold
`x0` is uncertain; treat it as a **hyperparameter tuned out-of-fold**, and report sensitivity.

---

### Decision D — Inclusion / HE rule and hard-vs-soft mechanism (§3b)

**Decision.** **Soft, evidence-graded spike-and-slab inclusion** as the mechanism, run in **two
disciplined modes side by side**: a **hard-HE "confirmatory" model** (trusted backbone) and a **soft
"discovery" model** (orphans admitted at low prior, must clear OOD to be believed). **Redefine HE** away
from the miRTarBase-only rule by folding in TarBase-v9 method/confidence, a PMID-deduped AGO-CLIP axis
(chimeric confirmatory), context++/K_D magnitude, and breast literature (Decision E builds the grade).

**How.** Inclusion indicator `z_{g,m} ~ Bernoulli(π_{g,m})`, `π` a monotone fusion of the re-graded
evidence (Decision E). The confirmatory model sets `π=1` on the HE set and `π=0` elsewhere (hard mask);
the discovery model sets `π∈(0,1)` graded, high for strong edges, low-but-nonzero for orphans with real
sequence/physical support — so a weak edge enters **only if the expression insists** while the prior
holds the noise floor down. The two are the same generative model at two prior settings, so they are
directly comparable and the discovery model **nests** the confirmatory one.

**Why.**
- *It is not fence-sitting — the two answer different §0 questions:* confirmatory = "reproduce curated
  knowledge more parsimoniously / test it"; discovery = "extend to edges curation never had." Both
  references are live per §0; forcing a single mode discards one question.
- *Biology (§1d):* orphans often carry strong context++/K_D/chimeric support (60% of our high-value
  orphans have direct duplex reads **[assert: Manakov, repo memory]**) — hard-restricting throws away
  the continuous signal and the discovery lane. But orphans are also where false positives live
  (predictors over-predict ≥7-fold **[assert: Pinzón 2017]**) — so they enter *soft*, gated by OOD, not
  as trusted coefficients.
- *Identifiability:* the hard-HE model is the most identifiable (cleanest support); it is the backbone
  whose weights you quote. The soft model trades identifiability for discovery under an explicit
  validation gate.

**Alternatives kept in play.**
- **Hard-restrict only.** *Wins* when the goal is a defensible, maximally-identifiable curated-weight
  refinement and you distrust discovery. *Loses* all orphan discovery and the continuous signal.
- **Soft only (no confirmatory backbone).** *Wins* for maximal discovery. *Loses* a trusted reference
  and invites over-fitting the orphan tail. Dominated by running both.
- **Grade granularity** is itself a choice: binary HE/not, vs a 3-tier (HE / evidence-tail / orphan),
  vs a fully continuous `π`. *Continuous wins* for the soft model (uses all information); *tiered wins*
  for reporting and for the hard confirmatory mask. Ship continuous `π` with tier labels for
  interpretability.

**Assumptions & risks.** Assumes the evidence grade is calibrated to true functionality. *Fails* if a
source is miscalibrated (e.g. CLIP depth admits binding-not-repression) — detect by the source-reliability
calibration in Decision E (does this evidence class predict held-out protein repression?); down-weight the
offending source. Spike-and-slab can be unstable across folds for near-threshold edges — report
**posterior inclusion probabilities**, not hard in/out, and treat PIP as the discovery score.

---

### Decision E — Multi-source prior fusion (§3b/§3e) **[core differentiator]**

**Decision.** A **two-layer evidence model** honouring the edge-level vs site-level split: **site-level
evidence → magnitude prior `μ` and interaction support**; **edge-level evidence → inclusion `π`, sign,
and prior strength `τ`**. Fuse sources by **calibrated log-odds/log-magnitude addition** with
source-reliability weights *learned against held-out functional repression*, not hand-set. Cap data
override by `τ`.

**How.**
- **Site-level (magnitude + interactions).** context++ / biochemical K_D are per-site, continuous,
  mechanistic, and defined for orphans. Aggregate sites in the gene's 3′UTR by the **occupancy sum**
  (not a mean), so multiplicity is encoded:
  `μ_{g,m} = log(1 + Σ_{sites s} occ_prior(K_{s}))`, `occ_prior` from the biochemical Langmuir at a
  reference abundance. This gives a magnitude prior for **every** predicted edge, orphans included, and
  is strictly better than the current `log(1+evidence)` (which is flat within a gene, max/min≈1.8×, and
  blind to multiplicity). Site **coordinates** prime interaction support: distinct sites 8–40 nt apart →
  a cooperative product candidate (Decision G); same/overlapping sites → sub-additive grouping. **Which
  predicted sites are actually occupied** is marked by the *site-level* physical layers — **POSTAR3 AGO2
  peaks** (dense), TarBase-v9 CLIP coordinates, and Manakov chimeric duplex coordinates — which up-weight
  `occ_prior` on covered sites and confirm the cooperative-pair support.
- **Edge-level (inclusion + sign + confidence).** A logit for `π`:
  `logit π_{g,m} = w0 + w_tier·tier(curated MTI: reporter/protein/perturbation > RNA > binding)
                       + w_meth·TarBase9_method_conf(PMID-deduped)  # NB regulation ≈ all-Negative → not a direction channel
                       + w_ago·log1p(AGO-CLIP: ENCORI∪POSTAR3∪TarBase-CLIP, PMID-deduped) + w_chim·log1p(chimeric duplex reads)
                       + w_lit·log1p(breast PMIDs) + w_seq·seq_support_indicator`   # POSTAR3 peaks → μ (site level).
  **Sign — with a big caveat** (checked on the live file): TarBase v9's `regulation` is **near-degenerate,
  4,724,236 Negative vs 301 Positive** — an almost-constant curation default (most rows are CLIP, which
  measures *binding, not direction*), **not** an informative per-edge sign, so it is *not* the rich
  directional prior my first draft claimed. Actionable use: route the **~301 explicit Positive** edges to
  the *activating/indirect* class or exclude them (a validated miRNA-up→target-up "edge" must not be forced
  to `M≥0`); **otherwise take sign from genuine functional assays** (miRTarBase functional tiers,
  luciferase, western) plus the `M≥0` prior — never from TarBase's blanket "Negative".
- **The two CLIP layers are different tools — lean on AGO-CLIP as primary, and use *several* AGO-CLIP
  sources** (feedback, 2026-07-04). **AGO-CLIP (HITS/PAR/eCLIP)** is dense, replicated — but **miRNA-
  *agnostic*** (it shows AGO occupies a region; *which arm* is inferred by seed-matching), so it is the
  **workhorse for region/edge inclusion `π`** and for marking bound sites in the cooperativity window
  (Decision G). **Do not rely on ENCORI alone — there are multiple in-repo AGO-CLIP sources, and ENCORI
  is sparse** (the method_dev site-ladder found ENCORI AGO-CLIP covered only ~3% of predicted sites):
  - **POSTAR3** carries **~2.36M AGO2 binding-site rows** (`data/RBP-RNA/POSTAR3.parquet`) — a far denser,
    uniformly-processed compendium; already in-repo (used for lncRNA-RBP, repurposable here).
  - **ENCORI/starBase** ships edge-level *interaction calls* (the current α=0.5 boost).
  - **breast-context AGO2-PAR-CLIP** (Farazi 2014, subtype-resolved) — a **tissue-matched** physical prior
    and a subtype signal (ties to Decision H); a literature add, not yet in-repo.
  **Do not double-count — and this critically includes TarBase v9** (feedback, 2026-07-04). Checked
  against the live file: **TarBase v9 is ~93% AGO-CLIP itself** (HITS-CLIP 2.93M + PAR-CLIP 1.49M of 4.72M
  human rows) plus chimeric (qCLASH/CLASH/Chimeric-fragments ~0.14M, the **same class as Manakov**), drawn
  from only **998 distinct PubMed IDs** — so ENCORI, POSTAR3, and TarBase's CLIP rows are largely the
  **same underlying experiments**. Therefore treat *physical AGO-CLIP binding as ONE axis* =
  **union/max over {ENCORI, POSTAR3, TarBase-CLIP rows}, de-duplicated by `article_pubmed_id`**, never a
  sum; and fold TarBase's qCLASH/CLASH rows into the **chimeric** axis with Manakov. Following the repo's
  site-ladder rule, ENCORI is the edge-level representative (`π`) and POSTAR3 the site-level representative
  (`μ`/occupancy). **TarBase v9's genuine value-add is NOT "more CLIP"** — it is the thin **functional
  layer** (Luciferase 3.3k, Western 1.3k, pSILAC 6.5k, qPCR 961 → functional-repression tier, PMID-deduped
  against miRTarBase), plus **confidence, tissue/cell_line, `microt_score`, and site coordinates**.
  **Chimeric eCLIP / CLASH (Manakov)** uniquely ligates the miRNA to its target, giving the direct
  **(arm, site) pair** — invaluable for *arm assignment* and interaction priming — **but it is sparse and
  its reproducibility is contested** (a view in the field is that chimeric maps are hard to replicate and
  capture many non-canonical, uncertain-function duplexes) **[concern raised; treat as a caveat, not
  settled]**. Therefore: **use AGO-CLIP depth as the primary physical `π` signal**; use chimeric only as a
  **confirmatory arm-assignment / site-coordinate** layer *where its depth is adequate*, and **fall back
  to AGO-CLIP region + seed inference when chimeric is too sparse** for an arm/gene. Neither sets sign or
  large `μ` on its own — proximity ≠ repression (the 0.5 CLIP discount, now mechanistic); functional /
  directional evidence is still required for magnitude and sign.
- **Provenance de-duplication is GLOBAL — one `(PMID × assay-class)` counts once, across *every*
  database and *both* the physical and the functional layers, before any fusion** (feedback, 2026-07-04).
  **This is not done today and is a real gap:** the current `evidence_score` is built from miRTarBase's own
  per-class **study counts** (`n_reporter__functional_mti_studies`, …) with the ENCORI AGO-CLIP depth
  **added on top** (α=0.5), and TarBase is not yet wired — so there is **no cross-database PMID dedup**, and
  the edges table the model consumes (`mirna_hallmark_edges.tsv.gz`) has **already collapsed to counts with
  no PMID column**. The same paper therefore double-counts: (i) a **CLIP** study surfaced by ENCORI *and*
  POSTAR3 *and* TarBase (TarBase = 998 PMIDs → millions of rows); (ii) a **functional** paper (a luciferase
  or western assay) curated by **both miRTarBase and TarBase** — the functional tier double-counts exactly
  like the CLIP tier. **Fix — build a unified evidence ledger at *edge-build* time** (upstream of the
  count-collapse), keyed by `(m, g, pmid, assay_class, direction, tissue)`, drawing PMIDs from
  miRTarBase `References (PMID)`, TarBase `article_pubmed_id`, and each CLIP/chimeric study; **dedupe to
  distinct `(m, g, pmid, assay_class)`** (dedupe on PMID×assay, *not* bare PMID — one paper with both a
  reporter *and* a western contributes two independent assay-classes), then aggregate per edge per class as
  a **count of distinct supporting PMIDs**, not rows/peaks/DB-hits. Only then apply the tier weights + log1p.
  Two provisos: (a) the **sequence/biochemical priors** (context++, K_D, TargetScan) have **no PMID and are
  a separate orthogonal channel** — they are *predictions*, correctly excluded from experimental dedup;
  (b) where a source exposes no PMID (some ENCORI cache columns → only `clipExpNum`; POSTAR3 → dataset IDs),
  dedupe by the best available experiment/GEO/cell-line identifier and **cap** the DB's own count as an
  upper bound, flagging the edge as lower-provenance.
- **The fusion is METHOD-centric, not source-centric** (the right frame for the model we are building).
  The evidence axis is the **assay method** — {reporter/luciferase, protein/western, qPCR/RNA, perturbation,
  AGO-CLIP, chimeric/CLASH, degradome} — and **every database pours its distinct PMIDs *into* those method
  bins**; a database is *never* its own additive term. This is the correct structure precisely because it is
  where the dedup lives (dedup = distinct PMIDs *per method, across sources*). *Diagnosis of today's score:*
  it is **method-centric *within* miRTarBase** (its cross-count classes `reporter/protein/rna/perturbation/
  binding` already *are* assay methods) **but source-centric *across* databases** — ENCORI is bolted on as
  its own composite `+ α·enc_depth`, so the **AGO-CLIP method is split** (miRTarBase `n_binding` at w=0.5
  *and* ENCORI `clipExpNum` at w=2.0·α, un-pooled, un-deduped), and worse, `enc_depth` **mixes sequence
  *predictions*** (`𝟙[TargetScan∧miRanda∧PITA]`) into what should be a physical-assay term. The build folds
  all of these into one method-centric, PMID-deduped set of bins, with predictions moved to the separate
  sequence/K_D channel (μ).
- **Fusion & calibration — the method weights are LEARNED, not hand-set, and *not* on TCGA coupling.**
  The per-method contributions combine additively on the logit/log scale (Bayesian evidence combination)
  **over the de-duplicated distinct-PMID counts** (never raw per-DB counts). The method weights `w_method`
  (reporter vs western vs CLIP …), the weak-discount, and the relative binding-vs-functional weight are
  **fit** — a logistic/Cox regression (or the outer Bayesian model's hyperprior) predicting an **external
  validated-repression outcome** the coupling test never sees (perturbation ΔlogFC from transfection/KO
  series; a held-out functional-MTI label), by **nested cross-validation, with uncertainty on each weight**.
  **This is a real correction to the current pipeline**, whose weights (`w_rep=3.0, w_prot=2.5, w_rna=1.5,
  w_bind=0.5, δ=0.3, α=0.5`) are **hand-set constants** "selected" once by a small grid
  (`dual_spine_comparison`) on **TCGA coupling** — i.e. tuned on the very outcome the model is later scored
  against (**circular**, and it violates the guardrail below). Everything non-prior (`λ, τ, x0`, pool
  params, low-rank rank) is likewise tuned **out-of-fold**, never in-sample. `τ_{g,m}` (prior strength)
  decreases with evidence: HE edges get tight `τ` (data barely moves them); orphans get loose `τ` (data free
  to speak, but low `π` holds the floor).
- **Override budget & non-circularity.** How far expression may move an edge is set by `τ`. All of
  `π,μ,τ` and the interaction support come from **sequence + curation + external-tissue assays only** —
  never TCGA `(X,Y)` — so the guardrail (§8 of the earlier note) holds.

**Why.**
- *It is the stated differentiator (§1 of the earlier note):* prior at **magnitude**, not just support;
  quality-weighted, assay-directness-graded, breast-audited. GenMiR++/LASSO use sequence for support
  only.
- *Data (§2c), corrected against the live file:* TarBase-v9 (4.72M rows) is **~93% AGO-CLIP** (HITS/PAR-
  CLIP) + chimeric (qCLASH/CLASH) from **998 PMIDs**, with **site coordinates** and `microt_score` but a
  **near-degenerate direction** (301 Positive) — so its first-class value is **site coordinates + method
  breadth/confidence + a thin functional layer**, fused into the *shared* CLIP/chimeric axes (PMID-deduped,
  not additive to ENCORI/POSTAR3/Manakov), *not* an independent directional layer. Manakov chimeric
  (direct duplex, **site-level**) is likewise partly shared with TarBase's qCLASH. The biochemical K_D
  gives the mechanistic magnitude for orphans. The genuinely under-used, non-redundant signals are the
  **site coordinates** (interaction priming) and **K_D magnitude**, not "more CLIP mass".
- *Edge vs site reconciliation:* the split resolves the brief's hardest fusion question — edge-level
  assertions set *whether and which sign*, site-level scores set *how much and which interactions*.

**Alternatives kept in play.**
- **Single fused scalar (no edge/site split).** *Wins* on simplicity. *Loses* the ability to prime
  interactions and to use orphan site magnitudes — dominated for the cooperativity and orphan goals.
- **Hand-set source weights (current).** *Wins* on transparency and needs no calibration data. *Loses*
  calibration; keep as the initialisation and as a robustness check (does calibration move conclusions?).
- **Feature/FM model for `μ`** (Decision A alternative) instead of the occupancy-sum. *Wins* when you
  want a *learned* sequence→magnitude map (more expressive than the biochemical prior). *Loses*
  circularity headroom. Use it only to fill `μ` where the biochemical model has no K_D (non-canonical
  sites).

**Assumptions & risks.** Assumes sources are conditionally informative and the calibration set is
representative. **The dominant risk is source correlation / counting the same paper more than once — in
*both* the CLIP tier and the functional tier**: ENCORI/POSTAR3/TarBase-v9 (~93% HITS/PAR-CLIP, 998 PMIDs)
share CLIP experiments, **and** miRTarBase ∩ TarBase share the same curated *functional* papers (luciferase,
western). **Mandatorily de-duplicate by `(PMID × assay-class)` globally** (the provenance-ledger step above),
or the prior fakes independent confirmation — and note this is a change from the current pipeline, which
does *not* dedup and *adds* ENCORI on top of miRTarBase counts. TarBase's `regulation` is near-constant
"Negative", so it must **not** be used as a directional/inclusion signal beyond the ~301 explicit Positives;
carry its `tissue`/`cell_line` to down-weight non-breast rows, and take sign from functional assays.

---

### Decision F — Family identifiability (§3d)

**Decision.** **Hierarchical partial pooling toward a seed-family latent**: `M_{g,m}=M_{g,F(m)}+δ_{g,m}`,
`δ~N(0,ρ_F²)`. The **identified estimand is family→gene** `M_{g,F}`; the per-arm split is a **nomination
reported with its posterior width**, never a point coefficient. Use the **real seed** (TargetScan
`miR_Family_Info`, positions 2–7/2–8), not the name-stem proxy.

**How.** Same-seed arms are near-collinear columns of `X` (within-family ρ≈2× between; hubs 20–37%
redundant **[assert: repo]**). Partial pooling makes `M_{g,F}` identified (the shared, well-conditioned
direction) and `δ` weakly identified: when members are perfectly collinear, `δ`'s posterior = its prior
("cannot split"); when a member's abundance diverges (arm-switch, differential expression, isomiR),
`δ` becomes identified and the model *can* resolve the member.

**Why.**
- *Biology (§1g):* seed-family redundancy is *mandated*, not a numerical nuisance — same-seed members
  hit overlapping sites and are substitutable (sub-additive). The group structure *is* the sub-additivity
  model (Decision G row 1).
- *Identifiability:* plain LASSO keeps one member arbitrarily and unstably across folds — reporting an
  unidentified quantity as if identified. Partial pooling reports the ridge honestly ("family carries W;
  member A = 0.6W±0.3, not distinguishable from B") — matching the repo's identity-vs-magnitude and
  "per-arm is a nomination" stance.

**Alternatives kept in play.**
- **Group-LASSO (family in/out together).** *Wins* as the frequentist middle ground (cheaper than
  Bayes, gives the group selection). *Loses* the calibrated within-family uncertainty (still needs a
  post-hoc split rule). Use in the frequentist surrogate.
- **Collapse to one family predictor (mean/max/sum).** *Wins* for maximum power and as the MVP (already
  built: `mir301_family_depth`); the pooled family abundance is the cleanest single predictor. *Loses*
  the ability to ever resolve a member. Use when the split is not the question.
- **No grouping (per-arm).** *Wins* only if members are genuinely independent (distinct seeds
  mislabelled as a family). *Loses* to collinearity otherwise — dominated for true families.

**Assumptions & risks.** Assumes the seed family is the right grouping. *Fails* for arms whose *isomiR*
distribution shifts the seed (5′ isomiRs change positions 2–7 → a different family) — detect via isomiR
data (Decision/Open-question); a 5′-heterogeneous arm should be split by isomiR before grouping. Also
fails if pooling is too strong (masks a real member-specific edge) — detect by posterior-predictive
checks on the member with divergent abundance.

---

### Decision G — Competition, cooperativity, ceRNA (§3e)

**Decision.** Three couplings, three mechanisms, three resolutions:
1. **Seed redundancy (miRNA↔miRNA, same site, sub-additive)** → the **family group structure**
   (Decision F); optional per-site saturation for overlapping-but-different-seed arms.
2. **Cooperativity (miRNA↔miRNA, distinct sites 8–40 nt, super-additive)** → **primed product terms**
   `+ψ_{g,(m,m')}·occ_{m,g}·occ_{m',g}`, support drawn from TargetScan/TarBase site coordinates, strong
   spike-and-slab with a low base rate.
3. **ceRNA / global saturation (target↔target, shared finite RISC)** → a **cohort/sample-level shared
   free-AGO pool**, *not* a per-edge ceRNA term. This makes the promiscuity denominator `D(m)`
   **obsolete**.

**How.**
- *Cooperativity:* only arm-pairs with **distinct predicted sites 8–40 nt apart on the same 3′UTR**
  get a product term (Saetrom/Grimson/Broderick window **[assert]**). The prior does double duty:
  main-effect support (edge-level) + interaction support (site coords). Products are collinear with
  their mains, so promote a product only if it beats its mains out-of-fold.
- *Shared pool (ceRNA + promiscuity):* free-loading scalar per sample
  `φ_s = A_total,s / (1 + Σ_{m,g} occ_{m,g,s})`, and realised weight `w_eff = φ_s·occ`. A promiscuous
  arm targeting many sites raises the denominator → **dilutes its own per-target occupancy**; a
  highly-expressed target sinks occupancy from others through the same pool. This *is* the target-site-
  abundance / AGO-depletion mechanism, per-arm and mechanistic, replacing the hand-set `D(m)`.

**Why.**
- *Biology (§1d/§1f):* cooperativity and redundancy have **opposite signs** and must not be conflated
  (the earlier note's §4). ceRNA is real but narrow — physiological only near equimolarity/low target
  abundance **[assert: Denzler 2014/16; Bosson 2014]** — so a **global** shared-pool term is the honest
  resolution, and a per-edge ceRNA term would over-fit a rare effect.
- *Retiring `D(m)`:* promiscuity as a *consequence* of finite-RISC competition is both more principled
  and testable — the model should reproduce the promiscuity penalty endogenously, and we can check that
  the fitted per-target occupancy of hub arms is depressed without any `D(m)`.

**Alternatives kept in play.**
- **No interaction terms (pure additive).** *Wins* for identifiability and as the MVP — additivity is a
  fine *local* first order (§1c). *Loses* all cooperativity/saturation. Start here; add terms only if
  they earn OOF.
- **Per-edge ceRNA terms.** *Wins* if a specific low-abundance family (miR-25-like) shows physiological
  competition. *Loses* generally (over-fits a rare effect). Reserve for named low-abundance cases with
  evidence.
- **Keep `D(m)` as a hand-set denominator.** *Wins* on simplicity if the shared-pool fit is unstable.
  *Loses* principle. Keep as fallback; report whether the pool reproduces `D(m)`'s effect.

**Assumptions & risks.** Product terms risk over-fitting (`A(A−1)/2` candidates) — mitigated by strict
site-priming + shrinkage + OOF promotion. The shared-pool `φ_s` may be under-identified without an
AGO-capacity measurement — anchor it to the RISC-capacity covariate already in the repo (`ago_gate`),
and detect failure by whether `φ_s` correlates with measured AGO/TNRC6 levels.

---

### Decision H — State axis: GTEx→NAT→tumour (§7)

**Decision.** **Nested/hierarchical state model** `M_tumour = M_healthy + Δ`, shared support,
partially-pooled magnitudes, with the **healthy (GTEx) fit as the prior for the tumour fit**. Separate
`M` per state kept as an audit; state-differenced `ΔY=−M·Δx` used on the **tumour-vs-subtype-centroid**
contrast (platform-constant) where differencing helps identification. **NAT is its own state, never a
healthy stand-in.**

**How.** Fit `M_H` on GTEx healthy breast; use it as the prior mean for `M_T`, so `M_T` shrinks toward
`M_H` and `Δ = M_T − M_H` is "what the tumour data insists changed." This **separates acquired abundance
(`Δx`, the robust cross-state signal) from acquired per-site weight (`Δ` in `M`)** — the §1i distinction.
Biological prediction: `Δ≠0` should concentrate where **APA shortens 3′UTRs** (distal site loss → weight
drops) and where **AGO/TNRC6/Dicer capacity** shifts (global rescale). Test `Δ=0` as a falsifiable null:
if `Δ≈0` everywhere, a single `M` + acquired abundance suffices (a clean result).

**Why.**
- *Gauge (Decision C):* two **separately-gauged** matrices are not comparable — nesting shares the gauge
  by construction, so `Δ` is meaningful. This is the decisive argument against separate-`M`-as-primary.
- *Power:* NAT `n=104`, GTEx `n=346` (cross-platform) are weak for independent fits; nesting borrows the
  tumour `n≈1,060` to anchor and estimates only the *deviation*.
- *ICP precedent (§2.5):* Meng 2019 uses PAM50 subtypes as environments and finds invariance robust to
  TF confounding but **misses subtype-specific edges**. So run the nested model **and** flag edges whose
  `Δ` (or whose weight) is subtype/state-specific as *findings* (state-dependent regulation), not
  failures — capturing exactly what pooling-for-invariance discards.
- *Biology (§1i, §2.7):* APA/site-loss and effector-capacity change make `M` state-dependent in
  principle — this must be **tested**, and the nested model is the test.

**Alternatives kept in play.**
- **Separate `M` per state.** *Wins* for cleanest per-state interpretation if each state had ample,
  same-platform `n`. *Loses* power and cross-gauge comparability. Audit only.
- **Pooled + state interaction (`M_0 + state offset`).** *Wins* for power and simplicity. *Loses* the
  APA-driven support contraction (assumes shared support). Near-equivalent to nested when `Δ` is small;
  prefer nested because it makes `Δ` explicit and testable.
- **State-differenced `ΔY=−M·Δx`.** *Wins* on identification (baselines cancel — addresses Decision B's
  transcription confound directly). *Loses* to GTEx cross-platform noise and NAT field-effect. Use on the
  **tumour-vs-subtype-centroid** contrast (platform constant), pairing with Decision B's Δ variant.
- **Subtype-as-environment instead of tissue-state.** *Wins* for causal robustness (ICP). Complementary,
  not exclusive — run subtype environments *within* the tumour state for the invariance check.

**Assumptions & risks.** Assumes GTEx healthy is a valid prior anchor despite platform differences —
mitigated by the existing quantile-normalisation and MIMAT-join machinery, and by the miTED-anchor
imputation for GTEx-collapsed arms (already built). *Fails* if GTEx-TCGA platform effects masquerade as
`Δ` — detect by whether `Δ` concentrates in the GTEx-imputed/low-confidence arms (an artifact signature)
vs in APA-shortened genes (a biology signature).

---

### Decision I — Aggregation on a gene, attribution to a miRNA (§4)

**Decision.** **Aggregate = the fitted model's predicted repression** (occupancy-summed then saturated),
not a raw sum of hand-weighted contributions. **Attribution = expose two distinct objects, never
conflated:** **magnitude/force** = realised contribution `w_eff(m,g,s)` and its gene sum; **identity** =
a **Shapley decomposition of the fitted aggregate** across regulators. **Retire the softmax share.**

**How.**
- *Aggregate:* sum in **occupancy** space (additive loading), then apply the global saturation/pool —
  `ρ_{g,s} = sat(Σ_m w_eff + Σ coop)`. This replaces the current `AGG = Σ contributions`.
- *Identity (who owns the gene):* the **Shapley value** of each regulator's contribution to the fitted
  aggregate — the unique fair split that handles collinearity (redundant same-seed arms *split* credit
  rather than one taking it all). For the linear-in-occupancy regime this reduces to the **LMG /
  averaging-over-orderings** R² decomposition; for the nonlinear occupancy aggregate use sampling/kernel
  Shapley. This **replaces the softmax share**, which (a) spends abundance twice (once in the exponential
  share, once in `logrpm`) and (b) is rank-invisible to Spearman — both documented repo defects
  (`FORMULAS.md` §5a). The softmax share is *nested* as the degenerate Shapley of a proportional-abundance
  model, so we can show what it approximated and where it failed.
- *Magnitude (how much total):* `w_eff` and `Σ_m w_eff` directly — the force/pressure question.

**Why.**
- *The identity≠magnitude trap (repo memory `attribution-identity-vs-magnitude`):* a share summing to 1
  answers "who," not "how much." Coupling/repression claims need **magnitude** + the coupling test; "which
  arm owns the gene" needs **identity** (Shapley). Report both; enforce the repo discipline that
  **share ⊥ coupling** (a high identity share is *never* evidence of repression — only the coupling test
  is; `FORMULAS.md` §5a caveat carries over to Shapley unchanged).
- *Gene↔miRNA spectrum:* all one `M`. Column of `M` over `g` = an arm's total targeting role; row over
  `m` = a gene's total incoming pressure (layer impact); a **seed family / polycistron (miR-17~92)** =
  the group latent (Decision F) or a defined arm-set aggregate — the model exposes these *for free*
  without extra machinery.

**Alternatives kept in play.**
- **Softmax share (current).** *Wins* on cheapness and it is already wired. *Loses* on the
  double-counting and rank-invisibility. Keep only as a documented baseline to show the improvement.
- **Fitted contribution directly, no Shapley.** *Wins* when regulators are near-orthogonal (Shapley≈raw
  contribution) and you want speed. *Loses* under collinearity (over-credits one of a redundant pair) —
  which is exactly the family case. Use raw contribution for magnitude, Shapley for identity.
- **Sum vs mean aggregation; signed/positive/absolute coherence; pruning; RISC-gating.** These are
  *reporting* choices on the fitted aggregate. *Sum wins* for total force; *mean* for per-regulator
  intensity; *signed coherence* as a diagnostic (do regulators agree in direction?). Ship sum as
  headline, the others as diagnostics, mirroring the current share-metric menu but on the *fitted* model.

**Assumptions & risks.** Shapley is exponential in regulators — use the linear/LMG closed form where
valid and sampling elsewhere; for genes with many regulators, cluster to families first (Decision F) then
Shapley over families (also the identified unit). Risk: users read identity as repression — mitigate by
never shipping identity without the coupling flag beside it (the repo's `driver_candidate` = coupled AND
high share discipline).

---

### Decision J — Evaluation & references (§5, §0)

**Decision.** A **falsification ladder** with pre-named criteria, each scored against the *right*
reference per §0, plus pre-named **abandon** conditions. Identification lives in the *fit* (the `C` block
and the instrument), not only in the checks — Q1 and Q3 are one problem.

(Full test specs, success thresholds, and abandon conditions in §5.)

**Why.** The references are not interchangeable (§0): sometimes the bar is **raw abundance** (does
structure beat just-abundance?), sometimes **curated fixed-M** (does learning beat/parsimoniously match
hand weights?), sometimes a **shuffled-evidence null** (is the prior doing work?). The **coupling-
invariance lemma** (repo) forces a key design point: a *single-edge* Spearman coupling **cannot**
distinguish edge weightings (it is invariant to monotone `edge_w` and to the gauge) — so to show learned
weights beat curated ones you **must** test at the **aggregate** (multi-regulator) level or by
**out-of-fold prediction**, never by single-edge coupling. This is why Bars 1–2 are aggregate/predictive.

**Alternatives kept in play.** Which reference is "the bar" depends on the claim being made — keep all
three live and report against each, rather than declaring "beat the heuristic" the sole goal (the brief
is explicit that it is not).

**Assumptions & risks.** The instrument may be weak (many loci have little CN variance) — report
first-stage F and restrict to focal segments; if no arm has an adequate instrument, the causal bar is
*inconclusive*, not *failed*. Orthogonal-layer coverage (CPTAC) is partial — pre-register the gene set
with adequate protein coverage.

---

## 4. Model specification — recommended default (and variants)

Every symbol defined. This is precise enough to implement.

### 4.1 Indices, data, gauge
- `s ∈ {1..S}` samples (tumour `S≈1,060` — primary tumour with PAM50 **Normal-like excluded**, `n=36`
  dropped via `config.EXCLUDE_NORMAL_LIKE`, from the ~1,096 primary-tumour miRNA set; NAT `104`; GTEx
  `346`; metastatic `7` — **dropped**, too few).
- `m ∈ A_P` candidate arms for program `P`; `g ∈ G_P` genes; `F(m)` = seed family of `m`.
- `X_{s,m} = log2(RPM_{s,m}+1)` — **the fixed abundance gauge** (Decision C); no per-arm variance
  standardisation for the magnitude/attribution model.
- `Y_{s,g}` = primary readout = tumour mRNA `log2(TPM+1)` (Decision B).
- `C_{s,·}` = **minimal-core** confounders (Decision B): **purity = CPE/ESTIMATE** (decided 2026-07-04;
  not the expression-derived Cancer-Epithelial fraction), target-gene CN, **one malignant-compartment
  proliferation** (Cancer-Epithelial CIBERSORTx HiRes GEP; the E2F/G2M/MKI67 metagenes and
  `thornsson_proliferation` as *sensitivity only* — not stacked), and a **transcription proxy**
  (intronic/pre-mRNA or TF-regulon or latent factors). **Program-conditional** add-ons for stroma/immune
  programs: the **breast single-cell-atlas deconvolution fractions** (CIBERSORTx Cancer/Normal Epithelial,
  CAFs, PVL, Endothelial, immune; or BayesPrism) — gold-standard vs the metagene (retest). **HRD and batch
  are off by default**
  (near-redundant / single-platform; include only if an ablation or batch diagnostic justifies them).

### 4.2 Occupancy link (abundance → repression), gauge-anchored
For arm `m`, gene `g`, sample `s`, with per-site effective dissociation `κ_{g,m}` (low = strong site):
```
occ_{m,g,s} = a_{m,s} / (a_{m,s} + κ_{g,m}),      a_{m,s} = (X_{s,m} − x0)_+      # threshold-gated abundance
```
`x0` = functional-threshold floor on the RPM scale **[assert: Mullokandov 2012]**, tuned out-of-fold.
Unsaturated approximation (`a ≪ κ`, the MVP): `occ ≈ a_{m,s}·θ_{g,m}`, `θ_{g,m} = 1/κ_{g,m}` = the
learned non-negative weight `M_{g,m}`.

Shared free-AGO pool (Decision G; ceRNA + promiscuity), per sample:
```
φ_s = A_total,s / (1 + Σ_{m,g} occ_{m,g,s}),      w_eff(m,g,s) = φ_s · occ_{m,g,s}
```
`A_total,s` anchored to the RISC-capacity covariate (`ago_gate`).

### 4.3 Repression aggregate and observation
```
ρ_{g,s} = Σ_{m ∈ R(g)} w_eff(m,g,s)  +  Σ_{(m,m') ∈ coop(g)} ψ_{g,mm'} · occ_{m,g,s} · occ_{m',g,s}
Y_{s,g} = α_g + C_{s,·} β_g − ρ_{g,s} + ε_{s,g},      ε ~ N(0, σ_g²)
```
`β_g` (confounder block) **unpenalised**; `coop(g)` = arm-pairs with distinct sites 8–40 nt apart.

### 4.4 Weights, priors, groups (the learned `M`)
```
θ_{g,m} = z_{g,m} · η_{g,m},        η_{g,m} ~ N⁺(μ_{g,m}, τ_{g,m}²)          # non-negative magnitude
z_{g,m} ~ Bernoulli(π_{g,m})                                                  # spike-and-slab inclusion
η_{g,m} = η_{g,F(m)} + δ_{g,m},     δ_{g,m} ~ N(0, ρ_F²)                       # seed-family partial pooling
ψ_{g,mm'} ~ z^coop_{g,mm'}·N⁺(μ^coop, τ²),   z^coop ~ Bernoulli(π^coop)       # primed cooperativity, low base rate
```
Priors from evidence fusion (Decision E), **from sequence/curation/external only**:
```
μ_{g,m}     = log(1 + Σ_{sites} occ_prior(K_site))                            # site-level magnitude (context++/K_D)
logit π_{g,m} = w0 + w_tier·tier + w_meth·TarBase9_method_conf(PMID-deduped)
                   + w_ago·log1p(AGO-CLIP: ENCORI∪POSTAR3∪TarBase-CLIP, PMID-deduped) + w_chim·log1p(chimeric: Manakov∪TarBase-qCLASH)
                   + w_lit·log1p(breast PMIDs) + w_seq·seq_ind                # NB TarBase regulation ≈ all-Negative → not a direction term
τ_{g,m}     = τ0 · decreasing(evidence_strength)                              # HE tight, orphan loose
sign        : from functional assays + M≥0 prior ; only the ~301 TarBase Positive edges → activating/excluded class
```
Source weights `w_*` **calibrated** against held-out functional repression (not hand-set).

### 4.5 Estimation
- **Gold:** mean-field VI or Gibbs over `(z, η, δ, ψ, β, σ)`; report **posterior inclusion probabilities**
  (discovery score) and posterior widths (family non-identifiability).
- **Surrogate / MVP:** non-negative **adaptive group-elastic-net**: minimise
  `Σ_{g,s}(Y − α − Cβ + ρ)² + Σ_{g,m} λ_{g,m}|θ_{g,m}| + group penalties`, with `λ_{g,m}=λ/w_prior(g,m)`
  (adaptive from `μ`), seed-family groups, `θ≥0`, `β` unpenalised. Tune `λ`, `x0`, pool params by **nested
  patient-level CV**, never in-sample.

### 4.6 Where it nests the current heuristic (§0 sanity anchor)
Freeze `z=1` on the HE support, `θ_{g,m}=edge_w` (curated, fixed), linear occupancy (`occ=a·θ`), `φ=1`
(no pool), no cooperative terms, `C=∅`, abundance transform = `expr_mult` → then
`ρ_{g,s} = Σ_m expr_mult(m,s)·edge_w(m,g)` = **exactly the current pressure** (`FORMULAS.md` §1). So the
heuristic is the **infinite-prior-strength, no-confounder, linear-occupancy, hard-HE limit**. Anything
learned reads as "the prior plus what the data insists on."

### 4.7 Variants kept in play (specified enough to tell apart)
| Variant | Change from default | When preferred |
|---|---|---|
| **Gene-focused adaptive-lasso** | one penalised regression per gene, no cross-gene pooling | MVP; deep single-driver dives; direct compare to per-edge coupling |
| **Low-rank RRR/NMF overlay** | `M ≈ UVᵀ`, rank `r` | explain residual covariance; comodule roll-up; borrow max strength |
| **Feature/FM prior generator** | `μ_{g,m}=φ(g)ᵀΘψ(m)` | fill `μ` for non-canonical orphans with no K_D |
| **Y = protein / discordance** | swap readout; MoPC-style conditioning | function-first genes; translational-class edges; OOD validation |
| **State: separate / differenced** | fit `M_H,M_T` apart, or `ΔY=−MΔx` | audit; platform-constant subtype-centroid contrast |
| **Subtype environments (ICP)** | invariance across PAM50 as a causal filter | causal robustness check; flag subtype-specific edges |

---

## 5. Evaluation & falsification plan

Pre-named tests, references, success criteria, and abandon conditions. All coupling tests are
**partial-Spearman adjusting `C`**, on **held-out patients** (patient-level CV folds; the prior never
touches the tested `(X,Y)`).

**Bar 1 — beat raw abundance (primary, §0 ref 1).** Does the learned **realised aggregate** `ρ_{g,s}`
anti-correlate with held-out target mRNA better than **raw arm abundance** (and than family-summed
abundance)? *Metric:* paired ΔSpearman across folds, partial on `C`. *Success (pre-named):* median
ΔρSpearman > 0 with 95% CI excluding 0 across the §1h driver panel **and** for ≥ a pre-registered
fraction of program genes. **Abandon-to-abundance** if the learned weights add nothing beyond abundance
here — this is the decisive bar.

**Bar 2 — match/extend curated fixed-M (§0 ref 2).** Learned-M vs curated-M on the same OOF coupling and
on **OOF prediction** of `Y` (because the coupling-invariance lemma says single-edge coupling can't
separate weightings — test at the aggregate/predictive level). *Success:* learned ≥ curated (non-
inferiority) while **nesting** it, **or** learned discovers ≥1 OOD-validated orphan edge curated-M lacks.
Beating is *not* required — parsimonious reproduction counts.

**Bar 3 — prior does real work (§0 ref 3).** Permute the evidence→(π,μ) map (shuffle within degree
strata), refit. *Success:* real prior beats shuffled OOF. Guards against "any sparse structure would do."

**Bar 4 — causal, not covariational (the differentiator).** miRNA-locus **CN instrument** 2SLS:
`x_m = γ·CN_locus(m) + Cδ + u`; `y_g = −M·x̂_m + Cβ + ε`, on **focal** segments (pleiotropy control).
*Report* first-stage F (>10 = not weak). *Success:* causal `M` concordant in sign/rough magnitude with the
conditioned `M` for arms with adequate instruments. *Inconclusive* (not failed) where instruments are
weak.

**Bar 5 — orthogonal layer & independent cohort (anti-circularity crux).** Fit on tumour mRNA; score the
learned edges, unrefit, on **CPTAC protein** and **protein-vs-mRNA discordance** (MoPC-style
partial-corr conditioned on the arm) — the translational residual the mRNA fit never saw; and on an
**independent mRNA cohort** (SCAN-B / METABRIC / Buffa). *Success:* learned edges show protein-level
repression **and** discordance signal beyond curated-M; independent-cohort coupling replicates in
direction. **Manakov chimeric-eCLIP** physical duplex as a binding check on discovered orphans.

**Face validity (not a quantitative bar).** Recover the §1h canonical edges: let-7→KRAS/MYC/HMGA2,
miR-200→ZEB1/2, miR-29→collagens/DNMT3, miR-21→PTEN/PDCD4, miR-17~92→E2F1/PTEN/BIM, miR-34→CDK4/6/BCL2,
miR-221/222→CDKN1B/ESR1. Failure to recover well-instrumented canonical edges is a red flag on the fit.

**Circularity guardrail (non-negotiable).** `π,μ,τ` and interaction support from sequence/curation/
external only; expression estimates magnitudes; never test coupling on the fit's `(X,Y)`; keep learned-M
and curated-M side by side and report disagreements as findings.

---

## 6. Phased build plan (cheapest, most diagnostic first; each phase gated on a pre-named win)

> **Companion:** `LEARNED_MODEL_BUILD_PLAN.md` maps these phases to the repo — **what to reuse vs rebuild**
> (the data/infra spine + eval machinery are reused; the model core + evidence ledger are new; the current
> heuristic is frozen as the baseline; ~100 one-off scripts are archived) and a **target `learned/` package
> architecture** with a low-risk *classify-first, build-alongside, archive-lazily* migration.

**Phase 0 — Observation-model audit (cheapest, highest-value).** On the existing per-edge coupling,
add the **transcription-control block** to `C` and measure how much canonical edges move (Decision B
ablation). *Gate:* the transcription confound is quantified; if adding the proxy flips/kills spurious
positive weights on known repressive edges, the block is mandatory downstream. No new estimator yet.

**Phase 1 — Gene-focused adaptive-lasso MVP.** On a hub panel (PTEN, GATA3, ESR1, ZEB1, IRF1): non-
negative adaptive-lasso, `λ_{g,m}=λ/w_prior`, seed-family groups, `C` (with transcription block)
unpenalised, occupancy in the unsaturated (linear) approximation. *Gate:* recovers curated drivers **and
beats raw-abundance coupling out-of-fold** (Bars 1–2 on the panel).

**Phase 2 — Evidence-fusion prior + soft discovery. Build the provenance ledger FIRST.** Before any
weight-tuning, construct the unified **`(m,g,PMID,assay-class)` ledger** (miRTarBase + TarBase-v9 + ENCORI
+ POSTAR3 + Manakov) and **de-duplicate globally by PMID×assay** — physical *and* functional — replacing
the current "miRTarBase counts + additive ENCORI" scoring. Then wire biochemical K_D / context++ into
`(π,μ,τ)` and run confirmatory (hard-HE) and discovery (soft) side by side. *Gate:* dedup materially changes
≥ some edges' evidence (sanity: hub CLIP edges collapse toward distinct-PMID counts); Bar 3 (beat
shuffled-evidence null); ≥1 discovery edge that clears Bar 5 (protein/independent-cohort).

**Phase 3 — Program-wise hierarchical model.** Full default (§4) on 3–5 Hallmarks (EMT, P53, G2M, IFNγ,
HYPOXIA): family partial pooling, occupancy link + threshold gauge, shared pool. *Gate:* family→gene
estimands stable across folds; retire `D(m)` (pool reproduces promiscuity penalty); OOF coupling ≥
curated-M at the aggregate level.

**Phase 4 — Cooperativity + ceRNA terms.** Add primed product terms (site-coord window) and the shared
free-AGO pool depletion. *Gate:* any product/pool term **improves held-out fit** and survives shrinkage;
otherwise drop (additive is fine).

**Phase 5 — Causal instrument.** CN-locus 2SLS on the hub panel and program genes with adequate
instruments. *Gate:* first-stage F reported; causal-vs-conditioned concordance for well-instrumented arms
(Bar 4).

**Phase 6 — State axis.** Nested `M_T = M_H + Δ` with GTEx prior; test `Δ=0`. *Gate:* `Δ` concentrates in
APA-shortened / capacity-shifted genes (biology signature) not in GTEx-imputed arms (artifact signature).

**Phase 7 — Universe overlay (only if 1–6 pay off).** Low-rank RRR/NMF for structural roll-up and
comodules; the feature/FM model for non-canonical orphan `μ`. *Gate:* comodules interpretable; couple OOF.

Each gate mirrors `PRESSURE_FUTURE_OPTIONS.md` discipline: promote a lane only on a pre-named held-out /
orthogonal-layer win, never on in-sample fit.

---

## 7. Open questions & what would change the answer

1. **Is `M` actually state-dependent, or only `x`?** (Decision H.) The nested model tests `Δ=0`. If `Δ≈0`,
   the whole state story collapses to acquired abundance — simpler, and it retires the separate/differenced
   variants. **What would resolve it:** wiring **APA/PDUI** (Xia tumour, APAatlas normal) so site loss is
   an explicit state covariate; if `Δ` tracks PDUI shortening, `M` is genuinely state-dependent.
2. **How much does isomiR resolution change the predictor unit?** 5′ isomiRs shift the **seed** (positions
   2–7) → different family, different targets. This could break some family collinearity (Decision F) and
   correct the prior's seed. **What would resolve it:** run the arm-level model, then re-run with isomiR-
   resolved arms on the arms with high 5′ heterogeneity; if edges move, isomiR is first-class, not a
   detail. *Inferred* to matter for a minority of arms **[infer]**.
3. **Does TarBase-v9's direction label change which edges are admissible? — largely answered: no.** I
   checked the live file: only **301 of 4.72M rows are "Positive"** (the field is an almost-constant
   "Negative" default, because ~93% of rows are CLIP = binding, not direction). So the activating-edge
   exclusion is a **tiny** effect, not a major prior channel; the non-negativity constraint is fine for
   almost all edges. The *open* residual: are those ~301 Positives real activating interactions worth a
   dedicated class, and does miRTarBase/breast-context add more? **What would resolve it:** map the 301
   Positives (and any curated up-regulating edges) onto the Hallmark panel and inspect.
4. **Which observation target wins (Decision B) ripples everywhere.** If protein/discordance turns out to
   be the only clean readout (transcription confound intractable in mRNA), then: primary `Y` shifts to
   CPTAC, `n` drops, the program-wise pooling (Decision A) becomes *more* important (borrowing strength
   under small `n`), and the state axis narrows to states with CPTAC coverage. Conversely if the
   transcription proxy fully de-confounds mRNA, mRNA stays primary and CPTAC is purely validation.
5. **Is the a·w threshold `x0` estimable, or must it be asserted?** (Decision C.) If OOF tuning finds a
   stable `x0`, the gauge is data-anchored; if not, it stays a biological assertion (Mullokandov) with a
   sensitivity band. **What would resolve it:** the CN instrument — a dose perturbation that crosses the
   threshold would reveal the nonlinearity directly.
6. **Does the shared-pool term truly retire `D(m)`?** (Decision G.) **What would resolve it:** fit with the
   pool and *without* `D(m)`; check that hub-arm per-target occupancy is depressed to the degree `D(m)`
   hand-imposed. If yes, `D(m)` is gone; if the pool is under-identified, `D(m)` stays as a fallback.
7. **Family vs arm as the shipped estimand.** If per-arm splits never identify (δ posteriors stay at
   prior), the honest deliverable is family→gene throughout, and the per-arm "driver nomination" language
   must be enforced repo-wide. **What would resolve it:** the posterior widths of `δ` across the driver
   panel.

**Things in the brief I would flag as needing correction or sharpening** (§8 item 7):
- **Chimeric eCLIP is site-level — but rely on AGO-CLIP as the *primary* CLIP layer** (feedback,
  2026-07-04). The brief's §2c lists Manakov under edge-level ("physical proof a pairing occurs"), but
  chimeric reads carry **duplex coordinates** — so they can **prime cooperative interactions** (Decision
  G) *and assign the arm*, not just gate inclusion. **However, chimeric is sparse and its reproducibility
  is contested**, so the design leans on **AGO-CLIP (ENCORI) read depth** — dense, replicated, already
  in-repo — as the primary physical `π` signal, with chimeric a *confirmatory* arm-assignment / coordinate
  layer used only where dense enough, falling back to AGO-CLIP region + seed inference when chimeric is
  too sparse (Decision E). Their key difference: **AGO-CLIP is miRNA-agnostic** (needs seed inference);
  **chimeric gives the direct (arm, site) pair**. The design underused AGO-CLIP in its first draft — this
  corrects it.
- **"Site exists" is a weak prior; graded K_D is the strong one.** The flanking-context ~100-fold affinity
  range **[assert: McGeary 2019]** means a binary site-type mask under-uses the sequence evidence. The
  magnitude prior should be the **continuous K_D/occupancy**, not the site-type tier.
- **The current heuristic's biggest defect is the missing transcription control**, not the softmax or
  `D(m)` per se. Fixing the observation model (Decision B) likely moves more edges than any weight change.
- **Composition adjustment is large-effect, not cosmetic — and the metagene under-corrects it** (feedback
  + repo `core_coupling_deconv_retest`). Tissue-matched **deconvolution compartment fractions** (the
  user's CIBERSORTx breast single-cell-atlas run; BayesPrism) are the gold-standard composition covariate;
  the concrete cost of getting this wrong is the **miR-29b→ECM/collagen** edges collapsing from −0.3/−0.4
  to ≈0 once CAF/stromal fraction is deconvolved. Use compartment fractions for stroma/immune programs and
  a **malignant-compartment (Cancer-Epithelial) proliferation**, not a bulk metagene (Decision B).
- **TarBase v9 is mostly CLIP, not a fresh directional layer** (feedback + live-file check). ~93% of its
  4.72M rows are HITS/PAR-CLIP (+ qCLASH chimeric) from **998 PMIDs**, and `regulation` is
  **4,724,236 Negative vs 301 Positive** — so it double-counts ENCORI/POSTAR3 CLIP and provides almost no
  direction. Fold it into the shared PMID-deduped CLIP/chimeric axes; its unique value is site coordinates,
  confidence, `microt_score`, and a thin functional layer — **not** direction or extra binding mass.
- **The current scores are NOT PMID-deduplicated across databases — and this must be fixed globally**
  (feedback). `evidence_score` today = miRTarBase per-class *study counts* with ENCORI CLIP depth **added
  on top** (no dedup), and the model's edges table has already collapsed to counts (no PMID column). So the
  same paper is counted twice both as **CLIP** (ENCORI/POSTAR3/TarBase) *and* as **functional** (a
  luciferase/western paper in both miRTarBase and TarBase). Fix = a `(m,g,PMID,assay-class)` ledger built at
  edge-build time, deduped by PMID×assay (physical **and** functional), counting distinct PMIDs per class
  (Decision E). Sequence/K_D priors are PMID-free and stay a separate channel.
- **Use several AGO-CLIP sources, not just ENCORI** (feedback). ENCORI AGO-CLIP is sparse (~3% site
  coverage in the site-ladder); **POSTAR3 carries ~2.36M AGO2 peaks** in-repo and is the dense site-level
  representative. Keep them non-double-counted (same assay class): ENCORI edge-level, POSTAR3 site-level.
- **Metastatic `n=7`** is too small to model; drop rather than pool.
- **The AGO gate is currently non-competitive** (cohort saturating, no pool depletion); ceRNA/promiscuity
  need the *shared-pool* form (Decision G).

---

## 8. References (refreshed July 2026)

Mechanism & magnitude — Guo, Ingolia, Weissman, Bartel (2010) *Nature* "Mammalian microRNAs
predominantly act to decrease target mRNA levels"; Eichhorn et al. (2014) *Mol Cell* "mRNA
destabilization is the dominant effect… by the time substantial repression ensues"
(https://pmc.ncbi.nlm.nih.gov/articles/PMC4292926/); Selbach et al. (2008) & Baek et al. (2008) *Nature*
pSILAC/proteomics.

Site efficacy & 3′UTR — Agarwal et al. (2015) *eLife* TargetScan7 context++
(https://www.targetscan.org/); McGeary, Lin, Shi, Bartel et al. (2019) *Science* "The biochemical basis
of microRNA targeting efficacy" (https://pmc.ncbi.nlm.nih.gov/articles/PMC7051167/); Soutschek et al.
(2022) *Bioinformatics* scanMiR (https://doi.org/10.1093/bioinformatics/btac110).

Cooperativity & competition — Saetrom et al. (2007) *NAR* distance constraints
(https://academic.oup.com/nar/article/35/7/2333/1094167); Grimson et al. (2007) *Mol Cell* determinants
beyond seed; Broderick et al. (2011) *RNA*; Rinck/Sætrom-line (2013) *RNA* transcriptome enriched for
cooperativity-permitting distance (https://pubmed.ncbi.nlm.nih.gov/23696004/); Denzler et al. (2014,
2016) *Mol Cell* ceRNA thresholds (https://pubmed.ncbi.nlm.nih.gov/24793693/,
https://pmc.ncbi.nlm.nih.gov/articles/PMC5048918/); Bosson et al. (2014) *Mol Cell* AGO-iCLIP
susceptibility (https://www.sciencedirect.com/science/article/pii/S1097276514007503); Mullokandov et al.
(2012) *Nat Methods* sensor/decoy functional threshold.

Expression modelling — Huang et al. GenMiR++ (2007); Muniategui et al. TaLasso
(https://www.researchgate.net/publication/51480609); Zhang et al. (2011) *Bioinformatics* SNMNMF
(https://dx.doi.org/10.1093/bioinformatics/btr206); PIMiM; Bayesian sparse reduced-rank regression
(https://pmc.ncbi.nlm.nih.gov/articles/PMC5628626/, https://pmc.ncbi.nlm.nih.gov/articles/PMC7333421/);
MoPC (2024) mRNA–protein partial correlation conditioned on miRNA
(https://pmc.ncbi.nlm.nih.gov/articles/PMC10956802/).

Causal / identification — Meng et al. (2019) *BMC Bioinformatics* invariant causal prediction with PAM50
environments (https://pmc.ncbi.nlm.nih.gov/articles/PMC6419852/); miRNA Mendelian randomisation (various,
2024–2025, e.g. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12885728/).

Protein & discordance — CPTAC breast proteogenomics (Mertins et al. 2016; Krug et al. 2020, as reviewed
https://tcr.amegroups.org/article/view/10135/html); MoPC (above).

State / field effect — epigenetic field cancerisation in breast (medRxiv 19002014); field cancerisation
review (https://pmc.ncbi.nlm.nih.gov/articles/PMC9322418/); mutational landscape of adjacent normal
breast (*Cell Reports Medicine* 2025).

Breast / Hallmark miRNA — subtype-specific miRNA & AGO2-PAR-CLIP subtype target regulation
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053773/); subtype miRNA networks review
(https://www.mdpi.com/2227-9059/10/3/651).

CLIP / chimeric — Helwak et al. (2013) *Cell* CLASH (https://pmc.ncbi.nlm.nih.gov/articles/PMC3650559/);
chimeric eCLIP / miR-eCLIP (bioRxiv 2022.02.13.480296; eclipsebio); mirTarCLASH
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11971479/).

Evaluation & overfitting — Pinzón et al. (2017) *Genome Res* "microRNA target prediction programs predict
many false positives" (https://pmc.ncbi.nlm.nih.gov/articles/PMC5287229/); Kern et al. (2021) *NAR*
validation of miRNA target pathways (https://academic.oup.com/nar/article/49/1/127/6030235).

---

*Prepared as the §9 deliverable for `LEARNED_MODEL_DISCUSSION_PROMPT.md`. Reasoning-first: every
decision carries its mechanism, its grounding in biology/data/identifiability, the alternatives kept in
play with when each wins, and its failure modes. The recommended default nests the current heuristic and
is falsifiable against raw abundance, curated-M, and shuffled-evidence nulls.*
