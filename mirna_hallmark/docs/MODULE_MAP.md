# mirna_hallmark — module map (post-reorg 2026-07-04, waves 1+2)

Navigation for the subproject after the sprawl reorg. **Top-level went 160 → 36 `.py`** (35 modules +
`__init__`). All one-off analyses/figures and their theme-base utilities moved into a themed `analyses/`
+ `eval/` tree; **the top level now holds only the data/infra spine, the frozen baseline heuristic, and
`run_all`'s live pipeline steps.**

- **Wave 1** moved the 74 in-degree-0 *leaves* (zero import rewrites — nothing imported them).
- **Wave 2** moved 50 *theme-base* modules (imported by the wave-1 analyses), **rewriting 234 import sites
  across 106 files**. Verified: `py_compile` all OK; no dangling old-path refs; spine + `run_all` + moved
  bases (incl. `outcome_survival` ×17 importers, `cross_state_landscape` ×19) + their importers import OK.

> The moved files were git-**untracked** (this subproject was never committed except `method_dev/`), so
> moves were filesystem-level; `git add` to put the new tree in history. A pre-wave-2 backup of all `.py`
> is at `scratchpad/mh_backup_prewave2.tgz` (session-local).

---

## The tree

```
mirna_hallmark/
  <36 top-level .py: SPINE + BASELINE + pipeline + run_all — see below>
  analyses/                    # finished one-off analyses + their theme bases (off the build path)
    dcis_ev/      (20)         # DCIS + EV arc: loaders, CAF/NF, mir29c, ev_mirna_* screens/deep/followup
    outcome/      (13)         # survival/prognostic: outcome_survival base + outcome_* + prognostic
    cnv_locus/    (12)         # miRNA-locus CN/SV, CCLE, dosage attribution
    pressure_dev/ (12)         # pressure heuristic sensitivity/comparison/validation/prognostic
    spatial/      (10)         # spatial_common base + visium/xenium/mibi/deconv
    edge_panels/   (7)         # edge pressure/corr panels + edge_prior_refinement + edge_breast_context
    cptac/         (6)         # CPTAC one-offs (orphan discovery/confound, acquired pressure, …)
    evidence_dev/  (5)         # ENCORI/evidence-scoring sensitivity + comparison + fetch
    cross_state/   (4)         # cross_state_landscape/coupling/deep_dive/expression_panels bases
    nmf/           (4)         # nmf_* + within_* signatures
    builders/      (3)         # _build_depmap_dependency / _build_tf_census / _build_meth_chunked (cache)
    mir301/        (2)         # mir301_family_depth / mir301_focus_genes
    figures/       (2)         # make_grant_figure / make_seed_grant_figures
    ops/           (2)         # rerun_* batch scripts
    architecture/  (1)         # geneset_architecture
    misc/         (14)         # cross-theme one-offs (comovement, promiscuity, methylation, escape, …)
  eval/            (7)         # reused OOD/validation: coupling_permutation, held_out_tuning,
                               #   core_coupling_{deconv,composition}_retest, buffa/cptac_validation,
                               #   targetscan_orphan_coupling
  method_dev/                  # (unchanged — already an organized subpackage)
```

Run any moved module by its new dotted path, e.g. `python -m mirna_hallmark.analyses.outcome.outcome_advanced`.

---

## The 35 modules that stayed top-level

Rule: a module stays only if it is **shared infra the learned model / eval will reuse**, the **frozen
baseline heuristic**, or a **live `run_all` pipeline step**. Nothing else.

**SPINE — data/infra, REUSE (21).** `config` `data_loaders` `stats` `hallmark_sets` `mirna_arm_resolve`
`healthy_anchor` `gtex_mirna_matrix` `arm_expression` `estimate_scores` `brca_deconvolution`
`mirna_locus_cnv` `gene_roles` `robustness_checks` `coupling_inference` `build_edges` `evidence_scoring`
`encori_edges` `seed_family_justification` `family_normal_reference` `tcga_batch` `cptac_batch`.
  *Evidence-refactor targets* (→ the PMID-deduped method-centric ledger, BUILD_PLAN §1b): `build_edges`
  (emit ledger not counts), `evidence_scoring` (superseded by learned fusion), `encori_edges`.

**BASELINE — frozen heuristic, do NOT extend (8).** `pressure_build` `pressure_engine` `ago_gate`
`hallmark_interaction` `mirna_state_class` `hybrid_pressure` `dual_spine_comparison`
`pressure_signature_curated_info`. The learned model must *nest and beat/match* these.

**PIPELINE — live `run_all` steps (5) + orchestrator.** `coupling_predictor_comparison`
`pam50_gene_resolution` `stratum_characterization` `subtype_contrasts` `visibility_archetype_contrasts`;
and **`run_all`** (the only true full-graph leaf; drives the pipeline).

---

## The new model
- **`learned/` package** (BUILD_PLAN §3) — **scaffolded 2026-07-04**, imports the SPINE in place.
  `data.py` (assemble Y/X/C/w_prior from spine loaders), `regression.py` (gene-focused NN adaptive-lasso,
  confounder-residualised, evidence-weighted penalty), `mvp.py` (Phase-1 OOF gate vs raw abundance);
  documented stubs for `priors`/`occupancy`/`families`/`cooperativity`/`state`/`attribution`/`evidence/`/`eval/`.
  **Phase-1 result:** `python -m mirna_hallmark.learned.mvp` — learned aggregate beats raw-abundance
  coupling OOF on **4/5** hub genes (PTEN, GATA3, ESR1, ZEB1; not CDKN1A), recovering canonical drivers
  (miR-200a/c→ZEB1, miR-206/miR-18a→ESR1, miR-106b/miR-29c→PTEN). Passes the Design §6 Phase-1 gate.
  **Phase-2 (evidence ledger) IMPLEMENTED:** `learned/evidence/{methods,ledger}.py` builds the
  PMID-deduped, method-centric ledger from miRTarBase + TarBase v9. `python -m
  mirna_hallmark.learned.evidence.ledger` → **1.92M raw (m,g,pmid,class) rows → 1.11M deduped (42%
  collapsed** — same paper in both DBs, functional tiers too), 6,313 distinct PMIDs; top fused edges
  canonical (miR-21→PDCD4/PTEN, miR-200c→ZEB1, miR-221→CDKN1B). Wired into the MVP
  (`--prior ledger`): holds the 4/5 gate and **improves driver recovery** (miR-21→PTEN surfaces).
  **Phase-3a (seed-family pooling) IMPLEMENTED:** `families.py` collapses regulators to TargetScan seed
  families (edge table `miRNA_family`; correctly splits miR-200a/141 vs miR-200b/c/429).
  `python -m mirna_hallmark.learned.mvp --prior ledger --family` → gate **4/5 → 5/5** (CDKN1A now passes —
  its collinear miR-15/16 & miR-17~92 families were split arbitrarily per-arm), stability rises where
  collinearity lives (ZEB1 0.70→0.81); the family→gene estimand is canonical (ZEB1 = the two miR-200
  subfamilies; CDKN1A = miR-15/16/195/424/497 + miR-17/20/93/106).
  **Composition confounder (Design §B) WIRED:** `data.assemble_gene(deconv=True)` / `mvp --deconv` adds the
  user's CIBERSORTx non-malignant compartment fractions (CAFs, immune, endothelial, PVL, normal-epi — NOT
  Cancer Epithelial, which would over-control) to `C`. Demonstration: coupling attenuates most on
  stroma-expressed miR-29→ECM targets (COL3A1 −0.35→−0.02, COL1A1 −0.42→−0.14, SPARC −0.31→−0.07) and
  barely on cell-intrinsic ESR1 (−0.38→−0.32) — the composition confound, separated (reproduces the repo
  deconvolution retest). Confirms Design §B "biggest hole" is real and now addressed.
  **CN-locus instrument (Design §5 Bar 4 — the causal test) IMPLEMENTED:** `instrument.py` (reuses
  `mirna_locus_cnv` ASCAT3 locus CN). Per edge: first-stage F (CN→arm), reduced-form (CN→target), vs
  observational coupling. `python -m mirna_hallmark.learned.instrument` → all 12 tested edges have valid
  instruments (F 24–174, first-stage ρ>0), **7/12 causal-concordant** (incl. canonical miR-17~92→BCL2/MEF2D,
  miR-320a→TFRC, let-7b→CDC34); flags miR-17→CDKN1A/RUNX1 as **possibly confounded** (obs<0 but reduced≈0) —
  the confounded-vs-causal discrimination only an instrument gives.
  **External data fetched (Tier 1, 2026-07-04) → `data/external_cache/`:** CollecTRI TF-regulon (64k signed
  edges + n_references strength) → `data._tf_regulon/_tf_proxy` STRENGTH-weighted transcription proxy
  (PTEN −0.037→−0.055 = instrument-causal level; ESR1 less over-correction); TargetScan v8.0 context++
  (1.4M sites, coords) → `data._targetscan_context` orphan magnitude, `assemble_gene(orphans=True)`
  **DISCOVERS miR-148a-3p→ESR1** (no curated edge, coupling −0.17); scanMiR K_D (2656 human KdModels, installed,
  for occupancy κ + orphan μ — not yet wired).
  **Evidence-weight calibration (does the fusion learn its hyperparameters non-circularly?) IMPLEMENTED:**
  `learned/evidence/calibrate.py` regresses the McGeary-2019 transfection readout (GSE140217, 17 miRNAs ×
  HeLa, batch-norm logTPM; leave-one-out mRNA repression) on the ledger's per-assay-class distinct-PMID
  counts, joined per **seed family** × gene (bridges: UCSC refGene NM_→symbol 99.7%; mirGeneDB-confirmed
  guide arms). `python -m mirna_hallmark.learned.evidence.calibrate` → (A) **curated edges repress more than
  non-curated** (within-family, p=2e-18) — curation captures function; (B) **functional-assay evidence
  predicts mRNA repression ~3× better than binding-only** (ρ 0.21 vs 0.08), confirming the coarse
  direct-functional > CLIP/chimeric ordering *non-circularly*, and the directness signal survives
  controlling total evidence (amount≠directness); (C) the fused `ledger_weight` tracks repression better
  than a raw PMID count (ρ 0.26 vs 0.20 — weighting is worth it). The *fine* order is readout-specific: on
  an mRNA readout qpcr_rna is the single best predictor and protein-only (western) under-predicts, exactly
  as expected → task-matched `methods.CLASS_WEIGHT_MRNA` (selectable `w_prior_source='ledger_mrna'`). **The
  OOF gate is robust to hand-set vs calibrated weights** (PTEN −0.349 vs −0.348, gate 5/5 both; the adaptive
  lasso uses the prior for ordering, which is stable) — so the hyperparameters are both data-groundable AND
  non-critical. Caveat: HeLa is heterologous, so this calibrates assay-*directness* (cell-agnostic), not
  breast target sets.
  **scanMiR biochemical K_D → occupancy κ + affinity prior IMPLEMENTED (2026-07-04):** `learned/kd.py`
  + `learned/scanmir_scan.R` (getKdModels("hsa") 2656 KdModels → findSeedMatches → aggregateMatches over
  MANE 3'UTRs of the 1432 HE genes [reuses `site_ladder.utr_seed_scan._gene_3utr`; no fetch] × 771 HE arms).
  `python -m mirna_hallmark.learned.kd` → **912,718 per-(arm,gene) predicted-repression rows** (biochemical
  logFC; miR-200c→ZEB1 −2.25, miR-21→PTEN −0.98 = canonical). Two wirings tested: (1) **affinity-aware
  per-edge κ** for the occupancy transform (`kd.edge_kappa`, arm+family; `mvp --occ --kd`): κ(m,g)=κ0·2^(β·rep)
  so strong edges saturate sooner. **NEGATIVE RESULT** — occ (naive *or* kd-κ) mostly *hurts* coupling vs the
  linear baseline (PTEN family −0.35→−0.11, ZEB1 −0.34→−0.12), kd-κ worse than naive, because saturating
  removes the *cross-sample abundance variance* the TCGA coupling lives on (the signal is sub-saturation,
  per Denzler/Bosson thresholds). occ's real effect is `abund_dom`→0 = retiring abundance-domination
  (attribution axis, ⊥ coupling) — the biochemistry is correct but for *within-sample* occupancy not
  *cross-sample* coupling; only the weak-coupling CDKN1A improves. (2) **affinity as adaptive prior**
  (`w_prior_source='scanmir'`): coupling ≈ ledger (wash) but **best selection stability yet, 0.90 vs 0.85**,
  and lowers PTEN abund_dom 0.23→0.14 — a principled continuous prior stabilizes family selection. So
  scanMiR K_D's productive slot is prior/orphan-magnitude/independent-M-check, NOT occ κ. Defaults unchanged
  (ledger+family canonical); scanMiR is a validated complement. Cache data/external_cache/scanmir/.
  **Fused prior (ledger existence + scanMiR magnitude) — tested, NON-ADDITIVE (2026-07-04):**
  `w_prior_source='fused'` (standardize each prior among a gene's regulators, `exp(½(z_ledger+z_scanmir))`).
  On the curated HE-only set it's *worse* than either component (gate 5/5→4/5, stability 0.86 < scanMiR's
  0.90) — because every HE regulator is already existence-supported, so ledger and scanMiR are redundant
  *magnitude* proxies there and any admixture breaks weak-coupling CDKN1A's ~0 abundance gate. The three-roles
  fusion's real slot is existence-gated *orphan* magnitude (where existence and magnitude diverge). Kept as a
  documented, un-tuned option; not default.
  **SV instrument for standalone miRNAs — tested, NOT VIABLE; CN already covers it (2026-07-04):** reused the
  SV pipeline (`analyses.cnv_locus.mirna_locus_sv_overlap` → `mirna_locus_sv_hits.tsv.gz`, 163,842 DEL/DUP ×
  hairpin rows). Focal (≤1Mb) + standalone (no cluster) leaves only ~12–21 carriers per locus, and the
  first stage ρ(focal-SV-sign → arm expr) ≈ 0 (miR-21 +0.05) — no power (SVs are subclonal; a heterozygous
  focal event barely moves post-transcriptionally-buffered mature miRNA). The **CN instrument (continuous
  ASCAT, all 946 samples) handles the same standalone arms** — miR-21 F=49, and causally confirms **PDCD4**
  (−0.127) & SPRY2 while flagging PTEN/BCL2/TIMP3 as observational-only. The SV motivation (cluster-pleiotropy)
  never applied to standalone arms. ⇒ SV data is *descriptive* (recurrent SV-disrupted loci), not an
  instrument; CN is the right causal tool. No new module added.
  **Bar 3 — shuffled-evidence null (Design §5 Bar 3) IMPLEMENTED:** `learned/eval/shuffled_null.py`
  (`permute_prior` hook on `mvp.oof_gate`). Two halves: **magnitude (μ)** — shuffle prior weights within a
  gene's HE regulators → **prior does NO coupling work** (Δρ +0.0003, Fisher p=0.39; the adaptive lasso
  fits coefficients to the data regardless of weights); **inclusion (π)** — curated HE arm-set vs random
  same-size arm-sets (uniform weights) → **prior DOES work** (mean Δρ −0.090, Fisher p=0.021; PTEN/ESR1/ZEB1
  Δ −0.13/−0.14/−0.15). ⇒ the finding that explains the session: *which* arms are curated carries coupling;
  *how* they're weighted does not — so every magnitude-prior refinement (ledger/scanMiR/fused/calibration)
  was a wash by construction. `--test inclusion|magnitude`.
  **Phase 6 — state axis Mₜ = M_H + Δ (Design §Decision H / §7 Q1) IMPLEMENTED:** `learned/state.py`. GTEx
  breast mRNA (`data/GTEx/gene_tpm_v10_breast.gct.gz`) + miRNA (`gtex_mirna_matrix`), **327 donor-paired**.
  Cross-state transfer with an **n-matched held-out control** (fit M on a 327-tumour subsample vs healthy
  M_H, both eval on held-out tumour — isolates state from estimation noise). Result: **M is state-dependent
  gene-specifically** — ZEB1 retention 0.93 (miR-200→ZEB1 wiring **conserved** healthy→tumour, only
  abundance shifts — and composition mismatch is a headwind, so this is robust), while PTEN 0.56 / ESR1 0.33
  / CDKN1A 0.65 remodel per-arm weights in tumour (ESR1: miR-22-3p stronger, miR-18-5p flips). Caveat:
  healthy-vs-tumour composition mismatch confounds the *state-dependent* calls (not the conserved ones).
  **Remaining plan:** Phase 3 cross-gene hierarchical Bayes; Phase 4 cooperativity/ceRNA; Decision I Shapley;
  fold ENCORI/POSTAR3/Manakov into ledger; Bar 5 independent-cohort on learned M.
  **Why the magnitude prior doesn't help — diagnosed (2026-07-04):** alpha sweep on the Bar-3 magnitude
  test shows higher shrinkage makes it *worse* (not underweighting); per-edge `Spearman(evidence, |coupling|)
  = +0.22` (weak) and among collinear regulators it washes out. Top- vs bottom-evidence HE edges couple
  about the same — and for BCL2/MYC/CDKN1B the *low*-evidence edges couple better. ⇒ **inclusion (who) ⊥
  magnitude (how much); evidence informs who, data informs how much.** The Bayesian implication: evidence →
  spike-and-slab inclusion π (works), data → magnitude μ, hierarchy → cross-gene pooling.
  **Discovery lane (where the prior pays off) IMPLEMENTED:** `learned/discovery.py` — TargetScan sequence
  inclusion nominates orphans; keep those adding repressive coupling *beyond* the curated aggregate
  (partial on C + HE-aggregate) with scanMiR support, then a built-in **deconv composition control**.
  Panel scan → 5/10 survive; headline **miR-135b-5p→GATA3** (−0.42→−0.36 deconv, fully novel, scanMiR-backed)
  and **miR-181d-5p→ESR1**. Constructive proof of the Bar-3 thesis.
  **miR-135b-5p→GATA3 VALIDATED (clears Phase-2 gate: a discovery clears Bar 5):** de novo TCGA coupling +
  scanMiR K_D −0.68 + deconv-robust −0.36 + **Bar 5 CPTAC protein** (independent prospective: ρ(miR-135b,
  GATA3 protein) −0.49 n=101, discordance −0.23; tcga105 −0.46) + **literature** (Cancer Cell Int 2021,
  10.1186/s12935-021-02015-6: miR-135b-5p↔GATA3 inverse in breast tumours, circ_0044234 ceRNA axis). Manakov
  HEK293T has no miR-135b–GATA3 duplex (miR-135b likely not expressed there — weak negative).
  **Phase 3 — hierarchical Bayes (cross-gene pooling + uncertainty) IMPLEMENTED:** `learned/hierarchical.py`
  (Gibbs sampler, β(m,g)~N(μ_m,τ²)). (1) **cross-gene magnitude pooling = negative** — pooled≈solo at n=987,
  worse at n=150 (Δ−0.036): a miRNA's per-target weight is *target-specific*, so the μ_m hierarchy is
  mis-specified (the right pooling is seed-family, Phase 3a). (2) **uncertainty payoff = positive** — 72/242
  edges identified (|β/sd|>2): miR-17~92→RB1 (z=−9.8), miR-34→MYC, miR-18/224→BCL2 … the Decision-F
  identifiability object a point-lasso can't give. Bayes' value = uncertainty + small-n regimes, not coupling.
  **Bar 5 — orthogonal protein layer (CPTAC), anti-circularity crux — IMPLEMENTED:** `eval/ood_protein.py`
  applies the TCGA-mRNA-fit M to CPTAC arms, correlating with protein + protein-vs-mRNA **discordance**.
  **TWO CPTAC cohorts, kept separate:** `prospective` = independent patients (non-circular OOD) vs
  `tcga105` = same TCGA patients (circular — layer-transfer only). Result: on prospective, **7/7 hub genes
  protein-coupled and 7/7 discordance-coupled** (translational residual), beats abundance 5/7; tcga105
  agrees in direction (robust, not a same-patient artifact). Discordance strongest on PDCD4/CDKN1B/PTEN.
  (tcga105 also lower-resolution: iTRAQ-4 MS vs prospective TMT-11 — doubly weaker.)
  **Bar 2 — learned vs curated fixed-M** (Design §0 ref 2): the MVP gate now reports the 3-way §0
  comparison (raw abundance / curated fixed-M = weights frozen at the prior / learned). Full-stack
  (family+ledger+deconv): **learned beats raw abundance 5/5, matches/beats curated fixed-M 4/5**, biggest
  win where the prior is weak (PTEN: abund +0.03, curated −0.02, learned −0.17, sparse 14/59); ESR1 the
  prior is already strong (−0.31→−0.32, "matches curated knowledge").
  **Transcription proxy (Design §B other half) — global latent factors OVER-CORRECT (negative result):**
  `--latent K` adds top-K genome-wide PCs to C; tested, it *eliminates real causally-confirmed edges*
  (miR-21→PTEN −0.17→+0.03, miR-17~92→BCL2 −0.29→−0.08, miR-206→ESR1 −0.32→−0.03) because the top PCs ARE
  the miRNA-regulated programs. ⇒ do NOT use global PCs; the transcription proxy must be orthogonal
  (intronic/pre-mRNA or gene-specific TF-regulon). Composition (deconv) is clean; transcription is not,
  via PCs. `_latent` kept as a documented negative control.
  **Phase-3b (partial-pooling identifiability) IMPLEMENTED:** `families.resolution_report(gene)` reads
  per-member resolvability from within-family abundance collinearity — family→gene is the identified unit;
  a member is nominable only where its abundance diverges. Result (hub panel): 4/8 selected multi-member
  families resolvable; the **miR-17~92 polycistron** on CDKN1A (corr 0.88) and **miR-200bc/429** on ZEB1
  (0.85) correctly held **family-only** (co-transcribed/co-varying), while miR-148/152→PTEN (0.09) and
  miR-18a/18b→ESR1 (0.40) resolve to a member. Honest "per-arm = nomination, not identified" (Design §F).

## What's left to do (not this reorg)
- **Optional `baselines/` packaging** — re-export shims for the 8 baseline modules (moving them would
  rewrite 10–44 importers each; lowest priority, cosmetic).
- **`git add`** the new tree if you want it version-controlled.
