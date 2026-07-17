# Learned model — VALIDATION DOSSIER (does it work, and where is it sharper)

The living evidence file for `mirna_hallmark/learned/`. Every claim here is a **real number, a recovered
positive control, or a specific sharpening/novelty example**, with the command to reproduce it. Companion to
`LEARNED_MODEL_ESTIMATOR_MAP.md` (what estimator does which job) — this doc answers *does each one actually
work, and what does it buy over the baseline*. Add to it as components are validated; keep the honest
negatives in §4 (they are part of the evidence).

Status legend: ✅ verified this arc · 📁 from a prior run (cited) · ⏳ to (re)generate under canonical M.

---

## §1 Component gates — the summary numbers

| Component | What it must show | Result | Reproduce |
|---|---|---|---|
| **Adaptive-lasso coupling** (the model) | learned X·M anti-correlates target OOF, FDR-controlled; beats raw abundance where weighting is applicable | ✅ **genome-wide 1119 POOLED-HE genes + BH/BY-FDR** (pooled-HE regen 2026-07-06): **597 (53%) FDR-significant coupling** (q_BY<0.05 & rho<0; q_BH 701). Learning **beats abundance 356/597 among FDR-sig genes**, 562/1119 (50%) overall — diluted by singletons where weighting is vacuous (§4). mean rho_model −0.123. *(per-regulator-count breakdown, HE-era 35%→88%, to re-tabulate.)* **⚠⚠ THE "BEATS ABUNDANCE" HALF OF THIS GATE IS RESTATED BY MH-127 — IT IS A *FITTING* RESULT, NOT A CURATION RESULT.** A **matched FAKE** regulator set (abundance-matched, site-free, non-seed-family, non-cluster, \|r\|<0.30 vs every real regulator), once **fitted**, reaches **78.9%** of the real model's 5-fold-OOF gene field (899 genes: real-fitted **−0.0465** vs fake-fitted **−0.0367**; **84.9%** under the structure-matched decoy, **91.3%** on multi-regulator genes), the paired real-vs-fake gap is **NOT significant** (Δ=−0.0098 [−0.0213,+0.0018], 52.5% of genes, **p=0.211**; FAKE2 p=0.577; MDE 0.0161 ⇒ **underpowered, not zero**), the **decoy gains from fitting as much as the real set does** (−0.0210 vs −0.0198), and **the fitted decoy BEATS real unfitted abundance** (Δ=−0.0100 **p=0.054** all genes; **−0.0357 p=3.6e-04** on multi-regulator genes). **⇒ an abundance baseline is the WRONG comparator: use a fitted matched decoy.** ✅ What DOES survive the decoy control is **out-of-cohort transfer** (next row). | `mvp.gate_fdr` → `output/learned/gate_fdr.tsv`; **decoy control:** `output/learned/mh127/` (`MH127_DECOY_MODEL_GENE_BUDGET.md`) |
| **⭐ DECOY-CONTROLLED TRANSFER** (MH-127 — the gate that replaces "beats abundance") | do the frozen TCGA β predict an INDEPENDENT cohort *better than a matched impossible design does*? | ✅ **YES, and it is the only curation-specific thing the model buys.** Frozen TCGA β applied to CPTAC (n=101, never refit; design z-scored **within** CPTAC because β are standardised): the **REAL** field predicts CPTAC **mRNA** (**ρ=−0.0189, p=8.9e-04**, 685 genes) and the **matched FAKE** field does **not** (−0.0056, **p=0.19**); paired REAL-vs-FAKE **p=0.030** (FAKE1) / **0.027** (FAKE2), and under FAKE2 it **survives composition adjustment** (Δ=−0.0126, **p=0.018**). Head-to-head partial is **asymmetric** (REAL\|FAKE −0.0198 p=3.6e-04; FAKE\|REAL −0.0041 p=0.29). **⚠ CPTAC PROTEIN CANNOT ARBITRATE:** under this from-scratch, correctly-pooled field the **REAL** model is itself **null** there (p=0.38 core / 0.75 comp) and the test is **underpowered** (MDE 0.018–0.020; power 0.56 for a ceiling-sized gap) — this is **not** "real = fake on protein". ⚠ **DISCREPANCY:** `eval/coupling_grid` (MH-125-fixed, re-run) still reports learned·protein·core −0.035 p=7.8e-12 over 1,162 genes; MH-127's 685-gene **trimmed** universe does not reproduce it ⇒ **diff the two field constructions.** | `output/learned/mh127/{p2_transfer_results.tsv.gz, mh127_cptac_rungs.tsv}` (§2b of the MH-127 doc) |
| **Prior magnitude** (does learning the magnitude matter, non-circularly) | functional evidence predicts repression more than binding-only, on external transfection data | 📁 functional-evidence ρ **+0.21** vs binding ρ **+0.08** (~3×), within-family, McGeary GSE140217 | `learned.evidence.calibrate` |
| **Prior robustness** (is the result hostage to prior hyperparameters) | gate stable hand-set vs data-calibrated class weights | 📁 PTEN −0.349 vs −0.348; gate 5/5 & 4/5; stability 0.85→0.87 — **non-critical** (prior used for *ordering*) | `w_prior_source="ledger"` vs `"ledger_mrna"` |
| **Spike-and-slab inclusion** (parked Bayesian alternative, now BUILT & benchmarked) | does an explicit evidence-graded inclusion prior `z~Bernoulli(π)` beat the lasso's soft-L1 selection on OOF coupling? | ✅ **it does NOT** — head-to-head on the SAME folds/C/family predictors. Hub panel: beats abundance 4/5, recovers curated drivers at PIP≈1.0, but matches/beats lasso only **2/5** (mean ρ −0.313 vs lasso −0.351). Genome-wide sample (107 genes): FDR-sig **35% vs lasso 45%**, mean ρ **−0.078 vs −0.125**, matches/beats lasso on **21%**. ⇒ inclusion prior buys a real per-edge PIP (interpretability) but not coupling; **confirms the design decision to ship the lasso and park the spike-and-slab** (rationale §2a) | `learned.spike_slab` (`run()`; `gate_fdr(limit=120)` → `output/learned/gate_fdr_spike_slab.tsv`) |
| **EB learned-shrinkage toward the prior** (the τ² question) | make the curated prior a *location* and let the data learn the anchor strength τ² — does adaptive anchoring beat the fixed-α lasso? | ✅ **coupling wash, but the learned τ² VALIDATES the design.** `β_m~N(c·ŵ_m,τ²)`, c/τ² learned by Gibbs. Hub panel mean ρ **−0.344 vs lasso −0.351**. Diagnostic `τ²/c²` (data dispersion ÷ prior-anchored part) = **71–57000 in 4/5 genes** ⇒ the data itself discounts the prior *magnitude* and shrinks little toward it (ESR1 the lone anchor-using gene, ratio 2.4, c=0.27). So "prior sets ordering not magnitude" is **derived from data**, not assumed. Side-finding: EB's ZEB1 win (−0.596 vs −0.438) is Gaussian-ridge>L1 on a cooperative family, **not** the prior (c≈0). | `learned.eb_shrink.run()` (rationale §2c) |
| **Learned-τ² NN-ridge vs lasso — GENOME-WIDE** (the coupling improvement that held) | does dense learned-τ² shrinkage beat the lasso's sparse selection on OOF coupling at scale, biology intact? | ✅ **YES, significant.** 655 multi-reg genes: `ebridge_nn` (L2, ν² learned, M≥0) mean ρ **−0.1683 vs lasso −0.1518, wins 58% (379/655), Wilcoxon p=9e-16.** `nnridge_cv`≈lasso ⇒ it's the LEARNED τ², not L2 (CV over-shrinks). Hyperprior-robust (InvGamma 0.01→20). CAVEATS: **tail-driven** (median Δ −0.003, mean −0.0165 — concentrated in many-regulator genes: CCL5, MET, E2F2, MMPs); **collinearity hypothesis FAILS** (Spearman(collin,Δ)=−0.074, non-monotonic) ⇒ mechanism is regulator-COUNT not cooperative-families. Dense (no selection) ⇒ better for the coupling job, not a lasso drop-in. | `learned.shrinkage_compare.full()` → `output/learned/shrinkage_compare.tsv` (rationale §2d) |
| **Student-t observation likelihood** (the dense mRNA likelihood's distributional form, MH-92) | are the dense residuals Gaussian, or heavy-tailed (→t) / count-like at low expr (→NB)? and does fixing it change any verdict? | ✅ **heavy-tailed → Student-t (NB rejected); a robustness upgrade, not a lever.** Diagnostic (85 genes × 1041 pts): pooled excess-kurt **3.5** (median 1.2/well-expressed), t-MLE **df≈7**, P(\|z\|>4)=33× Gaussian, from amplified-subset (+skew ERBB2/3) & near-floor (−skew, 3 genes); **NB rejected** (Spearman(expr,kurt)=+0.19, wrong sign). Built in `_gibbs_posterior(nu=)` (scale-mixture, Gibbs, `nu=None` bit-identical). VALIDATED: nests Gaussian (Δβ 0.001), **robust 1.76×** under 5% contamination, but **OOF coupling Δρ +0.001 n.s.** (Wilcoxon p=0.15), attribution corr **0.999** / 5-of-5 drivers ⇒ correctness/robustness, not performance, at n≈1000. CN-channel t (`t['nu']`) discounts collider-outlier segment (β 0.535→0.508); gate-SE bootstrap (`exclusion(n_boot=)`) confirms Sobel adequate (s²_π agrees 0.94–1.11×, δ²-floor-dominated). | `spike_slab._gibbs_posterior(nu=7)`; scratchpad `dense_resid_diag.py` / `t_oof.py`; METHODS §3b |
| **Coupling is PRIOR-ROBUST — no non-circular prior improves it (the unifying result, 2026-07-09)** | do evidence / seed-specificity / AGO-loading move OOF coupling, or only the inclusion LEVEL? | ✅ **inert to every prior; only p0 (→dense) moves it.** p0×spec_gain sweep: raising p0 climbs coupling toward dense (FDR-sig 45%→**64%**); spec_gain **FLAT at every p0** (FDR-sig 48/48·55/55). Seq-spec in the dense-model slab: ΔρOOF **+0.000** every gene/gain (likelihood dominates the slab at n≫p). AGO-loading·abundance: Δρ **−0.001**. Shuffled-evidence null: dense coupling prior-independent. ⇒ **coupling is carried by abundant arms + absorbed by the learned M; the priors' value is discovery/identity, not coupling.** | `genome_p0_spec_sweep`, `genome_seq_fdr`, `seq_spec_bench` (rationale §2b/2h) |
| **Learned-τ² for ATTRIBUTION** (EB posterior vs bagged NNLS) | can the same τ² model read the per-edge weight AND give calibrated uncertainty, matching bagged NNLS? | ✅ hub: reproducibility EB **0.996 ≥ NNLS 0.983** (single-lasso 0.03); agreement ρ 0.844; identified set |mean/sd|>2 recovers PTEN miR-21/103/181/182 (4/4), GATA3 miR-27, ESR1 miR-18, and **declines** §9-unidentified splits (ESR1 221/222, ZEB1 miR-200) — EB |z|>2 **reproduces §9's identifiability verdict**. CAVEAT: raw EB mean biased up for un-informed edges (half-normal prior mean>0) → read via |z|>2 or evidence-π. ⇒ deployable as NNLS point-estimate + EB posterior sd; one model could fuse coupling+attribution+§9. | `learned.attribution_eb.run()` (rationale §2e) |
| **Bagged-NNLS attribution** (stable coefficients) | per-edge weight reproducible under resampling (single lasso is not) | ✅ single-fit bootstrap corr **0.03 → bagged 0.99** (PTEN); ESR1 0.93, ZEB1 0.99 | `states.canonical_M(g,"01")` two seeds |
| **Family collapse + conditioned nomination** | classify family as single-driver / co-drivers / family-only by conditioned anti-corr | ✅ §3.2: 13:7:11 (PTEN miR-17~92 = miR-106b+miR-20a **co-drivers**; ESR1 miR-221/222 family-only) | `families.resolution_report(g)` |
| **Identity vs magnitude** (Decision I) | who EXPLAINS the target (Shapley R²) vs who exerts force (budget M·X̄) | ✅ diverge by abundance-vs-coupling: ESR1 miR-18 explains 0.54 on 0.20 budget (quiet on-target owner, indiv ρ −0.43); miR-22 budget 0.26 / identity 0.01 (loud passenger, indiv ρ −0.12) | `attribution.identity_vs_magnitude(g)` |
| **Sequence-designed identity model** (affinity-percentile, §16b — the SEQUENCE "who owns g", vs the row above's REALIZED Shapley identity) | does the abundance-removed, expression-free ownership model recover real ownership non-circularly? | 🟡 **face-valid + breast-binding-enriched, but seed-caveated & NO model value.** Recovers canonical specialists (miR-19/26→PTEN, miR-15/16→CCND1, miR-200/141→ZEB1) abundance-removed; **breast** MCF7 PAR-CLIP **p=7.5e-7** / MCF7+MDA HITS-CLIP p=5.9e-5 — non-circular through HE (CLIP∉HE, <2% overlap) but **seed-circular** (Manakov chimeric p=7e-8 did NOT independently replicate — Manakov-exclusive p=0.53; no breast chimeric exists); **arm-level L2 91% cross-source** (the robust part). MODEL-VALUE NULL: identifiability 0% top-driver flip genome-wide. | `attribution_identity`, `identity_payoff_genome` (SYNTHESIS §6c) |
| **CN-locus instrument** (causality) + cluster-exclusion | miRNA-locus CN → abundance → target (reduced<0), AND the effect is attributable to the specific arm, not its genomic cluster | ✅ 60 edges: **59/60 usable** (F>10), **46 causal-concordant**, but after the cluster check only **24 STRONG** (**17 cluster-clean** + **7 pleiotropic-but-arm-unique**). Clean-strong: miR-30d→GNAI2 (reduced −0.41), miR-200a→CTNNB1, let-7b→IFNB1; rescued: miR-96→ADCY6/FOXO1 (co-targeter miR-182 flips +), miR-141→PTEN (−0.28) | `instrument.run_clean(60)` → `output/learned/instrument/cn_instrument_clean.tsv` |
| **Cross-gene pooling** (honest negative) | pooling should help weak/small-n genes IF genes share regulators | ✅ **does NOT help** on p53 set (Δρ −0.002 full, −0.045 at n=150); neutral on miR-29 module (+0.001) — structure-bound, see §4 | `hierarchical.oof()`; `program_network.module_pooling_test` |
| **Subtype-specificity** | LEVEL vs common-M SLOPE vs WITHIN-subtype wiring (independent vs pooled) | ✅ pressure **LEVEL** differs 10/10 (Kruskal q≈0); common-M **SLOPE** differs 2/10 (RB1,CDKN1A). Independent two-fit `subtype_wiring` *suggested* PTEN wiring-remodels — but the **pooled Bayesian-nesting estimator** (`subtype_wiring_pooled`, H: ridge M_s=M_all+δ, regularised Δ vs n-matched null) shows that was a **noise-difference over-call: PTEN is borderline ~1.35× across λ=1–30 (not robust)**; only **RB1-Her2 modest (~1.5–1.66×)**. ⇒ **subtype wiring mostly CONSERVED**; robust subtype signals = LEVEL + common-M slope, not rewiring | `learned.subtype` (`subtype_coupling` / `subtype_wiring` / `subtype_wiring_pooled`) |

---

## §2 Positive controls — recovered known biology

The model recovers canonical, independently-established miRNA→target biology **without being told the answer**
(the prior is inclusion+evidence, not the in-cohort magnitudes):

- **miR-21 → PTEN** — canonical oncomiR/tumour-suppressor axis. Recovered as PTEN's **rank-1 budget driver**
  (share 0.31) and a top learned weight. ✅ `states.budget_shift("PTEN")`
- **miR-200 / miR-429 → ZEB1** — the canonical EMT double-negative switch. Recovered as ZEB1's drivers, and
  the wiring is **conserved healthy→tumour** (transfer retention 0.93 — the switch is intact, only abundance
  shifts). 📁 `states.stable_wiring("ZEB1")`, `state.cross_state_transfer`
- **miR-29 → collagen/ECM** — established ECM-repressor. Recovered as the **broadest EMT hub** (21 program
  targets, ADAM12/COL/LOX/MMP2). ✅ `program_network EMT`
- **miR-18 → ESR1**, **miR-34/449 → MET** — recovered as top drivers of their curated targets. ✅ §3

---

## §3 Sharpening & novelty — specific examples

The vignettes that show the method is *sharper* than the abundance/evidence baselines, and finds new biology.

### 3.1 Learning re-weights off literature-prominence onto in-cohort signal (MET)
The **frozen evidence weights give the wrong sign** for MET (rho_curated **+0.23**, i.e. "more miRNA → more
mRNA") because they are dominated by the *most-studied* edges. Learning:
- **keeps** miR-34/449 (evidence rank 1 — the one strong edge that also tracks MET),
- **data-selects** miR-128 (rank 12) and **miR-181a-3p (rank 21 — near the bottom)**,
- **drops** the well-studied-but-flat miR-1/206, miR-23, miR-199, miR-31, miR-144 (evidence ranks 2–6),

recovering correct repression (rho_model **−0.17**, gain +0.40). Same pattern on IGF1/ERBB2/E2F1/MYC. ✅
This is the core value: **literature prominence ≠ in-cohort predictive strength**, and the model corrects it.

### 3.2 Conditioned nomination resolves the driver(s) inside a collinear family — by distribution
Conditioning each family member on its family-mates classifies the family (§ METHODS 9) into
**single-driver / co-drivers / family-only** by how many members' anti-corr survives (part < −0.1).
Across the panel: **13 single-driver : 7 co-drivers : 11 family-only**. Key cases:
- **miR-17~92/106b → PTEN is CO-DRIVERS: miR-106b (−0.22) + miR-20a (−0.13)** — *not* miR-17 (flips to +0.09
  once conditioned = it was riding the aggregate). miR-106b is the *strongest* survivor but not the sole one.
  (An earlier "miR-106b is *the* driver" over-simplified; the co-driver-aware nomination is the correct read.)
- **ZEB1 = miR-429 + miR-200c** together (the canonical miR-200 family, co-drivers), miR-200b flips positive.
- **family-only** (unidentified): ESR1 miR-221/222 (both ~−0.05, collinear) — the model honestly says so.
The old abundance-divergence heuristic saw none of this (it never looked at the target).
✅ `families.resolution_report` → `output/learned/attribution/resolution_panel.tsv` (+ `budget_shift_panel.tsv`, 382 rows/20 genes)

### 3.3 The program network surfaces the regulatory architecture (all 5 Hallmarks)
Rolled up + persisted for all five priority programs (`output/learned/program_network/`). Each is a sparse
network (density ~0.02) **concentrated on a few broad hub families holding 73–82% of program weight**:

| Program | genes × families | edges | hubs (≥2) | hub weight | broadest hub (n_targets) |
|---|---|---|---|---|---|
| EMT | 110 × 153 | 342 | 60 | 76% | miR-29 (21, ADAM12) |
| P53 | 88 × 152 | 328 | 70 | 82% | miR-27 (10, BTG2); miR-25/32/92 (9, TP53) |
| G2M | 95 × 167 | 310 | 66 | 73% | let-7/98 (8, AURKB) |
| IFN | 67 × 141 | 247 | 47 | 73% | miR-17~92 (12, ZNFX1) |
| HYPOXIA | 92 × 160 | 328 | 68 | 81% | miR-155 (9, MXI1) |

The `focus` metric (max/sum) separates broad regulators (miR-29 = 0.13) from narrow ones (miR-21 = 0.68).
**Positive control — the recovered hubs are the canonical program regulators**: miR-17~92/106 (cell-cycle
/IFN), let-7/98 (G2M — AURKB), miR-155 (hypoxia), miR-29 & miR-200 (EMT), miR-15/16 (P53/G2M). And a
cross-program readout the per-gene view can't give: the **miR-17~92/106 cluster and the miR-15/16/195/424/497
cluster are top hubs in 3 programs each** — genuinely pleiotropic program regulators. ✅ `program_network EMT|P53|G2M|IFN|HYPOXIA`

### 3.4 Discovery — validated novel edges (full genome-wide run, POOLED-HE 2026-07-06)
✅ Full pass over **1571 pooled-HE genes**: **2955 candidate orphan edges vs 2 permutation-null hits → edge-level
FDR ≈ 0.001**; after composition control, **1579 deconv-robust across 460 genes, 289 fully novel**. *(Just below
the HE-era 3034/1589/291 because the ~640 TarBase functional edges are now KNOWN, not orphan discoveries.)*
- Top robust hits cross-validate the day's other findings: **miR-135b→GATA3** (partial_deconv −0.42, novel,
  retention 0.84) and the **miR-17~92/106b cluster → RABEP1, TGFBRAP1, AHNAK, IL6ST, AFF1, M6PR, WWP1, SYNE1**
  — nominated by **miR-106b/93/17**, the same driver the 4b conditioning picked inside the collinear family.
- **miR-135b→GATA3** also validated on CPTAC **protein −0.49** + literature ceRNA axis (spiker, 28%>floor
  ⇒ coupling in the miR-135b-high subset). 📁 111 orphans: 60% carry chimeric-eCLIP duplexes (Manakov).
- Output: `output/learned/discoveries.tsv` (1579 robust), `discovery_candidates.tsv` (2955). `discovery.run_all`
- **Program enrichment** (`learned.enrichment`, hypergeometric + BH): the 462 robust-target genes are
  FDR-enriched in **TGF_BETA_SIGNALING** (fold 1.7, q=0.043) and UV_RESPONSE_DN (q=0.028). Hub families
  (miR-106b/17) show **no single-program concentration** (nominal TNFα-NFκB p=0.015, q n.s.) — coherent with
  their cross-program pleiotropy (§3.3): broad regulators don't enrich one program.
- **TGF-β composition-robustness**: the TGF-β edges survive the finer **myCAF + within-malignant-EMT**
  controls (`assemble_gene(mycaf=True, mal_emt=True)`) — median coupling retention **0.90** (XIAP/ACVR1/
  SMURF1-2/ARID4B ≥0.9; only the EMT-linked KLF10 edges drop to ~0.65). Not a lumped-CAF or EMT-state artifact.
- ⏳ **FDR is currently edge-level (flat)** — should be recomputed at the **seed-family** unit (Simes-per-family
  BH + per-family min-statistic null, mirroring `coupling_inference.py`) so family-mate orphans aren't
  double-counted. See §5.

---

## §4 Honest negatives & caveats — also evidence

These bound the claims and show the gates are real (things fail them):

- **⭐⭐ THE DECOY MODEL (MH-127) — the deepest negative in this dossier: IN-COHORT, THE LEARNING IS NOT MORE
  THAN FITTING.** Learn β on a **matched fake** regulator set (abundance-matched, **site-free** — no
  TargetScan site *and* no scanMiR duplex — non-seed-family, non-genomic-cluster, |r|<0.30 vs every real
  regulator) and the fitted decoy reaches **79–91%** of the real model's OOF gene field; the paired
  real-vs-fake test is **n.s. (p=0.211 / 0.577)** and **no stratum rescues it** — not regulator count (1-3
  p=0.599 · 4-8 p=0.173 · 9-15 p=0.508; **26+ is not decoy-controllable at all**) and not Shapley class
  (SINGLE_DRIVER p=0.273 · CO_DRIVERS p=0.176 · DIFFUSE p=0.248), with the gap running **opposite** to the
  prediction (smallest where credit is concentrated) and the class ordering **indistinguishable from shuffled
  labels (perm p=0.92)**. Structural cause of the aggregate: **64% of trimmed genes have ONE family**, where
  `Xz@β` is a monotone rescaling of the abundance pool ⇒ **R3 ≡ R1 exactly** (99.5% exact ties) — no fitting
  happens there at all. **UNDERPOWERED ≠ NEGATIVE:** MDE 0.0161, and a leakage lens measured the per-gene SE
  is understated ~1.8× ⇒ *a real advantage >0.021 ρ is excluded; a smaller one is not*, and the point estimate
  favours real in 5/5 seeds and both decoys. ✅ **What survives:** out-of-cohort **transfer** (§1) and a
  **fit-free** stratifier (top tertile of curated evidence weight `w`: real > fake **fitted** q=0.0059 *and*
  **unfitted** q=0.0059 — but **not** under the structure-matched decoy, q=0.20 ⇒ a lead, not a result).
  ⇒ **RULE: an abundance baseline is not a control. Benchmark any aggregator against a FITTED matched decoy.**
  (`MH127_DECOY_MODEL_GENE_BUDGET.md`; `output/learned/mh127/`)
- **Cross-gene pooling doesn't earn its place** (§1). Fails on the p53/cell-cycle set (hurts at small-n), only
  neutral/marginal on the miR-29 ECM module — because the HE network is sparse (density 0.02) and effects are
  gene-specific. ⇒ per-gene is primary; `program_network` roll-up is the program view (aggregate *after*
  fitting, no pooling bias).
- **Occupancy saturating link failed the held-out gate** → linear a·M kept (`occupancy.py` parked); cooperativity/ceRNA gated behind it.
- **Prior hyperparameters are non-critical** — a feature, not a null: the result does not hinge on the class
  weights, because the model uses the prior for ordering/selection and estimates magnitudes from data.
- **Cross-state ΔM** compares an n≈1000 tumour fit to an n=104 NAT fit — the small-n fit is more diffuse, so
  the wiring sum partly reflects n, not only biology (read the top movers, not the tail).
- **Composition confounding** — GTEx-breast (stroma-rich) vs tumour (epithelial) mismatch is a headwind
  against cross-state transfer; "conserved" (ZEB1) is the robust call, "remodeled" (PTEN/ESR1) is
  composition-suspect without deconvolved M_H.
- **"Beats abundance" is only meaningful for multi-regulator genes** — 396/1033 HE genes are singletons
  (one regulator family), where X·M has the same ranking as abundance up to a per-fold scalar (mean |Δ|=0.03,
  47% exact ties); the OOF refit then only adds estimation variance (nothing to weight) → singletons "beat"
  just 19% (10% among FDR-sig). Not a model failure — the right tool for a singleton is the direct single-edge
  coupling (+FDR), not the aggregate. The model's value is specifically multi-regulator genes: the equal-weight
  abundance baseline degrades as more (noisy) regulators are averaged in, so the weighting advantage grows with
  count (beat-rate → ~95% at 8+). **Checked whether the singleton fit truly harms (2026-07-05): it does not**
  — mean harm rho_model−rho_abund = +0.026, and **0/104 real-coupling singletons lose their edge**; the harm
  is confined to the *gate metric*, while the attribution/card layer already reports singletons via direct
  partial-Spearman. So the penalty-bypass is not worth implementing (criterion: only if it really harms).
- **The CN instrument rejects some observationally-strong edges** — 14/60 are non-concordant; notably
  miR-17/20→RUNX1 and →CDKN1A have strong *observational* coupling (−0.30 to −0.40) but ~0 reduced-form
  (F>60, so not a weak-instrument artifact) — the causal test does not support them, a caution that
  observational anti-correlation ≠ causal repression for a subset of edges.
- **CN instrument is a CLUSTER perturbation** — of the 46 concordant edges, 22 are pleiotropic AND not
  uniquely attributable to the arm (a co-located co-targeter co-carries the effect), so the causal claim is
  cluster-level, not arm-specific. E.g. **miR-17→PTEN** is concordant but miR-18a co-carries it (−0.24), and
  the **miR-200/429→ZEB1** edges are irreducibly cluster-level (miR-429 always survives conditioning). The
  honest arm-specific causal set is the **24 "strong"** edges (§1), not all 46. The within-cluster attribution
  (does the arm uniquely keep its anti-corr conditioning on cluster-mates) is bounded-circular by design — CN
  supplies exogeneity, expression only resolves *which* cluster member carries it (`instrument.cluster_attribution`).

---

## §5 Multiple-testing / FDR status (hypothesis families, precursor-style)
**Principle** (mirrors the precursor `coupling_inference.py`): the **seed family is the hypothesis unit**
(Decision F — family→gene is identified, arm is a nomination), so family-mate arms are ONE test, not N.
Controls used: **BH** (independence), **Benjamini-Yekutieli** (arbitrary dependence — what shared regulators
induce), **Simes-per-family** (combine within-family p → BH across families), and a **per-family
min-statistic permutation null** (family-wise error).

| Test | hypothesis unit | FDR status |
|---|---|---|
| Gate (aggregate coupling) | gene | ⏳ genome-wide `mvp.gate_fdr` running — **BH + BY** across genes (BY for the shared-regulator dependence); family-collapse already removes arm redundancy *within* a gene |
| Discovery (orphan edges) | (family → gene) | ✅ edge-level permutation FDR ≈0.001; family collapse counted (3034→1878, robust→1100, novel→251). ⏳ a *proper* family-level FDR wants the **per-family min-statistic scan null** added to `discovery.scan_all` |
| Coupling (per-edge) | (family → gene) | ✅ `learned.coupling` — **family-aware permutation** partial-Spearman (Freedman-Lane, per-gene C incl. target_cn, seed-family dependence preserved) → p_perm / p_z / q_bh / q_by / **q_family** (Simes). `permutation_pvalues` confirmed generic (no precursor coupling); called per-gene for the per-gene target_cn. Validated (PTEN top = miR-106b −0.50). ✅ genome-wide (pooled-HE 2026-07-06): **5499 edges / 1507 genes → 1746 edge-sig (q_BY_global<0.05 & rho<0), 2222 family-sig → 1895 distinct significant seed-family→gene** (`output/learned/edge_coupling.tsv`) |
| Instrument (causal) | edge / genomic cluster | ✅ reduced-form **p + BY** added (`run_clean`): 24 strong-causal → **18 also FDR-sig** (q_BY<0.05); cluster-attribution in §1 |
| Identifiability (hierarchical) | edge | \|β/sd\|>2 per edge, ⏳ no cross-edge FDR |

Discovery is the one component with a real (permutation) null today; the rest need the p→BH/BY/Simes layer.
`gate_fdr` is the first, reusing `mirna_hallmark.stats.bh_fdr` + `coupling_inference.{benjamini_yekutieli,family_simes_fdr}`.

## §6 The CPTAC / PROTEIN layer — what validated, after the channel died

> ⛔ **THE PROTEIN CHANNEL'S CENTREPIECE IS FALSIFIED. `βᵗ` — a per-family translational-repression latent — is
> NOT SUPPORTED at n=101 (MH-103/MH-104).** The prior evidence for it was a **mediator LEAK**: every earlier βᵗ
> result built `protein_resid = P − a·M` with **`a` fit on ALL samples**, then OOF-ed `X → protein_resid` — the
> mediator was fit using the test folds. **Fit `a` inside the fold and the signal disappears.** Four independent
> leak-free framings agree, including the maximally-powered **AGGREGATE** one (TCGA-fixed weights, **1 df**, zero
> parameters fitted in CPTAC): **q<0.10 in 1/17 genes — only BCL2**; **PTEN, the doc's headline throughout,
> does NOT survive (z=0.53; d=−0.006, p=0.82).** ⇒ **NOT a leak, NOT an estimator choice, NOT a p>n artifact —
> the TRUE EFFECT SIZE meeting the cohort size.** Our own biology already said so: MH-34 found a translational
> residual in **11/1132 genes ≈ 1%**, which on a 17-gene panel predicts ~0.2 hits. **There was never a βᵗ field
> to fit.**
>
> ⚠ **What this does NOT overturn:** MH-34/35's marginal **association** results. The translational component is
> **real, weak, and widespread-in-direction** (13/17 repression-directed, binomial p≈0.025) but **not resolvable
> per gene at n=101**. **The falsification is of the MODELING OBJECT, not of the biology.**
>
> **The five results below are what SURVIVES that falsification — none of them depends on `βᵗ`.** Do not read
> any of them as reviving a protein channel; §6.5 is the measured reason one can never exist.

### 6.1 `locus_cn_cptac` — BUILT + VALIDATED, but fidelity ≠ power
An independent **second** CN-instrument cohort. **Median r = 0.997 vs ASCAT truth**; the coordinate hole is
patched (**482/482 loci**, via `mirna_mature_loci.csv` → `pre_gene_id`); matrix **CPTAC 122×414**, median CN 1.99.
(The gene-level→locus proxy itself hits **r=0.998** vs ASCAT with **no gene deserts** — median 28 kb, zero >1 Mb.)

⚠ **But high fidelity buys little power — the binding constraint is n, not measurement.** First stage (locus CN →
arm dose), F:

| cohort | median F | arms F>10 | families resolvable (≥2 distinct loci at F>10) |
|---|---|---|---|
| **TCGA** (n=961) | 3.16 | **234 (32%)** | **31** |
| **CPTAC** (n=101) | 1.10 | **59/685 (8.6%)** | **5** |
| both | — | 49 arms | **5** — let-7 · miR-15/16 · miR-17~92 · miR-181 · miR-25/92 |

First-stage γ replicates across cohorts at Spearman **0.349**. ⇒ the full CN-fused δ comparison is feasible but
**low-powered**; the abundance-based test (§6.2) is the well-powered one and is what we report. ✅ MH-104

### 6.2 V1 — δ-TRANSPORTABILITY VERIFIED (the cross-cohort assumption the program rests on)
Member **dose-share** (which member carries the family's dose = δ's primary input), TCGA n=1041 vs CPTAC n=101, on
the **58 multi-member families expressed in both** (canonical floor RPM≥10):
- **Same dominant member: 84.5%** vs **43.7% chance** ⇒ **1.93× enrichment**
- **median member-share correlation 0.991**; median total-variation distance **0.145**
- (all 260 families incl. noise-dominated: 73.5% / TV 0.230 — the degradation is entirely in undetectable arms)

⇒ **δ's abundance component transports.** ⚠ **SCOPE — read this before citing it:** this validates the
**abundance input** to `delta_pooling` **only**, *not* its full mRNA-width + CN + chimeric fusion.
✅ `output/learned/v1_dose_delivery_transport.tsv` (MH-104)

### 6.3 ⭐ BAR-5 PASSES (581 genes) — and it locates the programme's real limit
`β` fit on **TCGA mRNA only** (never sees CPTAC, never sees protein), scored on 581 genes. The essential control
is rung (i) — fit on TCGA half-A, score on TCGA half-B — the same-cohort out-of-sample **CEILING**, which is what
distinguishes "β doesn't transport" from "my aggregate is broken":

| rung | median ρ | % repression-directed | binomial p | **retention vs ceiling** |
|---|---|---|---|---|
| **(i) TCGA held-out mRNA** — same cohort, same layer (the CEILING) | **−0.117** | **88.8%** | **2e-88** | **1.000** |
| (ii-a) CPTAC mRNA — STAR | −0.017 | 54.6% | 1.5e-2 | 0.148 |
| **(ii-b) CPTAC mRNA — LinkedOmics** (the default) | **−0.023** | 56.1% | **1.8e-3** | **0.193** |
| **(iii) CPTAC protein** — cross-cohort **+** cross-layer | **−0.023** | **57.7%** | **1.3e-4** | **0.201** |

1. **The aggregate is NOT broken** (88.8% negative in TCGA held-out, p=2e-88) ⇒ the flat cross-cohort result is REAL.
2. **BAR-5 PASSES**: the learned M carries to an **independent cohort AND an independent layer** (p=1.3e-4) — the
   honest, FDR-controlled version of MH-83's hand-picked 7-gene "7/7".
3. ⭐ **THE COHORT JUMP COSTS ~80%; THE mRNA→PROTEIN JUMP COSTS ~0%.** Crossing cohorts: 0.117 → 0.023. Crossing
   mRNA→protein: 0.193 → 0.201 (protein ≈ mRNA). **This independently re-derives the βᵗ falsification from a
   different direction** (protein adds nothing beyond the transcript) **and locates the real barrier of the whole
   cross-cohort programme.** Corroborates the state session's 0.6% (MH-102d). ✅ MH-104

⚠ **The ~80% is a MEDIAN hiding a MIXTURE (MH-106).** Observed sd(ρ_protein) **0.136** vs **0.100** under a pure
null at n=101 ⇒ **1.36× excess dispersion**; corr(ceiling, ρ_CPTAC_protein) = **+0.176 (p=1.8e-5)** ⇒ a real
SUBSET transfers. **Corrected claim: β transports for the genes where β is WELL-DETERMINED** (many regulators ·
strong in-cohort coupling) — retention **0.22–0.39** — and not at all for the poorly-determined majority.
⛔ **RETRACTED (MH-109):** the `a_g`-stratified reading of that mixture ("the miRNA effect reaches protein only
where protein tracks mRNA") **does not hold** — an artifact of coarse binning + an unstable ratio-of-medians
retention statistic. The continuous test says nothing (a_OLS +0.010 p=0.88 · a_IV +0.032 p=0.65), and instrument
strength F — a pure nuisance — predicts equally well. **UNSUPPORTED, not disproven:** you cannot detect a
modulation of a signal that is itself at the noise floor.

### 6.4 `a_g` — VALIDATED as the mRNA→protein propagation slope (MH-105)
One coefficient per gene ⇒ well-powered (median **0.432**; median 0.397 over 9,196 genes; IQR 0.23–0.62).
⭐ **Falsification test PASSED — complex subunits are BUFFERED:** obligate-complex subunits (RPL/RPS/PSM/SNRP,
**n=139**) have median `a_g` **0.117** vs **0.437** for all other genes (**n=9,057**), **Mann-Whitney p=1.4e-29**
— protein-complex buffering, recovered genome-wide, exactly as predicted.

⚠ **REFINED (MH-106, official CORUM 5.3): buffering is a LARGE-OBLIGATE-COMPLEX phenomenon, NOT generic complex
membership.** CORUM subunits (n=3,519) median `a_g` **0.410** vs non-CORUM (n=5,677) **0.443** — real but **tiny**
(p=4.5e-4). The evidence is the **dose-response with complex size**: largest-complex **2–3 members → 0.407** vs
**≥31 members → 0.261**; **Spearman(complex size, a_g) = −0.079, p=2.6e-6.** The RPL/RPS/PSM/SNRP proxy above
over-states the *general* claim precisely because it **IS** the large-complex set.

⛔ **RETRACTED: "the unbiased `a` is up to 38% off the marginal one."** That came from a **9-gene hand-picked
panel**. Genome-wide the miRNA-confound correction is **IMMATERIAL** (paired median bias **−0.001**;
adjusted<marginal in **52.4%** = a coin-flip, 597 genes). ⚠ The *mechanism* is real — the miRNA is a genuine
confounder of the mRNA→protein slope (a common cause of M and P) and it does bite on strongly-coupled genes
(ZEB1 0.551→0.340 = −38%, BCL2 −25%, VIM −19%) — **but those are unrepresentative, not a global bias.**

⚠ **NAME COLLISION — two different `a`s, do not conflate:** `a_g` (here) is the **per-gene mRNA→protein
propagation slope**, and its per-gene variation *is* the biology. `gauge.Gauge.a` is **ONE GLOBAL** cross-cohort
scale **NUISANCE** (composition/platform/C), explicitly never claimed. See METHODS §18a.

### 6.5 ⭐ THE FISHER-INFORMATION BOUND — why protein can NEVER be a coupling lever
Pre-registered, and it is the structural reason §6.3's layer jump is free while the cohort jump is not:

| quantity | measured (under `build_C("cptac")`) |
|---|---|
| `a_g` (protein~mRNA slope) | **0.397** (median, 9,196 genes) |
| σ²_mRNA / σ²_protein | **0.79** (0.643 / 0.819) |
| **information ratio protein→β** | **≈4–6%** — **≤7.6% even at a_g = 1.0** (the pre-registered ceiling) |

⇒ **protein carries only ~4–6% of the mRNA channel's information about β, and ≤7.6% at ANY `a_g`. It cannot move
β. Robust to every plausible `a_g`.** ⭐ **Converges with the state channel's independently-derived 0.6%** for
CPTAC-mRNA→TCGA (MH-102d) — two sessions, two routes, same structural conclusion: **the value of an exogenous
source is a NEW LATENT, never a coupling gain on the same β.** (And the new latent βᵗ is itself falsified — see
the banner. Both doors are closed.)

⛔ **Do NOT cite the "1.7% / ≤8.8%" figures — they are VOID.** They were computed with **PAM50 in the confounder
block**, which `lineage_verdict` **prohibits**: PAM50 is computed *from* the mRNA, **27/50** of its classifier
genes are our targets, and it costs **−36% of |β|**. The values above were **re-measured, not rescaled**
(measured-only gate). ⭐ **The verdict is unchanged either way — at 1.7% or at 4–6%, protein cannot move β** — and
the correction makes the ceiling **tighter** (≤7.6% vs the void ≤8.8%). *(⚠ The **1.2%** once quoted for the ratio
used the ATTENUATED observational `a_g`=0.397; the CAUSAL CN-instrumented `a_IV`=0.893 gives ≈6%, and the direct
Bar-5 retention² check gives ≈4% — MH-108.)*

**The honest headline this axis delivers:** *the protein layer CONFIRMS the mRNA-mediated model and BOUNDS the
translational residual below the per-gene resolution of a 101-patient cohort.* A real, publishable, negative result.

## Reproduce index
```
mvp.oof_gate(gene, w_prior_source="ledger", family=True)     # coupling gate (beats abundance / curated)
mvp.gate_fdr()                                                # genome-wide gate + BH/BY FDR (54% sig; 77% beat among multi-reg)
learned.coupling.run() / coupling.edge_coupling(gene)         # per-edge family-aware permutation FDR (p_family)
learned.enrichment                                            # discovery/hub program enrichment (hypergeom+BH)
learned.subtype                                               # subtype-specific coupling (Kruskal level + interaction slope)
instrument.run_clean(60)                                      # CN causal + cluster-attribution (24 strong / 18 FDR-sig)
states.canonical_M(gene, "01")                                # canonical bagged attribution weight
attribution.identity_vs_magnitude(gene)                       # identity (Shapley coupling ownership) vs magnitude (budget)
states.budget_shift(gene)                                     # within-gene budget, family-apportioned
card.gene_card(gene) / card.decompose(gene)                   # per-edge attribution card / ΔM decomposition
families.resolution_report(gene)                              # conditioned driver nomination
program_network EMT|P53|G2M|IFN|HYPOXIA                       # program miRNA×gene network + summaries
program_network.module_pooling_test(anchor, genes)           # (b) pooling gate on a co-regulated module
learned.evidence.calibrate                                    # transfection calibration of class weights
```
