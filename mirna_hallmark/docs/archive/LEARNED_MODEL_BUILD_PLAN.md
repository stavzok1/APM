# Learned model — build plan: reuse-vs-rebuild & repo architecture

Companion to `LEARNED_MODEL_DESIGN_RESPONSE.md` (the design) and its §6 phased plan. This answers the
two "how do we start" questions: **(1) what we rebuild from scratch vs reuse**, and **(2) how to
organize the folder**. Written 2026-07-04.

> **STATUS (2026-07-04): the sprawl reorg is DONE (waves 1 + 2).** Top-level **160 → 36 `.py`**. Wave 1
> moved 74 in-degree-0 leaves (zero rewrites); wave 2 moved 50 theme-base modules (234 import-site rewrites
> across 106 files). Verified both times (`py_compile` all OK, no dangling refs, spine/`run_all`/moved
> bases+importers import OK). Top level now = SPINE + BASELINE + `run_all`'s pipeline steps only. See
> `MODULE_MAP.md`. **`learned/` package scaffolded + Phase-1 MVP passing** (2026-07-04): gene-focused NN
> adaptive-lasso (`learned/{data,regression,mvp}.py`, spine reused, evidence-weighted penalty, C
> residualised) beats raw-abundance coupling OOF on **4/5** hub genes and recovers canonical drivers
> (miR-200→ZEB1, miR-206/miR-18a→ESR1, miR-106b/miR-29c→PTEN). **Phase-2 evidence ledger IMPLEMENTED**
> (`learned/evidence/`): PMID-deduped method-centric fusion of miRTarBase + TarBase v9 — 1.92M→1.11M rows
> (**42% double-counting removed**), canonical top edges, wired into the MVP (`--prior ledger`; holds gate,
> surfaces miR-21→PTEN). **Phase-3a seed-family pooling IMPLEMENTED** (`families.py`, `--family`): gate
> **4/5→5/5** (CDKN1A passes), stability up where collinear (ZEB1 0.70→0.81), canonical family→gene
> estimand (miR-200 subfamilies on ZEB1; miR-15/16 + miR-17~92 on CDKN1A). **Composition confounder
> (Design §B) WIRED** (`--deconv`, CIBERSORTx non-malignant fractions): stroma-expressed miR-29→ECM edges
> attenuate to ~0 (COL3A1 −0.35→−0.02), cell-intrinsic ESR1 survives (−0.38→−0.32) — the confound
> separated. **CN-locus instrument (Bar 4, causal) IMPLEMENTED** (`instrument.py`): 12/12 valid (F 24–174),
> 7/12 causal-concordant (miR-17~92→BCL2/MEF2D, miR-320a→TFRC, let-7b→CDC34); flags miR-17→CDKN1A/RUNX1 as
> possibly confounded. **Bar-5 protein OOD (CPTAC) IMPLEMENTED** (`eval/ood_protein.py`): on the INDEPENDENT
> `prospective` cohort (NOT the same-patient `tcga105` — kept separate), TCGA-fit M is 7/7 protein-coupled +
> 7/7 discordance-coupled.
>
> **STATUS UPDATE (2026-07-04, late — reconciled against §5 Bars + §6 Phases):**
> **DONE:** Phase 0 (obs-model audit — TF-regulon strength-weighted transcription proxy [global PCs
> over-correct = negative control] + CIBERSORTx composition, both in `C`); Phase 1 (adaptive-lasso MVP, 5/5);
> Phase 2 (PMID-deduped method-centric ledger + **transfection-calibrated** method weights — functional≫
> binding ~3× non-circular, gate robust; `evidence/calibrate.py`, `CLASS_WEIGHT_MRNA`); Phase 3a/3b
> (seed-family pooling + identifiability `resolution_report`); Phase 5 (CN instrument + cluster/pleiotropy;
> **SV instrument tested → NOT VIABLE**, CN covers standalone arms). Bars 1,2,4 pass; Bar 5 (protein) passes.
> **scanMiR biochemical K_D** built (`kd.py`, 912k edges): as occ-κ = NEGATIVE (occ hurts cross-sample
> coupling — signal is sub-saturation), as prior = best stability 0.90; **fused prior** non-additive on HE set.
> **NOT BUILT (the remaining plan):** Phase 3 program-wise **hierarchical Bayes / cross-gene pooling** (still
> gene-by-gene frequentist surrogate — Decision A primary); Phase 4 **cooperativity + ceRNA shared-pool**
> (Decision G; `cooperativity.py` / `occupancy.shared_pool_scale` stubs); Phase 6 **state axis Mₜ=M_H+Δ**
> (Decision H; `state.py` stub; GTEx prior + APA/PDUI unwired); Phase 7 universe overlay (RRR/NMF/FM orphan
> μ); Decision I **Shapley attribution** (`attribution.py` stub); `priors.py` π/μ/τ spike-and-slab object.
> **Gaps in DONE phases:** Bar 3 (shuffled-evidence null) not run; ENCORI/POSTAR3/Manakov not folded into the
> ledger; Bar 5 independent-cohort (METABRIC/Buffa) not run on the *learned* M; orphan lane not Bar-5-validated.
> `baselines/` shims + `git add` still pending.
> **Meta-finding:** the simple linear adaptive-lasso is robust and clears every falsification bar; the
> saturation/competition machinery (occupancy κ, fusion) does **not** pay off in cross-sample TCGA coupling —
> caution before Phase 4 (cooperativity/ceRNA sit on the same saturation substrate that failed here).

## 0. The starting problem (measured)

`mirna_hallmark/` has **160 top-level `.py` files** plus `method_dev/` (14), `notebooks/`, `output/`,
`Grant/`, `MIRACLE/`, `presentation/`. Three very different concerns are tangled in one flat namespace:
(a) a **stable data/infra spine** (loaders, config, anchors, deconvolution, CN), (b) the **current
heuristic** (pressure/attribution/coupling), and (c) **~100 one-off analysis/figure scripts** that are
*finished findings* (the outcome arc, spatial, DCIS/EV, CCLE, NMF, grant figures). The learned model is a
*fourth* concern that does not exist yet. The mess is that (a)–(c) are not separated, so it is not obvious
what is load-bearing.

**Guiding principle: separate the four concerns, and do NOT big-bang-move 160 cross-imported files.** The
one-off analyses are *done results* referenced by the catalog/registry/ledger; physically relocating them
risks breaking reproducibility for no modelling gain. So the plan is **classify-first, build-alongside,
archive-lazily** (§4), not a mass `git mv`.

---

## 1. Reuse-vs-rebuild taxonomy

Verdicts: **REUSE** (import as-is, maybe extend config) · **REFACTOR** (keep the ingestion/logic, change
the interface) · **REBUILD** (new; nothing adequate exists) · **BASELINE** (freeze as the nested reference,
do *not* extend) · **ARCHIVE** (one-off; move off the build path eventually).

### 1a. Data / infra spine — **REUSE** (this is the win; it is solid)
| Module(s) | Role | Verdict | Note |
|---|---|---|---|
| `config.py` | paths, strata, cutoffs, `EXCLUDE_NORMAL_LIKE` | REUSE + extend | add a `learned_model` config block; do not fork |
| `data_loaders.py` | clinical, miRNA arms `X`, RNA `Y`, AGO/RISC, target-CN | REUSE | the data-access layer; already honours Normal-like exclusion |
| `hallmark_sets.py` | 50 Hallmark gene sets | REUSE | program scoping unit |
| `stats.py` | BH-FDR, kruskal, hypergeom, z-score | REUSE | |
| `mirna_arm_resolve.py` | arm-name/MIMAT resolution | REUSE | needed for every source join |
| `arm_expression.py` | detectability floor / expressed-arm universe | REUSE | candidate-arm gate |
| `healthy_anchor.py`, `gtex_mirna_matrix.py` | GTEx QN baseline + miTED imputation | REUSE | **required for the nested state model** (`M_H` prior) |
| `estimate_scores.py` | ESTIMATE ssGSEA immune/stromal | REUSE | purity=ESTIMATE + composition |
| `brca_deconvolution.py` | breast scRNA-atlas CIBERSORTx fractions (Cancer Epithelial …) | REUSE | the composition covariate + malignant-compartment proliferation source |
| `mirna_locus_cnv.py` | miRNA-locus CN by stratum | REUSE | **the CN instrument data** for 2SLS |
| `gene_roles.py`, `seed_family_justification.py` | DepMap/TF roles; family collinearity analysis | REUSE | interpretation + informs partial-pooling groups |

### 1b. Evidence / prior layer — **REFACTOR → the PMID-deduped, method-centric ledger** (Design §E)
| Module(s) | Verdict | Note |
|---|---|---|
| `build_edges.py` | REFACTOR | today it **collapses to per-class counts and drops PMIDs**; must instead emit the `(m,g,PMID,assay_class,direction,tissue)` **ledger** upstream of the collapse |
| `evidence_scoring.py` (`score_confidence_logclass`, ENCORI boost) | **SUPERSEDE** | the hand-set, source-centric, coupling-tuned score (Design §E diagnosis). Keep as a *baseline scorer* for comparison; the new fusion is method-centric + learned weights |
| `encori_edges.py`, `fetch_encori_mrna.py` | REUSE (ingestion) | feed the AGO-CLIP bin of the ledger (edge-level rep) |
| `edge_prior_refinement.py`, `edge_breast_context.py` | REUSE (ingestion) | breast-context PMIDs → literature bin |
| **NEW ingestion (REBUILD):** TarBase v9 (`data/miRNA/…TarBase-v9.tsv.gz`), POSTAR3-AGO (`data/RBP-RNA/POSTAR3.parquet`, repurpose `pipeline/lncRNA_interactions/postar3_*`), Manakov chimeric (`data/external_cache/manakov_chimeric/`), TargetScan context++/biochemical K_D | REBUILD | none are wired into the edge prior today |

### 1c. Evaluation / confounder machinery — **REUSE** (the eval half is largely built)
| Module(s) | Verdict | Note |
|---|---|---|
| `coupling_inference.py` | REUSE | the partial-Spearman ladder = OOF coupling test **and** the `C`-block adjustment engine |
| `coupling_permutation.py` | REUSE | shuffled-evidence null (Bar 3) |
| `held_out_tuning.py` | REUSE | patient-level CV folds (OOF discipline) |
| `robustness_checks.py` | REUSE | |
| `core_coupling_deconv_retest.py`, `core_coupling_composition_retest.py` | REUSE | deconvolution/composition covariates = the `C` block; already proved deconvolution > metagene |
| `cross_state_coupling.py` | REUSE | state-resolved coupling |
| `cptac_validation.py`, `buffa_validation.py`, `outcome_metabric.py`, `targetscan_orphan_coupling.py` | REUSE (as OOD) | protein / independent-cohort / orphan layers (Bar 5) |

### 1d. Current heuristic — **BASELINE** (freeze; the model must *nest* and *beat/​match* it)
`pressure_engine.py`, `pressure_build.py`, `hallmark_interaction.py` (CORE), `mirna_state_class.py`,
`geneset_architecture.py`. Also `ago_gate.py` → **partial precursor**: reuse its RISC-capacity estimate,
but the competitive **shared-pool** (ceRNA/promiscuity) is REBUILD. `denominator_attribution_sweep.py`,
`dual_spine_comparison.py`, `hybrid_pressure.py` = the `D(m)`/spine/coupling-tuning experiments →
BASELINE/reference only (they encode the circular tuning the design replaces).

### 1e. The learned model — **REBUILD from scratch** (nothing adequate exists)
occupancy link + threshold gauge + competitive pool · spike-and-slab evidence prior (π, μ, τ) ·
program-wise hierarchical Bayesian regression (+ frequentist adaptive group-elastic-net surrogate) ·
seed-family partial pooling (`mir301_family_depth.py` is a *single-family* precursor to generalize) ·
site-primed cooperative product terms · nested state model `M_T=M_H+Δ` · Shapley identity vs realized
magnitude attribution · learned method-weight calibration · CN 2SLS instrument estimator.

### 1f. The sprawl — **ARCHIVE** (finished one-offs; ~100 files)
`edge_*_panels.py`, `cross_state_*` analyses, `outcome_*.py` (18), `spatial_*.py` (12), `dcis_*.py` (10),
`ev_mirna_*.py` (10), most `cptac_*` one-offs, `mir301_*` deep dive, `nmf_*`/`within_*`, `make_*figure*`,
`mirna_locus_*` plots/galleries/SV, `ccle_*`, `mirna_cnv_*`, `visibility_*`, `subtype_contrasts`,
`stratum_characterization`, `spine_claim_audit`, `escape_mechanism`, `rerun_*`, `run_all.py`, `_build_*`.
These stay runnable but move **off the model build path** (they import the spine, not vice-versa).

---

## 2. What this means in one line
**Reuse the data/infra spine (1a) and the evaluation machinery (1c); refactor the evidence layer into the
ledger (1b); rebuild the model core (1e); freeze the heuristic as the baseline (1d); archive the ~100
one-offs (1f).** The genuinely new code is 1b-ingestion + 1e — a focused new package, not 160 files.

---

## 3. Target architecture (a clean sub-package, spine imported in place)

Create **`mirna_hallmark/learned/`** as the new home; it *imports* the reusable spine from the existing
top-level modules (no move needed to start). Physical relocation of the spine into `core/` is optional and
can follow once `learned/` proves out.

```
mirna_hallmark/
├── config.py  data_loaders.py  hallmark_sets.py  stats.py  mirna_arm_resolve.py     # SPINE (reuse in place)
│   arm_expression.py  healthy_anchor.py  gtex_mirna_matrix.py  estimate_scores.py
│   brca_deconvolution.py  mirna_locus_cnv.py  gene_roles.py
│
├── learned/                         # NEW package — the model
│   ├── evidence/
│   │   ├── ingest/                  #  miRTarBase, tarbase_v9, encori, postar3_ago, manakov, targetscan_kd
│   │   ├── ledger.py                #  (m,g,PMID,assay_class,dir,tissue) build + GLOBAL PMID×assay dedup
│   │   ├── methods.py               #  method-centric bins (reporter/western/qPCR/pert/AGO-CLIP/chimeric)
│   │   └── calibrate.py             #  LEARN method weights vs EXTERNAL repression label (not coupling)
│   ├── priors.py                    #  π (spike-and-slab), μ (site occupancy sum / K_D), τ
│   ├── occupancy.py                 #  a/(a+κ) link + functional-threshold gauge + competitive shared pool
│   ├── families.py                  #  seed-family partial pooling  (identified estimand = family→gene)
│   ├── cooperativity.py             #  site-primed product terms (8–40 nt window)
│   ├── regression.py                #  program-wise hierarchical Bayes  +  freq adaptive-group-EN surrogate
│   ├── state.py                     #  nested  M_T = M_H + Δ  (healthy prior)
│   ├── attribution.py               #  Shapley (identity)  vs  realized contribution (magnitude)
│   └── run.py                       #  orchestrator → learned_model/output/
│
├── learned/eval/                    # eval, mostly thin wrappers over reused machinery
│   ├── coupling.py   → coupling_inference        # OOF partial-Spearman + C block
│   ├── permutation.py → coupling_permutation     # shuffled-evidence null
│   ├── instrument.py                              # NEW: CN 2SLS (data ← mirna_locus_cnv)
│   ├── ood.py        → cptac_validation/buffa/metabric   # protein + independent-cohort
│   └── folds.py      → held_out_tuning
│
├── baselines/                       # FROZEN nested reference (thin re-exports; do NOT extend)
│   → pressure_engine, hallmark_interaction, mirna_state_class
│
└── legacy/                          # (lazy) the ~100 one-off analysis/figure scripts, archived
```

Everything under `learned/` is new or a thin adapter; `baselines/` and `legacy/` are just *classifications*
(can be re-export shims first, physical moves later).

---

## 4. Migration sequence (low-risk; maps to Design §6 phases)

**Step 0 — Classification map (cheapest, do first).** Add `docs/MODULE_MAP.md` tagging each of the 160
modules `spine | evidence | eval | baseline | analysis-oneoff`. This makes the folder navigable **without
moving anything** and is the prerequisite for safe archiving. *Gate:* every module tagged; the ~35
load-bearing ones (1a/1c) identified.

**Step 1 — Stand up `learned/` skeleton + reuse the spine.** Create the package; wire imports to the
existing spine modules; port the MVP path (Design Phase 1: gene-focused adaptive-lasso on the hub panel)
using `coupling_inference` + `held_out_tuning` for the OOF gate. *Gate:* Design Bars 1–2 on the panel.

**Step 2 — The evidence ledger (Design Phase 2, "build the ledger FIRST").** Implement `learned/evidence/`:
ingest miRTarBase (PMIDs) + TarBase v9 + ENCORI + POSTAR3 + Manakov + TargetScan/K_D → `ledger.py` →
global PMID×assay dedup → `methods.py` method-centric bins → `calibrate.py` learned weights vs external
repression label. Refactor `build_edges.py` to emit the ledger. *Gate:* dedup materially changes hub-edge
evidence; Bar 3 (shuffled null); weights beat the hand-set constants OOF on the external label.

**Step 3 — Model core** (`occupancy`, `families`, `regression`, then `cooperativity`, `state`,
`attribution`) per Design Phases 3–6, each behind its own gate.

**Step 4 — Instrument + OOD** (`learned/eval/instrument.py`, `ood.py`) per Design Phases 4–5.

**Step 5 — Archive lazily.** Once `learned/` + `baselines/` cover the build path, move 1f one-offs into
`legacy/` in batches, updating the catalog/registry as each batch moves. Never move a module still imported
on the build path (Step 0 map tells you which).

---

## 5. Two risks to name now
- **Do not fork the spine.** `config.py`/`data_loaders.py` must stay single-source; `learned/` extends
  config via a namespaced block, it does not copy loaders. A fork here re-creates the sprawl.
- **Baselines must stay runnable** — the whole evaluation (Design §5) compares the learned model *against*
  `pressure_engine`/raw-abundance, so freezing ≠ deleting. `baselines/` re-exports, it does not rewrite.

---

## 6. First concrete step
Write `docs/MODULE_MAP.md` (Step 0) and create the empty `learned/` package skeleton with the spine
imports wired — then port the Phase-1 gene-focused MVP. That is the cheapest move that (a) makes the folder
legible and (b) starts the build on the reusable half without touching the 100 one-offs.
