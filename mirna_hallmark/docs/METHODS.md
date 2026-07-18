# Methods — `mirna_hallmark`

> ⚠ **SCOPE (2026-07-18): the PRESSURE construction (§3 miRNA pressure, §3.1–3.5) is RETIRED — but §1 Hallmark
> universe, §2 miRTarBase edges, §2a evidence weight (feeds the ledger), §2b harmonization, and §4 AGO gate
> home methods still in use. Cite BY SECTION; verify against `learned/`. Not a wholly-dead doc.**

> ⚠ **THE PRESSURE-HEURISTIC ESTIMATOR SPECCED IN §3 IS RETIRED — read this before citing a §3 number (MH-115, 2026-07-13).**
> The evidence-weighted **pressure heuristic** (`compute_gene_pressure`, the `edge_w × expr_mult` construction
> below) was benchmarked against the learned β on the SAME design matrix, varying only the weight:
> on CPTAC protein — abundance −0.006 (p=0.20) · **heuristic −0.008 (p=0.11)** · **learned β −0.036 (p=5.3e-13)**;
> on CPTAC mRNA the heuristic is **p=0.83 (exactly zero)** and abundance has the **WRONG SIGN** (+0.014, p=1.9e-03);
> and even **in-sample on TCGA, raw abundance (−0.056) BEATS the heuristic (−0.046)** ⇒ **the evidence weighting
> adds nothing over raw abundance, ever.** Per-edge: β is entirely **additional** to abundance (β|abundance
> −0.140, p=1.5e-13) while abundance is **redundant** given β (abundance|β −0.012, p=0.52).
> ⇒ **This doc is the record of a BASELINE, not of the model.** The model is ONE dense learned-τ² non-negative
> Bayesian posterior per gene, two readouts — `MODELING_FRAMEWORK.md` → *THE CANONICAL MODEL STATEMENT*, spec in
> `LEARNED_MODEL_METHODS.md`. **Do not build on the construction below.**
> ✅ **Still valid and still cited:** the DESIGN decisions here are estimator-independent — the Hallmark universe,
> the edge universe, arm-name harmonization (→ `DATA_SOURCES`), the confounder ladder (→ `FORMULAS §7`), and
> **two analyses the learned model has never done**, now on the board **§Y**: the **AGO/RISC capacity gate as an
> INTERACTION** on β (§4 — never tested, by either arc), and the **Freedman–Lane permutation null** (§6a —
> which the learned model still lacks; MH-123 measured its t-null as 3–4× too narrow).



Formal description of the subproject's computations. Cohort and join
conventions follow the parent repo (`analysis/COHORT_AND_JOIN_CONVENTIONS.md`):
TCGA-BRCA primary tumors, 12-char participant key, multi-vial collapse by mean
(expression) / median (CNV).

## 1. Hallmark gene universe

The 50 MSigDB Hallmark sets are parsed from
`annotations/GSEA/h.all.v2026.1.Hs.txt` (`hallmark_sets.py`). The analysis
**universe** is the union of all member genes (~4,384). Membership is
many-to-many: a gene in *N* sets is represented in *N* contexts. The
`(miRNA, gene, hallmark_set)` edge table therefore repeats a gene's evidence
once per set it belongs to; cohort-level counts dedup on `(miRNA, gene)`.

## 2. miRTarBase edges (vectorized)

`build_edges.py` reuses the fast normalizer
`pipeline.genes.mirtarbase.load_mirtarbase` (gene-panel = Hallmark universe),
then computes the `(miRNA, gene)` interaction summary with **vectorized** groupby
aggregations (`compute_interaction_summary_fast`):

- explode to `(miRNA, gene, study, experiment_class, support_type)`;
- per-study boolean flags via `groupby(...).max()` (OR within a study);
- experiment × support **cross-flags** at study level;
- sum study-level flags to per-`(miRNA, gene)` counts (`n_*_studies`,
  `n_<exp>__<support>_studies`, `n_studies`);
- `evidence_score = 3·reporter + 3·binding + 2·protein + 1·rna + 1·perturbation`
  (study counts), matching the canonical formula.

The output schema matches the canonical `mirtar_interaction_summary.csv` columns
consumed by `analysis.expression.mirna_target_integration.load_mirtar_edges`, so
the pressure engine is drop-in compatible. The reference (slow) builder is
available via `--use-canonical-builder` for parity checks.

## 2a. S1 effective evidence weight (pressure spine default)

When loading edges for pressure computation (`pressure_build.load_mirtar_edges →
evidence_scoring.build_m0_edges`), the CSV `evidence_score` field is **replaced**
at runtime by the `confidence_logclass` scorer (`config.PRESSURE_EVIDENCE_SCORER`):

```
confidence_logclass(m, g)
  = w_rep  · log1p( n_reporter__functional_mti )
  + w_prot · log1p( n_protein__functional_mti  )
  + w_rna  · log1p( n_rna__functional_mti      )
  + w_pert · log1p( n_perturbation__functional_mti )
  + w_bind · log1p( n_binding__functional_mti  )
  + δ · [same sum over __functional_mti_weak columns]
```

`w_rep=3.0, w_prot=2.5, w_rna=1.5, w_pert=1.5, w_bind=0.5`; `δ=0.3`.
Binding-only CLIP earns 0.5 (proximity ≠ demonstrated repression); `log1p`
within each class gives diminishing returns, defusing publication-count bias.

**ENCORI AGO-CLIP depth boost** (`PRESSURE_ENCORI_ALPHA = 0.5`): for the ~25%
of M0 pairs with an ENCORI entry:

```
enc_depth(m, g) = 2·log1p(clipExpNum) + 1·log1p(degraExpNum)
                + 0.5·I(TargetScan ∧ miRanda ∧ PITA) + 0.25·log1p(pancancerNum)
evidence_score_eff = confidence_logclass + 0.5 · enc_depth   (0 if no ENCORI row)
```

Decided 2026-06-22 (`dual_spine_comparison` S1 grid; Basal coupling ≥ 42/50,
key 8/8 held). Source: `evidence_scoring.py`; applied in `pressure_build.py`.

**High-evidence flag** (`config.HIGH_EVIDENCE_*`): `n_functional_mti_studies ≥ 1`
AND (`n_reporter__functional_mti_studies ≥ 1` OR
`n_protein__functional_mti_studies ≥ 1`); fallback to `evidence_score ≥ 3` when
cross-count columns are unavailable.

## 2b. Mature-arm name harmonization (miRTarBase → GDC expression)

miRTarBase and TargetScan use miRBase-style mature-arm names. TCGA/GDC miRNA-seq
row labels often add **isoform suffixes** (e.g. `hsa-miR-383-5p` in miRTarBase vs
`hsa-miR-383-5p.1` in the expression matrix). Before 2026-06, unmapped arms were
**silently dropped** from pressure (`pressure_engine` filters to
`mirna.index`), which could zero out validated edges (e.g. miR-383-5p → IRF1).

Module: `mirna_arm_resolve.py`. Applied when loading miRTar edges for **pressure
and downstream miRNA→target tests** (`pressure_build.load_mirtar_edges`,
`resolve_edges_mirna`). Each resolution is tagged (`match` column) and written to
`output/logs/mirna_arm_resolve_audit.tsv`.

### Default policy (decided 2026-06-16)

| Tier | Rule | Rationale |
|------|------|-----------|
| **A — always** | Exact matrix name | No transform |
| **A** | Isoform suffixes `.1` / `.2` / `.3` on full name and on `-3p`/`-5p` terminus | Same mature product; GDC MIMAT isoform labels |
| **A** | Star-suffix strip (`*` → base) | miRBase star notation vs matrix |
| **A** | Loci-file map: `mirbase_mature_id` / MIMAT accession → expressible row (per ID only) | Curated annotation join |
| **A** | miRBase ``mature.fa`` (human ``hsa-*`` mature IDs) for MIMATs absent from loci | GDC rows left as raw MIMAT without loci entry (~1.4k arms) |
| **B — guarded** | Precursor stem without `-3p`/`-5p` → append mature arm **only when exactly one** candidate exists in the matrix (e.g. `544a` → `544a-3p`) | miRTar sometimes omits arm suffix; ambiguous if both arms measured |
| **C — excluded** | Opposite mature arm (`-3p` ↔ `-5p`) | Different molecules (different seeds/targetomes); only fixes ID errors at cost of wrong abundance |
| **C** | Letter slim (`151a` → `151`, `21a` → `21`) | Paralog / locus conflation risk |
| **C** | Pre-miRNA sibling collapse (map all mature IDs from one hairpin to one expressible row) | Same risk as opposite-arm flip |

**Hallmark M0 edge table impact** (miRTar evidence ≥ 2, permissive weights):
730 edge rows referenced miRTar arms absent from the expression matrix → **245
recovered** under Tier A+B → **485 rows** still missing (arms not measured in GDC).

### Where harmonization applies vs does not

**Applies (miRTarBase name → GDC expression row):**

- Combined gene pressure (`pressure_build`, `hallmark_interaction`)
- Hybrid modes M0/M7/M8/M11 (`hybrid_pressure`)
- Target-level anti-correlation (`target_combined_anticorr`)
- Hub-route OLS explainability (`gene_expression_explainability`)
- Robustness / hub partial Spearman (`robustness_checks` via `load_mirtar_edges`)

**Does not apply (separate join logic):**

- **CN→miRNA expression concordance** (`mirna_locus_cnv.py`): joins
  `entity_label` from the miRNA entity catalog directly to `load_mirna_arms().index`;
  arms not in the matrix are skipped (`if arm not in expr_arms: continue`). The
  catalog is built from the same loci/GDC ecosystem, so isoform mismatches are
  rare here; ~22 loci arms are genuinely absent from GDC and cannot be harmonized.
- **CCLE cell-line CN→expression** (`ccle_mirna_cn_concordance.py`): uses a
  **separate**, more permissive CCLE NanoString alias map (`151a→151-3p`, precursor
  stems, etc.) — not the miRTar harmonizer above.

**Indirect link to dosage work:** CN→expression establishes *which arms show dosage
propagation*; harmonization affects how miRTarBase *target edges* for those arms
enter **pressure → target RNA** tests (`concordance_target_join`,
`target_combined_anticorr`). It does **not** change the primary Spearman(CN, arm
expression) statistics in `mirna_cnv_expr_concordance.tsv`.

## 3. miRNA pressure

Per-(gene, sample) pressure uses ``pressure_build.compute_gene_pressure`` (wrapper
around ``pressure_engine``). **Expression weighting is identical across edge
universes**; modes differ only in which `(miRNA, gene)` pairs enter the sum and
how each pair's **structural edge weight** is set.

### 3.1 Expression weighting (all modes)

Default spine (``config.PRESSURE_*``, M0 modules):

```
edge_w(m→g) = log1p(evidence) / log1p(Σ_g log1p(evidence_m→g))   [target_norm=evidence_mass]
share(m,s)  = softmax( log2(RPM+1)_m,s − cohort_median_m ) among arms → g
pressure(g,s) = Σ_m  share(m,s) · z_m(s) · log2(RPM+1)_m,s · edge_w(m→g)   [expr_mode=softmax_z_logrpm]
```

(`evidence` = `evidence_score_eff` per §2a — i.e. `confidence_logclass + α·enc_depth`,
not the raw CSV study-count field.)

Default hybrid (``config.PRESSURE_HYBRID_*``, M7+):

```
edge_w from §3.3 fused weight / log1p(Σ_g [log1p(ev)+log1p(ts)])   [target_norm=combined_mass]
pressure(g,s) = Σ_m  share(m,s) · z_m(s) · edge_w   [expr_mode=softmax_z]
```

- ``structural_weight`` — set by edge universe (§3.2–3.3); passed as
  ``evidence_score`` into ``pressure_engine``.
- ``z_m(s)`` — cohort z-score of mature-arm log2(RPM+1).
- ``target_norm=evidence_mass`` — down-weights miRNAs with heavy total miRTar evidence mass.
- ``aggregate=sum`` — sum over targeting arms (``mean`` weakens signal; not default).
- Legacy z-only weighting preserved in ``pressure_sensitivity`` as ``z_sum_baseline``.

**Extended target normalization** (``pressure_sensitivity``, optional M7 runs in
``hybrid_param_sensitivity``):

| ``target_norm`` | Denominator per miRNA *m* |
|-----------------|---------------------------|
| ``evidence_mass`` | ``log1p( Σ_g log1p(evidence_{m→g}) )`` |
| ``ts_mass`` | ``log1p( Σ_g log1p(ts_weight_{m→g}) )`` (requires TS column on edges) |
| ``combined_mass`` | ``log1p( Σ_g [log1p(ev) + log1p(ts)] )`` — miRTar + TargetScan evidence mass |

**Absolute abundance anchoring** (multiply ``softmax_z`` share × z):

| ``expr_mode`` | Extra term |
|---------------|------------|
| ``softmax_z_logrpm`` | × ``log2(RPM+1)`` |
| ``softmax_z_absratio`` | × ``max(0.25, RPM/median_cohort)`` |
| ``softmax_z_blend`` | z multiplier uses ``0.5·share + 0.5·uniform`` among targeting arms |

M7 hybrid edges store **raw** ``evidence_score`` + ``ts_weight`` alongside fused
``edge_weight``; ``degree``/`none` use the fused weight, mass norms use raw columns.

**Sensitivity grids** (2026-06-16, post arm-resolve): ``pressure_sensitivity`` (M0
edges + TS merge) and ``hybrid_param_sensitivity`` (one-at-a-time ``β``, ``γ``,
``min_ts`` + weighting variants). Best Basal Hallmark counts on M0:
``softmax_z_logrpm`` + ``evidence_mass`` (42/50 neg-sig, 8/8 key). TS-only mass
norm collapses coupling (11/50). M7 ``combined_mass`` reaches **41/50** Basal in
default ``hybrid_pressure`` (param-grid peak **43/50**) and restores borderline
IRF1 hub partial (ρ≈−0.15–0.17, p≈0.025–0.045) vs M7 on ``degree`` norm
(ρ≈−0.11, p≈0.15) after isoform harmonization.

**AGO gate** (§4). TS-only / orphan arms additionally require cohort-median
log2(RPM+1) ≥ ``PRESSURE_TS_ABUNDANCE_FLOOR`` (default 1.0).

### 3.2 TargetScan sequence weights

TargetScan 7.1 context++ scores (`robustness_checks._load_targetscan_weights`):

```
ts_weight(m→g) = Σ |weighted context++ score|   (human hsa- rows; Hallmark genes)
```

Study-count independent (pure sequence/site model). Default floor for hybrid modes:
``min_ts = 0.25`` (`hybrid_pressure.HybridPressureSpec`).

On the Hallmark gene universe (~4,384 genes), TargetScan yields:

| Filter | (miRNA, gene) pairs | miRNAs | Genes |
|--------|--------------------:|-------:|------:|
| any TS score | 52,170 | 320 | 3,247 |
| \|TS\| ≥ 0.25 | 14,336 | 318 | 2,385 |
| \|TS\| ≥ 0.50 | 2,161 | 257 | 741 |
| \|TS\| ≥ 0.75 | 251 | 120 | 109 |
| \|TS\| ≥ 1.00 | 47 | 41 | 19 |

### 3.3 Edge universes and miRTar + TargetScan fusion

Modes are implemented in ``hybrid_pressure.build_hybrid_edges``; all resolved
through §2b before scoring.

| Mode | Edge set | Structural weight `w(m→g)` | Default use |
|------|----------|------------------------------|-------------|
| **M0** | miRTarBase, ``evidence_score ≥ 2`` | ``log1p(evidence)`` | **Spine**: `hallmark_interaction`, `target_combined_anticorr`, robustness, explainability |
| **M7** | HE miRTar + capped TS fill-in | see below | Hybrid primary (`config.PRESSURE_HYBRID_PRIMARY`); sensitivity vs M0 |
| **M8** | M0 pressure + λ·orphan track (separate sum) | dual-track | Sensitivity |
| **M11** | Top-k TS per gene (all TS pairs) | ``log1p(ts_weight)`` | Pure-sequence sensitivity |

**M0 spine (miRTar-only):** 5,315 `(miRNA, gene)` rows → 791 arms → 1,438 Hallmark
genes. **No TargetScan-only pairs enter pressure.**

**M7 tiered fusion** — for each selected `(m, g)`:

| Tier | Inclusion rule | `w(m→g)` |
|------|----------------|----------|
| `he_ts_boost` | high-evidence miRTar **and** `ts_weight ≥ min_ts` | ``log1p(evidence) · (1 + β · log1p(ts))`` (β=0.5) |
| `he` | high-evidence miRTar, no TS boost | ``log1p(evidence)`` |
| `gap` | not HE, **no** miRTar study (`n_studies=0`), TS orphan | ``γ · log1p(ts_weight)`` (γ=0.25) |
| `orphan` | not HE, miRTar row exists but below HE | ``γ · log1p(ts_weight)`` |

Orphan/gap arms are capped per gene: **top 5** TS arms for program genes, **top 15**
for the eight hub targets (`HUB_ROUTES`). M7 yields ~13,077 edges (827 arms, 2,791
genes) vs M0's 5,315.

**M11:** top-k TS arms per gene (same caps), ``w = log1p(ts_weight)`` — ~8,171 edges.

### 3.4 TS-only pairs excluded from the default spine

**Orphan definition** (`targetscan_orphan_coupling`): `ts_weight ≥ 0.25` and **not**
a high-evidence miRTarBase edge.

| Category | Pairs | Notes |
|----------|------:|-------|
| TS ≥ 0.25, not HE (full orphan pool) | **13,245** | 318 miRNAs, 2,352 genes — **all absent from M0** |
| … of which no miRTar row at all | 9,866 | pure sequence predictions |
| … of which miRTar exists but below HE | 3,379 | curated low-confidence |
| TS orphan, \|TS\| ≥ 0.50 | 1,867 | strong sequence support, still M0-excluded |
| TS orphan, \|TS\| ≥ 0.75 | 194 | |
| TS orphan, \|TS\| ≥ 1.00 | 26 | |

M7 adds **~7,840** capped orphan/gap pairs (after per-gene top-k); ~1,158 orphan
rows still fail GDC arm harmonization (§2b). Uncapped orphan pool remains in
``targetscan_orphan_coupling`` for gene-level coupling screens, not in M0 pressure.

### 3.5 Module → edge universe map

| Module | Edge universe |
|--------|---------------|
| `hallmark_interaction` | **M0** miRTar only |
| `target_combined_anticorr` | **M0** miRTar only |
| `robustness_checks` (main pressure) | **M0** miRTar only |
| `gene_expression_explainability` | **M0** route edges (+ IRF1 TS alt arms in OLS terms only) |
| `hybrid_pressure` | M0 / **M7** / M8 / M11 (config default modes: M7, M8, M11) |
| `targetscan_orphan_coupling` | TS orphans only (coupling tests, not full Hallmark pressure spine) |

``config.PRESSURE_HYBRID_PRIMARY = "M7"`` names the preferred **hybrid** mode for
extended runs; the cohort **headline spine** remains M0 unless a module explicitly
calls ``hybrid_pressure``.

## 4. AGO/RISC availability gate

`ago_gate.py`. **Biological rationale.** Seed-mediated repression is not free: a
miRNA only acts once it is (i) loaded into an Argonaute to form the miRISC core
and (ii) the loaded complex engages a TNRC6/GW182 effector that recruits the
degradation machinery. Both steps draw on a **finite, shared** pool that every
miRNA in the cell competes for, which is exactly why a per-sample,
miRNA-independent *availability* term is the right shape for a global gate.

We model the two steps as two both-required components, from RNA log2(TPM+1)
z-scored per gene across the cohort:

1. **Loadable miRISC core — `ago_load(s)`.** AGO1–4 all bind mature miRNAs and
   support canonical seed-mediated, TNRC6-dependent repression. **AGO2** is the
   predominant and most abundant Argonaute in most somatic cells and the only one
   with slicer (endonuclease) activity, so it carries the majority of the loaded
   pool; AGO4 is typically lowest-expressed. We therefore use an
   abundance/dominance-weighted mean rather than a flat mean:

   ```
   ago_load(s) = Σ_g w_g · z_g(s) / Σ_g w_g      over g in {AGO1..AGO4}
   w = {AGO1: 1.0, AGO2: 2.0, AGO3: 1.0, AGO4: 0.5}   (transparent prior, swept in §sensitivity)
   ```

   (AGO protein is stabilised by miRNA loading — unloaded AGO is degraded — so AGO
   abundance tracks loadable capacity. Slicing per se matters mainly for the rare
   extensively-complementary site, not the seed-match bulk; the AGO2 up-weight
   reflects abundance dominance, not catalysis.)

2. **Effector engagement — `effector(s)`.** AGO occupancy alone does **not**
   repress: the TNRC6A/B/C (GW182) scaffolds dock onto AGO via GW/WG motifs and
   recruit the PAN2–PAN3 and CCR4–NOT deadenylases → poly(A) shortening → DCP1/2
   decapping → XRN1 5′→3′ decay (plus translational repression). TNRC6 is thus the
   rate-limiting *effector*. `effector(s) = mean z of TNRC6A/B/C`.

Because both resources are required, capacity is **co-limited** — whichever is
scarcer rate-limits repression (Liebig's law of the minimum), not their sum (a
sum would let abundant AGO compensate for absent effector, which the biology
forbids):

```
capacity(s) = min( ago_load(s), effector(s) )        (effector_colimit, default)
            = 0.5·(ago_load + effector)               (compensatory blend — sensitivity variant)
            = ago_load(s)                              (--no-tnrc6 — sensitivity variant)
```

The bounded saturating gate maps capacity to a multiplier in `[gate_min, 1]` via
a logistic (sigmoid `σ(x) = 1 / (1 + e^-x)`):

```
gate(s) = gate_min + (1 - gate_min) · σ( k · (capacity(s) - c0) )
```

with defaults `gate_min = 0.5`, `k = 1`, `c0 = 0`. The floor `gate_min` encodes
that even a low-capacity tumour retains basal repression — the layer can attenuate,
never erase. At joint-typical capacity (both components at the cohort mean) the
gate is `(1 + gate_min) / 2 = 0.75`; the cohort-median gate sits slightly below
(≈0.74, mean ≈0.73 on TCGA-BRCA) because the *minimum* of two mean-zero components
is negatively biased — i.e. needing both resources high is genuinely harder than
needing one. Gated pressure = ungated pressure × gate(s).

> The gate is a **documented sensitivity layer, not a causal model**. AGO/TNRC6
> *mRNA* is an imperfect readout of functional RISC capacity (loading, AGO2
> phosphorylation/hydroxylation, and effector recruitment are all
> post-transcriptional), so the bounded form is deliberately conservative. Every
> interaction result is reported gated **and** ungated; empirically the gate
> rescales rather than reorders Hallmarks. Treat gated/ungated differences as
> sensitivity, not mechanism.

## 5. Hallmark roll-ups

- **Hallmark pressure** (hallmark × sample) = mean per-gene gated pressure over a
  Hallmark's member genes present in the pressure matrix.
- **Hallmark target expression** (hallmark × sample) = mean RNA log2(TPM+1) over
  present member genes.

## 6. Anti-correlation (interaction)

Per Hallmark, Spearman(pressure, target expression) across shared samples
(min n = 20); expected **negative**. A **partial Spearman** adjusting CPE and HRD
(`analysis.utils.common.loaders.partial_spearman`) is reported alongside the raw
coefficient. Tests are FDR-corrected (Benjamini-Hochberg) within each view
(gated / ungated). The same test is repeated within each PAM50 subtype.

## 6a. Dependence-aware FDR + permutation null (all resolutions)

`coupling_inference.py` / `coupling_permutation.py`. The §6 partial-Spearman is
re-tested at **every** rung (edge, gene, target-set/miRNA, program, universe) with
two additions that the plain BH headline lacks.

**Two senses of "family."** The word is overloaded; both are handled:

1. *Statistical testing family* — the set of tests sharing one BH correction. BH
   assumes independence, so we additionally report **Benjamini–Yekutieli** (`q_by`),
   valid under arbitrary dependence: BH with the threshold divided by the harmonic
   number `c(m) = Σ_{i≤m} 1/i`.
2. *miRNA seed family* — paralogous arms sharing a TargetScan seed (miR-29a/b/c, …).
   They co-vary and target the same sites, so their edge tests are **non-independent**.
   This is justified empirically, not assumed (`seed_family_justification`): within-
   seed-family arm pairs co-vary ≈2× random pairs (median Spearman 0.38 vs 0.17,
   Mann–Whitney p≈4e-14; shared-target Jaccard 0.05 vs 0.0), and the multi-arm hubs are
   20–37% redundant (Nyholt effective-test count). Where the test unit is a miRNA
   (edge / target-set) each seed family is collapsed to one **Simes** p, BH'd across
   families (`q_seedfamily`), so an `n`-paralogue hub spends one error slice; a per-family
   **min-statistic permutation FWER** (`p_seedfamily_fwer`) accompanies it. Discoveries
   are counted in *families*, not edges. Singleton-seed arms are unchanged.

**Permutation null.** Each coupling is re-tested by **Freedman–Lane** permutation: rank-
transform predictor, target, and covariates; OLS-residualize predictor and target ranks on
the covariate ranks (intercept + CPE + HRD + TCGA batch); the partial Spearman is the
correlation of the unit-scaled residuals. A null is built by permuting the predictor
residual's sample order, with **one shared permutation per draw across all rows of a
resolution**, so seed-family correlation is preserved (and a per-family min-statistic null
is available). Two p-values are emitted: the discrete exceedance `p_perm` (floor `1/(B+1)`,
the honest empirical check) and a smooth tail `p_z = Φ(z)` from the null moments (the
permuted statistic is ~Gaussian by CLT) — the FDR rides `p_z` so multiplicity is not capped
by the permutation floor. Cohort = primary tumours with PAM50 Normal-like excluded
(`EXCLUDE_NORMAL_LIKE`), n ≈ 1041. The per-gene target-CN rung of the §6 ladder stays in the
parametric path (a row-specific covariate breaks the shared-permutation vectorization).

**Universe.** The "how many of 50 programs couple" claim gets its own null: from the joint
program-level permutation, the per-draw count of programs crossing their own α-quantile
gives a null count distribution, and the observed coupled-program count is scored against it.

## 7. High-evidence target enrichment

Per Hallmark, a one-sided hypergeometric test asks whether the Hallmark's genes
are over-represented among **high-evidence miRNA target genes**, relative to the
Hallmark universe. With:

- `N` = universe size (genes in any Hallmark set),
- `K` = total high-evidence target genes in the universe,
- `n` = Hallmark size (member genes in the universe),
- `k` = Hallmark genes that are high-evidence targets,

the enrichment p-value is `P(X ≥ k)` for `X ~ Hypergeometric(N, K, n)`, and
`fold_enrichment = k / (n · K / N)` (observed / expected). FDR-corrected (BH).

## 8. miRNA-universe CNV + concordance

`mirna_locus_cnv.py` reuses `dosage_landscape_cnv` with `panel_arm_ids=None`
(full miRNA universe) to build a long sample × entity (locus / arm / family) CNV
table from ASCAT3 (primary tumors), summarized by the four strata via
`summarize_by_strata`.

**MIMAT arm CN:** aggregated from **paralog hairpin loci** (`cnv_miRNA.csv`
``mature_accessions`` → MIMAT). **Locus CN** uses ASCAT3 **allelic segments**
by default (`ANALYSIS.mirna_cnv_dosage.mirna_cnv_source=ascat3_segment`).

**Paralog weighting (expression_weighted; full rationale:**
`MIRNA_CNV_DOSAGE_REPORT.md` §2): cohort-median locus×MIMAT RPM → `weight_ref`;
normalize to `w_norm` per MIMAT. Tumor: `MIMAT_CN = Σ w_norm × CN_locus`.
Healthy reference: same weights, all loci diploid → composite **2** (not the
sum of paralog copies). Weights are cohort-fixed (not sample-specific) to avoid
CN→expression circularity.

Methods via ``ANALYSIS.mirna_cnv_dosage.mimat_cn_aggregate``:

| Method | Rule |
|--------|------|
| `median` (default) | Median locus CN across paralogs |
| `mean` / `max` | Mean or max locus CN |
| `expression_weighted` | RPM-weighted mean using **cohort-median** locus×MIMAT reference weights (same for all samples; avoids circularity in CN→expression concordance) |

**Focus-arm subtype tables** (`mirna_cnv_subtype_depth.py`; dosage report §5):
per arm × PAM50, median participant CN with a 95% bootstrap CI on that median
(2,000 percentile-bootstrap draws; all participants included, equal weight per
participant; paralog `w_norm` aggregation applied before summarization). See
`MIRNA_CNV_DOSAGE_REPORT.md` §5 for interpretation alongside Range / % amp.

**CN→expression concordance figures:** `mirna_locus_cnv.py` writes
`mirna_cnv_expr_concordance_scope_audit.tsv` and
`figures/mirna_cnv_expr_concordance_rho_boxplot.png` (marginal vs partial ρ;
full cohort, not gain-filtered). See dosage report §7.1.

Build isoform caches once:

```bash
.venv/bin/python3 scripts/mirna/build_mirna_locus_mimat_expression.py
# → data/miRNA/locus_mimat_rpm.tsv.gz
# → data/miRNA/locus_mimat_weight_reference.tsv
# → data/miRNA/locus_gdc_hairpin_map.tsv
```

**Locus mapping tables compared**

| Step | Table A | Table B | Overlap |
|------|---------|---------|---------|
| ASCAT3 CN at locus | MirGeneDB `MI*` (`cnv_miRNA.csv`) | ASCAT3 **allelic segments** | Max bp overlap → `cn_total` |
| ASCAT3 CN (legacy) | Same MirGeneDB `MI*` | ASCAT3 gene intervals | Max bp overlap → ENSG → copy number |
| Hairpin map (weights) | Same MirGeneDB `MI*` | GDC isoform hairpin coords (`hsa-mir-*`) | Max bp overlap on same chr → GDC hairpin id |
| MIMAT aggregation | Locus CN rows | `mature_accessions` on `cnv_miRNA.csv` | Logical map MI* → MIMAT (no geometry) |

~348/506 MirGeneDB loci map to a GDC hairpin; **158 unmapped** have a chromosome match but **zero bp** overlap with any GDC hairpin interval (annotation discordance). Those loci still receive locus-level CN; they contribute to MIMAT CN only via unweighted methods or with zero reference weight.

Run overlap audit: ``scripts/mirna/audit_mirna_locus_maps.py`` → ``data/miRNA/audit/``.

**CN→expression concordance**: per mature arm,
Spearman(median copy number, arm log2(RPM+1)) across participants (min n = 20),
with a CPE+HRD partial Spearman and BH FDR.

> **Arm naming:** this test joins catalog `entity_label` directly to GDC expression
> rows (§2b). It does **not** use the miRTarBase harmonizer. Isoform suffix
> mismatches are uncommon in the catalog; unmeasured arms are dropped. Target
> repression tests that consume miRTarBase edges **do** use §2b harmonization.

## 9. Stratum characterization

`stratum_characterization.py`. For each Hallmark × stratum level:
- **target-gene CNV burden**: mean ASCAT3 integer CN and % gain/% loss/% deep-del
  over member genes;
- **regulatory-miRNA expression**: mean/median log2(RPM+1) of the high-evidence
  miRNAs targeting the Hallmark.

## 10. Statistics helpers (`stats.py`)

`bh_fdr` (NaN-safe Benjamini-Hochberg), `kruskal_across_strata` (KW across
levels, min group size), `hypergeom_enrichment` (one-sided), `zscore_rows`.

## 11. CN state thresholds

Integer copy number → state: `deep_del` (0), `loss` (1), `neutral` (2),
`gain` (3), `amp` (≥4) — matching the parent CNV module.
