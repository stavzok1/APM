# Formulas reference — `mirna_hallmark`

> ⚠ **SCOPE (2026-07-18): the PRESSURE-HEURISTIC CONSTRUCTION here (§1 backbone · §3 expression · §4 modes ·
> §5 share) is RETIRED — but this doc is NOT wholly dead: §7 coupling ladder, §8 family aggregates, §11
> state-trajectory (the LIVE `mirna_state_class` module), §2 evidence weight, and §2b harmonization home formulas
> still in use. Cite BY SECTION and verify against `learned/` — do not treat the whole doc as a baseline record.**

> ⚠ **THE PRESSURE-CONSTRUCTION ESTIMATOR SPECCED IN §1–5 IS RETIRED — read this before citing a §1–5 number (MH-115, 2026-07-13).**
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



A single, detailed, *explainable* reference for every quantity the subproject
computes: the pressure backbone, every expression mode, every share / specificity
metric, the coupling (anti-correlation) ladder, the family-level aggregates, and
the orphan validation-priority composite. Each block gives the **formula** (in
plain ASCII so it renders in any viewer), the **plain-language meaning**, **what it
is used for**, and the **code location**.

This document is descriptive of the code as written — if you change
`pressure_engine.py` / `robustness_checks.py`, update here too.

---

## 0. Notation

| Symbol | Meaning |
|--------|---------|
| `m` | a miRNA **arm** (e.g. `hsa-miR-301a-3p`) |
| `g` | a target **gene** |
| `s` | a **sample** (participant) |
| `x[m,s]` | arm expression = `log2(RPM + 1)` of arm `m` in sample `s` |
| `R(g)` | the set of arms that **regulate** `g` (have a miRTarBase edge to it) |
| `mean_m`, `sd_m` | mean / SD of `x[m,:]` across the samples in scope |
| `med_m` | cohort **median** of `x[m,:]` |
| `targets(m)` | all genes arm `m` regulates (used in specificity denominators) |
| `mean_s(.)` | average over samples; `sum_m(.)` = sum over arms |

**"In scope" matters:** every mean / SD / median / softmax is recomputed on *only
the samples in the current scope* (whole cohort, or one PAM50 subtype). That is why
subtype tables are not just cohort numbers re-sliced.

---

## 1. Pressure backbone

The one master formula. Everything else is a piece of it or a summary of it.

```
pressure(g, s) = AGG over m in R(g) of [ expr_mult(m, s) * edge_w(m, g) ]

                 |________________|   |___________|   |___________|
                  aggregate (sum)      expression       evidence
                                       piece (S3-4)     piece (S2)
```

- `AGG` = **sum** (default) or **mean** over the gene's regulators (`config.PRESSURE_AGGREGATE`).
- The per-arm term inside the bracket is the **contribution**:

```
c(m, g, s) = expr_mult(m, s) * edge_w(m, g)
```

  **Every share, specificity, and family number in S5-S9 is built from `c`.**

*Code:* `pressure_engine.compute_gene_pressure` (sums over arms);
`compute_gene_pressure_contributions` (keeps the per-arm `c`).

---

## 2. Evidence piece — `edge_w(m, g)`

Edge weight starts from curated evidence, then is optionally divided by a per-arm
**normalization denominator** `D(m)` that controls promiscuity. Selected by
`config.PRESSURE_TARGET_NORM` (default `evidence_mass`).

**S1 spine (default, `config.PRESSURE_EVIDENCE_SCORER = confidence_logclass`).**
`evidence_score` passed to the engine is **not** the raw CSV study-count field
but a quality-weighted recompute from the per-experiment-class × support-type
cross-count columns:

```
confidence_logclass(m, g)
  = w_rep  · log1p( n_reporter__functional_mti )
  + w_prot · log1p( n_protein__functional_mti  )
  + w_rna  · log1p( n_rna__functional_mti      )
  + w_pert · log1p( n_perturbation__functional_mti )
  + w_bind · log1p( n_binding__functional_mti  )
  + δ · [same sum over __functional_mti_weak columns]
```

weights: `w_rep=3.0, w_prot=2.5, w_rna=1.5, w_pert=1.5, w_bind=0.5`;
`δ = 0.3` (curator-flagged Weak discount). Binding-only CLIP earns 0.5 vs
reporter/protein (proximity ≠ demonstrated repression); `log1p` within each
class gives diminishing returns — defuses publication-count inflation.

ENCORI AGO-CLIP depth boost (`PRESSURE_ENCORI_ALPHA = 0.5`): for the ~25% of
M0 pairs that have an ENCORI entry (`output/encori/encori_mirna_mrna_targets.parquet`):

```
enc_depth(m, g) = 2.0 · log1p(clipExpNum)
                + 1.0 · log1p(degraExpNum)
                + 0.5 · I(TargetScan ∧ miRanda ∧ PITA)
                + 0.25 · log1p(pancancerNum)

evidence_score_eff(m, g) = confidence_logclass(m, g)  +  α · enc_depth(m, g)
                                                            (α = 0.5; 0 if no ENCORI row)
```

Decided 2026-06-22 (`dual_spine_comparison` S1 grid; Basal coupling ≥ 42/50
held, 8/8 key Hallmarks).

*Code:* `evidence_scoring.score_confidence_logclass`,
`evidence_scoring.apply_encori_boost_to_edges`, `pressure_build.load_mirtar_edges`.

Base edge strength:

```
w0(m, g) = log(1 + evidence_score_eff(m, g))
```

(The `combined_mass` hybrid modes additionally add a TargetScan term
`log(1 + ts_weight)` — see METHODS.md §3.3.)

Per-arm denominator `D(m)` by mode:

```
none           D = 1                                              (raw evidence)
degree         D = log(1 + number_of_targets(m))                 (penalize wide targeting)
evidence_mass  D = log(1 + SUM over g' of log(1+evidence(m,g')))  (DEFAULT: penalize total curated mass)
ts_mass        D = log(1 + SUM over g' of log(1+ts_weight(m,g'))) (TargetScan-mass version)
combined_mass  D = log(1 + SUM over g' of [log_ev + log_ts])      (miRTar + TS combined)
```

Final edge weight:

```
edge_w(m, g) = w0(m, g) / max(D(m), 1e-9)
```

**Meaning:** a promiscuous hub miRNA (large total mass) gets each of its edges
down-weighted, so it cannot dominate a gene merely by targeting everything.

*Code:* `mirna_mass_denominators`, `_edge_weights`.

---

## 3. Expression piece — building blocks

Three independent blocks are combined into the modes in S4.

### 3a. z-score (dynamics)

```
z[m,s] = (x[m,s] - mean_m) / sd_m          (= 0 if sd_m = 0)
```

- **Meaning:** deviation of the arm from *its own* typical level. Centered at 0,
  unitless, signed.
- **Caveat:** an arm sitting at its mean -> `z ~ 0` -> its contribution is
  **zeroed**. A stably-high arm contributes nothing under z. This is exactly why
  the *attribution* step uses a no-z mode (see "Attribution policy" below).

### 3b. softmax share (relative abundance among a gene's regulators)

Logit = abundance centered on the cohort median, then a gene-local softmax over the
regulators of `g`:

```
                 exp( x[m,s] - med_m )
sm[m,s] = ---------------------------------------       and   SUM over m in R(g) of sm[m,s] = 1
           SUM over m' in R(g) of exp( x[m',s] - med_m' )
```

- **Meaning:** fraction of gene `g`'s regulator "abundance budget" held by arm `m`
  in sample `s`. Always positive; sums to 1 across `R(g)`.

### 3c. absolute-abundance anchors

```
logrpm[m,s]   = max( x[m,s], 0 )
absratio[m,s] = max( x[m,s] / med_m, 0.25 )
```

- **Meaning:** keep absolute level in the score, so a relative share between two
  barely-expressed arms is not treated like one between two highly-expressed arms.

*Code:* `expression_multiplier` (z / softmax base), `_softmax_rows`,
`_abs_ratio_matrix`.

---

## 4. Expression modes (`expr_mult`)

Each mode is a product of S3 blocks. Selected by `config.PRESSURE_EXPR_MODE` (spine)
or passed explicitly (attribution / sensitivity grid). Formulas:

```
z                   expr_mult = z[m,s]                              sign: +/-   (sensitivity baseline)
softmax             expr_mult = sm[m,s]                             sign: +     (pure relative abundance)
softmax_logrpm      expr_mult = sm[m,s] * logrpm[m,s]               sign: +     (ATTRIBUTION / shares, no-z)
softmax_absratio    expr_mult = sm[m,s] * absratio[m,s]             sign: +     (abundance-ratio variant)
z_softmax           expr_mult = softmax( z[m,s] )                   sign: +     (variant)
softmax_z           expr_mult = sm[m,s] * z[m,s]                    sign: +/-   (variant)
softmax_z_logrpm    expr_mult = sm[m,s] * z[m,s] * logrpm[m,s]      sign: +/-   (CANONICAL SPINE)
softmax_z_absratio  expr_mult = sm[m,s] * z[m,s] * absratio[m,s]    sign: +/-   (variant)
softmax_z_blend     expr_mult = (0.5*sm[m,s] + 0.5) * z[m,s]        sign: +/-   (partially-pooled share)
devhealthy          expr_mult = dev[m,s]                            sign: +/-   (healthy-anchored, raw)
softmax_devhealthy  expr_mult = sm[m,s] * dev[m,s]                  sign: +/-   (healthy-anchored share)
softmax_devhealthy_logrpm  expr_mult = sm[m,s] * dev[m,s] * logrpm[m,s]  sign: +/-  (HEALTHY-ANCHORED SPINE)
```

where the **healthy-anchored deviation** replaces the within-tumor z:

```
dev[m,s] = x[m,s] - median_healthy(m)          # log2 units, NOT divided by an SD
```

- `median_healthy(m)` = **GTEx breast** arm median (346 samples, correct tissue),
  **quantile-normalised onto the TCGA scale** (each GTEx arm gets the TCGA-tumor
  arm-median value at the rank it occupies in the GTEx distribution → rank-relative
  acquired elevation, cross-platform-safe). Joined by **MIMAT accession** (canonical;
  416 vs 395 arms vs name matching). Source: `healthy_anchor.gtex_qn_baseline()`.
  Gap arms (GTEx zeros/omits, via multi-mapping collapse of canonical members let-7a,
  miR-30a, miR-10a, miR-16, miR-17, miR-200b…) are filled by two layers, then floor0:
  - **Layer 1 (PRIMARY) — anchor-calibrated miTED rank-transfer** (`source="mited_anchor"`):
    arms measured in **both** GTEx and DIANA-miTED are anchors that define a monotonic
    miTED→GTEx map (quantile-quantile via `np.interp` on sorted anchor pairs); a gap arm's
    GTEx-equivalent level is predicted from its miTED value. Uses miTED only for **relative
    position among reference miRNAs** (platform-robust), sidestepping the compositional-RPM
    problem (raw miTED-vs-TCGA slope 0.74, not 1). **Arm-specific and validated** — leave-
    one-out against real GTEx values: Pearson 0.90, median error 0.86 log2. Preferred over
    seed-family, which assigns a degenerate **family-max to every collapsed member** (e.g.
    let-7a *and* let-7c → let-7b's level; miR-103a → its lone low seed-mate miR-107). Fills
    ~52 arms (incl. orphans like miR-182). Blood/platelet markers (miR-451a, miR-486,
    miR-92a, miR-144…) excluded from the fit and targets (miTED breast n=12 blood-contaminated).
  - **Layer 2 (FALLBACK) — GTEx-internal seed-family** (`source="gtex_family"`): only for
    gap arms absent from miTED but with a GTEx-reliable seed-mate → **max** seed-mate
    (collapsed = dominant member). Currently 0 arms.
  - Combined: tumor-relevant coverage rises 182→213 / 222 (96%) at ≥15 RPM; 87/88 at
    ≥255 RPM (see `expression_threshold_coverage.tsv`).
  - **NAT is rejected** (field-effect: sits tumour-ward → cancels acquired signal).
- Arms with **no GTEx, no seed-mate, and not in miTED** (or blood) → baseline 0 → `dev = x`
  (full = tumour-acquired; correct for genuine oncomiRs, conservative for blood). Flags in
  `gtex_qn_baseline.tsv`: `imputed` (seed-mate or miTED-anchor proxy → lower precision),
  `floor0_tumor_expressed`.

**Why three headliners:**

- `softmax_z_logrpm` (spine) = *relative share* x *within-tumor dynamics* x *absolute
  level*. `z` centres on the **tumor-cohort mean**, so it answers "who moves across
  tumours" (coupling, S7).
- `softmax_logrpm` (attribution) = *relative share* x *absolute level*, **no z**.
  Always `>= 0`: nothing zeroed for being stably expressed, `abs_share == signed_share`.
- `softmax_devhealthy_logrpm` (healthy-anchored) = *relative share* x *deviation from
  healthy* x *absolute level*. Centres on the **healthy median** instead of the tumor
  mean, so it answers "who carries pressure the tumour *acquired* over healthy
  tissue". No SD division -> a uniformly tumour-elevated arm is **not** gated to ~0
  the way z would gate it. Constitutively-high-in-healthy arms (miR-125b, miR-24,
  miR-145) shed share; tumour-acquired oncomiRs (miR-21, miR-182, miR-200c) gain it.

*Code:* `_expr_contribution`, `_dev_health_matrix`; baseline `healthy_anchor.py`.

---

## 5. Share metrics (per gene, built from contributions `c`)

Two axes of choice: **global** (ratio of means) vs **mean** (mean of ratios), and
**abs** (magnitude) vs **signed/positive** (direction-aware).

```
                       mean_s |c(m,:)|
global_abs_share(m) = ----------------------------------          <-- PRIMARY RANKING
                       SUM over m' of mean_s |c(m',:)|

                        mean_s  c(m,:)
global_signed_share(m)= ----------------------------------        <-- COHERENCE DIAGNOSTIC
                        SUM over m' of mean_s  c(m',:)

mean_abs_share(m)   = mean_s [  |c(m,s)|        / SUM_m' |c(m',s)|       ]

mean_pos_share(m)   = mean_s [  max(c(m,s),0)   / SUM_m' max(c(m',s),0)  ]
```

where `mean_s` averages over samples.

| Metric | Reads as | Notes |
|---|---|---|
| `global_abs_share` | share of the gene's total pressure **magnitude** | stable single ratio -> **rank on this** |
| `global_signed_share` | share of the gene's **net** pressure | can be `<0` or `>1` when arms cancel -> **diagnostic only** |
| `mean_abs_share` | per-sample magnitude share, averaged | sensitive when a sample's total ~ 0 |
| `mean_pos_share` | share of only the **repressive** signal | arms net-activating in a sample add 0 to that sample |

**Under a no-z mode** (`softmax_logrpm`) every `c >= 0`, so
`global_abs_share == global_signed_share` and `mean_abs == mean_pos`. That collapse
is the whole point of the attribution policy.

*Code:* `compute_gene_pressure_contributions`.

### 5a. Structural share — `gene_struct_share` (de-abundance-double-counted resolution)

The attribution contribution `c = edge_w · sm · logrpm` (mode `softmax_logrpm`) spends the
abundance variable `x[m]` **twice**: once exponentially inside the softmax share `sm`
(relative abundance among `R(g)`) and again linearly in `logrpm` (absolute level). Empirically
the within-gene `gene_share_tumor` ranks regulators ~2.4× more by abundance (within-gene
Spearman ≈ 0.60 vs arm mass) than by the curated evidence it is meant to weight (≈ 0.25 vs
`evidence_score`), because `edge_w` is near-flat within a gene (max/min ≈ 1.8×). So
`gene_share_tumor` answers *"who is most abundant on `g`"*, not *"who is `g`'s preferred
regulator."*

The **structural share** removes the second abundance factor — mode `softmax`
(`config.PRESSURE_STRUCT_EXPR_MODE`), `c_struct = edge_w · sm` — so relative abundance is spent
**once** (bounded, sums to 1 across `R(g)`) and `edge_w` regains proportional leverage:

```
                          mean_s |c_struct(m, g, :)|
gene_struct_share(m,g) = ----------------------------------------   (within-gene, mode = softmax)
                          SUM over m' in R(g) of mean_s |c_struct(m', g, :)|
```

This is the **"which arm" resolution axis**; the abundance-driven `gene_share_tumor` is retained
as a separate magnitude axis (and `mass_X` / cross-state deltas stay on `softmax_logrpm` for
cross-gene magnitude comparability — the structural mode is *not* cross-state comparable). The
coupling spine `softmax_z_logrpm` is unaffected. Columns: `gene_struct_share_{gtex,nat,tumor}`,
`gene_struct_share_rank_tumor`, `gene_struct_share_delta_tumor_gtex` (edge table); audit of how
often the dominant-regulator call flips vs the abundance axis →
`struct_vs_abundance_resolution.tsv`. On the **full cohort** (5,108 edges, 758 multi-regulator
genes, 2026-06-27 regen) the #1 regulator flips for **61%** of genes and the within-gene
abundance↔structural rank correlation is only **ρ ≈ 0.17** (median) — the two axes are nearly
orthogonal in crowded genes, confirming `logrpm` was dominating the ordering. Examples:
ABCA1 abundance→miR-93-5p but structural→**miR-33b-5p** (canonical validated regulator),
CDKN1A/p21 miR-10b-5p→miR-345-5p, CD274/PD-L1 let-7d-5p→miR-380-5p, EGFR miR-7-1-3p→miR-542-5p.
The coupling-driven `joint_edge_class` counts are **unchanged** by the split (spine untouched).

> **CAVEAT — share ⊥ coupling (do not read share as repression).** The structural split is a
> correctness fix on the *attribution/identity* axis **only**. Tested against realized
> anti-correlation at every level (both ladders, `struct_vs_abundance_coupling_assessment.tsv`):
> share — of *either* flavour — is essentially **orthogonal to realized coupling** (within-gene
> Spearman of share vs −`adj_rho` ≈ **0**: abundance +0.001, structural −0.021 on the abundance
> ladder; ±0 on the pressure ladder). The structural pick is **not** a better realized repressor:
> the #1-pick FDR-sig-repressor rate is *lower* for structural (34.7%→27.0% abundance ladder,
> 25.1%→19.9% pressure ladder), which is a **statistical-power artifact** — abundance favours
> higher-dynamic-range arms that more often clear FDR — not a truth gain (negative-coupling rate
> essentially unchanged, 60.4%→57.8%). On flip genes the structural pick is more repressive in
> only ~50% of cases (Wilcoxon p=0.90 / 0.30). **Therefore: a high `gene_struct_share` is evidence
> of structural preference, NEVER of repression.** Any "driver" claim must still require
> FDR-significant `adj_rho` (the §9 quadrant enforces this: `driver_candidate` = coupled AND
> top-quartile structural share; high share alone → `structural_only`). To improve *realized*
> resolution, sharpen the coupling estimator (`adj_rho`), not the share; to make preferred-regulator
> identity move *potential* pressure, sharpen `edge_w` (interaction quality) — not the abundance term.

*Code:* `mirna_state_class._struct_share_table`, `_struct_vs_abundance_audit`,
`_struct_coupling_assessment`; engine mode `pressure_engine` `softmax`.

---

## 6. Specificity — "global vs specific target"

How concentrated an arm's budget is on one gene. Two flavors, because orphans have
no realized (miRTar) edge.

**Realized** (curated edges; the network table):

```
                  mean_s |c(m, g, :)|
spec(m, g) = ----------------------------------------
              SUM over g' in targets(m) of mean_s |c(m, g', :)|
```

**TargetScan** (orphans; sequence-only):

```
                     ts_weight(m, g)
ts_spec(m, g) = ----------------------------------
                 SUM over g' of ts_weight(m, g')
```

`spec` is the fraction of the arm's *total* target-mass on one gene; the reference is
the even-spread baseline `1/|targets(m)|` (well under 1% for an arm with dozens-hundreds
of targets), so the thresholds read as a **relative concentration**, not an absolute
majority:

- **High** (`>= 0.05`) -> the gene commands a disproportionately large, several-fold-
  above-baseline share of the arm's budget = **specific / focused** target -> cleaner
  single-gene validation readout. (NOT "most" of the budget — even a focused arm rarely
  puts a majority on one gene.)
- **Low** (`<= 0.01`) -> at/below the even-spread floor; the arm spreads thin over many
  targets = **promiscuous hub**.

*Code:* `mir301_family_network._family_specificity`, `_ts_specificity_table`.

---

## 7. Coupling — the anti-correlation ladder (is the pressure *observed*?)

Partial **Spearman** correlation `rho` between an arm `x[m,:]` (or a family
predictor, S8) and target-gene RNA `y[g,:]`, removing confounders by regressing
them out of both variables and correlating the residual ranks. Built as a **ladder**
of nested covariate sets:

```
raw_rho                 covariates: none
rho_CPE_HRD             covariates: CPE (purity), HRD
rho_CPE_HRD_CN          covariates: CPE, HRD, target-gene copy number
rho_e2f_g2m             covariates: CPE, HRD, E2F+G2M proliferation metagene
rho_e2f_g2m_CN          covariates: CPE, HRD, target CN, E2F+G2M           <-- proliferation-surviving
rho_mki67 (_CN)         covariates: ... with MKI67-only proxy
rho_ortho_noE2F_MYC(_CN) covariates: ... with G2M+spindle metagene STRIPPED of E2F/MYC genes
```

Three proliferation proxies (differing in E2F/MYC dependence) test whether
repression survives proliferation adjustment without over-correcting through a
shared E2F path.

**Multiple testing:** Benjamini-Hochberg within each `(predictor_type, predictor,
scope)` family. Significance flags:

```
neg_sig_cn_fdr        = (rho_CPE_HRD_CN  < 0) AND (q_CPE_HRD_CN  < 0.05)    <-- headline
neg_sig_prolif_cn_fdr = (rho_e2f_g2m_CN  < 0) AND (q_e2f_g2m_CN  < 0.05)    <-- proliferation-surviving
```

**Variability gate** (so a flat arm's null coupling is not over-read), on
`log2(RPM+1)`:

```
low_variability = (sd_x < 0.5) AND (IQR_x < 0.75)
```

*Code:* `_partial_ladder`, `_proliferation_proxies`, `_build_cov`,
`mir301_family_depth._add_fdr_and_flags`, `_arm_variability`.

---

## 8. Family-level aggregates

`F` = the seed family `{301a-3p, 301b-3p, 130a-3p, 130b-3p, 454-3p}`.

### 8a. Family predictors (one value per sample, fed into the S7 ladder)

```
family_mean_log2rpm[s] = mean over m in F of x[m,s]
family_max_log2rpm[s]  = max  over m in F of x[m,s]
family_sum_log2rpm[s]  = log2( 1 + SUM over m in F of (2^x[m,s] - 1) )    (pool in linear RPM, re-log)
family_mean_z[s]       = mean over m in F of z[m,s]
```

- `sum_log2rpm` is true pooled abundance; `mean_z` is pooled *dynamics*. Each is
  correlated against the target like a single predictor
  (`predictor_type = family_aggregate`).

### 8b. Family-as-unit pressure share (resolution table)

```
family_abs_share(g, scope) = SUM over m in (F intersect R(g)) of global_abs_share(m)
```

recomputed **within each scope's samples** (cohort + each PAM50). This is the number
behind statements like "the family carries 30% of NR3C2's pressure in Basal vs 10%
in LumA" — a **relative composition** of the gene's regulatory crowd, *not*
expression and *not* coupling.

*Code:* `mir301_family_depth._family_predictors`,
`mir301_family_network._family_pressure_by_subtype`.

---

## 9. Regulator-resolution edge classification

For a family arm on a focus gene, position it **among all** of `R(g)`:

```
struct_share_rank_pct = rank_of(structural global_abs_share within gene) / n_regulators   <-- quadrant keys on this
share_rank_pct        = rank_of(global_abs_share within gene) / n_regulators   (abundance-driven realized share; reported)
abund_rank_pct        = rank_of(median_log2rpm   within gene) / n_regulators
```

`struct_share_rank_pct` is computed under the **structural** mode (`softmax`, §5a) so a merely-
abundant arm is not auto-labelled the gene's structural driver. The quadrant's `high_share` test
keys off it (falling back to the realized `share_rank_pct` only if unavailable). **The `coupled`
leg is independent and indispensable: share ⊥ coupling (§5a caveat), so top-quartile structural
share alone never implies repression — it yields `structural_only`, and only `coupled AND` high
structural share yields `driver_candidate`.**

Quadrant `edge_class = base | spec_tag`:

```
coupled  AND  struct top-quartile        ->  driver_candidate   (preferred regulator AND observed)
coupled  AND  not struct top-quartile    ->  real_but_minor     (observed but not a structural driver)
NOT coupled  AND  struct top-quartile    ->  structural_only    (preferred regulator but no observed coupling)
otherwise                                ->  weak_unsupported   (neither)
```

where `coupled = cohort_neg_sig_cn_fdr OR cohort_neg_sig_prolif_cn_fdr`, and

```
spec_tag = specific     if spec >= 0.05
           promiscuous  if spec <= 0.01
           moderate     otherwise
```

*Code:* `mir301_family_network._interpret_edge`, `_regulator_resolution`.

---

## 10. Orphan validation-priority composite

Ranks `(arm, gene)` edges for wet-lab follow-up. Only edges with any
FDR-significant negative coupling (cohort or best subtype) qualify. Each component
is z-scored across the qualifying candidate edges (`z(.)` below = standardize across
candidates), then combined:

```
priority = [ z(coupling) + z(ts_weight) + z(ts_spec) + 0.5 * z(subtype_conc_pos) ] / 3.5
           * novelty_weight
```

Component definitions:

```
coupling          = -rho   (proliferation-surviving CN-adj if it survives,
                            else CN-adj, else best-subtype CN-adj)
ts_weight         = TargetScan |context++| sequence prior
ts_spec           = S6 TargetScan specificity (focused target -> cleaner readout)
subtype_conc_pos  = max( 0, rho_cohort - rho_best_subtype )   (extra negativity in a subtype)
novelty_weight    = 1.0 for ts_orphan / mirtar_low ;  0.6 for already-published
```

`recommended_context` = the FDR-significant, more-negative PAM50 subtype (else
`cohort`), mapped to field-standard PAM50 cell lines for a suggested model.

> **Caveat baked into the ranking:** miRTarBase-HE absence is a *curation gap*, not
> novelty. Many top "orphan" edges are luciferase-validated elsewhere; the list
> ranks **in-tumor / breast-context confirmation** candidates, and the
> `literature_note` column flags edges with known external validation.

*Code:* `mir301_family_network._validation_priority`.

---

## 11. State-trajectory classification + gain/loss (`mirna_state_class`)

Formalises *how a miRNA's repressive pressure moves* across the three tissue states
(GTEx-healthy -> NAT -> tumor) at **two resolutions**. Module:
`mirna_hallmark.mirna_state_class`. Cross-state comparability uses the **no-z**
`softmax_logrpm` mode (the z-spine is standardised *within* each state, so it is not
comparable across states); GTEx is cross-platform, so contrasts use **within-state
percentile rank**, not absolute deltas.

### 11a. Global pressure mass + the gain/loss deltas (miRNA resolution)

```
mass_X(m)  = Σ_{g in HallmarkUniverse} c(m, g, X)        # X in {gtex, nat, tumor}
share_X(m) = mass_X(m) / Σ_m' mass_X(m')
rank_X(m)  = within-state percentile rank of mass_X(m)
```

`c(m, g, X) = edge_w(m, g) * expr_mult_softmax_logrpm(m, X)` (§1, §4 with `ATTR_MODE`).
The three legs (rank deltas, range −1..+1):

```
dHN(m) = rank_nat   - rank_gtex     # healthy -> NAT   (field effect)
dNT(m) = rank_tumor - rank_nat      # NAT -> tumor     (tumor-specific increment)
dHT(m) = rank_tumor - rank_gtex     # healthy -> tumor (ACQUIRED axis = headline)
```

A **gainer** has `dHT > 0` (gained repressive pressure over healthy); a **loser** has
`dHT < 0`. The healthy anchor here is the **GTEx state**, never NAT. (This is the
clean replacement for the old `delta_acquired_vs_zspine`, which reallocated share
between two scoring schemes *within tumor* and was neither a healthy->tumor nor a
NAT->tumor shift.)

**Magnitude axis (`log2fc_tumor_vs_healthy`) — dHT's complement.** `dHT` is a within-state
*rank* delta, so it **saturates**: an arm already near the top of the healthy distribution
cannot gain rank even when its absolute abundance rises several-fold (the max achievable
`dHT` falls to ~0.08 in the 80-95% healthy-abundance bin). 90% of >4-fold-up arms are
rank-invisible. The magnitude axis recovers them on the **same QN'd TCGA scale** (so it is
cross-platform-safe, unlike a raw GTEx-vs-TCGA FC):

```
log2fc_tumor_vs_healthy(m) = tcga_tumor_median(m) − healthy_baseline_tcga(m)   # NaN for floor0
acquired_magnitude_class   = strong_gain (>2) / moderate_gain (>1) / flat / moderate_loss / strong_loss
acquired_gainer            = (dHT > NET_THR)  OR  (log2fc > MAG_STRONG=2)    # rank OR magnitude
acquired_axis              = rank+magnitude / rank_only / magnitude_only / none
```

The two axes answer different questions: `dHT` = "gained **share** among regulators"
(saturating, rank); `log2fc` = "rose in **absolute abundance**" (the oncomiR surge —
miR-182/183/375/93/141/210 are `magnitude_only`, dHT≈0 but +9 to +13 log2). Use
`acquired_gainer` for the union; 236/543 arms gain on at least one axis (29 rank, 187
magnitude, 20 both).

**Primary class** on `(s1=dHN, s2=dNT, net=dHT)` (`NET_THR=0.15`, `MOVE_THR=0.10`):

```
healthy_unknown        net is NaN (no GTEx rank) -> fall back to nat_axis (tumor vs NAT)
stable                 |net| < 0.15
non_monotonic          |s1|,|s2| >= 0.10 AND sign(s1) != sign(s2)
progressive_{gain|loss}  both legs move with net   (healthy<NAT<tumor or reverse)
tumor_acquired_{..}    only the NAT->tumor leg moves (NAT not yet shifted)
field_effect_{..}      only the healthy->NAT leg moves (gain already in NAT)
net_{gain|loss}        net passes but a leg is NaN
```

### 11b. Subtype-specificity secondary class (Tau over PAM50)

On the per-PAM50 **acquired** share (`ANCHOR_MODE = softmax_devhealthy_logrpm`, §4;
abundance share if `--healthy-anchor` is off), with subtype share vector `a` over the
4 PAM50 subtypes:

```
tau = Σ_i (1 - a_i / max(a)) / (k - 1)        # k = 4 ; 0 = uniform, 1 = one-subtype
```

```
subtype_specific:<argmax>   tau >= 0.50
subtype_excluded:<argmin>   present in >= 3 subtypes AND a_min < 0.25 * median(a)
pan                         tau < 0.30
intermediate                otherwise
```

### 11c. Edge resolution — movement of pressure on a gene by a miRNA (m -> g)

Because `edge_w(m, g)` is **static across states**, the per-edge pressure-movement
*direction* is inherited from the miRNA's abundance trajectory (only the magnitude
localises onto `g`). So the meaningful edge class is the **join** of two axes:

1. **Pressure-movement (potential):** per-state `c(m, g, X)`, the same-platform
   increment `edge_dNT_abs = c_tumor - c_nat`, and
   `edge_share_tumor = c_tumor(m,g) / mass_tumor(m)` (this edge's fraction of the
   miRNA's tumor budget). Direction = the miRNA's `pressure_move` (gain/loss/flat).
   Complement: the **gene-side** share `gene_share_{X}(m,g) = c_X(m,g) / Σ_{m'∈reg(g)} c_X(m',g)`
   (this arm's fraction of the GENE's total incoming miRNA pressure; `edge_share` is arm-side,
   `gene_share` is gene-side — they answer opposite questions), with `gene_share_rank_tumor`
   (1 = the gene's dominant regulator) and `gene_share_delta_tumor_gtex` (grip tightening/loosening).
   **`gene_share_*` is abundance-driven** (built from `softmax_logrpm` `c`, abundance spent twice);
   the de-double-counted **`gene_struct_share_{X}`** (§5a, mode `softmax`) is shipped alongside as the
   *preferred-regulator* axis, with `gene_struct_share_rank_tumor` and
   `gene_struct_share_delta_tumor_gtex`. Use `gene_struct_share_*` for "which arm owns the gene";
   keep `gene_share_*` for realized magnitude.
2. **Coupling-movement (realised):** composition+prolif-adjusted partial Spearman
   `adj_rho_{gtex,nat,tumor}` (§7 ladder), with **NAT and GTEx computed separately**
   (unlike S7's `cross_state_coupling`, which pools them into "normal"). Yields the
   per-edge `coupling_archetype` (tumor_acquired_repressor, lost_repressor, …).

```
nat_decoupled : adj_rho_nat is FDR-sig (|.|>=RHO_FLOOR, q<Q_ALPHA) AND its sign differs
                from BOTH adj_rho_tumor and adj_rho_gtex   (the miR-9 +0.43-in-NAT case)

joint_edge_class:
  uncoupled            no partial-rho scored (n<MIN_N or capped out)
  nat_decoupled        NAT sign-flip vs both healthy(GTEx) and tumor
  non_monotonic        the miRNA's own trajectory is non-monotonic
  acquired_realized    arm acquired (rank OR magnitude) AND tumor-step realized repression
                       (`tumor_acquired_repressor` or `repression_replaces_coexpr`)
  field_established_realized
                       arm acquired AND `0RR` coupling (NAT+ tumor repressive, GTEx neutral)
  acquired_unrealized  arm acquired but repression flat/absent/positive in tumor
  constitutive         repressive in tumor AND in healthy, without arm acquired gate
  lost                 pressure loss   OR  lost_repressor
  stable               otherwise

Per-edge trajectory columns (strict `coupling_sig_pattern` R/P/0 vs directional
`coupling_dir_pattern` r/p/0 without FDR):
  nat_underpowered_repressive : NAT ρ<-RHO_FLOOR but NAT q>=Q_ALPHA and tumor sig R
  nat_tumor_repressive_concordant : both NAT and tumor ρ<-RHO_FLOOR (directional)
  edge_pressure_gain_gtex / edge_pressure_gain_nat : c_tumor > c_gtex / c_nat (sensitivity)

Arm acquired gate for joint classes uses `arm_acquired_gainer` when present
(`dHT>NET_THR` OR `log2fc_tumor_vs_healthy>MAG_STRONG`), else `pressure_move==gain`.
Default coupling covariates: `estimate_epi` (ESTIMATE immune+stromal + epithelial + prolif).
```

The per-miRNA `odd_nat_coupling` flag = "has any `nat_decoupled` edge."

### 11d. Special-case flags + secondary tables

```
odd_nat_abundance      sign(dHN)==sign(dHT) AND |dHT|>=0.10 AND |dHN|>=0.5*|dHT|
                       (NAT pressure already >= half-way to tumor: field-effect)
baseline_source        gtex | mited_anchor | gtex_family | floor0  (from healthy_anchor)
imputed_baseline       source in {mited_anchor, gtex_family}
blood                  arm in BLOOD_MIRNAS
low_healthy_confidence source != gtex OR blood  (the acquired/dev axis is less certain)
```

Outputs (under `output/tissue_reference/mirna_state_class/`):
`mirna_gene_edge_class.tsv` (edge), `mirna_gene_edge_class_by_pam50.tsv`
(per-edge tumor pressure + share + dominant subtype + edge-level Tau),
`mirna_state_class.tsv` (per miRNA), `mirna_state_class_by_hallmark.tsv`
(per-(miRNA, Hallmark set) mass/rank/`set_class`).

### 11e. Gene-centric acquired pressure + the floor0 artifact split

Per gene, total tumor pressure and the acquired gain over each healthy reference:

```
pressure_X(g)       = Σ_m c(m, g, X)                 # X in {tumor, nat, gtex}
acquired_vs_gtex(g) = pressure_tumor − pressure_gtex   # principal healthy anchor
acquired_vs_nat(g)  = pressure_tumor − pressure_nat    # collapse-free cross-check
```

**GTEx collapse-repair (DEFAULT, `impute_gtex=True`).** GTEx's uniquely-mappable pipeline
zeros canonical arms (let-7a, miR-30a, miR-182…) via multi-mapping collapse. Before any
contribution is computed, `_impute_gtex_state` repairs the GTEx **abundance** matrix: arms
that are collapsed there (state-matrix median <= `ABUND_FLOOR`) but carry a healthy reference
in the baseline (`source` in {gtex, mited_anchor, gtex_family}) get their `gtex_median`
(same log2(TPM+1) scale) broadcast across GTEx columns. So `c_gtex`, `mass_gtex`, `dHT` and
`acquired_vs_gtex` are **collapse-free at source** (278 arms repaired; this collapses the
let-7a/miR-30a/miR-182 artifact gainers from dHT≈+0.96 to ≈0). **Coupling is left on the raw
GTEx matrix** (it re-derives its own bundles), so partial-rho is unaffected. `--no-impute-gtex`
reverts to raw collapsed GTEx (audit only).

The split below is therefore now a **residual diagnostic** (post-repair `frac_gain_imputable_artifact ≈ 0`,
`acquired_vs_gtex_corrected ≈ acquired_vs_gtex`); it stays for transparency and for the
`--no-impute-gtex` audit path. The floor0 gain (`c_gtex<=0`) is split by the arm's
**healthy-anchor `baseline_source`**:

```
pos_gain(g)        = Σ_{m: gain>0} gain_gtex(m,g)             # gain_gtex = c_tumor − c_gtex
imputable_artifact : c_gtex<=0 AND source in {gtex, mited_anchor, gtex_family} AND not blood
true_acquired      : c_gtex<=0 AND source == floor0 AND not blood
blood              : c_gtex<=0 AND arm in BLOOD_MIRNAS

frac_gain_imputable_artifact = Σ imputable gain / pos_gain    # FALSE gain (correctable)
frac_gain_true_acquired      = Σ true gain      / pos_gain    # genuine acquisition
frac_gain_blood              = Σ blood gain     / pos_gain
acquired_vs_gtex_corrected   = acquired_vs_gtex − Σ imputable gain
```

`acquired_vs_gtex_corrected` treats imputable-artifact arms as healthy≈tumor (their
true gain ≈ 0), removing the collapse inflation while keeping genuinely-acquired arms.
Output: `gene_acquired_pressure.tsv` (sorted by `acquired_vs_gtex`).

### 11f. Card `shift_class` — the CALIBRATED two-axis progression annotation (MH-166)

The per-edge progression class emitted by `learned/card.py::_shift_class` (flows into
`canonical_card.tsv` → `progression_edge_card`). This is the **learned successor** to the §11c
heuristic `joint_edge_class` (which routes through the §6b-RETIRED pressure spine — do not read
the two as the same object). It is built from **two orthogonal axes**, deliberately separating
*when repression engages* from *when dose rises* (the old label thresholded coupling ONLY, at a
−0.1 cut on a null 3-4× too narrow — MH-123):

**Axis A — CALIBRATED coupling, per state.** For each state `X ∈ {hly, nat, tum}` the arm-level
partial-Spearman `coupling_X` (§7, learned-model residuals) is scored against a **per-state
site-free empirical null** (`eval/site_free_null.py::fit_state(X)`, Efron per-abundance-quintile
on pairs that *cannot* repress — no 3'UTR site), yielding `coupling_p_X` and `coupling_z_tum`:

```
(coupling_p_X, coupling_z_X) = NullLadder(X).pvalues([coupling_X], [abund_arm_X])
rep(X) := (coupling_p_X == coupling_p_X) and coupling_p_X < COUPLING_ALPHA   # COUPLING_ALPHA = 0.05
```

The null reproduces MH-123: tumour inflation **2.67×** (sd0 ≈ 0.08 vs theoretical 0.031), NAT 1.25×,
GTEx ~1.5–2× (lowest GTEx abundance bin unfit → sd0=1.0 guard). ⚠ set-level, NOT per-edge FDR.

**Axis B — SAME-PLATFORM dose (SOE).** The trusted paired NAT→tumour arm log-fold-change
`arm_lfc_NAT_TUM` (never the soft cross-platform QN healthy leg):

```
dose_gain := arm_lfc_NAT_TUM > 0.3     # dose_class = gain / loss (<-0.3) / flat
```

**Composite class** (`h,n,t = rep(hly),rep(nat),rep(tum)`; evaluated top-down):

```
# --- NOT calibrated-repressor in tumour (not t) ---
undetectable              spiker OR arm_pct_floor < 30
dose_acquired_uncoupled   dose_gain and not (h or n)             (⭐ NEW: dose up, NOT calibrated-coupled)
lost                      (h or n)                               (was coupled healthy/NAT, gone in tumour)
uncoupled                 otherwise                              (never coupled, no dose gain)
# --- calibrated-repressor in tumour (t) ---
composition_explained     retention < 0.4                        (coupling is composition-driven)
constitutive              h and n                                (calibrated-coupled healthy + NAT + tumour)
field_established_realized (not h) and n                         (coupled from NAT onward — field effect)
acquired_realized         (not h) and (not n) and dose_gain      (tumour-only coupling + dose acquired)
tumour_realized           (not h) and (not n) and not dose_gain  (tumour-only coupling, dose flat)
nat_decoupled             h and not n                            (coupled healthy, dropped by NAT, back in tumour)
```

Census (genome-wide `progression_edge_card`, MH-166): `undetectable` 1145 · `uncoupled` 1104 ·
`dose_acquired_uncoupled` **571** · `acquired_realized` 484 · `tumour_realized` 183 · `lost` 162 ·
`composition_explained` 161 · `field_established_realized` 77 · `nat_decoupled` 52 · `constitutive` **6**
(collapsed from 127 under calibration).

**Quantified companions** (emitted alongside the class, per the axiom-5 continuum rule — do not
hard-bin as a headline count): `realization_score = −coupling_z_tum` (higher = stronger calibrated
tumour repression); `dose_class`; the abundance-vs-wiring decomposition `Δ(M·a) = term_ABUND
(M_NAT·Δa) + term_WIRING (a_NAT·ΔM) + term_INTERACT (Δa·ΔM)` with `wiring_frac = |term_WIRING|/Σ|terms|`
(SOFT — per-family ΔM is n-limited at NAT; cross-check aggregate retention, `state.cross_state_transfer`).

⭐ **Calibration separates realization at EQUAL dose** (MH-166): `dose_acquired_uncoupled` realizes
+0.002 (≈0) vs `acquired_realized` −0.101 at matched dose (1.13 vs 1.19) — the exact case the old
−0.1 coupling-only label buried in `constitutive`/`tumour_realized`. MH-160 (class→within-patient
realization) RE-VALIDATED on the new class: core live-vs-acquired ordering holds/strengthens (−0.046
gene-clustered) but the decoy-arbitrated increment is now marginal (−0.026, p=0.098) as `constitutive`
collapsed 127→6 under calibration. Registry: **MH-166**.

### 11g. Healthy-leg provenance + the graded `healthy_uninformative` flag (MH-166 follow-up)

The healthy (GTEx) coupling leg `coupling_hly`/`coupling_p_hly` is **mechanically blind** where the GTEx
uniquely-mappable miRNA pipeline COLLAPSED the arm (multi-mapping paralogues → per-sample GTEx ≈ 0 → ρ≈0 by
construction, not biology). `card.py::_healthy_leg_map` cross-references the DIANA-miTED healthy-breast cohort
(`healthy_anchor`, n=12; positional rank stable — median rank_sd 2.8, the flag-driving oncomiRs rank_sd<1,
detected 12/12) to assign per-arm provenance and a graded flag:

```
healthy_leg  ∈  measured        GTEx arm above ABUND_FLOOR and <50% floored → coupling_hly trustworthy
                collapse_blind   GTEx floored BUT miTED shows expressed → coupling_hly BLIND (not absent)
                true_absent      floored in GTEx AND miTED → genuinely absent in healthy
healthy_potential = HA.gtex_qn_baseline()          # imputed (miTED/GTEx) healthy level = repression-POTENTIAL proxy
healthy_uninformative = (healthy_leg == collapse_blind) AND (healthy_potential > _FLOOR=log2(11))   # GRADED
```

**GRADED, not binary.** Repression is abundance-gated (AGO/RISC saturates), so a LOW-abundance collapse_blind arm
has abundance-bounded potential ⇒ its 'acquired' call is safe even blind; blindness only THREATENS the call when the
arm is functionally present in healthy (`healthy_potential > RPM≥10 floor`). This cuts the flag **826→305** genome-wide
(acquired-class contamination 24%→13%) and MH-160 is robust to excluding it (−0.046→−0.047, p=0.001).

**Semantics of the flag (MH-166 follow-up, measured):** a flagged edge is **dose-constitutive** (miTED-high) with
**coupling history blind** — NOT "acquired." The blind healthy *coupling* is recoverable for the ~7–30 arms with a
**same-seed co-transcription instrument** (a mappable cluster-mate, e.g. miR-20a-5p for miR-17-5p, tumour-validated
corr 0.88): the instrument's GTEx coupling estimates the seed's healthy repression. Demonstrated: miR-17→{ERBB3,
NCOA3, JAK1, MEF2D, ARID4B, TLR7} surrogate-GTEx ≈ −0.09 vs tumour ≈ −0.41 ⇒ **dose-constitutive but coupling-ACQUIRED**,
per-target-graded (TLR7 fully de-novo −0.001, MEF2D partially pre-existing −0.168). The genome-wide surrogate node +
per-target GTEx→NAT→TUM identity trajectory is a SCOPED follow-up (not yet built). Registry: **MH-166**.

---

## 12. Gene-set architecture: topology propagation + tumor-direction prior (`geneset_architecture`)

Quantifies the **net effect of a miRNA's pressure on a program's output**, accounting
for (a) where each gene sits in the set's regulatory topology and (b) whether the
program is pro- or anti-tumor in breast cancer. Network = OmniPath directed+signed,
induced on the set's members.

### 12a. Signed-influence propagation (finite-horizon Katz) + calibrated edge sign

**Calibrated continuous edge sign (replaces the binary flatten).** Each ordered gene pair
`(s→t)` aggregates *all* OmniPath records into a graded sign instead of a hard `±1/0`:

```
net_sign(s→t) = (n_stimulation − n_inhibition) / (n_stimulation + n_inhibition)   ∈ [−1, +1]
```

Pure activator → `+1`, pure inhibitor → `−1`, dual / conflicting (the *ambiguous* case that
used to flatten to 0) → a graded value. Per node, `effect_sign_soft = mean net_sign over
out-edges` (the binary `effect_sign` is kept for the legacy `arch_signed_pressure`).

Build signed adjacency `S` on the induced edges using `net_sign` (`S[i,j] = net_sign(i→j)`;
the legacy path used `+1/−1/0`). Net downstream influence of
*increasing* node `g` on the rest of the program:

```
T = Σ_{k=1..H} (α · S)^k                  # α = 0.25, H = 6
output_influence(g) = 1 + Σ_j T[g, j]      # self + net signed effect on all other nodes
```

A central activator scores `>0`; an inhibitor hub (represses many downstream) scores
`<0`. A miRNA represses `g` (magnitude `c_tumor(m,g) > 0`), so the program-output change
flips sign with the node's influence:

```
program_output_change(m) = Σ_g [ −c_tumor(m, g) · output_influence(g) ]
```

`>0` = the arm's repression nets to *activating* the program (it mostly hits inhibitor
hubs); `<0` = it nets to *repressing/damaging* the program.

### 12b. Tumor-direction prior -> `pro_tumor_pressure`

A coarse BC-context literature polarity per Hallmark set (`TUMOR_DIRECTION_PRIOR`):
`+1` pro-tumorigenic (EMT, HYPOXIA, GLYCOLYSIS, MTORC1, PI3K/AKT, E2F, G2M, MYC v1/v2,
MITOTIC_SPINDLE, ANGIOGENESIS, TGFβ, WNT, NOTCH, KRAS_UP, UPR, IL6/STAT3); `−1`
tumor-suppressive / anti-tumor immune (APOPTOSIS, P53, IFNα, IFNγ, INFLAMMATORY,
ALLOGRAFT_REJECTION); `0` ambiguous/metabolic (excluded from the call, **not** from the
topology score).

```
pro_tumor_pressure(m, set) = prior(set) · program_output_change(m)
```

`>0` = the arm's pressure pushes the program toward malignancy (activating a `+1`
program or **releasing a `−1` brake**); `<0` = toward suppression (damaging a `+1`
engine or de-repressing a `−1` program). Because miRNAs repress, the dominant pro-tumor
mechanism is **brake-release on `−1` programs** (miR-21-5p -> P53/APOPTOSIS), not
amplification of `+1` engines.

**Acquired-weighted parallel.** The same score weighted by the per-edge **acquired pressure**
`max(c_tumor − c_gtex, 0)` (collapse-repaired) instead of total tumor abundance isolates the
pro-tumor pressure that is **tumor-acquired** vs constitutive:

```
pro_tumor_acquired(m, set) = prior(set) · Σ_g [ −max(c_tumor − c_gtex, 0) · output_influence(g) ]
```

This runs alongside `pro_tumor_pressure` (not replacing it). Diagnostic finding: miR-375 and
miR-17-3p top *both* (genuinely acquired pro-tumor), while let-7d/let-7f/miR-125a/miR-130a
lead the abundance axis but drop out of the acquired axis (their pro-tumor pressure is
constitutive — they are abundant in healthy too).

Cross-set roll-up (`architecture_mirna_cross_set.tsv`, polarity sets only):
`n_sets_pro_tumor`, `n_sets_anti_tumor`, `sum_pro_tumor_pressure` (headline rank),
`top_pro_tumor_set`; parallel `sum_pro_tumor_acquired` / `n_sets_pro_tumor_acquired` /
`top_pro_tumor_acquired_set`. Also retains the topology-only `arch_signed_pressure`
(double-negative inhibitor-release) and PageRank/redundancy `w_arch` from the pilot.

### 12c. Family + total-layer rollups and brake-release annotation

- **Total-layer per set** (`architecture_all_set_summary.tsv`):
  `total_program_output_change = Σ_m program_output_change(m)` (sign = `net_program_direction`:
  `damaged` if `<0`); `total_pro_tumor_pressure`, `total_pro_tumor_acquired`. For **ambiguous
  (prior=0) sets** this total is the only directional readout kept (no malignancy sign); all 27
  are `damaged` (the miRNA layer net-reduces every program output).
- **Family rollup** (`architecture_family_set_pressure.tsv` long, `architecture_family_cross_set.tsv`
  wide): group arms by a **seed-approx family key** `_family_key` — arm name with the paralog
  letter + copy index dropped, numeric stem + 5p/3p kept (miR-30a/b/…-5p → `miR-30-5p`). A
  **name-stem proxy**, not the exact TargetScan seed (miR-200a vs 200b/c differ by one seed nt).
  Sum each family's per-edge pressures within every set, then roll up the pro-tumor breadth.
- **Brake-release annotation** (Q6, merged into `architecture_mirna_cross_set.tsv`): each arm's
  **own** Stage-1 trajectory is joined as `own_dHT` / `own_log2fc` / `own_pressure_move` /
  `own_acquired_axis`, so an anti-tumor (suppressive) arm can be read together with whether the
  brake itself is lost in tumor. Finding: the top anti-tumor TS-miRNAs are rank-flat but
  magnitude-**up** (not lost) — failure is downstream (realized/coupling), not arm loss.

### 12d. Gene-role malignancy score + per-(set,gene) hit rollup

**Gene-role malignancy score** (per set,miRNA; `malignancy_sign∈{−1 TSG,0,+1 onco}` from `gene_roles`):

```
mal_pro_tumor(m)      = Σ_g [ −malignancy_sign(g) · c_tumor(m,g) ]            # graph-free driver call
mal_pro_tumor_arch(m) = Σ_g [ −malignancy_sign(g) · c_tumor(m,g) · w_arch(g) ]  # × centrality MAGNITUDE only
```

`w_arch` is always ≥0 (PageRank/redundancy), so `mal_pro_tumor_arch` is **not** multiplied by signed
`output_influence` (an earlier version was — it double-counted direction and let one hub edge dominate).
This recovers miR-21-5p (topology `pro_tumor` +0.4 → `mal_pro_tumor` +34.5 #1) and the oncomiR roster.

**Per-set malignancy aggregate** (`architecture_all_set_summary.tsv`, the prior-FREE companion to
`total_pro_tumor_pressure = prior·total_program_output_change`):

```
total_mal_pro_tumor(set) = Σ_m mal_pro_tumor(m, set) = Σ_{edges in set} −malignancy_sign(g)·c_tumor(m,g)
malignancy_direction      = pro_tumor if total_mal_pro_tumor > 0 else anti_tumor
```

with `total_tsg_credit` / `total_onco_collateral` / `total_mal_pro_tumor_acquired`. Unlike the prior
total (one hand-assigned set sign, `ambig` for prior=0), this is defined for **all 50 sets** from the
curated roles and **resolves the 27 ambiguous sets** (TNFA_NFKB +20.7, DNA_REPAIR +1.6 → pro_tumor). It
validates against the prior (prior=−1 suppressive sets → malignancy-pro_tumor in 5/6) and ranks
P53_PATHWAY #1 (+53.7). **High-specificity / low-coverage**: only the ~232 curated drivers contribute, so
a small total means "few driver genes hit" (ALLOGRAFT +439 by prior, +6 by malignancy) — the two
aggregates are complementary, not interchangeable.

**Continuous, all-gene role weight** (the graded counterpart to the binary driver-only `malignancy_sign`;
`gene_roles.load_gene_dependency`, source DepMap 26Q1 CRISPR Chronos gene-effect over the 96-line breast
panel, builder `_build_depmap_dependency.py`):

```
dep_role_weight(g)        = gene_effect_breast(g)        # oriented like −malignancy_sign
                                                          #  <0 = essential/dependency (repress = anti-tumor)
                                                          #  ≳0 = TSG-like (repress = pro-tumor)
mal_pro_tumor_cont(m,set) = Σ_g dep_role_weight(g)·c_tumor(m,g)     # ALL genes, not just drivers
total_mal_pro_tumor_cont(set) = Σ_m mal_pro_tumor_cont(m)          # + _acquired (weight max(c_tumor−c_gtex,0))
```

Defined for ~18.5k genes ⇒ 530/543 arms scored (vs driver-only binary). **One-sided**: CRISPR resolves
the dependency (negative/anti-tumor) tail well but the TSG (positive/pro-tumor) tail weakly, so this is
the **all-gene graded onco-collateral / anti-tumor axis** (the curated `tsg_credit` stays the better
pro-tumor detector; corr cont↔binary ≈0.29). Resolves the binary coverage blind spot: miR-24-3p is +4.7
binary but **−112.6 continuous (#1 anti-tumor)** — it represses many non-driver breast dependency genes;
top dependency-repressed sets are the proliferation engines (G2M −362, E2F −330, MYC −190).

**Within-set role: topology + curated TF hierarchy.** The gene's *functional role inside the program*
(master regulator / effector / inhibitor) is the **topology**, distinct from its oncogenic role:
`w_arch` (reverse-PageRank), `output_influence` (signed-path propagation, §12a) and `effect_sign_soft`
are per-set, continuous, derived from the OmniPath-induced subgraph (Hallmark itself = membership only).
That misses curated TF identity (10% of `w_arch` hubs are TFs; ~48% of curated TFs are below median
`w_arch`), so a curated TF census is added (Lambert 2018, `gene_roles.load_tf_census`):

```
is_tf(g)   ∈ {True,False}        # curated human TF (1,639)
hub_tf(g)  = is_hub(g) ∧ is_tf(g)   # curated TF ∧ topology hub = high-confidence master regulator
w_hier(g)  = w_arch(g) · (1 + TF_HIER_BOOST·is_tf(g))   # MANUAL tunable boost (default 0.5); w_arch untouched
mal_pro_tumor_hier(m,set)        = Σ_g −malignancy_sign(g)·c_tumor(m,g)·w_hier(g)   # TF-up-weighted role score (+acquired)
program_output_change_hier(m)    = Σ_g −c_tumor(m,g)·output_influence(g)·(w_hier/w_arch)(g)  # TF-up-weighted topology damage
tf_pressure(m,set)               = Σ_{g:is_tf} c_tumor(m,g)        # raw pressure on curated TFs; n_tf_targets = count
```

These are **scored** at per-(set,miRNA), cross-set (`sum_mal_pro_tumor_hier`, `sum_tf_pressure`,
`n_tf_targets`) and per-set (`total_mal_pro_tumor_hier`, `total_program_output_change_hier`,
`total_tf_pressure`) levels. `mal_pro_tumor_hier` ≈ `mal_pro_tumor_arch` where no TF target is hit
(333/543 arms identical, corr 0.99) and exceeds it for the 95 TF-hitting arms (miR-10b/150/25/106b).
`hub_tf` isolates the transcriptional master regulators hit by miRNA pressure (P53→TP53, G2M→E2F1/2/3/
MYC/MYBL2, IFNγ→STAT1) — sharper than topology-hub alone (which mixes in signaling hubs RAC1/CASP3).
At set level the TF-up-weighted damage re-ranks engines: TNFA_NFKB −513→−589 (TF pressure 137, highest),
G2M −646→−701, P53 −455→−510.

**Per-(set,gene) hit rollup** (`set_gene_hits` → `architecture_set_gene_hits.tsv`): localizes each set's
pressure to genes, reading both lenses at gene resolution:

```
incoming_pressure_tumor(g) = Σ_{m∈reg(g)} c_tumor(m,g)             # how much miRNA pressure lands on g
prog_change_gene(g)        = incoming_pressure_tumor(g) · output_influence(g)   # topology/hub damage routed via g
mal_hit(g)                 = −malignancy_sign(g) · incoming_pressure_tumor(g)   # gene-role view (TSG pressure = pro-tumor)
is_hub(g)                  = output_influence(g) ≥ Q90  OR  w_arch(g) ≥ Q90   # per-set top-decile centrality
```

`incoming_pressure_tumor` is gene-level (set-agnostic); only `output_influence`/`w_arch`/`malignancy_sign`
are set-specific. Hub hits are usually program effectors (repressing them damages the program —
let-7→AURKA/PLK1 in E2F), malignancy hits are the TSGs whose repression is pro-tumor regardless of
network position; the striking convergence is **TP53 itself = the dominant P53-program hub hit, by
miR-375**. Built by `--hits-only` from the existing weight + edge tables (no OmniPath re-fetch).

---

## 13. Stage-2 target-module composite + signed gene-pressure NMF (`mirna_comovement`)

### 13a. Target-module composite repression (realized module pressure)
For miRNA `m` with target set `T(m)` (within a state `s` with sample set `S_s`), each target gene
`g` is standardised **within state** across that state's samples:

```
z(g, j) = ( x(g, j) − mean_{j'∈S_s} x(g, j') ) / sd_{j'∈S_s} x(g, j')
```

The **per-sample composite** is the (optionally edge-weighted) mean over targets:

```
comp(m, j)   = (1/|T'|) · Σ_{g∈T'} z(g, j)                  (equal-weight, legacy)
comp_w(m, j) = ( Σ_g w(m,g)·z(g, j) ) / Σ_g w(m,g)          (edge-weighted; w = evidence_score)
```

where `T' = T(m) ∩ expressed`. Realized module pressure = composition+proliferation-adjusted
**partial Spearman** of the composite against the arm abundance:

```
rho_composite_{s}   = partial_spearman( comp(m,·),   abundance_m(·) | ESTIMATE )
rho_composite_w_{s} = partial_spearman( comp_w(m,·), abundance_m(·) | ESTIMATE )
```

A strongly negative ρ = high miRNA ↔ coordinately-low target module (the whole module is repressed,
not just one gene). `delta_tumor_gtex = rho_composite_tumor − rho_composite_gtex` = **acquired**
module repression. Empirically `rho_composite_w ≈ rho_composite` (r=0.97–0.99): edge weighting
washes out when averaging a module of near-even-evidence targets. **De-coupling control:** comparing
the arm's module composite (this §) against a single lost edge (§11c) separates *arm-wide* realized
failure (module composite collapses too) from *target-specific escape* (module stays coupled, one
edge breaks).

### 13b-role. Role-stratified composite (gene role, *not* a reweight)
**Design choice (the "should the composite weight by gene role/TF?" question).** The composite
above answers a single, deliberately role-*agnostic* question: *do the arm's targets co-move down
when the arm is high* (realized coordinated repression). Whether that repression is pro- or
anti-tumor is a **different** axis, owned by the malignancy/architecture layer (§12). Folding gene
role into one weighted composite would re-merge the two axes (the same conflation that inflated the
old `mal_pro_tumor_arch`; §12). So instead of reweighting we **stratify**: recompute the *same*
partial-Spearman composite separately over the arm's tumour-suppressor targets and its oncogene
targets (`malignancy_sign` from `gene_roles`, ≥4 role-annotated targets and ≥25 shared samples):

```
rho_composite_tsg_{s}  = partial_spearman( comp(m, T_TSG(m)),  abundance_m | ESTIMATE )
rho_composite_onco_{s} = partial_spearman( comp(m, T_ONC(m)),  abundance_m | ESTIMATE )
```

Interpretation: a negative `rho_composite_tsg_tumor` = the arm coordinately represses its
*suppressor* targets in tumour (anti-tumour realized pressure); negative `rho_composite_onco_tumor`
= coordinated repression landing on *oncogene* targets (pro-tumour collateral). The unstratified
`rho_composite` remains the headline; the two strata are an interpretive overlay, available only for
arms with enough annotated targets per stratum (most arms have few oncogene/TSG targets, so coverage
is partial by construction). TF status is captured indirectly through the architecture layer
(`output_influence`/`w_arch`), where master-regulator centrality already up-weights TF hubs.

### 13b. Signed (gain+loss) gene-pressure NMF
Per-(gene, tumor-sample) **dev** pressure `P_dev(g, j)` (mode `softmax_devhealthy_logrpm` vs the
healthy baseline) is signed (`>0` more pressure than healthy, `<0` relieved). NMF requires
non-negativity, so split each gene into two channels and stack:

```
M = [ P_pos ; P_neg ],   P_pos(g,j)=max(P_dev,0) (rows "g|gain"),   P_neg(g,j)=max(−P_dev,0) (rows "g|loss")
W, H = NMF(M, k)          # 2·|genes| × samples  →  W (gene-channel × factor), H (factor × sample)
```

A factor can load on a gene's gain channel, loss channel, or both, so **direction is preserved**
(unlike NMF on `|P_dev|`). Per factor: `gain_mass_frac = Σ_{g|gain} W / (Σ_{g|gain} W + Σ_{g|loss} W)`
(>0.5 ⇒ a gain-dominated program), top gain-genes / loss-genes, dominant Hallmark, and per-PAM50
activity (Kruskal-Wallis on `H`). Empirically all programs have `gain_mass_frac ∈ [0.97, 0.99]`
⇒ the tumor gene-pressure landscape is near-purely *acquired gain*; there is no coherent
pressure-relief program, so the gain-only acquired NMF (§ comovement (7a)) loses nothing.

### 13c. Share NMF (`_nmf_share`) — regulatory-share programs
The raw-`c_tumor` NMF (§ comovement (3)) is dominated by the abundance backbone. Factorising the
**share** matrices removes that scale and exposes complementary structure:

```
gene_share_tumor(m,g) = c_tumor(m,g) / Σ_{m'∈reg(g)} c_tumor(m',g)   # column-stochastic (per gene)
nmf_share_dominant      = NMF(gene_share_tumor matrix, k)            # dominant-regulator programs
d_share(m,g)            = gene_share_tumor(m,g) − gene_share_gtex(m,g)
nmf_share_grip_change   = NMF( [max(d_share,0) ; max(−d_share,0)] , k)  # ±-channel signed grip change
```

`nmf_share_dominant`'s top factor stays the backbone (loading corr ≈1 with the c_tumor NMF) but the
**secondaries are uncorrelated (≈0)** — they isolate low-abundance arms that *own* their few targets'
budgets (miR-125b/24/27~29 modules; miR-200 family → ETV5/ITGB4). `nmf_share_grip_change` direction is
set by `gain_mass_frac` (as §13b) and reads which arms co-GAIN vs co-LOSE regulatory grip across genes
independent of pressure inflation (the share analogue of the signed gene-pressure NMF). The **level**
share matrices are deterministic renormalizations of `c`, so only the *secondary* dominant-regulator
programs and the *delta* (grip-change) carry signal beyond the raw-`c` NMF — that is why only these two
share NMFs are shipped (not an edge_share/level twin).

---

## 14. De-coupling validation — confound-safe conditioning (`decoupling_validation`)

A "lost"/`nat_decoupled` edge is a miRNA→gene partial-ρ that is repressive in NAT/healthy and ≈0 in
tumour. Lowering the bar to raw expression would break the whole confounder model, so de-coupling is
validated by **conditioning**, not re-correlating. For a lost edge `(m, g)` the tumour partial
Spearman of target `g` vs arm `m` (ESTIMATE-adjusted) is re-fit with one extra covariate `v` at a
time and the **re-emergence** Δ is measured:

```
rho_base          = partial_spearman( y_g, x_m | ESTIMATE )
rho_plus_v        = partial_spearman( y_g, x_m | ESTIMATE, v )   v ∈ {target_cn, target_meth, target_atac, p_others}
{v}_reemergence   = rho_base − rho_plus_v
{v}_masked        = (reemergence ≥ 0.10)  AND  (rho_plus_v < −0.10)
```

A masked flag = the repressive coupling **re-appears** once the extra factor is held constant, i.e.
that factor (copy number, promoter methylation, chromatin accessibility, or co-regulator competition)
was confounding the loss.

- **`target_cn`** — ASCAT3 integer copy number of `g` (per participant).
- **`target_meth`** — per-sample promoter (body-fallback) β for `g`, built over the **378
  decoupled-target Hallmark genes** (`hallmark_gene_methylation.tsv.gz`, ~785 tumours), *not* the APM
  panel.
- **`target_atac`** — per-sample promoter ATAC signal for `g`: **direct** value where the participant
  was assayed (n≈74), otherwise the **PAM50 subtype-representative** mean of same-subtype assayed
  participants. `atac_frac_imputed` records the fraction of conditioned samples carrying the imputed
  (subtype-mean) value — for those, the covariate only absorbs *between-subtype* accessibility
  differences, so a null ATAC result on heavily-imputed genes is coverage-limited, not conclusive.
  Because subtype imputation dilutes per-tumour signal, a **strict direct-only** variant is also
  emitted (`atac_direct_*`): `rho_base` and `rho_plus_atac` are **both recomputed on only the ~70
  directly-assayed tumours** (drop every imputed sample), giving an un-imputed proof-of-concept test
  (`atac_direct_masked` with the same Δ≥0.10 / ρ<−0.10 rule). This is the **preferred** ATAC readout
  (the imputed `atac_masked` is the cohort-wide sensitivity view): it requires `n_atac_direct ≥ 25`.
- **`p_others` (competition control)** — the gene's **competitive field**
  `P_others(g,s) = Σ_{m'∈R(g), m'≠m} c(m',g,s)` (total OTHER incoming realized pressure, `softmax_logrpm`).
  This is the **only** conditioning term that makes the rank-coupling edge-specific: raw abundance is
  arm-global (identical across all of m's targets), and `edge_w`/`z`/`logRPM` are rank-invisible or
  rank-identical to abundance for a Spearman (`coupling_predictor_comparison`), so the softmax is just
  a nonlinear rank-1 proxy for exactly this field — here it enters **linearly, as a covariate** (the
  principled form). `competition_masked` (reemergence ≥ 0.10 ∧ ρ_plus < −0.10) = a co-regulator was
  *hiding* m's repression → **not a true loss** (m is the real, masked regulator; e.g. miR-1246→NFIB,
  −0.29 → −0.40). `competition_explained` (reemergence ≤ −0.10) = the residual coupling **was** the
  field, not m (miR-221-3p→DVL2, −0.07 → +0.04). Empirically the lost set is **robust to competition**
  (1 masked / 1 explained of 210, median Δ≈0) → the losses are genuine, not co-regulator takeover.

  **Performance evaluation** (`evaluate()` → `decoupling_validation_eval.tsv`; bootstrap CI + permutation
  null + stability, n=200 each): observed `competition_masked` **0 robust** vs **0.24 expected under the
  P_others-permutation null** → the masking is **well-calibrated, no over-flagging** (the single point-
  estimate flag does not survive bootstrap). The competitive field *does* carry **real but negligible**
  information: **50% of edges (105/210)** show a permutation-significant competition effect (median
  perm-p 0.04), but the **effect size is ~13× below the masking bar** (median |Δ|=0.008, 95%CI width
  0.038). The bootstrap CI + robustness is applied **uniformly to all four conditionings**
  (`_boot_condition`, CN/meth/ATAC/competition): **0 bootstrap-robust masks of *any* kind** — the
  point-estimate cnv/meth/atac/competition flags (1/2/1/1) *all* fail to survive resampling ⇒ the
  per-edge masking layer is **entirely point-estimate noise**. The **loss call itself is only ~80%
  stable** under patient resampling. Edge-level decoupling is intrinsically **low-SNR**; aggregating
  edge→arm-module→gene does **not** magnify the effect (strong-subset ρ_gtex ≈ −0.15 at every level).
  **GO-UP test (`go_up_stability` → `decoupling_go_up_stability.tsv`) — REFUTED:** going to the
  arm-module level *dilutes*, it does not denoise. Paired per-arm (same 54 arms): **best single edge
  ρ_gtex −0.209 vs arm-module −0.074 (3× weaker), dilution in 93% of arms; best-edge call stability 97%
  vs module 32%** — the module composite pools the arm's 1–few functional targets with its many
  non-functional ones. **The realized repression signal is concentrated at specific edges, not diffuse.**
  So the within-state repression call is *stable at the EDGE level* (95–97%); the only genuinely low-SNR
  piece is the cross-state **transition** (tumor-flat ~80%, tiny Δ, 0 robust masks) — a biology/data
  limit (the healthy→tumor coupling change is small), not a resolution problem. Productive direction:
  **functional-edge-selected** reads, not module aggregation.

**GTEx-primary healthy anchor (§14, `healthy_anchor_class`).** The "was repressive in healthy"
reference is **GTEx** (true healthy, consistent with the §11 acquired axis, which *rejects* NAT for
field-effect); NAT is a same-platform **confirmation** tier, not an independent trigger. Replaces the
permissive `g_neg OR n_neg` in `_edge_archetype_row` for *interpretation*: `gtex_repressive` /
`nat_repressive` (each ρ<−`RHO_FLOOR` ∧ q<`Q_ALPHA`) → `healthy_anchor_class ∈ {gtex_anchored,
nat_only, neither}`. On the canonical lost/nat_decoupled set: **123 gtex_anchored, 13 nat_only
(field-effect-suspect — admitted only by the OR), 74 neither** (pressure-trajectory losses, not
coupling-anchored — a diagnostic the OR hid).

Module-control (`decoupling_mechanism`) and competing-regulator context (`gene_regulator_share`,
`competitors_carry_coupling`) are computed alongside (see `ANALYSES_CATALOG`). The
`focal_high_regulator_share` gate behind `competitors_carry_coupling` keys off the **structural**
rank (`struct_share_rank_pct`, §5a) — so "the focal lost edge wasn't the gene's main regulator
anyway, the coupled competitors carry the brake" is a statement about the *preferred* regulator,
not merely the most abundant one (falls back to the abundance rank only if structural is absent).

### 14a. Gene-level realized net-repression (`gene_corepression`) + `escape_scope`

Arm-side `target_corepression` (§13a) asks "does the arm's whole *target module* move down." The
gene-side mirror asks "is the *gene* genuinely under net miRNA repression." For gene `g`, build its
**aggregate incoming pressure** per sample `Pagg(g, j) = Σ_{m∈reg(g)} c(m,g)·f(abundance_m(j))`
(`compute_gene_pressure`, mode `softmax_logrpm`) and take the ESTIMATE-adjusted partial Spearman of
the gene's **own** expression against it, per state:

```
rho_gene_pressure_{s} = partial_spearman( expr_g(·), Pagg(g, ·) | ESTIMATE )      ( <0 = net-repressed )
gene_repression_class : constitutive_repressed (ref<−.1 & tumor<−.1) | lost_repression (ref<−.1, tumor≈0)
                        | gained_repression | never_repressed         (ref = NAT else GTEx)
```

A single **direction-aware `gene_repression_scope`** (`edge_gene_repression_scope`) reads each edge's
`joint_edge_class` against the target's `gene_repression_class`, covering **both** directions:

```
lost / nat_decoupled:   edge_local_escape (gene net-repressed in tumour ⇒ co-regulators carry the brake)
                      | gene_wide_derepression (gene lost_repression ⇒ gene lost miRNA control)
                      | gene_not_repressed (never realizedly repressed in bulk)
acquired_realized/_unrealized:  gene_wide_acquired_repression (gene gained_repression ⇒ the gene is
                                  newly under net repression and this edge is part of it)
                              | reinforces_constitutive (gene already constitutive_repressed)
                              | edge_specific_acquired   (gene not broadly newly-repressed)
```

For `lost`/`nat_decoupled` edges the lost-side labels are surfaced as **`escape_scope`** in
`decoupling_validation`. This separates a coupling change that *matters for the gene* (gene-wide
de-repression, or gene-wide acquired repression) from one absorbed by co-regulators or never decisive —
symmetrically for both the loss and the gain side.

---

## Attribution policy (the "Q1" decision)

Three decisions for three distinct questions:

1. **Rank** miRNA contribution with `global_abs_share` (stable single ratio).
2. **Diagnose** stack coherence with `global_signed_share` (sign retained; a
   divergence from `global_abs_share` = arms cancelling).
3. **Compute** both under `softmax_logrpm` (no-z), so `abs == signed` and stably
   expressed arms are not zeroed.

The z-spine `softmax_z_logrpm` stays canonical for **coupling** (where dynamics are
the signal). Config switches: `PRESSURE_EXPR_MODE` (spine) vs
`PRESSURE_ATTRIBUTION_EXPR_MODE` (no-z attribution).

---

## Config quick map

| Constant (`config.py`) | Controls |
|---|---|
| `PRESSURE_EXPR_MODE = softmax_z_logrpm` | spine expression mode (coupling) |
| `PRESSURE_ATTRIBUTION_EXPR_MODE = softmax_logrpm` | no-z mode for shares / specificity |
| `PRESSURE_TARGET_NORM = evidence_mass` | per-arm promiscuity denominator `D(m)` |
| `PRESSURE_AGGREGATE = sum` | how regulators combine into pressure |
| `PRESSURE_EVIDENCE_SCORER = confidence_logclass` | base evidence scorer for edges (assay-directness × log1p × Functional MTI cross-counts; replaces raw CSV study-count at edge-load time) |
| `PRESSURE_ENCORI_ALPHA = 0.5` | ENCORI AGO-CLIP depth boost on (miRNA, gene) pairs present in the ENCORI parquet (0 = disabled) |
| `PRESSURE_MIN_EVIDENCE`, `PRESSURE_WEIGHT_MODE`, `PRESSURE_ABUNDANCE_FLOOR` | edge filtering / weighting |
| `PRESSURE_HYBRID_*` | M7 / M8 / M11 miRTar + TS hybrid modes |
