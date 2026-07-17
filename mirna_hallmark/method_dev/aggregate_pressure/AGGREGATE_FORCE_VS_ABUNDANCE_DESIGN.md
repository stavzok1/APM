# Design note — aggregate "raw force" vs abundance-sum baseline (gene rung)

**Status:** design, not yet run. Parked plan.
**Rung:** gene-level aggregate incoming pressure (see skill `apm-gene-question`,
`docs/GENE_QUESTION_TAXONOMY.md`). NOT the edge rung (where abundance wins, MH-78)
and NOT the program rung.

## 1. Question

Does an **aggregated repressive "force"** on a gene — abundance × edge strength,
attenuated by miRNA promiscuity — detect net miRNA repression of that gene
**better than a naive abundance-sum baseline**, and is the **tumor→healthy change**
in that repression more evident in the force form than in the baseline?

This lives in the one regime where the answer can be yes: the **aggregation lemma**
(summing breaks edge-level monotonicity, so construction is consequential at the
gene aggregate — unlike the edge).

## 2. Predictor construction

For gene `g`, sample `s`, regulator set `R(g)`:

```
P_force(g, s) = Σ_{m ∈ R(g)}  a(m, s) · W(m, g)        W(m, g) = w_eff(m, g) / D(m)
P_base (g, s) = Σ_{m ∈ R(g)}  a(m, s)
```

- `a(m, s)` = **raw** miRNA abundance (NOT z-scored). This is a deliberate
  *loading-dominance* model: the aggregate is carried by the highest-mass
  contributors. Carry **log-sum** as a grid contrast (raw-sum vs log-sum are NOT
  rank-equivalent at the aggregate, so this is a real axis, not cosmetic).
- `w_eff(m, g)` = edge strength (evidence-weighted, per spine).
- `D(m)` = miRNA **promiscuity** as a *budget split* — the miRNA's output capacity
  shared across its targets. **NOT** a softmax/normalize-to-1 share (that is the
  IDENTITY axis, ⊥ coupling — see `[[attribution-identity-vs-magnitude]]`,
  FORMULAS §5a). Because `Σ` over a gene's regulators is unconstrained, magnitude
  survives — this is force/MAGNITUDE, not identity.
- **`D(m)` is static (v1)** — same target-set size in tumor and healthy. Then
  `W(m,g)` is a constant per-edge weight, so force-vs-base reduces cleanly to:
  **regulators weighted by strength/promiscuity vs weighted equally.** A
  context-specific `D(m)` (expressed-target load) is a v2 axis — deferred because
  it injects target expression into the predictor (circularity risk vs `y_g`).

## 3. Estimand & estimator

Per gene, **cross-sample within a cohort**:

```
ρ_force(g) = partial-Spearman( P_force(g, ·), y_g(·) | confounds )
ρ_base (g) = partial-Spearman( P_base (g, ·), y_g(·) | confounds )
```

- **Spearman**, not Pearson (raw-abundance tails would otherwise drive Pearson).
- Confound set (per APM guardrails + the gene skill): **CPE/purity + proliferation
  + target-gene CN + target methylation** (+ HRD on the tumor arm). Report **raw
  and adjusted**; headline is **adjusted**. Proliferation is the load-bearing one
  (it co-drives both miRNA abundance and target mRNA — the MH-78 confound).
- A net-repression call needs a **coherent** stack (G8: `|signed| ≫` cancel) **and**
  FDR-significant adjusted ρ — not a large aggregate alone. Split by `role` /
  `malignancy_sign` (`gene_roles`).

## 4. Tumor vs healthy — the Δρ test

```
Δρ(g) = ρ_tumor(g) − ρ_healthy(g)        (adjusted)
```

- **Tumor** = TCGA-BRCA primary tumors. **Healthy** = **GTEx breast** (NOT NAT —
  true healthy, no field effect). GTEx miRNA already ingested on TCGA arm names
  (`gtex_mirna_matrix.py`, `data/GTEx/miRNA_TPM_matrix_v10.txt.gz`): **346 breast
  SMLRNA samples**, 514 RNA-seq → paired unit ≈ **~300** (confirm exact
  SMLRNA∩RNASEQ intersection). Well-powered, unlike the NAT/outcome arc.
- **Cross-cohort batch:** diff **within-cohort correlations** (ρ is location/scale
  invariant per gene; Spearman ranks travel across platforms better than levels).
  Do **not** harmonize abundance levels — let each cohort's ranks stand, diff the
  ρ's. Requires identical `R(g)`, `w_eff`, `D(m)`, and confound model in both;
  proliferation computed consistently. Purity ≈ 1 in GTEx (non-issue).

### 4a. RESULT (2026-06-29) — the Δρ axis is REAL and validated; force still ties

`tumor_gtex_delta_rho()` → `gene_force_delta_rho.tsv`. TCGA tumor vs **GTEx breast (327
paired)**, symmetric composition+proliferation confounds (tumor-only purity/HRD/CN have no
GTEx analog), predictors {base=Σa, validated-α0 force=Σ a·log1p(ev)}, curated ± sets:

| predictor | set | n | med ρ_tumor | med ρ_GTEx | med Δρ | pos vs neg (MWU) |
|-----------|-----|---|-------------|------------|--------|------------------|
| base | positive (oncomiR-TSG) | 18 | −0.097 | −0.016 | **−0.082** | **p=0.0028** |
| base | negative (miRNA-poor) | 556 | −0.008 | −0.001 | −0.007 | |
| force_val_a0 | positive | 18 | −0.102 | −0.017 | **−0.086** | **p=0.0021** |
| force_val_a0 | negative | 556 | −0.008 | −0.001 | −0.007 | |

**Two findings:**
1. **The change-axis works.** Known oncomiR-repressed TSGs show **acquired repression in
   tumor** (ρ −0.10) that is ~absent in healthy (ρ −0.02), and Δρ **significantly separates
   them from miRNA-poor genes** (p≈0.002–0.003). The curated positive control fires — the
   tumor-vs-healthy *change* in aggregate miRNA coupling is real and lands on the right genes.
2. **Force still ties abundance** (Δρ −0.086 vs −0.082; p 0.0021 vs 0.0028) — a hair better,
   not meaningful, exactly like the level. Binding-weight adds nothing here either.

**Arc verdict (weighting):** abundance-sum is the right gene-rung dose on the *level* **and** the
*change*; the *weighting* construction ("force") is a settled negative. The win is the **reframe** —
the Δρ (change) is the publishable result, not the weighting. Caveats: Δρ conditions on
composition+proliferation (symmetric), not target-CN (no normal analog); positive n=18.

### 4b. Rectified ACQUIRED abundance (2026-06-29) — the FIRST construction to move magnitude

Replace raw abundance with **`max(a_m(s) − h_m, 0)`**, `h_m` = QN healthy baseline on TCGA scale
(`healthy_anchor.gtex_qn_baseline`) — i.e. only *supra-healthy* miRNA exerts force (nonlinear, so
NOT a centering no-op). Added as predictor `acq` to the Δρ:

| predictor | positive ρ_tumor | positive Δρ | pos-vs-neg separation | MWU p |
|-----------|------------------|-------------|-----------------------|-------|
| base (Σa) | −0.097 | −0.082 | −0.075 | 0.0028 |
| force val-α0 | −0.102 | −0.086 | −0.079 | 0.0021 |
| **acq** `max(a−h,0)` | **−0.124** | **−0.114** | **−0.104** | 0.0022 |

**Effect size jumps ~30–40%** (positive ρ_tumor −0.124 vs −0.097; separation −0.104 vs −0.075) —
the first construction to beat abundance on magnitude. **But significance is flat** (p 0.0022 vs
0.0028): rectification drops samples/genes (negatives 556→507) and adds variance. **Methodological
caveat:** the GTEx-side `max(a_gtex − h, 0)` is ill-posed (rectifying healthy against its *own*
baseline → near-degenerate), so the Δρ partly double-counts the healthy anchor. **`ρ_tumor` is the
cleaner readout** — and there `acq` clearly strengthens detection of repression of known oncomiR-TSGs
(−0.124 vs −0.097, ~28%). **Lead to pursue:** score `acq` by `ρ_tumor` (not Δρ — acquired already
folds in the healthy contrast), recover power, test `acq × evidence`. First promising signal in the arc.

### 4c. Acquired by ρ_tumor (2026-06-29) — clean WIN (`gene_acquired_tumor_coupling.tsv`)

Dropped the degenerate GTEx-side difference; score `ρ_tumor` directly under **canonical confounds
(CPE+HRD+batch+target-CN)**, predictors {base, force val-α0, `acq`=Σ max(a−h,0), `acq_ev`=acq×log1p(ev)}:

| predictor | med ρ_tumor (pos) | pos neg-FDR | med ρ (neg) | pos-vs-neg MWU | separation |
|-----------|-------------------|-------------|-------------|----------------|------------|
| base | −0.123 | 13/18 | −0.006 | 1.3e-04 | −0.116 |
| force val-α0 | −0.124 | 13/18 | −0.006 | 9.1e-05 | −0.118 |
| **acq** | **−0.145** | **14/18** | −0.007 | **5.1e-05** | **−0.138** |
| acq_ev | −0.141 | 14/18 | −0.007 | **4.3e-05** | −0.134 |

**Acquired abundance beats raw abundance on effect size (−0.145 vs −0.123, ~18%), separation, recovery
(14 vs 13/18), AND significance (~2.5× better p).** First clear construction win in the arc; holds under
the canonical model. Evidence weighting remains inert alone (force≈base); `acq_ev` a hair more significant
but the **rectified-acquired transform is the lever**, not the weight. The user's framing — "abundance as
distance-above-healthy", nonlinearly — is what works (plain weighting/promiscuity never did).

### 4d. Proliferation gate PASSED + four independent confirmations (2026-06-29)

Re-ran with **proliferation added** (CPE+HRD+batch+target-CN+prolif) and four orthogonal reads
(`gene_acquired_tumor_coupling.tsv`):

| predictor | gw net-FDR | med ρ_tumor (pos) | pos-vs-neg MWU | separation |
|-----------|-----------|-------------------|----------------|------------|
| base | 444 | −0.161 | 1.07e-05 | −0.157 |
| force val-α0 | 444 | −0.178 | 6.5e-06 | −0.174 |
| **acq** | **460** | **−0.182** | **4.96e-06** | −0.180 |
| acq_ev | 459 | −0.193 | 3.68e-06 | −0.191 |

**Paired per-gene Δρ = ρ_acq − ρ_base:** positives (18) median −0.018, 95%CI [−0.026,−0.006] (excludes
0), Wilcoxon p=8.2e-04; **orthogonal TSG label (gene_roles, n=44)** median −0.004, CI excludes 0,
p=4.7e-04; **genome-wide (n=1330) Wilcoxon p=1.35e-14** (acq helps almost every gene, tiny but
ultra-consistent; +16 net-repressed recovered).

**Verdict: acq is a ROBUST win** — survives proliferation; confirmed by (1) ±-set MWU, (2) paired
per-gene Wilcoxon with bootstrap CI excluding 0, (3) an orthogonal TSG label, (4) a genome-wide signed-rank
(p=1e-14). Magnitude is **modest but real and overwhelmingly directional**. Mechanism clean: only
supra-healthy miRNA drives acquired repression.

**Correction (paired head-to-head, §4e):** the apparent "evidence weighting beats base" from the n=18
*medians* above (−0.178 vs −0.161) was **noise** — the paired per-gene Wilcoxon refutes it. Weighting is
genuinely inert; the weighting rehabilitation does not survive a paired test. See §4e.

### 4e. Full paired head-to-head (2026-06-29) — the definitive read

All constructions scored by ρ_tumor (+prolif), each **paired per-gene vs base** (Wilcoxon + bootstrap
95%CI; * = CI excludes 0 ⇒ robustly beats abundance):

| construction | positives(18) | TSG gene_roles(44) | all genes | verdict |
|--------------|---------------|--------------------|-----------|---------|
| `val_a0` evidence wt | p=0.051 (CI touches 0) | ns | medΔ≈0 (p=2e-3, magnitude 0) | **inert** |
| `val_a1` +promiscuity | medΔ **+0.038** p=0.99 | ns | worse | **harmful** |
| `ts_a0` TargetScan wt | p=0.64 | ns | ns | **inert** |
| `ts_a1` +promiscuity | medΔ +0.030 p=0.95 | ns | worse | **harmful** |
| **`acq`** | **medΔ −0.018, p=8.2e-4 \*** | **p=4.7e-4 \*** | **p=1.3e-14 \*** | **WIN** |
| **`acq_ev`** | **medΔ −0.017, p=7e-4 \*** | **p=2e-4 \*** | **p=2.9e-13 \*** | **WIN** |
| `acq_ts` | p=0.24 | p=0.089 | p=0.04 | not robust (TS coverage) |

**Definitive verdict:** (1) **edge weighting (evidence OR TargetScan) is inert** — paired test kills the
apparent benefit even with proliferation conditioned; (2) **promiscuity (budget-split α=1) is actively
harmful** — confirmed in both currencies; (3) **acquired abundance `max(a−h,0)` is the only robust win**,
confirmed paired across positives, an orthogonal TSG label, and genome-wide (p=1e-14); (4) `acq_ev ≈ acq`
— the evidence weight adds nothing on top, the **acquired transform is the whole lever**. This vindicates
the cross-sample-vs-cross-regulator mechanism: weighting reshuffles regulators (inert for a correlation),
`acq` reshapes each arm's sample trajectory (moves it).

## 5. The moderator — contribution effective-N (NEW)

Whether force can beat base is **gene-specific** and measurable: it depends on
whether `W` has enough dispersion to reorder the `a·W` ranking vs the `a` ranking.
Formalize per gene from the contribution shares `p_m = a_m·W_m / Σ`:

```
n_eff(g) = 1 / Σ_m p_m²        (or exp(entropy); cohort-median for stability)
```

Used as a **moderator**, this turns the "is it cosmetic?" worry into a built-in
falsification:

- `n_eff ≈ 1` → aggregate ≈ one edge → edge lemma → abundance wins → force
  **cannot** help. Do not claim a gain here.
- `n_eff` high → many comparable contributors → construction matters → **headline
  the force-vs-base comparison on these genes.**
- If force "wins" uniformly across `n_eff`, something leaks (likely proliferation).

This is a *contribution-weighted* upgrade of the existing raw-count "crowding"
moderator (see §7).

### 5a. First result (2026-06-29) — built, and it forces a choice in §2

`gene_force_concentration.py` → `output/coupling_predictor_comparison/
gene_force_concentration.tsv` (1,408 genes, 759 multi-regulator). Per-gene `n_eff`
and `top1_share` under base (`c=a`) vs force (`c=a·evidence/D(m)`, static degree),
for log and raw-linear abundance:

| abundance | median n_eff (multi-reg) | frac n_eff≥2 | promiscuity reorders top arm |
|-----------|--------------------------|--------------|------------------------------|
| **log** (`max(log2RPM,0)`) | 2.53 | 60% | **31.6%** |
| **raw-linear** (RPM) | 1.37 | 21% | 8.7% |

**Implication for §2's raw-abundance commitment:** raw-linear abundance makes the
aggregate **near-single-driver for ~79% of genes** (loading dominance), where the
edge lemma says force *cannot* beat abundance — so the raw-linear force-vs-base test
has headroom only on the ~160-gene multi-regulator tail. Under **log** abundance the
aggregate is genuinely multi-contributor for 60% of genes and promiscuity reorders
the top regulator for ~1/3 — far more headroom. **Open decision:** run the headline
under **log** abundance, keep raw-linear as a (tail-restricted) contrast — not the
reverse. Either way, **stratify/gate the force-vs-base comparison by `n_eff≥2`**;
testing the full gene set dilutes it with single-driver genes by construction.
(Ex: PTEN 87 regs — n_eff 55 log vs 4.4 linear; top miR-21-5p→miR-29c-5p under force.)

### 5b. Step 2 result (2026-06-29) — force does NOT beat abundance on the tumor *level*

`force_vs_base_coupling()` → `gene_force_coupling.tsv`. Partial-ρ(predictor, expr |
CPE+HRD+batch+target-CN) on TCGA tumor, net repression = neg-FDR(q<0.05):

| construction | stratum | n | base neg-FDR | force neg-FDR | Δ median ρ |
|--------------|---------|---|--------------|---------------|------------|
| log  | all       | 1408 | 398 | **374** | 0.0000 |
| log  | n_eff≥2   | 450  | 159 | **145** | +0.0066 |
| lin  | n_eff≥2   | 186  | 73  | **70**  | −0.0024 |

**Force underperforms abundance**, even in the multi-driver tail where it geometrically
*can* act. Built-in correctness check: for single-driver genes (n_eff<2) force ≡ base
under Spearman (scale-invariant), Δρ=0 exactly — so all action is in multi-reg genes
(as the lemma predicts), and there the promiscuity budget-split modestly *hurts* (it
downweights high-abundance hubs that carry the coupling — consistent with the MH-78/MH-44
abundance dominance).

**Consequence for the program:** the within-tumor *level* test is not where force wins.
The remaining places it still could: (i) the tumor↔GTEx **Δρ** — force may track the
*change* even if it doesn't beat the level (the real §4 test, GTEx already powered);
(ii) a stronger `w_eff` (TargetScan context-score) and/or **genome-wide** `D(m)` instead
of within-universe degree; (iii) **AGO-gating** the aggregate. If none of these moves it,
the honest headline is that abundance-sum is the right gene-rung dose and "force" adds
nothing — itself a clean, MH-78-consistent result.

### 5c. Corrected construction (2026-06-29) — supersedes 5b; same verdict, now *fair*

§5b used within-universe target **degree** as `D(m)` — a strawman (under-penalizes hubs,
literature-biased). Rebuilt canonically: `D(m)` = genome-wide binding **budget** (evidence-mass,
**HE edges** only for validated), via `genome_wide_promiscuity.genome_wide_budget`;
`W=w_eff/D^α` = fraction of the arm's budget on g. Both currencies (validated `log1p(evidence)`
/ TargetScan `|context++|`) × α∈{0,0.5,1}. Budget rankings confirm the bias contrast: validated
top = miR-21-5p (most-studied), TargetScan top = miR-15/16/195 (biggest seed targetome).

Fair coupling (TCGA tumor, matched gene sets, **Δmedian ρ ≈ 0.0000 throughout**):

| construction | net-repressed (matched base → force) |
|--------------|--------------------------------------|
| validated α=0 (raw evidence wt) | 398 → **397** (tie) |
| validated α=0.5 | 398 → 386 |
| validated α=1 (full budget-split) | 398 → **367** (worse) |
| TargetScan α=0 (777-gene cov.) | 277 → **278** (tie) |
| TargetScan α=1 | 277 → **267** (worse) |

**Verdict (now robust to the construction objection):** *no* force construction beats abundance-sum;
binding weight (evidence **or** TargetScan context++) only ties; **promiscuity penalty monotonically
hurts** in α — strong support for catalytic/sublinear promiscuity (÷D over-penalizes the
high-abundance hubs that carry the coupling). TargetScan is **not** stronger than validated. Coverage
caveat: TS budget exists for only **320 arms**, so TS-force drops uncovered edges (concentration
H2 1.66 vs base 2.53). The tumor *level* test is settled: abundance-sum is the gene-rung dose. The
**only** remaining live test is the tumor↔GTEx **Δρ** (§4) — the change, not the level.

## 6. Validation set

Genes **known to be more severely miRNA-downregulated in cancer** = positive
control; expect a larger/stronger force `Δρ` than base.

- **Non-circularity:** the label source must be **independent** of the edge
  construction (not the same miRTarBase/TargetScan that defines the edges).
- **Negative set required:** genes with little/no miRNA regulation. The claim is
  "force *separates* positive from negative," not "positive shows an effect."

## 7. Build plan — reuse, don't fork

- **Per-gene table home:** `output/coupling_predictor_comparison/
  gene_construction_role_coherence.tsv` already has `gene, n_regulators,
  coherence, partial_rho_cn, role, malignancy_sign`. `n_eff(g)` lands here next
  to `n_regulators`.
- **Moderator slot:** `coupling_predictor_comparison.py:473` already stratifies by
  "crowding" = raw `n_regulators`. Wire `n_eff` into that same stratification
  (contrast: does weighted `n_eff` moderate the force-vs-base gap better than raw
  count?). Extend `gene_construction_grid`, per the skill's "don't fork" rule.
- **Concentration input exists:** `pressure_attribution_validation.py:386` already
  computes "top-1 regulator share (budget concentration)" — `n_eff` is its
  principled multi-term generalization; reuse that share computation.
- GTEx side: `gtex_mirna_matrix.py` (donor×arm) + `healthy_anchor` already there.

## 8. Guardrails (do not skip)

- **No circular tuning.** Do not tune any promiscuity exponent (`D(m)^α`) to make
  force diverge from base — that optimizes the wrong objective. If `α` is a knob,
  select it by **held-out detection** on the independent ±sets (CV), not by
  separation from baseline.
- **Functional-abundance floor.** A steeper promiscuity penalty up-weights
  *specific-but-rare* miRNAs; a miRNA below its repression copy-number floor cannot
  act. Keep a min-abundance gate; don't let `α` push the aggregate onto
  sub-threshold regulators.
- **Gate the aggregate, not edges.** AGO-gating `P_force` is a legitimate gene-rung
  axis; prefer the proliferation-purified gate; report gated AND ungated.
- **Coherence + role.** Report split by `malignancy_sign`; a canceling stack is not
  net pressure.

## 9. Prerequisite status (checked 2026-06-29)

1. **✓ RESOLVED — 327 paired GTEx breast donors** (SMLRNA ∩ RNASEQ; 346 miRNA / 514
   mRNA). Both modalities present: `data/GTEx/miRNA_TPM_matrix_v10.txt.gz` (donor×arm,
   TCGA-arm join via `gtex_mirna_matrix.py`) + `data/GTEx/gene_tpm_v10_breast.gct.gz`
   (59k genes × 514, loaders in `family_normal_reference.py`). Well-powered.
2. **✓ RESOLVED — proliferation (the shared confound) computes on GTEx** (`_prolif_metagene`
   on GTEx mRNA: 514 samples, z-range [−1.8, 2.4]). `cross_state_coupling._state_covariates`
   **already builds the GTEx covariate frame** (ESTIMATE composition + proliferation +
   GTEx SMNABTCH/SMGEBTCH batch) — reuse it. Tumor-only confounds (purity, HRD, somatic
   target-CN) correctly **don't apply** to normal tissue → GTEx ρ conditions on the
   relevant-there set; the comparable cross-cohort confound is **proliferation**.
   Methylation not in repo (was sensitivity-only). ⇒ **Δρ confound model is ready.**
3. **✓ RESOLVED — hand-curated ± set built.** `annotations/breast_oncomir_tsg_positive_set.tsv`
   (+ `.README.md` for provenance/non-circularity). **Positive** = 22 curated breast
   oncomiR-repressed TSGs from functional literature; **18 scored** (PTEN, PDCD4, BCL2L11,
   CDKN1A/B, FOXO1/3, TGFBR2, CDH1, BRCA1, TIMP3, THBS1, LATS2, SPRY2, BTG2, SERPINB5,
   SOCS1, TPM1); 4 dropped (HOXD10/DICER1/RECK/SMAD4 — not in Hallmark universe). Positives
   sit in the multi-regulator tail (n_eff≥2), where the construction can act. **Negative** =
   648 miRNA-poor scored genes (`n_regulators ≤ 1`). Label is literature-derived, not from
   any Δρ/coupling output → non-circular.

**Verdict: the primary Δρ (per-gene force vs base, tumor vs GTEx) is buildable NOW** —
prereqs 1 & 2 met. The ± validation is a confirmatory layer to add after, in parallel
with curating the positive set. Build should **reuse `cross_state_coupling`** (it already
does tumor↔GTEx coupling at the edge rung; the Δρ is its gene-aggregate sibling), not fork.
