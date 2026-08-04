# Stratification playbook — which axes to profile a result against, and how to read the tail

> **Scope: any "where does this effect live?" question** on an arm, edge, family or gene outcome.
> Companion to `EDGE_`/`GENE_`/`PATIENT_QUESTION_TAXONOMY.md`, which fix the *object*; this fixes the
> *profiling*. Written after §J + the coverage map (MH-228…MH-247), where six controls flipped a result
> and two of my own summary statistics were meaningless before they were fixed.

## 0. The tool is `learned/gene_axes.py` — do not hand-roll a tertile split

`scan(outcome, axes, circular=)` FDR-scans every numeric axis · `scan_categorical(outcome, cats,
circular=, controls=)` for class columns (⚠ `scan` is **numeric-only and silently skips them**) ·
`contrast(outcome, axes, k)` profiles the extremes · `sign_analysis` asks whether a gain is a **sign
correction** rather than a sharpening · `mask_degenerate` · `regulator_axes` / `weight_axes` / `self_axes`
build the axis blocks. `MIN_N_SCAN = 100` — below that `scan` returns empty, which looks like "no axes"
and is really "not enough units".

## 1. The axis inventory, by rung — with the actual column names

### ARM (`arm_card.tsv`, `learned/genomic_context.tsv`, `gene_axes._promiscuity_map`)

| block | columns |
|---|---|
| **abundance / detection** | `arm_med_rpm` · `arm_pct_floor` · `detection` · `spiker` · `arm_iqr` · `grank_{HLY,NAT,TUM}` |
| **cross-patient dynamics** | `wshift_sd_dGlobalRank` · `wshift_mean_dGlobalRank` · `wshift_own_specific_frac` · `wshift_mean_own_shift` |
| **genomic context** | `mir_class` (intergenic / sense_coding_host / sense_lncRNA_host) · `host` · `host_type` · `strand_rel` · `n_loci` |
| **state shift** | `arm_lfc_HLY_{NAT,TUM}_{raw,QN}` · `arm_lfc_NAT_TUM` · `dGlobal_{HLY_NAT,NAT_TUM,HLY_TUM}` |
| **loading / evidence** | `ago_loading` · `ago_loading_measured` · `chim_tarbase_{weight,n_genes}` · `chim_manakov_n_genes` · `he_degree*` |
| ⭐ **promiscuity** | `gene_axes._promiscuity_map()` — **the ONLY axis here that is abundance-ORTHOGONAL** (measured corr with abundance −0.011). See §4. |

### GENE (`realization/gene_card.tsv` — ⚠ **not** `learned/gene_card.tsv`, which does not exist)

| block | columns |
|---|---|
| **architecture** | `n_arms` · `n_fam` · `ctx_n_abund` · `ctx_n_fam_abund` · `ctx_frac_abund` · `ctx_n_fam_multi` |
| **weight / identity** | `top_beta_frac` · `concentration` · `frac_identified` · `n_identified` · `max_beta_frac_sd` · `identity_eq_magnitude` |
| **dose** | `ctx_dose_max` · `ctx_dose_med` · `ctx_gap_d_dose` |
| **categorical** | `dominant_edge_shift_class` (7 levels) · `gene_repression_class` (4) · `ctx_apriori_class` (4) |
| **the gene's own profile** | `gene_axes.self_axes` → `self_mean`, `self_sd`, `self_iqr`. ⚠ **sd and IQR predict; the MEAN does not** (MH-201). Always carry a dispersion term. |
| **role** | `gene_roles.load_gene_roles` → tsg / oncogene / dual / unknown |

⚠ **`w_max` is the max curated EVIDENCE weight over the gene's arms — NOT a β or dose share.** Values run
3.8–8.0, not 0–1. The actual share columns are **`top_beta_frac`** (median 0.937) and **`concentration`**.

### EDGE (`realization/edge_card.tsv`)

`coupling_z_{tum,nat,hly}` · `coupling_p_*` · `beta`, `beta_sd`, `z`, `beta_frac`, `pip_dense`,
`pip_discovery` · `share_{HLY,NAT,TUM}`, `rank_*`, `d_rank_*`, `dose_rank_*` · `shift_class` ·
`retention` (⚠ **sense (B)** — see the patient taxonomy §5) · `realization_score` · `orientation`
(SAME/OPPOSITE compartment, from ΔCAF geometry — **not** coupling-derived, so it is the rare
non-circular stratifier).

## 2. The mandatory controls — each of these flipped a result in §J

| # | control | why |
|---|---|---|
| 1 | **`mask_degenerate`** | one-unit genes make a weighted estimator ≡ unweighted **by construction**. **~50% of this universe** (MH-237: 531/1,070; the coverage map: 105/224). |
| 2 | **collinearity check first** | `n_fam ↔ n_arms` **ρ=0.97**; `top_beta_frac`/`concentration` are its inverse. Five "findings" were one axis (MH-236). |
| 3 | **conditional grid** when axes are collinear | equal-n `a × b`; MH-236's marginal `n_fam` gradient **reversed** conditionally on `pressure_sd` — and that reversal then turned out to be the decoy (§3). |
| 4 | **n-matched *and* EXPRESSION-matched** random control | n-matched alone gave 0.011 where expression-matched gives **0.275** at the same n ⇒ *any* well-expressed set looks like a 23× win (MH-234). |
| 5 | **circularity** | `scan_categorical(controls=…)`; compare `p_resid` to **`p_matched`**, never to `p` (different row sets: p=0.057 unmatched vs **0.0011** matched, same data). |
| 6 | ⭐ **an ORTHOGONAL control** | see §4 — the only thing that settles reliability-vs-biology. |

## 3. ⚠ VALIDATE THE CONTROL BEFORE READING THE CONTRAST IT MAKES

A control is a measurement too. Three ways they failed here:

- **The decoy LEAKS.** `corr(dX[real_arm], dX[fake_arm])` median **+0.211** against a **+0.124**
  random-arm baseline — the Hungarian matching makes it *more* correlated than chance. Decoy edges hit
  `p<0.05` at **9.2%** ⇒ **the decoy is not a null set**, `REAL − DECOY` always understates, and a
  conditional *reversal* it produced (MH-236) **vanished** when re-run REAL-only.
- **The matched control was matched on the wrong thing** (§2.4).
- **A "verify" function that asserts instead of testing.** `verify_gtex_constant` originally returned a
  hardcoded `True`. A check that cannot fail reads as coverage.

## 4. ⭐ The orthogonal control — the only clean reliability-vs-biology test

When every significant axis is abundance-flavoured, you cannot tell a measurement gradient from biology by
adding more abundance-flavoured axes. **Find a covariate provably orthogonal to the suspected confound and
test it.** Here: `promiscuity` (corr with abundance **−0.011**) vs retention = **−0.010, p=0.87**, while
abundance vs retention = **+0.339**. ⇒ retention is a **reliability** gradient — settled, not inferred.

⚠ **And a global gradient does NOT imply a within-unit weighting gain.** Abundance graded retention at
+0.339 across 267 arms, yet weighting a gene's 2–8 arms by abundance gained **+0.0000, p=0.628** — the
within-gene abundance range is narrow. (Same family as `pooled-effect-may-be-between-unit`.)

## 5. Characterising the TAIL — medians lie about dispersed quantities

**Always report `median · p90 · p95 · p99 · max`, not the median alone.** Measured here: adjusted
retention median **0.111** but p90 0.256, p99 0.404, **max 0.456** — the tail is **~4×** the median, and
every "≈0.13" summary in MH-243/244/245 concealed it.

Then, before naming anything from the top:

1. **Is the ranking n-contaminated?** 5 of the top-10 retention arms were excluded by complete-case
   (n=75–102): lower n ⇒ noisier r ⇒ extremes enriched. **Report the ranking restricted to full-n units.**
2. **Is the top a class or a list?** let-7 taking 4 of the top 12 is a class result and survived
   composition adjustment (0.219 vs 0.109); a single arm at 0.66 may be an outlier.
3. **Does an entity repeat across independent estimators?** miR-335-5p tops the arm rung, takes 3 of the
   top 10 edges, **and** is 4th on the pre-session estimator — that is what reproducibility looks like.
4. ⛔ **Never headline a ratio/threshold where the mass piles up** (axiom 5). `frac_budget_retained` had
   median **0.000** and p90 **0.979** — binary at 2 arms, therefore not reportable.
5. **Say what the sample was.** `_he_genes()[:N]` is **alphabetically sorted**; an alphabetical slice
   leaves distributions intact but makes every named top-entity an A–C artifact. **Sample randomly.**

## 6. Reporting checklist

- [ ] rung stated; `mask_degenerate` applied
- [ ] collinearity table before any marginal gradient; conditional grid if ρ>0.5
- [ ] random control **expression-matched**, not just n-matched
- [ ] categorical scans carry `controls=`; `p_resid` read against `p_matched`
- [ ] an orthogonal control if the axes are all confound-flavoured
- [ ] **median · p90 · p95 · p99 · max**, not the median alone
- [ ] top-N restricted to full-n units, and read as class-or-outlier
- [ ] the gene/arm sample is random, and said to be
- [ ] every control validated (does the decoy leak? can the check fail?)
