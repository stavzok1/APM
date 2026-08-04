# Patient-level question taxonomy — what "personalised" can mean, and which measure answers which

> **Scope: the PATIENT.** The 103 TCGA-BRCA patients with both a primary tumour (`01`) and a matched
> NAT (`11`). Question: *"what, if anything, is written in a patient's own normal tissue about their own
> tumour — and which of the many 'per-patient' quantities is actually being asked for?"* Companion to
> `EDGE_QUESTION_TAXONOMY.md` and `GENE_QUESTION_TAXONOMY.md`. Numbers live in `DISCOVERY_REGISTRY.md`
> (MH-158, MH-162, MH-228…MH-239); this file defines the **objects and their nulls**.

## 1. ⚠⚠ The rung trap — "per-patient" names TWO DIFFERENT RUNGS

The paired data is one `feature × patient` matrix. You can collapse it either way, and both results get
called "per-patient". **They are transposes and are not interchangeable** (axiom 7 — classify the
quantity by the call that produced it).

| rung | one value per… | computed ACROSS… | examples |
|---|---|---|---|
| **P-across** | **arm / edge / gene** | the 103 patients | `own_specific_frac`, `r_own`, `r_adj`, `b`, per-edge ρ |
| **P-each** | **patient** | the ~1,000–1,400 genes | `efficiency`, `rho_real`, `sim_own`, `own_minus_cohort` |

⛔ **The entire paired lane is P-across.** `_realize` correlates `pred` against `dy` **over the 103
patients**, so every ρ in it is a *cross-patient* quantity: the pairing cancels the baseline and then
discards the patient as the replicate. "Set-level only, never per-edge" is a consequence of that
sampling choice (per-edge SE ≈ 0.10 at n=103), **not** a fact about the data.

## 2. P-across measures (arm rung) — one value per arm

| name | definition | module | estimates |
|---|---|---|---|
| `own_specific_frac` | `var(cohort_NAT_mean − own_NAT) / var(own tumour−NAT shift)` ⇒ **≡ `var(N)/var(T−N)`, see §10** | `realization.dose_shift_arm` | how much of the patient's **dose shift** is their own baseline |
| `r_own` | `Spearman(tumour level, own-NAT level)` across patients | `mirna_nat_retention` | does the arm's **level persist** NAT→tumour within a patient |
| `r_perm` | same, patients permuted (15 draws); `retention = r_own − r_perm` | ″ | the pairing's contribution (`r_perm` ≈ 0.002) |
| `r_adj` | `r_own` with **each state residualised on its OWN CIBERSORTx block** | ″ | the composition-free version |
| `b` | OLS slope of tumour deviation on NAT deviation, composition-residualised | `offset_amplification` (MH-229) | **amplify (>1) / shrink (0–1) / replace (≈0)** |

`b = r · sd_T/sd_N` exactly — which is why `b` and `r_adj` correlate at **+0.636**, and why **only `b`
can separate amplification from replacement**; `r` alone cannot.

## 3. P-each measures (patient rung) — one value per patient

| name | definition | module | note |
|---|---|---|---|
| `efficiency` | `−corr_g(P_resid, Y_resid)`, **two-way centred** (gene AND patient main effects) | `patient_realization_efficiency` (MH-162A) | the canonical one — two-way centring is the correct form |
| `rho_real` / `gap` | within-patient ρ across genes, real minus site-free decoy | `patient_realization` (MH-235) | ⚠ **largely a rebuild of MH-162A**; one-way centred (weaker) |
| `sim_own`, `own_minus_cohort` | profile similarity of the tumour to own-NAT vs to the cohort-median NAT | `patient_nat_identity` | the cohort median wins for **93%** of patients |

## 4. ⭐ THE MEASURES ARE NOT ONE AXIS (MH-239)

Measured per arm:

```
corr(own_specific_frac, b     ) = −0.081      corr(r_own, r_adj) = +0.764
corr(own_specific_frac, r_adj ) = +0.224      corr(r_adj, b    ) = +0.636
```

⇒ **the DOSE axis and the RETENTION axis are orthogonal.** "How much of a patient's dose shift is their
own baseline" and "how much of their tumour level their NAT predicts" are **different properties**.
⛔ This does **not** support MH-158(C)'s "both axes converge on the patient's own-NAT baseline" at the
per-arm level. Exemplar: **miR-200c-3p**, `own_specific_frac` **0.85** (top decile) but `b` **0.028**.
Only **let-7f-5p** is high on both.

## 5. ⚠ "RETENTION" NAMES TWO UNRELATED QUANTITIES

| | definition | rung | reads as |
|---|---|---|---|
| **(A) patient-baseline retention** | `r_own` / `r_adj` / `r_own − r_perm` | arm | a **level persisting** across states within a patient |
| **(B) adjustment retention** | `rho_adj / rho_raw`; `composition_explained = <0.4` | edge | a **fraction of an effect surviving a covariate control** |

(B) is what `retention`, `median_retention`, `realized_retention`, `retention_rho`, `retention_beta`
mean on the cards, and what "within-PAM50 retention 0.81" / "Buffa LEVEL retention +0.358" mean. **The
word alone never disambiguates — read the producing call.**

## 6. The NULL taxonomy — three nulls, three different questions

| null | breaks | keeps | use it for |
|---|---|---|---|
| **arm substitution** (site-free decoy) | targeting | the pairing | is the signal site-mediated? |
| **patient permutation** | the pairing | targeting | is the signal patient-specific? |
| **split-half over genes** | — | — | is the per-patient value *reproducible*? |

⛔ **THE DECOY IS NOT A NULL SET (MH-238).** `corr(dX[real_arm], dX[fake_arm])` median **+0.211** over
1,722 pairs against a **+0.124** random-arm baseline — the Hungarian matching makes it **+0.087 more**
correlated than chance — and decoy edges hit `p_perm < 0.05` at **9.2%** (~1.8× nominal). ⇒ site-free
arms genuinely correlate with these targets, so **`REAL − DECOY` always understates** and the decoy
cannot serve as a false-positive floor. It is a *conservative* control, not a null.

✅ **The patient permutation is the clean one.** Use `coupling_inference.permutation_pvalues`
(Freedman–Lane), never a hand-rolled shuffle: it permutes the **residualised predictor** (a raw
permutation is only valid without covariates) and applies **one shared permutation across all rows per
draw**, preserving cross-edge/seed-family dependence that an independent per-edge shuffle would destroy.

⚠ Which vector is permuted differs by module and is **not** always the arm: `permutation_pvalues`
permutes the **ARM** (`null = (Xv[:, perm] * Yv).sum(axis=1)`); `patient_realization_efficiency` permutes
the **GENE** matrix columns; `mirna_nat_retention` permutes the **NAT arm**, per arm. For a plain
Spearman the choice is immaterial (symmetric); under covariate adjustment it is not.

## 7. What is settled

- **Amplification is EXCLUDED, not unobserved** — `b` = 0.183 arms / 0.137 genes, **>1 for 0.0% / 0.1%**;
  arms plateau at b=0.22 by the 4th NAT-expression quintile, so `b_true=1` needs 78% noise (MH-229).
- **No per-patient biomarker** — split-half reliability of `real − decoy` is **0.33** (MH-235), and
  MH-162A found a decoy *matching* the raw trait's reliability (+0.59 vs +0.50).
- **≈74% of the within-patient correlation is achieved by a STRANGER's pressure vector** — diag −0.113
  vs off-diagonal −0.084 on the 103×103 matrix (MH-235).
- **Hosted vs intergenic is the class axis** — `r_own` intergenic **0.028** vs coding-host **0.125** and
  lncRNA-host **0.116** ⇒ hosted-vs-intergenic, *not* coding-host co-transcription (MH-239, reproducing
  MH-158D to the third decimal). **let-7 takes 4 of the top 12** retention arms.
- ⚠ **"Imprinted" is NOT one class** — MEST-hosted miR-335-5p is 3rd in retention while the 14q32
  DLK1–DIO3 cluster is at the bottom (MH-232, p=0.0097 on the full 571-arm table).
- **The compartment-ORIENTATION stratum replicates under a second null** — MH-158's opposite-compartment
  gap (~0.017, p=2.7e-4, arm-substitution null) is **−0.022 at p=1.1e-06** under a patient-permutation
  null, same-compartment n.s. (p=0.154). Orientation is ΔCAF geometry, **not** coupling significance, so
  it is not circular. ⭐ **The single most load-bearing per-patient result.**
- **NAT is null POOLED, and every stratum that overturns it is CIRCULAR** — `field_established_realized`
  and `lost` are *defined* by `coupling_p_nat` significance while the outcome IS NAT coupling (MH-238).

## 8. ⛔ Do not rebuild on this cohort

| | why |
|---|---|
| the four-quantity per-patient vulnerability map | pre-loading dies against its decoy, 3.7e-06 → 0.061; no TSG enrichment (MH-233) |
| per-patient identity restricted to a "better" gene subset | the 1%→37% gradient is **entirely an expression artifact** — expression-matched controls give 0.275/0.269/0.285 (MH-234) |
| a per-patient realization biomarker | reliability 0.33; the ranking is ~2/3 noise (MH-235, MH-162A) |
| headroom | not established; the internal dose-response control fails (MH-230) |

## 9. Before adding a new per-patient measure

1. **State the rung** — P-across or P-each (§1). They are transposes; do not compare them.
2. **Write the covariance** if it touches both `N` and `T − N` — three of six §J items had a
   mechanically biased estimator, and J3's was biased *toward its own hypothesis* (memory
   `paired-design-shared-term-algebra`).
3. **Ship the decoy AND the patient permutation** — they answer different questions (§6), and the decoy
   alone understates.
4. **Report split-half reliability**, not only significance — a set-level effect with reliability 0.33
   is not a patient trait.
5. **Check circularity** — `gene_axes.scan_categorical(controls=…)` now emits `p_resid` against
   `p_matched`; compare to `p_matched`, never to `p` (different row sets).

## 10. ⚠ `own_specific_frac` IS NOT A FRACTION — and MH-158(A)'s wording does not follow

`nat_cohort` is a **per-arm constant**, so `var(cohort_mean − N) ≡ var(N)` and therefore

```
own_specific_frac  ≡  var(own NAT)  /  var(tumour − NAT)
```

verified identical to machine precision (max |diff| = 3.3e-16 over 571 arms; median 0.450, matching the
recorded 0.46). It is a **baseline-to-shift VARIANCE RATIO**, not a share:
**0.5% of arms exceed 1.0**, and **corr(own_NAT, own_shift) = −0.596** — the two components are far from
orthogonal (`cov(N, T−N) = cov(N,T) − var(N)`, negative whenever `cov(N,T) < var(N)`, i.e. almost always
at `r_own ≈ 0.11`). ⇒ **MH-158(A)'s "46% of the per-patient NAT→tumour dose-shift variance IS the
patient's own NAT baseline" is an invalid variance decomposition.** The quantity is real and usable;
the *share* reading is not. Same shared-term structure as the rest of this arc — `N` in both numerator
and denominator (memory `paired-design-shared-term-algebra`).

`r_own` for contrast is exactly what it looks like: raw `Spearman(tumour level, own-NAT level)` across
patients, **no covariates**; `r_adj` residualises each state on its own CIBERSORTx block.

## 11. Coverage map — what exists, at which rung, and the systematic gap

**Structural constraint first: COUPLING needs ≥2 samples, so coupling is inherently P-across.** There is
no within-one-patient coupling. Every "per-patient coupling" (`efficiency`, `rho_real`) is *across genes
within a patient* — a different rung, not a coupling. **Only LEVEL quantities can be genuinely per-patient.**

### LEVEL

| rung | within-state (NAT / tumour) | trans NAT→tumour | trans from constant GTEx ref |
|---|---|---|---|
| **arm** | `grank_*`, `arm_med_rpm`, `arm_iqr` — cohort | ✅ `own_specific_frac`, `mean_own_shift`, `r_own`/`r_adj`, `b`, `mean_dGlobalRank` | `dGlobal_HLY_NAT`, `dGlobal_HLY_TUM` — cohort |
| **family** | — | ✅ `_fam_own_specific_frac` only | ⬜ |
| **edge** | `share_NAT/TUM`, `rank_NAT/TUM`, `dose_rank_*` — cohort | ✅ `dShare_M_own`, `dShare_raw_own` — per-patient internally, **collapsed to a per-edge mean** | `share_HLY`, `rank_HLY`, `d_rank_HLY_NAT`, `d_rank_HLY_TUM` — cohort |
| **family_edge** | ⬜ | ⬜ | ⬜ |

### COUPLING (all P-across)

| rung | within-NAT | within-tumour | trans (Δ) | NAT→tumour cross |
|---|---|---|---|---|
| **edge** | ✅ `coupling_p_nat`/`z_nat`; §J9 `mode=nat` | ✅ `coupling_p_tum`; §J9 `mode=tumour` (**strongest, p=3.3e-16**) | ✅ `edge_rho_adj`; §J9 `mode=delta` | ✅ §J9 `mode=cross` — **null, 0 FDR edges** |
| **gene** | — | — | ✅ `realized_rho_adj` | ⬜ |
| **family** | — | — | ✅ `family_rho_adj`, `realize_family` | ⬜ |
| **family_edge** | ⬜ | ⬜ | ⬜ | ⬜ |
| **patient** (P-each) | — | — | ✅ `efficiency`, `rho_real` | — |

### ⬜ THE SYSTEMATIC GAP

**Every `share_*` / `rank_*` / `dose_rank_*` / `grank_*` / `d_rank_*` is a function of EXPRESSION
evaluated at the COHORT aggregate.** Each has a well-defined per-patient version that has never been
computed:

1. per-patient **budget share** (`β·exp` normalised within gene) per state
2. per-patient **budget rank** within gene, per state
3. per-patient **Δrank NAT→tumour** (only the cohort `d_rank_NAT_TUM` exists)
4. ⭐ per-patient **Δrank from the GTEx reference** — the cleanest of the four, because **GTEx is a
   CONSTANT reference for every sample**, so `rank_HLY` is a genuine fixed baseline and each patient's
   deviation from it is unambiguous — unlike the cohort NAT mean, which is estimated *from the same 103
   patients* and is precisely what makes `own_specific_frac` a shared-term quantity (§10)
5. **family and family_edge** versions of all of the above — the `family_edge` rung is empty throughout

⭐ **The infrastructure is already half-built:** `dose_shift_edge` and `dose_shift_arm` both
CONSTRUCT the per-patient vector (linear-RPM share matrices, `nat_dev`) and then **collapse it to a
per-edge/per-arm mean without persisting it**. Persisting those vectors is most of the work.

## 12. RETENTION — the consolidated result (MH-239, MH-243…MH-247)

**Definition.** `retention` = `r_own − r_perm`, where `r_own` = **Spearman(tumour LEVEL, NAT LEVEL) across
patients** — level-vs-level, **no differencing and no centring** (`r_perm` ≈ −0.002). The deviation-based
sibling is MH-229's slope, related exactly by `b = r · sd_T/sd_N`. `r_adj` = each state residualised on
its **own** CIBERSORTx block. ⚠ **Quote the ADJUSTED figure** — composition costs **16.2%** (share) /
**11.8%** (level), paired per cell. And `median(raw) − median(adj) ≠ median(raw − adj)`: for `r_level` the
medians move the *wrong* way while the paired cost is +0.0152. **Quote the paired cost.**

### The ladder (strongest first, composition-adjusted where marked)

| | value | note |
|---|---|---|
| **let-7 arms** | **+0.219** med, p90 **0.407** | adjusted; **2× background, SURVIVES adjustment** |
| **multi-arm-family cells** | **+0.141** med, p90 0.31 | adjusted; 1-arm cells 0.101 |
| P-each landscape diagonal | +0.143 | *different object* — across units within a patient |
| all-cell `r_share` / `r_budget_raw` / `r_level` | +0.134 / +0.133 / +0.129 | **raw** |
| **all-cell `r_share_adj`** | **+0.111** | ⭐ **the honest headline** |
| pre-session arm-rung `r_own` / `r_adj` | +0.110 / +0.105 | 571 arms |
| stranger-swap baseline | +0.035 | the floor |

**⚠ The tail is ~4× the median and every "≈0.13" above is a median:** `r_share_adj` p90 **0.256**, p95
**0.311**, p99 **0.404**, **max 0.456** (raw max 0.547).

### Where it is prominent — ⭐ persisted, not prose: `output/learned/realization/retention_atlas_top.tsv`

125 ranked entries across 3 rungs / 4 metrics. Highlights:

- **EDGE:** RB1|miR-335-5p **0.456** · NAMPT|miR-26b-5p 0.439 · NAMPT|miR-34a-5p 0.439 ·
  CCND2|let-7a-5p 0.436 · TGIF2|miR-34a 0.431 · ROCK1|miR-335-5p 0.430 · MYC|miR-335-5p 0.409
- **ARM** (median over the genes it regulates): miR-335-5p **0.419** · let-7e-5p 0.353 · miR-152-3p 0.333 ·
  miR-222-3p 0.303
- **ARM, pre-session independent estimator** (`mirna_nat_retention.r_adj`): miR-412-5p **0.665** ·
  miR-2355-5p 0.429 · let-7f-5p 0.416 · miR-335-5p 0.412
- **GENE** (β-weighted): NAMPT 0.439 · TRIB2 0.426 · VEGFC 0.346 · FASLG 0.338 · STAT5A 0.328 · CDKN1C
  0.325 · HMGA1 0.326 · BTG2 0.309
- ⭐ **miR-335-5p (host MEST, imprinted 7q32) takes 3 of the top 10 edges AND tops the arm rung AND is 4th
  on the independent pre-session estimator.** The most reproducible single entity on this axis.
- ⚠ **"Imprinted" is not one class** — MEST-hosted miR-335 is top; the 14q32 DLK1–DIO3 cluster is at the
  BOTTOM (MH-232, p=0.0097).

### ⛔ The ceiling — read before citing any of the above

**Retention is graded by how well an arm is MEASURED, not by what it does (MH-247).** The arm-covariate
scan is a pure detection signature (`detection` +0.356 q=5.8e-08, `arm_pct_floor` +0.375, `arm_med_rpm`
+0.362, `grank_TUM` +0.359, `wshift_sd_dGlobalRank` −0.348), and the **designed orthogonal control
settles it**: `promiscuity` — the one axis documented and here confirmed abundance-orthogonal (corr with
abundance −0.011) — correlates with retention at **−0.010, p=0.87**, while abundance does at **+0.339**.
Every strong predictor is abundance-correlated; the one that provably is not predicts nothing.

**And importance does NOT select retained arms (MH-247), under any weighting:** unweighted +0.1128 ·
β-only +0.1106 (p=0.116) · abundance-fraction +0.1096 (p=0.628) · budget-share +0.1154 (p=0.104).
⚠ I predicted the abundance weighting would gain mechanically from the +0.339 gradient; it did not —
**a global covariate gradient does not imply a within-unit weighting gain** (the weighting acts inside one
gene's 2–8 arms). Per-gene the two weightings *can* diverge sharply — **BTG2 0.205 unweighted / 0.307
β-weighted / −0.029 abundance-weighted** — but that is an exemplar, not a population effect.

### Artifacts

`retention_atlas_top.tsv` (ranked entities) · `retention_variants_edge.tsv` (1,032 cells × 4 variants +
fam_size/host/nat_lvl) · `gene_weighted_retention_3ways.tsv` (211 genes × 3 weightings) ·
`mirna_nat_retention{,_compadj}.tsv` (571 arms, pre-session) · `budget_vs_level_retention_edge.tsv`.
