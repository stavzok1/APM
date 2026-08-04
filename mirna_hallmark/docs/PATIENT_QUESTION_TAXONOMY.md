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
| `own_specific_frac` | `var(cohort_NAT_mean − own_NAT) / var(own tumour−NAT shift)` | `realization._within_patient_shift` | how much of the patient's **dose shift** is their own baseline |
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
