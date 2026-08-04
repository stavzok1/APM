# Edge-question taxonomy — what a single miRNA→target edge can be asked, and how

> **Scope: the EDGE.** A single `(miRNA arm m, target gene g)` pair, measured across samples in a
> cohort. This is the *classification layer* — it names the distinct questions, the design axes that
> parametrize them, and which are already built. **Numbers and claims live in the modules and
> `DISCOVERY_REGISTRY.md`**, not here. Companion to `FORMULAS.md` (§3–5, pressure construction) and
> `MODELING_FRAMEWORK.md` (§6.3, predictor-vs-resolution). For the **gene/program aggregate** the
> axes change (the softmax *helps* the aggregate while it *hurts* the edge — MH-44); that is a
> separate frame.

The recurring error this doc prevents: **conflating "which arm owns the gene" (attribution /
identity) with "does this arm actually repress the gene" (realized impact).** They are different
questions with different estimators, and the *share* answers the first, not the second
(`FORMULAS.md` §5a: share ⊥ coupling).

---

## 1. The coupling-invariance lemma (read this first)

**Claim.** For a within-sample-set **rank** coupling of a *single* edge, the predictor construction
matters **only through operations that are not per-arm-monotone in `x_m`** (`x_m = log2(RPM+1)` of
the arm).

- **Inert** (rank-identical → *identical* coupling): any centering reference (`x−μ_cohort`,
  `x−μ_subtype`, `x−μ_NAT`, `x−μ_healthy`), z-scaling (`/sd`), ratio (`x/ref`, ref>0), `log`,
  multiplication by a **constant** `edge_w`. Each is a strictly-increasing function of `x_m` for a
  fixed arm, so the per-sample ranks — and therefore the Spearman against any target — are unchanged.
- **Bites** (can change the coupling):
  1. **softmax / share** — mixes co-regulators into the denominator, so the share is *not* a
     monotone function of `x_m` across samples;
  2. **per-sample modulation** — e.g. the AGO/RISC gate `g(s)`, a sample-varying (not per-arm)
     factor that reorders samples;
  3. **changing the sample set** — subtype restriction, the paired NAT subset (a different rank pool);
  4. **paired-differencing** — `Δx = x_tumor − x_NAT` is a different estimand, not a transform of `x`;
  5. **the conditioning set Z** — confounders, `P_others` (this is where the action is);
  6. **the functional form** — rank vs linear (Pearson) vs interaction/spline.

**Why the softmax bites (point-1 correction).** Centering each arm by reference `c_m` inside a
softmax gives `sm(m,s) = e^{x_ms − c_m} / Σ_{m'∈R(g)} e^{x_m's − c_m'}`. Switching `c_m` (e.g.
cohort-median → healthy-median, `δ_m = med_m − h_m`) reweights arm `m` by `e^{δ_m}`:

```
sm_h(m,s) / sm_med(m,s) = e^{δ_m} · [Σ a_m's] / [Σ a_m's · e^{δ_m'}]      a_m's = e^{x_m's − med_m'}
```

The bracket **varies across samples** (the dominant co-regulators change sample to sample), so
`sm_h(m,·)` is *not* a monotone transform of `sm_med(m,·)` → the edge Spearman genuinely changes. So
a centering reference is inert for a *standalone* (linear) edge predictor but a **real,
sample-dependent lever once it feeds a softmax**. (Inert again in the degenerate cases:
single-regulator gene, or uniform `δ_m`.)

**Consequence.** *Reference* and *scaling* are **attribution / state axes** — they reshape potential,
occupancy, and the gene aggregate. They touch the *edge coupling* only by routing through (i) a
softmax or (iii)/(iv) a sample-design change. So: choose reference/scale on **biological-question**
grounds (absolute vs within-population vs acquired-over-normal), **not** to chase anti-correlation;
the only knobs that move a single-edge coupling are the **estimator** (Group II) and the
**sample/functional design**.

---

## 2. Axis Group I — predictor construction `P(m,s)`

One object: `P = comp( transform(x_ms; ref, scale); promisc, aff )`.

| Axis | Values | Biological question it sets | Inert for edge coupling? |
|---|---|---|---|
| **ref** (centering) | ∅ · cohort-μ · subtype-μ · NAT-cohort-μ · NAT-paired · healthy(GTEx)-μ | "relative to *what* baseline" — absolute / within-population / lineage / acquired-vs-normal | **Yes**, unless via softmax or paired/stratum design |
| **scale** | level (`x` / `x−ref`) · z (`/sd`) · ratio (`x/ref`) | magnitude-preserving vs purely-dynamical vs proportional | **Yes**, unless via softmax |
| **comp** (competition) | standalone · softmax-share over R(g) | dose vs occupancy/attribution | **No** — route (i) |
| **promisc** (in logit) | none · degree · evidence/ts/combined-mass | "dilution of the arm's molecules over its target set" → effective per-gene dose | only via softmax |
| **aff** (in logit) | none · `log edge_w` | affinity×abundance → fractional **occupancy** (Boltzmann) | only via softmax |

Notes:
- **`z` is not its own axis.** It is the cell (ref = cohort-μ) × (scale = sd). Reference ⊥ scaling:
  `(x−μ_healthy)/sd_tumor` and `(x−μ_subtype)/sd_subtype` are legal cells the engine does not ship
  (it ships the diagonal). `pressure_engine.ExprMode` enumerates the implemented cells.
- **Promiscuity in the logit** (`logit_m − log promisc_m`) is a *dilution-over-targets* model
  (molecules spread across many sites → less per gene). Currently `target_norm` divides the edge
  weight **outside** the softmax (`pressure_engine.mirna_mass_denominators`); putting it in the logit
  is the in-softmax variant (`_softmax_gene_logits` adds `log edge_w`, the affinity case).
- **`aff` in the logit = a competitive-binding (Boltzmann) occupancy model**:
  `occupancy_m ∝ edge_w·e^{x−ref}`, normalized over R(g). Implemented as `softmax_z_logrpm_spec` /
  `softmax_z_logrpm_enc`. Sharpens *who owns the gene* (E7), ⊥ realized coupling.

---

## 3. Axis Group II — coupling estimator `assoc(P, Y | Z, S, design, functional)`

| Axis | Values | What it controls |
|---|---|---|
| **predictor** | raw `x_m` · share-based | mostly collapses by the lemma; share is the only rank-changer |
| **confounders Z₀** | CPE · HRD · target-CN · batch · proliferation · composition · **host-gene** (intragenic) | honesty |
| **co-reg conditioning** | none · **precision** (resid Y only on `P_others`) · **attribution** (resid both) | power-gain vs honest-unique-effect |
| **stratum / set S** | whole cohort · subtype · state (GTEx/NAT/tumor) | "where the coupling lives" |
| **design** | cross-sectional · paired-delta (NAT subset) | association vs within-patient change |
| **functional form** | Spearman (rank) · Pearson · interaction/spline | monotone-marginal vs nonlinear/saturating |
| **sign handling** | signed · positive-part · abs | repression-only vs net |

**Precision vs attribution (the key Group-II distinction).** `P_others(m,g,s) = Σ_{m'≠m} c(m',g,s)`
is the gene's competitive field.
- **Precision (semi-partial):** residualize the **target** on `Z₀ ∪ {P_others}`, keep the raw focal
  predictor. Removes background *noise from the target*; **gains** power when the field is ~independent
  of the focal arm. This is the variance-reduction (ANCOVA) move.
- **Attribution (partial, both-sides):** residualize **both** sides on `Z₀ ∪ {P_others}`. Honest
  unique focal effect; **costs** power only where the field correlates with the focal arm — and that
  cost is honest (shared variance can't be attributed to `m`).

---

## 4. Axis Group III — object definition (what *is* the edge / target)

| Axis | Values | Why it's a real fork |
|---|---|---|
| **regulator unit** | arm (5p/3p) · mature · seed-family | family shares targets → the "edge" can be family-level |
| **target layer** | mRNA · protein (L2) · protein-beyond-mRNA (L1b) · APA-isoform | "anti-correlation vs *which* target" |
| **effective dose** | nominal `x_m` · AGO/RISC-gated `g(s)·x_m` | saturation modifies the dose per sample (route ii) |

---

## 5. Named estimands (E1–E10)

| ID | Question | Estimator | Reference / axis | Class |
|---|---|---|---|---|
| **E1** | How much repressive *dose* does `m` put on `g`? | `P = aff·gate(s)·x_m` (per-sample) | abs / cohort / healthy | potential |
| **E2** | Does `m`'s level track `g` down, net of confounders? | `partial-ρ(x_m, y_g \| CPE,HRD,CN)` | within-cohort | realized (spine) |
| **E3** | Does `m` repress `g` *uniquely*, after the other regulators? | attribution: resid **both** on `Z₀ + P_others` | within-cohort | realized |
| **E4** | Does `m` repress `g` once background *noise* is stripped from the target? | precision: resid **Y only** on `Z₀ + P_others`, raw `x` | within-cohort | realized (power) |
| **E5** | Does the tumor-*acquired* change in `m` track the acquired change in `g`? | paired-delta `ρ(Δx, Δy)` same-patient NAT | acquired / paired | realized |
| **E6** | Where does the coupling live? | E2 within PAM50 / state | subtype / state | realized |
| **E7** | Of `g`'s regulators, how much of the budget does `m` own? | softmax share, logit `(x−ref)[/sd] − log promisc + log edge_w` | per-sample, vs co-regs | attribution (⊥ coupling) |
| **E8** | Does the repression realize on protein, through vs beyond mRNA? | E2 with `Y` = protein (L2) / protein⊥mRNA (L1b) | within-cohort, CPTAC | realized |
| **E9** | Is the coupling independent of co-transcription with the host gene? | E2 + host-gene in `Z` (intragenic miRNAs) | within-cohort | realized (honesty) |
| **E10** | Does `m`'s effect *saturate / depend on* background load? | `+ x_m·P_others` interaction, or spline | within-cohort | realized (nonlinear) |

---

## 6. Done / not-done grid (living tracker)

| ID | Status | Home module(s) | Output anchor |
|---|---|---|---|
| E1 | **built** | `pressure_engine`, `pressure_build`, `ago_gate` | `output/edges/`, contribution maps |
| E2 | **built** | `hallmark_interaction`, `coupling_predictor_comparison` (edge), `mirna_state_class`, `edge_partial_corr_panels` | `output/coupling_predictor_comparison/edge_coupling_by_predictor.tsv.gz` (abundance: ~1,472 neg-sig, 2026-06-26) |
| E3 | **run (lost set only)** | `decoupling_validation` block 6 | `output/decoupling_validation/decoupling_validation.tsv` (Δρ≈0.008, ~210 edges) |
| **E4** | **run — NEGATIVE** (2026-06-27, MH-78): conditioning the target on `P_others` does **not** beat raw abundance. Precision 1,471→**1,369**; attribution →**1,390**; median Δρ≈0.0001; predicted precision-gain absent → the field is **not** a target-noise source | `coupling_predictor_comparison` (conditioning sweep) | `output/coupling_predictor_comparison/edge_conditioning_{comparison.tsv.gz,summary.tsv,mechanism.tsv}` |
| **E4′ logit grid** | **run — NEGATIVE** (MH-78): no bare share under any logit (reference cohort/healthy/NAT × measure rawx/dev/z × promiscuity) beats abundance (best ~906 vs 1,471); paired-NAT apparent edge = low-power noise (n=111) | `coupling_predictor_comparison --logit-grid` | `edge_logit_grid_{comparison.tsv.gz,summary.tsv}` |
| **E4″ function/occupancy** | **run — NEGATIVE** (MH-78): mass-action occupancy is closest (1,095/1,230 ≈90% under precision — reduces to saturating abundance at single-reg genes) but still loses; softmax/temp/sparsemax/linear all below; **AGO-gating HURTS** (median ρ flips +; full/partial/purified) | `coupling_predictor_comparison --function-grid` | `edge_function_grid_{comparison.tsv.gz,summary.tsv}` |
| E5 | **open** (only cohort-level cross-state exists, `cross_state_coupling`) | — | — |
| E6 | **partly built** | `pam50_gene_resolution`, `visibility_archetype_contrasts`, `subtype_contrasts`, `target_combined_anticorr` (by_pam50) | `output/target_combined_anticorr/target_combined_anticorr_by_pam50.tsv` |
| E7 | **built** (⊥ coupling) | `pressure_engine` modes, `encori_share_sensitivity`, struct-share (`mirna_state_class`) | `output/tissue_reference/mirna_state_class/struct_vs_abundance_*.tsv` |
| E8 | **built** | `cptac_validation`, `cptac_target_specificity`, `cptac_orphan_discovery` | `output/cptac_*/` |
| E9 | **open** | — (no host-gene/intragenic conditioning anywhere) | — |
| E10 | **open** | — (existing "interaction" hits are survival-model RSF, not coupling-interaction) | — |

---

## 7. Open levers + priority

**The predictor-construction space is now exhausted (MH-78).** Across conditioning (E4), logit
construction (E4′: reference × measurement × promiscuity), competition function (E4″:
softmax/temperature/sparsemax/linear/mass-action), AGO-gating (full/partial/purified), and a stricter
prune — **nothing beats raw ungated abundance**; mass-action occupancy is the closest (≈90% under
precision, by design) but never overtakes. Per the lemma, reference/scale are inert at the edge except
through the softmax, and the softmax/competition only *dilutes* the focal edge. So the lever is **not**
the predictor. Remaining untested moves are estimator/design-side and now lower-priority:

1. **E10 — interaction term** (`x_m·P_others`). Whether competition is nonlinear/saturating — the one
   form the linear conditioning could not capture. (Mass-action's saturation is a partial proxy and
   it didn't win, so the prior is weak.)
2. **E9 — host-gene conditioning**. Honesty fix for intronic hub arms (miR-25/93/106b in MCM7).
3. **E5 — paired-delta** over the matched-NAT subset (low power; n≈111 noise-dominated, per E4′).

The productive direction has moved to the **gene level**, where (unlike the edge) construction is
consequential — see `GENE_QUESTION_TAXONOMY.md` (and `PATIENT_QUESTION_TAXONOMY.md` for the per-patient rung).
