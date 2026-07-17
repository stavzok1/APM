# Aggregate miRNA pressure — definitions, tests, results

Reference card. Minimal prose. Chronological design log: `AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md`.
Code: `gene_force_concentration.py`, `genome_wide_promiscuity.py`. Outputs: `output/coupling_predictor_comparison/`.

---

## 1. Cohorts

| cohort | source | n | use |
|--------|--------|---|-----|
| tumor | TCGA-BRCA primary, Normal-like excluded | per-gene partial-corr **n=1077** (common miRNA∩mRNA∩covariate set; MIN_N=40) | ρ_tumor |
| healthy | GTEx breast | 346 miRNA, 514 mRNA, **327 paired (SMLRNA∩RNASEQ)** | ρ_GTEx, healthy baseline h |

## 2. Expression units & normalization

| matrix | unit |
|--------|------|
| miRNA tumor | log2(RPM+1), TCGA GDC |
| miRNA GTEx | log2(TPM+1), joined to TCGA arm names (MIMAT-primary) |
| mRNA tumor | log2 normalized (`D.load_rna`) |
| mRNA GTEx | log2(TPM+1) (`gene_tpm_v10_breast.gct.gz`) |
| abundance `a` | log2RPM clipped at 0 (`max(log2RPM,0)`) |
| **healthy baseline `h`** | **`gtex_qn_baseline`** = per-arm GTEx-breast median **quantile-normalized to the TCGA scale** (`healthy_anchor.py`); resolves cross-platform offset; used in acquired transform |

## 3. Evidence definitions (miRTarBase)

- **Canonical `evidence_score`** (config `PRESSURE_EVIDENCE_SCORER = "confidence_logclass"`, scorer
  `score_confidence_logclass` + ENCORI S1 boost). Formula:

  ```
  evidence_score(m,g) = Σ_c w_c·log1p(nf_c)  +  0.3·Σ_c w_c·log1p(nw_c)  +  0.5·enc_depth(m,g)
    c ∈ {reporter, protein, rna, perturbation, binding}
    w  = {reporter 3.0, protein 2.5, rna 1.5, perturbation 1.5, binding 0.5}   (assay directness)
    nf_c / nw_c = # Functional-MTI / # "Functional-MTI (Weak)" studies of class c   (weak ×0.3)
    enc_depth   = ENCORI CLIP read-depth for (m,g); α = PRESSURE_ENCORI_ALPHA = 0.5
  ```
  Range on analysis edges: [2.08, 15.01].
  **Justification (the run):** `log1p` *within class* so a pair with 10 reporter studies is not 10×
  a pair with 1 → defuses publication bias (the main reason fame, not biology, picks within-gene
  winners); CLIP/binding weighted 0.5 (physical proximity ≠ repression); weak-curator studies ×0.3.
  Validated by the **curation-bias robustness run (MH-18, `robustness_checks.py`)**: the hub ranking is
  *not* a study-count artifact — the top miRNA ranks #1 under evidence, binary/degree, **and** TargetScan
  context++ (rank-corr ρ=0.65 evidence-vs-TargetScan). *(This replaces the build_edges default
  `3·reporter+3·binding+2·protein+1·rna+1·perturbation` study-count scorer.)*
- **Analysis edge set** = miRTarBase Hallmark edges with **`confidence_logclass` score ≥ 2
  (`PRESSURE_MIN_EVIDENCE`) + ENCORI boost** — what `load_mirtar_edges` produces. (This is an
  evidence-score floor, *not* the older "Functional-MTI + low-throughput" HE mask.)
  **5,219 edges / 1,424 genes / 789 arms**; **5,101 / 1,408** after restricting to expressible arms
  (= the analysis set).
- weight currencies `w_eff(m,g)`: **validated** = `log1p(evidence_score)` (canonical, as above);
  **targetscan** = `abs(ctx)` where **`ctx` = TargetScan7 weighted context++ score** (hsa, summed over sites;
  more negative ctx = stronger predicted repression, so we use absolute value as a positive weight).
  TargetScan coverage is limited: its prediction file contains only **321 distinct hsa arms** (it scores
  conserved miRNA families, not all ~2,600 arms) → of the 789 analysis edge-arms, only **247 arms / 782
  genes / 2,222 edges** carry a ctx score. This is why all `ts_*` constructions are coverage-limited.
- The genome-wide budget `D(m)` (§4) uses the **same** canonical construction (confidence_logclass≥2
  + ENCORI), so the budget currency matches the edge w_eff currency (fixed 2026-06-29; conclusions unchanged).

## 4. Promiscuity definitions

- **`D(m)` = genome-wide binding BUDGET** = Σ over arm m's *entire* targetome of the per-edge weight (evidence-MASS, **not** a target count). Canonical: `genome_wide_promiscuity.genome_wide_budget`.
  - validated: Σ_{genome-wide, confidence_logclass≥2+ENCORI} `log1p(evidence)` — 1,036 arms, median D=3.55, top = **miR-21-5p** (study-biased)
  - targetscan: Σ_{genome-wide} `abs(ctx)` — 320 arms, median D=70, top = **miR-15/16/195** (seed-targetome)
- **budget-split weight** `W(m,g) = w_eff(m,g) / D(m)^α`; α∈{0,0.5,1}. α=0 = raw weight (no promiscuity); α=1 = fraction of budget.

## 5. Predictor constructions (gene aggregate `P(g,s)=Σ_{m∈R(g)} c_m(s)`)

| name | `c_m(s)` | axis |
|------|----------|------|
| base | `a` | abundance dose |
| val_a0 / val_a1 | `a·log1p(ev)` / `a·log1p(ev)/D_val` | weight / weight+promiscuity |
| ts_a0 / ts_a1 | `a·abs(ctx)` / `a·abs(ctx)/D_ts` | weight / weight+promiscuity |
| **acq** | `max(a − h, 0)` | acquired (supra-healthy) magnitude |
| acq_ev / acq_ts | `max(a−h,0)·log1p(ev)` / `max(a−h,0)·abs(ctx)` | acquired × weight |

## 6. N-effective abundance concentration (descriptor, not a test)

- **Computed on cohort-MEDIAN abundance per arm** (`a = median_s log2RPM`), i.e. **one value per gene
  from the median-abundance contribution mass — NOT a per-sample mean of n_eff.**
- per-gene contribution shares `p_m = c_m/Σc` (with `c_m = median-abundance × W`); **Hill spectrum**:
  H0=`n_regulators`, H1=`exp(−Σ p·lnp)`, **H2=`1/Σp²` (primary)**, H∞=`1/top1_share`.
- **Purpose:** moderator / falsification. `n_eff≈1` ⇒ aggregate ≈ one edge ⇒ edge lemma ⇒ construction
  *cannot* beat abundance; `n_eff` high ⇒ construction *can* act. Tells you WHERE a win is possible.

| construction (multi-reg genes) | H1 | H2 | H∞ | frac H2≥2 |
|---|---|---|---|---|
| base log | 2.81 | 2.53 | 2.04 | 60% |
| base raw-linear | 1.64 | 1.37 | 1.19 | 21% |
| force_val log | 2.79 | 2.45 | 1.90 | 59% |
| force_ts log | 1.87 | 1.66 | 1.36 | 27% |

Budget-split reorders top regulator vs abundance: validated 33%, ts 27%. **Use log abundance** (raw-linear collapses to single-driver). Single-driver genes: force ≡ base by Spearman scale-invariance.

## 7. Estimator & confounds

- coupling ρ = **partial Spearman(predictor, gene mRNA | confounds)**; more negative = stronger repression. MIN_N=40.
- **canonical tumor confounds = CPE + HRD + sequencing-batch + target-gene-CN + proliferation**
  (proliferation = E2F+G2M metagene; canonical here, the load-bearing confound per MH-78). The
  level-coupling head-to-head (Test 1) uses this full set.
- **Δρ confounds = "symmetric composition + proliferation"**, defined precisely:
  - *composition* = epithelial + immune + stromal **marker metagenes** + proliferation metagene. **Source:**
    hand-curated marker gene sets (`EPI_MARKERS` / `IMMUNE_MARKERS` / `STROMA_MARKERS` in
    `cross_state_coupling.py`; proliferation = E2F+G2M), each a mean-z over its markers from the cohort's mRNA.
  - *symmetric* = the **same** signatures computed **identically from each cohort's own mRNA**, so ρ_tumor
    and ρ_GTEx are adjusted by an identically-constructed frame (required for the difference to be comparable).
    Tumor-only CPE/HRD/target-CN/batch are **excluded** (no GTEx analog).
  - *Why not CIBERSORTx?* CIBERSORTx fractions exist **for TCGA only** (`annotations/cibersortx_brca_participant_features.tsv`;
    a prior core-coupling CIBERSORTx deconv retest is at `output/core_coupling_deconv_retest_cibersortx/`). They
    are a more principled composition estimate and **should be used as a sensitivity for the tumor-only Test 1**,
    but the Δρ needs *symmetric* covariates — using CIBERSORTx would require running it on **GTEx** with the same
    signature matrix (not done). The metagene was chosen because it is trivially symmetric across cohorts.

## 8. Controls

| set | definition | n (scored) |
|-----|-----------|------------|
| positive | curated breast oncomiR-repressed TSGs, literature-derived (non-circular); `annotations/breast_oncomir_tsg_positive_set.tsv` | 18 |
| negative | miRNA-poor genes (`n_regulators ≤ 1`) | 648 |
| orthogonal positive | `gene_roles` TSGs (`malignancy_sign=−1`) | ~44 |

---

## TEST 1 — level coupling ρ_tumor (full head-to-head; canonical confounds incl. proliferation)

- **Statistic:** ρ_tumor = partial Spearman(P, mRNA | CPE+HRD+batch+target-CN+**proliferation**). Three reads:
  (a) # net-repressed = ρ<0 **AND BH-FDR q<0.05** (FDR per predictor on the partial p); (b) positive-vs-negative
  **Mann-Whitney** (unpaired, across gene sets);
  (c) **paired per-gene test:** for each gene scored under *both* predictor and base, take the difference
  `d_g = ρ_pred(g) − ρ_base(g)`; **Wilcoxon signed-rank** on {`d_g`} (alternative pred<base, i.e. predictor
  more repressive); **bootstrap** = resample genes 2,000× → 95% percentile CI of median `d_g`; **\*** = CI
  entirely below 0 ⇒ predictor robustly beats base *gene-by-gene*. (Paired ⇒ more powerful than (b); run on
  positives, the orthogonal TSG set, and all genes.)
  *(Supersedes the earlier no-proliferation run; proliferation is canonical here.)*

| construction | gw net-FDR | med ρ pos | med ρ neg | pos<neg MWU | paired vs base: pos / TSG / all |
|--------------|-----------|-----------|-----------|-------------|----------------------------------|
| base | 444 | −0.161 | −0.004 | 1.1e-05 | — |
| val_a0 | 444 | −0.178 | −0.004 | 6.5e-06 | p=0.051 / ns / Δ≈0 → **inert** |
| val_a1 | 412 | −0.081 | −0.004 | 6.6e-04 | +0.038 (p=0.99) → **harmful** |
| ts_a0 | 304 | −0.163 | −0.005 | 3.5e-04 | p=0.64 / ns / ns → **inert** |
| ts_a1 | 295 | −0.088 | −0.005 | 6.2e-03 | worse → **harmful** |
| **acq** | **460** | **−0.182** | −0.002 | **5.0e-06** | **8.2e-4\* / 4.7e-4\* / 1.3e-14\*** |
| **acq_ev** | 459 | **−0.193** | −0.002 | **3.7e-06** | **7e-4\* / 2e-4\* / 2.9e-13\*** |
| acq_ts | 309 | −0.185 | +0.001 | 1.1e-04 | 0.24 / 0.089 / 0.041 → not robust |

paired median Δρ (acq−base): positives −0.018 [−0.026,−0.007]; TSG −0.004 [−0.023,−0.0001]; all −0.0001 (tiny, ultra-consistent).

- **Verdict:** weighting (evidence & TargetScan) **inert**; promiscuity **harmful**; **acq / acq_ev the only robust wins**; acq_ev≈acq ⇒ the acquired transform is the lever, weight adds nothing.
- **Limitation:** sees only *modulated* repression (needs cross-sample variance), **blind to constitutive repression** (low-variance-high-mean acquired arms invisible — burden axis untested); low sensitivity (miRNA <2-fold vs transcriptional variance); mRNA not protein; positive n=18 (orthogonal TSG + genome-wide back it); ts coverage 320 arms/777 genes.

## TEST 2 — Δρ = ρ_tumor − ρ_GTEx (the change)

- **Statistic:** per-gene Δρ; **NO per-gene FDR on Δρ.** Significance = **Mann-Whitney** (positive Δρ < negative Δρ). Confounds = symmetric composition+proliferation (§7). (Each cohort ρ has its own partial p; the difference is tested at the gene-set distribution level.)

| predictor | med ρ_tum (pos) | med ρ_GTEx (pos) | med Δρ (pos) | med Δρ (neg) | pos<neg MWU |
|-----------|-----------------|------------------|--------------|--------------|-------------|
| base | −0.097 | −0.016 | −0.082 | −0.007 | 0.0028 |
| force_val_a0 | −0.102 | −0.017 | −0.086 | −0.007 | 0.0021 |
| acq | −0.124 | −0.028 | −0.114 | −0.011 | 0.0022 |

- Validated positive control fires (positives separate from miRNA-poor, p≈0.002–0.003); force ties base; acq biggest effect.
- **Limitation:** GTEx-side `max(a_gtex−h,0)` is ill-posed (rectify healthy vs own baseline → near-degenerate) so Δρ double-counts the anchor for acq → **ρ_tumor (Test 1) is the cleaner metric for acq**; no per-gene FDR; GTEx n smaller (less power on ρ_GTEx); constitutive-repression blindness inherited.

---

## Verdict

`acq = Σ max(a−h,0)` is the one construction that beats summed abundance (robust, proliferation-survived, modest). Edge weighting and promiscuity do not. **Open gap:** coupling cannot see constitutive (low-variance) acquired pressure — the **burden** axis (mean acquired pressure, no correlation) is not yet tested against the controls.
