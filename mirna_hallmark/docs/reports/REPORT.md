# miRNA × Hallmark — detailed results report (TCGA-BRCA)

Run: orchestrated 2026-06-08; **pressure spine refreshed 2026-06-17**
(`softmax_z_logrpm` + `evidence_mass` M0; M7 hybrid `softmax_z` + `combined_mass`).
See `output/run_manifest.json` and `output/logs/downstream_refresh_20260617.log`.
Cohort: TCGA-BRCA primary tumors; participant-level (12-char) keys; multi-vial
collapse by mean (expression) / median (CNV). All inferential tests are FDR-corrected
(Benjamini-Hochberg) and, where marked, confounder-adjusted by **partial Spearman**
on CPE (purity) + Thornsson HRD. Methods: [`METHODS.md`](../METHODS.md). Tracked
claims: [`DISCOVERY_REGISTRY.md`](../DISCOVERY_REGISTRY.md).

Scope of evidence: 50 Hallmark sets, ~4,384-gene universe; 913,778
`(miRNA, gene, hallmark_set)` edges (14,147 high-evidence); per-gene pressure
computed for 1,324 target genes × 1,079 samples; miRNA-universe CNV over 1,419
mapped entities.

> **Reading guide.** Two effects recur and must be separated from biology:
> (1) a **global ploidy / aneuploidy** gradient (LumB/Her2/Basal high CN,
> Normal-like near-diploid) that inflates *all* Hallmark target-CN "spreads";
> (2) a **proliferation/immune co-expression** confound that turns some
> pressure↔expression correlations positive. Both are flagged inline.

---

## 1. Which programs are most miRNA-targeted? (high-evidence enrichment)

Per Hallmark, a hypergeometric test of whether its genes are over-represented
among high-evidence miRNA target genes (relative to the Hallmark universe).

| Hallmark | n genes | high-ev targets | frac | fold | q |
|----------|--------:|----------------:|-----:|-----:|--:|
| TGF_BETA_SIGNALING | 54 | 42 | 0.78 | **2.37** | 1.1e-10 |
| WNT_BETA_CATENIN_SIGNALING | 42 | 32 | 0.76 | **2.33** | 6.0e-08 |
| APOPTOSIS | 161 | 102 | 0.63 | 1.93 | 8.7e-15 |
| TNFA_SIGNALING_VIA_NFKB | 200 | 126 | 0.63 | 1.92 | 1.4e-17 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | 200 | 123 | 0.62 | 1.88 | 3.4e-16 |
| ANGIOGENESIS | 36 | 21 | 0.58 | 1.78 | 4.0e-03 |
| UV_RESPONSE_DN | 144 | 83 | 0.58 | 1.76 | 3.2e-09 |
| NOTCH_SIGNALING | 32 | 18 | 0.56 | 1.72 | 1.1e-02 |
| G2M_CHECKPOINT | 200 | 111 | 0.56 | 1.69 | 1.1e-10 |
| IL6_JAK_STAT3_SIGNALING | 87 | 48 | 0.55 | 1.68 | 5.1e-05 |
| PI3K_AKT_MTOR_SIGNALING | 105 | 57 | 0.54 | 1.66 | 1.7e-05 |
| **Least-targeted:** OXIDATIVE_PHOSPHORYLATION | 200 | 34 | 0.17 | 0.52 | 1.0 |
| BILE_ACID_METABOLISM | 112 | 23 | 0.21 | 0.63 | 1.0 |
| DNA_REPAIR | 150 | 34 | 0.23 | 0.69 | 1.0 |

**Observation.** The most miRNA-targeted programs are signal-transduction and
plasticity modules (TGF-β, WNT/β-catenin, NF-κB, EMT, Notch, PI3K-AKT-mTOR),
plus apoptosis and proliferation control. The least-targeted are
housekeeping/metabolic gene sets (OXPHOS, bile-acid, DNA repair).

**Interpretation.** This is exactly the layout expected if miRNAs preferentially
buffer **signaling and cell-state-decision** networks rather than constitutive
metabolic machinery. It also reflects curation/study bias in miRTarBase (heavily
studied oncogenic pathways), so the *ranking among signaling programs* is more
trustworthy than the absolute fold values.

**Literature.** Concordant with the canonical view that miRNAs are enriched
regulators of TGF-β/EMT, WNT, and NF-κB signaling, and with miRNA control of
apoptosis; OXPHOS/DNA-repair under-representation matches the weaker literature
on direct miRNA targeting of core metabolic enzymes.

---

## 2. Distinct miRNAs in context (which miRNAs drive which programs)

### 2.1 Broadest high-evidence regulators (most Hallmark sets touched)

| miRNA | # Hallmark sets | # target genes |
|-------|---------------:|---------------:|
| hsa-miR-27a-3p | 45 | 58 |
| hsa-miR-125b-5p | 45 | 63 |
| hsa-miR-21-5p | 44 | 75 |
| hsa-miR-34a-5p | 44 | 59 |
| hsa-miR-155-5p | 43 | 55 |
| hsa-miR-214-3p | 43 | 35 |
| hsa-miR-124-3p | 42 | 70 |
| hsa-miR-16-5p | 42 | 43 |
| hsa-miR-29a-3p | 41 | 55 |
| hsa-let-7a-5p | 41 | 35 |

These are the textbook breast-cancer "hub" miRNAs: **miR-21** and **miR-155**
(oncomiRs), **miR-34a** and **let-7** (p53-linked / tumor-suppressive),
**miR-125b**, **miR-16**, **miR-29a**, **miR-124**. Their breadth reflects both
true pleiotropy and the depth of their experimental validation.

### 2.2 Program-specific dominant regulators (high-evidence, by # target genes)

| Hallmark | dominant high-evidence miRNAs (n target genes) |
|----------|-----------------------------------------------|
| TGF_BETA_SIGNALING | miR-17-5p (5), miR-155-5p (4), miR-181a-5p (4), miR-27a-3p (4) |
| WNT_BETA_CATENIN_SIGNALING | **miR-34a-5p (11)**, miR-145-5p (5), miR-16-5p (5), miR-30c-5p (4) |
| EPITHELIAL_MESENCHYMAL_TRANSITION | **miR-29b/c/a-3p (23/17/17)**, miR-767-5p (11), miR-21-5p (10) |
| TNFA_SIGNALING_VIA_NFKB | miR-21-5p (12), miR-155-5p (11), miR-199a-5p (11), miR-34a-5p (11) |

**Observation / literature.** The program-specific leaders are mechanistically
"correct": **miR-34a** is an established WNT/β-catenin repressor (and p53
effector); the **miR-29 family** is the canonical EMT / ECM-remodeling brake
(targeting collagens and DNA-methyltransferases); **miR-21 + miR-155** anchor
the NF-κB/inflammatory axis; **miR-17-5p / miR-181a** are recognized TGF-β
modulators. This convergence of an unbiased target count onto the literature's
program-specific miRNAs is a strong internal validity check.

---

## 3. Does miRNA pressure suppress its target programs? (AGO-gated anti-correlation)

Per Hallmark, Spearman(AGO-gated Hallmark pressure, mean target-gene log2-TPM)
across samples; **negative = repression** as expected. `spearman_*` is the raw
coefficient; `partial_*` adjusts CPE+HRD.

**24 / 50 Hallmarks** show a significant **negative** gated anti-correlation
(q<0.05; raw Spearman). After CPE+HRD partial adjustment, **19 / 50** remain
negative at q<0.05. Strongest raw negatives:

| Hallmark | n | rho | q | partial_rho | partial_p |
|----------|--:|----:|--:|------------:|----------:|
| BILE_ACID_METABOLISM | 1077 | **-0.363** | 3.3e-33 | -0.266 | 6.0e-18 |
| ESTROGEN_RESPONSE_EARLY | 1077 | -0.265 | 2.0e-17 | -0.143 | 4.6e-06 |
| ESTROGEN_RESPONSE_LATE | 1077 | -0.242 | 1.1e-14 | -0.120 | 1.3e-04 |
| APICAL_SURFACE | 1077 | -0.237 | 3.3e-14 | -0.197 | 2.2e-10 |
| ADIPOGENESIS | 1077 | -0.226 | 5.9e-13 | -0.140 | 7.7e-06 |
| MYOGENESIS | 1077 | -0.215 | 7.6e-12 | -0.161 | 2.7e-07 |
| PEROXISOME | 1077 | -0.207 | 4.8e-11 | -0.135 | 1.6e-05 |
| ANGIOGENESIS | 1077 | -0.185 | 4.8e-09 | -0.191 | 9.0e-10 |
| ANDROGEN_RESPONSE | 1077 | -0.168 | 1.1e-07 | -0.117 | 1.9e-04 |
| OXIDATIVE_PHOSPHORYLATION | 1077 | -0.157 | 7.5e-07 | -0.086 | 6.0e-03 |
| COAGULATION | 1077 | -0.156 | 8.3e-07 | -0.167 | 9.0e-08 |

**Unexpected positives** (pressure↑ with expression↑), and what happens after adjustment:

| Hallmark | rho | partial_rho | comment |
|----------|----:|------------:|---------|
| G2M_CHECKPOINT | +0.254 | +0.152 | persists — proliferation co-expression |
| E2F_TARGETS | +0.191 | +0.107 | persists — proliferation co-expression |
| ALLOGRAFT_REJECTION | +0.203 | +0.034 | **collapses** after purity/PAM50 adjustment |
| IL6_JAK_STAT3_SIGNALING | +0.135 | -0.002 | **collapses** |
| INTERFERON_GAMMA_RESPONSE | +0.123 | -0.012 | **collapses** |

**Observations.**
- Differentiation / metabolic-identity programs (bile-acid, estrogen/androgen
  response, adipogenesis, myogenesis, peroxisome) show the cleanest repression
  signal that **survives confounder adjustment**.
- The positive correlations split into two kinds: **proliferation** programs
  (G2M, E2F) stay positive even after adjustment, while **immune** programs
  (allograft rejection, IL6/JAK/STAT3, IFN-γ) are **artifacts of tumor purity /
  immune infiltration** — they vanish under CPE+PAM50 adjustment.

**Interpretation.** Positive coupling for proliferation reflects co-regulation,
not de-repression: in fast-cycling tumors both proliferation-associated miRNAs
and proliferation genes rise together (a shared cell-cycle program), so the
per-sample pressure score cannot be read as causal repression there. The
adjusted negatives are the credible repression signals. This is the single most
important caveat for downstream use of the pressure score.

**AGO gating effect.** Gating with the bounded RISC capacity (gate ∈
[0.51, 0.93] across the cohort) **rescales but rarely reorders** Hallmarks:
gated coefficients track ungated almost 1:1. Treat the gate as a documented
sensitivity layer, not a mechanism (registry MH-6).

---

## 4. Stratum differences

### 4.1 Anti-correlation by PAM50 (interaction view)

Pressure↔expression coupling is **strongly subtype-dependent**, and the
cohort-wide negatives are disproportionately driven by the small **Normal-like**
group. Most-variable Hallmarks (rho range across PAM50):

| Hallmark | Basal | Her2 | LumA | LumB | Normal | range |
|----------|------:|-----:|-----:|-----:|-------:|------:|
| XENOBIOTIC_METABOLISM | -0.03 | +0.12 | +0.10 | +0.12 | **-0.49** | 0.61 |
| ALLOGRAFT_REJECTION | +0.16 | +0.39 | +0.12 | +0.19 | -0.13 | 0.53 |
| IL6_JAK_STAT3_SIGNALING | +0.04 | +0.32 | +0.09 | +0.15 | -0.19 | 0.51 |
| BILE_ACID_METABOLISM | -0.19 | -0.08 | -0.20 | -0.06 | **-0.56** | 0.50 |
| ESTROGEN_RESPONSE_EARLY | -0.03 | -0.21 | -0.01 | +0.05 | **-0.43** | 0.48 |

Example — ESTROGEN_RESPONSE_EARLY by subtype: LumA rho≈-0.01 (q=0.95),
LumB +0.05, Her2 -0.21, Basal -0.03, **Normal-like -0.43 (partial -0.56)**.

**Observation / caveat.** The cohort-level estrogen-response anti-correlation is
**not** a luminal phenomenon — within LumA/LumB it is ~0. It is carried by
Normal-like tumors (n=36; uncorrected p only), which also show the strongest
negatives for xenobiotic/bile-acid/apical programs. Normal-like BRCA samples are
known to be enriched for normal-tissue / stromal admixture, so this likely
reflects normal-epithelium miRNA–program relationships leaking into the tumor
cohort rather than a tumor-intrinsic luminal effect. **Do not** generalize the
pooled negatives to luminal tumors.

### 4.2 Target-gene copy number by PAM50

Mean target-gene CN is highest in LumB/Her2/Basal and near-diploid in
Normal-like, **uniformly across Hallmarks** (e.g. cholesterol 3.39 / WNT 3.30 /
estrogen-early 3.25 in LumB vs ~2.33 in Normal-like). The per-Hallmark "spread"
is therefore dominated by the genome-wide aneuploidy gradient, not by
program-specific dosage. Use these tables for stratum CN context, not as
evidence of Hallmark-selective amplification.

### 4.3 Regulatory-miRNA expression by stratum

**Basal** tumors express the highest levels of Hallmark-targeting miRNAs across
many programs (bile-acid, IFN-α response, apical-surface, peroxisome,
adipogenesis), consistent with a globally more active miRNA compartment in
basal-like disease.

- **Immune (Thornsson C1–C6):** modest spread (≤0.42 log2-RPM); C2 tends highest.
- **Stage (I–IV):** **minimal** variation (≤0.24) — Hallmark-directed miRNA
  expression is essentially **stage-independent**, i.e. a lineage/subtype feature
  rather than a progression feature.
- **TNBC (Lehmann-4):** **BL1** (basal-like-1, proliferative) highest for IFN-α
  / apical / bile-acid / peroxisome regulatory miRNAs; **LAR** (luminal-AR)
  consistently lowest — matching LAR's more luminal, less basal biology.

---

## 5. miRNA-universe copy number and CN→expression concordance

Across the **whole miRNA universe** (1,419 mapped locus/arm/family entities),
locus copy number propagates to mature-arm expression for a substantial fraction
of arms: **238 / 756** testable arms show a significant **positive** Spearman
(q<0.05); only **9** are significantly negative; median rho ≈ 0.06 (so the
*typical* arm is weakly coupled, but a third are clearly dosage-driven).

Top dosage-coupled arms:

| arm | n | rho | q | partial_rho |
|-----|--:|----:|--:|------------:|
| hsa-miR-151a-3p | 1043 | 0.48 | 7e-59 | 0.41 |
| hsa-miR-151a-5p | 1043 | 0.45 | 5e-51 | 0.35 |
| hsa-miR-30d-5p | 1043 | 0.41 | 8e-42 | 0.38 |
| hsa-miR-30b-5p | 1042 | 0.39 | 1e-36 | 0.38 |
| hsa-miR-320a | 1043 | 0.38 | 5e-35 | 0.37 |

**Observation / interpretation.** The most dosage-driven arm, **miR-151a**, sits
at **8q24.3** — the same amplicon neighborhood as MYC and AGO2 — so its tight
CN→expression coupling is the expected signature of recurrent 8q gain in breast
cancer. miR-30b/d (8q24.2) behave likewise. This qualifies the pressure model:
for these arms, miRNA abundance is partly a *structural CNV readout*, so their
contribution to "pressure" co-varies with aneuploidy and should be interpreted
in that light.

---

## 6. AGO / RISC capacity (rate-limiting layer)

- Per-sample gate spans **[0.512, 0.928]** (cohort capacity z-scores mapped
  through the bounded logistic).
- **AGO2 copy number** is broadly elevated in aggressive subtypes: LumB 5.52,
  Basal 4.90, Her2 4.34, LumA 3.70, **Normal-like 2.78** (near-diploid);
  % gain 80% (Basal) / 80% (LumB) / 31% (Normal). This is the expected footprint
  of **8q24 amplification** (AGO2/EIF2C2).
- **AGO2 expression** is highest in Basal (mean log2-TPM 2.45 vs ~2.19–2.27
  elsewhere; KW q≈6e-32) — but note LumB has the *highest* AGO2 CN yet only
  mid-range expression, i.e. AGO2 CN and mRNA are **not** tightly coupled at the
  cohort level.

**Interpretation.** RISC catalytic capacity (AGO2) is amplified in the
genomically unstable subtypes, which is *why* a capacity gate is biologically
motivated. But because AGO2 CN→mRNA coupling is loose and the bounded gate
barely reorders results (§3), the present data support AGO availability as a
**plausible modulator worth modeling**, not as a demonstrated rate-limiter in
this cohort. This is the honest boundary of the current evidence.

---

## 7. Comparison vs literature knowledge (summary)

| Finding (this study) | Literature status |
|----------------------|-------------------|
| TGF-β, WNT, NF-κB, EMT, Notch are the most miRNA-targeted Hallmarks | **Concordant** — canonical miRNA-regulated signaling/plasticity networks |
| miR-21, miR-155, miR-34a, let-7, miR-29, miR-125b are the hub regulators | **Concordant** — established BRCA oncomiRs / tumor-suppressor miRNAs |
| miR-29 family dominates EMT; miR-34a dominates WNT; miR-21/155 dominate NF-κB | **Concordant** — matches mechanism-specific literature |
| Differentiation/metabolic programs show clean pressure→repression | **Plausible/novel** — consistent with miRNA buffering of differentiation, not previously framed this way in BRCA |
| Proliferation (G2M/E2F) shows positive pressure–expression coupling | **Expected confound** — co-regulation within the cell-cycle program, not de-repression |
| Immune-program positives vanish after purity adjustment | **Expected confound** — immune signal tracks infiltration/purity |
| miR-151a/30b/30d CN→expression coupling at 8q24 | **Concordant** — recurrent 8q gain in breast cancer |
| AGO2 amplified across aggressive subtypes (8q24) | **Concordant** — known AGO2/8q24 amplification |

---

## 7b. Robustness of the basal hub (proliferation + curation-bias controls)

Module: `robustness_checks.py` (grant Aims 1–2). Basal n=197. Two pre-registered
objections, both tested on the real cohort.

**Aim 1 — proliferation (with over-adjustment control).** The hub (miR-17~92 /
miR-106b~25) is MYC/E2F-driven and Basal is the most proliferative subtype, so we
added a proliferation covariate to the CPE+HRD partial-Spearman ladder. Because the
first proxy (E2F+G2M metagene) shares the hub's upstream driver, we ran **three
proxies**: E2F+G2M, **MKI67 alone**, and an **E2F/MYC-stripped** G2M+MITOTIC_SPINDLE
metagene (271 genes). The negative coupling **strengthens under all three** —
proliferation is a suppressor confound, not over-adjustment:

| program (Basal) | CPE+HRD | +E2F/G2M | +MKI67 | +ortho(no-E2F/MYC) |
|-----------------|--------:|---------:|-------:|-------------------:|
| APOPTOSIS | −0.16 | −0.29 | −0.24 | **−0.33** |
| P53_PATHWAY | −0.18 | −0.26 | −0.23 | **−0.30** |
| EMT | −0.19 | −0.24 | −0.22 | −0.26 |
| TNFA_NFKB | −0.16 | −0.26 | −0.24 | −0.25 |
| PI3K_AKT_MTOR | −0.03 | −0.23 | −0.19 | **−0.22** |
| INTERFERON_GAMMA | −0.06 | −0.22 | −0.19 | −0.15 |
| NOTCH | −0.10 | −0.11 | −0.11 | −0.13 (NS) |

Key-program survival: **8/8** (E2F/G2M), **6/8** (MKI67), **8/8** (ortho). Independence claim
leads on **MKI67-alone** (the ortho metagene is *largely* E2F/MYC-independent, not
fully — G2M genes remain partly E2F-influenced). **Magnitude is illustrative, not an
effect size:** part of the deepening is variance deflation (high target↔proliferation
co-variance, e.g. PI3K +0.89, leaves a thin residual that inflates |partial-r|), so
the values are proxy-dependent (IFN-γ −0.22→−0.15 is the tell). Report **sign +
FDR-survival + count**; treat magnitude as direction only.

**Sign structure** (`prolif_sign_structure.tsv`): at **program** level proliferation
is positively associated with *both* pressure (≈+0.1) and target expression
(+0.36…+0.89). At **gene** level (the headline claim) the same sign holds but weaker:
`rho(prolif, target_expr)` ortho/MKI67 = PTEN +0.54/+0.24, TGFBR2 +0.37/+0.07, CDKN1A
+0.27/+0.08, VIM +0.24/−0.02; hub-arm side +0.16…+0.38. PTEN shows the clean
suppressor sign individually; p21/TGFBR2/VIM are weak-positive (VIM ~0 under MKI67) —
single TSGs track proliferation less than a growth-program score. Per-gene routes
(`hub_route_partial_corr.tsv`, ortho proxy): **PTEN** (all 3 arms), **CDKN1A**
(miR-106b/17), **TGFBR2** (both), **VIM** survive (8/14; 9/14 under E2F/G2M & MKI67;
13/14 cohort-wide); **IRF1** (miR-23a-3p p≈0.4–0.6) and **BCL2L11** fail under every
proxy.

**Aim 2 — curation bias.** Per-miRNA basal load recomputed under three weightings
of the same high-evidence edges (`mirna_basal_load_by_weighting.tsv`):

| miRNA | rank evidence | rank binary | rank TargetScan |
|-------|:-------------:|:-----------:|:---------------:|
| miR-17-5p | 1 | 1 | 1 |
| miR-20a-5p | 3 | 3 | 2 |
| miR-106b-5p | 9 | 10 | 5 |
| miR-23a-3p | 15 | 15 | 13 |

Load rank-correlation: evidence-vs-TargetScan ρ=0.65; binary≈evidence ρ=1.00. The
miR-17~92/106b~25 family tops all three, so the hub is **not** a study-count
artifact. Per-target (`per_target_top_regulators.tsv`), TargetScan reassigns the
*specific arm* (PTEN→miR-19a/b-3p + miR-29; **IRF1→miR-130b-3p/miR-301a/b**, not
miR-23a-3p/miR-106b-5p) — so cluster-level robustness is strong but single-edge
nominations, especially IRF1, are predictor-dependent.

**Read.** The proliferation- and weighting-robust core is a basal
miR-17~92/106b~25 program repressing PTEN / p21 / TGFBR2 / VIM. The IRF1
immune-priming route is exploratory.

---

## 8. Limitations

1. **Pressure is associative.** It is an evidence-weighted miRNA-abundance score,
   not a measured repression rate; positive proliferation/immune couplings show
   it must be read with confounder adjustment.
2. **miRTarBase curation bias** inflates well-studied oncogenic pathways; the
   *relative ranking among signaling programs* is more reliable than absolute folds.
3. **Normal-like (n=36)** drives several pooled anti-correlations and carries
   normal-tissue admixture — pooled negatives should not be generalized to
   luminal tumors.
4. **Target-CN spreads** are dominated by genome-wide aneuploidy, not
   Hallmark-selective dosage.
5. **AGO gate** is a bounded sensitivity layer; current data motivate but do not
   prove a rate-limiting role.
6. **Multi-set genes** are counted in each Hallmark context (dedup on (miRNA, gene)
   for cohort totals), so per-Hallmark results are not mutually independent.
   Enrichment **ranking** (item 2) is unaffected; cross-Hallmark significance
   counts should be read with this overlap in mind.

---

## Appendix — TargetScan orphan sensitivity (MH-23/24/25)

**Question:** Do sequence-predicted miRNA→target pairs *excluded* from miRTarBase
high-evidence show repression coupling independently?

**Program level:** Orphan TS-only pressure **inflated** neg-sig counts under the
prior softmax_z+degree panel (Basal 38/50 vs 24/50 under miRTar HE-only routes).
Under the **adopted logrpm+evidence_mass spine**, orphan and full M0 program
coupling both reach **42/50 Basal** (8/8 key) — no additive Basal key program
beyond the spine. Luminal breadth patterns from the prior run should not be
re-read as orphan-specific under the new weighting.

**Gene level (hub aggregation):** Within Basal, orphan arms cluster on **CDKN1A**
(8 CN-surviving arms; miR-301/130/454/499 + let-7/98), **PTEN** (miR-148b), and
**TGFBR2** (miR-520d). Seed-family coupling heatmaps show coherent blocks. IRF1/BIM
orphan routes remain negative.

**Subtype intersection (MH-25):** the curated HE hub is **constitutive** (4 active
miRNAs, 3/4 pan-subtype, Jaccard ≥0.75, zero Basal-only); the orphan TS layer carries
the **subtype-private** structure (52 active miRNAs; 5 Basal-only, 5 LumB-only, 3
pan-subtype; Jaccard 0.27–0.47, Basal most divergent). Orphan subtype specificity is
private low-evidence arms, not a coherent shared program.

**Combining HE + TS (recommendation):** do not flat-union; build an **agreement-weighted
/ confidence-tiered** set — Tier 1 HE (canonical), Tier 2 HE∩TS (cross-validated,
up-weighted, testable for deeper Basal coupling), Tier 3 TS-only orphans (separate,
exploratory). Inputs (`ts_weight`, evidence weights) already exist. See BIOLOGICAL_SYNTHESIS §10.5.

**miRTarBase-HE curation-gap coupling queue (MH-26):** orphans yield a *queue*, not a
pressure set — HE membership used as a **negative** filter on TS predictions:
`(TS-predicted) AND (zero miRTar HE studies) AND (CN+prolif-adjusted coupling)`. 38/62
orphan hub relations are absent from miRTarBase HE; 19 couple. Strongest sit outside
Basal (miR-183→FOXO1 ρ −0.48; miR-124→VIM ρ −0.46; miR-383→IRF1 ρ −0.42 Her2) plus the
Basal p21/PTEN block. **CORRECTION (all 19 edges PubMed-triaged):** miRTar-HE absence is
a *curation gap, not literature novelty* — **14/19 are luciferase-validated elsewhere**
(FOXO1, p21, TGFBR2, IRF1, PTEN, VIM, BIM), 1 is likely indirect (miR-124→VIM), and **only
4 are genuine direct-edge gaps** (miR-454→CDKN1A, miR-301b→CDKN1A, miR-24/miR-519d→IRF1).
The queue is for **breast-context confirmation** of under-curated edges, not new biology.
Output: `mirtar_gap_coupling_queue.tsv` (`literature_status`).

**Outputs:** `output/robustness/targetscan_orphan/`; notebook
`notebooks/inspect_robustness_tables.ipynb` §14.

---

*Regenerate inputs: `.venv/bin/python3 -m mirna_hallmark.run_all`. Source tables
under `output/{hallmark_interaction,stratum_characterization,mirna_locus_cnv,ago_gate}/`.*
