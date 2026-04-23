# 06 — Automation plan and hypothesis ranking

Which hypotheses can the engine run unattended, which need supervision, which are design-heavy. Plus a ranked execution order.

---

## Automation tiers

| Tier | Description | Throughput | Human-in-the-loop |
|---|---|---|---|
| **T1** | Feature → feature test (regression, mediation, interaction) on tidy tables with existing L1/L2 features | minutes–hours per hypothesis | Interpretation + figure caption + writing only |
| **T2** | Requires clustering, imputation, sensitivity sweeps, or moderate new feature engineering | days per hypothesis | Design review + parameter selection, then auto-run |
| **T3** | Requires new data access, new pipeline (HLA typing, external cohort harmonization), or novel inferential design (serial mediation with external validation) | weeks per hypothesis | Heavy — project-scoping |

---

## Tier assignments (full catalog)

| Hypothesis | Tier |
|---|---|
| H1, H2, H3, H7, H8, H9, H13, H14, H15, H17, H19, H20, H22, H26, H27, H28, H29, H30, H31, H32, H33, H35, H37, H38, H39, H40, H41, H42, H43, H44, H47, H49, H51, H52, H53, H54, H55, H58, H59, H60, H62 | **T1** |
| H5, H10, H11, H21, H23, H25, H34, H36, H45, H46, H48, H50, H56, H57, H61, H63, H64, H65 | **T2** |
| H6, H16, H18 (external replication), H24 | **T3** |

About **41 / 65 are T1, 18 are T2, 4 are T3**. H46 moves from T3 to T2 because HLA typing is available on the full cohort. H59–H62 and H63–H65 are added with HLA and TNBC-focus layers.

---

## Ranking for execution order

Criterion = `priority × tier_weight × literature_support × expected_novelty`. Manually ranked; the list below is the execution order for the first 12 months.

### Wave 1 — Foundation (month 1–2)

Runs first because everything else depends on their features.

1. **H51** — 13-HiChIP calibration of subtype-consensus prior (T1; unblocks half the catalog)
2. **H53** — Deviation from prior as feature (T1; shape of residuals informs all downstream regressions)
3. **H19** — HLA-locus co-methylation block (T1; defines a composite the rest of Family G depends on)

### Wave 2 — Easy wins, high literature support (month 2–3)

4. **H41** — IRF1/IRF2 ratio (T1; PMC6761035)
5. **H38** — Bottleneck principle min-rule (T1; yields clean plot)
6. **H44** — Pathway coherence (T1; yields a new sample phenotype)
7. **H8** — R_g rank explains B2M fragility (T1; one-figure result)
8. **H40** — STAT1/STAT3 antagonist ratio (T1)

### Wave 3 — Core novel claims (month 3–5)

9. **H1** ✱ — CRML beyond CNV + methylation (T1; C1+C2 anchor)
10. **H9** ✱ — Conservation-weighted BCD (T1; C2 anchor)
11. **H7** — R_g × burden interaction (T1)
12. **H14** ✱ — Four-quadrant escape map (T1)
13. **H22** — STAT1 methylation cascade (T1)

### Wave 4 — Proteostatic layer (month 5–7)

14. **H4** ✱ — mRNA-competent/protein-deficient class (T1; RPPA discordance)
15. **H15** — Proteasome-primed/MHC-lost (T1)
16. **H13** — NK ligand − shedding vs activated-NK (T1)
17. **H10** — RPPA archetype clustering (T2)
18. **H47** — Peptide-pipeline coherence (T1)

### Wave 5 — Non-coding layer (month 6–8)

19. **H26** — miR-148a/152 → HLA-G (T1)
20. **H27** — miR-125a/b → B2M/MHC-I (T1)
21. **H29** — Data-driven miR master regulators (T1)
22. **H31** — Coordinated miR-pressure phenotype (T1)
23. **H32** — HCP5 → PD-L2 ceRNA (T1)
24. **H35** — lncRNA–gene TAD colocalization (T1)

### Wave 6 — Subtype and within-panel (month 7–9)

25. **H3** — Subtype coding vs cis damage (T1)
26. **H39** — Paralog compensation (T1)
27. **H42** — SOCS1 negative-feedback leakiness (T1)
28. **H43** — EP300/CREBBP coactivator rate-limit (T1)
29. **H49** — NLRC5 vs CIITA coordination (T1)
30. **H54** — Core vs opportunistic enhancer effect sizes (T1)

### Wave 7 — HLA per-allele layer (month 4–6, runs parallel to Wave 3–5)

31. **H46** ✱ — HLA per-allele ASE loss (regulatory vs genomic) (T2; HLA typing available)
32. **H63** ✱ — Regulatory ASE attributable to per-allele CRML + β (T2)
33. **H64** — HLA supertype × structural-integrity interaction (T2; needs ICB cohort)
34. **H65** — Per-allele immunoediting signature (T2; needs neoantigen calls per allele)

### Wave 8 — TNBC hot/cold lead case (month 6–10)

35. **H59** ✱ — Four-class hot/cold decomposition (T1 compute, T2 validation)
36. **H60** — Hot-discordant class mechanism (T1)
37. **H61** — Cold-primed class STING-responsive (T2 + wet-lab W4)
38. **H62** — TNBC sub-subtype × four-class distribution (T1)

### Wave 9 — Late-stage core claims (month 8–12)

39. **H12** ✱ — CD8 × JSI moderation (T1; large-N interaction, do it right)
40. **H16** — Serial mediation CRML→CXCL→CD8→CYT→survival (T3)
41. **H18** ✱ — IVI predicts OS (T2; Cox + external cohort)
42. **H56** ✱ — Dual-axis IFN intactness (T2; requires extended panel)

### Wave 10 — Platform (month 10–12)

43. **H17** — Mechanism fingerprint across DDR, EMT (T1–T2)

### Waves 11+ — Remaining T2/T3 + exploratory T1

Order opportunistically based on Aim 2 / Aim 3 needs. The catalog has ~20 more hypotheses; completing 12–15 of them is enough for a rich thesis.

---

## Batch-test design

Run Waves 1–3 and 4–5 as a batch. Structure:

```
for hypothesis in wave:
    features = load_or_build(hypothesis.required_features)   # cache-aware
    result = engine.fit(hypothesis.design)                    # unified interface
    write_report(hypothesis.id, result)                       # fixed template
```

Hypothesis design is a declarative spec:

```yaml
id: H1
y: APM_expr_residual
x: CRML_snv + CRML_sv
covariates: [CNV_gene, promoter_beta, enhancer_beta, purity, PAM50, IFN_gamma]
moderator: null
method: elastic_net_with_bootstrap
bootstraps: 2000
strata: PAM50
report_figure: Fig_CRML_APM
```

Writing these YAMLs for Waves 1–7 is a one-time cost of about a day. After that, re-running any wave is hours, not weeks.

---

## Ranking signals (how the engine auto-prioritizes)

When new hypotheses are generated (e.g., after a lit-search pass), rank them by:

| Signal | Weight | Notes |
|---|---|---|
| Literature support (PubMed hits with mechanism) | 0.25 | Proxy for reviewer acceptance |
| Lead-set membership | 0.20 | Forced to the top |
| Tier (T1 > T2 > T3) | 0.20 | Efficiency |
| Claim coverage (covers C1–C5 gap) | 0.15 | Fills thesis gaps |
| Novelty (not covered by priors) | 0.10 | Publication value |
| Orthogonality to already-run hypotheses | 0.10 | Avoid redundant tests |

A `rank_hypotheses.py` script can emit a priority list from the catalog YAML.

---

## Human-in-the-loop checkpoints

Required for:

- Freezing the L1 feature contract (one-time, month 1).
- Approving clustering choices (H10, H21).
- Choosing external cohort(s) (G4, month 7–8).
- HLA typing pipeline decision (H46 go/no-go, month 3).
- Publication-figure selection per wave (one half-day per wave).
- Committee-facing narrative synthesis (one week per aim paper).

---

## Negative-result policy

A **null Tier-1 hypothesis is reported in the thesis**, not discarded. Pre-registering the lead set means every null carries weight. This is a feature, not a bug: it forces honest interpretation.

For T2/T3 hypotheses that fail, document the failure mode (underpowered, feature-engineering failure, external cohort unavailable) in `research_plan/negative_results.md` (create on first null result).
