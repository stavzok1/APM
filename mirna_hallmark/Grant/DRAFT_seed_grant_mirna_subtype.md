# HUCMED Seed Grant 2026 — Draft (Track A: Cancer)

> **Seed grant, not a full proposal.** This is scoped as a 12-month, USD 25,000
> proof-of-concept (POC) that turns preliminary computational results into a
> subtype-resolved miRNA biomarker / target shortlist and the preliminary data
> needed to compete for ISF / NIH / ERC. Format target (CFP §6): single PDF,
> **3 pages max incl. references & figures**, 11 pt, 1.5 spacing, ½″ margins;
> a concise budget justification and a 2-page CV per PI are appended.
> Bracketed `[PI INPUT]` items are deliberately left for the PIs to complete.

---

## Working title

**Subtype-resolved microRNA pressure on cancer Hallmark programs in breast
cancer: a proof-of-concept for stratified miRNA biomarkers and therapeutic
targets, focused on basal-like / triple-negative disease.**

---

## Scientific Abstract  *(≤250 words)*

Breast cancer is not one disease: its PAM50 subtypes (Luminal A, Luminal B,
HER2-enriched, basal-like / triple-negative, TNBC) differ in biology, prognosis
and treatment. microRNAs (miRNAs) are the dominant post-transcriptional layer
tuning oncogenic and tumor-suppressor programs, yet are still studied mostly
**cohort-wide**, averaging over the heterogeneity that defines the disease. As
preliminary data we built an evidence-weighted, AGO/RISC-gated model of miRNA
**pressure** on MSigDB Hallmark programs across TCGA-BRCA, **resolved by PAM50
subtype**. Even across 8 core programs, miRNA repression of Hallmark targets is
markedly subtype-specific: basal-like tumors carry by far the deepest, broadest
pressure — the strongest anti-correlations on nearly every program (**Fig. 1**) —
while luminal tumors are comparatively shielded. One miRNA family
(miR-17∼92 / 106b∼25) applies **monotonically stronger** pressure from Luminal A
to Basal on tumor-suppressor / cell-state programs (**Fig. 2**) — a subtype-graded
gradient, not on/off. This positions subtype-resolved miRNA pressure as both a
candidate **biomarker axis** and a **therapeutic-target map**, most actionable in
basal / TNBC, where targeted options are scarcest. This
seed will (i) expand this 8-core POC to all **50** MSigDB Hallmark programs; (ii)
harden it against **tumor-purity** confounding; (iii) **replicate** it in an
independent cohort (METABRIC) and at the **protein** level (CPTAC); and (iv)
distill a **minimal, most-informative miRNA set per subtype** — prioritizing
basal / TNBC and resolved **per tumor** — as a deployable biomarker / target
shortlist. The result is a subtype-specific, validated miRNA shortlist with POC
stratification value, positioning the team for ISF / NIH / ERC funding.

---

## Introduction and Hypothesis

**State of the art and the clinical challenge.** Breast cancer management is
already subtype-driven, but its molecular **regulators** are usually profiled
across the whole cohort. miRNAs are the central post-transcriptional rheostat of
the cancer cell — each can tune dozens of genes — and they are attractive both as
**non-invasive biomarkers** (stable in tissue and circulation) and as
**therapeutic targets** (mimics / antagomiRs). The unmet need is acute in
**basal-like / TNBC**, the most aggressive subtype with the fewest targeted
options. The gap this seed addresses: there is no **program-level, evidence-graded,
subtype-resolved** map of *which* miRNAs throttle *which* cancer programs in
*which* tumors — the prerequisite for using the miRNA layer for stratification or
intervention. Cohort-averaged analyses actively obscure this, because miRNA
regulation turns out to be strongly subtype-graded (below).

**Preliminary observation (motivating data).** In our TCGA-BRCA model, the
coupling between miRNA pressure and target-gene repression is **not uniform across
subtypes**. We show the **8 core** Hallmark programs (apoptosis, p53, PI3K–AKT–mTOR,
TGF-β, EMT, Notch, TNFα/NF-κB, IFN-γ — the tumor-suppressor, cell-state, signaling
and immune axes central to TNBC); across them, basal-like tumors show the
strongest, broadest repression and luminal tumors the weakest (**Fig. 1**), and a
single high-evidence miRNA family reproduces this as a clean LumA→LumB→Her2→Basal
gradient on tumor-suppressor / cell-state programs (**Fig. 2**). This 8-core set
is deliberately the displayed POC: it concentrates the tumor-suppressor,
cell-state, signaling and immune axes most relevant to basal / TNBC biology, rather
than asking reviewers to evaluate a broad discovery atlas. Studying miRNAs
**without** subtype stratification therefore averages a real, ordered biological
signal toward the null — the core motivation for this proposal.

**Central hypothesis.** *The strength of miRNA-mediated repression of cancer
Hallmark programs is subtype-specific and most intense in basal-like / TNBC;
consequently a small, subtype-resolved set of miRNAs carries biomarker
information and nominates therapeutic targets that a cohort-averaged analysis
misses.* The seed tests whether this 8-core POC generalizes across all 50 Hallmark
programs, survives tumor-purity adjustment, replicates across external cohorts and
protein-level readouts, and compresses to a minimal, informative miRNA panel per
subtype.

**Figure 1 — subtype heterogeneity of miRNA pressure (preliminary).**

![Spider: miRNA pressure to repression by subtype, 8 key Hallmark programs](../output/figures/seed_spider_keyhallmarks.png)

*8 key Hallmark programs; one trace per PAM50 subtype. Radius = anti-correlation
(−ρ) of miRNA Hallmark pressure vs target expression (within-subtype partial
Spearman, purity/HRD/proliferation/copy-number adjusted); **further from the
centre = stronger repression**. Pressure is built on the high-confidence edge
universe (high-evidence miRTarBase ∩ TargetScan-predicted; 235 miRNAs) so the
aggregate reflects sequence-supported targeting rather than indiscriminate
summation. Basal (red) is the outer envelope on nearly every program; luminal
tumors sit closer to the centre. n: LumA 503 / LumB 258 / Basal 188 / Her2 92
(Normal-like set aside).*

---

## Computational Medicine Approach

This seed integrates an existing computational platform (preliminary data in
hand) with technical extensions, independent-cohort and protein-level
replication, and a translational distillation step. It is purely computational
and reuses public, de-identified data — TCGA-BRCA, extended here with the
METABRIC and CPTAC cohorts.

### Universes and definitions already in place (preliminary data)

| Universe | Size / definition |
|---|---|
| **Samples** | TCGA-BRCA primary tumors with paired RNA-seq + miRNA-seq (pool n = 1,077). **PAM50 subtype-analysis cohort n = 1,041**: LumA 503 / LumB 258 / Basal 188 / Her2 92 (Normal-like, n = 36, set aside). |
| **8-core Hallmark programs** | The displayed POC uses 8 MSigDB Hallmark programs: apoptosis, p53 pathway, PI3K–AKT–mTOR, TGF-β, EMT, Notch, TNFα/NF-κB and IFN-γ response. Their gene universe is the **deduplicated union** of these sets (**964 genes**; 934 carry ≥1 miRNA edge). The seed aim expands this from 8 core programs to the full **50-program** MSigDB Hallmark collection. |
| **miRTarBase edges** | Within the 8-core universe: **135,800** unique (miRNA, gene) interactions over **2,930 miRNAs** and **934** genes (166,573 miRNA–gene–program rows). |
| **High-evidence subset** | Within the 8-core universe: **2,646** unique edges (594 miRNAs, 498 genes): ≥1 *functional MTI* study **and** a low-throughput validation (reporter or protein). |
| **TargetScan** | ~**1.40 M** weighted-context++ predictions, used to re-weight the same edges by **sequence** (curation-independent control). |
| **High-confidence pressure universe** (Fig. 1) | Within the 8-core universe, high-evidence ∩ TargetScan-predicted edges: **1,079** edges over **197 miRNAs** / 310 genes — the sequence-supported set on which Fig. 1 pressure is computed. |

**miRTarBase evidence score** (edge weight, study counts):
`evidence_score = 3·reporter + 3·binding + 2·protein + 1·rna + 1·perturbation`.

**miRNA pressure** (per gene *g*, sample *s*):
`pressure(g,s) = Σ_m z_m(s) · log(1 + evidence_{m→g})`, summed over miRNAs *m*
targeting *g* (edges with `evidence_score ≥ 2`), where `z_m(s)` is the cohort
z-score of mature-arm miRNA expression. Higher pressure ⇒ more predicted
repression. A Hallmark's pressure is the mean over its member genes.

**AGO/RISC gate** (documented sensitivity layer; every result reported gated
**and** ungated). RISC capacity `capacity(s) = mean z of AGO1–4`; gate
`g(s) = gate_min + (1 − gate_min)·σ(k·(capacity(s) − c0))` with
`gate_min = 0.5, k = 1, c0 = 0`, so `g ∈ [0.5, 1]`; gated pressure = pressure × g(s).

**Basal load** (miRNA ranking metric) = `Σ log1p(evidence_score)` over a miRNA's
high-evidence basal edges. The **miR-17∼92 / 106b∼25** family carries **22.1%**
of total basal repressive load (**28.2%** under the sequence-based **TargetScan
weight-score**, `Σ |weighted context++|` over the same edges) — a compact,
curation-robust hub.

### Aims

**Aim 1 — Expand the 8-core POC to all 50 Hallmark programs.** The current
preliminary result is intentionally shown on the 8 core programs most relevant to
basal / TNBC biology. The seed expands the same subtype-resolved, high-confidence
pressure framework to the full **50-program** MSigDB Hallmark collection, with
tumor purity, gated/ungated pressure, copy-number and subtype-aware FDR. This is
the technical scale-up that turns the POC figure into a systematic miRNA–program
atlas. *Deliverable:* a 50-program × subtype miRNA-pressure atlas identifying
which program axes are basal-specific, shared across subtypes, or absent.

**Aim 2 — Hardening against tumor-purity confounding.** Tumor purity varies by
subtype and can inflate or mask expression-based coupling. Extend the confounder
model to include **purity (CPE)** alongside the existing covariates, reporting
each coupling raw and purity-adjusted. *Deliverable:* a purity-robust subtype map
in which the basal-maximal signal is shown to be biology, not admixture.

**Aim 3 — Independent replication (METABRIC) and protein-level validation
(CPTAC).** Test whether the subtype-graded pressure map reproduces in an
independent cohort — **METABRIC** (RNA + PAM50) — and, because miRNAs act
post-transcriptionally, whether the predicted repression is even clearer at the
**protein** level in **CPTAC-BRCA** (the regime where the claim should be
strongest). This also serves as a power/stability check on the smaller subtypes.
*Deliverable:* cross-cohort and protein-level confirmation that the basal-maximal
miRNA axis is real and generalizes, not a single-cohort or mRNA-only artifact.

**Aim 4 — Minimal, most-informative miRNA set per subtype (biomarker / target
shortlist).** *Motivation (preliminary):* indiscriminate whole-miRNA pressure
**buries** real regulation. On the full aggregate Her2 looks almost unpressured,
yet restricting to the high-confidence universe roughly **doubles** its coupling
(mean ρ ≈ −0.05 → −0.12) and the miR-17∼92 / 106b∼25 cluster alone reaches
ρ = −0.37 — the diluting non-cluster bulk is **near-zero noise, not anti-signal**.
A compact, informative miRNA set therefore *recovers* signal that averaging
discards, which is exactly what a deployable assay needs. From the atlas we select
the smallest miRNA set that best separates a subtype's regulatory state —
**prioritizing basal / TNBC** — using load-ranked, sequence-corroborated hubs and a
sparse selection step. *Deliverable:* a short, ranked per-subtype miRNA panel
proposed as a stratification **biomarker** and a therapeutic-**target** list; the
basal panel is the headline POC output.

**Aim 5 — Subtype and per-sample deep dive (within-subtype heterogeneity).**
*(a, preliminary)* Resolving coupling to individual subtypes already shows the
effect is a **gradient**, not a basal-only switch: the miR-17∼92 / 106b∼25 family
deepens monotonically from Luminal A to Basal on tumor-suppressor / cell-state
programs (**Fig. 2**, e.g. p53 pathway −0.21 → −0.31 → −0.44 → −0.58).
*(b, extension)* Move from subtype **averages** to **per-sample** pressure
profiles to discern **within-subtype heterogeneity** — e.g. which individual
basal / TNBC tumors actually carry the hub program versus those that do not —
since a subtype label alone is too coarse for patient-level stratification. This
per-tumor resolution is what makes the Aim 4 panel clinically meaningful.
*Deliverable:* per-subtype, family/gene-resolved tables plus per-sample basal
pressure profiles quantifying intra-subtype heterogeneity.

**Figure 2 — one miRNA family, subtype-graded pressure (preliminary).**

![Heatmap: miR-17~92/106b~25 cluster coupling across subtypes, 8 key programs](../output/figures/seed_family_subtype_heatmap.png)

*miR-17∼92 / 106b∼25 cluster pressure → repression coupling (partial Spearman ρ;
negative = repression) for the 8 key programs across PAM50 subtypes. The same
family applies progressively deeper pressure LumA → LumB → Her2 → Basal — the
LumA/LumB deep-dive difference in one image.*

**Rigor.** Every coupling is reported under the confounder ladder (adding purity
in Aim 2); hubs are reported under evidence, degree, **and** sequence (TargetScan)
weighting; the AGO gate is reported gated and ungated. The cluster's
pressure→repression coupling is additionally robust to the **aggregation operator**
(summed vs count-normalized mean) and to partialling out a **global
miRNA-expression factor** (Basal ρ −0.50 → −0.42; Her2 −0.37 → −0.38), confirming
it is a specific module effect rather than a global-abundance artifact — whereas
the broad multi-miRNA aggregate, as expected, attenuates under the same controls
(staying negative and Basal-maximal). All aims are correlational and stated as
falsifiable, prioritized predictions for downstream experimental follow-up.

---

## Collaborative Plan  **[PI INPUT — scaffold]**

The CFP prioritizes ≥2 PIs bridging computation and clinic.
- **PI A — Computational / regulatory genomics ([lab/name]).** Owns the
  miRNA–Hallmark platform and Aims 1–4 (atlas, purity hardening, minimal-panel
  selection). *Provides the candidate subtype panels.*
- **PI B — Clinical / translational breast oncology ([name/institution]).**
  Grounds the biomarker shortlist in clinical questions and access to annotated
  (TNBC) material for downstream validation. *Provides clinical relevance and the
  path to a stratification assay.*  **[PI INPUT — confirm collaborator]**
- **Synergy.** Computation **nominates** small, subtype-specific miRNA panels and
  the exact tumors they should matter in; the clinical partner **prioritizes**
  the most translationally useful panel and frames the validation study.

---

## Path to Scale

Positive POC results target:
- **ISF** — mechanistic grant on subtype-specific miRNA control of cancer
  Hallmark programs in breast cancer.
- **NIH R01 / DoD BCRP** — biomarker + therapeutic-development program on the
  basal / TNBC miRNA panel (stratification assay + antagomiR/mimic intervention),
  leveraging the TCGA POC + independent-cohort and experimental validation.
- **ERC (Starting / Consolidator)** — broader "subtype-resolved regulatable
  cancer cell state" program, **[PI INPUT — eligibility]**.

---

## Budget Justification (USD 25,000, shared)  **[PI INPUT — indicative]**

| Item | Indicative | Notes |
|---|---:|---|
| Computational personnel (analyst, part-time) | ~$18,000 | Aims 1–4: 8-to-50 Hallmark expansion, purity hardening, replication, minimal-panel selection |
| HPC / cloud compute + storage | ~$3,000 | cohort-scale re-runs across 50 programs × subtypes |
| Data access / independent-cohort fees | ~$2,000 | external replication cohort access **[PI INPUT]** |
| Dissemination (publication / presentation) | ~$2,000 | acknowledging HUCMED + Azrieli OMICs Center |

*Computation-only seed; adjust to the two labs' real personnel rates and
institutional overhead.*

---

## Reporting & Affiliation (CFP compliance)

- Single summary report at the end of the 12-month period.
- All outputs acknowledge **HUCMED** and the **Azrieli OMICs Center**.
- Submit via the official online form by **15 June 2026**.

---

## Left for PI completion

- PI identities and the clinical/translational collaborator (≥2 PIs required).
- 2-page CV per PI; final budget tied to real rates.
- CFP page-count clarification (whether budget + CVs count against 3 pages).
- Choice of independent replication cohort for the Path-to-Scale framing.
- Figures regenerate with `.venv/bin/python3 -m mirna_hallmark.make_seed_grant_figures`
  (`output/figures/seed_spider_keyhallmarks.{png,pdf}`,
  `output/figures/seed_family_subtype_heatmap.{png,pdf}`).
