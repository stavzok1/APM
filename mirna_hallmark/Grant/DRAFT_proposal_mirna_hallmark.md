# HUCMED Seed Grant 2026 — Draft Proposal (Track A: Cancer)

> **Status:** Draft of the parts that can be authored from the existing
> `mirna_hallmark` subproject results. Sections needing human input are flagged
> **[NEEDS PI INPUT]**. Format target (from the CFP): single PDF, **3 pages max
> incl. references & figures**, 11 pt, 1.5 spacing, ½″ margins; budget
> justification + a 2-page CV per PI appended.
> **Resolve before compressing the body:** email the HUCMED program officer to
> confirm whether the budget justification and per-PI CVs count against the
> 3-page limit (the CFP says "the same PDF must additionally include," which is
> ambiguous). The answer changes how tightly the body must be cut.
> **Figure:** within 3 pages a single figure is likely the highest-value use of
> space — recommended panels: (a) basal multi-program pressure→repression coupling
> *before vs after proliferation adjustment* (showing it strengthens), and (b) the
> hub-miRNA rank under evidence/degree/TargetScan weighting (showing miR-17-5p #1
> in all three). Both are directly buildable from `output/robustness/` tables.
> **Robustness analyses (proliferation + curation-bias) are already run**
> (`mirna_hallmark/robustness_checks.py`); their numbers are folded in below.

---

## Working title

**A microRNA cluster that represses the tumor-suppressor and cell-state networks
of triple-negative breast cancer most strongly in basal-like tumors: a
regulatory-genomics proof of concept for miRNA-directed reprogramming.**

---

## Scientific Abstract  *(≤250 words)*

Triple-negative / basal-like breast cancer (TNBC) is the most aggressive breast
cancer subtype, with the fewest targeted options. Its biology is defined by loss
of tumor-suppressor/apoptotic control (PTEN, p21, BIM) and aberrant cell-state
programs — much of it through regulatable silencing, *without* genetic deletion.
microRNAs are the dominant post-transcriptional tuning layer, yet the field lacks a
program-level map of *which* miRNAs coordinate *which* programs in *which* tumors.
We built an evidence-weighted, AGO/RISC-gated model of miRNA **pressure** on all 50
MSigDB Hallmark programs across TCGA-BRCA, by **PAM50 subtype**. Basal-like tumors
carry the broadest miRNA pressure, coordinated by a compact hub — the
**miR-17∼92 / miR-106b∼25 cluster** — across tumor-suppressor and cell-state
programs (PTEN, p21, TGFBR2, VIM). The basal pressure→repression coupling
**survives proliferation adjustment** under E2F/G2M, MKI67-alone, and
E2F/MYC-stripped proxies (**8/8**, **6/8**, and **8/8** key programs negative and significant).
The hub is stable under study-count, degree, and sequence-based TargetScan
weighting (miR-17-5p ranks #1 in all three). The gene-resolved robust core is
**CDKN1A/TGFBR2/VIM** (with target copy-number adjustment); **PTEN** holds at
individual hub-arm level; apoptosis holds at *program* level but not at BIM when
tested with combined per-gene miRNA pressure. This
seed will (i) bound the per-patient miRNA-explained fraction and confirm it at the
**protein** (CPTAC) and **cohort** (METABRIC) level (also a stability check beyond
n≈197); (ii) pilot antagomiR de-repression of the tumor-suppressor axis in basal
cells. The computation yields a **falsifiable prediction the experiment can
refute**, positioning miRNA-directed reprogramming of TNBC for ISF / ERC / NIH
funding.

---

## Introduction and Hypothesis

**State of the art.** Basal-like / triple-negative breast cancer (TNBC) is the
most aggressive intrinsic subtype and the one with the fewest targeted therapies.
Its hallmark biology is the functional loss of tumor-suppressor and pro-apoptotic
control — **PTEN, FOXO1, CDKN1A (p21), BCL2L11 (BIM)** — together with aberrant
**EMT/cell-state and signaling** programs and, in many tumors, blunted
interferon/antigen-presentation output. Critically, much of this dysregulation
occurs **without genetic deletion**, i.e. through *regulatable*
transcriptional/post-transcriptional silencing that is, in principle,
reversible. microRNAs are the dominant post-transcriptional tuning layer of the
cancer cell, but the field still lacks a program-level, evidence-graded map of
*which* miRNAs throttle *which* oncogenic/suppressor programs in *which* tumors —
the prerequisite for treating the miRNA layer as a therapeutic dial.

**Our proof of concept (already generated).** Within the parent
regulatory-genomics platform we constructed an evidence-weighted,
AGO/RISC-availability–gated model of miRNA *pressure* on all 50 MSigDB Hallmark
gene programs across TCGA-BRCA (n ≈ 1,060 primary tumors), resolved to PAM50
subtypes (each subtype vs the complement of the other three; Normal-like
excluded). Three results motivate this proposal:
1. **The miRNA→repression coupling concentrates in Basal — measured within
   subtype.** We compute coupling **inside each PAM50 subtype** (partial Spearman,
   adjusted for purity/HRD/proliferation, mean program copy number, and — as a
   sensitivity check — **within-subtype z-scoring** of pressure and expression to
   control internal dynamic range). It is **densest and deepest in Basal** — negative
   and FDR-significant for **42/50 Hallmark programs (8/8 key)** — versus **LumB
   17/50 (4/8 key)**, **LumA 18/50 (4/8 key)**, and **Her2 0/50** (**Fig. 1A**). This
   is not a Basal-vs-rest contrast; LumB carries a weaker parallel signal.
   Proliferation is the largest basal-distinguishing axis and the Basal concentration
   **persists after it is removed** (point 2).
2. **A coherent, proliferation-robust tumor-suppressor / cell-state signature.**
   PTEN, CDKN1A (p21), TGFBR2 and VIM are the most consistently miRNA-pressured
   nodes in Basal — precisely the controls whose loss defines TNBC. The hub is the
   **dominant driver of the tumor-suppressor sub-axis**: **miR-17-5p / miR-20a-5p /
   miR-106b-5p** (and the cluster's miR-19 arm) are the top predicted regulators of
   **PTEN and p21**. **TGFBR2 and VIM are co-pressured cell-state nodes whose
   strongest *sequence*-predicted regulators lie outside the cluster** (miR-520/302,
   miR-30/138) — so we describe the cluster as the leading driver of the suppressor
   core, not the sole regulator of every node. Crucially, the
   basal pressure→target-repression coupling **survives — sign and FDR-significance
   intact — and does not weaken under proliferation adjustment** (8/8 key programs
   stay negative & significant). The same pattern holds with MKI67-alone and an
   E2F/MYC-stripped proliferation metagene. The coupling *deepens* on adjustment
   because proliferation is a *suppressor* confound — positively associated with
   **both** the hub miRNAs and the target genes, so it masks rather than creates the
   repression — but we report the stable quantities (sign, survival, FDR), not the
   deepened magnitudes, which are partly variance-deflation and proxy-dependent. Per gene,
   PTEN/p21/TGFBR2/VIM survive; **apoptosis is robust at the program level but does
   not localize to a single effector here** — neither individual hub miRNAs nor
   **combined per-gene miRNA pressure** (all high-evidence regulators summed) tracks
   BCL2L11/BIM after proliferation **and target copy-number** adjustment, so we do
   not over-claim a specific apoptotic node. At gene level with CN adjustment,
   **CDKN1A, TGFBR2 and VIM** remain negative and significant; **PTEN** holds via
   the hub arms (miR-106b/17/20a) though aggregate PTEN pressure is borderline once
   CN is included.
3. **A compact, druggable regulator hub that is not a curation artifact.** The
   **miR-17∼92 / miR-106b∼25 cluster** is the single most concentrated source of
   basal pressure: this one paralogous family (≈24 mature arms out of ~515 ranked
   miRNAs, <5% of the miRNA universe) carries **≈22% of the total evidence-weighted
   basal repressive load** (≈28% under sequence weighting) and supplies **4 of the
   top-10** miRNAs by load (miR-17-5p is #1) — that concentration is what makes the
   hub "compact."
   Because evidence-weighting could merely be upweighting famous oncomiRs, we
   recomputed the hub under (i) binary/degree and (ii) **sequence-based TargetScan
   weighted-context++** weighting: **miR-17-5p ranks #1 in all three** and the
   cluster family tops each (load rank-correlation evidence-vs-TargetScan ρ=0.65).
   Sequence weighting also reassigns the *specific* PTEN arm to **miR-19**, the
   established miR-17∼92 PTEN-repressing arm (Xiao/Mendell) — a recovery of known
   arm-level biology that *validates* the cluster claim (miR-19 is inside the
   cluster). A few inhibitors could therefore relieve several programs at once.
   The remaining basal pressure is distributed across **canonical basal/TNBC
   miRNAs** — a second polycistronic cluster (**miR-23∼27∼24**, ≈9%, which contains
   the IRF1-route candidate miR-23a-3p), plus miR-221/222, miR-9-5p, miR-135b-5p and
   miR-155-5p — all basal-enriched and recovered here as expected, an internal
   validation that the landscape is real. We lead on miR-17∼92/106b∼25 because it is
   both the single largest concentrated family (≈2.4× the next) **and the only
   top-load family whose pressure actually couples to target repression**: restricting
   the per-sample pressure to each family and re-testing within subtype, the cluster
   represses across subtypes but **most deeply in Basal** (median ρ −0.45 vs −0.19 in
   LumA), whereas miR-23∼27∼24 and miR-221/222 are basal-enriched in expression but
   show **no coherent coupling in any subtype**. They remain candidate secondary
   nodes only.

**What is and is not novel (positioning).** The *individual* edges are largely
established: miR-17∼92 → PTEN/p21 is textbook (Mendell/He). We therefore do
**not** claim these edges as discovery; we use the canonical map (miR-29→EMT,
miR-34a→WNT, miR-21/155→NF-κB) only as an internal-validity check. The contribution
is the **subtype-concentrated, multi-program coordination** — one compact cluster
simultaneously gating tumor-suppressor + cell-state + signaling programs, coupling
across subtypes but **most deeply in Basal**, beyond what proliferative state or
study-count bias explains.
A candidate basal **IRF1
immune-priming route** (miR-23a-3p/miR-106b-5p; edges validated by 3′UTR luciferase
in HCC/glioma) is the most exciting but **most fragile** claim: it does *not*
survive proliferation adjustment within Basal (p=0.43/0.62) and sequence prediction
favors miR-130b/301, so we carry it as an **exploratory, high-risk secondary** to
be tested, not assumed.

**Proliferation control.**
The hub *is* the proliferation machinery (miR-17∼92/106b∼25 are MYC/E2F-driven;
Basal is the most proliferative subtype), so elevated hub-miRNA **and** depressed
tumor suppressors could both be downstream of proliferation. In preliminary data,
the basal coupling remains negative and significant after adding proliferation.
This holds with **MKI67 alone** (a single histological proliferation marker) and an
E2F/MYC-stripped G2M+mitotic metagene, not only with the E2F+G2M score
(**Fig. 1B**). Aim 1 extends this to a per-patient explainability range.

**Central hypothesis.** *A compact miRNA cluster (miR-17∼92 / miR-106b∼25)
coordinately represses tumor-suppressor and cell-state programs (PTEN/p21/TGFBR2/VIM)
across breast-cancer subtypes but **most deeply in basal-like tumors** — beyond what
proliferative state or curation bias explains; inhibiting the cluster de-represses
several defining
programs together and reprograms the basal phenotype toward a more controllable
state. Whether the same cluster also de-represses immune-priming (IRF1/MHC-I) is a
higher-risk corollary to be tested, not assumed.*

---

## Computational Medicine Approach

The seed integrates the existing computational platform with targeted
experimental validation — the explicit HUCMED "computational + clinical/
translational" bridge.

**Preliminary data in hand.** A Hallmark-scoped high-evidence miRTarBase edge
graph; per-sample, evidence-weighted, AGO-gated pressure scores; per-gene PAM50
subtype resolution; partial-Spearman confounder control (**including target copy
number at gene level and mean program CN at Hallmark level**); and the robustness
analyses above — proliferation-adjusted coupling (sign and FDR survival under
**three** proliferation proxies), within-subtype dynamic-range sensitivity, and the
evidence/degree/TargetScan hub comparison (miR-17-5p ranks #1 in all three). These results are computed, documented and
directly extensible by the aims below.

**Aim 1 — Per-patient, bounded explainability of the basal hub (computational).**
Building on the completed proliferation control, model per-sample target expression
(PTEN, p21, TGFBR2, VIM and the broader basal node set) on the cluster miRNAs under
the **full confounder ladder** — purity (CPE), leukocyte fraction, copy-number,
methylation, **and proliferation** — and quantify, per PAM50 subtype, the miRNA
contribution as a **share of *total* target variance, decomposed against a stated
proliferation/CN baseline**, not as a residual partial correlation (a partial-r is
read off the variance left after the covariates and overstates the contribution).
The honest framing is explicit: proliferation typically owns the **majority** of a
target's variance and the cluster a **single-digit-percent, de-repressible
post-transcriptional layer on top of it** — that controllable fraction, not a large
share, is the POC quantity. Because the covariates are collinear, it is reported as a
**bounded sensitivity range across model specifications, not a point estimate**.
*Deliverable:* a per-patient, range-valued miRNA-*controllable* fraction for the
robust basal axis, with the proliferation baseline reported alongside. *This is
correlational and generates a falsifiable prediction that Aim 3 can refute.*

**Aim 2 — Protein + independent-cohort validation, and a power/stability check
(computational).** (a) Test whether the signature predicts **protein** suppression
in CPTAC-BRCA better than mRNA (miRNAs act post-transcriptionally — the regime where
the claim should be strongest). (b) Replicate the targeting map + basal cluster
signature in **METABRIC**. Because the within-Basal sample size (n≈197) is modest
against the five-plus-covariate ladder of Aim 1, METABRIC/CPTAC replication doubles
as a **power-and-stability check** — it directly answers "is n≈197 enough" by asking
whether the same axis reappears in independent cohorts. *Deliverable:* protein- and
cohort-independent confirmation that the proliferation-/weighting-robust basal axis
is real, generalizes, and is not a small-sample artifact.

**Aim 3 — Feasibility pilot for causality (clinical/experimental collaborator).**
Deliberately scoped to fit a seed budget, and **led by the robust axis**:
antagomiR/LNA inhibition of the cluster in **two basal lines (one basal-A, one
basal-B; luminal MCF-7 control)**, with the **core readout** being
tumor-suppressor/cell-state restoration — **PTEN/p21/TGFBR2 protein + an
apoptosis/viability assay**. Because sequence prediction maps individual targets to
specific cluster arms (e.g. PTEN → miR-19), we inhibit the **cluster broadly
(pooled/multi-arm LNAs)** rather than a single arm, so a target whose true arm
differs from miR-17-5p/20a/106b is still de-repressed; per-target arm-matching is a
secondary, mechanism-resolving readout. Surface **MHC-I by flow** is a *single,
exploratory long-shot* secondary readout: our own data say the cluster likely does
**not** drive IRF1 (route fails proliferation adjustment; TargetScan points to
miR-130b/301), so cluster antagomiR is **not predicted** to move MHC-I via that
route — we include it only as a cheap, hypothesis-generating add, not a secondary
endpoint. Organoids, full EMT/IFN-γ panels and additional lines are explicitly
deferred to the R01/ISF. *Deliverable:* a go/no-go causal signal on the core axis,
powering the larger program.

**Rigor.** Every claim is reported under the confounder ladder **including
proliferation**; the hub is reported under degree- and sequence-based weighting;
enrichment is reported as ranking; the AGO gate is a documented sensitivity layer
reported gated **and** ungated.
Aims 1–2 are explicitly correlational and stated as a falsifiable prediction for
Aim 3; the IRF1 route is carried as exploratory throughout.

---

## Path to Scale

Positive POC results target:
- **ISF (Israel Science Foundation)** — mechanistic grant on miRNA control of the
  basal tumor-suppressor/cell-state program (basic mechanism + patient-attribution
  model).
- **NIH R01 / DoD BCRP Breakthrough** — biomarker + therapeutic-development
  program on miRNA-directed reprogramming of TNBC (hub-miRNA inhibition to
  reactivate tumor-suppressor/apoptotic control, including immune sensitization as
  one readout), leveraging the TCGA/CPTAC/METABRIC + experimental package.
- **ERC (Starting/Consolidator)** — broader "regulatable cancer cell state"
  program if PI eligibility fits **[NEEDS PI INPUT]**.
- Translational step: a biomarker substudy of hub-miRNA-module status as a
  prognostic / response stratifier in TNBC, with the clinical collaborator.

---

## Collaborative Plan  **[NEEDS PI INPUT — scaffold below]**

The CFP prioritizes ≥2 PIs bridging basic research and clinical application.
- **PI A (Computational / regulatory genomics — [lab/name]).** Owns the
  miRNA–Hallmark platform, patient-level attribution (Aim 1), and protein/cohort
  validation (Aim 2). *Provides the candidate module and patient strata.*
- **PI B (Experimental / clinical breast oncology — [name/institution]).** Owns
  the antagomiR functional validation (Aim 3), TNBC cell-line/organoid models,
  and access to annotated TNBC tissue for downstream biomarker work.
  *Provides causal validation and clinical grounding.*
- **Synergy.** Computation **nominates** a small, testable miRNA hub and the
  exact tumors it should matter in; experiment **tests causality** and the
  protein/phenotype readouts computation cannot see; the loop refines the
  predictive model toward a TNBC stratification biomarker.

---

## Budget Justification (USD 25,000, shared)  **[NEEDS PI INPUT — indicative split]**

| Item | Indicative | Notes |
|------|-----------:|-------|
| Computational personnel (student/analyst, part-time) | ~$11,000 | Aim 1 (per-patient explainability ladder) + Aim 2 (CPTAC/METABRIC); core robustness already done |
| HPC / cloud compute + storage | ~$2,000 | cohort-scale re-runs |
| Wet-lab consumables — pilot (antagomiRs/LNAs, antibodies, flow reagents) | ~$9,000 | Aim 3: 2 lines + control, core readouts only |
| qPCR / apoptosis-viability assays | ~$3,000 | Aim 3 core readout |

*Scoped to the trimmed Aim 3 (two lines, core TSG + surface-MHC-I readout).
Adjust to the two labs' real personnel rates and institutional overhead.*

---

## Reporting & Affiliation (CFP compliance)

- Single summary report at the end of the 12-month period.
- All outputs acknowledge **HUCMED** and the **Azrieli OMICs Center**.
- Submit via the official online form by **15 June 2026**.

---

## Items I cannot author (require you)

- **PI identities & the clinical/experimental collaborator** (the CFP *requires*
  ≥2 PIs; the strongest framing needs a wet-lab/clinical breast-cancer PI).
- **2-page CV per PI.**
- **Final budget numbers** tied to real personnel costs / institutional rates.
- **CFP page-count clarification** — email the program officer (see Status note).
- **Figure (built):** `mirna_hallmark/output/figures/grant_figure_basal_heatmap.pdf`
  (one-column, B&W-safe). The affirmative thesis in one image — a **subtype × Hallmark
  heatmap of the proliferation-adjusted coupling**: the Basal column is broadly
  negative and concentrated on tumor-suppressor/cell-state programs (bracketed to the
  miR-17∼92/106b∼25 cluster), LumB carries a weaker version, and LumA/Her2 stay sparse
  — plus a compact Basal strip showing the signal holds under the E2F/MYC-independent
  proxy. The
  curation-bias result is carried in text (intro point 3), not the figure. Regenerate
  with `.venv/bin/python3 -m mirna_hallmark.make_grant_figure`.

---

## Open item for you to verify before submission

**The IRF1 immune route has been demoted by our own preliminary data, not just by
literature framing.** When we added a proliferation covariate to purity+HRD (under
all three proxies, including the E2F/MYC-independent one), the basal
IRF1↔miR-23a-3p/miR-106b-5p coupling **lost significance within Basal**
(p=0.43 / 0.62), and sequence-based TargetScan weighting nominates
**miR-130b-3p / miR-301** — not the proposed arms — as the top predicted IRF1
repressors. The route's edges remain validated in *other* tumors (3′UTR luciferase
in HCC/glioma; a published TGF-β→miR-23a-3p→IRF1→MHC-I axis), but in this cohort the
basal route is fragile. The draft therefore carries IRF1 as an **exploratory,
high-risk secondary** (a single surface-MHC-I flow readout in Aim 3), and leads on
the proliferation-/weighting-robust **PTEN/p21/TGFBR2/VIM** axis. Confirm you are
comfortable with this framing, and decide whether to keep IRF1 in the body
at all or move it entirely to "future directions."
