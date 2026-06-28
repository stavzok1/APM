---
title: "Porting the miRNA-Pressure Machine to Colorectal Cancer"
subtitle: "Literature landscape and where the machine wins"
date: "2026-06-28"
geometry: margin=1in
colorlinks: true
linkcolor: RoyalBlue
urlcolor: RoyalBlue
fontsize: 11pt
mainfont: "DejaVu Serif"
sansfont: "DejaVu Sans"
monofont: "DejaVu Sans Mono"
---

## 0. What the machine is (the lens that makes the literature gap legible)

The `mirna_hallmark` engine is a **potential × realized** framework, not a
biomarker-panel finder. It separates **pressure** (evidence-weighted,
expression-scaled *potential* repression over miRNA→target edges, modeled at
five nested resolutions and rolled up over the 50 MSigDB Hallmark programs) from
**coupling** (confounder-adjusted partial-Spearman *realization*), tests
realization at the **protein layer**, classifies movement along a
**healthy→NAT→tumour trajectory** into brake-release / acquired-realized classes,
and has a **compartment/EV extension** (the DCIS arc, MH-48..56). That arc's own
closing verdict was that it was blocked by exactly two things: *no within-patient
stage pairs* and *no compartment-resolved small-RNA*. Hold that thought — it is
the whole thesis.

The entire CRC literature below is, by contrast, **descriptive differential
expression plus correlational biomarker panels.** That mismatch is the
opportunity.

---

## 1. The CRC literature landscape, by axis

### 1.1 The premalignant ladder is real, dual, and exhaustively *described* — but never modeled as realization

CRC is the textbook multistep cancer, and unlike breast it has a genuinely
**sampled premalignant ladder** with two molecularly distinct routes:

- **Conventional adenoma–carcinoma sequence** (APC/WNT → KRAS → p53), and
- **Serrated pathway** (BRAF-mutant, CIMP-high, frequently MSI, right-sided)
  ([serrated pathway hallmarks](https://pmc.ncbi.nlm.nih.gov/articles/PMC6678087/);
  [proximal/BRAF/MSI link](https://www.nature.com/articles/modpathol201298)).

The miRNA dynamics across these steps are profiled in depth and, importantly,
**diverge by pathway**: in the adenoma lineage miR-20a/miR-21/miR-181b rise
progressively and miR-93 increases stepwise, whereas the serrated lineage shows
selective downregulation/plateauing
([Murakami 2025, serrated vs AD miRNA dynamics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12279426/)).
The "multi-hit" framing — different miRNAs switching at normal→adenoma vs
adenoma→carcinoma — is established
([Oberg, PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020465);
[normal-adenoma-carcinoma transition](https://pmc.ncbi.nlm.nih.gov/articles/PMC6607403/);
[polyp landscape](https://pubmed.ncbi.nlm.nih.gov/27925331/)). Specific drivers
are mechanistically nailed: the **miR-17-92 cluster** at 13q is MYC-driven, rises
across adenoma→carcinoma with 13q gain, and forces a carcinoma-like signature
when overexpressed in adenoma organoids, with **miR-92a** the dominant component
([13q/MYC](https://www.nature.com/articles/6605037);
[organoid](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9307940/)); **miR-31**
tracks the BRAF/serrated/CIMP axis
([EZH2/miR-31 serrated](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914316/));
**miR-34a** is a p53 target and is methylation-silenced, substituting for p53 loss
([miR-34a methylation/metastasis](https://aacrjournals.org/clincancerres/article/19/3/710/208961/));
**miR-143/145** is repressed by oncogenic KRAS in a feed-forward loop
([Ras feed-forward](https://genesdev.cshlp.org/content/24/24/2754.full)).

**The gap:** every one of these is a fold-change or a single-edge knockdown.
Nobody has asked, across the ladder, *which brakes are engaged in normal mucosa
and released at which step, after confounders, and confirmed at the protein
layer.* That is precisely the §6.4 join + §8 trajectory machinery — and CRC
offers a **paired, multi-state, real-lesion ladder** that breast (where the DCIS
arc had no within-patient stage pairs) structurally could not.

### 1.2 The single-cell atlases give the compartment map — but they are miRNA-*blind*

Two landmark atlases now chart polyp→CRC at cell-state resolution: **Chen et al.,
Cell 2021** (62 tumours / 128 datasets; conventional adenomas = WNT-driven
stem-cell expansion, serrated polyps = gastric metaplasia)
([Cell 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)01381-7)) and
**Becker, Curtis et al., Nat Genet 2022** (scATAC + scRNA over 48 polyps / 27
normal / 6 CRC; a stem-like continuum, FOX-regulated *pre-CAFs* maturing into
RUNX1-CAFs, expanding Tregs, and T-cell exhaustion at the cancer step)
([Nat Genet 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9279149/)).

These are **scRNA / scATAC / protein — there is no miRNA layer.** This is the
single most important structural fact: the field has built the
compartment/cell-state scaffold of CRC progression *with a miRNA-shaped hole in
it.* The machine lives in exactly that hole, and the atlases supply the
cell-state priors the composition covariates need (the pre-CAF→CAF and
myeloid/Treg programs are the CRC analog of the CAF stroma-reframe the DCIS arc
hit at MH-52/54).

### 1.3 Field cancerization: NAT is not neutral — and already carries prognostic miRNA signal

The "normal adjacent" tissue in CRC is a **field-cancerized** state, and miRNA
dysregulation extends into it
([field cancerization review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4385523/)).
Concretely, a **four-miRNA coordinate deregulation in tumour-adjacent mucosa
(miR-18a/miR-21/miR-182/miR-183) predicts relapse after curative resection**
([NAT relapse predictor](https://pmc.ncbi.nlm.nih.gov/articles/PMC5791208/)).
This validates the decision to model NAT as an *intermediate trajectory state*
rather than a clean baseline — and it means the acquired/field-effect classes
have a ready-made outcome anchor in CRC.

### 1.4 Tissue miRNA–mRNA networks & TCGA: anticorrelation done naively

TCGA-COAD/READ has been mined repeatedly for miRNA–mRNA anticorrelation networks
(hundreds of DE miRNAs, TargetScan/miRanda-filtered Pearson anticorrelation)
([integrative knowledge base](https://www.nature.com/articles/s41598-019-54358-w);
[regulatory network](https://pubmed.ncbi.nlm.nih.gov/33172374/)). These are
**unweighted, unconfounded, single-edge** correlations — no evidence weighting,
no promiscuity control, no purity/proliferation/CN adjustment, no permutation
null. This is the exact methodology the framework was built to supersede; running
the spine here is a direct, defensible upgrade with a known baseline to beat.

### 1.5 The protein layer: CPTAC-CRC exists, with the right shape, and is unused for pressure

[**Vasaikar et al., Cell 2019**](https://pmc.ncbi.nlm.nih.gov/articles/PMC6768830/)
is a **prospective 110-patient colon cohort with tumour, matched NAT, and
blood**, profiled by WXS, CNV array, RNA-seq, **miRNA-seq, and label-free
proteome (+phosphoproteome)** — the same data shape the §7 protein-realization
layer already consumes for breast, with the bonus of **matched blood**.
Proteogenomic discordance / post-transcriptional buffering in CRC is itself a
documented phenomenon (mRNA–protein decoupling, 3′-UTR/APA effects)
([post-transcriptional CRC review](https://pmc.ncbi.nlm.nih.gov/articles/PMC6693420/);
[APA dysregulation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7047281/)),
which is exactly the "beyond-measured-mRNA residual" the §7 layer isolates.
**Nobody has run evidence-weighted miRNA pressure → protein realization on
CPTAC-CRC.** This is the cleanest *runnable-today* paper section and the
non-circular positive control.

### 1.6 Circulating / EV biology: mechanistically rich, clinically incoherent

This is the loudest CRC literature and the messiest:

- **Selective EV sorting is a real, mechanistic phenomenon**, not random leakage:
  sumoylated **hnRNPA2B1** sorts specific miRNAs (GGAG/EXO-motif) into exosomes,
  **miR-1246** via APE1/hnRNPA2B1, and there is a **KRAS-dependent sorting
  program**
  ([hnRNPA2B1 sorting, Nat Commun](https://www.nature.com/articles/ncomms3980);
  [KRAS-dependent secreted miRNA](https://www.oncotarget.com/article/20117/text/)).
  So "which arm is exported" is *governed*, not noise — and the DCIS arc already
  built the transcriptional-loss vs vesicular-export distinction that this
  mechanism underwrites.
- **CAF/stromal exosomes** carry a defined CRC signature (miR-21 / miR-329 /
  miR-181a / miR-199b / miR-382 / miR-215; miR-92a-3p driving EMT/chemoresistance;
  miR-93-5p radioresistance) and concentrate in **CMS4**
  ([CAF exosomal miRNA](https://pubmed.ncbi.nlm.nih.gov/29283887/);
  [miR-93-5p CAF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7158087/)).
  **Tumour exosomes build the liver pre-metastatic niche** (miR-25-3p vascular
  permeability, miR-934 M2 polarization, miR-188-3p hepatic stellate activation)
  ([niche review](https://pmc.ncbi.nlm.nih.gov/articles/PMC8075282/)).
- **Diagnostic panels are abundant but mutually contradictory** — the field's own
  reviews attribute the non-replication to **sample-type dependence** (free vs EV
  vs EPCAM⁺)
  ([exosomal CRC diagnosis](https://pmc.ncbi.nlm.nih.gov/articles/PMC7795745/);
  [combinatorial profiles](https://link.springer.com/article/10.1007/s12672-024-01481-4)),
  and the **preanalytical literature explains why translation has failed**:
  hemolysis distorts key arms (miR-451a), there is no consensus endogenous
  normalizer, and spike-ins don't rescue it
  ([preanalytical challenges](https://academic.oup.com/clinchem/article/57/6/833/5621214);
  [hemolysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8584813/)). Fecal
  miR-92a/miR-223 complement FIT for screening
  ([fecal complement](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4891149/)).

**The gap is glaring:** the circulating field has *mechanism* (selective sorting)
and *clinical noise* (discordant panels) but **no model that links a circulating
arm back to its tissue compartment's pressure/coupling state** and adjudicates
*why* it is in blood (transcriptionally elevated in tumour vs selectively
exported vs stromal-CAF-derived). That linkage is the thing the machine uniquely
produces.

### 1.7 Subtype, immune, and outcome structure (the CRC-specific strata)

CRC's taxonomy is **CMS1-4** (MSI-immune / canonical-WNT-MYC / metabolic /
mesenchymal-TGFβ), replacing PAM50
([CMS, Nat Med](https://www.nature.com/articles/nm.3967)), with mature classifiers
including a **miRNA-based one (CMS-miRaCl)**
([CMS-miRaCl](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297751/);
[CMScaller](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5709354/)). miRNAs gate
immune evasion (miR-148a-3p and miR-15b-5p → PD-L1, MSI-linked)
([miRNA/PD-L1/MSI](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8551827/)). On
outcome, **tissue miR-21 is a robust prognostic marker** (OS HR ≈ 1.6–1.8 in
meta-analysis)
([miR-21 meta-analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC3827229/)), while
for **MRD/recurrence ctDNA dominates** (CIRCULATE-Japan / GALAXY) and circulating
miRNA remains a promising-but-underpowered, mechanism-free add-on
([ctDNA MRD](https://aacrjournals.org/clincancerres/article/31/2/328/751096/);
[miRNA clearance panel](https://onlinelibrary.wiley.com/doi/10.1002/jso.70166)).

### 1.8 One CRC-specific twist the gate should exploit: AGO2 is *lost*

Unlike the generic up-regulation of DROSHA/DICER in CRC, **AGO2 protein is reduced
in CRC and predicts worse DFS/OS, driving EMT and metastasis**
([AGO2/miR-185-3p/NRP1](https://www.nature.com/articles/s41419-021-03672-1);
[Drosha/Dicer up](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605964/)). The §4.2
**availability gate** assumes RISC capacity (AGO/TNRC6) bounds realized repression
— in CRC that assumption has a real, directional, prognostic substrate, so the
gate is *more* mechanistically motivated here than in breast, and AGO2 loss
becomes a finding to fold in rather than just a sensitivity layer.

---

## 2. The white space, stated precisely

Putting the landscape together, every CRC sub-literature is missing the same
thing, and it is the thing the machine *is*:

| CRC literature has… | …but lacks (what the machine adds) |
|---|---|
| Exhaustive DE across the adenoma/serrated ladder | Potential×realized separation; confounder-adjusted **brake-release classification** along the ladder |
| Single-cell compartment/cell-state atlases (Chen; Becker-Curtis) | **Any miRNA layer** at compartment resolution — the atlases are scRNA/scATAC/protein |
| Naive TCGA miRNA–mRNA anticorrelation networks | Evidence weighting, promiscuity control, purity/proliferation/CN adjustment, permutation null, **protein realization** |
| CPTAC-CRC multi-omics incl. proteome | Evidence-weighted **pressure→protein** realization; the discordance is documented but unmodeled |
| Mechanistic EV sorting (hnRNPA2B1, KRAS) + discordant clinical panels | A model **linking a circulating arm to its tissue pressure state** and attributing it (transcriptional vs exported vs CAF) |
| ctDNA-led MRD; tissue miR-21 prognosis | A **mechanism-anchored** circulating signal tied to the tumour's realized repression field |

And critically: the **two limits that closed the DCIS/EV arc** — no within-patient
stage pairs, no compartment-resolved small-RNA — are **removed by the cohort
described** (within-patient blood→polyp→NAT→tumour, with followup).

---

## 3. How the machine maps onto CRC (the substrate port)

The framework is substrate-agnostic by design; only the cancer-flavoured priors
change:

- **Strata:** CMS1-4 + MSI/MSS + sidedness + BRAF/KRAS, replacing PAM50.
  Classifiers (CMScaller, ColoType, CMS-miRaCl) are off-the-shelf.
- **Healthy anchor (§8.1):** GTEx **colon (sigmoid + transverse)** replaces GTEx
  breast.
- **Edge universe (§3):** miRTarBase / TargetScan / ENCORI are tissue-agnostic and
  transfer unchanged; CRC supplies *strong positive controls* (miR-21→PDCD4/PTEN,
  miR-143/145→KRAS, miR-17-92→targets, miR-34a→targets) to validate the spine
  recovers known biology unprompted, exactly as the breast oncomiR roster did.
- **Confounders (§6.2):** same stack; the **CMS4/desmoplasia stromal axis is the
  same stroma-reframe** the DCIS arc resolved (CAF activation), and **MSI/immune
  composition** is the CRC analog of the immune covariate. The
  Becker-Curtis/Chen cell-state programs feed composition deconvolution.
- **AGO gate (§4.2):** keep it, and additionally treat **AGO2 loss** as a
  directional CRC finding.
- **Program-polarity prior + gene-role layer (§5):** pan-cancer, transfer directly.

---

## 4. Leverage map (what to run, in order)

| Dataset | Role in CRC port | Status |
|---|---|---|
| **TCGA-COAD/READ** | Primary spine — miRNA+mRNA+CNV+methylation; CMS/MSI strata; beats the naive anticorrelation networks head-to-head | Runnable now |
| **GTEx colon (sigmoid+transverse)** | True-healthy anchor for the acquired axis | Runnable now |
| **CPTAC-CRC (Vasaikar, 110; tumour+NAT+blood)** | Protein-realization (§7) + NAT field state; the rigor anchor; matched blood is a bonus | Runnable now |
| **Public EV/circulating + fecal GEO CRC** | Replicate export-selectivity & the transcriptional-loss vs vesicular-export split before the new cohort | Runnable now |
| **Chen 2021 / Becker-Curtis 2022 atlases** | Cell-state composition priors + the miRNA-blind scaffold this layer complements | Reference now |
| **New blood/polyp/NAT/tumour + followup cohort** | The firsts in §5 — within-patient ladder, direct export mechanism, outcome-anchored pressure | When it lands |

A full **pre-arrival validation** of the entire logic is runnable *today* on the
first four rows — which de-risks the new cohort before a single sample exists.

---

## 5. The headline the new cohort uniquely unlocks

> **A within-patient, compartment-resolved potential×realized "brake-release
> ledger" across the blood→polyp→NAT→tumour ladder, outcome-anchored by followup,
> in which every circulating arm is mechanistically attributed (transcriptional
> elevation vs selective EV export vs CAF-stromal origin) back to its tissue
> pressure/coupling state.**

This single design sits at the convergence of three gaps no existing CRC study
closes simultaneously:

1. the premalignant-ladder DE literature has no **realization** model (the §6.4
   brake-release classes on a *real, paired, dual-pathway* ladder);
2. the EV/circulating literature has no **tissue-of-origin mechanism** (the
   transcriptional-vs-export split, now testable *directly within patient* rather
   than triangulated across isogenic lines as in DCIS); and
3. the framework's own binding limits (no within-patient stage pairs, no
   compartment-resolved small-RNA) are **dissolved by the cohort design**.

It is also serrated-vs-conventional-pathway aware (a structure breast lacks
entirely), giving a built-in second axis of biology the machine can resolve that
no prior pressure analysis could.

---

## 6. Four things to confirm before committing a roadmap

These change whether §5's components are "firsts" or "stronger replications":

1. **Blood modality** — EV / total small-RNA-seq vs a qPCR panel? (Determines
   whether the direct within-patient export-attribution test is possible at all.)
2. **Polyp↔tumour pairing** — same patients, or separate arms? (Determines
   within-patient ladder vs cross-sectional.)
3. **Pathway coverage** — are *both* conventional adenomas and serrated/SSL polyps
   sampled? (Unlocks the dual-pathway axis.)
4. **Followup endpoint** — recurrence/DFS vs OS, and is MRD/ctDNA co-collected?
   (Sets the outcome anchor and the ctDNA-complementarity framing.)
