---
title: "Modeling miRNA Pressure on Cancer Hallmarks — NotebookLM Slides Companion"
purpose: "Concrete numbers + worked example + per-slide build spec to sharpen the NotebookLM slide deck generated from MODELING_FRAMEWORK_EXTERNAL.md"
audience: "lab work-in-progress talk, ~45 min"
---

# HOW TO USE THIS FILE

This file has two parts.

- **PART A — Facts & worked example.** Upload this (or the whole file) to NotebookLM
  **as an additional source**, alongside `MODELING_FRAMEWORK_EXTERNAL.md`. It carries the
  exact counts, the single worked edge, the protein survivors, the trajectory examples,
  and the 17 orphan edges — so NotebookLM has real numbers to ground on instead of
  paraphrasing them away.
- **PART B — Per-slide build spec.** Paste this into the NotebookLM **slides "Customize"
  box** (the instruction field). It is the slide-by-slide ordering with the four
  insertions, the slide-7 fix, the slide-12 split, and the slide-13 sharpening.

Everything below is grounded in the repo (TCGA-BRCA spine; Normal-like set aside →
**1,065 retained tumours**, realized tests on **~1,041**). Numbers re-derive each rebuild;
treat ρ/q to 2 sig-figs.

---
---

# PART A — FACTS & WORKED EXAMPLE  (upload as a source)

## A0. The universe in numbers (say these out loud)

| Quantity | Value |
|---|---|
| Cohort | TCGA-BRCA primary tumours; **1,065 retained** (PAM50 Normal-like **36 set aside**); realized coupling on **~1,041** with complete miRNA+mRNA+covariates |
| miRNA arms (expressible) | **721** |
| Curated spine edges | **5,108** (miRTarBase functional core) |
| Target genes under pressure | **1,410** (of **4,384** genes in the union of the 50 programs) |
| Programs | **50** MSigDB Hallmark gene sets |
| PAM50 strata (retained n) | Luminal A **507** · Luminal B **264** · HER2-enriched **93** · Basal-like **197** |
| Healthy reference | **GTEx breast 346** (+ DIANA-miTED atlas for coverage repair) |
| Normal-adjacent (NAT) | **104** (same platform; field-effect intermediate state) |
| Protein cohorts | **CPTAC** same-patient **n=105** + independent **n≈101** |
| Independent replication | **Buffa 2011**, **210** paired tumours, no TCGA overlap |
| Priors | **232** curated drivers (114 oncogene / 106 TSG / 12 dual) · DepMap **96** breast lines · **1,639** TFs |

**The five resolutions:** `Universe > Program > Target-set > Gene > Edge`.

## A1. The master formula (the spine of the whole deck)

$$c(m,g,s) = \mathrm{expr}(m,s)\cdot w(m,g) \qquad \mathrm{pressure}(g,s)=\sum_{m\in R(g)} c(m,g,s)$$

- **`c(m,g,s)`** = the **edge-resolution pressure atom** (one miRNA `m`, one gene `g`, one sample `s`).
- Two high-level components: **`w(m,g)` = interaction weight** (static, evidence-based) and
  **`expr(m,s)` = expression component** (dynamic, sample-varying).
- A gene's pressure is the **sum** of `c` over its regulators `R(g)`.

**Component previews (one line each, before the detail slides):**
- `w(m,g) = ℓ(ev_eff(m,g)) / D(m)` — evidence numerator ÷ promiscuity denominator, `ℓ(n)=ln(1+n)`.
- `expr(m,s) = sm · z · max(x,0)` (spine mode) — **share × dynamics × abundance level**.

## A2. THE WORKED EDGE — `hsa-miR-17-5p → PTEN`  (use this ONE edge for slides 7–9 and again at protein)

This single edge threads the whole deck: it is a curated spine edge, it is mRNA-coupled,
it is **protein-validated**, and its target **PTEN** is a 91-regulator hub for the
gene-aggregation slide. Use it instead of a generic "collagen axis."

**Identity**
- Edge: `hsa-miR-17-5p → PTEN` (PTEN = PI3K/AKT tumour-suppressor phosphatase).
- Seed family: `miR-17-5p / 20-5p / 93-5p / 106-5p / 519-3p` (the **miR-17~92 oncomiR cluster**).
- Hallmark sets the edge sits in: `PI3K_AKT_MTOR_SIGNALING`, `APICAL_JUNCTION`, `UV_RESPONSE_DN`.
- Biology: a proliferation-cluster oncomiR repressing a canonical tumour suppressor — the
  textbook **proliferation-cluster → tumour-suppressor** hub.

**(1) Interaction weight `w(m,g)` — static evidence**
- Curated assay counts for this edge: **reporter assays n=1, protein assays n=1**;
  passes the high-evidence gate (`high_evidence = True`).
- Evidence numerator (with `ℓ(n)=ln(1+n)`, `ℓ(1)=0.693`):
  $$\mathrm{ev}=3.0\,\ell(n_\text{rep})+2.5\,\ell(n_\text{prot})+1.5\,\ell(n_\text{rna})+1.5\,\ell(n_\text{pert})+0.5\,\ell(n_\text{bind})$$
  Here reporter=1, protein=1 ⇒ leading terms `= (3.0+2.5)·0.693 ≈ 3.81` (other classes 0 for this edge).
- `ev_eff = ev + 0.5·enc` (ENCORI crosslink/degradome depth boost where present).
- **`w(m,g) = ℓ(ev_eff) / max(D(m), ε)`** — `D(m)` is miR-17-5p's per-arm **promiscuity
  scalar** (total credible evidence mass over all its targets), so a hub miRNA can't
  dominate any one gene just by being promiscuous.
- (Pipeline integer `evidence_score` for this edge = **9**; the miR-17~92 sibling
  `miR-20a-5p → PTEN` = **15**, two studies.)

**(2) Expression component `expr(m,s)` — dynamic, three subcomponents (spine mode)**
$$\mathrm{expr}(m,s)=\underbrace{\mathrm{sm}_{m,s}}_{\text{share among PTEN's regulators}}\times\underbrace{z_{m,s}}_{\text{within-tumour dynamics}}\times\underbrace{\max(x_{m,s},0)}_{\text{abundance level}}$$
with `x = log2(RPM+1)`, `z = (x−μ)/σ` (standardized within scope), and
`sm` = softmax share of miR-17-5p's abundance among the regulators converging on PTEN.
(Attribution mode drops `z`; acquired mode replaces `z` by deviation-from-healthy `dev = x − h`.)

**(3) Contribution `c` = `expr · w`** — the pressure this edge places on PTEN in sample `s`.

**(4) Realized coupling — does the pressure show up in tumours?**
Partial Spearman across **~1,041** tumours, miR-17-5p abundance vs **PTEN mRNA**,
confounders (purity, HRD/instability, PTEN copy number, proliferation, composition, batch) removed:

| Edge | partial-ρ | significance | reading |
|---|---:|---|---|
| **miR-17-5p → PTEN** | **−0.29** | q_by ≈ **1e-19**, perm p≈5e-4 | strong realized repression |
| miR-20a-5p → PTEN (sibling) | −0.32 | FDR-sig | cluster coheres |
| miR-106b-5p → PTEN (sibling) | −0.34 | FDR-sig | cluster coheres |
| seed-family aggregate → PTEN | — | q_seedfamily ≈ **1.5e-26** | family-level call |
| **miR-21-5p → PTEN** (contrast) | **−0.03** | **n.s.** (q≈0.30) | the *famous* pair is pressure-present but realization-weak — the platform's point |

**(5) Protein layer (CPTAC, independent cohort)** — PTEN protein residualized on its own
mRNA, widest covariate stack: **resid-ρ = −0.34, q = 0.008** (85 candidate cluster arms;
strongest in Basal). The proliferation-cluster→TS hub is recovered **at the protein level**.

> Teaching contrast to make explicit: **miR-21→PTEN is the textbook pair but does NOT
> realize** (ρ=−0.03); the miR-17~92 cluster does, at mRNA and protein. Realized ≠ famous.

## A3. PTEN AS A GENE — aggregation in practice (for the new gene slide before slide 10)

PTEN carries **91 regulators**. Aggregation answers two *different* questions, via two modes:

| Mode | Formula piece | "Who lands hardest on PTEN?" | Leader |
|---|---|---|---|
| **Attributive** (no-z: `sm·level`) | steady-state share | which arm carries the most **standing** pressure | **miR-9-5p**, gene-side share **0.080** (#1; rose 0.048→0.080, grip tightening) |
| **Dynamic** (spine: `sm·z·level`) | share × dynamics | which arm's **variation** tracks PTEN down | **miR-17~92 cluster** (17-5p/20a/106b, the coupled arms above) |

- The attribution leader (**miR-9-5p**) and the coupling leaders (**miR-17~92**) are
  *different arms* — exactly why both modes exist.
- Arm-side vs gene-side: PTEN is miR-9-5p's **#1** gene (gene-side share 0.080) but PTEN is
  a **trivial 0.006** slice of miR-9-5p's own portfolio (arm-side `edge_share`) — opposite
  questions, both kept.
- PTEN is a **heavily-contested** gene: it always keeps *some* coupled regulator while
  individual arms decouple, so no single arm "controls" it. Cohort context: **282 / 1,424**
  genes are net-repressed in tumour (gene-level partial-ρ of own expression vs aggregate
  incoming pressure).

## A4. PROTEIN REALIZATION — CPTAC (slide 12a, split out)

Three tests on mass-spec proteomes (standardized mRNA & protein):

| Test | Quantity | Expected |
|---|---|---|
| direct protein | ρ(pressure, protein) | negative |
| mRNA→protein gap | ρ(pressure, mRNA − protein) | positive |
| **protein beyond mRNA** | ρ(pressure, protein \| mRNA) | negative — the key test |

**Named FDR-significant protein-residual survivors** (prospective CPTAC, resid on own
mRNA, widest covariates):

| Gene ← driver arm | resid-ρ (q) | note |
|---|---:|---|
| **KRT8** ← miR-186-5p | −0.44 (6e-4) | keratin 8 |
| **DDAH1** ← miR-210-3p | −0.42 (9e-4) | hypoxia oncomiR |
| **BMP1** ← miR-194-5p | −0.40 (1e-3) | also a miR-29 target |
| **PEBP1** ← miR-155-5p | −0.39 (1e-3) | RKIP, MAPK brake |
| **ARHGDIA** ← miR-151a-5p | −0.39 (1e-3) | Rho-GDI |
| **MINPP1** ← miR-30b-5p | −0.35 (6e-3) | inositol phosphatase |
| **PTEN** ← miR-17-5p | **−0.34 (8e-3)** | **the cluster→TS hub (our worked edge)** |
| **HMGN2** ← miR-23a-3p | −0.33 (1e-2) | chromatin |
| **ZYX** ← miR-16-5p | −0.31 (2e-2) | zyxin |

**Headlines:** dozens of genes FDR-significant; repression is **predominantly mRNA-mediated**
(residualizing on mRNA collapses most of it, leaving a small beyond-mRNA residual);
survivors are **known, subtype-structured** edges incl. the PTEN hub; a **decoy control**
confirms target-specificity (each arm hits its cognate protein, not the proteins it does
not target); survives batch adjustment.
**Caveat to state:** prognostically the EMT signal runs through **LOX / ECM, not ZEB1/2**
(per the outcome deep-dive) — do **not** call ZEB1/2 the protein driver.

## A5. HEALTHY → NAT → TUMOUR TRAJECTORY (slide 12b, split out)

Three states: **GTEx healthy (346) → NAT (104) → tumour**. Two acquired axes (union call):
- **rank delta** `d_HT = r_tumour − r_healthy` (saturating, cross-platform-robust),
- **magnitude** `log2fc = x̃_tumour − b_healthy` (recovers the surge).
- BC oncomiR surge example: **miR-183 / 182 / 375 / 93 / 141, +9 to +13 log2**.

**Released brakes** — repression present in healthy, gone in tumour (`R00` pattern):
**189 edges / 109 arms.**

| Edge | ρ healthy → tumour | biology |
|---|---|---|
| **miR-193a-5p → ERBB2** | −0.28 → +0.07 | HER2 brake released |
| miR-224-5p → CXCR4 | −0.25 → +0.02 | metastatic-homing brake |
| miR-22-3p → MMP14 · miR-708-5p → MMP2 | −0.30 → +0.19 · −0.31 → +0.03 | invasion (MMP de-repression) |
| miR-326 / miR-23b-3p → NOTCH2 | −0.33 / −0.30 → ~0 | Notch de-repression |
| miR-363-3p → BCL2L11(BIM) | −0.28 → ~0 | apoptosis brake |

*Theme:* healthy breast actively holds oncogenic/invasion targets (ERBB2, CXCR4, MMP2/14,
NOTCH2) under miRNA repression; the tumour **releases the brake**.

**Acquired-and-realized** — pressure rises *and* coupling turns repressive: **640 edges / 129 arms.**

| Edge | ρ healthy → tumour | biology |
|---|---|---|
| **miR-744-3p → HLA-G** + **→ TGFB1** | +0.01 → −0.10 · +0.01 → −0.24 | **novel immune axis** (MHC tolerance ligand + TGFβ) |
| miR-590-5p → RBPJ | +0.02 → −0.11 | acquired Notch-effector repression |
| miR-301a-3p → BTG1 | −0.11 → −0.20 | strengthening anti-proliferative-TS repression |

(Field-established `0RR`: miR-21-5p→TGFBR3, miR-375→AHR, miR-25→PTEN.)

## A6. ORPHAN DISCOVERY SURFACE — and THE 17 (slide 13, sharpened)

**Orphan edge** = TargetScan/ENCORI-predicted, **absent from miRTarBase functional
curation** — deliberately kept out of the spine so it is an unbiased screen.

**The funnel:**
```
20,663 orphan edges screened (TargetScan ∪ ENCORI, no functional MTI)
   → 15,512 testable vs CPTAC protein
      → 539 protein anti-correlation candidates
         → 182 translational (beyond-mRNA residual)
            → 93 translational AND TCGA-replicated   (76/93 fully uncurated)
```
- **Triple-cohort (TCGA mRNA · CPTAC protein · Buffa mRNA):** of **492 genuine orphans**,
  **127** were CPTAC-sig + TCGA-replicated + Buffa-testable; **108 (85%)** keep a negative
  sign in Buffa (47 reach FDR) — far above the 50% chance rate; led by **miR-29→ECM (34/108)**.
- **De-novo anchor:** the **miR-29 → collagen/ECM** axis re-emerges *at uncurated targets*
  (miR-29a→COL11A1 ts_weight 0.80, 15 CLIP experiments; →COL6A2, →BMP1) — credibility check.

**The 17 composition-robust core (MH-39).** Re-test the prediction-only tier (TargetScan
sequence-only — **no CLIP, no curation**) after adding PAM50 + ESTIMATE immune/stromal +
endothelial + target-copy-number covariates: **17 of 24 survivors hold** (negative
beyond-mRNA `protein_resid` at q<0.1 under the widest covariate set; **11/17 at q<0.05**).
The 7 drops are exactly the composition artifacts (VCL×2, TPM4, and the hematopoietic
miR-142 cluster). **13 target genes, 17 edges:**

| Target (function) | TargetScan-only arm(s) | protein ρ | resid-ρ (q) |
|---|---|---|---|
| **PMM1** — mannose metabolism | miR-195-5p · miR-16-5p · miR-15a-5p | −0.40 · −0.39 · −0.39 | −0.40 (9e-4) · −0.31 (.012) · −0.30 (.012) |
| **PPM1A** — PP2C phosphatase (TGF-β/AKT brake) | miR-19b-3p · miR-19a-3p | −0.30 · −0.32 | −0.28 (.016) · −0.23 (.042) |
| **PURA** — transcription / replication | miR-532-5p · miR-452-5p | −0.45 · −0.33 | −0.19 (.093) · −0.18 (.093) |
| **KIF1B** — kinesin, pro-apoptotic | miR-15a-5p | −0.37 | −0.30 (.012) |
| **NCK2** — signalling adaptor | miR-142-3p | −0.07 | −0.31 (.012) |
| **VCAN** — versican (ECM) | miR-107 | −0.52 | −0.27 (.018) |
| **MAPRE1** — EB1, microtubule +end | miR-429 | −0.25 | −0.27 (.018) |
| **NUDT4** — Nudix hydrolase | miR-183-5p | −0.16 | −0.26 (.023) |
| **RAB1A** — ER→Golgi traffic | miR-101-3p | −0.33 | −0.24 (.041) |
| **CERCAM** — endothelial adhesion | let-7g-5p | −0.35 | −0.22 (.057) |
| **SKP1** — SCF ubiquitin-ligase core | miR-452-5p | −0.38 | −0.21 (.066) |
| **DR1** — NC2 transcriptional repressor | miR-182-5p | −0.34 | −0.19 (.093) |
| **GLUD1** — glutamate dehydrogenase | miR-30e-5p | −0.28 | −0.18 (.093) |

**What to emphasize about the 17:**
- **Counts:** 17 edges, 13 genes, **11/17 at q<0.05**; all from sequence-only predictions
  that survived a full composition + copy-number confounder stack.
- **Two seed-family convergences lead and are CAF-neutral / method-independent:** the
  triple-hit **miR-15/16/195 → PMM1** (a **metabolic convergence**, mRNA-flat but
  protein-strong = translational) and **miR-19a/b → PPM1A**.
- **NCK2 is the instructive single case:** direct protein ρ ≈ **0**, yet the **beyond-mRNA
  residual is strong (−0.31, q=.012)** — a translational-layer signal bulk mRNA cannot see.
- **General character of the genes:** mostly **metabolism** (PMM1, GLUD1, NUDT4),
  **signalling brakes** (PPM1A, NCK2, SKP1, RAB1A), **transcription/replication**
  (PURA, DR1, MAPRE1), plus single ECM (VCAN), apoptosis (KIF1B), endothelial (CERCAM) —
  i.e. NOT the curated collagen crowd; these are *new* metabolic/signalling nominations.
- **Honest status:** these are **wet-lab nominations** (a prioritized CLIP/luciferase
  queue), not validated edges — protein anti-correlation ≠ binding, co-expression confound
  persists, and seed-family arms (e.g. miR-15/16) cannot be separated.

## A7. VALIDATION & RIGOR — "why believe it"  (feeds the new slides 17–19, *before* the original ladder/limits)

### A7.1 Positive controls — known biology recovered *unprompted* (non-circular)

The priors — interaction weights, the high-evidence gate, the program-polarity prior,
gene-role annotations — **never see TCGA patient expression**. So when textbook biology
falls out of the patient data, it is genuine **internal validation**, not a fitted echo.

- **Tumour-suppressor-miRNA roster** emerges as the architecture **anti-tumour leaders**:
  **miR-30a, miR-24, miR-29b/c, let-7a**.
- **OncomiR roster** emerges from the **gene-role** malignancy score:
  **miR-21, miR-182, miR-10b, miR-93, miR-106b**.
- **Three canonical axes re-emerge de novo / as unsupervised sub-modules**, never hand-fed:
  **miR-29 → collagen/ECM**, **miR-17~92 → PTEN / replication** (our worked edge), and
  **let-7 / miR-30 → mitotic**. Their *absence* would have made the decomposition suspect.
- **High-evidence edges couple more strongly than sequence-only predictions** in **two**
  independent cohorts (TCGA + Buffa) — the curation gate is doing real work.

### A7.2 Negative controls — the signal is target-specific, not a backdrop

- **Decoy / target-specificity control** (`cptac_target_specificity`): compare each driver
  arm's coupling to its **cognate** protein vs **≥50 proteins it does NOT target** (decoys).
  Cognate couplings (median **ρ = −0.33**) are **~16× more negative** than the decoy null
  (**≈ −0.02**, centred at zero); **9/13 FDR drivers are target-specific.** The honest
  exceptions are **broad proliferation oncomiRs (PTEN←miR-17-5p)** whose decoy nulls are
  *themselves* somewhat negative — because they genuinely target many proliferation genes.
- **Empirical permutation null** (Freedman–Lane rank-residual; one shared sample
  permutation per draw, so seed-family dependence is preserved in the null) re-tests every
  coupling at **every resolution** — significance never rests on the asymptotic p alone.
- **Composition drop-out:** under the full immune/stromal/endothelial/CN stack, the **7
  composition artifacts drop out** (VCL ×2, TPM4, the hematopoietic miR-142 cluster) while
  the **17 hold** — confounded edges behave exactly as a negative control should.
- **Variability gate:** flat arms (low SD *and* low IQR) are flagged so their null coupling
  is not over-read.

### A7.3 Inference rigor

- **Dependence-aware FDR:** Benjamini–Yekutieli (valid under arbitrary dependence) **plus**
  a seed-family Simes collapse — co-varying paralogues spend **one** slice of the error
  budget, not several; discovery counts reported in **families**, not edges.
- **Out-of-sample, not overfit:** the pressure construction's candidate ranking reproduces
  on held-out patients — **5-fold CV, train→test rank-ρ = 1.0**.
- **Constants are a prior, not a fit:** ±25% random and ±50% systematic jitter of the
  evidence weights leaves the headline essentially unchanged (**CV ≈ 0.06**).

### A7.4 Independent-cohort replication (Buffa 2011, 210 tumours, no TCGA overlap)

- **Spine:** partial-ρ concordance **+0.32**; **67%** of TCGA FDR-negative edges keep a
  negative sign in Buffa — an independent-patient leg the same-patient proteome cannot give.
- **Orphans:** **108** triple-cohort-validated (TCGA mRNA · CPTAC protein · Buffa mRNA),
  **85%** sign-preserved — far above the 50% chance rate.

### A7.5 The validation ladder

| Rung | Evidence | Strength |
|---|---|---|
| Reproduces known biology | textbook rosters + axes recovered unprompted | established |
| Robustness & inference | permutation null + BY/seed-family FDR; 5-fold CV rank-ρ 1.0; constants CV 0.06 | statistical |
| mRNA coupling | FDR-negative, confounder-adjusted partial-ρ | statistical |
| State-resolved trajectory | brake-release / acquired-realized classes | descriptive / hypothesis |
| Protein layer | independent CPTAC proteome residual; decoy-specific | statistical |
| Independent cohort | Buffa +0.32, 67% sign-preserved | statistical / established |
| Wet-lab queue | 17 composition-robust orphans + 76 uncurated nominations | hypothesis |

---
---

# PART B — PER-SLIDE BUILD SPEC  (paste into the slides "Customize" box)

> Build a ~60-minute, ~21-slide work-in-progress lab talk from the two sources. Be
> quantitative: when a number, formula, or gene appears in the sources, **put it on the
> slide** — do not round to "several" or paraphrase it away. Use the companion file's
> worked edge (**miR-17-5p → PTEN**) as a through-line. Produce exactly these slides in
> this order:

1. **Title** — "Modeling miRNA Repressive Pressure on Cancer Hallmark Programs (breast cancer)."
2. **The three rungs** — pressure (potential) vs coupling (realized) vs causation; why they are kept distinct.
3. **Data & cohort, in numbers** — 1,065 retained tumours (36 Normal-like set aside), ~1,041 for coupling; 721 arms, 5,108 edges, 1,410 genes, 50 programs; PAM50 507/264/93/197; GTEx 346, NAT 104, CPTAC 105 + 101, Buffa 210.
4. **Five resolutions + the interaction universe** — Universe > Program > Target-set > Gene > Edge; spine = miRTarBase functional core; orphans held out.
5. **★ NEW — The master formula at a glance.** Show `c(m,g,s)=expr(m,s)·w(m,g)` and `pressure(g,s)=Σ c`. Name the **two high-level components**: `w(m,g)` (static evidence weight) and `expr(m,s)` (dynamic expression). Give the one-line previews `w = ℓ(ev_eff)/D(m)` and `expr = sm·z·level`. Say "the next slides open each component." (Insert this BETWEEN current slide 4 and current slide 5.)
6. **The evidence gate & weighting rationale** — high-evidence gate; evidence numerator (assay-directness weights 3.0/2.5/1.5/1.5/0.5, log-damping, weak-evidence discount) ÷ promiscuity denominator `D(m)`; weights are a transparent prior, not a fit.
7. **★ FIX — Worked edge weight: `miR-17-5p → PTEN` (ONE specific edge, not "the collagen axis").** An edge is one miRNA → one gene, so pick this one. Show its real assay counts (reporter=1, protein=1, high-evidence), the populated `ev` numerator (≈3.81 from the two direct-assay terms), and `w = ℓ(ev_eff)/D(miR-17-5p)`. Name the seed family (miR-17~92) and the biology (oncomiR cluster → tumour-suppressor PTEN).
8. **★ FIX — The expression component, whole then parts.** First show the **entire** spine expression formula `expr(m,s)=sm·z·max(x,0)`. THEN break it into its three subcomponents: `z` (within-tumour dynamics), `sm` (softmax share among a gene's regulators), `max(x,0)` (abundance level). Note the two sibling modes (no-z attribution; deviation-from-healthy acquired).
9. **★ NEW — Worked edge, expression → contribution.** Same `miR-17-5p → PTEN`: walk its three expression subcomponents, then combine with the weight from slide 7 to produce the **contribution `c = expr·w`** — the pressure atom this edge places on PTEN. (Insert AFTER slide 8.)
10. **★ NEW — Worked edge, realized coupling.** Same edge: partial Spearman across ~1,041 tumours of miR-17-5p abundance vs PTEN mRNA, confounders removed → **ρ=−0.29 (q_by≈1e-19)**; siblings miR-20a/106b also negative; seed-family q≈1.5e-26. **Contrast: miR-21-5p→PTEN ρ=−0.03, n.s.** — the famous pair does not realize. (Insert AFTER the slide above, BEFORE the gene-aggregation slide.)
11. **★ NEW — PTEN as a gene: aggregation in practice.** PTEN = **91 regulators**. Show pressure aggregated two ways: **attributive (no-z)** leader **miR-9-5p** (gene-side share 0.080, #1) vs **dynamic (spine-z)** coupling leaders **miR-17~92**; note arm-side vs gene-side (PTEN is miR-9-5p's #1 gene but a 0.006 slice of its portfolio). (Insert BEFORE the original gene-aggregation formula slide.)
12. **Gene aggregation formula** (original slide 10) — `pressure(g)=Σ_m c`; share & specificity; arm-side vs gene-side shares; sum (not mean) is the default roll-up.
13. **Architecture** — topology influence (OmniPath signed propagation) + program-polarity prior + gene-role malignancy; the result that miRNAs turn malignant mainly by **brake-release on suppressive programs**.
14. **★ SPLIT 12a — Protein realization (CPTAC).** Three tests (direct / mRNA-protein gap / **protein beyond mRNA**). Named FDR-significant survivors: **KRT8, DDAH1, BMP1, PEBP1, ARHGDIA, MINPP1, PTEN (the hub, −0.34 q=.008), HMGN2, ZYX**. Headlines: dozens FDR-sig, **predominantly mRNA-mediated**, target-specific (decoy control), batch-robust. **Do not call ZEB1/2 the driver — the EMT prognostic signal is LOX/ECM.**
15. **★ SPLIT 12b — Healthy→NAT→tumour trajectory.** Three states (GTEx 346 → NAT 104 → tumour); two acquired axes (rank delta + log2fc); oncomiR surge miR-183/182/375/93/141 (+9..+13 log2). **Released brakes (189 edges/109 arms): miR-193a→ERBB2, miR-224→CXCR4, miR-22/708→MMP14/MMP2, miR-326/23b→NOTCH2.** Acquired-realized (640/129): **miR-744-3p→HLA-G + TGFB1** (novel immune axis).
16. **★ SHARPEN 13 — Orphan discovery surface + the 17.** Show the funnel (20,663 → 15,512 → 539 → 182 → 93; triple-cohort 108/127). Then show **the 17 composition-robust edges explicitly** (the companion's 13-gene table) with anti-correlation values and gene annotations. Call out: **11/17 at q<0.05**; the metabolic lead **miR-15/16/195→PMM1**; **NCK2** as translational-only (direct ρ≈0, residual −0.31); general character = metabolism + signalling brakes, not the curated ECM crowd; miR-29→collagen as the de-novo anchor. Label all 17 **wet-lab nominations**.
17. **★ NEW — Positive controls: known biology recovered *unprompted*.** Priors never saw TCGA expression, so textbook recovery = internal validation. TS-miRNA roster (architecture leaders **miR-30a, miR-24, miR-29b/c, let-7a**); oncomiR roster (gene-role: **miR-21, miR-182, miR-10b, miR-93, miR-106b**); three axes re-emerge de novo (**miR-29→ECM, miR-17~92→PTEN, let-7/miR-30→mitotic**); high-evidence edges couple more than sequence-only in two cohorts.
18. **★ NEW — Negative controls + inference rigor.** Decoy/target-specificity: cognate coupling median **ρ=−0.33** vs decoy null **≈−0.02** (**~16×**), **9/13 drivers target-specific** (broad proliferation oncomiRs the honest exception). **Freedman–Lane permutation null at every rung**; **Benjamini–Yekutieli + seed-family** FDR; **5-fold CV rank-ρ=1.0** (not overfit); constants **±25%/±50% jitter, CV≈0.06**; the 7 composition artifacts drop out while the 17 hold.
19. **★ NEW — Independent replication (Buffa 2011, 210 tumours, no TCGA overlap).** Spine partial-ρ concordance **+0.32**; **67%** of TCGA neg-sig edges keep a negative sign; **108** triple-cohort-validated orphans (85% sign-preserved) — the independent-patient leg the same-patient proteome can't provide.
20. **Validation ladder & honest limits** *(original — KEEP)* — recovered (textbook edges unprompted) → robustness/inference → mRNA coupling → state-resolved trajectory → protein layer → independent cohort → wet-lab queue; limits: bulk-tissue ceiling, metabolite-mediated regulation invisible to topology, seed-family non-identifiability, pressure construction non-unique.
21. **Where the lab comes in** *(original — KEEP)* — 2–3 specific open questions (e.g. which of the 17 to CLIP/luciferase first; a matched-NAT-miRNA + proteome cohort to power acquired→protein; resolving seed-family attribution).

> Style: define every symbol in plain words the first time. Flag preliminary results as
> preliminary. End slide 21 by asking the lab the open questions.
