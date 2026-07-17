# Discussion prompt — building a learned miRNA→gene regulatory model (`mirna_hallmark`)

> **Purpose.** A self-contained brief to hand to a strong reasoning model (Fable 5) to *build with*.
> I want a rigorous, generative discussion of how to **learn** a miRNA→gene regulatory model from
> expression + evidence — the weights, how they aggregate/attribute, and how to validate them. I
> describe a hand-built heuristic we use today; treat it as **one baseline among several**, not a
> design to preserve and not the only thing to beat (see §0 on references). §1 is an extensive
> biology primer so you can reason about mechanism-grounded priors; §2 is the *exhaustive* data
> inventory (it corrects and extends our repo doc). §§3–8 pose the modeling questions with
> **developed candidate directions** — build them out, combine them, or replace them; they are
> starting points, not a menu to pick from. **§9 specifies the written deliverable I want back.**
>
> **Depth expected.** Treat this as a long, deep working session, not a quick answer. Reason at
> length; enumerate options and their tradeoffs *exhaustively*; work through the math *and* the
> mechanism; surface second-order effects and failure modes; prefer thoroughness over brevity —
> depth is the deliverable. **Favour presenting several viable directions and configurations, each
> with its tradeoffs and when it wins, plus a reasoned recommended default — over pruning to a single
> option.** **Refresh the literature from the web first, broadly** (see §1j) so the build reflects and
> improves on the state of the art rather than reinventing it.
>
> Repo companions (for me, not required reading for you): `docs/FORMULAS.md` (current heuristic as
> coded), `docs/LEARNED_REGULATORY_MATRIX_DESIGN.md` (an earlier regression sketch),
> `docs/DATA_SOURCES.md` (provenance — **known to be incomplete**, see §2).

---

## 0. Orientation, and what "better" means

We study **miRNA repression of the 50 MSigDB Hallmark programs in TCGA-BRCA**. A miRNA *arm* `m`
(e.g. `miR-29c-3p`) represses target genes `g`. I want a principled, *learned* model of (Q1) the
weight of each miRNA→gene interaction, (Q2) how weights aggregate on a gene and attribute back to a
miRNA, and (Q3) how to know the weights capture real targeting.

**Today's heuristic (a baseline, not scaffolding).** "Pressure" on a gene is a *fixed* score:

```
pressure(g,s) = Σ_{m∈R(g)}  expr_mult(m,s) · edge_w(m,g)

  edge_w(m,g) = log(1+evidence_score(m,g)) / D(m)     evidence weight ÷ promiscuity denominator
  D(m)        = hand-chosen per-arm "promiscuity" normaliser (total curated target mass)
  expr_mult   = softmax_z_logrpm = sm(m,s)·z(m,s)·logRPM(m,s)   hand-chosen; sm = gene-local softmax "share"
```

Every italicised piece (`D(m)`, the softmax share, the log-evidence weight) is a **hand-engineered
approximation** of effects — promiscuity, competition, relative abundance — that a fitted model
should produce endogenously. Known weaknesses we already observe: the softmax spends abundance twice
and is rank-invisible to the Spearman tests we score with; `edge_w` is near-flat within a gene so it
barely moves ranking; `D(m)` is an arbitrary correction. So this is a candidate to supersede — but:

**References the learned model is judged against — and which one matters depends on the question:**
1. **Raw miRNA abundance** — the null that weighting adds nothing beyond "how much miRNA is present."
   Often the *real* bar: does structure beat just-abundance?
2. **The curated fixed-M heuristic above** — does *learning* beat hand-set weights, or match them
   with fewer assumptions?
3. **Shuffled-evidence / permuted-support nulls** — is the prior doing real work?

It is **not** predetermined that beating the current heuristic is the goal. Sometimes the meaningful
comparison is against raw abundance; sometimes it is reproducing curated knowledge more parsimoniously
or extending it to edges it never had. Keep all references live; report against each.

One property worth keeping only as a sanity anchor: a good model should *nest* the heuristic (freeze
it at the prior → recover something like it), so anything learned reads as "the prior + what the data
insists on."

---

## 1. BIOLOGY PRIMER — mechanism you need to reason about priors and model form

This is deliberately thorough. Everything here is a lever on *what the weight means*, *what a sound
prior is*, and *what functional form is defensible*. Modeling implications are called out inline.

### 1a. Biogenesis, guide/passenger, and why "arm" is the unit
- Pol II transcribes a **pri-miRNA** (capped, polyadenylated). **Microprocessor** (Drosha + DGCR8)
  crops it to a ~70 nt **pre-miRNA** hairpin → Exportin-5 → cytoplasm → **Dicer** (+TRBP) removes the
  loop → a ~22 nt **duplex**. One strand (the **guide**) is loaded into an **Argonaute** (AGO1–4)
  to form **RISC**; the other (**passenger/star**) is usually degraded. Strand choice follows
  thermodynamic asymmetry (less-stable 5′ end) and a 5′-U preference, but **either arm (5p or 3p)
  can be the guide**, and the choice is tissue-dependent ("arm switching").
- ⇒ *Modeling:* the predictor unit is the **mature arm** (5p vs 3p resolved), not the hairpin. Guide
  vs passenger status is a prior on which arm is even competent to repress.

### 1b. The silencing mechanism and the seed
- The guide directs RISC to (mostly) the **3′UTR**. Specificity is dominated by the **seed** = guide
  nucleotides **2–7**. Canonical site types, in **decreasing efficacy**:
  **8mer** (seed match 2–8 + an A opposite position 1) > **7mer-m8** (2–8) > **7mer-A1** (2–7 + A1) >
  **6mer** (2–7, weak). **3′-supplementary** pairing (guide nt 13–16) adds modest affinity for some
  sites; fully non-seed "non-canonical" sites are mostly weak/ignorable.
- Effector step: AGO recruits **TNRC6/GW182**, which recruits the **CCR4–NOT** and **PAN2–PAN3**
  deadenylases → deadenylation → decapping → 5′→3′ decay, **plus** translational repression at
  initiation. **Guo et al. 2010:** at steady state the *majority* of repression is realised as
  **mRNA destabilisation** (visible at the mRNA level), with a translational-only residual.
- ⇒ *Modeling:* the **site-type hierarchy is a graded, quantitative magnitude prior** on `M`
  (8mer ≫ 6mer), independent of curated-MTI counts. And because most repression hits mRNA level,
  **mRNA is a legitimate primary readout** — but a translational-only component is invisible there
  and only appears in **protein / protein-vs-mRNA discordance** (ties to the observation-model fork, §6).

### 1c. Magnitude: miRNAs are fine-tuners, and combinatorial
- Individual target repression is typically **mild — most targets < 2-fold**. miRNAs buffer noise,
  sharpen thresholds, and confer robustness. Regulation is **many-to-many**: one 3′UTR carries sites
  for many miRNAs; one miRNA has hundreds of targets.
- ⇒ *Modeling:* a **linear-additive** main-effect term `Σ_m M_{g,m} x_m` is a reasonable *local*
  first order — but see 1e/1f for where additivity breaks (saturation, cooperativity). Mild effects
  also mean low SNR per edge → borrowing strength (priors, groups, low rank) is essential.

### 1d. Efficacy determinants (the TargetScan `context++` feature set — a real prior source)
Beyond site type, repression strength scales with, roughly in order:
- **Site multiplicity** — multiple sites in one 3′UTR repress more; sites **8–40 nt apart** act
  **cooperatively (super-additively)**.
- **Local AU-rich context** around the site; **site accessibility** (weak local RNA secondary
  structure / an "open" single-stranded region predicts efficacy — pairing competes with target
  self-structure); **positional bias along the 3′UTR** — functional sites are enriched **near the
  two ends**, depleted in the **~15 nt after the stop codon** and in the **interior of long UTRs**,
  with proximity to the poly-A signal also mattering. So *where* on the UTR a site sits, and the
  **UTR length/isoform** it sits in (§1f, APA), modulate its weight — not just its presence.
- **Seed-pairing stability (SPS)** and **3′-supplementary pairing**.
- **Target-site abundance (TA)** — a miRNA whose seed matches sites across *many* transcripts is
  spread thinner per target (a promiscuity/dilution effect).
- **Conservation** (branch-length score) — a prior on *functionality*, not affinity per se.
- ⇒ *Modeling:* `context++` gives a **continuous, sequence-based prior on entry magnitude** for
  *every* predicted edge, including ones with **no curated MTI** (orphans). Multiplicity/spacing is
  the biological trigger for **cooperative product terms**. TA is exactly what a promiscuity
  denominator crudely approximates — but here it is a per-*edge*, mechanistic quantity.

### 1e. Dose, the RISC-loading competition, and functional thresholds
- miRNA abundance spans **>4 orders of magnitude**; a handful of highly expressed miRNAs dominate
  **AGO loading**. Below a per-cell copy-number threshold (~hundreds of copies; Mullokandov 2012), a
  miRNA is **functionally inert regardless of site quality**. AGO/RISC is **limiting**.
- ⇒ *Modeling:* this is the mechanistic meaning of the **a·w scale gauge** (§3f): the effective
  weight is abundance-gated — low-abundance arms should be shrunk toward zero effect *even with a
  perfect site*. It also motivates a **global saturating nonlinearity** (finite RISC), not a purely
  linear map.

### 1f. Target-pool competition (ceRNA) and 3′UTR isoforms (APA)
- **Competing endogenous RNA (ceRNA)/sponge effect:** the total pool of available target sites
  dilutes a miRNA's per-target effect; a highly expressed target (or a sponge lncRNA/circRNA) can
  **relieve** repression on others. (Magnitude debated, but directionally real.)
- **Alternative polyadenylation (APA):** proliferating/cancer cells **shorten 3′UTRs**, deleting
  distal miRNA sites → targets **escape** repression. So **whether a site exists is sample-dependent.**
- ⇒ *Modeling:* ceRNA is competition among *targets* for a shared miRNA (opposite axis to seed
  competition among *miRNAs*). APA makes the **design-matrix support itself state/sample-varying** —
  the model can condition site presence on 3′UTR isoform.

### 1g. Seed families, redundancy, and context-dependence
- Arms sharing a seed hit **overlapping sites** → **functional redundancy**; they act as a group
  (TargetScan **families**). Same-seed members are near-substitutable at a shared site
  (**sub-additive**), which is the *opposite* of the cooperativity in 1d.
- A miRNA can only regulate **expressed** transcripts; its realised repertoire and even its
  onco/suppressor valence are **cell-type/context dependent**.
- ⇒ *Modeling:* group structure at the seed-family level is **biologically mandated**, not merely a
  statistical fix for collinearity (§3d). The model is inherently **cohort/tissue-specific**.

### 1h. The specific miRNAs that carry the Hallmark biology (for grounding priors)

| Arm / family | Key validated targets | Mechanism / role | Hallmark(s) |
|---|---|---|---|
| **let-7** family | KRAS/NRAS, MYC, HMGA2, LIN28B (double-neg loop) | tumour-suppressor; self-renewal, cell cycle | MYC/E2F proliferation, self-renewal |
| **miR-200** family + **miR-205** | ZEB1, ZEB2 (double-neg feedback → E-cadherin) | master epithelial-state switch | EMT (as suppressor) |
| **miR-29** (a/b/c) | COL1A1 & ECM collagens, DNMT3A/3B, TET, LOXL2 | matrix + DNA-methylation control | EMT/ECM, epigenetic |
| **miR-21** (oncomiR) | PTEN, PDCD4, TPM1, SPRY | proliferation, invasion, apoptosis resistance | PI3K/AKT, apoptosis (neg) |
| **miR-17~92 cluster** (17,18a,19a,20a,19b,92a) | E2F1 (feedback), PTEN, BIM/BCL2L11, CDKN1A/p21, THBS1 | MYC-induced "oncomiR-1"; **polycistron** | proliferation, angiogenesis, apoptosis evasion |
| **miR-155** | SOCS1, SHIP1/INPP5D | NF-κB/inflammation, immune | inflammatory response, IFN |
| **miR-34** (a/b/c) | CDK4/6, CCND1, BCL2, MET, SNAI1 | **p53 effector** | P53 pathway, apoptosis |
| **miR-221/222** | CDKN1B/p27, CDKN1C/p57, **ESR1** | ER-negativity, endocrine resistance | basal/TNBC, estrogen response (neg) |
| **miR-10b, miR-9** | HOXD10, CDH1/E-cadherin | pro-metastasis | EMT, invasion |
| **miR-145/143** | KRAS, MYC | tumour-suppressive cluster | proliferation (neg) |

Note the **miR-17~92 polycistron** and the **let-7 / miR-200 / miR-29 families** as the natural test
cases for family- and cluster-level aggregation (Q2) and for group structure (Q1).

### 1i. The state axis is biological, not just a covariate: healthy → NAT (field effect) → tumour
The regulatory landscape *itself* changes across states — abundance `x`, the available targets, and
plausibly the effective matrix `M` all move. This is the biology behind the §7 modeling question.
- **True-healthy breast (GTEx).** Baseline program of normal mammary epithelium/stroma: high
  tumour-suppressor miRNAs (let-7; miR-200/205 maintaining the epithelial state; miR-145/143;
  miR-125b), low oncomiRs. The only *bona fide* healthy reference.
- **NAT is NOT healthy — it is a field-effect / field-cancerisation tissue.** Histologically normal
  tumour-adjacent tissue carries subclinical molecular defects: shared stromal reprogramming,
  epigenetic drift, partial loss of tumour-suppressor miRNAs, altered immune/inflammatory tone. Using
  NAT as "normal" *understates* tumour change and can even invert directions. So GTEx is the healthy
  anchor and **NAT is its own state, never a stand-in for healthy**.
- **Tumour.** Global miRNA:RISC rewiring — **acquired** oncomiRs (miR-21, miR-155, miR-10b, miR-17~92),
  **lost** suppressors (miR-200 loss → EMT; let-7 loss → RAS/MYC de-repression), and shifts in
  RISC/AGO/TNRC6 machinery capacity that rescales *all* repression (§1e).
- **Why `M` itself may be state-dependent (not only `x`).** (i) **APA shortens 3′UTRs** in
  proliferating/tumour cells → site loss → the *support* of `M` contracts (§1f); (ii) the **expressed
  target repertoire** differs by state; (iii) **AGO/TNRC6/Dicer levels** change effector efficiency;
  (iv) ceRNA pools differ. A single fixed `M` across states is therefore an assumption to *test*, not
  a given — exactly the §7 question.
- **Acquired vs re-weighted.** The one robust cross-state signal we already have is *acquired*
  abundance (`max(tumour − healthy, 0)`); a state-aware model should separate "this arm is newly
  present" from "this arm's per-site weight changed."

### 1j. Literature to refresh from the web before building — go WIDE (do this as part of the work)
We seeded the current framework from a *partial* literature scan; **redo it broadly and deeply.** The
model you build should demonstrably reflect and improve on the state of the art across **every** facet
below, not only site scoring. Search, synthesise, and **cite**:
- **Silencing mechanism & magnitude** — current consensus on mRNA-destabilisation vs translational-
  repression split, effect-size distributions, fine-tuning/buffering, AGO/TNRC6 effector biology,
  functional dose thresholds.
- **Quantitative site efficacy & 3′UTR biology** — site-type hierarchy, **positional/location bias**,
  **site accessibility / secondary structure**, TargetScan `context++` and successors, APA /
  3′UTR-isoform effects, and modern ML site/efficacy predictors (deep-learning target models).
- **Seed families, cooperativity & competition** — cooperative spacing (8–40 nt) super-additivity,
  seed-redundancy sub-additivity, and the **ceRNA / sponge** literature (magnitude, when it matters,
  the criticisms) and AGO-pool saturation.
- **Expression-based regulatory *modeling*** — GenMiR++/GenMiR, sparse-regression / LASSO
  miRNA-target networks (Lu-type), comodule/factorisation methods (Theia, PIMiM, NMF lines), Bayesian
  and matrix-factorisation miRNA-target/disease work, reduced-rank / multivariate regression. For
  each: what to adopt, and what it *fails* at (magnitude priors, confounder control, family
  identifiability, cooperativity, an internal instrument).
- **Causal / identification methods** — miRNA / CNV **Mendelian-randomisation and instrumental**
  approaches, mediation, and confounder adjustment (purity / proliferation / composition) in tumour omics.
- **Protein-level & discordance** — miRNA effects at the proteome, protein-vs-mRNA discordance as a
  post-transcriptional readout (CPTAC-style), pSILAC / pulse-labelling studies.
- **State / cross-tissue & field-effect biology** — miRNA rewiring healthy→tumour, tumour-adjacent
  field cancerisation, cross-tissue targeting differences.
- **Breast-specific & Hallmark miRNA biology** — the §1h arms/families in breast cancer, subtype
  (PAM50 / TNBC) specificity, and anything that updates our target/role assignments.
- **Evaluation & overfitting** — how these networks are validated (held-out, orthogonal-layer,
  perturbation), and known failure modes / non-reproducibility.
Let this literature shape the object you commit to in §3a, the priors in §3b, and the evaluation in §5.

---

## 2. DATA & EVIDENCE — the exhaustive, corrected inventory

Any model is bounded by this. Our repo `DATA_SOURCES.md` is a good provenance map **but incomplete /
partly mischaracterised** — corrections flagged **⚠**. Exact counts are from the live matrices.

### 2a. Expression tensors (candidate inputs X / outputs Y), TCGA-BRCA spine

| Object | Shape / detail | Scale |
|---|---|---|
| **miRNA arm expression `X`** | **2,236 mature arms (MIMAT) × 1,207 samples**: **1,096 primary tumour, 104 NAT, 7 metastatic**. Detectable subset far smaller (detectability floor ~RPM ≥ 10 in ≥1 tumour) | `log2(RPM+1)` |
| **⚠ miRNA isomiR** | arm-resolved **isomiR** counts/RPM also available (`*mirna_isoform_*`) — sub-arm 5′/3′ heterogeneity; not in DATA_SOURCES | log2 RPM |
| **target mRNA `Y_mRNA`** | STAR counts + TPM, gene-level; **~4,384 Hallmark-union genes / 50 sets** (full transcriptome available) | log2 TPM |
| **⚠ target protein / phospho / acetyl** | CPTAC provides **proteome + phosphoproteome + acetylome** (not just proteome) on breast tumours | MS log-ratio |
| **target isoform / 3′UTR (APA)** | transcript-level expression + **PDUI** (Xia tumour, APAatlas normal) — site availability | — |
| **miRNA-locus copy number `CN_locus`** | ASCAT3 segment → locus overlap — a **natural perturbation of arm dose** (instrument, §5/§8) | integer CN |
| target-gene CNV, purity, HRD, proliferation, batch | confounders `C` | — |

### 2b. The state axis
`GTEx breast (true-healthy anchor) → NAT (adjacent, field-effect — NOT healthy) → tumour`, GTEx
quantile-normalised onto TCGA scale. Plus an in-situ→invasive (DCIS) axis from GEO. **TCGA is
cross-sectional** — no per-patient trajectory (rules out HMM/pseudotime unless external).

### 2c. Evidence / prior sources — **edge-level vs site-level** (this split drives §3e)

**Edge-level** (a whole miRNA→gene assertion):
- **miRTarBase** — curated MTIs graded by assay class (reporter / protein / RNA / perturbation /
  binding; functional vs weak). Our current edge table = **913,779 (miRNA, gene, hallmark) rows**;
  deduped high-evidence `(m,g)` = the model support.
- **⚠ TarBase v9** — **4,724,538 human interaction rows** (`data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz`),
  columns: experimental_method, **regulation direction (Positive/Negative)**, tissue, cell_line,
  confidence, PubMed id, **and 3′UTR site coordinates (start/end)**, microT score. **Absent from
  DATA_SOURCES.** A massive, method-graded, *directional*, partly site-resolved evidence layer — a
  first-class prior source we are not yet using.
- **ENCORI/starBase** — CLIP read-depth per edge (physical AGO-binding support).
- **⚠ Manakov chimeric eCLIP** — direct **miRNA–target duplex** reads (`data/external_cache/manakov_chimeric/`,
  many SRR runs). **Physical proof a specific pairing occurs**; validated ~60% of our high-value
  orphan edges. Absent from DATA_SOURCES.
- **PubMed breast-context** — literature PMID counts per edge (non-circular prior).

**Site-level** (sub-gene, 3′UTR coordinates — where cooperativity/competition live):
- **TargetScan** — full **`Predicted_Targets_Context_Scores`** + **`context++`** + Conserved/
  Nonconserved Family Info + the scoring scripts. Gives **site coordinates, seed type, conservation,
  and a continuous magnitude score for every predicted edge** (incl. orphans).
- **Seed families** (`miR_Family_Info`) — the grouping unit of §1g.
- **APA** (PolyASite, APAatlas normal-breast, Xia tumour ΔPDUI) — sample-varying site availability.
- **⚠ mirGeneDB** — **correction:** here it is **559 guide + 468 star mature *sequences* only**
  (`hsa_guide_mat.fas`, `hsa_star.fas`) — a curated guide/passenger **sequence catalog with NO
  experimental context**. Useful for arm assignment; do not treat as an assay/context source.

### 2d. External cohorts (held-out / orthogonal — never for fitting)
- **CPTAC** (prospective TMT + TCGA-105 iTRAQ): protein/phospho — an orthogonal *layer*.
- **⚠ SCAN-B** and **⚠ METABRIC** — large independent breast **mRNA** cohorts present in-repo (used
  in outcome modules); available as external target-side replication. Not in DATA_SOURCES' validation list.
- **Buffa / METABRIC-miRNA (GSE22216, n=210)** — independent invasive miRNA+mRNA.
- **EV & DCIS GEO corpus** — direction/state slotting (biomarker lane, kept separate).
- **CCLE breast lines** — miRNA/CN/expression for cell-line concordance.

---

## 3. Q1 — the weight of a miRNA→gene interaction (build the learned model)

**My framing:** move from static evidence weights + hand-coded promiscuity/competition → to
**learning** weights from expression under **evidence priors** (quantitative *and* qualitative), with
competition, cooperativity, and same-seed ambiguity emerging **axiomatically** inside the model.
Below are **developed candidate directions** — build on them, combine, or replace.

### 3a. Candidate object (react, then commit to something)
A non-negative repression matrix `M` in a multivariate multiple regression:
```
Y = − X·Mᵀ + C·B + E ,   M_{g,m} ≥ 0
```
`M_{g,m}` = repressive strength of arm `m` on gene `g` (our `edge_w`, learned). Not bilinear in
(x,y) (`y` is the outcome); bilinearity re-enters *legitimately* via low-rank `M ≈ U Vᵀ`
(reduced-rank regression / latent comodules — Theia/PIMiM/NMF family), via a **feature model for the
entries** `M_{g,m} = φ(g)ᵀ Θ ψ(m)` (to reach orphan edges from sequence features — a
factorization-machine reading of `context++`), and via **cooperative product terms** in `x` (§3e).
**Is regression even the right object?** Weigh RRR vs feature/FM models vs hierarchical Bayes vs a
mechanistic-link GLM, given `S≈1,100`, `A≈2,236`, `G≈4,384` (`p ≫ n`). Lay out the **viable objects
and configurations** with their tradeoffs and identifiability implications, and recommend a **default**
(saying when an alternative is preferable) — rather than collapsing to one prematurely.

### 3b. Where the prior enters — first the inclusion decision, then four magnitude channels

**The inclusion decision is itself a first-class modeling choice, not a preprocessing step — and the
grading rule is part of what you are designing, not a fixed input.** Our support is graded. Our
**current *working* definition of high-evidence (HE)** is: an edge with ≥1 *functional* MTI study
(reporter / protein / perturbation) **AND** an independent low-throughput validation; below HE is an
evidence-score tail; below *that* are **orphan** edges predicted by TargetScan/CLIP with **no curated
MTI at all** (often with strong `context++` / chimeric-eCLIP support). **Treat this HE definition as
provisional and very likely improvable** — it is a miRTarBase-tier rule that should probably be
redefined by folding in TarBase-v9 method + direction labels, CLIP / chimeric-eCLIP physical depth,
`context++` magnitude, and breast-context literature. "How much evidence, of what quality and
directness, counts as includable" is itself a design question. Given some grading, the mechanism is
still a fork:
- **Hard-restrict to a high-evidence set** — set `M_{g,m}=0` off it. Cleanest and most identifiable,
  but discards orphans and most of the matrix, and throws away the continuous `context++` /
  TarBase-v9 / CLIP signal on excluded edges (no discovery).
- **Soft, evidence-graded inclusion** — a **spike-and-slab** inclusion probability
  `π_{g,m}=f(evidence)`: high for strong edges, low-but-nonzero for orphans with real
  sequence/physical support, so a weak edge enters **only if the expression data insists** while the
  prior holds the noise floor down.
Lay out the viable options for **both** the support rule (what is includable, and how finely graded)
**and** the mechanism (hard gate vs soft prior vs a blend), with their tradeoffs and a recommended
default — keeping the strong alternatives in play rather than pruning to one; neither is settled.

**Then the magnitude** of an admitted edge enters in the remaining channels — the distinctive claim
being that our **quality-weighted, assay-directness-graded, breast-audited** evidence enters as
*magnitude*, not merely as a 0/1 mask:
1. **Penalty weight** — adaptive L1 `λ_{g,m}=λ/w_prior`: strong edge shrinks little, weak edge heavily.
2. **Prior mean** — `M_{g,m} ~ N⁺(μ_{g,m}, τ²)`, `μ` = curated/`context++` magnitude.
3. **Sign/box** — `M_{g,m} ≥ 0`.
(The inclusion prior above is the fourth — the *support* channel.)

The prior is now **multi-source and heterogeneous** (§2c): curated-MTI assay tiers, **TarBase v9
method+direction**, CLIP/chimeric physical depth, TargetScan `context++` (continuous, per-edge, incl.
orphans), breast literature. **How do you fuse these into a single inclusion probability + prior mean
+ prior strength per edge** — especially reconciling edge-level assertions with site-level
`context++`, and directional TarBase evidence with undirected CLIP binding? And how much may data
override the prior while staying non-circular (§5)?

### 3c. Sparsity & identifiability
`p ≫ n`; learnability comes from structure (evidence-sparsity, low rank, groups, primed interactions).
Specify the regularisation/structure and how its strength is chosen **out-of-fold**, not by in-sample fit.

### 3d. Same-seed / family problem — represent, don't fake-resolve
Paralogous same-seed arms co-vary (within-family ρ ≈ 2× between-family; hubs 20–37% redundant) →
**near-collinear columns of X**. Plain LASSO keeps one member arbitrarily, unstably across folds.
Directions: **group-LASSO at the seed-family level** (in/out together), **collapse to one family
predictor**, or **hierarchical partial pooling** toward a family latent that reports member weights
*with* their unidentifiability. The identified estimand is likely **family→gene**; the per-arm split
is a within-group attribution, not identified. (Biologically mandated, §1g — not just a numerical fix.)

### 3e. Competition & cooperativity — axiomatic, and needing the edge/site split
Three distinct couplings, different signs, different terms, different priors:

| Phenomenon | Site relationship | Effect | Model term | Prior trigger |
|---|---|---|---|---|
| **Seed redundancy / competition** (miRNA↔miRNA) | same/overlapping site | **sub-additive** | group structure (§3d) + optional per-site saturation | shared seed / overlapping **coords** |
| **Cooperativity** (miRNA↔miRNA) | *distinct* sites 8–40 nt apart | **super-additive** | product term `+ M_{g,(m,m')} x_m x_{m'}` (bilinear in x) | TargetScan/TarBase **site coords** in window |
| **ceRNA / pool competition** (target↔target) | shared finite RISC | targets relieve each other; saturating | global/cohort saturating nonlinearity + target-pool term | total load / AGO capacity; target abundance |

A plain linear `Σ_m M_{g,m}x_m` represents **none** of these. Key points to build:
- We can't fit all `A(A−1)/2` products — **site-level evidence primes which interactions exist**, so
  the prior does double duty: main-effect support (edge-level) **and** interaction support (site
  coords). Show the fusion concretely.
- Does modeling competition explicitly make the promiscuity denominator `D(m)` **obsolete** — i.e.
  does promiscuity become a *consequence* of one arm competing across many genes for finite RISC,
  rather than a hand-set denominator?

### 3f. The a·w scale gauge (central; you must resolve it)
`contribution = abundance · weight` is multiplicatively non-identifiable — rescale `x_m` up and
`M_{g,m}` down and predict the same `Y`. So **"the weight" is undefined without an anchoring gauge**;
`D(m)` is *implicitly* one. Biologically the gauge is the **RISC-loading share** and the
**functional-threshold** of §1e (low-abundance arms are inert regardless of site). Give a principled
gauge (per-arm abundance standardisation? row/column norm constraint on `M`? fix abundance scale, let
`M` carry magnitude? an abundance-gated *effective* weight?) and show what it makes "a strong edge" mean.

---

## 4. Q2 — aggregation on a gene, attribution to a miRNA

**My framing:** how do edges **sum on a gene**, and how do you **attribute** the aggregate back to a
single miRNA — with a gene↔miRNA-centric spectrum of extensions.
- **Gene aggregate** = row of `M` applied to `x`. Choices: sum vs mean; signed/positive/absolute
  coherence; pruning weak regulators; RISC-gating the aggregate.
- **Attribution to one arm.** Our current softmax "share" is a **leftover I want reconsidered.** In a
  *learned* model what is the right attribution object — a Shapley/LMG decomposition of the fitted
  aggregate, or the fitted `contribution(m,g,s)` directly? Beware the trap we've hit:
  **identity/attribution (who owns the repression, abundance-removed) ≠ magnitude/force (how much
  total)** — a share summing to 1 answers "who," not "how much." Which does each downstream use need,
  and should the model expose both?
- **The gene↔miRNA spectrum — all one `M`, different blocks:**
  - *miRNA→its target set*: a **column** of `M` aggregated over `g` (one arm's total role).
  - *gene-centric layer impact*: total pressure the **whole miRNA layer** applies to each gene.
  - *intermediates*: a **seed family**, a **polycistronic cluster** (miR-17~92), or any arm-set, onto
    a gene or a Hallmark. Does the §3d family/group grouping give this aggregation for free?

---

## 5. Q3 — evaluation: behaviors, tests, references, no overfitting

**My framing:** how do I know the weights capture real targeting, not covariation or overfitting?
- **References** — see §0: raw abundance, curated fixed-M, shuffled-evidence nulls; which is the
  meaningful bar depends on the claim.
- **Tests / expected behaviors:**
  - **Out-of-fold coupling** — does the learned aggregate anti-correlate with target expression on
    held-out patients better than raw abundance? (partial-Spearman adjusting `C`).
  - **CN "Mendelian randomisation" / instrument** — miRNA-locus CN naturally perturbs arm dose,
    largely not caused by the target. 2SLS: `x_m = γ·CN_locus + Cδ`; `y_g = −M·x̂ + Cβ`. Report
    first-stage F (weak instrument), restrict to focal segments (pleiotropy). Strongest "is it
    causal" test; almost no expression-network paper carries an internal instrument.
  - **Orthogonal-layer OOD** — validate at the **protein / discordance** layer (CPTAC) the fit never
    saw, and in an **independent cohort** (SCAN-B, METABRIC, Buffa). Optionally **chimeric-eCLIP**
    physical binding as an orthogonal support check on learned edges.
  - **Known-driver recovery** (PTEN, GATA3, ESR1, ZEB1, the §1h set) as face validity.
- **Held-out discipline & circularity guardrail** — the sequence/curation prior fixes support + prior
  mean; expression estimates magnitudes; **never** test coupling on the same (X,Y) used to fit. Keep
  learned-M and curated-M side by side; report disagreements as findings.
- **Identification is not only here.** The confounders that make correlation ≠ targeting live inside
  the *fit* (§3's `C·B` block and the instrument), not just the check. **What you can identify bounds
  what you can learn** — Q1 and Q3 are one problem.

---

## 6. Cross-cutting A — the observation model: predict *what*? (present all, decide with argument)

Most load-bearing choice, currently implicit. miRNA acts post-transcriptionally:
`target = transcription − repression(miRNAs) − other`. Regressing weights on steady-state **level**
makes transcription a confounder (a TF co-driving a miRNA and its target masquerades as a *positive*
weight). Options — all on the table, with tradeoffs; §1b says mRNA carries most (not all) of the effect:
1. **mRNA level** — most data; transcription-confounded; needs a transcription proxy in `C`.
2. **mRNA Δ vs reference** (NAT / GTEx / subtype centroid) — differences out shared baselines.
3. **mRNA residualised on a transcription proxy** — isolates the post-transcriptional part.
4. **Protein level** (CPTAC) — closer to function; smaller n; MS batch.
5. **Protein-vs-mRNA discordance** — arguably the *cleanest* post-transcriptional signal (mRNA
   present, protein down = the translational residual of §1b); needs paired CPTAC.
Lay out the viable choices with their tradeoffs, recommend a primary `Y` and a validation `Y` (noting
when another would be preferable), and say how the choice interacts with `C` (§3) and the state axis (§7).

---

## 7. Cross-cutting B — the H→NAT→Tumour axis: how does state enter? (open)

Three states (GTEx-healthy → NAT-field-effect → tumour), no trajectory. Argue among:
- **Separate `M` per state** (`M_H`, `M_T`) — cleanest interpretation, weakest per-state power; and
  how do you compare two separately-*gauged* matrices (§3f)?
- **Pooled fit, state as covariate/interaction** — one `M` + state offsets, or low-rank modulation
  `M_state = M_0 + Δ_state`.
- **State-differenced** `ΔY = −M·Δx` across the chain — better-identified (baselines cancel) or noisier?
- **Nest one into the next** — is tumour `M` a perturbation of healthy `M` (shared support, changed
  magnitudes)? Should the healthy fit be the *prior* for the tumour fit?
- Does the answer depend on the §6 observation-model choice (differencing pairs with "Δ vs reference")?

---

## 8. What I want out of this discussion

Please **build**, at length and exhaustively — treat each item below as a section to work through in
depth (math + mechanism + tradeoffs + failure modes), not a bullet to tick. For each, present the
**several viable directions / configurations** with their tradeoffs and *when each wins*, and give a
**reasoned recommended default**: this is not a shallow survey, but do **not** prune to a single option
prematurely — keeping strong alternatives alive is valuable. Start by **refreshing the literature
broadly (§1j)** and let it shape your answers.
1. **Lay out the viable model objects** (§3a) for `p≫n` with tradeoffs and a recommended default —
   plus what you'd *not* do, and why.
2. **Resolve the observation-model fork** (§6) — primary `Y` + validation `Y`, with the confounding argument.
3. **Give a principled a·w gauge** (§3f) and its consequences for cross-arm comparison.
4. **Specify the multi-source prior fusion** (§3b) and how **edge-level + site-level** evidence jointly
   prime main-effect *and* interaction (competition/cooperativity/ceRNA) support (§3e, §2c) — and
   whether explicit competition retires the promiscuity denominator.
5. **Handle family identifiability** (§3d) — group vs partial-pooling vs collapse; name the identified estimand.
6. **Design the evaluation** (§5) so a learned weight can beat *the right reference* out-of-fold and on
   an orthogonal layer, without circularity.
7. **Flag anything wrong, conflated, or missing** in this framing, including whether a richer source
   (e.g. TarBase v9's directional/tissue labels, isomiR resolution) changes the model you'd build.

Assume I'm comfortable with the statistics and the biology of §1 — be technical, mechanism-grounded,
and concrete. The current `softmax`/`D(m)`/`edge_w` machinery is **one baseline**, not a design to
preserve; the developed directions above are **starting points to build from**, not the only options.

---

## 9. The deliverable — what to produce, and in what form

Produce a **standalone written design report** I can read, audit, and act on: long, structured, and —
above all — **explaining what you decided, how it works, and why** (reasoning first, not just
conclusions). I should be able to reconstruct your logic and challenge any single step. Structure it:

1. **Executive summary** (≤1 page) — the model you'd build, the observation target `Y`, the prior
   scheme, the evaluation, and the 3–5 decisions that matter most, each in one line.
2. **Literature synthesis** (from §1j) — what the current state of the art establishes across the
   facets you searched, **with citations**, and explicitly what you **adopt** vs **improve on**. Flag
   anything that contradicts the biology (§1) or data (§2) assumptions here.
3. **Decision log — one entry per major choice**, in this exact shape so each is auditable:
   - **Decision:** what you chose.
   - **How:** the mechanism / math (equations where relevant).
   - **Why:** the reasoning from biology (§1), data (§2), and identifiability.
   - **Alternatives kept in play:** the other viable options/configurations, their tradeoffs, and
     **when each would be preferred** — prune only what is genuinely dominated, and say why.
   - **Assumptions & risks:** what must hold, how it fails, how you'd detect the failure.
   Cover at least: model object (§3a); observation target `Y` (§6); the inclusion/HE rule **and**
   hard-vs-soft mechanism (§3b); multi-source prior fusion (§3b/§3e); the a·w gauge (§3f); family /
   identifiability (§3d); competition / cooperativity / ceRNA terms (§3e); the state-axis treatment
   (§7); aggregation & attribution (§4); evaluation & references (§5, §0).
4. **The model specification(s)** — a full generative/estimation statement for the **recommended
   default** (objective, penalties, constraints, priors, interaction support, confounder block),
   precise enough to implement, every symbol defined; plus the **viable variants kept in play**,
   specified enough to tell apart and to know when each is preferred. Note explicitly where it
   **nests** the current heuristic (§0).
5. **Evaluation & falsification plan** — the exact out-of-fold / instrument / orthogonal-layer tests,
   the reference each is scored against (§0), the **pre-named** success criteria, and what result would
   make you **abandon** the learned model for raw abundance or curated-M.
6. **Phased build plan** — an ordered path from a minimal viable version to the full model, each phase
   with a decision gate (a named held-out / orthogonal win required to proceed). Cheapest, most
   diagnostic step first.
7. **Open questions & what would change the answer** — where you are uncertain, what data or analysis
   (e.g. wiring TarBase v9, isomiR, CPTAC protein) would resolve it, and how a different answer to §6
   or §7 would ripple through the rest.

Throughout, **make the reasoning legible**: show the tradeoffs, quantify where you can, and mark every
place you are *inferring* vs *asserting*. I value a well-argued "here is why this is hard and here is
the principled choice, here are the viable alternatives and when each wins, and here is my recommended
default" over a confident but shallow single pick — depth, breadth of options, and traceability are the deliverable.
