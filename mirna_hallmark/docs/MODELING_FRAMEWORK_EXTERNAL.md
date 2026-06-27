---
title: "Modeling microRNA Repressive Pressure on the Hallmark Programs of Cancer"
subtitle: "A conceptual and computational framework, applied to breast cancer"
---

## Abstract

A microRNA (miRNA) represses many target genes, and a gene is repressed by many
miRNAs, so the biologically meaningful object is not a single binding event but a
**field of repression** acting on a cell. We formalize this field as *pressure* — an
evidence-weighted, expression-scaled potential for repression — and separate it
cleanly from *coupling*, an observed, confounder-adjusted anti-correlation that asks
whether the pressure is actually **realized** in a tumour. We model both quantities at
five nested resolutions (the whole program universe, a single program, a miRNA's
target set, a gene, and a single interaction), across the healthy→tumour trajectory,
and we test realization at the messenger-RNA layer, the protein layer, and in
independent patient cohorts. The framework is substrate-agnostic; here it is
implemented on breast cancer using the fifty MSigDB Hallmark programs as the gene
universe. This document presents the concepts, the formulas with their rationale, and
the validation logic.

> **Central idea.** Pressure is a *potential* (what the wiring and abundance could
> deliver); coupling is *realization* (what is observed after confounders are removed);
> neither is causation (which requires a perturbation assay). Keeping these three rungs
> distinct — and joining potential with realization rather than conflating them — is the
> framework's organizing principle.

---

## 1. The modeling problem and stance

### 1.1 What "pressure" is, and why it is a potential

A miRNA loaded into the RISC/Argonaute complex base-pairs a target mRNA and represses
it (predominantly by destabilizing the message, with a smaller translational-block
component). Because the regulatory graph is densely many-to-many, we quantify the
aggregate downward force the miRNA layer places on a gene — its **pressure** — built
from three ingredients:

- **which interactions exist** (curated and predicted edges),
- **how strong each interaction is** (assay evidence),
- **how abundant the miRNA is** (expression).

Pressure is deliberately a *potential* quantity. It does not, by itself, prove that
repression occurred in a given tumour. That harder question is answered separately by
**coupling**: a confounder-adjusted anti-correlation between the miRNA *dose* (its
abundance at single-regulator resolution, or the aggregate pressure where regulators are
pooled — §6) and its target's expression. The decision to hold these two as distinct
measurements, and then to report their **join**, organizes everything that follows.

### 1.2 Why programs, why breast cancer

The fifty Hallmark gene sets are a curated, largely non-redundant atlas of the coherent
biological *programs* of cancer (proliferation, epithelial–mesenchymal transition,
hypoxia, apoptosis, p53, interferon response, and so on). Using them as the gene
universe makes every result interpretable as "pressure on a named program," and lets
pressure be rolled up and down a biological hierarchy.

Breast cancer is the implementation substrate because it offers the deepest matched
multi-omics — miRNA, messenger RNA, copy number, DNA methylation, chromatin
accessibility, and two independent mass-spectrometry proteomes — together with
well-defined molecular subtypes and an independent replication cohort. Nothing in the
pressure, coupling, or architecture machinery is breast-specific; only the
program-polarity prior and one or two confounders are cancer-flavoured.

### 1.3 Design principles

These recur throughout and are stated once.

1. **Static evidence × dynamic abundance.** The interaction weight encodes curated and
   sequence evidence, which does not change between samples; *all* sample- and
   state-to-state variation enters through miRNA expression. This separates "what is
   biologically possible" from "what this tumour is expressing."
2. **The prior is non-circular.** Interaction weights, the evidence filter, the
   program-polarity prior, and gene-role annotations never use the patient expression
   data they are later tested against, so they cannot leak into the findings.
3. **Adjusted and raw, gated and ungated, always.** Every coupling carries both its raw
   and its confounder-adjusted value; every availability-gated pressure carries its
   ungated twin. Differences are read as sensitivity, not as the headline.
4. **Honest rungs of evidence.** Potential is not realization is not causation; each
   claim is labelled by the strength of support behind it.

### 1.4 Data, cohort, and stratification

All quantities are computed on **TCGA breast cancer (TCGA-BRCA)**: primary tumours with
matched small-RNA sequencing and messenger-RNA sequencing, plus copy number, DNA methylation,
and chromatin accessibility where assayed; the 104 normal-adjacent samples provide the
intermediate tissue state (§8). The PAM50 **Normal-like** class (36 tumours; normal-tissue
admixture) is **set aside cohort-wide**, leaving **1,065 retained tumours**; the realized-coupling
tests run on the **~1,041** with complete miRNA + mRNA + covariate data. Over this cohort the
interaction universe is **721 miRNA arms, 5,108 curated edges, and 1,410 target genes** drawn
from the **50 MSigDB Hallmark programs** (4,384 genes in their union).

Tumours are stratified throughout by the **PAM50 molecular subtypes**, the field-standard
breast taxonomy — and most relative quantities are *recomputed within each subtype*, not
sliced from cohort numbers (§2):

| PAM50 subtype | n (retained) | Character |
|---|---:|---|
| Luminal A | 507 | hormone-receptor-positive, low-proliferation |
| Luminal B | 264 | hormone-receptor-positive, higher-proliferation |
| HER2-enriched | 93 | HER2-driven |
| Basal-like | 197 | largely triple-negative, high-proliferation |
| Normal-like | 0 (36 set aside) | low-alteration, normal-tissue-like — **set aside cohort-wide** (normal-tissue admixture); excluded from the realized analyses |

The external and reference datasets the framework draws on:

| Source | Role |
|---|---|
| **TCGA-BRCA** | primary cohort — miRNA, mRNA, copy number, methylation, accessibility, subtype |
| **GTEx** (breast, 346 samples) | true-healthy miRNA reference for the acquired axis (§8) |
| **DIANA-miTED** | independent healthy miRNA atlas, used to repair healthy-reference coverage |
| **CPTAC** breast proteomes — a same-patient set (n = 105) and an *independent* cohort (n about 101) | protein-layer validation (§7) |
| **Buffa 2011** (210 paired primary breast tumours, no TCGA overlap) | independent-cohort replication of coupling (§9) |
| **miRTarBase** | curated functional miRNA–target interactions — the experimental core |
| **TargetScan** | sequence (context++) target predictions |
| **ENCORI** | Argonaute-crosslink / degradome occupancy |
| **miRGeneDB** | guide-versus-passenger arm curation |
| **OmniPath** | directed, signed gene–gene interactions for program topology (§5) |
| **COSMIC Cancer Gene Census / OncoKB** | 232 curated drivers (114 oncogene, 106 tumour-suppressor, 12 dual) for gene-role |
| **DepMap** (CRISPR essentiality, 96 breast cell lines) | continuous, all-gene dependency weight (§5) |
| **Lambert 2018 TF census** (1,639 human transcription factors) | master-regulator identity within programs (§5) |

All sample counts are approximate and refer to the samples passing each analysis's join
and quality filters.

---

## 2. The hierarchy of resolution

The same pressure object is modeled at five nested resolutions, each contained in the one
before it:

> **Universe > Program > Target-set > Gene > Edge**

— that is, the whole program universe (all 50 programs, 4,384 genes of which 1,410 are
under pressure, 721 miRNA arms, 5,108 edges); one program; the genes a single miRNA arm
or seed family represses; one gene with its full crowd of regulators; and one
miRNA-to-gene interaction. Zooming in localizes a signal, zooming out contextualizes it;
each level is one row of the table below.

| Resolution | Unit | Question | Pressure quantity | Realization quantity |
|---|---|---|---|---|
| **Edge** | one interaction | Does this miRNA repress this gene? | contribution $c(m,g,s)$ | per-edge partial correlation |
| **Gene** | a gene + its regulators | How much pressure lands here; who carries it? | $\sum_m c$; gene-side share | gene-level net-repression |
| **Target-set** | one arm's / family's targets | Does the arm repress its whole module together? | per-arm mass / share | module composite correlation |
| **Program** | one program's genes | Net effect on the program's *output* | architecture score | fraction of program coupled |
| **Universe** | all programs | Cross-program structure, subtype programs | unsupervised factors | fraction of programs coupled |

Two features make this more than "sum the small things":

- **Scope recomputation.** Every mean, standard deviation, median, and softmax that
  feeds pressure is recomputed *within the samples of the current scope* (whole cohort,
  one subtype, or one tissue state). A subtype table is therefore not the cohort numbers
  re-sliced — relative abundances are re-derived inside the subtype.
- **Aggregation carries meaning.** A gene's pressure is the **sum** over its regulators;
  a program's effect is a **topology-weighted, sign-aware propagation** over its genes,
  not a sum. The roll-up rule is part of the model.

---

## 3. The interaction universe

The edges are the skeleton over which pressure is built. The universe is **primarily
experimental, with a controlled computational extension.**

### 3.1 Experimental core, computational extension

| Layer | Source | Contribution | Role |
|---|---|---|---|
| **Experimental core** | **miRTarBase** — curated functional miRNA–target interactions (reporter assays, protein blots, qPCR, perturbation, crosslinking) | literature-validated edges | the **spine** — only these enter the default pressure budget |
| Sequence prediction | **TargetScan** — seed-match + conservation (context++) scores | binding-site strength, study-count-independent | extension, discovery screens, hybrid sensitivity |
| Crosslink depth | **ENCORI** — Argonaute-crosslink / degradome occupancy | physical-binding depth | a depth boost on the subset of pairs with such data |

The default spine is the curated, experimentally supported set (5,108 edges over 721
expressible miRNA arms and 1,410 target genes). Sequence-only ("orphan") predictions do
**not** enter the default budget — they are screened separately for discovery (§9) so
that prediction noise can never inflate a headline.

### 3.2 The high-evidence gate

An interaction qualifies for the spine if it has functional support from a *direct*
assay — at least one functional study, including at least one reporter or protein
assay. This gate is a binary membership filter (which edges are credible enough to
model); the weight below is the continuous strength of those that pass.

### 3.3 Interaction weighting

The weight is an **evidence numerator** divided by a **promiscuity denominator**.

#### Evidence numerator

The numerator does not use a raw study count, which inflates with publication volume and
treats a proximity hit like a functional proof. Instead it is a quality-weighted,
log-damped sum over assay classes. Writing $\ell(n)=\log(1+n)$ for the study count of
one assay class,

$$
\mathrm{ev}(m,g) \;=\; 3.0\,\ell(n_{\text{reporter}}) + 2.5\,\ell(n_{\text{protein}})
+ 1.5\,\ell(n_{\text{rna}}) + 1.5\,\ell(n_{\text{pert}}) + 0.5\,\ell(n_{\text{bind}})
\;+\; 0.3 \!\!\sum_{\text{weak classes}}\!\! \ell(\cdot).
$$

Three rationales are built in:

- **Assay directness** — a reporter assay that *demonstrates* repression earns six times
  a crosslink hit that only shows *proximity*. This single correction stops
  physical-occupancy data from masquerading as functional proof.
- **Log damping within each class** — the tenth confirming paper adds little over the
  third; $\ell(\cdot)$ gives diminishing returns and defuses publication-count bias.
- **A weak-evidence discount** ($0.3$) — curator-flagged weak evidence still counts, at
  thirty percent.

These class weights are a *transparent prior, not a fit*. Perturbing all of them — random
±25% jitter and a systematic ±50% on each in turn — leaves the realized-coupling headline
essentially unchanged (the count of strong-negative genes stays within ~10% of baseline,
coefficient of variation ~0.06), so no result rests on the specific numbers.

Where physical-binding depth is available it is added (with overall weight $0.5$):

$$
\begin{aligned}
\mathrm{enc}(m,g) &= 2\,\ell(\text{crosslink depth}) + \ell(\text{degradome depth}) \\
&\quad + 0.5\,\mathbb{1}[\text{three predictors agree}] + 0.25\,\ell(\text{pan-cancer depth}), \\
\mathrm{ev}_{\text{eff}}(m,g) &= \mathrm{ev}(m,g) + 0.5\,\mathrm{enc}(m,g).
\end{aligned}
$$

#### Promiscuity denominator

A miRNA that targets everything must not dominate any one gene merely by being
promiscuous. So the base strength $\ell(\mathrm{ev}_{\text{eff}})$ is divided by the
arm's total evidence mass:

$$
D(m) = \ell\!\Big(\textstyle\sum_{g'} \ell\big(\mathrm{ev}(m,g')\big)\Big),
\qquad
w(m,g) = \frac{\ell\big(\mathrm{ev}_{\text{eff}}(m,g)\big)}{\max\big(D(m),\,\varepsilon\big)}.
$$

**What the denominator does, and does not, affect — measured, not asserted.** Because
$D(m)$ is a single per-arm scalar, it is invisible to two of the three things one might
think it tunes, and we verified both directly rather than claim a coupling benefit it does
not have. (i) *Realized coupling* is essentially unchanged: in patient cross-validation the
denominator is coupling-neutral — the un-normalized weighting is in fact marginally ahead on
out-of-sample coupling, and the spine's choice sits within ~0.001 of it (and the candidate
ranking reproduces perfectly across held-out folds, so nothing here is overfit). (ii)
*Specificity* is mathematically invariant to it: $\mathrm{spec}(m,g)$ is a within-arm ratio
in which $D(m)$ cancels exactly (verified: rank-correlation $1.0$). Where the denominator
*does* act is the one place it should — **cross-arm attribution**: it changes which arm is
credited as a gene's dominant regulator for about a tenth of genes and measurably
de-concentrates the hub structure (a handful of promiscuous high-evidence arms otherwise
dominate disproportionately many genes; the normalization spreads that dominance from ~43 to
~56 arms). So the denominator earns its place as **promiscuity control on the attribution
axis**, not as a coupling lever — and penalizing by total *credible evidence mass* rather
than raw target count is a deliberate middle (a pure degree penalty de-concentrates slightly
more), not a uniquely optimal choice.

### 3.4 The two directions of many-to-many

The two directions of the graph are handled by two different mechanisms; conflating them
is a common error.

- **One miRNA, many genes** is controlled by the denominator $D(m)$ above and by
  *specificity* (§5.4), which measures how concentrated an arm's budget is on a single
  gene (focused target versus promiscuous hub).
- **Many miRNAs, one gene** is handled by gene-level aggregation ($\sum_m c$) and by
  *gene-side share* — the fraction of a gene's total incoming pressure that one arm
  carries, with a rank (rank 1 = the gene's dominant regulator). A heavily targeted gene
  may have dozens of regulators; gene-side share lets us name the dominant one while
  recognizing that the same edge is a trivial slice of that miRNA's own portfolio. The
  arm-side and gene-side shares answer opposite questions and both are kept.

A gene that belongs to several programs generates one edge row per program for
program-level views; cohort-level counts deduplicate per interaction so a multi-program
gene is not over-counted.

### 3.5 The pan-context caveat and non-circular refinements

The curated evidence carries no tissue field — an interaction is "high-evidence" if it
was validated in *any* cell line. In this cohort only about an eighth of high-evidence
edges have any breast-specific validating paper. Rather than hide this, two refinements
quantify it, and both are kept strictly separate from the pressure budget so they remain
non-circular:

- a **tissue re-ranking** that re-orders a gene's regulators by adding breast-literature
  weight, recovering textbook tissue-specific assignments; and
- a **multi-axis prior** that composes evidence with sequence strength, tissue
  literature, *negative* (contradicting) evidence, and guide-versus-passenger arm
  status, each as a separate, visible term.

Both are **audit and nomination overlays, never folded into the default budget**, for
three reasons. First, *non-circularity*: the budget is later tested against patient
expression, so any term we might also want to validate against that same expression must
stay out of the budget, or the test becomes circular. Second, *unvalidated utility*: a
tissue re-ranking can only be shown to be an *improvement* against an independent
breast-specific miRNA dataset, and no suitable one exists (§10), so we do not trust it
enough to redefine the budget. Third, *they do not move the headline*: the refinements
reorder which specific arm is credited for a given gene but leave the hub- and
family-level structure unchanged — so their value is per-target nomination, not a change
to the aggregate field. Keeping them as overlays lets each be inspected and tested on its
own rather than silently baked in.

---

## 4. Expression and pressure

### 4.1 The dynamic ingredient

All sample- and state-to-state variation enters through expression. Let
$x_{m,s}=\log_2(\mathrm{RPM}+1)$ be arm $m$'s abundance in sample $s$, with in-scope
mean $\mu_m$, standard deviation $\sigma_m$, and median $\tilde x_m$, and let $h_m$ be
the arm's median in healthy tissue. Three building blocks are combined into modes:

$$
\underbrace{z_{m,s}=\frac{x_{m,s}-\mu_m}{\sigma_m}}_{\text{within-tumour dynamics}}
\qquad
\underbrace{\mathrm{sm}_{m,s}=\frac{e^{\,x_{m,s}-\tilde x_m}}{\displaystyle\sum_{m'\in R(g)} e^{\,x_{m',s}-\tilde x_{m'}}}}_{\text{relative share among a gene's regulators}}
\qquad
\underbrace{\mathrm{dev}_{m,s}=x_{m,s}-h_m}_{\text{deviation from healthy}}
$$

The standardized score $z$ measures how far an arm sits from *its own* typical level —
so a stably high arm sits at its own mean ($z$ near zero) and is invisible under $z$ (intended for coupling,
where dynamics are the signal; wrong for attribution). The softmax $\mathrm{sm}$ is the
fraction of a gene's regulator "abundance budget" an arm holds, and is where the *crowd
of regulators* enters. The healthy deviation $\mathrm{dev}$ has no standard-deviation
division, so a uniformly tumour-elevated oncomiR is not gated to zero — this drives the
acquired axis (§8). An absolute-level anchor $\max(x_{m,s},0)$ keeps the abundance scale
in the score.

These combine into three **pressure-construction modes**. The names denote *which
downstream analysis each construction feeds* — not a transform applied edge-by-edge inside
a correlation:

$$
\begin{aligned}
\textbf{spine (z):}\quad          & \mathrm{expr}(m,s) = \mathrm{sm}_{m,s}\cdot z_{m,s}\cdot \max(x_{m,s},0) && \text{(share × dynamics × level — feeds aggregate coupling)}\\
\textbf{attribution (no-}z\textbf{):}\quad & \mathrm{expr}(m,s) = \mathrm{sm}_{m,s}\cdot \max(x_{m,s},0) && \text{(nothing zeroed for being stable — feeds shares)}\\
\textbf{acquired:}\quad           & \mathrm{expr}(m,s) = \mathrm{sm}_{m,s}\cdot \mathrm{dev}_{m,s}\cdot \max(x_{m,s},0) && \text{(centred on healthy, not tumour)}
\end{aligned}
$$

Using the wrong construction silently corrupts the answer: dynamics for the coupling-facing
spine, no-dynamics for attribution (who carries the steady-state pressure), healthy-centred
for acquisition. **Important:** these build the *aggregate* pressure (over a gene or
program). A single-edge or single-module coupling does **not** apply this transform — it
uses the arm's raw abundance directly; §6.3 makes the predictor-by-resolution explicit.

### 4.2 The master formula

Everything downstream is a piece of, or a summary of, this:

$$
c(m,g,s) = \mathrm{expr}(m,s)\cdot w(m,g),
\qquad
\mathrm{pressure}(g,s) = \sum_{m\in R(g)} c(m,g,s).
$$

The contribution $c(m,g,s)$ is the **edge-resolution pressure atom**; a gene's pressure
is its sum over the gene's regulators $R(g)$. The default aggregation is the sum (a gene
with more, stronger regulators genuinely carries more pressure); the mean is a
sensitivity variant.

An optional **availability gate** scales pressure by per-sample RISC capacity, mapped
through a bounded sigmoid into $[g_{\min},1]$:

$$
\mathrm{gate}(s) = g_{\min} + (1-g_{\min})\,\sigma\!\big(k\,(\kappa(s)-c_0)\big).
$$

Repression is not free: a miRNA acts only once it is **loaded into an Argonaute** and the
loaded complex **engages a TNRC6/GW182 effector** that recruits the degradation machinery,
and both steps draw on a finite pool every miRNA competes for — which is why a per-sample,
miRNA-independent capacity term is the right shape for a global gate. Capacity therefore has
two both-required components, each a standardized abundance: a **loadable-core** term
$\mathrm{ago}(s)$ — an abundance-weighted mean over AGO1–4 that up-weights AGO2, the
predominant and only slicer-competent Argonaute — and an **effector** term $\mathrm{eff}(s)$
over TNRC6A/B/C. Because either resource can be rate-limiting, capacity is their
**co-limiting minimum** (Liebig's law of the minimum, not a sum — abundant AGO must not
paper over absent effector):

$$
\kappa(s) = \min\!\big(\mathrm{ago}(s),\,\mathrm{eff}(s)\big).
$$

Here $\sigma(\cdot)$ is the logistic sigmoid, $g_{\min}\in(0,1)$ the floor (the gate never
falls below it, so the layer can only attenuate, never erase — even a low-capacity tumour
keeps basal repression), $c_0$ the capacity midpoint, and $k$ the slope. The AGO weighting
reflects abundance dominance, not catalysis (slicing matters mainly for the rare
extensively-complementary site). It is reported as a sensitivity layer, not a mechanism —
AGO/TNRC6 *mRNA* is an imperfect readout of functional RISC capacity, so the bounded form is
deliberately conservative, and empirically it rescales rather than reorders programs.

### 4.3 Acquired and lost pressure, at every resolution

Acquired pressure is what the tumour gained over the healthy state; lost pressure is what
it shed. It is defined identically at each resolution by differencing the same
contribution across states (the cross-state comparison uses the no-dynamics mode, since
$z$ is standardized *within* a state and is not comparable across states):

$$
\begin{aligned}
\text{edge:}\quad     & \mathrm{gain}(m,g) = c_{\text{tumour}}(m,g) - c_{\text{healthy}}(m,g),\\
\text{gene:}\quad     & \mathrm{acquired}(g) = \mathrm{pressure}_{\text{tumour}}(g) - \mathrm{pressure}_{\text{healthy}}(g),\\
\text{program:}\quad  & \Pi_{\text{acq}}(m,\mathcal H) = \mathrm{prior}(\mathcal H)\cdot\!\sum_g \big[-\max(c_{\text{tumour}}-c_{\text{healthy}},0)\cdot \mathrm{infl}(g)\big].
\end{aligned}
$$

(The program-level form previews the architecture layer of §5: $\mathrm{infl}(g)$ is the
gene's net downstream influence *within the program* — how much repressing $g$ changes the
program's output, defined in §5.1 — and $\mathrm{prior}(\mathcal H)$ is the program's
tumour-direction polarity, §5.2.)

Two healthy anchors are complementary: a true-healthy reference (large baseline jumps,
the canonical oncomiR program) and a normal-adjacent reference on the *same* measurement
platform (no cross-platform artifact, ranking the final malignant step). Section 8 gives
the trajectory machinery in full.

### 4.4 Share and specificity

**Share** decomposes a gene's pressure among its regulators. Writing the over-bar for the
mean over samples and a dot for the sample index it runs over — so $\overline{\lvert c(m,\cdot)\rvert}$
is the across-sample mean of an edge's absolute pressure, and $\sum_{m'}$ runs over the gene's
regulators —

$$
\text{abs-share}(m) = \frac{\overline{\lvert c(m,\cdot)\rvert}}{\sum_{m'}\overline{\lvert c(m',\cdot)\rvert}}
\quad\text{(ranking)},
\qquad
\text{signed-share}(m) = \frac{\overline{c(m,\cdot)}}{\sum_{m'}\overline{c(m',\cdot)}}
\quad\text{(coherence diagnostic).}
$$

When the signed share diverges from the magnitude share, regulators are *cancelling*.
Under the no-dynamics attribution mode every contribution is non-negative, so the two
coincide — the reason attribution is computed there.

**Specificity** is the fraction of an arm's *total* target-mass — summed over all the genes
$T(m)$ it regulates — that lands on one gene:

$$
\mathrm{spec}(m,g) = \frac{\overline{\lvert c(m,g,\cdot)\rvert}}{\sum_{g'\in T(m)}\overline{\lvert c(m,g',\cdot)\rvert}}.
$$

The reference point is the even-spread baseline $1/\lvert T(m)\rvert$: an arm with dozens
to hundreds of targets that spread its budget evenly would put well under one percent on
each. So a value of a *few* percent already marks a gene that commands a disproportionately
large, focused share of the arm's budget — several-fold above the even-spread floor — and
hence a cleaner single-gene readout (here $\gtrsim 0.05$). A value at or below that floor
($\lesssim 0.01$) marks a gene the arm regulates only as one of a promiscuous crowd. (It is
a *relative concentration*, not an absolute majority — even a focused arm rarely puts most
of its budget on a single gene.)

> **Seed-family non-identifiability — a caveat on *every* single-arm number.** Share,
> gene-side rank, and specificity all credit a *named arm*, but paralogues that share a
> seed (for example the miR-29a/b/c or miR-302/372/373/520 families) target the same sites
> and co-vary in expression, so a per-arm attribution cannot say *which* family member
> actually carries the pressure — the credited arm may be standing in for its family. This
> is not a peripheral worry: across the cohort, within-seed-family arm pairs co-vary about
> twice as strongly as random pairs and the largest multi-arm hubs are 20–37% redundant
> (one effective regulator masquerading as several). We therefore read share/specificity at
> the **seed-family** level as the safe unit and treat the specific within-family arm
> assignment as a nomination, not a claim; the realized-coupling inference makes the same
> grouping explicit (the seed family is the unit for both the dependence-aware FDR and the
> permutation null — §6, §9.3). Singleton-seed arms are unaffected.

---

## 5. Architecture: modeling a set of genes

A program is not a bag of independent genes. A miRNA that represses an upstream activator
changes the whole downstream program; one that represses a tumour-suppressor is pro-tumour
regardless of network position. The architecture layer overlays three orthogonal lenses
on otherwise-flat pressure.

### 5.1 Topology — network position

Build the program's regulatory wiring from **OmniPath** — a curated meta-database of
directed, signed (activating / inhibiting) molecular interactions — induced on the
program's member genes. Each ordered gene pair gets a **calibrated continuous sign**, so a
conflicting (dual-mode) pair contributes a graded value instead of being flattened:

$$
\mathrm{ns}(s\to t) = \frac{n_{+} - n_{-}}{n_{+} + n_{-}} \in [-1,+1].
$$

A finite-horizon signed-path propagation then gives each gene its net downstream
influence:

$$
T = \sum_{k=1}^{H}(\alpha S)^{k},
\qquad
\mathrm{infl}(g) = 1 + \sum_j T_{g,j}
\qquad (\alpha = 0.25,\ H = 6),
$$

where $S$ is the signed adjacency. A central activator scores positive, an inhibitor hub
negative. Because a miRNA *represses* its target, the effect on program output flips with
the node's influence:

$$
\Delta_{\text{prog}}(m) = \sum_g \big[-c_{\text{tumour}}(m,g)\cdot\mathrm{infl}(g)\big].
$$

This is the only lens that captures *indirect*, network-mediated effects — everything
else is local.

### 5.2 Program-polarity prior — pro- or anti-tumour

A coarse, curated polarity per program ($+1$ pro-tumour engines such as EMT, hypoxia,
glycolysis, E2F, G2M, MYC; $-1$ tumour-suppressive or anti-tumour-immune programs such as
apoptosis, p53, interferon, inflammatory, immune-rejection; $0$ ambiguous, kept
directional but unpolarized):

$$
\Pi(m,\mathcal H) = \mathrm{prior}(\mathcal H)\cdot\Delta_{\text{prog}}(m).
$$

A positive value pushes toward malignancy — either by activating a $+1$ engine or by
**releasing a brake** on a $-1$ program. A clean mechanistic result follows directly:
because miRNAs repress, they *damage* the proliferation engines (negative pro-tumour on
E2F, G2M, EMT) and turn malignant mainly by **brake-release on the suppressive programs**
(for example, the canonical oncomiR releasing the p53 and apoptosis brakes, and another
releasing immune-evasion). Every program is, on net, *damaged* by the miRNA layer; only
the targeted sub-pattern is oncogenic.

### 5.3 Gene-role and malignancy — driver-level direction

Topology has a blind spot for canonical oncomiRs, because they mostly damage the engines
they broadly hit, and because some key targets act through metabolites rather than gene
nodes. A curated **gene-role** layer (high-confidence cancer drivers, each labelled
oncogene, tumour-suppressor, or dual) provides a per-gene malignancy sign
$s_{\mathrm{mal}}(g)\in\{-1,0,+1\}$:

$$
\mathrm{Mal}(m) = \sum_g \big[-s_{\mathrm{mal}}(g)\cdot c_{\text{tumour}}(m,g)\big],
\qquad
\mathrm{Mal}_{\text{arch}}(m) = \sum_g \big[-s_{\mathrm{mal}}(g)\cdot c_{\text{tumour}}(m,g)\cdot w_{\text{arch}}(g)\big],
$$

where the non-negative centrality weight $w_{\text{arch}}$ scales magnitude without
re-flipping the sign. This recovers the textbook oncomiR roster that topology missed —
external validation that the layer captures real driver biology.

> Topology and gene-role are kept side by side, not merged: they answer different
> questions — "what the wiring says" versus "what the driver annotation says" — and the
> canonical oncomiR is exactly the case where both must be read. For broad, abundant arms,
> the gene-role, coupling, and acquired axes are the reliable readouts, not the topology
> sign.

A continuous, **all-gene** complement uses functional-genomic dependency
(genome-wide CRISPR essentiality in breast cell lines): genes the tumour depends on
(strongly essential) versus tumour-suppressor-like genes,

$$
\mathrm{Mal}_{\text{cont}}(m) = \sum_g w_{\text{dep}}(g)\cdot c_{\text{tumour}}(m,g),
$$

which extends coverage from a few hundred curated drivers to nearly all genes and shows
the miRNA layer landing hardest on proliferation dependencies. A transcription-factor
census adds a master-regulator identity flag on top of pure topology centrality.

### 5.4 Roll-ups

Each architecture score rolls up across the universe (per arm across programs, per
program across arms, by seed family) and runs in an **acquired-weighted** parallel that
isolates the tumour-acquired pro-tumour pressure from the constitutive component. The
universe-level unsupervised structure — pressure programs and their subtype divergence —
is described in §9.2.

---

## 6. Coupling: realized pressure and confounding

Pressure is potential. **Coupling** asks whether repression is actually observed: does a
tumour with more of the miRNA carry less of the target, after the obvious confounders are
removed?

### 6.1 The core idea

Realized repression should appear as a **negative correlation** between the miRNA *dose*
and the target's expression. The dose is the arm's **raw abundance** at single-regulator
resolution (one edge, or one arm's target module) and the **aggregate pressure** — the
weighted, expression-transformed sum over regulators (§4) — where many regulators converge
on one gene or program (the only way to fold them into a single predictor). A naive
correlation, however, is dominated by confounders — a more cellular, more aneuploid, more
proliferative tumour drags everything with it. We therefore use a **partial Spearman
correlation**: regress the confounders out of *both* variables, then correlate the residual
ranks. A negative partial correlation that survives is realized repression. Rank
correlation (not linear) is used because expression is heavy-tailed and the expected
relationship is monotone, not strictly linear.

### 6.2 The confounders

Each canonical confounder is included for a specific mechanistic reason.

| Confounder | What it is | Why it confounds | 
|---|---|---|
| **Tumour purity** | fraction of malignant cells | stroma and immune infiltrate carry their own miRNAs and mRNAs, so a low-purity tumour looks "high miRNA / low target" compositionally |
| **Genomic instability** | homologous-recombination deficiency / scar burden | instability co-drives both arms; the most unstable subtype over-expresses many loci |
| **Target copy number** | the target gene's own copy number | an amplified target's mRNA rises for reasons unrelated to the miRNA, masking or faking repression |
| **Proliferation** | a cell-cycle expression signature | the cell cycle drags both the proliferation-miRNA clusters and their targets and must be removed |
| **Composition (immune/stromal/epithelial)** | deconvolved cell-content scores | the most complete composition control; the epithelial axis specifically de-confounds lineage co-expression |
| **Sequencing batch** | technical run/plate identity (the small-RNA and mRNA assay plates; per-state batch for the healthy and normal-adjacent references; the protein-assay plex for the proteome cohort) | a shared technical batch co-varies features on *both* arms, so it can fabricate **or** mask coupling; it is removed as **covariate dummies**, not by batch-effect *subtraction* (which over-corrects an unbalanced design and can manufacture the very feature–feature correlations that are the readout) |

Three subtleties make the proliferation handling careful rather than reflexive:

- **Three proliferation proxies** of differing cell-cycle-gene dependence are used in
  turn, to test whether repression survives proliferation adjustment *without
  over-correcting through a shared pathway* — the "we just removed the signal" objection.
- **Proliferation is a *suppressor* confound, not the source of the signal.** Within the
  most proliferative subtype it associates positively with *both* pressure and target
  expression, so adjusting for it makes the negative coupling *more* negative — it was
  masking, not faking. The coupling therefore strengthens under proliferation
  adjustment, a robustness result rather than a fragility.
- **A covariate ladder, not a single number.** Coupling is reported across a nested
  sequence of covariate sets (none → purity + instability → + target copy number → +
  proliferation), **with sequencing-batch dummies appended to every adjusted rung** (the raw
  rung is the only batch-free one), so a reader sees exactly what each adjustment costs. The headline
  significance flag requires a negative partial correlation with a false-discovery rate
  below five percent; a parallel flag requires survival of the proliferation-adjusted
  ladder. A variability gate (low standard deviation *and* low interquartile range) marks
  flat arms whose null coupling should not be over-read.

**Multiplicity, dependence, and an empirical null — at every rung.** Two precautions guard
the significance calls, and both are applied at *each* resolution below, not only the
headline one. First, the false-discovery control is made **dependence-aware**: because
seed-family paralogues are non-independent (within-family arm pairs co-vary about twice as
strongly as random pairs), plain Benjamini–Hochberg over edges would be anti-conservative,
so we add a Benjamini–Yekutieli correction (valid under arbitrary dependence) and, where the
test unit is a miRNA, collapse each seed family to a single test (Simes) so a multi-paralogue
hub spends one slice of the error budget rather than several — the de-duplicated discovery
count is reported in *families*, not edges. Second, every reported coupling is re-tested
against an **empirical permutation null** (Freedman–Lane residual permutation, one shared
sample permutation per draw so the family correlation structure is preserved), so the
significance does not rest on the asymptotic partial-correlation p alone; the same joint null
yields a count statistic for the universe-level "how many programs couple" claim. These run
at edge, gene, target-set, program, and universe resolution alike.

### 6.3 Coupling at every resolution

The same partial-correlation machinery is applied at each rung of the hierarchy, against
a different predictor and target.

| Resolution | Predictor | Target | The coupling object |
|---|---|---|---|
| **Edge** | arm abundance | target expression | per-edge partial correlation, per tissue state; a three-letter trajectory pattern across states |
| **Gene** | the gene's aggregate incoming pressure | the gene's own expression | whether the *gene* is net-repressed |
| **Target-set** | arm abundance | the mean standardized expression over the arm's targets (the *module composite*) | whether the arm's whole module moves down together |
| **Family** | a seed-family aggregate predictor | target | family-as-unit coupling |
| **Program** | per-program pressure | per-program mean target expression | the count of programs that couple |
| **Universe** | all programs | — | the fraction of programs and edges that couple, and cross-cohort concordance |

The **module composite** is the strongest confounder-independent signal: standardize each
target within a state, average across the arm's targets, then take one partial correlation
of the arm's abundance against that composite. A strongly negative value means the
*coordinated* module is repressed — something no single-edge correlation can see. A
role-stratified version (tumour-suppressor versus oncogene targets) is an interpretive
overlay, deliberately not folded into one number, since whether repression is pro- or
anti-tumour is the architecture layer's job (§5).

**Why the predictor switches with resolution — and what actually differs.** The edge and
aggregate predictors are *not* a reparametrization of one another, but most of the
difference is inert. Of pressure's pieces, the static evidence weight is a sample-constant
scalar (invisible to a rank test) and the within-arm $z$ is an affine transform of
abundance (rank-identical for a single arm), so neither changes a single-edge coupling; the
**softmax share** — which pulls in the co-regulators — is the one piece that does. A
head-to-head check (the same confounder-adjusted partial Spearman applied under every
predictor) confirms it:

- **Edge level:** the within-arm $z$ reproduces raw abundance (sign-concordance 0.96,
  recovering 94% of its FDR-significant edges), while the share *alone* recovers only ~49%
  and flips its median positive. So **raw abundance is the most sensitive single-edge
  predictor** — which is exactly why edge and module coupling use it.
- **Gene level:** the share alone is degenerate (per-gene shares sum to about 1 in every sample,
  so the predictor is near-constant), but the *full* aggregate pressure (share × dynamics ×
  level) recovers ~89% of a naive abundance-sum's hits **and finds ~2–3% more, with a deeper
  median** — so the pressure construction is a genuinely better way to fold many regulators
  into one dose than summing raw abundances.

Net: pressure differs from abundance **only** through the softmax share, and that term
*weakens* a single edge but *strengthens* the gene aggregate — which is why the predictor
is abundance at edge/module resolution and aggregate pressure at gene/program resolution.

### 6.4 The potential × realized join (pressure × coupling)

The purpose of measuring both potential and realization is to **join** them. In this join
*pressure* is the **potential** axis (at edge level its trajectory tracks the miRNA's
abundance, since the edge weight is static across states) and *coupling* is the
**realized** partial correlation — so the join is potential × realized, **not** a
correlation of one against the other. Each edge carries a pressure (potential) trajectory —
does the miRNA gain or lose abundance across states — and a coupling (realized) trajectory
— the per-state significance pattern. Their join is the headline class:

- **acquired-and-realized** — pressure rises *and* coupling turns repressive (a newly
  switched-on repressor);
- **field-established** — already repressive in normal-adjacent tissue, retained in
  tumour;
- **acquired-unrealized** — pressure rises but coupling never realizes;
- **released brake** — repressive in healthy tissue, gone in tumour;
- **constitutive** — repressive in both healthy and tumour.

This join lets us state, as a *realization* claim immune to the pressure-side caveats,
that healthy breast actively holds certain oncogenic and invasion targets under miRNA
repression and the tumour releases the brake. Whether such a loss reflects a genuine loss
of repression is then tested by conditioning on copy number, promoter methylation, and
chromatin accessibility one at a time (§9).

---

## 7. Protein realization: beyond the messenger-RNA layer

miRNA repression ultimately acts on the **protein** product. The messenger-RNA-level
coupling sees destabilization but cannot see a pure translational block, so realization
onto protein is tested directly, three ways, using mass-spectrometry proteomes. For a
target with standardized mRNA and protein levels,

$$
\begin{aligned}
\textbf{direct protein:}\quad & \rho\big(\mathrm{pressure},\ \mathrm{protein}\big) && \text{expected negative}\\
\textbf{messenger-to-protein gap:}\quad & \rho\big(\mathrm{pressure},\ \mathrm{mRNA}-\mathrm{protein}\big) && \text{expected positive}\\
\textbf{protein beyond mRNA:}\quad & \rho\big(\mathrm{pressure},\ \mathrm{protein}\mid \mathrm{mRNA}\big) && \text{expected negative}
\end{aligned}
$$

The third quantity — protein residualized on its *own* messenger RNA — is the key idea:
it removes everything the **measured** mRNA already explains, so a surviving negative
correlation is a protein-level effect that the bulk steady-state mRNA does **not** capture.
This is the distinction between reading protein *together with* the mRNA coupling (direct)
and *after subtracting* it (residual).

It is tempting to call that residual "translational repression," and it is consistent
with a translational block — but it is **not** exclusively that, so we deliberately label
it *beyond-(measured)-mRNA* rather than *translational*. At least four non-translational
routes produce the same residual: (i) **destabilization that bulk RNA-seq misses** —
miRNA-driven deadenylation / poly(A) shortening lowers a message's output before its total
count drops, and poly(A)-selected libraries register this differently from total-RNA
libraries; (ii) **3'-UTR isoform / alternative-polyadenylation effects**, where gene-level
mRNA quantification averages over isoforms that differ in whether the miRNA site is even
present; (iii) **protein-versus-mRNA turnover kinetics** — protein level is a time-integral
of past mRNA, so a transient or recently-onset mRNA dip can surface in protein while being
smoothed out of a single steady-state mRNA snapshot; and (iv) **indirect routes**, where
the miRNA represses an intermediate that changes the target's protein without changing the
target's own mRNA. Platform noise between sequencing and mass spectrometry adds further
slack. The residual is therefore best read as *regulation the mRNA snapshot cannot see*,
of which translational repression is one — biologically expected — component.

The headline findings: the spine validates at the protein layer in an *independent*
patient cohort (dozens of genes show significant negative pressure-to-protein coupling,
led by epithelial-keratin, hypoxia, and matrix targets — KRT8, DDAH1, BMP1 — *not* a
canonical EMT transcription factor, ZEB1/2 being absent from the protein-layer drivers);
the repression is **predominantly mRNA-mediated**
(residualizing on mRNA collapses most of it, leaving a small beyond-mRNA residual —
consistent with canonical mammalian miRNA biology, in which destabilization dominates and
a smaller component acts beyond the steady-state message); and the surviving residuals
are **known, subtype-structured edges**, including a protein-level recovery of the central
proliferation-cluster-to-tumour-suppressor hub. A decoy control confirms the residual
couplings are **target-specific** — each driver arm anti-correlates with its cognate
protein and not with the proteins it does not target — and they survive batch adjustment.

---

## 8. The healthy → normal-adjacent → tumour trajectory

How a miRNA's pressure *moves* across tissue states is itself a finding, so the acquired
concept gets dedicated machinery across three states: true-healthy tissue, normal-adjacent
tissue (which carries a pre-malignant field effect), and tumour.

### 8.1 The healthy anchor

The healthy reference is **GTEx breast** (346 normal-tissue samples), quantile-normalized
onto the tumour-cohort scale and joined by stable molecular accession. Because the GTEx
pipeline collapses a subset of canonical arms to zero, two imputation layers repair
coverage — a validated rank-transfer from an **independent healthy miRNA atlas
(DIANA-miTED)** (primary) and a seed-family fallback —
while arms with no reference anywhere are left explicitly un-imputed and flagged rather
than fabricated, since imputing a truly absent arm would erase a real acquired gain. The
healthy abundance matrix is repaired *at source* so the acquired axis is free of the
collapse artifact, while the realized-coupling layer is left on the raw matrix and so is
unaffected. The imputed baseline feeds only the acquired axis, never the coupling — so all
headline coupling findings are baseline-independent.

### 8.2 Two acquired axes

A single rank delta is insufficient, because rank *saturates*: an already-abundant arm
cannot gain rank even as its absolute abundance rises several-fold. Two complementary
axes are therefore used,

$$
\underbrace{d_{HT}(m) = r_{\text{tumour}}(m) - r_{\text{healthy}}(m)}_{\text{rank delta (saturating, cross-platform-robust)}}
\qquad
\underbrace{\mathrm{log2fc}(m) = \tilde x^{\,\text{tumour}}_m - b^{\,\text{healthy}}_m}_{\text{magnitude delta (recovers the surge)}},
$$

where $r$ is a within-state percentile rank and $b$ the healthy baseline on the
common scale. The magnitude axis recovers exactly the strongly up-regulated oncomiRs the
rank axis is blind to. The default acquired call is the union: an arm is an acquired
gainer if it gains on *either* axis.

### 8.3 Trajectory classification

The three legs are healthy→normal-adjacent (the field effect), normal-adjacent→tumour
(the tumour-specific increment), and healthy→tumour (the acquired headline). On these,
each miRNA is classed as stable, progressive gain or loss, tumour-acquired, field-effect,
or non-monotonic; each edge inherits a class from the join of the miRNA's movement with
its per-edge coupling trajectory (§6.4); and a subtype-specificity index marks whether an
acquired arm is pan-subtype or concentrated in one subtype. Because the acquired axis is
anchored to the true-healthy state — never the field-affected normal-adjacent tissue —
"acquired over healthy" means exactly that.

---

## 9. Validation

Findings are graded by the strength of evidence behind them: **statistically
supported** (a significant, confounder-adjusted effect), **descriptive** (a coherent
pattern), or **a nomination** (a hypothesis awaiting a perturbation assay). The three
validation questions are: did the model recover known biology, did it deepen it, and did
it find anything new?

Underpinning all three is an **inference-rigor layer** that guards the usual soft spots.
Every coupling claim is re-tested at *each* resolution against an empirical permutation
null (so significance does not rest on an asymptotic p alone) and a **dependence-aware**
false-discovery control — Benjamini–Yekutieli plus, where the unit is a miRNA, a seed-family
collapse, so co-varying paralogues cannot each spend an independent slice of the error
budget. The pressure construction itself is validated **out-of-sample**: in patient
cross-validation its candidate ranking reproduces on held-out patients, so it is not overfit
to the cohort it is tested against, and the weighting constants are stability-tested rather
than asserted.

### 9.1 Recovered — the positive controls

The model reproduces textbook biology *unprompted* — recovered from inputs that never saw
the patient expression data — which is evidence that the recovered structure is not an
artifact of the construction. The canonical tumour-suppressor-miRNA roster falls out of the
architecture anti-tumour leaders, and the canonical oncomiR roster out of the gene-role score. The
established miRNA-to-extracellular-matrix, proliferation-cluster-to-tumour-suppressor,
and let-7-to-mitotic axes all re-emerge as unsupervised sub-modules — internal
validation, since their *absence* would make the decomposition suspect. High-evidence
edges couple more strongly than sequence-only predictions in two independent cohorts, and
the miRNA-to-target anti-correlation **replicates in an independent patient cohort** (a
clear positive concordance, with two-thirds of significant repressive edges keeping a
negative sign) — an independent-patient leg the same-patient proteome cannot provide.

### 9.2 Deepened — the same players, newly resolved

The genuine methodological contribution is **state-resolved** and **within-program**
structure that standard enrichment and single-edge analyses cannot see:

- **State-resolved coupling** — the literature asserts a miRNA-target pair and the
  miRNA's down-regulation, but rarely measures the brake being *engaged in normal tissue
  and released in tumour*.
- **Within-program sub-modules** — a program resolves into a shared, pan-subtype backbone
  plus named subtype-specialist sub-modules (for example a basal immune-restriction
  module, and tumour-suppressor repression localized to a *minority* sub-module), where
  ordinary enrichment gives a single score per program.
- **Protein-layer confirmation** of the central hub in an independent proteome.
- **An end-to-end copy-number-to-network attribution** showing that the copy-number-driven
  component of miRNA dosage reproduces a large fraction of the repressive network and
  survives the full confounder stack.

### 9.3 Novel and candidate-novel — the discovery surface

The cleanest discovery surface is the **orphan edges**: sequence- and crosslink-predicted
interactions absent from functional curation, kept out of the spine precisely so they can
serve as an unbiased screen. Of roughly twenty thousand orphan edges screened against the
protein layer, a few thousand are testable, several hundred show protein anti-correlation,
a couple hundred survive residualization onto mRNA, and the strongest tier — translational
*and* replicated in a second cohort — numbers in the dozens, most of them fully
uncurated.

Three convergent strengthenings make these more than a prediction list:

- **Method self-validation** — the established miRNA-to-collagen axis emerges *de novo* at
  uncurated target genes, providing a credibility anchor.
- **Triple-cohort support** — about a hundred genuine orphan edges keep a negative sign
  across three independent datasets and two molecular layers (messenger RNA, protein, and
  an independent cohort), far above the chance rate.
- **A composition-robust core** — re-testing the prediction-only tier after adding
  subtype, immune, stromal, endothelial, and target-copy-number covariates leaves a clean
  core of about a dozen-and-a-half de-novo nominations (led by a metabolic convergence)
  while the composition-loaded artifacts drop out.

> These are **wet-lab nominations** — a prioritized queue for crosslinking and reporter
> assays — not validated edges. Protein anti-correlation is not direct binding,
> co-expression confounding persists (especially for matrix genes that co-vary with
> stroma), and seed-family arms cannot be separated — members of one seed family (for
> example the miR-29 paralogues a, b and c) share the same seed sequence, hence the same
> predicted target sites, and co-vary in expression, so an anti-correlation cannot say
> *which* member drives the repression of a shared target; a nominated arm may simply be
> standing in for its family. The de-novo recovery of the collagen axis is the anchor; the
> fully uncurated translational hits and the composition-robust core are the discovery
> surface.

### 9.4 The validation ladder

| Rung | Evidence | Strength |
|---|---|---|
| Reproduces known biology | recovers textbook edges unprompted | established |
| Robustness & inference | permutation null + dependence-aware (BY / seed-family) FDR at every rung; construction validated out-of-sample; constants stability-tested | statistical |
| Messenger-RNA coupling | significant, confounder-adjusted partial correlation | statistical |
| State-resolved trajectory | brake-release and acquired-realized classes | descriptive / hypothesis |
| Protein layer | independent proteome residual | statistical |
| Independent cohort | replication in new patients | statistical / established |
| Wet-lab queue | nominations for perturbation assays | hypothesis |

---

## 10. Design principles and honest limits

The stance, as it cashes out: static evidence times dynamic abundance keeps "what is
possible" separate from "what this tumour expresses"; a non-circular prior cannot leak into
the findings it is tested against — and this is now *checked, not just asserted*, since the
pressure construction reproduces its ranking on held-out patients (so it is not overfit to
the cohort it is validated on); two measurements (potential and realization) are joined,
never conflated; everything is reported adjusted *and* raw, gated *and* ungated; and every
claim is labelled by its rung of evidence.

The standing limits, stated plainly:

- **Gene-gene topology cannot see metabolite-mediated regulation** (a target acting
  through a lipid rather than a gene node) — a true ceiling, partly worked around by the
  gene-role layer.
- **The mechanism of a lost coupling is not fully pinned.** Copy number and promoter
  methylation are ruled out by conditioning, chromatin accessibility is largely null with
  a small minority at growth and cell-cycle loci, and the leading remaining candidate is
  3'-untranslated-region shortening and competitive (decoy-RNA / RNA-binding-protein)
  regulation.
- **Arms with no healthy reference** are a coverage gap, not acquisition — to be resolved
  by extending healthy-tissue coverage, never by imputing absent arms.
- **The program-polarity prior is coarse** (one polarity per program). The tissue-specific
  edge re-ranking (§3.5) is kept as an overlay because proving it actually *improves* on
  the pan-context weighting would require an independent breast-specific miRNA dataset to
  test the reassigned regulators against — and no suitable one exists (one major breast
  cohort lacks miRNA data; another's miRNA is re-derived from the same source rather than
  being independent). The refinement is built and inspectable, but its utility is not yet
  externally validated.
- **The pressure construction is one of several near-equivalent choices** — validated as
  non-overfit out-of-sample, but not uniquely optimal: on realized coupling the bare
  un-normalized weighting is marginally ahead, and the promiscuity denominator is
  coupling-neutral, earning its place on attribution (which arm is a gene's dominant
  regulator) rather than on coupling.
- **A bulk-tissue ceiling** — every quantity is an association, not a binding event;
  cell-composition is controlled but never fully eliminated.

---

## Glossary

| Term | Definition |
|---|---|
| **Pressure** | potential repression = expression × interaction weight, summed over a gene's regulators |
| **Contribution** $c(m,g,s)$ | the edge-resolution pressure atom; every other quantity is built from it |
| **Interaction weight** $w(m,g)$ | static evidence strength divided by a promiscuity denominator |
| **Acquired / lost pressure** | tumour pressure minus a healthy anchor, at edge, gene, or program level |
| **Share** | a regulator's fraction of a gene's total pressure (arm-side and gene-side) |
| **Specificity** | how concentrated an arm's budget is on one gene (focused target versus promiscuous hub) |
| **Coupling** | observed, confounder-adjusted anti-correlation = *realized* pressure (predictor = arm abundance at edge/module resolution, aggregate pressure at gene/program) |
| **Module composite** | partial correlation of an arm's abundance against the mean standardized expression of its whole target set |
| **Architecture** | topology influence + program-polarity prior + gene-role malignancy, applied to a program |
| **Released brake** | an interaction repressive in healthy tissue and lost in tumour |
| **Acquired gainer** | an arm that gains pressure over healthy, on a rank or a magnitude axis |
| **Orphan edge** | a sequence- or crosslink-predicted interaction absent from functional curation — the discovery screen |
| **Availability gate** | a per-sample bounded multiplier from RISC-machinery abundance — a sensitivity layer, not a mechanism |
