# Modeling miRNA repressive pressure on the Hallmark programs of cancer

**A conceptual *and* implemented framework.**
This document is the connective layer between the two reference docs: it explains
*what we are modeling and why each design choice is made*, then ties every concept
to its formula and its code. The two companions are:

- **`FORMULAS.md`** — the exact, code-faithful formula for every quantity (the *how*).
- **`LANDSCAPE_REPORT.md`** — the results: what the model found in TCGA-BRCA (the *what*).

Where this document gives a formula, the authoritative version is `FORMULAS.md`;
where it gives a number, the authoritative version is `LANDSCAPE_REPORT.md`. The
older `METHODS.md` / `BIOLOGICAL_SYNTHESIS.md` are consistent on mechanics but may
trail on the latest counts — defer to the two above.

> **One-sentence thesis.** A miRNA exerts *repressive pressure* on a gene; that
> pressure is a **potential** (how much repression the wiring + abundance *could*
> deliver) which may or may not be **realized** (an observed anti-correlation after
> confounders are removed); we model both, at five nested resolutions, across the
> healthy→tumor trajectory, and validate at the mRNA layer, the protein layer, and
> in independent cohorts.

---

## ⭐ THE CANONICAL MODEL STATEMENT — what *the* model IS

> **Scope note.** §§1–10 below describe the **pressure/coupling framework** — the resolutions,
> the edge universe, the confounders, the trajectory, the validation ladder — which stands.
> This block states **the estimator that is canonical for coupling, attribution,
> identifiability, and discovery**. Where an older section reads as though the heuristic
> pressure construction *is* the model, this block wins.
> *(Folded in from `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b — the canonical one-liner —
> which remains the capstone account of the whole architecture.)*

**One dense learned-τ² non-negative Bayesian posterior per gene, two readouts of it:** run
**π ≡ 1 (dense)** for **coupling + attribution + identifiability** (prior-inert by design), and
**evidence-π (inclusion)** for **discovery/selection**. Slab-scale carries μ·τ only. The
**seed-rarity** channel is **off by default** (a global rarity gain promotes obscure/unexpressed
seeds, not coupling specialists). Nothing else in the estimator space moved the needle:
dense-vs-inclusion is the sole axis, slab shape and CV-vs-EB are second-order. **This is the model
to cite as *the* learned model; the adaptive lasso is retired to a baseline.**

Per gene `g`, on confounder-residualised, z-scored **seed-family** predictors:

    r = X·β + ε ,  ε ~ N(0, σ²)                       # r = −resid(Y|C): repression → positive
    β_m = z_m · θ_m                                    # coefficient = inclusion × magnitude
    θ_m ~ N⁺(0, ν²)                                    # half-normal slab — β ≥ 0 ("miRNAs repress")
    z_m ~ Bernoulli(π_m)                               # inclusion indicator
    ν², σ² ~ InvGamma                                  # ν² (=τ²) LEARNED, not guessed

One knob, `π`, switches the two modes:
- **π ≡ 1 (DENSE)** → learned-τ² ridge. Every candidate in; ν² sets shrinkage. **The core-jobs fit.**
- **evidence/uniform/learned-base-rate π (INCLUSION)** → spike-and-slab. Genuine selection (PIP).
  **Discovery.**

**Why it wins (genome-wide):** the dense learned-τ² ridge beats the shipped adaptive lasso on OOF
coupling — mean ρ **−0.168 vs −0.152**, wins 58%, **Wilcoxon p = 9×10⁻¹⁶** — because at n≫p the
lasso's sparse selection over-shrinks the aggregate, and dense learned-τ² keeps the diffuse
many-arm signal. CV-ridge ≈ lasso, so **it is the *learned* τ², not L2 per se** — that, not
shrinkage, is the edge over the lasso.

**Two facts that constrain how you read any result from it:**
- **Engine = Gibbs** (`learned/spike_slab.py`: `_gibbs_posterior`, `_gibbs_ss`), not HMC.
- **The unit is the seed FAMILY.** Co-seed matures target identically, so the family→gene weight
  is the identified estimand; **an arm-level result is a nomination, not the estimand** (§F honesty
  rule — a member is nominable only where its abundance diverges from its family-mates).

Where the detail lives: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` (the whole architecture, four jobs
from one posterior, the discovery pipeline, the decoy control §6d) · `LEARNED_MODEL_METHODS{,_FORMAL}.md`
(formulas) · `LEARNED_MODEL_RATIONALE.md` (*why* per §) · `LEARNED_MODEL_VALIDATION.md` (numbers) ·
`docs/STATE_OF_PLAY.md` Axis 1 (current verdict).

---

## Contents

- [⭐ The canonical model statement — what *the* model IS](#-the-canonical-model-statement--what-the-model-is)
0. [The modeling problem and the stance](#0-the-modeling-problem-and-the-stance)
1. [The hierarchy of resolution](#1-the-hierarchy-of-resolution)
2. [The edge universe](#2-the-edge-universe)
3. [Expression modeling — the dynamic piece](#3-expression-modeling--the-dynamic-piece)
4. [Pressure — the master formula and its variants](#4-pressure--the-master-formula-and-its-variants)
5. [Architecture — modeling a set of genes](#5-architecture--modeling-a-set-of-genes)
6. [Coupling — realized pressure and confounding](#6-coupling--realized-pressure-and-confounding)
7. [Protein realization — beyond the mRNA layer](#7-protein-realization--beyond-the-mrna-layer)
8. [The healthy → NAT → tumor trajectory (acquired pressure)](#8-the-healthy--nat--tumor-trajectory-acquired-pressure)
9. [Validation — recovered, deepened, novel](#9-validation--recovered-deepened-novel)
10. [Cross-cutting design principles and honest limits](#10-cross-cutting-design-principles-and-honest-limits)
- [Appendix: terminology glossary](#appendix-terminology-glossary)

---

## 0. The modeling problem and the stance

### 0.1 What "miRNA pressure" is

A microRNA (miRNA) is a ~22-nt RNA that, loaded into the RISC/AGO complex, base-pairs
a target mRNA's 3′UTR and represses it (mRNA destabilization predominantly, some
translational block). One miRNA targets **many** genes; one gene is targeted by
**many** miRNAs. The biological object we want to quantify is therefore not a single
binding event but a **field of repression**: how much downward pressure the entire
miRNA layer places on a gene, a program, a tumor.

We call this **pressure**. It is deliberately a *potential* quantity — built from
**which edges exist** (curated/predicted interactions), **how strong each edge is**
(evidence), and **how abundant the miRNA is** (expression). It does *not* by itself
prove that repression happened in a given tumor. That second, harder question —
*was the pressure realized?* — is answered separately by **coupling** (an observed,
confounder-adjusted anti-correlation between the miRNA and its target).

**Keeping potential and realization as two separate measurements, then joining them,
is the central design decision of the whole framework.** Everything else hangs off it.

### 0.2 Why Hallmarks, why breast cancer, and why it generalizes

- **The 50 MSigDB Hallmark gene sets** are a curated, non-redundant atlas of the
  coherent biological *programs* of cancer (proliferation, EMT, hypoxia, apoptosis,
  p53, interferon response, …). Using them as the gene universe means every result
  is automatically interpretable as "pressure on *a named program*," and lets us roll
  pressure up and down a biological hierarchy (§1).
- **Breast cancer (TCGA-BRCA)** is the implementation substrate because it has the
  deepest matched multi-omics (miRNA-seq, RNA-seq, CNV, methylation, ATAC, and two
  independent CPTAC proteomes), plus well-defined molecular subtypes (PAM50) and
  independent replication cohorts (Buffa 2011).
- **The framework is substrate-agnostic.** Nothing in the pressure/coupling/
  architecture machinery is breast-specific; only the *prior* (the tumor-direction
  polarity, §5.2) and some confounders (HRD) are cancer-flavoured. The same model
  runs on any cohort with miRNA + mRNA and any gene-set collection.

### 0.3 The stance: a few non-negotiable principles

These recur everywhere; they are stated once here and assumed thereafter.

1. **Static edge weight × dynamic abundance.** The edge weight `w(m,g)` is fixed
   (it encodes curated/sequence evidence, which does not change between samples or
   states); all sample-to-sample and state-to-state variation enters through the
   miRNA's **expression**. This cleanly separates "what the literature says is
   possible" from "what this tumor is actually expressing."
2. **The prior is non-circular.** Edge weights, the HE filter, breast re-ranking,
   gene roles — none of these ever use TCGA expression. So they cannot leak into the
   coupling findings they are tested against. (This is why, e.g., the breast-context
   re-ranking in §2.6 is an *overlay*, never folded into the pressure budget.)
3. **Always report the adjusted *and* the raw, the gated *and* the ungated.** Every
   coupling has a raw Spearman beside its confounder-adjusted partial Spearman; every
   gated pressure has its ungated twin. Differences are read as *sensitivity*, not as
   the headline.
4. **Potential ≠ realization ≠ causation.** Pressure is potential; coupling is an
   observed association (still not binding); a luciferase/CLIP assay would be
   causation. The framework is honest about which rung a claim sits on, and tags every
   discovery `S / P / H / R` accordingly (§9).

---

## 1. The hierarchy of resolution

The single most important structural idea is that **the same pressure object is
modeled at five nested resolutions**, and each resolution answers a different
question with a different aggregation. Zooming in localizes a signal; zooming out
contextualizes it.

```
   ┌─────────────────────────────────────────────────────────────────┐
   │  UNIVERSE          all 50 Hallmark sets (4,384 genes; 1,410       │
   │                    under miRNA pressure), 721 arms, 5,108 edges   │
   │    └── HALLMARK     one program (e.g. EMT, P53_PATHWAY)           │
   │          └── TARGET-SET   the genes one miRNA arm / seed family   │
   │                            represses (T(m))                       │
   │                └── GENE     one target gene g and its full        │
   │                              crowd of regulators R(g)             │
   │                      └── EDGE   one (miRNA m → gene g) interaction │
   └─────────────────────────────────────────────────────────────────┘
```

| Resolution | Unit | The question it answers | Pressure metric | Coupling metric |
|---|---|---|---|---|
| **Edge** `(m → g)` | one interaction | Does *this* miRNA repress *this* gene? | `c(m,g,s)` (§4.1) | per-edge partial-ρ (§6.3) |
| **Gene** `g` | one gene + its regulators `R(g)` | How much total pressure lands on the gene; who carries it? | `pressure(g,s)=Σ_m c` (§4.1); gene-side `share` (§4.4) | `gene_corepression` (§6.3) |
| **Target-set** `T(m)` | one arm's / family's targets | Does the arm coordinately repress its whole module? | per-arm mass / share (§4.3) | module composite ρ (§6.3) |
| **Hallmark** | one program's genes | Net effect of the miRNA layer on a program's *output* | architecture `program_output_change` (§5) | per-Hallmark coupling count (§6.3) |
| **Universe** | all programs | Cross-program structure (subtype programs, NMF) | NMF factors, cross-set roll-ups (§5.6, LANDSCAPE §6) | fraction of programs coupled (§6.3) |

Two facts make this hierarchy non-trivial (not just "sum the small things"):

- **Scope recomputation.** Every mean / SD / median / softmax that feeds pressure is
  recomputed *within the samples in the current scope* (whole cohort, or one PAM50
  subtype, or one state). A subtype table is therefore **not** the cohort numbers
  re-sliced — the relative abundances are re-derived inside the subtype.
- **Different aggregations carry different meaning.** A gene's pressure is the **sum**
  over its regulators; a program's effect is a **topology-weighted, sign-aware
  propagation** over its genes (§5), not a sum. The roll-up rule is part of the model.

---

## 2. The edge universe

The edges are the skeleton: the set of `(miRNA → gene)` interactions over which
pressure is built. The universe is **primarily experimental, with a controlled
computational extension.**

### 2.1 Experimental core + computational extension

| Layer | Source | What it contributes | Role |
|---|---|---|---|
| **Experimental core** | **miRTarBase** functional MTI (reporter assay, Western/protein, qPCR, perturbation, CLIP) | curated, literature-validated interactions | **the spine** — only these enter default (M0) pressure |
| Sequence prediction | **TargetScan 7.1** (weighted context++ score) | site strength from seed match + conservation; study-count-independent | extension / orphan screens / hybrid modes (M7/M11) |
| Crosslink depth | **ENCORI** AGO-CLIP / degradome | physical-occupancy depth boost on existing pairs | depth boost (α) on ~25% of pairs |

The default landscape spine is **miRTarBase-only (M0)**: 5,108 edges over 721
expressible arms and 1,410 target genes (cross-state landscape figures;
`LANDSCAPE_REPORT.md` §0). Sequence-only ("orphan") pairs do **not** enter the
default budget — they are screened separately (§9.3) so that prediction noise can
never inflate the headline. The hybrid modes (M7 = HE miRTar + capped TargetScan
fill-in; M11 = pure TargetScan) exist as sensitivity layers, not the default.

### 2.2 The "high-evidence" (HE) gate

An edge is **high-evidence** if it has functional support from a *direct* assay:

```
HE(m,g)  ⟺  n_functional_mti_studies ≥ 1
            AND ( n_reporter__functional_mti ≥ 1  OR  n_protein__functional_mti ≥ 1 )
            (fallback: evidence_score ≥ 3 when cross-count columns unavailable)
```

The HE gate is a *binary membership filter* (which edges are credible enough to model
at all); the edge **weight** below is the *continuous strength* of those that pass.

### 2.3 Edge weighting — component by component, with rationale

The edge weight has two factors: an **evidence numerator** (how well-supported the
interaction is) divided by a **promiscuity denominator** (how thinly the miRNA
spreads itself). Both are deliberate.

#### 2.3a Evidence numerator — `confidence_logclass`

Crucially, the weight does **not** use miRTarBase's raw study-count field. That field
inflates with publication count and treats a CLIP proximity hit like a luciferase
proof. Instead the spine recomputes a **quality-weighted, log-damped** score from the
per-assay-class Functional-MTI counts:

```
confidence_logclass(m,g)
   = 3.0 · log1p(n_reporter__functional_mti)     ← luciferase: direct proof of repression
   + 2.5 · log1p(n_protein__functional_mti )     ← Western/protein: direct
   + 1.5 · log1p(n_rna__functional_mti     )     ← qPCR/mRNA level
   + 1.5 · log1p(n_perturbation__functional_mti)
   + 0.5 · log1p(n_binding__functional_mti )     ← CLIP only: proximity ≠ repression
   + 0.3 · [same sum over the curator-flagged "weak" columns]
```

Three rationales are baked in:

- **Assay-directness weighting** (`3.0 … 0.5`): a reporter assay that *demonstrates*
  repression earns 6× a binding-only CLIP hit that only shows *proximity*. This is the
  single most important correction — it stops physical-occupancy data from masquerading
  as functional proof.
- **`log1p` within each class**: diminishing returns. The 10th luciferase paper adds
  little over the 3rd; `log1p` defuses publication-count inflation so a heavily-studied
  textbook edge doesn't crush a solid but less-famous one.
- **`δ = 0.3` weak discount**: curator-flagged weak evidence still counts, but at 30%.

#### 2.3b ENCORI depth boost

For the ~25% of pairs with an ENCORI entry, physical AGO-CLIP depth is added at α=0.5:

```
enc_depth(m,g) = 2·log1p(clipExpNum) + 1·log1p(degraExpNum)
               + 0.5·I(TargetScan ∧ miRanda ∧ PITA) + 0.25·log1p(pancancerNum)

evidence_score_eff(m,g) = confidence_logclass(m,g) + 0.5 · enc_depth(m,g)   (0 if no ENCORI row)
```

The base edge strength is then `w0(m,g) = log(1 + evidence_score_eff)`.

#### 2.3c Promiscuity denominator — the "one miRNA, many genes" problem

A hub miRNA that targets *everything* should not be allowed to dominate any one gene
merely by being promiscuous. So `w0` is divided by a per-arm **mass denominator**
(default `evidence_mass`):

```
D(m) = log(1 + Σ_{g'} log(1 + evidence(m, g')))        ← total curated target mass of arm m
edge_w(m,g) = w0(m,g) / max(D(m), 1e-9)
```

| Mode | `D(m)` | Reading |
|---|---|---|
| `none` | 1 | raw evidence (no promiscuity control) |
| `degree` | `log(1 + #targets(m))` | penalize *number* of targets |
| **`evidence_mass`** (default) | `log(1 + Σ log(1+ev))` | penalize *total curated mass* |
| `combined_mass` | adds the TargetScan-mass term | miRTar + TS (hybrid modes) |

**Rationale (measured, `held_out_tuning` + `denominator_attribution_sweep`):** `D(m)` is
a single per-arm scalar, so it acts on attribution, **not** coupling — and we verified
this rather than assert a coupling benefit. (i) **Coupling-neutral:** in 5-fold patient
CV the denominator barely moves gene-level coupling — bare `none` is in fact marginally
*ahead* out-of-sample and the spine sits within ~0.001 of it; the candidate ranking
reproduces perfectly across folds (train→test rank-ρ 1.0), so the choice is **not overfit**
(non-circular) but is also **not coupling-maximal**. (ii) **Specificity-invariant:** `spec`
is a within-arm ratio in which `D(m)` cancels exactly (rank-ρ 1.0). (iii) **Where it acts —
cross-arm attribution:** it re-assigns a gene's dominant regulator for ~11% of genes
(dominant-arm concordance 0.89 vs `none`) and de-concentrates the hub structure (Gini
0.570→0.538; dominance spread from ~43 to ~56 arms) — exactly the promiscuity control it is
for. `evidence_mass` penalizes by *credible* evidence mass rather than raw `degree` (which
de-concentrates marginally more) — a principled middle, not a uniquely optimal pick.
**Note:** the old "sharpest/most stable coupling (Basal ≥42/50)" rationale was a
coupling-axis claim that the held-out refutes; the spine's value is on attribution.

### 2.4 Many-to-many, made explicit

The two directions of the many-to-many graph are handled by **two different
mechanisms**, and conflating them is the classic mistake.

- **One miRNA → many genes** is handled by the **denominator** `D(m)` (§2.3c) and by
  **specificity** (§4.4). The denominator stops a hub from dominating; specificity
  measures how concentrated an arm's budget is on a single gene (`spec(m,g)` high =
  focused/clean readout; low = promiscuous hub).
- **Many miRNAs → one gene** is handled by the **gene-level aggregation** (`pressure(g)
  = Σ_m c`, §4.1) and by **gene-side share** (`gene_share`, §4.4) — what fraction of a
  gene's *total incoming pressure* one arm carries, with a rank (1 = the gene's
  dominant regulator). A gene like PTEN has ~91 regulators; gene-side share is how we
  say "miR-9-5p is PTEN's #1 regulator" while that edge is a trivial slice of miR-9's
  own portfolio (arm-side `edge_share` 0.006). **`edge_share` is arm-side,
  `gene_share` is gene-side — they answer opposite questions** and both are kept.

### 2.5 Membership is many-to-many too (gene ∈ N Hallmarks)

A gene that belongs to *N* Hallmark sets generates *N* `(miRNA, gene, hallmark)` edge
rows. Per-Hallmark views keep all *N*; cohort-level counts **dedup on `(miRNA, gene)`**
so a multi-program gene is not counted *N* times.

### 2.6 The pan-context limitation and the non-circular priors (overlay, not budget)

miRTarBase carries **no tissue field** — an edge is "high-evidence" if validated in
*any* cell line. Only **12.3%** of HE edges have *any* breast-specific validating
paper. This is a real limitation of the prior, and it is quantified rather than hidden:

- **Breast re-ranking** (`edge_breast_context`): mine each edge's Functional-MTI PMIDs
  for breast MeSH / mammary cell lines, re-rank a gene's regulators by
  `log1p(evidence) + β·log1p(n_breast_pmids)`. This flips the *credited top regulator*
  for ~90/783 multi-regulator genes and recovers textbook breast edges (PTEN→miR-21-5p,
  ERBB2→miR-125b-5p). **It never touches the pressure budget** — it is an
  audit/nomination overlay, kept strictly non-circular (literature only, no TCGA
  expression).
- **Multi-axis prior** (`edge_prior_refinement`): composes evidence with four more
  auditable, TCGA-independent axes — TargetScan sequence strength, breast literature,
  *negative* evidence (Non-Functional MTI contradictions; ~26% of HE edges are
  contradicted ≥once), and miRGeneDB guide/passenger status — each a **separate visible
  term** so nothing is hidden. Also an overlay, not the budget.

The lesson encoded here: **the spine prior is pan-context and we say so loudly; the
breast-specific refinements are kept separate precisely so they can be tested rather
than assumed.**

---

## 3. Expression modeling — the dynamic piece

The edge weight is static; *all* sample/state variation enters through expression.
There are four building blocks, combined into modes. (`FORMULAS.md` §3–4 is canonical.)

Notation: `x[m,s] = log2(RPM+1)` of arm `m` in sample `s`; `mean_m, sd_m, med_m` are
in-scope statistics; `R(g)` = the arms regulating `g`.

### 3.1 The four building blocks

```
(a) z-score (dynamics)        z[m,s] = (x[m,s] − mean_m) / sd_m
(b) softmax share             sm[m,s] = exp(x[m,s]−med_m) / Σ_{m'∈R(g)} exp(x[m',s]−med_m')
(c) absolute anchors          logrpm[m,s] = max(x[m,s], 0)
(d) healthy-anchored dev      dev[m,s] = x[m,s] − median_healthy(m)        (log2 units, no SD)
```

- **z (a)** — deviation of the arm from *its own* typical level across the cohort.
  Signed, unitless. *Caveat:* an arm sitting at its mean has z≈0 and contributes
  nothing — a stably-high arm is invisible under z. This is intentional for coupling
  (where *dynamics* are the signal) but wrong for attribution (§4.3).
- **softmax share (b)** — the fraction of gene `g`'s regulator "abundance budget" held
  by arm `m` in sample `s`; positive, sums to 1 across `R(g)`. This is where *the crowd
  of regulators on one gene* enters: an arm only gets share to the extent it out-abounds
  its co-regulators on that gene.
- **absolute anchor (c)** — keeps absolute level in the score so a relative share
  between two barely-expressed arms isn't treated like one between two highly-expressed
  arms.
- **healthy dev (d)** — replaces the within-tumor z with deviation from the **healthy
  (GTEx-breast) median**. No SD division, so a uniformly tumor-elevated oncomiR is *not*
  gated to ~0 the way z would gate it. This is the engine of the acquired axis (§8).

### 3.2 The three headline modes (and why each exists)

| Mode | Formula | Sign | Used for | Why |
|---|---|---|---|---|
| **`softmax_z_logrpm`** | `sm · z · logrpm` | ± | **coupling spine** | relative share × within-tumor dynamics × absolute level; `z` centres on the tumor-cohort mean → answers "*who moves across tumors*" |
| **`softmax_logrpm`** | `sm · logrpm` | + | **attribution / shares** | no `z`, so nothing is zeroed for being stably expressed; `abs_share == signed_share`; cross-state comparable |
| **`softmax_devhealthy_logrpm`** | `sm · dev · logrpm` | ± | **acquired axis** | centres on the *healthy* median → "*who carries pressure the tumor acquired over healthy tissue*" |

The reason there are three (not one) is that they answer three *different* questions,
and using the wrong one silently corrupts the answer — z for coupling (dynamics), no-z
for attribution (who carries the steady-state pressure), healthy-dev for acquisition.
This is the "Q1 attribution policy" (`FORMULAS.md` §"Attribution policy"). **Note:** these
modes build the *aggregate* (gene/program) pressure; a single-edge or single-module
coupling does **not** apply the transform — it uses the arm's raw abundance directly (§6.3).

---

## 4. Pressure — the master formula and its variants

### 4.1 The one master formula

Everything in §4–5 is a piece of, or a summary of, this:

```
pressure(g, s) = Σ_{m ∈ R(g)}  expr_mult(m, s) · edge_w(m, g)
                 └───────────┘  └────────────┘   └──────────┘
                  aggregate       expression       evidence
                  (sum over        piece (§3)       piece (§2)
                   the gene's
                   regulators)

per-arm contribution:   c(m, g, s) = expr_mult(m, s) · edge_w(m, g)
```

`c(m,g,s)` is the **edge-resolution pressure** — the atom. The gene-resolution
pressure is its sum over `R(g)`. **Every share, specificity, architecture, and family
number is built from `c`.** Default aggregate = **sum** (a gene with more, stronger
regulators genuinely has more pressure on it); `mean` is a sensitivity variant.

### 4.2 Acquired / lost pressure — at every resolution

"Acquired" pressure is pressure the tumor *gained over the healthy state*; "lost" is
pressure it shed. It is defined identically at each resolution by differencing the same
`c` across states `X ∈ {gtex, nat, tumor}` (cross-state uses the no-z `softmax_logrpm`,
because z is standardized *within* a state and is not comparable across states):

```
EDGE     gain_gtex(m,g) = c_tumor(m,g) − c_gtex(m,g)          (the edge's acquired increment)
GENE     acquired_vs_gtex(g) = pressure_tumor(g) − pressure_gtex(g)     (principal healthy anchor)
         acquired_vs_nat(g)  = pressure_tumor(g) − pressure_nat(g)      (collapse-free same-platform cross-check)
HALLMARK pro_tumor_acquired(m,set) = prior(set) · Σ_g [ −max(c_tumor−c_gtex,0) · output_influence(g) ]   (§5)
```

The two healthy anchors are complementary: **GTEx→tumor** captures the large
healthy-baseline jump (the miR-21 program, top acquired genes PCBP1/PYCR1/PTPN14/…);
**NAT→tumor** is on the *same platform* (no cross-platform artifact) and ranks the final
malignant step (HUS1/ENG/PTEN/INHBA). See §8 for the full trajectory machinery, the
two acquired *axes* (rank vs magnitude), and the GTEx collapse-repair.

### 4.3 Share — who carries a gene's pressure

Share decomposes a gene's total pressure among its regulators. There are two axes of
choice (global vs per-sample, magnitude vs signed):

```
global_abs_share(m)    = mean_s|c(m,:)|        / Σ_{m'} mean_s|c(m',:)|      ← PRIMARY RANKING
global_signed_share(m) = mean_s c(m,:)         / Σ_{m'} mean_s c(m',:)        ← COHERENCE DIAGNOSTIC
mean_abs_share(m)      = mean_s[ |c(m,s)|      / Σ_{m'}|c(m',s)| ]
mean_pos_share(m)      = mean_s[ max(c,0)      / Σ_{m'}max(c',0) ]
```

- **Rank arms on `global_abs_share`** — a single stable ratio (magnitude share of the
  gene's total pressure).
- **Diagnose with `global_signed_share`** — when it diverges from the abs version, arms
  are *cancelling* (some net-activating in a sample). Can be <0 or >1; diagnostic only.
- Under the **no-z attribution mode** every `c ≥ 0`, so `abs == signed` and the two
  collapse — which is the whole point of attributing under `softmax_logrpm`.

`share` exists in **two complementary directions** (§2.4): arm-side `edge_share` (this
gene's fraction of the arm's own budget) and gene-side `gene_share` (this arm's fraction
of the gene's incoming pressure, with a rank).

### 4.4 Specificity — focused target vs promiscuous hub

How concentrated is an arm's budget on one gene? (Two flavours, because orphans have no
curated edge to measure realized share on.)

```
realized:    spec(m,g)    = mean_s|c(m,g,:)| / Σ_{g'∈targets(m)} mean_s|c(m,g',:)|
TargetScan:  ts_spec(m,g) = ts_weight(m,g)   / Σ_{g'} ts_weight(m,g')
```

`spec` is the fraction of the arm's *total* target-mass on one gene, read against the
even-spread baseline `1/|targets(m)|` (well under 1% for a promiscuous arm). So `spec ≥ 0.05`
→ a disproportionately large, several-fold-above-baseline share = **specific/focused** =
cleaner single-gene validation (a *relative concentration*, not an absolute majority);
`spec ≤ 0.01` → at/below the floor = **promiscuous hub**. Specificity feeds
the edge classification (driver vs minor vs structural-only) and the orphan validation
priority (§9.3) — a focused edge is a cleaner wet-lab readout.

> **Seed-family non-identifiability applies to every per-arm number here** (share, gene-side
> rank, `spec`). Paralogues sharing a seed (miR-29a/b/c, miR-302/372/373/520, …) hit the same
> sites and co-vary, so a single-arm attribution may simply credit one member as a stand-in
> for its family. This is measured, not assumed: within-seed-family arm pairs co-vary ≈2×
> random pairs (Mann–Whitney p≈4e-14) and the largest multi-arm hubs are 20–37% redundant
> (`seed_family_justification`). Read share/`spec` at the **seed-family** level as the safe
> unit; the specific within-family arm is a nomination. The realized-coupling layer uses the
> same seed-family grouping for its dependence-aware FDR and permutation null (§6.x,
> `coupling_inference` / `coupling_permutation`). Singleton-seed arms are unaffected.

---

## 5. Architecture — modeling a set of genes

This is the answer to "*how do we model a set of genes and their components*" — the
hardest part conceptually, because a program is **not** a bag of independent genes. A
miRNA repressing an upstream activator changes the whole downstream program; repressing
a tumor-suppressor is pro-tumor regardless of network position. The architecture layer
(`geneset_architecture`) overlays three orthogonal lenses on otherwise-flat pressure.
(`FORMULAS.md` §12; `LANDSCAPE_REPORT.md` §7.)

### 5.1 Lens 1 — topology (network position inside the program)

Build the program's regulatory wiring from **OmniPath** (directed + signed), induced on
the set's members. Each ordered gene pair gets a **calibrated continuous sign** instead
of a hard ±1/0:

```
net_sign(s→t) = (n_stimulation − n_inhibition) / (n_stimulation + n_inhibition)  ∈ [−1,+1]
```

(Pure activator → +1, pure inhibitor → −1, conflicting → a graded value — instead of
flattening dual edges to 0, which demonstrably over-counted.) Then a finite-horizon
signed-path propagation (Katz, α=0.25, 6 hops) gives each gene its **net downstream
influence**:

```
T = Σ_{k=1..6} (α·S)^k        output_influence(g) = 1 + Σ_j T[g,j]
```

A central activator scores `>0`; an inhibitor hub `<0`. Because a miRNA *represses* `g`,
the effect on program output flips sign with the node's influence:

```
program_output_change(m) = Σ_g [ −c_tumor(m,g) · output_influence(g) ]
```

`>0` = the arm's repression nets to *activating* the program (it mostly hit inhibitor
hubs); `<0` = it *damages* the program. **This is the only lens that captures indirect,
network-mediated effects** — everything else (flat pressure, coupling, gene-role) is
local.

### 5.2 Lens 2 — tumor-direction prior (is the program pro- or anti-tumor?)

A coarse, hand-curated BC-context polarity per Hallmark (`TUMOR_DIRECTION_PRIOR`):
`+1` pro-tumor engines (EMT, HYPOXIA, GLYCOLYSIS, E2F, G2M, MYC, …), `−1`
tumor-suppressive / anti-tumor-immune (APOPTOSIS, P53, IFNα/γ, INFLAMMATORY,
ALLOGRAFT), `0` ambiguous (kept directional but not polarized).

```
pro_tumor_pressure(m, set) = prior(set) · program_output_change(m)
```

`>0` = pushes toward malignancy (activating a `+1` engine, or **releasing a brake** on a
`−1` program). The headline mechanistic finding falls straight out: because miRNAs
*repress*, they **damage** the `+1` proliferation engines (negative pro-tumor on E2F/G2M/
EMT) and turn malignant mainly by **brake-release on the `−1` suppressive programs**
(miR-21-5p → P53/APOPTOSIS, miR-375 → ALLOGRAFT/immune-evasion). All 50 sets are net
*damaged* by the miRNA layer; only the *targeted sub-pattern* is oncogenic.

### 5.3 Lens 3 — gene-role / malignancy (driver-level direction)

Topology has a blind spot for canonical oncomiRs (it scores miR-21-5p ~neutral because
it mostly damages the engines it broadly hits, and PTEN acts on a *lipid*, not a gene
node). The fix is a curated **gene-role** layer (`gene_roles.py`: 232 high-confidence
drivers — 114 oncogene / 106 TSG / 12 dual, COSMIC-CGC/OncoKB consensus) giving each
gene `malignancy_sign ∈ {−1 TSG, 0, +1 onco}`:

```
mal_pro_tumor(m)      = Σ_g [ −malignancy_sign(g) · c_tumor(m,g) ]              (graph-free driver call)
mal_pro_tumor_arch(m) = Σ_g [ −malignancy_sign(g) · c_tumor(m,g) · w_arch(g) ]  (× centrality MAGNITUDE)
```

(`w_arch` = reverse-PageRank ÷ redundancy, always ≥0, so it scales magnitude without
re-flipping the sign — an earlier version that multiplied by *signed* influence
double-counted direction and was rejected.) This **recovers the textbook oncomiR
roster** (miR-21/182/10b/93/106b/25/221) that topology missed — external validation that
the layer captures real driver biology.

> **Topology and gene-role are kept side by side, not merged.** They answer different
> questions ("what the wiring says" vs "what the driver annotation says"), and miR-21 is
> the canonical case where you must read both. **Rule: for high-breadth abundant arms
> (n_targets ≫ median) read the gene-role / coupling / acquired axes, not the topology
> sign.**

### 5.4 The continuous, all-gene complement (DepMap dependency)

The binary `malignancy_sign` covers only 232 drivers. A graded, all-gene counterpart
uses **DepMap 26Q1 CRISPR (Chronos) gene-effect over the 96-line breast panel**:

```
dep_role_weight(g)   = gene_effect_breast(g)         (<0 = essential/dependency → repress = anti-tumor;
                                                       ≳0 = TSG-like → repress = pro-tumor)
mal_pro_tumor_cont(m) = Σ_g dep_role_weight(g) · c_tumor(m,g)      (ALL genes, 530/543 arms scored)
```

CRISPR resolves the *dependency* (anti-tumor) tail well but the TSG tail weakly, so this
is the **all-gene graded onco-collateral / anti-tumor axis**; the curated `tsg_credit`
stays the better *pro-tumor* detector — complementary (corr ≈0.29). It shows the miRNA
layer lands hardest on proliferation dependencies (G2M, E2F, MYC).

### 5.5 The within-program functional role (TF hierarchy)

Topology centrality misses curated **transcription-factor** identity (only ~10% of hubs
are curated TFs; ~48% of curated TFs are below-median centrality). A Lambert-2018 TF
census adds `is_tf` / `hub_tf` (curated TF ∧ topology hub = high-confidence master
regulator) and a tunable hierarchy weight `w_hier = w_arch·(1 + 0.5·is_tf)`, wired into
TF-up-weighted versions of the malignancy and topology-damage scores. This sharpens the
"which master regulators are hit" read: P53→TP53, G2M→E2F1/2/3/MYC/MYBL2, IFNγ→STAT1.

### 5.6 Roll-ups

Each architecture score rolls up across the universe: per-arm cross-set
(`sum_pro_tumor_pressure`, `n_sets_pro_tumor`, top set), per-set totals
(`total_program_output_change`, `total_mal_pro_tumor`), seed-family aggregates, and an
**acquired-weighted parallel** of every score (weight = `max(c_tumor−c_gtex,0)`) that
isolates the *tumor-acquired* pro-tumor pressure from the constitutive (e.g. miR-21's
TSG repression is 93% acquired; miR-374a's is constitutive). The universe-level
unsupervised structure (NMF programs, subtype-divergent factors, within-pathway
sub-modules) is the subject of `LANDSCAPE_REPORT.md` §6 and is summarized in the
glossary.

---

## 6. Coupling — realized pressure and confounding

Pressure is potential. **Coupling** asks whether the repression is actually *observed*
in the cohort: does a tumor with more of the miRNA have less of the target, **after the
obvious confounders are removed**?

### 6.1 The core idea

Realized repression should show up as a **negative correlation** between the miRNA *dose*
and target expression — the dose being the arm's raw abundance at single-regulator
resolution (edge/module) and the aggregate pressure where many regulators converge
(gene/program); see §6.3. But a naive correlation is dominated by confounders (a tumor that
is more cellular, more aneuploid, more proliferative drags *everything* with it). So we
use a **partial Spearman**: regress the confounders out of *both* the miRNA and the
target, then correlate the residual ranks. A negative partial-ρ that survives = realized
repression. Spearman (rank) rather than Pearson because expression is heavy-tailed and we
want monotone, not linear, association.

### 6.2 The confounders — the deep dive

This is where most of the model's defensiveness lives. Each canonical confounder is
included for a *specific* mechanistic reason, and each is mechanistically motivated, not
just thrown in.

| Confounder | What it is | Why it confounds miRNA↔target | Source |
|---|---|---|---|
| **Purity (CPE)** | Consensus tumor-purity estimate | Stroma/immune infiltrate carries its *own* miRNAs and mRNAs; a low-purity tumor looks "high miR / low target" for purely compositional reasons | `analysis.utils.common.loaders` |
| **HRD** | Homologous-recombination deficiency / genomic scar | Genome-instability burden co-drives both arms; HRD-high (basal) tumors over-express many loci | clinical |
| **Target copy number** | ASCAT3 integer CN of the *target* gene | If the target is amplified, its mRNA rises for a reason that has nothing to do with the miRNA — masking or faking repression | `load_cnv_target_genes` |
| **Proliferation** | E2F+G2M metagene (+ two alternatives) | Cell-cycle drags both miRNA clusters (miR-17~92) and their targets; *must* be removed to avoid "cell-cycle in disguise" | metagene |
| **ESTIMATE (+epithelial)** | immune + stromal (+ epithelial + prolif) composite — `estimate_epi`, the default | The most complete composition control; the epithelial axis specifically de-confounds lineage co-expression (miR-200 family) | ESTIMATE ssGSEA |
| **Sequencing batch** | TCGA analyte-plate dummies — miRNA-seq + RNA-seq plates (`plate_both`; rare plates pooled → ~14 levels/assay); per-state NAT/GTEx batch (barcode plate / GTEx SMNABTCH·SMGEBTCH) in cross-state; TMT-plex in CPTAC | Technical run/plate co-varies features on *both* arms, so it can fabricate **or** mask coupling; adjusted as **covariate dummies — not ComBat** (ComBat over-corrects the unbalanced TCGA design and can manufacture the feature–feature correlations that *are* the readout) | `data_loaders.augment_tcga_batch` / `tcga_batch` |

Three subtleties make the proliferation handling careful rather than reflexive:

- **Three proliferation proxies**, differing in their E2F/MYC dependence: the full
  E2F+G2M metagene, an **E2F/MYC-stripped** G2M+spindle metagene, and **MKI67 alone**.
  Reporting all three tests whether repression survives proliferation adjustment
  *without over-correcting through a shared E2F path* — the "we just removed the signal"
  objection.
- **Proliferation is a *suppressor* confound, not a driver of the signal.** Within Basal
  it is positively associated with *both* pressure and target expression, so adjusting
  for it makes the negative coupling *more* negative (it was masking, not faking). This
  is why the basal coupling *strengthens* under proliferation adjustment — a key
  robustness result, not a fragility.
- **The covariate ladder.** Coupling is reported as a *nested ladder* of covariate sets,
  not a single number, so a reader sees exactly what each adjustment costs:

```
raw_rho                  none
rho_CPE_HRD              CPE, HRD
rho_CPE_HRD_CN           + target CN
rho_e2f_g2m              CPE, HRD, E2F+G2M
rho_e2f_g2m_CN           + target CN          ← "proliferation-surviving" headline
rho_mki67(_CN)           MKI67-only proxy
rho_ortho_noE2F_MYC(_CN) E2F/MYC-stripped proxy
```

Since the **2026-06-26 Normal-excluded + batch refresh**, **TCGA sequencing-batch dummies are
appended to every non-raw rung** (`augment_tcga_batch`), so each confounder-adjusted ρ above is
*also* batch-corrected (the raw rung is the only batch-free one); cross-state coupling adds the
per-state NAT/GTEx batch at the same chokepoint, and CPTAC adds its TMT-plex.

Significance flags, **BH-FDR within each `(predictor_type, predictor, scope)` family**
(this "family" is the *statistical testing family* — the set of tests sharing a correction;
distinct from the *miRNA seed family* below):

```
neg_sig_cn_fdr        = (rho_CPE_HRD_CN < 0) AND (q_CPE_HRD_CN < 0.05)     ← headline
neg_sig_prolif_cn_fdr = (rho_e2f_g2m_CN < 0) AND (q_e2f_g2m_CN < 0.05)     ← proliferation-surviving
```

A **variability gate** (`sd_x < 0.5 AND IQR_x < 0.75` on log2RPM) marks flat arms whose
"null coupling" should not be over-read.

**Two senses of "family", both handled** (`coupling_inference` / `coupling_permutation`):

- *Statistical testing family* — the BH scope above. Made dependence-robust with
  **Benjamini–Yekutieli** (`q_by`), since correlated tests violate BH's independence.
- *miRNA seed family* — paralogues sharing a seed are non-independent (justified
  empirically: within-family arm ρ ≈2× between, multi-arm hubs 20–37% redundant,
  `seed_family_justification`). Where the test unit is a miRNA (edge / target-set), each
  seed family is collapsed to one **Simes** test (`q_seedfamily`) plus a min-statistic
  family-wise null (`p_seedfamily_fwer`); the de-duplicated discovery count is in *families*,
  not edges.

**Permutation null at every rung.** Beyond the asymptotic partial-ρ p, every coupling is
re-tested by **Freedman–Lane** rank-residual permutation, applied jointly (one shared sample
permutation per draw) so seed-family dependence is preserved in the null; a smooth tail-p
from the null moments feeds the FDR (so the empirical 1/(B+1) floor never caps it) while the
discrete exceedance p is kept as the honest check. The same joint null gives the universe
"how many programs couple" count its own null. Resolutions: edge, gene, target-set, program,
universe.

### 6.3 Coupling at every resolution

The partial-ρ machinery is applied at each rung of the hierarchy (§1), against a
different predictor/target pair:

| Resolution | Predictor | Target | Coupling object | Module |
|---|---|---|---|---|
| **Edge** | arm `x[m,:]` | target `y[g,:]` | per-edge partial-ρ, per state; the `coupling_sig_pattern` (R/P/0 over GTEx/NAT/tumor) | `mirna_state_class`, `edge_partial_corr` |
| **Gene** | aggregate incoming pressure `Pagg(g)` | gene's own expression `y[g,:]` | `gene_corepression`: is the *gene* net-repressed? (282/1,424 in tumor — GATA3, ESR1, AURKB) | `decoupling_validation` §14a |
| **Target-set** | arm abundance `x[m,:]` | mean within-state z over `T(m)` (the **module composite**) | does the arm's *whole module* move down together? (miR-18a/106b/17~92) | `mirna_comovement` §13a |
| **Family** | family predictor (`mean / max / sum / mean_z` over the seed family) | target | family-as-unit coupling | `mir301_family_depth` §8 |
| **Hallmark** | per-Hallmark pressure | per-Hallmark mean target expression | count of programs with FDR-neg coupling (Basal 42/50) | `hallmark_interaction` |
| **Universe** | all programs | — | fraction of programs/edges coupled; cross-cohort concordance | roll-ups |

The **module composite** deserves emphasis as the strongest confidence-independent
signal: for each arm, z-standardize each target within state, take the (optionally
edge-weighted) mean over targets, then one partial-ρ of arm abundance vs that composite.
A strongly negative composite-ρ means the *coordinated* module is repressed — something
no single-edge ρ can see. (Edge-weighted ≈ equal-weight, r=0.97–0.99, so equal-weight is
the default.) A **role-stratified** version (TSG vs oncogene targets) is an interpretive
*overlay*, deliberately not folded into one number — folding role into the composite would
re-conflate "is it repressed" with "is it pro/anti-tumor," which is the architecture
layer's job (§5).

**Predictor switches with resolution — and only the share matters.** The edge and
aggregate predictors are not a reparametrization, but most of the difference is inert: the
static edge weight is a sample-constant scalar (rank-invisible) and `z` is affine in
abundance (rank-identical for one arm), so neither moves a single-edge coupling; the
**softmax share** (which pulls in co-regulators) is the sole differentiator. A head-to-head
check — the same CPE+HRD+target-CN partial-ρ under every predictor (`coupling_predictor_comparison`) —
confirms it: at the **edge** level `z` reproduces abundance (sign-concordance 0.96, recovers
94% of FDR-sig edges) while share-alone recovers only ~49% and flips its median positive →
**raw abundance is the most sensitive single-edge predictor** (why edge/module coupling uses
it). At the **gene** level share-alone is degenerate (per-gene shares sum to ≈1 each sample
→ near-constant), but the full pressure recovers ~89% of the naive abundance-sum's hits and
finds **~2–3% more** with a deeper median → the construction folds many regulators into one
dose better than summing abundances. Net: pressure differs from abundance only through the
share, which *weakens* one edge but *strengthens* the gene aggregate (MH-44). The one remaining
*linear* lever — conditioning the target on the competitive field `P_others` to denoise the focal
edge — was tested over the full universe (E4, MH-78) and **does not beat raw abundance** (precision
1,471→1,369 neg-sig, median Δρ≈0.0001); the field is not a target-noise source, so the
abundance-is-best headline survives. The full edge-question frame (the coupling-invariance lemma:
*only* the softmax, per-sample modulation, sample-set/paired design, the conditioning set, and the
functional form can move a single-edge rank coupling — reference/scale are inert) is in
`EDGE_QUESTION_TAXONOMY.md`, with skill `apm-edge-question`.

At the **gene aggregate the lemma flips** — Σ breaks per-arm monotonicity, so construction *is*
consequential. The gene construction grid (MH-79; `GENE_QUESTION_TAXONOMY.md`, skill
`apm-gene-question`) finds a **NAT-(acquired-)referenced softmax share beats naive abundance-sum by
+12.6%** (448 vs 398 net-repressed genes) — the *reference* is the winning axis (the function barely
matters; z-sum ≈ baseline), while **AGO-gating the aggregate HURTS** (the edge proliferation-confound
generalizes) and pruning is neutral-to-harmful. (This supersedes the earlier "~2–3% more" gene-level
figure: that compared only `softmax_z_logrpm` to abundance-sum; the NAT-referenced share is the larger
gain.) Net-repression concentrates on known cancer genes (TSG 52% / onco 42% vs dual 27%), with ~35%
of net-repressed genes carrying an incoherent (canceling) stack.

### 6.4 The potential × realized join (pressure × coupling; the headline classes)

The whole point of measuring *both* potential and realization is to **join them**. Here
*pressure* is the **potential** axis (its edge-level trajectory tracks the miRNA's
abundance, the weight being static) and *coupling* is the **realized** partial-ρ — the join
is potential × realized, **not** a correlation of one against the other. Every edge carries
a pressure (potential) trajectory (does the miRNA's abundance gain/lose over states) and a
coupling (realized) R/P/0 trajectory. Their join is the headline class:

```
acquired_realized          pressure rises AND coupling turns repressive (00R / P0R)   ← genuinely switched on
field_established_realized  arm acquired + 0RR (repressive in NAT, kept in tumor)
acquired_unrealized         pressure up but coupling never realized
lost                        pressure flat + brake released (R00: repressive in healthy, gone in tumor)
constitutive               repressive in both healthy and tumor (R0R)
nat_decoupled               NAT-only sign-flip vs both states (essentially absent as an FDR class)
```

This join is what lets us say, e.g., *"healthy breast actively holds ERBB2/CXCR4/MMP2/
NOTCH2 under miRNA repression and the tumor releases the brake"* (the `R00` released-brake
biology) as a **coupling-derived** claim, immune to the pressure-side caveats. See §8 for
how the state axis is built and `LANDSCAPE_REPORT.md` §3–4 for the full class biology and
the de-coupling validation (CN/methylation/ATAC conditioning).

---

## 7. Protein realization — beyond the mRNA layer

miRNA repression ultimately acts on the **protein** product. The mRNA-level coupling
(§6) sees destabilization; it cannot see pure translational block. So the framework tests
whether pressure realizes onto **protein**, two ways, using CPTAC mass-spec proteomes.

### 7.1 Three layers

For a target gene with mRNA `rna_z` and protein `protein_z`:

```
L2  (direct protein)       partial-ρ( pressure , protein_z )                          expect < 0
L1  (slope-1 RNA→protein)  partial-ρ( pressure , rna_z − protein_z )                  expect > 0 (translational repression widens the gap)
L1b (protein_resid)        partial-ρ( pressure , protein residualized on its own mRNA) expect < 0 = genuine BEYOND-mRNA effect
```

**L1b is the key idea — "realization onto protein *through* vs *beyond* the mRNA."**
Residualizing protein on its own mRNA removes everything the *measured* mRNA explains; a
surviving negative ρ is a protein-level effect the bulk steady-state mRNA does not capture —
best read as *beyond-(measured)-mRNA*, **not** strictly *translational*. It is consistent
with a translational block but also with destabilization bulk RNA-seq misses (deadenylation /
poly(A) shortening; 3'-UTR-isoform / APA effects), protein-vs-mRNA turnover kinetics (protein
integrates past mRNA), or indirect routes — plus platform noise between sequencing and MS.
This is the "with or without the miRNA–mRNA coupling" distinction: L2 reads protein *together
with* whatever happened at the mRNA layer; L1b reads protein *after subtracting* the measured
mRNA.

### 7.2 The two CPTAC cohorts and the headline

- **CPTAC TCGA-105** (Mertins, iTRAQ): the *same* TCGA patients re-assayed. Protein
  coupling is essentially null (0/1006 genes; L1b 0) — the noisier same-patient iTRAQ
  cannot see it.
- **CPTAC prospective-122** (Krug, TMT): *independent* patients, with its own GDC
  miRNA-seq, so pressure is rebuilt **self-contained** (no cross-cohort matching). Here
  the spine **validates at protein**: 34/1143 genes FDR-neg on L2 (led by **ZEB1**, the
  canonical miR-200 target, ρ=−0.55), 6/50 Hallmarks (EMT top).

The mechanistic conclusion: **the repression is predominantly mRNA-mediated.**
Residualizing on mRNA collapses most of it (ZEB1 −0.55 → −0.24 NS), leaving only ~8
genes with a genuine protein-beyond-mRNA (L1b) component — canonical mammalian miRNA
biology (destabilization dominates, small translational residual). The L1b survivors are
**known, subtype-structured edges**, and critically they **recover the basal miR-17~92→
PTEN hub at the protein layer** (miR-17-5p ρ=−0.34, q=0.008) — the strongest external
corroboration of the subproject's central hub thesis.

### 7.3 Target specificity (the decoy control)

A reviewer worry: maybe the survivor arms are just *globally* anti-correlated with the
whole proteome (a confound). The **decoy control** (`cptac_target_specificity`) compares
each driver arm's coupling to its *cognate* protein against its coupling to every protein
it does **not** target (≥50 decoys). The cognate couplings (median ρ=−0.33) sit ~16× more
negative than the decoy null (≈−0.02, centred at zero): 9/13 FDR drivers are
target-specific (empirical p<0.05), and the 2 borderline cases are exactly the broad
proliferation oncomiRs (PTEN←miR-17-5p) whose decoy nulls *are* somewhat negative — i.e.
the control correctly discounts genuinely-promiscuous arms. The specificity survives the
literal TMT-plex batch adjustment.

---

## 8. The healthy → NAT → tumor trajectory (acquired pressure)

The "acquired" concept (§4.2) is given its own machinery (`mirna_state_class`) because
*how* a miRNA's pressure moves across the three tissue states is itself a finding.

### 8.1 The three states and the healthy anchor

```
GTEx-healthy breast   →   TCGA-NAT (normal-adjacent)   →   TCGA-tumor
       │                          │                              │
   true healthy            field-effect / pre-malignant      malignant
```

The **healthy anchor** is GTEx breast (346 samples), **quantile-normalized onto the TCGA
scale** and joined by **MIMAT accession** (not name). GTEx's uniquely-mappable pipeline
*collapses* ~278 canonical arms to zero (let-7a, miR-30a, miR-182…); two imputation
layers repair coverage:

1. **miTED anchor rank-transfer** (primary): arms measured in both GTEx and DIANA-miTED
   define a monotone map; a gap arm's GTEx-equivalent is predicted from its miTED value.
   Validated leave-one-out (Pearson 0.90, median error 0.86 log2).
2. **GTEx seed-family** (fallback). floor0 = no reference anywhere (left un-imputed and
   flagged, never fabricated — imputing a truly-absent arm would *erase* a real gain).

**GTEx collapse-repair** repairs the GTEx abundance matrix *at source* before any
contribution, so `c_gtex` / `dHT` / `acquired_vs_gtex` are collapse-free (the let-7a/
miR-30a/miR-182 collapse artifacts no longer register as fake gainers;
`frac_gain_imputable_artifact = 0/1,410`). **Coupling is left on the raw GTEx matrix**
(it re-derives its own bundles), so partial-ρ is untouched. The imputed baseline feeds
**only the acquired/dev axis** — never the z-spine or coupling — so all headline coupling
findings are baseline-independent.

### 8.2 Two acquired axes — rank and magnitude

A single rank delta is not enough, because rank **saturates**:

```
RANK      dHT(m) = rank_tumor(m) − rank_gtex(m)         (percentile-rank delta on the global pressure budget)
MAGNITUDE log2fc_tumor_vs_healthy(m) = tcga_tumor_median(m) − healthy_baseline_tcga(m)   (same QN'd scale)
```

`dHT` is cross-platform-robust but an already-abundant arm *cannot gain rank* (max dHT
≈0.08 in the 80–95% bin), so **81% of >4-fold-up arms are rank-invisible** — precisely the
BC oncomiR surge (miR-182/183/375/93/141, +9 to +13 log2) is `magnitude_only`. The
magnitude axis recovers them. **The default acquired readout is both axes:** `acquired_
gainer = rank-gain OR magnitude-gain`; 242/721 arms gain on ≥1 axis (69 rank-only, 141
magnitude-only, 32 both).

### 8.3 Trajectory classification, at three resolutions

```
miRNA (8.3a)   on (dHN, dNT, dHT):  stable / progressive_{gain,loss} / tumor_acquired / field_effect / non_monotonic
edge  (8.3c)   join of the miRNA's pressure-move with the per-edge coupling trajectory → joint_edge_class (§6.4)
gene  (8.3e)   acquired_vs_gtex / acquired_vs_nat with the floor0-artifact split as a residual diagnostic
```

The three legs are `dHN = nat−gtex` (field effect), `dNT = tumor−nat` (tumor-specific
increment), `dHT = tumor−gtex` (the **acquired headline**). A **subtype-specificity**
secondary class (Tau over PAM50) tags whether an acquired arm is pan or single-subtype
(§ glossary; `LANDSCAPE_REPORT.md` §5). The whole §4.2 acquired story is read off this
trajectory — and because it is anchored to the GTEx *state* (never NAT, which sits
tumor-ward by field effect), "acquired over healthy" means exactly that.

---

## 9. Validation — recovered, deepened, novel

Validation is framed by *what kind of claim* a result is, using the registry tags
(`S` statistical+FDR+adjusted / `P` pattern / `H` hypothesis / `R` reused-known). The
three questions are: *did we recover known biology?* *did we deepen it?* *did we find
anything new?*

### 9.1 Recovered — the positive controls

The model reproduces textbook biology *without being told to*, which is the evidence it
is real:

- **The TS-miRNA roster** falls out of the architecture anti-tumor leaders (miR-30a,
  miR-24, miR-29b/c, let-7a) and the gene-role oncomiR roster (miR-21/182/10b/93/106b).
- **miR-29→collagen/ECM**, **miR-17~92→PTEN/replication**, **let-7/miR-30→mitotic**,
  **miR-375 as immune-evasion oncomiR** all emerge as sub-modules and grip-change leaders
  in the within-pathway NMF — internal validation (if they had *not* appeared the
  decomposition would be suspect).
- **HE edges couple more than orphan edges** reproduces across cohorts (TCGA and Buffa).
- **Cross-cohort replication (Buffa 2011, 207 independent patients):** Buffa-vs-TCGA
  partial-ρ concordance +0.32; 67% of TCGA neg-sig edges keep a negative sign — an
  independent-patient leg the same-patient CPTAC cannot give.

### 9.2 Deepened — same players, new resolution

The platform's genuine methodological contribution is **state-resolved and
within-program** structure that standard GSEA / single-edge analyses cannot see:

- **State-resolved coupling** (`R00` brake-release, `00R` acquired-realized): the
  literature asserts the miRNA→target *pair* and the miRNA's down-regulation, but rarely
  measures *the brake being engaged in normal tissue and released in tumor* (strength
  `P/H`).
- **Within-program sub-modules**: a Hallmark resolves into a shared pan-subtype backbone
  + named subtype-specialist sub-modules (Basal CD69/ICAM1 immune-restriction in IFN-γ,
  PTEN-repression localized to a *minority* PI3K/P53 sub-module) — "same known players,
  newly resolved organization."
- **Protein-layer confirmation** of the basal miR-17~92→PTEN hub in an independent
  proteome (§7).
- **CN→regulatory-network attribution** (MH-43): the CN-driven component of miRNA dosage
  reproduces ~half-to-two-thirds of the bulk repressive network and survives the full
  aneuploidy+proliferation control stack.

### 9.3 Novel / candidate-novel — the orphan discovery surface

The cleanest discovery surface is the **orphan edges**: TargetScan/ENCORI-predicted
interactions **absent from miRTarBase functional curation** — kept out of the spine
precisely so they can serve as an unbiased discovery screen. The funnel:

```
20,663 orphan edges screened (TargetScan ∪ ENCORI, miRTarBase-functional-MTI-absent)
   → 15,512 testable against CPTAC protein
      → 539 protein anti-correlation candidates
         → 182 translational (protein_resid, beyond-mRNA)
            → 93 translational AND TCGA-replicated   ← strongest tier (76/93 fully uncurated · 29/93 carry both priors — overlapping subsets, not a partition)
```

Three convergent strengthenings make these more than a prediction list:

- **Method self-validation:** the **miR-29→collagen/ECM axis emerges de novo**
  (miR-29a→COL11A1 ts_weight 0.80, clipExpNum 15; →COL6A2, →BMP1) — a textbook axis
  uncurated *at these specific targets*. It is the credibility anchor.
- ⛔ **Triple-cohort validation (MH-38) — RETRACTED (MH-114, 2026-07-12). DO NOT CITE THE 108.**
  The screen that produced the 594 candidates had **no cell-composition block**; re-run adjusted, it
  gives **594 → 23 candidates, with ZERO collagen and ZERO ECM** — the miR-29→ECM axis that *was*
  this claim's headline (34 of the 108) does not survive. `output/buffa_validation/*` is still dated
  **2026-06-26 and has not been re-run**, so the figures below are computed from the retracted input.
  Corroborated independently by **MH-55**, 16 days earlier. Full retraction: registry **MH-38**.
  *Superseded text, kept for provenance:* of the **492 genuine orphans** (no miRTarBase row
  at all), **127** were CPTAC-protein-significant, TCGA-replicated, *and* testable in the
  independent Buffa cohort; of those **108 (85%) keep a negative sign in Buffa** — far above
  the 50% chance rate — giving **108 triple-cohort-validated orphan edges** (TCGA mRNA ·
  CPTAC protein · Buffa mRNA; 47 reach FDR in Buffa), led by miR-29→ECM (34/108 ECM targets).
- **Composition-robust core (MH-39) — the 17 TargetScan-only orphan edges:** re-testing the
  prediction-only tier (TS-sequence-only, no CLIP, no curation) after adding PAM50 + ESTIMATE
  immune/stromal + endothelial + target-CN covariates leaves **17 of 24 survivors holding**
  (negative beyond-mRNA `protein_resid` at q<0.1 under the widest covariate set; **11/17** at
  q<0.05). The 7 drops are exactly the composition artifacts — VCL ×2, TPM4, and the
  hematopoietic miR-142-3p/-5p cluster (WASL/STAM2/RNH1/HOOK3), which load on the CAF/immune
  axes. The 17 surviving edges, grouped by target (protein ρ = prospective CPTAC, base
  purity+CIN; resid ρ = protein residualized on its own mRNA under the *widest* covariate
  stack, with FDR q):

  | Target (function) | TS-only arm(s) (ts_weight) | protein ρ | resid ρ (q) |
  |---|---|---|---|
  | **PMM1** — mannose metabolism | miR-195-5p (0.27) · miR-16-5p (0.27) · miR-15a-5p (0.30) | −0.40 · −0.39 · −0.39 | −0.40 (9e-4) · −0.31 (.012) · −0.30 (.012) |
  | **PPM1A** — PP2C phosphatase (TGF-β/AKT brake) | miR-19b-3p (0.26) · miR-19a-3p (0.26) | −0.30 · −0.32 | −0.28 (.016) · −0.23 (.042) |
  | **PURA** — transcription / replication | miR-532-5p (0.26) · miR-452-5p (0.36) | −0.45 · −0.33 | −0.19 (.093) · −0.18 (.093) |
  | **KIF1B** — kinesin, pro-apoptotic | miR-15a-5p (0.45) | −0.37 | −0.30 (.012) |
  | **NCK2** — signalling adaptor | miR-142-3p (0.49) | −0.07 | −0.31 (.012) |
  | **VCAN** — versican (ECM) | miR-107 (0.27) | −0.52 | −0.27 (.018) |
  | **MAPRE1** — EB1, microtubule +end | miR-429 (0.27) | −0.25 | −0.27 (.018) |
  | **NUDT4** — Nudix hydrolase | miR-183-5p (0.82) | −0.16 | −0.26 (.023) |
  | **RAB1A** — ER→Golgi traffic | miR-101-3p (0.59) | −0.33 | −0.24 (.041) |
  | **CERCAM** — endothelial adhesion | let-7g-5p (0.34) | −0.35 | −0.22 (.057) |
  | **SKP1** — SCF ubiquitin-ligase core | miR-452-5p (0.29) | −0.38 | −0.21 (.066) |
  | **DR1** — NC2 transcriptional repressor | miR-182-5p (0.29) | −0.34 | −0.19 (.093) |
  | **GLUD1** — glutamate dehydrogenase | miR-30e-5p (0.27) | −0.28 | −0.18 (.093) |

  Two **seed-family convergences** lead and are CAF-neutral / method-independent (they hold
  even with target-gene CN partialled out): the triple-hit **miR-15/16/195→PMM1** (mRNA-flat,
  protein-strong = translational) and **miR-19a/b→PPM1A**. **NCK2** is the instructive single
  case — its *direct* protein ρ is ≈0 yet the **beyond-mRNA residual** is strong (−0.31,
  q=.012), a translational-layer signal bulk mRNA cannot see. All 17 remain **wet-lab
  nominations** (a prioritized CLIP/luciferase queue), not validated edges.

> **What these are and are NOT.** They are **wet-lab nominations** (a prioritized CLIP/
> luciferase queue), not validated edges — protein anti-correlation ≠ direct binding,
> co-expression confound persists (especially ECM↔stroma), and seed-family arms
> (miR-29a/b/c) cannot be separated. The miR-29→collagen recovery is the anchor; the
> fully-uncurated dual-prior translational hits and the 17 composition-robust prediction-
> only edges are the discovery surface.

### 9.4 The validation ladder, summarized

| Rung | Evidence | Example | Strongest tag |
|---|---|---|---|
| Reproduces known | recovers textbook edges unprompted | miR-29→collagen, miR-17~92→PTEN | `R` |
| mRNA coupling | FDR-neg partial-ρ, confounder-adjusted | Basal 42/50 Hallmarks | `S` |
| State-resolved | brake-release / acquired-realized trajectories | ERBB2/MMP/NOTCH2 R00 | `P`/`H` |
| Protein layer | independent CPTAC proteome | PTEN←miR-17-5p protein-resid | `S` |
| Independent cohort | Buffa replication | +0.32 concordance; ~~108 triple-validated~~ | ⛔ **RETRACTED** (MH-114 — orphan screen composition-blind, 594→23, zero ECM; not re-run. The **+0.32 concordance** is unaffected. See registry MH-38.) |
| Wet-lab queue | nominations for CLIP/luciferase | 17 composition-robust orphans, 76 uncurated | `H` |

---

## 10. Cross-cutting design principles and honest limits

Restating the stance (§0.3) as it cashes out, plus the standing caveats
(`LANDSCAPE_REPORT.md` §8):

- **Static weight × dynamic abundance** keeps "what's possible" separate from "what this
  tumor expresses."
- **Non-circular prior** — edge weights / HE filter / breast re-ranking / gene roles
  never use TCGA expression, so they cannot leak into the coupling they are tested
  against. **Checked, not just asserted:** the pressure construction's candidate ranking
  reproduces on held-out patients (5-fold CV, train→test rank-ρ 1.0), so it is not overfit
  to the cohort it is validated on (`held_out_tuning`, MH-47).
- **Two measurements, then their join** — pressure (potential) and coupling (realized)
  are never conflated; the joint class is the headline.
- **Adjusted *and* raw, gated *and* ungated, everywhere** — differences read as
  sensitivity. The AGO/RISC gate in particular is a *documented sensitivity layer, not a
  causal model* (it rescales but rarely reorders): capacity is the **co-limiting minimum**
  of a loadable-AGO term (AGO2-weighted) and a TNRC6/GW182 effector term — both required —
  and is floored so it can attenuate, never erase (`METHODS.md` §4).
- **Honest rungs** — potential ≠ realization ≠ causation; every claim is tagged.

Standing limits, stated plainly:

- **Gene-gene topology cannot see lipid/metabolite-mediated inhibition** (PTEN→PIP3) — a
  true ceiling; the gene-role layer is the partial workaround.
- **De-coupling mechanism is not fully pinned**: CN (0/220) and promoter methylation
  (1/188) are ruled out, ATAC is mostly null (a 4/134 accessibility-linked minority); the
  leading remaining candidate is 3′UTR APA/shortening + ceRNA/RBP competition.
- **floor0 arms** (no healthy reference) are a coverage gap, not acquisition — resolve by
  extending GTEx/miTED coverage, never by imputing absent arms.
- **The prior is pan-context and coarse** (one BC polarity per set); breast refinement is
  an overlay with no external breast-miRNA cohort yet to validate utility.
- **The pressure construction is one of several near-equivalent choices** — non-overfit
  out-of-sample but not uniquely optimal: bare un-normalized weighting is marginally ahead
  on coupling, and the promiscuity denominator is coupling-neutral, earning its place on
  attribution (dominant-regulator credit + hub de-concentration), not coupling (MH-47).
- **Bulk-RNA / bulk-protein ceiling** — everything is association, not binding;
  composition (stroma/immune) is controlled but never fully eliminated.

---

## Appendix: terminology glossary

| Term | One-line definition | Where |
|---|---|---|
| **Pressure** `c / pressure(g,s)` | potential repression = expression × edge weight, summed over a gene's regulators | §4.1 |
| **Edge weight** `edge_w(m,g)` | static evidence strength ÷ promiscuity denominator | §2.3 |
| **Contribution** `c(m,g,s)` | edge-resolution pressure atom; everything is built from it | §4.1 |
| **Acquired / lost pressure** | tumor pressure minus a healthy anchor (GTEx or NAT), at edge/gene/hallmark | §4.2, §8 |
| **Share** (`global_abs_share`) | a regulator's fraction of a gene's total pressure magnitude (arm-side vs gene-side) | §4.3 |
| **Specificity** `spec / ts_spec` | how concentrated an arm's budget is on one gene (focused vs promiscuous hub) | §4.4 |
| **Coupling** (partial-ρ) | observed, confounder-adjusted anti-correlation = *realized* pressure | §6 |
| **Module composite** | partial-ρ of an arm vs the mean-z of its whole target set (coordinated repression) | §6.3 |
| **Confounders** | CPE (purity), HRD, target CN, proliferation (3 proxies), ESTIMATE+epi | §6.2 |
| **Joint edge class** | the pressure × coupling join (acquired_realized, lost/R00, constitutive, …) | §6.4 |
| **Architecture** | topology (`output_influence`) + tumor-direction prior + gene-role malignancy on a set | §5 |
| **`pro_tumor_pressure`** | prior × propagated program-output change (brake-release vs engine-damage) | §5.2 |
| **`mal_pro_tumor`** | gene-role driver call: −malignancy_sign × pressure (curated; +DepMap continuous) | §5.3–5.4 |
| **L2 / L1 / L1b** | protein-direct / RNA-protein gap / protein-residualized-on-mRNA (beyond-mRNA) | §7.1 |
| **dHT / log2fc** | acquired rank axis (saturating) / acquired magnitude axis (recovers oncomiR surge) | §8.2 |
| **Orphan edge** | TargetScan/ENCORI-predicted, miRTarBase-uncurated — the discovery screen | §9.3 |
| **AGO/RISC gate** | per-sample bounded multiplier; capacity = min(AGO2-weighted load, TNRC6 effector) — co-limited, floored, sensitivity layer not mechanism | `METHODS.md` §4 |
| **NMF program** | unsupervised pressure factor (subtype-divergent programs vs shared acquired backbone) | `LANDSCAPE_REPORT.md` §6 |
| **Tau (subtype specificity)** | Yanai index of how single-subtype an arm's acquired pressure is | `LANDSCAPE_REPORT.md` §5 |
| **Strength tags** | `S` statistical / `P` pattern / `H` hypothesis / `R` reused-known | §9 |
```
