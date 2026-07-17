# Learning the miRNA→gene regulatory matrix M (design note)

**Status:** design / parked. Not in any spine. Companion to `MODELING_FRAMEWORK.md`
(the fixed-M pressure model), `FORMULAS.md`, `EDGE_QUESTION_TAXONOMY.md`,
`GENE_QUESTION_TAXONOMY.md`, and `PRESSURE_FUTURE_OPTIONS.md`.

Origin: a mathematician's framing of the edge as a matrix `M` mapping a miRNA
abundance vector `x` to a gene-expression vector `y`. This note records the
correct mathematical object, what it would add **over both prior literature and
our own fixed-M model**, and exactly where the prior, the confounders, and the
many-to-many structure (seed-family overlap, cooperative proximity) enter.

---

## 0. The object, stated precisely

```
X  ∈ ℝ^{S×A}   samples × miRNA-arms      (x[m,s] = log2(RPM+1), our spine input)
Y  ∈ ℝ^{S×G}   samples × target genes    (target expression)
C  ∈ ℝ^{S×K}   samples × confounders     (CPE, HRD, target-CN, prolif, batch)
M  ∈ ℝ^{G×A}   genes × arms              (the regulatory matrix = our edge_w)

        Y  =  −X·Mᵀ  +  C·B  +  E
```

`M_{g,m} ≥ 0` is the repressive strength of arm `m` on gene `g` (sign carried by
the leading `−`). **This is not new to us — it is our model written down.** Our
`pressure(g,s) = Σ_m expr_mult(m,s)·edge_w(m,g)` is row `g` of `−X·Mᵀ` evaluated
at sample `s` (`MODELING_FRAMEWORK.md` §4.1). The only difference the mathematician
proposes: **we currently *import* M from miRTarBase/TargetScan; he proposes to
*learn* it from the (X, Y) pairs.**

Two corrections to the framing he offered (recorded so they are not re-litigated):

- **Not bilinear in (x, y).** `y` is the outcome, not an input. The honest name is
  **multivariate multiple regression**; with a low-rank constraint on M,
  **reduced-rank regression (RRR)**. Bilinearity re-enters legitimately only in
  (a) a low-rank factorization `M ≈ U·Vᵀ` (bilinear in latent comodule factors —
  this is the Theia/PIMiM/NMF family) and (b) a feature model for M's *entries*
  `M_{g,m} = φ(g)ᵀ Θ ψ(m)` (bilinear in gene-features × miRNA-features — the
  factorization-machine reading of TargetScan-style scoring), and (c) **cooperative
  product terms** (§4), which are genuinely bilinear *in x*.
- **Not an HMM.** A Markov/HMM needs a state sequence (time or trajectory). TCGA is
  cross-sectional. The only Markov structure we have is the coarse 3-state
  GTEx→NAT→tumor chain, and we already model it as *state-differenced M·x*
  (`MODELING_FRAMEWORK.md` §8), not as learned transitions. HMM is the wrong tool
  here unless pseudotime/longitudinal data appears.

**Nesting property (the single strongest selling point):** freeze every coefficient
at its prior and the learned model collapses *exactly* to the current fixed-M
pressure. So the learned model is a strict generalization — fixed-M is the
infinite-penalty / point-mass-prior limit. Anything it finds is "the prior, plus
what the data insists on adding."

---

## 1. Where the prior goes — **not just non-zeroness**

This is the central question, and the answer is that a curated edge can inform the
estimate in **four** distinct places, only the first of which is "support":

| # | Channel | What enters | Mechanism |
|---|---|---|---|
| 1 | **Support / inclusion** | *which* `M_{g,m}` may be nonzero | hard mask to HE edges, or spike-and-slab inclusion prob `π_{g,m}` = f(evidence) |
| 2 | **Penalty weight (magnitude shrinkage)** | seed quality, assay-directness, evidence counts | adaptive/weighted L1: `λ_{g,m} = λ / w_prior(g,m)`. Strong edge → small penalty → free to be large; weak edge → heavy shrinkage |
| 3 | **Prior mean (location)** | the curated edge weight itself | Bayesian `M_{g,m} ~ N(μ_{g,m}, τ²)`, `μ_{g,m}` = our `edge_w`. Pulls toward the literature magnitude, not just toward 0 |
| 4 | **Sign / box constraint** | repression directionality | `M_{g,m} ≥ 0` (non-negativity), encoding "miRNAs repress" as a hard constraint |

So **yes — the prior on seed quality, Functional-MTI evidence counts/quality,
reporter/protein vs binding-only directness, and ENCORI CLIP depth all enter as
the magnitude prior (channels 2–3), not merely as a 0/1 mask.** This is exactly our
`confidence_logclass` numerator and `edge_w` denominator (`MODELING_FRAMEWORK.md`
§2.3) repurposed: the quantity we currently *use as the fixed weight* becomes the
*prior mean and inverse penalty* of the learned weight. TargetScan context++ and
miRGeneDB guide/passenger status enter the same way (and the bilinear feature model
in §0(b) is the way to let them inform entries we have no curated evidence for —
the orphan extension).

**This is the first differentiator over prior work.** The LASSO miRNA-target
papers (Lu 2011) use sequence prediction only for channel 1 (support). GenMiR++
uses sequence + a flat target indicator. None of them push a *quality-weighted,
assay-directness-graded* evidence score into channels 2–3. We already have that
score, non-circular and breast-audited.

---

## 2. Where the confounders go — the matrix view's cleanest gift

`C·B` is an additive block in the same regression (`MODELING_FRAMEWORK.md` §6.2
confounders: CPE, HRD, target-CN, the three proliferation proxies, sequencing
batch). **M is only identified after projecting C out of both X and Y** — i.e. M
is estimated in the residual space `(I − P_C)`. This is precisely our
partial-Spearman ladder, done jointly and on all entries at once instead of one
edge at a time on ranks.

Restated for the mathematician in his language: a confounder that drives both a
miRNA and a gene (proliferation does) creates **omitted-variable bias** in M unless
it is in the design. That single sentence is our entire §6.2, and it is the
**second differentiator**: essentially none of the prior expression-network methods
adjust for tumor purity / proliferation / target-CN / batch. Composition
confounding is the #1 failure mode of these networks and we control it by
construction.

Two ways to run it, both honest:
- **Conditioned (regression):** put `C` in as unpenalized covariates. Bias-corrects
  for *measured* confounders.
- **Instrumented (causal):** use the **miRNA-locus copy number as an instrument**
  for arm abundance (§5). Corrects for *unmeasured* confounders too. We already
  compute miRNA-locus CNV (`mirna_locus_cnv.py`) and have never used it as an
  instrument — only target-CN as a confounder. This is the one genuinely new
  estimator the matrix view unlocks on our existing data.

---

## 3. Seed-family overlap → collinearity → **group structure** (identifiability)

Paralogous arms sharing a seed (miR-29a/b/c, miR-200 family, miR-17~92 members)
hit the same sites and co-vary — we have measured this (within-family arm ρ ≈ 2×
between-family, hubs 20–37% redundant; `seed_family_justification`,
`MODELING_FRAMEWORK.md` §4.4). In a *joint* regression that documented caveat
becomes a hard numerical problem: **near-collinear columns of X**. Plain LASSO will
arbitrarily keep one family member and zero the rest, unstably across folds — the
non-identifiability re-surfacing as variable-selection instability.

The fix turns the caveat into a **modeling constraint**, which is the **third
differentiator**:

- **Group LASSO at the seed-family level** — all arms of a seed family form one
  group, in-or-out together; magnitude shared/distributed within the group. This
  literally implements "read share at the seed-family level as the safe unit; the
  within-family arm is a nomination."
- **Or collapse the family to one predictor** (mean/max/sum) — our existing
  family-coupling rung (`MODELING_FRAMEWORK.md` §6.3, `mir301_family_depth`).

Either way, the estimand is the **family→gene** weight; the per-arm split is
explicitly a within-group attribution, not an identified quantity. Per-edge methods
that ignore this report unstable single-arm coefficients as if identified.

---

## 4. Co-targeting vs cooperative proximity — **two opposite phenomena, don't conflate**

The question lumps "same/overlapping seed" with "synergy when in close proximity."
They pull in **opposite directions** and need different terms:

| Phenomenon | Site relationship | Effect on repression | Model term | Prior that triggers it |
|---|---|---|---|---|
| **Seed redundancy / competition** | same or overlapping site | **sub-additive** (substitutable; can't both bind) | group structure (§3) + optional per-site saturation | shared seed / overlapping site coords |
| **Cooperativity** | *distinct* sites in close proximity (~8–40 nt) | **super-additive / synergistic** | **product interaction** `+ M_{g,(m,m')}·x_m·x_{m'}` | TargetScan site coordinates within a window |

Key points:

- A plain linear `Σ_m M_{g,m}x_m` is **purely additive — it cannot represent
  either** synergy or saturation. Cooperativity needs an explicit **bilinear
  product term in x** (this is the only place a product of two *miRNA abundances*
  is biologically correct — and it answers the mathematician's "bilinear" instinct
  with a real mechanism, distinct from the latent-factor bilinearity of §0).
- **The prior decides the interaction support.** We do not fit all `A·(A−1)/2`
  products (hopeless, p ≫ n). We add a product term only for arm pairs with *distinct*
  predicted sites within the cooperativity window on the same 3′UTR — drawn from
  TargetScan site coordinates. So the prior is doing double duty: defining the
  *main-effect* support (§1) **and** the *interaction* support.
- **Global saturation we already model**, separately: the AGO/RISC capacity gate
  (`ago_gate.py`) is a cohort-level saturating nonlinearity. Local same-site
  competition (§3 sub-additivity) is a *different*, per-site effect and would need
  its own grouping/saturation — flag, do not assume the AGO gate covers it.

This pair (redundancy as group structure; cooperativity as primed product terms) is
the **fourth differentiator** and is, as far as the literature scan found, not done
in the expression-based miRNA-network methods.

---

## 5. Scope: gene-focused vs program-wise vs universe-wise

All three are the same M, restricted to different `(genes × arms)` blocks. They
trade interpretability against statistical strength:

| Lane | Block | Structure | Use it for | Cost |
|---|---|---|---|---|
| **Gene-focused** | one gene `g`, its regulators `R(g)` | one penalized regression (one row of M) | deep single-driver story (PTEN ~91 arms, GATA3, ESR1, ZEB1) — the cleanest validation targets | loses cross-gene strength |
| **Program-wise** (recommended primary) | one Hallmark's genes + their regulators | small M block, optional within-program low-rank (comodules) | native unit; interpretable as "regulatory matrix of EMT/P53"; ties to architecture layer (§5) | moderate |
| **Universe-wise** | all genes × all arms | RRR (low-rank) + group LASSO + interaction priors | structural roll-up, NMF comodules, borrows strength | heaviest; identifiability bites hardest; circularity risk highest |

**Recommendation: program-wise as the primary lane** (matches the Hallmark unit,
tractable p ≫ n, directly interpretable), **gene-focused** for the driver/orphan deep
dives, **universe-wise** only as the low-rank structural overlay (it is essentially
our NMF universe layer, `MODELING_FRAMEWORK.md` §5.6, applied to M rather than to
pressure).

---

## 6. The causal/perturbation lane (the instrument)

```
Stage 1:   x_m  =  γ · CN_locus(m)  +  C·δ  +  u        (miRNA-locus CN → arm abundance)
Stage 2:   y_g  =  −M_{g,m}·x̂_m     +  C·β  +  ε        (2SLS / MR estimate of the causal weight)
```

`CN_locus(m)` is a *natural perturbation* of the miRNA dose that is (largely) not
caused by the target's expression — so it breaks the confounding loop that
regression alone cannot. Prior art: CNV↔expression causal inference; miRNA
Mendelian-randomization. **We have the instrument already built** (`mirna_locus_cnv.py`,
`MIRNA_CNV_DOSAGE_REPORT.md`) and have only ever used CN as a *confounder on the
target side*. Caveats to honor: weak-instrument (many miRNA loci have little CN
variance — report first-stage F), and locus CN can hit neighboring genes
(pleiotropy) — restrict to focal segments. This is the strongest answer to "is the
weight causal," and almost no expression-network paper carries an internal
instrument.

---

## 7. What we have over prior work — summary

Relative to GenMiR++, the LASSO miRNA-target networks (Lu 2011), Theia/PIMiM
comodule methods, and the matrix-factorization miRNA-disease line:

1. **Prior at magnitude, not just support** (§1) — quality-weighted, assay-directness-graded,
   non-circular, breast-audited evidence enters as prior mean + adaptive penalty. The
   learned model **nests** our fixed-M model.
2. **Confounder-correct by construction** (§2) — purity/proliferation/target-CN/batch
   in the design; the #1 failure mode of these networks, controlled.
3. **Seed-family identifiability as a constraint** (§3), not an ignored collinearity.
4. **Redundancy vs cooperativity modeled with priced interaction terms** (§4).
5. **A built-in causal instrument** (miRNA-locus CN, §6).
6. **Program-wise lane on a curated program atlas** (§5) — every M block reads as a
   named Hallmark's regulatory matrix.
7. **Out-of-distribution validation at the protein layer** (CPTAC L1b,
   `MODELING_FRAMEWORK.md` §7) and in independent cohorts — most methods stop at mRNA.

---

## 8. The non-negotiable guardrail: circularity

Our central principle (`MODELING_FRAMEWORK.md` §0.3 #2) is that the prior never
touches the expression it is tested against. **If we learn M from the same (X, Y)
we later test coupling on, we destroy that.** The resolution is standard:

- **Sequence/curation prior fixes the support** (which entries can be nonzero) and
  the prior mean; **expression estimates the magnitudes**.
- **Validate out-of-sample**: patient-level CV folds (we already do 5-fold,
  `held_out_tuning`), and on an **orthogonal layer** the fit never saw — CPTAC
  protein (L1b) and/or the CN-instrument estimate.
- Keep the **learned M and the curated M side by side** and report where they
  disagree (the PTEN→miR-19 vs miR-17 kind of disagreement) as a *finding*, never
  silently let the data overwrite the prior.

---

## 9. Suggested build order

1. **Gene-focused, conditioned, group-LASSO** on a hub panel (PTEN, GATA3, ESR1,
   ZEB1, IRF1) — adaptive penalty from `edge_w`, `C` unpenalized, seed-family groups.
   Cheapest; directly comparable to the existing per-edge coupling and the
   `P_others` experiment (E4/MH-78). Decision gate: does the joint estimate recover
   the curated drivers and beat raw-abundance coupling out-of-fold?
2. **Add cooperative product terms** (§4) on that panel, support primed by TargetScan
   site proximity. Gate: do any survive the penalty + improve held-out fit?
3. **Program-wise** with within-program low-rank (comodules) on 3–5 Hallmarks
   (EMT, P53, G2M, IFNγ, HYPOXIA). Gate: comodules interpretable; couples on held-out.
4. **CN-instrument (2SLS)** on the same hub panel — report first-stage F, compare
   causal weight sign/magnitude to the conditioned estimate.
5. **Universe-wise RRR** only if 1–4 pay off — as the structural/NMF overlay.

Each gate mirrors `PRESSURE_FUTURE_OPTIONS.md` discipline: promote a lane only on a
pre-named held-out / orthogonal-layer win, never on in-sample fit.

---

## Related

| Path | Content |
|------|---------|
| `MODELING_FRAMEWORK.md` §2–6 | the fixed-M pressure + coupling this generalizes |
| `EDGE_QUESTION_TAXONOMY.md` | a single entry of M = the edge question (coupling-invariance lemma) |
| `GENE_QUESTION_TAXONOMY.md` | a row-sum of M = the gene-aggregate question ("construction is consequential") |
| `PRESSURE_FUTURE_OPTIONS.md` | upstream evidence-scoring lanes that feed the prior (§1) |
| `mirna_locus_cnv.py` / `MIRNA_CNV_DOSAGE_REPORT.md` | the instrument (§6) |
