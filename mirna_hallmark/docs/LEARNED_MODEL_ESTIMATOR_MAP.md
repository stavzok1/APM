# Learned model — ESTIMATOR MAP (what model for what, and why over alternatives)

Single source of truth for **which estimator does which job, and what it was chosen over** in
`mirna_hallmark/learned/`. For the **actual formulas/objectives** of each estimator see the companion
`LEARNED_MODEL_METHODS.md`; for the **conceptual reasoning behind each** (why this construction, what breaks
otherwise) see `LEARNED_MODEL_RATIONALE.md`. Companion to `LEARNED_MODEL_DESIGN_RESPONSE.md` (the design decision-log),
`LEARNED_MODEL_BUILD_PLAN.md` (what's built), `LEARNED_MODEL_WHATS_NEXT.md` (forward list). Last updated
2026-07-06 (pooled-HE inclusion migration; loss lens added; attribution un-stubbed).

**Inclusion policy (map-wide, 2026-07-06):** every estimator below draws a gene's candidate regulators from
**POOLED-HE** = miRTarBase-HE ∪ TarBase-v9 low-throughput functional (reporter/western/proteomics, non-weak) —
`ledger.pooled_he_edges()`, ~5,940 edges / 1,571 Hallmark genes. This replaced miRTarBase-HE-only
(`he_only=True` now returns pooled-HE); gate-stability verified on anchors (rho_model Δ≤0.003). Candidate arms
for the biochemical potential are scanMiR-modelled with a **TargetScan-context++ fallback** for arms scanMiR
lacks (recovers e.g. miR-137). Scope is the Hallmark universe by design; `he_only=False` → full hallmark set.

## The organizing principle: match the estimator to the estimand

> **UPDATE (MH-82/83/84):** the map below is the *original lasso-centric* production layout. It has since
> been consolidated into a **single dense learned-τ² Bayesian fit** that serves coupling + attribution +
> identifiability from one posterior (better/equal, genome-wide, CPTAC-validated), with an inclusion-mode
> variant for discovery. See **`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`** for the one-fit→4-uses architecture
> and `RATIONALE.md` §2a–§2g for why. The lasso rows below remain valid as the *shipped* pipeline + baseline.

> **UPDATE (MH-92, 2026-07-10) — observation likelihood:** the dense posterior's likelihood is now selectable
> **Gaussian (default) or Student-t** via `spike_slab._gibbs_posterior(nu=)` (scale-mixture of normals, per-obs
> λ_i, stays Gibbs-conjugate; `nu=None` = Gaussian, bit-identical). Chosen because the dense mRNA residuals are
> heavy-tailed (t-df≈7; NB rejected — no low-expr count gradient). It **robustifies** (down-weights amplified-subset /
> near-floor outlier tumours) without changing coupling/attribution verdicts at n≈1000 (OOF Δρ +0.001 n.s.). Same
> option on the CN channel term (`t['nu']`) and an optional bootstrap SE on the exclusion gate (`exclusion(n_boot=)`,
> confirms Sobel adequate). METHODS §3b (twins); VALIDATION §1; CHANNEL_FUSION §AB.

There is **one model** — the gene-focused non-negative adaptive lasso — and around it a small set of
purpose-built estimators, each matched to a specific question. The recurring mistake is to assume one
estimator must serve every question. It must not, because the questions have **different robustness**:

- **PREDICTION / SELECTION questions** (does coupling exist? does the model beat abundance? which new edges
  enter?) read only the *aggregate* X·M or an inclusion decision. Under collinearity the aggregate is
  **stable even when the individual coefficients are not** (classic result: collinearity inflates
  coefficient variance, not fitted-value variance). Verified: ρ_lasso ≈ ρ_baggedNNLS within 0.03, bootstrap
  SD ~0.01. → these use **the adaptive lasso** (the model, prior-bearing, selective).
- **COEFFICIENT-READ questions** (how much does *this* edge contribute? how does its weight shift across
  states?) read the individual M(m,g), which **is** the unstable object under collinearity (single-fit
  bootstrap corr 0.03). → these use **bagged NNLS on the fixed seed-family support** (stability by design).

So two estimators, by design, not by accident. Using bagged-NNLS for the gate would certify the wrong
object; using single-lasso coefficients for attribution would report noise.

## Map (production)

| # | Job / estimand | Estimator in use | Module | Chosen over |
|---|----------------|------------------|--------|-------------|
| 1 | **The model**: per-gene repression weights + aggregate pressure X·M | Gene-focused **non-negative adaptive lasso** (evidence-weighted penalty; C partialled out) | `regression.fit_gene` | OLS (no sign, unstable), plain lasso (no prior), ridge (no selection), curated fixed-M heuristic (the *nested* Bar-2 baseline it must beat) |
| 2 | **Validation gate**: is there real coupling; does learning beat raw abundance / match curated? | Out-of-fold **partial-Spearman** of X·M vs held-out mRNA \| C | `mvp.oof_gate` | raw-abundance mean (Bar-1 null), curated fixed-M (Bar-2 nested limit) — both computed alongside as the comparators |
| 3 | **Discovery**: which *new* (orphan) edges exist? | adaptive lasso **+ permutation-FDR null + scanMiR biochemical filter + deconv composition tag** | `discovery.py` | NNLS/fixed-support (structurally *cannot* select), naive correlation (no FDR / no confounders), lasso-alone (no biochemical or composition control) |
| 4 | **Attribution — budget/force**: per-edge/family budget share & rank, ΔM (wiring), NAT→T decomposition | **Bagged z-scored NNLS on the fixed pooled-HE seed-family support** (+ per-state variance floor). Budget contribution apportions the family weight to arms by **linear-RPM share of the family pool** (`Σ_arms = M_fam·X_fam`, no family-size inflation — §4a) | `states.canonical_M` / `_bagged_nnls`, `states.budget_shift` | single adaptive-lasso coefficient (unstable, corr 0.03), **level-scale NNLS** (ill-conditioned → recovers noise arms), summed M·log-abundance broadcast (inflates multi-arm families ~×size), hierarchical posterior mean (pooling mis-specified) |
| 4b | **Attribution — identity/ownership**: who *owns* the coupling (abundance-removed), collinearity-fair | **Shapley/LMG** decomposition of the fixed-M aggregate's R² across families | `attribution.shapley_identity`, `identity_vs_magnitude` | budget share (that's *force*, abundance-included — different question); raw per-family ρ (double-counts shared variance under collinearity). Discipline: identity ⊥ coupling (says *who*, not *whether*) |
| 5 | **Resolution**: which family member is the DRIVER? | **Conditioned partial anti-corr** — residualise the member arm & target on (other family members + C), driver = member whose anti-corr survives (<0); collinear members auto-fail (no residual variance) | `families.resolution_report` | old **abundance-divergence heuristic** (nominated the most-distinct arm, never looked at the target — replaced §4b); per-arm lasso (keeps one member arbitrarily). Estimand = family→gene; a member is a driver only if its *independent* variation tracks the target |
| 6 | **Coupling readout** (per-edge, per-state; the card's coupling_hly/nat/tum) | **Partial-Spearman** of arm vs gene \| C (estimator-free) | `card.py`, `states._aggregate_coupling` | any fitted M (unnecessary — a single edge needs no model; direct partial-ρ is assumption-light and robust) |
| 7 | **Causality**: is the edge causal, not reverse / confounded? | **CN-locus instrument** (miRNA-locus copy number → abundance → target); pleiotropy + cluster checks | `instrument.py` | partial correlation alone (cannot exclude reverse causation / shared upstream driver) |
| 8 | **Uncertainty / identifiability object**: how sure is each weight | **Hierarchical Bayes** (Gibbs) → posterior mean + **sd** (\|β/sd\|>2 = identified). Cross-gene μ_m pooling BUILT but **does not improve coupling** on p53/cell-cycle (gate 2026-07-05: Δρ −0.002 full-n, **−0.045 at n=150** — pooling *hurts* at small-n = mis-specified: effects gene-specific + genes barely share regulators, most miRNAs hit 1 program gene so μ_m is trivial). → used as the **uncertainty object only**, coupling stays per-gene | bootstrap SD (used for the point estimate); the pooling's intended weak-gene shrinkage empirically didn't materialize here |
| 9 | **The prior** (adaptive penalty): which edges, how strong | **Evidence ledger fusion** (PMID-deduped, method-centric, transfection-calibrated) for *existence*; **scanMiR K_D** for *magnitude* | `evidence/ledger.py`, `kd.py` | hand-weighted `evidence_score` (baseline, un-deduped), raw PMID count (no dedup / no directionality / no method weighting) |
| 10 | **Identification block C** (what's partialled out) | **CPE purity + target-gene CN + malignant-compartment proliferation** (+ CIBERSORTx non-malignant composition for cell-intrinsic). **NB `mal_prolif` is materially INCOMPLETE** — a residual-proliferation axis (\|ρ\|≤0.83 with the target after C) remains, yet a fuller *global* control OVER-controls the majority of genes ⇒ left as-is; the per-gene verdict is row 13 (MH-100) | `data.assemble_gene` | global PCs (over-correct — *are* the programs miRNAs regulate), TF-regulon proxy (over-controls −60–80%), naive plate dummies (strip biological heterogeneity), HRD (dropped, neutral), **a fuller global proliferation control** (over-controls 50% of genes → rejected, row 13) |
| 13 | **Proliferation confound-vs-mechanism** (per gene: is C's residual-proliferation crack a confound to control or the mechanism to keep) | `R = resid(broad cell-cycle metagene \| C)`, target-excluded; non-circular **OOF 2×2** train {C, C+R} × eval {C-space, C+R-space} → {**confound** (control helps both) / **over_control** (hurts both = mechanism) / **fragile** (flips with eval space) / neutral}; `rho_resid` = residual-prolif still in target. **Per-EDGE layer:** gene class × edge Δβ direction → {unmasked_real / deconfound / mechanism / prolif_robust / ambiguous}. **Arm-SOURCE (dose-delivery, §8):** per-arm abund-share × corr(arm, R) → which arm carries the family's proliferation confound (median concentration 0.85; systematically miR-93-5p). **Mechanism (`confound_mechanism`): HOST-GENE transcriptional, NOT CN** — miR-93/106b/25 intronic in MCM7; host mediates 65–97% vs locus CN 8–10% (axis BP) | `prolif_verdict.py` (`gene_verdict`/`run`; `run_genomewide` = 8-worker parallel 1507-gene; `arm_prolif_source`/`arm_source_scan`/`confound_mechanism`) | a GLOBAL fuller proliferation control (over-controls 50%, mean Δ_Cspace −0.032 → rejected; C row 10 unchanged); AL2 nascent-transcription latent (data-blocked, no intronic reads); reading Δβ alone (ambiguous — de-confound, over-control, and un-masking all shrink β) |

## Map (parked / scaffolded — deliberately NOT in production)

| Job | Estimator | Module | Why parked |
|-----|-----------|--------|-----------|
| Dose–response / retire D(m) | Saturating **occupancy link** a/(a+K) + threshold gauge | `occupancy.py` (STUB) | tested; the **saturation substrate failed** the held-out gate → kept linear a·M. Gate hard before reviving. |
| Cooperativity / ceRNA shared pool | Site-primed **product terms**; cohort shared free-AGO pool | `cooperativity.py` (STUB) | sits *on* the failed saturation substrate (Decision G) → not built until occupancy earns its place |
| Selection via explicit inclusion prior | **Spike-and-slab** `z~Bernoulli(π)` + half-normal slab, SSVS Gibbs (Design §D/E) | `spike_slab.py` (**BUILT + benchmarked**) | tested head-to-head vs the lasso on the SAME gate: recovers curated drivers with a real PIP but **loses on coupling** — genome-wide FDR-sig 35% vs lasso 45%, mean ρ −0.078 vs −0.125 (VALIDATION §1; rationale §2a). Kept as the documented Bayesian alternative; lasso ships. |

## The LOSS lens (built 2026-07-06 — the expression-free "designed & lost" axis)

Distinct from the coupling model above: it answers "who is *designed* to repress g" and "was a real repressor
*lost/gained*", which the coupling-gated estimators (M / budget / Shapley) structurally miss because a silenced
arm has no variance and no abundance → they read 0.

| # | Job / estimand | Estimator | Module | Chosen over |
|---|----------------|-----------|--------|-------------|
| 11 | **Structural identity**: who is DESIGNED to specifically repress g (expression-free) | `structural = potential · specificity · confidence`, family-pooled. potential = scanMiR/TargetScan biochem × (1+log1p evidence), direction-filtered; specificity = family-pooled evidence-mass fraction on g; confidence = min(1, pooled ev_depth/10) as a **multiplier** | `structural_identity.py` | M/budget/Shapley (all coupling-gated → 0 for a silenced arm); a biochemical-specificity floor (there is NO biochemical specialist — every arm hits ~700 targets ~uniformly, so specificity is purely an evidence construct); reading one representative arm's evidence (understudied-member artifact — pool the whole seed family) |
| 12 | **Loss/gain mechanism**: is a structural specialist epigenetically silenced (or de-repressed) in tumour | **Tumour−normal promoter Δβ** over CpG probes DIRECTLY overlapping the miRNA promoter (full HM450 manifest; no cCRE hop). Bidirectional: hyper→repression LOST, hypo→repression GAINED | `methylation.py` (`lost_specialist` integrates it) | cCRE-aggregated methylation (indirection; the gene-centric probe annotation misses miRNA loci); expression alone (can't distinguish epigenetic silencing from CN loss / transcriptional) |

**Program-level view (built, 2026-07-05):** since cross-gene pooling doesn't help (#8), the program view is a
**descriptive roll-up** of the canonical per-gene M — `program_network.py`: family×gene weight matrix +
per-family engagement (n targets, weight, focus=max/sum broad-vs-narrow) + per-gene incoming + network stats
(hub families co-targeting ≥2 genes). Aggregates AFTER fitting → no pooling bias. On EMT: 60 hub families
hold 76% of program weight; miR-29 broadest (21 targets), miR-200/17~92/30/103/21 the canonical hubs.

## The three "why over alternatives" that matter most

1. **Adaptive lasso over the curated fixed-M heuristic (#1 vs its Bar-2 baseline).** The heuristic *is* the
   lasso with the prior frozen (the nesting limit). The learned aggregate beats raw abundance OOF *and*
   ≥ the frozen-prior aggregate with fewer arms — the gate (#2) is exactly this test (40-gene scan
   2026-07-05: beats abundance **38/40**, ≥ curated **37/40**, mean sparsity 0.35). On well-studied hubs
   (PTEN/ESR1) curated ≈ learned, but on **non-hub genes the frozen evidence weight gives ~0 or WRONG-SIGN
   coupling** (MET rho_curated +0.23 → learned −0.17; IGF1/ERBB2/E2F1/MYC similar) because the evidence
   weight is dominated by the most-*studied* edges, not the most-*predictive* — e.g. MET keeps miR-34/449
   (evid rank 1) + data-selects miR-181a-3p (evid rank 21) and drops the well-studied-but-flat miR-1/206,
   miR-23, miR-199. So learning materially rescues coupling where literature prominence ≠ in-cohort strength;
   it is not merely sparsity. (Earlier "frozen weight nearly suffices" was a hub-panel artifact.)

2. **Bagged NNLS over the single lasso coefficient for attribution (#4).** PTEN has ~50 collinear weak
   regulators over ~1000 samples → any single sparse fit's per-edge weight is under-determined (bootstrap
   corr 0.03). Bagging over 40 resamples drives estimate variance ~1/B → corr 0.99 and recovers the curated
   drivers (miR-21/103a/181/182/96). z-scoring (not level) because level-NNLS is ill-conditioned and puts
   weight on noise arms; a per-state **variance floor** (std<0.1 → 0) because z-scoring an undetectable arm
   in small states (NAT n=104) divides by ~0 and manufactures weight (miR-105-3p M_NAT 0.157→0).

3. **Adaptive lasso, NOT bagged NNLS, for discovery (#3).** Discovery is a *selection* problem — which
   orphan edges cross an FDR-controlled null. Bagged NNLS has a *fixed* support and no penalty; it has no
   mechanism to admit or reject an edge. So bagging isn't "redundant here" — it's the wrong tool for the
   estimand. It is the *post-selection* attribution readout, downstream of the lasso's selection.

## One-line summary

> **Adaptive lasso = the model** (prediction, selection, validation, prior). **Bagged family-NNLS = the
> stable readout** of that model's coefficients (attribution). **Partial-Spearman = the assumption-light
> coupling measure** (single edges). **CN instrument = causality. Hierarchical Bayes = uncertainty.**
> Family collapse is the identifiability unit under all of them; the C block is the identification.
