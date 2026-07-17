# Learned model — METHODS (formula-level spec of every estimator)

The formal math for each estimator in `mirna_hallmark/learned/`. Plain-text notation (no LaTeX — renders in a
terminal). This doc is the **spec** — notation, objective, formula, one line of "what it computes" — kept
clean of prose. **Which estimator for which job, and why over the alternatives → §18 (in this doc).**
Companions: `LEARNED_MODEL_RATIONALE.md` (**the *why* behind every §here — section numbers
mirror this doc**), `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b (the converged model in one line), and
`LEARNED_MODEL_VALIDATION.md` (does it work). This is the "what is the actual objective/formula" reference.
(The precursor **pressure** model's math is separate, in `METHODS.md`.)

## §0 Notation & identification
Per gene g, aligned on tumour participants (PAM50 Normal-like dropped):
- `Y` = target mRNA, log2(TPM+1), vector over participants.
- `X` = regulator-arm abundance, log2(RPM+1), participant × arm (undetected arm → 0).
- `C` = confounder block (participant × conf), see §1.
- `w` = adaptive-penalty prior weight per arm, see §2.
- **Candidate regulators (inclusion, migrated 2026-07-06):** POOLED-HE = miRTarBase-HE ∪ TarBase-v9 low-throughput
  functional (reporter/western/proteomics, non-weak) — `ledger.pooled_he_edges()`, ~5,940 edges / 1,571 Hallmark
  genes; `assemble_gene(he_only=True)` returns this. Replaced miRTarBase-HE-only; gate-stable on anchors (Δρ≤0.003).
Residualisation (matrix-view identification, Design §2/§B): for a vector v and matrix C,
`resid(v | C) = v − C·(C⁺ v)` (OLS fit removed; intercept included). Every coupling/weight below is estimated
in the residual space `(I − P_C)`, i.e. C is partialled out of BOTH the predictor and the target.

## §1 Confounder block C  (`data.assemble_gene`)
Core (always): `C = [CPE, target_cn, mal_prolif]`.
- `CPE` = consensus tumour purity (per participant).
- `target_cn` = gene g's own ASCAT3 copy number (`_target_cn`) — the gene's genomic dosage.
- `mal_prolif` (`_malignant_prolif`) = malignant-compartment proliferation:
  `score = mean_g'∈(E2F∪G2M) zrow(Y_g')`  then  `mal_prolif = resid(score | CancerEpithelial_fraction)`
  (removes the composition-driven part, leaving per-malignant-cell proliferation).
Optional add-ons (joined as extra C columns):
- `deconv` → the 8 non-malignant CIBERSORTx fractions `[CAFs,T,Myeloid,B,Endothelial,PVL,NormalEpi,Plasmablasts]`
  (Cancer-Epithelial excluded — conditioning on the target's own compartment over-controls).
- `mycaf` (`_mycaf`) = `resid( mean z(ACTA2,TAGLN,POSTN,FAP,COL11A1,THBS2,MYH11,CNN1) | CAF_fraction )` — myCAF-vs-iCAF axis.
- `mal_emt` (`_mal_emt`) = `resid( mean z(VIM,FN1,CDH2,SNAI2,ZEB1,ZEB2,TWIST1,SPARC) | non-malignant fractions )` — within-cancer EMT state.
- HRD off; naive batch off (post-fit conditioning only, never in-fit).

**DEFAULTS (which C is actually used):** the **core `[CPE, target_cn, mal_prolif]` is always on.** The
non-malignant **CIBERSORTx fractions (`deconv`) are OPT-IN**: OFF for the gate/coupling (`assemble_gene`
`deconv=False`), but **forced ON for the canonical ATTRIBUTION** (`canonical_M`/`budget_shift`/`decompose` call
`deconv=True`) → those weights are *cell-intrinsic*. The **metagene-z composition** (`_state_metagene_cov`, a
per-program z-metagene) is only the **cross-state fallback** for GTEx (no CIBERSORTx there); NAT/tumour use the
CIBERSORTx fractions. So: tumour attribution = fractions; gate/coupling = core only unless deconv requested;
healthy = metagene-z until the pending GTEx CIBERSORTx run lands.
**Post-fit conditioning:** `batch` (TCGA plate dummies) is **never in the lasso** — an in-fit ComBat-style
correction removes biological heterogeneity along with the batch — so batch is appended to C **only for the
partial-Spearman ρ significance computation** (`oof_gate`), conditioning the *test*, not the *model*.
**Within-cancer vs non-malignant:** `mal_prolif` residualises the proliferation metagene on the **Cancer-
Epithelial** fraction → per-malignant-cell proliferation; the `deconv` block adds the **non-malignant**
fractions (CAFs/immune/…) to strip stromal/immune composition. **Cancer-Epithelial is deliberately excluded
from the `deconv` block** — conditioning on the compartment the target is expressed in over-controls the signal.

## §2 The prior w  (`evidence/ledger.py`, `evidence/methods.py`)
Evidence enters as MAGNITUDE, non-circularly:
- Per (edge, assay class) count **distinct PMIDs**; fused edge weight
  `w_m = Σ_class  w_class · log1p( distinct_PMID_count(m, class) )`.
- `w_class` (assay directness) either hand-set (`CLASS_WEIGHT`: reporter 3 > western/proteomics 2.5 > qpcr 1.5
  > degradome/chimeric 1 > ago_clip 0.5 > other 0) OR **learned** (`CLASS_WEIGHT_MRNA`): non-negative LS of
  measured transfection repression (McGeary GSE140217) on per-class log-PMID, within-family-centred.
- `w = clip(w_m, ≥1e-3)`. Alternatives: `scanmir` (K_D affinity), `fused` (existence×magnitude).
- **DEFAULT:** the `assemble_gene`/`oof_gate` signature default is `evidence_score` (hand-weighted), but the
  **canonical/reported fits force `ledger`** (the weighted, PMID-deduped, method-centric fusion) — so we use
  **weighted evidence, not binary inclusion**. Empirically the source is **non-critical** (hub panel, family
  gate): mean rho_model ledger −0.345 ≈ evidence −0.340 ≈ fused −0.342 ≈ scanmir −0.331; stability 0.86/0.86/
  0.85/0.84; ledger passes 5/5 gate, scanmir/fused 4/5 (scanmir breaks GATA3, fused breaks PTEN vs-curated).
  ⇒ ledger is the marginal best and is canonical; the prior sets ordering, the data sets magnitude (§3).

## §3 THE MODEL — gene-focused non-negative adaptive lasso  (`regression.fit_gene`)
> **SUPERSEDED FOR COUPLING (MH-82, genome-wide):** the **dense learned-τ² non-negative ridge**
> (`spike_slab._gibbs_posterior`, π≡1) BEATS this lasso on OOF coupling (mean ρ −0.168 vs −0.152, p=9e-16)
> and serves attribution (`mean±sd`) + identifiability (`|z|>2`) from the *same* fit. The lasso remains the
> *shipped-pipeline* estimator and the Bar-2 baseline the learned-τ² is validated against. Full architecture
> (one fit → 4 uses) in **`LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`**; estimator sweep in `RATIONALE.md` §2a–§2g.

Objective (one gene), on C-residualised ỹ=resid(Y|C), X̃=resid(X|C):

    M_hat = argmin_{M ≥ 0}  || ỹ − X̃·M ||²  +  α · Σ_m |M_m| / w_m

- Non-negativity encodes "miRNAs repress" (Design §Decision C).
- **Adaptive** penalty λ_m = α/w_m: strong-evidence edges resist shrinkage, weak edges shrink out. The prior is
  an INVERSE PENALTY (resistance to being zeroed), NOT a location — it sets ordering/selection, not magnitude.
- Implementation (Zou reparametrisation): scale column m by w_m → plain non-negative lasso (uniform α) → unscale
  `M_m = b_m·w_m`. Default α=0.005 (fixed; **non-critical** — global grid 0.001–0.02 gives flat mean rho_model
  −0.340, gate 8/8; α is a **sparsity knob** [nonzero 9→7], not a coupling knob; only 0.05 over-shrinks). α
  touches ONLY the lasso — the canonical attribution (§5) is **unpenalized** bagged NNLS.
- Aggregate repression pressure: `ρ_agg = X·M` (higher ⇒ more repression ⇒ lower Y).
- This is the MODEL used for coupling/selection/validation (§6) and discovery. Contribution ≈ a·M (unsaturated
  linear occupancy; the saturating a/(a+K) link was tried and failed the gate).

## §3b OBSERVATION LIKELIHOOD — Gaussian vs Student-t  (`spike_slab._gibbs_posterior(nu=)`, MH-92, 2026-07-10)
The dense §6b posterior fits `r = X·β + ε`. Its default likelihood is **Gaussian** `ε_i ~ N(0, σ²)` (voom-style, on
log-CPM residuals). A residual diagnostic (85 genes × 1041 patients) found the residuals **heavy-tailed** — pooled
excess-kurtosis 3.5 (median 1.2 across well-expressed genes), Student-t MLE **df ≈ 7**, tails P(|z|>4)=33× Gaussian —
from two mechanisms: amplified-subset tumours (+skew, ERBB2/ERBB3) and a few near-floor genes (−skew). **NB is rejected**
(no low-expression count gradient: Spearman(expr, kurtosis)=+0.19, *wrong* sign; the near-floor genes are a detectability-
floor issue, not a global count structure). So the honest upgrade is **Student-t** `ε_i ~ t_ν(0, σ²)`, implemented as a
**scale-mixture of normals** (stays Gibbs-conjugate, no HMC):

    ε_i ~ N(0, σ²/λ_i) ,   λ_i ~ Gamma(ν/2, ν/2)          # marginal ε_i ~ t_ν
    λ_i | ε_i,σ² ~ Gamma( (ν+1)/2 , (ν + ε_i²/σ²)/2 )      # extra per-obs Gibbs step

The per-coordinate conditional gains λ-weights: `A_m = Σ_i λ_i x_{im}²/σ²`, `B_m = Σ_i λ_i x_{im} r_{m,i}/σ²`, and
`σ² ~ InvGamma(a₀+n/2, b₀+½Σ_i λ_i ε_i²)`. A large-residual patient (amplification outlier, floor spike) gets small
`λ_i` ⇒ down-weighted. **`nu=None` ⇒ Gaussian, bit-identical to §6b** (no λ draws). **Verdict:** correctness/robustness
upgrade, **not** a performance lever at n≈1000 — nests Gaussian (Δβ≈0.001), robust 1.76× under 5% gross contamination,
but OOF coupling Δρ=+0.001 (n.s.) and attribution unchanged (corr 0.999, 5/5 drivers). The same scale-mixture applies
to the **CN channel term** (a t on the single `π̂_s` pseudo-obs, `t['nu']`; `s²_π` is estimated ⇒ t is honest + guards
collider-outlier segments) — verified to discount an outlier segment (β 0.535→0.508) while leaving a consistent one intact.

## §4 Seed-family collapse  (`families.collapse_by_family`, Design §F)
Same-seed arms are near-collinear ⇒ the identified estimand is family→gene (arm = nomination).
- Pooled predictor: `X_fam = log2( 1 + Σ_{m∈family} (2^{x_m} − 1) )`  (true RPM pool, re-logged).
- Pooled prior: `w_fam = max_{m∈family} w_m`.

## §5 CANONICAL ATTRIBUTION — bagged family NNLS  (`states._bagged_nnls`, `canonical_M`)

> ✅ **SETTLED 2026-07-17 (MH-138) — this section and §6b were never in conflict; ONE WORD was doing TWO JOBS.**
> **MAGNITUDE** (who delivers how much) = the dense Gibbs posterior mean `beta` ⇒ **§6b is right**.
> **IDENTITY** (who owns the gene) = `attribution.shapley_identity` (LMG on R²), which requires **fixed
> weights** — and those come from `canonical_M`, i.e. **the bagged NNLS specified below** ⇒ **§5 is right, for
> that job.** The code was already coherent: `attribution.identity_vs_magnitude:131` and
> `bayes_shapley_identity:91` both call `ST.canonical_M`, the latter exactly as `RATIONALE §2e` prescribed —
> **NNLS fixes the support and signs; the Gibbs draws supply the width.**
> **⭐ WHY identity cannot use the Gibbs mean (§2e's mechanism, now MEASURED on `readouts_edges.tsv`, n=5,117):**
> the half-normal slab `N⁺(0,ν²)` has a **strictly positive mean**, so an un-informed family **cannot be zeroed**
> — it relaxes to the prior. **`beta == 0` in 0/5,117 (0.00%); 100% are > 0**, where NNLS returns exact zeros.
> Feed that into a Shapley and you credit families the model cannot distinguish from noise. **Pooled, 73.0% of
> all β-mass sits on unidentified edges (|z|≤2)** (per-gene median 55.5%; 42.7% of genes at 100%, of which 452
> are single-family where β ≡ uniform anyway; `n_fam≥3` median 32.0%).
> Also settled: **bagged NNLS is not retired** — it is correct for the **cross-cohort gauge** (Gibbs's
> heavy-tailed posterior SDs break the errors-in-variables correction: a=4.1 vs a split-half truth of 1.0) ⇒
> ***bagged NNLS for the GAUGE, Gibbs for the MODEL***; on the model itself Gibbs reproduces better
> (split-half ρ **0.822** vs **0.729**).
> ⚠ The premise below (single-*lasso* coefficients are unstable) is pre-convergence framing — the lasso is a
> baseline now; the collinearity instability it describes is real regardless. See §18 and MH-138.

The stable readout of per-edge weight (single-lasso coefficients are unstable under collinearity, corr 0.03).
On C-residualised, z-scored family predictors:

    yr = − resid(Y | C)                                    # sign so repression → positive weight
    Xz = zscore_cols(X_fam)                                # (X_fam − mean)/sd
    Xz[:, sd < 0.1] = 0                                    # VARIANCE FLOOR: undetectable-in-state arm → no weight
    for b in 1..B:  idx_b = bootstrap(n);  c_b = NNLS(Xz[idx_b], yr[idx_b])
    M_fam = (1/B) Σ_b c_b                                  # BAGGED (B=40); variance ↓ ~1/B

- Reproducible (corr ~0.99 across seeds) and recovers curated drivers (PTEN miR-21/103a/181/182/96). z-scored
  because level-scale NNLS is ill-conditioned (recovers noise arms).
- **arm_level broadcast:** `M_arm = M_fam( family(arm) )` — family weight shared to member arms (per-arm split
  = nomination). Used everywhere COEFFICIENTS are read: budget, ΔM/wiring, decompose, card.
- ΔM / cross-state wiring: fit `M_fam` separately per state (tumour "01" / NAT "11" / GTEx) on the SAME fixed
  family support → `dM = M_state1 − M_state2` (fixed support ⇒ comparable estimand).

### §5a Budget contribution (family-size-correct)  (`states.budget_shift`)
Realized within-gene budget, apportioned so a multi-arm family isn't over-counted:

    contribution(arm, state) = M_fam · X_fam(state) · [ 2^{a_arm(state)} / Σ_{m∈fam} 2^{a_m(state)} ]
    share(arm)  = contribution / Σ_arms contribution ;   rank = descending rank of contribution

so Σ_{arm∈fam} contribution = M_fam · X_fam (no family-size inflation). A singleton's X_fam = its own abundance,
so contribution = M·a (unchanged).

## §6 VALIDATION GATE — out-of-fold aggregate coupling  (`mvp.oof_gate`, `mvp.gate_fdr`)
5-fold over participants; fit M on train, score held-out:

    oof_model[te] = X[te] · M_train ;   oof_abund[te] = mean_m X[te]         # per participant
    rho_model  = spearman( resid(oof_model | C), resid(Y | C) )              # expect < 0
    rho_abund  = spearman( resid(oof_abund | C), resid(Y | C) )              # Bar-1: equal-weight inclusion
    rho_curated= spearman( resid(X·w_prior | C), resid(Y | C) )              # Bar-2: FROZEN evidence weights (no fit)

Gate passes if `rho_model < rho_abund` (learning beats abundance) and `rho_model ≤ rho_curated` (≥ frozen prior).
Genome-wide FDR (`gate_fdr`): per gene one-sided partial-t p for rho_model<0 on df≈n−8, then BH and
Benjamini-Yekutieli across genes (BY for shared-regulator dependence). Only meaningful for ≥2-regulator genes
(singletons: X·M ≡ abundance up to scaling).

## §7 PER-EDGE COUPLING — family-aware permutation partial-Spearman  (`coupling.edge_coupling` → `coupling_inference.permutation_pvalues`)
Per gene, predictor rows = arms, target = gene mRNA (repeated), covariates = C, families = seed family:

    obs_m = Σ_s  unit(resid(rank(x_m) | rank(C)))_s · unit(resid(rank(Y) | rank(C)))_s   # partial Spearman
    null:  same permutation π of the predictor-residual sample order applied to EVERY row (preserves family dep.)
    p_perm = (1 + #{null ≤ obs}) / (n_perm+1) ;   p_z = Φ( (obs − μ_null)/σ_null )       # smooth tail from null moments
    q_bh = BH(p_z) ;  q_by = BenjaminiYekutieli(p_z) ;  q_family = Simes-per-family BH(p_z)   # §14

Called PER GENE (C carries per-gene target_cn). Genome-wide: second-level BY across all edges.

## §8 CAUSALITY — CN-locus instrument  (`instrument.edge_instrument`, `run_clean`)

> ⛔ **THERE IS NO LIVE "CAUSALITY" JOB — BOTH CN INSTRUMENTS ARE RETRACTED (banner added 2026-07-17).**
> This section still specs the machinery, and the machinery still runs — but **its causal interpretation is dead**,
> and the retraction is in the production code itself (`learned/instrument.py:1128–1148`).
> **(1) `pi_causal` is NOT an IV** (MH-124r/126, 2026-07-13). It is `γ_s·b_fam`, a **product-of-coefficients
> MEDIATION** estimator in which `b` is the **OBSERVATIONAL OLS slope** — the instrument-*orthogonal*
> (endogenous) variation, carrying every confound. It therefore **re-derived the anti-correlation it was built
> to validate.** `γ_s` — the only exogenous factor — is **site-blind** (HE +0.19931 vs decoy +0.19844, p=0.20).
> The clean reduced form, **clustered at the arm** (1,793 edges are only ~330 arms): **p=0.115, n.s.**
> **(2) The within-arm two-way-FE replacement is ALSO retracted** (MH-133, 2026-07-16): its control class was
> **71% real binders**; against a genuinely site-free control **τ=−0.0007, p=0.84**, and it fails the site-type
> efficacy ladder (8mer/7mer-m8/7mer-A1 not monotone, p=0.26) — copy number cannot see a site type, so genuine
> site-mediated repression *must* follow it. **These are refutations at n=216k–235k pairs, not power failures.**
> ⇒ **Edge existence rests on ONE observational line.** Do not cite this section as causal evidence; do not call
> `pi_causal` exogenous. The `F>10 ∧ T1-clean` admission rule remains a **validity** rule (exclusion fails for
> ~73% of well-instrumented edges) — it was never a reality test. Current state:
> [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) Axis 2.

miRNA-locus copy number CN(m) instruments arm dose:

    first stage:  x_m = γ·CN(m) + C·δ + u        F = ((RSS_reduced − RSS_full)/1)/(RSS_full/df)   usable if F>10
    reduced form: rho_reduced = spearman( resid(CN | C), resid(Y | C) )     expect < 0 (genetic dose ↓ target)
    observational: rho_obs   = spearman( resid(x_m | C), resid(Y | C) )
    causal_concordant = (rho_reduced < 0) AND (rho_obs < 0) AND F>10
    p_reduced = one-sided Spearman-p(rho_reduced<0) ;  q_reduced_by = BY(p_reduced)

**Cluster-exclusion** (`cluster_attribution`): CN moves the whole genomic cluster. For co-located miRNAs that
also target g (set S = {m} ∪ co-targeters), condition each on the others+C:
`part(k) = spearman( resid(x_k | {x_{j≠k}} ∪ C), resid(Y | {x_{j≠k}} ∪ C) )`. `arm_unique = part(m)<−0.1 AND
no co-targeter survives`. `strong_causal = concordant AND (cluster_clean OR arm_unique)`; `strong_causal_fdr`
adds `q_reduced_by<0.05`.

## §9 DRIVER RESOLUTION inside a family  (`families.resolution_report`)
For a model-selected multi-member family, each member's anti-corr is measured *net of its family-mates* (+C):
`part(arm) = spearman( resid(x_arm | {mates}∪C), resid(Y | {mates}∪C) )`. A member **survives** if
`part < −0.1` (real conditioned repression). The nomination reflects **how the surviving anti-corr
distributes** (not a forced single winner):
- **0 survivors → "family-only"** — the split is unidentified (collinear members lose their residual variance
  ⇒ part→0, auto-fail; e.g. ESR1 miR-221/222).
- **1 survivor → "single-driver"** (e.g. PTEN miR-181b, miR-429).
- **≥2 survivors → "co-drivers"** — the repression is carried by several members, reported as a set
  (e.g. PTEN miR-17~92 = **miR-106b −0.22 + miR-20a −0.13**; ZEB1 = **miR-429 + miR-200c**).
Columns: `resolution` (the class), `n_drivers`, `drivers` (the surviving set, strongest first), `driver`
(strongest survivor, compat), `member_partials` (every member's part). Empirically ~13:7:11 single:co:family-only.

## §10 UNCERTAINTY — hierarchical Bayes  (`hierarchical._gibbs`)
Per program (genes sharing regulators), on C-residualised z-scored X:

    y_g = X_g β_g + ε_g ,  ε_g ~ N(0, σ_g²)
    β(m,g) ~ N(μ_m, τ²)        # miRNA-level pooling ACROSS the genes m targets
    μ_m ~ N(0, ω²) ,  τ², σ_g² ~ InvGamma

Gibbs (conjugate): β_g|· ~ N(V(XᵀX/σ² + μ/τ²), V), V=(XᵀX/σ²+I/τ²)⁻¹; μ_m|· pooled across its genes; τ²,σ_g²|·
inverse-gamma. Output: posterior mean ± sd of β; **identified edge = |β/sd| > 2**. (Cross-gene magnitude pooling
does NOT improve OOF coupling here — used as the uncertainty object, coupling stays per-gene §6.)

## §11 SUBTYPE tests  (`subtype.py`)
- **Level** (Kruskal): `kruskal_across_strata( X·M , PAM50 )` — does pressure magnitude differ across subtypes.
- **Slope** (interaction): nested rank-OLS F-test `rank(Y) ~ rank(ρ_agg)·subtype + C` vs common-slope — does the
  COUPLING differ; BH across genes.
- **Wiring, independent** (`subtype_wiring`): fit canonical `M_fam` (§5) SEPARATELY per subtype; mean pairwise
  cross-subtype Spearman(M) vs an **n-matched noise floor** (random cohort subsamples). High variance (differences
  two small-n fits) — *over-calls* remodeling.
- **Wiring, POOLED (Decision H, `subtype_wiring_pooled`):** whole-cohort `M_all` (§5) is the **prior mean**;
  per subtype `M_s = M_all + δ_s` via ridge-toward-M_all, `δ_s = argmin_δ ‖(yr − Xz·M_all)[s] − Xz[s]·δ‖² +
  λ‖δ‖²` (closed form; λ=10). Δ is regularised → small-n noise shrinks to 0; remodeling = `‖δ_s‖₁` vs an
  n-matched null. **Corrects the independent version:** PTEN borderline (~1.35× across λ=1–30, not robust),
  only RB1-Her2 modest (~1.5×) ⇒ subtype wiring mostly conserved. This is the design's `M_T = M_H + Δ` nesting,
  subtype flavour (state flavour: same with `M_H` from GTEx as the prior mean).

## §12 DISCOVERY  (`discovery.py`)
Orphan edge (m,g) beyond the curated POOLED-HE regulators: partial-Spearman(resid(x_m|C∪pooledHE-aggregate), resid(Y|·)) <
min_partial, gated by scanMiR biochemical support and (deconv) composition retention. Edge-level FDR from a
permutation null-scan (real vs permuted hits). Family-level: collapse to seed-family→gene hypotheses (§14).

## §13 ENRICHMENT  (`enrichment.py`, `stats.hypergeom_enrichment`)
For a target set (discoveries / a hub miRNA's targets) vs Hallmark program P over the HE universe U:
`fold = (k/|P|)/(K/|U|)`, `p = P(X ≥ k)` hypergeometric with (|U|, K, |P|); BH across programs.

## §14 MULTIPLE-TESTING machinery  (`stats.bh_fdr`, `coupling_inference.{benjamini_yekutieli, simes_p, family_simes_fdr}`)
- **BH**: q_(i) = min_{j≥i} p_(j)·n/j (step-up).
- **Benjamini-Yekutieli** (arbitrary dependence): BH with α/H_n, H_n = Σ_{i=1}^n 1/i.
- **Simes per family**: family p = min_i p_(i)·k/i over the family's k tests; BH across families → per-test q.
- **Family-wise min-statistic null**: per draw take the family's most-extreme null statistic;
  p_family = (1+#{null_min ≤ obs_min})/(n_perm+1).
Hypothesis unit = seed family (Design §F): family-mate arms are ONE test, not N.

## §15 IDENTITY vs MAGNITUDE attribution  (`attribution.identity_vs_magnitude`, Design §Decision I)
Over the gene's nonzero seed-families, two objects — never conflated:
- **MAGNITUDE** (budget) = `M_fam·X_fam` normalised (§5a) — realized force (abundance-included).
- **IDENTITY** = **Shapley** credit for the FIXED-weight aggregate's explained variance:
  `φ_f = mean over orderings π of [ R²(S∪{f}) − R²(S) ]`, where `S` = f's predecessors in π,
  `pred_S = Σ_{f'∈S} M_{f'}·X̃_{f'}`, `R²(S) = corr(pred_S, ỹ)²`, X̃/ỹ = C-residualised (X̃ z-scored). Sampled
  orderings (n_perm=400; = **LMG** for the linear aggregate — exact = 2^p). `identity_share = φ_f / Σφ`.
- `gap = identity − magnitude`: **+** = explains more of Y than its budget (quiet on-target owner, e.g. ESR1
  miR-18: identity 0.54 / budget 0.20); **−** = high budget but doesn't explain Y (loud passenger, e.g. ESR1
  miR-22: budget 0.26 / identity 0.01). Budget rewards abundance (M·X̄); identity rewards on-target coupling
  — the gap is *primarily* abundance-vs-coupling, and *secondarily* collinearity (co-varying families split
  the shared R²). **Discipline:** identity ⊥ coupling (§6/§7) — the share distributes the coupling that
  *exists*; it is never evidence a driver exists. This RETIRES the precursor softmax-share identity.

## §16 STRUCTURAL IDENTITY — the loss lens  (`structural_identity.py`, Design §Decision I loss variant)
Expression-free "who is DESIGNED to specifically repress g", so a *silenced* specialist (M/budget/Shapley all 0
— no variance, no abundance) still surfaces. Candidate families = pooled-HE (§0) ∩ scanMiR-modelled arms, with a
**TargetScan-context++ fallback** for arms scanMiR lacks (recovers e.g. miR-137). Per seed-family f of g:
- `potential(f,g)` = biochem strength (scanMiR |repression| / TargetScan ts_mag) × (1 + log1p(evidence ledger)),
  direction-filtered (drop TarBase-v9 validated-activating). Family = its strongest member.
- `specificity(f,g)` = `Σ_{members} ev(·,g) / Σ_{members} Σ_{g'} ev(·,g')` — the family's POOLED evidence-mass
  fraction aimed at g (pool the WHOLE seed family, not one representative arm — else a thinly-studied member fakes
  concentration). There is NO biochemical specialist (every arm hits ~700 targets ~uniformly), so specificity is
  purely an evidence construct — no biochemical floor.
- `confidence(f)` = `min(1, pooled_ev_depth(f) / 10)`, ev_depth = # genes with any family evidence.
- `structural(f,g)` = `potential · specificity · confidence` (confidence a MULTIPLIER: no evidence ⇒ no ownership).
  `structural_share` = gene-local normalisation.

`lost_specialist(gene)` overlays baseline-active (GTEx/NAT > log2(11)) ∧ tumour-silenced (< floor) + §17.

## §16b SEQUENCE-DESIGNED OWNERSHIP — affinity-percentile attribution-identity  (`attribution_identity.py`, 2026-07-09)

A THIRD identity object: abundance-removed, expression-free SEQUENCE-designed owner. Distinct from §15 (Shapley
REALIZED coupling ownership, from TCGA) and §16 (structural = potential·specificity·confidence). Non-circular
(sequence + AGO chimeric, never TCGA X/Y).

- L1 FAMILY ownership (seed-family unit, §4; `seq_specificity`):
    aff(F,g)  = mean over co-seed matures m∈F of Σ_sites |weighted context++|(m→g)   # matures pooled by seed_family_map
    pct(F,g)  = rank-percentile of aff(F,g) within F's genome targetome                # "is g among F's strongest targets"
    own(F|g)  = softmax_F( z(pct(F,g)) )                                               # within-gene ownership distribution
  PERCENTILE, not concentration-share (share confounded by targetome size → crowned miR-96 over miR-200 on ZEB1);
  context++ > scanMiR K_D (checked). ⚠ **SUPERSEDED (MH-87):** that comparison was HE-restricted — on the GENOME-WIDE
  universe, all-site scanMiR K_D BEATS context++ at per-arm specialist recovery (0.89 vs 0.79 fair); the HE-scan was
  the artifact. `kd.genome_affinity` / `seq_specificity.affinity_percentile_kd`. (At family-L1 the swap is mixed — K_D's
  per-arm edge is muted at family resolution — so context++ stays the L1 default; K_D wins per-arm.) Recovers canonical specialists (miR-19/26→PTEN, miR-15/16→CCND1, miR-200/141→ZEB1).
- L2 ARM nomination within F (`chimeric_evidence`): distribution over member arms from pooled chimeric provenance
  (Manakov + TarBase qCLASH/CLASH, overlap-resolved) else expression×loading. Combined = own(F|g)·share(arm|F).

VERDICT: DESCRIPTIVE model (breast-binding-enriched p=7.5e-7 seed-caveated; arm-level 91% cross-source) with NO
regression value — inert on §6 coupling (p0×spec_gain sweep), §15 magnitude, §12 discovery, identifiability
(`identity_payoff_genome` 0% flip). Identity ⊥ what the model estimates; a hypothesis-generator, not a prior.
See VALIDATION §1, SYNTHESIS §6c.

## §17 LOSS/GAIN MECHANISM — promoter methylation gate  (`methylation.py`)
Tumour−normal promoter Δβ over CpG probes DIRECTLY overlapping the miRNA promoter (full sesame HM450 manifest,
485,577 probes; window [TSS−2000, TSS+500] around the pri-miRNA hairpin; **no cCRE hop** — the gene-centric probe
annotation misses miRNA loci, which is why the pipeline fell back to cCRE). Group means over 802 tumour vs 99
normal raw betas (no barcode map needed). **Bidirectional:**
- `hyper` = Δβ ≥ +0.15 ∧ tumour β ≥ 0.20 → arm SILENCED → repression LOST (target ↑).
- `hypo`  = Δβ ≤ −0.15 ∧ normal β ≥ 0.20 → arm DE-REPRESSED → repression GAINED (target ↓).

Validated: TSG-miRNAs mean Δβ +0.061 vs oncomiRs −0.070 (miR-124 +0.193, 48 CGI probes; miR-21 −0.176). Promoter-
probe-union betas cached at `output/learned/mirna_promoter_betas.parquet`.

## §18 ESTIMATOR SELECTION — what estimator for what job

**There is ONE model:** a **dense learned-τ² non-negative Bayesian posterior per gene**
(`spike_slab._gibbs_posterior`), read **two** ways — `π ≡ 1` (DENSE) for coupling + attribution +
identifiability, **evidence-π** (INCLUSION mode) for discovery. The half-normal slab `θ_m ~ N⁺(0, ν²)`
encodes "miRNAs repress"; **ν² is LEARNED, not guessed** — that, not L2, is the edge over the lasso
(`nnridge_cv ≈ lasso`). The **adaptive lasso (§3) is a BASELINE, not the model.** Canonical statement of the
converged model: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b; the *why* per estimator: `RATIONALE.md` §2a–§2g.

**Inclusion policy (section-wide):** every estimator below draws a gene's candidate regulators from
**POOLED-HE** (§0, `ledger.pooled_he_edges()`). Scope is the **Hallmark universe by design**;
`he_only=False` → the full hallmark set.

| # | Job / estimand | Estimator | Module | Chosen over |
|---|----------------|-----------|--------|-------------|
| 1 | **Coupling** — does repression exist; aggregate pressure | Aggregate `X·E[β]` of the **dense (π≡1)** posterior, on C-residualised z-scored family predictors | `spike_slab._gibbs_posterior` (π≡1) | **adaptive lasso** (§3) — OOF mean ρ **−0.168 vs −0.152**, Wilcoxon **p=9e-16**: at n≫p the lasso's sparse selection over-shrinks the aggregate, dense learned-τ² keeps the diffuse many-arm signal; **CV-tuned NN-ridge** (`nnridge_cv ≈ lasso` ⇒ the win is the *learned* ν², not L2); **inclusion mode** (selection discards the diffuse signal) |
| 2 | **Attribution — magnitude** (how much does this family contribute) | Posterior **mean ± sd** of the *same* dense fit | `spike_slab._gibbs_posterior`, `attribution_eb` | **bagged NNLS** (§5) — equal *and* unified: agree **0.80**, reproducibility **0.97 > 0.98·NNLS**, split-half ρ **0.822 vs 0.729**; the calibrated sd comes free, bagging only approximates it |
| 3 | **Attribution — identity/ownership** (who *owns* the coupling, abundance-removed, collinearity-fair) | **Shapley/LMG** decomposition of the fixed-weight aggregate's R² across families (§15) | `attribution.shapley_identity`, `identity_vs_magnitude` | **budget share** (§5a — that is *force*, abundance-included: a different question); raw per-family ρ (double-counts shared variance under collinearity). **Discipline: identity ⊥ coupling** — it says *who*, not *whether* |
| 4 | **Identifiability** (is this weight resolved) | **`\|z\| = \|mean/sd\| > 2`**, read off the same dense posterior | `spike_slab._gibbs_posterior`, `attribution_eb` | the **§9 conditioned-partial** — equal (precision **0.86**, recall **0.89**); `\|z\|>2` **IS** the full-conditional partial-Spearman, shrinkage-stabilised (Spearman **0.88**) ⇒ identifiability needs no separate estimator. Cross-gene hierarchical pooling (§10) does not improve coupling |
| 5 | **Discovery** (which *new* orphan edges exist) | **PIP** from the **inclusion mode** (evidence-π / learned base-rate spike-and-slab) + permutation-FDR + scanMiR duplex gate + deconv composition-robustness gate | `spike_slab._gibbs_ss`, `priors.inclusion_prior`, `discovery.py` (§12) | the dense **`\|z\|`** (over-conditions: discovery must nominate a *representative* of a collinear orphan group — model-averaging — which the full-conditional `\|z\|` drops; the two are **complementary**, SYNTHESIS §4); NNLS / fixed-support (structurally *cannot* select); the lasso (baseline) |
| 6 | **The cross-cohort GAUGE** (`β_source = a·β_target + c + δ`) | **Bagged NNLS** β/sd — **NOT** the Gibbs posterior | `gauge.py` (`_bagged_nnls_meansd` via `attribution_eb`) | the dense Gibbs posterior — see the **carve-out** below |
| 7 | **OOF prediction / the validation gate** | Out-of-fold **partial-Spearman** of the aggregate vs held-out mRNA \| C (§6), 5-fold over participants | `mvp.oof_gate`, `mvp.gate_fdr` | any in-sample fit statistic (the estimand is held-out prediction). *Which comparator the gate is scored against is a validation question, not an estimator one* → `STATE_OF_PLAY.md` Axis 4 |
| 8 | **Per-edge coupling readout** (single edge, per state) | **Partial-Spearman** of arm vs gene \| C (§7) — estimator-free | `coupling.edge_coupling`, `card.py` | any fitted M (a single edge needs no model; direct partial-ρ is assumption-light and robust) |

**Engine: Gibbs, not HMC — measured, not aesthetic.** Every live channel is **Gaussian-conjugate**; the
Student-t likelihood is a **scale mixture of normals** and stays conjugate (§3b); the protein discordance link
is **additive** ⇒ conjugate. HMC was warranted only for the A3 program-wise field, which the A4 probe **gated
out** (pooling is a wash, Δρ **−0.0006**). Only a future **non-conjugate** channel reopens this — binding =
Poisson/NB, methylation = Beta.

**⚠ The GAUGE carve-out (job 6) — fit it with bagged NNLS, keep Gibbs for the model.** Fed the dense Gibbs
posterior's β/sd, the gauge returns **a = 4.1** where the truth is **1.0**. **Gibbs is not the worse
estimator** (split-half reproducibility ρ **0.822** vs bagged NNLS's **0.729**) — the **errors-in-variables
correction** is what breaks: it divides by `Var(β̂) − mean(se²)`, and for Gibbs `mean(se²)` is dominated by a
heavy tail of a few enormous posterior SDs (`sqrt(mean(se²)) = 0.055` against a *typical* `se` of **0.015**),
collapsing the denominator (reliability **0.17** vs NNLS's **0.72**). `gauge._MIN_RELIABILITY` **refuses** that
gauge rather than silently returning 4.1.

**Where the prior acts.** The prior enters as **inclusion** (π, evidence-graded, §2) and as slab-**scale**
channels (μ = scanMiR biochemical magnitude, τ = evidence depth); the slab **location stays at 0**, so the
prior sets *ordering*, the data sets *magnitude*. The **dense (π≡1) fit is prior-inert by design** — the
prior's work is in inclusion/discovery, not in coupling.

---
*Every formula above has a runnable implementation at the cited module.function; see the reproduce index in
LEARNED_MODEL_VALIDATION.md.*
