# CN-locus instrument — design & identifiability

The causal layer of the learned model: use miRNA **copy number** as an instrumental variable to turn
"arm abundance correlates with target" (observational, confoundable) into "the arm's genetic dose **causes**
the target to move." This doc pins down the goal, the mechanism, and — the hard part — **what CN can and
cannot resolve**, before any refactor. Module: `learned/instrument.py`. Status: **active** — gap 1 (per-locus CN),
gap 2A (multi-IV + Hansen-J over-ID + locus-of-origin), **Ring-1 (within-family causal attribution: p_S ×
dose-delivery)**, the **gap 2B between-family exclusion FREQUENTIST INTERIM** (2026-07-09/10), the **coding-gene
pleiotropy gate** (`coding_pleiotropy`, MH-99) and its a-priori extreme the **intronic-HOST exclusion gate**
(`host_pleiotropy`, MH-101, 2026-07-11) all DONE; remaining = **gap 2B Bayes unification (CN-anchor on the cross-family
posterior covariance)** (§9).
⚠ **TWO NEW GAPS — READ `§9x` BEFORE USING THIS INSTRUMENT (MH-126, 2026-07-13):** (a) the confounder block is
**missing an ANEUPLOIDY/PLOIDY term**, which is why **86.4% of *impossible* (site-free) edges show a NEGATIVE
CN→target reduced form** — a derived per-sample ploidy fixes the calibration (**→ 50.0%, exactly the null**) but
is **not yet persisted**; and (b) **⛔ `pi_causal` is NOT exogenous** — it is `γ_s · b_fam` = first stage × the
**observational** OLS slope, and the only exogenous factor `γ_s` is **site-blind (p=0.20)**. This **downgrades
`MH124_ANTICOUPLING_VALIDITY.md` §5's CN existence leg.**
Companion to `COLLINEARITY_AND_IDENTIFIABILITY.md` (the CN overlay is that doc's Kind-B tool).

---

## 0. Goal

Observational coupling (arm abundance vs target mRNA) cannot separate **causation** from **confounding** — a
shared program (proliferation, EMT) can drive both the miRNA and the target. A miRNA's **locus copy number** is
a *genetic* perturbation of its dose that is (largely) not caused by the target → an **instrument**. If the
locus is amplified/deleted and the target then moves, that is causal, not shared-program covariation.

## 1. Mechanism — `edge_instrument(arm, gene)`

Three C-residualised regressions:
- **first stage**  `x_arm = γ·CN_locus + Cδ + u`  — does locus CN move the arm's dose? `F > 10` ⇒ usable.
- **reduced form** `y_g  = π·CN_locus + Cβ + e`  — does the genetic dose move the target? `π < 0` ⇒ repression.
- **observational** `partial-ρ(x_arm, y_g | C)` — for comparison.

**causal-concordant** = `reduced < 0 ∧ observational < 0 ∧ F > 10`. Weak instrument (`F ≤ 10`) = INCONCLUSIVE,
not failed. `run_clean` runs this over well-instrumented HE edges → the STRONG-causal set, BY-FDR on `p_reduced`.

## 2. The exclusion restriction — why the unit is everything

Validity needs: **locus CN affects the target ONLY through this arm.** But **a CN event moves an entire
co-amplified region at once** — every miRNA on it changes together. So the reduced form is the effect of the
**whole region**, not the single arm. The whole design problem is: *what is the region, and how do we get from
a region-level causal claim down to an arm/family?*

## 3. The unit hierarchy — what CN can and cannot resolve

Three genomic scales; CN's resolution stops at the **cluster**:

| scale | example | share CN? | CN can separate? | separator below CN |
|---|---|---|---|---|
| **two arms of one hairpin** | miR-19b **-3p** vs **-5p** | identical (one segment) | **No** | expression (arm selection) + sequence |
| **hairpins in one cluster (polycistron)** | miR-18a/19a/20a in 17~92 | identical in practice (cluster ≪ CN segment; empirically miR-18↔miR-20a CN ρ=1.00) | **No** | expression (mature abundances diverge, ρ≈0.55) + sequence |
| **across clusters / loci** (incl. paralogs on different chromosomes) | 17~92(chr13) vs 106a~363(chrX); the two miR-200 clusters | **independent** | **Yes** — multi-IV (§5) | — |

**Consequences:**
- CN's finest unit is the **cluster** (co-amplified region). It is assigned to the **combined** unit, never to
  a sub-arm/hairpin/family within it.
- **Complementarity:** across clusters, expression is confounded (co-transcription) but CN is exogenous → **CN
  wins**; within a cluster, CN is blind (identical) but mature abundances diverge → **abundance-coupling wins**.
- **Testable exception:** a CN breakpoint inside a cluster gives its hairpins slightly different CN → CN can
  then partially separate. Detect by same-cluster CN ρ < 1.0 (rare: clusters ~kb, CN segments ~Mb).

## 4. The four complication cases

- **(a) same seed family, same gene, same locus.** Same seed → identical sites, same CN → neither dose nor
  sequence splits them. → **collapse to family** (Kind A). Member "resolution" is at most a *dose-delivery*
  attribution (which copy carried the informative abundance variation), NEVER a repression-identity claim —
  same-seed members repress g identically per molecule (§8).
- **(b) co-targeters share the locus — different OR same family** (17~92 → PTEN). CN is the **cluster** effect;
  can't split by CN (same dose). Resolve *within* the cluster by expression conditioning (abundances diverge) +
  sequence, **distributing a portion** to each g-targeting member/family (§7). The *operation is the same* for
  different- and same-family co-targeters; only the **interpretation** differs: different-family portions attribute
  *different repressive mechanisms* (different sites); same-family portions attribute the *same* mechanism to
  *different dose-deliverers* (§8). NOT "one arm uniquely carries it."
- **(c) an arm from multiple loci** (miR-19b-3p = chr13 miR-19b-1 ⊕ chrX miR-19b-2). Its mature CN in `arm_cn`
  is a **mixture** → contaminated instrument. Fix = **multi-IV** (§5), not a sum.
- **(d) a seed family spanning many loci** (miR-19-3p = miR-19a chr13 + miR-19b chr13&chrX). One family = several
  independent-CN loci ⇒ the **family is the wrong unit** for CN (`_genomic_clusters` docstring says exactly this).

## 5. Multi-locus arms — multi-IV + over-identification, NOT summing

Summing/averaging an arm's multi-locus CN is **wrong**: the loci contribute **non-uniformly** to the mature
dose (`x = w₁·CN_L1 + w₂·CN_L2`, weights unequal and unknown), so `mean(CN)` doesn't track dose. Instead use
each locus's CN as a **separate instrument**:
- first stage `x ~ γ₁·CN_L1 + γ₂·CN_L2 + C` — the coefficients *estimate* each locus's contribution (handles
  non-uniformity rather than assuming it);
- reduced form `y ~ π₁·CN_L1 + π₂·CN_L2 + C` — a causal coefficient per locus.
- **Over-ID (Hansen-J, heteroskedasticity-robust — the built test; Sargan is its homoskedastic special case) bonus:**
  the two loci make the *same* mature product, so the Wald ratios `π₁/γ₁` and
  `π₂/γ₂` should be **equal** (both act through one molecule). Agreement confirms exclusion; disagreement is the
  **diagnosis that one locus is a SHARED cluster** (its neighbours also hit g) — **do NOT drop it, RECURSE**:
  run that locus through the whole tree (§7, SHARED → portion), extract A's portion of it, and reconcile with
  the clean locus. Multi-IV then *combines A's portions across its loci*, clean loci as cross-checks.

**AS-BUILT in `exclusion` (2026-07-11): instrument columns are DECOUPLED from the exposure.** X_fam stays per-**arm**
(`log2 Σ_m 2^{x_m}`, aggregated mature — §6), but each arm contributes one CN **instrument column per ACTIVE-SOURCE locus**
(`_arm_active_loci`), not just its focal locus. **Silent-source loci are excluded — they aren't valid instruments
(relevance fails: their CN doesn't move the mature dose), so a SINGLE-active-source arm's focal locus IS the complete
instrument set and the focal simplification is EXACT, not an approximation.** Empirically this is **783/785** HE arms; only
**2** (miR-550a-3p/5p) have ≥2 active sources, and they now get all their loci as separate instruments (each feeds the
segment grouping + the within-arm Hansen-J over-ID). So `exclusion`/`run_exclusion` are now multi_iv-complete for the CN
channel; `family_multi_iv` still uses focal-per-member (its own scope).

## 6. The CN / expression asymmetry (the crux)

- **CN is observed per locus** — each genomic segment has its own copy number ⇒ **separable** ⇒ multi-IV.
- **Expression is observed per mature arm** — the sequencer counts mature reads regardless of which locus made
  them ⇒ **one aggregated number, not splittable**.

So: the **repression is by the mature product** (locus-agnostic — the chr13 and chrX miR-19b-3p molecules are
identical and repress g identically); the locus CN only supplies *independent genetic perturbations of that one
product's dose*. **Corollary for attribution:** when distinguishing families by expression, a multi-locus member
contributes its **single aggregated mature abundance** — you do **not** (and cannot) sum per-locus, because
per-locus expression is unobserved.

## 7. The corrected architecture

1. **Instrument per cluster/locus** (multi-IV for multi-locus arms; over-ID check across an arm's loci) →
   a **cluster-level causal** claim (`reduced < 0`, `F > 10`).
2. **Map** cluster → the mature arms it produces (many-to-many).
3. **Attribute within/across** those arms on the **aggregated mature abundances**, **distributing a portion** —
   a **two-stage** split, because Shapley alone over-credits collinear passengers (it *divides* variance fairly
   but cannot *exonerate*: a co-transcribed non-targeter B, ρ≈0.8 with a real driver A, still gets ≈½R²(B) from
   the *shared* variance). So:
   - **between families:** `conditioned partial anti-corr EXONERATES` (prune families whose *unique* residual
     variation doesn't survive → 0 — this is the thing Shapley can't do), then `SHAPLEY` distributes among the
     survivors (`attribution.shapley_identity`, shares sum to the whole). **portion = CN π × Shapley(survivors)**.
   - **within a surviving family (= the built Ring-1 L2, §8):** split by **measured chimeric binding**
     (`chimeric_evidence`, arm-resolved physical occupancy — which member actually binds g) where present, else
     **abundance** (level; members' mature abundances diverge — real data median member-abundance |ρ|=**0.44**,
     62% of pairs <0.8, only 8% ≥0.95 — so level dose separates them, collapsing to `family-only` only for the ~8%
     truly co-transcribed pairs). **No exoneration here** — every same-seed member is a real regulator (identical
     sites), so the split *distributes* dose-delivery (a low-abundance member → small share, not zero); the
     conditioned-partial's *binary pruning* is what's banned within a family, not the split. **Three orthogonal
     signals, three questions — NOT a supersession:** abundance **level** = the dose-delivery *share* (occupancy ∝
     abundance here, since same-seed ⇒ identical affinity/K_D and ~equal loading); **chimeric** = that same occupancy
     *measured directly* (gold); abundance-**variation**/Shapley = *coupling* (whose fluctuation tracked g). The
     within-family split is a **share/occupancy** question **and family-level causality is already settled** (the CN
     instrument + identical sites ⇒ every member is a real regulator), so nothing remains to attribute by covariation
     → **level/measurement is the right functional and variation-Shapley is the wrong one here** (a high-abundance,
     low-variance member delivers most of the dose yet earns ≈0 variation-share). Variation-Shapley is **not
     discarded** — it's the **between-family** tool, where causality is *not* settled and unique variation must
     exonerate. Sequence can't help either (K_D constant, 3' non-predictive, §9). Hence **L2 = the occupancy share,
     estimated CO-PRIMARILY by chimeric and abundance** — chimeric = the *direct measurement* (gold **where it exists**,
     but sparse + out-of-tissue), abundance-level = the *principled in-context estimate* (occupancy ∝ abundance for
     same-seed ⇒ **not** a weak fallback); use the measurement when present, the abundance estimate otherwise. Variation
     lives one level up. Interpretation: within-family shares =
     *dose-deliverers of one mechanism* vs between-family shares = *different mechanisms*.
   Perfect within-family collinearity stays honestly unresolved (no chimeric + tied levels ⇒ equal split; partial
   → NaN). Only the **between-family** level uses conditioned-exoneration + Shapley(survivors); the
   **within-family** level uses the occupancy share (measured chimeric, else its abundance-level proxy) — the
   **exoneration gate is between-family only**. **Never** a single-arm "unique carrier"; **never** per-locus expression (§6).

This **replaced** (built 2026-07-10 as Ring-1 / `family_causal_attribution`) `cluster_attribution`'s
`arm_unique_anticorr` (a special case: exactly one survivor) with a **portion distribution** (the general case),
and moved the CN off `arm_cn` (mature-aggregated) onto the **precursor/cluster CN** so paralog arms stop
contaminating.

**The whole thing is recursive, built from TWO primitives.** Every case — clean, shared (same or different
family), multi-locus arm, pleiotropic locus — is these two composed:
1. **CN separates *clusters*** — attribute a causal effect to a co-amplified region (across independent-CN loci,
   incl. paralog loci via multi-IV).
2. **Within a cluster, distribute a *portion*** — **between families** by expression + sequence (different seeds →
   different mechanisms; sequence/K_D *do* separate them) and, **among same-seed members**, as a *dose-delivery*
   **occupancy share** — chimeric (measured) or abundance-**level** (its proxy, occupancy ∝ abundance); sequence/K_D
   are constant within a seed family so they **cannot** split it, and *variation*-coupling belongs one level up
   (between families) (§8/§9).
A multi-locus arm = several instruments, each run through the tree; a pleiotropic locus = the SHARED branch
applied to that locus, then A's portion combined across loci. No special cases — cluster-separation outside,
portion-distribution inside, all the way down.

## 8. Kind-A note — the family is the repression unit

Same-seed members hit identical sites ⇒ they repress g identically per molecule. So the seed family is the
**biochemically identified** unit; any member split is a *dose-delivery* nomination (which member delivered the
identical-per-molecule dose), never a repression-identity claim. Ring-1 resolves it as a single **occupancy share**
(L2) at two fidelities of the *same* quantity: **measured** directly (chimeric — `chimeric_evidence`, arm-resolved
occupancy) or, absent that, **proxied by abundance-level** (occupancy ∝ abundance for identical seeds; `l2_src`
flags which). These are **not** the abundance-*variation*/Shapley functional — that is a *coupling* question that
belongs between-family (where causality is unsettled), not to a within-family share where the CN instrument has
already settled causality. Sequence cannot help either — K_D is constant within a seed family (corr 1.00) and the
3' supplement is non-predictive (§9). Abundance-level separates members because mature abundances diverge (median
member |ρ|=0.44, 62% of pairs <0.8); only where neither chimeric nor divergent abundance exists does the split
honestly collapse to `family-only` (the ~8% near-collinear |ρ|≥0.95 co-transcribed pairs). The split is never a claim that one member out-represses its family-mate — only
which delivered the (identical) dose.

## 9. What's built vs the gaps

**Built** (`instrument.py`): `edge_instrument` (per-edge causal, now on **per-locus** CN — takes an optional
`locus=` to instrument a specific hairpin, the gap-2 multi-IV hook), `pleiotropy` (co-located co-targeters =
exclusion check), `cluster_attribution` (within-cluster expression resolution — LEGACY `unique-carrier`, superseded
by Ring-1 `family_causal_attribution`; §7),
`_genomic_clusters` + `cluster_instrument` (**the locus-level instrument** — groups by genomic position, per-
cluster reduced form + `members_resolvable`; CN now sourced per-locus), `run_clean` (60 HE edges → 24 STRONG
causal: 17 cluster-clean + 7 pleiotropic-but-arm-unique). **CN source primitives (gap 1):** `locus_cn`
(participant × `MI*` hairpin CN), `arm_loci_map` (arm → contributing loci + w_norm, from the MIMAT paralog map;
entity `pre_gene_id` fallback), `_arm_focal_locus` (dominant-w_norm hairpin; flags `multilocus`).

**Gap 1 — CN source — DONE (2026-07-09).** `edge_instrument` / `cluster_instrument` now source **per-locus**
hairpin CN (`locus_cn`), each arm keyed to its **focal locus** (highest-w_norm hairpin) instead of the
mature-aggregated `arm_cn` composite, so paralog arms (c) no longer contaminate. Single-locus arms (all 42
well-instrumented arms; arm CN ≡ locus CN exactly) reproduce the run **bit-identically** (24 STRONG / 18
FDR-sig unchanged); the change decontaminates multi-locus arms (e.g. miR-16-5p focal chr3 locus F=16.1 vs the
paralog-mixed arm-CN F=15.0) and exposes per-locus instruments via `edge_instrument(..., locus=)`. Output gains
`locus_id` / `n_loci` / `multilocus` provenance columns. `arm_cn` retained as LEGACY (comparison only).

**Gap 2A — multi-IV + over-ID — DONE (2026-07-09).** `multi_iv(arm, gene)` = 2SLS over an arm's active-source loci
+ **Hansen-J** over-ID (heteroskedasticity-robust); `locus_origin(arm)` = the gene-independent **locus-of-origin
assay** (separates AVAILABILITY = CN variance from ORIGIN = first-stage F/ρ → active_source / silent_source /
cn_invariant — a silent-source locus is CN-varying but doesn't make the mature arm, e.g. miR-92a/19b **chrX
106a~363** F≈2 vs **chr13 17~92** F≈94, so it's excluded as an instrument, its target-associations flagged as
horizontal pleiotropy). Side-finding: `focal_locus_concordance` (per-arm CN→expr on the focal locus) recovers
paralog arms (miR-92a/19b) the mature-aggregated concordance diluted below threshold; `well_instrumented_arms`
now sources it.

**Ring-1 — within-family causal attribution (§7) — BUILT (2026-07-10).** `family_multi_iv(gene, mem)` = the
family-level 2SLS (exposure = family pool `X_fam = log2 Σ 2^x`, instruments = member CNs, per-segment reduced
form + Hansen-J over-ID on independent-CN directions) — **replaced the mean-CN collapse**. `family_causal_attribution`
composes it into the hierarchical portion: **`portion(arm) = p_S(segment) × within-segment dose-delivery`**, where
`p_S` = the exogenous CN reduced-form share per independent-CN segment (F>10-gated), and the within-segment split
(§8, L2) = the members' **occupancy share** at two fidelities of the same quantity — **measured** (chimeric,
`chimeric_evidence`, arm-resolved) or **proxied by abundance-level** (`l2_src` flags which); *not* the
abundance-variation/Shapley functional, which is a coupling question that lives between-family (§7/§8). Biologically validated: canonical top deliverers (miR-16→BCL2/CCND1, miR-200c→ZEB1,
miR-93/106b→PTEN). **CN-specific design consequence of the general leg results** (the results themselves — K_D
discard-overturn, loading-inert, 3'-dead-end — are GENERAL, see the pointer below + `LEARNED_MODEL_WHATS_NEXT §K_D`):
because **K_D is CONSTANT within a family** (same seed, corr 1.00) it carries NO within-family information, so it lives
in the CROSS-family identity layer, **not** the within-segment split; and because sequence (K_D + 3') cannot separate
same-seed arms while **loading does not transfer** across contexts, the within-segment L2 split is **dose + measured
chimeric only** (loading dropped). `arm_cn` retained LEGACY.

**Gap 2B — exclusion / between-family — FREQUENTIST INTERIM DONE (2026-07-10); Bayes-anchor + coding-gene REMAINING.**
`p_S` is confounded when co-located **other-seed** families also target `gene` (the cluster CN moves them too).
`between_family_exclusion(gene, mem)` (CLI `--gap2b`) wires the interim: (1) enumerate co-located other-seed
co-targeters via `pleiotropy()`; (2) **condition** each segment's reduced form on their pooled abundance →
`seg_pi_cond` (blocks the CN_seg → other-family-abundance → target mediation), report `pi_shift` + renormalised
`seg_pS_cond`; (3) **attribute** ownership among {our family} ∪ {co-located families} with the between-family
exoneration gate + Shapley (`attribution._exonerate_between` + `shapley_identity` wired to the CN reduced form) →
`between_family_share` + `family_exonerated`. Verdicts: `between_clean` / `robust` / `partly_attenuated` /
`confounded` / `confounded_rider` / `weak_instrument`. **Result (full HE universe, 4,922 (gene,family)):** `p_S` is
exclusion-robust for **~99%** — **4,489 (91%) `between_clean`** (no co-located other-seed co-targeter), and of the
**433 confoundable**, 236 robust / 25 partly / 13 confounded (+63 rider, 96 weak); **median |pi_shift| over
confoundable F>10 segments = 0.007** (90th pctile 0.053). The confound concentrates on the polycistronic oncomiR
clusters — miR-17~92/106b~25 on **CCND1** (pi −0.156→−0.073, family keeps 56% Shapley ownership), **PTEN** (52%),
**TGFBR2/MDM2/BCL2L11/MYC**, and passenger **miR-18a→CCND1** (`confounded`: its CN effect is really the cluster's).
Interim caveat: conditions on the other families' *abundance* (mediation adjustment), NOT their exogenous CN.
**REMAINING:** the **Bayes unification** — the CN reduced form as a **causal anchor on the cross-family posterior
covariance** (the exogenous-CN form of step 2/3, the "want both" item) + within-family δ-pooling (`M_member = M_fam + δ_m`);
family-collapse is the seam (within = frequentist here, between = Bayes covariance).

**Coding-gene + HOST pleiotropy — ONE JOINT GATE (`coding_pleiotropy` MH-99 + host fold-in MH-101, 2026-07-11).** The other
exclusion threat: a **co-amplified CODING gene** (incl. the intronic **HOST**) that regulates `gene` directly (path
Z_s→CN_c→X_c→Y bypassing the miRNA). Both are handled by ONE gate that conditions δ_s on {co-located coding} ∪ {host}
**JOINTLY** and emits a single `pleio_down_weight`∈[0,1] that **inflates s²_π** (weaker CN pull, NOT a subtraction). This
retired the earlier separate-gates-combined-via-**MAX** heuristic: MAX assumed the two explained the *same* slice of δ_s,
whereas co-amplified confounders are collinear — the joint correctly **de-double-counts** (ESR1/miR-26a joint 0.50 vs MAX
0.72; PIK3CG/miR-7 0.67 vs 0.90) and captures the case where a host and a coding neighbour explain *different* slices.
- **Candidates:** co-located coding genes = **CN co-amplified** with the miRNA-locus CN (`coamp`, recovers RB1@13q14 at
  1.7 Mb) AND host-gate-aligned **partial Spearman(coding,target\|C)≥0.3**; a **proliferation-survival guard** drops prolif
  *markers* on prolif *effectors*. Condition δ_s, **ACT ONLY ON REDUCTIONS** (inflation = collider → ignore). RB1 anchor.
- **HOST = a-priori GUARD-EXEMPT candidate (MH-101, axis BP):** a SENSE-intronic miRNA's locus CN co-amplifies its HOST by
  construction, so if the host RELATES to the target (`host_target_relatedness`, partial ρ(host,target\|CPE,**miRNA**)≥0.3
  — miRNA partialled OUT to isolate the host's OWN target path) it is added to the conditioning set **bypassing** the
  coamp/min_corr/surv filters (it has its own validated gate + annotation). This is why the scan's prolif-survival guard
  can drop **MCM7** (for miR-17~92/106b~25) yet the host gate retains it — E2F3/RBL1 show `MCM7*` alone.
- **δ_s-SIGN CONFIDENCE (`pleio_confidence`):** a residual-miR leak into δ_s has sign −sign(γ) (the miR represses), so δ_s
  of the SAME sign as the first-stage γ can't be lost-miR signal ⇒ `clean_antirepression` (unambiguous pleiotropy);
  opposite ⇒ `repression_direction` (could be residual-miR OR a co-targeting other family → grade via miR-partial +
  `prolif_verdict`). ⚠ For proliferation-gene hosts the flag overlaps MH-100's confound-vs-mechanism — don't blind-apply.
- **✅ MULTI-LOCUS host-locus matching — RESOLVED (per-locus host, 2026-07-11).** Previously the host was annotated per-arm
  and assigned to the arm's focal locus, so a context-mixed multi-locus arm's focal could miss the coding-host copy. Now
  **`genomic_context.locus_host_map()`** classifies each CN locus at ITS OWN hg38 coordinate — bridge: `_precursor_loci`
  maps the CN `locus_id` (`MI*`, = `locus_cn` columns / `arm_loci_map`) → (chrom,start,end,strand) → `_classify_locus`
  (verified 506/506) — persisted to **`locus_context.tsv`** (506 CN loci: 211 coding-host / 235 lncRNA-host / 32 intergenic
  / 28 antisense). `exclusion` (via `_locus_coding_host` + `_arm_active_loci`), `host_pleiotropy`, `dose_attribute_host_downweights`,
  and `host_target_relatedness` (now every distinct coding host per arm + a `host_locus` col; 403→**444 related**) all condition
  each **ACTIVE** locus on **its own** host. **Correctness gain:** a coding host at a SILENT (non-instrumented) locus is no
  longer applied to the active instrument — miR-30c→NFYC (active = LINC00472 lncRNA) and miR-194→BPNT1 (active = MIR194-2HG)
  correctly DROP; miR-16→SMC4, miR-26a→CTDSP2, miR-93→MCM7 (correctly located) retained. The 6 multi-coding-host arms
  (miR-26a→CTDSP2+CTDSPL, miR-218→SLIT2+SLIT3, miR-103a→PANK2+PANK3, …) now score both hosts. Well-instrumented host
  down-weights 37→**36**.
**Batch (full HE universe, `run_exclusion(coding=True, host=True)`, 8-worker BLAS-pinned + optimized loop — precomputed
lstsq projectors, vectorised co-amp, scan restricted to host families — now ~2.4 min after the hot-path caching speedup,
see §perf):** base table **45,187 seg-rows / 17,053 strong / 4,342→4,343 clean / 12,698→12,710 pleiotropic** (was 45,156;
the +31 is the `.N` universe refresh finally persisting — verified via a base-only `host=False` run = identical 45,187, NOT
a code regression); each strong segment now carries `seg_dose_share` (its share of the family's CN identification). **Per-locus
36 host down-weights** (per-locus dropped the 2 silent-locus spurious ones — NFYC/BPNT1 — see the multi-locus RESOLVED note
above), each graded by `grade_host_downweights` (the ambiguous `repression_direction` subset resolved via the OOF lens: 6
confound-confirmed / 6 over-control / 11 host-immaterial). **FULL coding scan (2026-07-12): the host-family `do_scan`
restriction is LIFTED** (affordable after the 17× scan speedup), so `coding=True` scans EVERY family → **3059 pure-coding
down-weights** (was 11 under the restricted scan; mean 0.405, 1086 ≥0.5) — ~18% of strong segments carry a co-amplified
coding confounder (the coding subset of the MH-89 "73% exclusion-failing" rate). Top host: IL6/miR-26a@CTDSP2*
0.87, IGF1/miR-28@LPP* 0.84, HMGCS2/miR-425@DALRD3* 0.72, SLC7A5/miR-744@MAP2K4* 0.71, PIK3CG/miR-7@HNRNPK* 0.67,
CCNE1/miR-16@SMC4* 0.51, ESR1/miR-26a@CTDSP2* 0.50, FBP1/miR-21-3p@VMP1* 0.43, the MCM7 cluster 0.30–0.36. Columns
`pleio_down_weight`/`s2_pi_pleio`/`beta_weight_pleio`/`pleio_confidence`/`pleio_sources` (host tagged `*`) in `exclusion.tsv`.
(`host_pleiotropy` remains as a standalone host-only diagnostic.)
**BUILT** (`LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md`): the CN reduced form enters the converged one-posterior
(`DISCOVERY_SYNTHESIS §6b`) as a **conjugate Gaussian channel** on β — `π̂_s ~ N(Σ_{m∈s} γ_m β_m, s²_π)`. The
**between-family** cross-family covariance anchor is `channel_cn.multi_family_terms` / `between_family_bayes`
(**MH-93/94**, run — correct/validated but WEAK: the between split is mRNA-likelihood-dominated, CN adds only its
~0.7%-weight exogenous slice), and the **within-family** δ-pooling `M_member = M_fam + δ_m` is
`within_family.delta_pooling` (**MH-98**, run — *"the last open Gap-2B Bayes item, now CLOSED"*). **δ-pooling is
DELIBERATELY a standalone inverse-variance estimator, NOT an unfinished "fold-in":** it is a *separate, orthogonal* latent
(F2, design §6-F — member-level *dose-delivery* δ_m; same-seed members SHARE β_f per §8, so there is no per-member repression
β to add as a β-channel term), and the **fully-joint hierarchical Gibbs form was TESTED and found WORSE in MH-98** (at n≈1000
the mRNA collinearity dominates any CN prior and buries the CN-resolved member + gives unphysical negative M). Folding
member-level latents into ONE joint posterior only becomes meaningful at the **A3 de-collapsed (arm×gene) field** — i.e. it is
subsumed by the A-axis fit-unit continuation, NOT a standalone CN item. Gap 2B is the first instance of that broader
**channel-fusion** frame (M as a latent fused from many channels, not just mRNA+prior).

**GENERAL infrastructure Ring-1 CONSUMES (not CN-specific — own homes, do not duplicate here):** the affinity leg is
**genome-wide scanMiR K_D** — validated/discard-overturned, ⊥ context++, plus its scan-scope + context++×K_D-fusion
forward work — which is GENERAL identity/discovery infra (serves `discovery.py` + the selection prior + `attribution_identity`
as much as Ring-1) and lives in **`LEARNED_MODEL_WHATS_NEXT.md §K_D`**. The L2 within-segment signals — **chimeric**
measured binding + **loading** (inert) — are the general arm-resolution tools (`chimeric_evidence`, `ago_loading`,
`DATA_SOURCES.md`). **CN-SPECIFIC here** = the instrument (gap 1/2A) + Ring-1 + the three gaps above (2B between-family
exclusion, coding-gene pleiotropy, Bayes CN-anchor); Ring-1 merely *uses* K_D/chimeric as its dose legs.

### 9x. ⛔ GAP — the ANEUPLOIDY covariate, and the exogeneity of `pi_causal` (MH-126, 2026-07-13)

Two findings that change how this instrument must be **used** and **reported**. Neither has been persisted into
`instrument.py` yet — this section is the record of what was measured and what must happen before it is.
Full detail: **`MH126_ANEUPLOIDY_AND_GRADED_ATTRIBUTION.md` §1**.

**(A) THE CONFOUNDER BLOCK IS MISSING A PLOIDY TERM — and that is why the RAW instrument fails.**
The instrument's `C` (`assemble_gene`) is `[CPE, target_cn, mal_prolif]` — **verified to contain no aneuploidy
term**. Consequence, measured on 3,111 edges / 482 genes (HE + abundance-matched **site-free decoys**):
**86.4% of IMPOSSIBLE (decoy) edges show a NEGATIVE CN→target reduced form.** The mechanism is **not**
co-amplification (which predicts a *positive* reduced form) but a **PLOIDY SIGN FLIP**:
`corr(resid(CN_locus | C), ploidy) = +0.437` (100% positive) while `corr(resid(Y | C), ploidy) = −0.135`
(negative in **94.6% of 482 genes** — `target_cn` is already in `C`, and TPM dilutes) ⇒ a negative reduced form
for **any** pair, site or no site. **One per-gene scalar explains R² = 0.554** of the reduced form.
**FIX (measured, not yet persisted):** add a **derived per-sample `ploidy`** = mean ASCAT3 CN over the 506
miRNA hairpin loci in `locus_cn()`. Decoy %neg **86.4% → 50.0% (exactly the null)**; HE 89.3% → 53.6%;
HE exclusion-clean 34.6% → **55.6%**; and the **first stage STRENGTHENS** (mean ρ +0.096 → +0.118) ⇒ it removes
nuisance, not signal. The **official** measures (Taylor arm-level aneuploidy score, FGA) do almost nothing
(81.5% / 82.9% alone); the derived ploidy does it all (49.9%).
**BEFORE PERSISTING:** (i) **leave-one-out ploidy** — the focal locus currently contributes 1/506 to its own
control; (ii) a **VIF / conditioning check** (`ploidy` is r=0.828 with `target_cn`, already in `C`);
(iii) **trace the ripple** — `exclusion.tsv`, `cn_edge_validation.tsv`, `cn_family_*.tsv`, `channel_cn` all go
stale. ⚠ **CONTRACT CHANGE: under the fixed `C` the exclusion GATE stops discriminating** (57.7% of decoys
"clean" vs 55.6% of HE) — it is no longer a validity filter.

**(B) ⛔ `pi_causal` IS NOT AN EXOGENOUS ESTIMATE — STOP REPORTING IT AS ONE.**
`exclusion()` returns `pi_causal = γ_s · b_fam`, where `γ_s` is the first stage but **`b_fam` is the OLS partial
slope of `Y` on `X_fam` given `[C, Z_s]`** — i.e. it is fitted on exactly the variation in `X_fam` that is
**orthogonal to the instrument** (the **endogenous** variation). Identity re-verified by re-running
`INS.exclusion`; `b_fam` reproduces the plain **no-copy-number-anywhere** partial correlation (ρ≈0.93).
* The **only** exogenous factor, **`γ_s`, is SITE-BLIND**: HE +0.19931 vs DECOY +0.19844, **ratio 1.004×,
  p=0.20**. All of the HE-vs-decoy "separation" lives in the observational `b_fam`.
* An independent verifier **reproduced the whole "3× exogenous CN validation" from a simulation with ZERO
  causal effect.**
* The one genuinely exogenous statistic — the **CN→target reduced form** — shows **no HE-vs-decoy effect** once
  the unit is the **LOCUS** (τ = −0.0015, **p=0.67**) or once arms are matched on their own site-free
  background (p=0.038 → **p=0.71**).
⇒ **`MH124_ANTICOUPLING_VALIDITY.md` §5's "exogenous CN instrument" existence leg is DOWNGRADED.** Any claim
that needs exogeneity must use the **reduced form** or a proper 2SLS/LATE ratio. **Edge existence currently has
NO surviving exogenous validation** (its other two legs are observational) — restoring one is the top open item
on this arc.

**(C) UNIT OF ANALYSIS.** 3,111 edges = only **482 genes / 409 arms / 277 loci** (ICC of `pi_causal`: gene 0.372
/ arm 0.352 / locus 0.300), and per-arm abundance matching is a **per-arm constant** so it cannot remove a
locus-level confound. Under **gene + arm two-way FE, 59–72% of any HE-vs-decoy gap is composition.** Every
HE-vs-decoy test on this instrument must cluster at the gene/arm/locus, not the edge.

## 10. Code map

| object | function |
|--------|----------|
| per-edge causal test | `instrument.edge_instrument` |
| first-stage F | `instrument._first_stage_F` |
| pleiotropy / exclusion check | `instrument.pleiotropy` |
| within-cluster expression resolution (LEGACY unique-carrier — superseded by `family_causal_attribution`/Ring-1) | `instrument.cluster_attribution` |
| genomic-cluster grouping (the CN unit) | `instrument._genomic_clusters` |
| **locus-level instrument** | `instrument.cluster_instrument` |
| well-instrumented arm filter | `instrument.well_instrumented_arms` |
| pipeline → STRONG causal set + FDR | `instrument.run_clean` |
| **per-locus CN matrix** (instrument source, gap 1) | `instrument.locus_cn` |
| arm → contributing loci (+w_norm) | `instrument.arm_loci_map` |
| arm → focal (dominant) hairpin locus | `instrument._arm_focal_locus` |
| arm → ALL active-source instrument loci (multi_iv-complete; focal-only for single-active = 783/785 arms) | `instrument._arm_active_loci` |
| **multi-IV + Hansen-J over-ID** (gap 2A) | `instrument.multi_iv` |
| **locus-of-origin assay** (active/silent source) | `instrument.locus_origin` |
| focal-locus CN→expr concordance (side-finding) | `instrument.focal_locus_concordance` |
| **family-level 2SLS** (per-segment p_S, over-ID) | `instrument.family_multi_iv` |
| **Ring-1 within-family attribution** (p_S × dose-delivery; +`portion_abs`=share×β, `cn_driven_fraction`) | `instrument.family_causal_attribution` |
| batch Ring-1 (all families×genes; CLI `--ring1`) | `instrument.run_family_causal_attribution` |
| **Gap 2B between-family exclusion** (condition p_S on co-located other-seed co-targeters; `seg_pi_cond`/`pi_shift`/`between_family_share`) | `instrument.between_family_exclusion` |
| batch Gap 2B (all families×genes; CLI `--gap2b`) | `instrument.run_between_family_exclusion` |
| **CN-channel exclusion gate** (over-ID + T1 on X_fam → CONTINUOUS `pi_causal`+`s²_π`; scope-3 for CHANNEL_FUSION; `coding=`/`host=` add the pleiotropy down-weights) | `instrument.exclusion` |
| batch exclusion gate (all families×genes; `workers` parallel; CLI `--exclusion`) | `instrument.run_exclusion` |
| **JOINT coding+HOST pleiotropy gate** (co-amplified coding ∪ guard-exempt host, conditioned together → `pleio_down_weight`; `host_genes=`/`scan=`; MH-99+101) | `instrument.coding_pleiotropy` |
| standalone host-only diagnostic (per-edge host attribution + δ_s-sign confidence; now per-locus) | `instrument.host_pleiotropy` |
| **PER-LOCUS host map** (`MI*` locus → host classified at its OWN coordinate; the multi-locus fix, MH-101 per-locus) | `genomic_context.locus_host_map` → `locus_context.tsv` |
| CN locus → its own coding host (exclusion consumes this instead of the arm representative) | `instrument._locus_coding_host` |
| estimator-grade DOSE-ATTRIBUTION of host down-weights (`host_beta_share` = host locus's share of family CN identification; abundance vs δ-pooling dominant-arm) | `instrument.dose_attribute_host_downweights` → `host_dose_attribution.tsv` |
| GRADE the ambiguous `repression_direction` host down-weights via the OOF lens `dC` (confound_confirmed / over_control / host_immaterial / unresolved) | `instrument.grade_host_downweights` → `host_downweight_grades.tsv` |
| gene×participant frame → cached (values, name→i, name→j) for numpy positional pulls (the 17× coding-scan fix) | `instrument._as_np` |
| per-segment `seg_dose_share` (= `beta_weight` / family-total; the estimator-grade dose-leverage annotation) | column in `exclusion.tsv` |

**Perf (2026-07-11/12, PURE CACHING/VECTORIZATION — outputs proven identical):** the exclusion/learned hot path is **6.4× faster** (3.63→0.57 s/gene base; ~24→2.4 min) AND the **coding scan is 17× faster (117→7 ms/family)**, so the full `coding=True host=True` batch runs in **3.8 min** (was ~48 min). All cProfile-driven; the costs were repeated data loads / per-item scipy, not the model numerics. **Base loads:** `families.family_of` rebuilt the arm→family map on ~1030 calls/run; `instrument.locus_origin` re-read the full miRNA matrix + clinical per new arm; `data._target_cn` re-read the CN matrix per gene; `data.assemble_gene` copied the whole edge table per gene — all now cached. **Coding scan:** the proliferation PC was re-SVD'd per FAMILY (it's gene-level) → cached per gene; the per-candidate scipy `spearmanr` → a vectorized rank+matmul (verified 6e-16); pandas `.loc[genes, parts]` block reindex (77% of the scan) → numpy positional indexing on cached arrays (`_as_np`; verified `array_equal`, no dup index). Verified base-identical (a `host=False,coding=False` run gives the same 45,187 seg-rows; the RB1/coding down-weights match).
| between-family exoneration gate + Shapley (wired to the CN reduced form) | `attribution._exonerate_between` + `attribution.shapley_identity` |
| genome-wide scanMiR K_D (all-site, validated) | `kd.genome_affinity` / `kd.genome_affinity_pct` |
| per-arm global K_D specificity | `seq_specificity.affinity_percentile_kd` |
| locus/arm/family CN builder | `mirna_locus_cnv` |
