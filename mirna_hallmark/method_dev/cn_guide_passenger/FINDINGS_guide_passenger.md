# FINDINGS — guide/passenger CN reduced-form contrast (powered up), 2026-07-18

## PREDICTION OUTCOME (axiom 1a) — I WAS WRONG.
I predicted the contrast would DIE under proper hairpin clustering + a matched null (predicted p 0.05–0.30,
"underpowered not refuted"). It SURVIVED. My "change-my-mind" criterion (survives hairpin clustering AND
beats a matched null at p<0.05) was MET. ⇒ from here every downstream claim gets extra checks, not arguments.

## THE DESIGN (why it is confound-immune)
Both arms of a hairpin share the IDENTICAL locus CN ⇒ the reduced form ρ(CN,target|C) is identical for the
two arms for a given gene. The contrast τ_G − τ_P is identified ONLY by which arm's SEED targets which gene
(sequence = exogenous to CN). Ploidy AND co-amplification difference out. C=[CPE,target_cn,mal_prolif,PLOIDY]
(derived ploidy = mean ASCAT3 CN over 506 hairpin loci; LOO negligible, corr≥0.99997, verified by me).
Estimator: reduced ~ α_arm + γ_gene + τ_G·site_G + τ_P·site_P. Site = TargetScan context++; baseline =
genuinely site-free (no TargetScan AND no scanMiR duplex).

## RESULTS (all re-derived by me; role-swap null B=2000 is the primary inference)
| spec | n | gap τ_G−τ_P | role-swap null p | note |
|---|---|---|---|---|
| unweighted ALL | 126,384 | −0.01302 | **0.0015** | hairpin+gene clustered analytic p=0.0096 |
| trans-only (diff chrom) | 120,266 | −0.01273 | **0.0035** | co-amp impossible — SURVIVES |
| cis-only (same chrom) | 6,118 | −0.01575 | 0.037 | small n, noisier |
| confident guide |Δab|≥2, unwtd | 76,588 | −0.02337 | **0.0005** | clearest pairs = strongest |
| trans + confident | 72,869 | −0.02400 | **0.0005** | |
| F-weighted | 126,384 | −0.01041 | 0.083 | weaker — but F-wt down-weights passengers (low-F by defn) |

Leverage: leave-one-hairpin-out gap range [−0.0110, −0.0159] (159 hairpins) ⇒ NOT one oncomiR cluster.

## DECOMPOSITION + the honest wrinkle
τ_G ≈ −0.006 (guide targets repressed), τ_P ≈ +0.006 to +0.018 (passenger POSITIVE).
- GUIDE-only (trans): τ=−0.00346, hairpin+gene-clustered analytic p=0.197 (n.s.); site-perm p=0.0005
  BUT the site-perm null is ANTI-CONSERVATIVE (shuffles sites across genes ⇒ destroys the gene-level
  confound structure) — do NOT trust it. Analytic (clustered) is the honest one ⇒ guide-alone n.s.
- PASSENGER-only (trans): τ=+0.00593, n.s. (analytic p=0.84); the +sign is CONSISTENT across specs and
  is NOT mechanistically explained (plausibly the ~6× fewer passenger predicted-sites are noise, but a
  noisy 0 should not be systematically +). FLAGGED as not-understood.
- WHY the contrast is significant while neither arm is: the within-hairpin role-swap DIFFERENCES OUT the
  shared gene-level confounds (same genes, same hairpins, same CN) ⇒ the contrast is identified far more
  precisely than either arm alone (null sd 0.0044). This is the intended behaviour of the design, not a
  paradox — but it does mean the signal is a RELATIVE (guide-vs-passenger) claim, not "guide represses"
  as an absolute.

## RECORD CORRECTION (independent of the positive result)
Reproduced MH-133's site-type efficacy ladder exactly: 8mer −0.00174 / 7mer-m8 −0.00007 / 7mer-A1 −0.00057,
ladder diff −0.00117, p=0.256, per-bin SE ≈0.0025 — NOT ONE bin distinguishable from 0. The ladder is
UNDERPOWERED (MDE on the gradient ≈0.0044 > the whole site effect), NOT a refutation. MH-133/129's
"genuine REFUTATION, not power failure" language for the ladder + the F>10-gated site-free spec is
OVER-STATED. (The two NAMED estimators — pi_causal mediation, within-arm-FE headline — remain dead.)

## VERDICT (to be gated by rigor-auditor)
The guide/passenger within-hairpin contrast is a SMALL, honestly-bounded, but ROBUST exogenous signal that
curated guide-arm targets are more CN-repressed than passenger-arm targets — surviving hairpin clustering,
a role-swap permutation null, trans-only (co-amp-free), and leave-one-hairpin-out, and strengthening on the
clearest pairs. It is the FIRST CN-based exogenous signal in the program that survives proper inference.
It does NOT resurrect per-edge existence (effect ~−0.013, ~5× attenuated from structural; relative not
absolute) and carries an UNEXPLAINED positive-passenger term. Next: (a) rigor-auditor; (b) if it passes,
CCLE cross-line replication (composition-free second cohort) is the natural confirm.

## QC / PROFILING (user-requested, 2026-07-18) — all TRANS, role-swap null B=2000
- QC-1 both-arms-target: only 57 (hp,gene) both-targeted (11,305 guide-only / 1,829 passenger-only).
  Excluding both-targeted genes: gap −0.01268, p=0.0035 — UNCHANGED. Not a driver. ✓
- QC-2 guide-definition: signal is ABUNDANCE-SPECIFIC. guide=abundance −0.0127 p=0.0015; guide=VARIANCE
  +0.0007 p=0.57 (NULL); guide=detectability −0.0023 p=0.30. abundance↔variance concordance 0.72.
  ⚠ PUZZLE: a clean CN-repression signal should track the axis that lets CN ACT (variance), not abundance.
  Cuts against the "repression" reading; consistent with "functional capacity (abundance) gates it."
- QC-3 dose-response on |ΔabundE|: [0,1) +0.004 p=0.62 → [1,2) −0.010 p=0.22 → [2,4) −0.022 p=0.008 →
  [4,∞) −0.025 p=0.0005. Clean MONOTONE gradient (no effect when arms are similar; strong when guide
  dominates). Strengthens the abundance-asymmetry interpretation. Positive passenger grows with |ΔabundE|
  too (tP=+0.019 at |ΔabundE|≥4) — still unexplained.
- NET: strengthens "real abundance-asymmetry effect" (dose-response, both-targeting ruled out); complicates
  the causal-repression reading (variance-null). Stays PROVISIONAL; CCLE confirm more important.

## CONTEXT-AWARE READOUTS (user-requested, 2026-07-18) — TRANS, role-swap null B=2000
### Readout A — guide's DOMINANCE in the gene's regulatory profile (rank by abundance among the gene's
detectable site-bearing regulators). guide is the gene's TOP regulator in only 7% of edges (median share 0.7%).
| guide rank in gene | n_guide_sites | tau_G | gap | p |
|---|---|---|---|---|
| TOP (rank 1) | 787 | −0.0132 | −0.0203 | 0.013 |
| rank 2 | 772 | −0.0167 | −0.0239 | 0.0045 |
| MINOR (rank≥3) | 8,063 | −0.0073 | −0.0142 | 0.0030 |
⭐ The GUIDE-side effect (τ_G) SCALES with dominance (~2× stronger for a top-2 regulator vs minor), while
τ_P is FLAT (~+0.007) across strata ⇒ the dominance stratification ISOLATES the guide biology from the
dominance-invariant positive-passenger offset. Consistent with "a dominant functional guide represses its
targets"; also explains the QC-2 variance-null (dominance is an abundance-LEVEL property, not a variance one).

### Readout B — family-level aggregation: STRUCTURALLY POWER-LIMITED, does not add power.
The guide/passenger contrast is a PER-HAIRPIN construct (the two arms of a hairpin are different seed
families). Only ~7 guide-families contain ≥3 both-arm hairpins; per-family contrasts are barely identified
(only 1 family met ≥2 hp AND ≥5 guide&passenger sites, gap≈0 degenerate). ⇒ family aggregation is not a
viable power gain here; the productive enrichment is the per-edge dominance lens (Readout A).

## UPDATED CHARACTERIZATION (post-enrichment)
The finding is best stated as: within a hairpin, the DOMINANT (high-abundance) guide arm's TargetScan
targets show more negative locus-CN↔expression coupling than the passenger arm's targets, and this
guide-side effect SCALES with how dominant the guide is in the gene's own regulatory profile — a graded,
abundance-driven, exogenous (ploidy- & co-amp-differenced) RELATIVE signal. Still not per-edge existence;
still carries a dominance-invariant positive-passenger offset (unexplained); CCLE is the independent confirm.

## ⛔ CO-CLUSTER LEAK — THE GUIDE-REPRESSION READING IS LARGELY FALSIFIED (user Q2, 2026-07-18)
The within-hairpin differencing removes the SHARED CN symmetrically, but a co-cluster (polycistron) member
that shares the same CN and CO-TARGETS the guide's genes more than the passenger's leaks ASYMMETRICALLY into
tau_G. Measured:
- Guide target genes are co-targeted by a detectable co-cluster member **43.2%** of the time vs passenger
  targets **16.6%** — a large asymmetry, exactly the confounding channel.
- **SOLO hairpins (80/159; no detectable co-cluster arm ⇒ leak impossible): tau_G = +0.0014 (GONE, was
  −0.0058), gap −0.0079, p=0.089 (n.s.).** No guide repression on clean hairpins.
- **CLUSTERED hairpins (79/159): tau_G = −0.0113, gap −0.0160, p=0.005.** The signal lives ENTIRELY here.
- Excluding co-cluster-co-targeted genes: tau_G collapses −0.0058 → −0.0019.
⇒ The guide-side "repression" signal is SUBSTANTIALLY CO-CLUSTER LEAK (polycistron neighbours co-repressing
shared targets through the shared cluster CN), NOT the guide's own CN-driven repression. On clean solo
hairpins there is NO guide effect. The residual gap is the still-unexplained positive passenger.
⚠ Readout A's dominance gradient is likely co-confounded (dominant guides = oncomiR-cluster members).

## STRONGEST CASES + CANONICAL CONTROLS (face validity of the reduced form itself)
Top guide arms by mean CN reduced form are canonical breast regulators: miR-30d/b (−0.093), miR-16/15b
(−0.075), miR-200c (−0.068), miR-93 (−0.062), miR-182/96, miR-21 — ALL cluster members. Canonical edges:
miR-16→BCL2 −0.227, miR-16→CCND1 −0.184, miR-93→PTEN −0.152, miR-200c→ZEB1 −0.130 (strong negatives ✓);
but miR-155→SOCS1 +0.066, miR-17→CDKN1A +0.059 (positive). miR-93(guide)→PTEN and miR-106b(passenger)→PTEN
share the IDENTICAL reduced form −0.152 (same cluster CN) — the co-cluster situation in one canonical edge.
⇒ The CN reduced form DOES register strong canonical edges (face validity), but cannot attribute them below
the CLUSTER — consistent with the whole CN_INSTRUMENT §3 unit-hierarchy doctrine.

## ⭐ FINAL VERDICT (revised after the co-cluster control)
The guide/passenger contrast does NOT survive the co-cluster control as guide-repression: it is largely
polycistron co-targeting leak (clean solo-hairpin contrast is n.s., tau_G≈0). The CN instrument is NOT
revived by this design. What stands: (1) the RECORD CORRECTION (ladder underpowered; MH-129/133 over-stated)
— independent and intact; (2) face validity that CN detects strong canonical edges but only at CLUSTER
resolution (per §3). Both NAMED instruments stay dead; per-edge existence stays dead. No new positive
registry finding is warranted for the guide/passenger contrast.

## ⚠ CORRECTION to the "co-cluster FALSIFIED" section above (abundance-aware, user-caught 2026-07-18)
The solo-vs-clustered comparison is CONFOUNDED WITH ABUNDANCE: solo guides median 4.97 log2RPM vs clustered
7.37 (~5×), |dAB| 2.50 vs 3.80. So "solo tau_G≈0, n.s." does NOT by itself prove co-cluster leak — solo
guides are simply lower-abundance/less-dominant, and the effect is abundance-graded. The CLEAN test holds the
arm set fixed (clustered hairpins, arm FE absorbs abundance level) and splits guide sites by co-targeting:
- guide, CO-targeted by a co-cluster member: tau=−0.0147 (p=0.002, n=4,914)
- guide, NON-co-targeted (its OWN targets):   tau=−0.0062 (p=0.116, n=2,507)  <- SAME magnitude as baseline tau_G
- passenger: +0.0044 (p=0.62)
⇒ CORRECTED VERDICT: co-cluster leak is REAL and PARTIAL (~2.4× inflation on co-targeted genes ⇒ ~half the
apparent guide effect is polycistron co-targeting), but the guide's OWN effect does NOT vanish — it is
−0.0062 on non-co-targeted genes, marginally n.s. (p=0.116). NOT "falsified"; rather the signal is SMALL,
partly co-cluster-inflated, and its clean own-component is not individually significant. This REINFORCES the
rigor-auditor's PROVISIONAL/relative verdict and adds a demonstrated partial confound. Still not a clean CN
revival; still needs an independent cohort (CCLE) to say anything stronger.

## ⭐ DECISIVE SHARPENING — the clean guide-own effect FAILS the site-efficacy ladder (2026-07-18)
Readout = one number per (guide,gene): partial-Spearman(locus CN | C, gene expr | C) over ~1065 tumors; same
CN, different gene; arm FE absorbs the guide's mean, tau = site vs site-free within-arm.
On the CLEAN guide-own subset (non-co-targeted, leak-free, n=6,448 sites):
- SITE-TYPE ladder: 8mer −0.0002 / 7mer-m8 −0.0027 / 7mer-A1 −0.0026 — NON-monotone (8mer weakest, backwards).
- CONTINUOUS context++ slope: beta=−0.00074, p=0.53 (slightly WRONG sign).
⇒ The clean guide-own effect does NOT scale with site strength — the SAME ladder that killed the within-arm-FE
instrument (MH-133). Either not repression, or the CN reduced form is too weak to resolve site strength; either
way NO site-discriminating exogenous signal.
Passenger-positive EXPLAINED: passenger arms ~32× lower abundance (2.95 vs 7.97 log2RPM); tau +0.005..+0.008,
site-type-independent, never significant = LOW-ABUNDANCE NOISE offset, not anti-repression.

## ⭐⭐ CONVERGED FINAL VERDICT (CN revival, fully investigated)
The guide/passenger within-hairpin contrast does NOT revive the CN instrument. Fully decomposed it is:
(a) ~half CO-CLUSTER LEAK (polycistron neighbours co-repressing shared targets via shared CN, ~2.4× on
co-targeted genes); (b) a LOW-ABUNDANCE PASSENGER NOISE offset (+0.007); (c) a small guide-own component
(−0.006, p=0.12) that FAILS the site-efficacy ladder ⇒ not demonstrably site-mediated repression.
No new positive registry finding is warranted. The two named instruments stay dead; per-edge existence stays
dead; edge existence still rests on one observational line. CN reduced form registers strong canonical edges
(face validity) but only at CLUSTER resolution (CN_INSTRUMENT §3).
WHAT STANDS (clean, independent, rigor-auditor PASS): the RECORD CORRECTION — MH-133/129 over-stated the
retraction ("genuine refutation" of an underpowered ladder); downgrade to "underpowered / can't-tell".
METHOD LESSON: always profile abundance/variance before interpreting a stratified contrast (the solo/clustered
abundance confound); most apparent structure here was abundance- or cluster-driven.

## ⭐⭐⭐ DEFINITIVE CLOSE — the contrast is an ABUNDANCE confound (user's matched-design idea, 2026-07-18)
Tried to make the comparison fair (comparable gene-sets + comparable HIGH position). Two decisive facts:
1. The passenger is STRUCTURALLY never high-position: guide sites rank1/2 = 787/772; passenger rank1/2 = 1/14
   (passenger = low-abundance arm ⇒ almost never a dominant regulator). Matched-high-position is IMPOSSIBLE
   (15 qualifying passenger sites). Not a data gap — intrinsic to "passenger".
2. COMPARABLE-ABUNDANCE hairpins ⇒ contrast VANISHES: |dAB|<1.0 gap +0.0043 p=0.65; <1.5 +0.0023; <2.0
   +0.0006 p=0.53. The signal exists ONLY at large |dAB| and is gone (even +) when arms are comparable.
⇒ The contrast is an ABUNDANCE-ASYMMETRY effect ("the higher-abundance arm's targets anti-correlate with the
locus CN more"), NOT seed-mediated repression. ⚠ This OVERTURNS the earlier read of the |dAB| dose-response
as a strength — it is the CONFOUND's fingerprint (a seed effect would NOT collapse at comparable abundance).
DEFINITIVE: guide/passenger is foreclosed for CN edge-existence (CN can't separate arms, §3) AND empirically
an abundance confound. No positive finding. Record correction (MH-133/129 over-stated) is the clean output.
