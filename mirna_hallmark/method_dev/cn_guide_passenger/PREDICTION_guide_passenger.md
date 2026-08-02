# PREDICTION — guide/passenger powered-up test (written BEFORE running; axiom 1a)

Date: 2026-07-18. Author: this session.

## The design (confirmed by reading cn_falsify.py F3)
- Hairpins with BOTH arms present; guide = higher-median-abundance arm, passenger = lower.
- Both arms share the IDENTICAL locus CN column ⇒ reduced form CN→Y is identical for the two arms for a
  given gene. The contrast tau_G − tau_P is identified ONLY by which arm's seed targets which gene
  (sequence = exogenous to CN). This is a genuine "asymmetry the confound is blind to" (axiom 4).
- Estimator: reduced ~ alpha_m + gamma_g + tau_G*ts_G + tau_P*ts_P; C=[CPE,target_cn,mal_prolif,PLOIDY].
- Prior result (agent-reported, NOT yet reproduced by me): unweighted diff −0.013, one-sided p=0.010;
  F-weighted p=0.17; confident-guide p=0.29. Fragile (1/3 specs).

## My priors going in
- The pooled reduced form is NULL (MH-126 locus-unit τ=−0.0015 p=0.67).
- Memory ⛔ exogenous-identification-or-nothing: CN existence rests on ONE observational line; both named
  instruments dead. Strong prior that CN carries little/no exogenous site-mediated signal in bulk TCGA.
- p=0.010 in 1 of 3 specs, dying under F-weighting, is axiom-5 territory (near-sig, spec-dependent).

## PREDICTIONS (falsifiable)
1. Reproduction: I will reproduce the F3 unweighted one-sided p ≈ 0.010 ± 0.01 and F-weighted p ≈ 0.17.
2. Under PROPER hairpin-level clustering (the guide/passenger rows are perfectly CN-dependent within a
   hairpin), the contrast SE will INFLATE and the unweighted p will WEAKEN to > 0.05 (likely 0.05–0.20).
3. Under a MATCHED NULL (shuffle which arm-of-the-pair is "guide" / permute site assignment within
   hairpin×gene preserving abundance+CN), the observed gap will fall inside the null 90% CI ⇒ NOT
   significant. Predicted null-based p ~ 0.10–0.30.
4. On the STRUCTURAL (Wald) scale (reduced / first-stage per arm) the direction holds (tau_G < tau_P) but
   the SE grows (division by a small, noisy first stage) ⇒ p WEAKER than the reduced-form scale.

## Bottom line prediction
The guide/passenger signal does NOT survive proper clustering + a matched null as p<0.05. It is a real
DIRECTION (guide represses more than passenger — biologically expected) but UNDERPOWERED, not a resolved
exogenous existence leg. Verdict I expect to write: "directionally consistent, does not reach significance
under honest inference ⇒ CN still cannot validate edges; retraction stands, but as 'underpowered' not
'refuted'."

## What would change my mind (a genuine surprise / real lead)
If (2) AND (3) both hold p<0.05 — i.e. the gap survives hairpin clustering AND beats a matched
site-permutation null — then guide/passenger IS the program's first surviving exogenous existence leg,
and it should go to CCLE for a second-cohort check. I will NOT soften this criterion after seeing results.
