# State of play — `mirna_hallmark`

> **Goal:** the ONE doc you read to know where the project stands. Per axis: what **stands**, what is
> **dead**, what is **open** — each with its MH id and date, so you can go deeper without guessing.
> **What belongs here:** one short block per axis, verdicts only. **NOT** methods (→ `FORMULAS.md`,
> `LEARNED_MODEL_METHODS.md`), **NOT** the finding-by-finding record (→ `DISCOVERY_REGISTRY.md`),
> **NOT** run history (→ `ANALYSIS_RUN_LEDGER.md`), **NOT** forward plans (→ `PROGRAM_FORWARD_BOARD.md`).
> **Update trigger:** any arc that changes a verdict. If you write an `MH-` registry row that supersedes
> or retracts something, update the matching block here **in the same pass**.
> **Sync-partner:** `DISCOVERY_REGISTRY.md` (the record of record — this doc is its executive summary).

**Last updated: 2026-08-03** · covering the registry through **MH-207**.
⚠ This header read *"through MH-204"* while the registry already carried **MH-205** — the same
registry-runs-ahead lag the "How to read this" section warns about, caught 2026-08-03.

---

> ⭐ **The PAIRED / per-patient axis ran in full 2026-08-04 (§J, MH-228…MH-239). Before asking or building ANY per-patient question, read [PATIENT_QUESTION_TAXONOMY.md](PATIENT_QUESTION_TAXONOMY.md)** — it defines the P-across vs P-each rung split, which measures are orthogonal (dose vs retention, ρ=−0.081), the TWO unrelated senses of "retention", why the site-free decoy is **not** a null set (p<0.05 at 9.2%), and the four things not to rebuild on this cohort.

## How to read this

**The `DISCOVERY_REGISTRY.md` is the source of truth.** It is append-only and current. Every other
doc — including the ones labelled "canonical" — may lag it. Where this doc and any other doc
disagree, the registry wins, then this doc, then everything else.

⚠ **The registry runs ahead of the docs.** MH-131…153 (the most consequential fortnight in the
project) exist **only as registry rows**. There is no `MH13x_*.md` for them, and that is now policy —
the per-finding-doc pattern is **banned** (`../CLAUDE.md` §3). Do not conclude from a doc's silence
that nothing happened.

---

## The one-paragraph version

The last two weeks were a **self-falsification campaign**, and it worked. Most of the June headlines
are dead. What survives is narrower, better-founded, and honestly bounded. The model is a good
**coupling estimator** and is **not** an attribution engine. Edge existence rests on **one
observational line** — both copy-number instruments are retracted. The curated edge set beats a
properly-built fake by **≈0.012–0.013 ρ** — a value two independent designs now agree on (MH-137/139),
and which *shrank* at every control improvement before converging there. The
single strongest result in the program is a **negative control that correctly reports zero**
(MH-136). Two whole channels — state and protein — were measured and cancelled, for the same
structural reason: information, not effort.

---

## Axis 1 — The learned model (what it *is*)

**STANDS.** One **dense learned-τ² non-negative Bayesian posterior per gene, two readouts** —
π≡1 dense → coupling/attribution/identifiability; evidence-π → discovery. Half-normal slab
(`θ_m ~ N⁺(0,ν²)`) encodes "miRNAs repress". **ν² is learned, not guessed** — that, not L2, is the
whole edge over the lasso (`nnridge_cv ≈ lasso`; OOF −0.168 vs −0.152, p=9e-16). **Adaptive lasso is
RETIRED to a baseline.** Canonical: `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` §6b (2026-07-13).

**Engine: Gibbs, not HMC** — measured, not aesthetic. Every live channel is Gaussian-conjugate;
Student-t is a scale mixture; the protein discordance link is *additive* ⇒ conjugate; HMC was only
warranted for the A3 field, which the A4 probe gated out (pooling is a wash, Δρ −0.0006). Gibbs is
also the better estimator (split-half ρ=0.822 vs bagged-NNLS 0.729). **Carve-out:** fit the *gauge*
with NNLS — Gibbs's heavy-tailed posterior SDs break the errors-in-variables correction.

**What it is FOR (MH-126b):** a **de-confounded, jointly-apportioned, transferable coupling
estimator**. **What it is NOT:** an attribution engine. *It must not be sold as an answer to "which
miRNA represses this gene".*

**OPEN / known limits:**
- **Posterior widths are TOO NARROW — but only ONE of the three numbers is stable (RE-DERIVED with
  error bars 2026-08-01, MH-185).** Measured over **4 seeds** (200 genes, 3 splits), ratio =
  reported_se / true_sampling_SD from INDEPENDENT half-cohorts (CALIBRATED = 0.9–1.1):
  **bagged NNLS 0.745 ± 0.009** (inflation 1.34 — confirms the recorded 1.37×) ·
  **Gibbs Gaussian 0.815 ± 0.084** (1.23×; range 0.74–0.93) ·
  **Student-t ν=7 0.849 ± 0.051** (1.18×).
  ⚠ **The Gibbs figures are IMPRECISE, not just different** — the seed SD is 0.084, which is why two
  runs of the same harness flipped the verdict OVERCONFIDENT↔CALIBRATED. Quote them with the spread;
  the old point estimates (1.29× / 1.13×) were unstable statistics, not wrong readings.
  ⛔ **"ν=7 substantially FIXES the calibration (0.77→0.89)" is TOO STRONG.** It improves the mean
  (0.815→0.849) **and halves the seed variance** (0.084→0.051), so its direction is right — but it does
  **not reach the calibrated band**, and 0.89 sits above the observed spread. Do not turn it on expecting
  a fix.
  ⛔ **A SUBTYPE RANDOM EFFECT WILL NOT FIX IT EITHER — tested, refuted (MH-185).** Subtype-STRATIFIED
  half-splits did **not** shrink the true sampling SD (Δ ratio **−0.026 ± 0.051**, p=0.63; true SD rose
  0.01742→0.01822) ⇒ the unmodelled heterogeneity is **finer-grained than PAM50**, so the fix is
  **REPORTING, not modelling**.
  ✅ **SHIPPED: honest widths are now ON THE EDGE CARD** — `cal_beta_sd` / `cal_z` / `cal_identified`
  (+`_lo`/`_hi` spanning the constant's own ±1 SD). **Identifiability falls 24.8% → 19.9% [18.1–22.4%]**
  on the 5,648 arm-edge card. ⚠ that is the ARM-EDGE rate; the **27.8% of 5,117 UNITS** quoted below is a
  different universe — do not compare them directly.
  ⚠ **`LEARNED_MODEL_STATE_CHANNEL_PLAN.md` §11(a) still says "~40% (NNLS 28%, Gibbs 39%)" and "SBC
  promoted" — SUPERSEDED.** The earlier figure was measured on a biased gene subset (the
  scale-dependent sd-floor bug); **SBC is RETIRED as the wrong tool.**
  ⭐ **MH-92's `nu=7` substantially FIXES the calibration (0.77 → 0.89)** — so Student-t is
  "robustify-not-lever" *for coupling*, but it **is** the calibration fix. Turn it on for any power calc.
  ⚠ **Two things do NOT inherit the inflation:** **Shapley shares are RATIOS, so a common β inflation
  cancels entirely** (re-derived: MH-94's width is **1.8× too WIDE**, not too narrow) — and note
  `readouts.py` does not pass `nu`, so the shipped table is Gaussian.
- **Identifiability ceiling:** median posterior SD/|β| = **0.799**; only **27.8% of 5,117 units
  identified** (|z|>2). Not a fitting failure — same-seed arms share the binding site.
- **The n≈1000 ceiling:** every internal lever (CN channel, δ-pooling, isomiR refit, t-likelihood,
  cross-gene pooling) lands "immaterial at n≈1000".
- ✅ **E6 "self-partialling bug (live)" — NOT LIVE; the tag is RETRACTED (MH-142, 2026-07-17).**
  Verified: every call site scores **orphans**, which are disjoint from `he_agg` by construction —
  **0 of dossier's 6,744 candidate pairs are HE edges** (vs 5,938 HE edges). Nothing is partialled on
  a covariate containing itself in production. E6 **stands as a trap for CLASS-MATCHED comparisons**
  (its evidence column: *"class-matched estimator comparison"*) — push HE arms through this estimator,
  as MH-124's own experiment did, and you get a ~36× depressed HE rate **as an artifact of your design**.
  ⛔ **Do not "fix" it by removing `he_agg`** — that term is what makes an orphan's coupling *partial on
  what curation already captures*, i.e. the point of the lane. Guard: `eval/_e6_production_check.py`.
- The **shipped readouts do not include the CN channel** (`readouts.py` calls `_gibbs_posterior`
  without `channels=`) — the production table is mRNA-likelihood + evidence-prior only.

**FIXED, don't re-litigate:** the `_rtnorm_pos` **sampler bug** (inverse-CDF saturated for
`mu/s < −7.0345` ⇒ deterministic *negative* draws, impossible under a half-normal slab; contaminated
3.15% of β). Replaced with Robert (1995) rejection sampling, KS-verified. ⚠ It **overturns MH-119's
recorded cause** — the `share`=43.7 explosion was logged as "anti-repressive biology"; it was this bug.

---

## Axis 2 — Edge existence

**⛔ It rests on ONE observational line:** *curated edges anti-correlate more than abundance-matched
site-free ones.* The "three independent validations" headline is **FALSE** — the other two legs were
both observational and never independent of each other.

**Both CN instruments are RETRACTED:**
- `pi_causal = γ·b` where **b is an observational OLS slope** ⇒ a product-of-coefficients mediation
  estimator, not an IV. It re-derived the anti-correlation it was meant to validate. γ — the only
  exogenous factor — is **site-blind** (HE +0.19931 vs decoy +0.19844, p=0.20). The clean reduced
  form, arm-clustered: **p=0.115, n.s.** (MH-124r/126, 2026-07-13). The retraction is in the
  production code at `learned/instrument.py:1128–1148`.
- The within-arm two-way-FE replacement: its control class was **71% real binders**. Against a
  genuinely site-free control **τ=−0.0007, p=0.84**, and it fails the site-type efficacy ladder
  (8mer/7mer-m8/7mer-A1 not monotone, p=0.26) — copy number cannot see a site type, so genuine
  site-mediated repression *must* follow it (MH-133, 2026-07-16). ⚠ **LANGUAGE CORRECTED (MH-157,
  2026-07-18): these are UNDERPOWERED / can't-tell, NOT "refutations."** The ladder's MDE@80% (0.00442)
  exceeds the whole 8mer effect (0.00174); no bin is distinguishable from 0. The "τ=−0.0007, p=0.84" was
  the F>10-**gated** spec — F-weighted on the same control gives p=0.057. **The two named instruments stay
  dead regardless** (their retraction is unaffected); only the framing is downgraded.

**⛔ AND A THIRD DESIGN — guide/passenger within-hairpin — IS ALSO FORECLOSED (MH-157, 2026-07-18).** The
one confound-immune-*looking* CN design (two arms of a hairpin share the identical locus CN ⇒ ploidy &
co-amplification difference out; contrast identified only by which seed targets which gene) does **not**
revive the instrument. It is (a) **conceptually foreclosed** — CN cannot separate hairpin arms (`CN_INSTRUMENT
§3`), so the contrast is a *gene-set* comparison, not edge identification; and (b) **empirically an abundance
confound** — it vanishes when arms are comparable-abundance (|ΔabundE|<1 gap +0.004, p=0.65) and is ~half
co-cluster leak (guide targets co-targeted by a polycistron neighbour 43% vs passenger 16%); the clean
guide-own residual (−0.006, p=0.12) fails the site-efficacy ladder (context++ slope p=0.53). **Edge existence
still has NO surviving exogenous validation.**

**⚠ Also bounding this axis:** the per-edge null is **3–4× too narrow**. Site-free (cannot-repress)
pairs pass "FDR q<0.05" at **25–35%**; under an honest empirical null **0.0%** of curated HE edges
survive per-edge (was a nominal 31.3%). **This does not say the edges are fake** — a typical HE edge
is ρ≈−0.15 against σ₀≈0.09. It says the signal is a **set-level distributional shift, not per-edge
significance** (MH-123/124). Any FDR count in any doc predating this rests on the uncalibrated null.
**⛔ And a permutation null does NOT rescue it (MH-154):** `coupling_permutation`'s 2000 Freedman–Lane
draws were measured to reproduce the **theoretical** null (σ₀=0.0309 vs 0.0311, R²=0.9997) — **permuting
destroys the very structure that inflates the real null.** Nor does **BY**, which corrects *dependence*,
never *scale*: MH-45's BY count (**33.3%**) sits inside the **35.3%** rate for *impossible* edges.
**The calibration standard is `eval/site_free_null.py`** (rescued from `/tmp`, reproduces MH-123 exactly).
It is **conservative** — the site-free class holds some real non-canonical targets.

**⭐ THE DISCOVERY LANE, RESOLVED (MH-155/156).** Per-edge and per-family discovery in bulk are **EMPTY**
under the honest FDR (empirical heavy-tailed null → Simes within family → BH across families; `q_family_emp<0.05
= 0` everywhere). BY was demoted to a worst-case sensitivity bound — it corrected the wrong defect (dependence,
not the heavy tail; it even kept 2 tail false-positives the empirical FDR caught). **The defensible claim is
CONVERGENT-EVIDENCE, not per-edge:** coupling co-varies with independent biochemistry, abundance-controlled, on
multiple axes — **site-level chimeric overlap partial ρ=−0.096 (p=2.7e−15), scanMiR K_D −0.091 (p=7e−14),
site-confidence rung −0.050 (p=3.9e−5); site COUNT n.s.** (quality not quantity) — and the site–duplex
coincidence climbs the site-confidence ladder monotonically (9%→27%). Deliverable = **A1∩site-level-chimeric,
157 edges** (miR-18a/miR-17~92 cluster). Effects are SMALL by construction (bulk is weak per-edge); the strength
is the convergence. **⚠ NOT a per-edge licence.** Home: `learned/discovery.py`, `method_dev/site_ladder/`.

⭐ **QUANTIFIED + GIVEN A PRODUCER 2026-08-02 (MH-198)** — `eval/discovery_gold_set.py` →
`output/learned/discovery_gold_set.tsv` (sha256-stamped). It had existed only as a hand-merge of two files,
the same shape that rotted the literature sets (MH-196). **157 edges · 90 genes · 18 arms · ⭐ 11 SEED
FAMILIES**, and the concentration is now measured rather than gestured at: **96/157 (61%) are ONE family
(miR-17/20/93/106/519) and 139/157 (88%) are the miR-17~92 polycistron + its paralogue clusters
(miR-106b~25, miR-106a~363), which are CO-TRANSCRIBED.** ⇒ **the honest unit is ~one oncogenic cluster's
realized target set, not 157 independent findings — always quote the 11 families beside the 157 edges.**
✅ What is strong: **94/157 sit on genes with NO curated regulator at all**, and only **3/157** are same-seed
paralogues of a regulator already curated for that gene ⇒ not curation leaking back in.
⛔ **`KIF3B` is NOT in the site-level gold set** (edge-level CLIP only) — the handoff's
miR-18a→{STAM2, KIF3B, MAP3K1, NEDD4} quartet is **edge-level; at site level it is a TRIO.**
⚠ **DO NOT compress this axis to "discovery is empty"** — that is true ONLY of the per-edge FDR and
misstates the lane's actual deliverable.
✅ **FRAMING DECIDED 2026-08-02 (user: "just write both") — BOTH RUNGS WRITTEN, MH-199 + MH-200.**
**MH-199 (cluster rung, BIOLOGICAL):** the cluster's share climbs 25.3% → 66.4% → 88.5% through the filters,
and ⭐ **the enrichment SURVIVES both confounds I expected to kill it** — conditioned on Manakov detection
cluster A1-rate **20.1% vs 2.5% (OR=9.76, p=3.5e−28)**; abundance-quintile matched **54.5% vs 13.2%**, gap in
all 5 quintiles. ⚠⚠ **BUT QUOTE THE OR, NEVER THE 88%** — the `site_manakov` rate is 40.1% vs 14.1% and stays
elevated at every abundance quintile ⇒ the final share is enrichment **plus HEK293T assay-coverage bias.**
Biological reading is *"a canonical oncomiR polycistron is UNDER-CURATED relative to its realized activity"*,
**not** a novel regulator.
**MH-200 (edge rung, CANDIDATE QUEUE):** the 157 edges are a **ranked validation queue**, not 157 claims —
per-edge FDR is empty by construction, and the set is 11 families. ⚠ **No breast chimeric exists in any
source** (Manakov=HEK293T, TarBase CLASH=other tissues) ⇒ hurdle 3 is cross-tissue corroboration, never
breast-context validation.

**OPEN:** restoring an exogenous existence validation is, by both docs' own assessment, **the single
highest-value open item in the program** — and it is now **TESTED, not aspirational (MH-159, 2026-07-18).**
**⛔ GERMLINE cis-eQTL MR built + run + refuted.** Germline genotype is genuinely exogenous (randomized at
meiosis; no reverse causation / tumour-state confound), so a germline cis-eQTL instruments a miRNA's dose.
Built at n=1,075 (near-complete overlap): SNP6 birdseed → hg38 cis-SNP matrix (1,075×128k). But the
instruments are **honestly WEAK** (cross-fitted OOF: only **10 arms F>10**; the single-SNP "38" was
winner's-curse inflation) and **NON-SPECIFIC** — the reduced form is repression-directed for the miRNA's
curated targets **no more than for random non-target genes** (curated +0.0027 vs control +0.0049, paired
p=0.66; still null with the strong multi-SNP instruments, p=0.42). ⇒ **germline fixes EXOGENEITY but not
EXCLUSION** (cis-SNPs still hit host genes/clusters/paralogs) and the miRNA cis-heritability is small; SNP
aggregation doesn't rescue it. **Edge existence still rests on ONE observational line.** ⬜ Only untested
extension: full genotype imputation (unlikely to change it — low OOF-F + specificity null). Home:
`method_dev/cn_guide_passenger/GERMLINE_FINDINGS.md`.
For the discovery lane specifically: ✅ ~~correlation-matched family null~~ **DONE (MH-197)** — pseudo-families
now match the ρ̄ of the specific candidate cell they stand in for; the old uniform draw was **20–25% too
narrow** (null sd ×1.196, robust-sd ×1.249, Levene p=1.5e−21 on a paired A/B). **Verdict unchanged at 0
survivors — as it had to be, since widening an anti-conservative null cannot create one — but the headline
set-level shift attenuates, median `null_z` −1.6 → −1.32.** ⚠ The `scrambled-seed` variant was scoped and
NOT built: shuffling a seed changes which genes an arm is *predicted* to hit, not its expression, so it
matches site statistics but does **not** remove MH-123's arm-lane leak (a separate defect, opposite
direction — still open). ⚠ Verified 2026-08-03: the scrambled-seed null is **still unbuilt** (no code).
⭐ **NEW 2026-08-03 (MH-207) — THE FAMILY LANE NOW CARRIES EVIDENCE, AND THE ARM LANE'S CONVERGENCE DOES
NOT TRANSFER TO IT.** `discovery.attach_evidence_family` pools each family's member arms (chimeric / ledger /
context++) — coverage on 3,119 families: ledger 76.5%, chimeric 45.3%.
- ⭐ **PMIDs are UNION-DEDUPED, never summed — summing inflates depth 2.55×** (33,012 → 12,935; worst row
  NFAT5/miR-17~92 **159 → 34 = 4.7×**), because same-seed members are routinely assayed in one paper.
- ⚠⚠ **Every pooled statistic is MONOTONE IN FAMILY SIZE** (chimeric support **25.2% at 1 arm → 89.1% at
  7+**), and **42.2% of rows are single-arm families where the family lane IS the arm lane** by construction.
- ⛔ **So the naive result is an artifact.** Pooled MWU p=0.0074 on a **Δ of 0.017** in null_z; conditioned,
  the **size-stratified Stouffer is Z=−0.31, p=0.76 with incoherent signs (4/7; the two largest strata point
  OPPOSITE ways)** and the rank-partial on (size, abundance) is **ρ=−0.0101, p=0.571**. ⇒ **MH-155's
  "chimeric-present edges couple stronger, p=4.6e−8" is an ARM-RUNG claim — do not broadcast it to families**
  (MH-111: the unit at which a property acts predicts which tool can see it).
- ⛔⛔ **AND LEDGER PMID DEPTH IS DEAD TOO — rigor-gated the same session.** The pooled ρ=−0.0578
  (p=0.00124) was **entirely BETWEEN-gene**: the **within-gene** contrast (2,566 families / 663 genes,
  gene-demeaned, size+abundance-controlled) gives **ρ=+0.0032, p=0.889** — zero and nominally the wrong
  sign; the within-gene fame null centres at −0.0010 ± 0.0331 (**empirical p=0.55**); the study-bias channel
  is **larger here than MH-196 measured** (+0.330 mean / +0.500 median, 70.7% of 604 genes); and the internal
  null control **fires the wrong way** (`ev_n_arms_supported` +0.0607, p=0.0081).
⇒ ⭐ **NET: NEITHER evidence axis carries convergence at the FAMILY rung, while the ARM rung does
  (MH-155, p=4.6e−8). The lane's deliverable is the ATTACH + its two design rules, not a finding — and the
  general lesson is to run an evidence-vs-coupling test at the rung the evidence is RECORDED at** (MH-111).
Remaining: unify the two lanes by family size (`family_size_degenerate` now marks the degenerate half) ·
subtype-stratified (PAM50) run — ⚠ but note **the board records MH-165 as having closed the subtype-discovery
hope at chance** (337 vs 287±35, z=+1.4); reconcile before spending on it.

---

## Axis 3 — Attribution / Shapley

**STANDS: the model RANKS the canonical family high; it does NOT NAME it.**
**⭐ UPDATED 2026-08-02 (MH-196) — on a VERSIONED ground truth 10× larger, β is NO LONGER at chance, and
abundance's apparent superiority does NOT survive a fame control. Cite MH-196, not MH-126 Test 1.**

⛔ **First, why MH-126's numbers cannot be re-checked.** Its n=32 set **is not on disk**, and the four
neighbouring literature sets that are (n=16/19/21/23, plus a 173-row power test) **all have zero producers
in the repo** — every one was written 2026-07-13 by a `/tmp` scratchpad since deleted. Five mutually
inconsistent hand-curated lists, no definition, no version. The old row is kept below for history only.

**THE SET IS NOW MECHANICAL AND AUDIT-STAMPED** (`eval/lit_ground_truth.py`): canonical family = argmax of
**distinct low-throughput functional PMIDs** (reporter/western/proteomics, weak excluded) from the
PMID-deduped ledger, admitted only when the gene is in the design, the family is in that gene's design, ≥2
families compete, and the top is unambiguous. **tier2 329 genes / 92 families · tier3 200/63 · tier4 132/48**,
with a sha256 + source-file stamp in `lit_ground_truth_provenance.json`. Positive control: it recovers
PTEN→miR-21, ZEB1→miR-200, BCL2→miR-15/16, CDKN1A→miR-17/20/93 with no hand curation.

**RAW RANK** (0 = top, 0.5 = chance; cluster-bootstrapped on the 92 families):

| readout | MH-126 (n=32, unreproducible) | **MH-196 (n=329 / 92 fams)** | 95% CI | p |
|---|---|---|---|---|
| **abundance** | 0.240 | **0.359** | [0.291, 0.437] | 8e-4 |
| **`identity`** | 0.370 | **0.416** | [0.361, 0.473] | 0.004 |
| **β** | 0.518 — *chance* | **0.436** | [0.382, 0.491] | **0.021** |
| **β_deconv** | — | 0.453 | [0.393, 0.510] | 0.11 n.s. |

⇒ **β is significantly better than chance, and MH-126's 0.518 falls OUTSIDE the new CI.** Raw ordering is
still abundance > identity > β (paired β−abundance +0.077, p=0.005).

**⭐⭐ BUT THE RAW ORDERING IS A STUDY-BIAS ARTEFACT, AND THE FAME NULL SHOWS IT.** The ground truth is
*defined* by PMID depth and abundant miRNAs get studied more — measured here at the family level, within
gene: **spearman(abundance, LT-func depth) = +0.187 mean / +0.244 median, positive in 67.3% of genes,
p=3.2e-10**. Substituting **another gene's canonical family that also sits in this gene's design** holds
fame constant and varies only gene-specificity:

| readout | real | fame-null | **Δ gene-specific** | 95% CI | p |
|---|---|---|---|---|---|
| **β** | 0.439 | 0.525 | **−0.085** | [−0.164, −0.005] | **0.038** ✅ |
| **`identity`** | 0.424 | 0.501 | −0.077 | [−0.156, +0.005] | 0.066 |
| **abundance** | 0.377 | **0.429** | −0.052 | [−0.165, +0.066] | **0.40 ✗** |

⇒ **abundance ranks ANY famous family at 0.429 — its raw win is largely publication bias, and its
gene-specific component is not distinguishable from zero.** β and identity sit at chance on the null
(0.50–0.53), so what they have is gene-specific. Controls fire hard (`lit_depth` Δ=−0.647, `w_ledger`
−0.540), so the null has ample power.

**⚠ THE HEAD-TO-HEAD IS TIED, AND IS UNRESOLVABLE IN PRINCIPLE.** Δ(β)−Δ(abundance) = **−0.033, 95% CI
[−0.165, +0.085]**. ⭐ **This CLOSES the board's "scale the ground-truth set" item with arithmetic rather
than a deferral:** the cluster unit is the canonical family, only **~330 families are the top-cited
regulator of any gene in the entire curated literature**, and exhausting all of them still leaves a
half-width of 0.065 > |−0.033|. Resolving it needs **~1,241 clusters = 3.8× every canonical family ever
published.** ⇒ **the ground truth is NOT the binding constraint; scaling genes cannot fix this.**
(Recoverable headroom does exist — the 480 ambiguous ties would add 262 families — but 426 of them rest on
a single paper, and even taking all of them the contrast stays inside the CI.)

— *historical, MH-126 Test 1, kept because the numbers are cited elsewhere; the set is unreproducible* —
The 32 genes held only **13 distinct families**, and abundance is family-constant ⇒ the independent
unit is the **family**. Cluster-bootstrapped, shapley survived. **Argmax is at chance at every k.**

**MH-124 §4b and MH-126 do NOT conflict** — both say β is at chance; §4b's positive claim was about
`shapley_identity` and dose, which MH-126 reproduces. §4b survives as the **estimand argument**: the
literature's "miR-X represses gene Y" is a **level** claim (knockdown/overexpression); β is a
**covariance** functional. Asking β to recover it is an estimand mismatch, not a model failure.

**WITHDRAWN:** shapley's **evidence gradient** — 85% is family fame (family-preserving null p=0.297).
Only β's gradient survives, barely (b≈+0.031, p=0.074–0.097 under family FE). Shapley wins the
**rank** test and loses the **gradient** test; β vice versa. Dissociation, not inconsistency.

**The w-circularity gate PASSES, bit-identically:** β is unchanged (max|Δβ| = **0.0**, 0/48 genes)
under shuffled *and* constant `w`, with the positive control firing (`w_fam` moved 42/48). Mechanism,
not luck: `X_fam` is pure abundance and the dense Gibbs is called with `π = ones(p)`.
⚠ **But `pip_discovery` and `prior_pi` ARE w-contaminated by construction** — any "do canonical
regulators score higher?" test on those columns is **circular**.

✅ **FIXED 2026-07-17 (MH-140) — the table now has a TRUE identity column, and it disagrees with
magnitude in 1 of every 4 multi-family genes.** `readouts.py` emits **`identity`** = `shapley_identity`
(Shapley/LMG on **R²(Xz·m, yr)** — a NON-additive value function, so collinear overlap is genuinely split),
with **m = bagged NNLS**. The old `share` (= β_f/Σβ) is renamed **`beta_frac`** and demoted to MAGNITUDE:
for the ADDITIVE aggregate Shapley(f) ≡ β_f, so it split nothing (MH-122), and `top_family_identity` was
`argmax β` in disguise — **identical to `top_family_magnitude` in 99.35% of genes**.
**Now: WHO ≠ HOW-MUCH.** Agreement **85.3%** genome-wide, **75.2% among the 819 multi-family genes**.
**The mechanism, at genome scale:** `m_nnls` is exactly 0 in **1,623/5,117** families (31.7%); `beta` is
exactly 0 in **0/5,117** — the half-normal slab cannot zero an un-informed family, NNLS can.
⭐ **Identity can also say UNDEFINED:** 213/1,549 genes (13.8%) have NaN identity, and in **100%** of them
NNLS zeroed every family ⇒ no coupling to distribute. `beta` structurally cannot make that statement.
**Control: purely additive** — `beta`/`beta_sd`/`z`/`pip_dense`/`pip_discovery` bit-identical (max|Δ|=0.0).
⚡ Three independent wins, all controlled: `shapley_identity` **8.2–9.7×** (output-identical); **1 of 4 Gibbs
chains per gene was DEAD** (the evidence-π chain ran on both C blocks; `dec["pip_disc"]` was never read) —
gated, output-identical; and **N_ITER/BURN 2000/700 → 1000/350** (MH-141), licensed against the sampler's own
reseed jitter and verified NOISE-not-BIAS. **Run 5.9 → 1.6 min (3.7×) while emitting MORE.** ⭐ `identity` is
**bit-identical** under the chain cut — it rides the NNLS weights, not the Gibbs mean.

**✅ SETTLED (MH-138): "§5 bagged NNLS vs §6b dense posterior" was never a conflict — one word, two
jobs.** MAGNITUDE = the Gibbs posterior mean (`beta`) ⇒ §6b. IDENTITY = `shapley_identity`, which
needs **fixed weights** from `canonical_M` = **bagged NNLS** ⇒ §5, for that job. And NNLS is
*required* there, not vestigial: the half-normal slab has a **strictly positive mean**, so an
un-informed family **cannot be zeroed** — measured, `beta == 0` in **0/5,117**, 100% positive —
whereas NNLS returns exact zeros. The built `bayes_shapley_identity` already does the right thing:
**NNLS fixes support+signs, Gibbs draws supply the width.** Rule of thumb elsewhere: ***bagged NNLS
for the GAUGE, Gibbs for the MODEL.***

✅ **CLOSED 2026-08-02 (MH-196) — ~~scale the literature set~~.** The versioned pull exists
(`eval/lit_ground_truth.py`, **329 genes / 92 families**, sha256-stamped) and replaced five producer-less
lists. **β is no longer at chance (0.436, p=0.021)** and abundance's apparent superiority does not survive
a fame control (Δ=−0.052, p=0.40 vs β's −0.085, p=0.038). ⭐ **And the item is CLOSED BY ARITHMETIC, not by
doing it again: only ~330 canonical families exist in the whole curated literature, and exhausting them
still leaves a CI half-width of 0.065 > the 0.033 effect — resolving β-vs-abundance needs ~1,241 clusters,
3.8× every canonical family ever published.** Do NOT re-open this as a scaling task.

---

## Axis 4 — Decoy validation (does curation beat a matched fake?)

**MH-137 is the row of record, ✅ CORROBORATED by MH-139 (both 2026-07-17).**

**Core gap −0.0119 (p=3.4e-05); deconv −0.0090; retention 0.76; 1,322 genes** (dose-matched by a hard
|Δdose|<1 caliper). **MH-139 re-measured it the other way** — caliper OFF, all 1,349 genes kept, dose
residual removed post hoc by regressing gap on Δdose with **b re-derived on this decoy (+0.01122,
95% CI [+0.0018,+0.0198])**: **gap at Δdose=0 = −0.0129 [−0.0205, −0.0052]**. **−0.0119 lies inside
that CI** ⇒ the caliper's 14.4% universe restriction (which drops precisely the HIGH-DOSE regulators,
1,164 vs 31 RPM median, p=2e-179) **cost nothing material.** Two designs, one answer: **≈−0.012 to −0.013.**

**✅ OOF POOLING DEFECT — MEASURED, 4%, MH-137/139 STAND (MH-181, 2026-08-01).** `oof_budget` pooled an
arbitrary-scale budget across folds with a per-fold-refit β (scale-non-invariant; it broke an exact null
elsewhere). Fixed and re-run: each ARM moved appreciably (real −0.0860→−0.0901, decoy −0.0690→−0.0737)
but they moved TOGETHER, so the **gap shifted only +0.0007** — dose-corrected **−0.0129 → −0.0124
[−0.0200,−0.0048]**, inside MH-139's own CI. ⚠ single-arm users (not the gap) were off by ~0.004–0.005 ρ.
Blast radius audited to ONE live site: the defect needs an arbitrary-scale index, so `gene_atlas`'s
y-unit OOF ceiling is unaffected.

**✅ DESIGN-MATCH CONTROL PASSES (MH-167, 2026-08-01).** The decoy matches ARMS and collapses to families
only afterward, so it is **not** matched on design **width** or within-design **collinearity** — both
asymmetries are real (p=4.6e-42 / p=6.7e-03; curation bundles same-seed mates, fakes do not) — but
**neither moves the gap**: contamination bounded at **≤~0.0025 ρ** (point estimate ~3%). **Δdose remains the
only design axis with a detectable effect, and it is already corrected.** ⚠ The pre-registered prediction
(that these make −0.0129 *conservative*) was **wrong** — the effect is undetected, not directional; carry the
**bound**, not a conservative framing. ⚠ MH-139's evidence artifact had been clobbered by a truncated
38-gene run and was **regenerated + verified identical** on the same date (see its row).

⚠ **The trajectory runs DOWN as controls improve** — this is the single most important thing to
internalise on this axis:

> **−0.045** (v3, loading-handicapped) → **−0.0306** (6mer-contaminated) → **−0.0147** (dose-weighted)
> → **−0.0119** (evidence-excluded + calipered)

Three defects fixed in MH-137: a **126× evidence hole** (only 2% of miRTarBase was excluded, so
CLIP/CLASH-supported *real* edges served as "decoys"; now 1,790,439 pairs excluded); the **6mer gate**
(Poisson-gated, not presence — a single 6mer is 24% chance noise at the median UTR); the **dose
caliper** (Δdose −3.05 → −0.105).

**The honest statement is a BOUND, not a win:** *real beats fake by ≤0.012 ρ, and every control
improvement so far has shrunk it.*
⭐⭐ **BUT THE "NEXT CONTROL FIX" IS NO LONGER OPEN — MEASURED AND CLOSED 2026-08-03 (MH-206).** This block
used to end *"whether the next fix (Manakov chimeric eCLIP + POSTAR3, still uncovered) takes it to zero is
genuinely open."* Manakov has been **in** the exclusion since MH-137 (448,330 pairs), and the POSTAR3 leg is
now settled on both the fact and the mechanism:
- **POSTAR3 IS on disk** (`data/external/POSTAR/human (1).txt.gz`, 676 MB, 2026-06-30) — but it is the **RBP
  binding-site table**: 221 RBPs, **2,360,006 AGO2 records, ZERO miRNA-named entries** ⇒ **no pair-level
  `(arm,gene)` call is derivable from it.**
- ⭐ **An AGO-peak ∩ seed-site layer is a NO-OP BY CONSTRUCTION** — `build_decoys` already requires
  `(a,g) ∉ sites`, and the strong-site + Poisson-6mer layers exclude every seeded arm. Only an **arm-resolving
  source on SEEDLESS pairs** can add anything; that is what Manakov is, and it moved **2.7%**.
- **ENCORI is 97.0% union-redundant**: 20,404 local pairs vs 3,861 in the labelled layer, but 19,799 already
  covered; of the 605 new, 578 are in the site map anyway and **0 of the 4,937 assigned decoys carry one.**
⇒ **the residual ≈−0.012 is NOT an evidence-hole artifact.** ⚠ Quote "97% covered / 0 decoys affected", not
"5.3× under-ingested" (axiom 5). Genuinely closing further needs POSTAR3's **separate miRNA-target/degradome
module** — not what is on disk — and layer (4) bounds its scale at a few percent.

**⭐ STANDS — MH-136, the arc's first real positive control, and its strongest result:** using the
genome-wide 6mer map, **187 genes whose curated edges are ALL seedless show gap +0.0006 — exactly
zero**; **810 all-seeded genes show −0.0325**. A priori, sequence-only, non-circular, and the model
behaves correctly on **both** sides. *No site → no repression → no coupling → no gap.*
It also closes the **21% positive-coupling mystery** (MH-132): controlled for regulator count,
`frac_seedless` predicts positive coupling (+0.042, p=2.2e-04).

**The internal null is back:** the **1-family stratum** (where β ≡ uniform exactly) gives
**−0.0007, p=0.15** on the clean control. MH-135's retraction of it was an artifact of the broken decoy.

**STANDS: the advantage over *unfitted abundance* is a FITTING effect, not curation** (MH-115/127).
A fitted matched fake also beats real unfitted abundance. **Rule: an abundance baseline is not a
control. Benchmark any aggregator against a fitted matched decoy.**

**Domain of validity — ⛔ THE "27% DOMAIN" IS RETIRED AS STATED (MH-147, 2026-07-17).** Re-tested on the
canonical decoy (1,349 genes):
- ✅ **Width axis (`n_fam ≥ 3`) HOLDS** — n=519, gap **−0.0264, p=4.2e-06**.
- ⚠ **Evidence axis (`w_max > median`) — "ADDS NOTHING GIVEN WIDTH" is TOO STRONG; it is UNDERPOWERED
  (softened MH-169, 2026-08-01).** MH-147 measured, within `n_fam ≥ 3`, high (−0.0298) vs low (−0.0130),
  **MWU p=0.293**. Re-measured on the restored canonical decoy: **A_COMPETENT −0.0335 (n=392) vs C_WEAK
  −0.0073 (n=69), MWU p=0.0618** — same VERDICT (does not clear the bar) but far closer than recorded,
  on a **69-gene** low arm. Per MH-157's language rule: **can't-tell, not "adds nothing".**
  **The operating rule still reduces to `n_fam ≥ 3`.**
- ⛔ **The complement is NOT zero** — MH-130's load-bearing "outside: +0.0002, p=0.72" becomes
  **−0.0114, p=6.3e-03 (n=936)**. The model works outside its claimed domain, just weaker.
- ⇒ **There is no partition. Competence is a GRADIENT in design width** — consistent with
  `spearman(ceiling, n_fam)=+0.551`. **The defensible claim:** the gap SCALES WITH REGULATOR-DESIGN
  WIDTH, is ~2.3× larger at `n_fam ≥ 3` than in the complement, and is **nowhere exactly zero**.
  ✅ **Reproduced 2026-08-01 (MH-169): −0.0295 vs −0.0105 (p=3.6e-03), ~2.8×; internal null holds
  (`n_fam==1` → adj −0.0033).** ⛔ **But WIDTH CANNOT BE SEPARATED FROM ABUNDANCE** — `spearman(n_fam,
  n_abund)=+0.780`, only **20** genes are narrow-but-abundant, and **both partials vanish**
  (width|abundance +0.005 p=0.85; abundance|width −0.016 p=0.55). ⇒ they are ONE axis here; the gap
  tracks whichever you measure and is **not attributable to either**. The sharpest separator is
  MEASURABILITY (ceiling ≤0 → +0.0094, a clean zero; >0.05 → −0.0621) — but that is **partly
  definitional** (ceiling≤0 forces ρ_real≈0). Home: **MH-169**, `eval/decoy_by_gene_strata.py`.

**The gene lens — ✅ REPRODUCIBLE 2026-07-17 (MH-144), and one headline RETIRED:**
- ⛔ **"17.6% NOT MEASURABLE" is FRAGILE — do not cite it.** It is a threshold count evaluated exactly
  where the distribution **piles up**: 39.8% of genes sit within ±0.01 of zero ceiling. A faithful
  re-implementation gives **25.5% on MH-130's own genes** while the per-gene ceilings agree at
  **corr 0.9956** — a systematic shift of **−0.00136** (0.06% of the R² scale) reclassifies ~8% of the
  universe. **The ROBUST statement: `ceiling ≤ 0.02` → ~52–53%**, which reproduces across implementations.
  ⚠ Same failure mode as MH-138's per-gene median (see Axis 3) — **twice in one day**.
- ✅ **What reproduces** (current 1,549-gene universe): ONE seed family **47.1%** · `A_COMPETENT`
  **27.6%** · spearman(ceiling, n_fam) **+0.551** ⇒ the lens's STRUCTURE is sound.
- ⚠ **MH-130's universe was stale + undocumented** — 1,421 genes, a strict subset of today's 1,549; the
  150 it never saw are 96.9% single-family. It almost certainly built on the decoy-eligible set (cf.
  MH-127's 1,436), which its doc never records.
- The ceiling is **design capacity, not biology** — target detection fails only 7.9% of genes;
  **regulator-design width is the binding constraint.**
- **48.1% of genes have ONE seed family, where β ≡ uniform exactly** ⇒ *"does the learning add
  anything?"* is **UNDEFINED for half the universe — a non-question, not a null.**

---

⭐⭐ **2026-08-03 (MH-203/204) — THE FITTED-FAKE QUESTION IS NOW ANSWERED ON THREE LAYERS, AND THE
ANSWER SPLITS BY LAYER, NOT BY COHORT.** Everything above compares β against an *unweighted* reference. The
harder test — a **FITTED** fake (β refit on the site-free decoy design, identical folds/C/estimator) —
now exists in all three places:

| test | real | decoy | gap | p |
|---|---|---|---|---|
| **mRNA, TCGA** (n=1,354) | −0.0807 | −0.0579 | **−0.0189** | **6.69e-09** |
| **mRNA, transported to BUFFA** (multi-family, n=243) | −0.0334 | +0.0035 | **−0.0374** | **3.54e-06** |
| protein, RPPA (n=180) | −0.0186 | −0.0264 | −0.0086 | 0.554 |

⇒ **β beats a fitted, matched, site-free fake on mRNA in BOTH cohorts** — which retires the standing caveat *"every win is fitted β vs
UNFITTED abundance"*. ⚠ **It also survives the subtype control:** within-PAM50 gap **−0.0154 (p=4.2e-08)**,
**retention 0.81** (real −0.0528 vs decoy −0.0363, n=1,354) ⇒ the mRNA effect is **not** lineage-scale.

⛔⛔ **THE PROTEIN ROW IS *NOT* A LAYER BOUNDARY — CORRECTED SAME-DAY (MH-205).** This block first read
*"fails only on protein ⇒ a layer boundary, not a general negative"*. **That was wrong, and the error was
comparing a 180-gene protein test against a 1,354-gene mRNA test.** Three things were confounded — gene
selection, power, and layer:
* ⭐ **on the SAME 180 genes the *mRNA* gap is −0.0103, p=0.0565 — also n.s.** Those genes are weak on
  **both** layers; RPPA's antibody panel is signalling-biased, not a random draw. Matched, protein retains
  **84%** of the mRNA gap, not the 45% the naive comparison implies.
* ⭐ **the MDE at n=180 genes is 0.0309 / 0.0232 — both observed values sit INSIDE their own noise floor.**
  Detecting the mRNA-sized gap needs **~514 genes**.
* ⭐ **NOT the phospho contamination I had named** — 13/180 genes, and dropping them moves −0.0086 → −0.0089.
* ⭐⭐ **paired on gene at 5× the genes, the layers are indistinguishable:** CPTAC `tcga105` (849 genes, the
  SAME patients ⇒ layer-only) **+0.0018, p=0.287**; `prospective` (897, TMT-11) **+0.0056, p=0.341**.
* ⚠ **What survives:** the direction is consistent in all three instruments (Stouffer z=+2.44, p=0.015 —
  **anti-conservative**, 151/180 genes recur), so a **small** protein penalty of order 0.002–0.006 ρ stays
  plausible — an order of magnitude below what was claimed.
⇒ **the protein half is UNRESOLVED, not a measured negative. The mRNA results are untouched.**

✅ **AND IT CARRIES A ZERO-CROSSING INTERNAL CONTROL — the strongest this arc has produced.** Stratifying the
gap by `gene_repression_class` crosses zero exactly where repression stops: constitutive_repressed −0.1464
(n=115) · gained −0.0881 (n=165) · lost −0.0126 (n=250) · **never_repressed +0.0034 (n=814)**. The
specificity signal is **absent exactly where the biology says it should be**.

⚠⚠ **STALE-ARTIFACT WARNING, and it inverted a verdict.** `decoy_family_betas.tsv` / `decoy_betas.tsv`
(2026-07-13) overlap the CURRENT matched pairs by **1.8%** and predate MH-137's three fixes ⇒ their decoys
are **too weak and flatter real**. Reading them made MH-201's Buffa decoy control read *"n=120, p=0.335,
underpowered"* when the true answer is **n=243, p=3.5e-06**. Both files are now **archived** at
`output/learned/ARCHIVE_PRE_MH204/`; the live artifacts are `decoy_oof_budgets.parquet` +
`decoy_family_betas_oof.tsv` (**4,889 cells / 1,395 genes**). **Zero live consumers of the stale tables.**

---

## Axis 5 — CPTAC / protein

**⛔ DEAD — do not rebuild:** `βᵗ` (the translational-repression latent, the plan's centrepiece) is
**not supported at n=101**, and all prior positive evidence was a **LEAK** (mediator slope `a` fitted
on all samples, then OOF-ed; fit `a` inside the fold and it vanishes). Four leak-free framings agree;
decisively the aggregate (1 df, not 60): **BH q<0.10 in 1/17 genes — BCL2 alone; PTEN d=−0.006,
p=0.82** (MH-103).

**⚠ MH-34 is NOT overturned** — 13/17 stay repression-directed (binomial p≈0.025). The falsification
is of the **modelling object**, not the biology. The translational residual is ~1% of genes; 17 genes
predicts ~0.2 hits; we got 1. **There was never a βᵗ field. Only n fixes this.**

**STANDS — protein can never be a coupling lever:** it carries **4–6%** of the mRNA channel's Fisher
information about β, and **≤7.6% at any `a_g`** (pre-registered ceiling; a 5× correction to the
number moved the conclusion zero). Converges with the state channel's independent 0.6%.

**STANDS — MH-114, the best-designed test in the program:** ~**57%** of the protein validation is
**compartment arithmetic**. Degree-preserving shuffled null, stratified by compartment orientation:
**shuffled non-edges reproduce the real gradient exactly (−0.1290 vs −0.1272).** Pairs that *cannot*
repress produce the same bulk anti-correlation. The true compartment-blind effect is **≈−0.011**,
surviving adjustment at **≈−0.0065** (p=2.0e-03). Damage: gene FDR-neg **34→3**, hallmark **9→0**,
EMT **−0.330 → +0.162 (sign flip)**, orphan candidates **594→23**.

⚠ **The protein transfer is NOT decoy-controlled** — a site-free fitted fake reaches **99%** of it
(MH-130d). It is real arithmetic, but **not evidence the curated edge set is right.** The only
decoy-controlled transfer in the arc is the **CPTAC mRNA** cell inside class A.

✅ **RETENTION MEASURED 2026-08-01 (MH-171) — the channel is ATTENUATED, not erased.** MH-107's open
item is closed for `ood_protein` (its code fix had landed but was never run/recorded): **prospective —
the INDEPENDENT cohort, n=101 — median retention 0.797, 7/7 genes still protein-coupled after
composition**; TCGA-105 (circular) 0.656, 7/7. ⭐ **ZEB1 retains WORST in both (0.246 / 0.308 ⇒
`composition_explained`)** — a per-gene confirmation of the compartment-arithmetic claim above, from a
different design than MH-114's shuffled null, agreeing on which gene. ESR1 / BCL2 / CDKN1B are
cell-intrinsic in both. ⚠ 7 hub genes only; ⚠ GATA3 (`n_reg=1`) disagrees across cohorts — do not cite;
⛔ `model_beats_abund` is FITTED-vs-UNFITTED and is **not** a control.
✅ **`MH-83` / `MH-84` — RESOLVED 2026-08-01 (MH-182): a REGISTRY-FILING gap, not missing science.**
The work is recorded in `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md` (capstone of MH-82→MH-84) with numbers
in `LEARNED_MODEL_VALIDATION.md` §1. Nothing to re-run; no rows to back-fill.


⚠ **MH-114's "34→3" was measured on the retired heuristic** (`compute_gene_pressure`), not the
learned posterior (MH-115) — and its FDR counts rest on the uncalibrated null (MH-123). The
*conclusion* is unharmed; **do not reprint the counts.**

**The confounder-ARCHITECTURE principle (MH-111) — the most generalisable thing here:** proliferation
and host localise to an **ARM** (LOCUS properties); composition localises to the **FAMILY** (an
EXPRESSION property, because same-seed members sit at different loci but are co-expressed).
**The unit at which a confound acts predicts which tool can fix it.**

---

⭐⭐ **2026-08-03 (MH-203) — THE RPPA CHANNEL IS OPEN AT n=866, AND IT IS NOT A REPEAT OF THE CPTAC
ARC.** `data/rppa/` was never touched by this subproject: **881 × 461 participant×antibody, 866 with
miRNA — the SAME patients, so there is no cohort boundary.** The power argument is the whole point:
detectable |ρ| at 80% power is **0.276 at CPTAC's n=101** vs **0.095 at n=866**, against a typical
gene-level coupling of **0.07–0.12** ⇒ **CPTAC's detection floor was ~3× the effect it was asked to
measure.** ⚠ This does **not** rescue MH-103 — that died of a **mediator leak**, a design flaw, and
`protein-βᵗ` stays dead.

**⭐ THE RESULT — DESTABILISATION, NOT TRANSLATION.** Residualising protein on its own mRNA isolates the
translational channel, which the project has never had the n to separate: **`rho_mrna` carries the signal;
`rho_disc` (protein-beyond-mRNA) does not.** ⇒ what the model measures is transcript destabilisation.
⚠ **The fitted-decoy control on protein reads −0.0086, p=0.554 (n=180) — but that is UNDERPOWERED AND
GENE-SELECTED, not a measured negative (MH-205; see Axis 4).** On the same 180 genes mRNA is *also* n.s.
(p=0.0565), both sit inside a 0.023–0.031 MDE, and paired on gene at 849–897 CPTAC genes the two layers are
indistinguishable (+0.0018 / +0.0056). ⇒ **do not cite a "layer split".**
✅ **What DOES survive, and is now measured on 849/897 genes instead of 180: the protein-beyond-mRNA
channel is flat everywhere (−0.0011, p=0.657; +0.0034, p=0.659) ⇒ DESTABILISATION, NOT TRANSLATION.**

---

## Axis 6 — Progression / state

**Two different objects. Do not conflate them.**

**⛔ The state CHANNEL is CLOSED — measured and cancelled** (MH-102d, 2026-07-12).
`learned/channel_state.py` **was never built** (verified absent). τ — the acquired-vs-constitutive
payload, the axis's entire deliverable — is **indistinguishable from zero** (GTEx τ=0.0009, info
**0.6%**; NAT 0.7%; CPTAC 0.7%). Structural: channel precision ∝ **a²**, and a≈0.11 ⇒ **even
unlimited GTEx donors add ~1%**. *Composition attenuation does not merely bias the healthy β — it
destroys its information content.* The additive `M_T = M_H + Δ` form **would have faked the
headline** (Δ>0 on nearly every edge, manufactured out of attenuation). **Do not rebuild it.**

**✅ The abundance-level cross-state LANDSCAPE is built and run** — GTEx-healthy → TCGA-NAT → tumor,
721 arms / 5,108 edges. Verified against the persisted 2026-06-27 artifact: **acquired_realized 640
(129 arms) · lost 209 (89 arms) · stable 2,141 · acquired_unrealized 1,775**. ⚠ Its **per-edge FDR
class labels** rest on the null MH-124 measured as 3–4× too narrow; `LANDSCAPE_REPORT.md` carries no
banner about this. ⚠ NAT is structurally underpowered (n≈104): `nat_decoupled` = **1 edge**.

**STANDS — the panel question is SETTLED, zero new CIBERSORTx runs needed:** β is **ρ=0.94** across
entirely different atlases (Wu-9 vs HBCA-11) and ρ=0.99 across reference resamples, but only ρ=0.81
vs **no C at all** (mean |β| halves). ⇒ ***which* composition control you use moves β ~10× less than
*whether* you control composition.** The "reference reconciliation blocker" is **RETIRED**.

**What survives the cancellation:** `learned/gauge.py` and `learned/confounders.py` — *they are what
killed the axis, which is what falsifiable infrastructure is for.* **Standing rule:** any future
cross-cohort channel runs `gauge.calibrate()` → read `info_ratio` + τ **BEFORE** it is built.

**✅ NEW — the WITHIN-PATIENT PAIRED axis is a THIRD object (MH-158, 2026-07-18), and it works where the
cross-cohort channel couldn't:** differencing tumour − **own-NAT** (n=103 matched pairs) removes the patient
baseline, bypassing the composition attenuation that killed the channel. Engine `learned/realization.py`
(M = `M_complement`, fit on the 975 non-matched tumours — leak-free). **DOSE/potential:** the acquired dose is
~46% patient-specific (`own_specific_frac` median **0.46**) ⇒ own-NAT ≠ cohort/subtype anchor; within-gene
share shifts recover canonical grip (miR-21→PTEN). **COUPLING/realized:** a WEAK set-level signal — real edges
beat matched site-free decoys by **~0.017 ρ in the OPPOSITE-compartment stratum after composition-Δ**
(gene-clustered p=2.7e-4; SAME n.s.). Own-NAT gives a stronger raw within-patient ρ than cohort (p=2.8e-6) **but a decoy control DEMOTES that advantage to SHARED-baseline** (decoys show the same own−cohort gap, MWU p=0.92) — the repression-specific signal is the real-vs-decoy gap alone, reference-independent. ⚠ Tiny
(~1/6 single-gene SE), **NOT per-edge**, single cohort; **NAT≠healthy** ⇒ final malignant step only.
rigor-auditor **A** (CONCERNS→re-scoped: null-referenced gap not −0.050/64%; retention median 0.61 not 38%).
✅ **Phase 2/3 DONE (2026-07-18, rigor-auditor PASS, all 9 predictions pre-registered+confirmed).** Res-3 family:
nonlinear pooling adds no coupling (`family_agg` −0.046 ≈ arm-level −0.050). **Owner-convergence** (Res-4): the
realization-Shapley owner tracks the static-identity owner above a within-gene random-owner null (z=6.1) **but the
coincidence is partly MECHANICAL (shared M vector), realization~dose is n.s. (z=1.6), and realization~static is not
distinguishable from static~dose ⇒ the realization-SPECIFIC ownership content is NULL.** Retention×realization
near-null. Genomic-context (P): host-coupled retain+realize more than intergenic (antisense exception noted).
Hallmark (P, descriptive): metabolic/differentiation programs most realized. Res-5 within-family at chance (prereg
null). Own-retaining 7% not clinically distinct (fragile purity hit, n=5). Details in MH-158's registry row.
✅ **CLASS → REALIZATION (MH-160, rigor-auditor A):** the cross-state `shift_class` PREDICTS within-patient
realization, ordered by repression-LIVENESS not dose — live-in-NAT edges {constitutive, field_established}
realize more than acquired-only ones (Δ −0.039, gene-clustered) despite LOWER dose; anti-dose (dose+coupling_tum
controlled −0.032), decoy-arbitrated (direct paired diff −0.024, p=0.01); partial circularity (ρ=0.234 with the
defining coupling_nat) quantified + defused. Enabled by the integrated `edge_card`.
✅ **`shift_class` REBUILT → CALIBRATED TWO-AXIS (MH-166, 2026-07-19):** the old coupling-only −0.1 label is
SUPERSEDED by a two-axis annotation = per-state **calibrated** coupling (`site_free_null` per state, replaces −0.1;
tumour null 2.67× reproduces MH-123) × **same-platform** `arm_lfc_NAT_TUM` dose (replaces the soft QN healthy leg).
New quantified cols (`coupling_p_*`, `coupling_z_tum`, `realization_score`, `dose_class`, `wiring_frac`) + NEW class
`dose_acquired_uncoupled` (571 edges: dose up, no calibrated coupling). ⭐ **At EQUAL dose the calibrated coupling
axis alone separates realization** — `dose_acquired_uncoupled` +0.002 vs `acquired_realized` −0.101 (dose 1.13 vs
1.19). MH-160 RE-VALIDATED: core ordering holds/strengthens (−0.046) but the decoy-arbitrated increment is now
MARGINAL (−0.026, p=0.098) as `constitutive` collapsed 127→6 under calibration. Home: `card.py::_shift_class`.
⚠ **Healthy-leg fallacy (MH-166 follow-up, 2026-07-19):** the GTEx coupling leg is BLIND where the uniquely-mappable
pipeline collapsed multi-mapping arms. Added `healthy_leg`/`healthy_potential`/graded `healthy_uninformative`
(miTED-aware, 305 flagged; contamination 24%→13%; MH-160 robust). Flagged arms = dose-constitutive, coupling-blind;
recoverable via same-seed co-transcription surrogate (miR-17 via miR-20a ⇒ coupling-ACQUIRED). Genome-wide surrogate
node + per-target GTEx→NAT→TUM identity trajectory = NEXT ARC. See FORMULAS §11g + DISCOVERY_REGISTRY MH-166.
⭐⭐ **EXTERNAL VALIDATION ARC (MH-166, 2026-07-19) — the big current-state update:**
• **PRECISION validated externally, TWO orthogonal axes:** McGeary&Bartel functional (19 arms, Fisher **p=7.96e-33**,
  HE **5.6×** enriched vs null) + breast chimeric binding GSE263552 (**27×**). The edges we use are real.
• **Family-collapse validated FUNCTIONALLY** (miR-29a≈29b, per-target ρ 0.87–0.95) — but only the seed-HOMOGENEOUS
  case; within-family TCGA arm-corr is HETEROGENEOUS (med 0.47, 6 tight/6 loose; miR-29 raw 0.44 ≈ partial 0.43,
  so co-expression ≠ the 0.87 seed-driven interchangeability). ≈25% seed-heterogeneous families still gated per-family.
• **⭐ COUPLING IS REALIZATION, NOT EXISTENCE (supersedes the "weak set-level signal" framing above):** breast TCGA
  coupling is INDEPENDENT of edge reality (Fisher OR=0.99 vs McGeary-functional) and ~13% SENSITIVE — it measures
  realized DOMINANT regulation across patients (observational, variance-driven), a harder/different question than
  existence. High-specificity (2.9–5× null) / low-sensitivity. ⇒ judge EXISTENCE on functional/binding; use coupling
  as REALIZATION (presence = dominant realized regulator; absence ≠ not-real). Do NOT use coupling as an existence arbiter.
• **RECALL → DISCOVERY (expansion pending):** HE captures ~6% of functional targets (precision-by-design). Built the
  tiered candidacy + `model_expansion_list.tsv` (functional × seed-type/Kd × context × chimeric × family-mate) tied to
  the project's coupling-discovery lane → **21 gold∩external crown jewels** (miR-17~92→TMBIM6/AHNAK/XIAP; miR-182→MET
  is literature-VALIDATED in breast, the rest novel-but-assay-backed). ⚠ **abundance is NOT a "tilt"/artifact — abundant
  arms being dominant realized regulators IS the biology; a "real edge" = structure(seed) + physical(binding) +
  impact(coupling), all three. The abundance association only limits RECALL (the low-abundance tail), not precision —
  do NOT down-weight coupling.** ⬜ **OPEN: MODEL EXPANSION + Gibbs REFIT (board item 0) — fold candidates in, re-fit,
  read new β + gene pressure.** Files in `output/learned/realization/{model_expansion_list,missed_edges_*,breast_chimeric_*}.tsv`.

---

## Axis 7 — Buffa / independent cohort / outcome

**⚠ "METABRIC" in MH-73/74/75/76 means Buffa, n=207 — verified in code**
(`pressure_prognostic_signature.py:166`, `pressure_prognostic_gene_centric.py:97` both import
`load_buffa`). ⛔ **BUT "Buffa's miRNA arm IS the METABRIC invasive miRNA subset" — WRITTEN HERE UNTIL
2026-08-03 — IS WRONG (MH-201; the full disproof is in the METABRIC block below).** Buffa = **Oxford**,
surgery **1989–92**; METABRIC = Cambridge/Vancouver; **0/207** arrays match at r>0.90. Buffa is GSE22216
miRNA + GSE22219 mRNA, **207 paired tumours** (not 210), **DRFS only, 77 events**, no TCGA overlap — a
**genuinely independent cohort**, which is what makes MH-201/202/204 worth anything. ⇒ the
+0.056 / +0.060 / +0.028 / +0.019 are **four views of ONE 207-patient result**, not four validations.

**STANDS (the well-earned negatives):**
- **Pressure magnitude is exhaustively non-prognostic** — 0 recurrence FDR across 3,368 features.
- **Functional > magnitude** — realized +0.060 vs received-pressure −0.016 (a *within-cohort*
  contrast that cohort-specificity cannot explain).
- **miR-210 → DRFS reproduces Buffa 2011** — a working **positive control**, which is what makes the
  pressure nulls credible rather than merely underpowered.

⭐⭐ **2026-08-03 (MH-202) — THE *LEARNED β* IS NOW ON OUTCOME TOO, AND THE VERDICT SURVIVES THE UPGRADE.**
Everything above tested the **heuristic** pressure, so "pressure is non-prognostic" was open to the reply
*"you tested the weaker object"*. `learned/eval/beta_prognostic.py` closes that: Buffa DRFS + TCGA OS/DFI on
the learned coupling model. **Per-gene YES, aggregate NO — and the WEIGHTS earn nothing** (β-weighted
budgets do not beat unweighted abundance on outcome). ⇒ **a good coupling estimator is not thereby a
prognostic one; the two orderings are different objects.**
✅ **Credible because GATED, not merely null:** `gate()` requires miR-210→DRFS to reproduce Buffa 2011 and
**the module refuses to report if it fails** — which is exactly what separates a well-earned negative from
an underpowered one. ⚠ TCGA is the **weaker** cohort here (OS only, no DRFS) — the same direction, not an
independent confirmation.

**⛔ DEAD — the "TCGA is sparse so it cannot arbitrate" defense.** MH-76 ran the **frozen-panel**
test (Buffa-trained, frozen, scored on TCGA), which removes overfitting as an explanation entirely:
**DFI +0.002, OS +0.006; panel-alone C 0.48–0.53 ≈ random.** MH-73/74 still assert the superseded
framing and need annotating.

**⛔ SETTLED BY MEASUREMENT 2026-07-17 — "108 orphans triple-validated" IS 3.** Re-ran
`eval/buffa_validation` against the composition-adjusted screen (594 → 23 candidates since 2026-07-12,
MH-114): **`triple_validated` 108 → 3; ECM/collagen among them 30 → 0** (all rows 64 → 0). The three
survivors are all miRGeneDB **guide** arms — **miR-195-5p→PSMD7 · miR-181b-5p→IRS2 · miR-22-3p→STK39**
— and only **1** is neg-sig in Buffa. **The miR-29→collagen/ECM headline is gone.** The pre-registered
expectation ("single digits, possibly zero ECM") was correct.
⚠ **The module had been DEAD, not merely stale.** `GEO_DIR` was a `__file__`-relative hop
(`parent.parent/"data"`), correct when the module sat at the subproject top level — but it was moved
into `eval/`, which added a directory level and pointed it at a non-existent `mirna_hallmark/data/`.
Every input raised `FileNotFoundError`; the GEO cache was at the **repo root** all along. Fixed to
`C.REPO_ROOT`. **That is why this lane was never re-run after MH-114: it could not run.** Ripple
checked — every other `__file__`-relative path in the subproject is correct (`learned/readouts.py`,
`learned/cptac_data.py`, `presentation/make_figures.py` use `parents[2]`; `config.py` correctly uses
`parent.parent` at the top level).
✅ **CONTROL — the HE-edge legs reproduce EXACTLY** (concordance **+0.319** · sign-concordance
**0.593** · TCGA neg-sig same-sign **0.673** · neg-sig-in-Buffa **0.128**), so only the orphan lane
moved. Pre-fix outputs preserved at `output/buffa_validation_PRE_MH114/`.
⚠⚠ **THOSE FOUR NUMBERS ARE UNDER-CONTROLLED — RE-RUN 2026-08-02 (MH-201 stage 3b).** They condition on
**[prolif, ER] only**; the lane predates any Buffa composition block, which is exactly the setup MH-114
showed produces a 90% artifact gradient. Re-scored under the composition-aware C (reproduction gate on the
`raw` block returns all four to the digit, so the comparison is valid):

| metric | recorded (raw) | +ESTIMATE/CAF | +NNLS lineages | retention |
|---|---|---|---|---|
| Buffa↔TCGA concordance | 0.319 | **0.215** | 0.217 | **0.672** |
| sign-concordance | 0.593 | 0.567 | 0.572 | 0.956 |
| TCGA neg-sig same-sign | 0.673 | 0.643 | 0.597 | 0.955 |
| **TCGA neg-sig SIGNIFICANT** | **0.128** | **0.052** | **0.038** | **0.402** |
| Buffa neg-sig fraction | 0.052 | 0.026 | 0.016 | — |

⇒ **~⅓ of the headline concordance and ~60% of the significance-based replication is composition.** ⭐ **But
the DIRECTION metrics barely move (retention 0.955–0.956) — the same dissociation as the gene-level transfer:
direction/level survives adjustment, significance/magnitude does not.** Quote the adjusted column.

✅ **CLOSED 2026-08-02 (MH-201) — ~~β has never touched Buffa~~.** `learned/eval/ood_cohort.py` built;
1,364 genes, β FROZEN (zero params fit in Buffa), FAMILY rung. Stage-0 gate reproduced all six recorded
HE-edge controls exactly before any β was involved.
**⭐ THE RESULT IS A DISSOCIATION, AND IT INVERTS MH-174.** Against a **C-ABLATION ceiling** (the TCGA
ceiling recomputed under the reduced C Buffa can support — without it, "cohort boundary" is confounded with
"weaker confounder block"; dropping `target_cn`+CPE alone costs **33%** of TCGA |ρ|, −0.0781 → −0.0521):
* **LEVEL transports** — β **−0.0186** (58.5% repression-directed, p=1.5e-05) vs unweighted abundance
  −0.0022, **Δ −0.0072, p=6.2e-13**; vs a within-gene β-permutation **p=7.8e-17**; **retention +0.358**.
  Restricted to genes with any identified (|z|>2) β mass: **−0.0297 vs −0.0091, p=3.9e-10.**
* **RANK does NOT separate β from abundance** — +0.222 vs +0.183, increment **+0.040 [−0.013, +0.093],
  p=0.13 n.s.** ⇒ **MH-174's "transport the RANK, never the LEVEL" is a TCGA→CPTAC statement, not a law.**
* ✅ **Survives the pre-registered MH-114 null**: the orientation gradient is arithmetic (real −0.0181 vs
  shuffled −0.0209) but the **compartment-blind contrast is −0.0204 / −0.0232, SAME SIGN in both
  orientations** — genuine repression, ~2× CPTAC's ≈−0.011.
⚠⚠ **43% of genes are SINGLE-family, where β is MATHEMATICALLY INERT** (`M = β·Z` vs `Z` ⇒ identical
Spearman; verified max|Δ| = 0.00e+00 over 587 genes). That stratum alone gives rank +0.234 **from abundance**.
Never pool it with multi-family genes.
⚠ The in-cohort refit (−0.0099) is **NOT a ceiling** — transporting beats refitting at n=207.
⭐ **WHERE THE TRANSFER LIVES — the 2×2 gene profile (concentration × identifiability):**
`peaked/determined` (n=244) is strongest on every axis (level **−0.0355**, rank **+0.238**, Δ vs abundance
**+0.063**, p=1e-08); `peaked/UNDETERMINED` (n=87) has **no effect at all** (p=0.33). ⇒ **β concentration
buys nothing until the peak sits on a family the model could determine.**
⛔ Do NOT cite "corr(n_fam, HHI) = −0.667" as evidence that big designs are diffuse — **that is the HHI
floor (1/k); floor-corrected it is −0.075.**
⭐ **WHERE IT WORKS (82-modifier FDR scan on Δ = ρ_β − ρ_abund):** ① **more regulators ⇒ bigger β advantage**
(`ctx_ceiling` −0.185, `n_fam` −0.162, `n_arms` −0.156); ② ⭐ **Buffa DYNAMIC RANGE matters, mean level does
not** (`buffa_sd` −0.144 q=0.0016, `buffa_iqr` −0.123 — but `buffa_mean` q=0.45 n.s.: a flat gene cannot
couple however highly expressed); ③ **composition-driven genes transfer WORSE** (`comp_tcga_mrna_driver_share`
+0.127); ④ ⭐ **CPTAC-transferable genes are Buffa-transferable** (`cptac_t105_agg_rho_rna_raw` → Buffa level
+0.234, q=5e-06) — two independent target cohorts agreeing on which genes transport.
⭐⭐ **THE STRONGEST MODIFIER OF ALL IS ON THE REGULATOR SIDE: `reg_dose_hhi` +0.186 (q=2.1e-05)** — when
the regulators' ABUNDANCE is concentrated in one family, **β LOSES to unweighted abundance**, because the
sum already IS that family. With `reg_frac_flat` +0.180 and `ident_n_def` −0.172. ⇒ **β earns its keep only
where the regulator ensemble is SPREAD and MEASURABLE; what matters is how the DATA distributes, not how β
does.** Shapley enters as COVERAGE (`ident_n_def`) — every distributional form of it is null.
⚠ **THE HARM: β helps 55.7% of genes but HURTS 34.7%** (mean gain −0.0614 vs loss +0.0397, 1.55×; tails
26 vs 4 at |Δ|>0.15; net −0.0204). ⭐ **BUT IT IS GATEABLE — and the two axis families do DIFFERENT jobs:
COMPOSITION axes predict the SIGN (does β help or hurt), ENSEMBLE axes predict the MAGNITUDE.** Gating on
the ensemble axis alone leaves harm FLAT (33.9%); **BOTH gates together → hurt 34.7% → 27.6%, β-sign-right
58.5% → 69.3%, mean Δ −0.0468, on 19% of genes — and retains every showcase gene (MYC, CCND1, HIF1A,
FOXC1, TXNIP, EGFR).** FN1 is not a surprise: 89th percentile on `comp_tcga_mrna_driver_share`.
⭐ **BUT THE GATES REMOVE COIN-FLIPS, NOT DAMAGE.** At n=207, SE(ρ)=0.070; only **8 of the 35 post-gate
hurt genes exceed 1 SE** ⇒ genuine harm is **6.3%**, not 27.6%. Down the gate ladder the **nominal** hurt
falls 34.7 → 27.6 → 21.6 → **18.8%** while the **genuine-harm rate is PINNED at 6.8 / 6.3 / 6.8 / 6.2%**.
An 81-axis FDR scan inside the gated set gives ONE hit (`median_retention` q=0.04 — a third
composition-flavoured axis) **which does not separate the harmed genes** (p=0.80). ⇒ **the ~6–7% where β
does real damage (MET, SMAD7, IL6, CDKN1A, ETS1, BCL2L2, XIAP, CDK4) is irreducible with the current
axes.** ⚠ Always separate a nominal sign-flip rate from a magnitude-qualified one.
⛔ **The intuitive "they're immune genes" explanation is REFUTED** — immune-hallmark membership (p=0.78)
and the gene's own immune-compartment loading (ρ=−0.023, p=0.56) both null, though IL6 really is the most
immune-loaded of the eight (+0.32). ⭐ **The composition axis that predicts harm is STROMAL/VASCULAR**
(driver = CAFs 178 · Endothelial 67 · PVL 58, vs 81 for all immune compartments combined) ⇒ **mesenchymal
content is what fails to transport between two bulk breast cohorts; say WHICH compartment, never just
"composition".**

**⚠ "METABRIC-full (EGA-pending)" is doing unearned work across five docs.** There is **no accession,
no DAR id, no submission date, no owner** anywhere in the repo. The evidence supports *"never applied
for"*, not *"pending"*.
⛔ **CORRECTED 2026-08-02 (MH-201): the rest of this paragraph was WRONG — Buffa is NOT a METABRIC subset.**
Verified three ways: **0/207** Buffa samples match any of the 1,980 METABRIC arrays at r>0.90 (best-match r
median **0.270**, margin **0.0094**); **49%** have no clinical candidate on (age±1, nodes, ER, size±2mm,
grade); and provenance disagrees outright — **Buffa = Oxford, Weatherall Institute, surgery 1989–1992**
(GSE22220), METABRIC = Cambridge/Vancouver. ⇒ **METABRIC-full would be a DIFFERENT independent cohort, not
merely a power upgrade**, and the request should be re-assessed on that basis. ⚠ But note METABRIC public
has **no miRNA at all**, so it cannot run this test either — the EGA request is specifically for the miRNA
layer. ⛔ Related: `target_cn` is **unobtainable** for Buffa (GSE22220 has exactly two subseries, miRNA +
mRNA; Buffa 2011 generates no CNV; no later study profiles those tumours). Do not re-search it.

---

## The four things worth doing next

1. ✅ **DONE 2026-08-02/03 (MH-201 + MH-204) — `β(TCGA) → Buffa mRNA`.** `learned/eval/ood_cohort.py`;
   **LEVEL transports (−0.0186, retention +0.358), RANK does not separate β from abundance**, and — after
   the stale-decoy fix — **the transported β beats a FITTED fake, −0.0374, p=3.5e-06.** The pre-registered
   MH-114 compartment stratification ran: the gradient is near-identical real-vs-shuffled, as predicted, so
   the compartment-blind contrast is the readout.
2. ✅ **DONE 2026-07-17 — `buffa_validation` re-run; MH-38 / MH-55 / MH-73 / MH-74 all annotated.**
3. ✅ **DONE 2026-08-01 (MH-168/169) — competence map rebuilt against the canonical decoy.**
   ⚠ **It was BLOCKED, not neglected:** `gene_atlas.py` had been **dead for 15 days** — the same
   `__file__`-hop-after-a-move bug as MH-38, introduced by commit `e5d5d84` relocating it one level
   deeper; it ran its full compute and died on the final write. Fixing it exposed a **second** bug:
   `load_cnv_target_genes` clobbered its own cache with each single-gene request and, under 8 workers,
   corrupted the gzip — making **`assemble_gene` fail for every gene**. Both fixed (MH-168).
   Strata result in **MH-169**; the atlas reproduces MH-144's structural numbers exactly.
4. **Restore an exogenous existence validation** — both `CN_INSTRUMENT.md` and MH-126 name this as the
   highest-value open item, and nothing else changes the program's standing as much.

---

## Doc traps — read this before trusting a doc

| doc | trap |
|---|---|
| `LEARNED_MODEL_ESTIMATOR_MAP.md` | Billed as the "canonical map", but self-dates to **2026-07-06** and its main table is **pre-convergence** (lasso-as-primary, "Bayes = uncertainty layer"). Only its two banner UPDATEs and the pooled-HE policy are current. |
| `MH127_*` / `MH128_*` / `MH130_*` | Read as current; their central magnitudes are **superseded by MH-137**. MH-130's *atlas* half stands; its *gap* half does not. |
| `MH124_*` §5 / §5b | **VOID** (`pi_causal` not exogenous). **Do not cite.** §4b is untouched and stands. |
| `CPTAC_PROTEIN_CHANNEL_PLAN.md` | ✅ **FIXED + ARCHIVED 2026-07-17.** §6's VOID PAM50 numbers (1.7% / ≤8.8%) corrected to **4–6% / ≤7.6%**, the self-nested correction strings de-nested, and a third live citation of the void figure found at §0a. The doc is now in [`archive/`](archive/README.md) — its channel is dead (`βᵗ` falsified at n=101) and its **survivors were extracted to `LEARNED_MODEL_VALIDATION.md` §6.1–§6.5** (`locus_cn_cptac` · δ-transportability · Bar-5 · `a_g` · the Fisher bound). |
| `handoffs/HANDOFF_PROTEIN_AND_COMPOSITION.md` | Dated 2026-07-12; **12+ rows have landed since**. Superseded on the estimator (MH-115), the decoy control (MH-127/130/132), and the null calibration (MH-123/124). |
| `PROGRAM_FORWARD_BOARD.md`, `WHATS_NEXT.md`, `LEARNED_MODEL_WHATS_NEXT.md` | **Three competing forward docs, all predating MH-133…137.** All stale on the CN instrument. Consolidation pending. |
| `LANDSCAPE_REPORT.md` | Its per-edge FDR class counts rest on the **uncalibrated null** (MH-123/124). No banner yet. |
| `presentation/` F25 + `talk.pptx` | `talk.qmd` prose was fixed 2026-07-12 and is honest; **the F25 figure beside it is from 2026-06-24** and shows the old numbers. Deck needs re-rendering. |
