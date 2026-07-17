# isomiR-aware modeling — plan (per-sample seed attribution → corrected family predictors)

> **Goal:** specify how to move from the pooled isomiR QC (MH-96) to a **per-sample, seed-resolved family predictor**
> that routes each arm's reads to the family its seed actually targets, and propagate that correction through every
> submodule that builds a family pool `X_fam` (coupling, attribution, CN instrument, pressure, discovery).
> **What belongs here:** the per-sample estimation, the corrected-predictor math, the propagation map, the validation/
> gate, and the scope. NOT results (→ VALIDATION/registry), NOT the QC taxonomy itself (→ `within_family.py` / MH-96).
> **Update trigger:** when a phase lands or the design changes.
> **Sync-partner:** `learned/within_family.py` (the QC + `seed_shift`/`shifted_seed_family`), `families.collapse_by_family`
> (the predictor choke point), `CN_INSTRUMENT.md` / `LEARNED_MODEL_CHANNEL_FUSION_DESIGN.md` (the `X_fam` consumers).

> **Status (2026-07-12): Phase 0 + 1 + 2 BUILT; Phase 4 general gate RUN → the refit stays DEFAULT-OFF for coupling.**
> `seed_composition` (`c_{s,m,f}`) is cached; `collapse_by_family(isomir=True)` is the opt-in corrected predictor; the
> general coupling gate over the seed-heterogeneous families confirms the miR-29 verdict. **Scored where it ACTS** (user-flagged
> — the 656-edge mean −0.001 diluted it): on the **59 families that materially jump to another known family**, the RE-ORDERED
> edges move **real & material (mean |Δρ| 0.047–0.057, up to 0.14/edge) but BIDIRECTIONALLY** (signed mean ≈0 — tightens some,
> loosens others) ⇒ NOT a systematic coupling improvement. So Phase 3 (propagate) is correctly NOT triggered (adopt only if
> coupling improves); the deliverable is a genuine per-edge **re-attribution** (attribution/discovery), not a coupling lever.

---

## 0. Why (the finding that forces this)

The seed-family collapse (§8: same-seed arms → one `β_f`) assumes every read of an arm carries the family's seed. The
isomiR QC (MH-96, raw TCGA/CPTAC isoform-seq) shows this fails for a real minority: **~a quarter of expressed multi-
member families have a member whose dominant 5'-isomiR moves the seed off the family** (`shifted_seed_family`: 26+/122
arms' shifted seed lands in a *different* known family — e.g. **miR-29b-3p ~19% (per-sample) carries the miR-767-5p
seed**, a headline family). So a fraction of a heterogeneous arm's reads target a *different* genome-wide targetome and
should not count toward `X_fam,f`.

**Two facts from the per-sample diagnostic (683 TCGA samples) that shape the design:**
1. **The shift fraction VARIES per sample** — miR-29b-3p 0.03–0.50 (mean 0.19, sd 0.08); miR-29c-5p 0.20–0.67. A single
   static number cannot represent it.
2. **Pooling is BIASED** — miR-181b-3p / miR-374a-3p read ~0.48 pooled but ~0.01 per-sample (pooling conflates inter-
   sample 5'-dominant disagreement into false intra-arm heterogeneity). ⇒ **per-sample is more correct for the QC itself,
   not only for the model.** (Phase 0 fixes this; some MH-96 flags are pooled artifacts and drop out.)

**When the correction actually moves the model (the gate-relevant subtlety):** a *constant* per-arm fraction is absorbed
by the model's per-arm z-scoring/`β` (a pure rescale). The correction changes coupling only because the fraction **varies
per sample** (fact 1) — so per-sample is not a nicety, it is the whole reason the correction is non-trivial. For multi-
member families it also re-weights *which member* dominates the pool (member-specific fractions), which shifts the pool's
sample pattern even before per-sample variation.

---

## 1. Phases

### Phase 0 — per-sample QC (BUILT, `within_family.seed_shift`)
Recompute `seed_shift` per sample then average: (1) pooled reads → per-arm **global-mode 5'** (the canonical-seed
reference, robust to miRBase-version offset); (2) per file, the fraction of that sample's within-locus (±6 nt) mature
reads NOT at the mode; arm summary = **mean ± sd of per-sample fractions** over samples with ≥ `min_reads` (20). Fixes
the pooling bias. Re-run `family_structure` / `shifted_seed_family` on the full TCGA cohort → the corrected taxonomy.
Perf: `_parse_mature` is `usecols` + vectorised, files parsed in **parallel** (`ProcessPoolExecutor`, 8 cores) → full
cohort ~seconds.

**COHORT MATCHING (cross-cutting — the manifest is not the model cohort).** The 1207 isoform files are **1096 Primary
Tumor + 104 NAT + 7 Metastatic**, and among primaries some are **Normal-like PAM50** which the model EXCLUDES. So:
- The manifest carries `id`→`cases.submitter_id` (participant)→`sample_type`; `_file_meta()` maps it and `_iso_files`
  filters to **Primary Tumor** by default (NAT/Metastatic dropped) and tags each file with its participant.
- **Normal-like subtype is inherited, not re-filtered:** the per-sample matrix (Phase 1) is keyed by participant and
  **intersected with the model's `X` columns** (already primary-tumor + normal-like-excluded via `data_loaders`
  EXCLUDE_NORMAL_LIKE) — so matching on participant ID automatically enforces the model's exact cohort.
- (State-stratified seed-shift — does the isomiR fraction differ NAT vs tumor? — is a downstream extension, using the
  `sample_type` tag; not the default.)

### Phase 1 — per-sample seed-composition matrix (BUILT, `within_family.seed_composition`)
Per (participant × arm) `canon_frac = c_{s,m}` = fraction of arm m's reads in participant s carrying the arm's
canonical (global-mode) seed — the **per-sample** object that corrects `x(m,s)`, NOT the `seed_shift` average (a
constant per-arm fraction is absorbed by z-scoring; the per-sample variation is the whole point — e.g. miR-29b-3p
`canon_frac` ranges 0.49–0.98 across tumors). Low-read `(s,m)` cells are **EB-shrunk** toward the arm's cohort
canonical fraction (pseudo-count, weight ∝ read depth), so noise isn't read as a seed switch. **Keyed by participant
→ intersect with the model's `X` columns** (verified: 1078 of the model's 1079 participants match).

**JUMP leg BUILT** — `seed_composition` returns the full `c_{s,m,f}` (long: participant, arm, target_family, frac,
is_canonical): each arm's reads are split over destination families by their 5'-offset seed (canonical at offset 0;
JUMP families where a shifted seed = another family's, via `_seed_family_maps`; `orphan` where the shifted seed matches
no known family → dropped). So a family **gains** an arm's jumped reads and **loses** the reads that jump away.
Verified: **miR-29b-3p → {miR-29-3p 0.81, miR-767-5p 0.16, orphan 0.02}**. Low-read cells Dirichlet-shrunk toward the
arm's cohort family-distribution. **Cleanups DONE:** per-cell fracs **renormalised to Σ=1** (mean 0.9999); and a real
bug fixed — `_iso_files` was keyed by (participant, path) so **multi-aliquot participants were double-counted** (per-cell
sum hit 3.0), now **one file per participant** (1078, matching the model's mean-collapse). The obscure low-read `MIMAT*`
arms' 0.5 splits are correct shrinkage (34% are below-floor) → excluded downstream by the expression filter, not a bug.

**PHANTOM-REGULATOR SCAN (`phantom_scan`, bounds MH-97).** Systematic generalisation: an HE regulator that is
**unexpressed** (median RPM < floor) **whose seed is manufactured as a 5'-isomiR by an abundant arm** (a jump into its
family) — its curated edges can only be that arm's isomiRs here. **Result: only miR-767-5p qualifies** (1 RPM; donors
miR-29a-3p 4105 RPM / miR-29b-3p 500 RPM; 10/11 shared ECM targets). So the phenomenon is **rare and specific**, not
pervasive — it validates *and* bounds the miR-767 finding.

### Phase 2 — isomiR-corrected family predictor (the core change) — ✅ BUILT (opt-in, 2026-07-12)
`families.collapse_by_family(..., isomir=True, comp=None)` — DEFAULT OFF (the 22 callers are unaffected):
- now:      `X_fam,f[s] = log2( 1 + Σ_{m∈f} (2^{x_{s,m}} − 1) )`
- corrected: `X_fam,f[s] = log2( 1 + Σ_{m: c_{s,m,f}>0} c_{s,m,f} · (2^{x_{s,m}} − 1) )`
An arm contributes to its canonical family (weight `c_canon`) **and** the family its isomiR jumps to (weight `c_jump`);
novel-seed reads are dropped. **Arms below the isomiR read-floor (no `c` row) keep the canonical assumption** (full RPM to
their own family), so a family with no isomiR coverage reduces EXACTLY to the plain pool (verified). `c_{s,m,f}` from
`families._seed_comp()` — cached from `within_family.seed_composition` → `output/learned/seed_composition.tsv.gz` (1.89 M
rows, 1078 participants; parsed once). Verified: COL1A2's miR-29 family = only `miR-29c-3p` (canon 0.99, low jumper) ⇒ tiny
correction, as designed (member-dependent).

### Phase 3 — propagate to the `X_fam` consumers  ⚠ THIS CHANGES EVERYTHING
`X_fam` is a **foundational predictor**, so correcting it is **not localized** — it propagates to **coupling (OOF ρ),
attribution (`canonical_M`, identity/magnitude), discovery (which orphans enter), the CN instrument (`X_fam` is the
exposure/mediator in `exclusion` / `family_multi_iv` / `channel_cn`), the Bayesian Shapley (MH-93/94/95), and pressure —
the entire spine.** `collapse_by_family` is the choke point (most consumers inherit it), but **audit the modules that
build `X_fam` themselves** — chiefly `instrument.py` — and route them through the corrected pool. Grep target: `2.0 **`
/ `log2(` over family members.

**Because it changes everything, the discipline is:**
- **Default OFF, opt-in (`isomir=True`), gated** — adopt only if Phase 4 shows it *improves* held-out coupling, not
  merely *changes* it. The QC knowledge (which arms are mis-collapsed) is banked regardless of whether we refit.
- **Headline re-validation is mandatory** — re-run the established verdicts under corrected `X_fam` and diff them:
  coupling OOF, the CN exclusion/β, the Shapley identities, discovery FDR, and especially **miR-29→ECM** (a project
  headline AND the most-affected arm — miR-29b ~19%→miR-767); check whether the correction *sharpens* miR-29's coupling
  (cleaner instrument) and whether miR-767's targets *gain*.
- **Expected targeted, not global (a hypothesis to confirm):** median per-sample shift 0.056 ⇒ most families untouched;
  the refit should move only the ~17 seed-heterogeneous + ~34 dose-partitioned families. The gate must verify this
  (a global shift would signal a bug or a scale error, not biology).

### Phase 4 — validation & gate (does it earn its place?)

**GATE RESULT — miR-29 (the headline family), 2026-07-10.** The correction is REAL knowledge but a COUPLING WASH on the
donor side, and its value on the recipient side is ATTRIBUTION, not coupling magnitude:
- **Donor (miR-29 loses ~20%):** both miR-29a-3p (~23%) and miR-29b-3p (~16%) jump to miR-767-5p. On miR-29→ECM
  coupling: **mean Δρ +0.004, tighter on 3/7 targets = a WASH** — the per-sample down-weight is ~absorbed by z-scoring
  (same "real-but-coupling-immaterial-at-this-depth" pattern as the t-likelihood / CN channel / gate SEs).
- **Recipient (miR-767 gains):** miR-767-5p is a **PHANTOM regulator — median 1 RPM** (vs miR-29a 4105 / miR-29b 500), so
  its own signal is meaningless noise (`rho_unc` is the rank-corr of a ~constant near-zero variable). Its seed **GCACCAU
  = miR-29's +1 isomiR seed exactly**, so its 11 curated targets ARE miR-29's ECM set (**10/11 shared**; only PMP22 is
  miR-767-only), and the jumped-in miR-29 reads (~1000 RPM) *dominate* the corrected pool ⇒ **`rho_corr` ≈ miR-29's own
  abundance coupling with each gene.** That coupling is heterogeneous after composition adjustment (per MH-51): it
  **tightens** where miR-29 repression is still detectable (**SERPINH1 +0.11 → −0.24**, COL4A2/MMP2/COL3A1/SPARC) and
  **loosens** where it is composition-attenuated (**FBN1** — a real miR-29 target, not "noise") or genuinely non-miR-29
  (**PMP22**, the one miR-767-only gene). mean Δ −0.01, 5/11 tighter — the mix is miR-29's own composition-sensitive
  coupling, not random noise.
- **VERDICT (tightened to the defensible claim):** the **strong, clean claim is ATTRIBUTION** — miR-767-5p is
  **unexpressed in BRCA (1 RPM)**, its seed is miR-29's +1 isomiR seed, so its curated seed-edges **can only be
  miR-29-5'-isomiR-supplied here** (the model credits a phantom; the molecules are miR-29). This is a sequence+expression
  fact, robust. The coupling test only *illustrates* it (noisy, n=1-RPM-dominated, composition-sensitive) — it is **not**
  a coupling-improvement claim. Net: the isomiR correction's payoff is **attribution/discovery, not a coupling-magnitude
  lever** (same place the CN channel's value landed) → stays **default-OFF for the coupling model**; the discovery
  (phantom regulator, ~20% miR-29 mis-collapse) is the deliverable. Reproduce: `within_family.seed_composition` + the
  miR-29/miR-767 gate (scratchpad).

**GENERAL GATE RESULT — all seed-heterogeneous families, 2026-07-12 (`isomir_coupling_gate.tsv`).** Ran the coupling
Δρ (partial Spearman(`X_fam,f`, target | C), plain vs `isomir=True`) over **253 families with an expressed member of mean
canonical frac < 0.85 · 656 family-target edges**: **mean Δρ −0.0014, mean |Δρ| 0.016, only 41/656 (6%) edges exceed
|Δρ|>0.05 and those scatter BOTH directions** (no systematic tightening). No family shows a robust coupling gain (the
largest per-family means are n=3 noise). ⇒ **the miR-29 wash GENERALIZES — the correction is a systematic coupling wash**
(the recurring "real-but-immaterial-at-n≈1000" pattern), so the refit stays **default-OFF for coupling**; Phase 3 is not
triggered. The mechanism (`collapse_by_family(isomir=True)`) + `seed_composition` are banked for the attribution/discovery use.

**⚠ AGGREGATION FIX (2026-07-12, user-flagged): the mean over ALL edges DILUTES the effect — score it where it ACTS.**
Averaging Δρ over 656 broad edges (mean −0.0014) understated it, because most families barely re-order `X_fam` (and a
rank correlation only moves where the correction changes sample RANK order — a constant re-weight is absorbed). Refocused on
the **59 families with a member that materially jumps to another KNOWN family (mean jump ≥10%) · 211 edges**, stratified by
`rank_change = 1−Spearman(X_plain, X_corr)` (`isomir_gate_jumpfamilies.tsv`): on the RE-ORDERED edges the effect is **real and
material — mean |Δρ| 0.047 (rank>0.05) → 0.057 (rank>0.1), up to 0.14 per edge** (miR-520→CDKN1A −0.14 tighter, miR-520g→SMAD7
−0.06; but miR-133a→FUS +0.12 looser, miR-520g→VEGFA +0.08). **BUT the signed mean stays ≈0 (−0.003 → −0.001) — the moves are
BIDIRECTIONAL** (tightens where miR-repression survives composition adjustment, loosens where attenuated/non-target), so it is
**NOT a systematic coupling improvement**. ⇒ the verdict holds but is sharper: the correction is a genuine per-edge
**re-attribution** (real, not a wash-to-zero where it acts) that does NOT net-improve held-out coupling ⇒ **default-OFF for the
coupling model, but a real attribution/discovery tool** — exactly the miR-29 pattern, now shown across 59 jump families.

**Gate spec (for other families):**
- **Coupling:** does the isomiR-corrected `X_fam,f` couple *better* with family f's targets (removing non-targeting reads
  should tighten it), and does the recipient family (miR-767) *gain* coupling with its targets? **Gate on held-out OOF
  improvement** (patient CV), as every change in this program. If no improvement → keep as a QC/attribution flag, not a
  model change.
- **Null:** shuffle the per-sample fractions across samples → the correction's effect must collapse (else it's leakage).
- **Scope:** only the genuinely per-sample-heterogeneous arms (a couple dozen) are affected; the bulk (~median shift
  0.04) is unchanged → a targeted correction, cheap to gate.

---

## 2. Caveats to carry
- **Read-fraction ≠ loaded-fraction.** Repression depends on which isomiRs are RISC-loaded; 5'-isomiRs may load
  differently. The read-fraction split is first-order; the loaded refinement ties to `ago_loading` (a `c_{s,m,f}` scaled
  by loaded fraction). Flag, don't block.
- **Functionality.** The correction assumes the shifted isomiR is active on its seed's targets; the Phase-4 coupling test
  is what confirms that per arm.
- **Intra- vs inter-sample.** The per-sample `sd` distinguishes an arm that is a *mixture within* every sample from one
  that *switches* dominant 5' across samples; both reduce `c_canon`, but only the latter creates cross-sample targeting
  differences. Report both.
- **Direction knownness.** Only ~21% of shifted seeds map to a *known* family; the rest are orphan seeds (dropped from
  the canonical pool but not re-attributed). The orphan fraction is a pure down-weight of the canonical family.

---

## 3. Estimand honesty
The corrected `X_fam,f` estimates *the abundance of reads carrying family f's seed*, per sample — a cleaner instrument
for "family-f regulatory input" than total arm abundance. It does **not** change the §8 principle (same-seed → identical
targeting); it enforces it, by ensuring only same-seed reads enter the same-seed pool.
