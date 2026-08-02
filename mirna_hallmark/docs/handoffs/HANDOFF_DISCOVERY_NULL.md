# Handoff — Discovery lane, the site-free null, and the FDR (MH-154 / MH-155)

> **Goal:** hand off the 2026-07-17 discovery-calibration arc — what was built, what was measured, what is
> WRONG in the record, and what to do next.
> **What belongs here:** the discovery lane's null/FDR/lane state and its open work. NOT a registry finding
> (deferred by decision — no MH discovery *finding* row until the convergent-evidence framing is decided).
> **Update trigger:** when the correlation-matched null, the size-routed test, or the A1∩chimeric table lands.
> **Sync-partner:** `PROGRAM_FORWARD_BOARD.md` (topic D Discovery) — keep the one-line status there in sync.

Read first: `docs/STATE_OF_PLAY.md` (edge-existence + attribution axes), registry rows **MH-123, MH-154**,
and `eval/site_free_null.py`. Code home: `learned/discovery.py`, `eval/site_free_null.py`.

---

## 0. The one-paragraph state

Per-edge and per-family miRNA→target **discovery in bulk TCGA is empty under an honest null** — computed the
right way now (empirical heavy-tailed null → Simes within family → BH across families), not via BY's misdirected
worst-case. The real signal is **set-level** (candidate set shifted to median null_z ≈ −1.6; ⭐ the FAMILY lane is
**−1.32** on the corrected co-expression-matched null, MH-197 — the ARM lane is unchanged) and **convergent-evidence** (coupling tracks independent chimeric/K_D evidence). The defensible
deliverable is **A1∩chimeric**: 296 edges with bulk coupling + composition-robustness + orthogonal physical duplex,
led by miR-18a→{STAM2, KIF3B, MAP3K1, NEDD4}. Two correctness fixes are owed (below).

---

## 1. What MH-154 settled (committed, in the registry)

- **The permutation null IS the theoretical null.** `coupling_permutation`'s 2000 Freedman–Lane draws reproduce
  the t-null to 3 s.f. (σ₀=0.0309 vs 0.0311, R²=0.9997). A permutation destroys the *unmodelled* structure that
  inflates the real null, so it cannot rescale it. **Porting permutation to the learned model buys nothing.**
- **BY corrected the wrong defect.** MH-45's dependence-robust count (BY, ρ<0) = 33.3%, inside the 35.3% rate at
  which *impossible* edges clear the same gate. BY handles dependence; the defect was SCALE.
- **The calibration standard is `eval/site_free_null.py`** — MH-123's Efron site-free null, rescued from a dead
  `/tmp` scratchpad into a repo module; reproduces MH-123 to the digit.

## 2. What MH-155 built in `learned/discovery.py` (committed)

- **Site-free null fitted IN-LOOP per gene** on that gene's exact `Cext` (so the null's estimand = the
  candidates': same C, same he_agg, same sample support). Verified bit-equal to the real per-arm path (5.7e-12).
- **he_agg = canonical Gibbs β** (`readouts_arm_edges.tsv` `beta`/`beta_deconv` — the ATTRIBUTION block of `canonical_card`; was the orphan `arm_attribution_card.tsv`, deleted 2026-07-18), applied on β's z-scored
  `_prep` scale — NOT the §6b-retired adaptive lasso (user-caught). A stronger, canonical control.
- **batch OFF** in `Cext` (user-asked ×2): canonical C = [CPE, target_cn, mal_prolif] carries no batch, and
  batch is shared structure already IN the site-free null (σ₀ 0.094→0.087, ~7%, no calibration gain).
- **FDR hierarchy** (calibrate): `q_family_emp` **PRIMARY** = empirical arm-p (heavy-tail) → Simes within
  (gene, seed-family) → BH across families. `q_empirical` (arm) secondary. `q_by_arm` = worst-case SENSITIVITY
  bound only. **All lanes/universes: q_family_emp<0.05 = 0.**
- **Flags:** `same_seed_he` (paralogue of a curated regulator), `he_max_corr` (near-duplicate of a specific
  known regulator), `no_he_gene` (lane 2).
- **Evidence attach:** `ts_mag`, `chimeric_wt`/`chimeric_src` (Manakov eCLIP + TarBase CLASH — arm-resolving;
  ⚠ NO breast chimeric in any source), `ev_classes`/`ev_npmid` (ledger, PMID-deduped across mirTarBase+TarBase).
  scanMiR/TargetScan are PRECONDITIONS, not attachments.
- **Two new lanes:** LANE 2 (`universe='hallmark'`) — the 2,792 Hallmark genes with no curated HE arm; LANE 1
  (`scan_families`/`run_families`) — whole seed families uncurated for the gene, as a dose aggregate.
- **6.5× speedup** (measured): the epicenter was the candidate loop's per-arm sklearn-residualise + spearmanr,
  NOT `assemble_gene` (0.02s warm). One QR([1,Cext]) residualiser/gene applied as a matmul; `OMP_NUM_THREADS=1`
  kills BLAS oversubscription. Genome-wide 20→3.1 min (arm HE), ~10 min (arm+family Hallmark, co-running).

Outputs: `output/learned/{discoveries.tsv (Hallmark arm), discoveries_family.tsv, discovery_sitefree_null.tsv,
discovery_family_null.tsv, discovery_fall_diagnosis.tsv}`.

## 3. The measured findings (honest readouts — NOT per-edge FDR)

- **The site-free null is HEAVY-TAILED, not Gaussian:** 2.7× the Gaussian tail at z=−3, **8.4× at z=−4**. So
  BH/BY on a fitted-Gaussian p under-counts tail false-positives — BH's 654 "survivors" are that artifact. The
  empirical FDR (heavy-tailed null-draw reference) gives 0, and on the FAMILY lane it caught **2 tail
  false-positives BY had missed** (BY on Gaussian p kept 2; empirical kept 0). BY was giving *wrong* answers, not
  just conservative ones.
- **PRDS is plausible:** 74.4% of candidate-arm pairs are positively correlated (mean +0.117) — so BH is valid
  IF the null were Gaussian; it isn't, so the empirical FDR is the honest primary regardless.
- **Set-level shift is stable and real:** median null_z ≈ −1.6 on every lane and universe. ⭐ **CORRECTED
  2026-08-02 (MH-197): the FAMILY lane is −1.32** once the null's pseudo-families match the candidates'
  internal co-expression; the ARM lane is untouched (single arms have no internal ρ̄).
- **Independent evidence concentrates in stronger coupling** (the concept-level claim): chimeric-duplex-present
  edges couple stronger (MWU p=4.6e−8), ledger PMIDs (p=9e−5), scanMiR strength (ρ=−0.078, p=2.6e−5);
  TargetScan magnitude does NOT (p=0.27) — a real asymmetry (K_D tracks coupling, site-count doesn't).

## 4. ⛔ CORRECTIONS TO THE RECORD (do not propagate the wrong versions)

- **Commit `3bc2bdc`'s message is WRONG where it says "real seed families are NOT more internally co-expressed
  than random (0.153 vs 0.115, p=0.82)."** That comparison was diluted (averaged over hundreds of
  undetected-arm families) and used a null that was NOT dosage-matched. **REDONE correctly** (dosage-matched
  null, detected multi-arm families, stratified): real within-family ρ̄ **0.42** vs matched null **0.11**,
  Wilcoxon **p=1.2e−6**, strongest in high-dosage families. ⇒ **the family-lane null IS mis-calibrated** (it sums
  ~independent arms; real families sum arms correlated at 0.42). The mis-calibration is anti-conservative (null
  too clean), so the family lane's 0 survivors understates how null it is — but it must be fixed.
  ✅ **FIXED 2026-08-02 (MH-197).** Independently re-derived (0.421 vs 0.118, MWU p=5.6e−22) and then measured
  **on the discovery lane's OWN candidate population, where it is WORSE than recorded here: real ρ̄ +0.488 mean
  / +0.543 median vs random draw +0.118.** `_corr_matched_group` now draws each pseudo-family to match the ρ̄ of
  the **specific candidate cell it stands in for** (not a global average). Paired A/B on identical genes, cells
  and seed, multi-arm draws only (n=2,946): **null sd 0.1299 → 0.1554 (1.196×), robust-sd ×1.249, Levene
  p=1.5e−21** ⇒ the old null was **20–25% too narrow**. Reversible via `scan_families(corr_matched=False)`.
- **My within-family "refutation" was a false negative from a bad measurement.** The user caught it. Lesson:
  when comparing real vs a matched null, match on the confound axis (dosage) AND restrict to the relevant
  (detected, high-dosage) stratum — an unweighted mean over a diluted population hides a concentrated effect.

## 5. Simes adequacy + BH-across-families (the reasoning, for whoever extends the FDR)

- **Simes is VALID (conservative) under the positive within-family dependence** (Sarkar, PRDS/MTP2) — it will not
  inflate. But it is **UNDERPOWERED on exactly the important families:** a single strong arm in an m-arm family
  is penalised ~m (this is the A2 "dilution" class, 400 candidates), and *distributed*-dosage signal is caught
  better by the **dose-aggregate** than by combining weak arm-p's. Because the high-dosage families are both
  multi-arm and internally correlated (§4), the aggregate is the right unit there.
- **BH-across-families > old arm-BH** (654→0): fixes both the paralogue double-counting (collapse) and the scale
  (empirical p). ⚠ Residual dependence remains: same-cluster seed families (miR-17~92 → miR-18, miR-19,
  miR-17/20/93, miR-25/92 are co-transcribed) are correlated, so "across families" isn't fully independent — BH
  stays valid only under PRDS (positive co-transcription supports it).

## 6. THE DELIVERABLE — A1∩chimeric, now SITE-LEVEL (convergent evidence, not per-edge FDR)

The `discovery_fall_diagnosis.tsv` "A1" class = deconv-robust candidates that pass a *standard* BH but fail the
honest FDR. **69% of A1 (296/429) carry an independent physical duplex** (edge-level chimeric). These have
triply-convergent support — bulk coupling + composition-robustness + orthogonal chimeric assay — a claim *against
the null by concept*, since the honest per-edge FDR is unreachable. Flagship: **miR-18a→{STAM2, KIF3B, MAP3K1,
NEDD4}**, all with Manakov duplexes; **NEDD4 confirmed by three chimeric sources** (Manakov + TarBase chimeric +
qCLASH).
> ⛔ **CORRECTION 2026-08-02 (MH-198): the flagship quartet is EDGE-level and one member does NOT survive the
> site-level sharpening below. `KIF3B` has `site_manakov = False`** (edge-level CLIP only — `site_clip_any=True`,
> `best_type=7mer-A1`, `partial_coupling=−0.588`), so it is in the 296-edge edge-level set and **NOT in the
> 157-edge site-level gold set.** STAM2, MAP3K1 and NEDD4 hold. Quote the quartet only at edge level; at site
> level it is a **trio**.

**⭐ SHARPENED TO SITE-LEVEL (MH-155, `discovery_site_evidence.py`).** The site ladder now runs on the discovery
orphans: each candidate carries whether its *predicted 3'UTR site* physically coincides with a chimeric duplex,
not just whether a duplex exists somewhere on the gene. **Site-level Manakov overlap tracks coupling at MWU
p=1.9e−20** (vs edge-level 4.6e−8 — 12 orders sharper), and the coincidence rate climbs monotonically up the
site-confidence ladder (7mer-A1 Manakov 9% → 8mer+conserved+3'-supp 27%) — a gradient a site-free arm cannot
produce. **A1 sharpens 296 (edge-level) → 157 (predicted site coincides with a duplex)** — the gold set. Output:
`output/learned/discoveries_sitelevel.tsv` (`site_manakov`, `site_clip_any`, `best_type`, `n_sites`).

### 6b. ⭐ THE GOLD SET NOW HAS A PRODUCER, AND TWO CAVEATS THAT MUST TRAVEL WITH IT (MH-198, 2026-08-02)

Until now the 157 existed **only as a hand-merge** of `discovery_fall_diagnosis.tsv` × `discoveries_sitelevel.tsv`
— no producer, no provenance, no attached caveats. *That is precisely how the five literature ground-truth sets
rotted (MH-196).* Built now by **`eval/discovery_gold_set.py`** → `output/learned/discovery_gold_set.tsv`
(+ sha256/provenance JSON). Shape: **157 edges · 90 genes · 18 arms · ⭐ 11 SEED FAMILIES.**

⚠ **(1) IT IS NOT 157 INDEPENDENT FINDINGS — the set is one cluster.**

| edges | seed family |
|---|---|
| **96 (61%)** | miR-17-5p/20-5p/93-5p/106-5p/519-3p |
| 24 | miR-19-3p |
| 10 | miR-25-3p/32-5p/92-3p/363-3p/367-3p |
| 9 | miR-18-5p |
| 9 | miR-130-3p/301-3p/454-3p |
| ≤3 each | 6 further families |

**139/157 = 88% belong to the miR-17~92 polycistron and its paralogue clusters** (miR-106b~25, miR-106a~363),
which are **co-transcribed** — their doses are not independent draws. ⇒ **the honest unit is ~one oncogenic
cluster's realized target set, not 157 discoveries. Always quote the 11 families beside the 157 edges.**

✅ **(2) WHAT IS GENUINELY STRONG:** **94/157 sit on genes with NO curated regulator at all** (`no_he_gene`) —
novel gene territory, not curation re-derived. And only **3/157** are same-seed paralogues of a regulator already
curated for that gene, so the set is not curation leaking back in via family membership.

⚠ **Re-derived the concept-level statistic:** site-coincident edges couple **−0.2619 (n=1,409) vs −0.2358
(n=5,408), MWU p=2.28e−27** over the full site-level table. §6's **1.9e−20 was a narrower scope** — both are
real; cite the scope with the number.

## 7. OPEN WORK (suggested next steps, roughly ordered)

1. ✅ **DONE 2026-08-02 (MH-197) — ~~Correlation-matched family null~~.** Implemented as the "lighter version"
   below (greedy co-expressed site-free groups), but matched **per candidate cell** rather than to a global ρ̄.
   See §4. **The verdict did NOT change and could not have** — the old null was anti-conservative, so widening
   it can only keep 0 survivors at 0 (confirmed: 0 → 0 under `q_family_emp`, `q_empirical` and `q_by_arm`).
   ⭐ **What DID change is the headline robust readout: the set-level shift attenuates, median `null_z`
   −1.6 → −1.32.** The lane's deliverable is now measured against a null that matches the candidates'
   co-expression structure.
   ⚠ **The `scrambled-seed` gold standard was scoped and NOT built, deliberately — it does not do what §7 said.**
   Scrambling a seed changes which GENES an arm is predicted to hit, not the arm's expression; you still score a
   REAL arm's expression against the gene, and that arm can still be a genuine non-canonical regulator. So it
   matches SITE statistics but **does not remove MH-123's leak.** The arm-lane leak is a separate defect running
   in the OPPOSITE direction (real signal inside the null ⇒ null too WIDE ⇒ conservative) and remains open; it
   needs an object with real expression that provably cannot bind, which the scrambled seed is not.
   ⚠ Also rejected on measurement: **real multi-arm families as donors** (structurally exact ρ̄) — only **68**
   exist detected, **9 of size ≥4 and 1 of size ≥7**, against 177 candidate cells needing size 7.
2. **Unify the two lanes by family size.** Route multi-arm families through the dose-aggregate (family lane),
   singletons through the arm test — because Simes underpowers the correlated multi-arm families (§5). One test,
   right unit per family.
3. **A1∩chimeric convergent-evidence table** → the candidate registry finding (framed as convergent evidence,
   NOT per-edge significance). Decide the framing with the user before writing an MH row.
4. **⭐ Subtype-stratified (PAM50) run (user-asked).** Pooled-cohort washout can hide subtype-specific coupling;
   `scan_all`/`scan_families` accept a gene list but need a subtype sample mask threaded through `assemble_gene`.
5. **Family-lane evidence attach** — pool members' chimeric/ledger (currently arm-only).
6. ~~Dose-response test — coupling vs site-count / K_D quality.~~ **DONE (MH-155):** `discovery_site_evidence.py`
   — site-level chimeric overlap tracks coupling (p=1.9e−20) and climbs the site-confidence ladder monotonically.
   Remaining: formalise the site-count → coupling dose-response as a single test and add K_D (scanMiR) as a
   parallel axis (scanMiR already tracks coupling, ρ=−0.078; TargetScan magnitude does not).
7. **Revive the rest of the `site_ladder`** — `utr_seed_scan` (site counts), `validate_ladder_experimental`
   (the ladder is now runnable end-to-end; the L5 rung's POSTAR is optional, MH-149). It is the natural home for
   the dose-response and the site-level convergent-evidence deliverable.

## 8. Data / infra notes for the next session

- **Canonical Gibbs β lives in `output/learned/readouts_arm_edges.tsv`** (`beta` core, `beta_deconv`; = the ATTRIBUTION block of `canonical_card.tsv`. Was the orphan `arm_attribution_card.tsv`, deleted 2026-07-18 — no writer, superseded by `readouts.run(level="arm")`);
  100% HE coverage on the discovery universe. Apply on the `attribution_eb._prep` scale (z-scored, C-residualised).
- **`OMP_NUM_THREADS=1` for every parallel discovery run** — workers were BLAS-oversubscribing (~130% CPU each).
- **Smoke tests are broken on a pre-existing dead GENCODE path** (memory `smoke-test-stale-paths`), unrelated to
  this work — verify discovery changes with targeted end-to-end runs, not `tests.run_smoke_tests`.
- `chimeric_evidence.py` (Jul 9, parquet cache) is current; ⚠ NO breast chimeric exists in any source
  (Manakov=HEK293T, TarBase CLASH=other tissues) — cross-tissue corroboration only.
