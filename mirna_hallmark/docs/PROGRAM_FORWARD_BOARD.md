# Program forward board — the whole open-work map (`mirna_hallmark`)

> **Goal:** ONE consolidated view of every open forward thread across the whole subproject — edge existence,
> the decoy/claim-1 arc, attribution, cross-cohort validation, the discovery lane, the DCIS/EV lane, and infra —
> so nothing is lost between the detail docs. This is the *index of what's next*; the *depth* lives in the
> registry row named on each line.
> **What belongs here:** a one-to-three-line entry per open thread, with a status marker and a pointer to its
> `MH-##` row or home doc. NOT methods/rationale (→ their home docs), NOT results (→ `DISCOVERY_REGISTRY.md`),
> NOT done-and-closed items except as compact ✅/⛔ anchors that stop someone rebuilding a dead thing.
> **Update trigger:** when a thread is started/finished, or a new thread is opened anywhere.
> **Sync-partner:** `docs/STATE_OF_PLAY.md` (the current-state verdicts this board plans against; if the two
> disagree, STATE_OF_PLAY wins, and the registry wins over both).

**Last updated: 2026-08-03** · planning against the registry through **MH-207**.
⚠ Several 2026-07-17 items were found STALE on 2026-08-01 (already done) and are now marked ✅ —
verify a 🔨 before building against it.
⚠⚠ **AND AGAIN ON 2026-08-03: a sweep of §B + the discovery lane found FIVE more stale entries** — three
already DONE (`site_free_null` wiring · the per-GENE roll-up, hidden by MH-181's **rename** · the ENCORI/
POSTAR3 hole, closed by measurement) and two mis-stated in a way that would have sent work in the wrong
direction (POSTAR3 "not in this repo" — it is; "acquire it" — acquiring it would not have helped).
**The base rate of stale entries on this board is high. Verify before building, always.**
*Supersedes and replaces `WHATS_NEXT.md` and `LEARNED_MODEL_WHATS_NEXT.md` (both archived; both predated MH-133…137).*

Status: ✅ done (anchor only) · 🔨 immediate/flagged · ⬜ open · 🔬 investigate · ⛔ blocked/dead · 🚫 deliberate skip.

---

## 🔨 DO THESE FIRST

Highest value ÷ cost. ⚠ Two of the original four turned out to be claims that did not survive checking
(MH-38's 108 orphans were 3; E6's "live bug" was never live) — verify the claim before building the fix.

0. 🔨 **⭐ MODEL EXPANSION + REFIT (user-requested 2026-07-19, MH-166 external arc) — NOT YET RUN.** Fold the
   externally-validated missed edges into the HE universe and **re-run the Gibbs** (β is a joint per-gene
   posterior — a new edge re-decomposes the whole gene, so it needs a refit, not an appended weight). Priority
   candidates, in order: **(a)** the 21 gold∩external crown-jewels (internal coupling-discovery `dossier_gold`
   ∩ external functional/binding) — miR-17~92→TMBIM6/AHNAK/NEDD4/SERTAD3; **(b)** the 134 multi-evidence
   (`model_expansion_list.tsv`, ≥2 of functional/binding/coupling, not family-covered); **(c)** the double-breast
   set. Note some targets (TMBIM6, NEDD4, SERTAD3) have **0 current edges** → new *genes*, not just regulators.
   ⚠ prioritize by evidence QUALITY (seed type / functional lfc), not raw count — the list tilts to abundant
   arms. Files: `output/learned/realization/{model_expansion_list,missed_edges_*,breast_chimeric_discovery_*}.tsv`.
   After refit: read the new edges' β + updated gene pressure; treat coupling as REALIZATION not existence. **(→ MH-166)**

1. ✅ **DONE 2026-07-17 — `eval/buffa_validation` re-run: the "108 triple-validated orphans" ARE 3.**
   `triple_validated` **108 → 3**; **ECM/collagen among them 30 → 0**. Survivors (all guide arms):
   **miR-195-5p→PSMD7 · miR-181b-5p→IRS2 · miR-22-3p→STK39**; 1 neg-sig in Buffa. ⚠ **It required a code fix —
   the module had been DEAD since the `eval/` reorg** (`GEO_DIR` was a `__file__`-relative hop pointing at a
   non-existent `mirna_hallmark/data/`; the cache was at the repo root). **That is why it was never re-run.**
   HE legs reproduce exactly (+0.319 / 0.593 / 0.673) = the control. `talk.qmd` updated. **(→ MH-38)**
2. ✅ **DONE 2026-07-17 (verified 2026-08-01) — MH-38 / MH-55 / MH-73 / MH-74 ALL carry their staleness banners.**
   ~~all four still carry **S/R** strength tags with **no stale~~
   marker**. MH-38/55: input collapsed 594→23. MH-73/74: superseded by MH-76's frozen-panel test. Per the
   one-home rule, the retraction banner goes in *their* rows. **(→ CLAUDE.md doc protocol §4)**
3. ✅ **NOT A BUG — E6's `(live)` tag RETRACTED 2026-07-17 (MH-142).** Verified in code and data: every call
   site scores **orphans**, disjoint from `he_agg` by construction — **0 of dossier's 6,744 candidate pairs are
   HE edges** (vs 5,938 HE). Nothing self-partials in production. E6 **stands as a trap for CLASS-MATCHED
   comparisons** (its evidence: *"class-matched estimator comparison"*) — MH-124's own experiment pushed HE arms
   through this estimator and got a ~36× depressed HE rate **from its own design**. ⛔ **Do NOT remove `he_agg`**:
   it is what makes an orphan's coupling *partial on what curation captures*. Guard: `eval/_e6_production_check.py`.
   **(→ MH-142)**
4. ✅ **DONE 2026-07-17 (MH-140) — `readouts.share` relabelled AND the table gained a TRUE identity column.**
   `identity` = `shapley_identity` (LMG on R², m = bagged NNLS); `share`→`beta_frac` (demoted to MAGNITUDE);
   + `identity_deconv`, `identity_eq_magnitude`. **WHO ≠ HOW-MUCH in ~25% of multi-family genes** (agree 75.2%,
   n=819; was 99.35% when "identity" was argmax β). `m_nnls`=0 in 31.7% of families vs `beta`=0 in 0/5,117.
   213 genes (13.8%) now honestly report identity = UNDEFINED. Change is purely additive (5 columns
   bit-identical); `shapley_identity` optimized 8.2–9.7× output-identical. **(→ MH-140, MH-138)**
5. ✅ **DONE 2026-08-02 (MH-201) — `β(TCGA) → Buffa mRNA`.** `learned/eval/ood_cohort.py`, 1,364 genes,
   β frozen. **LEVEL transports** (β −0.0186 vs abundance −0.0022, p=6.2e-13; retention **+0.358** vs the
   **C-ABLATION** ceiling, which is the only honest denominator — dropping `target_cn`+CPE costs 33% of
   TCGA |ρ|). **RANK does NOT separate β from abundance** (+0.040 [−0.013,+0.093], p=0.13) ⇒ **inverts
   MH-174 for this boundary.** Survives the pre-registered MH-114 null (compartment-blind Δ −0.020/−0.023,
   same sign both orientations). ⚠ 43% of genes are single-family where β is **mathematically inert** —
   never pool them. ⛔ Corrected en route: Buffa is **not** a METABRIC subset; `target_cn` unobtainable;
   a plain symbol join drops 16.7% of β genes incl. **ZEB1** (=`TCF8`). — original item follows —
   ⬜ ~~**`β(TCGA) → Buffa mRNA`**~~ — the learned model's **only untested boundary is the cohort boundary**
   (MH-104: **~80% loss at the cohort boundary, ~0% at the layer boundary**), and Buffa (n=207, miRNA+mRNA) is
   the only independent-patient cohort in the repo that can test it. **β has never touched Buffa** — one `STUB`
   docstring in `learned/eval/__init__.py`; no `ood_cohort.py`. Inputs are on disk.
   ⚠ **Pre-register MH-114's compartment-orientation stratification** — both cohorts are bulk breast and share
   the CAF confound, so a clean replication proves nothing on its own. **(→ axiom 4; MH-104, MH-114)**
6. ✅ **DONE 2026-07-18 (MH-161) — REFUTED: admissibility is NOT a decoy-gap sharpener.** The pre-registered
   prediction (gap grows) failed — admissible vs all gap_core is **n.s. at every width** (pooled p=0.39); the
   complement guardrail fails (inadmissible not ~0); the only firmed signal is that admissible edges are more
   **cell-intrinsic** than inadmissible in **single-family genes** (n_fam=1 gap_deconv −0.0154 vs −0.0027, MWU
   p=0.023, held K=1→3 fold-seeds). **Keep as a diagnostic flag, not a filter.** Census: has_site 89.9%,
   expressed 63.9%, admissible 60.4%. `eval/admissibility_bench.py`; `output/learned/admissibility_{edge_tags,bench}.tsv`.
   Method gotcha: `decoy_bench --seeds` is order-only (fold seed is the real MC axis). **Pre-registration kept below for provenance.**
   ⬛ ~~**CURATED-EDGE ADMISSIBILITY — can this HE edge repress IN BREAST AT ALL? (user-proposed 2026-07-17)**~~
   An HE edge earns its place by a luciferase assay somewhere — not necessarily in **this tissue**. Two a-priori,
   Y-blind, model-blind gates ask whether it *could* operate here:
   * **SITE:** does the arm have a seed site in the target's 3'UTR? (TargetScan context++ ∪ scanMiR duplex ∪
     Poisson-gated 6mer — `decoy_bench._site_maps`, already built.)
   * **EXPRESSION:** is the arm expressed in breast — **tumour (TCGA) OR healthy (GTEx)**? (`_FLOOR = log2(11)`,
     RPM≥10, `card.py`; GTEx breast is ingested via `gtex_mirna_matrix.py`. **Healthy matters**: an arm silent in
     tumour but present in healthy breast is a *lost* regulator, not an inadmissible one.)
   **An edge failing BOTH cannot repress here by any modelled mechanism** ⇒ it is curation TRANSFERRED from
   another tissue/cell line, and it is sitting in the design diluting the real arm.
   **⭐ HALF OF THIS IS ALREADY MEASURED, and it lands:** MH-136 — **187 genes whose curated edges are ALL
   seedless give gap +0.0006, EXACTLY ZERO**, vs −0.0325 for 810 all-seeded genes. *No site → no repression →
   no coupling → no gap.* MH-132 independently found the 21% positive-coupling genes are **~9/13 positive in
   GTEx NORMAL BREAST too** ⇒ *"a **curation-transferability** problem, not a model defect"* — the same claim
   from the other side. **The expression half has never been tested.**
   **THE TEST (cheap — both artifacts exist):** tag every HE edge `admissible = has_site ∧ expressed_breast`;
   re-run `decoy_bench` on the admissible subset. **PREDICTION (register before running): the gap GROWS** —
   inadmissible edges are noise in the REAL arm, so removing them should sharpen the contrast beyond MH-137's
   −0.0119 / MH-139's −0.0129. **⛔ REFUTED 2026-07-18 (MH-161) — it did NOT grow; no sharpening at any width.
   Guardrail (d) mattered: the gate shrinks n_fam 2→1, and the residual signal is cell-intrinsic-vs-compartment
   in single-family genes only, not a magnitude gain.**
   ⚠ **Guardrails.** (a) This is **DIAGNOSTIC first, a filter second** — axiom 2a: flag, don't delete. A gate
   that removes 30% of the universe needs its own null before it changes a headline. (b) It is **NOT** the
   decoy's site-free rule (that defines the FAKE; this classifies the REAL). (c) Report on the gated set AND
   the complement — if the complement is not ~0, the gate is not what we think. (d) `n_fam` shrinks under the
   gate, and MH-147 showed the gap **scales with design width** — so re-check the width axis inside the gated
   set, or the gain is just width. **(→ MH-136, MH-132, MH-147; `ARM_EXPRESSION_FLOOR.md`, `SILENT_ARM_REMOVAL.md`)`**

6b. 🔨 **⭐ THE EXCLUDED SAME-SEED MATES — NEW 2026-08-19 (MH-248), AND IT IS THE CHEAPEST ROUTE INTO ITEM 0.**
   A gene's design holds a **median 0.67** of its seed family, and since same-seed mates bind the SAME site
   the absences are curation, not biology. Measured: **~187 excluded mates are above the expression floor,
   couple below −0.10, AND survive conditioning on their in-design partner** (median −0.1795 → −0.1382, so
   only 23% was the mate's signal; mates correlate a median 0.479, >0.8 in just 6%). ⭐ **The reason they are
   absent is the ASSAY LADDER, not a coverage hole: 86.6% of the actionable tail is `HT_ONLY`** — non-weak
   `ago_clip`/`chimeric`/`qpcr_rna` evidence with **no low-throughput functional assay** — and **LT_FUNC = 0**,
   so nothing properly evidenced is missing and there is no design-construction defect to fix.
   ⛔ **NOT explained by fame (75 PMIDs vs 60 for in-design arms), NOR by abundance.**
   **THE THREE OPEN LEGS, in order of cost:**
   * ⬜ **(a) COHORT TRANSFER — untested and the cheapest.** Everything above is TCGA-only. Repeat in
     **Buffa** (n=207, genuinely independent) and **CPTAC prospective** / **GTEx-NAT**. A mate that couples
     in two cohorts is a different object from one that couples in one.
   * ⬜ **(b) PUSH THE EVIDENCE AXIS.** Score each excluded mate on chimeric / degradome / CLIP **depth** and
     **site type**, and test whether HT_ONLY depth predicts its coupling. If it does, the ladder is
     informative and the exclusion rule can be re-drawn on evidence rather than on assay *kind*.
   * ⬜ **(c) FOLD THE SURVIVORS IN AND REFIT.** β is a joint per-gene posterior, so a new edge
     **re-decomposes the whole gene** — this needs a Gibbs refit, not an appended weight. **This is item 0
     (MH-166) with a concrete, pre-screened candidate list attached**, which is what that item has lacked.
   ⚠ **Pooled, adding the mates back is a WASH on the family aggregate** (median Δ −0.0002; 56% improve /
   41% dilute) — the gain is confined to NARROW designs and REVERSES at width ≥3. Do NOT sell this as a
   uniform improvement. The mean-reversion-immune extremes are CDH1/miR-17~92 (−0.059→−0.315), BCL10
   (−0.180→−0.342), BRCA1 (−0.012→−0.150), DLC1/miR-141-200a (−0.273→−0.408).
   * ⬜ **(d) THE DILUTION LEG — NEW, AND IT CONSTRAINS (c).** Adding mates back dilutes the family
     aggregate materially (Δ > +0.05) in **10.2%** of cells — not the 40.9% the raw sign split suggests,
     most of which is noise. ⭐ **The diluters are identified by ONE axis: their own coupling**
     (−0.0020 vs −0.0870, p=1.4e−22). **Abundance (p=0.61), fame (p=0.37) and EVIDENCE STATE (p=0.72) all
     fail to separate them** — the 2×2 is flat in abundance (uncoupled dilute 55.2% sparse / 52.5%
     abundant; coupled 31.8% / 32.0%). ⛔⛔ **So there is no evidence-based pre-screen for the refit**, and
     screening on coupling is selection-on-outcome. **(c) must therefore add all mates and accept the wash,
     or screen on an OUT-OF-FOLD coupling estimate — decide this BEFORE running the refit.**
     Dilution is family-structured and worth a per-family rule: miR-196-5p (4/4 cells dilute), miR-302
     (58%, median 8 added), miR-29-3p (58%), miR-34/449 (65%) vs gainers miR-99/100 (Δ −0.099),
     miR-146-5p (12.5%), miR-141/200a (25%), miR-200bc/429 (7.7%).
     ⭐⭐ **THE NO-PRE-SCREEN RESULT IS NOW WELL-SUPPORTED: 19 axes scanned, BH-corrected, ONE survives.**
     Dispersion (`arm_iqr` q=0.47, `abund_sd` q=0.75, `abund_iqr` q=0.75, spiker q=0.47) and EVERY
     cross-state shift (lfc HLY→TUM q=0.88, field leg q=0.92, `dGlobal_*` q=0.92, family-share shifts
     q=0.80–0.92) are null. ⚠ Two nominal hits die under multiplicity — lfc NAT→TUM (p=0.016, **q=0.12**)
     and `detection` (p=0.020, **q=0.12**). Only `own coupling ρ` survives (q=2.6e−21), and it is the
     outcome. ⇒ **stop looking for an a-priori screen; design (c) around out-of-fold instead.**

   * ⬜ **(e) THE GENE-RUNG QUESTION — NOT YET ASKED (user-queued 2026-08-19).** Everything above is the
     FAMILY-cell aggregate. The operational quantity is the gene's **total incoming pressure** summed over
     ALL its families, so the open test is what adding the excluded mates does to **per-gene aggregate
     anti-coupling** — the rung `gene_card.total_pressure` / `realized_rho_adj` live on. A family-cell wash
     does NOT imply a gene-level wash: the mates are unevenly distributed across a gene's families, so
     gene-level effects can concentrate in the wide-design genes even where the per-cell median is zero.
     ⚠ Use `apm-gene-question` first (aggregation function, prune, gating are all live choices there), and
     mind the degeneracy — 48% of genes have ONE family, where the gene rung IS the family rung.
   ⚠⚠ **CANDIDATES, NOT EDGES** — no decoy, no null, one cohort. **(→ MH-248, MH-166)**

7. ⬜ **Restore an exogenous existence validation** — see §A. Not cheap, but both `CN_INSTRUMENT.md` and MH-126
   name it **the single highest-value open item in the program**, and nothing else changes the program's
   standing as much.

---

## Y. Mined from the PRESSURE arc — analyses the LEARNED model can do and never has (2026-07-17)

`METHODS.md` specs the pressure arc. The estimator is retired (MH-115) — but **two of its analyses were never
ported to the learned model, and both target live gaps.** Found by asking what METHODS does that §6b doesn't,
rather than by merging the docs. (Its other content is not unique: **coupling** → `FORMULAS §7`, the
anti-correlation ladder; **data/edges/harmonization** → `DATA_SOURCES §0/§1/§2`.)

- ✅ **DONE 2026-07-18 (MH-164) — UNDETECTED, gate retired to rationale-only.** Ran the never-before-run test:
  per-gene `r ~ M + ago + prolif + M:ago + M:prolif`, ago ⊥ confounders, shuffled-ago null. The ago-specific
  `M:ago` interaction is **+0.00656, perm p=0.116, ~40% power (MDE 1.7× observed)** — a sub-significant trend
  in the predicted (positive) direction, NOT "AGO does not modify β". The naive Stouffer Z=+7.8 was
  pseudo-replication (shared patients). Over-control excluded (raw +0.00118 vs ⊥C +0.00656 — residualisation
  un-masks, ago's 40% prolif component carries the opposite modification). Per the rule below (*null retires
  the gate*): the AGO gate sits on biological rationale, not demonstrated coupling-modification, at n≈1000.
  `learned/ago_interaction.py`. rigor-auditor CONCERNS→admissible (P). **Details kept below for provenance.**
  ⬛ ~~**⭐ DOES RISC AVAILABILITY MODIFY β? — the learned model has NEVER used the AGO/RISC capacity gate.**~~
  **VERIFIED:** `learned/` contains only `ago_loading.py` (*which arm* is loaded, from Manakov chimeric —
  already measured **coupling-INERT**). **`ago_gate.py` DOES exist and ~25 pressure-arc modules import it** —
  it is absent specifically from `learned/`: **no per-sample RISC-capacity term anywhere in the learned tree.**
  *(Corrected 2026-07-17 — this item previously said "no `ago_gate` … anywhere", which is false repo-wide; the
  claim is true only of the learned lane.)* METHODS §4's gate is a different object: a **per-sample, miRNA-independent AVAILABILITY**
  term (`ago_load(s)` from AGO1–4 + a required TNRC6/GW182 effector, both-required, z-scored per gene) — the
  finite shared pool every miRNA competes for.
  **AND THE PRESSURE ARC NEVER TESTED IT EITHER.** It applies the gate as a **multiplier** and reports every
  result "gated AND ungated", calling it *"a documented sensitivity layer, **not a causal model**"*. The old
  `WHATS_NEXT §2.3` posed the real test — *fit `expression ~ pressure × AGO_capacity` and test the INTERACTION*
  — and it was never run. **The learned model makes it sharp for the first time:** β is the *de-confounded*
  coupling, so `r ~ (X·β) × ago_load(s)` asks whether RISC availability modifies the **real** slope, not the
  heuristic's. **A positive interaction upgrades AGO from "motivated" to "demonstrated"; a null retires the
  gate.** Either is worth having — the gate currently sits in the spine on biological rationale alone.
  ⚠ Pre-register: `ago_load(s)` is built from RNA of AGO1–4/TNRC6, so it correlates with global transcription
  and with proliferation — the interaction needs `mal_prolif` and composition in C, and a **shuffled-`ago_load`
  null** (axiom 4: state what the artifact would do to the test). **(→ METHODS §4; `ago_gate.py`)**

- ✅ **RESOLVED BY MEASUREMENT (MH-154) — AND THE ANSWER IS THE OPPOSITE OF WHAT THIS ITEM SAID.**
  This item used to read *"the pressure arc already built the null the learned model lacks … so the port is
  the **permutation** leg, not the BY leg."* **Both halves are wrong. Measured, 2026-07-17:**
  - **The permutation null IS the theoretical null.** Backing the scale out of `coupling_permutation`'s own
    `p_z` column (MH-45, **2000** Freedman–Lane draws): **σ₀ = 0.0309, μ₀ = +0.0001** (fit R²=0.9997) vs the
    theoretical **0.0311, 0** — agreement to **three significant figures**. **Porting it would have bought
    nothing.** Structural reason, and it generalises: **a permutation destroys the very structure that
    inflates the real null.** Freedman–Lane residualises on C ⇒ removes the *modelled* confound only; MH-123
    showed composition adjustment does **not** close the gap ⇒ the culprit is **UNMODELLED**, and no
    permutation can preserve an unmodelled confound. A real-data impossible-edge class can.
  - **BY did not protect MH-45 either.** Its dependence-robust count is **1,013/3,040 = 33.3%** (BY, ρ<0) —
    **inside** the **35.3%** rate at which *impossible* edges clear the same gate under the same C
    (CPE+HRD+batch = MH-123's core C). BY corrects **dependence among tests**; it cannot rescue a **mis-scaled**
    null. ⇒ **Do not read a BY q here as protection.** *(Both senses of "family" remain correctly implemented —
    the defect is the null they sit on, not the correction.)*
  - **The learned model's permutation port already existed and was never run:** `learned/coupling.py`
    (family-aware FL + `q_bh`/`q_by`/Simes `q_family`, per-gene C adaptation) is **imported by nothing**, with
    **no output, no ledger row, no registry row**. Given the above it is **not worth wiring** — keep it as a
    baseline, do not promote it.
  - ⇒ **THE CALIBRATION STANDARD IS THE SITE-FREE NULL, and it is now a module:** `eval/site_free_null.py`
    (rescued 2026-07-17 from a dead session's `/tmp` scratchpad — an MH-149 casualty in waiting; reproduces
    MH-123 to the digit). Efron form (**fit** location+scale per arm-abundance quintile, score N(μ₀,σ₀));
    exceedance counting is resolution-limited and can never fire BH. ⚠ It is **conservative** (the site-free
    class holds some real non-canonical targets) — the truth lies between it and the theoretical null.
  - ✅ **DONE — VERIFIED IN CODE 2026-08-03; this item was STALE.** ~~wire `site_free_null` into
    `discovery`/`dossier` as the gate~~. It is wired in **three** places: `learned/discovery.py` fits the null
    **IN-LOOP per gene** (`calibrate` imports `site_free_null._fit_bins`; its docstring explains *why* in-loop
    rather than calling `site_free_null.fit()` — that module residualises on **C only**), `learned/card.py`
    uses `SFN.fit_state()` for the per-state **calibrated** coupling that replaced the −0.1 cut (MH-166), and
    `learned/eval/dossier.py` replaced its `realized_coupling < −0.1` rule with a **sweep** on the calibrated
    quantity. ⚠ The `permute` path the item describes is still present but is **explicitly marked DEPRECATED
    in the source** (`discovery.py:497`), not the live gate. **Reserve "discovery" language for the AGGREGATE lane.**
  **(→ MH-123, MH-154; `eval/site_free_null.py`, `output/site_free_null/`)**

- ⬜ **DISCOVERY LANE — calibrated + expanded (MH-155), finding still DEFERRED.** `learned/discovery.py` now: fits
  the site-free null IN-LOOP per gene (canonical Gibbs he_agg, batch off), reports `q_by_arm` (BY) + `q_seedfamily`
  (Simes→BY, the de-duplicated hypothesis), attaches lower evidence (`chimeric_wt`, `ev_classes/ev_npmid` deduped
  via the ledger, `ts_mag`), flags `same_seed_he` / `he_max_corr` / `no_he_gene`. **Two new lanes:** LANE 2
  (`universe='hallmark'`) adds the ~2,792 Hallmark genes with no curated HE arm; LANE 1 (`scan_families`/`run_families`)
  tests whole seed families uncurated for the gene as a dose aggregate. **Honest state:** per-edge ~nothing survives
  BY (2 arm / 1 family on the HE universe) — the signal is the SET-LEVEL shift (median null_z −1.6 ⭐ **now
  −1.32 on the corrected family null, MH-197**) and the
  concordance of INDEPENDENT evidence with coupling (chimeric-present edges couple stronger, MWU p=4.6e−8).
  ✅ **FAMILY-NULL CORRECTNESS FIX LANDED (MH-197):** pseudo-families are now co-expression-matched per
  candidate cell; the old uniform draw was **20–25% too narrow**. Survivor count is unchanged at 0 (a wider
  null cannot create one) — what moved is the set-level readout above. ⬜ **NEXT:**
  ① ✅ **CLOSED 2026-08-02 — framing decided by the user ("just write both"): BOTH rungs written, MH-199
  (cluster, biological) + MH-200 (edge, candidate queue).** ⚠ Quote MH-199's assay/abundance-controlled OR
  (≈7–10), **never the 88% share** — that figure carries a HEK293T detection bias. — original item follows — The table is BUILT and now has a producer + provenance
  (`eval/discovery_gold_set.py` → `discovery_gold_set.tsv`, MH-198): **157 edges · 90 genes · 11 seed
  families**, each clearing bulk coupling + composition-robustness + a **site-coincident** Manakov duplex.
  **The framing choice is real because the set is 88% ONE co-transcribed cluster** (miR-17~92 + miR-106b~25 +
  miR-106a~363) — headline it as *the miR-17~92 cluster's realized target set in breast* (honest about the
  concentration, a biological claim) **or** as *a 157-edge convergent-evidence candidate table* (honest about
  the method, quieter about the biology). ⚠ **Do not compress this lane to "discovery is empty"** — that holds
  only for the per-edge FDR; ② ⛔ **DONE/NEGATIVE 2026-07-18 (MH-165) — the HE surfacing hope is at chance** (337 vs 287±35 shuffle, z=+1.4;
   83% family-lone), so **non-HE orphan discovery per subtype is NOT worth running** (weaker/noisier at n=93–197
   than HE, which already failed). ⬛ ~~**SUBTYPE-STRATIFIED run (user-asked)**~~ — PAM50-stratified
  discovery, since coupling washed out in the pooled cohort can surface within a subtype;
  ③ ✅ **DONE 2026-08-03 (MH-207) — family-lane evidence-attach BUILT, and it delivers a NEGATIVE that
  matters.** `discovery.attach_evidence_family` (wired into `run_families`) pools each family's member arms:
  `ev_npmid_fam` / `ev_classes_fam` / `ev_n_arms_supported` · `chimeric_wt_max` / `chimeric_n_arms` /
  `chimeric_src_fam` · `ts_mag_max` / `ts_n_arms` · `family_size_degenerate`. Purely additive (all 21
  pre-existing columns bit-identical). Coverage on 3,119 families: ledger **76.5%**, chimeric **45.3%**,
  context++ 98.9%, none-at-all 19.
  * ⭐ **PMIDs are UNION-DEDUPED, NEVER SUMMED — and it matters: summing inflates depth 2.55×**
    (33,012 → 12,935 over 1,590 multi-arm families; **20,077 double-counted PMID-slots**). Worst single row:
    NFAT5 / miR-17-5p/20-5p/93-5p/106-5p/519-3p, **sum 159 vs union 34 = 4.7×**. Same-seed members are
    routinely assayed in the SAME paper, so a summed depth is largely a family-SIZE proxy.
  * ⚠⚠ **EVERY pooled statistic is MONOTONE IN FAMILY SIZE** — chimeric support climbs **25.2% (1 arm) →
    89.1% (7+ arms)**, ledger **60.5% → 99.5%**. Condition on `n_family_arms` or the comparison is circular.
  * ⛔ **AND THE ARM LANE'S CHIMERIC CONVERGENCE DOES *NOT* TRANSFER TO THE FAMILY RUNG.** Naive pooled MWU
    looks real (**p=0.0074**) but is a size artifact and a trivial effect anyway (median null_z −1.329 vs
    −1.312, **Δ=−0.017**). Conditioned: **size-stratified Stouffer Z=−0.31, p=0.76** with **incoherent signs
    (4/7 strata negative; the two LARGEST strata point OPPOSITE ways)**, and the rank-partial on
    (family size, arm abundance) gives **ρ=−0.0101, p=0.571**. Two independent methods agree. ⇒ **MH-155's
    "chimeric-present edges couple stronger (p=4.6e-8)" is an ARM-rung result; do NOT broadcast it to families**
    (MH-111's confounder-architecture principle: the unit predicts the tool).
  * ⛔⛔ **AND THE LEDGER-DEPTH SIGNAL IS DEAD TOO — RIGOR-GATED SAME SESSION.** It briefly read *"what
    survives conditioning: partial ρ=−0.0578, p=0.00124"*. **The decisive control is WITHIN-GENE** (depth is
    (arm,gene)-specific, so a within-gene contrast removes gene-level study intensity, expression level,
    dynamic range and composition at once): on **2,566 families / 663 genes**, gene-demeaned and
    size+abundance-controlled, **ρ = +0.0032, p = 0.889 — zero, and nominally the wrong sign.** The pooled
    effect was **entirely BETWEEN-gene.** Fame null (depth shuffled within gene, 500 draws) centres at
    −0.0010 ± 0.0331 with the observation at **empirical p=0.55**; the study-bias channel is **larger here
    than MH-196 measured** (within-gene spearman(depth, abundance) **+0.330 mean / +0.500 median, 70.7% of
    604 genes** vs +0.187/+0.244/67.3%); and an **internal null control fires the WRONG way**
    (`ev_n_arms_supported` ρ=+0.0607, p=0.0081). ⚠ The degeneracy split is **incoherent, not a rescue**
    (multi-arm −0.0636 p=0.059 vs single-arm +0.0651 p=0.133 — opposite signs, neither significant).
  ⇒ ⭐ **NET: the family lane's deliverable is the ATTACH and its two design rules — NOT a finding. Neither
    evidence axis carries convergence at the family rung, while the ARM rung does (MH-155, p=4.6e−8) ⇒ run an
    evidence-vs-coupling test at the rung the evidence is RECORDED at** (generalises MH-111).
  **(→ MH-155, MH-207; `output/learned/{discoveries,discoveries_family}.tsv`)**

- 🚫 **Do NOT merge `METHODS.md` into `FORMULAS.md`.** ~87 KB of spec for an estimator MH-115 **retired**;
  consolidating it buys tidiness on a baseline. Both keep a retired-estimator banner and stay. The value in
  METHODS is the two analyses above, and they are now mined.

### Y.2 — PRESSURE-ARC → LEARNED PORT MAP (2026-07-18, full audit of the 12 pressure modules + pressure_dev)

The retired estimator is dead, but the arc's **auxiliary machinery and biological controls** are reusable and the
learned model lacks them. **Do NOT archive these files — extract from them.** Verified by reading, not docstrings.

**⬜ THREE NEW ANALYSES the learned model wants and has never done:**
1. ✅ **DONE 2026-07-18 (MH-165) — subtype coupling HETEROGENEITY is REAL (distributional), but the discovery
   hope FAILS.** Per-edge realized coupling is heterogeneous in strength AND sign across PAM50 (16% vs a
   marginal-preserving null 5%, 3×; net of composition-slope-modification; NOT a detectability artifact —
   best_sub doesn't track n, ~half sign-flip). ⚠ a **distributional** shift (tag P), not validated edge lists;
   no new betas (fixed pooled β re-scored). The "coupling surfaces in a subtype" hope is **at chance** (337 vs
   287±35, z=+1.4; 83% family-lone), and subtype β-remodeling (wiring) is weak (+12% over placebo). `learned/
   subtype.py::genomewide_edge_heterogeneity`. Method: use the within-subtype y-permutation null (a label-shuffle
   is anti-conservative, MH-154). **Details/pre-registration below.**
   ⬛ ~~**SUBTYPE-STRATIFIED LEARNED LANE**~~ — unifies the discovery §E ask (PAM50-stratified discovery) with the
   pressure arc's subtype contrasts (`subtype_contrasts` / `visibility_archetype_contrasts` / `pam50_gene_resolution`).
   `learned/subtype.py` does only the per-PAM50 **coupling test** — NOT "who is pressured, by which miRNAs, the
   immune axis, cold-Basal/hot-Luminal archetypes" on learned β. Build once, serves both. **HIGHEST value.**
2. **LEARNED-β PROGNOSTIC** — the outcome arc (`analyses/outcome/`, 6 modules) is pressure-based, and pressure
   *magnitude* was measured non-prognostic (`outcome-prognostic-arc-verdict`). Untested: does **learned β** predict
   outcome where the heuristic didn't? Source: `analyses/pressure_dev/pressure_prognostic_{signature,gene_centric,improve}`.
3. **EDGE-SOURCE SENSITIVITY** — does learned coupling change with the edge *source* (ENCORI vs miRTarBase)?
   Only ever asked of the pressure spine (`dual_spine_comparison`). The learned model fixes on pooled-HE/ledger.

**⬜ REUSABLE MACHINERY / CONTROLS to extract into the learned model (the "covered" files are NOT redundant):**
- **`mirna_state_class`**: `_mirna_trajectory` = per-miRNA **dosage/rank trajectory** GTEx→NAT→tumor (states.py has
  coupling-per-state, NOT the dosage trajectory) + `_tau` **subtype-specificity** + struct-vs-abundance audit.
- **`robustness_checks`**: `_proliferation_proxies` = **three** proxies differing in E2F/MYC dependence incl. the
  **E2F/MYC-INDEPENDENT over-control** (the reviewer answer that miR-17~92 hubs *are* E2F/MYC targets) + Basal
  scoping + partial-corr ladder. Control constructions `prolif_verdict` may not have.
- **`pressure_attribution_validation`**: **biological control sets** — positive controls (canonical breast
  regulators abundance-sum masks, MH-31/32) + `gene_corepression` + `per_target_arm_reranking`. A source of
  control genes/sets for the learned decoy/validation lane.
- **`healthy_anchor`**: GTEx↔TCGA **arm harmonization** + field-effect-aware healthy baseline (anchors GTEx, not NAT).
- **`coupling_predictor_comparison`**: abundance-sum baseline + edge-conditioning-mechanism (partly covered by `decoy_bench`).

**N/A — genuinely retired (the Gibbs has no such knobs; do not port):** `hybrid_pressure`, `pressure_engine`, and the
`pressure_dev/` weighting sweeps (`hybrid_param_sensitivity`, `pressure_constant_sensitivity`, `pressure_evidence_sensitivity`,
`pressure_sensitivity`, `pressure_layer_comparison`, `denominator_attribution_sweep`, `build_highconf_pressure`).

---

## Z. Absorbed from the parked design docs (archived 2026-07-17)

Five parked docs were archived; these are the only threads in them that were still alive.

- ⬜ **The DEFINITIVE orphan retest** (from `PRIMARY_TUMOR_ROADMAP`, and MH-55's own stated next step):
  *"re-run orphan coupling on the TCGA spine with composition metagenes"* — **n=1,041, not MH-55's n=40 DCIS**.
  MH-55 screened orphans in within-patient DCIS and called it *screening*; the fixed CPTAC screen (594→23,
  MH-114) and the re-run `buffa_validation` (108→3, MH-38) both now exist, so the powered version is the
  natural close. **(→ MH-55, MH-38, MH-114)**
- 🚫 **Pressure-formula options A–E** (`PRESSURE_FUTURE_OPTIONS`) — **DEAD BY SUPERSESSION, do not revive.**
  They tune the evidence-weighted pressure heuristic, which **MH-115 RETIRED**: on CPTAC protein the heuristic
  is p=0.11 and on CPTAC mRNA p=0.83, while abundance has the WRONG SIGN and the learned β gives −0.036
  (p=5.3e-13). *"The evidence weighting adds nothing over raw abundance, ever."* **Tuning a retired estimator
  is not an option, it is a category error.** **(→ MH-115)**
- 🚫 **`LEARNED_REGULATORY_MATRIX_DESIGN`** — the design note for learning M, self-declared *"parked, not in
  any spine"*. **The learned model is production** (§6b). Historical; archived.
- ⛔ **The CRC port** (`CRC_PORT_LITERATURE_SCAN`, 2026-06-28) — a literature landscape for porting the machine
  to colorectal cancer. **Blocked on the obvious: the breast machine's own existence claim rests on ONE
  observational line** (§A). Porting an unvalidated machine to a second tissue multiplies the claim, not the
  evidence. Revisit only after §A. **(→ §A)**
- ⛔ **METABRIC-full / the prognostic design** (`PRESSURE_PROGNOSTIC_DESIGN`) — see §D. Blocked on a DAR that
  **was never filed**, and MH-76's frozen-panel test already answered the generalization question (**no**).

---

## A. Edge existence — the one open foundation

**⛔ Both CN instruments are RETRACTED. Edge existence rests on ONE observational line** (curated edges
anti-correlate more than abundance-matched site-free ones). The "three independent validations" headline is
**FALSE** — the other two legs were both observational and never independent of each other.
**(→ STATE_OF_PLAY Axis 2; MH-124r/126, MH-133)**

- ⛔ **`pi_causal` — DEAD, do not rebuild.** `pi_causal = γ·b` with **b an observational OLS slope** ⇒ a
  product-of-coefficients mediation estimator, not an IV; it re-derived the anti-correlation it was meant to
  validate. Clean reduced form, arm-clustered: **p=0.115, n.s.** Retraction is in production code at
  `learned/instrument.py:1128–1148`. `MH124_*` §5/§5b are **VOID — do not cite** (§4b is untouched).
- ⛔ **Within-arm two-way-FE replacement — DEAD.** Its control class was **71% real binders**; against a
  genuinely site-free control **τ=−0.0007, p=0.84**, and it fails the site-type efficacy ladder
  (8mer/7mer-m8/7mer-A1 not monotone, p=0.26). **Refutations at n=216k–235k pairs, not power failures.** (MH-133)
- ⬜ **THE OPEN ITEM: a new exogenous handle.** The identification strategy must be an **asymmetry the confound
  is blind to** (axiom 4). Copy number cannot see a site type — that is what killed both CN designs; whatever
  replaces them must pass the same ladder.
- 🔬 **AI2 — second instrument: germline-eQTL / methylation-as-instrument.** Mapped, unbuilt. Over-identifies
  the same β and breaks reverse-causation differently than somatic CN (a stronger Sargan/Hansen jointly).
  Needs an ancestry-PC C block (axis BQ). **The leading candidate for the item above.** **(→ CHANNEL_FUSION §9 AI)**
- ⚠ **Bounding this axis regardless:** the per-edge null is **3–4× too narrow**; site-free pairs pass "FDR
  q<0.05" at **25–35%**; under an honest empirical null **0.0%** of curated HE edges survive per-edge. The
  signal is a **set-level distributional shift, not per-edge significance** (MH-123/124). **Any FDR count in any
  doc predating this rests on the uncalibrated null** — see §G for the two docs that still need banners.
- 🚫 **Do not carry forward** the mapped-but-unbuilt CN-instrument refinements (BL subclonal/CCF-weighted CN,
  AD3 allele-specific CN on imprinted loci). They sharpen a design whose exclusion restriction is dead; they are
  not a route back to identification.
- ⚠ **Note on the CN *channel* (distinct from the CN *instrument*):** the fusion channel was only ever
  "real but weak" (~0.7% of β's information), and **the shipped readouts do not include it at all**
  (`readouts.py` calls `_gibbs_posterior` without `channels=`). The soft-F-weight gate (`_ADMIT_SOFT_C=10`) and
  the coding/host pleiotropy gate are built and correct but currently **do not touch the production table**.
  Decide deliberately whether to wire them or retire them — do not leave this ambiguous. **(→ MH-91, MH-99/101)**

## B. Decoy / claim-1 — the surviving positive result

**MH-137 is the row of record and supersedes every earlier magnitude: gap −0.0119 (p=3.4e-05), deconv −0.0090,
retention 0.76, 1,322 genes.** ⚠ **The magnitude has SHRUNK at every control fix:**
**−0.045 → −0.0306 → −0.0147 → −0.0119.** The honest statement is a **BOUND, not a win**.

- ✅ **DONE 2026-08-01 (MH-168/169) — competence map rebuilt on the canonical decoy.** ⚠ **It was BLOCKED, not
  neglected:** `gene_atlas.py` had been **dead for 15 days** (the MH-38 `__file__`-hop-after-a-move bug, from
  commit `e5d5d84`; it ran the full compute then died on the write), and fixing it exposed a **self-destroying
  CNV cache** that was making `assemble_gene` fail for **every** gene. Both fixed (MH-168). Atlas reproduces
  MH-144 structurally (47.1% / 27.6% / +0.551). Strata: **MH-169** — width reproduces MH-147 (−0.0295 vs
  −0.0105) but **width and abundance are inseparable** (ρ=+0.780; both partials vanish). **(→ MH-168, MH-169)**
- ✅✅ **CLOSED BY MEASUREMENT 2026-08-03 (MH-206) — THE EVIDENCE HOLE CANNOT MOVE THE GAP, AND BOTH PRIOR
  STATEMENTS OF THIS ITEM WERE WRONG.** The board and `decoy_bench._site_maps`'s own docstring said *"POSTAR3
  miRNA-target is NOT in this repo ⇒ cannot be closed with local data; what remains is ACQUIRING POSTAR3"*.
  Wrong on the fact **and** on the mechanism:
  * ⭐ **POSTAR3 IS on disk** — `data/external/POSTAR/human (1).txt.gz`, **676 MB, downloaded 2026-06-30**,
    one directory from where the docstring looked. But it is the **RBP binding-site table**: **221 distinct
    RBPs, 2,360,006 AGO2 records**, TNRC6A/B/C ~3k, DICER1 8.5k — and **ZERO miRNA-named entries** ⇒ it
    yields no pair-level `(arm, gene)` call, only **miRNA-ANONYMOUS occupancy intervals**.
  * ⭐⭐ **AND AN AGO-PEAK ∩ SEED-SITE LAYER IS A NO-OP BY CONSTRUCTION** — the structural point neither
    doc had: `build_decoys` (`decoy_bench.py:280`) already requires `(a,g) ∉ ctx["sites"]`, and layers
    (1)+(2) exclude **every** arm with a strong site or a Poisson-significant 6mer in that gene. An arm whose
    seed sits under an AGO peak is **already ineligible**. Only an **arm-resolving source on SEEDLESS pairs**
    can add anything — which is precisely what Manakov chimeric is, and it moved **2.7%**.
  * **ENCORI: 5.3× under-ingested as a LABELLED layer, 97.0% redundant as a SET.** `data/external_cache/
    encori/miRNATarget/` (8,647 per-gene files, 6,431 non-empty) holds **20,404 distinct (arm,gene) pairs**
    vs the shipped `ENCORI_starBase` **3,861** — but **19,799/20,404 (97.0%)** are already in the union via
    TarBase v9 / miRTarBase / Manakov. Of the **605 new**, **578 are already in the full site map** via
    layers (1)/(2), and ⭐ **0 of the 4,937 assigned decoys carry one (0.000%)**.
  ⇒ **The arc's "decisive remaining test" is answered: it is not a live test.** ⚠ Read the pair of numbers
  together per axiom 5 — "5.3× under-ingested" is the fragile framing, "97% union-covered, 0 decoys affected"
  is the robust one. What would still be needed is POSTAR3's **separate miRNA-target/degradome module**
  (not what is on disk), and layer (4) bounds its plausible scale at a few percent.
  **(→ MH-206; `eval/decoy_bench._site_maps` docstring corrected in code)**
- ✅ **DONE 2026-08-01 (MH-169) — tested on the restored canonical decoy: A_COMPETENT −0.0335 (n=392) vs C_WEAK −0.0073 (n=69), MWU p=0.0618.** Same verdict (does not clear the bar) but far closer than MH-147's recorded p=0.293, on a 69-gene arm ⇒ **underpowered, not "adds nothing"**. ~~Its only support fails~~
  the arc's own both-fake-sets rule (q=0.006 FAKE1, q=0.20 FAKE2). ⇒ carry MH-130's **"27% domain" as a live
  hypothesis, not an established partition**. The **width axis (`n_fam ≥ 3`) is CONFIRMED** and sharpened
  (3–4 families −0.0314, p=0.0006) — ⚠ but **not monotone under MH-137: it peaks at 3–4 and falls at 5+**.
- ✅ **Anchors — do not re-litigate:** MH-136's **seedless positive control** (187 all-seedless genes, gap
  **+0.0006 = exactly zero**; 810 all-seeded **−0.0325**) is the arc's first real positive control and stands ·
  the **1-family internal null is back** (−0.0007, p=0.15 on the clean control; MH-135's retraction of it was an
  artifact of the broken decoy) · **an abundance baseline is NOT a control — benchmark any aggregator against a
  fitted matched decoy** (MH-115/127) · MH-130b's whole **SNR-strata arc is dead** (double-log bug);
  ✅ `mh127_snr_strata.tsv` **ARCHIVED OUT OF THE LIVE TREE 2026-08-03** →
  `output/learned/ARCHIVE_POISONED_MH130b/` with the bug and the corrected re-run in its README.
  Zero code consumers at archive time (verified: only doc references).
- 🚫 **DROPPED 2026-08-03 — Aggregate "raw force" vs abundance-sum.** Design in
  `method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md` (design-only; **no module was ever
  built** — verified, `grep aggregate_force|raw_force` returns no code). It fails **both** of its own gates
  and is now superseded twice over: its `w_eff` **IS the evidence-weighted pressure heuristic MH-115 RETIRED**
  (*"adds nothing over raw abundance, ever"*), and its actual question — the doc's own words, *"regulators
  weighted by strength/promiscuity vs weighted equally"* — **has been answered in a strictly better frame by
  MH-201/204**: learned β vs unweighted abundance-sum (Δ −0.0072, p=6.2e-13 on Buffa) and vs a **fitted
  matched decoy** on three layers. Tuning a retired estimator against a non-control comparator is a category
  error (cf. the §Z ruling on PRESSURE_FUTURE_OPTIONS).
  ✅ **PURSUED AND CLOSED 2026-08-03 (MH-208) — and my note here was itself a CONTAINER ERROR.** It read
  *"no promiscuity / target-count axis exists in `learned/gene_axes.py` (verified)"*. True of that file, and
  **misleading**: a promiscuity annotation **had existed since 2026-06-28**
  (`analyses/misc/genomewide_promiscuity.py`, cached, one consumer — the retired-arc
  `coupling_predictor_comparison`). I checked the CONTAINER I expected it in, not the QUANTITY (axiom 7 —
  the same error I had just written up in MH-206). ⭐ The user caught the structural reason: promiscuity is an
  **ARM-level property and there is no ARM CARD** (`card_rungs.py` defines edge/family/gene only), so it had
  nowhere to live and never entered the learned model's axis system.
  * ⛔ **The annotation was a FAME axis:** curated `he_degree` ↔ #distinct-PMIDs **ρ=+0.736**, ↔ abundance
    +0.556, but ↔ a curation-free sequence targetome only **+0.124** — top-10 lists **disjoint**, medians
    **2 vs 3,634 targets**. The retired heuristic's `w/D(m)` was dividing by fame.
  * ✅ **FIXED + WIRED:** `build_sequence_promiscuity` (scanMiR RBNS K_D, 746 arms, cached, producer'd);
    `load_promiscuity` defaults to it with `fill="nan"` (missing = UNSCANNED); `reg_promisc_*` added to
    `gene_axes.regulator_axes(dose, arms=)`.
  * ⛔ **But the axis does NOT earn its keep:** marginal `reg_promisc_max` −0.119 (q=0.013) / `sd` −0.097
    (q=0.035) on MH-201's β-over-abundance margin **vanish** under abundance-concentration control
    (**−0.037 / −0.039, n.s.**), and the curated axis that fired was `reg_dose_hhi` in disguise (ρ=+0.381
    with it; +0.045 n.s. once controlled). ⭐⭐ **The transferable lesson: arm-level abundance-orthogonality
    (ρ=−0.004) DOES NOT SURVIVE AGGREGATION to the regulator ensemble (−0.167).**
  ⇒ **Do not revive the aggregate-force design. The promiscuity infrastructure is built and correct; the
  axis is a null.** **(→ MH-208)**
- ✅✅ **AND THE ROOT CAUSE IS FIXED — THE ARM CARD EXISTS (MH-209, user-diagnosed).** The user's read of
  MH-208 — *"we don't have a miRNA card"* — was the correct structural diagnosis: `card_rungs` knew only
  `edge`/`family`/`gene`, so per-arm facts had **no rung to live on** and six modules each re-derived their
  own. `learned/arm_card.py` → `output/learned/arm_card.tsv`, registered as the 4th card.
  ⭐⭐ **v2 SHIPPED 2026-08-03 (MH-211): 2,450 arms × 227 cols across 25 blocks** (`--check` CLEAN on all
  four cards; registry 328 → **559** rows). ⛔ **v1's "3,241 × 47" is RETRACTED** — its locus block was
  keyed by `MI*` locus id, so 506 non-arm rows inflated the denominator and `clustered`/`mirgenedb` were
  100% NaN. Corrected universe **2,450**, of which **2,236 = the entire TCGA matrix** (no expressed arm lost).
  * ⭐ **Provenance-prefixed by design**, because the targetome universes are NOT interchangeable:
    `seq_` genome-wide (746) · `site_` ⚠ **Hallmark-scoped, 1,432 genes** (771) · `ts_` TargetScan
    per-site (⚠ 321) · ⭐ `sdsz_` a **FOURTH** universe (MANE 3′UTR seed scan, 2,370) · `cur_`/`fame_`
    ⛔ fame, for CONTROL only · `cov_*` — **missing = UNSCANNED, never zero**.
    Taxonomy now has **three facets**: provenance · `bc_` (broadcast DOWN from a coarser unit — the
    opposite relation to `AGG_OF`) · `hx_` (heuristic / retired-estimator provenance).
  * ⚠ **The definitional control was REFUTED IN THE LETTER and is reported, not buried:** four of five
    MH-208 correlations reproduce to the digit, but curated↔fame moved **+0.736 → +0.752** because v1's
    universe held 76 non-arm bare stems (which alone carry only +0.370) ⇒ **MH-208 is STRENGTHENED**.
  * ✅ Other controls: no-impute 0 leaks · key-join cross-check `foc` vs `cnvc` **ρ=+0.982** ·
    broadcast constancy exact. ⛔ Excluded by audit: the whole §6b-retired `mirna_state_class` table,
    `comovement.dHT` (MH-210's fabricated quantity), `arm_survival.gse_*` (bare-stem matched),
    `gtex floor0`, and `ago_loading.loading()`'s `.fillna(1.0)`.
  * ⚠ Two of my own design assumptions died on contact: **sites-vs-genes is not an independent axis**
    (ρ=+0.927) and the **8mer share is TOOL-DEPENDENT** (ρ=+0.620, medians 2.2× apart) — never compare it
    across sources.
  ⬜ **Still to fold onto the card** (all arm-level, all already built elsewhere): `arm_rung` β +
  identifiability · `rarity_bench` seed rarity · `kd` per-arm κ0. **(→ MH-209)**
- ✅ **DONE 2026-08-04 — the card system is now READABLE without code.** `analyses/ops/gen_cards.py --build`
  → `docs/derived/cards/` : six standalone HTML pages (search any arm/gene/edge/family, or a column name to
  find which card carries it). Four audits ship with it and stand alone — `--check` (registry↔card, two
  independent signals), `--roundtrip` (lossless encoding, NaN mask exact), `--links` (cross-card referential
  integrity), `--verify-js` (runs the pages' OWN javascript under node vs the source TSVs). Details in
  `ANALYSES_CATALOG.md` + the ledger; **no numbers restated here.**
  ⬜ Follow-ups: rename `CARDS["family"]` → `gene_family` and its file (13 modules / ~30 refs) once no
  session is live in `realization.py`; re-run `--build` after any card rebuild (the footer stamps every
  input's sha, and `--check` fails the build if a card drifts from the registry).
- 🔨 **NEW, SURFACED BY MH-208 — give MH-201's axis scan a PRODUCER.**
  `output/learned/ood_cohort/{ood_cohort_regulator_features,ood_cohort_modifier_scan}.tsv` carry MH-201's
  headline `reg_dose_hhi` (q=2.1e-05) and its entire 82-modifier FDR scan, and **nothing in the repo writes
  them** (verified: `ood_cohort.py` writes `ood_cohort_genes.tsv`, `edge_leg_*` and the manifest). This is the
  **MH-196 shape** that rotted the literature ground-truth sets. The outcome column is reconstructible from
  produced data (`delta = rho_buffa_metagene − rho_abund_metagene`, verified to max|diff| 1.0e-16), so this is
  a write-the-producer job, not a re-derivation. **(→ MH-208, MH-196, MH-201)**

## C. Attribution / Shapley

**STANDS: the model RANKS the canonical family high; it does NOT NAME it.** ⭐ **Cite MH-196** (n=329 / 92
families, versioned + stamped) — **not** MH-126 Test 1 (n=32) or MH-124 §4b (n=21), whose sets are
**unreproducible** (no producer in the repo; MH-126c's is not on disk at all). β **0.436, p=0.021 — no longer
at chance**; `identity` 0.416 (p=0.004); abundance 0.359. Under the **family-fame null** only β's
gene-specific component is significant (Δ=−0.085, p=0.038) while abundance's is not (Δ=−0.052, p=0.40).
**Argmax is still at chance at every k.** **(→ STATE_OF_PLAY Axis 3)**

- ✅ **DONE + CLOSED (MH-196, 2026-08-02) — ~~Scale the literature ground-truth set~~.** The versioned pull
  exists (`eval/lit_ground_truth.py`, **329 genes / 92 families**, sha256-stamped) and replaces five
  producer-less lists plus MH-126c's n=32, **which was never on disk**. Two results and a closure:
  **(1) β is NO LONGER at chance** — 0.436 [0.382, 0.491], p=0.021; MH-126c's 0.518 is outside the CI.
  **(2) Abundance's win is largely PUBLICATION BIAS** — under a family-fame null its gene-specific component
  is Δ=−0.052, **p=0.40 (n.s.)**, while β's is Δ=−0.085, **p=0.038**. The study-bias channel is +0.187 within
  gene (p=3.2e-10).
  **(3) ⭐ THE ITEM IS CLOSED BY ARITHMETIC, NOT BY DOING IT AGAIN.** The head-to-head is TIED
  (−0.033 [−0.165, +0.085]) and **unresolvable in principle**: the cluster unit is the canonical family,
  only **~330 exist in the whole curated literature**, and exhausting them leaves half-width 0.065 > 0.033
  — resolving it needs ~1,241 clusters, **3.8× every canonical family ever published.**
  ⇒ **the ground truth is NOT the binding constraint and scaling genes cannot fix this.** Do not re-open
  this as a scaling task. Revised standing sentence: *"β carries a MEASURED, gene-specific gradient toward the
  canonical family; whether it beats abundance is beyond what curated literature can ever resolve."*
- ⬜ **One canonical rerun of the persisted attribution outputs** (`output/learned/programs/*`, cards) now the
  estimator is settled. Optional: canonicalize the `shift_vs_weight` per-edge weight diagnostic.
- ✅ **DONE (DO-FIRST #4, MH-140)** — `share`→`beta_frac`, true `identity` column shipped.
- ⚠ **Do not build a "do canonical regulators score higher?" test on `pip_discovery` or `prior_pi`** — both are
  **w-contaminated by construction**; the test is circular. (β itself passes the w-circularity gate
  bit-identically: max|Δβ| = 0.0 under shuffled *and* constant w.)
- 🚫 **Do not carry "dose weighting is dead" from MH-130e into the attribution arc** — it refuted the extension
  of §4b to the **coupling** estimand only. §4b's attribution claim is a different estimand and is not overturned.

## D. Cross-cohort & external validation

- ⬜ **`β(TCGA) → Buffa mRNA`** — see DO-FIRST #5. The program's only untested boundary.
- ⚠ **"METABRIC" in MH-73/74/75/76 means Buffa, n=207** — verified in code. The +0.056 / +0.060 / +0.028 /
  +0.019 are **four views of ONE 207-patient result, not four validations.** Do not recount them.
- ⬜ **METABRIC-full: file the DAR or stop citing it.** "EGA-pending" appears in **five docs** with **no
  accession, no DAR id, no submission date, no owner** anywhere in the repo. The evidence supports *"never
  applied for"*, not *"pending"*. ⚠ It would also be a **power** upgrade, not an **independence** one — same
  platform, same consortium; Buffa is a subset of it.
- 🔨 **CPTAC composition re-runs (MH-107) — PARTLY DONE, and the item needed re-scoping (2026-08-01, MH-171).**
  ✅ **`ood_protein` DONE:** its code fix had landed but had **never been run or recorded**. Prospective
  (independent, n=101) **median retention 0.797, 7/7 still protein-coupled**; TCGA-105 (circular) 0.656, 7/7.
  **ZEB1 worst in both (0.246/0.308 ⇒ `composition_explained`)** — confirms the compartment-arithmetic claim
  per-gene, from a different design than MH-114's shuffled null. **(→ MH-171)**
  ⬜ **`dossier.tier3_protein` — code already runs BOTH blocks and emits `coupling_retention` /
  `composition_class`; still needs a RUN + a record.**
  ✅ **`MH-83` / `MH-84` — RESOLVED 2026-08-01 (MH-182): a REGISTRY-FILING gap, not missing science.**
  The ids were never written, but the work is fully recorded — `LEARNED_MODEL_DISCOVERY_SYNTHESIS.md`
  declares itself the capstone of **MH-82→MH-84** and carries their provenance; numbers in
  `LEARNED_MODEL_VALIDATION.md` §1, decisions in `LEARNED_MODEL_RATIONALE.md` §2a–§2g.
  ⇒ **nothing to re-run and no rows to back-fill** (reconstructing them from a summary doc would be
  arithmetic on someone else's numbers — axiom 1a). The one number that WAS docstring-only is homed
  in **MH-172**.

  🚫 **DO NOT re-run `eval/cptac_validation.py` as a composition fix** — it imports `compute_gene_pressure`,
  the estimator **MH-115 RETIRED**; a re-run would mint citable-looking numbers from a dead estimator, which
  is exactly why STATE_OF_PLAY forbids reprinting MH-114's counts. It needs a **PORT to the learned
  posterior**, not a re-run. That port is the real open work here (MH-34/35/36/39/40).
  ⚠ MH-114 still supersedes MH-107's magnitudes (~57% compartment arithmetic; compartment-blind effect
  ≈−0.011, surviving adjustment at ≈−0.0065, p=2.0e-03).
- ⬜ **The CPTAC protein transfer is NOT decoy-controlled** — a site-free fitted fake reaches **99%** of it
  (MH-130d). It is real arithmetic but **not evidence the curated edge set is right**. *"The protein channel has
  never faced a decoy that could lose"* — the arc's own words. The only decoy-controlled transfer is the
  **CPTAC/mRNA** cell inside the competent class.
- ⬜ **V1 δ-transportability, remaining scope.** The *abundance* input transports (same member dominates 84.5%
  vs 43.7% chance; member-share correlation 0.991, MH-104). The **full CN/chimeric-fused δ** is untested and
  low-powered — the CPTAC CN instrument is F>10 for only **59/685 arms (8.6%)** vs TCGA's 32%; **49 arms / 5
  families** usable in both.
- 🔬 **Manakov chimeric-duplex validation of the discovery shortlist** — also the input to §B's evidence-hole fix.
- ⚠ **Standing rule for any future cross-cohort channel:** run `gauge.calibrate()` → read `info_ratio` + τ
  **BEFORE** building it. This is what killed the state channel; it is what falsifiable infrastructure is for.

## E. Discovery lane

- ✅ **DONE 2026-07-18 — the per-EDGE card is built as the CANONICAL two-block card** (`learned/canonical_card.py`
  → `output/learned/edge_card_base.tsv`, 5,650 (gene,arm) rows / 1,421 genes / 3,941 both-block). One `(gene,arm)`
  key, two provenance-tagged estimator blocks joined: **ATTRIBUTION** (Gibbs β/`identity`/`beta_frac`/`pip`, from
  `readouts_arm_edges.tsv`) + **PROGRESSION** (regime `arm_med_rpm`/`arm_pct_floor`/`spiker`; budget
  `states.budget_shift` rank/share + GTEx→NAT→tumour Δ; composition `retention_rho`; 3-state coupling
  `coupling_{hly,nat,tum}`; **`shift_class`**). ⚠ **`shift_class` is LEARNED-NATIVE** (`card._shift_class`,
  partial-ρ), **NOT** `mirna_state_class.joint_edge_class` (retired-pressure) — same names, different label
  system; never map one onto the other. Companion `canonical_card_provenance.tsv` maps every column→estimator.
  Also consolidated the card DUPLICATION (deleted 3 orphans, repointed `discovery.py`). Ledger 2026-07-18.
  ✅ **DONE — VERIFIED 2026-08-03; this item was STALE, and a RENAME is why.** ~~REMAINING: the per-GENE
  roll-up (G-series…)~~ — `output/learned/realization/gene_card.tsv` is built: **1,549 genes × 106 columns**
  (2026-08-02), and **all five named G-series items are present and populated**: net repression
  (`gene_repression_class`, `gene_net_repressed_tumor`, 1,409 non-null) · coherence+role (`realization_owner`,
  `static_owner_family`, `owner_agrees`, `regulatory_handoff`) · budget concentration (`concentration`,
  `top_beta_frac`, 1,549) · composition fraction (`median_retention`, `n_composition_explained`,
  `n_cell_intrinsic`) · shift-class (`dominant_edge_shift_class`, 1,260). Registered in
  `output/learned/card_registry.tsv` (edge 161 / gene 106 / family 61 columns).
  ⚠ **Why the board missed it — axiom 6, a naming collision:** the artifacts were **RENAMED on 2026-08-01
  (MH-181)** — `progression_gene_card` → `gene_card`, `progression_edge_card` → `edge_card`,
  `canonical_card` → `edge_card_base` (`learned/card_context.py:72`). This item was written against the old
  name and became invisible. **(→ MH-181; `card_context.py`, `card_rungs.py`)**
- 🔬 **Top-discovery deep-dives.** (a) the **miR-17~92 / miR-106b-25 target cluster** (miR-106b/93/17/19a →
  RABEP1, AHNAK, TGFBRAP1, IL6ST, AFF1, WWP1, M6PR — ubiquitous arms, cell-intrinsic, strong, weakly-curated):
  lit scan each, then CPTAC protein + Manakov binding + within-patient realization; is it a coordinated
  cluster→program hit? ⚠ note miR-93-5p is also the systematic **proliferation-confound arm** (MH-100) — check
  the `prolif_verdict` flag per edge before reading these as coupling. (b) **miR-135b-5p → GATA3** (fully novel,
  a **spiker**, 28%>floor ⇒ coupling lives in the miR-135b-high subset): which tumours, is the subset a clinical
  stratum? (c) the **568 fully-novel candidates** → triage (cell-intrinsic + expressed/spiker + biochemically
  supported), lit-novelty check, shortlist for protein + Manakov.
  ⚠ Every FDR-based candidate count here inherits the **uncalibrated per-edge null** (§A).
- ⬜ **Epigenetic-rewiring discovery — the bidirectional methylation lane.** Detector built
  (`structural_identity` + `methylation`; positive controls pass: miR-124/9/129/137 hyper→LOST, miR-21/155
  hypo→GAINED); needs genome-wide scale-up. **Axis A** = LOST specialist (baseline-active, tumour-silenced,
  promoter HYPER-methylated Δβ≥+0.15 → target de-repressed). **Axis B** = GAINED specialist (baseline-silent,
  tumour-active, HYPO-methylated Δβ≤−0.15 → target acquires repression). Nominate on `structural_identity`
  (`structural_share ≥ 0.10` ∧ `confidence ≥ 0.5` — expression-free, so it survives silencing); the ONE
  inferential claim is the **target consequence** (group Δ + within-tumour g-mRNA vs promoter β; permutation
  null, BY-FDR). **Must condition:** purity/composition (global tumour hypomethylation is the known Axis-B trap
  ⇒ demand CGI-promoter-local, not open-sea) · CN-vs-methylation (cross `mirna_locus_cnv`) · imprinted loci
  (DLK1-DIO3 — flag & exclude; miR-127 baseline β 0.92 is exactly this). **Prereq: the UUID→TCGA-barcode map.**
  Power-bounded (normal n=97); the tumour arm is deep.
- ⬜ **Target-side mechanism legs — the reference-blind quadrant.** `identity ≫ structural` is not only
  "opportunistic": it can be a coupling via a site the **reference genome does not encode**. Build the
  **seed-SNV leg** (APM has SNV/VEP): scan each gene's 3′UTR tumour SNVs for seed-match **gain** (acquired
  site) or **loss** (the target-side analog of the silenced specialist). Same family: **A-to-I editing** and
  **APA**. ⇒ the identity/structural quadrant is a **mechanism router** — loss side {methylation, seed-SNV-loss,
  CN-loss}; gain side {hypomethylation, seed-SNV-gain, editing/APA-exposed, uncurated-real, residual-confound}.
  Subsumes the old "cross EMT/NF-κB genes with APM SV/3′UTR disruption" item.
- 🔬 **context++ × K_D fusion.** The two are **nearly orthogonal — Spearman +0.065** over 78k (family,gene)
  pairs (K_D = scanMiR RBNS *binding*; context++ = TargetScan *repression-outcome* + conservation + site
  context) ⇒ exactly where a fusion can beat either. MH-86's "don't combine" was HE-conditioned and used the
  worse (HE) K_D. Test for specificity + discovery candidacy. **(→ MH-87)**
- ⬜ **K_D scan scope + discovery floor sweep.** Current genome-wide K_D = **746 detectable arms** (RPM≥10).
  Optional: (1) a complete **~2,600-arm reference scan** (one `scanmir_genomewide_par.R` run) as a reusable
  `(arm,gene)` resource; (2) a floor sensitivity sweep (RPM ≥ 50/10/1) — a re-*filter*, not a re-scan; the
  floor is the sensitivity/precision lever, not an FDR lever.
- ⬜ **CLIP-source candidate expansion** (ENCORI / POSTAR3 physical sites) beyond TargetScan nomination.
- ⬜ **Regulation-pattern archetypes + seed-rarity.** (a) Tabulate the per-gene pairwise-ranking concordance
  matrix over the four rankings (abundance · budget · identity · structural) → cluster genes into archetypes
  (abundance-controlled / weight-controlled / designed-but-silenced / redundant-collinear) — the systematic
  version of the ad-hoc PTEN (structural anti-abundance ρ=−0.65) / ESR1 reads. (b) Ground **structural
  specificity in targetome RARITY**: how many genes carry the arm's 7–8mer genome-wide, tied to scanMiR K_D +
  context++ (rare strong site = genuine specialist; ubiquitous seed = pseudo-specialist) — makes "designed to
  *specifically* repress g" a sequence claim. (c) **GLOBAL** rank concordance via `program_network` engagement
  (Σ_g M) — loud-but-unwired vs scarce-but-hub arms at genome scale.
- ⬜ **State-stratified discovery** — run the discovery lane **per state** (own permutation-FDR null + scanMiR +
  deconv tag), then compare: found-in-tumour-not-NAT = candidate **acquired** regulator · found-in-all =
  constitutive · found-only-NAT = field/lost. **Distinct from the wiring axis** (which uses a *fixed* support to
  *compare* known edges). Power-bounded: NAT n=104 / GTEx n=327 ≪ tumour ~1000; and *undetectability* is worst
  in healthy ⇒ "not found in healthy" ≠ "not a healthy edge". NAT is the more feasible end.

## F. Progression / state — the landscape, not the channel

**Two different objects. Do not conflate them.** **(→ STATE_OF_PLAY Axis 6)**

- ⛔ **The state CHANNEL is CLOSED — measured and cancelled. DO NOT REBUILD.** `learned/channel_state.py` was
  never built (verified absent). τ — the acquired-vs-constitutive payload, the axis's **entire deliverable** —
  is indistinguishable from zero (GTEx τ=0.0009, info **0.6%**; NAT 0.7%; CPTAC 0.7%). Structural: channel
  precision ∝ **a²**, and a≈0.11 ⇒ **even unlimited GTEx donors add ~1%**. The additive `M_T = M_H + Δ` form
  **would have faked the headline** (Δ>0 on nearly every edge, manufactured out of composition attenuation).
  **(→ MH-102d, `LEARNED_MODEL_STATE_CHANNEL_PLAN.md §10`)**
- ✅ **The abundance-level cross-state LANDSCAPE is built and run** (GTEx-healthy → TCGA-NAT → tumour, 721 arms
  / 5,108 edges; acquired_realized 640 · lost 209 · stable 2,141 · acquired_unrealized 1,775).
- ✅ **DONE (verified 2026-08-01) — `docs/reports/LANDSCAPE_REPORT.md` carries the banner.** ~~needs a banner — its per-edge FDR class labels rest on the null MH-124 measured~~
  as **3–4× too narrow**, and it carries no warning. ⚠ NAT is structurally underpowered (n≈104):
  `nat_decoupled` = **1 edge**.
- ⬜ **Abundance-vs-wiring decomposition** (Δx·M vs x·ΔM) as a posterior readout + **subtype interaction tests
  (Gelman–Stern)** — a readout over the landscape, not a channel.
- ⬜ **GTEx logFC gauge check** — `parallel_view.change_trajectory` reports both the naive RPM−TPM (`rawFC_GN`)
  and the QN-into-TCGA (`relFC_GN`) GTEx→NAT change; validate the miRNA RPM≈TPM assumption against the QN
  version to find where they diverge (housekeeping/spike arms).

## G. Correctness debt, docs, figures, infra

- ✅ **DONE (verified 2026-08-01) — the simplex-illusion correction is in `METHODS §1` and its twin.** ~~Fix the rationale. It justifies~~
  the C block by *"Cancer-Epithelial is deliberately EXCLUDED — conditioning on the compartment the target is
  expressed in over-controls the signal"*. **That is VOID:** the 9 Wu-major fractions sum to exactly 1.000000
  and **R²(Cancer-Epithelial ~ the 8 held-out) = 1.0000 in both cohorts** ⇒ the hold-out is a **simplex
  illusion**. The practice is nevertheless **safe for a different, measured reason**: purity contributes only
  **0.9%** of the retention loss while the stromal MIX contributes **33.1%**. **(→ MH-111/114)**
- ⬜ **Regenerate presentation F25 + re-render the deck.** `talk.qmd` prose was fixed 2026-07-12 and is honest;
  **the F25 figure beside it is from 2026-06-24** and shows the old numbers.
- ⬜ **Posterior calibration — ✅ DOC CONFLICT RESOLVED 2026-07-17; the fix itself is still open.**
  **SETTLED VALUE: widths are ~25% too narrow** — bagged NNLS **0.73×** · Gibbs Gaussian **0.77×** ·
  **Student-t ν=7 0.89×** (MH-102e, `learned/calibration.py`). The two docs never actually disagreed on the
  measurement — this board reports **reliability** factors and `CPTAC_PROTEIN_CHANNEL_PLAN.md` reports their
  **reciprocals** as **inflation** factors (1.37× / 1.29× / 1.13×). The lone outlier was
  `LEARNED_MODEL_STATE_CHANNEL_PLAN.md` §11(a)'s **"~40% (NNLS 28%, Gibbs 39%) / SBC promoted"**, measured on a
  **biased gene subset** (the scale-dependent sd-floor bug) — now **bannered SUPERSEDED at source**.
  **SBC stays RETIRED as the wrong tool** (the fault is the MODEL — unmodeled between-participant heterogeneity —
  not the sampler, so SBC would PASS and mislead).
  **To apply:** multiply reported widths by the inflation factor; **turning `nu=7` on is a live recommendation**
  (⚠ `readouts.py` does not pass `nu` — the shipped table is Gaussian, i.e. the 0.77× rung).
  **Triage:** inverse-variance **fusion weights are SAFE** (SE *ratios* cancel); **all OOF/permutation results
  are UNTOUCHED**.
  ✅ **MH-94 WIDTH SUB-CONFLICT — DIRECTION SETTLED 2026-07-17 (MH-146); the flagship is NOT.** MEASURED
  (`calibration.shapley_resampled_width`, 64 family-rows): the posterior width **UNDERSTATES by ~1.4×** on
  families that carry identity (median 1.37 at identity>0.05) ⇒ **the board's DIRECTION was right; the CPTAC
  plan's "1.8× too WIDE" is RETIRED.** ⛔ But the tool's own printed headline ("understates by 112.86×") is a
  **zero-denominator artifact** — 48% of rows have `posterior_sd ≈ 0` (zero-identity families) where the ratio
  explodes to 3881–7655×. ⬜ **Open:** gate that summary on `identity_mean > 0.05`, and re-run MH-94's
  flagship (PTEN miR-141/200a 0.77±0.41 ⛔ **[RETRACTED — MH-150, 2026-07-17: re-run at MH-94's OWN config gives **0.003 ± 0.012, PIP 0.042**. The SHARE collapsed 250×, not just its width; the exemplar is void. PTEN's identity is now miR-17~92 (0.501) + miR-21/590 (0.223) — the canonical regulators. Prime suspect: the `_rtnorm_pos` sampler bug, which MH-94 predates. DO NOT CITE.]** vs MH-102e's re-derived 0.243) at ITS OWN config — the configs here
  are not comparable, and 0.41→0.243 points the OTHER way. **(→ MH-146)**
- ✅ **CLOSED (MH-151/152) — the isomiR cache needs NOTHING.** Coverage of the model cohort is **99.9%**
  (1,078/1,079 participants; ZERO missing). The manifest's 18 extra Primary Tumor files are **duplicate
  vials**, not missing people — a rebuild and a download both gain nothing.
- ✅ **ANSWERED (MH-153) — YES, it survives, and the risk is elsewhere.** miR-17~92's dose is **corr 0.95**
  under `isomir=True` and **PTEN's identity argmax is UNCHANGED** (miR-18-5p; identity corr 0.92). **Donor
  families are BUFFERED BY POOLING** — losing ⅔ of one member's reads is absorbed by the rest.
- ⬜ **THE REAL EXPOSURE, newly surfaced (MH-153) — audit RECEIVERS and the ORPHAN MASS, not donors:**
  (a) small families have no pool to dilute a donation — `miR-543` gains **+351%** from jumped reads;
  (b) the `orphan` bucket **DELETES 4.1% of total RPM** (654.3M → 627.4M) — reads whose shifted seed
  matches no family are dropped from EVERY family. **`X_fam` does not conserve mass and nothing flags it.**
- 🔨 **BEFORE any OOF isomiR run (MH-153): fit the shrinkage target on TRAIN ONLY.** `seed_composition`
  Dirichlet-shrinks low-read cells toward the **cohort** family-distribution — a cross-sample X-leak (no
  Y-leak). Confined to low-read cells and the path is default-OFF, so nothing shipped is affected.
  **(→ MH-153, MH-152)**
- ⬜ **Infra:** ✅ the `learned/` tree IS tracked (verified 2026-08-01; the 2026-08-01 commit added the rest) · fold
  ENCORI/POSTAR3/Manakov into the ledger (**union, not summed**) · `baselines/` re-export shims · vectorize the
  per-sample segment×locus overlap in `mirna_locus_cnv` (interval join per chromosome — the main slow step in
  `run_all`) · covariate-protected ComBat batch **only if** a channel reintroduces batch (never naive dummies).
- ⬜ **Doc consolidation, continued.** This board absorbed `WHATS_NEXT.md` + `LEARNED_MODEL_WHATS_NEXT.md`.
  Still queued: merge `LEARNED_MODEL_METHODS` → `FORMULAS.md`; merge `MODULE_MAP` → `ANALYSES_CATALOG`; archive
  `LEARNED_MODEL_{DESIGN_RESPONSE,DISCUSSION_PROMPT,BUILD_PLAN}`. See `STATE_OF_PLAY.md`'s **doc-traps table**
  for which docs read as current but are not.
- ⬜ **Low-priority coordinate/coverage debt** (deferred, cheap, no claim depends on it): extend MirGeneDB coords
  to the high-abundance expr-only MIMAT tail (~1,418 GDC MIMATs lack loci; ~24 have mean log2RPM>2; lead gap
  **miR-625-3p**) · the MI*→GDC hairpin bridge (348/506 mapped, 158 zero-bp) + per-locus `pct_mapped` reporting
  · multi-set governance (set-size-weighted / Jaccard-deduplicated Hallmark scoring).

## H. DCIS / pre-malignant / EV lane

- ⛔ **HTAN Duke (Strand) DCIS atlas — the RNA + outcome instrument, ACCESS-BLOCKED.** Two uses: gene-level
  miR-29 targets (COL1A1/3A1/4A1/5A1, LOX/LOXL2, SPARC, FBN1) across normal→DCIS→IBC **and by recurrence
  outcome**; and the MH-48 abundance-asymmetry signature + **miR-185-5p** as a DCIS→recurrence predictor — the
  one question the public cohorts (3 events) cannot answer. Patient bulk sequencing is dbGaP/EGA-controlled
  (DAR), but **derived L3/L4 are often open and are sufficient** (L1/L2 not needed) — confirm the specific open
  files on the portal before committing. **Check periodically.**
- ⬜ **EV recipient-cell delivery — GSE297447** (the EV **mRNA** companion to GSE297448): are the targets of the
  *exported* miRNAs (miR-140-3p, miR-205) de-repressed in recipient cells / present as EV cargo? Closes the
  export→delivery loop the MCF10A EVs only half-answer.
- ⬜ **Rab27A / secretion-machinery axis** — does export-selectivity (cellular-down→EV-up) track
  Rab27A/secretory-pathway expression across progression? Separates selective packaging from bulk shedding.
- ⬜ **miR-185-5p** — the one arm surviving composition adjustment as an epithelial grade-decliner (MH-51,
  adj ρ=−0.75, q_adj=0.026). Cheap in existing data: its targets, Hallmark links, and whether it is a bona fide
  epithelial progression suppressor (vs the stroma-driven miR-29c).
- ⬜ **Stroma-aware compartment extension.** The framework (MH-1..47) is epithelial-centric, but miR-29c→ECM
  showed the invasion biology lives in non-epithelial compartments; the **myoepithelial-loss program
  (miR-145/-143/-126)** is the LATE invasion leader (MH-48) and deserves its own module rather than being
  adjusted away.
- ⬜ **The desmoplastic driver, if not miR-29** (MH-54's open end). ⚠ **The remaining clean test is IN-SITU**
  (spatial / single-cell CAF), **not another cultured-CAF set** — GSE196354 was inconclusive due to **culture
  drift** (FAP/ACTA2 only +0.15 NS, collagens flat), so cultured CAFs under-represent the in-vivo myCAF state.
- ⬜ **A DCIS-STAGED sorted/LCM small-RNA-seq accession** — the one gap the field does not fill. Scan persists.
- 🚫 **Do NOT** run a broad small-RNA-seq sweep, or acquire another plasma-EV cohort (the tissue-export
  signature does **not** translate to a robust circulating marker — miR-205-5p replicated in only 1 of 6
  cohorts, miR-140-3p in 1; plasma EV is multi-source). 🚫 Do **not** use the miRNA-inference tools (SiCmiR,
  STmiR) as evidence for the miR-29c→ECM axis — they infer miRNA activity *from target expression* and are
  **circular** for exactly this question; landscape/hypothesis only.


---

## J. ⭐ WHAT IS ALREADY WRITTEN IN A PATIENT'S NAT — the paired lane's unoccupied axis (planned 2026-08-03)

> ## ✅ §J RAN IN FULL 2026-08-04 — all six items closed (MH-228 … MH-233). **Verdict below; the item text
> that follows is the ORIGINAL PLAN, kept for provenance. Read this box first.**
>
> **The answer to the section's title: LESS than it assumed — at the LEVEL axis. What the paired design
> actually buys is the CHANGE axis, and there it holds up.**
>
> | item | verdict |
> |---|---|
> | **J1** | ⛔ premise **INVERTED** — tumours do not converge, **both layers DIVERGE** (mRNA +0.670, miRNA +0.329 log2), and the output diverges ~2× harder than the input ⇒ **the miRNA layer is not the dominant source of tumour individuality**. PAM50 explains ~nothing. |
> | **J2** | ⛔ the tumour does **not** exaggerate the patient — `b` = 0.183/0.137, **>1 for 0.0%/0.1%**, r² ≈ 2%/0.8%. Amplification **excluded**, not merely unobserved. The tumour's individuality is **NEW, not inherited**. |
> | **J3** | ⬜ headroom **not established** (p=0.058; the internal dose-response control fails). MH-162B stands. Power-limited; the design is the deliverable. |
> | **J4** | ⛔ pre-loading **does not survive its decoy** (REAL p=3.7e-06 → REAL−DECOY p=0.061). **No TSG enrichment.** Centrepiece closed negative. |
> | **J5** | ✅✅ **losses DO de-repress** — dy_gain −0.017 (p=9e-18), dy_loss +0.023 (p=4e-28), separation p=1.2e-24. **The lane's central quantity survives its own sign test.** The one clear win. |
> | **J6** | ✅ the (small) personal signal is **not** field/reverse causation (no gradient in `Cancer Epithelial`); ⚠ but 14q32 is where the constitutional setpoint is most **overwritten**, so that lever cannot adjudicate. |
>
> **⚠⚠ THE METHODOLOGICAL FINDING, and the most reusable thing here: THREE of the six items had a
> proposed estimator that was MECHANICALLY BIASED, and each would have shipped a false headline.**
> J2's `corr(nat_dev, own_shift)` shares `N` across both sides (+0.556 real vs +0.626 null, 100% of arms
> positive — reads as "the tumour amplifies the person in every arm"); J3's `corr(dy, head)` shares `N`
> too and is biased **toward** its own hypothesis, so it **could not fail**; J4's un-decoyed ρ is
> p=3.7e-06 and dies at p=0.061 against its matched decoy. ⇒ **on a paired design, write out the algebra
> of any statistic that touches both `N` and `T−N` before running it, and never report a per-patient
> quantity without its decoy.**
>
> **⛔ AND: FIVE of six items show NO role dependence** (J1 p=0.573 · J2 p=0.71 · J5 p=0.110 · J6 n/a ·
> J4 p=0.333). The "read it through the gene's cancer role" lens is **consistently null in this arc** —
> treat a role split here as a control that has already failed to fire, not as a live hypothesis.
>
> **➕ FOLLOW-UPS RUN SAME DAY (user-directed), MH-234 / MH-235:**
> **J7 — per-patient identity by LEVEL × PROFILE.** `patient_nat_identity` pools **all arms and only arms**;
> generalised to 3 levels × 11 profiles. `frac_own_wins` climbs 0.010 → 0.369 as the set narrows, which
> against a plain n-matched control looks like a 23×/4.2×/2.5× win — ⛔ **and dies completely against an
> EXPRESSION-MATCHED control** (0.252 vs **0.275**, 0.350 vs 0.269, 0.369 vs 0.285). The unmatched control
> gives 0.011 where the matched one gives 0.275 at the same n. **There is no hidden per-patient identity at
> any level or profile tested.** The control is now expression-matched inside the function.
> **J8 — the TRANSPOSE (within-patient, across 1,070 genes).** `rho_real` **−0.0580** vs `rho_decoy`
> **−0.0180**, gap **−0.0320, p=4.3e-11, real beats decoy in 78.6% of patients** — **the first per-patient
> quantity in this program to survive a matched decoy.** ⛔ But split-half over genes gives Spearman–Brown
> **r = 0.329**, so ~11% of the between-patient variance is reproducible ⇒ **a set-level quantity computed
> patient-wise, not a biomarker** — which is why the gap tracks nothing (field p=0.76, ΔC p=0.87, subtype
> flat). ⭐ **This closes MH-162A's explicitly untested object: the reliability of `real − decoy` is 0.33.**
> ⇒ **do not build a per-patient miRNA-realization biomarker on this cohort.**
>
> **⬜ What is genuinely still open:** the per-PATIENT lost set (J5 did the gene-level sign test only);
> Δpromoter-methylation as J6's third lever (70/103, staged); and ~~MH-162A's untested object~~ — **✅ CLOSED by MH-235: the reliability of `eff_real − eff_decoy` is 0.33.**

**The structural fact that motivates this section.** NAT enters this codebase in exactly THREE roles: a
differencing reference (`dX = Xt − Xn`), a nuisance (`dC`, residualised out of every ρ), and a cohort
median. **Never as a description of the patient.** And `_realize` correlates `pred` against `dy` over the
**103 patients** — so every ρ in the lane is a *cross-patient* correlation. The pairing cancels the
baseline and then discards the patient as a replicate. **`n=103` → per-edge SE≈0.10 → "set-level only" is
a consequence of that sampling choice, not a fact about the data**, which holds 103 patients × ~1,200
genes × ~5,600 edges.

**The reframe — a patient's NAT carries FOUR separable kinds of information the lane collapses into one
subtraction:**

| fold | what it is | status |
|---|---|---|
| **Constitution** | germline/imprinted setpoint — 14q32 DLK1–DIO3 = **56 HE arms → 169 HE genes**, NAT-detected *better* than tumour (50 vs 49 arms at >50% of patients) | cohort-level only |
| **Field** | how far the neighbourhood already moved — CIBERSORTx gives NAT a `Cancer Epithelial` fraction of **median 0.024, max 0.213** | computed, then residualised away |
| **Headroom** | the target's own NAT level bounds what the malignant step can do | ⬜ **zero code — `grep headroom` = 0 hits in 337 files** |
| **Architecture** | the tissue the tumour arose in — 113 NAT CIBERSORTx rows, **103/103 paired**, 0 missing | used only as a nuisance |

**The organising question: what does the malignant transformation PRESERVE, AMPLIFY, and ERASE?** — read
through the target gene's **role in cancer** (TSG/oncogene/dual), because that is where the impact is: a
patient's field determines *which route to cancer is already open in them*.

- ✅ **J1. Convergence vs divergence, n-matched — DONE, and the PREMISE INVERTED (MH-228, 2026-08-04).**
  There is **no convergence**. On the matched design (same 103 patients, same features, complete-case)
  **BOTH layers DIVERGE in tumour**, and the *output* diverges ~2× harder than the *input*:
  `log2(sd_tum/sd_nat)` = **mRNA +0.670** (1.59×, 1,436 genes) vs **miRNA +0.329** (1.26×, 221 arms), both
  outside all 500 within-patient label permutations. Only 8.1% of genes / 16.7% of arms are more variable
  in NAT. **PAM50 explains ~nothing** (within-subtype +0.640 / +0.308) ⇒ it is *within-subtype* dispersion.
  **No role dependence** (TSG +0.741 vs oncogene +0.807, p=0.573). ⛔ The old "sd 0.600→0.237" was a
  **naming collision + a per-cohort OBJECT mismatch** in `gauge.cohort_matrices`, not biology — the NAT leg
  reproduces, the TCGA leg does not, and the ordering is inverted (see MH-228; MH-102d's *verdict* is
  untouched because the state channel runs on z-scored β). ⇒ **retire "routes funnelling to one state"
  wherever it appears**, and note the consequence for the rest of §J: **the miRNA input cannot be the
  dominant source of the tumour's expression individuality**, so J4 must ask how *much* of the personal
  signal it carries, not merely whether it survives.
- ✅ **J2. Does the tumour exaggerate who the patient already was? — DONE. NO (MH-229, 2026-08-04).**
  ⛔ **And `corr(nat_dev, own_shift)`, the test this board asked for, is MECHANICALLY BIASED** — `N` sits on
  both sides with opposite sign, so it is positive by construction. Measured **+0.556 arms / +0.484 genes,
  100% / 99.9% of features positive — with a broken-pairing null of +0.626 / +0.537 and p_perm = 1.000.**
  It would have shipped "the tumour amplifies the person in every arm". **Do not use it anywhere.**
  ✅ The unbiased form is the LEVEL-on-LEVEL slope `T_dev = b·N_dev + e`: **b = 0.183 arms / 0.137 genes**
  (null ≈0, p_perm 0.005), **>1 for 0.0% / 0.1%**, median r 0.141/0.092 ⇒ the person's own NAT explains
  **r² ≈ 2% / 0.8%** of where their tumour sits. ⭐ Amplification is **excluded, not unobserved**: b plateaus
  at 0.22 by the 4th NAT-expression quintile, so `b_true=1` would need 78% noise in the best-measured arms,
  and genes show **no** dilution signature at all. No role dependence (p=0.71). ⇒ with J1: **the tumour's
  extra dispersion is NEW variance, not amplified personality — partially preserve + REPLACE.**
  ⚠ This is a level-axis result and does **not** touch the *change* axis: `b<1` means high-NAT patients are
  pulled DOWN toward a common tumour level, which is **J3's headroom mechanism in disguise.**
- ⬜➡️ **J3. Headroom — RUN, NOT ESTABLISHED (MH-230, 2026-08-04). MH-162B stands; do not re-run on this
  cohort.** ⚠ The test as stated below **could not fail** (headroom IS `N`, the outcome is `T−N`);
  replaced with ANCOVA `T ~ N + pred + pred×N_c + ΔC`. The discriminator (power is sign-SYMMETRIC,
  headroom is sign-ASYMMETRIC) points the right way — `a3_gain −0.148 < 0 < a3_loss +0.269`, `disc`
  −0.614 vs a clean split-permutation null of −0.059 — but the **pre-registered test gives p=0.058**, pooled
  `a3` is **+0.176 (wrong sign)**, and ⛔ **the internal dose-response control FAILS**: the effect does not
  scale with `nat_spread` (narrow −1.288 · mid −0.067 · wide −0.532, MW p=0.753) and narrow is strongest.
  Power-limited (an interaction on ~50 patients/half/gene; needs ≈970 genes with both halves, have 638).
  **The design is the deliverable.** Original text: **J3. Headroom.** A repressor gained where the target is already at floor **cannot** realize;
  repression needs the target present, de-repression needs transcription to resume — yet the lane treats Δ
  symmetrically everywhere. ⚠ `ctx_ceiling` is a FALSE FRIEND (an OOF R², not an expression level). ⭐ This
  **re-opens MH-162B's "acquired-but-unrealized = detection power" null with a competing mechanism** — keep
  MH-162B's power control in the design so the two explanations compete fairly.
- ✅ **J4. The per-patient vulnerability map — CLOSED, NEGATIVE (MH-233, 2026-08-04). Do not rebuild it
  on this cohort.** Sharpened onto the one surviving question (J2 voided the level axis; J5 validated the
  change axis): does a patient's own NAT pre-loading predict how that gene moves in THEIR tumour?
  3,298 matched REAL/DECOY edge pairs: **REAL +0.0103 (p=3.7e-06) · DECOY +0.0040 · REAL−DECOY +0.0063,
  p=0.061, frac<0 0.491.** ⇒ **the un-decoyed statistic would have shipped a p=3.7e-06 headline the matched
  control removes** — MH-162A's failure mode re-run, and why that condition governs §J. ⛔ **No TSG
  enrichment** (board question 1 answered NO, and backwards: TSG gap −0.0040, n=380, real edges *below*
  their own decoys, vs oncogene +0.0060, n=555, MW p=0.333).
  Original text: **J4. The per-patient cancer-gene vulnerability map (the centrepiece).** Per patient × cancer gene:
  **pre-loading** (NAT pressure) · **headroom** (NAT level) · **acquisition** · **realization**, split by
  `gene_roles`. Is the patient-specific component of NAT pressure enriched on **TSGs**? Does NAT
  pre-loading predict that TSG's realized repression in that patient's tumour? ⛔ **MH-162A's condition
  governs this whole item: every per-patient metric ships with its own site-free decoy or it does not
  ship** (a decoy matched the per-patient trait's reliability, +0.59 vs +0.50). ⚠ Note MH-162A never tested
  the reliability of `eff_real − eff_decoy` — the decoy-corrected quantity is a *different* object from the
  one that died.
- ✅ **J5. The opposite-sign control — DONE, and it PASSES (MH-231, 2026-08-04).** `dy_gain` **−0.0170**
  (p=9.1e-18) and `dy_loss` **+0.0230** (p=3.8e-28) across 1,003 HE genes — opposite signs, both
  significant, **de-repression the STRONGER leg**. Separation +0.0570 vs a broken-pairing null of
  **+0.0000**, paired Wilcoxon **p=1.2e-24**, opposite-sign in **57.0%** of genes (chance 25%). Replicates
  independently within TSG (n=46) and oncogene (n=72). **No role dependence** (cancer vs unknown p=0.110 —
  recorded as null, not as the 1.7× it superficially looks like). ⇒ **the pooled realization ρ is not a
  one-directional "gain ⇒ down" artifact; the paired lane's central quantity survives its own sign test.**
  ⬜ Still open from the original item: the **per-patient** lost set (only the gene-level sign test is
  done). Original text: **J5. The LOST repertoire + the opposite-sign control.** What a patient lost is only visible with
  paired data. Cohort `lost` exists (161 edges, ρ_adj −0.043); **no patient-level lost set, and no
  de-repression test on TCGA paired data** (the only one in the repo is GEO DCIS, ECM-only). ⭐ Losing a
  repressor predicts **de-repression — the opposite sign**. If gains repress but losses do not de-repress,
  the realization signal is not doing what its name says. **The strongest internal control the paired
  design can supply, never run.**
- ✅ **J6. Controls — DONE (MH-232, 2026-08-04). Reverse causation is NOT the explanation — via lever (b).**
  **(b) FIELD GRADIENT passes:** `b` within NAT `Cancer Epithelial` tertiles = low **0.171** · mid **0.221**
  · HIGH **0.134** — no gradient, high-field lowest. Field contamination would have to RISE with the field.
  ⚠ `Cancer Epithelial` is **not** in the C block (it is the malignant fraction) — read it from the raw
  deconvolution table. **(a) CONSTITUTIONAL RESTRICTION INVERTS and so cannot do its job:** imprinted 14q32
  arms have **LOWER** b (+0.101, n=13) than the rest (+0.186, n=208), MW **p=0.0082**, and it is **not** the
  abundance artifact (level-matched −0.0769, p=0.017, 85% below same-decile peers; matching removes only 9%
  of the raw gap). ⇒ **14q32 is where the malignant step most completely OVERWRITES the constitutional
  setpoint** — coherent with its near-uniform silencing (68% lose) — so the set is the place the signal is
  *weakest*, not where it survives, and the verdict rests entirely on (b). Coverage 13/52 (complete-case;
  depletion not differential, 25% vs 28%). ⬜ Δpromoter-methylation lever not run.
  Original text: **J6. Controls, because "field" and "reverse causation" are the same observation.** (a)
  **Constitutional restriction** — imprinted 14q32 arms are germline-set and cannot be tumour-induced the
  way a field effect can; if the patient-specific signal survives inside those 56 arms / 169 genes, reverse
  causation is weaker. ⚠ no curated imprinting list exists; `genomic_context.host ∈ {MEG3, MEG8, MEG9,
  DIO3OS, MIR493HG}` covers 52/56 and is a *proxy*. (b) **Field gradient** — `Cancer Epithelial` in NAT
  (max 0.213) as a covariate and stratifier, not a nuisance. Δpromoter-methylation (70/103, `N::` columns
  already staged in `mirna_promoter_betas.parquet`) folds in here as a constitution-vs-field lever.

🚫 **Explicitly NOT this section:** ΔRISC-as-moderator, isomiR-Δ, `M_reference` sweep — each is "an
existing per-arm axis with a Δ in front of it" (module application, not a paired question). And nothing on
the §F/§I dead list.

⚠⚠ **TRAPS THAT WILL SILENTLY CORRUPT ANY NAT ANALYSIS — read before writing code.**
`data_loaders.load_mirna_arms()`/`load_rna()` **average tumour and NAT** (see MH-223; now take an opt-in
`sample_type`) · `ago_gate`'s per-sample RISC capacity is that same average with **zero NAT keys**, so it
looks available for all 103 and is state-blind · `CPE`/`abs_purity`/`target_cn` are tumour quantities that
join happily onto NAT rows · missing miRNA is **NaN, not 0** (`.fillna(0)` turns an abstention into a
measurement) · **no technical replicate exists anywhere** (0 blood-normal, 0 NAT duplicates) ⇒ measurement
noise cannot be separated from biological retention, only bounded by `r_perm`.

## I. ⛔ Dead anchors — do not rebuild

One line each. The detail is in the registry row named.

- ⛔ **`βᵗ`, the CPTAC translational-repression latent** — the plan's centrepiece, **not supported at n=101**;
  all prior positive evidence was a **LEAK**. Four leak-free framings agree; the aggregate (1 df) gives BH
  q<0.10 in **1/17 genes — BCL2 alone**; **PTEN d=−0.006, p=0.82**. ⚠ **MH-34 is NOT overturned** (13/17 stay
  repression-directed, binomial p≈0.025) — the falsification is of the **modelling object**, not the biology.
  **There was never a βᵗ field. Only n fixes this.** **(→ MH-103)**
- ⛔ **Protein as a coupling lever** — it carries **4–6%** of the mRNA channel's Fisher information about β, and
  **≤7.6% at any `a_g`**. It must justify itself as a **different question**, never as a coupling gain. **(→ MH-108)**
- ⛔ **The state channel** — see §F. **(→ MH-102d)**
- ⛔ **"TCGA is sparse so it cannot arbitrate"** — MH-76's **frozen-panel** test removes overfitting as an
  explanation entirely: **DFI +0.002, OS +0.006; panel-alone C 0.48–0.53 ≈ random.**
- ✅ **Pressure magnitude is exhaustively non-prognostic** — 0 recurrence FDR across 3,368 features. A
  **well-earned negative**, credible because **miR-210 → DRFS reproduces Buffa 2011** as a working positive
  control. Functional > magnitude (realized +0.060 vs received-pressure −0.016, a *within-cohort* contrast).
- ✅ **CIBERSORTx / the panel question is SETTLED — ZERO new runs.** β is **ρ=0.94** across entirely different
  atlases (Wu-9 vs HBCA-11), ρ=0.99 across reference resamples, but only ρ≈0.58–0.81 vs **no C at all** (mean
  |β| halves) ⇒ ***which* composition control you use moves β ~10× less than *whether* you control it.** Wu-9
  stays the panel for all four cohorts. **Do NOT run** a fine 42–49-way panel, TCGA×HBCA, or refsample
  "harmonisation". Two gaps are **not** deconvolution problems: **adipose** (no reference has an adipocyte class
  and GTEx breast is adipose-dominant ⇒ needs a marker metagene) and **immune C-completeness** (Wu-9 lymphoid
  tracks Thorsson r=0.31 vs LM22's 0.70 — ⚠ LM22 + Thorsson are **validators, never covariates**).
- ✅ **The n≈1000 ceiling** — every internal lever (CN channel, δ-pooling, isomiR refit, Student-t, cross-gene
  pooling, the A3/A4 field) landed **"immaterial at n≈1000"** ⇒ **the frontier is exogenous, not more internal
  refinement.** Don't re-open them.
- ✅ **Engine: Gibbs, not HMC** — measured, not aesthetic. ⚠ **Carve-out:** fit the *gauge* with NNLS (Gibbs's
  heavy-tailed posterior SDs break the errors-in-variables correction). **Only a future non-conjugate channel
  reopens this** — binding = Poisson/NB, methylation = Beta.
- ✅ **`_rtnorm_pos` sampler bug — FIXED, don't re-litigate** (contaminated 3.15% of β; Robert (1995) rejection
  sampling, KS-verified). ⚠ It **overturns MH-119's recorded cause**: `share`=43.7 was this bug, not
  "anti-repressive biology".
- 🚫 **Occupancy / ceRNA saturation substrate** (Decision C/G) — failed its held-out gate; Phase-4 cooperativity
  sits on it (gate HARD on held-out if ever revisited). Likewise the `priors.py` spike-and-slab π/μ/τ object
  (permutation-FDR discovery already *is* the inclusion inference), the A/D/F pooling partials (verified no
  gain), and any channel whose link isn't identifiable.
- 🚫 **Fusion honesty-layer readouts** (axes U/AO/AP) and the **binding channel** (already a prior, `priors.μ` —
  the fusion frame only renames it): free readouts, but they **materialize only if a 2nd non-redundant channel
  lands**. Per §A and this section, none is on the board.
- ⬜ **The one design partial still worth having: B, the transcription-rate proxy** (the design's "most
  load-bearing"). An **intronic / pre-mRNA-read** nascent-transcription proxy as a C column, so a positive M
  isn't a TF co-driving miRNA and target. **HIGH value but RISKY** (over-control — the TF-regulon proxy failed
  exactly this way) and **data-blocked** on intronic reads. Check availability first.

---

## Open questions worth a dedicated study (not scheduled)

- 🔬 Is the basal **immune-pressure tone** a *cause* (miRNA-driven IFN silencing) or a *consequence* (low
  infiltration → low immune gene/miRNA co-program)? ⚠ post-MH-107/114, any version of this must run under the
  composition C and report retention, or it is compartment arithmetic by construction.
- 🔬 Does **8q24 amplification** (miR-151a, AGO2) create a coordinated "amplicon-driven RISC + guide" state that
  amplifies pressure non-linearly?
- 🔬 Are the **differentiation-program** repressions (estrogen / adipogenesis / myogenesis) a generic
  loss-of-identity signal shared with other carcinomas (pan-cancer test)?
- ⬜ **Per-sample hub-route OLS explainability** (`gene_expression_explainability.py`, in progress): nested OLS
  on `log2(TPM+1) ~ CN + CPE + PAM50 + promoter meth [+ route miRNA block]`. Early Basal read: CDKN1A ΔR²≈0.16
  (all) / 0.05 (rna_low); PTEN ≈0.11 / 0.05; **IRF1 curated route ≈0.006 but the TargetScan alt (miR-130b/301)
  ≈0.086 / 0.30** ⇒ the IRF1 grant route should use `curated_plus_targetscan` (combined Basal ΔR²≈0.10 / 0.32).
  ⚠ The IRF1 route **failed proliferation adjustment within Basal** (MH-19) and is exploratory until this
  resolves it.

---

## ✅ Closed 2026-08-03 (MH-202 / MH-203 / MH-204) — three board items, and one caveat retired

- ✅ **`β(TCGA) → Buffa`** (was the #1 next-thing on `STATE_OF_PLAY`). Dissociation: **LEVEL transports,
  RANK does not separate β from abundance.** After the MH-204 stale-decoy fix it also **beats a fitted
  fake, −0.0374, p=3.5e-06.** → `learned/eval/ood_cohort.py`.
- ✅ **The learned-β OUTCOME arm** — per-gene yes, **aggregate no, and the weights earn nothing.**
  Positive-control gated on miR-210→DRFS. TCGA is the weaker cohort (OS only). → `beta_prognostic.py`.
- ✅ **The PROTEIN layer, properly powered** — RPPA at **n=866**, same patients, no cohort boundary.
  **Destabilisation, not translation.** → `rppa_protein.py`.
- ⭐⭐ **RETIRED CAVEAT: *"every win is FITTED β vs UNFITTED abundance"*.** A fitted, matched, site-free
  decoy now exists on all three layers. **β wins on mRNA in TCGA (−0.0189, p=6.7e-09) and in Buffa
  (−0.0374, p=3.5e-06); it does NOT win on protein (p=0.554).** ⇒ **a LAYER split, not a general
  negative** — and the mRNA win survives the subtype control (within-PAM50 gap −0.0154, retention 0.81).

✅ **AND THAT QUESTION IS NOW CLOSED — SAME DAY, AND IT CORRECTS THE ROW ABOVE (MH-205).** *"Is the protein
null biology or instrument?"* → **INSTRUMENT: gene selection and power — and NOT the phospho contamination
I had named as the suspect** (13/180 genes; dropping them moves −0.0086 → −0.0089).
**On the same 180 genes the mRNA gap is also n.s. (−0.0103, p=0.0565), both sit inside their own MDE, and
paired on gene at 5× the gene count the layers are indistinguishable — CPTAC `tcga105` +0.0018 (p=0.287),
`prospective` +0.0056 (p=0.341).** ⇒ **the "layer split" is WITHDRAWN; the protein half is UNRESOLVED.**
A consistent direction across all three instruments (Stouffer z=+2.44, anti-conservative) leaves a **small**
penalty ~0.002–0.006 ρ plausible. ✅ Destabilisation-not-translation survives, now on 849/897 genes.
→ `learned/eval/layer_boundary.py`.

**⬜ What THIS opens.** The binding constraint is now explicit, and it is **genes with a protein readout, not
patients**: ~514 genes are needed for the mRNA-sized effect, and no cohort supplies both a large gene count
and a large patient count (RPPA 180×866; CPTAC ~880×105). The honest options are (a) a gene-powered proteome
cohort, or (b) accepting the protein layer is **bounded** — |gap| ≲ 0.006 — which is itself a real constraint
and far cheaper than chasing it.

---

## ⬜ Ingest the versioned COSMIC Cancer Gene Census as a real annotation (user-directed, 2026-08-19)

**Why now.** Column-review unit 14 (MH-250) measured what `role` actually covers: **`unknown` on 1,262 of
1,420 genes — 11.1% coverage**. The source is `gene_roles._builtin_table()`, **232 genes hardcoded as Python
sets**. Per `gene_roles.py` it is a **COSMIC-CGC / OncoKB consensus, breast-cancer-prioritised** — so the
problem is *not* that the list is arbitrary. The problem is that it is **frozen, unversioned and subset**:
no census version, no download date, no refresh path, and no way to tell what a newer census would add.

**What to do.**
1. **Download the CGC** (`cancer_gene_census.csv`) from the COSMIC download portal. ⚠ **Requires a free
   academic registration and carries a non-commercial licence** — record the **census version and download
   date** in `docs/DATA_SOURCES.md`, which is half the point of the exercise.
2. **Ingest to `annotations/cosmic_cgc.tsv`** with a builder under `scripts/annotations/`, following the
   `_build_tf_census.py` / `humantfs_lambert2018_tf.tsv` pattern already used for the Lambert TF census —
   that is the in-repo precedent for a licensed curated list, so do not invent a new shape.
3. **Map CGC's `Role in Cancer` field** (`oncogene` / `TSG` / `fusion`, comma-separated and often multiple)
   onto the existing 3-value vocabulary. ⚠ **`fusion` has no home in the current scheme** and must not be
   silently folded into `oncogene`. **Keep Tier 1 and Tier 2 distinguishable** — Tier 2 is weaker evidence
   and a reader must be able to gate on it.
4. **Wire it as the OVERRIDE, not a code edit** — `config.GENE_ROLES_OVERRIDE` already exists for exactly
   this and `load_gene_roles` merges it on top (override wins). No change to `gene_roles.py` is needed.

**⛔ Before adopting it, run these two checks — a bigger list is not automatically a better one.**
- **CONCORDANCE ON THE OVERLAP.** Where the built-in and CGC both call a gene, do they agree? Report the
  disagreement rate and the disagreeing genes **by name**. The built-in is BC-prioritised and CGC is
  pan-cancer, so a systematic split is plausible and would be a finding, not a nuisance.
- ⭐⭐ **AND RE-RUN EVERY `role`-STRATIFIED RESULT ON BOTH ANNOTATIONS.** Coverage goes from 158 genes to
  whatever CGC supplies, which changes the *universe* of every oncogene-vs-TSG contrast — so a moved result
  is confounded between "better annotation" and "different genes". **Report both, and treat a result that
  appears under only one of them as unresolved.** This is axiom 4: state what the change would do to the
  test before running it.

**Downstream to re-check** (grep `malignancy_sign` / `load_gene_roles`): `geneset_architecture`
(`sum_mal_pro_tumor`, `sum_mal_pro_tumor_cont`), `mal_pro_tumor_hier`, and the edge card's `role` column.
⚠ The **continuous** counterpart `load_gene_dependency` (DepMap, ~18.5k genes) is unaffected and is the
better instrument where a binary role is not required — **CGC does not replace it.**

---

## ⬜ EXAMINE: LEFTY2 / miR-302 — the minimal case where β, identity and sd_arm tell three different stories

**User-flagged 2026-08-19** ("interesting, should be noted for future examination"), surfaced by
`ratio_blowup_audit` as the largest `sd_arm` on the edge card. It is worth a look because it is **the whole
magnitude-vs-identity question compressed into three rows**, on a biologically coherent pair.

| arm | β | `beta_arm` | `sd_arm` | z | `identity` | `coupling_tum` |
|---|---|---|---|---|---|---|
| hsa-miR-302a-3p | **0.91840** | 1.093 | **1.302** | 0.94 | **0.0** | −0.154 |
| hsa-miR-302d-3p | **0.91885** | 1.056 | **1.097** | 1.01 | **0.0** | −0.104 |
| hsa-miR-373-3p | **0.03470** | 0.022 | 0.014 | 1.53 | **1.0** | +0.103 |

**What makes it a clean test case.** LEFTY2 has exactly **3 edges, all in ONE seed family**
(`miR-302-3p/372-3p/373-3p/520-3p`, `n_arm_in_cell = 3`) — so there is no cross-family confounding to
unpick. Within it:
- **β says** miR-302a and miR-302d dominate, at values agreeing to **four decimals** (0.91840 / 0.91885).
- **`sd_arm` says** neither is identified: the posterior SD (1.302, 1.097) **exceeds the estimate**, and
  1.302 is close to the largest β anywhere on the card (1.521). z ≈ 1 for both.
- **`identity` says** miR-373-3p is the *entire* story (1.0) while the two high-β arms get **0.0** — the arm
  with a **26× smaller** β takes all the attribution credit.

**The questions to answer, in order.**
1. ⭐ **Is `identity = 0.0` here CORRECT or DEGENERATE?** Shapley/LMG is collinearity-fair by design, so two
   indistinguishable arms *should* split credit — but splitting should give ~0.5 each, not 0.0 each. A
   Shapley value of exactly 0 for the two largest coefficients wants explaining before it is trusted.
   ⇒ read `identity_coherence` and the NNLS weights for this gene; check whether the pair cancels rather
   than splits (memory `ratio-readouts-need-a-denominator-gate` — a value outside a model's own support is
   a bug report, not a finding).
2. **Does the sign disagreement matter?** `coupling_tum` is **negative for both miR-302 arms** (−0.154,
   −0.104) and **positive for miR-373-3p** (+0.103) — the arm receiving all the identity credit is the one
   whose observed coupling runs the *wrong way* for repression.
3. **Is the biology real?** miR-302 is the ESC/pluripotency cluster and LEFTY2 is a Nodal antagonist in
   TGF-β signalling — a coherent pairing, which is exactly why an unidentified β here is worth resolving
   rather than filtering away.

⚠ **Do not "fix" this by gating `sd_arm`** — it is not a ratio and its tail is real (MH-257). The wide
posterior is the *finding*: it says the design cannot separate these two arms, which is a statement about
identifiability, not noise. ⇒ this is a candidate for the **arm-resolution** machinery
(`armres_*`, `arm_sep_z`) and for the within-family refit that is default-OFF.

---

## ⬜ DIAGNOSE: `cptac_t105` breaks the edge→gene no-op identity on 29% of single-arm genes

**Found 2026-08-19** by the no-op test (`CARD_RUNG_DOCTRINE.md` §2b), which exists because of this question.

On a gene with exactly one arm, `cptac_{cohort}_agg_rho_prot` (gene card, a β-weighted sum over the gene's
arms) is a weighted sum of ONE term and must equal `cptac_{cohort}_rho_prot` (edge card) **exactly**.

- **prospective: 0 of 466 disagree** — the identity holds perfectly, so it is the correct expectation.
- **tcga105: 112 of 389 disagree (28.8%)**, median |diff| **0.0587**, p90 0.2374, **max 0.6446**.

**Already ruled out** — do not re-test these:
- Not a sample-count difference: edge-side n is **105** on both the matching and mismatching sets.
- Not a mis-specified no-op set: all 112 have `n_arms == 1` **and** exactly one scored edge.
- Not a general cptac defect: the prospective cohort, same code path, is exact on all 466.

**Where to look:** `learned/cptac_card.py`'s tcga105 branch — specifically whether the gene-level and
edge-level scorers see the same patient set, the same C block, and the same layer matrix. The asymmetry
against `prospective` is the strongest clue: whatever differs is **cohort-specific**, not estimand-level.

⚠ **Until resolved: do not cross-reference edge-rung and gene-rung `cptac_*` values on tcga105.** The
prospective cohort is unaffected. ⚠ Note this is *separate* from MH-253 (the `NaN < NaN` mask on
`agg_beats_abund_prot`) and from the 2026-08-04 card staleness — both were checked and neither explains it.

---

## ⬜ STAMP THE BASE CARD — a freshness contract, not a merge (user-asked 2026-08-19)

**The question was whether the base card can be retired and everything collapsed into one card design.**
Measured first: it cannot, and the reason is worth recording so nobody re-opens it.

| layer | cols | owns | cost to rebuild |
|---|---|---|---|
| **base** `edge_card_base.tsv` | **84 of 186 (45%)** | the **FITTED** quantities — `coupling_` 12, `beta_` 7, `dose_` 7, `share_`/`rank_`/`d_` 12, `retention_` 3 | multi-worker per-gene build over 1,549 genes (warm-then-fork) |
| **annotation** `card_ladders` / `card_context` | **102 (55%)** | context and joins — `cptac_` 20, `ctx_` 18, `esub_` 8, `cal_` 7, `adm_` 4 | seconds; additive, aborts if a pre-existing column changes |

⇒ **one card design would not simplify anything** — it would move the expensive per-gene fit into the layer
that currently runs in seconds on every annotate. **The base is a cache of the expensive thing, and that is
the right shape.**

⛔⛔ **THE REAL DEFECT IS ONE LEVEL DOWN: it is a cache with NO FRESHNESS CONTRACT.**
`edge_card_base_provenance.tsv` records what each column *is* (column → block → estimator, 93 rows) but
**nothing about WHEN it was built or FROM WHICH INPUTS**, so nothing can invalidate it. That single gap
produced four separate defects already recorded:

- **MH-252** — 129 genes carry model output with zero edge rows (gene card 1,549 vs edge card 1,420).
- **MH-266** — prunes that printed `✅ DROPPED` and did not take, because `_annotate` re-adds base blocks.
- **MH-269** — the renames needed a post-annotation mechanism for the same reason.
- The delivered card silently mixes a **2026-08-04** base with same-day annotations.

**THE FIX — a stamp, not a merge:**
1. On `canonical_card.build()`, write `edge_card_base_manifest.json`: the build timestamp, the mtime + size
   of every input it read, the universe size, and a hash of the gene list.
2. In `_annotate`, compare that manifest against the inputs' current state. **Refuse (or warn loudly and
   stamp the output) when the base is older than an input it was built from.**
3. Surface the base's vintage on the card browser and in the dossier, so a reader never has to infer it from
   an mtime.

⚠ **Do this BEFORE the rebuild, not after** — the rebuild is the moment the manifest is cheapest to write,
and without it the next drift is silent again.

---

## ⬜ WRITE THE 79 MISSING COLUMN DESCRIPTIONS (surfaced by the reference, MH-273)

`gen_column_reference --gaps` lists **79 columns with no column-specific glossary entry and no decodable
token** — they fall back to their block's description, which is why the reference badges them
`block-level`. This is the writing queue, and it is short enough to clear.

⭐ **Start with `ctx_ceiling`.** It is the measurability ceiling that MH-144 turns on, it appears on the
edge, gene and gene_family cards, and it has **never had an entry of its own** — every reader has been
getting the `ctx_` block text. Its siblings `ctx_measurable` and `ctx_collapse` are in the same state.

Others worth doing early because they are quoted: `healthy_leg` / `healthy_potential` /
`healthy_uninformative` (the healthy-leg triage), `cal_identified_hi` / `_lo` (the calibration band),
`echim_manakov_w` / `echim_tarbase_w` (evidence weights), `gene_net_repr` / `gene_dominated`,
`term_WIRING` / `term_INTERACT`, `realization_score`, `frac_identified`.

⚠ **Write the LEAD first and keep it short** — the reference splits at the first ⛔/⚠/⭐ marker, so a lead
that opens with a caveat renders as a caveat. One plain sentence saying what the column *is*, then the
warnings.

⚠ **Do not batch-generate these.** The decoder already covers everything mechanically decodable; what is
left needs someone to read the producer. A generated sentence here would be indistinguishable from an
authored one in the badge, which is precisely the confusion the badges exist to prevent.

