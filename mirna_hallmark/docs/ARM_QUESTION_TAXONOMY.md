# Arm-level question taxonomy — the arm, the arm-within-its-family, and the family itself

> **Scope: the ARM rung and its two neighbours.** A miRNA arm `a`, what can be asked about it alone,
> what can be asked about it *relative to its seed-family mates*, and what belongs to the family as an
> entity. Companion to `EDGE_QUESTION_TAXONOMY.md` (one miRNA→gene pair) and
> `GENE_QUESTION_TAXONOMY.md` (a gene's total incoming pressure); the unit bookkeeping is in
> `CARD_RUNG_DOCTRINE.md`. Numbers live in `DISCOVERY_REGISTRY.md`.

## 1. ⭐ The lemma — **the arm is a MEASURED unit but not an ESTIMATED unit**

`canonical_M` **always** fits on the **seed-family** support: `_state_family_data` →
`FAM.collapse_by_family` → `_bagged_nnls`. `arm_level=True` merely **broadcasts** each family's weight to
its member arms, and `states.py` says so explicitly — the per-arm split is

> *"a NOMINATION (abundance apportions the realized contribution), **not a resolved coefficient**"*

That one fact partitions every question at this rung, and it is why the edge and gene lemmas have no
counterpart here: **there is no arm-level estimator to characterise.** What varies is whether the thing
you are asking about is *measured on the arm* or *allocated to it*.

## 2. The A / B / C partition — decide this first

| | the question | rung | machinery | default answer |
|---|---|---|---|---|
| **A** | anything **measured from the arm alone** — detection, abundance, floor, spread, genomic context, AGO loading, isomiR composition, promiscuity | true **arm** | `arm_card.tsv` (key `[arm]`) | answerable |
| **B** | **"which arm in this family does the work?"** | **`arm-in-family`** — *not* the arm rung | `arm_resolvable`, `arm_sep_z`, `arm_dbeta`, `arm_credit_share`, `oof_rho_arm` vs `oof_rho_fam`, `oof_drho`, `n_arm_in_cell`, `family_dose_share`, `family_role` | ⚠ **often "not resolvable"** — treat that as the answer, not a failure |
| **C** | properties of the **family as an entity** — size, dose concentration, seed heterogeneity, collapse validity, redundancy | **seed_family** | `seed_family_card.tsv` (key `[seed_family]`) | answerable |

⛔ **The B columns were caught empirically** (MH-166): they vary *inside* a `(gene, family)` cell, so they
describe an arm **within its family**, not the family as a whole. Putting them on the family rung was the
bug.

## 3. Axis Group I — **can this arm carry a claim at all?** (lead with this)

⭐ **This is the only rung where the MEASUREMENT axis is primary rather than a caveat.** MH-247 measured
the NAT→tumour retention signal to be a **reliability gradient**: every strong predictor was
abundance-flavoured (`detection` ρ=+0.356 · `arm_pct_floor` +0.375 · `arm_med_rpm` +0.362 ·
`grank_TUM` +0.359), and the one axis provably orthogonal to abundance predicted **nothing**.

| column | reads as |
|---|---|
| `detection` · `arm_pct_floor` · `spiker` | is it seen, and is it seen as a spike in a few samples? |
| `arm_med_rpm` · `arm_iqr` · `grank_{HLY,NAT,TUM}` | level and spread, absolute and as a within-sample percentile |
| `wshift_sd_dGlobalRank` | how much its global rank swings **across patients** |

⇒ **Establish this before any biological reading.** A ranking of arms by almost any outcome will
reproduce this axis unless you control for it.

## 4. Axis Group II — arm-alone properties

- **Genomic context** (`MIRNA_GENOMIC_CONTEXT_AXIS.md`): `mir_class` ∈ {intergenic, sense_coding_host,
  sense_lncRNA_host} · `host` · `host_type` · `strand_rel` · `n_loci`. ⚠ Measured MH-246:
  hosted-vs-intergenic retention is 0.111 vs 0.094 **after composition adjustment** — much weaker than the
  raw 0.125 vs 0.028, so most of that split was composition and/or rung.
- **AGO loading**: `ago_loading`, `ago_loading_measured`.
- **isomiR composition** (`ISOMIR_AWARE_MODELING.md`): 5′-dominant fraction, per-sample QC.
- ⭐ **Promiscuity** (`gene_axes._promiscuity_map`) — **the one arm axis that is abundance-ORTHOGONAL**
  (measured corr with abundance −0.011). It is therefore the *only* clean way to separate a biological
  arm effect from the reliability gradient of §3. **It uses the SEQUENCE targetome, never the curated
  one** — the curated degree (`he_degree*`) is a FAME axis (ρ=+0.736 with PMID count).

## 5. Axis Group III — arm-in-family resolution (rung B)

The question is not *"what does this arm do"* but *"can this arm be distinguished from its family-mates
at all?"* — and the machinery exists because usually it cannot.

`arm_resolvable` / `arm_sep_z` (is the separation real?) · `oof_rho_arm` vs `oof_rho_fam`, `oof_drho`
(does modelling the arm beat modelling the family out-of-fold?) · `arm_credit_share`, `arm_dbeta`,
`family_dose_share`, `family_role` (how the family's realized pressure is apportioned) · `n_arm_in_cell`.

⚠ **One-unit degeneracy.** Where a family has a single arm in the cell, arm ≡ family **by construction**
and every weighted-vs-unweighted comparison is inert. That is ~43–50% of the universe;
`gene_axes.mask_degenerate` before any such comparison.

## 6. Axis Group IV — the family as an entity (rung C)

`seed_family_card.tsv` is gene-free and carries: membership (`fam_n_members`, `fam_n_expressed`), dose
shape (`fam_dose_{total,med,max}_log2`, `fam_dose_hhi`, `fam_dominant_arm`, `fam_dominant_share`),
variability (`fam_var_{med,max,hhi}`), **seed heterogeneity** (`fam_seed_shift_max_gated`,
`fam_n_shift_gt20_gated`), targetome (`fam_targetome_{seq,seed}_med`), evidence/fame
(`fam_fame_npmid_sum`, `fam_cur_degree_max`), guide/passenger (`fam_n_guide`, `fam_n_passenger`),
context homogeneity (`fam_ctx_*`) and `fam_redund_m_eff`.

⚠ **`gene_family_card.tsv` is NOT this rung** — it is keyed `[gene, family]` and cannot express a property of
the family itself. That distinction is why `seed_family_card` exists (see `CARD_RUNG_DOCTRINE.md` §1).

## 7. Traps with a recorded cost

- ⛔ **Within-gene shares are EDGE rung, never arm.** `share_HLY_meas`, `rank_HLY_meas`,
  `d_rank_HLY_TUM_meas`, `n_HLY_meas`, `share_*`, `rank_*` come from `states.budget_shift`, which is
  **per gene**. They cannot be read as arm properties.
- ⛔ **NaN there means "GTEx cannot see this arm" (abstention, excluded from the denominator), NOT "owns
  none of the healthy budget".** Confusing the two *was* the bug (MH-210).
- **An arm-rung column must be constant within arm ACROSS genes** — if it varies by gene it is not an arm
  property whatever its name suggests (MH-214), and `card_rungs --check` will say so.
- **Bare-stem arms conflate 5p/3p.** `arm_card._arm_key` drops them deliberately; matching a bare stem
  onto a suffixed family is exactly the conflation the key exists to refuse.
- ⬜ **Pending: `.N` arm-name normalisation** recovers ~18 HE arms + ~300 edges. Apply at the next full
  refresh (memory `universe-redefinition-pending-refresh`).
- **Complete-case / detection filtering selects the ranking.** 5 of the top-10 retention arms were
  excluded by complete-case (n=75–102): lower n ⇒ noisier estimate ⇒ extremes enriched. Report rankings
  restricted to full-n arms (`STRATIFICATION_PLAYBOOK.md` §5).

## 8. ⛔ Settled — do not rebuild

| | |
|---|---|
| **AGO loading (arm)** | built and validated, and **coupling-INERT**. Its value is identity/QC, not coupling. |
| **isomiR-aware refit** | **BUILT, default-OFF for coupling.** Bidirectional, not a systematic improvement — the recurring "real-but-immaterial-at-n≈1000" pattern. Should move only the ~17 seed-heterogeneous + ~34 dose-partitioned families. |
| **family collapse** | valid for bulk **except** those ~17 seed-heterogeneous families, where a fraction of reads target a different genome-wide targetome. |
| **arm-rung retention** | a **measurement-reliability gradient**, settled by the promiscuity orthogonal control (MH-247) — not a biological ranking. |
