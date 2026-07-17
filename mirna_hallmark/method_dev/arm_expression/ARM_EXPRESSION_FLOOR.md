# Arm expression-detectability floor (canonical)

**Module:** `mirna_hallmark/arm_expression.py` · **Config:** `ARM_EXPRESSED_MIN_RPM=10`,
`ARM_EXPRESSED_MIN_FRAC=0.01` · **Table:** `output/matrices/arm_expression_tiers.tsv`.

## Problem

Having a curated HE edge ≠ being expressed in TCGA-BRCA. Of 789 HE-edge arms: **73 aren't in the TCGA
miRNA matrix** (no data); of the 716 with data, **345 (48%) are virtually unexpressed** (cohort-median
RPM<1) and only **220 (28%) are robustly expressed** (median RPM≥10) — though those 220 carry **63% of
HE-edge degree** (the hubs are the expressed ones). With `PRESSURE_ABUNDANCE_FLOOR=0`, all silent arms were
nominally included, inflating regulator counts (`n_eff`) and any non-abundance-weighted construction.

## Are the silent arms ever functional? (two rescue mechanisms, tested)

1. **A co-induced silent group occupying RISC — NO.** Median-silent arms are **≤2% of the total miRNA pool
   even in the single most extreme tumor** (per-tumor median = 0 arms at RPM≥10; one outlier reaches 28 arms,
   still ~2%). RISC is dominated by the high-abundance hubs in every sample; a silent group cannot co-rise
   enough to compete.
2. **A low miRNA controlling an even-lower target (stoichiometry) — NOT MEASURABLE here, and invisible
   anyway.** miRNA-seq RPM and mRNA-seq TPM are separate libraries with no copy-number calibration, so the
   stoichiometric ratio can't be computed. And the only cohort-level test of impact is coupling
   (anti-correlation), which a near-flat low-abundance arm cannot produce (≈zero cross-sample variance).

So the silent arms are inert at the cohort level — but the cut must not pretend to judge *function*.

## The floor — DETECTABILITY (max-based), not cohort-median

A cohort-median cut would wrongly discard context-specifically induced arms. Canonical rule:

> **KEEP** an arm if it reaches **RPM ≥ `ARM_EXPRESSED_MIN_RPM` (10)** in **≥ 1 tumor** (max RPM);
> **REMOVE** only arms that **never** reach the floor in any tumor.

So 111 median-silent arms that *spike* >10 RPM in a few tumors (<1% of the cohort) are **kept** (induction in
a tumor subset can be real there); only the never-detected are dropped. `frac_expressed` flags the weakest.
This is a **noise filter, not a functional verdict** — surviving low arms self-attenuate in the
abundance-weighted aggregate, and their function is decided by coupling. See `SILENT_ARM_REMOVAL.md`.

**Tiers** (`arm_expression_tiers()`): **robust** (median RPM≥10) > **conditional** (max RPM≥10 but median<10
— context-specific/inducible) > **silent** (max RPM<10 — never reaches the floor; removed). Over all 2,236
arms in the matrix: robust 242 · conditional 476 · **silent 1,518** → 718 expressed. Over the 789 HE-edge
arms: **482 expressed (220 robust + 262 conditional) / 307 silent** (234 never-detected + 73 absent from the
matrix).

## Sensitivity — headline UNCHANGED on the expressed-arm universe

`tumor_acquired_coupling(expressed_only=True)` (canonical confounds + prolif), expressed vs full:

| predictor | full med ρ(pos) | expressed-only med ρ(pos) |
|-----------|-----------------|---------------------------|
| base | −0.161 | −0.162 |
| **acq** | **−0.182** | **−0.182** |
| **acq_ev** | **−0.193** | **−0.193** |
| val_a0 (evidence) | inert | inert |
| val_a1 (promiscuity) | harmful | harmful |

`acq`/`acq_ev` still win, promiscuity still harmful, weighting still inert; base is unchanged (kept low arms
self-attenuate). **The floor changes no conclusion; it just removes the never-detected noise arms.**

## API / usage

- `arm_expression_tiers()` → per-arm DataFrame (median/p90/max log2RPM, frac_expressed, tier, expressed).
- `expressed_arms()` → frozenset of arms passing the floor.
- `filter_expressed_edges(edges)` → drop silent-arm edges (locus-suffix tolerant).
- Landscape (`he_edge_arm_landscape.py`) tags every arm's tier; silent arms drawn as grey rings.

**Stance:** apply as a default noise filter for non-abundance-weighted constructions and for `n_eff`;
keep it optional/reported for abundance-weighted analyses (where silent arms self-attenuate). Never use it
to claim a low arm is non-functional — that's a coupling question, not an abundance one.
