"""CARD GLOSSARY — what each column MEANS, as opposed to what unit it lives on.

    .venv/bin/python3 -m mirna_hallmark.learned.card_glossary            # coverage report
    .venv/bin/python3 -m mirna_hallmark.learned.card_glossary --emit     # write card_glossary.tsv
    .venv/bin/python3 -m mirna_hallmark.learned.card_glossary --col beta # one column, every card

⭐ WHY THIS EXISTS. `card_rungs` answers *what UNIT is this number on*; `domain` answers *which ROWS is it
defined on*. Neither answers **what the number IS**, and until now that lived only in the producing
modules' docstrings — `CARD_RUNG_DOCTRINE §2` records the gap explicitly. A reader of
`docs/derived/cards/*.html` could see a value, its rung and its coverage, and still not know what it meant.

⛔ **LOOKUP IS CARD-SCOPED, ALWAYS.** A `(card, prefix)` entry beats a global `("", prefix)` one. Global
prefix matching is the exact defect MH-220 fixed in `card_rungs.DOMAIN` and `card_ladders.PREFIXES` after
it produced THREE separate bugs (`sub_`, `chim_`, `comp_`), because two cards independently grow blocks
with the same prefix and one card's metadata then silently describes the other's columns.

⚠ **AN ABSENT DESCRIPTION IS RECORDED AS ABSENT, NEVER GUESSED.** `describe()` returns None and the
coverage report counts it. A plausible-sounding invented description is worse than a blank one: it reads
as authority. Fill gaps from the producing module, not from the column name.

⚠ **THE RUNG STILL TRAVELS SEPARATELY.** A description here never states the unit — read `card_registry`
for that. `beta` is described once, and it is edge-rung on `edge` and family-rung on `gene_family`.

⛔⛔ **DUPLICATE KEYS ARE CHECKED, BECAUSE A DICT LITERAL SWALLOWS THEM SILENTLY.** A repeated key is
resolved at construction with the LAST one winning — no error, no warning. On 2026-08-19 that killed four
freshly-written `gene_family` entries (superseded stubs sat later in `COLUMNS`), and the guard written in
response IMMEDIATELY found a second, older instance nobody had noticed: a circular stub
*"Shapley/LMG identity attribution — see the `identity_` block"* had been overwriting the real `identity_`
description for every `identity_*` column. **A silent overwrite in a documentation table is worse than a
missing entry, because the stale text still reads as authoritative.** `duplicate_keys()` parses this
file's own AST and `main()` refuses to stay quiet about it — validated by injecting a known duplicate.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

OUT = pathlib.Path(__file__).resolve().parents[1] / "output" / "learned"

# ---------------------------------------------------------------- block descriptions
# key = (card, prefix); card "" applies to every card without a more specific entry.
BLOCKS: dict[tuple[str, str], str] = {
    # ---- design / context, annotated onto several cards by card_context.annotate() ----
    ("", "ctx_"): "Gene-design context from the decoy/atlas layer (`learned/analyses/gene_atlas.py` via "
                  "`card_context`): how wide and how measurable this gene's regulator design is, and the "
                  "real-vs-decoy gap measured on it. Gene-rung wherever it appears — on an edge or family "
                  "row it REPEATS by construction. ⛔⛔ **AND THE COST OF FORGETTING THAT IS LARGE — "
                  "MEASURED 2026-08-19, not asserted.** An EDGE-row statistic over any `ctx_` column "
                  "weights each gene by its DESIGN WIDTH (5,649 rows over 1,420 genes; up to **91 edges on "
                  "one gene**), and ceiling / gap / abundance all TRACK width — so the weighting is not "
                  "neutral, **it IS the confound**. Same column, two rungs, different answers:\n"
                  "        ctx_ceiling     edge-row median 0.0809  vs gene-row 0.0138   (~6x)\n"
                  "        ctx_n_abund                     4.00           vs      1.00   (4x)\n"
                  "        ctx_dose_max                   13.06           vs      8.97\n"
                  "        ctx_gap_deconv                 -0.0196         vs     -0.0085 (2.3x)\n"
                  "        ctx_gap_core                   -0.0186         vs     -0.0148\n"
                  "⇒ **compute every `ctx_` statistic on the GENE card, or de-duplicate to one row per "
                  "gene first.** Two of these are load-bearing: the decoy gap and the measurability "
                  "ceiling.",
    ("edge", "ctx_arm_dose"): "The ARM's dose in this gene's design context. ⚠ **ARM rung, not gene** — "
                              "unlike the rest of the `ctx_` block it is constant within ARM across genes "
                              "(verified 0/737 violations), so it repeats down the arm's rows instead of "
                              "the gene's.",
    ("edge", "ctx_arm_abundant"): "Is this ARM above the design's abundance threshold. ⚠ **ARM rung, not "
                                  "gene** — the other exception in this block (verified 0/737).",
    ("edge", "ctx_collapse"): "Mean arms collapsed per family in this gene's design (median 1.07) — how "
                              "much the seed-family collapse actually merged. Near 1 means the collapse "
                              "was almost a no-op for this gene.",
    ("edge", "ctx_d_fam"): "Design-WIDTH asymmetry between the real design and its matched decoy "
                           "(real minus fake families). ⚠ Median 1 — curation bundles same-seed mates and "
                           "fakes do not; MH-167 bounded the contamination this causes at <=0.0025 rho.",
    ("edge", "ctx_d_collin"): "Within-design COLLINEARITY asymmetry, real minus decoy. The second "
                              "design-match axis MH-167 tested; also bounded and not moving the gap.",
    ("edge", "ctx_gap_d_dose"): "Dose asymmetry between the real and decoy designs (median −0.43). ⚠ This "
                                "is the ONE design axis with a detectable effect on the gap, and it is "
                                "already corrected for in the headline number.",
    ("", "gctx_"): "Genomic context of the arm's precursor locus (`genomic_context.classify_he_arms`): "
                   "host gene, sense/antisense orientation, intron-vs-exon, promoter class. Strand-aware.",
    ("", "comp_"): "Cell-composition drivers of this gene's coupling (`compartment_*`): which deconvolved "
                   "fraction explains the correlation, its retention under adjustment, and its share. "
                   "⚠ This is DRIVER-SHARE, not the gene's own expression LOADING on a compartment — "
                   "`gene_axes` measured only the former to predict anything.",
    ("", "cov_"): "Coverage flags — was this row ever SCANNED by the named source. ⛔ Read before reading "
                  "the block it gates: a blank in a gated block means `— not scanned`, which is NOT a zero. "
                  "⭐ **THE FLAG→BLOCK MAP IS MEASURED, NEVER CURATED** (`gen_cards.cov_map`), requiring "
                  "`flag False ⟺ block entirely blank` in BOTH directions at ≥0.995. **All 16 bound "
                  "mappings hold at precision 1.000 / recall 1.000** — exact, not approximate. "
                  "⛔⛔ **AND THE NAMES LIE, WHICH IS WHY IT IS MEASURED (verified 2026-08-19).** A "
                  "hand-written table would bind by name and be WRONG for most of them: `cov_famrole` vs "
                  "`famrole_` recall **0.911** · `cov_attr` vs `attr_` **0.750** · `cov_ago_dom` vs `ago_` "
                  "**0.461** · `cov_comp` vs `cload_` **0.139** · `cov_meth` vs `bc_` **0.354** — all fail "
                  "the bar. `attr_` is in fact governed by **`cov_beta`** (1.000/1.000), not by `cov_attr`. "
                  "⚠⚠ **COVERAGE OF THE GUARANTEE IS PARTIAL: only 16 of the arm card's 42 blocks have a "
                  "governing flag, and 8 of the 21 flags govern nothing** (`cov_ago_dom`, `cov_attr`, "
                  "`cov_comp`, `cov_cptac`, `cov_expr`, `cov_famrole`, `cov_meth`, `cov_subtype`). ⇒ **for "
                  "the other 26 blocks — including `fame_`, `adm_`, `chim_`, and the new `fst_` and "
                  "`disc_` — a blank CANNOT be distinguished from a measured zero.** Do not read the "
                  "three-empty-states promise as covering the whole card.",

    # ---- the model ----
    ("", "beta_"): "Posterior summaries of the coupling coefficient from the dense Gibbs fit "
                   "(`learned/readouts.py`). `_sd` is the posterior SD, `_frac` the share of the gene's "
                   "total beta. ⚠ `beta_frac` is MAGNITUDE, not identity — it was renamed from `share` "
                   "precisely because it splits nothing (MH-140).",
    ("", "cal_"): "CALIBRATED posterior width (`learned/calibration.py`): the reported SD corrected for "
                  "the measured 1.18–1.34x narrowness, with `_lo`/`_hi` spanning the correction "
                  "constant's own +/-1 SD. Identifiability falls 24.8% -> 19.9% [18.1–22.4] under it (MH-185).",
    ("", "pip_"): "Posterior inclusion probability. `pip_dense` comes from the pi=1 coupling readout; "
                  "`pip_discovery` from the evidence-pi readout. ⛔ `pip_discovery` and `prior_pi` are "
                  "w-CONTAMINATED by construction — any 'do canonical regulators score higher' test on "
                  "them is circular.",
    ("", "oof_"): "Out-of-fold predictive correlation, arm-level vs family-level, and their difference "
                  "(`oof_drho` = arm minus family). Negative = the arm-resolved fit predicts repression "
                  "better. ⚠ Defined only in multi-arm family cells; elsewhere arm IS family.",
    ("", "aid_"): "Arm-IDENTIFIABILITY rollup for the arm: how often this arm is separable from its "
                  "same-seed mates (`arm_resolvable`, `arm_sep_z`, `oof_drho` medians).",

    # ---- attribution ----
    ("", "attr_"): "Attribution rollup over the arm's edges — how much of the gene's coupling this arm is "
                   "credited with, and over how many reliable edges.",
    ("", "identity_"): "Shapley/LMG attribution on R^2 with bagged-NNLS weights — WHO regulates the gene, "
                       "collinearity-fair. ⭐ Genuinely different from magnitude: they disagree in ~24% of "
                       "multi-family genes, and identity can honestly return UNDEFINED where NNLS zeroed "
                       "every family (`beta` structurally cannot).",
    ("", "lit_"): "⛔⛔ **A STUDY-COUNT AXIS, NOT A BIOLOGICAL GROUND TRUTH — read the name with "
                  "suspicion.** `eval/lit_ground_truth.py` defines the 'canonical' family of a gene as the "
                  "**argmax of DISTINCT PubMed IDs** carrying low-throughput functional evidence (reporter / "
                  "western / proteomics; miRTarBase-Weak and TarBase-Positive excluded). That measures **how "
                  "much a pair was PUBLISHED ON**, not how much it represses. Its virtue over what it replaced "
                  "is AUDITABILITY — mechanical, versioned, sha256-stamped, and it displaced five hand lists "
                  "with no producer in the repo — **not** that it is true. ⚠ It is the SAME KIND of quantity "
                  "as the `fame_` block, which the doctrine says to control for and never read as biology; "
                  "the difference is only that `fame_` is per-arm and this is an argmax per gene. ⚠ Defined "
                  "on 329 genes (21%) / 92 families, and the independent unit is the FAMILY, not the gene. "
                  "⚠ PMID depth tracks abundance, so this makes ABUNDANCE the baseline to beat and scaling n "
                  "cannot fix it. ⇒ 'agrees with lit' means 'picks the most-published family', and any "
                  "comparison against it needs a FAME NULL to carry meaning.",
    ("gene", "lit_family"): "The most-published seed family for this gene by distinct low-throughput "
                            "functional PMIDs. ⛔ A publication-depth argmax, not an established regulator.",
    ("gene", "lit_n_pmid"): "Distinct low-throughput functional PMIDs behind `lit_family`. Median **3** — the "
                            "'ground truth' rests on a handful of papers per gene.",
    ("gene", "lit_margin"): "PMID count of `lit_family` minus the runner-up. Median **1** — for half the "
                            "admitted genes the canonical family is decided by a SINGLE paper.",
    ("arm", "fame_assay_total_studies"): "`study_id.nunique()` — DISTINCT curated studies for this arm. "
                                         "⚠⚠ **NOT the denominator of the class marginals**: those are "
                                         "`has_<class>.sum()`, i.e. ROW counts, so a study spanning k "
                                         "target genes contributes k to a marginal and 1 here. The "
                                         "marginals therefore EXCEED this total on 25% of arms (max excess "
                                         "188). Do not compute class/total as a fraction.",
    ("arm", "fame_assay_unique_targets"): "`gene.nunique()` — DISTINCT curated TARGET GENES for this arm. "
                                          "⚠ It correlates with `fame_assay_total_studies` at r=0.999996 "
                                          "but is a DIFFERENT quantity, differing on **293 of 3,039 arms**; "
                                          "they track only because most (miRNA, gene) pairs carry exactly "
                                          "ONE study. ⛔ It was briefly pruned as 'redundant' on that "
                                          "correlation and RESTORED once the ingestion was read — a "
                                          "correlation cannot tell identity-by-definition from "
                                          "identity-in-this-pull.",
    ("arm", "fame_assay_binding_studies"): "`has_binding.sum()` — curated ROWS with binding evidence. "
                                           "⚠ Differs from its `binding__functional_mti_weak` cell on only "
                                           "12 of 3,039 arms (nearly all binding evidence is "
                                           "functional-weak), but the two are different quantities and can "
                                           "diverge as curation changes.",
    ("arm", "fame_npmid"): "Distinct PMIDs citing this arm. \u26d4 A FAME axis \u2014 control for it, never "
                           "read it as biology. (Its bit-identical twin `fame_led_n_pmid` was PRUNED "
                           "2026-08-19; this is the surviving name.)",
    ("arm", "arb_n_edges"): "Edges this arm regulates \u2014 equivalently the genes, since the source holds "
                            "one row per (arm, gene). (`arb_n_genes` was PRUNED 2026-08-19 as bit-identical "
                            "for that structural reason.)",
    ("arm", "arb_n_identity_reliable"): "Edges of this arm whose `identity` share is usable "
                                        "(`identity_reliable`). \u26d4 Until 2026-08-19 this was gated on "
                                        "`beta_frac_reliable`, which admits 100% of rows since MH-124 made "
                                        "\u03b2 strictly positive \u2014 so it equalled `arb_n_edges` and the "
                                        "gate did nothing. Now gated on identity's own coherence.",
    ("arm", "arb_max_identity"): "Largest usable `identity` share among this arm's edges. \u26d4 Read ONLY "
                                 "because it is gated: ungated it reached **+740.0** (identity is a SIGNED "
                                 "share summing to 1, 9.9% negative), now \u22120.255 \u2026 +1.000.",
    ("", "fame_"): "Study-attention axis for the arm — distinct PMIDs and distinct curated genes. "
                   "⛔ NOT a targetome measure: it correlates +0.75 with the curated degree and only +0.12 "
                   "with a sequence targetome (MH-208). Use it as a CONFOUND to control, never as biology.",
    ("", "cur_"): "CURATED degree — how many genes this arm is curated against (miRTarBase functional MTI). "
                  "⛔ A fame axis (see `fame_`), not a measure of how many genes the arm can actually hit.",
    ("", "seq_"): "SEQUENCE targetome degree — genes with a strong predicted duplex (scanMiR RBNS K_D), "
                  "breast-expressed. ⭐ This is the abundance-ORTHOGONAL promiscuity direction (rho=-0.004 "
                  "with dose); every other regulator axis is entangled with dose.",

    # ---- evidence / sites ----
    ("", "site_"): "Seed-site counts in the target 3'UTR by site type (8mer / 7mer / 6mer / non-canonical) "
                   "from the genome-wide scan. ⚠ A 6mer alone is ~24% chance noise at the median UTR — the "
                   "decoy layer gates them by Poisson, not presence.",
    ("", "ts_"): "TargetScan context++ scores for the arm's sites — number of sites, mean and best context "
                 "score. More negative context++ = stronger predicted repression.",
    ("", "sdsz_"): "MANE 3'UTR seed scan. ⚠ A FOURTH targetome universe — never blend it with the curated, "
                   "sequence-K_D or TargetScan universes; they get four separate axes by policy.",
    ("", "adm_"): "ADMISSIBILITY — could this edge repress in breast AT ALL? `has_site` (seed match in the "
                  "3'UTR) AND `expressed` (arm present in TCGA tumour OR GTEx healthy breast). ⭐ The "
                  "project's most load-bearing conditioning variable, and it was missing from every card "
                  "until MH-216. ⚠ Keep it a diagnostic FLAG — as a filter it did not sharpen the decoy gap.",
    ("", "echim_"): "CHIMERIC ligation evidence at its NATIVE (gene, arm) rung — a physical duplex read, "
                    "per source, never pooled. ⭐ Chimeric is the one evidence type that RESOLVES the mature "
                    "arm. ⚠ No breast chimeric source exists (Manakov = HEK293T), so this is cross-tissue "
                    "corroboration, never breast-context validation.",
    ("", "chim_"): "Chimeric-evidence ROLLUP for the arm (aggregated over its edges). ⚠ Distinct from "
                   "`echim_`, which is the per-edge form — they were one prefix until MH-220 and the "
                   "collision destroyed a 5-column block.",
    ("", "kd_"): "scanMiR RBNS binding affinity. \u2b50 `kd_affinity_pct` answers *is this gene among "
                 "THIS ARM's strongest targets* \u2014 a WITHIN-ARM percentile over the arm's genome-wide "
                 "K_D targetome, so **HIGHER = a stronger target for this arm**. \u26a0 Being a within-arm "
                 "rank it says nothing about ABSOLUTE affinity: a promiscuous arm's 99th percentile and a "
                 "specialist's 99th percentile are not comparable in K_D \u2014 `kd_repression` (the raw "
                 "scanMiR value) ships beside it for exactly that reason. Sequence-only, so a-priori and "
                 "non-circular with respect to expression.",

    # ---- expression / abundance ----
    ("", "abund_"): "The arm's own abundance profile across the cohort — median, SD, IQR, detection "
                    "fraction. ⭐ The DYNAMIC RANGE (sd/iqr) is the informative half; the MEAN predicts "
                    "much less (measured in MH-201).",
    ("", "tier_"): "Expression tier of the arm against the detection floor — class, fraction of samples "
                   "expressed, p90 and max log2 RPM.",
    ("", "hly_"): "HEALTHY baseline levels for the arm — GTEx breast median/percentile and TCGA-NAT median.",
    ("", "hx_"): "⛔ IMPUTED healthy baseline (miTED rank-transfer or seed-mate substitution) — LABELLED as "
                 "imputed, never presented as a measurement. Read `hx_hly_baseline_imputed` before using "
                 "any healthy-leg quantity.",
    ("", "dose_"): "Dose (abundance) of the regulator and how it shifts across states, plus host-gene "
                   "co-transcription fractions where the arm is intragenic.",

    # ---- cross-state / progression ----
    ("", "share_"): "The arm's share of the gene's total regulatory dose in one state (HLY / NAT / TUM). "
                    "Within-gene, sums to 1 across the gene's arms.",
    ("", "rank_"): "The arm's WITHIN-GENE rank by dose in one state. ⚠ Distinct from `grank_`, which is "
                   "the global (cohort-wide) rank — a within-gene rank cannot see cohort-level shifts.",
    ("", "grank_"): "The arm's GLOBAL (cohort-wide) dose rank in one state. This is what the dominant-"
                    "regulator and handoff columns are computed from.",
    ("", "field_"): "Within-patient NAT field-effect fit — is the regulator already established in the "
                    "patient's adjacent normal before the tumour. Defined on the paired subset.",
    ("", "wshift_"): "Within-patient paired quantities on the ~103 matched tumour/NAT participants. "
                     "⚠ P-ACROSS rung: one value per arm computed ACROSS patients, not per patient.",
    ("", "real_"): "Realization rollup for the ARM — how many of its edges are scored, and its median and "
                   "best realized coupling.",
    ("", "armres_"): "⭐ IS THIS GENE MODELABLE AT ARM RESOLUTION — the gene-level roll-up of arm "
                   "identifiability (`card_ladders.gene_arm_resolution`, 2026-08-19). ⚠ THE DENOMINATOR IS "
                   "THE POINT: only the **273 of 1,549 genes (17.6%)** with at least one MULTI-ARM cell can "
                   "be asked; elsewhere the arm rung IS the family rung by construction, so `armres_class` is "
                   "`not_applicable` — a design fact, not a negative. ⚠ Strongly BIMODAL among those 273 — "
                   "ALL cells resolvable in 68 genes, NONE in 168 — so never read its mean.",
    ("gene", "armres_class"): "`arm_modelable` (every multi-arm cell resolvable, 68 genes) · `partial` (37) · "
                            "`family_only` (no cell resolvable, 168). Absent = the gene has no multi-arm "
                            "cell at all, so the question is undefined.",
    ("gene", "armres_med_drho"): "Median `oof_drho` over the gene's multi-arm cells — negative means arm "
                               "resolution predicts repression better out of fold. Median across genes "
                               "−0.0012, best −0.206: the gain is a TAIL, concentrated in the miR-29 family.",
    ("", "fst_"): "⭐ THE ARM'S SHARE OF ITS OWN SEED FAMILY, ACROSS STATES — gene-free "
                  "(`card_ladders.family_state_shift`, 2026-08-19). Fills a real gap: every other "
                  "cross-state share is anchored to a GENE (`share_*`, `rank_*`, `d_rank_*`) or to the whole "
                  "cohort (`grank_*`, `dGlobal_*`), while the arm's standing inside its own family existed "
                  "only as a single-state snapshot (`famrole_*`). ⭐ Unlike the gene-level "
                  "`regulatory_handoff` it CANNOT be a design-width artifact — the denominator is the "
                  "family, which is fixed. ⚠ MASK `fst_n_members == 1` first: 66.9% of arms are singletons "
                  "where share ≡ 1 by construction. ⛔⛔ USE THE NAT LEG. The `_HLY` leg crosses GTEx→TCGA "
                  "and carries a platform artifact — near-total flips (|Δshare|>0.90) run 7.1%, and "
                  "requiring every sibling to be measured makes it WORSE (11.6%), versus **1.4% on the "
                  "same-platform NAT→TUM leg**. ✅ Defensible: `fst_d_share_NAT_TUM` median +0.032, "
                  "|Δ|>0.10 for 43.8%, dominant member changes in **26.7%** (56/210). DESCRIPTIVE — a "
                  "registry candidate, no decoy and no null yet.",
    ("arm", "cov_fst"): "Was the family-share block computed for this arm at all — i.e. did it have the "
                        "abundance input. True for **485 of 2,450** arms; the rest were NEVER COMPUTED, "
                        "which is genuinely UNSCANNED and not a zero share. ⚠⚠ **One flag CANNOT govern "
                        "the whole `fst_` block, because the block has THREE coverage strata**: TUM-only, "
                        "TUM+NAT, and TUM+NAT+HLY (the healthy leg reaches only ~17% of arms). So "
                        "`gen_cards.cov_map` correctly refuses to bind it — the map assumes a block is "
                        "homogeneous, and this one is not. Read `cov_fst` for *was it computed*, then "
                        "`fst_hly_measured` for *does it have a healthy leg*. Two flags, two questions.",
    ("arm", "fst_hly_family_measured"): "Do ALL members of this arm's family have a measured GTEx level. "
                                        "Emitted because it is INFORMATIVE, not because it repairs the "
                                        "healthy leg: the stricter guard RAISES the near-total-flip rate "
                                        "(7.1% → 11.6%), which is what localised the defect to the "
                                        "GTEx→TCGA platform boundary rather than to missing members.",
    ("", "greal_"): "Realization ladder for the GENE — the mirror of `real_`: how many of THIS gene's "
                    "regulators realize at each calibrated-coupling depth (c05/c10/c15/c20/c30). "
                    "✅ Cross-rung control: the gene-side totals reproduce the arm-side exactly.",
    ("", "acquired_"): "Acquisition of regulatory dose relative to a healthy reference (GTEx or NAT). "
                       "⚠ The GTEx leg passed through a multi-mapping collapse that FABRICATED acquisition "
                       "for 115 arms until MH-210; read the `_meas` (measured-only) variants beside it.",

    # ---- protein ----
    ("", "cptac_"): "CPTAC protein channel. `prosp` = the PROSPECTIVE cohort (independent patients, TMT-11, "
                    "n~101); `t105` = TCGA-105 (the SAME patients as the mRNA fit, so layer-only, and "
                    "circular for anything the mRNA fit used). `rho_rna`/`rho_prot` are coupling against "
                    "each layer; `rho_disc` is protein-BEYOND-mRNA (protein residualised on its own mRNA); "
                    "`ret_` is the fraction surviving composition adjustment. ⭐ `rho_disc` is flat "
                    "everywhere => the action is transcript DESTABILISATION, not translational repression.",
    ("", "tcga_"): "The TCGA mRNA leg carried beside the CPTAC columns so the two are comparable on the "
                   "same gene rows.",

    # ---- subtype / clinical ----
    ("", "esub_"): "PAM50 subtype-stratified coupling at the EDGE rung — one fit per subtype plus the "
                   "pooled fit and an adjusted q. ⚠ Named `esub_` and NOT `sub_` deliberately: the arm "
                   "card owns `sub_` for its per-arm AGGREGATE of the same source (MH-213/220). "
                   "⚠ Subtype heterogeneity is DISTRIBUTIONAL, not a per-edge label.",
    ("", "sub_"): "PAM50 subtype block aggregated to the ARM. See `esub_` for the per-edge form.",
    ("", "surv_"): "TCGA outcome panel legs. ⛔ TCGA only — the GSE/Buffa leg is bare-stem matched and is "
                   "not in this block. Pressure magnitude is exhaustively non-prognostic (0 FDR hits).",

    # ---- copy number ----
    ("", "cnv_"): "Copy number at the arm's own locus — mean CN and amplification/deletion summaries.",
    ("", "cnvc_"): "CN-to-expression concordance fit for the arm (756 arms): does locus CN actually move "
                   "this arm's abundance. ⛔ DESCRIPTIVE — both CN instruments are retracted, so this is "
                   "not causal evidence.",
    ("", "foc_"): "Paralog-decontaminated FOCAL-locus concordance. ⛔ DESCRIPTIVE, explicitly not an "
                  "instrument.",

    # ---- isomiR / family structure ----
    ("", "iso_"): "isomiR profile for the arm — 5' shift fraction and seed-shift occurrence. Drives the "
                  "isomiR-corrected X_fam, which is BUILT but default-OFF (coupling wash, bidirectional).",
    ("", "isoc_"): "isomiR COMPOSITION — what fraction of the arm's reads are canonical, orphan, or "
                   "donated from a neighbouring arm.",
    ("", "fam_"): "Seed-family structure the arm sits in — membership, affinity and dominance within it. "
                  "⭐ In a typical multi-member family ONE member carries ~62% of the dose, which is the "
                  "structural fact that justifies the arm/family split.",
    ("seed_family", "fam_single_member"): "⭐⭐ **THE DOMINANT FACT ABOUT THIS ENTIRE CARD: 1,639 of "
                                          "1,959 families (83.7%) are SINGLETONS.** For those rows every "
                                          "'family-level' statistic is TRIVIALLY that one arm's value, and "
                                          "`fam_dominant_share` is forced to exactly 1.0. ⛔ **Mask on this "
                                          "before ANY family-level comparison** — pooling the singletons is "
                                          "the degeneracy trap (`gene_axes.mask_degenerate`), and it is "
                                          "worse here than anywhere else in the card system.",
    ("seed_family", "fam_n_members"): "Arms in the seed family (median 1, max 22). ⚠ Only **320 families "
                                      "have more than one**. This is the family's TOTAL membership — not "
                                      "how many appear in any one gene's design, which is "
                                      "`n_arm_in_cell` and is a median 0.67 fraction of this (MH-248).",
    ("seed_family", "fam_dominant_share"): "The top member's share of the family's total dose. ⛔ Median "
                                           "1.0000 over ALL families is meaningless — singletons force it. "
                                           "⭐ **On the 320 MULTI-member families the median is 0.6326, "
                                           "independently reproducing MH-215's recorded 0.619**, and one "
                                           "member carries >80% of the dose in **28.1%** of them. That "
                                           "concentration is the structural fact justifying the arm/family "
                                           "split.",
    ("", "fam_dose_hhi"): "Concentration of DOSE across the family's members. ✅ **FLOOR-CORRECTED — "
                          "verified 2026-08-19**: its minimum sits far below the raw 1/k floor at every "
                          "member count (0.0000 at k=2 where the raw floor is 0.500), and "
                          "spearman(n_members, hhi) = **+0.1896** rather than the strongly negative value "
                          "a raw index would give. ⚠⚠ **CONTRAST `concentration` ON THE GENE CARD, which is "
                          "NOT corrected** and reads −0.9399 against design width. Same trap, two different "
                          "handlings in one codebase — check which form you have before comparing them.",
    ("", "fam_targetome_"): "Family-level targetome size. ⛔⛔ **TWO DIFFERENT UNIVERSES, NEVER TO BE "
                            "BLENDED** (the doctrine gives the four targetome universes four separate "
                            "axes): `_seed_med` is the MANE 3′UTR seed scan — fill **96.9%**, median "
                            "**2,040**; `_seq_med` is the scanMiR K_D sequence targetome — fill **30.5%**, "
                            "median **3,547**. They correlate at spearman 0.828 where both exist but sit "
                            "**1.78× apart in scale** and differ 3× in coverage, so a mixed statistic over "
                            "them is meaningless.",
    ("", "fam_redund_"): "Within-family redundancy — effective member count, over-count fraction, and the "
                         "median within-family correlation. ⚠ **Fill 4.3% (84 families)** — it requires a "
                         "correlation fit across members, so it exists only for well-measured multi-member "
                         "families. Read it as a case series, not a distribution.",
    ("", "fam_ctx_"): "The family's GENOMIC context, pooled over its members: host-coding fraction, "
                      "co-transcribed fraction, distinct coding hosts, and whether the members share a "
                      "context. Fill 29.1%. ⚠ `fam_ctx_context_homogeneous` False for only 47 of 571 — most "
                      "families that can be assessed are context-homogeneous, so this rarely discriminates.",
    ("", "famrole_"): "The arm's ROLE inside its seed family — member count, its share of family abundance, "
                      "and its rank among members.",
    ("", "family_"): "Family-level properties carried onto this row (size, the arm's dose share of it, and "
                     "its role). ⚠ On the gene_family card the family IS the grain; on the edge card these "
                     "REPEAT across the family's arms by construction.",

    # ---- AGO / machinery ----
    ("", "ago_"): "AGO/RISC loading gate for the arm. ⛔ `ago_loading` ends in `.fillna(1.0)`, so a MEASURED "
                  "guide (1.0) is indistinguishable from an arm never assayed (also 1.0) — 32.4% of the "
                  "column is that default. Use `ago_loading_measured` to tell them apart. The axis is "
                  "coupling-INERT: its value is identity/QC, not coupling.",

    # ---- co-expression / methylation / misc ----
    ("", "comv_"): "Co-expression module membership and top co-expressed partner. Pure co-expression — "
                   "the cross-state class columns are deliberately excluded from this block.",
    ("", "bc_"): "Breast-cancer locus annotation for the arm — CpG probes, CGI overlap and related "
                 "epigenetic context at the precursor locus.",
    ("", "arb_"): "Per-ARM rollup of its edge set — how many edges and genes it regulates and its mean "
                  "|beta| over them. ⚠ An aggregate: it inherits the edge rung's caveats.",
    # ── UNIT 12 (2026-08-19) — ⛔ THE `arm_` PREFIX IS AMBIGUOUS *WITHIN* THE EDGE CARD.
    # It carries BOTH the arm-resolution FIT quantities (arm_sep_z / arm_dbeta / arm_resolvable /
    # arm_credit_share / arm_id_status) AND seven lifted ABUNDANCE-and-TRAJECTORY columns. A prefix block
    # cannot serve both, so the seven get exact entries here — which beat any prefix, on every card.
    ("", "arm_med_rpm"): "The arm's MEDIAN RPM across tumour samples — its abundance level, gene-free. "
                         "Defined on the 19.8% of arms with tumour expression data.",
    ("", "arm_pct_floor"): "⛔ **THE NAME READS BACKWARDS.** Computed as `100·mean(x > FLOOR)` — the percent "
                           "of samples ABOVE the detection floor, i.e. **HIGH = WELL DETECTED**, not "
                           "'percent at the floor'. Verified: ρ=**+0.917** with `arm_med_rpm`. "
                           "Median arm sits at **5%**; 162 arms at 0%, 91 at 100%. Feeds `arm_spiker` "
                           "(`pct_floor < 40 AND arm_iqr > 1.5`).",
    ("", "arm_iqr"): "Interquartile range of the arm's tumour RPM — its DYNAMIC RANGE. ⭐ Axiom 8: "
                     "dispersion predicts where coupling is estimable; the MEAN level does not. Pairs with "
                     "`arm_pct_floor` to define `arm_spiker`.",
    ("", "arm_lfc_NAT_TUM"): "log-FC of the arm's median abundance, matched NAT → tumour. ⭐ **THE ONLY "
                             "GAUGE-CLEAN LEG** — same platform, same patients; median **+0.080**, i.e. "
                             "centred on zero as a fold-change should be. This is the trajectory column "
                             "`card.py` uses for the dose axis, and the one to prefer.",
    ("", "arm_lfc_HLY_TUM_QN"): "log-FC GTEx-healthy → TCGA-tumour, on quantile-normalised values. "
                                "⚠ **CROSS-PLATFORM: median +0.910, NOT zero.** The offset is a gauge "
                                "constant `c` (memory `cross-cohort-gauge`), not biology. ⇒ **its LEVEL is "
                                "uninterpretable; only the CONTRAST between arms is.**",
    ("", "arm_lfc_HLY_TUM_raw"): "As `arm_lfc_HLY_TUM_QN` but against the RAW GTEx median. ⛔⛔ **QN DOES "
                                 "NOT FIX THE OFFSET — measured 2026-08-19: raw median +1.060 vs QN +0.910, "
                                 "so quantile normalisation removes only 0.14 of a ~1.0 gauge shift.** And "
                                 "the two rank-agree at only **ρ=+0.632** on the same 335 arms, so they are "
                                 "not interchangeable either. Treat BOTH healthy legs as offset.",
    ("", "arm_lfc_HLY_NAT_raw"): "log-FC GTEx-healthy → TCGA-NAT, raw. ⚠ Same ~+1 cross-platform gauge "
                                 "offset (median **+0.970**) — and since NAT is *not* tumour, a nonzero "
                                 "median here is the offset's own fingerprint: two tissues that should "
                                 "differ least show the same +1 as the tumour leg.",
    # ── UNIT 13 — the two `site_frac_*` columns have DIFFERENT DENOMINATORS. Nothing in the block said so.
    ("", "site_frac_8mer"): "8mer share of the arm's **CANONICAL** sites (`site_n_8mer / site_n_canonical`; "
                            "median 0.150). ⛔ **NOT the same denominator as `site_frac_canonical`** — "
                            "these two neighbours divide by different things, which is why 8mer share "
                            "(0.150) reads *larger* than canonical share (0.054). ⚠ And it is TOOL-"
                            "DEPENDENT: agrees with `ts_frac_8mer` at only ρ=+0.620.",
    ("", "site_frac_canonical"): "Canonical share of the arm's **TOTAL** sites (`site_n_canonical / "
                                 "site_n_total`; median 0.054) — a different denominator from "
                                 "`site_frac_8mer`. Low because non-canonical sites outnumber canonical "
                                 "~12:1 (medians 5,388 vs 439).",
    ("", "site_sites_per_gene"): "Canonical sites per targeted gene (median 1.34) — site DENSITY, "
                                 "separating an arm with many sites on few genes from one spread thin. "
                                 "`site_n_canonical / site_n_genes_canonical`.",
    ("", "site_repression_med"): "Median predicted repression score over the arm's sites (all negative by "
                                 "construction; median −0.316). A PREDICTION from sequence, not a "
                                 "measurement — never cite it as evidence of observed repression.",
    ("", "site_repression_min"): "Strongest (most negative) predicted repression among the arm's sites "
                                 "(median −1.964). ⚠ A MIN over a set whose size varies 450→40,511 — it "
                                 "partly measures how many sites were drawn, not how strong the best is.",
    ("", "arm_"): "Arm-resolved fit quantities inside a family cell — the arm's own beta, its separation "
                  "z from its same-seed mates, and whether it is `arm_resolvable`. Defined only where a "
                  "family cell holds more than one arm (20.3% of edges).",
    ("", "term_"): "Decomposition of the acquisition score into ABUNDANCE, WIRING and INTERACTION terms.",
    ("", "coupling_"): "Partial-Spearman coupling of arm dose against target mRNA in one state, "
                       "conditioned on C (composition + target CN + proliferation), with its p and z. "
                       "⚠ The per-edge null is 3–4× too narrow — read these as a DISTRIBUTION, not as "
                       "per-edge significance. ⭐ **THE HONEST SCALE is `coupling_z_tum`, the effect in "
                       "NULL SDs: median −0.32, with |z|>2 on only 15.5% of edges** (`coupling_p_tum` "
                       "<0.05 on 16.2%). ⛔⛔ **`coupling_tum` IS MISSING ON 27.9% OF EDGES, AND THAT IS "
                       "STALENESS, NOT A PROPERTY OF THE DATA — diagnosed 2026-08-19.** Of 25 sampled "
                       "missing edges, **24 produce a finite coupling when recomputed now** (including "
                       "**PTEN→miR-19b-3p = −0.3436**, a canonical edge simply absent). The cause is that "
                       "`edge_card_base.tsv` is dated **2026-08-04** and has not been rebuilt, while its "
                       "inputs moved on — so the delivered card layers fresh annotations over a 15-day-old "
                       "base. ⛔ **I first reported the missingness as 'not at random, concentrated on "
                       "wide/well-measured genes and seedless edges' and drew a measurability moral from "
                       "it. THAT IS RETRACTED**: those correlates describe what CHANGED between the base "
                       "build and now, not a measurement property. ⚠ Four mechanisms were tested and "
                       "refuted first (arm absent from the design, zero variance, whole-gene failure, NaN "
                       "propagation) — a structured missingness is a PROVENANCE question before it is a "
                       "data question. ⇒ **treat every `coupling_*`, `share_*`, `rank_*`, `shift_class` "
                       "and `arm_lfc_*` value as 2026-08-04 vintage until `canonical_card.build()` is "
                       "re-run.** ⚠ `coupling_fam` (99.9% complete, spearman +0.901) is a COMPLETENESS "
                       "control only — it comes from the SAME stale base and controls nothing for freshness.",
    ("edge", "coupling_fam"): "Family-level coupling for the edge's (gene, seed_family) cell. ⭐ **99.9% "
                              "complete against `coupling_tum`'s 72.1%, and spearman +0.901 where both "
                              "exist** — so it is the coverage-robust stand-in whenever `coupling_tum`'s "
                              "non-random missingness could bias a comparison. 1,570 edges have this and "
                              "no arm-level value.",
    ("edge", "coupling_z_tum"): "The tumour coupling expressed in NULL SDs — the CALIBRATED effect size, "
                                "and the scale to quote. Median **−0.32**; |z|>2 on 15.5%, |z|>3 on 6.4%. "
                                "⚠ Use this rather than the raw ρ when comparing across states, since the "
                                "null's width differs between them.",
    ("edge", "coupling_p_tum"): "p against the CALIBRATED site-free null (not a theoretical t-null). "
                                "Median 0.375; <0.05 on 16.2%. ⛔ Do NOT read a per-edge FDR off this — "
                                "site-free pairs that CANNOT repress pass at 25–35% under the older "
                                "uncalibrated null, which is why the defensible claim is a set-level "
                                "distributional shift.",
    ("edge", "coupling_hly_surrogate"): "Healthy coupling via a same-seed SURROGATE arm, used where the "
                                        "arm itself has no GTEx signal. ⚠ Fill only **3.4%** — it rescues "
                                        "few edges. `coupling_p_hly_resolved` is the direct|surrogate "
                                        "resolution and reaches 59.3%.",
    ("", "retention_"): "Fraction of an effect surviving covariate adjustment (adjusted / raw). "
                        "⚠⚠ 'Retention' names TWO unrelated quantities in this project — on the CARDS it "
                        "always means this one, NOT the patient-baseline NAT->tumour persistence.",
    ("", "realized_"): "Realized (as opposed to potential) coupling for the gene — the raw and adjusted "
                       "rho, its retention, and how many regulators it was computed over.",
    ("", "realization_"): "Within-patient realization quantities — the score, its identity allocation, and "
                          "who owns it.",
    ("", "dominant_"): "The top-ranked regulator of this gene in one state (HLY / NAT / TUM), by GLOBAL "
                       "dose rank. The inputs to `regulatory_handoff`.",
    ("", "top_"): "The gene's leading regulator and its value, by magnitude or by identity. "
                  "⚠ `top_family_magnitude` and `top_family_identity` are DIFFERENT questions and disagree "
                  "in ~24% of multi-family genes.",
    ("", "prior_"): "The evidence prior fed to the discovery readout. ⛔ w-contaminated by construction.",
    ("", "model_"): "Which model variant produced this row.",
    ("", "net_"): "Net (signed) pressure summed over the gene's regulators.",
    ("", "frac_"): "A fraction — read the suffix for the numerator and `domain` for the denominator's rows.",
    ("", "mean_"): "A per-row mean of a lower-rung quantity; check `agg_of` for which rung it summarises.",
    ("", "median_"): "A per-row median of a lower-rung quantity; check `agg_of`.",
    ("", "composition_"): "Composition-adjustment verdict for this row (e.g. `composition_explained` when "
                          "the effect does not survive the deconvolution block).",
    ("", "cload_"): "Compartment LOADING of the arm. ⚠ Loading is NOT driver-share — see `comp_`. "
                    "All 11 loading axes were null under FDR where driver-share separated the harm tails.",
    ("", "seed_"): "Seed sequence and the seed family it defines.",
    ("", "gene_"): "Gene-level properties carried onto this row (its regulator count, repression class, "
                   "net-repression flag). On the edge card these REPEAT within gene by construction.",
    ("", "surrogate_"): "Surrogate healthy instrument used where a direct healthy leg is unavailable, and "
                        "its correlation with the quantity it stands in for.",
    ("", "healthy_"): "Healthy-state legs and flags — including `healthy_uninformative`, which marks rows "
                      "where the healthy baseline cannot carry a claim.",
    ("", "wiring_"): "Fraction of the acquisition score attributable to WIRING rather than abundance.",
    ("", "acquisition_"): "Composite acquisition score across states.",
    ("", "shift_"): "Cross-state shift class for the edge (acquired / lost / constitutive / uncoupled / ...).",
    ("", "own_"): "Within-patient own-NAT referenced quantity. ⚠ P-across rung.",
    ("", "is_"): "A boolean flag; read the suffix.",
    ("", "static_"): "The static (dose-based) owner of the gene, for comparison against the realization owner.",
    ("", "owner_"): "Whether the realization owner and the static owner agree. Defined only where both exist.",
    ("", "regulatory_"): "`regulatory_handoff`: does the gene's globally-dominant regulator DIFFER between "
                         "healthy and tumour, OR between NAT and tumour. ⚠ It fires on EITHER leg — 27 of "
                         "328 genes fire only via the NAT leg and are NOT healthy->tumour switches.",
    ("", "delta_"): "A tumour-minus-normal difference for this row.",
    ("", "pressure_"): "Aggregate incoming pressure on the gene in the named state.",
    ("", "total_"): "Sum over the gene's regulators.",
    ("", "max_"): "Maximum over the row's constituent units.",
    ("", "sd_"): "Standard deviation of the named quantity.",
    ("", "rho_"): "A Spearman correlation; read the suffix for what against what.",
    ("", "edge_"): "An edge-rung quantity carried onto this row.",
    ("", "spiker"): "Flag: the arm's abundance is carried by a few spiking samples rather than a broad "
                    "distribution \u2014 a measurement-reliability warning, not biology.",
    ("", "detection"): "Detection status of the arm against the expression floor.",
    ("", "dGlobal_"): "Change in the arm's GLOBAL (cohort-wide) dose rank between two states. "
                      "HLY_NAT is the FIELD delta, NAT_TUM the malignant step.",
    ("", "d_rank_"): "Change in the arm's WITHIN-GENE dose rank between two states. \u26a0 Distinct from "
                     "`dGlobal_`: a within-gene rank cannot see a cohort-wide shift.",
    ("", "dShare_"): "Within-patient change in the arm's share of the gene's dose \u2014 `_M_own` referenced "
                     "to the model's M, `_raw_own` to the raw budget. P-across rung.",
}

# ---------------------------------------------------------------- per-column overrides
COLUMNS: dict[tuple[str, str], str] = {
    ("edge", "beta"): "EDGE-rung coupling coefficient — `readouts.run(level='arm')`, the same-seed collapse "
                      "REMOVED. Posterior mean of the dense Gibbs fit. Never exactly 0: the half-normal "
                      "slab has a strictly positive mean and cannot zero an un-informed regulator."
                          " ⭐⭐ **AND THE TWO CARDS' `beta` DIFFER EVEN WHERE THE COLLAPSE IS A NO-OP.** Measured 2026-08-19: on **SINGLETON** cells — one arm, nothing to collapse — `beta_edge == beta_family` on only **21.6%** (median |diff| 0.0004, max 0.2444); on multi-arm cells, **0.2%** (median 0.0100, max 0.9680). The family fit re-designs the WHOLE gene, so collapsing a gene's OTHER families changes this one's competitors and moves its beta. ⇒ **never mix the two cards' beta, not even for singleton families.** ✅ Edge-rung verified: beta varies across arms in **467 of 467** multi-arm cells.",
    ("gene_family", "beta"): "FAMILY-rung coupling coefficient — `readouts.run(level='family')`, the "
                             "same-seed collapse APPLIED. Same estimator as the edge card's `beta`, "
                             "DIFFERENT UNIT. Do not place the two side by side without saying so."
                          " ⭐⭐ **AND THE TWO CARDS' `beta` DIFFER EVEN WHERE THE COLLAPSE IS A NO-OP.** Measured 2026-08-19: on **SINGLETON** cells — one arm, nothing to collapse — `beta_edge == beta_family` on only **21.6%** (median |diff| 0.0004, max 0.2444); on multi-arm cells, **0.2%** (median 0.0100, max 0.9680). The family fit re-designs the WHOLE gene, so collapsing a gene's OTHER families changes this one's competitors and moves its beta. ⇒ **never mix the two cards' beta, not even for singleton families.** ✅ Edge-rung verified: beta varies across arms in **467 of 467** multi-arm cells.",
    ("edge", "identified"): "`|z| > 2` on the UNCALIBRATED posterior SD, where `z = beta / beta_sd` "
                            "(MH-83 records precision 0.86 / recall 0.89 against a held-out standard). "
                            "⚠ The posterior SD is measured 1.18–1.34× too NARROW, so this over-admits: "
                            "`cal_identified` is the honest version (19.9% [18.1–22.4] vs 24.8%). "
                            "⚠ The `abs()` is vestigial — the half-normal slab makes β strictly positive, "
                            "so `z` never goes below 0 (observed range 0.47…22.33).",
    ("edge", "z"): "`beta / beta_sd` for the EDGE-level fit (`readouts.run(level='arm')`, the whole gene's "
                   "arms as separate predictors). ⚠⚠ NOT the same as `z_arm`, which is a DIFFERENT fit — "
                   "restricted to one (gene, seed_family) cell (`arm_rung.py`). They correlate 0.971 on the "
                   "1,147 edges carrying both, but they answer different questions: `z` is 'is this arm "
                   "identified against the gene's whole design', `z_arm` is 'is it separable from its "
                   "same-seed mates'. Both are strictly positive (β cannot be negative).",
    ("edge", "arm_med_rpm"): "Median linear RPM of this ARM across the cohort. Arm-rung, so it REPEATS "
                             "across the arm's genes. Fill 0.686 — the gap is arms with no expression "
                             "record, not zeros.",
    ("edge", "arm_iqr"): "Interquartile range of the arm's abundance — its DYNAMIC RANGE. ⭐ Dispersion is "
                         "usually the informative half against level (axiom 8), though it did NOT separate "
                         "diluting from improving mates in MH-248 (q=0.47).",
    ("edge", "arm_pct_floor"): "Percent of samples in which the arm sits above the detection floor "
                               "(median 99). Feeds `arm_spiker` together with `arm_iqr`.",
    ("", "arm_lfc_"): "Log fold-change of the ARM's abundance between two states. ⚠⚠ **ONLY "
                      "`arm_lfc_NAT_TUM` IS SAME-PLATFORM (TCGA→TCGA). Every `HLY` leg crosses "
                      "GTEx→TCGA**, and that boundary carries a measured artifact — in the `fst_` work "
                      "near-total share flips were 5–8× more common across it than within TCGA. ⛔⛔ **AND "
                      "THE TWO HLY→TUM VARIANTS DISAGREE MATERIALLY: `_QN` vs `_raw` correlate at only "
                      "spearman +0.544, median |diff| 1.49, and they DISAGREE ON SIGN FOR 23.9% OF "
                      "EDGES.** Which leg you pick changes the direction for ~1 edge in 4, so state which "
                      "one you used. ⚠ The QN bridge is CORRECT and settled but carries a standing "
                      "trust-weighting — prefer rank measures over QN magnitudes for cross-platform work.",
    ("edge", "hly_leg_concordant"): "Do the two HLY→TUM legs (`_QN` and `_raw`) at least AGREE ON SIGN. "
                                    "**Concordant 76.1%, DISCORDANT 23.9%** of the 3,433 edges carrying "
                                    "both. ⛔ **NOT a correctness flag** — agreement does not make either "
                                    "leg right; it is a DISAGREEMENT flag, which is what the data "
                                    "supports. Gate any cross-state claim on it rather than resting on an "
                                    "arbitrary choice of leg, and remember both legs cross GTEx→TCGA "
                                    "while only `arm_lfc_NAT_TUM` is same-platform.",
    ("edge", "arm_credit_share"): "This arm's share of its (gene, seed_family) cell's credit. ✅ Verified "
                                  "to sum to 1.0000 within every one of the 467 multi-arm cells, and it "
                                  "genuinely varies PER ARM inside the cell (467/467).",
    ("edge", "arm_id_status"): "How the arm's identity within its cell was resolved: `resolved_by_dose` "
                               "2,344 · `singleton` 2,260 · `unvalidated_balanced` 1,045. ⚠ **Fill is 1.0 "
                               "while its arm-in-family siblings are 0.203** — because `singleton` is a "
                               "VALID answer for a one-arm cell (95.8% of singletons are), where the "
                               "siblings are NaN. Same rung, different domain. ⛔ `unvalidated_balanced` "
                               "means dose CANNOT resolve which arm acts — it needs a breast chimeric-eCLIP "
                               "or a per-arm knockdown, not more modelling.",
    ("edge", "arm_dbeta"): "max−min β across the cell's member arms. ⚠⚠ **A CELL quantity REPEATED across "
                           "the cell's arms — verified constant within all 467 multi-arm cells.** Reading "
                           "it as a per-arm fact double-counts; the per-arm quantities are `beta_arm` / "
                           "`z_arm` / `arm_credit_share`, which do vary (467/467).",
    ("edge", "arm_sep_z"): "`arm_dbeta` in units of the cell's own pooled posterior SD — can these "
                           "same-seed arms be told apart at all. ⚠⚠ **Also a CELL quantity, constant "
                           "within all 467 cells.** Median 1.59.",
    ("edge", "arm_resolvable"): "Are the cell's arms separable AND does splitting them help out-of-fold. "
                                "True for 31.7% of edges in multi-arm cells. ⚠⚠ **A CELL verdict, constant "
                                "within all 467 cells** — not a property of the individual arm.",
    ("edge", "z_arm"): "`beta_arm / sd_arm` from the WITHIN-CELL fit (`arm_rung.py`) — the arm's z inside "
                       "its (gene, seed_family) cell only, defined on the 20.3% of edges in multi-arm cells. "
                       "Asks *which COPY of this seed carries the signal*, which is near-collinear by "
                       "construction. ⛔ It does NOT replace the family β; the honest deliverable is "
                       "`arm_resolvable` (a verdict on separability), not a confident per-arm β.",
    ("edge", "n_arm_in_cell"): "How many arms sit in this edge's (gene, seed_family) CELL — i.e. how many "
                               "of that family's members appear in THIS gene's design. ⚠ Not the family's "
                               "total membership: that is `famrole_n_members` on the arm card / the "
                               "seed_family card. >1 for 20.3% of edges; outside those the arm rung IS the "
                               "family rung by construction.",
    ("edge", "arm_detection"): "Fraction of patients in which this ARM is measured at all (`readouts._detection`). "
                           "⭐ MH-117: the ONE measurement proxy that arm-level β tracks within a family "
                           "(+0.305, p=4e-3) — abundance (p=0.65) and variance (p=0.95) do not. Emitted so a "
                           "reader can GATE on it; deliberately NOT used to drop arms (axiom 2a: flag, "
                           "don't delete). Arm-rung, so it repeats across the arm's genes.",
    ("edge", "arm_spiker"): "`arm_pct_floor < 40 AND arm_iqr > 1.5` — the arm sits at/below the detection floor "
                        "in most samples yet swings widely when present. Measured: spiker arms have median "
                        "`arm_pct_floor` 22 vs 99 for the rest; True on 171 edges. ⚠ A MEASUREMENT-RELIABILITY "
                        "WARNING, NOT BIOLOGY — a coupling carried by a few spiking samples is fragile: it "
                        "rests on a handful of points, so it is high-variance and easily sign-flipped. Gate "
                        "on it before reading any coupling on a spiker arm.",
    ("edge", "n_HLY_meas"): "How many of this gene's arms have a MEASURED healthy (GTEx) level — as opposed "
                            "to one manufactured by the multi-mapping collapse. ⛔ Read this BEFORE any "
                            "`share_HLY` / `rank_HLY` / `d_rank_HLY_*`: GTEx v10's uniquely-mappable pipeline "
                            "ZEROES the canonical member of several seed families (let-7a, miR-30a, miR-16, "
                            "miR-17, miR-200b…), which forced share=0 → rank LAST → dHT maximally positive "
                            "and FABRICATED healthy→tumour acquisition for 115 arms (MH-210). The `_meas` "
                            "variants recompute over measured arms only.",
    ("edge", "d_rank_HLY_TUM_meas"): "`d_rank_HLY_TUM` recomputed over arms with a MEASURED healthy level "
                                     "only — the collapse-safe version. Prefer it over the bare column.",
    ("edge", "dShare_M_own"): "Mean over the 103 paired patients of this arm's change in M-WEIGHTED share of "
                              "its gene's total regulator dose, tumour minus that patient's OWN NAT. Built "
                              "from linear RPM (2^x−1), column-normalised per patient so the gene's arms sum "
                              "to 1, then weighted by the model's `M`. ⇒ shares sum to 1 within gene, so the "
                              "dShares sum to ~0 — it is a REALLOCATION, not a level change. Defined only "
                              "where a gene has ≥2 regulators. P-across rung (one value per edge, computed "
                              "ACROSS patients — never a per-patient quantity).",
    ("edge", "dShare_raw_own"): "As `dShare_M_own` but on RAW abundance share, with the model's `M` weights "
                                "removed. ⭐ Read the PAIR: agreement means the reallocation is pure "
                                "abundance; divergence means the model's weighting is doing the work.",
    ("edge", "identity"): "Shapley/LMG share of the gene's R^2 credited to this edge, with bagged-NNLS "
                          "weights. ⭐ WHO, not how much. NaN where NNLS zeroed every family — that is an "
                          "honest 'undefined', not missing data.",
    ("edge", "m_nnls"): "Bagged-NNLS weight. Exactly 0 in 31.7% of families — this is what lets `identity` "
                        "say a regulator contributes nothing, which `beta` structurally cannot.",
    ("edge", "beta_frac_abs"): "|E[β_f]| / Σ|E[β]| — the share computed from the POINT ESTIMATES. "
                               "⚠⚠ **NOT the same as `beta_frac`, and the difference is JENSEN, not a "
                               "bug.** `beta_frac` is the mean over Gibbs DRAWS of (β_f/Σβ) — E[ratio] — "
                               "while this is the ratio of E. They differ on **4,985 of 5,644 rows**, and "
                               "the gap IS the posterior width: **spearman(beta_frac_sd, |gap|) = +0.809**, "
                               "median |gap| rising 0.00005 (tight) → 0.00844 (very wide), max **0.112**. "
                               "⇒ use `beta_frac` + `beta_frac_sd` when uncertainty matters; use this when "
                               "you need a plain bounded share of the point estimates. ✅ Verified: it "
                               "reconstructs β/Σβ over the gene's ARMS to 5e-4 (rounding), and both sum to "
                               "1 within gene on 100% of genes.",
    ("edge", "beta_frac_sd"): "Posterior SD of the per-draw fraction β_f/Σβ. ⭐ It is what makes "
                              "`beta_frac` an honest share rather than a point ratio — and it is the axis "
                              "that predicts the `beta_frac` vs `beta_frac_abs` gap (ρ=+0.809). "
                              "⚠ MH-94's PTEN miR-141/200a case (0.77 ± 0.41) is the reason this column "
                              "exists: a point Shapley hid a genuinely unidentified split.",
    ("_retired", "beta_frac_reliable"): "`net_pressure ≥ 0.5 AND |beta_frac_sd| ≤ 1`. ⛔⛔ **CONSTANT TRUE on "
                                    "all 5,644 rows — and that is the CORRECT answer, not a broken gate.** "
                                    "MH-119 built it because `beta_frac` exploded to 999% when a gene's βs "
                                    "cancelled; MH-124 then found those negative βs were a SAMPLER BUG and "
                                    "fixed it, so β is now strictly positive (min +0.00096) and `beta_frac` "
                                    "is properly bounded [0.0015, 1.0000] with 0 rows above 1. **The defect "
                                    "it guards can no longer occur.** ⚠ Vestigial — it carries no "
                                    "information and is a prune candidate; do NOT read it as evidence that "
                                    "a gate was applied to something that needed one. (Contrast "
                                    "`identity_reliable`, which guards the cancellation that DOES still "
                                    "occur — in `identity`, not in β.)",
    ("edge", "beta_deconv"): "β refit on the DECONV design (core C + the 8 non-malignant Wu-major "
                             "lineages, Cancer-Epithelial held out). ⚠ Not 'more controlled = better': a "
                             "miRNA acting through a cell-STATE shift PRODUCES a composition change, so "
                             "the composition is partly the MECHANISM. Read `retention_beta` beside it.",
    ("edge", "beta_sd_deconv"): "Posterior SD on the deconv design. ⚠ Genuinely different from `beta_sd` — "
                                "they differ on 5,597 of 5,644 rows (corr 0.999), so do not treat one as a "
                                "stand-in for the other.",
    ("edge", "beta_frac_deconv"): "`beta_frac` recomputed on the deconv design — the composition-adjusted "
                                  "share of the gene's budget.",
    ("edge", "beta_arm"): "The arm's OWN β from the within-cell fit (`arm_rung.py`), defined only in "
                          "multi-arm cells (fill 0.203). ⚠ Distinct from `beta`, which is the edge-level "
                          "readout over the gene's whole design. ⭐ Genuinely varies per arm inside a cell "
                          "(467/467 cells), unlike `arm_dbeta`/`arm_sep_z`, which are cell constants.",
    ("edge", "beta_frac"): "MAGNITUDE share (beta_f / sum beta). ⛔ Renamed from `share`: for an additive "
                           "aggregate Shapley(f) == beta_f, so this splits NOTHING. It is not identity.",
    ("gene", "concentration"): "`beta.max() / beta.sum()` — the TOP FAMILY's share of the gene's total β; median **0.954**, i.e. most genes are dominated by one family. "
                               "⛔⛔ **BOUNDED BELOW BY 1/n_fam, so it falls as designs widen WHATEVER the "
                               "biology does** — the moving-support trap (axiom 8) at its most extreme here: "
                               "**spearman(concentration, n_fam) = −0.9399**, versus **−0.2028** once the "
                               "floor is removed. The data sits ON the floor: at `n_fam==1` it is exactly "
                               "1.0000 for all **730 genes (47%)**, and at `n_fam==2` the minimum is 0.5005 "
                               "against a floor of 0.5000. ⇒ **never correlate it with anything that tracks "
                               "design width; use `concentration_adj`, or report it with `n_fam` beside it.**",
    ("gene", "concentration_adj"): "`(concentration − 1/n_fam) / (1 − 1/n_fam)` ∈ [0,1] — how concentrated "
                                   "the gene is RELATIVE to the most diffuse its own design allows. Removes "
                                   "the 1/k floor: spearman vs `n_fam` falls from **−0.9399 to −0.2028**. "
                                   "⚠ **NaN at `n_fam == 1` (730 genes) and that is CORRECT** — with one "
                                   "family there is nothing to concentrate, so the question does not exist; "
                                   "the NaN is the honest answer, not a gap (`gene_axes.mask_degenerate`).",
    ("", "disc_"): "⭐ THE 157-EDGE DISCOVERY QUEUE, made visible from the cards. ⛔ It CANNOT be an edge "
                   "flag: the gold set and the curated edge card are **DISJOINT — 0 of 157 pairs are on "
                   "the edge card** — because they are candidates on pairs curation never admitted. So it "
                   "is surfaced as a COUNT at the two rungs that overlap. ⚠⚠ **COVERAGE IS ASYMMETRIC:** "
                   "all **18/18** gold arms are on the arm card (accounting for all 157 edges), but only "
                   "**39 of 90** gold genes are on the gene card (63 edges) — the rest have **no curated "
                   "regulator at all**, which is what made them candidates. ⇒ **a 0 on the GENE card means "
                   "'no gold edge OR invisible here', never 'the queue is empty for this gene'.** Read the "
                   "queue from `discovery_gold_set.tsv`; treat these as a pointer, never a denominator. "
                   "⚠ 11 seed families, **61% is ONE family** — quote families beside edges. ⚠ A ranked "
                   "VALIDATION QUEUE, not 157 findings: per-edge FDR is empty by construction.",
    ("arm", "disc_n_gold_edges"): "How many of the 157 discovery-queue edges name this arm. Concentrated: "
                                  "miR-106b-5p **39** · miR-20a-5p 25 · miR-93-5p 25 · miR-19a-3p 19 · "
                                  "miR-18a-5p 9. ⚠ 18 arms carry the whole queue.",
    ("arm", "disc_gold_families"): "Which seed families this row's discovery-queue edges belong to, semicolon-separated. Usually ONE at the arm rung (an arm has one seed). \u26a0 The whole queue is 11 families and **61% is miR-17/20/93/106/519** \u2014 the honest unit is ~one oncogenic polycistron's realized target set, not 157 independent findings.",
    ("gene", "disc_gold_families"): "Which seed families this row's discovery-queue edges belong to, semicolon-separated. Usually ONE at the arm rung (an arm has one seed). \u26a0 The whole queue is 11 families and **61% is miR-17/20/93/106/519** \u2014 the honest unit is ~one oncogenic polycistron's realized target set, not 157 independent findings.",
    ("gene", "disc_n_gold_edges"): "How many of the 157 discovery-queue edges name this gene (max 4). "
                                   "⛔ Only **39 of the 90** gold genes appear on this card, so this "
                                   "accounts for **63 of 157** edges — a 0 is NOT evidence of absence.",
    ("gene", "n_arms"): "Arms in the LEARNED MODEL's design — `X.shape[1]` from `assemble_gene`, carried "
                        "via `gene_atlas`. This is the width the Gibbs actually fit. ⚠⚠ **NOT the same as "
                        "`n_regulators`**, which counts distinct arms in the RETIRED heuristic lane's edge "
                        "table: they agree on 1,103 genes and differ on 306, in BOTH directions (ABCA1 13 "
                        "vs 12, ABCC1 7 vs 8). Both read as 'how many regulators'; only this one describes "
                        "the fit.",
    ("", "heur_"): "⛔⛔ **THE RETIRED-HEURISTIC BLOCK — provenance, not model.** All five columns come from "
                   "`analyses/misc/mirna_comovement.py` → `tissue_reference/mirna_comovement/"
                   "gene_corepression.tsv`, i.e. the **§6b-RETIRED pressure heuristic**, NOT the learned "
                   "Gibbs. They were renamed as a BLOCK on 2026-08-19 so the provenance is visible at a "
                   "glance and they read as one unit. ⚠ Renamed **on the card only** — the source TSV and "
                   "its direct readers keep the original names. ⚠ Do not mix them with model columns in one "
                   "statement without saying which lane each came from.",
    ("gene", "heur_repression_class"): "Cross-state repression class from the RETIRED heuristic lane "
                                       "(never_repressed / lost_repression / gained_repression / "
                                       "constitutive_repressed). 867 of 1,409 genes are `never_repressed`.",
    ("gene", "heur_net_repressed_tumor"): "Is the gene net-repressed in tumour, per the RETIRED heuristic "
                                          "lane. True for 282 of 1,409 genes.",
    ("gene", "heur_rho_pressure_tumor"): "Spearman of the gene's HEURISTIC aggregate pressure against its "
                                         "mRNA in tumour. ⛔ The heuristic is retired — the learned "
                                         "equivalent is `realized_rho_adj`. Do not present this as the "
                                         "model's coupling.",
    ("gene", "heur_delta_tumor_nat"): "Tumour-minus-NAT change in the heuristic pressure correlation. "
                                      "⛔ Retired lane; the learned paired-design equivalent lives in the "
                                      "`realization` columns.",
    ("gene", "heur_n_regulators"): "⚠⚠ **A DIFFERENT PIPELINE'S COUNT — not the model's.** It is "
                              "`edges.groupby('gene')['miRNA'].nunique()` inside "
                              "`analyses/misc/mirna_comovement.py`, i.e. distinct arms in the **heuristic "
                              "pressure lane's edge table** — and that lane is the **§6b-RETIRED** pressure "
                              "heuristic. `n_arms` is the LEARNED MODEL's design width (`X.shape[1]` from "
                              "`assemble_gene`, via `gene_atlas`). ⇒ two counts, two universes, two "
                              "admission rules. **Measured: they agree on 1,103 genes, the model has MORE "
                              "on 256, FEWER on 50, range −15…+3**, and 140 gene-card genes have NO "
                              "`n_regulators` at all (co-repression covers 1,424 of 1,549). **Use `n_arms` "
                              "for anything about the FIT.** Read this only beside the co-repression "
                              "columns it arrived with (`gene_repression_class`, "
                              "`gene_net_repressed_tumor`, `rho_gene_pressure_tumor`, `delta_tumor_nat`).",
    ("gene", "n_identified"): "Families reaching |z| > 2 on the UNCALIBRATED SD. **Zero for 691 genes "
                              "(44.6%)** — meaning the model cannot distinguish ANY of that gene's families "
                              "from zero. ⛔ That is *cannot resolve WHO*, **not** *no regulation*. "
                              "⭐ It is a MEASURABILITY fact, and both gradients are steep: by design width "
                              "**66.6% (1 family) → 43.4% → 28.9% → 15.2% → 3.4% (6+)**, and by the "
                              "a-priori `ctx_ceiling` **84.9% (≤0) → 47.7% → 26.2% → 10.7% → 1.7% (>0.15)**. "
                              "The ceiling gradient is the cleaner one: it is decoy-independent and never "
                              "touches β. ⚠ A `total_pressure` split looks even sharper (96.1% → 12.6%) but "
                              "is PARTLY CIRCULAR — `total_pressure` is Σβ and z is β/β_sd. "
                              "⚠ Under the CALIBRATED width it would be lower still.",
    ("gene", "n_pip_disc_gt50"): "Count of the gene's families with `pip_discovery > 0.5` under the "
                                 "evidence-π readout. ⛔ **RENAMED from `n_discovered` 2026-08-19 because "
                                 "the old name read as a RESULT**: it is >0 for **725 genes (46.8%)**, while "
                                 "per-edge and per-family discovery are **EMPTY** under the honest "
                                 "empirical FDR — the lane's deliverable is a convergent-evidence QUEUE "
                                 "(157 edges / 11 seed families), not per-gene discoveries. ⛔ Also "
                                 "**w-CONTAMINATED** by construction (`pip_discovery` rides the evidence "
                                 "prior), so any 'do canonical regulators score higher' test on it is "
                                 "CIRCULAR.",
    ("gene", "n_cell_intrinsic"): "Edges whose coupling SURVIVES the composition block (retention ≥ 0.7). "
                                  "⚠ With `n_composition_explained` it does NOT partition the gene's edges — "
                                  "verified: the two sum to ≤ `n_arms` on 100% of genes (median sum 1 vs "
                                  "n_arms 2), because the middle class `partial` belongs to neither.",
    ("gene", "n_composition_explained"): "Edges whose coupling does NOT survive composition adjustment "
                                         "(retention < 0.4). ⚠ See `n_cell_intrinsic`: the two are the ENDS "
                                         "of a three-way classification, not a partition. ⚠⚠ And a "
                                         "composition-explained edge is not necessarily an artifact — a "
                                         "miRNA acting in CAFs on a CAF gene is removed by conditioning on "
                                         "CAF fraction and then labelled this (axiom 8's stratified-retention "
                                         "warning). Pair it with a specificity control.",
    ("gene", "n_hallmark_sets"): "How many MSigDB Hallmark programs contain this gene (1–10, median 2). "
                                 "Membership only — it says nothing about regulation.",
    ("gene", "w_max"): "The MAXIMUM curated evidence weight `w` over the gene's regulators "
                       "(`gene_atlas`: `nanmax(w)`), where `w` comes from the same PMID-deduped ledger that "
                       "feeds `lit_` and `fame_`. ⛔ NOT a β and NOT a dose share — `top_beta_frac` and "
                       "`concentration` are the shares. ⚠⚠ It is therefore a STUDY-DEPTH quantity, and "
                       "`w_max > median` is HALF the definition of `A_COMPETENT` — so the competence class "
                       "is part fame axis. ⚠ It carries THREE rungs across the cards: gene here, gene on the "
                       "edge card (repeated), and family on `gene_family` — where the label is WRONG "
                       "(measured 2026-08-19: 0 of 1,549 genes have >1 distinct value across their "
                       "families, so it is gene-rung repeated there too).",
    ("gene_family", "w_max"): "Maximum curated EVIDENCE weight for this (gene, family) cell. Family-rung "
                              "here; gene-rung on the gene card.",
    ("gene", "identity_eq_magnitude"): "Do the identity and magnitude answers name the SAME top family. "
                                       "Exactly 100% agreement at n_fam==1 by construction; 76% at "
                                       "n_fam>=2. Read it only on multi-family genes.",
    ("gene", "regulatory_handoff"): "Does the gene's globally-dominant regulator differ between healthy and "
                                    "tumour OR between NAT and tumour. Fires on 328/1200 genes where "
                                    "defined; 301 are genuine healthy->tumour switches, 27 fire only via "
                                    "the NAT leg. ⚠ Mechanically requires >=2 regulators — 1.6% at "
                                    "n_fam==1 vs 38.6% at n_fam>=2, so ALWAYS condition on design width.",
    ("gene", "ctx_measurable"): "Is the gene's coupling measurable at all given its design (positive OOF "
                                "ceiling). False for 877/1549 genes. ⚠ Partly definitional — a "
                                "non-positive ceiling forces the real correlation to ~0.",
    ("gene", "ctx_apriori_class"): "A-priori competence class from the gene atlas (A_COMPETENT / B_PARTIAL "
                                   "/ C_WEAK / D_UNDEFINED). ⛔ There is NO competence partition — the "
                                   "decoy gap is a GRADIENT in design width and is nowhere exactly zero. "
                                   "Use the class as a stratifier, never as a filter.",
    ("gene", "ctx_gap_core"): "Real-minus-decoy coupling gap for this gene. Pooled it is ~-0.012, but the "
                              "per-gene IQR STRADDLES ZERO — this is a set-level distributional shift, not "
                              "a per-gene effect. Do not rank genes by it.",
    ("edge", "oof_drho"): "Arm-level OOF rho MINUS family-level. Negative = arm resolution predicts better. "
                          "Median only -0.0038, but the tail reaches -0.21 — the median understates it "
                          "~50x. Characterise the tail, never the median.",
    ("edge", "ago_loading"): "AGO/RISC loading discount. ⛔ FABRICATED for 32.4% of edges: the producer "
                             "ends `.fillna(1.0)`, so an unmeasured arm is indistinguishable from a "
                             "measured guide. Use `ago_loading_measured` to separate them.",
    ("edge", "ago_loading_measured"): "Was `ago_loading` actually MEASURED for this edge (vs defaulted to "
                                      "1.0). True for 3,820/5,648 edges (67.6%).",
    ("edge", "adm_has_site"): "Does the arm have a seed match in this target's 3'UTR. ⭐ The project's most "
                              "load-bearing conditioning variable: genes whose curated edges are ALL "
                              "seedless show a decoy gap of exactly zero (+0.0006) vs -0.0325 for "
                              "all-seeded genes.",
    ("edge", "cal_identified"): "Is the edge identified (|z| > 2) under the CALIBRATED posterior width. "
                                "19.9% [18.1–22.4] of arm-edges — down from 24.8% on the uncalibrated SD.",
    ("gene", "top_identity"): "The largest `identity` share among the gene's families. ⛔⛔ **UNGATED AND "
                              "UNBOUNDED — it still reaches +740.007, with 80 genes above 1.** `identity` "
                              "is a SIGNED share (~10% negative), so its max has no upper bound; the "
                              "`identity_reliable` gate shipped to the edge and gene_family cards but "
                              "**never reached this DERIVED column**, because gating the inputs does not "
                              "gate a statistic already computed from them. ⇒ **use `top_identity_gated`.** "
                              "⚠ The blow-ups concentrate where the denominator vanishes: 75% of the "
                              "affected genes have `ctx_ceiling` ≤ 0.02 vs 31.5% of the rest.",
    ("gene", "top_identity_gated"): "`top_identity` recomputed over `identity_reliable` edges only — "
                                    "**max 1.0000, 0 genes above 1** (vs +740.007 ungated). Defined on "
                                    "1,203 genes. ⭐ Ships BESIDE the raw column rather than replacing it: "
                                    "`readouts` owns `top_identity`, and silently rewriting another "
                                    "module's value is how provenance rots.",
    ("gene", "top_identity_n_reliable"): "How many of the gene's families survive the identity gate — the "
                                         "DENOMINATOR behind `top_identity_gated`. **Median 2.** ⚠ A max "
                                         "over 2 families is not the same statement as a max over 12; "
                                         "axiom 5 says print the denominator, so it is a column.",
    ("arm", "arb_max_identity"): "Largest `identity` among this arm's edges, over reliable rows only. "
                                 "✅ Now runs **−0.255 … +1.000** — it reached **+740.0** until the inert "
                                 "`beta_frac_reliable` gate was replaced (2026-08-19). Read `arb_n_identity_"
                                 "reliable` beside it for the denominator.",
    ("gene", "top_family_identity"): "The family Shapley/LMG credits with the gene's regulation. NaN for "
                                     "197 genes (12.7%) where NNLS zeroed everything — an honest "
                                     "'undefined', which the magnitude answer cannot express.",
    ("gene", "top_family_magnitude"): "The family with the largest beta. ⚠ This is argmax over a quantity "
                                      "that is never zero, so it always returns SOMETHING.",
    ("gene", "lit_agrees_identity"): "Does `identity` name the same family as the **most-published** one for "
                                     "this gene. True for 100/323. ⛔ NOT an accuracy measure — the reference "
                                     "is a study count (see the `lit_` block). Against `magnitude` (96/329) "
                                     "the contrast is **14 vs 9 discordant genes, McNemar exact p=0.405 — "
                                     "a tie**.",
    ("edge", "shift_class"): "Cross-state class for the edge. ⚠ Its per-edge FDR labels rest on the "
                             "uncalibrated null (3–4x too narrow); read the classes, not their counts.",
    ("gene", "gene_repression_class"): "Cross-state repression class for the gene (never_repressed / "
                                       "lost_repression / gained_repression / constitutive_repressed). "
                                       "867 of 1409 genes are never_repressed.",
}

# ---- keys, counts and the remaining bare-name columns (filled 2026-08-18) ----
_KEYS = {"gene": "KEY — HGNC gene symbol.",
         "arm": "KEY — mature miRNA arm (e.g. hsa-miR-21-5p). Normalised through `_arm_key()`, which "
                "gates against the canonical universe and LOGS drops rather than silently dropping.",
         "family": "KEY — seed family of the regulator, as used by the family-rung fit.",
         "seed_family": "KEY — the seed family itself, one row per family."}
for _c in ("arm", "edge", "gene", "gene_family", "seed_family"):
    for _k, _t in _KEYS.items():
        COLUMNS.setdefault((_c, _k), _t)

COLUMNS.update({
    ("_retired_edge", "p_fam"): "\u26d4\u26d4 **NOT A P-VALUE.** `p` here is the model's DESIGN DIMENSION \u2014 the NUMBER "
                       "of family predictors in this gene's fit (integers 1..90). The name reads as a "
                       "p-value and is a genuine collision; verified 2026-08-19 against the values.",
    ("gene_family", "p_fam"): "\u26d4\u26d4 **NOT A P-VALUE** \u2014 the number of family predictors in the gene's "
                              "fit. **Bit-identical to `n_fam` on this card** (integers 1..12), so it is also "
                              "a redundant column.",
    ("edge", "z"): "beta / beta_sd on the UNCALIBRATED posterior SD. Use `cal_z` for the honest width.",
    ("edge", "identified"): "|z| > 2 on the UNCALIBRATED SD (24.8%). `cal_identified` is the honest "
                            "version (19.9%).",
    ("edge", "pip_dense"): "\u26a0 CONSTANT by construction \u2014 the coupling readout is called with \u03c0\u22611, so "
                           "every edge has the same value. Carries no information; `pip_discovery` is the "
                           "one that varies. (Bit-identical to `net_pressure` on this card.)",
    ("gene_family", "pip_dense"): "\u26a0 CONSTANT by construction (\u03c0\u22611 dense readout). Bit-identical to "
                                  "`net_pressure`.",
    ("_retired", "net_pressure"): "\u26a0 CONSTANT on this card and bit-identical to `pip_dense`.",
    ("edge", "n_samples"): "Samples the edge was fit on \u2014 CONSTANT (1,040) across the card.",
    ("gene_family", "n"): "Samples the cell was fit on \u2014 CONSTANT (1,040) across the card.",
    ("gene", "n_dense_included"): "Families entering the dense (\u03c0\u22611) readout. \u26a0 **Bit-identical to "
                                  "`n_fam`** \u2014 redundant; do not treat as an independent axis.",
    ("", "identity_coherence"): "1 / Σ|identity| per GENE, in (0,1] — identity's own sign coherence. 1.0 "
                                "when nothing cancels, → 0 as suppressor cancellation grows. ⭐ This is the "
                                "quantity `beta_frac_reliable` was supposed to watch and could not: it "
                                "guards **β** cancellation, which MH-124's sampler fix made impossible, "
                                "while the cancellation that actually occurs is in `identity` (~10% of "
                                "values are negative).",
    ("", "identity_abs"): "|identity| / Σ|identity| per gene — the BOUNDED companion to `identity`, always "
                          "in [0,1] (verified: 0 non-NaN values outside it on either card). Use it when a "
                          "share is wanted and the sign is not; NaN exactly where `identity` is NaN.",
    ("", "identity_reliable"): "`identity_coherence ≥ 0.5` AND `|identity| ≤ 1`. Admits **93.3%** of edge "
                               "rows and **92.6%** of gene_family rows, and excludes **every** row whose "
                               "|identity| exceeds 1 (119/119 and 118/118). ⛔ `identity`, `top_identity` "
                               "and `arb_max_identity` are only interpretable on these rows — ungated they "
                               "reach +166 and +740.",
    ("gene_family", "n_arms"): "How many arms of THIS family sit in THIS gene's design. ✅ FAMILY rung, "
                               "VERIFIED (2026-08-19): it varies within gene for 241 of 1,549 genes, and "
                               "the per-gene SUM of these equals the gene card's `n_arms` exactly — a clean "
                               "cross-rung consistency check. ⚠ Do NOT confuse with the family's TOTAL "
                               "membership, which is `famrole_n_members` / the seed_family card; a gene "
                               "sees a median 0.67 of its family (MH-248).",
    ("gene_family", "arms"): "The arms collapsed into this (gene, family) cell, **semicolon-separated** "
                             "(`hsa-miR-17-5p;hsa-miR-20a-5p`). ⚠ Verified: splitting on `;` reproduces "
                             "`n_arms` for EVERY row. It is a plain string, NOT a Python list literal and "
                             "NOT comma-separated — both mis-parses were made while writing this entry.",
    ("gene_family", "retention"): "`beta_deconv / beta_core` at the family grain — the ADJUSTMENT sense of "
                                  "retention (what fraction of the coupling survives the composition "
                                  "block), never the patient-baseline sense. ⚠⚠ **41.4% of rows exceed 1** "
                                  "(adjustment INCREASED the effect) and it reaches **9.24** — a ratio whose "
                                  "denominator can vanish. ✅ **`retention_reliable` is the gate and it "
                                  "WORKS**: it admits 27.1% of rows, and among those only **8.9% exceed 1** "
                                  "versus **53.5%** of the rest. Read `retention` ONLY where "
                                  "`retention_reliable` is true.",
    ("gene_family", "m_nnls"): "Bagged-NNLS weight for this (gene, family) cell. **Exactly 0 in 30.1%** of "
                               "rows — this is what lets `identity` say a family contributes nothing, which "
                               "`beta` structurally cannot (the half-normal slab has a strictly positive "
                               "mean, so `beta` is exactly 0 in 0 rows).",
    ("gene_family", "identity"): "Shapley/LMG identity at the FAMILY rung — same estimator as the edge "
                                 "card's `identity`, different unit. ⛔⛔ **CARRIES THE SAME UNBOUNDED "
                                 "SIGNED-SHARE PROBLEM, AND IS NOT YET GATED HERE.** Measured 2026-08-19: "
                                 "range **−739.01 … +740.01**, **9.5% negative** (legitimate suppressor "
                                 "contributions), 83 rows above 1, summing to exactly 1 per gene. "
                                 "⏳ `identity_reliable` / `identity_coherence` / `identity_abs` were added "
                                 "✅ **GATED SINCE 2026-08-19 — and NO refit was needed.** The gate is a "
                                 "PURE FUNCTION of `identity` and its per-gene sum, so "
                                 "`card_ladders.family_identity_gate()` derives it in place and re-derives "
                                 "it after every rebuild. Verified: it admits 92.6% of rows and excludes "
                                 "**118 of 118** with |identity|>1; gated `top_identity` maxes at exactly "
                                 "+1.0000 (ungated: +740.0). **Read `identity` only where "
                                 "`identity_reliable`, or read the bounded `identity_abs` instead.**",
    ("gene_family", "n_fam"): "How many families compete at this gene. beta is mathematically INERT at 1.",
    ("gene", "n_fam"): "How many seed families regulate this gene. \u2b50 48% of genes have exactly ONE, where "
                       "beta \u2261 uniform by construction \u2014 'does the learning add anything' is UNDEFINED "
                       "for half the universe, which is a non-question, not a null.",
    ("", "disc_"): "⭐ THE 157-EDGE DISCOVERY QUEUE, made visible from the cards. ⛔ It CANNOT be an edge "
                   "flag: the gold set and the curated edge card are **DISJOINT — 0 of 157 pairs are on "
                   "the edge card** — because they are candidates on pairs curation never admitted. So it "
                   "is surfaced as a COUNT at the two rungs that overlap. ⚠⚠ **COVERAGE IS ASYMMETRIC:** "
                   "all **18/18** gold arms are on the arm card (accounting for all 157 edges), but only "
                   "**39 of 90** gold genes are on the gene card (63 edges) — the rest have **no curated "
                   "regulator at all**, which is what made them candidates. ⇒ **a 0 on the GENE card means "
                   "'no gold edge OR invisible here', never 'the queue is empty for this gene'.** Read the "
                   "queue from `discovery_gold_set.tsv`; treat these as a pointer, never a denominator. "
                   "⚠ 11 seed families, **61% is ONE family** — quote families beside edges. ⚠ A ranked "
                   "VALIDATION QUEUE, not 157 findings: per-edge FDR is empty by construction.",
    ("arm", "disc_n_gold_edges"): "How many of the 157 discovery-queue edges name this arm. Concentrated: "
                                  "miR-106b-5p **39** · miR-20a-5p 25 · miR-93-5p 25 · miR-19a-3p 19 · "
                                  "miR-18a-5p 9. ⚠ 18 arms carry the whole queue.",
    ("arm", "disc_gold_families"): "Which seed families this row's discovery-queue edges belong to, semicolon-separated. Usually ONE at the arm rung (an arm has one seed). \u26a0 The whole queue is 11 families and **61% is miR-17/20/93/106/519** \u2014 the honest unit is ~one oncogenic polycistron's realized target set, not 157 independent findings.",
    ("gene", "disc_gold_families"): "Which seed families this row's discovery-queue edges belong to, semicolon-separated. Usually ONE at the arm rung (an arm has one seed). \u26a0 The whole queue is 11 families and **61% is miR-17/20/93/106/519** \u2014 the honest unit is ~one oncogenic polycistron's realized target set, not 157 independent findings.",
    ("gene", "disc_n_gold_edges"): "How many of the 157 discovery-queue edges name this gene (max 4). "
                                   "⛔ Only **39 of the 90** gold genes appear on this card, so this "
                                   "accounts for **63 of 157** edges — a 0 is NOT evidence of absence.",
    ("gene", "n_arms"): "How many miRNA arms regulate this gene (median 2).",
    ("gene", "n_regulators"): "Regulator count used by the realization lane.",
    ("gene", "n_identified"): "How many of the gene's families are identified (|z|>2). Median 1.",
    ("gene", "n_discovered"): "Families passing the evidence-pi discovery readout. \u26d4 Per-edge and "
                              "per-family discovery are EMPTY under the honest FDR \u2014 read this as a "
                              "ranked queue, never as a set of findings.",
    ("gene", "n_composition_explained"): "Edges whose coupling does NOT survive the composition block.",
    ("gene", "n_cell_intrinsic"): "Edges whose coupling DOES survive composition adjustment.",
    ("gene", "n_hallmark_sets"): "How many MSigDB Hallmark programs this gene belongs to.",
    ("edge", "n_HLY_meas"): "How many of the gene's arms have a MEASURED (not imputed, not "
                            "multi-mapping-collapsed) healthy level. Read before any `share_HLY`/`rank_HLY`.",
    # ── UNIT 19 (2026-08-19) — edge `identity_` (5) + arm `ago_` (4) + arm `cur_` (2).
    ("edge", "identity_reliable"): "Is this edge's Shapley `identity` trustworthy — `identity_coherence` "
        "≥ 0.5 AND |identity| ≤ 1? ⛔⛔ **MASKED 2026-08-19: it could not say 'unknown'.** A bare `&` of "
        "two comparisons returns False on NaN, so an edge whose identity was NEVER COMPUTED read as "
        "*unreliable* — **218 of the 377 False rows (57.8%)**. Anyone filtering `== False` to study "
        "unreliable attribution got a set that was 58% *not computed*, and the unreliable RATE was "
        "overstated **2.3×: 6.7% → 2.9%** on the measurable set. Now three-state: **True 5,272 · False 159 "
        "· NaN 218.**",
    ("edge", "identity_abs"): "|identity| as a fraction of the gene's total |identity|, **clipped to "
        "[0,1]**. ⚠⚠ **THE CLIP IS LOSSY AND HIDES THE TAIL: 550 rows (9.7%) sit at exactly 1.000**, which "
        "may be a genuine 1.0 or a clipped 740. ⇒ read `identity_reliable` beside it; a clipped row is by "
        "construction one the gate rejects.",
    ("edge", "identity_deconv"): "Shapley identity recomputed under the composition-adjusted design. "
        "⛔⛔ **UNGATED AND UNBOUNDED — range [−540.3, +541.3].** `identity` is a SIGNED share, so it "
        "explodes wherever the gene's Σ|identity| approaches 0 (axiom 5). **Never quote it raw**; gate on "
        "`identity_reliable` or use `identity_abs`.",
    ("edge", "identity_allocated"): "**Family identity with phantom minor arms forced to 0 — the "
        "allocation actually used downstream.** ⛔ **Ungated — "
        "range [−180.4, +166.0]**, same signed-share pathology as `identity_deconv`. Gate before use.",
    ("edge", "identity_coherence"): "`1 / Σ|identity|` over the gene, clipped to (0,1] — how much the "
        "gene's identity shares CANCEL. Median **0.996**, i.e. the shares nearly sum to 1 for most genes; "
        "the low tail is where the signed shares fight and the ratio columns blow up. The `≥ 0.5` half of "
        "`identity_reliable`.",
    ("arm", "ago_dom_src"): "Provenance of `ago_dom` — which CLIP dataset the AGO-loading dominance came "
        "from. ⚠ **CONSTANT: `manakov` on all 678 rows**, because only one source is wired. KEPT rather "
        "than pruned for the same reason as `pip_dense`: it records a PROVENANCE choice, and it would "
        "start varying the moment a second dataset is added. A constant column is not automatically dead — "
        "ask whether it names a design decision or a defunct check.",
    ("arm", "ago_dom"): "Fraction of the arm's Manakov AGO-CLIP reads assigned to it rather than its "
        "sibling arm — the loading dominance. Median 0.658 on the 27.7% of arms with CLIP coverage. "
        "⚠ MH-record: this axis is **coupling-INERT** — its value is identity/QC, not prediction.",
    ("arm", "ago_reads"): "Raw AGO-CLIP read count (median 6, max 210,722). ⚠ Fill **60.0%** against "
        "`ago_dom`'s 27.7% — a read count exists for many arms whose dominance is undefined, because "
        "dominance needs BOTH arms of the hairpin. Do not infer coverage of one from the other.",
    ("arm", "cur_he_degree"): "How many HE genes this arm regulates in the curated design — the arm's "
        "curation degree (median 3, max 125).",
    ("arm", "cur_he_degree_expr"): "As `cur_he_degree` but counting only targets above the expression "
        "floor. ⛔ **NOT redundant with it, though it looks it** — identical on only **57.4%** of the 875 "
        "arms, differing on **373** (max 18 targets), spearman 0.974. The gap IS the unexpressed-target "
        "fraction; a degree that shrinks under the expression filter is an arm whose curated targets are "
        "largely silent in this cohort.",
    # ── UNIT 18 (2026-08-19) — the `comp_` block (16 cols, IDENTICAL on the gene and gene_family cards).
    ("", "comp_tcga_mrna_driver_share"): "Share of the gene's mRNA coupling attributable to the single "
        "compartment that most changes it — **axiom 8's composition GATE**, the axis that separates the "
        "harm tails at p=0.00999 where the ensemble axis gives p=0.10. ⛔⛔ **EXACTLY COMPLEMENTARY TO "
        "`comp_tcga_mrna_driver_ret`: share ≡ 1 − ret, verified to 4.4e-16 on all 1,027 rows.** They are "
        "ONE quantity in two orientations — a Spearman scan sees ρ and −ρ, so **scanning both doubles the "
        "apparent evidence and wastes a test**. ✅ `gene_axes.py` builds only this one; keep it that way.",
    ("", "comp_tcga_mrna_driver_ret"): "`1 − comp_tcga_mrna_driver_share`, exactly (see there). The "
        "RETENTION orientation: what fraction of the coupling survives removing the driver compartment. "
        "⚠ **Denominator gated at `RHO_GATE = 0.05`, and that is too low for a stable readout: 11.3% of "
        "rows still fall outside [0,1].** At |ρ_raw| ≥ 0.10 that drops to **4.6%** and the median moves "
        "+0.646 → **+0.720** ⇒ quote the 0.10-gated value (axiom 5c: sweep the threshold, report where the "
        "arms agree).",
    ("", "comp_cptac_prot_driver_share"): "Protein-layer twin of `comp_tcga_mrna_driver_share`. Same exact "
        "complementarity with `_driver_ret` (8.9e-16 on 918 rows). ⚠ Noisier layer: **24.2%** of rows fall "
        "outside [0,1] even at the 0.05 gate, dropping to **16.2%** at |ρ_raw| ≥ 0.10.",
    ("", "comp_cptac_prot_driver_ret"): "`1 − comp_cptac_prot_driver_share`, exactly. Median +0.398 raw, "
        "**+0.469** gated to |ρ_raw| ≥ 0.10 — quote the gated value. See the mRNA twin.",
    ("", "comp_tcga_mrna_driver"): "WHICH compartment drives the gene's mRNA coupling — the one whose "
        "removal moves it most (`CAFs` 279 · `purity` 174 · `prolif` 148 · `Endothelial` 122 · `Myeloid` "
        "109). ⭐⭐ **NOT the same as either `top_` column, and the gap is large: the driver is neither the "
        "budget's nor the target's top compartment on 36.1% of genes** (protein: 41.9%), and it matches "
        "`top_budget` on only 39.2% / `top_target` on 41.9%. This is the measured form of `gene_axes.py`'s "
        "loading-vs-driver trap: *'is this gene expressed by stroma'* is not *'is this gene's "
        "miRNA-target relationship an artifact of stroma'*.",
    ("", "comp_tcga_mrna_top_budget"): "The compartment the gene's miRNA BUDGET loads on most — a property "
        "of the REGULATORS. ⚠ Agrees with `comp_tcga_mrna_top_target` on only **20.4%** of genes (protein "
        "17.2%), so the two sides of the edge load on different compartments most of the time. Neither is "
        "the driver — see `comp_tcga_mrna_driver`.",
    ("", "comp_tcga_mrna_top_target"): "The compartment the TARGET GENE's expression loads on most. ⚠ All "
        "11 loading axes were NULL under FDR in MH-201 (best q=0.18; CAFs p=0.59 with the sign running the "
        "opposite way) ⇒ **never report a loading null as if it refuted a composition effect** — the "
        "driver-share axis is the one that predicts.",
    ("", "comp_tcga_mrna_rho_raw"): "The gene's UNADJUSTED aggregate mRNA coupling — the DENOMINATOR of "
        "`comp_tcga_mrna_driver_ret`. Fill 90.2%, median −0.056. Read it beside any retention: a ratio "
        "whose denominator sits near zero is a coin-flip with decimals (axiom 5).",
    ("", "comp_cptac_prot_rho_raw"): "Unadjusted aggregate PROTEIN coupling — denominator of "
        "`comp_cptac_prot_driver_ret`. Fill 73.2%, median −0.018, i.e. **centred much closer to zero than "
        "the mRNA twin**, which is why the protein retention is the noisier of the two.",
    # ── UNIT 17 (2026-08-19) — the `cptac_` block (24 cols, two cohorts × three layers).
    ("", "cptac_prosp_agg_beats_abund_prot"): "Does the learned β-weighted aggregate couple more negatively "
        "to PROTEIN than a plain abundance sum, in the CPTAC prospective cohort? ⛔⛔ **THIS COLUMN WAS "
        "WRONG UNTIL 2026-08-19.** It was a bare `rho_agg < rho_abund`, and **`NaN < NaN` is False, not "
        "NaN** — so it was defined on all **1,420** genes while its inputs cover **1,127**, and **all 293 "
        "unmeasurable rows read False, none True.** ⇒ the headline rate was deflated **25.1% → 31.7%** once "
        "the denominator is the measurable set. ✅ The comparison itself was always correct where both "
        "inputs exist (100% match). ⚠⚠ **AND IT IS A REFERENCE, NOT A CONTROL** — a FITTED β against an "
        "UNFITTED abundance sum; the fitted-fake decoy is the control (`eval/decoy_bench.py`).",
    ("", "cptac_t105_agg_beats_abund_prot"): "As the prospective twin, on the TCGA-105 overlap cohort. Same "
        "`NaN < NaN → False` defect, fixed 2026-08-19: **488** unmeasurable rows had read False, deflating "
        "the rate **21.1% → 32.1%** (299 of 932 measurable). ⭐ **Note what the fix did to the cohort "
        "comparison: 25.1% vs 21.1% became 31.7% vs 32.1% — the two cohorts AGREE, and the 4-point gap was "
        "manufactured entirely by unequal missingness.**",
    ("", "cptac_prosp_agg_ret_prot"): "`rho_prot(composition-adjusted) / rho_prot_raw` — the ADJUSTMENT "
        "sense of retention, on the prospective cohort. ✅ **The denominator is ALREADY GATED** at "
        "`RHO_GATE` in `cptac_card.py:322`, so values outside [0,1] are NOT the vanishing-denominator "
        "artifact of axiom 5 (verified: 0 rows with |raw| < 0.02). They are real. Of **912** defined rows: "
        "**56.7% attenuate (0–1, the expected case)**, **21.8% FLIP SIGN (<0)**, **21.5% AMPLIFY (>1)**. "
        "⚠ The flips sit on WEAKER raw coupling (|ρ_raw| median **0.127 vs 0.164**), so read them as "
        "noise-driven sign instability, not systematic reversal. Gated to |ρ_raw| ≥ 0.10 the median is "
        "**+0.450**. ⚠ Runs WITH a composition block — MH-107's defect (no composition variant at all) is "
        "fixed here; see memory `cptac-protein-composition-confound`.",
    ("", "cptac_t105_agg_ret_prot"): "As the prospective twin, TCGA-105. Of **649** defined rows: **59.8% "
        "attenuate**, **10.3% flip sign**, **29.9% amplify**; median gated to |ρ_raw| ≥ 0.10 is **+0.748** "
        "— markedly higher retention than prospective's +0.450, on half the sample. ⚠ Do not read that as "
        "a cohort biology difference without an n-matched comparison.",
    ("", "cptac_prosp_n_arms"): "Arms contributing to this gene's prospective-cohort aggregate (median 2, "
        "max 61). ⚠ The DENOMINATOR behind every `prosp_agg_*` column — a ρ built on 2 arms is not the "
        "same statement as one built on 61.",
    ("", "cptac_t105_n_arms"): "Arms behind the TCGA-105 aggregate (median 2, max 57). See the prospective "
        "twin.",
    # ── UNIT 16 (2026-08-19) — the `arb_` block. Internally clean; its UNIVERSE is the story.
    ("arm", "arb_n_edges"): "Edges this arm carries in the LEARNED MODEL's own readouts. ⛔⛔ **DOES NOT "
        "RECONCILE WITH THE EDGE CARD, BY DESIGN AND BY VINTAGE — do not treat a mismatch as a bug.** The "
        "whole `arb_` block rolls up **`readouts_arm_edges.tsv` (2026-08-03)**, not "
        "`realization/edge_card.tsv`. Measured: Σ`arb_n_edges` = **5,792** vs **5,649** edge-card rows, "
        "disagreeing on **58 arms** (miR-124-3p 108 vs 87). ⭐ **The cause is a GENE-universe gap, not an "
        "edge filter: the 158 extra pairs sit on 129 genes that have ZERO rows on the edge card** — and all "
        "129 are on the GENE card. ⇒ `arb_` is the wider, model-native universe; the edge card is the "
        "narrower delivered one. ✅ Internally the block is consistent: threshold nesting holds on **0 "
        "violations**, every `_ident` count is a proper subset, and `arb_frac_identified` = "
        "`arb_n_identified / arb_n_edges` on **100%** of rows.",
    ("arm", "arb_frac_identified"): "`arb_n_identified / arb_n_edges` (verified exact on 100% of rows). "
        "⚠ `identified` is `|z| > 2` on the **UNCALIBRATED** posterior SD, which over-admits (24.8% vs "
        "`cal_identified`'s 19.9%) ⇒ this fraction inherits the loose bar. ⚠ Its denominator is the "
        "model-native universe, not the edge card's — see `arb_n_edges`.",
    ("arm", "arb_max_abs_beta"): "Largest |β| among the arm's edges. ⛔ **Identical to `arb_mean_abs_beta` "
        "on the 229 arms (9.3%) that carry exactly ONE edge** — a max over one unit is that unit (axiom 8's "
        "degeneracy trap). Gate on `arb_n_edges > 1` before comparing the two.",
    ("arm", "arb_mean_abs_beta"): "Mean |β| over the arm's edges. See `arb_max_abs_beta` for the one-edge "
        "degeneracy (9.3% of arms, where the two are equal by construction).",
    # ── UNIT 15 (2026-08-19) — the `fst_` block, audited against itself. A defect I shipped.
    ("arm", "fst_d_share_NAT_TUM"): "Change in the arm's share of its seed family's linear dose, matched "
        "NAT → tumour, **re-normalised over the members measured in BOTH states**. ⛔⛔ **THIS COLUMN WAS "
        "WRONG UNTIL 2026-08-19 AND THE ERROR WAS THE WHOLE SIGNAL.** The three `fst_share_*` come from "
        "three arm-card abundance columns with very different fill (TUM `arm_med_rpm` **19.8%**, NAT "
        "`hly_nat_median` **74.0%**, HLY `hly_gtex_median` **17.0%**), so the raw difference subtracted two "
        "fractions with DIFFERENT denominators — the TUM denominator held a mean **0.66** members against "
        "NAT's **1.88**, and only **21.2%** of families matched. It read **mean +0.156, 67.6% positive, "
        "p=1.2e-12** — *'arms gain family share in tumour'* — which was **entirely the coverage**. On a "
        "common support: **mean +0.0000, 52.6% positive, p=0.75.** ⚠⚠ **AND THE MEAN IS NOW STRUCTURALLY "
        "ZERO** (shares over a fixed support sum to 1 in both states, so the deltas sum to 0 within family) "
        "⇒ **never read the level; read the DISPERSION**: |Δ| median **0.0449**, 90th pct **0.2234**, max "
        "**0.6273**; **23.0%** of arms shift >0.10 of family share, 12 shift >0.25. ⛔ **AND 86 of 135 rows "
        "sit in 2-member families, where the two arms' deltas are EXACT MIRROR IMAGES by construction** — "
        "the row count overstates the independent information ~2× in the dominant stratum. Largest genuine "
        "switches: miR-190b vs miR-190a-5p (±0.627), miR-96-5p vs miR-1271-5p (±0.557), miR-196a vs "
        "miR-196b-5p (±0.435). NaN where <2 members are measured in both states (a one-member common set "
        "forces both shares to 1.0 and the delta to a fake 0).",
    ("arm", "fst_d_share_HLY_TUM"): "As `fst_d_share_NAT_TUM` but GTEx-healthy → tumour, same "
        "common-support correction. n=68 (HLY is the sparsest leg at 17.0%). Post-fix **mean −0.0000, "
        "47.1% positive, p=0.93**. ⚠ Cross-platform: it inherits the gauge offset of MH-249, so read it as "
        "a within-family CONTRAST only. Same structural-zero and mirror-image caveats as its NAT twin.",
    ("arm", "fst_n_common_NAT_TUM"): "How many of the family's arms are measured in BOTH NAT and tumour — "
        "the DENOMINATOR behind `fst_d_share_NAT_TUM`, shipped because the column was wrong for want of "
        "it. **Mode 2 (86 of 135 rows)**; at exactly 2 the family's two deltas are mirror images.",
    ("arm", "fst_n_common_HLY_TUM"): "Members measured in BOTH GTEx-healthy and tumour — the denominator "
        "behind `fst_d_share_HLY_TUM`. See `fst_n_common_NAT_TUM`.",
    ("arm", "fst_share_TUM"): "The arm's fraction of its seed family's total linear dose in tumour, from "
        "`arm_med_rpm` (fill **19.8%**). ⛔ **Median 1.000 — because `fst_n_members` is 1 for most arms and "
        "a singleton's share is 1.0 BY CONSTRUCTION** (axiom 8's degeneracy trap). **Mask on "
        "`fst_n_members > 1` before any family comparison.** ⚠ Its denominator sums only over MEASURED "
        "members, which is why the deltas need a common support — see `fst_d_share_NAT_TUM`.",
    ("arm", "fst_share_NAT"): "Share of family dose in matched NAT, from `hly_nat_median` (fill **74.0%** — "
        "the best-covered leg, 3.7× TUM). Same singleton degeneracy: median 1.000, mask on "
        "`fst_n_members > 1`.",
    ("arm", "fst_share_HLY"): "Share of family dose in GTEx-healthy, from `hly_gtex_median` (fill "
        "**17.0%** — the sparsest leg). Same singleton degeneracy; and cross-platform, so contrast only.",
    # ── UNIT 14 (2026-08-19) — the BARE-name columns, the six that collide across cards.
    ("edge", "role"): "⛔ **NOT the arm's role — the GENE's cancer role** (`oncogene` / `tsg` / "
                      "`oncogene/tsg` / `unknown`), from `gene_roles.load_gene_roles`. Rung is **gene**, "
                      "verified constant within gene (max 1 distinct value). ⚠⚠ **COVERAGE IS 11.1%: "
                      "`unknown` on 1,262 of 1,420 genes.** The source is a **232-gene dictionary "
                      "hardcoded as Python sets** (114 oncogene / 106 tsg / 12 both). ⚠ Per "
                      "`gene_roles.py` it IS a **COSMIC-CGC / OncoKB consensus**, breast-cancer-"
                      "prioritised and deliberately conservative (unlisted ⇒ sign 0, never a guessed "
                      "polarity) — so the defect is not that the source is arbitrary, it is that the "
                      "source is **FROZEN, UNVERSIONED and SUBSET**: no census version, no download "
                      "date, no refresh path. A stratification by `role` therefore runs on 158 of "
                      "1,420 genes. Ingesting the versioned CGC is queued on `PROGRAM_FORWARD_BOARD.md`. "
                      "Its purpose is `malignancy_sign`, i.e. scoring "
                      "the DIRECTION of repressing a node — not a general cancer-gene annotation.",
    ("gene_family", "n_samples"): "Samples the family cell was fit on — **CONSTANT 1,040** across the card. "
                                  "Renamed from bare `n` (2026-08-19): unit C applied that rename through "
                                  "`realization._normalise_edge_names`, which funnels the EDGE card only, "
                                  "so this card kept `n` for two weeks after I recorded the change as done "
                                  "(doctrine §4.5 — this card has no `_finish_card` call site). "
                                  "`card_ladders.normalise_family_card()` is now its funnel.",
    ("gene_family", "z"): "`beta / beta_sd` for the family cell, on the **UNCALIBRATED** posterior SD. "
                          "⛔ The per-edge theoretical null is **3-4× too narrow** (memory "
                          "`edge-null-is-miscalibrated`) ⇒ use `cal_z`. Median 1.261.",
    ("gene_family", "identified"): "`|z| > 2` on the UNCALIBRATED SD — **27.1%** of cells. `cal_identified` "
                                   "is the honest version. ⚠ Not a discovery: per-edge discovery is EMPTY "
                                   "under the empirical FDR (MH-154/155); this flags fit precision only.",

    ("edge", "n_fam"): "How many families compete at this edge's GENE. Gene-rung — it repeats across the "
                       "gene's edges by construction. beta is mathematically inert where it is 1.",
    ("edge", "w_max"): "Maximum curated EVIDENCE weight over the gene's regulators. Gene-rung repeat. "
                       "\u26d4 Not a beta and not a dose share.",
    ("edge", "n_arm_in_cell"): "How many arms were collapsed into this edge's family cell. >1 for 20.3% of "
                               "edges \u2014 outside those, the arm rung IS the family rung by construction.",
    ("edge", "kd_affinity_pct"): "Is this gene among the arm's strongest predicted targets \u2014 a "
                                 "WITHIN-ARM percentile, **higher = stronger target**. Climbs with coupling "
                                 "(rho=\u22120.090 within seeded edges, p=1.3e\u221207), which is one rung of the "
                                 "convergent-evidence ladder. \u26d4 Not an absolute affinity and not "
                                 "comparable across arms of different promiscuity.",
    ("edge", "echim_any"): "Does this edge carry ANY chimeric duplex evidence. \u26a0\u26a0 It is never False \u2014 "
                           "only True (1,065 edges) or NaN. A blank means 'no chimeric record', which "
                           "conflates *not scanned* with *scanned and absent*; test `!= 1`, never `== 0`.",
    ("edge", "arm_med_rpm"): "Median linear RPM of this ARM across the cohort. Arm-rung, so it REPEATS "
                             "across the arm's genes. Fill 0.686 — the gap is arms with no expression "
                             "record, not zeros.",
    ("edge", "arm_iqr"): "Interquartile range of the arm's abundance — its DYNAMIC RANGE. ⭐ Dispersion is "
                         "usually the informative half against level (axiom 8), though it did NOT separate "
                         "diluting from improving mates in MH-248 (q=0.47).",
    ("edge", "arm_pct_floor"): "Percent of samples in which the arm sits above the detection floor "
                               "(median 99). Feeds `arm_spiker` together with `arm_iqr`.",
    ("", "arm_lfc_"): "Log fold-change of the ARM's abundance between two states. ⚠⚠ **ONLY "
                      "`arm_lfc_NAT_TUM` IS SAME-PLATFORM (TCGA→TCGA). Every `HLY` leg crosses "
                      "GTEx→TCGA**, and that boundary carries a measured artifact — in the `fst_` work "
                      "near-total share flips were 5–8× more common across it than within TCGA. ⛔⛔ **AND "
                      "THE TWO HLY→TUM VARIANTS DISAGREE MATERIALLY: `_QN` vs `_raw` correlate at only "
                      "spearman +0.544, median |diff| 1.49, and they DISAGREE ON SIGN FOR 23.9% OF "
                      "EDGES.** Which leg you pick changes the direction for ~1 edge in 4, so state which "
                      "one you used. ⚠ The QN bridge is CORRECT and settled but carries a standing "
                      "trust-weighting — prefer rank measures over QN magnitudes for cross-platform work.",
    ("edge", "hly_leg_concordant"): "Do the two HLY→TUM legs (`_QN` and `_raw`) at least AGREE ON SIGN. "
                                    "**Concordant 76.1%, DISCORDANT 23.9%** of the 3,433 edges carrying "
                                    "both. ⛔ **NOT a correctness flag** — agreement does not make either "
                                    "leg right; it is a DISAGREEMENT flag, which is what the data "
                                    "supports. Gate any cross-state claim on it rather than resting on an "
                                    "arbitrary choice of leg, and remember both legs cross GTEx→TCGA "
                                    "while only `arm_lfc_NAT_TUM` is same-platform.",
    ("edge", "arm_credit_share"): "This arm's share of its (gene, seed_family) cell's credit. ✅ Verified "
                                  "to sum to 1.0000 within every one of the 467 multi-arm cells, and it "
                                  "genuinely varies PER ARM inside the cell (467/467).",
    ("edge", "arm_id_status"): "How the arm's identity within its cell was resolved: `resolved_by_dose` "
                               "2,344 · `singleton` 2,260 · `unvalidated_balanced` 1,045. ⚠ **Fill is 1.0 "
                               "while its arm-in-family siblings are 0.203** — because `singleton` is a "
                               "VALID answer for a one-arm cell (95.8% of singletons are), where the "
                               "siblings are NaN. Same rung, different domain. ⛔ `unvalidated_balanced` "
                               "means dose CANNOT resolve which arm acts — it needs a breast chimeric-eCLIP "
                               "or a per-arm knockdown, not more modelling.",
    ("edge", "arm_dbeta"): "max−min β across the cell's member arms. ⚠⚠ **A CELL quantity REPEATED across "
                           "the cell's arms — verified constant within all 467 multi-arm cells.** Reading "
                           "it as a per-arm fact double-counts; the per-arm quantities are `beta_arm` / "
                           "`z_arm` / `arm_credit_share`, which do vary (467/467).",
    ("edge", "arm_sep_z"): "`arm_dbeta` in units of the cell's own pooled posterior SD — can these "
                           "same-seed arms be told apart at all. ⚠⚠ **Also a CELL quantity, constant "
                           "within all 467 cells.** Median 1.59.",
    ("edge", "arm_resolvable"): "Are the cell's arms separable AND does splitting them help out-of-fold. "
                                "True for 31.7% of edges in multi-arm cells. ⚠⚠ **A CELL verdict, constant "
                                "within all 467 cells** — not a property of the individual arm.",
    ("edge", "z_arm"): "beta_arm / sd_arm \u2014 the ARM-resolved z inside the family cell, as opposed to the "
                       "family-level `z`.",
})


def describe(col: str, card: str = "") -> str | None:
    """Description for one column, CARD-SCOPED. Resolution order, most specific first:

        COLUMNS[(card, col)]   this card's exact override
        COLUMNS[("", col)]     an ALL-CARDS exact entry            <- added 2026-08-19
        BLOCKS[(card, prefix)] this card's prefix block, longest prefix wins
        BLOCKS[("", prefix)]   the all-cards prefix block
        None                   nothing recorded — never a guess

    ⭐ **THE `("", col)` STEP EXISTS BECAUSE ITS ABSENCE WAS A TRAP (user-directed).** `BLOCKS` is
    PREFIX-matched and `COLUMNS` is EXACT-matched, so an `("", name)` key was meaningful in one table and
    silently inert in the other: written into `COLUMNS`, it could never be found, because lookup asks for
    `("arm", name)` / `("gene", name)`. **That mistake was made three times in one session**
    (`detection`/`spiker`, the `n`/`p_fam` batch, `disc_gold_families`) and each time the only symptom was
    a coverage number one or two short of 100%. Accepting `("", col)` here makes an all-cards exact entry
    mean what it looks like it means, in BOTH tables.
    ⚠ A per-card entry still wins over an all-cards one — that is the point of card-scoping (MH-220).
    """
    for key in ((card, col), ("", col)):
        if key in COLUMNS:
            return COLUMNS[key]
    for key_card in (card, ""):
        hits = [(p, t) for (c, p), t in BLOCKS.items() if c == key_card and col.startswith(p)]
        if hits:
            return max(hits, key=lambda kv: len(kv[0]))[1]      # longest prefix wins
    return None


def annotate(registry: pd.DataFrame) -> pd.DataFrame:
    """Add a `description` column to a card_registry-shaped frame."""
    out = registry.copy()
    out["description"] = [describe(c, k) for k, c in zip(out["card"], out["column"])]
    return out


def coverage() -> pd.DataFrame:
    reg = pd.read_csv(OUT / "card_registry.tsv", sep="\t", dtype=str).fillna("")
    return annotate(reg)


def duplicate_keys() -> list[tuple[str, str, int, int]]:
    """⛔ THE DUPLICATE-KEY GUARD — a dict literal SILENTLY accepts a repeated key and the LAST one wins.

    This cost a real error on 2026-08-19: four freshly-written `gene_family` entries were dead on arrival
    because superseded stubs for the same keys sat LATER in `COLUMNS`. Nothing complained — the module
    imported, `--emit` reported 100% coverage, and the only symptom was that the new text did not appear
    in `--col`. **A silent overwrite in a documentation table is worse than a missing entry**, because the
    stale text still reads as authoritative.

    Python gives no runtime hook for this (the duplicate is resolved at literal construction), so the
    check is on the SOURCE: parse this file's AST and report any key repeated inside `BLOCKS` or
    `COLUMNS`, with both line numbers. Runs in `main()`, so it can never recur unnoticed.
    """
    import ast

    src = pathlib.Path(__file__).read_text()
    tree = ast.parse(src)
    dups: list[tuple[str, str, int, int]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.AnnAssign | ast.Assign):
            continue
        targets = [node.target] if isinstance(node, ast.AnnAssign) else node.targets
        names = [t.id for t in targets if isinstance(t, ast.Name)]
        if not ({"BLOCKS", "COLUMNS"} & set(names)) or not isinstance(node.value, ast.Dict):
            continue
        seen: dict[str, int] = {}
        for k in node.value.keys:
            try:
                lit = ast.literal_eval(k)
            except Exception:
                continue
            tag = repr(lit)
            if tag in seen:
                dups.append((names[0], tag, seen[tag], k.lineno))
            else:
                seen[tag] = k.lineno
    return dups


def main() -> None:
    args = sys.argv[1:]

    dups = duplicate_keys()
    if dups:
        print(f"⛔ {len(dups)} DUPLICATE KEY(S) — the LATER one silently wins and the earlier text is DEAD:")
        for where, key, first, second in dups:
            print(f"   {where}[{key}]  defined at line {first}, OVERWRITTEN at line {second}")
        print("   fix: delete the superseded entry, then re-run.\n")
    if "--col" in args:
        col = args[args.index("--col") + 1]
        reg = pd.read_csv(OUT / "card_registry.tsv", sep="\t", dtype=str).fillna("")
        hit = reg[reg["column"] == col]
        if hit.empty:
            print(f"{col!r} is on no card."); return
        for r in hit.itertuples():
            print(f"\n[{r.card}] {col}   rung={r.rung or '-'}  agg_of={r.agg_of or '-'}")
            print(f"  {describe(col, r.card) or '(NO DESCRIPTION RECORDED)'}")
        return

    d = coverage()
    have = d["description"].notna()
    print(f"described: {have.sum()}/{len(d)} ({have.mean():.1%})")
    print("\nby card:")
    for card, sub in d.groupby("card"):
        h = sub["description"].notna()
        print(f"  {card:12s} {h.sum():3d}/{len(sub):3d}  ({h.mean():5.1%})")
    miss = d[~have]
    if len(miss):
        print(f"\n{len(miss)} columns with NO description — fill from the producing module, do not guess:")
        for card, sub in miss.groupby("card"):
            print(f"  [{card}] " + ", ".join(sub["column"].head(18)) + (" …" if len(sub) > 18 else ""))
    if "--emit" in args:
        d.to_csv(OUT / "card_glossary.tsv", sep="\t", index=False)
        print(f"\nwrote {OUT/'card_glossary.tsv'}")


if __name__ == "__main__":
    main()
