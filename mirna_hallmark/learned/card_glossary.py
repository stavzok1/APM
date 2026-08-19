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
                  "row it REPEATS by construction.",
    ("", "gctx_"): "Genomic context of the arm's precursor locus (`genomic_context.classify_he_arms`): "
                   "host gene, sense/antisense orientation, intron-vs-exon, promoter class. Strand-aware.",
    ("", "comp_"): "Cell-composition drivers of this gene's coupling (`compartment_*`): which deconvolved "
                   "fraction explains the correlation, its retention under adjustment, and its share. "
                   "⚠ This is DRIVER-SHARE, not the gene's own expression LOADING on a compartment — "
                   "`gene_axes` measured only the former to predict anything.",
    ("", "cov_"): "Coverage flags — was this row ever SCANNED by the named source. ⛔ Read before reading "
                  "the block it gates: a blank in a gated block means `— not scanned`, which is not a zero.",

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
    ("", "arm_"): "Arm-resolved fit quantities inside a family cell — the arm's own beta, its separation "
                  "z from its same-seed mates, and whether it is `arm_resolvable`. Defined only where a "
                  "family cell holds more than one arm (20.3% of edges).",
    ("", "term_"): "Decomposition of the acquisition score into ABUNDANCE, WIRING and INTERACTION terms.",
    ("", "coupling_"): "Partial-Spearman coupling of arm dose against target mRNA in one state, "
                       "conditioned on C (composition + target CN + proliferation), with its p and z. "
                       "⚠ The per-edge null is 3–4x too narrow — read these as a distribution, not as "
                       "per-edge significance.",
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
                      "slab has a strictly positive mean and cannot zero an un-informed regulator.",
    ("gene_family", "beta"): "FAMILY-rung coupling coefficient — `readouts.run(level='family')`, the "
                             "same-seed collapse APPLIED. Same estimator as the edge card's `beta`, "
                             "DIFFERENT UNIT. Do not place the two side by side without saying so.",
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
    ("gene", "n_arms"): "Arms in the gene's DESIGN (the model's own count, from `gene_atlas`). "
                        "⚠⚠ **NOT the same as `n_regulators`** — the two come from different sources and "
                        "**differ on 306 of 1,409 genes, in BOTH directions** (ABCA1 13 vs 12, ABCC1 7 vs 8). "
                        "Both read as 'how many regulators'; only this one matches the design the model fit.",
    ("gene", "n_regulators"): "Regulator count from the CO-REPRESSION lane "
                              "(`tissue_reference/mirna_comovement/gene_corepression.tsv`), not from the "
                              "model design. ⚠⚠ **Differs from `n_arms` on 306 of 1,409 genes in both "
                              "directions** — a different pipeline counted a different universe. Use "
                              "`n_arms` for anything about the FIT; use this only alongside the other "
                              "co-repression columns it arrived with (`gene_repression_class`, "
                              "`gene_net_repressed_tumor`, `rho_gene_pressure_tumor`, `delta_tumor_nat`).",
    ("gene", "n_dense_included"): "Families entering the dense (π≡1) readout. ⛔ **Bit-identical to `n_fam` "
                                  "on every row** — verified 2026-08-19. Redundant; do not treat as an "
                                  "independent axis (it silently produced duplicate scan results once).",
    ("gene", "n_identified"): "Families with |z| > 2 on the UNCALIBRATED SD. ⚠ **691 genes (44.6%) have "
                              "ZERO** — the identifiability ceiling is the dominant fact about this column, "
                              "not a data gap. Under the CALIBRATED width it would be lower still.",
    ("gene", "n_discovered"): "Families passing the evidence-π discovery readout; **>0 for 725 genes "
                              "(46.8%)**. ⛔⛔ **This is NOT a count of discoveries.** Per-edge and "
                              "per-family discovery are EMPTY under the honest empirical FDR — the "
                              "defensible deliverable is a convergent-evidence QUEUE (157 edges / 11 "
                              "families). Reading 46.8% of genes as 'having discoveries' inverts the axis.",
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
    ("edge", "arm_resolvable"): "Is this arm statistically separable from its same-seed mates in this "
                                "gene's fit. True for 31.7% of edges in multi-arm cells, rising 20% -> 48% "
                                "from 2-arm to 3+-arm cells.",
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
    ("edge", "n_samples"): "Number of samples the edge's coupling was fit on.",
    ("gene_family", "n"): "Number of samples the (gene, family) cell was fit on.",
    ("_retired_edge", "p_fam"): "\u26d4\u26d4 **NOT A P-VALUE.** `p` here is the model's DESIGN DIMENSION \u2014 the NUMBER "
                       "of family predictors in this gene's fit (integers 1..90). The name reads as a "
                       "p-value and is a genuine collision; verified 2026-08-19 against the values.",
    ("gene_family", "p_fam"): "\u26d4\u26d4 **NOT A P-VALUE** \u2014 the number of family predictors in the gene's "
                              "fit. **Bit-identical to `n_fam` on this card** (integers 1..12), so it is also "
                              "a redundant column.",
    ("edge", "z"): "beta / beta_sd on the UNCALIBRATED posterior SD. Use `cal_z` for the honest width.",
    ("gene_family", "z"): "beta / beta_sd for the family cell, uncalibrated.",
    ("edge", "identified"): "|z| > 2 on the UNCALIBRATED SD (24.8%). `cal_identified` is the honest "
                            "version (19.9%).",
    ("edge", "pip_dense"): "\u26a0 CONSTANT by construction \u2014 the coupling readout is called with \u03c0\u22611, so "
                           "every edge has the same value. Carries no information; `pip_discovery` is the "
                           "one that varies. (Bit-identical to `net_pressure` on this card.)",
    ("gene_family", "pip_dense"): "\u26a0 CONSTANT by construction (\u03c0\u22611 dense readout). Bit-identical to "
                                  "`net_pressure`.",
    ("edge", "net_pressure"): "\u26a0 CONSTANT on this card and bit-identical to `pip_dense`.",
    ("gene_family", "net_pressure"): "\u26a0 CONSTANT on this card and bit-identical to `pip_dense`.",
    ("edge", "n_samples"): "Samples the edge was fit on \u2014 CONSTANT (1,040) across the card.",
    ("gene_family", "n"): "Samples the cell was fit on \u2014 CONSTANT (1,040) across the card.",
    ("gene", "n_dense_included"): "Families entering the dense (\u03c0\u22611) readout. \u26a0 **Bit-identical to "
                                  "`n_fam`** \u2014 redundant; do not treat as an independent axis.",
    ("gene_family", "identified"): "|z| > 2 for the family cell, uncalibrated.",
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
    ("gene", "n_arms"): "How many miRNA arms regulate this gene (median 2).",
    ("gene", "n_regulators"): "Regulator count used by the realization lane.",
    ("gene", "n_dense_included"): "Families entering the dense (pi=1) coupling readout.",
    ("gene", "n_identified"): "How many of the gene's families are identified (|z|>2). Median 1.",
    ("gene", "n_discovered"): "Families passing the evidence-pi discovery readout. \u26d4 Per-edge and "
                              "per-family discovery are EMPTY under the honest FDR \u2014 read this as a "
                              "ranked queue, never as a set of findings.",
    ("gene", "n_composition_explained"): "Edges whose coupling does NOT survive the composition block.",
    ("gene", "n_cell_intrinsic"): "Edges whose coupling DOES survive composition adjustment.",
    ("gene", "n_hallmark_sets"): "How many MSigDB Hallmark programs this gene belongs to.",
    ("edge", "n_HLY_meas"): "How many of the gene's arms have a MEASURED (not imputed, not "
                            "multi-mapping-collapsed) healthy level. Read before any `share_HLY`/`rank_HLY`.",
    ("edge", "identity_allocated"): "Family identity with phantom minor arms forced to 0 \u2014 the allocation "
                                    "actually used downstream.",
    ("edge", "arm_id_status"): "Whether the arm's identity could be resolved within its family cell.",
    ("edge", "role"): "The arm's role in this gene's regulation, as assigned by the realization lane.",
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
    ("edge", "z_arm"): "beta_arm / sd_arm \u2014 the ARM-resolved z inside the family cell, as opposed to the "
                       "family-level `z`.",
})


def describe(col: str, card: str = "") -> str | None:
    """Description for one column, CARD-SCOPED. Exact override first, then the card's own prefix
    entry, then the global prefix entry. Returns None when nothing is recorded — never a guess."""
    if (card, col) in COLUMNS:
        return COLUMNS[(card, col)]
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
