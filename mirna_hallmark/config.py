"""Configuration for the miRNA x Hallmark subproject.

Single place for input paths (mostly delegated to ``pipeline.config.PATHS``),
output locations, stratum definitions, evidence cutoffs, and the AGO/RISC gate
parameters. Importing this module also ensures the repo root is on ``sys.path``
so the subproject can be run with ``python -m mirna_hallmark.<module>`` from the
repository root.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

# --- repo root on sys.path (so `pipeline` / `analysis` import cleanly) ------- #
# mirna_hallmark/config.py -> parent = subproject dir -> parent.parent = repo root
REPO_ROOT = Path(__file__).resolve().parent.parent
SUBPROJECT_DIR = REPO_ROOT / "mirna_hallmark"

import sys  # noqa: E402

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pipeline.config import PATHS  # noqa: E402

# --------------------------------------------------------------------------- #
# Output tree (self-contained under the subproject)
# --------------------------------------------------------------------------- #
OUTPUT_ROOT = SUBPROJECT_DIR / "output"
EDGES_DIR = OUTPUT_ROOT / "edges"
MIRNA_LOCUS_CNV_DIR = OUTPUT_ROOT / "mirna_locus_cnv"
AGO_GATE_DIR = OUTPUT_ROOT / "ago_gate"
STRATUM_DIR = OUTPUT_ROOT / "stratum_characterization"
INTERACTION_DIR = OUTPUT_ROOT / "hallmark_interaction"
CONTRASTS_DIR = OUTPUT_ROOT / "subtype_contrasts"
VISIBILITY_ARCHETYPE_DIR = OUTPUT_ROOT / "visibility_archetype_contrasts"
PAM50_GENE_RESOLUTION_DIR = OUTPUT_ROOT / "pam50_gene_resolution"
FIGURES_DIR = OUTPUT_ROOT / "figures"
HYBRID_PRESSURE_DIR = OUTPUT_ROOT / "hybrid_pressure"
TARGET_COMBINED_ANTICORR_DIR = OUTPUT_ROOT / "target_combined_anticorr"
# End-to-end CN -> regulatory-network attribution (CN-driven vs CN-residual
# pressure) + combined multi-locus CN-load granularity contest.
CN_DOSAGE_ATTRIBUTION_DIR = OUTPUT_ROOT / "cn_dosage_attribution"
CCLE_CN_CONCORDANCE_DIR = OUTPUT_ROOT / "ccle_mirna_cn_concordance"
# TCGA tissue + NAT + GTEx cross-state reference (not plasma EV).
TISSUE_REFERENCE_DIR = OUTPUT_ROOT / "tissue_reference"
# Optional translational lane: circulating / exosome miRNA GEO cohorts only.
# Not part of ``run_all``; kept in this subproject for shared pressure helpers.
PLASMA_EV_DIR = OUTPUT_ROOT / "plasma_ev"
LOGS_DIR = OUTPUT_ROOT / "logs"

# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
HALLMARK_JSON = REPO_ROOT / "annotations" / "GSEA" / "h.all.v2026.1.Hs.txt"

# miRTarBase (raw download + regenerated Hallmark-scoped summary)
MIRTARBASE_CSV = PATHS.mirtarbase_csv
MIRTARBASE_SUMMARY = PATHS.mirtar_interaction_summary_csv
MIR_FAMILY_INFO = PATHS.mir_family_info
# Regenerated for the Hallmark universe (build_edges.py writes here):
MIRTARBASE_HALLMARK_DIR = EDGES_DIR / "mirtarbase"
MIRTAR_HALLMARK_SUMMARY = MIRTARBASE_HALLMARK_DIR / "mirtar_interaction_summary.csv"
HALLMARK_EDGES_TSV = EDGES_DIR / "mirna_hallmark_edges.tsv.gz"

# Expression matrices
MIRNA_EXPRESSION = PATHS.mirna_expression          # log2(RPM+1), MIMAT rows
MIRNA_MATURE_LOCI = PATHS.mirna_mature_loci_csv      # MIMAT -> arm map
MIRNA_PRECURSOR_LOCI = PATHS.mirna_precursor_loci_csv  # MI* precursor loci (CNV)
RNA_EXPRESSION = PATHS.rna_expression              # log2(TPM+1), gene-symbol rows

# CNV (ASCAT3: segments default; gene tables legacy fallback)
CNV_DIR = PATHS.cnv_dir
CNV_MANIFEST = PATHS.cnv_annotations_path
CNV_GENE_TABLES_DIR = PATHS.cnv_gene_tables_dir

# SV (Manta strict set + pipeline mir_hits; 07 is sufficient for miRNA overlap)
SV_OUTPUT_ROOT = PATHS.sv_output_root
SV_FINAL_DIR = SV_OUTPUT_ROOT / "07_final_sv_with_fimo"
MIRNA_LOCUS_SV_DIR = MIRNA_LOCUS_CNV_DIR / "sv_overlap"
# Optional; may be stale vs current SV dir — miRNA overlap computes VAF from Manta reads directly.
SV_DISTRIBUTIONS_FEATURES = REPO_ROOT / "analysis/output/sv_distributions/primary_all/per_sample_sv_features.tsv.gz"

# Clinical / strata
CLINICAL_UNIFIED = PATHS.brca_clinical_immune_unified
TNBC_SUBTYPES = REPO_ROOT / "annotations" / "BRCA_TNBC_Lehmann2021_subtypes.tsv.gz"

# CCLE / DepMap (cell-line CN→expression replication; see ccle_mirna_cn_concordance.py)
CCLE_DATA_DIR = REPO_ROOT / "data" / "CCLE"
CCLE_MIRNA_GCT = CCLE_DATA_DIR / "CCLE_miRNA_20181103.gct"
CCLE_MODEL = CCLE_DATA_DIR / "Model.csv"
CCLE_CN_SEGMENTS = CCLE_DATA_DIR / "OmicsAbsoluteCNSegmentsProfile.csv"
CCLE_OMICS_PROFILES = CCLE_DATA_DIR / "OmicsProfiles.csv"
CCLE_DEFAULT_PROFILES = CCLE_DATA_DIR / "OmicsDefaultModelProfiles.csv"
CCLE_DEPMAP_FILES_INDEX = CCLE_DATA_DIR / "depmap_files.csv"
CCLE_MIRNA_EXPR_PSEUDOCOUNT = 50.0  # Ghandi et al. 2019 CCLE miRNA normalization
CCLE_CONCORDANCE_MIN_N = 15
# CCLE partial Spearman covariates (OLS-residualize then Spearman; see partial_spearman).
# TCGA partial uses CONFOUNDER_NUMERIC (CPE + thornsson_hrd_score).
CCLE_PARTIAL_COVARIATE = "lineage"  # all-stratum default: OncotreeLineage dummies
CCLE_BREAST_PARTIAL_COVARIATE = "pam50"  # breast strata: PAM50 group dummies
# Oncotree lineages with enough CCLE lines for stratified concordance (min n enforced at runtime).
CCLE_LINEAGE_STRATA = (
    "Breast",
    "Lung",
    "Bowel",
    "Ovary/Fallopian Tube",
    "Skin",
    "CNS/Brain",
    "Esophagus/Stomach",
)
# CCLE NanoString uses legacy / combined probe names not matching mirBase mature ids.
CCLE_COMBINED_MIRNA_PROBES: Dict[str, str] = {
    "hsa-miR-17-5p": "hsa-miR-106a+hsa-miR-17",
    "hsa-miR-106a-5p": "hsa-miR-106a+hsa-miR-17",
}
CCLE_FOCUS_8Q_ARMS_PATH = MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_focus_subtype_comparison.tsv"

# Plasma exosome miRNA (GEO GSE270497; independent of TCGA tissue spine)
EV_MIRNA_GSE270497_DIR = PATHS.mirna_ev_gse270497_dir
EV_MIRNA_GSE270497_SAMPLES = PATHS.mirna_ev_gse270497_samples
EV_MIRNA_GSE270497_LOG2_TPM = PATHS.mirna_ev_gse270497_log2_tpm
EV_MIRNA_GSE197020_DIR = PATHS.mirna_ev_gse197020_dir
EV_MIRNA_GSE197020_SAMPLES = PATHS.mirna_ev_gse197020_samples
EV_MIRNA_GSE197020_LOG2_CPM = PATHS.mirna_ev_gse197020_log2_cpm
EV_MIRNA_GSE255660_DIR = PATHS.mirna_ev_gse255660_dir
EV_MIRNA_GSE255660_SAMPLES = PATHS.mirna_ev_gse255660_samples
EV_MIRNA_GSE255660_LOG2_EXPR = PATHS.mirna_ev_gse255660_log2_expr
EV_MIRNA_GSE241784_DIR = PATHS.mirna_ev_gse241784_dir
EV_MIRNA_GSE241784_SAMPLES = PATHS.mirna_ev_gse241784_samples
EV_MIRNA_GSE241784_LOG2_NORM = PATHS.mirna_ev_gse241784_log2_norm
EV_MIRNA_GSE241785_DIR = PATHS.mirna_ev_gse241785_dir
EV_MIRNA_GSE241785_SAMPLES = PATHS.mirna_ev_gse241785_samples
EV_MIRNA_GSE241785_LOG2_NORM = PATHS.mirna_ev_gse241785_log2_norm
EV_MIRNA_GSE301416_DIR = PATHS.mirna_ev_gse301416_dir
EV_MIRNA_GSE301416_SAMPLES = PATHS.mirna_ev_gse301416_samples
EV_MIRNA_GSE301416_LOG2_DCQ = PATHS.mirna_ev_gse301416_log2_dcq

# --------------------------------------------------------------------------- #
# AGO / RISC rate-limiting gate
# --------------------------------------------------------------------------- #
# miRNA repression needs two both-required resources (see METHODS.md §4):
#   1. loading into Argonaute (miRISC core). AGO1-4 all bind mature miRNAs and
#      support canonical seed-mediated repression; AGO2 is the predominant, most
#      abundant, and only slicer-competent Argonaute in most somatic cells, so it
#      carries the majority of the loaded miRNA pool and is up-weighted below.
#   2. effector engagement via TNRC6/GW182 (TNRC6A/B/C), which dock onto AGO and
#      recruit the CCR4-NOT / PAN2-PAN3 deadenylases -> decapping -> decay. AGO
#      occupancy without TNRC6 does not repress, so TNRC6 is the rate-limiting
#      *effector*, modelled as a co-limiting term (Liebig's law of the minimum).
AGO_GENES: Tuple[str, ...] = ("AGO1", "AGO2", "AGO3", "AGO4")
RISC_EXTRA_GENES: Tuple[str, ...] = ("TNRC6A", "TNRC6B", "TNRC6C")

# Relative contribution of each Argonaute to loadable miRISC capacity. A transparent
# abundance/dominance prior (AGO2 carries most of the load; AGO4 is typically lowest),
# NOT a fitted value -- it is one of the constants swept in pressure_constant_sensitivity.
AGO_LOAD_WEIGHTS: Tuple[Tuple[str, float], ...] = (
    ("AGO1", 1.0), ("AGO2", 2.0), ("AGO3", 1.0), ("AGO4", 0.5),
)


@dataclass(frozen=True)
class AgoGateParams:
    """Parameters for the saturating AGO/RISC availability gate.

    Per sample ``s`` the gate is built from two both-required components:

    - ``ago_load(s)``  = AGO-weighted mean per-gene z of AGO1-4 (weights
      ``ago_load_weights``) -- loadable miRISC core capacity.
    - ``effector(s)``  = mean per-gene z of TNRC6A/B/C -- the GW182 effector that
      converts AGO occupancy into deadenylation/decay.

    With ``effector_colimit`` the capacity is the *minimum* of the two (whichever
    resource is scarcer rate-limits repression); otherwise it is their mean (a
    compensatory blend, kept only as a sensitivity variant). Capacity maps to a
    multiplier in ``[gate_min, 1]`` via a logistic, so very-low-capacity samples
    damp miRNA pressure toward ``gate_min`` and high-capacity samples are ~1.

    This is a documented sensitivity layer, NOT a causal model: every analysis
    reports gated AND ungated results.
    """

    include_tnrc6: bool = True   # TNRC6/GW182 effector is biologically required
    effector_colimit: bool = True  # capacity = min(ago_load, effector) (Liebig minimum)
    ago_load_weights: Tuple[Tuple[str, float], ...] = AGO_LOAD_WEIGHTS
    gate_min: float = 0.5     # floor of the multiplicative gate
    gate_k: float = 1.0       # logistic steepness on the capacity z-score
    gate_midpoint: float = 0.0  # capacity z at which gate = (1 + gate_min) / 2


AGO_GATE = AgoGateParams()

# --------------------------------------------------------------------------- #
# Evidence / pressure thresholds
# --------------------------------------------------------------------------- #
# "High-evidence" miRNA->target edge: at least one Functional MTI study AND a
# low-throughput validation (luciferase reporter or protein/western). Falls back
# to an evidence_score floor when the cross-count columns are unavailable.
HIGH_EVIDENCE_MIN_FUNCTIONAL_MTI = 1
HIGH_EVIDENCE_REQUIRE_LOWTHROUGHPUT = True   # reporter or protein
HIGH_EVIDENCE_MIN_SCORE = 3                   # fallback when columns missing

# miRNA-pressure edge filter (mirrors analysis.expression.mirna_target_integration)
PRESSURE_MIN_EVIDENCE = 2
PRESSURE_WEIGHT_MODE = "permissive"           # legacy tiered_study_weights; superseded by scorer below
# Unified S1 spine (dual_spine_comparison 2026-06-22): confidence_logclass + ENCORI depth boost.
PRESSURE_EVIDENCE_SCORER = "confidence_logclass"
PRESSURE_ENCORI_ALPHA = 0.5
PRESSURE_CUTOFF_Z = 0.84                      # ~top 20% under normal approx

# Expression multiplier for pressure (see pressure_engine.py)
# Spine (M0): absolute logrpm anchor + evidence-mass norm (pressure_sensitivity 2026-06-16).
PRESSURE_EXPR_MODE = "softmax_z_logrpm"
# Gene-level *attribution* mode (who carries pressure among a gene's regulators).
# No-z on purpose: under softmax_logrpm every contribution is positive, so the
# realized share is unambiguous (abs == signed) and avoids the z-term zeroing an
# arm at its cohort mean. The z-spine above stays canonical for Hallmark coupling.
PRESSURE_ATTRIBUTION_EXPR_MODE = "softmax_logrpm"
# Healthy-anchored mode: replaces the within-tumor z-term with deviation from the
# GTEx true-healthy median, x(m,s) - median_healthy(m) (log2 units, no SD division).
# Requires a per-arm baseline from `healthy_anchor.gtex_qn_baseline()` (GTEx medians
# quantile-normalised onto the TCGA scale). Arms absent in GTEx are floored to 0, so
# their full expression counts as tumour-acquired. Use this when the question is
# "acquired-vs-healthy pressure"; the z-spine stays canonical for within-cohort coupling.
PRESSURE_HEALTHY_ANCHOR_MODE = "softmax_devhealthy_logrpm"
# Structural-resolution mode (who is the gene's *preferred* regulator, de-confounded
# from absolute abundance). The attribution mode above multiplies the softmax share
# `sm` (relative abundance, once) by `logrpm` (absolute abundance, again), so the
# within-gene resolution share is governed by abundance² and `edge_w` (evidence) is
# swamped. The structural mode drops the second abundance factor: c_struct = edge_w·sm,
# so relative abundance is spent once and evidence regains proportional leverage. Used
# ONLY for the within-gene edge-resolution / decoupling "which arm" axis; the spine
# (coupling) and cross-state magnitude mass stay on the modes above. See FORMULAS §5/§9.
PRESSURE_STRUCT_EXPR_MODE = "softmax"
PRESSURE_TARGET_NORM = "evidence_mass"   # none | degree | evidence_mass | ts_mass | combined_mass
PRESSURE_AGGREGATE = "sum"         # sum | mean (mean weakens signal in sensitivity grid)
# Hybrid primary (M7): relative share×z + miRTar+TS combined mass (IRF1 + Hallmark tradeoff).
PRESSURE_HYBRID_EXPR_MODE = "softmax_z"
PRESSURE_HYBRID_TARGET_NORM = "combined_mass"
# Cohort-median log2(RPM+1) floor — drop arms below this before pressure (TS tracks
# use PRESSURE_TS_ABUNDANCE_FLOOR when stricter).
PRESSURE_ABUNDANCE_FLOOR = 0.0
PRESSURE_TS_ABUNDANCE_FLOOR = 1.0

# Arm expression-DETECTABILITY floor (canonical; see mirna_hallmark/arm_expression.py).
# KEEP an arm if it reaches >= ARM_EXPRESSED_MIN_RPM in at least ONE tumor (max RPM; respects context-specific
# induction); REMOVE only arms that NEVER reach the floor in any tumor (truly silent — cannot occupy RISC or
# show coupling). A noise filter, not a functional verdict. ARM_EXPRESSED_MIN_FRAC is reported (frac_expressed),
# not the cut.
ARM_EXPRESSED_MIN_RPM = 10.0
ARM_EXPRESSED_MIN_FRAC = 0.01

# Hybrid edge modes for extended claims (hybrid uses PRESSURE_HYBRID_* weighting).
PRESSURE_HYBRID_PRIMARY = "M7"       # tiered miRTar + TS fusion + combined_mass
PRESSURE_HYBRID_SENSITIVITY = ("M8", "M11")  # sum-tracks; pure TS capped

# --------------------------------------------------------------------------- #
# Strata + confounders
# --------------------------------------------------------------------------- #
# (clinical column, output layer name)
STRATUM_SPECS: List[Tuple[str, str]] = [
    ("PAM50_final", "PAM50_final"),
    ("tnbc_subtype_4", "tnbc_subtype_4"),
    ("thornsson_immune_subtype", "thornsson_immune_subtype"),
    ("pathologic_stage_collapsed", "pathologic_stage_collapsed"),
]

# Tier-A confounders for inferential tests (per analysis conduct rule).
CONFOUNDER_NUMERIC = ("CPE", "thornsson_hrd_score")
CONFOUNDER_CATEGORICAL = ("PAM50_final",)

# --- Global cohort + batch switches (Normal-like discarded + batch-corrected rerun) ---
# Discard the PAM50 Normal-like class cohort-wide (n=36; normal-tissue admixture; MH-8
# shows it drives spurious cohort-level correlations). Applied centrally in
# data_loaders.load_clinical_strata, so every module inherits the filtered cohort.
EXCLUDE_NORMAL_LIKE = True
NORMAL_LIKE_LABEL = "Normal"

# TCGA-BRCA sequencing-batch adjustment for the partial-correlation spine: add the
# analyte-plate dummies as COVARIATES (not ComBat) to every confounder-adjusted coupling
# test, via tcga_batch.augment_cov. "plate_both" = miRNA-seq + RNA-seq plates.
# kind in {none, tss, plate, mirna_plate, plate_both}.
TCGA_BATCH_KIND = "plate_both"
# Per-state batch for the cross-state (GTEx/NAT/tumor) coupling: tumor/NAT use the analyte
# plate parsed from the TCGA barcode (NAT = the sample-11 aliquot's plate), GTEx uses the
# sample-attribute nucleic-acid (SMNABTCH) + sequencing (SMGEBTCH) batches.
GTEX_SAMPLE_ATTR = REPO_ROOT / "data" / "GTEx" / "GTEx_v10_SampleAttributesDS.txt"

PAM50_ORDER = (("LumA", "LumB", "Her2", "Basal")
               if EXCLUDE_NORMAL_LIKE else ("LumA", "LumB", "Her2", "Basal", "Normal"))

# Subtype groupings for contrast analyses. Normal-like is deliberately SET ASIDE
# (small n + normal-tissue admixture confounds tumor-intrinsic comparisons).
SUBTYPE_GROUPS = {
    "luminal": ("LumA", "LumB"),
    "nonluminal": ("Her2", "Basal"),      # HER2 + Basal grouped by similarity
}
SUBTYPE_SINGLETONS = ("LumA", "LumB", "Her2", "Basal")  # one-vs-rest, Normal excluded
SUBTYPE_EXCLUDE_FROM_CONTRAST = ("Normal",)

# Immune / inflammatory Hallmark programs (for the cold-Basal / hot-Luminal axis).
IMMUNE_HALLMARKS = (
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
)
# APM cold-Basal IRF1 miRNA route (D76): IRF1-low is largely miRNA-explained in
# cold-Basal via these arms. We check for their high-evidence edges here.
IRF1_ROUTE_MIRNAS = ("hsa-miR-23a-3p", "hsa-miR-106b-5p")

# FDR level for multiple testing
FDR_ALPHA = 0.05
CONCORDANCE_TARGET_TOP_N_ARMS = 25


def ensure_output_dirs() -> None:
    """Create the subproject output tree."""
    for d in (
        OUTPUT_ROOT,
        EDGES_DIR,
        MIRTARBASE_HALLMARK_DIR,
        MIRNA_LOCUS_CNV_DIR,
        AGO_GATE_DIR,
        STRATUM_DIR,
        INTERACTION_DIR,
        CONTRASTS_DIR,
        VISIBILITY_ARCHETYPE_DIR,
        PAM50_GENE_RESOLUTION_DIR,
        FIGURES_DIR,
        HYBRID_PRESSURE_DIR,
        CCLE_CN_CONCORDANCE_DIR,
        TISSUE_REFERENCE_DIR,
        PLASMA_EV_DIR,
        LOGS_DIR,
    ):
        d.mkdir(parents=True, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# ⭐ RHO_GATE — THE ONE DENOMINATOR GATE FOR EVERY RETENTION RATIO (MH-257, 2026-08-19).
#
# ⛔ **WHY THIS MOVED HERE.** `retention = rho_adjusted / rho_raw` is computed in SEVEN places. Four gated
# the denominator at 0.05 — each declaring its OWN `RHO_GATE = 0.05` constant with the comment *"axiom 5:
# don't divide by a vanishing raw coupling"* (`cptac_card`, `card_context`, `compartment_card`,
# `eval/dossier`). The other three guarded only against division by EXACTLY zero (`1e-6`, `1e-9`), which is
# not a gate at all: `realization.py` produced `realized_retention` reaching **720.8** and `retention_rho`
# reaching **1169.0**, and `readouts.py` the same for `retention`.
#
# ⇒ the same trap, handled correctly in four places and not in three, inside one codebase — the same shape
# as the floor-corrected `fam_dose_hhi` vs the uncorrected `concentration` (column review unit 10).
# A threshold that is a POLICY must have ONE home, or it is not a policy.
#
# **MEASURED IMPACT** of applying it where it was missing (`realized_retention`, n=1,259):
#   max|x| **720.8 → 4.55** (158×) · rows dropped **281 (22.3%)** · median **+0.610 → +0.625** — i.e. it
#   removes the pathological tail and leaves the bulk alone, which is what a correct gate looks like.
# ⚠ **RETENTION IS NOT BOUNDED TO [0,1]** and gating does not make it so — **41.1%** still sits outside
#   after the gate. Sign flips and amplification are real (MH-253); only the blow-ups are the artifact.
RHO_GATE = 0.05
