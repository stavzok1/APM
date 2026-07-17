"""Assay-method taxonomy for the method-centric evidence ledger (Design §Decision E).

Maps raw experimental-method strings (miRTarBase `Experiments`, TarBase `experimental_method`) to a
small set of canonical **assay classes**, each with an assay-directness weight (does the assay show
*repression* or only *proximity*?). The ledger counts DISTINCT PMIDs per (edge, assay_class); the fused
weight is Σ_class w_class · log1p(distinct_PMID_count). This is method-centric (the axis is the assay,
pooling all databases), not source-centric.
"""
from __future__ import annotations

# raw method substring (lowercased) -> canonical assay class. First containing-match wins (ordered).
_RULES: list[tuple[str, str]] = [
    ("luciferase", "reporter"), ("gfp reporter", "reporter"), ("reporter assay", "reporter"),
    ("qclash", "chimeric"), ("clash", "chimeric"), ("chimeric", "chimeric"), ("impact-seq", "chimeric"),
    ("degradome", "degradome"), ("pare", "degradome"),
    ("hits-clip", "ago_clip"), ("par-clip", "ago_clip"), ("clip-seq", "ago_clip"), ("iclip", "ago_clip"),
    ("clip", "ago_clip"), ("ago-ip", "ago_clip"), ("rip-seq", "ago_clip"), ("immunoprecipita", "ago_clip"),
    ("biotin", "ago_clip"),
    ("western", "western"), ("immunoblot", "western"),
    ("psilac", "proteomics"), ("silac", "proteomics"), ("proteomic", "proteomics"), ("2d-dige", "proteomics"),
    ("qrt-pcr", "qpcr_rna"), ("qpcr", "qpcr_rna"), ("rna-seq", "qpcr_rna"), ("rnaseq", "qpcr_rna"),
    ("srna-seq", "qpcr_rna"), ("microarray", "qpcr_rna"), ("rpf-seq", "qpcr_rna"), ("northern", "qpcr_rna"),
    ("sequencing", "qpcr_rna"),
]

# assay-directness weights (reporter/protein > RNA > binding); 'other' excluded from the fused weight.
CLASS_WEIGHT: dict[str, float] = {
    "reporter": 3.0,      # luciferase / seed-mutant: direct functional readout
    "western": 2.5,       # protein blot: functional, one step removed
    "proteomics": 2.5,    # (p)SILAC / MS: protein-level functional
    "qpcr_rna": 1.5,      # qPCR / RNA-seq / microarray: functional but indirect
    "degradome": 1.0,     # cleavage evidence
    "chimeric": 1.0,      # direct miRNA-target duplex (physical, arm-assigning; repro-contested)
    "ago_clip": 0.5,      # AGO occupancy: proximity, not demonstrated repression
    "other": 0.0,         # IHC/IF/ELISA/flow/ISH/review/ChIP: kept for provenance, no weight
}

# TASK-MATCHED (mRNA) weights — DATA-GROUNDED from the transfection calibration (learned.evidence.calibrate,
# McGeary 2019 GSE140217; 8021 within-family curated edges). The coarse decision (direct-functional >
# binding-only) is confirmed non-circularly: functional-evidence ρ=+0.21 vs binding ρ=+0.08 (~3×). But on an
# *mRNA* readout the fine order differs from the protein-informed hand-set in an interpretable way — RNA-level
# evidence (qpcr_rna) is the single best predictor and protein-only (western) under-predicts mRNA change.
# Since the model's Y is TCGA mRNA, these are the task-matched weights. Provenance: learned NNLS coefficients
# (within-family, potency-de-confounded) normalized to qpcr=2.5, with functional classes floored at 1.0 so
# NNLS collinearity zeros (western/degradome) don't discard genuine functional evidence. Not the default;
# selectable via w_prior_source='ledger_mrna' for the robustness check that the gate doesn't depend on it.
CLASS_WEIGHT_MRNA: dict[str, float] = {
    "reporter": 1.86, "western": 1.00, "proteomics": 1.96, "qpcr_rna": 2.50,
    "degradome": 1.00, "chimeric": 0.36, "ago_clip": 0.62, "other": 0.0,
}


def classify(method: str) -> str:
    m = str(method).strip().lower()
    for sub, cls in _RULES:
        if sub in m:
            return cls
    return "other"
