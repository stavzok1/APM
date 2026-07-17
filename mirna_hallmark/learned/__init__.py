"""Learned miRNA→gene regulatory model (see docs/LEARNED_MODEL_DESIGN_RESPONSE.md).

This package is the NEW model. It **imports the spine in place** (config, data_loaders,
coupling_inference, evidence_scoring, …) and compares against the frozen BASELINE heuristic
(pressure_engine, hallmark_interaction). Nothing here rewrites the spine.

Layout (BUILD_PLAN §3):
    data.py         assemble aligned (Y, X, C, w_prior) for a gene/program from the spine
    regression.py   gene-focused non-negative adaptive-lasso fit (Phase-1 estimator)
    mvp.py          Phase-1 driver: fit hub panel + OOF gate vs raw abundance (Design §6 Phase 1)
    priors.py       π/μ/τ from the PMID-deduped method-centric ledger + K_D   (stub → Design §E)
    occupancy.py    a/(a+κ) link + threshold gauge + shared free-AGO pool      (stub → Design §C/§G)
    families.py     seed-family partial pooling (identified estimand=family→gene)(stub → Design §F)
    cooperativity.py site-primed product terms (8–40 nt window)                (stub → Design §G)
    state.py        nested M_tumour = M_healthy + Δ                            (stub → Design §H)
    attribution.py  Shapley identity vs realized magnitude                      (stub → Design §I)
    evidence/       the (m,g,PMID,assay) ledger + method-centric dedup + learned weights (stub → Design §E)
    eval/           thin wrappers over reused coupling/permutation/instrument/OOD (stub → Design §5)
"""
