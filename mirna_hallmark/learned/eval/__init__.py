"""Evaluation = thin wrappers over REUSED machinery (Design §5, BUILD_PLAN §1c). STUB.

The eval half is largely built — reuse, don't rebuild:
    coupling      → coupling_inference.partial_spearman_batch   (OOF coupling + C block)   [Bar 1/2]
    permutation   → analyses/../eval/coupling_permutation        (shuffled-evidence null)   [Bar 3]
    instrument    → NEW 2SLS on mirna_locus_cnv (CN dose instrument)                        [Bar 4]
    ood_protein   → cptac_validation (protein + protein-vs-mRNA discordance, MoPC-style)    [Bar 5]
    ood_cohort    → buffa_validation / outcome_metabric / SCAN-B                            [Bar 5]
    folds         → held_out_tuning (patient-level CV)
References scored against (Design §0): raw abundance (primary), curated fixed-M, shuffled-evidence null.
Coupling-invariance lemma ⇒ test weightings at the AGGREGATE / predictive level (mvp.oof_gate), not
single-edge coupling.
"""
