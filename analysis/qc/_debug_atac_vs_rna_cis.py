#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    from analysis.qc._context import load_context, load_sample_ids
    from analysis.qc._report_utils import ensure_out_dir
    from analysis.qc.run_cross_modal_qc import _qc_atac_vs_rna_cis
    from pipeline.RNA_exp.signatures import GeneSets

    scratch = Path("analysis/trial_max_coverage_cohort/processing_scratch_20260421_192526/cohort_processing_outputs.json")
    ctx = load_context(scratch)
    sids = load_sample_ids(ctx)
    out = ensure_out_dir("tmp_atac_debug")
    gs = GeneSets()
    per_peak, per_gene, summ = _qc_atac_vs_rna_cis(out, sample_ids=sids, genes=list(gs.apm_class_i))
    print("summary:", summ)
    print("per_peak:", per_peak.shape)
    print("per_gene:", per_gene.shape)


if __name__ == "__main__":
    main()

