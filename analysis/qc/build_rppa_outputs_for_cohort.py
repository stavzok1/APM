#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.qc.context import load_context, load_sample_ids  # noqa: E402
from pipeline.rppa.rppa_main import run_rppa_pipeline  # noqa: E402
from pipeline.rppa.rppa_config import RPPAPathConfig  # noqa: E402
from pipeline.config import PATHS  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description="Build RPPA outputs (panel_scores, blocks, etc.) for the cohort.")
    ap.add_argument("--scratch-json", type=str, required=True)
    ap.add_argument("--out-dir", type=str, default="", help="Override output dir (defaults to ctx.rppa_output_dir).")
    args = ap.parse_args()

    ctx = load_context(Path(args.scratch_json))
    sample_ids = load_sample_ids(ctx)
    if not sample_ids:
        raise SystemExit("No sample IDs found for cohort.")

    out_dir = Path(args.out_dir) if args.out_dir else Path(ctx.rppa_output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Run full RPPA pipeline; it loads per-sample RPPA files under PATHS.working_dir/rppa/samples by default.
    # We set output_dir to the cohort-specific directory so QC can find panel_scores.csv there.
    run_rppa_pipeline(
        sample_dir=RPPAPathConfig.sample_data_dir,
        # Config default path may not exist; use the repo's shipped annotation TSV.
        annotation_path=(Path(PATHS.annotations_dir) / "rppa" / "TCGA_antibodies_descriptions.gencode.v36.tsv"),
        output_dir=out_dir,
        metadata_path=RPPAPathConfig.sample_metadata_tumor,
        rna_expr=None,
        save_outputs=True,
    )
    print(f"[OK] wrote RPPA outputs to: {out_dir}")


if __name__ == "__main__":
    main()

