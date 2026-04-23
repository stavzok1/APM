from __future__ import annotations

import argparse
from pathlib import Path

from pipeline.config import OUTPUT_SUBDIRS, PATHS
from pipeline.SNV.vcf_loader import load_mutect_snv_batch


def main() -> None:
    parser = argparse.ArgumentParser(description="Create per-sample SNV outputs from VEP VCFs")
    parser.add_argument(
        "--pattern",
        default="*.vep.vcf",
        help="Glob pattern within PATHS.snv_vcf_dir (default: *.vep.vcf)",
    )
    args = parser.parse_args()

    out_dir = Path(PATHS.working_dir) / OUTPUT_SUBDIRS["snv"]
    out_dir.mkdir(parents=True, exist_ok=True)

    load_mutect_snv_batch(
        vcf_dir=PATHS.snv_vcf_dir,
        pattern=args.pattern,
        output_dir=out_dir,
        samples_tsv=PATHS.snv_samples_tsv,
        save_per_sample=True,
        save_combined=False,
        save_outputs=False,
    )


if __name__ == "__main__":
    main()

