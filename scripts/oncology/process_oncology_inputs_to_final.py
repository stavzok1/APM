#!/usr/bin/env python3
"""
Drive each oncology submodule from its **input** dirs in ``pipeline.config.PATHS``
and write **final** tables to the canonical output locations documented in
``pipeline/config.py`` (and ``OUTPUT_SUBDIRS`` where applicable).

Inputs (defaults):
  SNV         → PATHS.snv_vcf_dir
  SV          → PATHS.sv_vcf_dir
  CNV         → PATHS.cnv_dir (segment files)
  Methylation → PATHS.methylation_samples_dir (+ manifest)
  RPPA        → PATHS.rppa_samples_dir

Outputs (defaults):
  SNV         → PATHS.snv_output_dir (``OUTPUT_SUBDIRS['snv']`` / per_sample/…)
  SV          → PATHS.sv_output_root (final ChIP tables in ``PATHS.sv_chip_enriched_dir`` when run with motifs+ChIP)
  CNV         → PATHS.cnv_output_dir + PATHS.cnv_gene_tables_dir
  Methylation → ``<working_dir>/Methylation/{reference,per_sample,cohort}`` (OUTPUT_SUBDIRS)
  RPPA        → PATHS.rppa_processed_dir

Examples::

    ./.venv/bin/python3 scripts/oncology/process_oncology_inputs_to_final.py --modules snv,cnv
    ./.venv/bin/python3 scripts/oncology/process_oncology_inputs_to_final.py --modules sv --sv-skip-vep --sv-skip-motifs --sv-skip-chip
    ./.venv/bin/python3 scripts/oncology/process_oncology_inputs_to_final.py --modules sv --sv-with-chip
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))


def _repo_root() -> Path:
    return _REPO


def cmd_snv(args: argparse.Namespace) -> None:
    from pipeline.config import PATHS, PIPELINE_GENE_PANEL
    from pipeline.SNV.vcf_loader import load_mutect_snv_batch

    vcf_dir = Path(args.snv_vcf_dir or PATHS.snv_vcf_dir)
    out_dir = Path(args.snv_out_dir or PATHS.snv_output_dir)
    samples_tsv = Path(args.snv_samples_tsv or (PATHS.annotations_dir / "SNV" / "samples.tsv"))
    print(f"[SNV] vcf_dir={vcf_dir}\n      out_dir={out_dir}")
    load_mutect_snv_batch(
        vcf_dir,
        primary_genes=PIPELINE_GENE_PANEL,
        output_dir=out_dir,
        samples_tsv=samples_tsv if samples_tsv.exists() else None,
        save_per_sample=True,
        save_combined=True,
        pattern="*.vcf*",
        run_fimo=args.snv_fimo,
    )


def cmd_sv(args: argparse.Namespace) -> None:
    from pipeline.SV.pipeline import run_sv_pipeline
    from pipeline.config import PATHS

    print(
        f"[SV] vcf_dir={PATHS.sv_vcf_dir}\n"
        f"     output_root={PATHS.sv_output_root}\n"
        f"     chip_enriched (when enabled)={PATHS.sv_chip_enriched_dir}"
    )
    run_sv_pipeline(
        vcf_dir=PATHS.sv_vcf_dir,
        output_root=PATHS.sv_output_root,
        skip_vep=args.sv_skip_vep,
        skip_motifs=args.sv_skip_motifs,
        skip_chip=args.sv_skip_chip,
    )


def cmd_cnv(args: argparse.Namespace) -> None:
    from pipeline.CNV.runner import process_cnv_directory
    from pipeline.config import PATHS

    genes_path = Path(args.cnv_genes_path or PATHS.cnv_genes)
    if not genes_path.exists():
        genes_path = PATHS.genes_all_features
    out_dir = Path(args.cnv_out_dir or PATHS.cnv_output_dir)
    gene_root = Path(args.cnv_gene_tables_dir or PATHS.cnv_gene_tables_dir)
    gene_root.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[CNV] cnv_dir={PATHS.cnv_dir}\n      annotated={out_dir}\n      gene_tables={gene_root}")
    process_cnv_directory(
        str(PATHS.cnv_dir),
        str(genes_path),
        str(PATHS.lncrnas_genes_centric),
        str(PATHS.regulatory_elements_table),
        str(PATHS.cnv_annotations_path),
        str(out_dir),
        mirna_path=str(PATHS.mirna_mature_loci_csv)
        if PATHS.mirna_mature_loci_csv.exists()
        else str(PATHS.mirna_path),
        gene_tables_root=str(gene_root),
    )


def cmd_methylation(args: argparse.Namespace) -> None:
    from pipeline.Methylation.methylation_table import run_methylation_pipeline
    from pipeline.config import PATHS

    print(
        f"[Methylation] betas_dir={PATHS.methylation_samples_dir}\n"
        f"              bundle_dir={PATHS.methylation_output_dir}"
    )
    run_methylation_pipeline(
        probe_reference_path=PATHS.methylation_probe_reference,
        sample_manifest_path=PATHS.methylation_sample_manifest,
        sample_beta_dir=PATHS.methylation_samples_dir,
        working_dir=PATHS.methylation_output_dir,
        build_reference=not args.methylation_skip_reference,
        build_per_sample=True,
        build_cohort=args.methylation_cohort,
    )


def cmd_rppa(args: argparse.Namespace) -> None:
    from pipeline.rppa.rppa_main import run_rppa_pipeline
    from pipeline.config import PATHS

    sample_dir = Path(args.rppa_samples_dir or PATHS.rppa_samples_dir)
    out_dir = Path(args.rppa_out_dir or PATHS.rppa_processed_dir)
    ann = Path(args.rppa_annotation or PATHS.rppa_antibody_annotation_csv)
    print(f"[RPPA] samples={sample_dir}\n       out={out_dir}")
    run_rppa_pipeline(
        sample_dir=sample_dir,
        annotation_path=ann,
        output_dir=out_dir,
        metadata_path=PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv"
        if (PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv").exists()
        else None,
        save_outputs=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--modules",
        type=str,
        default="snv,cnv,methylation,rppa",
        help="Comma-separated: snv, sv, cnv, methylation, rppa",
    )
    # SNV
    ap.add_argument("--snv-vcf-dir", type=Path, default=None)
    ap.add_argument("--snv-out-dir", type=Path, default=None)
    ap.add_argument("--snv-samples-tsv", type=Path, default=None)
    ap.add_argument(
        "--snv-fimo",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="SNV: run MEME FIMO on reference windows (needs bedtools + fimo on PATH).",
    )
    # SV
    ap.add_argument("--sv-skip-vep", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--sv-skip-motifs", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--sv-skip-chip", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument(
        "--sv-with-chip",
        action="store_true",
        help="Run motif + ChIP steps (needs bedtools, FIMO, unified ChIP); implies --no-sv-skip-motifs --no-sv-skip-chip",
    )
    # CNV
    ap.add_argument("--cnv-genes-path", type=Path, default=None)
    ap.add_argument("--cnv-out-dir", type=Path, default=None)
    ap.add_argument("--cnv-gene-tables-dir", type=Path, default=None)
    # Methylation
    ap.add_argument("--methylation-cohort", action="store_true", help="Also build cohort matrices (slow).")
    ap.add_argument(
        "--methylation-skip-reference",
        action="store_true",
        help="Reuse existing probe reference under Methylation/reference if present.",
    )
    # RPPA
    ap.add_argument("--rppa-samples-dir", type=Path, default=None)
    ap.add_argument("--rppa-out-dir", type=Path, default=None)
    ap.add_argument("--rppa-annotation", type=Path, default=None)

    args = ap.parse_args()

    mods = [m.strip().lower() for m in args.modules.split(",") if m.strip()]
    if args.sv_with_chip:
        args.sv_skip_motifs = False
        args.sv_skip_chip = False

    for m in mods:
        if m == "snv":
            cmd_snv(args)
        elif m == "sv":
            cmd_sv(args)
        elif m == "cnv":
            cmd_cnv(args)
        elif m in ("methylation", "meth"):
            cmd_methylation(args)
        elif m == "rppa":
            cmd_rppa(args)
        else:
            print(f"Unknown module: {m}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
