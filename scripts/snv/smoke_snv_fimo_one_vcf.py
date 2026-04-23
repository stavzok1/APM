#!/usr/bin/env python3
"""
Smoke test: load one small VEP Mutect2 slice with ``run_fimo=True``.

``fimo`` is resolved like the pipeline (``APM_FIMO_BIN``, PATH, or ``VEP_ENV``/``vep_env``
conda roots). ``bedtools`` is required too.
Run from repo root: ``.venv/bin/python3 scripts/snv/smoke_snv_fimo_one_vcf.py``
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import pandas as pd

from pipeline.config import PATHS, PIPELINE_GENE_PANEL  # noqa: E402
from pipeline.SNV.snv_fimo import resolve_fimo_argv  # noqa: E402
from pipeline.SNV.vcf_loader import load_mutect_snv_vcf  # noqa: E402


def main() -> None:
    vcf = ROOT / (
        "data/SNV/vep_vcfs/TCGA-BRCA.001ccb2b-4d30-4d95-aa0a-20e999095361."
        "wgs.tumor_normal.gatk4_mutect2.raw_somatic_mutation.APM_1Mb.vep.vcf"
    )
    if not vcf.is_file():
        raise SystemExit(f"Missing VCF: {vcf}")
    fimo0 = resolve_fimo_argv()[0]
    fimo_ok = Path(fimo0).is_file() or shutil.which(fimo0) is not None
    run_fimo = True
    if not fimo_ok:
        print("[SKIP] fimo not found (APM_FIMO_BIN, PATH, VEP_ENV/vep_env). Install MEME to exercise motif scanning.")
        print("       Pipeline still loads VCF; run with fimo available to validate hits.")
        run_fimo = False
    out = Path("/tmp/apm_snv_fimo_smoke")
    out.mkdir(parents=True, exist_ok=True)
    mirna_df = None
    if PATHS.mirna_mature_loci_csv.is_file():
        mirna_df = pd.read_csv(PATHS.mirna_mature_loci_csv)
    df, _n, _t = load_mutect_snv_vcf(
        vcf,
        primary_genes=PIPELINE_GENE_PANEL,
        elements_path=PATHS.regulatory_elements_table,
        mirna_df=mirna_df,
        apply_filter=True,
        save_outputs=False,
        run_fimo=run_fimo,
        fimo_work_dir=out / "fimo_work",
        ref_fasta=PATHS.sv_reference_fasta,
        meme_file=PATHS.sv_meme_file,
    )
    print("rows", len(df))
    if "fimo_hits" not in df.columns:
        raise SystemExit("missing fimo_hits column")
    nh = int((df["fimo_hits"].apply(len) > 0).sum())
    print("with_fimo_hits", nh)
    if nh:
        row = df[df["fimo_hits"].apply(len) > 0].iloc[0]
        print("example_motif", row["fimo_hits"][0].get("motif_id"))


if __name__ == "__main__":
    main()
