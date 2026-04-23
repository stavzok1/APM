#!/usr/bin/env python3
"""
Run one representative existing input per oncology module (SNV, SV, CNV, Methylation, RPPA)
into an isolated tree under ``analysis/smoke_oncology_modules/run_<UTC>/`` to verify wiring.

Defaults favor a **fast** SV pass (VEP/motifs/ChIP skipped). Use ``--full-sv`` for the full SV stack.

Usage (repo root)::

    .venv/bin/python3 scripts/oncology/smoke_oncology_one_per_module.py
    .venv/bin/python3 scripts/oncology/smoke_oncology_one_per_module.py --full-sv
    .venv/bin/python3 scripts/oncology/smoke_oncology_one_per_module.py --skip-snv --skip-sv   # CNV + Methylation + RPPA only
    .venv/bin/python3 scripts/oncology/smoke_oncology_one_per_module.py --only-rppa            # RPPA only
    .venv/bin/python3 scripts/oncology/smoke_oncology_one_per_module.py --skip-snv --skip-sv --skip-cnv --skip-methylation
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _symlink(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if not src.exists():
        return
    if os.name == "posix":
        os.symlink(src.resolve(), dst, target_is_directory=False)
    else:
        shutil.copy2(src, dst)


def _first_glob(root: Path, patterns: Tuple[str, ...]) -> Optional[Path]:
    if not root.is_dir():
        return None
    for pat in patterns:
        found = sorted(root.glob(pat))
        if found:
            return found[0]
    return None


def _first_methylation_row(manifest: Path, beta_dir: Path) -> Optional[Tuple[str, Path]]:
    if not manifest.exists():
        return None
    import pandas as pd

    df = pd.read_csv(manifest, sep="\t", low_memory=False)
    if "Sample ID" not in df.columns or "File Name" not in df.columns:
        return None
    for _, row in df.iterrows():
        fn = str(row["File Name"]).strip()
        sid = str(row["Sample ID"]).strip()
        p = beta_dir / fn
        if fn and sid and p.is_file():
            return sid, p
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--full-sv",
        action="store_true",
        help="Run SV with motifs + ChIP (slow; needs tools and ChIP cache).",
    )
    ap.add_argument(
        "--rppa-n-samples",
        type=int,
        default=25,
        help="Symlink this many RPPA TSVs into the smoke dir for QC (default 25).",
    )
    ap.add_argument("--skip-snv", action="store_true", help="Do not run SNV smoke.")
    ap.add_argument("--skip-sv", action="store_true", help="Do not run SV smoke.")
    ap.add_argument("--skip-cnv", action="store_true", help="Do not run CNV smoke.")
    ap.add_argument(
        "--skip-methylation",
        action="store_true",
        help="Do not run methylation smoke.",
    )
    ap.add_argument("--skip-rppa", action="store_true", help="Do not run RPPA smoke.")
    ap.add_argument(
        "--only-rppa",
        action="store_true",
        help="Shorthand: skip SNV, SV, CNV, and Methylation (RPPA only).",
    )
    args = ap.parse_args()

    if args.only_rppa:
        args.skip_snv = True
        args.skip_sv = True
        args.skip_cnv = True
        args.skip_methylation = True

    from pipeline.config import PATHS, PIPELINE_GENE_PANEL

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out = ROOT / "analysis" / "smoke_oncology_modules" / f"run_{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    report: Dict[str, str] = {"output_root": str(out)}

    # --- SNV: first VCF under raw or vep dirs ---
    if args.skip_snv:
        print("[SNV] skip (--skip-snv)")
        report["snv"] = "skipped_cli"
    else:
        from pipeline.SNV.vcf_loader import load_mutect_snv_batch

        snv_src = _first_glob(PATHS.snv_raw_vcf_dir, ("*.vcf.gz", "*.vcf")) or _first_glob(
            PATHS.snv_vcf_dir, ("*.vcf.gz", "*.vcf")
        )
        snv_in = out / "snv_inputs"
        if snv_src:
            snv_in.mkdir(parents=True, exist_ok=True)
            _symlink(snv_src, snv_in / snv_src.name)
            print(f"[SNV] {snv_src.name}")
            load_mutect_snv_batch(
                snv_in,
                primary_genes=PIPELINE_GENE_PANEL,
                output_dir=out / "snv_output",
                samples_tsv=PATHS.annotations_dir / "SNV" / "samples.tsv"
                if (PATHS.annotations_dir / "SNV" / "samples.tsv").exists()
                else None,
                save_per_sample=True,
                save_combined=True,
                pattern="*.vcf*",
            )
            report["snv"] = "ok"
        else:
            print("[SNV] skip (no *.vcf / *.vcf.gz under PATHS.snv_raw_vcf_dir or snv_vcf_dir)")
            report["snv"] = "skip_no_input"

    # --- SV: first VCF under sv_vcf_dir ---
    if args.skip_sv:
        print("[SV] skip (--skip-sv)")
        report["sv"] = "skipped_cli"
    else:
        from pipeline.SV.pipeline import run_sv_pipeline

        sv_src = _first_glob(PATHS.sv_vcf_dir, ("*.vcf.gz", "*.vcf")) or _first_glob(
            PATHS.working_dir / "SV" / "raw_somatic_sv_test", ("*.vcf.gz", "*.vcf")
        )
        sv_in = out / "sv_inputs"
        if sv_src:
            sv_in.mkdir(parents=True, exist_ok=True)
            _symlink(sv_src, sv_in / sv_src.name)
            print(f"[SV] {sv_src.name} (full={args.full_sv})")
            run_sv_pipeline(
                vcf_dir=sv_in,
                output_root=out / "sv_pipeline_output",
                skip_vep=True,
                skip_motifs=not args.full_sv,
                skip_chip=not args.full_sv,
            )
            report["sv"] = "ok_full" if args.full_sv else "ok_fast"
        else:
            print("[SV] skip (no VCF under PATHS.sv_vcf_dir or SV/raw_somatic_sv_test)")
            report["sv"] = "skip_no_input"

    # --- CNV: first segment file ---
    if args.skip_cnv:
        print("[CNV] skip (--skip-cnv)")
        report["cnv"] = "skipped_cli"
    else:
        from pipeline.CNV.runner import process_cnv_directory

        cnv_src = _first_glob(PATHS.cnv_dir, ("*.seg.txt", "*.seg", "*.txt", "*.tsv"))
        cnv_in = out / "cnv_inputs"
        if cnv_src:
            cnv_in.mkdir(parents=True, exist_ok=True)
            _symlink(cnv_src, cnv_in / cnv_src.name)
            print(f"[CNV] {cnv_src.name}")
            genes_path = PATHS.cnv_genes if PATHS.cnv_genes.exists() else PATHS.genes_all_features
            process_cnv_directory(
                str(cnv_in),
                str(genes_path),
                str(PATHS.lncrnas_genes_centric),
                str(PATHS.regulatory_elements_table),
                str(PATHS.cnv_annotations_path),
                str(out / "cnv_annotated"),
                mirna_path=str(PATHS.mirna_mature_loci_csv)
                if PATHS.mirna_mature_loci_csv.exists()
                else str(PATHS.mirna_path),
                gene_tables_root=str(out / "cnv_gene_tables"),
            )
            report["cnv"] = "ok"
        else:
            print("[CNV] skip (no segment file under PATHS.cnv_dir)")
            report["cnv"] = "skip_no_input"

    # --- Methylation: first manifest row with existing beta file; isolated bundle ---
    if args.skip_methylation:
        print("[Methylation] skip (--skip-methylation)")
        report["methylation"] = "skipped_cli"
    else:
        from pipeline.Methylation.methylation_table import run_methylation_pipeline

        meth_bundle = out / "methylation_smoke"
        meth_row = _first_methylation_row(PATHS.methylation_sample_manifest, PATHS.methylation_samples_dir)
        if meth_row:
            sid, beta_p = meth_row
            import pandas as pd

            sub = pd.DataFrame([{"Sample ID": sid, "File Name": beta_p.name}])
            man = out / "methylation_one_sample.tsv"
            sub.to_csv(man, sep="\t", index=False)
            print(f"[Methylation] {sid} / {beta_p.name}")
            run_methylation_pipeline(
                probe_reference_path=PATHS.methylation_probe_reference,
                sample_manifest_path=man,
                sample_beta_dir=PATHS.methylation_samples_dir,
                working_dir=meth_bundle,
                build_reference=True,
                build_per_sample=True,
                build_cohort=False,
            )
            report["methylation"] = "ok"
        else:
            print("[Methylation] skip (no manifest row with existing beta under methylation_samples_dir)")
            report["methylation"] = "skip_no_input"

    # --- RPPA: symlink N tumor TSVs, relaxed QC ---
    if args.skip_rppa:
        print("[RPPA] skip (--skip-rppa)")
        report["rppa"] = "skipped_cli"
    else:
        from pipeline.rppa.rppa_main import run_rppa_pipeline

        rppa_dir = PATHS.rppa_samples_dir
        rppa_in = out / "rppa_inputs"
        if rppa_dir.is_dir():
            files = sorted(rppa_dir.glob("*.tsv"))[: max(1, args.rppa_n_samples)]
            if files:
                rppa_in.mkdir(parents=True, exist_ok=True)
                for p in files:
                    _symlink(p, rppa_in / p.name)
                print(f"[RPPA] {len(files)} samples → {rppa_in}")
                run_rppa_pipeline(
                    sample_dir=rppa_in,
                    annotation_path=PATHS.rppa_antibody_annotation_csv,
                    output_dir=out / "rppa_processed",
                    metadata_path=PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv"
                    if (PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv").exists()
                    else None,
                    min_samples_per_target=2,
                    min_targets_per_sample=5,
                    save_outputs=True,
                )
                report["rppa"] = "ok"
            else:
                print("[RPPA] skip (no *.tsv under rppa_samples_dir)")
                report["rppa"] = "skip_no_input"
        else:
            print("[RPPA] skip (rppa_samples_dir missing)")
            report["rppa"] = "skip_no_dir"

    summary_path = out / "smoke_report.json"
    summary_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nWrote {summary_path}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
