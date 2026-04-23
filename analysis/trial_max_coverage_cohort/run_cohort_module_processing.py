#!/usr/bin/env python3
"""
Run SNV, SV (motifs + ChIP → 09_chip_enriched), CNV, methylation, and RPPA pipelines
for the representative tumor vials listed in ``cohort_summary.tsv`` (from
``build_trial_cohort.py``).

Creates a **scratch** directory of symlinks so each module only sees the cohort inputs,
then runs the same processors as ``scripts/oncology/process_oncology_inputs_to_final.py``.

By default, **outputs** go to canonical ``PATHS`` locations (SNV/SV/RPPA/Methylation/CNV),
not under the scratch tree. Use ``--isolate-module-outputs`` to keep SNV/SV/RPPA outputs
isolated under scratch (legacy trial behavior).

**50-sample max-coverage, canonical outputs (bash; copy as a block)** — do not type
literal placeholders; use the directory ``build_trial_cohort`` prints (or capture it in ``$OUT``)::

    OUT=analysis/trial_max_coverage_cohort/run_$(date -u +%Y%m%d_%H%M%S)
    .venv/bin/python3 analysis/trial_max_coverage_cohort/build_trial_cohort.py --n 50 --out "$OUT"
    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary "$OUT/cohort_summary.tsv"

If you already built once and know the folder (example: ``run_20260418_074615``)::

    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary analysis/trial_max_coverage_cohort/run_20260418_074615/cohort_summary.tsv

Methylation only, reusing an existing annotated probe table at
``PATHS.methylation_output_dir/reference/probe_annotations.parquet``::

    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary path/to/cohort_summary.tsv \\
        --methylation-only --methylation-reuse-probe-reference

Run a single module (examples)::

    # CNV only
    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary path/to/cohort_summary.tsv \\
        --cnv-only

    # SV only
    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary path/to/cohort_summary.tsv \\
        --sv-only

    # RPPA only
    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary path/to/cohort_summary.tsv \\
        --rppa-only

Optional: ``--max-samples 50`` to cap rows if ``cohort_summary.tsv`` lists more specimens
than you want to process in one pass.

Usage (from repo root)::

    .venv/bin/python3 analysis/trial_max_coverage_cohort/run_cohort_module_processing.py \\
        --cohort-summary analysis/trial_max_coverage_cohort/run_20260417_170605/cohort_summary.tsv
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
import json
from typing import Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pandas as pd

from pipeline.sample_ids import normalize_tcga_id


def _load_build_trial_cohort():
    p = ROOT / "analysis" / "trial_max_coverage_cohort" / "build_trial_cohort.py"
    spec = importlib.util.spec_from_file_location("build_trial_cohort", p)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _fid_norm(fid: str) -> str:
    return str(fid).replace("-", "").lower()


def _resolve_manifest_file(fn: str, file_id: str, search_dirs: List[Path]) -> Optional[Path]:
    """
    Locate a GDC manifest ``File Name`` on disk. Tries exact name under each search dir,
    common suffix variants, then any file whose name contains the normalized File ID.
    """
    fn = str(fn).strip()
    fid = str(file_id).strip()
    seen: Set[str] = set()
    dirs: List[Path] = []
    for d in search_dirs:
        if d is None:
            continue
        d = Path(d)
        if not d.is_dir():
            continue
        key = str(d.resolve())
        if key in seen:
            continue
        seen.add(key)
        dirs.append(d)

    if not fn:
        return None

    name = Path(fn).name
    alts = {fn, name}
    if fn.endswith(".vcf.gz"):
        alts.add(fn[:-3])  # .vcf
    elif fn.endswith(".vcf"):
        alts.add(fn + ".gz")
    if fn.endswith(".txt"):
        alts.add(fn[:-4] + ".seg")
    elif fn.endswith(".seg"):
        alts.add(str(Path(fn).with_suffix(".txt")))

    for d in dirs:
        for cand in alts:
            p = d / cand
            if p.is_file():
                return p

    if fid and fid.lower() != "nan":
        fnorm = _fid_norm(fid)
        for d in dirs:
            try:
                for p in d.iterdir():
                    if not p.is_file():
                        continue
                    pn = p.name.replace("-", "").lower()
                    if fnorm in pn or name == p.name:
                        return p
            except OSError:
                continue

    return None


def _lookup_gdc_manifest_row(
    tumor_map: Dict[str, Tuple[str, str]],
    vial: str,
) -> Optional[Tuple[str, str]]:
    """Match cohort ``vial`` to a manifest row when dict keys differ slightly (sample vs vial, aliquot)."""
    if vial in tumor_map:
        return tumor_map[vial]
    tv = normalize_tcga_id(vial)
    for mk_raw, row in tumor_map.items():
        mk = str(mk_raw).strip()
        mv = normalize_tcga_id(mk)
        if tv.aliquot and mv.aliquot and tv.aliquot == mv.aliquot:
            return row
        if tv.sample_vial and mv.sample_vial and tv.sample_vial == mv.sample_vial:
            return row
        if tv.sample and mv.sample and tv.sample == mv.sample:
            return row
    return None


def _manifest_tumor_map_multi(tsv: Path) -> Dict[str, List[Tuple[str, str]]]:
    """
    tumor_sample_vial -> list[(file_name, file_id)]

    CNV (and some other modalities) can have multiple rows per tumor sample (e.g.
    ASCAT3 segment-level + ASCAT3 gene-level). The original single-row mapping
    drops one of them; use this to preserve all rows.
    """
    btc = _load_build_trial_cohort()
    if not tsv.exists():
        return {}
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    tumor_row = btc._load_sample_module_coverage().tumor_sample_vial_from_gdc_manifest_row
    out: Dict[str, List[Tuple[str, str]]] = {}
    for _, row in df.iterrows():
        tv = tumor_row(row)
        if not tv:
            continue
        fn = str(row.get("File Name", "")).strip()
        fid = str(row.get("File ID", "")).strip()
        if not fn:
            continue
        out.setdefault(str(tv), []).append((fn, fid))
    # de-dup while keeping order
    for k, xs in list(out.items()):
        seen = set()
        uniq: List[Tuple[str, str]] = []
        for fn, fid in xs:
            key = (fn, fid)
            if key in seen:
                continue
            seen.add(key)
            uniq.append((fn, fid))
        out[k] = uniq
    return out


def _lookup_gdc_manifest_rows(
    tumor_map_multi: Dict[str, List[Tuple[str, str]]],
    vial: str,
) -> List[Tuple[str, str]]:
    """Multi-row version of `_lookup_gdc_manifest_row` (returns all matches, possibly empty)."""
    if vial in tumor_map_multi:
        return tumor_map_multi[vial]
    tv = normalize_tcga_id(vial)
    out: List[Tuple[str, str]] = []
    for mk_raw, rows in tumor_map_multi.items():
        mk = str(mk_raw).strip()
        mv = normalize_tcga_id(mk)
        if tv.aliquot and mv.aliquot and tv.aliquot == mv.aliquot:
            out.extend(rows)
            continue
        if tv.sample_vial and mv.sample_vial and tv.sample_vial == mv.sample_vial:
            out.extend(rows)
            continue
        if tv.sample and mv.sample and tv.sample == mv.sample:
            out.extend(rows)
            continue
    # de-dup while keeping order
    seen = set()
    uniq: List[Tuple[str, str]] = []
    for fn, fid in out:
        key = (fn, fid)
        if key in seen:
            continue
        seen.add(key)
        uniq.append((fn, fid))
    return uniq


def _resolve_snv_vcf_deep(
    fn: str,
    file_id: str,
    flat_dirs: List[Path],
    snv_tree_root: Path,
    *,
    max_rglob_files: int = 80_000,
) -> Optional[Path]:
    """Flat dirs first, then recursive ``*.vcf*`` under ``snv_tree_root`` (e.g. ``data/SNV``)."""
    hit = _resolve_manifest_file(fn, file_id, flat_dirs)
    if hit:
        return hit
    if not snv_tree_root.is_dir():
        return None
    name = Path(fn).name
    fid = str(file_id).strip()
    fnorm = _fid_norm(fid) if fid and fid.lower() != "nan" else ""
    n = 0
    for p in snv_tree_root.rglob("*.vcf*"):
        n += 1
        if n > max_rglob_files:
            break
        if not p.is_file():
            continue
        if p.name == name:
            return p
        if fnorm and fnorm in p.name.replace("-", "").lower():
            return p
    return None


def _symlink(src: Path, dst: Path) -> bool:
    if not src or not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        if os.name == "posix":
            os.symlink(src.resolve(), dst, target_is_directory=False)
        else:
            import shutil

            shutil.copy2(src, dst)
        return True
    except OSError:
        import shutil

        shutil.copy2(src, dst)
        return True


def _subset_methylation_manifest(
    full_manifest: Path,
    vials: List[str],
    out_path: Path,
) -> int:
    df = pd.read_csv(full_manifest, sep="\t", low_memory=False)
    if "Sample ID" not in df.columns:
        raise ValueError("Methylation manifest missing Sample ID")
    vset = set(vials)
    sub = df[df["Sample ID"].astype(str).str.strip().isin(vset)]
    if sub.empty:
        # try short sample ids (…-01A) vs full vials
        short = {v[:12] if len(v) > 12 else v for v in vials}
        sub = df[df["Sample ID"].astype(str).str.strip().isin(short | vset)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(out_path, sep="\t", index=False)
    return len(sub)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--cohort-summary",
        type=Path,
        default=ROOT / "analysis/trial_max_coverage_cohort/run_20260417_170605/cohort_summary.tsv",
    )
    ap.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Process at most this many vials from cohort_summary (head order; default: all rows).",
    )
    ap.add_argument(
        "--isolate-module-outputs",
        action="store_true",
        help="Write SNV/SV outputs under scratch and RPPA under PATHS.rppa_processed_dir/trial_cohort_<UTC> (legacy). "
        "Default: canonical PATHS outputs for SNV, SV, RPPA (same as process_oncology_inputs_to_final.py).",
    )
    ap.add_argument(
        "--scratch",
        type=Path,
        default=None,
        help="Working dir for input symlinks (default: trial_max_coverage_cohort/processing_scratch_<UTC>)",
    )
    ap.add_argument(
        "--sv-skip-vep",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="SV: skip VEP by default (faster). Pass --no-sv-skip-vep to run VEP on the cohort VCFs.",
    )
    ap.add_argument(
        "--sv-resume-motifs-from-fimo",
        action="store_true",
        help="SV: skip VCF→CSV; run motif steps 4–7 only (reuse 02_processed_sv_csv + 03_sv_bed + "
        "05_fimo_tsv; refresh 06/07/08/09). Does not re-run VEP, FIMO, or SV VCF parsing.",
    )
    ap.add_argument(
        "--snv-fimo",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="SNV: run MEME FIMO on a reference window per variant (needs bedtools + fimo on PATH). "
        "Use --no-snv-fimo to skip.",
    )
    ap.add_argument(
        "--snv-chip",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="SNV: overlap variants with unified ChIP peaks (strict POS overlap; ``PATHS.chip_unified``). "
        "Use --no-snv-chip to skip.",
    )
    ap.add_argument(
        "--extra-snv-dir",
        action="append",
        default=[],
        type=Path,
        help="Extra directory to search for VEP Mutect2 outputs (…APM_1Mb…; repeatable).",
    )
    ap.add_argument(
        "--extra-sv-dir",
        action="append",
        default=[],
        type=Path,
        help="Extra directory to search for SV VCFs (repeatable).",
    )
    ap.add_argument(
        "--extra-cnv-dir",
        action="append",
        default=[],
        type=Path,
        help="Extra directory to search for CNV segment files (repeatable).",
    )
    ap.add_argument(
        "--methylation-only",
        action="store_true",
        help="Only subset the methylation manifest and run the methylation pipeline (skip SNV, SV, CNV, RPPA).",
    )
    ap.add_argument(
        "--snv-only",
        action="store_true",
        help="Only link inputs and run SNV (skip SV, CNV, Methylation, RPPA).",
    )
    ap.add_argument(
        "--sv-only",
        action="store_true",
        help="Only link inputs and run SV (skip SNV, CNV, Methylation, RPPA).",
    )
    ap.add_argument(
        "--cnv-only",
        action="store_true",
        help="Only link inputs and run CNV (skip SNV, SV, Methylation, RPPA).",
    )
    ap.add_argument(
        "--rppa-only",
        action="store_true",
        help="Only link inputs and run RPPA (skip SNV, SV, CNV, Methylation).",
    )
    ap.add_argument(
        "--methylation-reuse-probe-reference",
        action="store_true",
        help="Do not rebuild the annotated probe reference: load existing "
        "``PATHS.methylation_output_dir/reference/probe_annotations.parquet`` if present "
        "(same as ``run_methylation_pipeline(build_reference=False)``).",
    )
    args = ap.parse_args()

    only_flags = [args.methylation_only, args.snv_only, args.sv_only, args.cnv_only, args.rppa_only]
    if sum(bool(x) for x in only_flags) > 1:
        raise SystemExit(
            "Only one of --methylation-only/--snv-only/--sv-only/--cnv-only/--rppa-only can be set."
        )

    run_snv = not any(only_flags) or args.snv_only
    run_sv = not any(only_flags) or args.sv_only
    run_cnv = not any(only_flags) or args.cnv_only
    run_meth = not any(only_flags) or args.methylation_only
    run_rppa = not any(only_flags) or args.rppa_only

    btc = _load_build_trial_cohort()
    from pipeline.config import PATHS, PIPELINE_GENE_PANEL
    from pipeline.CNV.runner import process_cnv_directory
    from pipeline.Methylation.methylation_table import run_methylation_pipeline
    from pipeline.SNV.mutect_manifest_paths import (
        resolve_vep_mutect_vcf,
        vep_apm1mb_filename_candidates,
    )
    from pipeline.SNV.vcf_loader import load_mutect_snv_batch
    from pipeline.SV.pipeline import run_sv_pipeline
    from pipeline.rppa.rppa_main import run_rppa_pipeline

    cohort = pd.read_csv(args.cohort_summary, sep="\t")
    if "representative_tumor_sample_vial" not in cohort.columns:
        raise SystemExit("cohort_summary.tsv must contain representative_tumor_sample_vial")
    rows = cohort.dropna(subset=["representative_tumor_sample_vial"])
    vials = [str(v).strip() for v in rows["representative_tumor_sample_vial"].tolist() if str(v).strip()]
    if args.max_samples is not None:
        cap = max(0, int(args.max_samples))
        vials = vials[:cap]
        print(f"Capped cohort to --max-samples={cap} → {len(vials)} vials")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    scratch = args.scratch or (ROOT / "analysis" / "trial_max_coverage_cohort" / f"processing_scratch_{stamp}")
    scratch.mkdir(parents=True, exist_ok=True)
    snv_in = scratch / "snv_inputs"
    sv_in = scratch / "sv_inputs"
    cnv_in = scratch / "cnv_inputs"
    rppa_in = scratch / "rppa_inputs"
    for d in (snv_in, sv_in, cnv_in, rppa_in):
        d.mkdir(parents=True, exist_ok=True)

    snv_tree = PATHS.working_dir / "SNV"
    # Manifest names match raw Mutect2 under vcfs_snv; ``load_mutect_snv_batch`` expects VEP ``…APM_1Mb…`` under vep_vcfs.
    snv_vep_dirs = [
        PATHS.snv_vcf_dir,
        snv_tree / "vep_vcfs",
        *args.extra_snv_dir,
    ]
    snv_search_dirs = [
        *snv_vep_dirs,
        PATHS.snv_raw_vcf_dir,
        snv_tree / "vcfs_snv",
        snv_tree / "raw_vcfs",
        snv_tree,
        *args.extra_snv_dir,
    ]
    sv_search_dirs = [
        PATHS.sv_vcf_dir,
        PATHS.working_dir / "SV" / "raw_somatic_sv",
        PATHS.working_dir / "SV" / "raw_somatic_sv_test",
        *args.extra_sv_dir,
    ]
    cnv_search_dirs = [
        PATHS.cnv_dir,
        PATHS.working_dir / "CNV_TCGA" / "CNV_extracted",
        *args.extra_cnv_dir,
    ]

    snv_m = btc._manifest_tumor_map(PATHS.annotations_dir / "SNV" / "samples.tsv")
    sv_m = btc._manifest_tumor_map(PATHS.annotations_dir / "SV" / "samples.tsv")
    cnv_m_multi = _manifest_tumor_map_multi(PATHS.annotations_dir / "CNV" / "samples.tsv")
    meth_df = (
        pd.read_csv(PATHS.methylation_sample_manifest, sep="\t", low_memory=False)
        if PATHS.methylation_sample_manifest.exists()
        else pd.DataFrame()
    )
    rppa_df = (
        pd.read_csv(PATHS.annotations_dir / "rppa" / "samples.tsv", sep="\t", low_memory=False)
        if (PATHS.annotations_dir / "rppa" / "samples.tsv").exists()
        else pd.DataFrame()
    )

    def _sample_manifest_lookup(df: pd.DataFrame) -> Dict[str, pd.Series]:
        if df.empty or "Sample ID" not in df.columns:
            return {}
        out: Dict[str, pd.Series] = {}
        for _, row in df.iterrows():
            sid = str(row["Sample ID"]).strip()
            if sid not in out:
                out[sid] = row
        return out

    meth_by_sid = _sample_manifest_lookup(meth_df)
    rppa_by_sid = _sample_manifest_lookup(rppa_df)

    log_lines: List[str] = []
    linked: Dict[str, Dict[str, bool]] = {v: {} for v in vials}
    # SNV: distinguish “row in GDC manifest” (used for max-coverage cohort rank) vs VCF on disk
    snv_prep = {"manifest_rows": 0, "vcf_linked": 0}

    if not (run_snv or run_sv or run_cnv or run_rppa):
        log_lines.append("no SNV/SV/CNV/RPPA enabled: skipped input linking")

    def meth_row(vial: str) -> Optional[pd.Series]:
        if not meth_by_sid:
            return None
        v = str(vial).strip()
        r = meth_by_sid.get(v)
        if r is not None:
            return r
        tid = normalize_tcga_id(vial)
        samp = (tid.sample or vial).strip()
        return meth_by_sid.get(samp)

    def rppa_row(vial: str) -> Optional[pd.Series]:
        if not rppa_by_sid:
            return None
        v = str(vial).strip()
        r = rppa_by_sid.get(v)
        if r is not None:
            return r
        tid = normalize_tcga_id(vial)
        samp = (tid.sample or vial).strip()
        return rppa_by_sid.get(samp)

    if run_snv or run_sv or run_cnv or run_rppa:
        for vial in vials:
            if run_snv:
                snv_t = _lookup_gdc_manifest_row(snv_m, vial)
                if snv_t:
                    snv_prep["manifest_rows"] += 1
                    fn, fid = snv_t[0], snv_t[1]
                    vep_src = resolve_vep_mutect_vcf(fn, vep_dirs=snv_vep_dirs)
                    if vep_src is None:
                        deep_dirs = list(snv_vep_dirs) + [
                            PATHS.snv_raw_vcf_dir,
                            snv_tree / "vcfs_snv",
                            snv_tree / "raw_vcfs",
                        ]
                        for cand in vep_apm1mb_filename_candidates(fn):
                            hit = _resolve_snv_vcf_deep(cand, fid, deep_dirs, snv_tree)
                            if hit is not None and Path(hit).is_file():
                                vep_src = Path(hit)
                                break
                    if vep_src is None or not vep_src.is_file():
                        linked[vial]["snv"] = False
                        log_lines.append(
                            f"MISS snv {vial} manifest={fn!r} file_id={fid!r}: "
                            f"no VEP APM_1Mb VCF (expected names like {vep_apm1mb_filename_candidates(fn)[:2]!r} under vep dirs)"
                        )
                    else:
                        ok = _symlink(vep_src, snv_in / vep_src.name)
                        linked[vial]["snv"] = ok
                        if ok:
                            snv_prep["vcf_linked"] += 1
                        if not ok:
                            log_lines.append(
                                f"MISS snv {vial} manifest={fn!r} found={vep_src!r} symlink_failed"
                            )
                else:
                    linked[vial]["snv"] = False
                    log_lines.append(f"MISS snv manifest key for vial={vial!r} (check annotations/SNV/samples.tsv)")

            if run_sv:
                sv_t = _lookup_gdc_manifest_row(sv_m, vial)
                if sv_t:
                    fn, fid = sv_t[0], sv_t[1]
                    src = _resolve_manifest_file(fn, fid, sv_search_dirs) or btc._sv_path(fn)
                    ok = _symlink(src, sv_in / Path(fn).name)
                    linked[vial]["sv"] = ok
                    if not ok:
                        log_lines.append(f"MISS sv {vial} manifest={fn!r} tried={sv_search_dirs!r} last={src}")
                else:
                    linked[vial]["sv"] = False
                    log_lines.append(f"MISS sv manifest key for vial={vial!r}")

            if run_cnv:
                cnv_rows = _lookup_gdc_manifest_rows(cnv_m_multi, vial)
                if cnv_rows:
                    ok_any = False
                    ok_n = 0
                    for fn, fid in cnv_rows:
                        src = _resolve_manifest_file(fn, fid, cnv_search_dirs) or btc._cnv_path(fn)
                        ok = _symlink(src, cnv_in / Path(fn).name)
                        ok_any = ok_any or ok
                        ok_n += int(bool(ok))
                        if not ok:
                            log_lines.append(
                                f"MISS cnv {vial} manifest={fn!r} tried={cnv_search_dirs!r} last={src}"
                            )
                    linked[vial]["cnv"] = bool(ok_any)
                    linked[vial]["cnv_linked_files"] = int(ok_n)
                    linked[vial]["cnv_manifest_files"] = int(len(cnv_rows))
                else:
                    linked[vial]["cnv"] = False
                    linked[vial]["cnv_linked_files"] = 0
                    linked[vial]["cnv_manifest_files"] = 0
                    log_lines.append(f"MISS cnv manifest {vial}")

            mr = meth_row(vial)
            if mr is not None and "File Name" in mr.index:
                src = btc._meth_path(str(mr["File Name"]))
                linked[vial]["methylation_manifest"] = True
            else:
                linked[vial]["methylation_manifest"] = False
                log_lines.append(f"MISS meth manifest row {vial}")

            if run_rppa:
                rr = rppa_row(vial)
                if rr is not None and "File Name" in rr.index:
                    src = btc._rppa_path(str(rr["File Name"]))
                    ok = _symlink(src, rppa_in / Path(str(rr["File Name"])).name)
                    linked[vial]["rppa"] = ok
                    if not ok:
                        log_lines.append(f"MISS rppa {vial} -> {src}")
                else:
                    linked[vial]["rppa"] = False
                    log_lines.append(f"MISS rppa manifest {vial}")

    meth_manifest = scratch / "methylation_cohort_manifest.tsv"
    if run_meth:
        n_meth = _subset_methylation_manifest(PATHS.methylation_sample_manifest, vials, meth_manifest)
        log_lines.append(f"Methylation manifest subset rows: {n_meth}")
    else:
        n_meth = 0
        log_lines.append("Methylation manifest subset rows: 0 (methylation disabled)")

    log_path = scratch / "cohort_processing_prepare.log"
    log_lines.append(
        f"SNV prep summary: manifest_hits={snv_prep['manifest_rows']} among {len(vials)} vials, "
        f"vcfs_linked_to_scratch={snv_prep['vcf_linked']}"
    )
    log_path.write_text("\n".join(log_lines) + "\n" + json.dumps(linked, indent=2) + "\n", encoding="utf-8")
    print(f"Scratch: {scratch}")
    print(f"Prepare log: {log_path}")
    if snv_prep["manifest_rows"] and snv_prep["vcf_linked"] == 0:
        print(
            "[prep] SNV: manifest matched cohort vials, but **no** VEP ``…APM_1Mb…`` Mutect2 VCF was resolved "
            "(see MISS snv lines in the prepare log). Raw names live under PATHS.snv_raw_vcf_dir (vcfs_snv); "
            "this step symlinks the VEP subset from PATHS.snv_vcf_dir (vep_vcfs). Use ``--extra-snv-dir`` for "
            "extra VEP trees."
        )
    elif snv_prep["manifest_rows"] == 0:
        print(
            "[prep] SNV: no cohort vial matched a tumor row in annotations/SNV/samples.tsv "
            "(TCGA id mismatch vs manifest keys)."
        )

    if args.isolate_module_outputs:
        snv_out = scratch / "snv_output"
        sv_out_root = scratch / "sv_pipeline_output"
        rppa_out = PATHS.rppa_processed_dir / f"trial_cohort_{stamp}"
        print("Module outputs: ISOLATED under scratch (SNV/SV) + trial RPPA subdir (--isolate-module-outputs).")
    else:
        snv_out = PATHS.snv_output_dir
        sv_out_root = PATHS.sv_output_root
        rppa_out = PATHS.rppa_processed_dir
        print("Module outputs: CANONICAL PATHS (PATHS.snv_output_dir, sv_output_root, rppa_processed_dir, …).")
    snv_out.mkdir(parents=True, exist_ok=True)
    sv_out_root.mkdir(parents=True, exist_ok=True)
    rppa_out.mkdir(parents=True, exist_ok=True)

    # --- SNV ---
    if run_snv:
        if any(snv_in.iterdir()):
            print("\n[SNV] batch …")
            load_mutect_snv_batch(
                snv_in,
                primary_genes=PIPELINE_GENE_PANEL,
                output_dir=snv_out,
                samples_tsv=PATHS.annotations_dir / "SNV" / "samples.tsv",
                save_per_sample=True,
                save_combined=True,
                # Require VEP-annotated Mutect2 VCFs (filename contains ``.vep.vcf``)
                pattern="*.vep.vcf*",
                run_fimo=args.snv_fimo,
                run_chip=args.snv_chip,
            )
        else:
            print("[SNV] skip (no inputs linked under snv_inputs/)")
            print(
                f"      Details: open {log_path} — look for MISS snv lines "
                "(vials must match annotations/SNV/samples.tsv; each row needs a VEP "
                "``…APM_1Mb…`` Mutect2 file under PATHS.snv_vcf_dir / data/SNV/vep_vcfs or ``--extra-snv-dir``)."
            )

    # --- SV ---
    if run_sv:
        if any(sv_in.iterdir()):
            print(f"\n[SV] pipeline (ChIP on) → {sv_out_root} …")
            run_sv_pipeline(
                vcf_dir=sv_in,
                output_root=sv_out_root,
                skip_vep=args.sv_skip_vep,
                skip_motifs=False,
                skip_chip=False,
                resume_motifs_from_fimo=args.sv_resume_motifs_from_fimo,
            )
        else:
            print("[SV] skip (no inputs linked)")

    # --- CNV ---
    if run_cnv:
        if any(cnv_in.iterdir()):
            print("\n[CNV] annotate …")
            genes_path = PATHS.cnv_genes if PATHS.cnv_genes.exists() else PATHS.genes_all_features
            process_cnv_directory(
                str(cnv_in),
                str(genes_path),
                str(PATHS.lncrnas_genes_centric),
                str(PATHS.regulatory_elements_table),
                str(PATHS.cnv_annotations_path),
                str(PATHS.cnv_output_dir),
                mirna_path=str(PATHS.mirna_mature_loci_csv)
                if PATHS.mirna_mature_loci_csv.exists()
                else str(PATHS.mirna_path),
                gene_tables_root=str(PATHS.cnv_gene_tables_dir),
            )
        else:
            print("[CNV] skip (no inputs linked)")

    # --- Methylation ---
    if run_meth and n_meth > 0 and meth_manifest.exists():
        ref_dir = PATHS.methylation_output_dir / "reference"
        ref_parquet = ref_dir / "probe_annotations.parquet"
        ref_csv = ref_dir / "probe_annotations.csv"
        if args.methylation_reuse_probe_reference:
            if ref_parquet.is_file():
                print(
                    f"\n[Methylation] reusing annotated probe reference ({ref_parquet}) "
                    "(--methylation-reuse-probe-reference)"
                )
            elif ref_csv.is_file():
                print(
                    f"\n[Methylation] reusing annotated probe reference ({ref_csv}) "
                    "(legacy CSV; prefer regenerating once to write probe_annotations.parquet)"
                )
            else:
                print(
                    f"\n[Methylation] --methylation-reuse-probe-reference set but neither "
                    f"{ref_parquet} nor {ref_csv} exists; will build reference from GDC probe table."
                )
        print("\n[Methylation] per-sample + cohort (subset) …")
        run_methylation_pipeline(
            probe_reference_path=PATHS.methylation_probe_reference,
            sample_manifest_path=meth_manifest,
            sample_beta_dir=PATHS.methylation_samples_dir,
            working_dir=PATHS.methylation_output_dir,
            build_reference=not args.methylation_reuse_probe_reference,
            build_per_sample=True,
            build_cohort=True,
        )
    else:
        if run_meth:
            print("[Methylation] skip (no manifest rows)")

    # --- RPPA ---
    if run_rppa:
        if any(rppa_in.iterdir()):
            print("\n[RPPA] cohort subset …")
            run_rppa_pipeline(
                sample_dir=rppa_in,
                annotation_path=PATHS.rppa_antibody_annotation_csv,
                output_dir=rppa_out,
                metadata_path=PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv"
                if (PATHS.annotations_dir / "rppa" / "metadata" / "TUMOR" / "samples.tsv").exists()
                else None,
                min_samples_per_target=3,
                min_targets_per_sample=15,
                save_outputs=True,
            )
        else:
            print("[RPPA] skip (no inputs linked)")

    summary = {
        "scratch": str(scratch),
        "methylation_only": bool(args.methylation_only),
        "methylation_reuse_probe_reference": bool(args.methylation_reuse_probe_reference),
        "isolate_module_outputs": bool(args.isolate_module_outputs),
        "snv_output": str(snv_out),
        "sv_output_root": str(sv_out_root),
        "sv_chip_enriched": str(sv_out_root / "09_chip_enriched"),
        "rppa_output_dir": str(rppa_out),
        "cnv_annotated_dir": str(PATHS.cnv_output_dir),
        "cnv_gene_tables_dir": str(PATHS.cnv_gene_tables_dir),
        "methylation_working_dir": str(PATHS.methylation_output_dir),
        "methylation_subset_manifest": str(meth_manifest),
    }
    (scratch / "cohort_processing_outputs.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\nDone. Output index:", scratch / "cohort_processing_outputs.json")


if __name__ == "__main__":
    main()
