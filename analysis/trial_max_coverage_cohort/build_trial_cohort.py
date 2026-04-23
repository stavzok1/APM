#!/usr/bin/env python3
"""
Select ~N primary-tumor (01) specimens with maximal omics module coverage, then
materialize a small trial folder: cohort tables, wide-matrix subsets, and
per-sample path manifests (symlinks for small files when possible).

Core modules (must maximize first): SNV, SV, CNV, Methylation, RPPA
Secondary (wide tables / sparse assays): RNA, miRNA, ATAC, HLA, HiCHIP

Uses the same presence logic as ``analysis/sample_module_coverage.py``.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import normalize_tcga_id, tcga_sample_type_two_digit  # noqa: E402

CORE_MODULES = ["SNV", "SV", "CNV", "Methylation", "RPPA"]
WIDE_MODULES = ["RNA", "miRNA", "ATAC", "HLA"]


def _load_sample_module_coverage():
    p = ROOT / "analysis" / "sample_module_coverage.py"
    spec = importlib.util.spec_from_file_location("sample_module_coverage", p)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _header_tcga_columns(path: Path) -> List[str]:
    if not path.exists():
        return []
    best: List[str] = []
    best_n = -1
    for sep in ("\t", ","):
        for kwargs in ({}, {"index_col": 0}):
            try:
                hdr = pd.read_csv(path, sep=sep, nrows=0, low_memory=False, **kwargs)
            except Exception:
                continue
            cols = [str(c).strip() for c in hdr.columns]
            n = sum(1 for c in cols if c.startswith("TCGA"))
            if n > best_n:
                best_n = n
                best = cols
    return [c for c in best if c.startswith("TCGA")]


def _participant_from_tcga_col(col: str) -> Optional[str]:
    tid = normalize_tcga_id(col)
    return tid.participant


def _is_primary_tumor_barcode(col: str) -> bool:
    return tcga_sample_type_two_digit(col) == "01"


def _group_tumor_presence_by_participant(
    vial_sets: Dict[str, Set[str]],
) -> Dict[str, Set[str]]:
    """participant -> set of module names with at least one primary-tumor (01) hit."""
    out: Dict[str, Set[str]] = {}
    for mod, vs in vial_sets.items():
        for raw in vs:
            if tcga_sample_type_two_digit(raw) != "01":
                continue
            tid = normalize_tcga_id(raw)
            p = tid.participant
            if not p:
                continue
            out.setdefault(p, set()).add(mod)
    return out


def _all_tumor_participants(
    by_p: Dict[str, Set[str]],
    module_to_cols: Dict[str, List[str]],
) -> Set[str]:
    parts = set(by_p.keys())
    for cols in module_to_cols.values():
        for c in cols:
            if _is_primary_tumor_barcode(c):
                pp = _participant_from_tcga_col(c)
                if pp:
                    parts.add(pp)
    return parts


def _wide_presence_for_participant(
    participant: str,
    module_to_cols: Dict[str, List[str]],
) -> Set[str]:
    """Which wide-matrix modules have a TCGA column for this participant (tumor 01)."""
    have: Set[str] = set()
    for mod, cols in module_to_cols.items():
        for c in cols:
            if not _is_primary_tumor_barcode(c):
                continue
            if _participant_from_tcga_col(c) == participant:
                have.add(mod)
                break
    return have


def _pick_representative_tumor_vial(participant: str, vial_sets: Dict[str, Set[str]]) -> Optional[str]:
    """Prefer a tumor vial present in SNV, else any tumor 01 vial for this participant."""
    cands: List[str] = []
    for vs in vial_sets.values():
        for raw in vs:
            if tcga_sample_type_two_digit(raw) != "01":
                continue
            tid = normalize_tcga_id(raw)
            if tid.participant == participant and (tid.sample_vial or tid.sample):
                cands.append(str(tid.sample_vial or tid.sample))
    if not cands:
        return None
    snv = vial_sets.get("SNV", set())
    for v in sorted(set(cands)):
        if v in snv:
            return v
    return sorted(set(cands))[0]


def _manifest_tumor_map(tsv: Path) -> Dict[str, Tuple[str, str]]:
    """
    tumor_sample_vial -> (file_name, file_id)
    """
    smc = _load_sample_module_coverage()
    tumor_row = smc.tumor_sample_vial_from_gdc_manifest_row
    if not tsv.exists():
        return {}
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    out: Dict[str, Tuple[str, str]] = {}
    for _, row in df.iterrows():
        tv = tumor_row(row)
        if not tv:
            continue
        fn = str(row.get("File Name", "")).strip()
        fid = str(row.get("File ID", "")).strip()
        if fn:
            out[str(tv)] = (fn, fid)
    return out


def _snv_path(fn: str) -> Path:
    """Prefer VEP ``…APM_1Mb…`` under snv_vcf_dir; else raw manifest name under snv_raw_vcf_dir."""
    from pipeline.SNV.mutect_manifest_paths import (
        resolve_vep_mutect_vcf,
        vep_apm1mb_filename_candidates,
    )

    vep_roots = [
        PATHS.snv_vcf_dir,
        PATHS.working_dir / "SNV" / "vep_vcfs",
    ]
    hit = resolve_vep_mutect_vcf(fn, vep_dirs=vep_roots)
    if hit is not None:
        return hit
    raw = PATHS.snv_raw_vcf_dir / fn
    if raw.exists():
        return raw
    cands = vep_apm1mb_filename_candidates(fn)
    return PATHS.snv_vcf_dir / cands[0] if cands else PATHS.snv_vcf_dir / fn


def _sv_path(fn: str) -> Path:
    return PATHS.sv_vcf_dir / fn


def _cnv_path(fn: str) -> Path:
    return PATHS.cnv_dir / fn


def _meth_path(fn: str) -> Path:
    p = PATHS.methylation_samples_dir / fn
    if p.exists():
        return p
    return PATHS.methylation_output_dir / fn


def _rppa_path(fn: str) -> Path:
    for base in (
        getattr(PATHS, "rppa_samples_dir", PATHS.working_dir / "rppa" / "samples"),
        PATHS.working_dir / "rppa" / "samples",
        PATHS.working_dir / "rppa",
        PATHS.annotations_dir / "rppa",
    ):
        cand = base / fn
        if cand.exists():
            return cand
    return PATHS.working_dir / "rppa" / "samples" / fn


def _snv_variants_csv(vial: str) -> Path:
    """Pipeline-processed SNV table (per-sample), if present on disk."""
    snv_root = getattr(PATHS, "snv_output_dir", PATHS.working_dir / "snv_somatic_annotated")
    for sub in (
        snv_root / "per_sample" / vial,
        snv_root / vial,
        PATHS.working_dir / "SNV" / "pipeline_output" / vial,
        PATHS.working_dir / "SNV" / "per_sample_outputs" / vial,
    ):
        cand = sub / "snv_variants.csv"
        if cand.exists():
            return cand
    return Path()


def _sv_final_with_fimo_csv(vial: str) -> Path:
    """SV pipeline CSV for this tumor vial: prefer 09_chip_enriched, else 07_final_sv_with_fimo."""
    chip_d = getattr(PATHS, "sv_chip_enriched_dir", PATHS.sv_output_root / "09_chip_enriched")
    for d in (chip_d, PATHS.sv_output_root / "07_final_sv_with_fimo"):
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.csv")):
            name = p.name
            stem = p.stem
            if name.startswith(vial) or stem.startswith(vial) or vial in stem:
                return p
    return Path()


def _maybe_symlink(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if not src.exists():
        return
    try:
        if os.name == "posix":
            os.symlink(src.resolve(), dst, target_is_directory=False)
        else:
            shutil.copy2(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _subset_wide_tsv(
    in_path: Path,
    out_path: Path,
    id_col: str,
    sample_cols: List[str],
    *,
    sep: str = "\t",
) -> None:
    if not in_path.exists() or not sample_cols:
        return
    usecols = [id_col] + [c for c in sample_cols if c != id_col]
    df = pd.read_csv(in_path, sep=sep, usecols=lambda c: c in set(usecols), low_memory=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=sep, index=False)


def _subset_atac_csv(in_path: Path, out_path: Path, sample_cols: List[str]) -> None:
    if not in_path.exists() or not sample_cols:
        return
    hdr = pd.read_csv(in_path, nrows=0, low_memory=False)
    cols = list(hdr.columns)
    lead = [c for c in cols if not str(c).startswith("TCGA")]
    # Keep all non-TCGA leading columns (coordinates) + selected participants
    want = lead + [c for c in sample_cols if c in cols]
    df = pd.read_csv(in_path, usecols=want, low_memory=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10, help="Number of tumor participants to select")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: analysis/trial_max_coverage_cohort/run_<UTC>)",
    )
    args = ap.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out_dir = args.out or (ROOT / "analysis" / "trial_max_coverage_cohort" / f"run_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    subsets = out_dir / "subsets"
    per_sample = out_dir / "per_sample"
    subsets.mkdir(exist_ok=True)
    per_sample.mkdir(exist_ok=True)

    smc = _load_sample_module_coverage()
    vial_sets = smc.collect_vial_sets()
    participant_sets = smc.collect_participant_sets(vial_sets)

    by_p = _group_tumor_presence_by_participant(vial_sets)

    module_to_cols = {
        "RNA": _header_tcga_columns(PATHS.rna_expression)
        or _header_tcga_columns(PATHS.rna_expression_raw),
        "miRNA": _header_tcga_columns(PATHS.mirna_expression_tsv),
        "ATAC": _header_tcga_columns(PATHS.atac_case_level_matrix),
        "HLA": _header_tcga_columns(PATHS.hla_types_normalized)
        if PATHS.hla_types_normalized.exists()
        else _header_tcga_columns(PATHS.hla_samples_tsv),
    }

    hichip_hdr = _header_tcga_columns(PATHS.hichip_tcga_processed_csv)
    # HiCHIP columns are participants; treat as present if participant in header list
    hichip_parts = {normalize_tcga_id(c).participant for c in hichip_hdr if c.startswith("TCGA")}

    all_parts = _all_tumor_participants(by_p, module_to_cols)
    rows = []
    for p in sorted(all_parts):
        mods = by_p.get(p, set())
        core_hit = sum(1 for m in CORE_MODULES if m in mods)
        wide_have = _wide_presence_for_participant(p, module_to_cols)
        wide_hit = sum(1 for m in WIDE_MODULES if m in wide_have)
        hichip = int(p in hichip_parts or p in participant_sets.get("HiCHIP", set()))
        rows.append(
            {
                "participant": p,
                "n_core_01": core_hit,
                "n_wide_01": wide_hit,
                "hichip": hichip,
                "mods_vial_level": ",".join(sorted(mods)),
                "mods_wide_01": ",".join(sorted(wide_have)),
            }
        )

    df_rank = pd.DataFrame(rows)
    df_rank = df_rank.sort_values(
        ["n_core_01", "n_wide_01", "hichip", "participant"],
        ascending=[False, False, False, True],
    )
    df_rank.insert(0, "rank", range(1, len(df_rank) + 1))
    top = df_rank.head(args.n).copy()

    # Representative vial + resolved paths
    snv_m = _manifest_tumor_map(PATHS.annotations_dir / "SNV" / "samples.tsv")
    sv_m = _manifest_tumor_map(PATHS.annotations_dir / "SV" / "samples.tsv")
    cnv_m = _manifest_tumor_map(PATHS.annotations_dir / "CNV" / "samples.tsv")
    meth_df = pd.read_csv(PATHS.methylation_sample_manifest, sep="\t", low_memory=False) if PATHS.methylation_sample_manifest.exists() else pd.DataFrame()
    rppa_df = pd.read_csv(PATHS.annotations_dir / "rppa" / "samples.tsv", sep="\t", low_memory=False) if (PATHS.annotations_dir / "rppa" / "samples.tsv").exists() else pd.DataFrame()

    def meth_row_for_vial(vial: str) -> Optional[pd.Series]:
        if meth_df.empty or "Sample ID" not in meth_df.columns:
            return None
        hit = meth_df[meth_df["Sample ID"].astype(str).str.strip() == vial]
        if hit.empty:
            tid = normalize_tcga_id(vial)
            samp = tid.sample or vial
            hit = meth_df[meth_df["Sample ID"].astype(str).str.strip() == samp]
        if hit.empty:
            return None
        return hit.iloc[0]

    def rppa_row_for_vial(vial: str) -> Optional[pd.Series]:
        if rppa_df.empty or "Sample ID" not in rppa_df.columns:
            return None
        hit = rppa_df[rppa_df["Sample ID"].astype(str).str.strip() == vial]
        if hit.empty:
            tid = normalize_tcga_id(vial)
            samp = tid.sample or vial
            hit = rppa_df[rppa_df["Sample ID"].astype(str).str.strip() == samp]
        if hit.empty:
            return None
        return hit.iloc[0]

    long_rows: List[dict] = []
    picked_cols: Dict[str, List[str]] = {m: [] for m in WIDE_MODULES}

    rep_vial_by_p: Dict[str, str] = {}
    for _, r in top.iterrows():
        p = str(r["participant"])
        vial = _pick_representative_tumor_vial(p, vial_sets) or ""
        rep_vial_by_p[p] = vial

        # Wide column names (exact header strings) for subsetting
        for mod in WIDE_MODULES:
            for c in module_to_cols.get(mod, []):
                if _is_primary_tumor_barcode(c) and _participant_from_tcga_col(c) == p:
                    picked_cols[mod].append(c)
                    break

        if not vial:
            continue

        # Paths
        snv_t = snv_m.get(vial)
        sv_t = sv_m.get(vial)
        cnv_t = cnv_m.get(vial)
        paths = {
            "SNV": _snv_path(snv_t[0]) if snv_t else Path(),
            "SV": _sv_path(sv_t[0]) if sv_t else Path(),
            "CNV_seg": _cnv_path(cnv_t[0]) if cnv_t else Path(),
            "CNV_annotated": PATHS.cnv_output_dir / f"{vial}_cnv_annotated.csv",
        }
        paths["SNV_variants_csv"] = _snv_variants_csv(vial)
        paths["SV_final_csv"] = _sv_final_with_fimo_csv(vial)
        mr = meth_row_for_vial(vial)
        paths["Methylation"] = _meth_path(str(mr["File Name"])) if mr is not None else Path()
        rr = rppa_row_for_vial(vial)
        paths["RPPA"] = _rppa_path(str(rr["File Name"])) if rr is not None else Path()

        for mod, pt in paths.items():
            long_rows.append(
                {
                    "participant": p,
                    "representative_tumor_sample_vial": vial,
                    "table": mod,
                    "path": str(pt) if pt else "",
                    "exists": bool(pt and Path(pt).exists()),
                }
            )

        # Symlink small tables under per_sample/<participant>/
        psdir = per_sample / p.replace("-", "_")
        if rr is not None:
            rp = paths["RPPA"]
            if rp and Path(rp).exists():
                _maybe_symlink(Path(rp), psdir / "RPPA" / Path(rp).name)
        if mr is not None:
            mp = paths["Methylation"]
            if mp and Path(mp).exists():
                _maybe_symlink(Path(mp), psdir / "Methylation" / Path(mp).name)
        snv_csv = paths.get("SNV_variants_csv")
        if snv_csv and Path(snv_csv).exists():
            _maybe_symlink(Path(snv_csv), psdir / "SNV" / Path(snv_csv).name)
        sv_csv = paths.get("SV_final_csv")
        if sv_csv and Path(sv_csv).exists():
            _maybe_symlink(Path(sv_csv), psdir / "SV" / Path(sv_csv).name)
        cnv_ann = paths.get("CNV_annotated")
        if cnv_ann and Path(cnv_ann).exists():
            _maybe_symlink(Path(cnv_ann), psdir / "CNV" / Path(cnv_ann).name)

    top = top.copy()
    top["representative_tumor_sample_vial"] = top["participant"].map(rep_vial_by_p)
    top.to_csv(out_dir / "cohort_summary.tsv", sep="\t", index=False)
    paths_df = pd.DataFrame(long_rows)
    paths_df.to_csv(out_dir / "per_sample_input_paths.tsv", sep="\t", index=False)
    # One row per sample: key modalities for quick joins (includes processed tables when found).
    final_tbl = []
    for _, r in paths_df.iterrows():
        if r["table"] not in (
            "SNV",
            "SV",
            "SNV_variants_csv",
            "SV_final_csv",
            "CNV_annotated",
            "Methylation",
            "RPPA",
        ):
            continue
        final_tbl.append(dict(r))
    pd.DataFrame(final_tbl).to_csv(out_dir / "per_sample_final_tables.tsv", sep="\t", index=False)

    # Wide subsets (one combined list for RNA/miRNA column naming)
    rna_cols = picked_cols.get("RNA", [])
    if rna_cols and (PATHS.rna_expression.exists() or PATHS.rna_expression_raw.exists()):
        src = PATHS.rna_expression if PATHS.rna_expression.exists() else PATHS.rna_expression_raw
        sep = "\t" if str(src).endswith(".tsv") else ","
        idc = "sample" if sep == "\t" else str(pd.read_csv(src, sep=sep, nrows=0).columns[0])
        _subset_wide_tsv(src, subsets / f"{src.stem}_cohort{args.n}.tsv", idc, rna_cols, sep=sep)

    mir_cols = picked_cols.get("miRNA", [])
    if mir_cols and PATHS.mirna_expression_tsv.exists():
        hdr = pd.read_csv(PATHS.mirna_expression_tsv, sep="\t", nrows=0)
        idc = "sample" if "sample" in hdr.columns else str(hdr.columns[0])
        _subset_wide_tsv(PATHS.mirna_expression_tsv, subsets / f"mirna_cohort{args.n}.tsv", idc, mir_cols, sep="\t")

    atac_cols = picked_cols.get("ATAC", [])
    if atac_cols and PATHS.atac_case_level_matrix.exists():
        _subset_atac_csv(PATHS.atac_case_level_matrix, subsets / f"atac_case_level_cohort{args.n}.csv", atac_cols)

    # HiCHIP: subset participant columns (keep first 6 coord cols)
    if hichip_hdr and PATHS.hichip_tcga_processed_csv.exists():
        want_parts = [c for c in hichip_hdr if normalize_tcga_id(c).participant in set(top["participant"])]
        if want_parts:
            _subset_atac_csv(PATHS.hichip_tcga_processed_csv, subsets / f"hichip_TCGA_BRCA_cohort{args.n}.csv", want_parts)

    pd.DataFrame(
        [{"module": m, "column": c} for m, cols in picked_cols.items() for c in cols]
    ).to_csv(out_dir / "wide_matrix_columns_used.tsv", sep="\t", index=False)

    print(top.to_string(index=False))
    print(f"Wrote -> {out_dir}")


if __name__ == "__main__":
    main()
