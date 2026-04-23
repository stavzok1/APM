#!/usr/bin/env python3
"""
Build annotations/RNA/samples.tsv, annotations/ATAC/samples.tsv, annotations/HLA/samples.tsv,
and annotations/HiCHIP/samples.tsv from wide matrices (RNA/ATAC/HiCHIP) and from HLA aliquot_id values.

Rows use the same core columns as GDC manifests under annotations/SNV|SV|CNV/ so downstream
tools can treat ``Sample ID`` as the primary tumor barcode (sample_vial when parseable).

Run from repo root:
    python scripts/annotations/build_annotation_sample_manifests.py

Requires: pandas (project venv).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import normalize_tcga_id  # noqa: E402


GDC_COLUMNS = [
    "File ID",
    "File Name",
    "Data Category",
    "Data Type",
    "Project ID",
    "Case ID",
    "Sample ID",
    "Tissue Type",
    "Tumor Descriptor",
    "Specimen Type",
    "Preservation Method",
    "source_column",
]


def _read_wide_matrix_header(matrix_path: Path) -> Optional[pd.DataFrame]:
    """
    Header only; try tab and comma. Prefer the parse with the most ``TCGA*`` columns
    (comma CSVs wrongly read as tab collapse into one column).
    """
    if not matrix_path.exists():
        return None
    best_hdr: Optional[pd.DataFrame] = None
    best_score = 0
    for sep in ("\t", ","):
        for kwargs in ({}, {"index_col": 0}):
            try:
                hdr = pd.read_csv(matrix_path, sep=sep, nrows=0, low_memory=False, **kwargs)
            except Exception:
                continue
            score = sum(1 for c in hdr.columns if str(c).strip().startswith("TCGA"))
            if score > best_score:
                best_score = score
                best_hdr = hdr
    return best_hdr


def _wide_matrix_tcga_headers(matrix_path: Path) -> List[str]:
    hdr = _read_wide_matrix_header(matrix_path)
    if hdr is None:
        return []
    out: List[str] = []
    for c in hdr.columns:
        s = str(c).strip()
        if s.startswith("TCGA"):
            out.append(s)
    return out


def _manifest_rows_from_headers(
    headers: List[str],
    *,
    source_file: Path,
    data_category: str,
    data_type: str,
    prefix: str,
    allow_participant: bool = False,
) -> pd.DataFrame:
    rows = []
    for i, raw_col in enumerate(headers):
        tid = normalize_tcga_id(raw_col)
        sid = tid.sample_vial or tid.sample or (tid.participant if allow_participant else None)
        if not sid or not str(sid).startswith("TCGA"):
            continue
        rows.append(
            {
                "File ID": f"{prefix}-{i:05d}",
                "File Name": source_file.name,
                "Data Category": data_category,
                "Data Type": data_type,
                "Project ID": "TCGA-BRCA",
                "Case ID": tid.participant or "",
                "Sample ID": sid,
                "Tissue Type": "Tumor",
                "Tumor Descriptor": "Primary",
                "Specimen Type": "Unknown",
                "Preservation Method": "Unknown",
                "source_column": raw_col,
            }
        )
    return pd.DataFrame(rows, columns=GDC_COLUMNS)


def _csv_headers_from_nth_col(path: Path, *, start_col_idx0: int) -> List[str]:
    """
    Read first line of comma CSV and return headers from a 0-based start index.
    HiChIP processed matrices in this repo have 6 coordinate columns, then TCGA columns.
    """
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
    cols = [c.strip() for c in header.split(",")]
    if len(cols) <= start_col_idx0:
        return []
    return cols[start_col_idx0:]


def build_hichip_manifest() -> pd.DataFrame:
    p = PATHS.hichip_tcga_processed_csv
    headers = _csv_headers_from_nth_col(p, start_col_idx0=6)
    headers = [h for h in headers if str(h).strip().startswith("TCGA")]
    if not headers:
        return pd.DataFrame(columns=GDC_COLUMNS)
    return _manifest_rows_from_headers(
        headers,
        source_file=p,
        data_category="Chromatin Interaction",
        data_type="HiChIP processed matrix",
        prefix="synthetic-HiCHIP",
        allow_participant=True,
    )


def build_mirna_manifest() -> pd.DataFrame:
    p = PATHS.mirna_expression_tsv
    headers = _wide_matrix_tcga_headers(p)
    if not headers:
        return pd.DataFrame(columns=GDC_COLUMNS)
    return _manifest_rows_from_headers(
        headers,
        source_file=p,
        data_category="Transcriptome Profiling",
        data_type="miRNA Expression Quantification (Xena arm-specific)",
        prefix="synthetic-miRNA",
    )


def build_rna_manifest() -> pd.DataFrame:
    headers: List[str] = []
    src: Optional[Path] = None
    for p in (PATHS.rna_expression, PATHS.rna_expression_raw):
        if p.exists():
            h = _wide_matrix_tcga_headers(p)
            if h:
                headers = h
                src = p
                break
    if not headers or src is None:
        return pd.DataFrame(columns=GDC_COLUMNS)
    return _manifest_rows_from_headers(
        headers,
        source_file=src,
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        prefix="synthetic-RNA",
    )


def build_atac_manifest() -> pd.DataFrame:
    p = PATHS.atac_case_level_matrix
    headers = _wide_matrix_tcga_headers(p)
    if not headers:
        return pd.DataFrame(columns=GDC_COLUMNS)
    return _manifest_rows_from_headers(
        headers,
        source_file=p,
        data_category="ATAC-seq",
        data_type="Chromatin Accessibility Case Matrix",
        prefix="synthetic-ATAC",
    )


def build_hla_manifest() -> pd.DataFrame:
    p = PATHS.hla_types_tsv
    if not p.exists():
        return pd.DataFrame(columns=GDC_COLUMNS)
    df = pd.read_csv(p, sep=",", low_memory=False)
    if "aliquot_id" not in df.columns:
        return pd.DataFrame(columns=GDC_COLUMNS)
    uniq = df["aliquot_id"].dropna().astype(str).unique()
    rows = []
    for i, raw in enumerate(sorted(uniq)):
        tid = normalize_tcga_id(raw)
        sid = tid.sample_vial or tid.sample
        if not sid or not str(sid).startswith("TCGA"):
            continue
        rows.append(
            {
                "File ID": f"synthetic-HLA-{i:05d}",
                "File Name": p.name,
                "Data Category": "Biospecimen",
                "Data Type": "HLA Types",
                "Project ID": "TCGA-BRCA",
                "Case ID": tid.participant or "",
                "Sample ID": sid,
                "Tissue Type": "Tumor",
                "Tumor Descriptor": "Primary",
                "Specimen Type": "Unknown",
                "Preservation Method": "Unknown",
                "source_column": raw,
            }
        )
    return pd.DataFrame(rows, columns=GDC_COLUMNS)


def main() -> None:
    ap = argparse.ArgumentParser(description="Write RNA/ATAC/HLA/HiCHIP/miRNA samples.tsv under annotations/.")
    ap.add_argument("--dry-run", action="store_true", help="Print counts only; do not write files.")
    args = ap.parse_args()

    rna_df = build_rna_manifest()
    atac_df = build_atac_manifest()
    hla_df = build_hla_manifest()
    hichip_df = build_hichip_manifest()
    mirna_df = build_mirna_manifest()

    print(f"RNA rows: {len(rna_df)} (from matrices under PATHS.rna_expression / rna_expression_raw)")
    print(f"ATAC rows: {len(atac_df)} (from PATHS.atac_case_level_matrix)")
    print(f"HLA rows: {len(hla_df)} (unique aliquot_id from PATHS.hla_types_tsv)")
    print(f"HiCHIP rows: {len(hichip_df)} (from PATHS.hichip_tcga_processed_csv columns starting at 7th col)")
    print(f"miRNA rows: {len(mirna_df)} (from PATHS.mirna_expression_tsv columns)")

    if args.dry_run:
        return

    out_rna = PATHS.rna_samples_tsv
    out_atac = PATHS.atac_samples_tsv
    out_hla = PATHS.hla_samples_tsv
    out_hichip = PATHS.hichip_samples_tsv
    out_mirna = PATHS.mirna_samples_tsv
    for p in (out_rna.parent, out_atac.parent, out_hla.parent, out_hichip.parent, out_mirna.parent):
        p.mkdir(parents=True, exist_ok=True)

    if len(rna_df):
        rna_df.to_csv(out_rna, sep="\t", index=False)
        print(f"Wrote {out_rna}")
    else:
        print(f"Skipped {out_rna} (no TCGA columns found in RNA matrices)")

    if len(atac_df):
        atac_df.to_csv(out_atac, sep="\t", index=False)
        print(f"Wrote {out_atac}")
    else:
        print(f"Skipped {out_atac} (matrix missing or no TCGA columns): {PATHS.atac_case_level_matrix}")

    if len(hla_df):
        hla_df.to_csv(out_hla, sep="\t", index=False)
        print(f"Wrote {out_hla}")
    else:
        print(f"Skipped {out_hla} (HLA table missing or empty): {PATHS.hla_types_tsv}")

    if len(hichip_df):
        hichip_df.to_csv(out_hichip, sep="\t", index=False)
        print(f"Wrote {out_hichip}")
    else:
        print(f"Skipped {out_hichip} (HiCHIP matrix missing or no TCGA columns): {PATHS.hichip_tcga_processed_csv}")

    if len(mirna_df):
        mirna_df.to_csv(out_mirna, sep="\t", index=False)
        print(f"Wrote {out_mirna}")
    else:
        print(f"Skipped {out_mirna} (miRNA matrix missing or no TCGA columns): {PATHS.mirna_expression_tsv}")


if __name__ == "__main__":
    main()
