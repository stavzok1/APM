#!/usr/bin/env python3
"""
Count TCGA sample types per module: primary tumor (01), blood-derived normal (10),
solid-tissue normal (11), and other/unknown.

Interpretation:
- For vial/sample-keyed modules, classify from the two-digit sample type in the barcode.
- For participant-only modules, assume primary tumor (01) as requested.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import count_unique_tcga_sample_types, normalize_tcga_id  # noqa: E402


def _read_manifest_sample_ids(tsv: Path) -> List[str]:
    if not tsv.exists():
        return []
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if "Sample ID" not in df.columns:
        return []
    out: List[str] = []
    for raw in df["Sample ID"].dropna().astype(str):
        out.extend([p.strip() for p in raw.split(",") if p.strip()])
    return out


def _read_wide_matrix_tcga_headers(path: Path) -> List[str]:
    if not path.exists():
        return []
    best: List[str] = []
    best_score = -1
    for sep in ("\t", ","):
        for kwargs in ({}, {"index_col": 0}):
            try:
                hdr = pd.read_csv(path, sep=sep, nrows=0, low_memory=False, **kwargs)
            except Exception:
                continue
            cols = [str(c).strip() for c in hdr.columns]
            tcga = [c for c in cols if c.startswith("TCGA")]
            if len(tcga) > best_score:
                best_score = len(tcga)
                best = tcga
    return best


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Primary (01) vs blood normal (10) vs solid normal (11) per module."
    )
    ap.add_argument("--out", type=Path, required=True, help="Output TSV path")
    args = ap.parse_args()

    ann = PATHS.annotations_dir

    module_to_ids: Dict[str, List[str]] = {
        "SNV": _read_manifest_sample_ids(ann / "SNV" / "samples.tsv"),
        "SV": _read_manifest_sample_ids(ann / "SV" / "samples.tsv"),
        "CNV": _read_manifest_sample_ids(ann / "CNV" / "samples.tsv"),
        "Methylation": _read_manifest_sample_ids(PATHS.methylation_sample_manifest),
        "RPPA": _read_manifest_sample_ids(ann / "rppa" / "samples.tsv"),
        "RNA": _read_manifest_sample_ids(PATHS.rna_samples_tsv)
        or _read_wide_matrix_tcga_headers(PATHS.rna_expression),
        "ATAC": _read_manifest_sample_ids(PATHS.atac_samples_tsv)
        or _read_wide_matrix_tcga_headers(PATHS.atac_case_level_matrix),
        "HLA": _read_manifest_sample_ids(PATHS.hla_samples_tsv),
        "miRNA": _read_manifest_sample_ids(PATHS.mirna_samples_tsv)
        or _read_wide_matrix_tcga_headers(PATHS.mirna_expression_tsv),
    }

    rows = []
    for mod, ids in module_to_ids.items():
        n_total, n01, n10, n11, nother = count_unique_tcga_sample_types(ids)
        n_norm = n10 + n11
        rows.append(
            {
                "module": mod,
                "n_unique_samples_total": n_total,
                "n_primary_tumor_01": n01,
                "n_normal_blood_derived_10": n10,
                "n_normal_solid_tissue_11": n11,
                "n_normal_total_10_plus_11": n_norm,
                "n_other_or_unknown": nother,
                "pct_primary_01": round(100.0 * n01 / max(n_total, 1), 3) if n_total else 0.0,
                "pct_normal_blood_10_of_normals": round(100.0 * n10 / max(n_norm, 1), 3) if n_norm else 0.0,
                "pct_normal_solid_11_of_normals": round(100.0 * n11 / max(n_norm, 1), 3) if n_norm else 0.0,
                "notes": "",
            }
        )

    participant_only = {
        "HiCHIP": _read_manifest_sample_ids(PATHS.hichip_samples_tsv),
        "Immune_advanced": [
            str(x)
            for x in pd.read_csv(PATHS.immune_subtype_annotations_normalized, sep="\t", low_memory=False)
            .iloc[:, 0]
            .dropna()
            .astype(str)
            .tolist()
        ]
        if PATHS.immune_subtype_annotations_normalized.exists()
        else [],
        "Immune_thornsson": [
            str(x)
            for x in pd.read_csv(PATHS.thornsson_immune_table_normalized, sep="\t", low_memory=False)
            .iloc[:, 0]
            .dropna()
            .astype(str)
            .tolist()
        ]
        if PATHS.thornsson_immune_table_normalized.exists()
        else [],
        "Clinical_unified": [
            str(x)
            for x in pd.read_csv(PATHS.brca_clinical_immune_unified, sep="\t", low_memory=False)["participant"]
            .dropna()
            .astype(str)
            .tolist()
        ]
        if PATHS.brca_clinical_immune_unified.exists()
        else [],
    }
    for mod, ids in participant_only.items():
        uniq_p = {normalize_tcga_id(x).participant for x in ids if normalize_tcga_id(x).participant}
        n_p = len(uniq_p)
        rows.append(
            {
                "module": mod,
                "n_unique_samples_total": n_p,
                "n_primary_tumor_01": n_p,
                "n_normal_blood_derived_10": 0,
                "n_normal_solid_tissue_11": 0,
                "n_normal_total_10_plus_11": 0,
                "n_other_or_unknown": 0,
                "pct_primary_01": 100.0 if n_p else 0.0,
                "pct_normal_blood_10_of_normals": 0.0,
                "pct_normal_solid_11_of_normals": 0.0,
                "notes": "participant_only_assumed_primary_tumor",
            }
        )

    out = pd.DataFrame(rows).sort_values("module")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
