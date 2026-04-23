#!/usr/bin/env python3
"""
Harmonize TCGA ATAC raw counts with healthy (normal) ATAC vectors.

- Loads TCGA raw ATAC counts matrix (tsv)
- Filters TCGA samples by _primary_disease using TCGA samples table + ATAC identifier mapping
- Appends healthy sample columns from a folder of *.vector files
- Writes a combined raw counts matrix (tsv)

This script is designed to receive PATHS as an input (dict-like or object with attributes).




USAGE:

I. from another file (programmatically):
from atac.harmonize_atac_counts import harmonize_atac_counts

PATHS = {
    "TCGA_RAW_ATAC_COUNTS_PATH": "data/TCGA_ATAC/TCGA-ATAC_PanCan_Raw_Counts.txt",
    "HEALTHY_COUNTS_DIR": "data/TCGA_ATAC/normal_counts",
    "TCGA_SAMPLES_TABLE_PATH": "data/TCGA.tsv",
    "ATAC_META_PATH": "data/TCGA_ATAC/TCGA_identifier_mapping.txt",
    "OUT_PATH": "data/TCGA_ATAC/TCGA_BRCA_ATAC_with_normals_raw.txt",
}

df = harmonize_atac_counts(
    paths=PATHS,
    first_sample_col=7,
    primary_disease="breast invasive carcinoma",
    tcga_cancer_type="BRCA",
)



II. from command line:
1. first create a json file with the paths
{
  "TCGA_RAW_ATAC_COUNTS_PATH": "data/TCGA_ATAC/TCGA-ATAC_PanCan_Raw_Counts.txt",
  "HEALTHY_COUNTS_DIR": "data/TCGA_ATAC/normal_counts",
  "TCGA_SAMPLES_TABLE_PATH": "data/TCGA.tsv",
  "ATAC_META_PATH": "data/TCGA_ATAC/TCGA_identifier_mapping.txt",
  "OUT_PATH": "data/TCGA_ATAC/TCGA_BRCA_ATAC_with_normals_raw.txt"
}

2. then run the script from the command line:
python harmonize_atac_counts.py \
  --paths-json atac_paths.json \
  --first-sample-col 7 \
  --primary-disease "breast invasive carcinoma" \
  --tcga-cancer-type BRCA

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Union

import os
import pandas as pd


# -------------------------
# PATHS handling
# -------------------------
PathLike = Union[str, Path]


@dataclass(frozen=True)
class AtacPaths:
    tcga_raw_atac_counts: Path
    healthy_counts_dir: Path
    tcga_samples_table: Path
    atac_meta: Path
    out_file: Path

    @staticmethod
    def from_paths(paths: Any, *, out_file: Optional[PathLike] = None) -> "AtacPaths":
        """
        Build AtacPaths from a PATHS input.

        Supported PATHS formats:
        1) dict-like with keys:
           - TCGA_RAW_ATAC_COUNTS_PATH
           - HEALTHY_COUNTS_DIR
           - TCGA_SAMPLES_TABLE_PATH
           - ATAC_META_PATH
           - OUT_PATH (optional; can be overridden by out_file)
        2) object-like with the same attribute names.

        If out_file is provided, it overrides OUT_PATH.
        """
        def _get(name: str) -> Any:
            if isinstance(paths, Mapping) and name in paths:
                return paths[name]
            if hasattr(paths, name):
                return getattr(paths, name)
            return None

        tcga_raw = _get("TCGA_RAW_ATAC_COUNTS_PATH")
        healthy_dir = _get("HEALTHY_COUNTS_DIR")
        tcga_tbl = _get("TCGA_SAMPLES_TABLE_PATH")
        atac_meta = _get("ATAC_META_PATH")

        if tcga_raw is None or healthy_dir is None or tcga_tbl is None or atac_meta is None:
            raise ValueError(
                "PATHS must provide TCGA_RAW_ATAC_COUNTS_PATH, HEALTHY_COUNTS_DIR, "
                "TCGA_SAMPLES_TABLE_PATH, and ATAC_META_PATH (as dict keys or attributes)."
            )

        out_from_paths = _get("OUT_PATH")

        chosen_out = out_file if out_file is not None else out_from_paths
        if chosen_out is None:
            raise ValueError(
                "Output path not provided. Provide out_file=... or PATHS['OUT_PATH'] / PATHS.OUT_PATH."
            )

        return AtacPaths(
            tcga_raw_atac_counts=Path(tcga_raw),
            healthy_counts_dir=Path(healthy_dir),
            tcga_samples_table=Path(tcga_tbl),
            atac_meta=Path(atac_meta),
            out_file=Path(chosen_out),
        )


# -------------------------
# Core logic
# -------------------------
def filter_tcga_samples(
    *,
    counts: pd.DataFrame,
    tcga_samples_table: pd.DataFrame,
    atac_meta: pd.DataFrame,
    first_sample_col: int,
    primary_disease: str,
) -> pd.DataFrame:
    """
    Keep only TCGA ATAC sample columns belonging to cases whose _primary_disease matches primary_disease.

    counts:
      TCGA ATAC raw counts table, with metadata columns first, then sample columns.
    tcga_samples_table:
      must include columns: "_primary_disease", "sample"
    atac_meta:
      must include columns: "Case_ID", "bam_prefix"
    first_sample_col:
      index (0-based) of first sample column in counts
    """
    # TCGA sample case IDs (short: first 3 fields)
    tcga_samples = (
        tcga_samples_table.loc[tcga_samples_table["_primary_disease"] == primary_disease, "sample"]
        .astype(str)
        .tolist()
    )
    tcga_samples_short = ["-".join(s.split("-")[:3]) for s in tcga_samples]
    print(f"Found {len(tcga_samples_short)} TCGA cases matching _primary_disease='{primary_disease}'")

    atac_meta = atac_meta.copy()
    atac_meta["Case_ID_short"] = atac_meta["Case_ID"].astype(str).apply(lambda x: "-".join(x.split("-")[:3]))

    # bam_prefix in mapping, normalize to match counts column naming (replace '-' with '_')
    atac_bam_prefixes = atac_meta.loc[atac_meta["Case_ID_short"].isin(tcga_samples_short), "bam_prefix"].astype(str)
    atac_bam_prefixes = [s.replace("-", "_") for s in atac_bam_prefixes.tolist()]

    # Keep only columns present in counts to avoid KeyError
    present = [c for c in atac_bam_prefixes if c in counts.columns]
    missing = [c for c in atac_bam_prefixes if c not in counts.columns]

    print(f"ATAC mapped sample columns present in counts: {len(present)}")
    if missing:
        print(f"Warning: {len(missing)} mapped bam_prefix columns not found in counts (showing up to 10): {missing[:10]}")

    meta_part = counts.iloc[:, :first_sample_col].copy()
    sample_part = counts.loc[:, present].copy()
    return pd.concat([meta_part, sample_part], axis=1)


def add_healthy_counts(*, counts: pd.DataFrame, vector_folder: PathLike) -> pd.DataFrame:
    """
    Add healthy sample columns from a folder containing *.vector files.
    Each .vector is assumed to be a single-column list of counts, with no header.
    """
    vector_folder = Path(vector_folder)
    if not vector_folder.exists():
        raise FileNotFoundError(f"Healthy counts dir not found: {vector_folder}")
    if not vector_folder.is_dir():
        raise NotADirectoryError(f"Healthy counts dir is not a directory: {vector_folder}")

    counts = counts.copy()
    n_rows_expected = counts.shape[0]

    added = 0
    for fpath in sorted(vector_folder.iterdir()):
        if fpath.is_file() and fpath.name.endswith(".vector"):
            sample_name = fpath.name[:-7]  # remove ".vector"
            vec = pd.read_csv(fpath, header=None).iloc[:, 0]

            if len(vec) != n_rows_expected:
                raise ValueError(
                    f"Row mismatch in {fpath.name}: expected {n_rows_expected}, got {len(vec)}"
                )

            counts[sample_name] = vec.values
            added += 1
            print(f"Added healthy sample: {sample_name}")

    print(f"Total healthy samples added: {added}")
    return counts


def harmonize_atac_counts(
    *,
    paths: Any,
    first_sample_col: int,
    primary_disease: str,
    tcga_cancer_type: str,
    out_file: Optional[PathLike] = None,
) -> pd.DataFrame:
    """
    Main entry point (callable from notebooks/pipeline).

    paths: dict-like or object-like PATHS container (see AtacPaths.from_paths)
    tcga_cancer_type: used for default naming conventions / logs (no filtering logic here)
    out_file: optional override for output filepath
    """
    ap = AtacPaths.from_paths(paths, out_file=out_file)

    # Load inputs
    tcga_raw = pd.read_csv(ap.tcga_raw_atac_counts, sep="\t")
    tcga_samples_table = pd.read_csv(ap.tcga_samples_table, sep="\t")
    atac_meta = pd.read_csv(ap.atac_meta, sep="\t")

    # Filter TCGA samples by disease
    tcga_filtered = filter_tcga_samples(
        counts=tcga_raw,
        tcga_samples_table=tcga_samples_table,
        atac_meta=atac_meta,
        first_sample_col=first_sample_col,
        primary_disease=primary_disease,
    )

    # Add healthy counts
    all_samples = add_healthy_counts(counts=tcga_filtered, vector_folder=ap.healthy_counts_dir)

    # Write output
    ap.out_file.parent.mkdir(parents=True, exist_ok=True)
    all_samples.to_csv(ap.out_file, sep="\t", index=False)
    print(f"Wrote combined raw ATAC counts for {tcga_cancer_type} to: {ap.out_file}")

    return all_samples


# -------------------------
# Optional CLI (nice to have)
# -------------------------
if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Create a combined TCGA+healthy ATAC raw counts matrix for a given disease."
    )
    parser.add_argument(
        "--paths-json",
        required=True,
        help=(
            "JSON string or path to JSON file containing PATHS with keys: "
            "TCGA_RAW_ATAC_COUNTS_PATH, HEALTHY_COUNTS_DIR, TCGA_SAMPLES_TABLE_PATH, ATAC_META_PATH, "
            "and either OUT_PATH or provide --out-file."
        ),
    )
    parser.add_argument("--out-file", default=None, help="Optional override output path.")
    parser.add_argument("--first-sample-col", type=int, default=7, help="0-based first sample column index.")
    parser.add_argument("--primary-disease", required=True, help="Value of _primary_disease to keep.")
    parser.add_argument("--tcga-cancer-type", required=True, help="Label used for logs / naming.")

    args = parser.parse_args()

    # Load PATHS dict from JSON string or JSON file
    pj = args.paths_json
    if os.path.exists(pj):
        with open(pj, "r", encoding="utf-8") as f:
            paths_dict = json.load(f)
    else:
        paths_dict = json.loads(pj)

    harmonize_atac_counts(
        paths=paths_dict,
        first_sample_col=args.first_sample_col,
        primary_disease=args.primary_disease,
        tcga_cancer_type=args.tcga_cancer_type,
        out_file=args.out_file,
    )
