#!/usr/bin/env python3
"""
Collapse TCGA ATAC technical replicates to case-level by averaging
quantile-normalized log2-CPM values.

- TCGA replicates → averaged per Case_ID
- Non-TCGA samples (e.g. normals) are left untouched

USAGE:
I. from another file (programmatically):
from atac.collapse_atac_replicates import run_collapse_atac_replicates

PATHS = {
    "ATAC_LOGCPM_QN_PATH": "data/TCGA_ATAC/TCGA_BRCA_ATAC_with_normals_logCPM_QN.txt",
    "ATAC_META_PATH": "data/TCGA_ATAC/TCGA_identifier_mapping.txt",
    "OUT_PATH": "data/TCGA_ATAC/TCGA_BRCA_ATAC_case_level_logCPM_QN.txt",
}

collapsed = run_collapse_atac_replicates(
    paths=PATHS,
    first_sample_col=7,
)

II. from command line:
python collapse_atac_replicates.py \
  --paths-json atac_paths.json \
  --first-sample-col 7
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Union

import pandas as pd


PathLike = Union[str, Path]


@dataclass(frozen=True)
class AtacCollapsePaths:
    input_matrix: Path
    atac_meta: Path
    out_file: Path

    @staticmethod
    def from_paths(paths: Any, *, out_file: Optional[PathLike] = None) -> "AtacCollapsePaths":
        def _get(name: str):
            if isinstance(paths, Mapping) and name in paths:
                return paths[name]
            if hasattr(paths, name):
                return getattr(paths, name)
            return None

        inp = _get("ATAC_LOGCPM_QN_PATH")
        meta = _get("ATAC_META_PATH")
        outp = out_file if out_file is not None else _get("OUT_PATH")

        if inp is None or meta is None or outp is None:
            raise ValueError(
                "PATHS must contain ATAC_LOGCPM_QN_PATH, ATAC_META_PATH, and OUT_PATH"
            )

        return AtacCollapsePaths(
            input_matrix=Path(inp),
            atac_meta=Path(meta),
            out_file=Path(outp),
        )


def collapse_atac_replicates(
    *,
    log_cpm_qn: pd.DataFrame,
    atac_meta: pd.DataFrame,
    first_sample_col: int,
    bam_prefix_col: str = "bam_prefix",
    case_id_col: str = "Case_ID",
) -> pd.DataFrame:
    """
    Collapse ATAC technical replicates to TCGA case-level.
    """

    meta_qn = log_cpm_qn.iloc[:, :first_sample_col]
    expr_qn = log_cpm_qn.iloc[:, first_sample_col:]

    atac_meta = atac_meta.copy()
    atac_meta["bam_prefix_norm"] = atac_meta[bam_prefix_col].str.replace("-", "_")

    atac_meta = atac_meta[atac_meta["bam_prefix_norm"].isin(expr_qn.columns)]

    col_to_case = dict(zip(atac_meta["bam_prefix_norm"], atac_meta[case_id_col]))

    group_keys = [col_to_case.get(c, c) for c in expr_qn.columns]

    collapsed_expr = expr_qn.T.groupby(group_keys).mean().T

    return pd.concat([meta_qn, collapsed_expr], axis=1)


def run_collapse_atac_replicates(
    *,
    paths: Any,
    first_sample_col: int,
    out_file: Optional[PathLike] = None,
) -> pd.DataFrame:
    ap = AtacCollapsePaths.from_paths(paths, out_file=out_file)

    log_cpm_qn = pd.read_csv(ap.input_matrix, sep="\t")
    atac_meta = pd.read_csv(ap.atac_meta, sep="\t")

    collapsed = collapse_atac_replicates(
        log_cpm_qn=log_cpm_qn,
        atac_meta=atac_meta,
        first_sample_col=first_sample_col,
    )

    ap.out_file.parent.mkdir(parents=True, exist_ok=True)
    collapsed.to_csv(ap.out_file, index=False)

    print(f"Wrote case-level ATAC matrix to: {ap.out_file}")
    return collapsed


if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(
        description="Collapse TCGA ATAC technical replicates to case level"
    )
    parser.add_argument("--paths-json", required=True)
    parser.add_argument("--first-sample-col", type=int, default=7)
    parser.add_argument("--out-file", default=None)

    args = parser.parse_args()

    if os.path.exists(args.paths_json):
        with open(args.paths_json) as f:
            paths = json.load(f)
    else:
        paths = json.loads(args.paths_json)

    run_collapse_atac_replicates(
        paths=paths,
        first_sample_col=args.first_sample_col,
        out_file=args.out_file,
    )
