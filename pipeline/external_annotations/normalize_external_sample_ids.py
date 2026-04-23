from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ..sample_ids import add_tcga_id_columns_inplace, normalize_tcga_id


def normalize_immune_subtype_annotations(
    path: Path,
    *,
    id_col: str = "TCGA Participant Barcode",
) -> pd.DataFrame:
    """
    Loads BRCA immune subtype annotations and appends normalized TCGA keys.

    Default file (config): annotations/BRCA_immune_subtypes_advanced.tsv
    ID column observed: "TCGA Participant Barcode" (typically TCGA-..-..-01).
    """
    df = pd.read_csv(path, sep="\t")
    if id_col in df.columns:
        add_tcga_id_columns_inplace(df, raw_id_col=id_col)
    return df


def normalize_thorsson_immune_table(
    path: Path,
    *,
    id_col: str = "sample_id",
) -> pd.DataFrame:
    """
    Loads the Thorsson/Thornsson immune table and appends normalized TCGA keys.

    Default file (config): annotations/Thornsson_immune_table.tsv
    ID column observed: "sample_id" (TCGA-..-..-01).
    """
    df = pd.read_csv(path, sep="\t")
    if id_col in df.columns:
        add_tcga_id_columns_inplace(df, raw_id_col=id_col)
    return df


def normalize_hla_types(
    path: Path,
    *,
    id_col: str = "aliquot_id",
) -> pd.DataFrame:
    """
    Loads HLA calls and appends normalized TCGA keys.

    File: annotations/HLA_types.tsv
    ID column observed: "aliquot_id" (full aliquot barcode).
    """
    # HLA_types.tsv in this project is comma-separated despite the .tsv suffix.
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
    sep = "," if header.count(",") > header.count("\t") else "\t"
    df = pd.read_csv(path, sep=sep)
    if id_col in df.columns:
        add_tcga_id_columns_inplace(df, raw_id_col=id_col)
    return df


def build_matrix_sample_metadata(sample_ids, *, source_name: Optional[str] = None) -> pd.DataFrame:
    """
    For wide matrices where samples are column headers (RNA/ATAC),
    create a tidy sample metadata table with normalized TCGA keys.
    """
    rows = []
    for s in sample_ids:
        ids = normalize_tcga_id(str(s))
        rows.append(
            {
                "raw_sample_id": ids.raw,
                "participant": ids.participant,
                "sample": ids.sample,
                "sample_vial": ids.sample_vial,
                "aliquot": ids.aliquot,
                "source": source_name,
            }
        )
    return pd.DataFrame(rows)

