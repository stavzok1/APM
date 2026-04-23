from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS
from pipeline.external_annotations import (
    build_matrix_sample_metadata,
    normalize_hla_types,
    normalize_immune_subtype_annotations,
    normalize_thorsson_immune_table,
)


def _read_matrix_header_columns(path: Path, *, sep: str = "\t") -> list[str]:
    """
    Read only the first line of a delimited text file and return column names.
    Works for very wide matrices.
    """
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
    return header.split(sep)


def main() -> None:
    out_dir = Path(PATHS.annotations_dir) / "_normalized"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Immune tables (tidy)
    immune_adv = normalize_immune_subtype_annotations(PATHS.immune_subtype_annotations)
    immune_adv.to_csv(out_dir / "BRCA_immune_subtypes_advanced.normalized.tsv", sep="\t", index=False)

    thorsson = normalize_thorsson_immune_table(PATHS.thornsson_immune_table)
    thorsson.to_csv(out_dir / "Thornsson_immune_table.normalized.tsv", sep="\t", index=False)

    # HLA (tidy; comma-separated file lives at PATHS.hla_types_tsv)
    hla = normalize_hla_types(Path(PATHS.hla_types_tsv), id_col="aliquot_id")
    hla.to_csv(out_dir / "HLA_types.normalized.tsv", sep="\t", index=False)

    # RNA matrix sample metadata (columns = samples)
    # Treat only headers that look like TCGA barcodes as samples.
    rna_cols = _read_matrix_header_columns(PATHS.rna_expression_raw, sep="\t")
    rna_samples = [c for c in rna_cols if str(c).startswith("TCGA-")]
    rna_meta = build_matrix_sample_metadata(rna_samples, source_name="RNAexp_TCGA")
    rna_meta.to_csv(out_dir / "RNA_expression.sample_metadata.tsv", sep="\t", index=False)

    # ATAC matrix sample metadata (columns = samples)
    # Treat only headers that look like TCGA barcodes as samples.
    atac_cols = _read_matrix_header_columns(PATHS.atac_case_level_matrix, sep=",")
    atac_samples = [c for c in atac_cols if str(c).startswith("TCGA-")]
    atac_meta = build_matrix_sample_metadata(atac_samples, source_name="TCGA_ATAC_case_level")
    atac_meta.to_csv(out_dir / "ATAC_case_level.sample_metadata.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

