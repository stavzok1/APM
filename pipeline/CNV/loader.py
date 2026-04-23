"""
CNV file loading and sample identification.

Reads raw CNV segment tables and resolves the real TCGA sample ID
via the GDC annotation manifest.
"""

import pandas as pd


def extract_sample_id_from_annotations(annotations, file_name):
    """Extract real sample ID using tumor-descriptor logic."""
    row = annotations[annotations["File Name"] == file_name]
    if row.empty:
        raise ValueError(f"No annotation row found for file_name={file_name!r}")

    row = row.iloc[0]
    tumor_desc = str(row["Tumor Descriptor"])
    sample_id_field = str(row["Sample ID"])
    parts = [p.strip() for p in sample_id_field.split(",")]

    if tumor_desc.startswith("Not"):
        if len(parts) < 2:
            raise ValueError(
                f"Expected 2 sample IDs in {sample_id_field!r} "
                f"for tumor_desc={tumor_desc!r}"
            )
        return parts[1]
    return parts[0]


def load_cnv_file(path, sep="\t"):
    """Read a single raw CNV segment file (ASCAT / GDC format)."""
    df = pd.read_csv(path, sep=sep)
    required = {"Chromosome", "Start", "End", "Copy_Number"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CNV file {path} missing required columns: {missing}")
    return df
