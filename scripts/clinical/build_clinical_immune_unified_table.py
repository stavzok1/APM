from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import add_tcga_id_columns_inplace  # noqa: E402


def main() -> Path:
    out_path = Path(PATHS.annotations_dir) / "BRCA_clinical_immune_unified.tsv"

    clinical = pd.read_csv(PATHS.brca_clinical, sep="\t")
    immune = pd.read_csv(PATHS.immune_subtype_annotations, sep="\t")

    clinical_keep = clinical[
        [
            "sampleID",
            "pathologic_stage",
            "pathologic_T",
            "pathologic_N",
            "pathologic_M",
        ]
    ].copy()
    clinical_keep = clinical_keep.rename(columns={"sampleID": "sample_id"})

    immune_keep = immune[["sample_id", "PAM50_final", "CPE"]].copy()

    merged = clinical_keep.merge(immune_keep, on="sample_id", how="left")
    add_tcga_id_columns_inplace(merged, raw_id_col="sample_id")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    return out_path


if __name__ == "__main__":
    p = main()
    print(str(p))

