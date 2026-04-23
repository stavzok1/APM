from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class DdrScoresProvider(CovariateProvider):
    """
    Provider for `data/covariates/DDRscores.txt`.

    This is treated as an external-but-canonical cohort covariate table:
    it includes HRD metrics, mutational signature exposures, purity/ploidy,
    and TP53-related scores.
    """

    name = "ddr_scores"

    def __init__(
        self,
        *,
        path: Path,
        cancer_acronym: str = "BRCA",
    ) -> None:
        self.path = Path(path)
        self.cancer_acronym = str(cancer_acronym)

    def build(self) -> ProviderResult:
        if not self.path.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample",
                diagnostics={"error": f"missing {self.path}", **_diagnose_index(df)},
            )

        # File is TSV with quoted header names.
        df = pd.read_csv(self.path, sep="\t", low_memory=False)
        # Unquote column names if needed.
        df.columns = [str(c).strip().strip('"') for c in df.columns]

        if "acronym" in df.columns:
            df = df[df["acronym"].astype(str) == self.cancer_acronym].copy()

        id_col = "patient_id" if "patient_id" in df.columns else None
        if id_col is None:
            out = pd.DataFrame(index=pd.Index([], name="sample"))
            return ProviderResult(
                name=self.name,
                df=out,
                index_key="sample",
                diagnostics={"error": "missing patient_id column", **_diagnose_index(out)},
            )

        # patient_id entries look like TCGA-..-..-01 (sample-level).
        add_tcga_id_columns_inplace(df, raw_id_col=id_col, overwrite=False)

        out = df.copy()
        out["sample"] = out["sample"].astype("string")
        out.index = out["sample"].where(out["sample"].notna(), out[id_col].astype(str)).astype(str)
        out.index.name = "sample"

        # Keep numeric columns as-is; caller prefixes columns.
        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample",
            diagnostics={
                "path": str(self.path),
                "acronym": self.cancer_acronym,
                **_diagnose_index(out),
            },
        )

