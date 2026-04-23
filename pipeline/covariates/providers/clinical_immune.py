from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class ClinicalImmuneProvider(CovariateProvider):
    name = "clinical_immune"

    def __init__(self, *, unified_tsv: Path) -> None:
        self.unified_tsv = Path(unified_tsv)

    def build(self) -> ProviderResult:
        if not self.unified_tsv.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample",
                diagnostics={"error": f"missing {self.unified_tsv}", **_diagnose_index(df)},
            )

        df = pd.read_csv(self.unified_tsv, sep="\t", low_memory=False)
        if "sample_id" in df.columns:
            raw_col = "sample_id"
        elif "sample" in df.columns:
            raw_col = "sample"
        else:
            raw_col = None

        if raw_col:
            add_tcga_id_columns_inplace(df, raw_id_col=raw_col, overwrite=False)

        # Prefer 'sample' as index; fall back to raw.
        idx = df.get("sample")
        if idx is None:
            idx = df[raw_col].astype(str) if raw_col else pd.Series([], dtype="string")

        out = df.copy()
        out.index = idx.astype(str)
        out.index.name = "sample"

        # Drop the raw index column if it duplicates the index.
        for c in ("sample_id",):
            if c in out.columns:
                # keep the explicit sample_id column if user wants; but avoid duplication explosion
                pass

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample",
            diagnostics={"unified_tsv": str(self.unified_tsv), **_diagnose_index(out)},
        )

