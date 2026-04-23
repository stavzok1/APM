from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


@dataclass(frozen=True)
class ExternalCovariateSpec:
    name: str
    path: Path
    sep: str = "\t"
    id_col: str = "sample_id"  # raw tcga id column in the external table


class ExternalTableProvider(CovariateProvider):
    """
    Load a user-provided covariate table and normalize TCGA identifiers.

    Intended uses:
    - HRD scores (CHORD/HRDetect/HRD-LOH/etc) exported from elsewhere
    - curated mutation status tables
    """

    def __init__(self, spec: ExternalCovariateSpec) -> None:
        self.spec = spec
        self.name = f"external:{spec.name}"

    def build(self) -> ProviderResult:
        p = Path(self.spec.path)
        if not p.exists():
            df = pd.DataFrame(index=pd.Index([], name=self.spec.id_col))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="raw",
                diagnostics={"error": f"missing {p}", **_diagnose_index(df)},
            )

        df = pd.read_csv(p, sep=self.spec.sep, low_memory=False)
        if self.spec.id_col in df.columns:
            add_tcga_id_columns_inplace(df, raw_id_col=self.spec.id_col, overwrite=False)
            df = df.set_index(self.spec.id_col, drop=True)
        df.index = df.index.astype(str)
        df.index.name = df.index.name or "raw_id"

        return ProviderResult(
            name=self.name,
            df=df,
            index_key="raw",
            diagnostics={"path": str(p), "id_col": self.spec.id_col, **_diagnose_index(df)},
        )

