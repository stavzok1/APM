from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import pandas as pd


@dataclass(frozen=True)
class ProviderResult:
    name: str
    df: pd.DataFrame
    index_key: str  # one of: sample_vial, sample, participant, aliquot, raw
    diagnostics: Dict[str, object]


class CovariateProvider:
    """
    A covariate provider yields a tidy table keyed by some TCGA id granularity.
    Providers should include normalized id columns whenever possible.
    """

    name: str = "provider"

    def build(self) -> ProviderResult:
        raise NotImplementedError


def _diagnose_index(df: pd.DataFrame) -> Dict[str, object]:
    return {
        "n_rows": int(df.shape[0]),
        "n_cols": int(df.shape[1]),
        "index_name": df.index.name,
        "n_unique_index": int(df.index.nunique()) if df.index is not None else None,
        "has_duplicates": bool(df.index.has_duplicates) if df.index is not None else None,
    }

