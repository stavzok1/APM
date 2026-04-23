from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_from_index_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class RppaProvider(CovariateProvider):
    name = "rppa"

    def __init__(self, *, combined_path: Path) -> None:
        self.combined_path = Path(combined_path)

    def build(self) -> ProviderResult:
        if not self.combined_path.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={"error": f"missing {self.combined_path}", **_diagnose_index(df)},
            )

        if self.combined_path.suffix.lower() == ".parquet":
            df = pd.read_parquet(self.combined_path)
        else:
            df = pd.read_csv(self.combined_path)

        # rppa_main writes index as sample_id; keep as index and add normalized columns.
        if df.index.name is None or df.index.name == "":
            # If CSV was read without index, try to recover
            if "sample_id" in df.columns:
                df = df.set_index("sample_id", drop=True)
        df.index = df.index.astype(str)
        df.index.name = df.index.name or "sample_id"

        add_tcga_id_columns_from_index_inplace(df, overwrite=False)

        # Choose index granularity: prefer sample_vial if derivable.
        if "sample_vial" in df.columns:
            out = df.copy()
            out["sample_vial"] = out["sample_vial"].astype("string")
            out.index = out["sample_vial"].where(out["sample_vial"].notna(), out.index.astype(str)).astype(str)
            out.index.name = "sample_vial"
        else:
            out = df.copy()
            out.index.name = "sample_id"

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample_vial" if out.index.name == "sample_vial" else "raw",
            diagnostics={"combined_path": str(self.combined_path), **_diagnose_index(out)},
        )

