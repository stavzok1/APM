from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_from_index_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


def _latest_run_dir(root: Path, prefix: str = "run_") -> Optional[Path]:
    if not root.exists():
        return None
    runs = sorted(
        [p for p in root.iterdir() if p.is_dir() and p.name.startswith(prefix)],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return runs[0] if runs else None


class CoverageProvider(CovariateProvider):
    name = "coverage"

    def __init__(
        self,
        *,
        analysis_root: Path,
        run_id: Optional[str] = None,
    ) -> None:
        self.analysis_root = Path(analysis_root)
        self.run_id = run_id

    def build(self) -> ProviderResult:
        base = self.analysis_root / "analysis" / "sample_coverage" / "output"
        if self.run_id:
            run_dir = base / self.run_id
        else:
            # Canonical location: output/current (preferred), else latest run_*
            current = base / "current"
            run_dir = current if current.exists() else _latest_run_dir(base)
        if run_dir is None:
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={"error": "no sample_coverage run dir found", **_diagnose_index(df)},
            )

        presence = run_dir / "tables" / "omics" / "sample_vial_presence.tsv"
        if not presence.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={
                    "error": f"missing {presence}",
                    "run_dir": str(run_dir),
                    **_diagnose_index(df),
                },
            )

        df = pd.read_csv(presence, sep="\t", index_col=0)
        df.index = df.index.astype(str)
        df.index.name = "sample_vial"

        # Add normalized columns from index; do not overwrite if present.
        add_tcga_id_columns_from_index_inplace(df, overwrite=False)

        # Summary covariates.
        module_cols = [c for c in df.columns if c not in {"participant", "sample", "sample_vial", "aliquot"}]
        bool_cols = [c for c in module_cols if df[c].dropna().isin([0, 1, True, False]).all()]
        if bool_cols:
            df["n_modules_present"] = df[bool_cols].astype(int).sum(axis=1)

        return ProviderResult(
            name=self.name,
            df=df,
            index_key="sample_vial",
            diagnostics={"run_dir": str(run_dir), "presence_path": str(presence), **_diagnose_index(df)},
        )

