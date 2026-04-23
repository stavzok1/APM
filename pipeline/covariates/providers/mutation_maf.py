from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_from_index_inplace

from ..genomics import classify_tp53_from_maf, flag_hotspot_mutations_from_maf
from ._base import CovariateProvider, ProviderResult, _diagnose_index


class MutationMafProvider(CovariateProvider):
    """
    Optional provider: build sample-level mutation covariates from a MAF-like file.

    Designed for cohort post-processing (not part of the main pipeline).
    """

    name = "mutations_maf"

    def __init__(
        self,
        *,
        maf_path: Path,
        sep: str = "\t",
        include_pik3ca_hotspots: bool = True,
    ) -> None:
        self.maf_path = Path(maf_path)
        self.sep = sep
        self.include_pik3ca_hotspots = bool(include_pik3ca_hotspots)

    def build(self) -> ProviderResult:
        if not self.maf_path.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={"error": f"missing {self.maf_path}", **_diagnose_index(df)},
            )

        maf = pd.read_csv(self.maf_path, sep=self.sep, low_memory=False)
        parts = []

        tp53 = classify_tp53_from_maf(maf)
        parts.append(tp53)

        if self.include_pik3ca_hotspots:
            # Common BRCA hotspots.
            pik3ca = flag_hotspot_mutations_from_maf(
                maf,
                gene="PIK3CA",
                aa_hotspots=["p.E542K", "p.E545K", "p.H1047R", "p.H1047L"],
            )
            parts.append(pik3ca)

        out = pd.concat(parts, axis=1) if parts else pd.DataFrame(index=pd.Index([], name="sample_vial"))
        out.index = out.index.astype(str)
        out.index.name = "sample_vial"
        add_tcga_id_columns_from_index_inplace(out, overwrite=False)

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample_vial",
            diagnostics={"maf_path": str(self.maf_path), **_diagnose_index(out)},
        )

