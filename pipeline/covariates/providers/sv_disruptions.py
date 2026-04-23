from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from ._base import CovariateProvider, ProviderResult, _diagnose_index
from pipeline.SV.gene_mirroring import SvGeneSummaryConfig, summarize_sv_gene_hits


class SvDisruptionProvider(CovariateProvider):
    """
    Summarize SV `gene_hits` into per-sample, per-gene disruption covariates.

    Uses processed SV CSVs that already include `gene_hits` with region flags and `hit_side`.
    This is BND-aware via `hit_side` (bp1/bp2/point vs span).

    Important limitation:
    - There is no explicit `CDS_flag` in current SV `gene_hits`. We approximate coding disruption
      using `exon_flag` (+ start/stop flags) and (when present) `transcript_type == protein_coding`.
    """

    name = "sv_disruptions"

    def __init__(
        self,
        *,
        sv_output_root: Path,
        stage_dir: str = "07_final_sv_with_fimo",
        genes_of_interest: Optional[Iterable[str]] = None,
    ) -> None:
        self.sv_output_root = Path(sv_output_root)
        self.stage_dir = str(stage_dir)
        self.genes_of_interest = set(g.upper() for g in (genes_of_interest or [])) or None

    def build(self) -> ProviderResult:
        d = self.sv_output_root / self.stage_dir
        if not d.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={"error": f"missing {d}", **_diagnose_index(df)},
            )

        out = summarize_sv_gene_hits(
            config=SvGeneSummaryConfig(sv_output_root=self.sv_output_root, stage_dir=self.stage_dir),
            genes_of_interest=sorted(self.genes_of_interest) if self.genes_of_interest else None,
        )

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample_vial",
            diagnostics={
                "sv_dir": str(d),
                "genes_of_interest": sorted(self.genes_of_interest) if self.genes_of_interest else None,
                **_diagnose_index(out),
            },
        )

