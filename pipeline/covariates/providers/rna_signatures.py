from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence

import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_from_index_inplace
from pipeline.RNA_exp.signatures import compute_gene_set_scores_from_tpm, default_gene_sets

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class RnaSignatureProvider(CovariateProvider):
    name = "rna_signatures"

    def __init__(
        self,
        *,
        tpm_path: Path,
        gene_col: str = "gene_symbol",
        sample_ids: Optional[Sequence[str]] = None,
    ) -> None:
        self.tpm_path = Path(tpm_path)
        self.gene_col = str(gene_col)
        self.sample_ids = list(sample_ids) if sample_ids is not None else None

    def build(self) -> ProviderResult:
        if not self.tpm_path.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample",
                diagnostics={"error": f"missing {self.tpm_path}", **_diagnose_index(df)},
            )

        # If sample_ids not provided, read the header to get TCGA columns.
        if self.sample_ids is None:
            hdr = pd.read_csv(self.tpm_path, sep="\t", nrows=0)
            cols = [c for c in hdr.columns if str(c).strip().startswith("TCGA-")]
            self.sample_ids = cols

        gene_sets = default_gene_sets()
        scores = compute_gene_set_scores_from_tpm(
            self.tpm_path,
            sample_ids=self.sample_ids,
            gene_col=self.gene_col,
            gene_sets=gene_sets,
        )
        scores.index = scores.index.astype(str)
        scores.index.name = "sample_id"

        add_tcga_id_columns_from_index_inplace(scores, overwrite=False)

        # Canonical for RNA is usually sample (…-01), even if index includes vial; normalize if possible.
        if "sample" in scores.columns:
            out = scores.copy()
            out["sample"] = out["sample"].astype("string")
            out.index = out["sample"].where(out["sample"].notna(), out.index.astype(str)).astype(str)
            out.index.name = "sample"
        else:
            out = scores
            out.index.name = "sample_id"

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample",
            diagnostics={"tpm_path": str(self.tpm_path), "n_input_samples": len(self.sample_ids), **_diagnose_index(out)},
        )

