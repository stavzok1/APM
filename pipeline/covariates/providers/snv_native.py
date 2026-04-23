from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from pipeline.sample_ids import normalize_tcga_id
from pipeline.config import PATHS
from pipeline.SNV.vcf_loader import _extract_primary_sample_id_from_samples_tsv

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class SnvNativeProvider(CovariateProvider):
    """
    Pipeline-native SNV covariates from `PATHS.snv_output_dir/per_sample/*_snv_summary.json`.

    Note: your current SNV pipeline is APM-region-focused; this provider is best-effort and
    may not capture genome-wide drivers (TP53/PIK3CA) unless present in your SNV region.
    """

    name = "snv_native"

    def __init__(self, *, snv_output_dir: Path) -> None:
        self.snv_output_dir = Path(snv_output_dir)

    def build(self) -> ProviderResult:
        per_sample = self.snv_output_dir / "per_sample"
        if not per_sample.exists():
            df = pd.DataFrame(index=pd.Index([], name="sample_vial"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="sample_vial",
                diagnostics={"error": f"missing {per_sample}", **_diagnose_index(df)},
            )

        rows: List[Dict[str, object]] = []
        for p in sorted(per_sample.glob("*_snv_summary.json")):
            try:
                d = json.loads(p.read_text())
            except Exception:
                continue

            # Map this per-sample artifact to a tumor sample id:
            # read 1-row preview from the paired variants CSV to get source_file (VCF),
            # then look up the tumor sample in annotations/SNV/samples.tsv.
            variants_csv = Path(str(p).replace("_snv_summary.json", "_snv_variants.csv"))
            sample_vial: Optional[str] = None
            if variants_csv.exists():
                try:
                    preview = pd.read_csv(variants_csv, nrows=1)
                    if "source_file" in preview.columns and len(preview) > 0:
                        vcf_name = str(preview["source_file"].iloc[0]).strip()
                        if vcf_name:
                            tumor_id = _extract_primary_sample_id_from_samples_tsv(
                                PATHS.annotations_dir / "SNV" / "samples.tsv",
                                Path(vcf_name).name,
                            )
                            ids = normalize_tcga_id(str(tumor_id))
                            sample_vial = ids.sample_vial or ids.sample or ids.participant or ids.raw
                except Exception:
                    sample_vial = None

            if not sample_vial:
                # Fallback: if the filename includes a full TCGA barcode token, use it; else skip.
                stem = p.name
                tok = next((t for t in stem.split(".") if t.startswith("TCGA-") and len(t.split("-")) >= 4), "")
                if tok:
                    ids = normalize_tcga_id(tok)
                    sample_vial = ids.sample_vial or ids.sample or ids.participant or ids.raw

            if not sample_vial:
                continue

            row = {"sample_vial": str(sample_vial)}
            # Keep a few generic SNV burden/quality covariates.
            for k in (
                "n_total_variants",
                "n_unique_positions",
                "n_variants_in_ccres",
                "tumor_vaf_median",
                "tumor_vaf_mean",
            ):
                if k in d:
                    row[f"snv_{k}"] = d[k]
            rows.append(row)

        out = pd.DataFrame(rows).set_index("sample_vial") if rows else pd.DataFrame(index=pd.Index([], name="sample_vial"))
        out.index = out.index.astype(str)
        out.index.name = "sample_vial"

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="sample_vial",
            diagnostics={"snv_output_dir": str(self.snv_output_dir), **_diagnose_index(out)},
        )

