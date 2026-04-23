from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from pipeline.sample_ids import add_tcga_id_columns_from_index_inplace

from ._base import CovariateProvider, ProviderResult, _diagnose_index


class HiChipParticipantProvider(CovariateProvider):
    """
    Participant-level covariates from the TCGA-wide HiChIP processed matrix.

    The file shape in this repo (PATHS.hichip_tcga_processed_csv):
      - first columns are loop anchors (chr1,s1,e1,chr2,s2,e2)
      - remaining columns are TCGA participants (TCGA-XX-YYYY)

    We treat values as continuous interaction strengths and compute simple summaries per participant.
    """

    name = "hichip_participant"

    def __init__(self, *, path: Path, chunksize: int = 50_000) -> None:
        self.path = Path(path)
        self.chunksize = int(chunksize)

    def build(self) -> ProviderResult:
        if not self.path.exists():
            df = pd.DataFrame(index=pd.Index([], name="participant"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="participant",
                diagnostics={"error": f"missing {self.path}", **_diagnose_index(df)},
            )

        # Read header to get participant columns.
        hdr = pd.read_csv(self.path, nrows=0)
        part_cols = [c for c in hdr.columns if str(c).strip().startswith("TCGA-")]
        if not part_cols:
            df = pd.DataFrame(index=pd.Index([], name="participant"))
            return ProviderResult(
                name=self.name,
                df=df,
                index_key="participant",
                diagnostics={"error": "no TCGA participant columns found", **_diagnose_index(df)},
            )

        sums = pd.Series(0.0, index=part_cols, dtype=float)
        nnz = pd.Series(0, index=part_cols, dtype="int64")
        counts = pd.Series(0, index=part_cols, dtype="int64")

        usecols = part_cols  # we can skip the anchor columns entirely
        for chunk in pd.read_csv(self.path, usecols=usecols, chunksize=self.chunksize, low_memory=False):
            # Ensure numeric
            for c in part_cols:
                chunk[c] = pd.to_numeric(chunk[c], errors="coerce")
            vals = chunk[part_cols]
            sums += vals.sum(axis=0, skipna=True)
            counts += vals.notna().sum(axis=0).astype("int64")
            nnz += (vals.fillna(0.0) > 0).sum(axis=0).astype("int64")

        mean = sums / counts.replace(0, np.nan)
        frac_nz = nnz / counts.replace(0, np.nan)

        out = pd.DataFrame(
            {
                "hichip_sum": sums,
                "hichip_mean": mean,
                "hichip_frac_nonzero": frac_nz,
                "hichip_n_rows": counts,
            }
        )
        out.index = out.index.astype(str)
        out.index.name = "participant"
        add_tcga_id_columns_from_index_inplace(out, overwrite=False)

        return ProviderResult(
            name=self.name,
            df=out,
            index_key="participant",
            diagnostics={
                "path": str(self.path),
                "n_participants": int(len(part_cols)),
                "chunksize": self.chunksize,
                **_diagnose_index(out),
            },
        )

