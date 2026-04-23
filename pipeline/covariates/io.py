from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

import json
import pandas as pd


@dataclass(frozen=True)
class CovariatesRunPaths:
    out_dir: Path
    table_parquet: Path
    table_csv: Path
    metadata_json: Path


def default_run_id(prefix: str = "run") -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"{prefix}_{ts}"


def resolve_output_paths(*, repo_data_dir: Path, run_id: str) -> CovariatesRunPaths:
    out_dir = Path(repo_data_dir) / "covariates" / run_id
    return CovariatesRunPaths(
        out_dir=out_dir,
        table_parquet=out_dir / "cohort_covariates.parquet",
        table_csv=out_dir / "cohort_covariates.csv",
        metadata_json=out_dir / "metadata.json",
    )


def write_covariates_outputs(
    df: pd.DataFrame,
    *,
    paths: CovariatesRunPaths,
    metadata: Optional[Dict[str, Any]] = None,
    write_csv: bool = True,
) -> None:
    paths.out_dir.mkdir(parents=True, exist_ok=True)

    # Parquet is the canonical output (preserves types + is fast).
    df.to_parquet(paths.table_parquet, index=True)

    if write_csv:
        df.to_csv(paths.table_csv, index=True)

    meta = dict(metadata or {})
    meta.setdefault("n_rows", int(df.shape[0]))
    meta.setdefault("n_cols", int(df.shape[1]))
    meta.setdefault("index_name", df.index.name)
    meta.setdefault("written_utc", datetime.now(timezone.utc).isoformat())

    paths.metadata_json.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")

