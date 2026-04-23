from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd

from pipeline.sample_ids import normalize_tcga_id


def utc_run_id(prefix: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"{prefix}_{ts}"


@dataclass(frozen=True)
class CohortRunContext:
    scratch_json: Path
    scratch_dir: Path

    methylation_working_dir: Path
    methylation_subset_manifest: Path

    cnv_annotated_dir: Path
    cnv_gene_tables_dir: Path

    snv_output: Path
    sv_output_root: Path
    rppa_output_dir: Path


def load_context(scratch_json: Path) -> CohortRunContext:
    scratch_json = Path(scratch_json)
    d = json.loads(scratch_json.read_text(encoding="utf-8"))
    scratch_dir = Path(d.get("scratch") or scratch_json.parent)
    return CohortRunContext(
        scratch_json=scratch_json,
        scratch_dir=scratch_dir,
        methylation_working_dir=Path(d["methylation_working_dir"]),
        methylation_subset_manifest=Path(d["methylation_subset_manifest"]),
        cnv_annotated_dir=Path(d["cnv_annotated_dir"]),
        cnv_gene_tables_dir=Path(d["cnv_gene_tables_dir"]),
        snv_output=Path(d["snv_output"]),
        sv_output_root=Path(d["sv_output_root"]),
        rppa_output_dir=Path(d.get("rppa_output_dir") or (Path("data") / "rppa" / "processed")),
    )


def _load_linked_json_from_prepare_log(prepare_log: Path) -> Dict[str, object]:
    txt = Path(prepare_log).read_text(encoding="utf-8", errors="replace")
    start = txt.find("{")
    if start < 0:
        return {}
    blob = txt[start:].strip()
    try:
        return json.loads(blob)
    except Exception:
        start2 = txt.rfind("\n{")
        if start2 >= 0:
            try:
                return json.loads(txt[start2 + 1 :].strip())
            except Exception:
                return {}
        return {}


def load_sample_ids(ctx: CohortRunContext) -> List[str]:
    """
    Prefer methylation manifest sample ids if present; else fall back to prepare-log vials.
    Returns TCGA vial ids (e.g. TCGA-XX-YYYY-01A) when available.
    """
    man = Path(ctx.methylation_subset_manifest)
    if man.exists():
        df = pd.read_csv(man, sep="\t")
        if "Sample ID" in df.columns:
            return [str(x).strip() for x in df["Sample ID"].tolist() if str(x).strip()]

    prep = Path(ctx.scratch_dir) / "cohort_processing_prepare.log"
    if prep.exists():
        linked = _load_linked_json_from_prepare_log(prep)
        out = [str(k).strip() for k in linked.keys() if str(k).strip()]
        out = sorted(set(out))
        return out
    return []


def to_tcga_sample_id(x: str) -> str:
    """Normalize any TCGA-ish id to a TCGA sample id (e.g. TCGA-XX-YYYY-01)."""
    tid = normalize_tcga_id(str(x))
    return str(tid.sample or "").strip()

