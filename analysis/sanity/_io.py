from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


def utc_run_id(prefix: str) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"{prefix}_{ts}"


@dataclass(frozen=True)
class ScratchContext:
    scratch_json: Path
    methylation_working_dir: Path
    methylation_subset_manifest: Path
    cnv_annotated_dir: Path
    cnv_gene_tables_dir: Path
    snv_output: Path
    sv_output_root: Path
    rppa_output_dir: Path


def load_scratch_context(scratch_json: Path) -> ScratchContext:
    scratch_json = Path(scratch_json)
    d = json.loads(scratch_json.read_text(encoding="utf-8"))
    return ScratchContext(
        scratch_json=scratch_json,
        methylation_working_dir=Path(d["methylation_working_dir"]),
        methylation_subset_manifest=Path(d["methylation_subset_manifest"]),
        cnv_annotated_dir=Path(d["cnv_annotated_dir"]),
        cnv_gene_tables_dir=Path(d["cnv_gene_tables_dir"]),
        snv_output=Path(d["snv_output"]),
        sv_output_root=Path(d["sv_output_root"]),
        rppa_output_dir=Path(d.get("rppa_output_dir") or (Path("data") / "rppa" / "processed")),
    )


def load_sample_ids_from_methylation_manifest(path: Path) -> List[str]:
    """
    Methylation subset manifest written by cohort runner.
    Columns include: Sample ID, File Name.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    if "Sample ID" not in df.columns:
        raise KeyError(f"Manifest missing 'Sample ID': {path} cols={list(df.columns)}")
    return [str(x).strip() for x in df["Sample ID"].tolist() if str(x).strip()]


def _load_linked_json_from_prepare_log(prepare_log: Path) -> Dict[str, object]:
    """
    `run_cohort_module_processing.py` writes:
      <text lines>\n{ ...linked dict... }\n
    We recover the last JSON object from the file.
    """
    prepare_log = Path(prepare_log)
    txt = prepare_log.read_text(encoding="utf-8", errors="replace")
    start = txt.find("{")
    if start < 0:
        return {}
    # The first "{" is the linked dict in current format.
    blob = txt[start:].strip()
    try:
        return json.loads(blob)
    except Exception:
        # Fallback: take from the last "\n{" occurrence
        start2 = txt.rfind("\n{")
        if start2 >= 0:
            try:
                return json.loads(txt[start2 + 1 :].strip())
            except Exception:
                return {}
        return {}


def load_sample_ids(ctx: ScratchContext) -> List[str]:
    """
    Preferred source is the methylation subset manifest (it enumerates cohort sample_ids).
    If methylation was disabled (or the manifest is missing), fall back to the vial keys
    recorded in `cohort_processing_prepare.log`.
    """
    man = Path(ctx.methylation_subset_manifest)
    if man.exists():
        return load_sample_ids_from_methylation_manifest(man)

    prepare_log = ctx.scratch_json.parent / "cohort_processing_prepare.log"
    if prepare_log.exists():
        linked = _load_linked_json_from_prepare_log(prepare_log)
        out = [str(k).strip() for k in linked.keys() if str(k).strip()]
        out = sorted(set(out))
        if out:
            return out
    return []


def ensure_out_dir(run_id: str) -> Path:
    out = Path("analysis") / "sanity" / "output" / run_id
    out.mkdir(parents=True, exist_ok=True)
    return out

