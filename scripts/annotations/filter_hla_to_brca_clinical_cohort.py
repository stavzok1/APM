#!/usr/bin/env python3
"""
Keep HLA rows whose aliquot maps to a TCGA participant that appears in the BRCA clinical table
(with at least one sample row in that table).

Writes:
  annotations/HLA/HLA_types.unfiltered.comma.csv  — one-time copy of the full comma table (source of truth for reruns)
  annotations/HLA_types.tsv                       — filtered comma-separated table (same columns as before)

If the unfiltered archive already exists, filtering is always applied from that archive so reruns stay stable.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import normalize_tcga_id  # noqa: E402


def _clinical_participants_with_samples(clinical_path: Path) -> set[str]:
    if not clinical_path.exists():
        raise FileNotFoundError(f"BRCA clinical file not found: {clinical_path}")
    df = pd.read_csv(clinical_path, sep="\t", low_memory=False)
    col = None
    for c in ("sampleID", "sample_id", "bcr_sample_barcode"):
        if c in df.columns:
            col = c
            break
    if col is None:
        raise KeyError(f"No sample id column in {clinical_path}; expected sampleID / sample_id / bcr_sample_barcode")
    out: set[str] = set()
    for raw in df[col].dropna().astype(str):
        raw = raw.strip()
        if not raw or raw.lower() == "nan":
            continue
        tid = normalize_tcga_id(raw)
        if tid.participant:
            out.add(tid.participant)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Filter HLA_types to BRCA clinical cohort participants.")
    ap.add_argument(
        "--clinical",
        type=Path,
        default=None,
        help="BRCA clinical wide table (default: PATHS.brca_clinical)",
    )
    ap.add_argument(
        "--hla-out",
        type=Path,
        default=None,
        help="Filtered HLA output (default: PATHS.hla_types_tsv)",
    )
    ap.add_argument(
        "--unfiltered-archive",
        type=Path,
        default=None,
        help="Archive path for full HLA table (default: annotations/HLA/HLA_types.unfiltered.comma.csv)",
    )
    args = ap.parse_args()

    clinical_path = args.clinical or PATHS.brca_clinical
    hla_out = args.hla_out or PATHS.hla_types_tsv
    archive = args.unfiltered_archive or (PATHS.annotations_dir / "HLA" / "HLA_types.unfiltered.comma.csv")

    allowed = _clinical_participants_with_samples(clinical_path)
    print(f"BRCA clinical participants with ≥1 sample row: {len(allowed)}")

    archive.parent.mkdir(parents=True, exist_ok=True)

    if not archive.exists():
        if not hla_out.exists():
            raise FileNotFoundError(f"Neither archive nor HLA output exists: {archive} / {hla_out}")
        shutil.copy2(hla_out, archive)
        print(f"Created unfiltered archive: {archive}")
    else:
        print(f"Using unfiltered archive: {archive}")

    df = pd.read_csv(archive, sep=",", low_memory=False)
    if "aliquot_id" not in df.columns:
        raise KeyError("HLA table must contain aliquot_id")

    def participant_for_row(aid: object) -> str:
        tid = normalize_tcga_id(str(aid).strip())
        return tid.participant or ""

    parts = df["aliquot_id"].map(participant_for_row)
    mask = parts.isin(allowed)
    filtered = df.loc[mask].copy()
    print(f"HLA rows: {len(df)} -> {len(filtered)} (dropped {len(df) - len(filtered)})")

    hla_out.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(hla_out, sep=",", index=False)
    print(f"Wrote filtered HLA: {hla_out}")


if __name__ == "__main__":
    main()
