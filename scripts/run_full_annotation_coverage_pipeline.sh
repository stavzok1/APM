#!/usr/bin/env bash
# HLA clinical cohort filter → unified stage column → normalize → manifests →
# sample_module_coverage → clinical omics stratification (writes under latest run_*).
#
# Note: this script must run under bash. If invoked via `sh`, `pipefail` may not exist.
set -eu
set -o pipefail 2>/dev/null || true

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="${ROOT}/.venv/bin/python3"
if [[ ! -x "$PY" ]]; then
  echo "Missing venv python: $PY" >&2
  exit 1
fi

"$PY" scripts/annotations/filter_hla_to_brca_clinical_cohort.py
"$PY" scripts/annotations/update_brca_clinical_unified_stage.py
"$PY" scripts/annotations/normalize_external_annotations.py
"$PY" scripts/annotations/build_annotation_sample_manifests.py
"$PY" analysis/sample_module_coverage.py

RUN="$(ls -td "${ROOT}/analysis/sample_coverage/output/run_"* 2>/dev/null | head -1 || true)"
if [[ -z "$RUN" ]]; then
  echo "No analysis/sample_coverage/output/run_* directory found." >&2
  exit 1
fi

PRES="${RUN}/tables/omics/participant_presence_omics.tsv"
if [[ ! -f "$PRES" ]]; then
  echo "Missing $PRES" >&2
  exit 1
fi

"$PY" analysis/clinical_omics_stratification.py --presence "$PRES" --out "${RUN}/stratification"

echo "RUN=$RUN"
ls "$RUN" | head -40
echo "--- stratification ---"
ls "${RUN}/stratification"
