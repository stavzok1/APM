#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# Defaults / inputs
# -------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATHS_JSON="${HERE}/atac_paths.json"
MANIFEST="${HERE}/normals_manifest.tsv"

PEAKSET=""               # required
FIRST_SAMPLE_COL="7"
PRIMARY_DISEASE="breast invasive carcinoma"
TCGA_CANCER_TYPE="BRCA"

FORCE=0
FROM_STAGE="recount"     # recount|harmonize|normalize|collapse
TO_STAGE="collapse"

usage() {
  cat <<EOF
Usage:
  $0 --peakset <TCGA_ATAC_PeakSet.bed> [options]

Options:
  --paths-json <file>          (default: atac_paths.json)
  --manifest <file>            (default: normals_manifest.tsv)
  --first-sample-col <int>     (default: 7)
  --primary-disease <string>   (default: breast invasive carcinoma)
  --tcga-cancer-type <string>  (default: BRCA)
  --from <stage>               (default: recount)
  --to <stage>                 (default: collapse)
  --force                      re-run steps even if outputs exist

Stages: recount, harmonize, normalize, collapse

Example:
  $0 --peakset data/TCGA_ATAC/TCGA_ATAC_PeakSet.180502.bed
  $0 --peakset data/TCGA_ATAC/TCGA_ATAC_PeakSet.180502.bed --from harmonize
EOF
}

# -------------------------
# Parse args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --peakset) PEAKSET="$2"; shift 2 ;;
    --paths-json) PATHS_JSON="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --first-sample-col) FIRST_SAMPLE_COL="$2"; shift 2 ;;
    --primary-disease) PRIMARY_DISEASE="$2"; shift 2 ;;
    --tcga-cancer-type) TCGA_CANCER_TYPE="$2"; shift 2 ;;
    --from) FROM_STAGE="$2"; shift 2 ;;
    --to) TO_STAGE="$2"; shift 2 ;;
    --force) FORCE=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$PEAKSET" ]]; then
  echo "Error: --peakset is required"
  usage
  exit 1
fi

[[ -f "$PATHS_JSON" ]] || { echo "Error: PATHS JSON not found: $PATHS_JSON"; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "Error: manifest not found: $MANIFEST"; exit 1; }
[[ -f "$PEAKSET" ]] || { echo "Error: peakset not found: $PEAKSET"; exit 1; }

# -------------------------
# Helpers
# -------------------------
stage_index() {
  case "$1" in
    recount) echo 1 ;;
    harmonize) echo 2 ;;
    normalize) echo 3 ;;
    collapse) echo 4 ;;
    *) echo 999 ;;
  esac
}
should_run_stage() {
  local s="$1"
  local si fi ti
  si="$(stage_index "$s")"
  fi="$(stage_index "$FROM_STAGE")"
  ti="$(stage_index "$TO_STAGE")"
  [[ "$si" -ge "$fi" && "$si" -le "$ti" ]]
}

# Read value from JSON
get_json_value() {
  python - <<PY
import json
p = json.load(open("$PATHS_JSON"))
key = "$1"
if key not in p:
    raise KeyError(f"Missing key in JSON: {key}")
print(p[key])
PY
}

# Pull paths from JSON
OUT_RAW="$(get_json_value OUT_PATH)"
OUT_LOGCPM_QN="$(get_json_value ATAC_LOGCPM_QN_PATH)"
OUT_CASE_LEVEL="$(get_json_value ATAC_CASE_LEVEL_PATH)"
HEALTHY_COUNTS_DIR="$(get_json_value HEALTHY_COUNTS_DIR)"

mkdir -p "$HEALTHY_COUNTS_DIR"

# -------------------------
# Stage 1: recount normals
# -------------------------
if should_run_stage recount; then
  echo "=== STAGE: recount normals ==="

  # Expect header: sample_id<TAB>bam_path
  tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample_id bam_path; do
    [[ -n "${sample_id:-}" ]] || continue
    [[ -n "${bam_path:-}" ]] || { echo "Error: missing bam_path for sample_id=$sample_id"; exit 1; }
    [[ -f "$bam_path" ]] || { echo "Error: BAM not found: $bam_path"; exit 1; }

    out_vec="${HEALTHY_COUNTS_DIR}/${sample_id}.vector"
    if [[ -f "$out_vec" && "$FORCE" -eq 0 ]]; then
      echo "[skip] $sample_id vector exists: $out_vec"
      continue
    fi

    echo "[run] recount: $sample_id"
    "${HERE}/recount_atac.sh" "$PEAKSET" "$bam_path" "$HEALTHY_COUNTS_DIR"

    # recount outputs PREFIX.vector where PREFIX = BAM basename without .bam
    bam_base="$(basename "$bam_path" .bam)"
    produced_vec="${HEALTHY_COUNTS_DIR}/${bam_base}.vector"

    if [[ ! -f "$produced_vec" ]]; then
      echo "Error: recount did not produce expected vector: $produced_vec"
      exit 1
    fi

    # Rename to sample_id if needed
    if [[ "$bam_base" != "$sample_id" ]]; then
      mv "$produced_vec" "$out_vec"
    else
      # Ensure naming consistent
      mv "$produced_vec" "$out_vec" 2>/dev/null || true
    fi

    [[ -f "$out_vec" ]] || { echo "Error: final vector missing for $sample_id: $out_vec"; exit 1; }
    echo "[ok] wrote: $out_vec"
  done
fi

# -------------------------
# Stage 2: harmonize TCGA + normals
# -------------------------
if should_run_stage harmonize; then
  echo "=== STAGE: harmonize TCGA + normals ==="
  if [[ -f "$OUT_RAW" && "$FORCE" -eq 0 ]]; then
    echo "[skip] already exists: $OUT_RAW"
  else
    python "${HERE}/harmonize_atac_counts.py" \
      --paths-json "$PATHS_JSON" \
      --first-sample-col "$FIRST_SAMPLE_COL" \
      --primary-disease "$PRIMARY_DISEASE" \
      --tcga-cancer-type "$TCGA_CANCER_TYPE"
  fi
fi

# -------------------------
# Stage 3: normalize (edgeR logCPM + QN)
# -------------------------
if should_run_stage normalize; then
  echo "=== STAGE: normalize logCPM + QN ==="
  if [[ -f "$OUT_LOGCPM_QN" && "$FORCE" -eq 0 ]]; then
    echo "[skip] already exists: $OUT_LOGCPM_QN"
  else
    Rscript "${HERE}/normalize_atac_edger.R" \
      "$OUT_RAW" \
      "$OUT_LOGCPM_QN" \
      "$FIRST_SAMPLE_COL"
  fi
fi

# -------------------------
# Stage 4: collapse replicates to Case_ID
# -------------------------
if should_run_stage collapse; then
  echo "=== STAGE: collapse replicates ==="
  if [[ -f "$OUT_CASE_LEVEL" && "$FORCE" -eq 0 ]]; then
    echo "[skip] already exists: $OUT_CASE_LEVEL"
  else
    # collapse script reads ATAC_LOGCPM_QN_PATH + ATAC_META_PATH from JSON,
    # and we override output to ATAC_CASE_LEVEL_PATH explicitly to avoid ambiguity.
    python "${HERE}/collapse_atac_replicates.py" \
      --paths-json "$PATHS_JSON" \
      --first-sample-col "$FIRST_SAMPLE_COL" \
      --out-file "$OUT_CASE_LEVEL"
  fi
fi

echo "Pipeline complete."
echo "Raw combined:     $OUT_RAW"
echo "LogCPM + QN:      $OUT_LOGCPM_QN"
echo "Case-level final: $OUT_CASE_LEVEL"
