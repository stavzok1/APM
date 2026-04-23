#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Master SV → motif pipeline (resumable + logs + preflight + JSON config)
#
# Usage (JSON mode):
#   ./run_sv_motif_pipeline.sh --config sv_paths.json [--threads N] [--force] [--from STEP] [--to STEP]
#
# Usage (positional mode, backwards compatible):
#   ./run_sv_motif_pipeline.sh [--threads N] [--force] [--from STEP] [--to STEP] \
#       SCRIPTS_DIR RAW_VCF_DIR OUTPUT_ROOT GENE_FEATURES SAMPLES_TSV \
#       REG_ELEM_FILE LNCRNA_FEATURES REF_FASTA MEME_FILE FLANK_SIZE
#
# Steps:
#   vep | process | bed | fasta | fimo | intersect | recombine

FORCE=0
FROM_STEP="vep"
TO_STEP="recombine"
THREADS=4
CONFIG_JSON=""

usage() {
  echo "JSON mode:" >&2
  echo "  $0 --config sv_paths.json [--threads N] [--force] [--from STEP] [--to STEP]" >&2
  echo "" >&2
  echo "Positional mode:" >&2
  echo "  $0 [--threads N] [--force] [--from STEP] [--to STEP] SCRIPTS_DIR RAW_VCF_DIR OUTPUT_ROOT GENE_FEATURES SAMPLES_TSV REG_ELEM_FILE LNCRNA_FEATURES REF_FASTA MEME_FILE FLANK_SIZE" >&2
  echo "" >&2
  echo "STEP names: vep | process | bed | fasta | fimo | intersect | recombine" >&2
}

# -------------------------
# Parse optional flags
# -------------------------
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=1; shift ;;
    --from) FROM_STEP="$2"; shift 2 ;;
    --to) TO_STEP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --config) CONFIG_JSON="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

if [[ $# -gt 0 ]]; then
  POSITIONAL+=("$@")
fi

# -------------------------
# Helpers
# -------------------------
need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing command in PATH: $1" >&2; exit 1; }
}
need_file() {
  [[ -f "$1" ]] || { echo "Missing file: $1" >&2; exit 1; }
}
need_dir() {
  [[ -d "$1" ]] || { echo "Missing directory: $1" >&2; exit 1; }
}

json_get() {
  local key="$1"
  python - <<PY
import json, sys
p = json.load(open("${CONFIG_JSON}"))
k = "${key}"
if k not in p:
    sys.stderr.write(f"Missing key in config JSON: {k}\\n")
    sys.exit(2)
print(p[k])
PY
}

step_index() {
  case "$1" in
    vep) echo 1 ;;
    process) echo 2 ;;
    bed) echo 3 ;;
    fasta) echo 4 ;;
    fimo) echo 5 ;;
    intersect) echo 6 ;;
    recombine) echo 7 ;;
    *) echo 999 ;;
  esac
}
should_run_step() {
  local s="$1"
  local si fi ti
  si="$(step_index "$s")"
  fi="$(step_index "$FROM_STEP")"
  ti="$(step_index "$TO_STEP")"
  [[ "$si" -ge "$fi" && "$si" -le "$ti" ]]
}

run_step() {
  local name="$1"; shift
  local marker="${OUTPUT_ROOT}/.${name}.done"
  local step_log="${LOG_DIR}/${name}.log"

  if [[ -f "$marker" && "$FORCE" -eq 0 ]]; then
    echo "[SKIP] $name (marker exists: $marker)"
    return 0
  fi

  echo "============================================================"
  echo "[RUN ] STEP: $name"
  echo "[LOG ] $step_log"
  echo "============================================================"

  ( "$@" ) 2>&1 | tee "$step_log"
  touch "$marker"

  echo "[DONE] STEP: $name"
  echo
}

on_error() {
  local exit_code=$?
  echo
  echo "[ERROR] Pipeline failed (exit code $exit_code). Last command: ${BASH_COMMAND}" >&2
  echo "[ERROR] See logs in: ${LOG_DIR}" >&2
  exit "$exit_code"
}
trap on_error ERR

# -------------------------
# Resolve inputs: JSON mode or positional mode
# -------------------------
if [[ -n "$CONFIG_JSON" ]]; then
  need_file "$CONFIG_JSON"
  need_cmd python

  SCRIPTS_DIR="$(json_get SCRIPTS_DIR)"
  RAW_VCF_DIR="$(json_get RAW_VCF_DIR)"
  OUTPUT_ROOT="$(json_get OUTPUT_ROOT)"
  GENE_FEATURES="$(json_get GENE_FEATURES)"
  SAMPLES_TSV="$(json_get SAMPLES_TSV)"
  REG_ELEM_FILE="$(json_get REG_ELEM_FILE)"
  LNCRNA_FEATURES="$(json_get LNCRNA_FEATURES)"
  REF_FASTA="$(json_get REF_FASTA)"
  MEME_FILE="$(json_get MEME_FILE)"
  FLANK_SIZE="$(json_get FLANK_SIZE)"

  # Optional overrides from JSON if present (otherwise keep CLI defaults)
  THREADS_JSON="$(python - <<PY
import json
p=json.load(open("${CONFIG_JSON}"))
print(p.get("THREADS",""))
PY
)"
  if [[ -n "$THREADS_JSON" ]]; then
    THREADS="$THREADS_JSON"
  fi

else
  if [[ ${#POSITIONAL[@]} -ne 10 ]]; then
    usage
    echo "Got ${#POSITIONAL[@]} positional args; expected 10 (or use --config)." >&2
    exit 1
  fi

  SCRIPTS_DIR="${POSITIONAL[0]}"
  RAW_VCF_DIR="${POSITIONAL[1]}"
  OUTPUT_ROOT="${POSITIONAL[2]}"
  GENE_FEATURES="${POSITIONAL[3]}"
  SAMPLES_TSV="${POSITIONAL[4]}"
  REG_ELEM_FILE="${POSITIONAL[5]}"
  LNCRNA_FEATURES="${POSITIONAL[6]}"
  REF_FASTA="${POSITIONAL[7]}"
  MEME_FILE="${POSITIONAL[8]}"
  FLANK_SIZE="${POSITIONAL[9]}"
fi

# -------------------------
# Normalize / validate
# -------------------------
SCRIPTS_DIR="$(realpath "$SCRIPTS_DIR")"
RAW_VCF_DIR="$(realpath "$RAW_VCF_DIR")"
OUTPUT_ROOT="$(realpath -m "$OUTPUT_ROOT")"
GENE_FEATURES="$(realpath "$GENE_FEATURES")"
SAMPLES_TSV="$(realpath "$SAMPLES_TSV")"
REG_ELEM_FILE="$(realpath "$REG_ELEM_FILE")"
LNCRNA_FEATURES="$(realpath "$LNCRNA_FEATURES")"
REF_FASTA="$(realpath "$REF_FASTA")"
MEME_FILE="$(realpath "$MEME_FILE")"

if ! [[ "$FLANK_SIZE" =~ ^[0-9]+$ ]]; then
  echo "FLANK_SIZE must be a non-negative integer. Got: $FLANK_SIZE" >&2
  exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "THREADS must be a positive integer. Got: $THREADS" >&2
  exit 1
fi

need_dir "$SCRIPTS_DIR"
need_dir "$RAW_VCF_DIR"
need_file "$GENE_FEATURES"
need_file "$SAMPLES_TSV"
need_file "$REG_ELEM_FILE"
need_file "$LNCRNA_FEATURES"
need_file "$REF_FASTA"
need_file "$MEME_FILE"

need_cmd bash
need_cmd python
need_cmd realpath
need_cmd tee

# Tools likely used downstream
need_cmd vep || true
need_cmd bedtools || true
need_cmd samtools || true
need_cmd fimo || true

# -------------------------
# Directory structure
# -------------------------
VEP_VCF_DIR="${OUTPUT_ROOT}/01_vep_vcfs"
PROC_CSV_DIR="${OUTPUT_ROOT}/02_processed_sv_csv"
SV_BED_DIR="${OUTPUT_ROOT}/03_sv_bed"
SV_FASTA_DIR="${OUTPUT_ROOT}/04_sv_fasta"
FIMO_TSV_DIR="${OUTPUT_ROOT}/05_fimo_tsv"
FIMO_MERGED_DIR="${OUTPUT_ROOT}/06_sv_fimo_merged"
FINAL_CSV_DIR="${OUTPUT_ROOT}/07_final_sv_with_fimo"
LOG_DIR="${OUTPUT_ROOT}/logs"

mkdir -p \
  "$VEP_VCF_DIR" "$PROC_CSV_DIR" "$SV_BED_DIR" "$SV_FASTA_DIR" \
  "$FIMO_TSV_DIR" "$FIMO_MERGED_DIR" "$FINAL_CSV_DIR" "$LOG_DIR"

PIPELINE_LOG="${LOG_DIR}/pipeline.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$PIPELINE_LOG") 2>&1

echo "[INFO] Using scripts from: $SCRIPTS_DIR"
echo "[INFO] RAW_VCF_DIR:        $RAW_VCF_DIR"
echo "[INFO] OUTPUT_ROOT:        $OUTPUT_ROOT"
echo "[INFO] FROM_STEP:          $FROM_STEP"
echo "[INFO] TO_STEP:            $TO_STEP"
echo "[INFO] FORCE:              $FORCE"
echo "[INFO] THREADS:            $THREADS"
echo "[INFO] FLANK_SIZE:         $FLANK_SIZE"
if [[ -n "$CONFIG_JSON" ]]; then
  echo "[INFO] CONFIG_JSON:        $(realpath "$CONFIG_JSON")"
fi
echo "[INFO] Pipeline log:       $PIPELINE_LOG"
echo

# Propagate force into sub-scripts when supported
VEP_FORCE_ARGS=()
FIMO_FORCE_ARGS=()
if [[ "$FORCE" -eq 1 ]]; then
  VEP_FORCE_ARGS+=(--force)
  FIMO_FORCE_ARGS+=(--force)
fi

################################################################################
# STEP 1 — VEP annotation
################################################################################
if should_run_step vep; then
  # Assumes you updated run_vep_batch.sh to accept --fork/--force (recommended)
  run_step vep \
    bash "${SCRIPTS_DIR}/run_vep_batch.sh" \
      "${VEP_FORCE_ARGS[@]}" \
      --fork "$THREADS" \
      "$RAW_VCF_DIR" \
      "$VEP_VCF_DIR"
fi

################################################################################
# STEP 2 — Process VEP VCFs → strict_sv_set.csv
################################################################################
if should_run_step process; then
  run_step process \
    python "${SCRIPTS_DIR}/SV_vcf_handling.py" \
      "$VEP_VCF_DIR" \
      "$GENE_FEATURES" \
      "$SAMPLES_TSV" \
      "$PROC_CSV_DIR" \
      "$REG_ELEM_FILE" \
      "$LNCRNA_FEATURES"
fi

################################################################################
# STEP 3 — Produce SV BED intervals
################################################################################
if should_run_step bed; then
  run_step bed \
    python "${SCRIPTS_DIR}/SV_create_bed.py" \
      "$PROC_CSV_DIR" \
      "$SV_BED_DIR" \
      "$FLANK_SIZE"
fi

################################################################################
# STEP 4 — FASTA extraction for FIMO
################################################################################
if should_run_step fasta; then
  run_step fasta \
    bash "${SCRIPTS_DIR}/create_fasta_for_sv_bed.sh" \
      "$REF_FASTA" \
      "$SV_BED_DIR" \
      "$SV_FASTA_DIR"
fi

################################################################################
# STEP 5 — Run FIMO on FASTA
################################################################################
if should_run_step fimo; then
  # Assumes you updated run_fimo_on_fasta.sh to accept --jobs/--force (recommended)
  run_step fimo \
    bash "${SCRIPTS_DIR}/run_fimo_on_fasta.sh" \
      "${FIMO_FORCE_ARGS[@]}" \
      --jobs "$THREADS" \
      "$MEME_FILE" \
      "$SV_FASTA_DIR" \
      "$FIMO_TSV_DIR"
fi

################################################################################
# STEP 6 — Intersect FIMO hits with SV BEDs
################################################################################
if should_run_step intersect; then
  run_step intersect \
    bash "${SCRIPTS_DIR}/run_fimo_sv_intersect.sh" \
      "$FIMO_TSV_DIR" \
      "$SV_BED_DIR" \
      "$FIMO_MERGED_DIR"
fi

################################################################################
# STEP 7 — Recombine motif hits into processed SV CSVs
################################################################################
if should_run_step recombine; then
  run_step recombine \
    python "${SCRIPTS_DIR}/recombine_sv_fimo.py" \
      "$FIMO_MERGED_DIR" \
      "$PROC_CSV_DIR" \
      "$FINAL_CSV_DIR"
fi

echo "[INFO] Pipeline complete."
echo "[INFO] Final annotated files: $FINAL_CSV_DIR"
echo "[INFO] Logs directory: $LOG_DIR"
echo
