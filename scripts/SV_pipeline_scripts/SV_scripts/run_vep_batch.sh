#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Usage:
#   ./run_vep_batch.sh [--fork N] [--dir-cache PATH] [--fasta PATH] [--force] <input_vcf_dir> <output_dir>
#
# Examples:
#   ./run_vep_batch.sh --fork 8 input_vcfs out_vep
#   ./run_vep_batch.sh --fork 4 --fasta /path/GRCh38.fa input_vcfs out_vep

FORKS=4
DIR_CACHE="${HOME}/.vep"
FASTA="${HOME}/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FORCE=0

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fork) FORKS="$2"; shift 2 ;;
    --dir-cache) DIR_CACHE="$2"; shift 2 ;;
    --fasta) FASTA="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--fork N] [--dir-cache PATH] [--fasta PATH] [--force] <input_vcf_dir> <output_dir>" >&2
      exit 0
      ;;
    *)
      POSITIONAL+=("$1"); shift ;;
  esac
done

if [[ ${#POSITIONAL[@]} -ne 2 ]]; then
  echo "Usage: $0 [--fork N] [--dir-cache PATH] [--fasta PATH] [--force] <input_vcf_dir> <output_dir>" >&2
  exit 1
fi

INPUT_DIR="${POSITIONAL[0]}"
OUTPUT_DIR="${POSITIONAL[1]}"

mkdir -p "$OUTPUT_DIR"

# Basic checks
command -v vep >/dev/null 2>&1 || { echo "Error: vep not in PATH" >&2; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input dir not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "Error: FASTA not found: $FASTA" >&2; exit 1; }
[[ -d "$DIR_CACHE" ]] || { echo "Error: VEP cache dir not found: $DIR_CACHE" >&2; exit 1; }
[[ "$FORKS" =~ ^[0-9]+$ ]] || { echo "Error: --fork must be integer" >&2; exit 1; }

shopt -s nullglob
vcfs=( "$INPUT_DIR"/*.vcf )
shopt -u nullglob

if [[ ${#vcfs[@]} -eq 0 ]]; then
  echo "No VCF files found in: $INPUT_DIR"
  exit 0
fi

for vcf in "${vcfs[@]}"; do
  base="$(basename "$vcf")"
  out_file="$OUTPUT_DIR/$base"

  if [[ -f "$out_file" && "$FORCE" -eq 0 ]]; then
    echo "[skip] VEP output exists: $out_file"
    continue
  fi

  echo "[run ] VEP on: $base (fork=$FORKS)"

  vep \
    --input_file "$vcf" \
    --output_file "$out_file" \
    --offline \
    --cache \
    --dir_cache "$DIR_CACHE" \
    --species homo_sapiens \
    --assembly GRCh38 \
    --fasta "$FASTA" \
    --format vcf \
    --vcf \
    --force_overwrite \
    --everything \
    --symbol \
    --variant_class \
    --check_ref \
    --fork "$FORKS"
done

echo "All VCFs processed. Output written to: $OUTPUT_DIR"
