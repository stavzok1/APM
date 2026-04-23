#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Usage:
#   ./run_fimo_on_fasta.sh [--jobs N] [--thresh X] [--force] <meme_file> <fasta_dir> <out_dir>
#
# Example:
#   ./run_fimo_on_fasta.sh --jobs 8 core.meme sv_fastas fimo_out

JOBS=1
THRESH="1e-4"
FORCE=0

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2 ;;
    --thresh) THRESH="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help)
      echo "Usage: $0 [--jobs N] [--thresh X] [--force] <meme_file> <fasta_dir> <out_dir>" >&2
      exit 0
      ;;
    *)
      POSITIONAL+=("$1"); shift ;;
  esac
done

if [[ ${#POSITIONAL[@]} -ne 3 ]]; then
  echo "Usage: $0 [--jobs N] [--thresh X] [--force] <meme_file> <fasta_dir> <out_dir>" >&2
  exit 1
fi

MEME_FILE="${POSITIONAL[0]}"
FASTA_DIR="${POSITIONAL[1]}"
OUT_DIR="${POSITIONAL[2]}"

mkdir -p "$OUT_DIR"

command -v fimo >/dev/null 2>&1 || { echo "Error: fimo not found in PATH." >&2; exit 1; }
[[ -f "$MEME_FILE" ]] || { echo "Error: MEME file not found: $MEME_FILE" >&2; exit 1; }
[[ -d "$FASTA_DIR" ]] || { echo "Error: FASTA directory not found: $FASTA_DIR" >&2; exit 1; }
[[ "$JOBS" =~ ^[0-9]+$ ]] || { echo "Error: --jobs must be integer" >&2; exit 1; }

# Build a list of fasta files
shopt -s nullglob
fastas=( "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.fasta )
shopt -u nullglob

if [[ ${#fastas[@]} -eq 0 ]]; then
  echo "No FASTA files found in: $FASTA_DIR"
  exit 0
fi

export MEME_FILE OUT_DIR THRESH FORCE

run_one_fa() {
  local fa="$1"
  local base
  base="$(basename "$fa")"
  base="${base%.fa}"
  base="${base%.fasta}"

  local out_tsv="$OUT_DIR/${base}_fimo.tsv"
  if [[ -f "$out_tsv" && "$FORCE" -eq 0 ]]; then
    echo "[skip] $out_tsv exists"
    return 0
  fi

  local tmp_oc="$OUT_DIR/${base}_fimo_run"
  rm -rf "$tmp_oc"
  mkdir -p "$tmp_oc"

  echo "[run ] FIMO: $fa → $out_tsv"
  fimo --oc "$tmp_oc" --thresh "$THRESH" "$MEME_FILE" "$fa"

  if [[ -f "$tmp_oc/fimo.tsv" ]]; then
    mv "$tmp_oc/fimo.tsv" "$out_tsv"
    echo "[ok  ] $out_tsv"
  else
    echo "[warn] fimo.tsv not found for $fa" >&2
  fi

  # Clean temp directory
  rm -rf "$tmp_oc"
}
export -f run_one_fa

# Parallelize across FASTA files
printf '%s\0' "${fastas[@]}" | xargs -0 -n 1 -P "$JOBS" bash -c 'run_one_fa "$0"' 

echo "Done. All FIMO TSVs collected in: $OUT_DIR"
