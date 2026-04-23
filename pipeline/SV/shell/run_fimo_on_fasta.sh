#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_fimo_on_fasta.sh <meme_file> <fasta_dir> <out_dir> [--thresh X] [--jobs N]

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <meme_file> <fasta_dir> <out_dir> [--thresh X] [--jobs N]" >&2
    exit 1
fi

MEME_FILE="$1"
FASTA_DIR="$2"
OUT_DIR="$3"
shift 3

# Defaults
THRESH="1e-4"
JOBS=1

# Parse optional args
while [[ $# -gt 0 ]]; do
    case $1 in
        --thresh)
            THRESH="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

mkdir -p "$OUT_DIR"

# Check FIMO availability
if ! command -v fimo >/dev/null 2>&1; then
    echo "Error: fimo not found in PATH." >&2
    exit 1
fi

if [ ! -f "$MEME_FILE" ]; then
    echo "Error: MEME file not found: $MEME_FILE" >&2
    exit 1
fi

if [ ! -d "$FASTA_DIR" ]; then
    echo "Error: FASTA directory not found: $FASTA_DIR" >&2
    exit 1
fi

shopt -s nullglob

for fa in "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.fasta; do
    [ -e "$fa" ] || continue

    base=$(basename "$fa")
    base="${base%.fa}"
    base="${base%.fasta}"

    tmp_oc="$OUT_DIR/${base}_fimo_run"

    # Ensure clean tmp dir
    rm -rf "$tmp_oc"
    mkdir -p "$tmp_oc"

    echo "Running FIMO for $fa → $tmp_oc"

    fimo \
        --oc "$tmp_oc" \
        --thresh "$THRESH" \
        "$MEME_FILE" \
        "$fa"

    if [ -f "$tmp_oc/fimo.tsv" ]; then
        mv "$tmp_oc/fimo.tsv" "$OUT_DIR/${base}_fimo.tsv"
        echo "  Saved: $OUT_DIR/${base}_fimo.tsv"
    else
        echo "  WARNING: fimo.tsv not found for $fa" >&2
    fi

    # Cleanup
    find "$tmp_oc" -type f ! -name '*.tsv' -delete || true
    rmdir "$tmp_oc" 2>/dev/null || true
done

echo "Done. All FIMO TSVs collected in: $OUT_DIR"
