#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./create_fasta_for_sv_bed.sh <reference.fa> <bed_dir> <output_dir>

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fa> <bed_dir> <output_dir>" >&2
    exit 1
fi

FASTA="$1"
BED_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

for bed in "$BED_DIR"/*.bed; do
    [ -e "$bed" ] || continue
    
    base=$(basename "$bed" .bed)
    out_fa="$OUT_DIR/$base.fa"

    echo "Processing $bed → $out_fa"

    bedtools getfasta \
        -fi "$FASTA" \
        -bed "$bed" \
        -name \
        > "$out_fa"
done

echo "Done. FASTA files saved to $OUT_DIR"
