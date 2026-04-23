#!/bin/bash

# Usage:
#   ./create_fasta_for_sv.sh reference.fa bed_dir output_dir

FASTA="$1"
BED_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

for bed in "$BED_DIR"/*.bed; do
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
