#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_vep_batch.sh <input_vcf_dir> <output_dir>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf_dir> <output_dir>" >&2
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

# Loop through all VCF files in the input directory
for vcf in "$INPUT_DIR"/*.vcf; do
    [ -e "$vcf" ] || continue  # skip if none found

    base=$(basename "$vcf")
    out_file="$OUTPUT_DIR/$base"

    echo "Running VEP on: $base"

    vep \
      --input_file "$vcf" \
      --output_file "$out_file" \
      --offline \
      --cache \
      --dir_cache ~/.vep \
      --species homo_sapiens \
      --assembly GRCh38 \
      --fasta ~/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      --format vcf \
      --vcf \
      --force_overwrite \
      --everything \
      --symbol \
      --variant_class \
      --check_ref \
      --fork 4

done

echo "All VCFs processed. Output written to: $OUTPUT_DIR"
