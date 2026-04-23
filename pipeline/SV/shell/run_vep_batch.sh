#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_vep_batch.sh <input_vcf_dir> <output_dir> [--fork N] [--force]

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_vcf_dir> <output_dir> [--fork N] [--force]" >&2
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
shift 2

# Defaults
FORK=4
FORCE=""

# Parse optional args
while [[ $# -gt 0 ]]; do
    case $1 in
        --fork)
            FORK="$2"
            shift 2
            ;;
        --force)
            FORCE="--force_overwrite"
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

mkdir -p "$OUTPUT_DIR"

# Loop through all VCF files
for vcf in "$INPUT_DIR"/*.vcf; do
    [ -e "$vcf" ] || continue

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
      $FORCE \
      --everything \
      --symbol \
      --variant_class \
      --check_ref \
      --fork "$FORK"

done

echo "All VCFs processed. Output written to: $OUTPUT_DIR"
