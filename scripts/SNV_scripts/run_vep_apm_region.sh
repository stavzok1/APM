#!/usr/bin/env bash
set -euo pipefail

# Paths (relative to repo root)
VCF_DIR="data/SNV/vcfs_SNV"
OUT_DIR="data/SNV/vep_vcfs"
REGION_BED="data/SNV/apm_genes_1Mb.bed"
FASTA="$HOME/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Create output dir if it doesn't exist
mkdir -p "${OUT_DIR}"

# Build REGION_BED if missing (uses PIPELINE_GENE_PANEL + SNV/CNV stratifiers from pipeline/config.py)
if [ ! -s "${REGION_BED}" ]; then
    echo "Building region bed: ${REGION_BED}"
    python -m pipeline.SNV.apm_region_bed --out "${REGION_BED}" --flank-bp 1000000
fi

# Loop over all VCFs in VCF_DIR
for vcf in "${VCF_DIR}"/*.vcf; do
    # Skip if no files match
    [ -e "$vcf" ] || { echo "No VCF files found in ${VCF_DIR}"; break; }

    # Get base name without path and extension
    filename=$(basename "$vcf")                 # e.g. TCGA-BRCA.xxx.vcf
    base="${filename%.vcf}"                     # e.g. TCGA-BRCA.xxx

    echo "Processing: $filename"

    # 1) Intersect with ±1Mb APM regions
    #    We keep only variants that overlap the regions in REGION_BED.
    #    -header keeps the VCF header intact.
    intersect_vcf="${OUT_DIR}/${base}.APM_1Mb.vcf"
    bedtools intersect \
        -header \
        -a "$vcf" \
        -b "$REGION_BED" \
        > "$intersect_vcf"

    # 2) Run VEP on the intersected VCF
    out_vep="${OUT_DIR}/${base}.APM_1Mb.vep.vcf"

    echo "  Running VEP on: $intersect_vcf"
    vep \
      --offline \
      --cache \
      --assembly GRCh38 \
      --everything \
      --fasta "$FASTA" \
      --vcf \
      --input_file "$intersect_vcf" \
      --output_file "$out_vep"

    echo "  Done: $out_vep"
    echo
done

echo "All VCFs processed."
