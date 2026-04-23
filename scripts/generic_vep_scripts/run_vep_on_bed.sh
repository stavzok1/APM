#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $0 -i <vcf_dir> -o <out_dir> -b <regions.bed> [-f <fasta>] [-p <glob>] [-- <vep_args...>]

Required:
  -i    Input directory containing VCFs
  -o    Output directory
  -b    BED file of regions to keep (bedtools intersect)

Optional:
  -f    FASTA for VEP (default: \$HOME/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa)
  -p    Glob pattern for VCFs (default: *.vcf)
  --    Everything after -- is passed directly to vep

Examples:
  $0 -i vcfs_SNV -o vep_vcfs -b apm_genes_1Mb.bed
  $0 -i vcfs_SNV -o vep_vcfs -b apm_genes_1Mb.bed -f /path/genome.fa -- --fork 4 --buffer_size 5000
EOF
}

VCF_DIR=""
OUT_DIR=""
REGION_BED=""
FASTA="${HOME}/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
PATTERN="*.vcf"

# Parse flags (supports "--" passthrough)
VEP_EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) VCF_DIR="$2"; shift 2 ;;
    -o) OUT_DIR="$2"; shift 2 ;;
    -b) REGION_BED="$2"; shift 2 ;;
    -f) FASTA="$2"; shift 2 ;;
    -p) PATTERN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --)
      shift
      VEP_EXTRA_ARGS=("$@")
      break
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

# Validate required args
if [[ -z "$VCF_DIR" || -z "$OUT_DIR" || -z "$REGION_BED" ]]; then
  echo "Error: missing required arguments."
  usage
  exit 1
fi

# Validate inputs exist
if [[ ! -d "$VCF_DIR" ]]; then
  echo "Error: input VCF directory does not exist: $VCF_DIR"
  exit 1
fi
if [[ ! -f "$REGION_BED" ]]; then
  echo "Error: BED file does not exist: $REGION_BED"
  exit 1
fi
if [[ ! -f "$FASTA" ]]; then
  echo "Error: FASTA file does not exist: $FASTA"
  exit 1
fi

# Check dependencies
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found in PATH"; exit 1; }
command -v vep      >/dev/null 2>&1 || { echo "Error: vep not found in PATH"; exit 1; }

mkdir -p "$OUT_DIR"

shopt -s nullglob
vcfs=( "$VCF_DIR"/$PATTERN )
shopt -u nullglob

if [[ ${#vcfs[@]} -eq 0 ]]; then
  echo "No files matched '$PATTERN' in $VCF_DIR"
  exit 0
fi

for vcf in "${vcfs[@]}"; do
  filename="$(basename "$vcf")"
  base="${filename%.vcf}"

  echo "Processing: $filename"

  # 1) Filter by regions
  intersect_vcf="${OUT_DIR}/${base}.filtered.vcf"
  bedtools intersect -header -a "$vcf" -b "$REGION_BED" > "$intersect_vcf"

  # 2) Run VEP
  out_vep="${OUT_DIR}/${base}.filtered.vep.vcf"
  echo "  Running VEP on: $(basename "$intersect_vcf")"

  vep \
    --offline \
    --cache \
    --assembly GRCh38 \
    --everything \
    --fasta "$FASTA" \
    --vcf \
    --input_file "$intersect_vcf" \
    --output_file "$out_vep" \
    "${VEP_EXTRA_ARGS[@]}"

  echo "  Done: $out_vep"
  echo
done

echo "All VCFs processed."
