#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <TCGA_ATAC_PeakSet.bed> <ENCODE_sample.bam> <out_dir>"
    exit 1
fi

PEAKSET="$1"
BAM="$2"
OUTDIR="$3"

mkdir -p "$OUTDIR"

PEAK_BASENAME="$(basename "$PEAKSET")"
PEAK_NOHEADER="${OUTDIR}/${PEAK_BASENAME%.bed}.noheader.bed"
SORTED_PEAK="${OUTDIR}/${PEAK_BASENAME%.bed}.sorted.bed"

BAM_BASENAME="$(basename "$BAM")"
PREFIX="${BAM_BASENAME%.bam}"

SORTED_BAM="${OUTDIR}/${PREFIX}.sorted.bam"
GENOME="${OUTDIR}/${PREFIX}.from_bam.genome"
COUNTS="${OUTDIR}/${PREFIX}.counts.txt"
VECTOR="${OUTDIR}/${PREFIX}.vector"

echo ">>> Using peakset:        $PEAKSET"
echo ">>> Using BAM:            $BAM"
echo ">>> Output dir:           $OUTDIR"
echo ">>> Vector output:        $VECTOR"
echo

# 1) Remove header from peakset
# (Optional improvement: cache this once; see below)
echo "[1/7] Making no-header BED..."
if [[ ! -f "$PEAK_NOHEADER" ]]; then
  echo "[1/7] Making no-header BED..."
  tail -n +2 "$PEAKSET" > "$PEAK_NOHEADER"
else
  echo "[1/7] Using cached no-header BED: $PEAK_NOHEADER"
fi

# 2) Sort + index BAM
echo "[2/7] Sorting BAM..."
samtools sort -o "$SORTED_BAM" "$BAM"

echo "[3/7] Indexing BAM..."
samtools index "$SORTED_BAM"

# 3) Create genome file from BAM header
echo "[4/7] Creating genome file from BAM header..."
samtools view -H "$SORTED_BAM" \
  | grep '^@SQ' \
  | awk '{
      chr=""; len="";
      for (i=1; i<=NF; i++) {
          if ($i ~ /^SN:/) { split($i,a,":"); chr=a[2] }
          if ($i ~ /^LN:/) { split($i,a,":"); len=a[2] }
      }
      if (chr != "" && len != "") { print chr "\t" len }
  }' > "$GENOME"

# 4) Sort peakset BED according to genome order
echo "[5/7] Sorting peakset BED according to genome..."
bedtools sort -i "$PEAK_NOHEADER" -g "$GENOME" > "$SORTED_PEAK"

# 5) Run coverage
echo "[6/7] Running bedtools coverage..."
bedtools coverage -sorted -g "$GENOME" -a "$SORTED_PEAK" -b "$SORTED_BAM" -counts > "$COUNTS"

# 6) Extract vector
echo "[7/7] Extracting vector..."
cut -f5 "$COUNTS" > "$VECTOR"

echo "Done. Vector file: $VECTOR"
