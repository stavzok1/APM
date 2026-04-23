#!/usr/bin/env bash
set -euo pipefail

# Master SV → motif pipeline
#
# Usage:
#   ./run_sv_motif_pipeline.sh \
#       SCRIPTS_DIR \
#       RAW_VCF_DIR \
#       OUTPUT_ROOT \
#       GENE_FEATURES \
#       SAMPLES_TSV \
#       REG_ELEM_FILE \
#       LNCRNA_FEATURES \
#       REF_FASTA \
#       MEME_FILE \
#       FLANK_SIZE
#
# Example:
#   ./run_sv_motif_pipeline.sh \
#       sv_pipeline/scripts \
#       SV_vcfs \
#       SV_output \
#       genes.csv \
#       samples.tsv \
#       elements.csv \
#       lncrnas.csv \
#       ~/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#       motifs.meme \
#       150

if [ "$#" -ne 10 ]; then
    echo "Usage: $0 SCRIPTS_DIR RAW_VCF_DIR OUTPUT_ROOT GENE_FEATURES SAMPLES_TSV REG_ELEM_FILE LNCRNA_FEATURES REF_FASTA MEME_FILE FLANK_SIZE" >&2
    exit 1
fi

SCRIPTS_DIR="$1"
RAW_VCF_DIR="$2"
OUTPUT_ROOT="$3"
GENE_FEATURES="$4"
SAMPLES_TSV="$5"
REG_ELEM_FILE="$6"
LNCRNA_FEATURES="$7"
REF_FASTA="$8"
MEME_FILE="$9"
FLANK_SIZE="${10}"

# Normalize SCRIPTS_DIR
SCRIPTS_DIR=$(realpath "$SCRIPTS_DIR")

echo "[INFO] Using scripts from: $SCRIPTS_DIR"

# Internal directory structure
VEP_VCF_DIR="${OUTPUT_ROOT}/01_vep_vcfs"
PROC_CSV_DIR="${OUTPUT_ROOT}/02_processed_sv_csv"
SV_BED_DIR="${OUTPUT_ROOT}/03_sv_bed"
SV_FASTA_DIR="${OUTPUT_ROOT}/04_sv_fasta"
FIMO_TSV_DIR="${OUTPUT_ROOT}/05_fimo_tsv"
FIMO_MERGED_DIR="${OUTPUT_ROOT}/06_sv_fimo_merged"
FINAL_CSV_DIR="${OUTPUT_ROOT}/07_final_sv_with_fimo"
LOG_DIR="${OUTPUT_ROOT}/logs"

mkdir -p "$VEP_VCF_DIR" "$PROC_CSV_DIR" "$SV_BED_DIR" "$SV_FASTA_DIR" "$FIMO_TSV_DIR" "$FIMO_MERGED_DIR" "$FINAL_CSV_DIR" "$LOG_DIR"

################################################################################
# STEP 1 — VEP annotation
################################################################################
echo "[STEP 1] Running VEP batch..."
bash "${SCRIPTS_DIR}/run_vep_batch.sh" \
    "$RAW_VCF_DIR" \
    "$VEP_VCF_DIR" \
    2>&1 | tee "${LOG_DIR}/01_vep.log"

echo "[STEP 1] Completed"
echo

################################################################################
# STEP 2 — Process VEP VCFs → strict_sv_set.csv
################################################################################
echo "[STEP 2] Processing VEP VCFs..."

# Find the sv_pipeline Python module directory
PYTHON_MODULE_DIR=$(dirname "$SCRIPTS_DIR")

python3 -c "
import sys
sys.path.insert(0, '$(dirname $PYTHON_MODULE_DIR)')
from sv_pipeline.pipeline import run_sv_vcf_processing
from pathlib import Path

run_sv_vcf_processing(
    vcf_dir=Path('$VEP_VCF_DIR'),
    output_dir=Path('$PROC_CSV_DIR'),
    genes_path=Path('$GENE_FEATURES'),
    elements_path=Path('$REG_ELEM_FILE'),
    samples_tsv_path=Path('$SAMPLES_TSV'),
    lncrnas_path=Path('$LNCRNA_FEATURES'),
)
" 2>&1 | tee "${LOG_DIR}/02_process.log"

echo "[STEP 2] Completed"
echo

################################################################################
# STEP 3 — Produce SV BED intervals
################################################################################
echo "[STEP 3] Creating BED intervals..."

python3 -c "
import sys
sys.path.insert(0, '$(dirname $PYTHON_MODULE_DIR)')
from sv_pipeline.bed_intervals import create_beds_from_directory
from pathlib import Path

create_beds_from_directory(
    csv_dir=Path('$PROC_CSV_DIR'),
    output_dir=Path('$SV_BED_DIR'),
    flank=$FLANK_SIZE,
)
" 2>&1 | tee "${LOG_DIR}/03_bed.log"

echo "[STEP 3] Completed"
echo

################################################################################
# STEP 4 — FASTA extraction for FIMO
################################################################################
echo "[STEP 4] Extracting FASTA from SV BEDs..."
bash "${SCRIPTS_DIR}/create_fasta_for_sv_bed.sh" \
    "$REF_FASTA" \
    "$SV_BED_DIR" \
    "$SV_FASTA_DIR" \
    2>&1 | tee "${LOG_DIR}/04_fasta.log"

echo "[STEP 4] Completed"
echo

################################################################################
# STEP 5 — Run FIMO on FASTA
################################################################################
echo "[STEP 5] Running FIMO..."
bash "${SCRIPTS_DIR}/run_fimo_on_fasta.sh" \
    "$MEME_FILE" \
    "$SV_FASTA_DIR" \
    "$FIMO_TSV_DIR" \
    2>&1 | tee "${LOG_DIR}/05_fimo.log"

echo "[STEP 5] Completed"
echo

################################################################################
# STEP 6 — Intersect FIMO hits with SV BEDs
################################################################################
echo "[STEP 6] Intersecting FIMO hits with SV BED intervals..."
bash "${SCRIPTS_DIR}/run_fimo_sv_intersect.sh" \
    "$FIMO_TSV_DIR" \
    "$SV_BED_DIR" \
    "$FIMO_MERGED_DIR" \
    2>&1 | tee "${LOG_DIR}/06_intersect.log"

echo "[STEP 6] Completed"
echo

################################################################################
# STEP 7 — Recombine motif hits into processed SV CSVs
################################################################################
echo "[STEP 7] Recombining FIMO hits into SV CSVs..."

python3 -c "
import sys
sys.path.insert(0, '$(dirname $PYTHON_MODULE_DIR)')
from sv_pipeline.motif_scanning import recombine_all_sv_fimo
from pathlib import Path

recombine_all_sv_fimo(
    merged_bed_dir=Path('$FIMO_MERGED_DIR'),
    sv_csv_dir=Path('$PROC_CSV_DIR'),
    output_csv_dir=Path('$FINAL_CSV_DIR'),
)
" 2>&1 | tee "${LOG_DIR}/07_recombine.log"

echo "[STEP 7] Completed"
echo
echo "[INFO] Final annotated files: $FINAL_CSV_DIR"
echo
echo "[INFO] Pipeline complete!"
