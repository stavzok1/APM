#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_fimo_sv_intersect.sh fimo_out_sv/fimo.tsv attempt2_strict_sv_set_flanks_and_overlaps.bed fimo_out_sv

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fimo_tsv> <sv_bed> <out_dir>" >&2
    exit 1
fi

FIMO_TSV="$1"
SV_BED="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

echo "[INFO] FIMO TSV: $FIMO_TSV"
echo "[INFO] SV BED  : $SV_BED"
echo "[INFO] OUT DIR : $OUT_DIR"

# --------------------------------------------------------------
# Step 1: Python – convert fimo.tsv → fimo_hits.bed
# --------------------------------------------------------------

python3 << EOF
import pandas as pd

fimo_tsv_path = "${FIMO_TSV}"
out_bed_path  = "${OUT_DIR}/fimo_hits.bed"

print(f"[PY] Reading FIMO TSV: {fimo_tsv_path}")

df = pd.read_csv(fimo_tsv_path, sep="\\t")
print("[PY] Columns:", list(df.columns))
print("[PY] N rows before:", len(df))

# Coerce start/stop
df["start_num"] = pd.to_numeric(df["start"], errors="coerce")
df["stop_num"]  = pd.to_numeric(df["stop"],  errors="coerce")

df = df.dropna(subset=["start_num", "stop_num"]).copy()
print("[PY] N rows after dropping bad coords:", len(df))

bed = pd.DataFrame({
    "chrom": df["sequence_name"].astype(str),
    "start": (df["start_num"].astype(int) - 1).clip(lower=0),
    "end":   df["stop_num"].astype(int),
    "motif_id": df["motif_id"],
    "p_value": df["p-value"],
    "strand": df["strand"],
})

bed.to_csv(out_bed_path, sep="\\t", header=False, index=False)
print(f"[PY] Wrote BED file: {out_bed_path}")
EOF

echo "[INFO] Head of fimo_hits.bed:"
head -5 "$OUT_DIR/fimo_hits.bed" | cat -t || true

# --------------------------------------------------------------
# Step 2: bedtools intersect – Join SV intervals with motif hits
# --------------------------------------------------------------

echo "[INFO] Running bedtools intersect..."

bedtools intersect \
    -a "$SV_BED" \
    -b "$OUT_DIR/fimo_hits.bed" \
    -wa -wb \
    > "$OUT_DIR/sv_fimo_merged.bed"

echo "[INFO] Wrote: $OUT_DIR/sv_fimo_merged.bed"
echo "[INFO] Done."
