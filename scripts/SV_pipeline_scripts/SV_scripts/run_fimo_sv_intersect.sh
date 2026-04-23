#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_fimo_sv_intersect.sh fimo_tsv_dir sv_bed_dir out_dir
#
# Assumes:
#   FIMO TSVs:  fimo_tsv_dir/{base}_fimo.tsv
#   SV BEDs:    sv_bed_dir/{base}.bed
#
# Outputs (in out_dir):
#   {base}_fimo_hits.bed
#   {base}_fimo_merged.bed

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fimo_tsv_dir> <sv_bed_dir> <out_dir>" >&2
    exit 1
fi

FIMO_DIR="$1"
SV_BED_DIR="$2"
OUT_DIR="$3"

mkdir -p "$OUT_DIR"

echo "[INFO] FIMO TSV DIR: $FIMO_DIR"
echo "[INFO] SV BED DIR  : $SV_BED_DIR"
echo "[INFO] OUT DIR     : $OUT_DIR"

shopt -s nullglob

for FIMO_TSV in "$FIMO_DIR"/*_fimo.tsv; do
    [ -e "$FIMO_TSV" ] || { echo "[WARN] No *_fimo.tsv files in $FIMO_DIR"; break; }

    base=$(basename "$FIMO_TSV")
    base="${base%_fimo.tsv}"

    SV_BED="$SV_BED_DIR/${base}.bed"

    if [ ! -f "$SV_BED" ]; then
        echo "[WARN] No matching BED for base '${base}' at: $SV_BED — skipping."
        continue
    fi

    FIMO_HITS_BED="$OUT_DIR/${base}_fimo_hits.bed"
    MERGED_BED="$OUT_DIR/${base}_fimo_merged.bed"

    echo
    echo "[INFO] Processing base: $base"
    echo "[INFO]   FIMO TSV: $FIMO_TSV"
    echo "[INFO]   SV BED  : $SV_BED"
    echo "[INFO]   Hits BED: $FIMO_HITS_BED"
    echo "[INFO]   Merged  : $MERGED_BED"

    # ----------------------------------------------------------
    # Step 1: Python – convert {base}_fimo.tsv → {base}_fimo_hits.bed
    # ----------------------------------------------------------
    python3 << EOF
import pandas as pd

fimo_tsv_path = "${FIMO_TSV}"
out_bed_path  = "${FIMO_HITS_BED}"

print(f"[PY] Reading FIMO TSV: {fimo_tsv_path}")

df = pd.read_csv(fimo_tsv_path, sep="\\t")
print("[PY] Columns:", list(df.columns))
print("[PY] N rows before:", len(df))

# Coerce start/stop
df["start_num"] = pd.to_numeric(df["start"], errors="coerce")
df["stop_num"]  = pd.to_numeric(df["stop"],  errors="coerce")

df = df.dropna(subset=["start_num", "stop_num"]).copy()
print("[PY] N rows after dropping bad coords:", len(df))

# Clean up chromosome / sequence_name
chrom = df["sequence_name"].astype(str)

# strip trailing '.0' (e.g. '17.0' -> '17')
chrom = chrom.str.replace(r"\\.0$", "", regex=True)

# strip 'chr' prefix if present, to match your SV BEDs
chrom = chrom.str.replace(r"^chr", "", regex=True)

bed = pd.DataFrame({
    "chrom": chrom,
    "start": (df["start_num"].astype(int) - 1).clip(lower=0),
    "end":   df["stop_num"].astype(int),
    "motif_id": df["motif_id"],
    "p_value": df["p-value"],
    "q_value": df["q-value"],
    "score" : df["score"],
    "strand": df["strand"],
})

bed.to_csv(out_bed_path, sep="\\t", header=False, index=False)
print(f"[PY] Wrote BED file: {out_bed_path}")
EOF

    echo "[INFO] Head of ${FIMO_HITS_BED}:"
    head -5 "$FIMO_HITS_BED" | cat -t || true

    # ----------------------------------------------------------
    # Step 2: bedtools intersect – Join SV intervals with motif hits
    # ----------------------------------------------------------
    echo "[INFO] Running bedtools intersect for base: $base"

    bedtools intersect \
        -a "$SV_BED" \
        -b "$FIMO_HITS_BED" \
        -wa -wb \
        > "$MERGED_BED"

    echo "[INFO] Wrote merged BED: $MERGED_BED"

done

echo
echo "[INFO] Done processing all FIMO TSVs."
