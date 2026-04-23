#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# Usage check
# ----------------------------
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input_folder> <output_folder>"
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# ----------------------------
# Safety checks
# ----------------------------
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input folder does not exist: $INPUT_DIR"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ----------------------------
# Extract all .gz files
# ----------------------------
find "$INPUT_DIR" -type f -name '*.gz' -print0 |
while IFS= read -r -d '' f; do
  base="$(basename "${f%.gz}")"
  out="$OUTPUT_DIR/$base"

  echo "Extracting: $f → $out"
  gunzip -c "$f" > "$out"
done

echo "Done."
