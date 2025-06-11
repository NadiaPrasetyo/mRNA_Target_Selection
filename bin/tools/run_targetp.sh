#!/bin/bash

# Usage: ./run_targetp.sh <input_fasta> <output_dir> [batch_size]
# Example: ./run_targetp.sh input.fasta results 100

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input_fasta> <output_dir> [batch_size]"
  exit 1
fi

FASTA="$1"
OUTPUT_DIR="$2"
BATCH_SIZE="${3:-100}"

mkdir -p "$OUTPUT_DIR/tmp"

# Use first header to derive prefix
FIRST_HEADER=$(grep '^>' "$FASTA" | head -n1 | cut -d '|' -f1 | sed 's/^>//')
PREFIX="${FIRST_HEADER}"

# Locate TargetP binary
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGETP_BIN="${SCRIPT_DIR}/targetp-2.0/bin/targetp"

if [ ! -x "$TARGETP_BIN" ]; then
  if ! command -v targetp &> /dev/null; then
    echo "❌ Error: TargetP executable not found at $TARGETP_BIN or in PATH."
    exit 1
  fi
  TARGETP_BIN="targetp"
fi

echo "[INFO] Running TargetP on: $FASTA"
echo "[INFO] Output prefix: $PREFIX"
echo "[INFO] Output dir: $OUTPUT_DIR"

"$TARGETP_BIN" \
  -fasta "$FASTA" \
  -org non-pl \
  -format short \
  -batch "$BATCH_SIZE" \
  -gff3 \
  -mature \
  -prefix "$PREFIX" \
  -tmp "$OUTPUT_DIR/tmp" \
  -stdout \
  > "$OUTPUT_DIR/${PREFIX}_targetp.txt"

# Rename plot if generated
PLOT_FILE="${OUTPUT_DIR}/${PREFIX}_targetp_plot.png"
if compgen -G "${OUTPUT_DIR}/*_plot.png" > /dev/null; then
  mv "${OUTPUT_DIR}"/*_plot.png "$PLOT_FILE"
  echo "[INFO] Plot saved to: $PLOT_FILE"
fi

echo "[INFO] TargetP results written to: $OUTPUT_DIR/${PREFIX}_targetp.txt"
