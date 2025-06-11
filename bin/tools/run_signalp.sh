#!/bin/bash

# Usage: ./run_signalp.sh <input_fasta> <output_dir> [batch_size]
# Example: ./run_signalp.sh input.fasta results 10000

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <input_fasta> <output_dir> [batch_size]"
  exit 1
fi

FASTA="$1"
OUTPUT_DIR="$2"
BATCH_SIZE="${3:-10000}"

mkdir -p "$OUTPUT_DIR/tmp"

# Extract the first FASTA header to use as a prefix (e.g., antigen_77)
FIRST_HEADER=$(grep '^>' "$FASTA" | head -n1 | cut -d '|' -f1 | sed 's/^>//')
PREFIX="${FIRST_HEADER}"

# Path setup for SignalP
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIGNALP_DIR="${SCRIPT_DIR}/signalp-5.0b"
SIGNALP_BIN="${SIGNALP_DIR}/bin/signalp"
EXPECTED_BIN="${SIGNALP_DIR}/bin/bin/signalp"

if [ ! -f "$EXPECTED_BIN" ]; then
  mkdir -p "$(dirname "$EXPECTED_BIN")"
  cp "$SIGNALP_BIN" "$EXPECTED_BIN"
fi

echo "[INFO] Running SignalP on: $FASTA"
echo "[INFO] Output directory: $OUTPUT_DIR"
echo "[INFO] Output prefix: $PREFIX"

"$EXPECTED_BIN" \
  -fasta "$FASTA" \
  -format long \
  -mature \
  -prefix "$PREFIX" \
  -batch "$BATCH_SIZE" \
  -stdout \
  -tmp "$OUTPUT_DIR/tmp" \
  > "$OUTPUT_DIR/${PREFIX}_signalp_phred.txt"

# If a plot is generated, rename it
PLOT_FILE="${OUTPUT_DIR}/${PREFIX}_plot.png"
if compgen -G "${OUTPUT_DIR}/*_plot.png" > /dev/null; then
  mv "${OUTPUT_DIR}"/*_plot.png "$PLOT_FILE"
  echo "[INFO] Plot saved to: $PLOT_FILE"
fi

echo "[INFO] SignalP results written to: $OUTPUT_DIR/${PREFIX}_signalp_phred.txt"
