#!/bin/bash

# Resolve the directory where this script is located (to handle relative paths correctly)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIGNALP_DIR="${SCRIPT_DIR}/signalp-5.0b"
SIGNALP_BIN="${SIGNALP_DIR}/bin/signalp"
SIGNALP_EXPECTED_BIN="${SIGNALP_DIR}/bin/bin/signalp"

# Step 1: Ensure the expected bin/bin/signalp path exists
if [ ! -f "$SIGNALP_EXPECTED_BIN" ]; then
    mkdir -p "${SIGNALP_DIR}/bin/bin"
    cp "$SIGNALP_BIN" "$SIGNALP_EXPECTED_BIN"
fi

# Step 2: Run SignalP with passed arguments
"$SIGNALP_EXPECTED_BIN" "$@"
#!/bin/bash

# Usage: ./run_signalp.sh <input_fasta> <output_dir> <organism> <prefix> [batch_size]
# Example: ./run_signalp.sh input.fasta results euk myrun 10000

set -euo pipefail

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <input_fasta> <output_dir> <organism: euk|gram+|gram-|arch> <prefix> [batch_size]"
  exit 1
fi

FASTA=$1
OUTPUT_DIR=$2
ORGANISM=$3
PREFIX=$4
BATCH_SIZE=${5:-10000}

mkdir -p "$OUTPUT_DIR"

signalp \
  -fasta "$FASTA" \
  -org "$ORGANISM" \
  -format short \
  -gff \
  -mature \
  -prefix "$PREFIX" \
  -batch "$BATCH_SIZE" \
  -stdout \
  -tmp "$OUTPUT_DIR/tmp" \
  > "$OUTPUT_DIR/${PREFIX}_summary.txt"
