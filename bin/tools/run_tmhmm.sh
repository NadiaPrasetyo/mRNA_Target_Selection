#!/bin/bash

set -euo pipefail

# --- Input check and setup patch ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 /path/to/TMHMM2.0a /path/to/input.fasta /path/to/output_dir"
    exit 1
fi

TMHMM_DIR=$(realpath "$1")
INPUT_FASTA=$(realpath "$2")
OUTPUT_DIR=$(realpath "$3")

TMHMM_SCRIPT="$TMHMM_DIR/bin/tmhmm"

# Check existence
if [ ! -f "$TMHMM_SCRIPT" ]; then
    echo "Error: tmhmm script not found in $TMHMM_DIR/bin/"
    exit 1
fi

PERL=$(which perl)

# Optional: validate required model and helper script files exist
REQUIRED_FILES=(
    "$TMHMM_DIR/lib/TMHMM2.0.model"
    "$TMHMM_DIR/bin/tmhmmformat.pl"
)

for f in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Missing required file: $f"
        exit 1
    fi
done

# Patch the script shebang (only if not already patched)
if ! head -1 "$TMHMM_SCRIPT" | grep -q "$PERL"; then
    cp "$TMHMM_SCRIPT" "$TMHMM_SCRIPT.bak"
    sed -i "s|^#!/usr/bin/perl|#!$PERL|" "$TMHMM_SCRIPT"
fi

mkdir -p "$OUTPUT_DIR"

BASENAME=$(basename "$INPUT_FASTA")
BASENAME="${BASENAME%.*}"

OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_tmhmm_result.txt"

# --- Run TMHMM in long format ---
"$TMHMM_SCRIPT" -long "$INPUT_FASTA" | \
  awk '
  /^#/ { print; next }
  {
    $1 = $1 "_tmhmm"
    print
  }' > "$OUTPUT_FILE"

echo "✔ TMHMM run complete. Output saved to $OUTPUT_FILE"
