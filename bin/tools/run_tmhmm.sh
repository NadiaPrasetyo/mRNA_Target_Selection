#!/bin/bash

# Usage: ./configure_tmhmm.sh /path/to/TMHMM2.0a

set -e

# Input check
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/TMHMM2.0a"
    exit 1
fi

TMHMM_DIR=$(realpath "$1")
TMHMM_SCRIPT="$TMHMM_DIR/bin/tmhmm"

# Check existence
if [ ! -f "$TMHMM_SCRIPT" ]; then
    echo "Error: tmhmm script not found in $TMHMM_DIR/bin/"
    exit 1
fi

# Find required executables
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

# Backup original script
cp "$TMHMM_SCRIPT" "$TMHMM_SCRIPT.bak"

# Patch the script with correct paths
sed -i "s|^#!/usr/bin/perl|#!$PERL|" "$TMHMM_SCRIPT"

# Export necessary environment variables in script header
sed -i "1 a\
# Auto-configured paths\n\
TMHMM_LIB=\"$TMHMM_DIR/lib\"\n\
TMHMM_FORMAT=\"$TMHMM_DIR/bin/tmhmmformat.pl\"\n\
export TMHMM_LIB TMHMM_FORMAT\n" "$TMHMM_SCRIPT"

# Replace any hardcoded tmhmmformat.pl or lib path references in the script
sed -i "s|tmhmmformat.pl|\"$TMHMM_DIR/bin/tmhmmformat.pl\"|g" "$TMHMM_SCRIPT"
sed -i "s|TMHMM2.0.model|\"$TMHMM_DIR/lib/TMHMM2.0.model\"|g" "$TMHMM_SCRIPT"

chmod +x "$TMHMM_SCRIPT"

echo "✔ TMHMM configuration complete. Backup created at: $TMHMM_SCRIPT.bak"

#!/bin/bash

# Usage: ./run_tmhmm.sh <input_fasta> <output_dir> <prefix>
# Example: ./run_tmhmm.sh input.fasta results myrun

set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <input_fasta> <output_dir> <prefix>"
  exit 1
fi

FASTA=$1
OUTPUT_DIR=$2
PREFIX=$3

mkdir -p "$OUTPUT_DIR"

# You can switch this to decodeanhmm + model if needed
tmhmm "$FASTA" > "$OUTPUT_DIR/${PREFIX}_tmhmm.txt"
