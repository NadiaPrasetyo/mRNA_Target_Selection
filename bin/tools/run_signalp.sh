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
