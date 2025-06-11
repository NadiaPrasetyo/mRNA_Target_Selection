#!/bin/bash

# Usage: ./configure_targetp.sh /path/to/targetp-1.1

set -e

# Input check
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/targetp-1.1"
    exit 1
fi

TARGETP_DIR=$(realpath "$1")
TARGETP_SCRIPT="$TARGETP_DIR/targetp"
TMP_DIR="/tmp"  # You can change this to another default if preferred

# Check targetp file
if [ ! -f "$TARGETP_SCRIPT" ]; then
    echo "Error: targetp file not found in $TARGETP_DIR"
    exit 1
fi

# Locate executables
PASTE=$(which paste)
PERL=$(which perl)
AWK=$(which awk)
SH=$(which sh)
ECHO=$(which echo)

# Handle Linux echo quirk
if [[ "$(uname)" == "Linux" ]]; then
    ECHO="/bin/echo -e"
fi

# Ensure TMP has sticky bit
if [ ! -d "$TMP_DIR" ]; then
    echo "Creating TMP directory: $TMP_DIR"
    mkdir -p "$TMP_DIR"
fi

if [[ $(stat -c "%A" "$TMP_DIR") != drwxrwxrwt ]]; then
    echo "Setting sticky bit on $TMP_DIR"
    chmod 1777 "$TMP_DIR"
fi

# Create a backup of the targetp file
cp "$TARGETP_SCRIPT" "$TARGETP_SCRIPT.bak"

# Update configuration block
awk -v tp="$TARGETP_DIR" -v tmp="$TMP_DIR" -v paste="$PASTE" -v perl="$PERL" \
    -v awk="$AWK" -v sh="$SH" -v echo="$ECHO" '
BEGIN { section = 0 }
{
    if (/^### GENERAL SETTINGS, CUSTOMIZE/) section = 1
    if (section == 1 && /^# Substitute your chosen location for TargetP software:/) {
        getline; print "TARGETP=" tp; next
    }
    if (section == 1 && /^# determine where to store temporary files/) {
        getline; print "TMP=" tmp; next
    }
    if (section == 1 && /^# Substitute paste:/) {
        getline; print "PASTE=" paste; next
    }
    if (section == 1 && /^# Substitute perl:/) {
        while ($0 ~ /^#?PERL=/) getline
        print "PERL=" perl; next
    }
    if (section == 1 && /^# Substitute nawk, gawk or equivalent:/) {
        getline; print "AWK=" awk; next
    }
    if (section == 1 && /^# Substitute POSIX-compliant shell:/) {
        getline; print "SH=" sh; next
    }
    if (section == 1 && /^# Substitute echo:/) {
        getline; print "ECHO=" echo; next
        section = 2  # Done editing config
        next
    }
    print
}
' "$TARGETP_SCRIPT.bak" > "$TARGETP_SCRIPT"

chmod +x "$TARGETP_SCRIPT"

echo "✔ Configuration complete. Original backed up to: $TARGETP_SCRIPT.bak"
