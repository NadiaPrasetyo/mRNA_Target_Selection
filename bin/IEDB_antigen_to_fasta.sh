#!/bin/bash

# === USAGE CHECK ===
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 [curl|wget] input.[xlsx|csv] -o output.fasta"
    exit 1
fi

# === CONFIGURATION ===
DOWNLOAD_TOOL="wget"
if [[ "$1" == "curl" || "$1" == "wget" ]]; then
    DOWNLOAD_TOOL="$1"
    shift
fi

# === INPUT PARSING ===
INPUT_FILE="$1"
shift
if [[ "$1" != "-o" || -z "$2" ]]; then
    echo "Error: Missing or invalid output option. Usage: $0 [curl|wget] input.xlsx -o output.fasta"
    exit 1
fi
FASTA_OUTPUT="$2"
shift 2

PARSE_SCRIPT="bin/parse_input_txt.py"
TEMP_TXT_FILE=$(mktemp)

# === CHECK ENTREZ DIRECT ===
if ! command -v efetch &> /dev/null; then
    echo "[INFO] Entrez Direct (efetch) is not installed."
    read -p "Do you want to install it now? [y/N] " INSTALL_CHOICE

    if [[ "$INSTALL_CHOICE" =~ ^[Yy]$ ]]; then
        echo "[INFO] Installing Entrez Direct using $DOWNLOAD_TOOL..."
        if [[ "$DOWNLOAD_TOOL" == "curl" ]]; then
            sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
        else
            sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
        fi
        export PATH=$PATH:$HOME/edirect
    else
        echo "[ABORTED] Please install Entrez Direct manually and try again."
        exit 1
    fi
fi

# === PARSE INPUT FILE TO TEXT FILE ===
echo "[INFO] Extracting UniProt IDs from '$INPUT_FILE'..."
python3 "$PARSE_SCRIPT" --input "$INPUT_FILE" --output "$TEMP_TXT_FILE"

if [[ ! -f "$TEMP_TXT_FILE" ]]; then
    echo "[ERROR] Temporary UniProt ID file not created."
    exit 1
fi

# === FETCH FASTA SEQUENCES ===
echo "[INFO] Fetching FASTA sequences from NCBI..."
grep -v '^#' "$TEMP_TXT_FILE" | efetch -db protein -format fasta > "$FASTA_OUTPUT"

if [[ $? -eq 0 ]]; then
    echo "[DONE] FASTA sequences saved to: $FASTA_OUTPUT"
else
    echo "[ERROR] Failed to fetch FASTA sequences."
    exit 1
fi

# === CLEAN UP ===
rm -f "$TEMP_TXT_FILE"
echo "[INFO] Temporary files cleaned up."
echo "[INFO] Script completed successfully."
exit 0
