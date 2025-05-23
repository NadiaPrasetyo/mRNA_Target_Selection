#!/bin/bash

# === USAGE CHECK ===
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 [curl|wget] input.xlsx -o output.txt"
    exit 1
fi

# === CONFIGURATION ===
DOWNLOAD_TOOL="wget"
if [[ "$1" == "curl" || "$1" == "wget" ]]; then
    DOWNLOAD_TOOL="$1"
    shift
fi

# === INPUT PARSING ===
INPUT_EXCEL="$1"
shift
if [[ "$1" != "-o" || -z "$2" ]]; then
    echo "Error: Missing or invalid output option. Usage: $0 [curl|wget] input.xlsx -o output.txt"
    exit 1
fi
OUTPUT_TXT="$2"
shift 2

PARSE_SCRIPT="bin/parse_xlsx_txt.py"
FASTA_OUTPUT="antigens.fasta"

# === 1. CHECK ENTREZ DIRECT ===
if ! command -v efetch &> /dev/null; then
    echo "Entrez Direct (efetch) is not installed."
    read -p "Do you want to install it now? [y/N] " INSTALL_CHOICE

    if [[ "$INSTALL_CHOICE" =~ ^[Yy]$ ]]; then
        echo "Installing Entrez Direct using $DOWNLOAD_TOOL..."
        if [[ "$DOWNLOAD_TOOL" == "curl" ]]; then
            sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
        else
            sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
        fi
        export PATH=$PATH:$HOME/edirect
    else
        echo "Aborting. Please install Entrez Direct manually and try again."
        exit 1
    fi
fi

# === 2. PARSE EXCEL TO TEXT ===
echo "Running Python script to extract UniProt IDs..."
python3 "$PARSE_SCRIPT" --input "$INPUT_EXCEL" --output "$OUTPUT_TXT"

# === 3. FETCH FASTA SEQUENCES ===
if [[ ! -f "$OUTPUT_TXT" ]]; then
    echo "Error: $OUTPUT_TXT not found. Make sure the parsing step completed successfully."
    exit 1
fi

echo "Fetching FASTA sequences for UniProt IDs..."
grep -v '^#' "$OUTPUT_TXT" | efetch -db protein -format fasta > "$FASTA_OUTPUT"

echo "Saved FASTA sequences to $FASTA_OUTPUT"
