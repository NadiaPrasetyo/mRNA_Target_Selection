#!/bin/bash

set -euo pipefail

INPUT_TSV="$1"
OUTPUT_FILE="$2"
TMP_DIR="intermediate/tmp_fasta_$(date +%s)"

# === SETUP ===
mkdir -p "$TMP_DIR"
> "$OUTPUT_FILE"

# === READ TSV & FETCH SEQUENCES ===
mapfile -t lines < <(tail -n +2 "$INPUT_TSV")

for line in "${lines[@]}"; do
    # Extract fields by splitting on tabs
    IFS=$'\t' read -r antigen query1 query2 <<< "$line"

    # Clean carriage returns
    antigen=$(echo "$antigen" | tr -d '\r\n')
    query1=$(echo "$query1" | tr -d '\r\n')
    query2=$(echo "$query2" | tr -d '\r\n')

    echo "🔍 Fetching: $antigen"

    # Primary query (UniProt/SwissProt)
    result=$(esearch -db protein -query "$query1" | efetch -format fasta 2>/dev/null)

    if [[ -n "$result" ]]; then
        echo "$result" > "$TMP_DIR/${antigen// /_}.fasta"
        echo " ✅ UniProt/SwissProt sequence found."
    else
        echo " ⚠️ UniProt/SwissProt not found, trying fallback..."
        result=$(esearch -db protein -query "$query2" | efetch -format fasta 2>/dev/null)
        if [[ -n "$result" ]]; then
            echo "$result" > "$TMP_DIR/${antigen// /_}.fasta"
            echo " ✅ Fallback sequence retrieved."
        else
            echo " ❌ No sequence found for: $antigen"
        fi
    fi

done

# === COMBINE RESULTS ===
echo "🧬 Concatenating individual FASTA files..."
cat "$TMP_DIR"/*.fasta > "$OUTPUT_FILE"
echo "✅ Final combined FASTA written to $OUTPUT_FILE"

# === CLEANUP ===
rm -rf "$TMP_DIR"
echo "🧹 Temporary files removed."
