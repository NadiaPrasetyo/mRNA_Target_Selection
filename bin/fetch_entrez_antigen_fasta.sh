#!/bin/bash

INPUT_TSV="data/S.aureus/entrez_queries.tsv"
OUTPUT_FILE="data/S.aureus/fetched_antigens.fasta"
TMP_DIR="data/S.aureus/tmp_fasta"

# Prepare output and temp folder
> "$OUTPUT_FILE"
mkdir -p "$TMP_DIR"

# Read all lines except header into an array
mapfile -t lines < <(tail -n +2 "$INPUT_TSV")

for line in "${lines[@]}"; do
    # Extract fields by splitting on tabs
    IFS=$'\t' read -r antigen query1 query2 <<< "$line"

    # Strip carriage returns/newlines from variables (fix Windows line endings)
    antigen=$(echo "$antigen" | tr -d '\r\n')
    query1=$(echo "$query1" | tr -d '\r\n')
    query2=$(echo "$query2" | tr -d '\r\n')

    echo "🔍 Fetching: $antigen"

    # Primary query (SwissProt)
    result=$(esearch -db protein -query "$query1" | efetch -format fasta 2>/dev/null)

    if [[ -n "$result" ]]; then
        echo "$result" > "$TMP_DIR/${antigen// /_}.fasta"
        echo " ✅ SwissProt sequence found."
    else
        echo " ⚠️ SwissProt not found, trying fallback..."
        result=$(esearch -db protein -query "$query2" | efetch -format fasta 2>/dev/null)
        if [[ -n "$result" ]]; then
            echo "$result" > "$TMP_DIR/${antigen// /_}.fasta"
            echo " ✅ Fallback sequence retrieved."
        else
            echo " ❌ No sequence found for: $antigen"
        fi
    fi

    sleep 0.3
done

# Concatenate all individual fasta files into the final output
cat "$TMP_DIR"/*.fasta > "$OUTPUT_FILE"

echo "🧬 Concatenating individual FASTA files..."
echo "✅ Final combined FASTA written to $OUTPUT_FILE"

