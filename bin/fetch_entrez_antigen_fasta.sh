#!/bin/bash

INPUT_CSV="data/S.aureus/entrez_queries.tsv"
OUTPUT_FILE="data/S.aureus/fetched_antigens.fasta"

# Ensure output file is fresh
> "$OUTPUT_FILE"

# Skip header, process TSV
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r antigen query1 query2; do
    echo "🔍 Fetching: $antigen"
    
    # Primary query (SwissProt)
    result=$(esearch -db protein -query "$query1" | efetch -format fasta 2>/dev/null)
    
    if [[ -n "$result" ]]; then
        echo "$result" >> "$OUTPUT_FILE"
        echo " ✅ SwissProt sequence found."
    else
        echo " ⚠️ SwissProt not found, trying fallback..."
        result=$(esearch -db protein -query "$query2" | efetch -format fasta 2>/dev/null)
        if [[ -n "$result" ]]; then
            echo "$result" >> "$OUTPUT_FILE"
            echo " ✅ Fallback sequence retrieved."
        else
            echo " ❌ No sequence found for: $antigen"
        fi
    fi

    # Respect EDirect rate limits
    sleep 0.3
done