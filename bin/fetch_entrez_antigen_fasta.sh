#!/bin/bash

INPUT_CSV="data/S.aureus/entrez_queries.csv"
OUTPUT_FILE="data/S.aureus/fetched_antigens.fasta"

# Skip header line and loop through each row
tail -n +2 "$INPUT_CSV" | while IFS=',' read -r antigen query1 query2; do
    # Remove any stray carriage returns or quotes
    antigen=$(echo "$antigen" | tr -d '\r"')
    query1=$(echo "$query1" | tr -d '\r"')
    query2=$(echo "$query2" | tr -d '\r"')

    echo "🔍 Fetching: $antigen"
    
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
    sleep 0.3
done
