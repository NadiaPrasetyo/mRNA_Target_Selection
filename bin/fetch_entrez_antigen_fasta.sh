#!/bin/bash

input="data/S.aureus/entrez_queries.txt"
output_dir="data/S.aureus/fasta_results"
mkdir -p "$output_dir"

# Clear the output file if it exists
> "$OUTPUT_FILE"

while IFS='|||' read -r antigen query swissprot_query fallback_query; do
    echo "Fetching: $antigen"

    # Try SwissProt query
    result=$(esearch -db protein -query "$swissprot_query" | \
             efetch -format fasta 2>/dev/null)

    if [[ -n "$result" ]]; then
        echo "$result" >> "$OUTPUT_FILE"
        echo " - SwissProt sequence retrieved."
    else
        echo " - SwissProt not found, trying general S. aureus proteins..."
        result=$(esearch -db protein -query "$fallback_query" | \
                 efetch -format fasta 2>/dev/null)

        if [[ -n "$result" ]]; then
            echo "$result" >> "$OUTPUT_FILE"
            echo " - Fallback sequence retrieved."
        else
            echo " WARNING:  FAILURE ( $(date) )"
            echo " - No sequence found for: $antigen"
        fi
    fi

    # Slight delay for safety
    sleep 0.3
done < "$INPUT_FILE"
