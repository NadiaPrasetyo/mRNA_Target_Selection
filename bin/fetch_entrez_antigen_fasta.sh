#!/bin/bash

input="data/S.aureus/entrez_queries.txt"
output_dir="data/S.aureus/fasta_results"
mkdir -p "$output_dir"

while IFS="|||" read -r antigen main_query fallback_query; do
    echo "Fetching: $antigen"

    # First try SwissProt
    result=$(esearch -db protein -query "$main_query" | efetch -format fasta)
    if [ -z "$result" ]; then
        echo " - SwissProt not found, trying general S. aureus proteins..."
        result=$(esearch -db protein -query "$fallback_query" | efetch -format fasta)
    fi

    if [ -n "$result" ]; then
        echo "$result" > "${output_dir}/${antigen// /_}.fasta"
        echo " - Success"
    else
        echo " - No sequence found"
    fi
done < "$input"
