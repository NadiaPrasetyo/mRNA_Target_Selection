#!/bin/bash

INPUT_CSV="data/S.aureus/entrez_queries.csv"
OUTPUT_FILE="data/S.aureus/fetched_antigens.fasta"

tail -n +2 "$INPUT_CSV" | awk -F',' '
    BEGIN { OFS=","; }
    {
        # Reconstruct fields safely, assuming 3 total fields
        line = $0
        gsub(/\"/, "", line)
        split(line, fields, /,/)
        antigen = fields[1]
        query1 = fields[2]
        query2 = fields[3]

        # Print a marker to separate logic blocks
        print "===FETCH===" antigen "---" query1 "---" query2
    }
' | while IFS='---' read -r marker antigen query1 query2; do
    if [[ $marker != "===FETCH===" ]]; then continue; fi

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
