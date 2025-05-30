#!/bin/bash

set -e

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pathogen_directory> <file_name.csv>"
    exit 1
fi

PATHOGEN_DIR="$1"
CSV_FILE="$2"
CSV_PATH="data/${PATHOGEN_DIR}/${CSV_FILE}"
OUTPUT_DIR="data/${PATHOGEN_DIR}/strain_genomes"

mkdir -p "$OUTPUT_DIR"

# Max number of parallel jobs
MAX_JOBS=3

# Function to fetch a strain's sequence
fetch_sequence() {
    local strain="$1"
    local embl_id="$2"
    local output_file="$3"

    if [ -n "$embl_id" ]; then
        echo "Fetching $strain ($embl_id)..."
        esearch -db nucleotide -query "${embl_id}[Accession]" | efetch -format fasta > "$output_file" 2>/dev/null
        if [ $? -ne 0 ]; then
            echo "Error fetching $strain ($embl_id)"
        fi
    else
        echo "Skipping $strain: missing EMBL_ID"
    fi
}

# Limit parallel jobs
parallel_jobs() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        wait -n
    done
}

# Read and process each line (skip header)
tail -n +2 "$CSV_PATH" | while IFS=, read -r Strain EMBL_ID CC LocusTag
do
    OUTPUT_FILE="${OUTPUT_DIR}/${Strain}.fasta"
    parallel_jobs
    fetch_sequence "$Strain" "$EMBL_ID" "$OUTPUT_FILE" &
done

# Wait for all background jobs to finish
wait

echo "All FASTA sequences saved to $OUTPUT_DIR"
