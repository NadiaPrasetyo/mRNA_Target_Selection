#!/bin/bash

set -e

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pathogen_directory> <file_name.csv>"
    exit 1
fi

# Arguments
PATHOGEN_DIR="$1"
CSV_FILE="$2"
CSV_PATH="data/${PATHOGEN_DIR}/${CSV_FILE}"
OUTPUT_DIR="data/${PATHOGEN_DIR}/strain_genomes"

# Temp files
EMBL_LIST_FILE=$(mktemp)
FASTA_BULK_FILE=$(mktemp)

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Create a mapping of EMBL_ID to Strain
declare -A strain_map

# Read CSV and build EMBL list and map
tail -n +2 "$CSV_PATH" | while IFS=, read -r Strain EMBL_ID CC LocusTag
do
    if [ -n "$EMBL_ID" ]; then
        echo "$EMBL_ID" >> "$EMBL_LIST_FILE"
        strain_map["$EMBL_ID"]="$Strain"
    else
        echo "Skipping $Strain: missing EMBL_ID"
    fi
done

# Fetch all sequences in bulk
echo "Fetching all EMBL IDs in bulk..."
efetch -db nucleotide -format fasta -id "$(paste -sd, "$EMBL_LIST_FILE")" > "$FASTA_BULK_FILE"

# Split bulk FASTA into individual files by header match
echo "Splitting bulk FASTA..."
awk -v outdir="$OUTPUT_DIR" -v csv="$CSV_PATH" '
    BEGIN {
        FS=",";
        while ((getline < csv) > 0) {
            if (NR == 1) continue; # skip header
            strain[$2] = $1;
        }
    }
    /^>/ {
        if (seq != "") {
            print seq > file;
            close(file);
        }
        id = substr($1, 2);  # remove ">"
        gsub(/ .*/, "", id); # keep accession only
        name = strain[id];
        if (name == "") name = id; # fallback
        file = outdir "/" name ".fasta";
        print $0 > file;
        seq = "";
        next;
    }
    {
        seq = seq $0 "\n";
    }
    END {
        if (seq != "") {
            print seq > file;
            close(file);
        }
    }
' "$FASTA_BULK_FILE"

# Cleanup temp files
rm -f "$EMBL_LIST_FILE" "$FASTA_BULK_FILE"

echo "Done. Sequences saved to $OUTPUT_DIR"
