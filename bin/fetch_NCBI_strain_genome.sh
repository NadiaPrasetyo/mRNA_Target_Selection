#!/bin/bash
set -e

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pathogen_directory> <file_name.csv>"
    exit 1
fi

PATHOGEN_DIR="$1"
CSV_FILE="$2"
CSV_PATH="data/${PATHOGEN_DIR}/${CSV_FILE}"
OUTPUT_DIR="data/${PATHOGEN_DIR}/strain_genomes"
LOG_FILE="${OUTPUT_DIR}/unmatched_ids.log"

# Keep temp files persistent for inspection
EMBL_LIST_FILE="data/${PATHOGEN_DIR}/embl_ids.txt"
FASTA_BULK_FILE="data/${PATHOGEN_DIR}/fasta_bulk.fa"

mkdir -p "$OUTPUT_DIR"
> "$LOG_FILE"   # clear log file

# Extract EMBL IDs (strip header) and save to file without "[Accession]"
tail -n +2 "$CSV_PATH" | cut -d',' -f2 | sed '/^$/d' > "$EMBL_LIST_FILE"

echo "Fetching all EMBL IDs in bulk..."
efetch -db nucleotide -format fasta -id "$(paste -sd, "$EMBL_LIST_FILE")" > "$FASTA_BULK_FILE"

echo "Splitting bulk FASTA..."

awk -v outdir="$OUTPUT_DIR" -v csv="$CSV_PATH" -v logfile="$LOG_FILE" '
    BEGIN {
        FS = ",";
        # Build EMBL_ID -> strain mapping from CSV
        while ((getline < csv) > 0) {
            if (NR == 1) continue;
            strain = $1;
            id = $2;
            gsub(/\r/, "", id);  # remove Windows CR if present
            base_id = id;
            sub(/\..*/, "", base_id);
            # Map both versioned and base IDs
            strain_name[id] = strain;
            strain_name[base_id] = strain;
        }
    }
    /^>/ {
        if (seq != "" && matched) {
            print seq > file;
            close(file);
        }
        matched = 0;
        id = substr($1, 2);
        sub(/ .*/, "", id);     # cut at first space
        base_id = id;
        sub(/\..*/, "", base_id);

        strain = (id in strain_name) ? strain_name[id] : ((base_id in strain_name) ? strain_name[base_id] : "");

        if (strain != "") {
            gsub(/ /, "_", strain);     # replace spaces with underscores
            file = outdir "/" strain ".fasta";
            print $0 > file;
            matched = 1;
        } else {
            print "Unmatched EMBL ID: " id >> logfile;
        }
        seq = "";
        next;
    }
    {
        if (matched) seq = seq $0 "\n";
    }
    END {
        if (seq != "" && matched) {
            print seq > file;
            close(file);
        }
    }
' "$FASTA_BULK_FILE"

echo "Done."
echo "FASTA files are in: $OUTPUT_DIR"
echo "EMBL IDs list kept in: $EMBL_LIST_FILE"
echo "Bulk FASTA file kept in: $FASTA_BULK_FILE"
echo "Unmatched IDs logged in: $LOG_FILE"
