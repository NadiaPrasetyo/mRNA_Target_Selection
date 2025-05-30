: '
/**
 * fetch_NCBI_strain_genome.sh
 *
 * Purpose:
 *   This script automates the retrieval of nucleotide FASTA sequences from NCBI for a list of EMBL IDs
 *   specified in a CSV file. It maps each EMBL ID to a strain name, downloads all sequences in bulk,
 *   and splits them into individual FASTA files named after their respective strains. Unmatched EMBL IDs
 *   are logged for review.
 *
 * Usage:
 *   ./fetch_NCBI_strain_genome.sh <pathogen_directory> <file_name.csv>
 *
 *   - <pathogen_directory>: Subdirectory under "data/" where the CSV and output will reside.
 *   - <file_name.csv>: CSV file (in the specified directory) with strain names and EMBL IDs (columns: strain,embl_id).
 *
 * Example:
 *   ./fetch_NCBI_strain_genome.sh ecoli ecoli_strains.csv
 *
 * Output:
 *   - Individual FASTA files for each strain in data/<pathogen_directory>/strain_genomes/
 *   - Log file of unmatched EMBL IDs (unmatched_ids.log) in the same output directory.
 *
 * Author: Nadia
 */
'
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

mkdir -p "$OUTPUT_DIR"
> "$LOG_FILE"   # clear log file

# Create temp files for EMBL ID list and bulk FASTA
EMBL_LIST_FILE=$(mktemp)
FASTA_BULK_FILE=$(mktemp)

# Ensure temp files are deleted on exit
trap 'rm -f "$EMBL_LIST_FILE" "$FASTA_BULK_FILE"' EXIT

echo "Extracting EMBL IDs..."
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
            gsub(/\r/, "", id);
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
        sub(/ .*/, "", id);
        base_id = id;
        sub(/\..*/, "", base_id);

        strain = (id in strain_name) ? strain_name[id] : ((base_id in strain_name) ? strain_name[base_id] : "");

        if (strain != "") {
            gsub(/ /, "_", strain);
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

if [ ! -s "$LOG_FILE" ]; then
    rm -f "$LOG_FILE"
fi

echo "Done."
echo "FASTA files written to: $OUTPUT_DIR"
echo "Unmatched EMBL IDs logged in: $LOG_FILE"
