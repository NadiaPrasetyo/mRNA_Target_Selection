: '
/**
 * fetch_NCBI_strain_genome.sh
 *
 * Purpose:
 *   This script automates the retrieval of nucleotide FASTA sequences from NCBI for a list of EMBL IDs
 *   specified in a CSV file. It maps each EMBL ID to a strain name, downloads all sequences in bulk,
 *   and splits them into individual FASTA files named after their respective strains.
 *
 *   After nucleotide sequences are fetched, the script uses `seqkit translate` to convert them into
 *   amino acid sequences. Translation is performed in parallel using multiple threads.
 *
 * Features:
 *   - Avoids re-downloading genome sequences if output files already exist.
 *   - Avoids re-translating sequences if translated files already exist.
 *   - Allows configurable parallel translation via --threads option.
 *   - Verifies that `seqkit` is installed before proceeding.
 *
 * Usage:
 *   ./fetch_NCBI_strain_genome.sh [--threads N] <pathogen_directory> <file_name.csv>
 *
 *   Options:
 *     --threads N           Number of threads to use for parallel translation (default: 4)
 *     <pathogen_directory>  Subdirectory under "data/" where the CSV and output will reside
 *     <file_name.csv>       CSV file (in the specified directory) with strain names and EMBL IDs
 *                           (columns: strain,embl_id)
 *
 * Example:
 *   ./fetch_NCBI_strain_genome.sh --threads 6 ecoli ecoli_strains.csv
 *
 * Output:
 *   - Individual FASTA files for each strain in:
 *       data/<pathogen_directory>/strain_genomes/
 *   - Translated amino acid FASTA files (one per strain), with suffix "_translated.fasta"
 *   - Log file of unmatched EMBL IDs: unmatched_ids.log (in the same output directory)
 *
 * Author: Nadia
 */
'
#!/bin/bash
set -e

# Default thread count
NUM_THREADS=4

# Argument parsing
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            NUM_THREADS="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# Restore positional args
set -- "${POSITIONAL[@]}"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [--threads N] <pathogen_directory> <strain_file_name.csv>"
    exit 1
fi

# Check if Entrez Direct (efetch) is installed
if ! command -v efetch &> /dev/null; then
    echo "Error: 'efetch' (Entrez Direct) is not installed. Please install it before running this script."
    echo "Visit https://www.ncbi.nlm.nih.gov/books/NBK179288/ for installation instructions."
    exit 1
fi

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed. Please install it before running this script."
    echo "Visit https://bioinf.shenwei.me/seqkit/ for installation instructions."
    exit 1
fi

# Check if mmseqs2 is installed
if ! command -v mmseqs &> /dev/null; then
    echo "Error: 'mmseqs' is not installed. Please install it before running this script."
    echo "Visit https://github.com/soedinglab/MMseqs2 for installation instructions."
    exit 1
fi


PATHOGEN_DIR="$1"
CSV_FILE="$2"
CSV_PATH="data/${PATHOGEN_DIR}/${CSV_FILE}"
OUTPUT_DIR="data/${PATHOGEN_DIR}/strain_genomes"
LOG_FILE="${OUTPUT_DIR}/unmatched_ids.log"

mkdir -p "$OUTPUT_DIR"
> "$LOG_FILE"   # clear log file

# Check if FASTA files already exist
FASTA_COUNT=$(find "$OUTPUT_DIR" -maxdepth 1 -name "*.fasta" ! -name "*_translated.fasta" | wc -l)

if [ "$FASTA_COUNT" -eq 0 ]; then
    echo "No genome FASTA files found. Proceeding to fetch..."

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
            while ((getline < csv) > 0) {
                if (NR == 1) continue;
                strain = $1;
                id = $2;
                gsub(/\r/, "", id);
                base_id = id;
                sub(/\..*/, "", base_id);
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

    echo "Genome FASTA files saved in: $OUTPUT_DIR"
else
    echo "Genome FASTA files already exist in $OUTPUT_DIR. Skipping fetch."
fi

echo "Translating nucleotide sequences to amino acid sequences in parallel using $NUM_THREADS threads..."

export OUTPUT_DIR  # make available to subshells
export -f

find "$OUTPUT_DIR" -maxdepth 1 -name "*.fasta" ! -name "*_translated.fasta" | \
xargs -P "$NUM_THREADS" -I {} bash -c '
    fasta="{}"
    translated="${fasta%.fasta}_translated.fasta"
    if [ ! -f "$translated" ]; then
        echo "Translating: $(basename "$fasta")"
        seqkit translate -f 1 -M "$fasta" > "$translated"
    else
        echo "Translation already exists for: $(basename "$fasta")"
    fi
'

echo "Translation complete."

echo "Done."
