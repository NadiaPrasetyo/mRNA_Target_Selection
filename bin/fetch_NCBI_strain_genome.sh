###############################################################################
# fetch_NCBI_strain_genome.sh
#
# Updated to include --random and --random-num options for fetching random
# complete genomes for a specified pathogen while maintaining the original usage.
###############################################################################
#!/bin/bash
set -e

# Default thread count and random genome count
NUM_THREADS=4
RANDOM_NUM=5

# Function to display help
print_help() {
    echo "Usage: $0 [OPTIONS] <pathogen_directory> <strain_file_name.csv>"
    echo ""
    echo "Options:"
    echo "  --threads N          Number of threads to use for parallel processing (default: 4)"
    echo "  --random \"<pathogen_name>\""
    echo "                       Fetch random complete genomes for the specified pathogen."
    echo "  --random-num N       Number of random genomes to fetch (default: 5)"
    echo "  --help               Display this help message and exit."
    echo ""
    echo "Usages:"
    echo "  $0 --random \"<pathogen_name>\" [--random-num N] <pathogen_directory>"
    echo "      Fetch random complete genomes for a pathogen."
    echo ""
    echo "  $0 [--threads N] <pathogen_directory> <strain_file_name.csv>"
    echo "      Fetch genomes based on a CSV file containing strain names and EMBL IDs."
    echo ""
    exit 0
}

# Argument parsing
POSITIONAL=()
RANDOM_MODE=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            NUM_THREADS="$2"
            shift 2
            ;;
        --random)
            RANDOM_MODE=true
            PATHOGEN_NAME="$2"
            shift 2
            ;;
        --random-num)
            RANDOM_NUM="$2"
            shift 2
            ;;
        --help | -h)
            print_help
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

if [ "$RANDOM_MODE" = true ]; then
    if [ -z "$PATHOGEN_NAME" ]; then
        echo "Error: --random requires a pathogen name in quotes."
        exit 1
    fi

    if [ "$#" -ne 1 ]; then
        echo "Usage: $0 --random \"<pathogen_name>\" [--random-num N] <pathogen_directory>"
        exit 1
    fi

    PATHOGEN_DIR="$1"
    OUTPUT_DIR="data/${PATHOGEN_DIR}/strain_genomes"
    LOG_FILE="${OUTPUT_DIR}/unmatched_ids.log"

    mkdir -p "$OUTPUT_DIR"
    > "$LOG_FILE"   # clear log file

    echo "Fetching $RANDOM_NUM random complete genomes for pathogen: $PATHOGEN_NAME"

    # get the reference genome for the pathogen
    REFERENCE_GENOME=$(esearch -db nucleotide -query "\"${PATHOGEN_NAME}\"[Organism] AND \"complete genome\"[All Fields]" | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary -element Caption | \
        head -n 1) # Get the first result as the reference genome

    if [ -z "$REFERENCE_GENOME" ]; then
        echo "Error: No reference genome found for pathogen: $PATHOGEN_NAME"
        exit 1
    fi

    # Extract the sequence length from the reference genome
    REFERENCE_LENGTH=$(esearch -db nucleotide -query "\"${REFERENCE_GENOME}\"[Caption]" | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary -element Slen)

    # Fetch random genome IDs that have sequence length within the range of 500,000 bp above or under the reference genome of the pathogen
    RANDOM_GENOME_IDS=$(esearch -db nucleotide -query "\"${PATHOGEN_NAME}\"[Organism] AND \"complete genome\"[All Fields] AND (\"$((REFERENCE_LENGTH - 500000))\"[SLEN] : \"$((REFERENCE_LENGTH + 500000))\"[SLEN])" | \
        efetch -format docsum | \
        xtract -pattern DocumentSummary -element Caption | \
        shuf -n "$RANDOM_NUM") # Get random IDs by shuffling the results and selecting the top N

    if [ -z "$RANDOM_GENOME_IDS" ]; then
        echo "Error: No genomes found for pathogen: $PATHOGEN_NAME"
        exit 1
    fi

    # Fetch genomes in bulk
    FASTA_BULK_FILE=$(mktemp)
    trap 'rm -f "$FASTA_BULK_FILE"' EXIT

    echo "Fetching genomes..."
    efetch -db nucleotide -format fasta -id "$RANDOM_GENOME_IDS" > "$FASTA_BULK_FILE"

    echo "Splitting bulk FASTA..."
    awk -v outdir="$OUTPUT_DIR" '
        /^>/ {
            if (seq != "") {
                print seq > file;
                close(file);
            }
            seq = "";
            id = substr($1, 2);
            gsub(/ /, "_", id);
            file = outdir "/" id ".fasta";
            print $0 > file;
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

    echo "Genome FASTA files saved in: $OUTPUT_DIR"
else
    if [ "$#" -ne 2 ]; then
        echo "Usage: $0 [--threads N] <pathogen_directory> <strain_file_name.csv>"
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

# Function to remove asterisks from translated FASTA files
remove_asterisks_from_translations() {
    echo "Removing asterisks from translated FASTA files..."
    for file in "$OUTPUT_DIR"/*_translated.fasta; do
        if [ -f "$file" ]; then
            sed -i 's/\*//g' "$file"
        fi
    done
    echo "Asterisks removed from translated files."
}

# Call the function after translation
remove_asterisks_from_translations

echo "Done."
