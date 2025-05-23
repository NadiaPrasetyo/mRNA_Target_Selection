#!/bin/bash

set -e

# === USAGE CHECK ===
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 -i input1.xlsx input2.xlsx ... -o output.fasta"
    exit 1
fi

# === PARSE ARGUMENTS ===
INPUT_FILES=()
OUTPUT_FILE=""
PARSE_INPUTS=false

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i)
            PARSE_INPUTS=true
            ;;
        -o)
            PARSE_INPUTS=false
            shift
            OUTPUT_FILE="$1"
            ;;
        *)
            if $PARSE_INPUTS; then
                INPUT_FILES+=("$1")
            fi
            ;;
    esac
    shift
done

if [[ "${#INPUT_FILES[@]}" -lt 2 || -z "$OUTPUT_FILE" ]]; then
    echo "❌ You must provide at least two input Excel files and an output file."
    exit 1
fi

# === PREPARE TEMP WORKSPACE ===
TMP_DIR="data/tmp_run_$(date +%s)"
mkdir -p "$TMP_DIR"
ENTREZ_TSV="$TMP_DIR/entrez_queries.tsv"

# === ITERATIVELY REMOVE DUPLICATES ===
echo "🔄 Filtering duplicates from input files..."

# Copy the first input to start the chain
cp "${INPUT_FILES[0]}" "$TMP_DIR/base.xlsx"

for ((i = 1; i < ${#INPUT_FILES[@]}; i++)); do
    echo "🔎 Comparing: ${INPUT_FILES[$i]}"

    python3 bin/remove_duplicates_xlsx.py \
        --file_with_ids "$TMP_DIR/base.xlsx" \
        --file_new_antigens "${INPUT_FILES[$i]}" \
        --tsv_output "$ENTREZ_TSV"

    # Replace base for next iteration
    cp "${INPUT_FILES[$i]}" "$TMP_DIR/base.xlsx"
done

# === FETCH FASTA SEQUENCES ===
echo "🚀 Fetching FASTA sequences from NCBI..."

bin/fetch_entrez_antigen_fasta.sh "$ENTREZ_TSV" "$OUTPUT_FILE"

echo "✅ Completed: Output FASTA written to $OUTPUT_FILE"
echo "🗑️ Cleaning up temporary files..."
rm -rf "$TMP_DIR"
echo "🧹 Temporary files removed."
echo "✅ All done!"
