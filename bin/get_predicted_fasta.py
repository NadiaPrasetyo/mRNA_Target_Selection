import os
import argparse
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences based on accession numbers from a CSV.")
    parser.add_argument("fasta_folder", help="Folder containing *_matched_antigens.fasta files")
    parser.add_argument("csv_file", help="CSV file containing accession numbers")
    parser.add_argument("-o", "--output", help="Output FASTA file path (default: parent of FASTA folder/predicted_sequences.fasta)")
    args = parser.parse_args()

    fasta_folder = args.fasta_folder
    csv_file = args.csv_file
    output_fasta = args.output

    if output_fasta is None:
        parent_folder = os.path.dirname(os.path.abspath(fasta_folder))
        output_fasta = os.path.join(parent_folder, "predicted_sequences.fasta")

    # Read accessions from CSV
    df = pd.read_csv(csv_file)
    accessions_set = set(df['accession'].astype(str))

    # Filter sequences
    filtered_sequences = []

    for filename in os.listdir(fasta_folder):
        if filename.endswith("_matched_antigens.fasta"):
            fasta_path = os.path.join(fasta_folder, filename)
            for record in SeqIO.parse(fasta_path, "fasta"):
                # Extract accession from header: >antigen_96|A0A391AB27||WP_000669040.1|tpos:1-208
                header_parts = record.id.split('|')
                if len(header_parts) > 1:
                    accession = header_parts[1]
                    if accession in accessions_set:
                        filtered_sequences.append(record)

    if filtered_sequences:
        SeqIO.write(filtered_sequences, output_fasta, "fasta")
        print(f"Saved {len(filtered_sequences)} sequences to {output_fasta}")
    else:
        print("No sequences matched the accessions in the CSV.")

if __name__ == "__main__":
    main()

