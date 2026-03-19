import argparse
import os
import pandas as pd
from Bio import SeqIO

def get_proteome_fastas(proteome_dir):
    proteome_fasta_files = [
        f for f in os.listdir(proteome_dir) if f.endswith(".fasta")
    ]
    return proteome_fasta_files

def get_filtered_accessions(input_dir):
    # accession column of files with _filtered_features.csv
    filtered_accessions = []
    for file in os.listdir(input_dir):
        if "_filtered_features.csv" in file:
            df = pd.read_csv(os.path.join(input_dir, file), usecols=["accession"])
            filtered_accessions.extend(df["accession"].tolist())
    return filtered_accessions

def filter_proteome(proteome_fasta_files, filtered_accessions, output_dir, proteome_dir):
    for fasta_file in proteome_fasta_files:
        filtered_fasta_file = os.path.join(output_dir, fasta_file.replace(".fasta", "_filtered.fasta"))
        with open(os.path.join(proteome_dir, fasta_file)) as fin:  
            with open(filtered_fasta_file, "w") as fout:
                for record in SeqIO.parse(fin, "fasta"):
                    for filtered_accession in filtered_accessions:
                        if filtered_accession in record.id:
                            SeqIO.write(record, fout, "fasta")

def summarize_to_csv(output_dir):
    # combine all the filtered fasta files into a single csv with fields: accession, sequence
    filtered_fasta_files = [
        f for f in os.listdir(output_dir) if f.endswith("_filtered.fasta")
    ]
    combined_csv = os.path.join(output_dir, "combined_filtered_proteome.csv")
    with open(combined_csv, "w") as fout:
        fout.write("accession,sequence\n")
        for fasta_file in filtered_fasta_files:
            with open(os.path.join(output_dir, fasta_file)) as fin:
                for record in SeqIO.parse(fin, "fasta"):
                    fout.write(f"{record.id},{record.seq}\n")
def main(input_dir, proteome_dir, output_dir):
    proteome_fasta_files = get_proteome_fastas(proteome_dir)
    filtered_accessions = get_filtered_accessions(input_dir)
    filter_proteome(proteome_fasta_files, filtered_accessions, output_dir, proteome_dir)
    summarize_to_csv(output_dir)

    print(f"✅ Finished! Saved {len(filtered_accessions)} accessions to {output_dir}")
    print(f"Next step: get the uniparc and uniprot accessions using match_uniparc_uniprotkb.py")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter proteome data based on subcellular location data (output of combine_filter_subcellular_loc.py)."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Path to directory containing prediction output files."
    )
    parser.add_argument(
        "-p", "--proteome_dir",
        required=True,
        help="Path to directory containing proteome fasta files."
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=False,
        default="data/",
        help="Optional directory to save per-stem FASTA and combined CSV outputs."
    )
    args = parser.parse_args()

    main(args.input_dir, args.proteome_dir, args.output_dir)

    