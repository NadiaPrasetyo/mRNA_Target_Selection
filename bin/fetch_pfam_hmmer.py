#!/usr/bin/env python3
"""
fetch_pfam_hmmer.py
Runner for fetching Pfam HMM profiles for proteins listed in an antigen CSV.
Overview:
    - Parses a Pfam-A.hmm file to create a mapping of Pfam accessions.
    - Reads a CSV file containing protein data and extracts Pfam IDs.
    - Fetches HMM profiles for the specified Pfam IDs using the HMMER tool.
    - Saves the fetched HMM profiles in a structured directory format for downstream analysis.
Arguments:
    pathogen_directory (str): Directory name under the `data/` folder.
    --pathogen_name (str): Prefix used in filenames (e.g., "staphylococcus_aureus").
    --input (str): Path to the input CSV file with antigen data (default: <pathogen_directory>/<pathogen_name>_compiled_proteins.csv).
    --output-dir (str): Directory where Pfam HMM profiles will be saved (default: <pathogen_directory>/pfam_hmms).
    --pfam_hmm (str): Path to the Pfam-A.hmm file (required).
Requirements:
    - Python packages: argparse, csv, os, subprocess, sys, re, tempfile, shutil.
    - HMMER must be installed and accessible in the system PATH (e.g., `conda install -c bioconda hmmer`).
    - A valid Pfam-A.hmm file must be provided.
Outputs:
    <output_dir>/                                      # Directory containing fetched HMM profiles.
    <output_dir>/<uniprot_id>_<pfam_id>.hmm            # Individual HMM files for each Pfam ID.
Notes:
    - The script ensures that the Pfam-A.hmm file is indexed before fetching profiles.
    - Logs warnings for missing or invalid Pfam IDs in the input CSV file.
    - Temporary files are used to store intermediate data and are automatically cleaned up.
    - Ensure that the HMMER tool is installed and accessible before running the script.
Author: Nadia
"""
import argparse
import csv
import os
import subprocess
import sys
import re
import tempfile
import shutil

csv.field_size_limit(sys.maxsize)

def check_hmmer_installed():
    """Ensure that hmmfetch from HMMER is available in PATH."""
    if shutil.which("hmmfetch") is None:
        print("❌ Error: hmmfetch not found. Please install HMMER (e.g., `conda install -c bioconda hmmer`).")
        sys.exit(1)


def load_pfam_accession_map(hmmfile):
    """
    Parse Pfam-A.hmm and create a mapping {PFxxxxx: PFxxxxx.yy}.
    Args:
        hmmfile (str): Path to the Pfam-A.hmm file.
    Returns:
        dict: Mapping of base Pfam accessions to full accessions.
    """
    accession_map = {}
    with open(hmmfile, "r") as f:
        for line in f:
            if line.startswith("ACC"):
                acc_full = line.split()[1].strip()   # PFxxxxx.yy
                acc_base = acc_full.split(".")[0]   # PFxxxxx
                accession_map[acc_base] = acc_full
    return accession_map


def main():
    """
    Main function to execute the pipeline for fetching Pfam HMM profiles.
    """
    parser = argparse.ArgumentParser(
        description="Fetch Pfam HMM profiles for proteins listed in an antigen CSV.",
        usage="python fetch_pfam_hmmer.py <pathogen_directory> [--pathogen_name <name>] or [--input <input_csv>] --pfam_hmm <Pfam-A.hmm>"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("--pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus_aureus")')
    parser.add_argument("--input", help="Input CSV file with antigen data (default: <pathogen_directory>/<pathogen_name>_compiled_proteins.csv)")
    parser.add_argument("--output-dir", help="Output directory for Pfam HMM profiles (default: <pathogen_directory>/pfam_hmms)")
    parser.add_argument("--pfam_hmm", required=True, help="Path to Pfam-A.hmm file")

    args = parser.parse_args()

    check_hmmer_installed()

    base_dir = os.path.join("data/", args.pathogen_directory)
    pathogen_name = args.pathogen_name.replace(" ", "_")

    input_csv = args.input or os.path.join(base_dir, f"{pathogen_name}_compiled_proteins.csv")
    output_dir = os.path.join(base_dir, args.output_dir if args.output_dir else "pfam_hmms")
    os.makedirs(output_dir, exist_ok=True)

    # Load accession mapping from Pfam-A.hmm
    accession_map = load_pfam_accession_map(args.pfam_hmm)

    # index the Pfam-A.hmm
    if not os.path.isfile(args.pfam_hmm):
        print(f"❌ Error: Pfam-A.hmm file not found at {args.pfam_hmm}. Please provide a valid file.", file=sys.stderr)
        sys.exit(1)

    index_file = f"{args.pfam_hmm}.ssi"
    if not os.path.isfile(index_file):
        try:
            subprocess.run(["hmmfetch", "--index", args.pfam_hmm], check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ Error: Failed to index Pfam-A.hmm file. Ensure the file is valid and not corrupted.\n{e}", file=sys.stderr)
            sys.exit(1)

    # Temporary keys file
    with tempfile.NamedTemporaryFile(mode="w", delete=True) as keys_tmp:
        pfam_keys = []

        # Collect Pfam accessions from CSV
        with open(input_csv, newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                uniprot_id = row["uniprot_accession"].strip()
                pfam_field = row["pfam"].strip()
                if not pfam_field:
                    continue

                # Pfam column may have multiple IDs separated by commas/semicolons
                pfam_ids = re.split(r"[;, ]+", pfam_field)
                for pfam_id in pfam_ids:
                    pfam_id = pfam_id.strip()
                    if not pfam_id:
                        continue

                    pfam_base = pfam_id.split(".")[0]  # normalize
                    if pfam_base in accession_map:
                        pfam_full = accession_map[pfam_base]
                        pfam_keys.append(pfam_full)
                        keys_tmp.write(pfam_full + "\n")
                        keys_tmp.flush()

                        # Fetch HMM for this accession
                        hmm_outfile = os.path.join(output_dir, f"{uniprot_id}_{pfam_base}.hmm")
                        try:
                            subprocess.run(
                                ["hmmfetch", args.pfam_hmm, pfam_full],
                                check=True,
                                stdout=open(hmm_outfile, "w"),
                            )
                        except subprocess.CalledProcessError:
                            print(f"⚠️ Warning: Failed to fetch {pfam_full} for {uniprot_id}", file=sys.stderr)
                    else:
                        print(f"⚠️ Warning: Pfam ID {pfam_id} not found in Pfam-A.hmm", file=sys.stderr)

        # keys_tmp will auto-delete on exit

    print(f"✅ Done. HMMs saved in {output_dir}")


if __name__ == "__main__":
    """Main entry point for the script."""
    main()
