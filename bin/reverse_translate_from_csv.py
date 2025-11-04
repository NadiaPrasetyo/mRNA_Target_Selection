#!/usr/bin/env python3
import csv
import subprocess
import tempfile
import os
import sys
import logging
from pathlib import Path
import argparse


def setup_logging(args):
    """
    Set up logging depending on --verbose flag.
    """
    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    output_dir = pathogen_path / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        log_file = output_dir / "reverse_translate.log"
        handlers = [
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode="a")
        ]
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=handlers,
            force=True
        )
        logging.info(f"Verbose mode on. Logging to {log_file}")
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s: %(message)s",
            force=True
        )


def run_backtranseq(protein_fasta, output_fasta, emboss_exe, codon_file=None):
    """
    Run EMBOSS backtranseq command on a given protein FASTA.
    """
    cmd = [emboss_exe, "-sequence", protein_fasta, "-outfile", output_fasta, "-auto"]
    if codon_file:
        cmd += ["-cfile", codon_file]

    logging.debug(f"Running backtranseq: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"backtranseq failed for {protein_fasta}: {e.stderr.decode().strip()}")
        raise


def fasta_from_row(row):
    """
    Convert a CSV row to a FASTA formatted string.
    """
    header_parts = [row["uniprot_accession"]]
    if row.get("protein_name"):
        header_parts.append(row["protein_name"])
    if row.get("organism_name"):
        header_parts.append(f"[{row['organism_name']}]")

    header = " ".join(header_parts)
    sequence = row["sequence"].strip().replace(" ", "")
    return f">{header}\n{sequence}\n"


def process_csv(input_csv, output_csv, emboss_exe, codon_file=None):
    """
    Read input CSV, run backtranseq for each protein, and save the updated CSV.
    """
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)

    logging.info(f"Reading input CSV: {input_csv}")
    logging.info(f"Writing output CSV: {output_csv}")

    output_rows = []

    with open(input_csv, newline='') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        if "nucleotide_sequence" not in fieldnames:
            fieldnames.append("nucleotide_sequence")

        for row in reader:
            acc = row.get("uniprot_accession", "UNKNOWN")
            seq = row.get("sequence", "").strip()

            if not seq:
                logging.warning(f"{acc}: Missing protein sequence, skipping.")
                output_rows.append(row)
                continue

            logging.debug(f"Processing {acc}...")

            with tempfile.TemporaryDirectory() as tmpdir:
                fasta_path = Path(tmpdir) / f"{acc}.fasta"
                output_fasta = Path(tmpdir) / f"{acc}_nt.fasta"

                with open(fasta_path, "w") as f:
                    f.write(fasta_from_row(row))

                try:
                    run_backtranseq(fasta_path, output_fasta, emboss_exe, codon_file)
                except Exception:
                    row["nucleotide_sequence"] = ""
                    output_rows.append(row)
                    continue

                with open(output_fasta) as f:
                    nt_seq = "".join([l.strip() for l in f if not l.startswith(">")])
                    row["nucleotide_sequence"] = nt_seq

            output_rows.append(row)
            logging.debug(f"{acc}: Reverse translation complete.")

    with open(output_csv, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    logging.info(f"✅ Processing complete. Output written to {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Reverse translate protein sequences using EMBOSS backtranseq and update CSV."
    )
    parser.add_argument("--input_csv", required=True, help="Path to input CSV file.")
    parser.add_argument("--output_csv", required=True, help="Path to output CSV file.")
    parser.add_argument("--pathogen_dir", required=True, help="Pathogen subdirectory name under data/")
    parser.add_argument("--sequence_dir", required=False, help="Sequence directory (not used directly).")
    parser.add_argument("--output_dir", required=True, help="Subdirectory for output/log files.")
    parser.add_argument("--tool_root", required=True, help="Path to EMBOSS tool root directory.")
    parser.add_argument("--codon_file", required=False, help="Path to codon usage file (e.g., Ecoli.cut).")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to file and console.")
    args = parser.parse_args()

    setup_logging(args)

    emboss_exe = str(Path(args.tool_root) / "backtranseq")
    logging.debug(f"Using EMBOSS binary: {emboss_exe}")

    process_csv(args.input_csv, args.output_csv, emboss_exe, args.codon_file)


if __name__ == "__main__":
    main()
