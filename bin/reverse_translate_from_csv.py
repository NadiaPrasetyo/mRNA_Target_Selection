#!/usr/bin/env python3
import csv
import subprocess
import tempfile
import sys
import logging
from pathlib import Path
import argparse
import shutil


def setup_logging(pathogen_dir: Path, verbose: bool):
    """
    Set up logging configuration.
    """
    pathogen_dir.mkdir(parents=True, exist_ok=True)
    log_file = pathogen_dir / "reverse_translate.log"

    if verbose:
        handlers = [
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode="a")
        ]
        level = logging.DEBUG
    else:
        handlers = [logging.StreamHandler(sys.stdout)]
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=handlers,
        force=True
    )

    if verbose:
        logging.info(f"Verbose logging enabled. Log file: {log_file}")
    else:
        logging.info("Logging set to INFO level.")


def check_backtranseq():
    """
    Ensure the EMBOSS backtranseq command is available in PATH.
    """
    path = shutil.which("backtranseq")
    if path is None:
        logging.error("❌ EMBOSS 'backtranseq' not found in PATH. Please install EMBOSS or update your PATH.")
        sys.exit(1)
    logging.info(f"Found EMBOSS backtranseq at: {path}")
    return path


def run_backtranseq(input_fasta: Path, output_fasta: Path):
    """
    Run EMBOSS backtranseq on a single FASTA file (handles multiple sequences).
    """
    cmd = ["backtranseq", "-sequence", str(input_fasta), "-outfile", str(output_fasta), "-auto"]
    logging.debug(f"Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"backtranseq failed: {e.stderr.decode().strip()}")
        raise


def parse_fasta(filepath: Path):
    """
    Parse a FASTA file into a dict of {header: sequence}.
    Header = first token on '>' line.
    """
    sequences = {}
    header = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]  # UniProt accession assumed first
                sequences[header] = ""
            else:
                if header:
                    sequences[header] += line
    return sequences


def process_csv(pathogen_dir: Path, input_csv: Path, output_csv: Path):
    """
    Read input CSV, run backtranseq once for all proteins, and save the updated CSV.
    """
    logging.info(f"Reading input CSV: {input_csv}")
    logging.info(f"Writing output CSV: {output_csv}")

    with open(input_csv, newline='') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames or []
        if "nucleotide_sequence" not in fieldnames:
            fieldnames.append("nucleotide_sequence")

        rows = list(reader)

    # Create combined FASTA file for all valid sequences
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_input = Path(tmpdir) / "all_proteins.fasta"
        fasta_output = Path(tmpdir) / "all_proteins_nt.fasta"

        valid_count = 0
        with open(fasta_input, "w") as fasta:
            for row in rows:
                acc = row.get("uniprot_accession", "UNKNOWN")
                seq = (row.get("sequence") or "").strip()
                if not seq:
                    logging.warning(f"{acc}: Missing protein sequence, skipping.")
                    continue
                fasta.write(f">{acc}\n{seq}\n")
                valid_count += 1

        if valid_count == 0:
            logging.error("No valid protein sequences found in the input CSV.")
            sys.exit(1)

        logging.info(f"Wrote {valid_count} protein sequences to {fasta_input}")

        # Run backtranseq once
        run_backtranseq(fasta_input, fasta_output)

        # Parse output FASTA (nucleotide sequences)
        nt_sequences = parse_fasta(fasta_output)
        logging.info(f"Parsed {len(nt_sequences)} nucleotide sequences from EMBOSS output.")

    # Merge results back into CSV
    for row in rows:
        acc = row.get("uniprot_accession", "").strip()
        if acc in nt_sequences:
            row["nucleotide_sequence"] = nt_sequences[acc]
        else:
            row["nucleotide_sequence"] = ""
            logging.debug(f"{acc}: No nucleotide sequence found in EMBOSS output.")

    with open(output_csv, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    logging.info(f"✅ Processing complete. Output written to: {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Reverse translate all protein sequences in a CSV using EMBOSS backtranseq."
    )
    parser.add_argument(
        "--pathogen_dir",
        required=True,
        help="Pathogen directory under data/ where logs and outputs are stored."
    )
    parser.add_argument(
        "--input_csv",
        required=True,
        help="Path to input CSV file containing protein sequences."
    )
    parser.add_argument(
        "--output_csv",
        required=True,
        help="Path to output CSV file to save results."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging to file and console."
    )
    args = parser.parse_args()

    pathogen_dir = Path("data") / args.pathogen_dir
    input_file = pathogen_dir / args.input_csv
    output_file = pathogen_dir / args.output_csv
    setup_logging(pathogen_dir, args.verbose)

    # check that the files and directories exist
    if not input_file.exists():
        logging.error(f"Input CSV file does not exist: {input_file}")
        sys.exit(1)

    if not pathogen_dir.exists():
        logging.error(f"Pathogen directory does not exist: {pathogen_dir}")
        sys.exit(1)


    check_backtranseq()
    process_csv(pathogen_dir, input_file, output_file)


if __name__ == "__main__":
    main()
