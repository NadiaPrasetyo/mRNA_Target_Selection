"""
/**
 * @file align_antigens.py
 * @brief Aligns antigen protein sequences against strain genome sequences using MMseqs2.
 *
 * This script processes compiled antigen protein data and aligns them against translated
 * strain genome sequences using the MMseqs2 tool. It extracts best matches and saves
 * aligned antigen sequences and metadata to structured output files.
 *
 * General Flow:
 *   1. Converts antigen CSV to FASTA format.
 *   2. Searches each strain’s translated genome for matching antigens using MMseqs2.
 *   3. Extracts best hits and saves results to TSV and FASTA files.
 *   4. Optionally includes query sequences in output.
 *
 * Parameters:
 *   pathogen_directory (str): Folder under `data/` where the antigen and strain data are located.
 *   pathogen_name (str): Name of the pathogen (used to infer antigen file names).
 *   threads (int): Number of parallel workers to use.
 *   output_dir (str): Name of subdirectory for results under `data/<pathogen_directory>/`.
 *   fetch_qseq (bool): Whether to include query sequences in MMseqs2 output.
 *
 * Usage:
 *   python align_antigens.py <pathogen_directory> <pathogen_name> [--threads N] [--output-dir X] [--fetch-qseq]
 *
 * Example:
 *   python align_antigens.py sars_cov_2 "SARS-CoV-2" --threads 8 --fetch-qseq
 *
 * Output:
 *   - <output_dir>/<strain>_alignment.tsv
 *   - <output_dir>/<strain>_best_hits.tsv
 *   - <output_dir>/<strain>_matched_antigens.fasta
 *
 * Requires:
 *   - MMseqs2 installed and accessible in PATH
 *   - Antigen protein CSV with uniprot_accession, protein_name, sequence fields
 *   - Translated strain genome FASTA files in the format *_translated.fasta
 *
 * Author: [Your Name]
 */
"""

#!/usr/bin/env python3
import os
import csv
import sys
import subprocess
import tempfile
import shutil
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

"""
/**
 * @brief Converts antigen CSV records into a FASTA file.
 *
 * Each row in the input CSV is formatted into a FASTA record using accession,
 * protein name, and sequence fields.
 *
 * @param csv_path (str): Path to the antigen CSV file.
 * @param fasta_path (str): Output path for the generated FASTA file.
 * @return: None
 */
"""
def extract_antigens_to_fasta(csv_path, fasta_path):
    with open(csv_path, newline='') as csvfile, open(fasta_path, 'w') as f_out:
        reader = csv.DictReader(csvfile)
        for idx, row in enumerate(reader):
            acc = row['uniprot_accession']
            name = row['protein_name']
            seq = row['sequence'].replace('\r', '').replace('\n', '')
            f_out.write(f">antigen_{idx}|{acc}|{name}\n{seq}\n")

"""
/**
 * @brief Runs MMseqs2 easy-search and processes alignment results.
 *
 * Executes alignment between antigens and a strain genome. Then extracts
 * the best hits and sequences from MMseqs2 output.
 *
 * @param strain_fasta_path (str): Path to strain's translated genome FASTA.
 * @param antigen_fasta (str): Path to antigen FASTA file.
 * @param results_dir (str): Directory to save output.
 * @param fetch_qseq (bool): If True, includes query sequences in the output.
 * @return: strain_name (str)
 */
"""
def run_mmseqs2_and_process(strain_fasta_path, antigen_fasta, results_dir, fetch_qseq=False):
    strain_fasta = Path(strain_fasta_path)
    strain_name = strain_fasta.stem.replace("_translated", "")
    raw_result = results_dir / f"{strain_name}_alignment.tsv"
    best_result = results_dir / f"{strain_name}_best_hits.tsv"
    antigen_seqs_out = results_dir / f"{strain_name}_matched_antigens.fasta"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_fields = [
            "query", "target", "pident", "nident", "alnlen",
            "evalue", "bits", "mismatch", "qcov", "tcov", "tstart", "tend", "taln"
        ]
        if fetch_qseq:
            output_fields.append("qseq")

        subprocess.run([
            "mmseqs", "easy-search",
            antigen_fasta, str(strain_fasta),
            str(raw_result), tmpdir,
            "--format-mode", "4",
            "--format-output", ",".join(output_fields)
        ], check=True)

        extract_best_hits_with_sequences(raw_result, best_result, antigen_seqs_out, fetch_qseq)

    print(f"[✓] {strain_name} aligned. Hits + sequences saved.")
    return strain_name

"""
/**
 * @brief Extracts best alignment hits and writes them to TSV and FASTA.
 *
 * For each query, the best hit (based on percent identity) is retained and written to
 * a TSV summary and a FASTA with matched target sequences.
 *
 * @param raw_tsv_path (str): Path to raw MMseqs2 alignment TSV file.
 * @param output_tsv_path (str): Output path for filtered best hits TSV.
 * @param fasta_out_path (str): Output path for matched sequences FASTA.
 * @param fetch_qseq (bool): Whether to include query sequences in output.
 * @return: None
 */
"""
def extract_best_hits_with_sequences(raw_tsv_path, output_tsv_path, fasta_out_path, fetch_qseq):
    best_hits = {}
    with open(raw_tsv_path) as f:
        lines = f.readlines()
        if lines and lines[0].lower().startswith("query"):
            lines = lines[1:]

        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue

            try:
                query, target = parts[0], parts[1]
                pident = float(parts[2])
                evalue = parts[5]
                mismatch = parts[7]
                qcov = parts[8]
                tcov = parts[9]
                tstart, tend = parts[10], parts[11]
                taln = parts[12]
                qseq = parts[13] if fetch_qseq and len(parts) > 13 else ""
            except (IndexError, ValueError):
                continue

            if query not in best_hits or pident > best_hits[query]['pident']:
                best_hits[query] = {
                    'query': query,
                    'target': target,
                    'pident': pident,
                    'evalue': evalue,
                    'mismatch': mismatch,
                    'qcov': qcov,
                    'tcov': tcov,
                    'tstart': tstart,
                    'tend': tend,
                    'taln': taln,
                    'qseq': qseq,
                }

    headers = ["query", "target", "pident", "evalue", "mismatch", "qcov", "tcov", "tstart", "tend", "taln"]
    if fetch_qseq:
        headers.append("qseq")

    with open(output_tsv_path, 'w') as f_out:
        f_out.write('\t'.join(headers) + '\n')
        for hit in best_hits.values():
            row = [str(hit[h]) for h in headers]
            f_out.write('\t'.join(row) + '\n')

    with open(fasta_out_path, 'w') as fasta_out:
        for hit in best_hits.values():
            header = f"{hit['query']}|{hit['target']}|tpos:{hit['tstart']}-{hit['tend']}"
            fasta_out.write(f">{header}\n{hit['taln']}\n")

"""
/**
 * @brief Entry point to execute the antigen-to-strain alignment workflow.
 *
 * Validates input files and directories, prepares antigen FASTA, then runs
 * MMseqs2 alignments in parallel for each strain.
 *
 * @param pathogen_dir (str): Directory name under `data/`.
 * @param pathogen_name (str): Scientific name of the organism.
 * @param num_threads (int): Number of parallel processes to use.
 * @param output_dir (str): Output subdirectory name.
 * @param fetch_qseq (bool): If True, include query sequences in results.
 * @return: None
 */
"""
def main(pathogen_dir, pathogen_name, num_threads, output_dir, fetch_qseq):
    base_dir = Path(f"data/{pathogen_dir}")
    pathogen_tag = pathogen_name.replace(" ", "_").lower()
    antigen_csv = base_dir / f"{pathogen_tag}_compiled_proteins.csv"
    strain_dir = base_dir / "strain_genomes"
    results_dir = Path(base_dir/output_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    if not antigen_csv.exists():
        print(f"Error: Antigen CSV file {antigen_csv} does not exist.")
        sys.exit(1)
    if not strain_dir.exists() or not any(strain_dir.glob("*_translated.fasta")):
        print(f"Error: No strain FASTA files found in {strain_dir}.")
        sys.exit(1)
    if not shutil.which("mmseqs"):
        print("Error: MMseqs2 is not installed or not found in PATH.")
        sys.exit(1)
    if not base_dir.exists() or not base_dir.is_dir():
        print(f"Error: Base directory {base_dir} does not exist.")
        sys.exit(1)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        antigen_fasta = tmp.name
    extract_antigens_to_fasta(antigen_csv, antigen_fasta)

    strain_fastas = list(strain_dir.glob("*_translated.fasta"))
    print(f"Running MMseqs2 on {len(strain_fastas)} strains with {num_threads} workers...")

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {
            executor.submit(run_mmseqs2_and_process, f, antigen_fasta, results_dir, fetch_qseq): f
            for f in strain_fastas
        }
        for future in as_completed(futures):
            future.result()

    os.remove(antigen_fasta)
    print("All alignments complete.")

"""
/**
 * @brief CLI wrapper for align_antigens.py.
 *
 * Parses command-line arguments and invokes the main workflow.
 */
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Align antigens to strain genomes using MMseqs2"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus aureus")')
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Custom output directory (default: data/<pathogen_directory>/mmseqs_results)"
    )
    parser.add_argument(
        "--fetch-qseq",
        action="store_true",
        help="Include qseq (query sequence) in MMseqs2 output"
    )

    args = parser.parse_args()

    if args.threads < 2:
        print("Please specify at least 2 threads with --threads.")
        sys.exit(1)

    # Set default output dir if not provided
    if args.output_dir is None:
        args.output_dir = f"mmseqs_results"

    main(args.pathogen_directory, args.pathogen_name, args.threads, args.output_dir, args.fetch_qseq)
