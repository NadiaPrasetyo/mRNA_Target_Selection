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

def extract_antigens_to_fasta(csv_path, fasta_path):
    with open(csv_path, newline='') as csvfile, open(fasta_path, 'w') as f_out:
        reader = csv.DictReader(csvfile)
        for idx, row in enumerate(reader):
            acc = row['uniprot_accession']
            name = row['protein_name']
            seq = row['sequence'].replace('\r', '').replace('\n', '')
            f_out.write(f">antigen_{idx}|{acc}|{name}\n{seq}\n")

def run_mmseqs2_and_process(strain_fasta_path, antigen_fasta, results_dir):
    strain_fasta = Path(strain_fasta_path)
    strain_name = strain_fasta.stem.replace("_translated", "")
    raw_result = results_dir / f"{strain_name}_alignment.csv"
    best_result = results_dir / f"{strain_name}_best_hits.csv"
    antigen_seqs_out = results_dir / f"{strain_name}_matched_antigens.fasta"

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run([
            "mmseqs", "easy-search",
            antigen_fasta, str(strain_fasta),
            str(raw_result), tmpdir,
            "--format-mode", "4",
            "--format-output", "query,target,pident,nident,alnlen,evalue,bits,tstart,tend,qseq,tseq"
        ], check=True)

        extract_best_hits_with_sequences(raw_result, best_result, antigen_seqs_out)
        print(f"[✓] {strain_name} aligned. Hits + sequences saved.")
    return strain_name

def extract_best_hits_with_sequences(csv_path, output_path, fasta_out_path):
    best_hits = {}

    with open(csv_path) as f:
        print("First line:", line.strip())
        header_skipped = False
        for line in f:
            if not header_skipped:
                header_skipped = True
                continue  # skip header row
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            query, target, pident = parts[0], parts[1], float(parts[2])
            tstart, tend = parts[7], parts[8]
            qseq, tseq = parts[9], parts[10]

            if query not in best_hits or pident > best_hits[query]['pident']:
                best_hits[query] = {
                    'query': query,
                    'target': target,
                    'pident': pident,
                    'tstart': tstart,
                    'tend': tend,
                    'qseq': qseq,
                    'tseq': tseq,
                    'line': line.strip()
                }

    # Save best hits TSV
    with open(output_path, 'w') as f_out:
        f_out.write("query\ttarget\tpident\tnident\talnlen\tevalue\tbits\ttstart\ttend\tqseq\ttseq\n")
        for hit in best_hits.values():
            f_out.write(hit['line'] + "\n")

    # Write target matched sequences to FASTA
    with open(fasta_out_path, 'w') as fasta_out:
        for hit in best_hits.values():
            header = f"{hit['query']}|{hit['target']}|tpos:{hit['tstart']}-{hit['tend']}"
            fasta_out.write(f">{header}\n{hit['tseq']}\n")


def main(pathogen_dir, pathogen_name, num_threads, output_dir):
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
            executor.submit(run_mmseqs2_and_process, f, antigen_fasta, results_dir): f
            for f in strain_fastas
        }
        for future in as_completed(futures):
            future.result()

    os.remove(antigen_fasta)
    print("All alignments complete.")

import argparse

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

    args = parser.parse_args()

    if args.threads < 2:
        print("Please specify at least 2 threads with --threads.")
        sys.exit(1)

    # Set default output dir if not provided
    if args.output_dir is None:
        args.output_dir = f"mmseqs_results"

    main(args.pathogen_directory, args.pathogen_name, args.threads, args.output_dir)

