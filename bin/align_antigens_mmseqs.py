#!/usr/bin/env python3
import os
import csv
import sys
import subprocess
import tempfile
import shutil
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
    raw_result = results_dir / f"{strain_name}_alignment.tsv"
    best_result = results_dir / f"{strain_name}_best_hits.tsv"

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run([
            "mmseqs", "easy-search",
            antigen_fasta, str(strain_fasta),
            str(raw_result), tmpdir,
            "--format-mode", "4",
            "--format-output", "query,target,pident,nident,alnlen,evalue,bits"
        ], check=True)
        extract_best_hits(raw_result, best_result)
        print(f"[✓] {strain_name} aligned and best hits saved.")
    return strain_name

def extract_best_hits(tsv_path, output_path):
    best_hits = {}

    with open(tsv_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            query, target, pident = parts[0], parts[1], float(parts[2])

            if query not in best_hits or pident > best_hits[query]['pident']:
                best_hits[query] = {
                    'query': query,
                    'target': target,
                    'pident': pident,
                    'line': line.strip()
                }

    with open(output_path, 'w') as f_out:
        f_out.write("query\ttarget\tpident\tnident\talnlen\tevalue\tbits\n")
        for hit in best_hits.values():
            f_out.write(hit['line'] + "\n")

def main(pathogen_dir, num_threads):
    base_dir = Path(f"data/{pathogen_dir}")
    antigen_csv = base_dir / f"{pathogen_dir}_compiled_proteins.csv"
    strain_dir = base_dir / "strain_genomes"
    results_dir = base_dir / "mmseqs_results"
    results_dir.mkdir(exist_ok=True) # Ensure results directory exists

    #check if the files and directories exist
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

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: align_antigens_mmseqs.py <pathogen_directory> [--threads N]")
        sys.exit(1)
    pathogen_directory = sys.argv[1]
    threads = 4 # Default to 4 threads
    if "--threads" in sys.argv:
        idx = sys.argv.index("--threads")
        try:
            threads = int(sys.argv[idx + 1])
        except (IndexError, ValueError):
            print("Invalid value for --threads. Must be an integer.")
            sys.exit(1)
    if threads < 2:
        print("Please specify at least 2 threads with --threads.")
        sys.exit(1)

    main(pathogen_directory, threads)
