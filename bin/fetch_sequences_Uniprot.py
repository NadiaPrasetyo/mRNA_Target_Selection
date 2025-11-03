"""
fetch_sequences_Uniprot.py
Command-line tool to fetch protein sequence and metadata from UniProt based on antigen data.

Overview:
    - Loads antigen records (antigen names, gene names, UniProt IDs) from a compiled CSV file.
    - Queries the UniProt API to retrieve full protein information for each antigen.
    - Fetches RefSeq nucleotide or protein sequences for matching entries, if available.
    - Parses and standardizes protein metadata including sequence, organism, domains, and features.
    - Compiles and saves the protein data into a new CSV file for downstream analysis.

Arguments:
    pathogen_directory (str): Subdirectory under `data/` containing pathogen data and antigen CSV.
    pathogen_name (str): Full name of the organism (used in UniProt queries).
    --output (str): Optional output CSV file path for compiled protein data. Defaults to
                    `data/<pathogen>/<organism>_compiled_proteins.csv`.
    --input (str): Optional input CSV file path with antigen data. Defaults to
                   `data/<pathogen>/<organism>_compiled_antigens.csv`.
    --fasta (bool): If specified, outputs protein data in FASTA format instead of CSV.

Requirements:
    - Input CSV file with antigen data present in the specified pathogen directory.
    - Python packages: argparse, csv, requests, os, re, unicodedata, time.

Usage Example:
    python fetch_sequences_Uniprot.py sars_cov_2 "SARS-CoV-2" --output proteins.csv --input antigens.csv

Outputs:
    - A CSV file (or FASTA file if --fasta is specified) containing compiled protein metadata for the specified organism, including:
        - UniProt accession
        - Protein name
        - Protein sequence
        - Organism name
        - Pfam domains
        - RefSeq nucleotide or protein sequences (if available)

Author: Nadia
"""
import argparse
import csv
import logging
import os
import random
import shutil
import tempfile
import zipfile
import time
import requests
from tqdm import tqdm
from functools import wraps

# ---------------------------
# Configuration
# ---------------------------
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
DEFAULT_THREADS = 4
DEFAULT_RANDOM_NUM = 5
RATE_LIMIT_RPS = 5          # 5 requests per second
MIN_REQUEST_INTERVAL = 1.0 / RATE_LIMIT_RPS


# ---------------------------
# Rate Limiting Decorator
# ---------------------------

_last_request_time = 0.0

def rate_limited_request(func):
    """Decorator to enforce rate limit and retry on 429 responses."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        global _last_request_time
        elapsed = time.time() - _last_request_time
        if elapsed < MIN_REQUEST_INTERVAL:
            time.sleep(MIN_REQUEST_INTERVAL - elapsed)

        for attempt in range(5):  # up to 5 retries
            response = func(*args, **kwargs)
            if response.status_code != 429:
                _last_request_time = time.time()
                return response

            retry_after = int(response.headers.get("Retry-After", 2))
            logging.warning(f"Rate limited (HTTP 429). Retrying after {retry_after}s...")
            time.sleep(retry_after)

        raise RuntimeError("Exceeded retry limit after repeated 429 responses.")
    return wrapper


@rate_limited_request
def safe_get(url, **kwargs):
    """Wrapper around requests.get with rate limiting and retries."""
    return requests.get(url, **kwargs)


# ---------------------------
# Utilities
# ---------------------------

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def download_and_extract_zip(url: str, accession: str, output_dir: str):
    """Download a genome ZIP and extract FASTA files."""
    r = safe_get(url, stream=True, headers={"accept": "application/zip"})
    r.raise_for_status()

    with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as tmp:
        for chunk in r.iter_content(chunk_size=8192):
            tmp.write(chunk)
        tmp_path = tmp.name

    with zipfile.ZipFile(tmp_path, "r") as zf:
        for member in zf.namelist():
            if member.endswith((".fna", ".faa", ".fasta")):
                zf.extract(member, output_dir)
                if member.endswith(".fna"):
                    shutil.move(
                        os.path.join(output_dir, member),
                        os.path.join(output_dir, f"{accession}.fasta")
                    )
                elif member.endswith(".faa"):
                    shutil.move(
                        os.path.join(output_dir, member),
                        os.path.join(output_dir, f"{accession}_proteins.fasta")
                    )

    os.remove(tmp_path)


# ---------------------------
# Fetching genomes
# ---------------------------

def fetch_complete_genomes_for_taxon(taxon: str) -> list:
    """Return a list of dicts with assembly info for a taxon at complete genome level."""
    assemblies = []
    next_page_token = None
    page_size = 500  # smaller requests to avoid 504 timeouts

    logging.info(f"Querying NCBI Datasets API for '{taxon}' complete genomes...")

    while True:
        url = f"{API_BASE}/taxon/{requests.utils.quote(taxon)}/dataset_report"
        params = {
            "filters.assembly_level": "complete_genome",
            "page_size": page_size,
        }
        if next_page_token:
            params["page_token"] = next_page_token

        # retry loop for 5xx errors
        for attempt in range(5):
            try:
                r = safe_get(url, params=params, headers={"accept": "application/json"})
                if r.status_code >= 500:
                    raise requests.exceptions.HTTPError(f"{r.status_code} Server Error")
                r.raise_for_status()
                break
            except requests.exceptions.RequestException as e:
                wait_time = 2 ** attempt
                logging.warning(f"Error fetching page ({e}). Retrying in {wait_time}s...")
                time.sleep(wait_time)
        else:
            raise RuntimeError(f"Failed after multiple retries for taxon: {taxon}")

        data = r.json()
        total = data.get("total_count", 0)
        logging.info(f"Found {total} complete genomes for {taxon}")

        for record in data.get("reports", []):
            acc = record.get("accession")
            if acc:
                assemblies.append({"accession": acc})

        next_page_token = data.get("next_page_token")
        if not next_page_token:
            break

        # pause slightly between page requests (stay within rate limit)
        time.sleep(MIN_REQUEST_INTERVAL)

    logging.info(f"Retrieved {len(assemblies)} total assemblies for {taxon}")
    return assemblies


def fetch_random_genomes(taxon: str, n: int, output_dir: str):
    logging.info(f"Fetching list of complete genomes for {taxon}...")
    assemblies = fetch_complete_genomes_for_taxon(taxon)

    if not assemblies:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    chosen = random.sample(assemblies, min(n, len(assemblies)))
    logging.info(f"Selected {len(chosen)} genomes for download.")
    ensure_dir(output_dir)

    for asm in tqdm(chosen, desc="Downloading"):
        acc = asm["accession"]
        url = f"{API_BASE}/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        download_and_extract_zip(url, acc, output_dir)

    logging.info(f"Saved genomes in {output_dir}")


def fetch_from_csv(csv_path: str, output_dir: str):
    ensure_dir(output_dir)

    ids = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            embl_id = row.get("RefSeq_ID")
            if embl_id:
                ids.append(embl_id)

    if not ids:
        raise RuntimeError("No valid RefSeq IDs found in CSV.")

    logging.info(f"Fetching {len(ids)} genomes from CSV...")
    for acc in tqdm(ids, desc="Downloading"):
        url = f"{API_BASE}/accession/{acc}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA&hydrated=FULLY_HYDRATED"
        try:
            download_and_extract_zip(url, acc, output_dir)
        except Exception as e:
            logging.warning(f"Failed to fetch {acc}: {e}")

    logging.info(f"Saved genomes in {output_dir}")


# ---------------------------
# Main
# ---------------------------

def main():
    parser = argparse.ArgumentParser(description="Fetch NCBI genomes (random or CSV-based).")
    parser.add_argument("pathogen_dir", help="Output directory under data/")
    parser.add_argument("--csv_file", nargs="?", help="CSV file with strain names and IDs")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--random", dest="random_pathogen", help="Fetch random genomes for pathogen")
    parser.add_argument("--random-num", type=int, default=DEFAULT_RANDOM_NUM)
    parser.add_argument("--output-dir", help="Output directory (default: pathogen_dir/strain_genomes)", default="strain_genomes")

    args = parser.parse_args()
    setup_logging()

    base_dir = os.path.join("data", args.pathogen_dir)
    output_dir = os.path.join(base_dir, args.output_dir)
    ensure_dir(output_dir)

    if args.random_pathogen:
        fetch_random_genomes(args.random_pathogen, args.random_num, output_dir)
    else:
        if not args.csv_file:
            parser.error("CSV file required when not using --random")
        csv_path = os.path.join(base_dir, args.csv_file)
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        fetch_from_csv(csv_path, output_dir)

    # Cleanup
    ncbi_tmp_dir = os.path.join(output_dir, "ncbi_dataset")
    if os.path.exists(ncbi_tmp_dir):
        shutil.rmtree(ncbi_tmp_dir)

    logging.info("Done.")


if __name__ == "__main__":
    main()