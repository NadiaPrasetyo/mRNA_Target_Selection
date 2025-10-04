#!/usr/bin/env python3
"""
fetch_ncbi_strain_genome.py

Rewritten from fetch_NCBI_strain_genome.sh
Uses NCBI Datasets v2 API for fetching complete genomes.

Dependencies:
    pip install requests tqdm biopython seqkit
    # seqkit still used for fast translation

Usage:
    python fetch_ncbi_strain_genome.py --random "Staphylococcus aureus" \
        --random-num 5 data/staph_aureus

    python fetch_ncbi_strain_genome.py data/staph_aureus strains.csv
"""

import argparse
import csv
import logging
import os
import random
import shutil
import tempfile
import zipfile
import logging
import requests
from tqdm import tqdm

# ---------------------------
# Configuration
# ---------------------------
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
DEFAULT_THREADS = 4
DEFAULT_RANDOM_NUM = 5

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


# ---------------------------
# Utilities
# ---------------------------

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def download_and_extract_zip(url: str, accession: str, output_dir: str):
    """Download a genome ZIP and extract FASTA files."""
    r = requests.get(url, stream=True, headers={"accept": "application/zip"})
    r.raise_for_status()

    with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as tmp:
        for chunk in r.iter_content(chunk_size=8192):
            tmp.write(chunk)
        tmp_path = tmp.name

    with zipfile.ZipFile(tmp_path, "r") as zf:
        # Extract only .fna/.fasta
        for member in zf.namelist():
            if member.endswith((".fna", ".faa", ".fasta")):
                zf.extract(member, output_dir)
                # Rename .fna to .fasta for consistency
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
    """
    Return a list of dicts with assembly info for a taxon at complete genome level.
    Handles pagination to retrieve all available assemblies.
    """
    assemblies = []
    next_page_token = None

    while True:
        url = f"{API_BASE}/taxon/{requests.utils.quote(taxon)}/dataset_report"
        params = {
            "filters.assembly_level": "complete_genome",
            "page_size": 1000,  # maximum allowed page size
        }
        if next_page_token:
            params["page_token"] = next_page_token

        r = requests.get(url, params=params, headers={"accept": "application/json"})
        r.raise_for_status()
        data = r.json()
        total = data.get("total_count", 0)
        logging.info(f"Found {total} complete genomes for {taxon}")

        for record in data.get("reports", []):
            acc = record.get("accession")
            if acc:
                logging.info(f"Processing assembly {len(assemblies) + 1}/{total}: {acc}")
                assemblies.append({
                    "accession": acc
                })

        next_page_token = data.get("next_page_token")
        if not next_page_token:  # No more pages
            break

    return assemblies


def fetch_random_genomes(taxon: str, n: int, output_dir: str):
    logging.info(f"Fetching list of complete genomes for {taxon}...")
    assemblies = fetch_complete_genomes_for_taxon(taxon)

    if not assemblies:
        raise RuntimeError(f"No complete genomes found for {taxon}")

    chosen = random.sample(assemblies, min(n, len(assemblies)))
    logging.info(f"Selected {len(chosen)} genomes for download.")
    logging.info(f"Chosen assemblies: {chosen}")
    ensure_dir(output_dir)

    logging.info(f"Downloading {len(chosen)} random genomes...")
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
        raise RuntimeError("No valid EMBL/GenBank IDs found in CSV.")

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
    parser = argparse.ArgumentParser(
        description="Fetch NCBI genomes (random or CSV-based) and translate them."
    )
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

    shutil.rmtree(os.path.join(output_dir, "ncbi_dataset")) if os.path.exists(os.path.join(output_dir, "ncbi_dataset")) else None  # remove zip extraction dir if exists
    logging.info("Done.")


if __name__ == "__main__":
    main()
