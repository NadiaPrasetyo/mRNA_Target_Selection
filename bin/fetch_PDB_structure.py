"""
fetch_PDB_structure.py
Fetches related PDB sequences for protein FASTA files using UniProt accession or sequence similarity.

This script automates the retrieval of experimentally determined PDB sequences corresponding to protein sequences provided in FASTA format. It first attempts to map each sequence to PDB entries using UniProt accession numbers parsed from FASTA headers. If no match is found, it performs a sequence-based search with high identity cutoff. Retrieved PDB sequences are saved in a specified output directory.

General Function:
    - Parses all FASTA files in a given directory.
    - Extracts UniProt accession numbers from FASTA headers using robust pattern matching.
    - Queries the RCSB PDB database for related structures using either UniProt accession or sequence similarity (>90% identity).
    - Downloads canonical PDB sequences for matched entries and saves them as text files.
    - Logs progress, warnings, and errors to console and optionally to a log file.

Arguments:
    - pathogen_dir: Path to the pathogen-specific directory under data/.
    - sequence_dir: Subdirectory containing FASTA files within pathogen_dir.
    - --threads: Number of threads to use (currently not utilized).
    - --output-dir: Directory to save fetched PDB sequences (default: pdb_sequences).
    - --verbose: Enables verbose logging and writes logs to file.

Dependencies:
    - Biopython (for FASTA parsing)
    - rcsbapi (for RCSB PDB queries)
    - Python standard libraries: os, argparse, re, logging, pathlib

Author: Nadia
"""
import os
import re
import json
import logging
import argparse
import subprocess
from pathlib import Path
from typing import List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from Bio import SeqIO


SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def setup_logger(verbose: bool, log_file: Path):
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if verbose:
        log_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure log directory exists
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=log_format, handlers=handlers)


def extract_uniprot_accession(header: str) -> Optional[str]:
    pipe_parts = header.split("|")
    for part in pipe_parts:
        if re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]", part):
            return part

    first_word = header.split()[0].replace(">", "")
    if re.fullmatch(r"[A-NR-Z][0-9][A-Z0-9]{3}[0-9]", first_word):
        return first_word

    return None


def search_by_uniprot(accession: str) -> List[str]:
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": accession,
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "operator": "exact_match",
                        "value": "UniProt",
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"
                    }
                }
            ]
        },
        "return_type": "entry"
    }

    response = requests.post(SEARCH_API_URL, json=payload)

    if response.status_code != 200 or not response.json():
        logging.error(f"Failed UniProt search: {response.status_code}")
        return []  # 🔧 always return list

    data = response.json()
    return [item["identifier"] for item in data.get("result_set", [])]

def fetch_alphafold_structure(accession: str, output_dir: Path) -> bool:
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    logging.info(f"🧠 Attempting AlphaFold fetch for: {accession}")

    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            logging.warning(f"AlphaFold fetch failed for {accession} (status {response.status_code})")
            return False

        predictions = response.json()
        if not predictions:
            logging.warning(f"No AlphaFold prediction found for {accession}")
            return False

        best_model = predictions[0]
        pdb_url = best_model.get("pdbUrl")
        if not pdb_url:
            logging.warning(f"No PDB URL in AlphaFold result for {accession}")
            return False

        filename = f"{accession}_AF.pdb"
        dest = output_dir / filename
        logging.info(f"🔗 AlphaFold model found: {pdb_url}")
        logging.info(f"⬇️  Downloading AlphaFold PDB to: {dest.name}")

        r = requests.get(pdb_url)
        if r.status_code == 200:
            with open(dest, "wb") as f:
                f.write(r.content)
            return True
        else:
            logging.error(f"❌ Failed to download AlphaFold PDB: HTTP {r.status_code}")
            return False

    except Exception as e:
        logging.error(f"AlphaFold error for {accession}: {e}")
        return False


def search_by_sequence(sequence: str) -> List[str]:
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": 0.9,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined"
        }
    }

    response = requests.post(SEARCH_API_URL, json=payload)
    if response.status_code != 200:
        logging.error(f"Failed sequence search: {response.status_code}")
        return []

    data = response.json()
    return [item["identifier"] for item in data.get("result_set", [])]

def download_pdb(pdb_id: str, output_dir: Path, accession: Optional[str]):
    suffix = accession if accession else "NOACCN"
    filename = f"{pdb_id}_{suffix}.pdb"
    dest = output_dir / filename

    if dest.exists():
        logging.debug(f"{filename} already exists, skipping.")
        return

    url_standard = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    logging.info(f"⬇️  Attempting canonical download for {pdb_id}...")

    try:
        subprocess.run(["wget", "-q", "-O", str(dest), url_standard], check=True)
        logging.info(f"✅ Downloaded {filename} successfully.")
        return
    except subprocess.CalledProcessError:
        logging.warning(f"⚠️ Canonical PDB download failed for {pdb_id}. Trying biological assembly...")

    # Fall back to assembly ID 1 (most common case)
    assembly_url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
    gz_filename = f"{pdb_id}_assembly1_{suffix}.pdb.gz"
    gz_dest = output_dir / gz_filename

    try:
        subprocess.run(["wget", "-q", "-O", str(gz_dest), assembly_url], check=True)
        logging.info(f"✅ Downloaded biological assembly (pdb1.gz) for {pdb_id} as {gz_filename}")
    except subprocess.CalledProcessError:
        logging.error(f"❌ Failed to download biological assembly for {pdb_id}")



def process_fasta_dir(sequence_dir: Path, output_dir: Path, threads: int):
    fasta_files = list(sequence_dir.glob("*.fasta"))
    if not fasta_files:
        logging.warning(f"No FASTA files found in: {sequence_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    all_records = []

    for fasta_file in fasta_files:
        all_records.extend(list(SeqIO.parse(fasta_file, "fasta")))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_record, record, output_dir) for record in all_records]

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"❌ Error during record processing: {e}")


def process_record(record, output_dir: Path):
    header = record.description
    sequence = str(record.seq)
    logging.info(f"🧬 Processing: {header}")

    accession = extract_uniprot_accession(header)
    pdb_ids = []
    alphafold_downloaded = False

    if accession:
        logging.info(f"🔍 Searching by UniProt accession: {accession}")
        pdb_ids = search_by_uniprot(accession)

        if not pdb_ids:
            alphafold_downloaded = fetch_alphafold_structure(accession, output_dir)

    if not pdb_ids and not alphafold_downloaded:
        logging.info(f"🔁 Falling back to sequence-based search...")
        pdb_ids = search_by_sequence(sequence)

    # ❌ Only warn if no structure at all (PDB nor AlphaFold)
    if not pdb_ids and not alphafold_downloaded:
        logging.warning(f"❌ No PDB entries found for: {header}")

    for pdb_id in pdb_ids:
        download_pdb(pdb_id, output_dir, accession)


def main():
    parser = argparse.ArgumentParser(description="Fetch all related PDB entries from FASTA files using RCSB Web API.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--output-dir", type=Path, default=Path("pdb_sequences"), help="Directory to save fetched PDB sequences (default: pdb_sequences)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output and log to file")
    args = parser.parse_args()

    pathogen_path = Path("data") / args.pathogen_dir
    output_path = pathogen_path / args.output_dir
    full_sequence_path = pathogen_path / args.sequence_dir
    log_file = output_path / "fetch_pdb_sequences.log"


    setup_logger(args.verbose, log_file)

    logging.info(f"🚀 Starting PDB fetch from {full_sequence_path}")
    process_fasta_dir(full_sequence_path, output_path, args.threads)
    logging.info("✅ Finished.")


if __name__ == "__main__":
    main()
