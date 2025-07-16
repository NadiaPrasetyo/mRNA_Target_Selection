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

def extract_accession(header: str) -> Optional[str]:
    header_pattern = re.compile(r'(?:[^|]*\|)?(?P<accession>[A-Z0-9]+)(?:\||\s)')
    match = header_pattern.match(header)
    if match:
        return match.group("accession")
    return None

def search_by_uniprot(accession: str) -> List[str]:
    logging.info(f"🔍 Searching PDB by UniProt accession: {accession}")
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

    try:
        response = requests.post(SEARCH_API_URL, json=payload)
        if response.status_code != 200:
            logging.error(f"⚠️ UniProt search failed for {accession} - Status code: {response.status_code}")
            return []

        data = response.json()
        result_set = data.get("result_set", [])

        if not result_set:
            logging.warning(f"⚠️ No PDB entries found for UniProt accession: {accession}")
            return []

        ids = [item["identifier"] for item in result_set]
        logging.info(f"📦 PDB IDs for {accession}: {ids}")
        return ids

    except Exception as e:
        logging.exception(f"💥 Exception during UniProt search for {accession}: {e}")
        return []


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
            logging.error(f"⚠️ Failed to download AlphaFold PDB: HTTP {r.status_code}")
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

    # Skip if good file already exists
    if dest.exists():
        if dest.stat().st_size > 0:
            logging.debug(f"{filename} already exists, skipping.")
            return
        else:
            logging.warning(f"⚠️ Removing empty file: {dest}")
            dest.unlink()

    # Primary download attempt: canonical PDB
    url_standard = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    logging.info(f"⬇️  Attempting canonical PDB download for {pdb_id}...")

    try:
        subprocess.run(["wget", "-q", "-O", str(dest), url_standard], check=True)
        if dest.stat().st_size > 0:
            logging.info(f"✅ Downloaded {filename} successfully.")
            return
        else:
            logging.warning(f"⚠️ Canonical PDB file is empty: {dest}. Removing.")
            dest.unlink()
    except subprocess.CalledProcessError:
        if dest.exists():
            logging.warning(f"⚠️ Removing partial file from failed canonical download: {dest}")
            dest.unlink()

        logging.warning(f"⚠️ Canonical PDB download failed for {pdb_id}. Trying biological assembly (.cif.gz)...")

    # Fallback: biological assembly in .cif.gz format
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    cif_filename = f"{pdb_id}_assembly1_{suffix}.cif.gz"
    cif_dest = output_dir / cif_filename

    try:
        subprocess.run(["wget", "-q", "-O", str(cif_dest), cif_url], check=True)
        if cif_dest.exists() and cif_dest.stat().st_size > 0:
            logging.info(f"✅ Downloaded biological assembly (.cif.gz) for {pdb_id} as {cif_filename}")
        else:
            logging.error(f"⚠️ Biological assembly .cif.gz is empty for {pdb_id}. Removing.")
            if cif_dest.exists():
                cif_dest.unlink()
    except subprocess.CalledProcessError:
        logging.error(f"⚠️ Failed to download biological assembly (.cif.gz) for {pdb_id}")
        if cif_dest.exists():
            cif_dest.unlink()



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

def was_new_file_created(before: set, after: set) -> bool:
    return any(f.stat().st_size > 0 for f in after - before)


def process_record(record, output_dir: Path):
    header = record.description
    sequence = str(record.seq)
    logging.info(f"🧬 Processing: {header}")

    accession = extract_accession(header)
    logging.info(f"🔍 Extracted UniProt accession: {accession if accession else 'None'}")
    pdb_ids = []
    pdb_downloaded = False
    alphafold_downloaded = False

    if accession:
        pdb_ids = search_by_uniprot(accession)

        for pdb_id in pdb_ids:
            before = set(output_dir.glob("*"))
            download_pdb(pdb_id, output_dir, accession)
            after = set(output_dir.glob("*"))
            if was_new_file_created(before, after):
                pdb_downloaded = True

        # 💥 Always try AlphaFold if accession exists (brute-force mode)
        alphafold_downloaded = fetch_alphafold_structure(accession, output_dir)

    # 🪂 Fallback to sequence-based search only if neither method yielded results
    if not pdb_downloaded and not alphafold_downloaded:
        logging.info("🔁 Falling back to sequence-based search...")
        pdb_ids = search_by_sequence(sequence)

        for pdb_id in pdb_ids:
            before = set(output_dir.glob("*"))
            download_pdb(pdb_id, output_dir, accession)
            after = set(output_dir.glob("*"))
            if was_new_file_created(before, after):
                pdb_downloaded = True

    # 🧾 Final logging
    if pdb_downloaded and alphafold_downloaded:
        logging.info(f"✅ PDB and AlphaFold models retrieved for: {header}")
        # delete the AlphaFold file if PDB was downloaded
        af_file = output_dir / f"{accession}_AF.pdb"
        if af_file.exists():
            logging.info(f"🗑️ Deleting AlphaFold file: {af_file.name} (PDB found)")
            af_file.unlink()
    elif pdb_downloaded:
        logging.info(f"✅ PDB structure(s) retrieved for: {header} (AlphaFold not found)")
    elif alphafold_downloaded:
        logging.info(f"✅ AlphaFold model retrieved for: {header} (PDB not found)")
    else:
        logging.warning(f"❌ All structure retrieval attempts failed for: {header}")


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
