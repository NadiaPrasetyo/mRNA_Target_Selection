#!/usr/bin/env python3
"""
fetch_PDB_structure.py
Command-line tool to fetch related PDB and AlphaFold structures for protein FASTA files using UniProt accession or sequence similarity.

Overview:
    - Parses all FASTA files in a specified directory.
    - Extracts UniProt accession numbers from FASTA headers using pattern matching.
    - Searches the RCSB PDB database for related structures using UniProt accession or, if unavailable, by sequence similarity (>90% identity).
    - Downloads canonical PDB files. If download fails, retrieves AlphaFold models instead.
    - Saves downloaded structures in a specified output directory.
    - Logs progress, warnings, and errors to the console and optionally to a log file.

Features:
    - Supports multithreading for faster processing of multiple FASTA records.
    - Automatically resolves redundancy by preferring PDB structures over AlphaFold models.
    - Provides detailed logging for debugging and tracking progress.
    - Handles incomplete or missing data gracefully, ensuring robust execution.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    sequence_dir (str): Subdirectory within pathogen_dir containing FASTA files.
    --threads (int): Number of threads to use for parallel processing (default: 4).
    --output-dir (Path): Directory to save fetched PDB/AlphaFold structures (default: pdb_sequences).
    --verbose: Enables verbose logging and writes logs to file.

Requirements:
    - FASTA files present in the specified sequence directory.
    - Python packages: biopython, requests.
    - wget command-line tool (for downloading PDB files).

Usage Example:
    python fetch_PDB_structure.py sars_cov_2 protein_fastas --threads 8 --output-dir pdbs --verbose

Outputs:
    - Exactly ONE structure file per accession:
        <pdb_id>_<accession>.pdb  (preferred)
        OR
        <accession>_AF.pdb        (if no PDB available)
    - Log file (if verbose mode is enabled):
        fetch_pdb_sequences.log

Author: Nadia
"""

import re
import logging
import argparse
import subprocess
from pathlib import Path
from typing import List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from Bio import SeqIO


SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


# ------------------- Logging ------------------- #
def setup_logger(verbose: bool, log_file: Path):
    """Set up logging configuration.
    Args:        
        verbose (bool): If True, enable verbose logging to file.
        log_file (Path): Path to the log file.
    """
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler()]
    if verbose:
        log_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure log directory exists
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(level=level, format=log_format, handlers=handlers)


# ------------------- Accession Extraction ------------------- #
def extract_accession(header: str) -> Optional[str]:
    """Extract UniProt accession from FASTA header.
    Args:
        header (str): FASTA header string.
    Returns:
        str: Extracted UniProt accession or None if not found.
    This function uses a regex pattern to match the UniProt accession format.
    If the header does not contain a valid accession, it returns None.
    """
    header_pattern = re.compile(r'(?:[^|]*\|)?(?P<accession>[A-Z0-9]+)(?:\||\s)')
    match = header_pattern.match(header)
    if match:
        return match.group("accession")
    return None


# ------------------- RCSB Search ------------------- #
def search_by_uniprot(accession: str) -> List[str]:
    """Search PDB entries by UniProt accession.
    Args:
        accession (str): UniProt accession number.
    Returns:
        List[str]: List of PDB IDs matching the UniProt accession.
    This function constructs a search query to the RCSB PDB API to find entries
    that match the given UniProt accession. It returns a list of identifiers for the matching PDB entries.
    If no entries are found, it returns an empty list.
    """
    logging.info(f"🔍 Searching PDB by UniProt accession: {accession}")

    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal",
                 "service": "text",
                 "parameters": {
                     "operator": "exact_match",
                     "value": accession,
                     "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"}},
                {"type": "terminal",
                 "service": "text",
                 "parameters": {
                     "operator": "exact_match",
                     "value": "UniProt",
                     "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"}}
            ]
        },
        "return_type": "entry"
    }

    try:
        r = requests.post(SEARCH_API_URL, json=payload, timeout=20)
        r.raise_for_status()
        data = r.json()
        hits = [x["identifier"] for x in data.get("result_set", [])]
        if not hits:
            logging.warning(f"⚠️ No PDB entries found for {accession}")
        return hits
    except Exception as e:
        logging.error(f"💥 UniProt search failed for {accession}: {e}")
        return []


def fetch_alphafold_structure(accession: str, output_dir: Path) -> bool:
    """Fetch AlphaFold structure for a given UniProt accession.
    Args:
        accession (str): UniProt accession number.
        output_dir (Path): Directory to save the AlphaFold PDB file.
    Returns:
        bool: True if AlphaFold model was successfully downloaded, False otherwise.
    This function queries the AlphaFold API for the specified UniProt accession.
    If a model is found, it downloads the PDB file and saves it in the specified output directory.
    If the model is not found or an error occurs, it logs a warning and returns False.
    """
    af_path = output_dir / f"{accession}_AF.pdb"
    if af_path.exists() and af_path.stat().st_size > 0:
        logging.debug(f"✅ AlphaFold already exists: {af_path.name}")
        return af_path
    
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    logging.info(f"🧠 Fetching AlphaFold for {accession}...")

    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            logging.warning(f"⚠️ AlphaFold fetch failed (status {r.status_code})")
            return None
        preds = r.json()
        if not preds or "pdbUrl" not in preds[0]:
            logging.warning(f"⚠️ No AlphaFold prediction found for {accession}")
            return None

        pdb_url = preds[0]["pdbUrl"]
        r_pdb = requests.get(pdb_url, timeout=30)
        if r_pdb.status_code == 200:
            with open(af_path, "wb") as f:
                f.write(r_pdb.content)
            logging.info(f"✅ AlphaFold model downloaded: {af_path.name}")
            return af_path
        else:
            logging.error(f"⚠️ Failed to download AlphaFold PDB: HTTP {r_pdb.status_code}")
            return None
    except Exception as e:
        logging.error(f"💥 AlphaFold error for {accession}: {e}")
        return None


def search_by_sequence(sequence: str) -> List[str]:
    """Search PDB entries by protein sequence similarity.
    Args:
        sequence (str): Protein sequence to search for.
    Returns:
        List[str]: List of PDB IDs matching the sequence.
    This function sends a sequence search request to the RCSB PDB API.
    It uses a POST request with the sequence as the search parameter.
    If the search is successful, it returns a list of PDB IDs that match the sequence.
    If the search fails or no results are found, it logs an error and returns an empty list.
    """
    logging.info("🔎 Performing sequence-based PDB search...")
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

    try:
        r = requests.post(SEARCH_API_URL, json=payload, timeout=30)
        r.raise_for_status()
        return [x["identifier"] for x in r.json().get("result_set", [])]
    except Exception as e:
        logging.error(f"💥 Sequence search failed {r.status_code}: {e} ")
        return []


# ------------------- Downloaders ------------------- #
def download_pdb(pdb_id: str, output_dir: Path, accession: str) -> Optional[Path]:
    """Download PDB file for a given PDB ID using wget.
    Args:
        pdb_id (str): PDB ID to download.
        output_dir (Path): Directory to save the downloaded PDB file.
        accession (Optional[str]): UniProt accession number for naming the file.
    If the accession is provided, it will be used to name the file.
    If not, a default suffix "NOACCN" will be used.
    The function first attempts to download the canonical PDB file.
    If that fails, it tries to download the biological assembly in .cif.gz format.
    If the canonical PDB file already exists and is non-empty, it skips the download.
    """    
    dest = output_dir / f"{pdb_id}_{accession}.pdb"

    if dest.exists() and dest.stat().st_size > 0:
        logging.debug(f"✅ PDB already exists: {dest.name}")
        return dest

    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    logging.info(f"⬇️ Attempting PDB download for {pdb_id}...")

    try:
        subprocess.run(["wget", "-q", "-O", str(dest), url], check=True)
        if dest.stat().st_size > 0:
            logging.info(f"✅ Downloaded PDB: {dest.name}")
            return dest
        else:
            logging.warning(f"⚠️ Empty PDB file for {pdb_id}, removing.")
            dest.unlink(missing_ok=True)
            return None
    except subprocess.CalledProcessError:
        logging.warning(f"⚠️ Failed to download canonical PDB: {pdb_id}")
        if dest.exists():
            logging.debug(f"Removing incomplete file: {dest.name}")
            dest.unlink()
        return None




def process_fasta_dir(sequence_dir: Path, output_dir: Path, threads: int):
    """Process all FASTA files in the specified directory.
    Args:
        sequence_dir (Path): Directory containing FASTA files.
        output_dir (Path): Directory to save fetched PDB sequences.
        threads (int): Number of threads to use for parallel processing.
    This function reads all FASTA files in the given directory, extracts UniProt accessions,
    searches for related PDB entries, and downloads the structures.
    It uses multithreading to speed up the processing of multiple records.
    """
    fasta_files = list(sequence_dir.glob("*.fasta"))
    if not fasta_files:
        logging.warning(f"⚠️ No FASTA files found in {sequence_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    records = [rec for f in fasta_files for rec in SeqIO.parse(f, "fasta")]

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [pool.submit(process_record, rec, output_dir) for rec in records]
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:
                logging.error(f"❌ Error processing record: {e}")

# ------------------- Core Processing ------------------- #
def process_record(record, output_dir: Path):
    """Process a single FASTA record to fetch PDB and AlphaFold structures.
    Args:
        record (SeqRecord): Biopython SeqRecord object containing the FASTA sequence.
        output_dir (Path): Directory to save the fetched PDB sequences.
    This function extracts the header and sequence from the record, searches for related PDB entries
    using the UniProt accession, and attempts to download the structures.
    If the accession is not found, it falls back to sequence-based search.
    It logs the progress and results of the operations.
    """
    header = record.description
    seq = str(record.seq)
    logging.info(f"\n🧬 Processing: {header}")

    accession = extract_accession(header)
    if not accession:
        logging.warning(f"⚠️ No accession found in: {header}")
        return

    logging.info(f"🔑 Accession: {accession}")

    # Step 1: Try PDB via accession
    pdb_file = None
    hits = search_by_uniprot(accession)
    if hits:
        pdb_file = download_pdb(hits[0], output_dir, accession)

    # Step 2: Always fetch AlphaFold
    af_file = fetch_alphafold_structure(accession, output_dir)

    # Step 3: Resolve redundancy
    if pdb_file and af_file:
        logging.info(f"🗑️ Removing AlphaFold for {accession} (PDB available)")
        af_file.unlink(missing_ok=True)

    # Step 4: Fallback to sequence if neither structure found
    if not pdb_file and not af_file:
        logging.info(f"🔁 Falling back to sequence-based search for {accession}")
        seq_hits = search_by_sequence(seq)
        if seq_hits:
            pdb_file = download_pdb(seq_hits[0], output_dir, accession)

    # Step 5: Final reporting
    if pdb_file:
        logging.info(f"✅ Final structure for {accession}: {pdb_file.name}")
    elif af_file:
        logging.info(f"✅ Final structure for {accession}: {af_file.name}")
    else:
        logging.warning(f"❌ No structure found for {accession}")


# ------------------- CLI ------------------- #
def main():
    """Main function to parse arguments and initiate PDB fetching."""
    parser = argparse.ArgumentParser(description="Fetch all related PDB entries from FASTA files using RCSB Web API. usage: python fetch_PDB_structure.py <pathogen_dir> <sequence_dir>")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("--output-dir", type=Path, default=Path("pdb_sequences"), help="Directory to save fetched PDB sequences (default: pdb_sequences)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output and log to file")
    args = parser.parse_args()

    pathogen_path = Path("data") / args.pathogen_dir
    seq_path = pathogen_path / args.sequence_dir
    out_path = pathogen_path / args.output_dir
    log_file = out_path / "fetch_pdb_sequences.log"

    setup_logger(args.verbose, log_file)
    logging.info(f"🚀 Starting structure fetch from {seq_path}")

    process_fasta_dir(seq_path, out_path, args.threads)

    logging.info("✅ Finished fetching structures.")


if __name__ == "__main__":
    main()
