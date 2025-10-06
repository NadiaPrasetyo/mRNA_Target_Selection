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
    data/<pathogen_dir>/<output-dir>/*.pdb      # Downloaded PDB structures
    data/<pathogen_dir>/<output-dir>/*_AF.pdb   # Downloaded AlphaFold models (if PDB unavailable)

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

        return [item["identifier"] for item in result_set]

    except Exception as e:
        logging.exception(f"💥 Exception during UniProt search for {accession}: {e}")
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
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    logging.info(f"🧠 Attempting AlphaFold fetch for: {accession}")

    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            logging.warning(f"⚠️ AlphaFold fetch failed for {accession} (status {response.status_code})")
            return False

        predictions = response.json()
        if not predictions:
            logging.warning(f"⚠️ No AlphaFold prediction found for {accession}")
            return False

        pdb_url = predictions[0].get("pdbUrl")
        if not pdb_url:
            logging.warning(f"⚠️ No PDB URL in AlphaFold result for {accession}")
            return False

        dest = output_dir / f"{accession}_AF.pdb"
        logging.info(f"🔗 AlphaFold model found: {pdb_url}")
        logging.info(f"⬇️ Downloading AlphaFold PDB to: {dest.name}")

        r = requests.get(pdb_url)
        if r.status_code == 200:
            with open(dest, "wb") as f:
                f.write(r.content)
            logging.info(f"✅ AlphaFold model downloaded for {accession}")
            return True
        else:
            logging.error(f"⚠️ Failed to download AlphaFold PDB: HTTP {r.status_code}")
            return False

    except Exception as e:
        logging.error(f"💥 AlphaFold error for {accession}: {e}")
        return False


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
    suffix = accession if accession else "NOACCN"
    dest = output_dir / f"{pdb_id}_{suffix}.pdb"

    # Skip if good file already exists
    if dest.exists():
        if dest.stat().st_size > 0:
            logging.debug(f"{dest.name} already exists, skipping.")
            return
        else:
            logging.warning(f"⚠️ Removing empty file: {dest}")
            dest.unlink()

    # Primary download attempt: canonical PDB
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    logging.info(f"⬇️ Attempting PDB download for {pdb_id}...")

    try:
        subprocess.run(["wget", "-q", "-O", str(dest), url], check=True)
        if dest.stat().st_size > 0:
            logging.info(f"✅ Downloaded PDB: {dest.name}")
            return
        else:
            logging.warning(f"⚠️ Empty PDB file for {pdb_id}, removing.")
            dest.unlink(missing_ok=True)
    except subprocess.CalledProcessError:
        if dest.exists():
            logging.warning(f"⚠️ Removing partial file from failed canonical download: {dest}")
            dest.unlink()

        logging.warning(f"⚠️ Canonical PDB download failed for {pdb_id}. Trying biological assembly (.cif.gz)...")


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
        logging.warning(f"⚠️ No FASTA files found in: {sequence_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    all_records = [record for fasta in fasta_files for record in SeqIO.parse(fasta, "fasta")]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_record, record, output_dir) for record in all_records]

        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"❌ Error during record processing: {e}")

def was_new_file_created(before: set, after: set) -> bool:
    """Check if any new file was created in the output directory.
    Args:
        before (set): Set of files before the operation.
        after (set): Set of files after the operation.
    Returns:
        bool: True if a new file was created, False otherwise.
    This function compares the sets of files before and after an operation.
    It returns True if any new file was created (i.e., if the size of the after set is greater than the before set).
    If no new files were created, it returns False.
    """
    return any(f.stat().st_size > 0 for f in after - before)


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
                break

        # 💥 Always try AlphaFold if accession exists (brute-force mode)
        alphafold_downloaded = fetch_alphafold_structure(accession, output_dir)

    # Fallback: sequence search if no UniProt accession or no hits
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
        logging.info(f"✅ AlphaFold model retrieved for: {header}")
    else:
        logging.warning(f"❌ No structure available for: {header}")

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
    sequence_path = pathogen_path / args.sequence_dir
    output_path = pathogen_path / args.output_dir
    log_file = output_path / "fetch_pdb_sequences.log"

    setup_logger(args.verbose, log_file)
    logging.info(f"🚀 Starting structure fetch from {sequence_path}")

    process_fasta_dir(sequence_path, output_path, args.threads)
    logging.info("✅ Finished fetching structures.")


if __name__ == "__main__":
    main()
