"""
/**
 * @file generate_random_proteins_per_strain.py
 * @brief Fetches random reviewed UniProt protein entries for a given organism, excluding known antigens.
 *
 * This script loads known antigen protein names and their sequence length bounds from a compiled CSV file,
 * then queries the UniProt API to randomly sample reviewed protein entries matching the organism and
 * sequence criteria. It excludes any protein names that match known antigens. Parsed entries are saved
 * in a new CSV file for further analysis or use.
 *
 * General Flow:
 *   1. Load antigen names and sequence length bounds from a CSV file.
 *   2. Query the UniProt API for random reviewed proteins within the same sequence length range.
 *   3. Parse the UniProt entries using a shared parsing function.
 *   4. Save the parsed protein information to a new CSV file.
 *
 * Parameters:
 *   pathogen (str): Short identifier for the organism's data folder (used as a subdirectory under data/).
 *   organism (str): Full organism name (used for querying UniProt and naming output files).
 *
 * Usage:
 *   python generate_random_proteins_per_strain.py <pathogen_subfolder> <organism_name>
 *
 * Example:
 *   python generate_random_proteins_per_strain.py s_aureus "Staphylococcus aureus"
 *
 * @author Nadia
 */
"""

import os
import csv
import random
import requests
import argparse
from bin.fetch_sequences_Uniprot import parse_uniprot_entry  # ✅ using shared parser

"""
/**
 * @brief Extracts and lowercases antigen protein names from a CSV file.
 *
 * Reads the input file and collects a set of all known antigen protein names,
 * converted to lowercase and stripped of whitespace.
 *
 * @param antigen_file (str): Path to the compiled antigens CSV file.
 * @return (set): Set of known antigen protein names (lowercased).
 */
"""
def get_antigen_protein_names(antigen_file):
    names = set()
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row.get("protein_name", "").strip().lower()
            if name:
                names.add(name)
    return names

"""
/**
 * @brief Determines the minimum and maximum lengths of antigen sequences.
 *
 * Reads the antigen CSV file and calculates sequence length bounds to use for filtering
 * UniProt proteins.
 *
 * @param antigen_file (str): Path to the compiled antigens CSV file.
 * @return (tuple): Minimum and maximum sequence lengths (int, int), or (None, None) if unavailable.
 */
"""
def get_antigen_length_bounds(antigen_file):
    lengths = []
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            seq = row.get("sequence", "").strip()
            if seq:
                lengths.append(len(seq))
    if not lengths:
        return None, None
    return min(lengths), max(lengths)

"""
/**
 * @brief Randomly fetches reviewed UniProt protein entries for a specific organism.
 *
 * Sends queries to the UniProt API to retrieve protein entries for a given organism,
 * filtering by sequence length and excluding names found in the known antigen set.
 *
 * @param n (int): Number of protein entries to fetch.
 * @param organism (str): Full organism name.
 * @param antigen_names (set): Set of antigen protein names to exclude.
 * @param min_len (int): Minimum protein sequence length.
 * @param max_len (int): Maximum protein sequence length.
 * @return (list): List of randomly sampled UniProt entry JSONs.
 */
"""
def fetch_random_uniprot_protein_entries(n=200, organism="Staphylococcus aureus", antigen_names=set(), min_len=None, max_len=None):
    headers = {"Accept": "application/json"}
    encoded_org = requests.utils.quote(organism)

    seq_length_param = f"&seqLength={min_len}-{max_len}" if min_len and max_len else ""
    base_url = f"https://www.ebi.ac.uk/proteins/api/proteins?reviewed=true&organism={encoded_org}{seq_length_param}"

    print(f"[INFO] Fetching record count for organism: {organism} with seqLength={min_len}-{max_len}")

    try:
        test_url = f"{base_url}&offset=0&size=1"
        response = requests.get(test_url, headers=headers)
        response.raise_for_status()
        max_records = int(response.headers.get("x-pagination-totalrecords", 0))
        if max_records == 0:
            print(f"[FATAL] No proteins found for organism with sequence length bounds.")
            return []
        print(f"[INFO] Found {max_records} total reviewed protein records within sequence range.")
    except Exception as e:
        print(f"[ERROR] Failed to get protein count: {e}")
        return []

    selected_entries = []
    tried_offsets = set()

    print(f"[INFO] Sampling {n} random reviewed UniProt protein entries...")
    while len(selected_entries) < n and len(tried_offsets) < max_records:
        offset = random.randint(0, max_records - 1)
        if offset in tried_offsets:
            continue
        tried_offsets.add(offset)

        try:
            url = f"{base_url}&offset={offset}&size=1"
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            entries = r.json()

            if not isinstance(entries, list):
                continue

            for entry in entries:
                name = entry.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "").lower()
                if name and name not in antigen_names:
                    selected_entries.append(entry)
                    print(f"  [INFO] Sampled: {name} ({len(selected_entries)}/{n})")

        except Exception as e:
            print(f"[WARN] Offset {offset} failed: {e}")
            continue

    if len(selected_entries) < n:
        print(f"[WARN] Only collected {len(selected_entries)} proteins out of {n} requested.")

    return selected_entries

"""
/**
 * @brief Main pipeline for generating non-antigen protein candidates for a given organism.
 *
 * Loads known antigen names and sequence bounds, fetches random non-antigen UniProt entries,
 * parses the entries using a shared parser, and writes the results to a CSV file.
 *
 * @param pathogen (str): Pathogen short name (used as subdirectory in `data/`).
 * @param organism (str): Full organism name.
 * @return: None
 */
"""
def main(pathogen, organism):
    organism_tag = organism.lower().replace(" ", "_")
    pathogen_dir = os.path.join("data", pathogen)
    antigens_file = os.path.join(pathogen_dir, f"{organism_tag}_compiled_proteins.csv")
    output_file = os.path.join(pathogen_dir, f"{organism_tag}_random_proteins.csv")

    if not os.path.exists(antigens_file):
        print("[FATAL] Missing antigen protein file.")
        return

    # Step 1: Get antigen names and length bounds
    antigen_names = get_antigen_protein_names(antigens_file)
    min_len, max_len = get_antigen_length_bounds(antigens_file)

    if not min_len or not max_len:
        print("[FATAL] Could not determine antigen sequence length bounds.")
        return

    # Step 2: Fetch random UniProt proteins
    entries = fetch_random_uniprot_protein_entries(
        n=200,
        organism=organism,
        antigen_names=antigen_names,
        min_len=min_len,
        max_len=max_len
    )

    # Step 3: Parse UniProt entries
    protein_data = []
    for entry in entries:
        protein_name = entry.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        parsed = parse_uniprot_entry(entry)
        if parsed:
            protein_data.append(parsed)

    if not protein_data:
        print("[WARN] No protein entries were successfully parsed.")
        return

    # Step 4: Write to CSV
    with open(output_file, "w", newline='') as outfile:
        fieldnames = [
            "uniprot_accession",
            "protein_name",
            "short_name",
            "function",
            "domains",
            "features",
            "sequence",
            "organism_name"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {len(protein_data)} protein sequences to {output_file}")

"""
/**
 * @brief Entry point of the script. Parses command-line arguments and starts the pipeline.
 *
 * Expects two arguments: <pathogen_subfolder> and <organism_name>.
 */
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch random reviewed UniProt protein entries for a given organism, excluding known antigens.",
        usage="python generate_random_proteins_per_strain.py <pathogen_directory> <pathogen_name>"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Full organism name (e.g., "staphylococcus aureus")')
    args = parser.parse_args()

    main(args.pathogen_directory, args.pathogen_name)
