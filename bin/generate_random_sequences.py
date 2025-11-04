#!/usr/bin/env python3
"""
generate_random_sequences.py
Command-line tool to fetch random reviewed UniProt protein entries for a given organism, excluding known antigens.

Overview:
    - Loads known antigen protein names and their sequence length bounds from a compiled CSV file.
    - Queries the UniProt API to randomly sample reviewed protein entries matching the organism and sequence criteria.
    - Excludes any protein names that match known antigens or have already been processed.
    - Parses UniProt entries using a shared parsing function.
    - Optionally trims human protein sequences to match antigen length bounds.
    - Saves the parsed protein information to a new CSV file for further analysis.

Arguments:
    pathogen_directory (str): Subdirectory under `data/` containing pathogen data.
    pathogen_name (str): Full organism name (used for querying UniProt and naming output files).
    --human (optional): If specified, fetches and processes human proteins instead of pathogen proteins.

Requirements:
    - Compiled antigen protein CSV file present in the pathogen data directory.
    - Python packages: argparse, csv, os, random, requests.

Usage Example:
    python generate_random_sequences.py s_aureus "Staphylococcus aureus"
    python generate_random_sequences.py human "Homo sapiens" --human

Outputs:
    data/<pathogen_directory>/random_compiled_proteins.csv   # Randomly sampled non-antigen protein entries
    data/<pathogen_directory>/human_compiled_proteins.csv    # Trimmed human protein entries (if --human is used)

Author: Nadia
"""

import os
import csv
import random
import requests
import argparse
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))  # Add parent directory to sys.path
from bin.fetch_sequences_Uniprot import parse_uniprot_response, fetch_refseq_nucleotide

# Increase the CSV field size limit to handle very large sequences
csv.field_size_limit(sys.maxsize)

def get_antigen_protein_names(antigen_file):
    """Extract antigen protein names from the CSV file.
    Args:
        antigen_file (str): Path to the CSV file containing antigen data.
    """
    names = set()
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row.get("protein_name", "").strip().lower()
            if name:
                names.add(name)
    return names


def get_antigen_length_bounds(antigen_file):
    """Calculate the minimum and maximum lengths of antigen sequences from the CSV file.
    Args:
        antigen_file (str): Path to the CSV file containing antigen data.
    Returns:
        tuple: (min_length, max_length) of antigen sequences, or (None, None) if no valid sequences found.
    """
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


def trim_human_proteins_to_length(protein_data, min_length, max_length):
    """Randomly trim human protein sequences (and nucleotide sequences) to fit within given bounds."""
    for protein in protein_data:
        seq = protein.get("sequence", "").strip()
        nuc_seq = protein.get("nucleotide_sequence", "").strip()

        if not seq or not nuc_seq:
            continue

        seq_len = len(seq)
        if seq_len < min_length:
            print(f"Skipping protein {protein.get('id')} - length {seq_len} < min_length {min_length}")
            continue

        target_len = random.randint(min_length, min(max_length, seq_len))
        max_offset = seq_len - target_len
        offset = random.randint(0, max_offset) if max_offset > 0 else 0

        protein["sequence"] = seq[offset:offset + target_len]
        protein["nucleotide_sequence"] = nuc_seq[offset * 3: offset * 3 + target_len * 3]

    return protein_data

def get_entry_name(entry):
    """Extract a usable protein name from UniProt entry.
    Args:
        entry (dict): UniProt entry as returned by the API.
    Returns:
        str: Protein name in lowercase, or primary accession if no name found.
    """
    protein_desc = entry.get("proteinDescription", {})

    # Recommended name
    rec_name = protein_desc.get("recommendedName", {})
    full_name = rec_name.get("fullName", {}).get("value", "")
    if full_name:
        return full_name.lower()

    # Alternative names
    for alt in protein_desc.get("alternativeNames", []):
        full_name = alt.get("fullName", {}).get("value", "")
        if full_name:
            return full_name.lower()

    # Fallback
    return entry.get("primaryAccession", "").lower()


# Extract next link from headers
def get_next_link(headers):
    """
    Extract the 'next' link from the response headers for pagination.
    Args:
        headers (dict): Response headers from the UniProt API.
    Returns:
        str or None: URL for the next page of results, or None if not present.
    """
    link_header = headers.get("link")
    if not link_header:
        return None
    # Example format: <url>; rel="next"
    parts = link_header.split(";")
    if len(parts) >= 2 and 'rel="next"' in parts[1]:
        url = parts[0].strip()[1:-1]  # remove < and >
        return url
    return None

def fetch_random_uniprot_protein_entries(n=200, organism="Staphylococcus aureus", antigen_names=set(), min_len=None, max_len=None):
    """
    Fetch random reviewed UniProt protein entries for a given organism, excluding known antigens.
    Uses the new UniProt REST API (https://rest.uniprot.org/uniprotkb/search).
    Args:
        n (int): Number of random protein entries to fetch.
        organism (str): Organism name for the UniProt query.
        antigen_names (set): Set of antigen protein names to exclude.
        min_len (int or None): Minimum sequence length to filter proteins (inclusive).
        max_len (int or None): Maximum sequence length to filter proteins (inclusive).
    Returns:
        dict: JSON response from UniProt API containing the sampled protein entries.
    """
    headers = {"Accept": "application/json"}
    query_parts = [f'organism_name:"{organism}"', "reviewed:true"]
    if min_len and max_len:
        query_parts.append(f"length:[{min_len} TO {max_len}]")
    query = " AND ".join(query_parts)

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "size": 500, # Initial request to get total count
        "query": query,
        "fields": "accession,protein_name,organism_name,length,sequence,xref_pfam,xref_refseq",
    }

    print(f"[INFO] Querying UniProt: {query}")
    try:
        response = requests.get(base_url, headers=headers, params=params)
        response.raise_for_status()
        results = response.json()
        # get the response header to get total results
        response_headers = response.headers
        total_records = int(response_headers.get("x-total-results", "0"))
        if total_records == 0:
            print("[FATAL] No proteins found for organism with sequence length bounds.")
            return {"results": []}
        print(f"[INFO] Found {total_records} reviewed protein records within sequence range.")
    except Exception as e:
        print(f"[ERROR] Failed to query UniProt API: {e}")
        return {"results": []}

    # Collect initial page of results
    all_entries = results.get("results", [])
    cursor = get_next_link(response_headers)
    
    # Paginate until no more pages
    while cursor:
        try:
            r = requests.get(cursor, headers=headers)
            r.raise_for_status()
            data = r.json()
            all_entries.extend(data.get("results", []))
            cursor = get_next_link(r.headers)
        except Exception as e:
            print(f"[WARN] Pagination fetch failed: {e}")
            break

    print(f"[INFO] Collected {len(all_entries)} total protein entries from UniProt.")
    
    # Filter and deduplicate
    filtered_entries = []
    seen_names = set()

    for entry in all_entries:
        name = get_entry_name(entry)
        if not name or name in antigen_names or name in seen_names:
            continue
        seen_names.add(name)
        filtered_entries.append(entry)

    print(f"[INFO] {len(filtered_entries)} entries remain after antigen filtering and deduplication.")

    # Random sample from the filtered pool
    if not filtered_entries:
        print("[FATAL] No eligible proteins left after filtering.")
        return {"results": []}

    sampled = random.sample(filtered_entries, k=min(n, len(filtered_entries)))

    return {"results": sampled}


def main(pathogen, organism, include_human=False):
    """Main function to generate non-antigen protein candidates.
    Args:
        pathogen (str): Directory name under data/ for the pathogen.
        organism (str): Full organism name (e.g., "Staphylococcus aureus").
        include_human (bool): If True, fetch human proteins instead of pathogen proteins.
    """
    organism_tag = organism.lower().replace(" ", "_")
    pathogen_dir = os.path.join("data", pathogen)
    antigens_file = os.path.join(pathogen_dir, f"{organism_tag}_compiled_proteins.csv")
    output_file = os.path.join(pathogen_dir, "human_compiled_proteins.csv" if include_human else "random_compiled_proteins.csv")

    if not os.path.exists(antigens_file):
        print("[FATAL] Missing antigen protein file.")
        return

    antigen_names = get_antigen_protein_names(antigens_file)
    min_len, max_len = get_antigen_length_bounds(antigens_file)
    if not min_len or not max_len:
        print("[FATAL] Could not determine antigen sequence length bounds.")
        return

    entries = fetch_random_uniprot_protein_entries(
        n=200,
        organism="Homo sapiens" if include_human else organism,
        antigen_names=antigen_names,
        min_len=min_len if not include_human else None,
        max_len=max_len if not include_human else None
    )

    protein_data = parse_uniprot_response(entries)

    if not protein_data:
        print("[WARN] No protein entries were successfully parsed.")
        return

    if include_human:
        protein_data = trim_human_proteins_to_length(protein_data, min_len, max_len)

    with open(output_file, "w", newline='') as outfile:
        fieldnames = [
            "uniprot_accession",
            "protein_name",
            "short_name",
            "function",
            "domains",
            "features",
            "sequence",
            "organism_name",
            "pfam",
            "nucleotide_sequence"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {len(protein_data)} protein sequences to {output_file}")


if __name__ == "__main__":
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Fetch random UniProt protein entries for a given organism, excluding known antigens.",
        usage="python generate_random_sequences.py <pathogen_directory> <pathogen_name>"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Full organism name (e.g., "staphylococcus aureus")')
    parser.add_argument("--human", action="store_true", help="Take human proteins")
    args = parser.parse_args()

    main(args.pathogen_directory, args.pathogen_name, include_human=args.human)
