import os
import sys
import csv
import json
import random
import requests
import subprocess
import re
from io import StringIO
from Bio import SeqIO
from fetch_sequences_Uniprot_NCBI import (
    clean_antigen_name,
    load_antigen_keywords,
    protein_matches,
    parse_protein_entry,
    fetch_protein_data,
    fetch_protein_data_ncbi,
    parse_genpept_entries,
    extract_keywords
)

UNIPROT_URL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&organism=Staphylococcus%20aureus"

def fetch_all_uniprot_proteins():
    print("[INFO] Fetching protein names from UniProt...")
    headers = {"Accept": "application/json"}
    try:
        r = requests.get(UNIPROT_URL, headers=headers)
        r.raise_for_status()
        proteins = r.json()
        names = []
        for entry in proteins:
            name = entry.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
            if name:
                names.append(name)
        print(f"[INFO] Retrieved {len(names)} protein names from UniProt")
        return names
    except Exception as e:
        print(f"[ERROR] Failed to fetch from UniProt: {e}")
        return []

def get_random_non_antigen_proteins(all_proteins, antigen_keywords, n=200):
    clean_set = []
    for name in all_proteins:
        if not any(kw in name.lower() for kw in antigen_keywords):
            clean_set.append(name)
    random.shuffle(clean_set)
    selected = clean_set[:n]
    print(f"[INFO] Selected {len(selected)} random non-antigen proteins")
    return selected

def main(pathogen):
    pathogen_dir = os.path.join("data", pathogen)
    strains_file = os.path.join(pathogen_dir, f"{pathogen}_strains.csv")
    antigens_file = os.path.join(pathogen_dir, f"{pathogen}_compiled_antigens.csv")
    output_file = os.path.join(pathogen_dir, f"{pathogen}_random_proteins.csv")

    if not os.path.exists(strains_file) or not os.path.exists(antigens_file):
        print("[FATAL] Missing input files.")
        return

    antigen_keywords = load_antigen_keywords(antigens_file)
    all_protein_names = fetch_all_uniprot_proteins()
    random_proteins = get_random_non_antigen_proteins(all_protein_names, antigen_keywords)

    with open(strains_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        strain_list = list(reader)

    protein_data = []
    total_matched = 0

    for protein_name in random_proteins:
        keywords = {kw.lower() for kw in extract_keywords(protein_name).split()}
        print(f"[INFO] Processing protein: {protein_name}")

        for strain_row in strain_list:
            strain = strain_row.get("Strain", "")
            embl_id = strain_row.get("EMBL_ID", "")

            entries = []
            parsed_entries = []

            if embl_id:
                entries = fetch_protein_data(embl_id)
                for entry in entries:
                    parsed = parse_protein_entry(entry, strain, keywords)
                    if parsed:
                        parsed["random_protein_name"] = protein_name
                        parsed_entries.append(parsed)

            if not parsed_entries:
                gp_text = fetch_protein_data_ncbi(strain, protein_name)
                if gp_text:
                    ncbi_parsed = parse_genpept_entries(gp_text, strain, keywords)
                    for item in ncbi_parsed:
                        item["random_protein_name"] = protein_name
                        parsed_entries.append(item)

            protein_data.extend(parsed_entries)
            total_matched += len(parsed_entries)

        print(f"  [INFO] Matched {total_matched} entries so far.")

    if not protein_data:
        print("[WARN] No matches found.")
        return

    with open(output_file, "w", newline='') as outfile:
        fieldnames = ["random_protein_name", "strain", "uniprot_accession", "protein_name", "short_name", "function", "domains", "sequence"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {total_matched} matched sequences to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_random_proteins_per_strain.py <pathogen_subfolder>")
    else:
        main(sys.argv[1])
