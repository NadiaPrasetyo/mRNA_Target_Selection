import os
import sys
import csv
import random
import requests
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

def fetch_random_uniprot_protein_names(n=200, organism="Staphylococcus aureus", antigen_keywords=set()):
    headers = {"Accept": "application/json"}
    encoded_org = requests.utils.quote(organism)

    # Step 1: Get max number of records
    print(f"[INFO] Fetching record count for organism: {organism}")
    try:
        test_url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=1&organism={encoded_org}"
        response = requests.get(test_url, headers=headers)
        response.raise_for_status()
        max_records = int(response.headers.get("x-pagination-totalrecords", 0))
        if max_records == 0:
            print(f"[FATAL] No proteins found for organism: {organism}")
            return []
        print(f"[INFO] Found {max_records} total protein records for organism.")
    except Exception as e:
        print(f"[ERROR] Failed to get protein count: {e}")
        return []

    # Step 2: Sample random protein names
    selected_names = set()
    tried_offsets = set()

    print(f"[INFO] Sampling {n} random proteins (excluding antigens)...")
    while len(selected_names) < n and len(tried_offsets) < max_records:
        offset = random.randint(0, max_records - 1)
        if offset in tried_offsets:
            continue
        tried_offsets.add(offset)

        try:
            url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset={offset}&size=1&organism={encoded_org}"
            r = requests.get(url, headers=headers)
            r.raise_for_status()
            entries = r.json()

            if not isinstance(entries, list):
                continue

            for entry in entries:
                name = entry.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
                if name and not any(kw in name.lower() for kw in antigen_keywords):
                    selected_names.add(name)
                    print(f"  [INFO] Sampled: {name} ({len(selected_names)}/{n})")

        except Exception as e:
            print(f"[WARN] Offset {offset} failed: {e}")
            continue

    if len(selected_names) < n:
        print(f"[WARN] Only collected {len(selected_names)} proteins out of {n} requested.")

    return list(selected_names)

def main(pathogen, organism):
    pathogen_dir = os.path.join("data", pathogen)
    strains_file = os.path.join(pathogen_dir, f"{pathogen}_strains.csv")
    antigens_file = os.path.join(pathogen_dir, f"{pathogen}_compiled_antigens.csv")
    output_file = os.path.join(pathogen_dir, f"{pathogen}_random_proteins.csv")

    if not os.path.exists(strains_file) or not os.path.exists(antigens_file):
        print("[FATAL] Missing input files.")
        return

    antigen_keywords = load_antigen_keywords(antigens_file)
    random_proteins = fetch_random_uniprot_protein_names(n=200, organism=organism, antigen_keywords=antigen_keywords)

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
    if len(sys.argv) < 3:
        print("Usage: python generate_random_proteins_per_strain.py <pathogen_subfolder> <organism_name>")
    else:
        main(sys.argv[1], sys.argv[2])
