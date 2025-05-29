import os
import sys
import csv
import random
import requests
from fetch_sequences_Uniprot_NCBI import parse_uniprot_entry  # ✅ using shared parser

def get_antigen_protein_names(antigen_file):
    names = set()
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row.get("protein_name", "").strip().lower()
            if name:
                names.add(name)
    return names

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

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python generate_random_proteins_per_strain.py <pathogen_subfolder> <organism_name>")
    else:
        main(sys.argv[1], sys.argv[2])
