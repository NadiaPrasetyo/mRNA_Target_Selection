import csv
import requests
import os
import re
import unicodedata
import sys

UNIPROT_API_BASE = "https://www.ebi.ac.uk/proteins/api/proteins"

def clean_antigen_name(name):
    if not isinstance(name, str):
        return ""
    name = name.strip()
    if name.startswith("[") and name.endswith("]"):
        name = name[1:-1].strip()
    if name.startswith("'") and name.endswith("'"):
        name = name[1:-1].strip()

    name = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore').decode('ascii')
    name = name.lower()
    name = re.sub(r'\(.*?\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)
    name = re.sub(r'\b(alpha|beta|gamma|delta|epsilon|zeta|theta|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega)[ -]?', '', name, flags=re.IGNORECASE)
    name = re.split(r'[,/;]', name)[0]
    name = re.sub(r'\s+', ' ', name).strip()
    name = re.sub(r'^\W+|\W+$', '', name)
    return name


def fetch_uniprot_by_accession(accession, organism):
    url = f"{UNIPROT_API_BASE}?offset=0&size=-1&accession={accession}&organism={organism.replace(' ', '%20')}"
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        if response.ok:
            data = response.json()
            if isinstance(data, list) and data:
                return data[0]  # Return the first result
            elif isinstance(data, dict):  # Just in case a single dict is returned
                return data
    except requests.RequestException as e:
        print(f"[ERROR] Request failed for accession {accession}: {e}")
    return None


def fetch_uniprot_by_gene(gene_name, organism):
    headers = {"Accept": "application/json"}
    query = f"{UNIPROT_API_BASE}?offset=0&size=-1&gene={gene_name}&organism={organism}"
    try:
        response = requests.get(query, headers=headers)
        if response.ok:
            data = response.json()
            return data[0] if data else None
    except requests.RequestException:
        pass
    return None

def parse_uniprot_entry(entry):
    try:
        accession = entry.get("accession", "")
        organism = entry.get("organism", {}).get("names", [{}])[0].get("value", "").lower()
        name = entry.get("id", "")

        # Handle protein name (fallback to submittedName if recommendedName is missing)
        protein_data = entry.get("protein", {})
        protein_name = ""
        if "recommendedName" in protein_data:
            protein_name = protein_data["recommendedName"].get("fullName", {}).get("value", "")
        elif "submittedName" in protein_data:
            protein_name = protein_data["submittedName"][0].get("fullName", {}).get("value", "")

        # Extract function comment
        function = ""
        for comment in entry.get("comments", []):
            if comment.get("type") == "FUNCTION":
                function = comment.get("text", [{}])[0].get("value", "")
                break

        # Extract domains
        domains = [
            ref.get("properties", {}).get("entry name", "")
            for ref in entry.get("dbReferences", [])
            if ref["type"] in ("InterPro", "Pfam")
        ]

        # Extract features: include type, description, and location
        features = []
        for feat in entry.get("features", []):
            feat_type = feat.get("type", "")
            desc = feat.get("description", "")
            begin = feat.get("begin", "")
            end = feat.get("end", "")
            if begin and end:
                features.append(f"{feat_type}:{desc}({begin}-{end})")

        sequence = entry.get("sequence", {}).get("sequence", "")

        return {
            "uniprot_accession": accession,
            "organism_name": organism,
            "protein_name": protein_name,
            "short_name": name,
            "function": function,
            "domains": ";".join(domains),
            "features": ";".join(features),
            "sequence": sequence
        }
    except Exception as e:
        print(f"[ERROR] Failed to parse entry for {entry.get('accession', 'unknown')}: {e}")
        return None


def load_antigen_records(file_path):
    records = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            antigen = {
                "antigen_name": clean_antigen_name(row.get("antigen_name", "")),
                "gene_name": row.get("gene_name", "").strip(),
                "uniprot_id": row.get("Uniprot_ID", "").strip()
            }
            records.append(antigen)
    return records

def main(pathogen, organism):
    pathogen_dir = os.path.join("data", pathogen)
    organism_tag = organism.replace(" ", "_").lower()
    antigens_file = os.path.join(pathogen_dir, f"{organism_tag}_compiled_antigens.csv")
    output_file = os.path.join(pathogen_dir, f"{organism_tag}_compiled_proteins.csv")

    if not os.path.exists(antigens_file):
        print(f"[FATAL] Antigen file not found: {antigens_file}")
        return

    antigen_records = load_antigen_records(antigens_file)
    protein_data = []
    seen_accessions = set()

    for record in antigen_records:
        antigen_name = record["antigen_name"]
        gene_name = record["gene_name"]
        uniprot_id = record["uniprot_id"]

        print(f"[INFO] Processing antigen: {antigen_name}")

        entry = None
        used_method = None

        if uniprot_id:
            print(f"  [DEBUG] Attempting fetch by accession: {uniprot_id}")
            entry = fetch_uniprot_by_accession(uniprot_id, organism)
            used_method = "accession"
        elif gene_name:
            print(f"  [DEBUG] Attempting fetch by gene: {gene_name}")
            entry = fetch_uniprot_by_gene(gene_name, organism)
            used_method = "gene"

        if entry:
            parsed = parse_uniprot_entry(entry)

            if parsed and parsed["uniprot_accession"] not in seen_accessions:
                protein_data.append(parsed)
                seen_accessions.add(parsed["uniprot_accession"])
                print(f"  [✓] {used_method.upper()} match: {parsed['uniprot_accession']}")
            else:
                print("  [WARN] Entry found but could not be parsed or is duplicate")
        else:
            print("  [WARN] No entry found for antigen")


    if not protein_data:
        print("[WARN] No matching proteins found.")
        return

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ["uniprot_accession", "protein_name", "short_name", "function", "domains", "features", "sequence", "organism_name"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {len(protein_data)} proteins to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fetch_protein_data.py <pathogen_subfolder> <organism_name>")
    else:
        main(sys.argv[1], sys.argv[2])
