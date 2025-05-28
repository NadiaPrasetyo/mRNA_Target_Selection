import csv
import json
import requests
import sys
import os
import re
import subprocess
import xml.etree.ElementTree as ET

def fetch_protein_data_ncbi(strain_name, antigen_name):
    strain_full = f'Staphylococcus aureus subsp. aureus {strain_name}'
    antigen_escaped = antigen_name.replace('"', '')  # Avoid breaking inner quotes

    # Full query string with properly quoted parts
    query = f'"{strain_full}"[All Fields] AND ("{antigen_escaped}"[Protein Name] OR "{antigen_escaped}"[All Fields])'

    try:
        # Wrap entire query in single quotes so shell sees it as one argument
        cmd = f'esearch -db protein -query \'{query}\' | efetch -format xml'
        output = subprocess.check_output(cmd, shell=True, text=True)
        root = ET.fromstring(output)
        return root.findall(".//GBSeq")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] NCBI fetch failed for strain: {strain_name}, antigen: {antigen_name}\n{e}")
        return []


def fetch_protein_data(embl_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/EMBL:{embl_id}?offset=0&size=-1&reviewed=true"
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"[ERROR] API failed for EMBL ID {embl_id}: HTTP {response.status_code}")
            return []
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Network exception for EMBL ID {embl_id}: {e}")
        return []

def extract_keywords(antigen_name, max_words=3):
    clean = re.sub(r"\(.*?\)", "", antigen_name)
    clean = re.sub(r"\[|\]|UniProt:[A-Z0-9]+", "", clean)
    clean = clean.replace("'", "").strip()
    words = clean.split()
    return " ".join(words[:max_words]).lower()

def load_antigen_keywords(antigen_file):
    keywords = set()
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            raw_name = row.get("antigen_name", "")
            if raw_name:
                keywords.add(extract_keywords(raw_name))
    print(f"[INFO] Loaded {len(keywords)} antigen keyword patterns (UniProt)")
    return keywords

def load_antigen_names(antigen_file):
    names = []
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row.get("antigen_name", "").strip()
            if name:
                names.append(name)
    print(f"[INFO] Loaded {len(names)} antigen names (NCBI)")
    return names

def protein_matches(protein_name, short_name, keywords):
    haystack = f"{protein_name} {short_name}".lower()
    return any(kw in haystack for kw in keywords)

def parse_protein_entry(entry, strain, keywords):
    try:
        accession = entry.get("accession", "")
        name = entry.get("id", "")
        protein_name = entry.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        short_name_joined = entry.get("id", "")

        if not protein_matches(protein_name, short_name_joined, keywords):
            return None

        function = ""
        for comment in entry.get("comments", []):
            if comment.get("type") == "FUNCTION":
                function = comment.get("text", [{}])[0].get("value", "")
                break

        domains = []
        for ref in entry.get("dbReferences", []):
            if ref["type"] in ("InterPro", "Pfam"):
                domain_name = ref.get("properties", {}).get("entry name", "")
                if domain_name:
                    domains.append(domain_name)

        features = []
        for f in entry.get("features", []):
            if f.get("type") == "VARIANT":
                original = f.get("original", "")
                variations = f.get("variation", [])
                features.append(f"{original}:{','.join(variations)}")
            elif f.get("type") == "PEPTIDE":
                desc = f.get("description", "")
                begin = f.get("begin", "")
                end = f.get("end", "")
                features.append(f"{desc} ({begin}-{end})")

        sequence = entry.get("sequence", {}).get("sequence", "")

        return {
            "strain": strain,
            "uniprot_accession": accession,
            "protein_name": protein_name,
            "short_name": short_name_joined,
            "function": function,
            "domains": ";".join(domains),
            "features": ";".join(features),
            "sequence": sequence
        }

    except Exception as e:
        print(f"[ERROR] Failed to parse protein entry for strain {strain} ({entry.get('accession', 'N/A')}): {e}")
        return None

def parse_ncbi_protein_entry(entry, strain, keywords):
    try:
        accession = entry.findtext("GBSeq_primary-accession", "")
        protein_name = entry.findtext("GBSeq_definition", "")
        sequence = entry.findtext("GBSeq_sequence", "").upper()
        short_name_joined = entry.findtext("GBSeq_locus", "")

        if not protein_matches(protein_name, short_name_joined, keywords):
            return None

        features = []
        for feat in entry.findall(".//GBFeature"):
            key = feat.findtext("GBFeature_key", "")
            if key == "Region":
                desc = ""
                for qual in feat.findall("GBFeature_quals/GBQualifier"):
                    if qual.findtext("GBQualifier_name") == "note":
                        desc = qual.findtext("GBQualifier_value")
                loc = feat.findtext("GBFeature_location", "")
                features.append(f"{desc} ({loc})")

        return {
            "strain": strain,
            "uniprot_accession": accession,
            "protein_name": protein_name,
            "short_name": short_name_joined,
            "function": "",
            "domains": "",
            "features": ";".join(features),
            "sequence": sequence
        }

    except Exception as e:
        print(f"[ERROR] Failed to parse NCBI protein entry for strain {strain}: {e}")
        return None

def main(pathogen):
    pathogen_dir = os.path.join("data", pathogen)
    strains_file = os.path.join(pathogen_dir, f"{pathogen}_strains.csv")
    antigens_file = os.path.join(pathogen_dir, f"{pathogen}_compiled_antigens.csv")
    output_file = os.path.join(pathogen_dir, f"{pathogen}_compiled_proteins.csv")

    if not os.path.exists(strains_file):
        print(f"[FATAL] Strains file not found: {strains_file}")
        return
    if not os.path.exists(antigens_file):
        print(f"[FATAL] Antigen file not found: {antigens_file}")
        return

    antigen_keywords = load_antigen_keywords(antigens_file)  # for UniProt
    antigen_names = load_antigen_names(antigens_file)        # for NCBI

    protein_data = []
    strain_count = 0
    matched_count = 0

    with open(strains_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            strain_count += 1
            strain = row.get("Strain", "")
            embl_id = row.get("EMBL_ID", "")

            if not embl_id:
                print(f"[WARN] Missing EMBL ID for strain: {strain}")
                continue

            print(f"[INFO] Processing strain {strain} (EMBL: {embl_id}) [{strain_count}]")

            entries = fetch_protein_data(embl_id)
            if not entries:
                print(f"[WARN] No UniProt entries found for EMBL ID {embl_id}, trying NCBI...")
                for antigen in antigen_names:
                    ncbi_entries = fetch_protein_data_ncbi(strain, antigen)
                    for entry in ncbi_entries:
                        parsed = parse_ncbi_protein_entry(entry, strain, {antigen.lower()})
                        if parsed:
                            protein_data.append(parsed)
                            matched_count += 1

            print(f"  [INFO] Found {len(entries)} protein entries")

            for entry in entries:
                parsed = parse_protein_entry(entry, strain, antigen_keywords)
                if parsed:
                    protein_data.append(parsed)
                    matched_count += 1

            print(f"  [INFO] Matched {matched_count} proteins so far")

    if not protein_data:
        print("[WARN] No matching proteins found. CSV not written.")
        return

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ["strain", "uniprot_accession", "protein_name", "short_name", "function", "domains", "features", "sequence"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(protein_data)

    print(f"[DONE] Wrote {matched_count} proteins to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fetch_protein_data.py <pathogen_subfolder>")
    else:
        main(sys.argv[1])
