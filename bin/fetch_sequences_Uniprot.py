"""
/**
 * @file fetch_sequences_Uniprot.py
 * @brief Fetches protein sequence and metadata from UniProt based on antigen data.
 *
 * This script reads a compiled antigen CSV file containing antigen names, gene names,
 * and UniProt IDs. It queries the UniProt API to fetch full protein information
 * including sequence, function, domains, and features, then compiles and saves
 * the protein data into a new CSV file.
 *
 * General Flow:
 *   1. Load antigen records from a compiled CSV.
 *   2. For each antigen, fetch corresponding protein data via UniProt API.
 *   3. Parse and standardize protein metadata.
 *   4. Save compiled protein data to a CSV.
 *
 * Parameters:
 *   pathogen (str): Folder name for the pathogen (used for file paths).
 *   organism (str): Full name of the organism (used in UniProt queries).
 *
 * Usage:
 *   python fetch_sequences_Uniprot.py <pathogen> <organism>
 *
 * Example:
 *   python fetch_sequences_Uniprot.py sars_cov_2 "SARS-CoV-2"
 *
 * @author Nadia
 */
"""
import csv
import requests
import os
import re
import unicodedata
import argparse

UNIPROT_API_BASE = "https://www.ebi.ac.uk/proteins/api/proteins"

"""
/**
 * @brief Cleans and normalizes an antigen name.
 *
 * Removes special characters, normalizes Unicode to ASCII, strips Greek letters
 * and other extraneous tokens to return a clean, comparable name string.
 *
 * @param name (str): The antigen name to clean.
 * @return (str): The cleaned antigen name.
 */
"""
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

"""
/**
 * @brief Fetches UniProt data using a UniProt accession number.
 *
 * Queries the UniProt API for a given accession and organism. Returns
 * the first matching protein entry if found.
 *
 * @param accession (str): UniProt accession number.
 * @param organism (str): Organism name for the query.
 * @return (dict or None): JSON response from UniProt or None if failed.
 */
"""
def fetch_uniprot_by_accession(accession, organism):
    url = f"{UNIPROT_API_BASE}?offset=0&size=-1&accession={accession}&organism={organism.replace(' ', '%20')}"
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        if response.ok:
            data = response.json()
            if isinstance(data, list) and data:
                return data[0]
            elif isinstance(data, dict):
                return data
    except requests.RequestException as e:
        print(f"[ERROR] Request failed for accession {accession}: {e}")
    return None

"""
/**
 * @brief Fetches UniProt data using a gene name.
 *
 * Queries the UniProt API by gene name and organism. Returns the first result if available.
 *
 * @param gene_name (str): Gene name to query.
 * @param organism (str): Organism name for the query.
 * @return (dict or None): JSON response from UniProt or None if no match found.
 */
"""
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

"""
/**
 * @brief Parses a UniProt protein entry into a structured dictionary.
 *
 * Extracts relevant metadata fields such as function, domains, features, and sequence.
 *
 * @param entry (dict): Raw JSON entry from UniProt API.
 * @return (dict or None): Structured dictionary of protein data or None if parsing fails.
 */
"""
def parse_uniprot_entry(entry):
    try:
        accession = entry.get("accession", "")
        organism = entry.get("organism", {}).get("names", [{}])[0].get("value", "").lower()
        name = entry.get("id", "")

        protein_data = entry.get("protein", {})
        protein_name = ""
        if "recommendedName" in protein_data:
            protein_name = protein_data["recommendedName"].get("fullName", {}).get("value", "")
        elif "submittedName" in protein_data:
            protein_name = protein_data["submittedName"][0].get("fullName", {}).get("value", "")

        function = ""
        for comment in entry.get("comments", []):
            if comment.get("type") == "FUNCTION":
                function = comment.get("text", [{}])[0].get("value", "")
                break

        domains = [
            ref.get("properties", {}).get("entry name", "")
            for ref in entry.get("dbReferences", [])
            if ref["type"] in ("InterPro", "Pfam")
        ]

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

"""
/**
 * @brief Loads antigen data from a CSV file.
 *
 * Reads the compiled antigen file and standardizes fields for downstream processing.
 *
 * @param file_path (str): Path to the antigen CSV file.
 * @return (list): List of antigen dictionaries with cleaned fields.
 */
"""
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

"""
/**
 * @brief Main function for fetching and compiling UniProt protein data.
 *
 * Reads antigens, fetches UniProt entries by ID or gene, extracts relevant metadata,
 * and writes the compiled protein data to a CSV file.
 *
 * @param pathogen (str): Folder name of the pathogen.
 * @param organism (str): Full name of the organism to search against.
 */
"""
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

"""
/**
 * @brief Entry point for the script. Parses CLI arguments and calls main().
 *
 * Usage:
 *   python fetch_sequences_Uniprot.py <pathogen_subfolder> <organism_name>
 */
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch protein sequence and metadata from UniProt based on antigen data.",
        usage="python fetch_sequences_Uniprot.py <pathogen_directory> <pathogen_name>"
    )
    parser.add_argument("pathogen_directory", help="Directory name under data/")
    parser.add_argument("pathogen_name", help='Prefix used in filenames (e.g., "staphylococcus aureus")')
    args = parser.parse_args()

    main(args.pathogen_directory, args.pathogen_name)
