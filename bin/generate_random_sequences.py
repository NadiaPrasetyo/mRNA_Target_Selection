import csv
import os
import re
import sys
import random
import string
import requests
from io import StringIO
from Bio import SeqIO
import subprocess

from fetch_protein_data import (
    load_antigen_keywords,
    protein_matches,
    parse_protein_entry,
    fetch_protein_data,
    fetch_protein_data_ncbi,
    parse_genpept_entries
)

def random_dna_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def get_antigen_lengths(antigen_file):
    lengths = []
    with open(antigen_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            seq = row.get("sequence", "")
            if seq:
                lengths.append(len(seq))
    return lengths

def generate_non_antigens(pathogen, target_count=200):
    pathogen_dir = os.path.join("data", pathogen)
    strains_file = os.path.join(pathogen_dir, f"{pathogen}_strains.csv")
    antigens_file = os.path.join(pathogen_dir, f"{pathogen}_compiled_proteins.csv")
    output_file = os.path.join(pathogen_dir, f"{pathogen}_non_antigens.csv")

    if not os.path.exists(strains_file) or not os.path.exists(antigens_file):
        print("[FATAL] Required files missing.")
        return

    antigen_keywords = load_antigen_keywords(antigens_file)
    antigen_lengths = get_antigen_lengths(antigens_file)
    length_range = (min(antigen_lengths), max(antigen_lengths)) if antigen_lengths else (150, 450)

    collected = []
    strain_count = 0

    with open(strains_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if len(collected) >= target_count:
                break

            strain = row.get("Strain", "")
            embl_id = row.get("EMBL_ID", "")
            if not embl_id:
                continue

            print(f"[INFO] Trying strain {strain} (EMBL: {embl_id})")

            entries = fetch_protein_data(embl_id)
            if not entries:
                print(f"[INFO] Falling back to NCBI for {strain}")
                dummy_antigen = "hypothetical protein"  # Generic fallback
                gp_text = fetch_protein_data_ncbi(strain, dummy_antigen)
                ncbi_entries = parse_genpept_entries(gp_text, strain, set())

                for entry in ncbi_entries:
                    if len(collected) >= target_count:
                        break
                    if not protein_matches(entry["protein_name"], entry["short_name"], antigen_keywords):
                        collected.append(entry)
                continue

            for entry in entries:
                if len(collected) >= target_count:
                    break
                parsed = parse_protein_entry(entry, strain, set())  # Don't use keywords here
                if parsed and not protein_matches(parsed["protein_name"], parsed["short_name"], antigen_keywords):
                    collected.append(parsed)

            strain_count += 1

    # Fill with random sequences if we’re short
    while len(collected) < target_count:
        rand_len = random.randint(*length_range)
        fake_seq = random_dna_sequence(rand_len)
        collected.append({
            "strain": "synthetic",
            "uniprot_accession": f"SYN{len(collected)+1}",
            "protein_name": "Random Sequence",
            "short_name": f"random_seq_{len(collected)+1}",
            "function": "",
            "domains": "",
            "sequence": fake_seq
        })

    print(f"[DONE] Collected {len(collected)} non-antigen entries")

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ["strain", "uniprot_accession", "protein_name", "short_name", "function", "domains", "sequence"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(collected)

    print(f"[OUTPUT] Wrote non-antigens to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_non_antigens.py <pathogen_subfolder>")
    else:
        generate_non_antigens(sys.argv[1])
