#!/usr/bin/env python3
import pandas as pd
import re
import argparse
import os

# === ARGUMENT PARSING ===
parser = argparse.ArgumentParser(description="Extract UniProt IDs from a CSV or Excel file.")
parser.add_argument('--input', '-i', required=True, help="Path to input file (.csv or .xlsx)")
parser.add_argument('--output', '-o', required=True, help="Path to output UniProt ID text file")
args = parser.parse_args()

input_file = args.input
output_txt = args.output

# === DETERMINE FILE TYPE ===
file_ext = os.path.splitext(input_file)[-1].lower()

if file_ext == '.csv':
    df = pd.read_csv(input_file)
elif file_ext in ['.xlsx', '.xls']:
    df = pd.read_excel(input_file)
else:
    raise ValueError(f"Unsupported file type: {file_ext}. Use .csv or .xlsx")

# === EXTRACT FIRST COLUMN ===
antigen_iris = df.iloc[:, 0].dropna().astype(str).str.strip()

# === PROCESS UNIProt IDs ===
uniprot_ids = []
missing_uniprot = []

for iri in antigen_iris:
    match = re.search(r'uniprot:([A-Z0-9]+)', iri, re.IGNORECASE)
    if match:
        uniprot_ids.append(match.group(1))
    else:
        missing_uniprot.append(iri)

# === WRITE TO OUTPUT TEXT FILE ===
with open(output_txt, "w") as f:
    for uid in sorted(set(uniprot_ids)):
        f.write(uid + "\n")

    if missing_uniprot:
        f.write("\n# Antigens without UniProt IDs\n")
        for item in missing_uniprot:
            f.write("# " + item + "\n")

print(f"[INFO] Wrote {len(uniprot_ids)} UniProt IDs to '{output_txt}'")
