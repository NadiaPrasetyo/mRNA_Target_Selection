import pandas as pd
import re
import csv
import unicodedata
import argparse

# === PARSE ARGS ===
parser = argparse.ArgumentParser()
parser.add_argument('--file_with_ids', required=True)
parser.add_argument('--file_new_antigens', required=True)
parser.add_argument('--tsv_output', required=True)
args = parser.parse_args()

# === CONFIGURATION ===
file_with_ids = args.file_with_ids
file_new_antigens = args.file_new_antigens
tsv_output = args.tsv_output

# === CLEANING FUNCTION ===
def clean_antigen_name(name):
    if not isinstance(name, str):
        return ""
    
    # Normalize to ASCII and lowercase
    name = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore').decode('ascii')
    name = name.lower()

    # Remove contents inside parentheses or brackets
    name = re.sub(r'\(.*?\)', '', name)
    name = re.sub(r'\[.*?\]', '', name)

    # Remove known Greek letter prefixes (e.g. alpha-, gamma-)
    name = re.sub(r'\b(alpha|beta|gamma|delta|epsilon|zeta|theta|kappa|lambda|mu|nu|xi|omicron|pi|rho|sigma|tau|upsilon|phi|chi|psi|omega)[ -]?', '', name, flags=re.IGNORECASE)

    # Remove after comma, slash, or semicolon
    name = re.split(r'[,/;]', name)[0]

    # Remove extra whitespace, non-word characters on ends
    name = re.sub(r'\s+', ' ', name)
    name = name.strip()
    name = re.sub(r'^\W+|\W+$', '', name)

    return name

# === LOAD AND CLEAN DATA ===
df_with_ids = pd.read_excel(file_with_ids)
df_new = pd.read_excel(file_new_antigens)

# Clean known antigens
known_antigens = df_with_ids.iloc[:, 0].dropna()
known_antigens_cleaned = known_antigens.str.replace(r"\s*\(UniProt:[^)]+\)", "", regex=True)
known_antigens_cleaned = known_antigens_cleaned.map(clean_antigen_name)

# Clean new antigens
new_antigens = df_new.iloc[:, 0].dropna().map(clean_antigen_name)

# === FIND UNIQUE NEW ANTIGENS ===
new_only = sorted(set(new_antigens) - set(known_antigens_cleaned))

# === GENERATE QUERIES ===
query_rows = []
for antigen in new_only:
    query1 = f'{antigen} [All Fields] NOT partial[All Fields] AND "Staphylococcus aureus"[Organism] AND swissprot[filter]'
    query2 = f'{antigen} [All Fields] NOT partial[All Fields] AND "Staphylococcus aureus"[Organism]'
    query_rows.append((antigen, query1, query2))

# === SAVE TO TSV ===
with open(tsv_output, "w", newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='\\')
    writer.writerow(["Antigen", "SwissProtQuery", "FallbackQuery"])
    for row in query_rows:
        writer.writerow(row)

print(f"✅ Wrote {len(query_rows)} cleaned antigen queries to {tsv_output}")
