import pandas as pd
import re
import csv

# === CONFIGURATION ===
file_with_ids = "data/S.aureus/IEDB_S.aureus_antigen_table.xlsx"     # First file with UniProt IDs
file_new_antigens = "data/S.aureus/Howden2023_S.aureus_virulence-factor-list.xlsx" # Second file with just names
query_csv = "data/S.aureus/entrez_queries.csv"


# === CLEANING FUNCTION ===
def clean_antigen_name(name):
    # Remove parentheses and their contents, and truncate at commas or slashes
    name = re.sub(r'\(.*?\)', '', name)
    name = re.split(r'[,/]', name)[0]
    return name.strip()

# === LOAD AND CLEAN DATA ===
df_with_ids = pd.read_excel(file_with_ids)
df_new = pd.read_excel(file_new_antigens)

known_antigens = df_with_ids.iloc[:, 0].dropna()
known_antigens_cleaned = known_antigens.str.replace(r"\s*\(UniProt:[^)]+\)", "", regex=True).map(clean_antigen_name)

new_antigens = df_new.iloc[:, 0].dropna().map(clean_antigen_name)

# === FIND UNIQUE TO NEW LIST ===
new_only = sorted(set(new_antigens) - set(known_antigens_cleaned))

with open(query_csv, mode="w", newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(["Antigen", "SwissProtQuery", "FallbackQuery"])
    for antigen in new_only:
        swissprot = f'{antigen} [All Fields] NOT partial[All Fields] AND "Staphylococcus aureus"[Organism] AND swissprot[filter]'
        fallback = f'{antigen} [All Fields] NOT partial[All Fields] AND "Staphylococcus aureus"[Organism]'
        writer.writerow([antigen, swissprot, fallback])

print(f"✅ Wrote {len(new_only)} cleaned antigen queries to {query_csv}")
