import pandas as pd
import re

# === CONFIGURATION ===
file_with_ids = "data/S.aureus/IEDB_S.aureus_antigen_table.xlsx"     # First file with UniProt IDs
file_new_antigens = "data/S.aureus/Howden2023_S.aureus_virulence-factor-list.xlsx" # Second file with just names
query_list_output = "data/S.aureus/entrez_queries.txt"

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

# === GENERATE QUERY FILE ===
with open(query_list_output, "w") as f:
    for antigen in new_only:
        query = f'"{antigen} [All Fields] NOT partial[All Fields] AND \\"Staphylococcus aureus\\"[Organism] AND swissprot[filter]"'
        fallback_query = f'"{antigen} [All Fields] NOT partial[All Fields] AND \\"Staphylococcus aureus\\"[Organism]"'
        f.write(f"{antigen}|||{query}|||{fallback_query}\n")

print(f"Prepared {len(new_only)} Entrez queries (excluding partials) in '{query_list_output}'")