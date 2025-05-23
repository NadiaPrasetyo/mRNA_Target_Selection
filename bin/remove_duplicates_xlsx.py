import pandas as pd
import re

# === CONFIGURATION ===
file_with_ids = "data/S.aureus/IEDB_S.aureus_antigen_table.xlsx"     # First file with UniProt IDs
file_new_antigens = "data/S.aureus/Howden2023_S.aureus_virulence-factor-list.xlsx" # Second file with just names
query_list_output = "data/S.aureus/entrez_queries.txt"

# === READ FILES ===
df_with_ids = pd.read_excel(file_with_ids)
df_new = pd.read_excel(file_new_antigens)

# Extract antigen names (first column)
known_antigens = df_with_ids.iloc[:, 0].dropna().str.replace(r"\s*\(UniProt:[^)]+\)", "", regex=True).str.strip()
new_antigens = df_new.iloc[:, 0].dropna().str.strip()

# === FIND NEW (NON-DUPLICATE) ANTIGENS ===
new_only = sorted(set(new_antigens) - set(known_antigens))

# === GENERATE QUERY STRINGS FOR SHELL SCRIPT ===
with open(query_list_output, "w") as f:
    for antigen in new_only:
        query = f'"{antigen} [All Fields] AND \\"Staphylococcus aureus\\"[Organism] AND swissprot[filter]"'
        fallback_query = f'"{antigen} [All Fields] AND \\"Staphylococcus aureus\\"[Organism]"'
        f.write(f"{antigen}|||{query}|||{fallback_query}\n")

print(f"Prepared {len(new_only)} Entrez queries in '{query_list_output}'")
