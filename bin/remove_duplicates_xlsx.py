import pandas as pd
import re


# === CONFIGURATION ===
file_with_ids = "data/S.aureus/IEDB_S.aureus_antigen_table.xlsx"     # First file with UniProt IDs
file_new_antigens = "data/S.aureus/Howden2023_S.aureus_virulence-factor-list.xlsx" # Second file with just names
tsv_output = "data/S.aureus/entrez_queries.tsv"

# === CLEANING FUNCTION ===
def clean_antigen_name(name):
    # Remove parentheses and their contents, truncate at commas or slashes
    name = re.sub(r'\(.*?\)', '', name)
    name = re.split(r'[,/]', name)[0]
    return name.strip()

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
pd.DataFrame(query_rows, columns=["Antigen", "SwissProtQuery", "FallbackQuery"]).to_csv(tsv_output, sep='\t', index=False)
print(f"✅ Wrote {len(query_rows)} cleaned antigen queries to {tsv_output}")
