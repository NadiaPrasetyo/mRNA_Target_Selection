import pandas as pd
import re

# === CONFIGURATION ===
input_excel = "data/S.aureus/IEDB_S.aureus_antigen_table.xlsx"  # Replace with your actual file name
output_txt = "data/S.aureus/uniprot_ids.txt"

# === READ EXCEL FILE ===
df = pd.read_excel(input_excel)
antigen_names = df.iloc[:, 0].dropna()

# === PROCESS ANTIGEN NAMES ===
uniprot_ids = []
missing_uniprot = []

for name in antigen_names:
    match = re.search(r'UniProt:([A-Z0-9]+)', name)
    if match:
        uniprot_ids.append(match.group(1))
    else:
        missing_uniprot.append(name.strip())

# === WRITE TO TEXT FILE ===
with open(output_txt, "w") as f:
    # Write UniProt IDs
    for uid in sorted(set(uniprot_ids)):
        f.write(uid + "\n")
    
    # Write fallback list if any are missing IDs
    if missing_uniprot:
        f.write("# Antigens without UniProt IDs")
        for item in missing_uniprot:
            f.write("\n# "+ item )

print(f"Wrote {len(uniprot_ids)} UniProt IDs and {len(missing_uniprot)} unmatched antigen names to '{output_txt}'")