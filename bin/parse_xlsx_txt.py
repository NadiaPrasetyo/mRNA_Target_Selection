import pandas as pd
import re
import argparse

# === ARGUMENT PARSING ===
parser = argparse.ArgumentParser(description="Extract UniProt IDs from an IEDB Excel table.")
parser.add_argument('--input', '-i', required=True, help="Path to input Excel file")
parser.add_argument('--output', '-o', required=True, help="Path to intermediate UniProt ID text file")
args = parser.parse_args()

input_excel = args.input
output_txt = args.output

# === READ EXCEL FILE ===
df = pd.read_excel(input_excel)
antigen_names = df.iloc[:, 0].dropna()

# === PROCESS ANTIGEN NAMES ===
uniprot_ids = []
missing_uniprot = []

for name in antigen_names:
    match = re.search(r'UniProt:([A-Z0-9]+)', str(name))
    if match:
        uniprot_ids.append(match.group(1))
    else:
        missing_uniprot.append(str(name).strip())

# === WRITE TO TEXT FILE ===
with open(output_txt, "w") as f:
    # Write UniProt IDs
    for uid in sorted(set(uniprot_ids)):
        f.write(uid + "\n")
    
    # Write fallback list if any are missing IDs
    if missing_uniprot:
        f.write("# Antigens without UniProt IDs")
        for item in missing_uniprot:
            f.write("\n# " + item)

print(f"[INFO] Wrote {len(uniprot_ids)} UniProt IDs to temporary file '{output_txt}'")
