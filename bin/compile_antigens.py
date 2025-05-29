import sys
import pandas as pd
import os
from glob import glob
import re

def extract_uniprot_id(iri):
    if pd.isna(iri):
        return None
    match = re.search(r"(?:UNIPROT:)?([A-Z0-9]+)$", iri)
    if match:
        return f"{match.group(1)}"
    return None

def load_iedb_antigens(file_path):
    if not os.path.exists(file_path):
        print(f"Warning: IEDB file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_csv(file_path)

    df_out = df.rename(columns={
        "source_organism_names": "source_organism",
        "host_organism_names": "host_organisms",
        "parent_source_antigen_names": "antigen_name"
    })

    # Extract and format Uniprot_ID
    if "parent_source_antigen_iri" in df.columns:
        df_out["Uniprot_ID"] = df["parent_source_antigen_iri"].apply(extract_uniprot_id)
    else:
        df_out["Uniprot_ID"] = None
    df_out["gene_name"] = None  # IEDB does not provide gene names
    df_out["source"] = "IEDB"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]

def load_literature_antigens(file_path, source_organism):
    if not os.path.exists(file_path):
        print(f"Warning: Literature file not found: {file_path}. Skipping.")
        return pd.DataFrame()
    df = pd.read_excel(file_path)
    df_out = pd.DataFrame()
    df_out["antigen_name"] = df["Name"]
    df_out["gene_name"] = df["Gene"]
    df_out["source_organism"] = source_organism
    df_out["host_organisms"] = "Homo sapiens"
    df_out["Uniprot_ID"] = None
    df_out["source"] = "literature"
    return df_out[["source_organism", "host_organisms", "antigen_name", "gene_name", "Uniprot_ID", "source"]]

def main(short_name, long_name):
    organism_tag = str(long_name).replace(" ", "_").lower()
    base_path = f"data/{short_name}"

    # Load IEDB
    iedb_file = os.path.join(base_path, f"{organism_tag}_IEDB_antigens.csv")
    iedb_df = load_iedb_antigens(iedb_file)

    # Load literature/patents (XLSX files)
    literature_dfs = []
    for xlsx_file in glob(os.path.join(base_path, "*.xlsx")):
        lit_df = load_literature_antigens(xlsx_file, long_name)
        if not lit_df.empty:
            literature_dfs.append(lit_df)
    literature_df = pd.concat(literature_dfs, ignore_index=True) if literature_dfs else pd.DataFrame()

    # Combine all sources
    combined_df = pd.concat([iedb_df, literature_df], ignore_index=True)
    if combined_df.empty:
        print("No antigen data found. Exiting.")
        return

    # Output
    output_file = os.path.join(base_path, f"{organism_tag}_compiled_antigens.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"Compiled antigen data saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compile_antigens.py <short_name> <long_name>")
        sys.exit(1)
    short_name = sys.argv[1]
    long_name = sys.argv[2]
    if not os.path.exists(f"data/{short_name}"):
        os.makedirs(f"data/{short_name}")
        print(f"Created directory: data/{short_name}")
    main(short_name, long_name)
