import argparse
import pandas as pd
import os
from glob import glob

def load_iedb_antigens(file_path):
    df = pd.read_csv(file_path)
    df_out = df.rename(columns={
        "source_organism_names": "source_organism",
        "host_organism_names": "host_organisms",
        "parent_source_antigen_names": "antigen_name"
    })
    df_out["source"] = "IEDB"
    return df_out[["source_organism", "host_organisms", "antigen_name", "source"]]

def load_literature_antigens(file_path, source_organism):
    df = pd.read_excel(file_path)
    df_out = pd.DataFrame()
    df_out["antigen_name"] = df["Name"]
    df_out["source_organism"] = source_organism
    df_out["host_organisms"] = "Homo sapiens"
    df_out["source"] = "literature"
    return df_out[["source_organism", "host_organisms", "antigen_name", "source"]]

def main(short_name, long_name):
    base_path = f"data/{short_name}"

    # Load IEDB
    iedb_file = os.path.join(base_path, "IEDB_antigens.csv")
    iedb_df = load_iedb_antigens(iedb_file)

    # Load literature/patents (XLSX files)
    literature_dfs = []
    for xlsx_file in glob(os.path.join(base_path, "*.xlsx")):
        lit_df = load_literature_antigens(xlsx_file, long_name)
        literature_dfs.append(lit_df)
    literature_df = pd.concat(literature_dfs, ignore_index=True) if literature_dfs else pd.DataFrame()

    # Combine all sources
    combined_df = pd.concat([iedb_df, literature_df], ignore_index=True)

    # Output
    output_file = os.path.join(base_path, f"{short_name}_compiled_antigens.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"Compiled antigen data saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile antigen data for a pathogen")
    parser.add_argument("short_name", help="Short name for the pathogen (e.g. S.aureus)")
    parser.add_argument("long_name", help="Full name for the pathogen (e.g. Staphylococcus aureus)")
    args = parser.parse_args()
    main(args.short_name, args.long_name)
