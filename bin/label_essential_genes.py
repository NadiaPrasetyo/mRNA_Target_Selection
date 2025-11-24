import argparse
import pandas as pd
from pathlib import Path
import re

def main():
    # --- Parse arguments ---
    parser = argparse.ArgumentParser(
        description="Add essential/non-essential labels to a gene/protein CSV."
    )
    parser.add_argument(
        "--input-csv",
        required=True,
        help="Input CSV file with proteins/genes to be labeled.",
    )
    parser.add_argument(
        "--essential-csv",
        required=True,
        help="CSV file listing essential genes (must include column 'Gene').",
    )
    parser.add_argument(
        "--output-file",
        required=False,
        help="Optional: Output CSV file (defaults to {input-file-stem}_essential.csv).",
    )

    args = parser.parse_args()

    # --- Resolve paths ---
    input_csv = args.input_csv
    essential_csv = args.essential_csv

    if args.output_file:
        output_csv = args.output_file
    else:
        output_csv = Path("results") / f"{Path(args.input_csv).stem}_essential.csv"

    # --- Load CSVs ---
    print(f"📖 Reading input file: {input_csv}")
    main_df = pd.read_csv(input_csv)

    print(f"📖 Reading essential gene list: {essential_csv}")
    ess_df = pd.read_csv(essential_csv)

    # --- Prepare essential gene list ---
    essential_genes = (
        ess_df["Gene"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.lower()
        .drop_duplicates()
        .tolist()
    )

    # --- Function to check if any gene in a row is essential ---
    def label_essential(gene_names_str):
        if pd.isna(gene_names_str):
            return "non-essential"
        
        # Split on commas OR semicolons OR whitespace
        genes = re.split(r"[;,]", gene_names_str.lower())
        genes = [g.strip() for g in genes if g.strip()]
        
        return "essential" if any(g in essential_genes for g in genes) else "non-essential"

    # --- Apply labeling ---
    main_df["essential_label"] = main_df["gene_names"].apply(label_essential)

    # --- Write output ---
    main_df.to_csv(output_csv, index=False)
    print(f"✅ Done! Saved labeled CSV to: {output_csv}")

if __name__ == "__main__":
    main()
