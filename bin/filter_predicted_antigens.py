import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Filter raw CSV by allergenicity/binder criteria and cross-check with predictions."
    )
    parser.add_argument(
        "--input-raw",
        required=True,
        help="Path to the main raw CSV file (with feature data)."
    )
    parser.add_argument(
        "--input-pred",
        required=True,
        help="Path to the prediction CSV file (with accession, prob_antigen, etc.)."
    )
    parser.add_argument(
        "-o", "--output-file",
        default=None,
        help="Optional path for the output CSV file. "
             "Defaults to results/filtered_<input_pred_filename_stem>.csv"
    )

    args = parser.parse_args()

    # Define paths
    raw_path = Path(args.input_raw)
    pred_path = Path(args.input_pred)
    output_path = (
        Path(args.output_file)
        if args.output_file
        else Path("results") / f"filtered_{pred_path.stem}.csv"
    )

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # --- Load the CSV files ---
    df_raw = pd.read_csv(raw_path)
    df_pred = pd.read_csv(pred_path)

    # --- Filter raw data ---
    filtered_raw = df_raw[
        (df_raw["allergenicity_hybrid_score"] < 0.3) &
        (df_raw["mhci_num_strong_binders"] > 0) &
        (df_raw["mhcii_num_strong_binders"] > 0)
    ]

    removed = df_raw[~df_raw.index.isin(filtered_raw.index)]
    print(f"ℹ️  Removed {len(removed)} rows that did not meet allergenicity/binder criteria.")

    # --- Merge with predictions ---
    merged = pd.merge(
        filtered_raw,
        df_pred,
        on="accession",
        how="inner",  # only keep matching accessions
        suffixes=("_raw", "_pred")  # rename duplicate columns clearly
    )

    merged_removed = pd.merge(
        removed,
        df_pred,
        on="accession",
        how="outer",  # keep all accessions to see which were missing
        indicator=True,
        suffixes=("_raw", "_pred")  # rename duplicate columns clearly
    )
    
    # --- Determine which prob_antigen column to use ---
    prob_col = "prob_antigen_raw" if "prob_antigen_raw" in merged.columns else "prob_antigen"


    # --- Select only desired columns ---
    final_df = merged[[
        "accession",
        prob_col,
        "pred_label",
        "protein_names",
        "gene_names"
    ]]

    removed_df = merged_removed[[
        "accession",
        prob_col,
        "pred_label",
        "protein_names",
        "gene_names",
        "_merge"
    ]]

    # --- Save results ---
    final_df.to_csv(output_path, index=False)
    removed_output_path = output_path.parent / f"removed_{pred_path.stem}.csv"
    removed_df.to_csv(removed_output_path, index=False)

    print(f"✅ Filtered and merged data saved to: {output_path}")
    print(f"Rows before filtering: {len(df_raw)}")
    print(f"Rows after filtering: {len(filtered_raw)}")
    print(f"Rows after merge: {len(final_df)}")

if __name__ == "__main__":
    main()

