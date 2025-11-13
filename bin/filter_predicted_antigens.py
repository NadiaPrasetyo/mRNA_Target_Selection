import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Filter CSV by allergenicity and binder criteria."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input CSV file."
    )
    parser.add_argument(
        "-o", "--output_file",
        type=str,
        default=None,
        help="Optional path for the output CSV file. "
             "Defaults to results/filtered_<input_filename_stem>.csv"
    )

    args = parser.parse_args()

    input_path = Path(args.input_file)
    output_path = (
        Path(args.output_file)
        if args.output_file
        else Path("results") / f"filtered_{input_path.stem}.csv"
    )

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load the CSV file
    df = pd.read_csv(input_path)

    # Apply filters
    filtered_df = df[
        (df["allergenicity_hybrid_score"] < 0.3) &
        (df["mhci_num_strong_binders"] > 0) &
        (df["mhcii_num_strong_binders"] > 0)
    ]

    # Save the filtered results
    filtered_df.to_csv(output_path, index=False)

    print(f"✅ Filtered data saved to: {output_path}")
    print(f"Rows before filtering: {len(df)}")
    print(f"Rows after filtering: {len(filtered_df)}")

if __name__ == "__main__":
    main()
