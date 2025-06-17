"""
extract_epitopes.py
Command-line tool to extract unique MHCI and MHCII epitopes from JSON output files.

Overview:
    - Scans specified directories for MHCI and MHCII epitope prediction JSON files.
    - Extracts peptide sequences passing user-defined IC50 and percentile thresholds.
    - Aggregates and deduplicates epitopes across all input files.
    - Optionally writes the list of unique epitopes to an output file.

Arguments:
    epitope_dir (str): Directory containing 'mhci' and 'mhcii' subdirectories with JSON files.
    --ic50-threshold (float, optional): IC50 threshold for MHCI epitopes (default: 500.0).
    --mhci-percentile (float, optional): Percentile threshold for MHCI epitopes (default: 2.0).
    --mhcii-percentile (float, optional): Percentile threshold for MHCII epitopes (default: 10.0).
    --output-file (str, optional): Output file to save extracted epitopes.

Requirements:
    - JSON files with peptide prediction results in the specified directory structure.
    - Python packages: argparse, pathlib, json.

Usage Example:
    python extract_epitopes.py results/epitopes --ic50-threshold 250 --mhci-percentile 1.5 --mhcii-percentile 5 --output-file selected_epitopes.txt

Outputs:
    Prints the number of peptides extracted from each file and the total unique epitopes.
    Optionally writes sorted unique epitopes to the specified output file.

Author: Nadia
"""
import json
from pathlib import Path
import argparse

def extract_mhci_epitopes(file_path, ic50_threshold=500.0, percentile_threshold=2.0):
    """
    Extract MHCI epitopes from a JSON file based on IC50 and percentile thresholds.
    Args:
        file_path (str or Path): Path to the JSON file containing MHCI results.
        ic50_threshold (float): IC50 threshold for filtering epitopes.
        percentile_threshold (float): Percentile threshold for filtering epitopes.
    Returns:
        Set[str]: A set of unique peptide sequences that meet the thresholds.
    """
    with open(file_path) as f:
        data = json.load(f)

    epitopes = set()
    for result in data.get("results", []):
        if result.get("type") != "peptide_table":
            continue
        cols = result["table_columns"]
        rows = result["table_data"]
        for row in rows:
            row_data = dict(zip(cols, row))
            ic50 = float(row_data.get("ic50", 99999))
            percentile = float(row_data.get("percentile", 100.0))
            if ic50 <= ic50_threshold and percentile <= percentile_threshold:
                epitopes.add(row_data["peptide"])
    return epitopes

def extract_mhcii_epitopes(file_path, percentile_threshold=10.0):
    """
    Extract MHCII epitopes from a JSON file based on percentile threshold.
    Args:
        file_path (str or Path): Path to the JSON file containing MHCII results.
        percentile_threshold (float): Percentile threshold for filtering epitopes.
    Returns:
        Set[str]: A set of unique peptide sequences that meet the percentile threshold.
    """
    with open(file_path) as f:
        data = json.load(f)

    epitopes = set()
    for result in data.get("results", []):
        if result.get("type") != "peptide_table":
            continue
        cols = result["table_columns"]
        rows = result["table_data"]
        for row in rows:
            row_data = dict(zip(cols, row))
            percentile = float(row_data.get("percentile", 100.0))
            if percentile <= percentile_threshold:
                epitopes.add(row_data["peptide"])
    return epitopes

def extract_all_epitopes(epitope_dir, ic50_threshold=500.0, mhci_percentile=2.0, mhcii_percentile=10.0):
    """
    Collect unique epitopes from MHCI and MHCII JSON files under `epitope_dir/mhci/` and `epitope_dir/mhcii/`.
    Args:
        epitope_dir (str or Path): Directory containing 'mhci' and 'mhcii' subdirectories with JSON files.
        ic50_threshold (float): IC50 threshold for MHCI epitopes.
        mhci_percentile (float): Percentile threshold for MHCI epitopes.
        mhcii_percentile (float): Percentile threshold for MHCII epitopes.
    Returns: 
        Set[str]
    """
    epitope_dir = Path(epitope_dir)
    mhci_dir = epitope_dir / "mhci"
    mhcii_dir = epitope_dir / "mhcii"

    all_epitopes = set()

    for file in mhci_dir.glob("*.json"):
        all_epitopes.update(extract_mhci_epitopes(file, ic50_threshold, mhci_percentile))

    for file in mhcii_dir.glob("*.json"):
        all_epitopes.update(extract_mhcii_epitopes(file, mhcii_percentile))

    return all_epitopes


def main(epitope_dir, ic50_threshold, mhci_percentile, mhcii_percentile, output_file):
    """
    Main function to extract and optionally save unique epitopes from specified directories.
    Args:
        epitope_dir (str or Path): Directory containing 'mhci' and 'mhcii' subdirectories with JSON files.
        ic50_threshold (float): IC50 threshold for MHCI epitopes.
        mhci_percentile (float): Percentile threshold for MHCI epitopes.
        mhcii_percentile (float): Percentile threshold for MHCII epitopes.
        output_file (str, optional): Output file to save extracted epitopes.
    """
    epitope_dir = Path(epitope_dir)
    mhci_dir = epitope_dir / "mhci"
    mhcii_dir = epitope_dir / "mhcii"

    all_epitopes = set()

    print("🔍 Scanning MHCI JSON files...")
    for file in mhci_dir.glob("*.json"):
        peptides = extract_mhci_epitopes(file, ic50_threshold, mhci_percentile)
        print(f"  - {file.name}: {len(peptides)} peptides")
        all_epitopes.update(peptides)

    print("🔍 Scanning MHCII JSON files...")
    for file in mhcii_dir.glob("*.json"):
        peptides = extract_mhcii_epitopes(file, mhcii_percentile)
        print(f"  - {file.name}: {len(peptides)} peptides")
        all_epitopes.update(peptides)

    print(f"\n✅ Total unique epitopes extracted: {len(all_epitopes)}")

    if output_file:
        with open(output_file, "w") as out:
            for epitope in sorted(all_epitopes):
                out.write(epitope + "\n")
        print(f"💾 Epitopes written to: {output_file}")

if __name__ == "__main__":
    """
    Main entry point for the script to enable command-line execution.
    Parses command-line arguments and calls the main function to extract epitopes.
    """
    parser = argparse.ArgumentParser(description="Extract MHCI and MHCII epitopes from JSON output")
    parser.add_argument("epitope_dir", help="Directory containing 'mhci' and 'mhcii' subdirectories with JSON files")
    parser.add_argument("--ic50-threshold", type=float, default=500.0, help="IC50 threshold for MHCI epitopes")
    parser.add_argument("--mhci-percentile", type=float, default=2.0, help="Percentile threshold for MHCI epitopes")
    parser.add_argument("--mhcii-percentile", type=float, default=10.0, help="Percentile threshold for MHCII epitopes")
    parser.add_argument("--output-file", help="Output file to save extracted epitopes (optional)")

    args = parser.parse_args()
    main(args.epitope_dir, args.ic50_threshold, args.mhci_percentile, args.mhcii_percentile, args.output_file)
