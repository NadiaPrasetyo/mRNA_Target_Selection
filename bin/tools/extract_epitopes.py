import json
from pathlib import Path
import argparse

def extract_mhci_epitopes(file_path, ic50_threshold=500.0, percentile_threshold=2.0):
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
    Returns: Set[str]
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
    parser = argparse.ArgumentParser(description="Extract MHCI and MHCII epitopes from JSON output")
    parser.add_argument("epitope_dir", help="Directory containing 'mhci' and 'mhcii' subdirectories with JSON files")
    parser.add_argument("--ic50-threshold", type=float, default=500.0, help="IC50 threshold for MHCI epitopes")
    parser.add_argument("--mhci-percentile", type=float, default=2.0, help="Percentile threshold for MHCI epitopes")
    parser.add_argument("--mhcii-percentile", type=float, default=10.0, help="Percentile threshold for MHCII epitopes")
    parser.add_argument("--output-file", help="Output file to save extracted epitopes (optional)")

    args = parser.parse_args()
    main(args.epitope_dir, args.ic50_threshold, args.mhci_percentile, args.mhcii_percentile, args.output_file)
