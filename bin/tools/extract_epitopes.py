"""
extract_epitopes.py

Command-line tool to extract unique MHCI and MHCII epitopes from JSON output files.

Author: Nadia
"""

import argparse
import json
from pathlib import Path
from typing import Set, Tuple


def extract_mhci_epitopes(file_path: Path, ic50_threshold: float, percentile_threshold: float) -> Set[Tuple[str, str]]:
    """Extracts MHCI peptides from a single JSON file."""
    return _extract_epitopes_from_file(
        file_path,
        ic50_threshold=ic50_threshold,
        percentile_threshold=percentile_threshold,
        mhc_class="mhci"
    )


def extract_mhcii_epitopes(file_path: Path, percentile_threshold: float) -> Set[Tuple[str, str]]:
    """Extracts MHCII peptides from a single JSON file."""
    return _extract_epitopes_from_file(
        file_path,
        ic50_threshold=None,
        percentile_threshold=percentile_threshold,
        mhc_class="mhcii"
    )


def _extract_epitopes_from_file(
    file_path: Path,
    ic50_threshold: float = None,
    percentile_threshold: float = 100.0,
    mhc_class: str = "mhci"
) -> Set[Tuple[str, str]]:
    """
    Generalized parser for MHC epitope JSON files.

    Returns:
        A set of (peptide, allele) pairs passing the thresholds.
    """
    try:
        with file_path.open() as f:
            data = json.load(f)
    except Exception as e:
        print(f"⚠️ Failed to read {file_path.name}: {e}")
        return set()

    epitopes = set()
    for result in data.get("results", []):
        if result.get("type") != "peptide_table":
            continue

        columns = result.get("table_columns", [])
        for row in result.get("table_data", []):
            row_data = dict(zip(columns, row))
            peptide = row_data.get("peptide")
            allele = row_data.get("allele", "unknown")
            if not peptide:
                continue

            percentile = float(row_data.get("percentile", 100.0))
            if percentile > percentile_threshold:
                continue

            if mhc_class == "mhci":
                ic50 = float(row_data.get("ic50", 99999))
                if ic50 > ic50_threshold:
                    continue

            epitopes.add((peptide, allele))

    return epitopes


def extract_all_epitopes(
    epitope_dir: Path,
    ic50_threshold: float = 500.0,
    mhci_percentile: float = 2.0,
    mhcii_percentile: float = 10.0
) -> Set[Tuple[str, str]]:
    """
    Extracts all epitopes from both MHCI and MHCII subdirectories.
    """
    all_epitopes = set()

    mhci_dir = epitope_dir / "mhci"
    if mhci_dir.exists():
        for file in mhci_dir.glob("*.json"):
            all_epitopes.update(extract_mhci_epitopes(file, ic50_threshold, mhci_percentile))

    mhcii_dir = epitope_dir / "mhcii"
    if mhcii_dir.exists():
        for file in mhcii_dir.glob("*.json"):
            all_epitopes.update(extract_mhcii_epitopes(file, mhcii_percentile))

    return all_epitopes


def write_epitopes_to_file(epitopes: Set[Tuple[str, str]], output_file: Path):
    """
    Writes sorted epitopes to the output file.
    """
    with output_file.open("w") as f:
        for peptide, allele in sorted(epitopes):
            f.write(f"{peptide}\t{allele}\n")


def main(epitope_dir, ic50_threshold, mhci_percentile, mhcii_percentile, output_file):
    """
    Extracts epitopes and optionally writes them to a file.
    """
    epitope_dir = Path(epitope_dir)
    all_epitopes = extract_all_epitopes(
        epitope_dir,
        ic50_threshold,
        mhci_percentile,
        mhcii_percentile
    )

    print(f"✅ Total unique (peptide, allele) pairs: {len(all_epitopes)}")

    if output_file:
        output_file = Path(output_file)
        write_epitopes_to_file(all_epitopes, output_file)
        print(f"💾 Epitopes written to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract unique MHCI and MHCII epitopes from prediction JSON files."
    )
    parser.add_argument("epitope_dir", help="Directory with 'mhci' and 'mhcii' subfolders")
    parser.add_argument("--ic50-threshold", type=float, default=500.0,
                        help="IC50 threshold for MHCI (default: 500)")
    parser.add_argument("--mhci-percentile", type=float, default=2.0,
                        help="Percentile threshold for MHCI (default: 2.0)")
    parser.add_argument("--mhcii-percentile", type=float, default=10.0,
                        help="Percentile threshold for MHCII (default: 10.0)")
    parser.add_argument("--output-file", help="Optional output file to save extracted epitopes")

    args = parser.parse_args()
    main(
        args.epitope_dir,
        args.ic50_threshold,
        args.mhci_percentile,
        args.mhcii_percentile,
        args.output_file
    )
