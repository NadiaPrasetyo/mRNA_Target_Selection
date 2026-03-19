#!/usr/bin/env python3
"""
Combine DeeplocPro, SignalP, and TargetP prediction outputs into unified feature matrices
for each genome (matching by stem, e.g. GCF_000009645.1).

Author: [Your Name]
Date: 2025-10-24
"""

import os
import re
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


def parse_deeplocpro_file(path):
    """Parse a single DeeplocPro CSV file."""
    results = []
    try:
        df = pd.read_csv(path)
        prob_cols = [
            ("cell_wall_surface", "Cell wall & surface"),
            ("extracellular", "Extracellular"),
            ("cytoplasmic", "Cytoplasmic"),
            ("cytoplasmic_membrane", "Cytoplasmic Membrane"),
            ("outer_membrane", "Outer Membrane"),
            ("periplasmic", "Periplasmic"),
        ]
        for i, row in df.iterrows():
            try:
                accession = f"{row['ACC'].split('|')[1]}" if '|' in row['ACC'] else row['ACC']
            except IndexError:
                accession = "unknown"
            for loc, col in prob_cols:
                try:
                    prob = float(row[col])
                    results.append({
                        "accession": accession,
                        "feature": "deeplocpro",
                        "subfeature": f"prob_{loc}",
                        "value": prob
                    })
                except (ValueError, KeyError):
                    continue
    except Exception as e:
        logging.error(f"DEEPLOCPRO: Failed parsing {path}: {e}")
    return results


def parse_signalp_file(path):
    """Parse a single SignalP .txt file."""
    results = []
    try:
        with open(path) as f:
            for i, line in enumerate(f):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    try:
                        accession = parts[0].split('|')[1] if '|' in parts[0] else parts[0]
                        results.append({"accession": accession, "feature": "signalp", "subfeature": "prob_signalp", "value": float(parts[2])})
                        results.append({"accession": accession, "feature": "signalp", "subfeature": "prob_other", "value": float(parts[3])})
                    except Exception:
                        continue
    except Exception as e:
        logging.error(f"SIGNALP: Failed parsing {path}: {e}")
    return results


def parse_targetp_file(path):
    """Parse a single TargetP .txt file."""
    results = []
    try:
        with open(path) as f:
            for i, line in enumerate(f):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    try:
                        accession = parts[0].split('|')[1] if '|' in parts[0] else parts[0]
                        results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_noTP", "value": float(parts[2])})
                        results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_SP", "value": float(parts[3])})
                        results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_mTP", "value": float(parts[4])})
                    except Exception:
                        continue
    except Exception as e:
        logging.error(f"TARGETP: Failed parsing {path}: {e}")
    return results


def extract_stem(filename):
    """
    Extract a genome stem (e.g. GCF_000009645.1) from a file name.
    """
    match = re.search(r"(GCF_\d+\.\d+)", filename)
    return match.group(1) if match else None


def combine_predictions_by_stem(directory, output_dir=None):
    """
    Combine DeeplocPro, SignalP, and TargetP predictions for each genome stem.
    Missing feature columns are filled with NA.
    """
    files = os.listdir(directory)
    stems = sorted({extract_stem(f) for f in files if extract_stem(f)})

    all_matrices = {}

    # define the complete set of expected columns
    expected_columns = [
        # DeeplocPro
        "deeplocpro_localisation",
        "deeplocpro_prob_cell_wall_surface",
        "deeplocpro_prob_extracellular",
        "deeplocpro_prob_cytoplasmic",
        "deeplocpro_prob_cytoplasmic_membrane",
        "deeplocpro_prob_outer_membrane",
        "deeplocpro_prob_periplasmic",
        # SignalP
        "signalp_prob_signalp",
        "signalp_prob_other",
        # TargetP
        "targetp_prob_noTP",
        "targetp_prob_SP",
        "targetp_prob_mTP",
    ]

    for stem in stems:
        logging.info(f"Processing genome stem: {stem}")

        deeplocpro_file = next((os.path.join(directory, f) for f in files if f.endswith(".csv") and stem in f), None)
        signalp_file = next((os.path.join(directory, f) for f in files if "signalp" in f and stem in f and f.endswith(".txt")), None)
        targetp_file = next((os.path.join(directory, f) for f in files if "targetp" in f and stem in f and f.endswith(".txt")), None)

        results = []
        if deeplocpro_file:
            results += parse_deeplocpro_file(deeplocpro_file)
        if signalp_file:
            results += parse_signalp_file(signalp_file)
        if targetp_file:
            results += parse_targetp_file(targetp_file)

        if not results:
            logging.warning(f"No parsed results found for {stem}")
            continue

        df = pd.DataFrame(results)
        matrix = df.pivot_table(
            index="accession",
            columns=["feature", "subfeature"],
            values="value",
            aggfunc="first"
        )

        # flatten and fill missing features with NA
        matrix.columns = [f"{feat}_{sub}" for feat, sub in matrix.columns]
        matrix = matrix.reindex(columns=expected_columns, fill_value=pd.NA)
        matrix.reset_index(inplace=True)

        all_matrices[stem] = matrix

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"{stem}_combined_features.csv")
            matrix.to_csv(output_path, index=False, na_rep="NA")
            logging.info(f"Saved combined feature matrix for {stem} → {output_path}")

    return all_matrices

def filter_remove_cytoplasmic_no_signalp(results):
    """
    Filter out rows where the highest scoring DeeplocPro feature is "cytoplasmic"
    and the highest scoring SignalP feature is "other" and the highest scoring
    TargetP feature is "noTP".

    Parameters
    ----------
    results : pandas.DataFrame
        DataFrame containing the combined feature matrices.

    Returns
    -------
    filtered_results : pandas.DataFrame
        DataFrame containing the filtered feature matrices.
    """
    deeplocpro_cols = [c for c in results.columns if c.startswith("deeplocpro_prob_")]  # bug fix: was `df`, also scoped to prob_ cols
    signalp_cols    = [c for c in results.columns if c.startswith("signalp_")]
    targetp_cols    = [c for c in results.columns if c.startswith("targetp_")]

    deeplocpro_idxmax = results[deeplocpro_cols].idxmax(axis=1)
    signalp_idxmax    = results[signalp_cols].idxmax(axis=1)
    targetp_idxmax    = results[targetp_cols].idxmax(axis=1)

    # Rows to REMOVE
    mask_remove = (
        (deeplocpro_idxmax == "deeplocpro_prob_cytoplasmic")  # condition 1
        |
        (                                                      # condition 2
            (signalp_idxmax == "signalp_prob_other") &
            (targetp_idxmax == "targetp_prob_noTP")
        )
    )

    return results[~mask_remove]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combine DeeplocPro, SignalP, and TargetP outputs by genome stem.")
    parser.add_argument(
        "-d", "--directory",
        required=True,
        help="Path to directory containing prediction output files."
    )
    parser.add_argument(
        "-o", "--output_dir",
        required=False,
        default="data/"
        help="Optional directory to save per-stem CSV outputs."
    )
    parser.add_argument(
        "--filter",
        action="store_true",
        help="Filter out proteins that were predicted to be localised in the cytoplasm and did not have any signal peptide were excluded according to the scores predicted by the tools."
    )
    args = parser.parse_args()

    all_matrices = combine_predictions_by_stem(args.directory, args.output_dir)

    for stem, matrix in all_matrices.items():                 
        if args.filter:
            matrix = filter_remove_cytoplasmic_no_signalp(matrix)

        output_path = os.path.join(args.output_dir, f"{stem}_combined_features.csv")
        matrix.to_csv(output_path, index=False, na_rep="NA") # overwrites the file in-place
        logging.info(f"Saved filtered matrix for {stem} → {output_path}")

    logging.info(f"Processed {len(results)} genomes successfully.")
