"""
Feature Extraction for Protein Sequences.
Overview:
    - Extracts a comprehensive set of features from a given protein sequence using the protlearn library.
    - Supports multiple feature extraction methods, including aaindex1 and others defined in FEATURE_FUNCTIONS.
    - Handles exceptions gracefully to ensure robust feature extraction.

Arguments:
    seq (str): A protein sequence represented as a string of amino acid characters.

Process:
    - Iterates through all feature extraction functions defined in FEATURE_FUNCTIONS.
    - For the "aaindex1" feature, filters the extracted indices to include only a predefined set of indices.
    - For other features, extracts all available descriptors and their corresponding values.
    - Handles cases where the feature extraction function returns either a tuple (array, descriptors) or a single output.

Error Handling:
    - If a feature extraction function raises an exception, attempts to process the output in an alternative way.
    - If the alternative processing also fails, records the error message for the corresponding feature.

Outputs:
    - A dictionary where keys are feature names (or feature_name_descriptor for detailed features) and values are the extracted feature values.
    - For features that encounter errors, the value will be an error message string.

Requirements:
    - Python packages: protlearn.
    - FEATURE_FUNCTIONS must be a dictionary mapping feature names to their corresponding extraction functions.
    - The "aaindex1" feature requires a predefined set of indices to filter the extracted descriptors.

Notes:
    - This function is designed to handle protein sequences only.
    - Ensure that the FEATURE_qFUNCTIONS dictionary is properly configured before calling this function.
    - Logs or additional debugging information can be added to enhance traceability during execution.

Author: Nadia
"""
import os
import csv
import logging
from Bio import PDB
from Bio.PDB import MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1

from protlearn.features import (
    length,
    aaindex1
)

FEATURE_FUNCTIONS = {
    'length': length,
    'aaindex1': aaindex1
}

def extract_sequence(file_path):
    """
    Extract amino acid sequence (1-letter) from a PDB or CIF file.
    Skips invalid/unsupported/compressed files instead of crashing.
    Args:
        file_path (str or Path): Path to the input PDB or CIF file.
    """
    file_path_str = str(file_path).lower()

    if file_path_str.endswith((".pdb.gz", ".cif.gz")):
        logging.warning(f"Skipping compressed file (decompress first): {file_path}")
        return None

    if not (file_path_str.endswith(".pdb") or file_path_str.endswith(".cif")):
        logging.warning(f"Skipping unsupported file (expected .pdb or .cif): {file_path}")
        return None

    seq = []
    structure = None
    try:
        if file_path_str.endswith(".cif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDB.PDBParser(QUIET=True)

        structure = parser.get_structure("protein", file_path)
    except Exception as e:
        logging.warning(f"Failed structured parse of {file_path}: {e}")
        return extract_sequence_loose(file_path)

    # Extract residues if structure parsed
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    resname = residue.resname.capitalize()
                    if resname in protein_letters_3to1:
                        seq.append(protein_letters_3to1[resname])

    return "".join(seq) if seq else None

def extract_sequence_loose(file_path):
    """
    Fallback: extract sequence by scanning ATOM records directly.
    Only uses CA atoms to avoid duplicates.
    Args:
        file_path (str or Path): Path to the input PDB or CIF file.
    """
    seq = []
    try:
        with open(file_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[13:15].strip() == "CA":
                    resname = line[17:20].strip().capitalize()
                    if resname in protein_letters_3to1:
                        seq.append(protein_letters_3to1[resname])
    except Exception as e:
        logging.warning(f"Loose parse also failed for {file_path}: {e}")
        return None

    return "".join(seq) if seq else None

def extract_all_features(seq):
    """Run all protlearn feature extractors on a protein sequence.
    Args:
        seq (str): Protein sequence as a string of amino acid characters.
    Returns:
        dict: Mapping of feature names to their extracted values.
    """
    features = {}
    for name, func in FEATURE_FUNCTIONS.items():
        try:
            if name == "aaindex1":
                # Filter aaindex1 to include only specific indices
                selected_indices = {
                    "ARGP820101", "JOND750101", "BHAR880101", "CHOC750101",
                    "DAYM780101", "DAYM780201", "GRAR740101", "GRAR740102",
                    "GRAR740103", "JOND750102", "KYTJ820101"
                }
                arr, desc = func([seq])
                values = arr[0]
                for d, v in zip(desc, values):
                    if d in selected_indices:
                        features[f"{name}_{d}"] = v
            else:
                arr, desc = func([seq])
                values = arr[0]
                for d, v in zip(desc, values):
                    features[f"{name}_{d}"] = v
        except Exception as e:
            try:
                out = func([seq])
                if isinstance(out, (list, tuple)) and len(out) == 2:
                    arr, desc = out
                    values = arr[0]
                    for d, v in zip(desc, values):
                        features[f"{name}_{d}"] = v
                else:
                    features[name] = out
            except Exception as inner_e:
                features[name] = f"Error: {inner_e}"
    return features

def run(input_file, tool_root, output_dir):
    """Main entry point: extract features from PDB or CIF and write to CSV.\
    Args:
        input_file : Path to the input PDB or CIF file.
        tool_root (unused): Path to the directory containing external tools kept for interface consistency.
        output_dir : Directory where output files will be saved. Results are placed in a 'protlearn' subfolder.
    Returns:
        Path to the output CSV file containing extracted features.
    """
    output_dir = os.path.join(output_dir, "protlearn")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{input_file.stem}.csv")

    sequence = extract_sequence(input_file)
    if not sequence:
        logging.warning(f"No sequence extracted from {input_file}, skipping.")
        return None

    features = extract_all_features(sequence)

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["feature", "value"])
        for k, v in features.items():
            writer.writerow([k, v])

    logging.info(f"Extracted {len(features)} features from {input_file} → {output_file}")
    return output_file
