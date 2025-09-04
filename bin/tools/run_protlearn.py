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
    """Run all protlearn feature extractors on a protein sequence."""
    features = {}
    for name, func in FEATURE_FUNCTIONS.items():
        try:
            if name == "aaindex1":
                # Filter aaindex1 to include only specific indices
                selected_indices = [
                    "ARGP820101", "JOND750101", "BHAR880101", "CHOC750101",
                    "DAYM780101", "DAYM780201", "GRAR740101", "GRAR740102",
                    "GRAR740103", "JOND750102", "KYTJ820101"
                ]
                arr, desc = func([seq], props=selected_indices)
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
    """Main entry point: extract features from PDB or CIF and write to CSV."""
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
