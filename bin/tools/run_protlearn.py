import os
import csv
from Bio import PDB
from Bio.PDB import MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1

from protlearn.features import (
    length,
    aac,
    aaindex1
)

FEATURE_FUNCTIONS = {
    'length': length,
    'aac': aac,
    'aaindex1': aaindex1
}

def extract_sequence(file_path):
    """
    Extract amino acid sequence (1-letter) from a PDB or CIF file.
    Falls back to text-based parsing if Biopython chokes.
    """
    seq = []
    structure = None

    try:
        if str(file_path).lower().endswith((".cif", ".mmcif")):
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("protein", file_path)
        else:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("protein", file_path)
    except Exception as e:
        print(f"[WARN] Failed structured parse of {file_path}: {e}")
        return extract_sequence_loose(file_path)

    # Extract residues if structure parsed
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    resname = residue.resname.capitalize()
                    if resname in protein_letters_3to1:
                        seq.append(protein_letters_3to1[resname])
    return "".join(seq)

def extract_sequence_loose(file_path):
    """
    Fallback: extract sequence by scanning ATOM records directly.
    Only uses CA atoms to avoid duplicates.
    """
    seq = []
    with open(file_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                resname = line[17:20].strip().capitalize()
                if resname in protein_letters_3to1:
                    seq.append(protein_letters_3to1[resname])
    return "".join(seq)

def extract_all_features(seq):
    """Run all protlearn feature extractors on a protein sequence."""
    features = {}
    for name, func in FEATURE_FUNCTIONS.items():
        try:
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
    features = extract_all_features(sequence)

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["feature", "value"])
        for k, v in features.items():
            writer.writerow([k, v])

    print(f"Extracted {len(features)} features from {input_file} → {output_file}")
    return output_file
