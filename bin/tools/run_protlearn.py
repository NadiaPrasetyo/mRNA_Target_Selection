import os
import csv
from Bio import PDB
from Bio.Data.IUPACData import protein_letters_3to1

from protlearn.features import (
    length,
    aac,
    aaindex1,
    ngram,
    entropy,
    posrich,
    motif,
    atc,
    binary,
    cksaap,
    ctd,
    ctdc,
    ctdt,
    ctdd,
    moreau_broto,
    moran,
    geary,
    paac,
    apaac,
    socn,
    qso,
)

FEATURE_FUNCTIONS = {
    'length': length,
    'aac': aac,
    'aaindex1': aaindex1,
    'ngram': ngram,
    'entropy': entropy,
    'posrich': posrich,
    'motif': motif,
    'atc': atc,
    'binary': binary,
    'cksaap': cksaap,
    'ctd': ctd,
    'ctdc': ctdc,
    'ctdt': ctdt,
    'ctdd': ctdd,
    'moreau_broto': moreau_broto,
    'moran': moran,
    'geary': geary,
    'paac': paac,
    'apaac': apaac,
    'socn': socn,
    'qso': qso,
}

def extract_sequence_from_pdb(pdb_file):
    """Extract amino acid sequence (1-letter) from a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    resname = residue.resname.capitalize()
                    if resname in protein_letters_3to1:
                        seq.append(protein_letters_3to1[resname])
                    else:
                        # Skip unknown/non-standard residues
                        continue
    return ''.join(seq)

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
            # Fallback for non-standard outputs
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
    """Main entry point: extract features from PDB and write to CSV."""
    output_dir = os.path.join(output_dir, "protlearn")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{input.stem}.csv")

    sequence = extract_sequence_from_pdb(input_file)
    features = extract_all_features(sequence)

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["feature", "value"])
        for k, v in features.items():
            writer.writerow([k, v])

    print(f"Extracted {len(features)} features from {input_file} → {output_file}")
    return output_file
