"""
run_dssp.py
Use PyDSSP to assign secondary structure from a PDB file.

Author: Nadia
"""
import logging
from pathlib import Path
import torch
import gemmi
import numpy as np
import pydssp
from tools import common

def _extract_backbone_coords(pdb_file: Path) -> np.ndarray:
    """Extract backbone coordinates (N, CA, C, O) as a NumPy array for PyDSSP.
    Shape: (residues, atoms=4, xyz=3)
    """
    structure = gemmi.read_structure(str(pdb_file))
    model = structure[0]  # first model
    coords = []

    for chain in model:
        for res in chain:
            atoms = []
            for atom_name in ["N", "CA", "C", "O"]:
                atom = res[atom_name] if atom_name in res else None
                if atom is None:
                    break
                pos = atom.pos  # gemmi.Position object
                atoms.append([pos.x, pos.y, pos.z])
            if len(atoms) == 4:
                coords.append(atoms)
    coords = np.array(coords, dtype=np.float32)
    return coords

def run(input_file, tool_root, output_dir):
    """
    Run PyDSSP on input PDB file and write DSSP output to file.
    
    Args:
        input_file (str or Path): input PDB path
        tool_root (str or Path): unused, API compatibility
        output_dir (str or Path): directory to save DSSP output
    """

    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_file.stem}.dssp"

    try:
        coords = _extract_backbone_coords(input_file)

        if len(coords) == 0:
            logging.error(f"No complete backbone residues found in {input_file}")
            return False

        # PyDSSP expects (batch, residues, atoms, xyz); add batch dim
        coords_torch = torch.from_numpy(coords[None, ...])
        dssp_assignment = pydssp.assign(coords_torch, out_type='c3')[0]  # take first batch

        # Save DSSP assignment to text file, one line
        with open(output_file, 'w') as f:
            f.write("".join(dssp_assignment) + f"  {input_file.name}\n")

        logging.info(f"✅ DSSP completed: {output_file}")
        return True

    except Exception as e:
        logging.error(f"❌ PyDSSP failed for {input_file}")
        logging.exception(e)
        return False
