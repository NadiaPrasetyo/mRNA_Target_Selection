"""
run_dssp.py
Utility to execute PyDSSP secondary structure assignment on a given PDB file.

Overview:
    - Supports input: .pdb
    - Extracts backbone coordinates (N, CA, C, O) for secondary structure assignment
    - Runs PyDSSP to assign secondary structure
    - Stores output in specified directory as a .dssp file

Arguments:
    input_file (str or Path): PDB file to analyze.
    tool_root (str or Path): Unused, included for API compatibility.
    output_dir (str or Path): Output directory for DSSP results.

Requirements:
    - PyDSSP must be installed and available in the Python environment
    - Input PDB file must contain residues with complete backbone atoms (N, CA, C, O)

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
    """
    Extract backbone coordinates (N, CA, C, O) for PyDSSP.
    Skips residues missing any backbone atom.
    Args:
        pdb_file (Path): Path to the input PDB file.

    Returns:
        np.ndarray: Array of backbone coordinates (N, CA, C, O) for each residue.
    """
    structure = gemmi.read_structure(str(pdb_file))
    model = structure[0]  # first model
    coords = []

    for chain in model:
        for res in chain:
            atoms = []
            missing_atom = False
            for atom_name in ["N", "CA", "C", "O"]:
                try:
                    atom_group = res[atom_name]  # may raise RuntimeError if missing
                except RuntimeError:
                    missing_atom = True
                    logging.warning(
                        "Skipping residue %s%d in chain %s: missing atom %s",
                        res.name, res.seqid.num, chain.name, atom_name
                    )
                    break

                if not atom_group:
                    missing_atom = True
                    logging.warning(
                        "Skipping residue %s%d in chain %s: empty AtomGroup for %s",
                        res.name, res.seqid.num, chain.name, atom_name
                    )
                    break

                atom = atom_group[0]  # pick first atom if multiple altlocs
                pos = atom.pos
                atoms.append([pos.x, pos.y, pos.z])

            if not missing_atom and len(atoms) == 4:
                coords.append(atoms)

    coords = np.array(coords, dtype=np.float32)
    if coords.size == 0:
        raise RuntimeError(f"No residues with complete backbone found in {pdb_file}")
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
