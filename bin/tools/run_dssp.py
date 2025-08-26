"""
run_dssp.py
Parse PDB with Biopython, preserve headers, then run mkdssp with --mmcif-dictionary.

Author: Nadia
"""
import logging
from pathlib import Path
from tools import common
from Bio.PDB import PDBParser, PDBIO
import os
import subprocess
import tempfile

def write_clean_pdb_with_header(input_file):
    """
    Parse a PDB with Biopython, preserve header lines, and write a temporary PDB
    suitable for mkdssp.

    Args:
        input_file (str or Path): Original PDB file path.

    Returns:
        str: Path to the temporary PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(input_file.stem, str(input_file))

    # Read header lines from original PDB
    header_lines = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith(("HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR")):
                header_lines.append(line)
            else:
                break  # stop at first ATOM line

    # Write temporary PDB combining header + cleaned coordinates
    temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    # write header
    for line in header_lines:
        temp_pdb.write(line)

    # write ATOM/HETATM via PDBIO
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb.name, write_end=True)

    temp_pdb.close()
    return temp_pdb.name

def run(input_file, tool_root, output_dir):
    """
    Run mkdssp on a cleaned PDB file prepared by Biopython.

    Args:
        input_file (str or Path): Path to the input PDB structure file.
        tool_root (str or Path): Unused, kept for API consistency.
        output_dir (str or Path): Directory to save the DSSP output.
    """
    common.create_conda_env_if_needed()

    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_file.stem}.dssp"

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        logging.error("❌ CONDA_PREFIX not set. Cannot locate DSSP dictionary.")
        return False
    dic_path = Path(conda_prefix) / "share/libcifpp/mmcif_pdbx.dic"

    temp_pdb_path = None
    try:
        # 1. Write cleaned PDB with headers
        temp_pdb_path = write_clean_pdb_with_header(input_file)

        # 2. Call mkdssp directly on the cleaned PDB
        cmd = [
            "mkdssp",
            "--mmcif-dictionary", str(dic_path),
            temp_pdb_path,
            str(output_file), "--verbose"
        ]
        subprocess.run(cmd, check=True)

        logging.info(f"✅ DSSP completed: {output_file.name}")
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ mkdssp failed: {input_file.name}")
        logging.error(e)
        return False
    except Exception as e:
        logging.error(f"❌ DSSP processing failed: {input_file.name}")
        logging.error(e)
        return False
    finally:
        # Clean up temporary file
        if temp_pdb_path and Path(temp_pdb_path).exists():
            try:
                os.remove(temp_pdb_path)
            except Exception:
                pass
