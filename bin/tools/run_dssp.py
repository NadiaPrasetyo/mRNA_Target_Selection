"""
run_dssp.py
Parse PDB with Biopython, then run mkdssp with --mmcif-dictionary.

Author: Nadia
"""
import logging
from pathlib import Path
from tools import common
from Bio.PDB import PDBParser, PDBIO
import os
import subprocess
import tempfile

def run(input_file, tool_root, output_dir):
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

    try:
        # 1. Parse PDB and write a clean temporary PDB
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(input_file.stem, str(input_file))
        temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        io = PDBIO()
        io.set_structure(structure)
        io.save(temp_pdb.name)
        temp_pdb.close()  # ensure it's written

        # 2. Call mkdssp directly on the cleaned PDB
        cmd = [
            "mkdssp",
            "--mmcif-dictionary", str(dic_path),
            temp_pdb.name,
            str(output_file)
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
        try:
            if temp_pdb:
                os.remove(temp_pdb.name)
        except Exception:
            pass
