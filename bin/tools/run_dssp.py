"""
run_dssp.py
Parse PDB with Biopython, rewrite in strict PDB format using Biopandas, then run mkdssp.

Author: Nadia
"""
import logging
from pathlib import Path
from tools import common
from Bio.PDB import PDBParser
import os
import subprocess
import tempfile
from biopandas.pdb import PandasPdb

def write_strict_pdb(input_file):
    """
    Parse a PDB with Biopython, then rewrite it in strict PDB format using Biopandas.
    This ensures compatibility with mkdssp.

    Args:
        input_file (str or Path): Original PDB file path.

    Returns:
        str: Path to the temporary strictly formatted PDB file.
    """
    # Parse with Biopython to sanitize structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(input_file))

    # Save coordinates to temporary PDB using PDBIO
    temp_io = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    temp_io.close()
    from Bio.PDB import PDBIO
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_io.name)

    # Read with Biopandas and rewrite in strict PDB format
    ppdb = PandasPdb().read_pdb(temp_io.name)
    strict_temp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    strict_temp.close()
    ppdb.to_pdb(path=str(strict_temp.name), records=['ATOM', 'HETATM'], gz=False, append_newline=True)

    # Remove intermediate file
    os.remove(temp_io.name)
    return strict_temp.name

def run(input_file, tool_root, output_dir):
    """
    Run mkdssp on a cleaned, strictly formatted PDB.

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

    strict_pdb_path = None
    try:
        # 1. Write strictly formatted PDB
        strict_pdb_path = write_strict_pdb(input_file)

        # 2. Run mkdssp directly
        cmd = [
            "mkdssp",
            "--mmcif-dictionary", str(dic_path),
            strict_pdb_path,
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
        # Clean up temporary strict PDB
        if strict_pdb_path and Path(strict_pdb_path).exists():
            try:
                os.remove(strict_pdb_path)
            except Exception:
                pass
