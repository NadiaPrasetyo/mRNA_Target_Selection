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

def write_mkdssp_compatible_pdb(input_file):
    """
    Parse input PDB and write a temporary file with a mkdssp-compatible header.
    """
    import tempfile
    from Bio.PDB import PDBParser, PDBIO

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(input_file))

    temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")

    # mkdssp-compatible header
    header_lines = [
        "HEADER    DSSP GENERATED\n",
        "TITLE     DSSP CLEAN PDB\n",
        "COMPND    MOL_ID: 1;\n",
        "SOURCE    MOL_ID: 1; ORGANISM_SCIENTIFIC: UNKNOWN;\n",
        "KEYWDS    DSSP\n",
        "EXPDTA    X-RAY DIFFRACTION\n",
        "AUTHOR    GENERATED\n"
    ]
    temp_pdb.writelines(header_lines)

    # Write ATOM/HETATM
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb.name, write_end=True)

    temp_pdb.close()
    return temp_pdb.name


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
        strict_pdb_path = write_mkdssp_compatible_pdb(input_file)

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
