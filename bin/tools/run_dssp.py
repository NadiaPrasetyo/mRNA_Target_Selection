"""
run_dssp.py
Run DSSP on a given structure file using Biopython's DSSP wrapper.

Author: Nadia
"""
import logging
from pathlib import Path
from tools import common
from Bio.PDB import PDBParser, DSSP
import os
import subprocess

def run(input_file, tool_root, output_dir):
    """
    Run DSSP on the specified PDB file using Biopython's DSSP wrapper.

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

    # Locate the mmcif dictionary inside conda
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        logging.error("❌ CONDA_PREFIX not set. Cannot locate DSSP dictionary.")
        return False
    dic_path = Path(conda_prefix) / "share/libcifpp/mmcif_pdbx.dic"

    # Monkey patch subprocess.run to include --mmcif-dictionary safely
    original_run = subprocess.run

    def patched_run(*popenargs, **kwargs):
        if isinstance(popenargs[0], list) and popenargs[0][0].endswith("mkdssp"):
            # Copy the command list and append the dictionary argument
            cmd = list(popenargs[0]) + ["--mmcif-dictionary", str(dic_path)]
            popenargs = (cmd, *popenargs[1:])
        return original_run(*popenargs, **kwargs)

    subprocess.run = patched_run

    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(input_file.stem, str(input_file))
        model = structure[0]  # first model

        # Run DSSP (mkdssp will receive --mmcif-dictionary)
        dssp = DSSP(model, str(input_file), dssp="mkdssp")

        # Write DSSP data
        with open(output_file, "w") as f:
            for key in dssp.keys():
                res_data = dssp[key]
                # key is a tuple (chain_id, (hetfield, resseq, icode))
                chain_id, res_info = key
                resseq = res_info[1]
                icode = res_info[2].strip()
                res_str = f"{chain_id}:{resseq}{icode}"
                line = res_str + "\t" + "\t".join(map(str, res_data)) + "\n"
                f.write(line)

        logging.info(f"✅ DSSP completed: {output_file.name}")
        return True

    except Exception as e:
        logging.error(f"❌ DSSP failed: {input_file.name}")
        logging.error(e)
        return False

    finally:
        # Restore original subprocess.run
        subprocess.run = original_run
