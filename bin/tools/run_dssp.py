"""
run_dssp.py
Run DSSP on a given structure file using the local DSSP binary (conda hcc::dssp).

Author: Nadia
"""
import logging
import subprocess
from pathlib import Path
from tools import common


def run(input_file, tool_root, output_dir):
    """
    Run DSSP on the specified input file using the local mkdssp binary.

    Args:
        input_file (str or Path): Path to the input structure file (PDB, CIF, or CIF.GZ).
        tool_root (str or Path): Unused, kept for API consistency.
        output_dir (str or Path): Directory to save the DSSP output.
    """

    common.create_conda_env_if_needed()
    
    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{input_file.stem}.dssp"

    try:
        # Run mkdssp with positional input/output arguments
        subprocess.run(
            ["mkdssp", str(input_file), str(output_file)],
            check=True,
            capture_output=True,
            text=True
        )
        logging.info(f"✅ DSSP completed: {output_file.name}")
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ DSSP failed: {input_file.name}")
        logging.error(e.stderr)
        return False

    except FileNotFoundError:
        logging.error("❌ mkdssp binary not found. Make sure your conda env with hcc::dssp is activated.")
        return False
