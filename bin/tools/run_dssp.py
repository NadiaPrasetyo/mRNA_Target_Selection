"""
run_dssp.py
Run DSSP via ssbio wrapper (no --mmcif-dictionary needed).

Author: Nadia
"""
import logging
from pathlib import Path
import pandas as pd
from tools import common
from ssbio.protein.structure import properties as dssp_props


def run(input_file, tool_root, output_dir):
    """
    Run DSSP using ssbio wrapper (no mmCIF conversion or dictionary needed).

    Args:
        input_file (str or Path): Path to the input PDB or mmCIF file.
        tool_root (str or Path): Unused, kept for API consistency.
        output_dir (str or Path): Directory to save the DSSP output.

    Returns:
        bool: True if DSSP succeeded and output was written, False otherwise.
    """
    common.create_conda_env_if_needed()

    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_file.stem}.dssp.csv"

    try:
        logging.info(f"▶ Running DSSP via ssbio on {input_file.name}")

        # Get DSSP dataframe
        dssp_df = dssp_props.dssp.get_dssp_df_on_file(
            str(input_file),
            outfile=None,
            outdir=str(output_dir),
            force_rerun=True,
        )

        if dssp_df is None or len(dssp_df) == 0:
            logging.error(f"❌ DSSP returned no results for {input_file.name}")
            return False

        # Save as CSV for downstream tools
        dssp_df.to_csv(output_file)
        logging.info(f"✅ DSSP completed: {output_file.name}")
        return True

    except Exception as e:
        logging.error(f"❌ DSSP processing failed: {input_file.name}")
        logging.exception(e)
        return False
