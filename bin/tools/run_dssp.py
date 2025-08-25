"""
run_dssp.py
Run DSSP on a given structure file using the DSSP API.


Author: Nadia
"""
import logging
import requests
import time
from pathlib import Path

def run(input_file, tool_root, output_dir):
    """
    Run DSSP on the specified input file after necessary conversions.
    Via the API of DSSP
    Args:
        input_file (str or Path): Path to the input structure file (PDB, CIF, or CIF.GZ).
        tool_root (str or Path): Unused
        output_dir (str or Path): Directory to save the DSSP output.
    """
    input_file = Path(input_file)
    output_dir = Path(output_dir)/"dssp"
    output_dir.mkdir(parents=True, exist_ok=True)

    dssp_URL = "https://pdb-redo.eu/dssp/do"
    """
    Using the DSSP API

    Using the DSSP API is as simple as doing a REST call to https://pdb-redo.eu/dssp/do using POST containing a parameter called data containing the mmCIF or PDB formatted structure and optionally a parameter called format that can contain either mmcif or dssp. The result will be an HTTP reply containing only the output from mkdssp.

    """

    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = {"data": input_file.read_text(), "format": "dssp"}

    try:
        response = requests.post(dssp_URL, headers=headers, data=data)
        response.raise_for_status()

        output_file = output_dir / f"{input_file.stem}.dssp"
        output_file.write_text(response.text)
        logging.info(f"✅ DSSP completed: {output_file.name}")
        wait_time = 1  # Set a wait time (in seconds) before the next request
        time.sleep(wait_time)
    except requests.RequestException as e:
        logging.error(f"❌ DSSP failed: {input_file.name}")
        logging.error(e)
