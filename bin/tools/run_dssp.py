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
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)

    dssp_URL = "https://pdb-redo.eu/dssp/do"

    # Read file content
    pdb_text = input_file.read_text()

    # Do NOT set headers manually – requests will handle it
    data = {
        "data": pdb_text,   # structure file as text
        "format": "dssp"    # optional, defaults to dssp
    }

    try:
        response = requests.post(dssp_URL, data=data)
        response.raise_for_status()

        output_file = output_dir / f"{input_file.stem}.dssp"
        output_file.write_text(response.text)
        logging.info(f"✅ DSSP completed: {output_file.name}")

        time.sleep(1)  # polite delay to avoid hammering server
    except requests.RequestException as e:
        logging.error(f"❌ DSSP failed: {input_file.name}")
        logging.error(e)
