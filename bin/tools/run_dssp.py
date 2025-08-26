"""
run_dssp.py
Run DSSP on a given structure file using the DSSP API.


Author: Nadia
"""
import logging
import requests
import time
from pathlib import Path

def run(input_file, tool_root, output_dir, max_retries=5):
    
    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)

    dssp_URL = "https://pdb-redo.eu/dssp/do"

    # Read file content
    pdb_text = input_file.read_text()

    data = {
        "data": pdb_text,
        "format": "dssp"
    }

    retries = 0
    wait_time = 2  # start with 2s
    while retries < max_retries:
        try:
            response = requests.post(dssp_URL, data=data, timeout=60)
            response.raise_for_status()

            output_file = output_dir / f"{input_file.stem}.dssp"
            output_file.write_text(response.text)
            logging.info(f"✅ DSSP completed: {output_file.name}")
            time.sleep(1)  # polite delay
            return True
        except requests.RequestException as e:
            retries += 1
            logging.warning(
                f"⚠️ DSSP attempt {retries} failed for {input_file.name}: {e}"
            )
            time.sleep(wait_time)
            wait_time *= 2  # exponential backoff

    logging.error(f"❌ DSSP failed after {max_retries} retries: {input_file.name}")
    return False
