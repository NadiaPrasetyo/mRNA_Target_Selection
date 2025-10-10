"""
run_discotope.py
Runner for DiscoTope epitope prediction tool.
Overview:
    - Applies the DiscoTope algorithm to predict discontinuous B-cell epitopes from a given PDB structure.
    - Uses a Conda environment ('discotope_env') to ensure all dependencies are met.
    - Saves prediction results in a structured directory format for downstream analysis.
Arguments:
    input_file (Path): Path to the input PDB file (e.g., /path/to/7c4s.pdb).
    tool_path (Path): Path to the DiscoTope CLI script (e.g., /base/discotope/src/predict_webserver.py).
    output_dir (Path): Directory where prediction results will be saved.
Requirements:
    - Python packages: subprocess, pathlib, logging, warnings, os.
    - A Conda environment named 'discotope_env' must be created and configured with the required dependencies.
    - The DiscoTope tool must be installed and accessible.
Outputs:
    <output_dir>/discotope/<input_file_stem>/           # Directory containing prediction results.
Notes:
    - This script automatically determines the structure type ('alphafold' or 'solved') based on the input file name.
    - Logs are generated to provide detailed information about the execution process.
    - Ensure that the Conda environment is activated and accessible before running the script.
Author: Nadia
"""
import os
import subprocess
import logging
from pathlib import Path
from tools import common
import warnings

def run(input_file, tool_path, output_dir):
    """
    Run DiscoTope on a single PDB structure.

    Args:
        input_file : Path to the input PDB file (e.g., /path/to/7c4s.pdb)
        tool_root : Path to the DiscoTope installation directory (where src/predict_webserver.py is located)
        output_dir : Directory to save the prediction results
    """
    
    # 🔇 Suppress the noisy pkg_resources deprecation warning from XGBoost
    warnings.filterwarnings(
        "ignore",
        message="pkg_resources is deprecated as an API",
        category=UserWarning,
        module="xgboost.compat"
    )

    input_file = Path(input_file)
    output_dir = Path(output_dir)/"discotope"/input_file.stem
    tool_path = Path(tool_path)  # e.g., /base/discotope/src/predict_webserver.py

    # Ensure paths exist
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input PDB file not found: {input_file}")

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    parts = input_file.stem.split("_")
    if len(parts) > 1 and parts[1].lower() == "af":
        struct_type = "alphafold"
    else:
        struct_type = "solved"

    # Build the command
    cmd = [
        "conda", "run", "-n", common.DISCOTOPE_ENV_NAME,
        "python", str(tool_path),
        "--models_dir", str(tool_path.parent.parent/"models"),  # e.g., /base/discotope/models
        "--cpu_only",
        "--pdb_or_zip_file", str(input_file),
        "--struc_type", struct_type,
        "--out_dir", output_dir
    ]

    logging.info(f"Running DiscoTope on file: {input_file} with struc_type: {struct_type}")

    # Execute the command
        # Execute the command
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"✅ DiscoTope finished successfully. Results saved in: {output_dir}")

    except subprocess.CalledProcessError as e:
        # Capture stderr for diagnosis
        err_msg = e.stderr if hasattr(e, "stderr") and e.stderr else str(e)
        logging.warning(f"⚠️ DiscoTope failed on {input_file.name} (exit code {e.returncode}). Skipping.\nDetails:\n{err_msg}")

        # Mark the output directory as failed for traceability
        fail_marker = Path(output_dir) / "FAILED.txt"
        with open(fail_marker, "w") as f:
            f.write(f"DiscoTope failed on {input_file}\nError code: {e.returncode}\n\n{err_msg}")

        # Do NOT raise the error — just skip
        return None
