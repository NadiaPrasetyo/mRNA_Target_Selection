import os
import subprocess
import shutil
import logging
from pathlib import Path
from tools import common
import warnings

def run(input_file, tool_path, output_dir):
    """
    Run DiscoTope on a single PDB structure.

    Params:
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

    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed(common.DISCOTOPE_ENV_NAME, common.DISCOTOPE_ENV_YML)

    input_file = Path(input_file)
    output_dir = Path(output_dir)/"discotope"/input_file.stem
    tool_path = Path(tool_path)  # e.g., /base/discotope/src/predict_webserver.py

    # Ensure paths exist
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input PDB file not found: {input_file}")

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    struct_type = "alphafold"  if input_file.stem.split("_")[1].lower() == "AF" else "solved"

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
    try:
        result = subprocess.run(cmd, check=True)
        logging.debug(f"DiscoTope command output:\n{result.stdout}")
        logging.debug(f"DiscoTope command error output:\n{result.stderr}")

        logging.info(f"✅ DiscoTope finished successfully. Results saved in: {output_dir}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"❌ DiscoTope failed with error code {e.returncode}") from e
