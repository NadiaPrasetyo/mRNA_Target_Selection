import os
import subprocess
import shutil
import logging
from pathlib import Path
from tools import common

def run(input_file, tool_path, output_dir):
    """
    Run DiscoTope on a single PDB structure.

    Params:
    input_file : Path to the input PDB file (e.g., /path/to/7c4s.pdb)
    tool_root : Path to the DiscoTope installation directory (where src/predict_webserver.py is located)
    output_dir : Directory to save the prediction results
    """
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed()

    input_file = Path(input_file)
    output_dir = Path(output_dir)/"discotope"/input_file.stem
    tool_path = Path(tool_path)  # e.g., /path/to/discotope

    # Ensure paths exist
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input PDB file not found: {input_file}")

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    struct_type = "alphafold"  if input_file.stem.split("_")[1].lower() == "AF" else "solved"

    # Build the command
    cmd = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "python", str(tool_path),
        "--cpu_mode",
        "--pdb_or_zip_file", str(input_file),
        "--struc_type", struct_type,
        "--out_dir", output_dir
    ]

    print(f"Running DiscoTope:\n{' '.join(cmd)}")

    # Execute the command
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ DiscoTope finished successfully. Results saved in: {output_dir}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"❌ DiscoTope failed with error code {e.returncode}") from e
