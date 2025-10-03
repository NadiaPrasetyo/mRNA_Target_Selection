"""
run_bcell.py
Runner for B-cell epitope prediction tools.

Overview:
    - Applies B-cell epitope prediction algorithms (currently Bepipred) to a given FASTA file.
    - Automatically patches deprecated imports and code in third-party tool dependencies for compatibility.
    - Parses and saves prediction results as CSV files for downstream analysis.
    - Optionally generates and saves plots for each prediction method.

Arguments:
    fasta_file (Path): Path to the input FASTA file containing protein sequence(s).
    tool_path (str): Path to the main B-cell prediction tool script (e.g., bcell.py).
    output_dir (Path): Directory where results and plots will be saved.
    plot (bool, optional): Whether to generate and save plots for each method (default: True).

Requirements:
    - Python packages: subprocess, pathlib, csv.
    - The B-cell prediction tool and its dependencies must be installed and accessible.
    - The script will attempt to patch deprecated code in 'configure.py' and 'src/util.py' if needed.

Outputs:
    <output_dir>/bcell/<fasta_file_stem>_<method>.csv   # Prediction results for each method
    <output_dir>/bcell/plots/                           # Plots (if enabled)

Author: Nadia
"""
import subprocess
from pathlib import Path
import shutil
from tools import common
import logging

def run(fasta_file: Path, tool_path: str, output_dir: Path):
    """
    Run BepiPred-3.0 predictor on a given FASTA file using the specified conda environment.

    Parameters
    ----------
    fasta_file : Path
        Path to input FASTA file containing protein sequence(s).
    tool_path : str
        Path to the BepiPred-3.0 installation directory (must contain bepipred3_CLI.py).
    output_dir : Path
        Path to the output directory where prediction results will be stored.

    Notes
    -----
    - This script assumes that a Conda environment named 'external_tools_env'
      is already created and contains all dependencies listed in requirements.txt.
    - It uses default prediction mode (vt_pred) and default thresholds.
    """
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed(common.EXT_TOOLS_ENV_NAME, common.EXT_TOOLS_ENV_YML)

    # Validate input
    fasta_file = Path(fasta_file)
    output_dir = Path(output_dir)/"bcell"/fasta_file.stem
    tool_path = Path(tool_path)
    # tool path would be like: /path/to/BepiPred3_src/bepipred3_CLI.py

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "python", str(tool_path),
        "-i", str(fasta_file),
        "-o", str(output_dir),
        "-pred", "vt_pred"
    ]

    logging.info(f"[INFO] Running BepiPred-3.0:\n{cmd}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logging.info("STDOUT:\n%s", result.stdout)
        logging.info("STDERR:\n%s", result.stderr)
        
        logging.info(f"[INFO] BepiPred-3.0 run completed successfully. results saved to {output_dir}")
        logging.info(result.stdout)

    except subprocess.CalledProcessError as e:
        logging.error("BepiPred-3.0 failed with exit code %s", e.returncode)
        logging.error("STDOUT:\n%s", e.stdout)
        logging.error("STDERR:\n%s", e.stderr)
        raise      # re-raise or handle as needed


