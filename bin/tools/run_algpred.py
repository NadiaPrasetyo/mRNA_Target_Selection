"""
run_algpred.py
Run the AlgPred2.0 allergenicity prediction tool via Conda environment.

Compatible with tool_runners[tool] pattern and concurrent.futures.

Requirements:
    - algpred2_dependencies.yml (defines Conda env)
    - pip-installable `algpred2` package inside that environment

Author: Nadia
"""

import subprocess
import logging
from pathlib import Path
import shutil

logging.basicConfig(level=logging.INFO)

CONDA_ENV_NAME = "algpred2_env"
CONDA_ENV_YML = Path("algpred2_dependencies.yml")

def create_conda_env_if_needed():
    """Create Conda environment if it doesn't exist."""
    logging.info(f"🔍 Checking for Conda environment '{CONDA_ENV_NAME}'...")
    result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if CONDA_ENV_NAME not in result.stdout:
        logging.info("📦 Conda environment not found. Creating from YAML...")
        subprocess.run(["conda", "env", "create", "-f", str(CONDA_ENV_YML)], check=True)
    else:
        logging.info("✅ Conda environment already exists.")

def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    """
    Main runner function compatible with pipeline:
    - tool_path: directory containing tools (unused here but kept for interface consistency)
    - input_fasta: input FASTA file path
    - output_dir: base output directory (tool-specific subdir will be created)
    """
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    create_conda_env_if_needed()

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve() / "algpred"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "result.csv"

    cmd = [
        "conda", "run", "-n", CONDA_ENV_NAME,
        "algpred2",
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "2",  # hybrid-based model
        "-d", "1"   # database features enabled
    ]

    logging.info(f"🚀 Running AlgPred2.0 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ AlgPred2.0 finished: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ AlgPred2.0 failed: {e}")
        raise
