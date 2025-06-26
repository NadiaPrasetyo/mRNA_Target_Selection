"""
run_algpred.py
Run the AlgPred2.0 allergenicity prediction tool via a specified conda environment.

Requirements:
    - Conda environment defined in `algpred2_dependencies.yml`
    - `algpred2` installed via pip inside that environment

Author: Nadia (refined by ChatGPT)
"""

import subprocess
import logging
from pathlib import Path
import shutil
import sys

CONDA_ENV_NAME = "algpred2_env"
CONDA_ENV_YML = Path("algpred2_dependencies.yml")

def create_conda_env():
    """Create the Conda environment if it doesn't exist."""
    logging.info(f"🔍 Checking for Conda environment '{CONDA_ENV_NAME}'...")
    result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if CONDA_ENV_NAME not in result.stdout:
        logging.info("📦 Creating Conda environment...")
        subprocess.run(["conda", "env", "create", "-f", str(CONDA_ENV_YML)], check=True)
    else:
        logging.info("✅ Conda environment already exists.")

def run_algpred(input_fasta: Path, output_dir: Path):
    """Run the algpred2 tool using the Conda environment."""
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve() / "Allergenicity"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}.csv"

    # Command to be run inside the Conda environment
    # Use method 2 (hybrid method) with -m 2
    algpred_cmd = (
        f"conda run -n {CONDA_ENV_NAME} "
        f"algpred2 "
        f"-i \"{input_fasta}\" -o \"{output_file}\" -m 2 -d 1"
    )

    logging.info("🚀 Running AlgPred2.0 inside Conda environment...")
    try:
        subprocess.run(algpred_cmd, shell=True, check=True)
        logging.info(f"✅ Finished. Output saved to: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ AlgPred2.0 failed with return code {e.returncode}")
        raise

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_algpred.py <input_fasta> <output_dir>")
        sys.exit(1)

    input_fasta = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not shutil.which("conda"):
        logging.error("❌ Conda not found in PATH. Please install Conda first.")
        sys.exit(1)

    create_conda_env()
    run_algpred(input_fasta, output_dir)

if __name__ == "__main__":
    main()
