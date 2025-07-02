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

CONDA_ENV_NAME = "algpred2_env"
CONDA_ENV_YML = Path("algpred2_dependencies.yml")

def patch_algpred_concat_bug():
    """
    Monkey-patch AlgPred2.0 source to fix .concat misuse.
    Only runs once if needed.
    """
    import sys
    import re

    logging.info("🩹 Checking AlgPred2.0 for known concat bug...")

    try:
        env_prefix = subprocess.run(
            ["conda", "run", "-n", CONDA_ENV_NAME, "python", "-c", "import sys; print(sys.prefix)"],
            capture_output=True, check=True, text=True
        ).stdout.strip()

        algpred_path = Path(env_prefix) / "lib/python3.10/site-packages/algpred2/python_scripts/algpred2.py"

        if not algpred_path.exists():
            logging.warning(f"⚠️ Cannot patch: {algpred_path} not found.")
            return

        with open(algpred_path, "r") as f:
            contents = f.read()

        if "df3.concat" not in contents:
            logging.info("✅ No patch needed; bug not present.")
            return

        patched_contents = re.sub(
            r"df3\.concat\((.*?)\)", 
            r"pd.concat([df3, \1])", 
            contents
        )

        if patched_contents != contents:
            logging.info("🔧 Patching algpred2.py to fix concat bug...")
            with open(algpred_path, "w") as f:
                f.write(patched_contents)
            logging.info("✅ Patch applied.")
        else:
            logging.info("✅ Patch not needed.")

    except Exception as e:
        logging.warning(f"⚠️ Failed to check or patch AlgPred2.0: {e}")

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
    patch_algpred_concat_bug()


    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_algpred.csv"

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
