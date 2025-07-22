
import logging
from pathlib import Path
import subprocess
import shutil

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

def run(tool_path: Path, input_fasta: Path, output_dir: Path, job_type: int = 1):
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
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_ifnepitope2.csv"

    cmd = [
        "conda", "run", "-n", CONDA_ENV_NAME,
        "ifnepitope2",
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-s", "1",  # host human
        "-j", str(job_type),  # job type: 1 for prediction
        "-d", "2"   # display mode 2: all peptides (not just allergens)
        # use default threshold of 0.49 and window lenght of 8
    ]

    logging.info(f"🚀 Running IfNePitope2 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ IfNePitope2 finished: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ IfNePitope2 failed: {e}")
        raise