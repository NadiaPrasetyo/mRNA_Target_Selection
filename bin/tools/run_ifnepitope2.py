import logging
from pathlib import Path
import subprocess
import shutil
import re

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


def patch_ifnepitope2_composition_bug():
    """
    Locate the ifnepitope2.py source file inside the conda env,
    and patch the line with 'cc.append(composition)' to avoid UnboundLocalError
    by adding an 'if' guard with proper indentation.
    """
    try:
        # Get the path of the python source file inside the conda env
        target_py_str = subprocess.check_output([
            "conda", "run", "-n", CONDA_ENV_NAME,
            "python", "-c",
            "import ifnepitope2.python_scripts.ifnepitope2 as m; print(m.__file__)"
        ], text=True).strip()
        target_py = Path(target_py_str)

        if not target_py.exists():
            logging.warning(f"⚠️ ifnepitope2 source file not found at {target_py}")
            return

        with open(target_py) as f:
            lines = f.readlines()

        patched_lines = []
        patched = False

        for line in lines:
            if re.match(r"^\s*cc\.append\(composition\)", line) and not patched:
                indent = re.match(r"^(\s*)", line).group(1)
                # Insert guarded append with proper indentation
                patched_lines.append(f"{indent}if 'composition' in locals():\n")
                patched_lines.append(f"{indent}    cc.append(composition)\n")
                patched = True
            else:
                patched_lines.append(line)

        if patched:
            backup_path = target_py.with_suffix(".bak")
            shutil.copy(target_py, backup_path)
            with open(target_py, "w") as f:
                f.writelines(patched_lines)
            logging.info(f"✅ Patched {target_py.name} for 'composition' bug.")
        else:
            logging.info("🔁 Patch not applied: target line not found or already patched.")

    except Exception as e:
        logging.warning(f"⚠️ Could not patch ifnepitope2.py: {e}")


def run(tool_path: Path, input_fasta: Path, output_dir: Path, job_type: int = 1):
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    # Assume create_conda_env_if_needed() is defined elsewhere
    create_conda_env_if_needed()

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_ifnepitope2.csv"

    # Only patch if job_type is 1
    if job_type == 1:
        patch_ifnepitope2_composition_bug()

    cmd = [
        "conda", "run", "-n", CONDA_ENV_NAME,
        "ifnepitope2",
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-s", "1",              # host human
        "-j", str(job_type),    # job type
        "-d", "2"               # display mode: all peptides
    ]

    logging.info(f"🚀 Running IfNePitope2 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ IfNePitope2 finished: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ IfNePitope2 failed: {e}")
        raise