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


def patch_ifnepitope2_bugfixes():
    """
    Patches the ifnepitope2.py script inside the conda environment:
    - Fixes 'composition' UnboundLocalError.
    - Adds filtering of NaN values before model prediction (for job_type == 1 only).
    """
    try:
        result = subprocess.run([
            "conda", "run", "-n", CONDA_ENV_NAME,
            "python", "-c",
            "import os, ifnepitope2.python_scripts; "
            "print(os.path.join(os.path.dirname(ifnepitope2.python_scripts.__file__), 'ifnepitope2.py'))"
        ], capture_output=True, text=True, check=True)

        target_py = Path(result.stdout.strip())
        if not target_py.exists():
            logging.warning(f"⚠️ Resolved ifnepitope2.py does not exist: {target_py}")
            return

        with open(target_py) as f:
            lines = f.readlines()

        patched_lines = []
        composition_patched = False
        nan_patched = False

        for i, line in enumerate(lines):
            # Patch 1: Prevent UnboundLocalError for 'composition'
            if re.match(r"^\s*cc\.append\(composition\)", line) and not composition_patched:
                indent = re.match(r"^(\s*)", line).group(1)
                patched_lines.append(f"{indent}if 'composition' in locals():\n")
                patched_lines.append(f"{indent}    cc.append(composition)\n")
                composition_patched = True
                continue

            # Patch 2: Skip NaN rows before prediction
            if "y_p_score1 = clf.predict_proba(data_test)" in line and not nan_patched:
                indent = re.match(r"^(\s*)", line).group(1)
                patched_lines.append(f"{indent}import numpy as np\n")
                patched_lines.append(f"{indent}if np.isnan(data_test).any():\n")
                patched_lines.append(f"{indent}    print('⚠️ Warning: Found NaNs in feature matrix, skipping those peptides.')\n")
                patched_lines.append(f"{indent}    data_test = data_test[~np.isnan(data_test).any(axis=1)]\n")
                patched_lines.append(f"{indent}    if len(data_test) == 0:\n")
                patched_lines.append(f"{indent}        print('❌ All peptides led to invalid feature vectors. Exiting.')\n")
                patched_lines.append(f"{indent}        return []\n")
                nan_patched = True

            patched_lines.append(line)

        if composition_patched or nan_patched:
            backup_path = target_py.with_suffix(".bak")
            shutil.copy(target_py, backup_path)
            with open(target_py, "w") as f:
                f.writelines(patched_lines)
            logging.info("✅ Applied patch(es) to ifnepitope2.py:")
            if composition_patched:
                logging.info("  - ✔️ Fixed 'composition' bug.")
            if nan_patched:
                logging.info("  - ✔️ Added NaN feature filtering.")
        else:
            logging.info("🔁 Patch not applied: target lines not found or already patched.")

    except subprocess.CalledProcessError as e:
        logging.warning(f"⚠️ Could not locate ifnepitope2.py in conda env: {e}")
    except Exception as e:
        logging.warning(f"⚠️ Could not patch ifnepitope2.py: {e}")


def run(tool_path: Path, input_fasta: Path, output_dir: Path, job_type: int = 1):
    """Run ifnepitope2 prediction tool from within the conda environment."""
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    create_conda_env_if_needed()

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_ifnepitope2.csv"

    if job_type == 1:
        patch_ifnepitope2_bugfixes()

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
