
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

def patch_ifnepitope2_if_needed(env_path: Path):
    """
    Patch the installed ifnepitope2.py script to fix the 'composition' UnboundLocalError bug.
    Only applies if not already patched.
    """
    script_path = env_path / "lib/python3.10/site-packages/ifnepitope2/python_scripts/ifnepitope2.py"
    if not script_path.exists():
        logging.warning("⚠️ Could not find ifnepitope2.py to patch.")
        return

    with open(script_path) as f:
        content = f.read()

    patch_marker = "# === PATCHED FOR composition bug ==="
    if patch_marker in content:
        logging.info("🛠️ ifnepitope2 already patched.")
        return

    patched_lines = []
    for line in content.splitlines():
        # Patch the part where 'composition' might be undefined
        if "cc.append(composition)" in line:
            patched_lines.append(f"{patch_marker}")
            patched_lines.append("        if 'composition' not in locals(): continue")
        patched_lines.append(line)

    with open(script_path, "w") as f:
        f.write("\n".join(patched_lines))

    logging.info("✅ Patched ifnepitope2.py for unbound 'composition' bug.")


def run(tool_path: Path, input_fasta: Path, output_dir: Path, job_type: int = 1):
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    create_conda_env_if_needed()

    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_fasta.stem}_ifnepitope2.csv"

    # Apply patch only for job_type=1
    # PATCH: Fix UnboundLocalError in ifnepitope2.py (only for job_type == 1)
    if job_type == 1:
        try:
            # Get the actual file path inside the conda environment
            target_py_str = subprocess.check_output([
                "conda", "run", "-n", CONDA_ENV_NAME,
                "python", "-c",
                "import ifnepitope2.python_scripts.ifnepitope2 as m; print(m.__file__)"
            ], text=True).strip()

            target_py = Path(target_py_str)

            if not target_py.exists():
                raise FileNotFoundError(f"Resolved path {target_py} does not exist.")

            with open(target_py) as f:
                lines = f.readlines()

            already_patched = any("if 'composition' in locals()" in line for line in lines)

            if not already_patched:
                new_lines = []
                for line in lines:
                    if line.strip() == "cc.append(composition)":
                        indent = " " * (len(line) - len(line.lstrip()))
                        new_lines.append(f"{indent}if 'composition' in locals():\n")
                        new_lines.append(f"{indent}    cc.append(composition)\n")
                    else:
                        new_lines.append(line)

                with open(target_py, "w") as f:
                    f.writelines(new_lines)

                logging.info("✅ Patched ifnepitope2.py for unbound 'composition' bug.")
            else:
                logging.info("🔁 ifnepitope2.py already patched.")

        except Exception as e:
            logging.warning(f"⚠️ Could not patch ifnepitope2.py: {e}")

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
