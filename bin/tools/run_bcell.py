import subprocess
from pathlib import Path

BCELL_METHODS = [
    "Chou-Fasman",
    "Emini",
    "Karplus-Schulz",
    "Kolaskar-Tongaonkar",
    "Parker",
    "Bepipred",
    "Bepipred-2.0"
]

def patch_configure(configure_path: Path):
    """
    Patch configure.py to replace deprecated pip import and usage.
    This replaces:
    - from pip._internal.utils.misc import get_installed_distributions
      with
      import pkg_resources
    - get_installed_distributions() calls with pkg_resources.working_set
    """
    if not configure_path.exists():
        print(f"configure.py not found at {configure_path}, skipping patch.")
        return

    # Read content
    text = configure_path.read_text()

    if "from pip._internal.utils.misc import get_installed_distributions" not in text:
        print("No deprecated import found in configure.py, skipping patch.")
        return

    print("Patching configure.py to fix pip import issue...")

    # Replace import line
    text = text.replace(
        "from pip._internal.utils.misc import get_installed_distributions",
        "import pkg_resources"
    ).replace(
        "get_installed_distributions()",
        "pkg_resources.working_set"
    )

    backup_path = configure_path.with_suffix(".py.bak")
    if not backup_path.exists():
        configure_path.rename(backup_path)
        print(f"Backed up original configure.py to {backup_path}")

    configure_path.write_text(text)
    print("Patch applied successfully.")

def run(fasta_file: Path, tool_path: str, output_dir: Path, plot: bool = True):
    fasta_stem = fasta_file.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    # Patch deprecated imports if needed
    configure_py_path = Path(tool_path).parent / "configure.py"
    patch_configure(configure_py_path)

    for method in BCELL_METHODS:
        print(f"Running BCell method: {method}...")
        cmd = [
            "python3",
            tool_path,
            "-m", method,
            "-f", str(fasta_file)
        ]

        if plot:
            # Specify plot directory
            cmd += ["--plot", str(output_dir)]

        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            print(f"✅ BCell [{method}] completed for {fasta_file.name}")
        except subprocess.CalledProcessError as e:
            print(f"❌ BCell [{method}] failed for {fasta_file.name}")
            print("Error:", e.stderr)
