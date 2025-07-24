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
import csv

# Constants for B-cell epitope prediction methods
BCELL_METHODS = [
    "Bepipred"
    # exclude "BepiPred-2.0" as it requires additional dependencies that are not supported anymore
]

def patch_numpy_float(util_path: Path):
    """
    Patch util.py to replace deprecated np.float with float.
    This is necessary for compatibility with newer versions of NumPy.
    Args:
        util_path (Path): Path to the util.py file in the B-cell prediction tool.
    """
    if not util_path.exists():
        print(f"util.py not found at {util_path}, skipping patch.")
        return

    content = util_path.read_text()
    if "np.float(" not in content:
        print("No deprecated np.float usage found in util.py, skipping patch.")
        return

    print("Patching util.py to replace np.float with float...")

    patched = content.replace("np.float(", "float(")

    backup_path = util_path.with_suffix(".py.bak")
    if not backup_path.exists():
        util_path.rename(backup_path)
        print(f"Backed up original util.py to {backup_path}")

    util_path.write_text(patched)
    print("Patch applied successfully.")

def patch_configure(configure_path: Path):
    """
    Patch configure.py to fix deprecated pip import.
    This is necessary for compatibility with newer versions of pip.
    Args:
        configure_path (Path): Path to the configure.py file in the B-cell prediction tool.
    """
    if not configure_path.exists():
        print(f"configure.py not found at {configure_path}, skipping patch.")
        return

    text = configure_path.read_text()

    if "from pip._internal.utils.misc import get_installed_distributions" not in text:
        print("No deprecated import found in configure.py, skipping patch.")
        return

    print("Patching configure.py to fix pip import issue...")

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

def parse_and_save_to_csv(output: str, output_file: Path):
    """
    Parse the output of the B-cell prediction tool and save it to a CSV file.
    Args:
        output (str): The raw output from the B-cell prediction tool.
        output_file (Path): The path where the CSV file will be saved.
    """
    lines = output.strip().splitlines()
    
    # Find the header (might not always be on the same line)
    for i, line in enumerate(lines):
        if line.startswith("Position"):
            header = line.split()
            data_lines = lines[i+1:]
            break
    else:
        print(f"No result header found in output.")
        return

    with output_file.open("w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for line in data_lines:
            if line.strip():
                writer.writerow(line.split())

def run(fasta_file: Path, tool_path: str, output_dir: Path, plot: bool = True):
    """
    Run B-cell epitope prediction methods on a given FASTA file.
    Args:
        fasta_file (Path): Path to the input FASTA file containing protein sequences.
        tool_path (str): Path to the main B-cell prediction tool script (e.g., bcell.py).
        output_dir (Path): Directory where results and plots will be saved.
        plot (bool): Whether to generate and save plots for each method (default: True).
    """
    fasta_stem = fasta_file.stem
    output_dir = output_dir / "bcell"
    output_dir.mkdir(parents=True, exist_ok=True)


    # Patch deprecated imports
    configure_py_path = Path(tool_path).parent / "configure.py"
    patch_configure(configure_py_path)

    util_py_path = Path(tool_path).parent / "src" / "util.py"
    patch_numpy_float(util_py_path)

    for method in BCELL_METHODS:
        print(f"Running BCell method: {method}...")
        cmd = [
            "python3",
            tool_path,
            "-m", method,
            "-f", str(fasta_file)
        ]

        if plot:
            plot_dir = output_dir / "plots"
            plot_dir.mkdir(parents=True, exist_ok=True)
            cmd += ["--plot", str(plot_dir)]


        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            print(f"✅ BCell [{method}] completed for {fasta_file.name}")
            
            # Save output to CSV
            csv_path = output_dir / f"{fasta_stem}_{method.replace(' ', '_')}.csv"
            parse_and_save_to_csv(result.stdout, csv_path)
            print(f"📄 Results saved to: {csv_path}")

        except subprocess.CalledProcessError as e:
            print(f"❌ BCell [{method}] failed for {fasta_file.name}")
            print("Error:", e.stderr)
