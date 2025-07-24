"""
run_mhci.py
Utility to execute MHCI epitope prediction tool on a given JSON input file.

Overview:
    - Runs an external MHCI prediction script using subprocess.
    - Organizes output files into a dedicated 'mhci' subdirectory within the specified output directory.
    - Handles process execution and error reporting.

Arguments:
    json_file (str or Path): Path to the input JSON file containing peptide data for MHCI prediction.
    tool_path (str or Path): Path to the MHCI prediction tool script to be executed.
    output_dir (str or Path): Directory where output files will be stored. Results are placed in an 'mhci' subfolder.

Requirements:
    - The specified MHCI tool script must be executable with Python 3 and accept '-j', '-o', and '-f' arguments.
    - Python packages: subprocess, pathlib.

Usage Example:
    run("input/peptides.json", "/tools/mhci_predictor.py", "results/epitopes")

Outputs:
    Prints a success message if the MHCI prediction completes successfully.
    Prints an error message and stderr output if the prediction fails.

Author: Nadia
"""
import subprocess
from pathlib import Path

def run(json_file, tool_path, output_dir):
    """
    Run MHCI prediction on the specified JSON file using the provided tool script.
    Args:
        json_file (str or Path): Path to the input JSON file containing peptide data.
        tool_path (str or Path): Path to the MHCI prediction tool script.
        output_dir (str or Path): Directory where output files will be stored.
    """
    output_dir = Path(output_dir) / "mhci"
    output_dir.mkdir(parents=True, exist_ok=True)

    out_base = Path(json_file).stem
    output_prefix = output_dir / out_base

    cmd = ["python3", tool_path, "-j", str(json_file), "-o", str(output_prefix), "-f", "json"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        print(f"✅ MHCI done: {json_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MHCI error: {json_file.name}")
        print(e.stderr)
