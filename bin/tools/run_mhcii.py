"""
run_mhcii.py
Utility to execute MHCII epitope prediction tool on a given JSON input file.

Overview:
    - Runs the specified MHCII prediction tool as a subprocess.
    - Organizes output files in a dedicated 'mhcii' subdirectory within the given output directory.
    - Handles subprocess errors and provides user-friendly status messages.

Arguments:
    json_file (str or Path): Path to the input JSON file containing peptide data for prediction.
    tool_path (str or Path): Path to the MHCII prediction tool script to be executed.
    output_dir (str or Path): Directory where output files will be stored (creates 'mhcii' subdirectory).

Requirements:
    - The specified tool must be executable with Python 3 and accept '-j', '-o', and '-f' arguments.
    - Python packages: subprocess, pathlib.

Usage Example:
    run("input/peptides.json", "/tools/mhcii_predictor.py", "results/epitopes")

Outputs:
    Prints a success message if the tool runs successfully, or an error message with details if it fails.

Author: Nadia
"""
import subprocess
from pathlib import Path

def run(json_file, tool_path, output_dir):
    """
    Run the MHCII prediction tool on the provided JSON file.
    Args:
        json_file (str or Path): Path to the input JSON file.
        tool_path (str or Path): Path to the MHCII prediction tool script.
        output_dir (str or Path): Directory where output files will be stored.
    """
    output_dir = Path(output_dir) / "mhcii"
    output_dir.mkdir(parents=True, exist_ok=True)

    out_base = Path(json_file).stem
    output_prefix = output_dir / out_base

    cmd = ["python3", tool_path, "-j", str(json_file), "-o", str(output_prefix), "-f", "json"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        print(f"✅ MHCII done: {json_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MHCII error: {json_file.name}")
        print(e.stderr)
