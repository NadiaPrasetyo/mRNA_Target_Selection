"""

Command-line tool to run DeepLocPro for protein subcellular localization prediction using the BioLib CLI.

Overview:
    - Executes DeepLocPro on a given FASTA input file.
    - Organizes output files into structured directories for results and plots.
    - Handles errors and logs output for debugging.

Arguments:
    tool_path (str): Working directory to run DeepLocPro in.
    input_file (str): Path to the FASTA input file.
    output_dir (str): Directory to store output files.
    group (str): Sample group name (e.g., 'positive').

Requirements:
    - DeepLocPro installed and available in PATH as 'deeplocpro'.
    - Input FASTA file containing protein sequences.

Outputs:
    <output_dir>/results*         # Result files from DeepLocPro.
    <output_dir>/plots/*          # Plot images and non-results files.

Author: Nadia
run_deeplocpro.py

"""
import subprocess
import logging
from pathlib import Path
import shutil

def run_deeplocpro(tool_path: str, input_file: str, output_dir: str, group: str):
    """
    Run DeepLocPro using the BioLib CLI and store results in a structured directory.

    Args:
        tool_path (str): Working directory to run DeepLocPro in
        input_file (str): Path to FASTA input file
        output_dir (str): Where to store output files
        group (str): Sample group name (e.g., 'positive')
    """
    working_dir = Path(tool_path).resolve()
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    plots_dir = output_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Running DeepLocPro in {working_dir} for {input_file.name} with group: {group}")

    try:
        cmd = [
            "deeplocpro",
            "-f", str(input_file),
            "-o", "output",
            "-p",
            "-d", "cpu",
            "-g", group
        ]

        result = subprocess.run(cmd, cwd=working_dir, capture_output=True, text=True)

        logging.debug(f"STDOUT:\n{result.stdout}")
        logging.debug(f"STDERR:\n{result.stderr}")

        # Check output directory inside working directory
        bio_output_dir = working_dir / "output"
        if not bio_output_dir.exists():
            raise FileNotFoundError("❌ No output directory created by DeepLocPro.")

        for f in bio_output_dir.iterdir():
            # Send plots and non-results-prefixed files to plots_dir
            if f.suffix.lower() in {".png", ".jpg", ".svg"} or not f.name.startswith("results"):
                shutil.move(str(f), plots_dir / f.name)
            else:
                shutil.move(str(f), output_dir / f.name)

        print(f"✅ DeepLocPro completed for {input_file.name}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
