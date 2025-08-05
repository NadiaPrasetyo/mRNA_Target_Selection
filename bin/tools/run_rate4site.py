"""
run_rate4site.py

Runner for the Rate4Site evolutionary rate estimation tool.

Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Executes Rate4Site to estimate evolutionary rates for each site in a protein or DNA sequence.
    - Outputs the results and the processed tree file to the specified directory.

Arguments:
    input_fasta (Path): Path to the input FASTA file containing the sequence data.
    input_tree (Path): Path to the input tree file in Newick format.
    output_dir (Path): Directory where output files will be saved.

Requirements:
    - A Conda environment with Rate4Site installed.
    - Conda available in PATH.

Outputs:
    <output_dir>/<input_fasta_stem>.out   # File containing evolutionary rate estimates.
    <output_dir>/<input_fasta_stem>.tree  # Processed tree file.
Author: Nadia
"""
import argparse
import subprocess
import logging
import shutil
from pathlib import Path
from tools import common

def run(input_fasta: Path, input_tree: Path, output_dir: Path):
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    # Ensure conda env is ready (assumes this function is defined elsewhere)
    common.create_conda_env_if_needed()

    input_stem = input_fasta.stem

    # Prepare output paths
    output_file = output_dir / f"{input_stem}.out"
    tree_out_file = output_dir / f"{input_stem}.tree"

    # Make sure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Command construction
    cmd = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "rate4site",      # Assuming rate4site is in the conda environment
        "-s", str(input_fasta),
        "-t", str(input_tree),
        "-o", str(output_file),
        "-x", str(tree_out_file),
        "-im",                  # Rate inference method flag:
                                # -im = rates are inferred using the maximum likelihood method
                                # -ib = rates are inferred using the empirical Bayes method 
        "-Mw"                   # The following amino-acids models are supported:
                                # DAY (-md), JTT (-mj), REV (-mr), aaJC (-ma), LG (-Ml), WAG (-Mw)
    ]

    logging.info(f"🚀 Running Rate4Site for {input_fasta.name} with tree {input_tree.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info("✅ Rate4Site completed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Rate4Site failed with error: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(
        description="Run Rate4Site with specified input files and output directory.",
        epilog="Rate4Site is a tool for estimating evolutionary rates at each site of a protein or DNA sequence."
    )
    parser.add_argument("-f", "--input_fasta", type=Path, required=True, help="Path to the input FASTA file.")
    parser.add_argument("-t", "--input_tree", type=Path, required=True, help="Path to the input tree file.")
    parser.add_argument("-o", "--output_dir", type=Path, required=True, help="Path to the output directory.")

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    try:
        run(args.input_fasta, args.input_tree, args.output_dir)
    except Exception as e:
        logging.error(f"❌ An error occurred: {e}")
        exit(1)

if __name__ == "__main__":
    main()
