"""
Runner for MAFFT multiple sequence alignment tool.
Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Runs MAFFT (L-INS-i algorithm) on the provided input FASTA file using the external_tools_env Conda environment.
    - Outputs the aligned sequences in Clustal format, along with a guide tree.
Arguments:
    tool_path (Path): Path to the MAFFT executable (ignored, required by interface).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
    batch_size (int, optional): Unused, present for interface compatibility.
Requirements:
    - ext_tools_dependencies.yml (defines Conda environment).
    - MAFFT installed in the specified Conda environment.
    - Conda available in PATH.
Outputs:
    <output_dir>/<input_fasta_stem>_aligned.fasta   # Aligned sequences in Clustal format
Author: Nadia
run_mafft.py

"""
import subprocess
from pathlib import Path
import logging
import shutil
import common

def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int):
    """
    Runs MAFFT using the external_tools_env conda environment.
    
    Parameters:
    - tool_path: Path to the MAFFT executable (ignored, required by interface).
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    - batch_size: Batch size for processing (not used here, kept for compatibility).
    """
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed()

    # Output file path
    output_file = output_dir / f"{input_fasta.stem}_aligned.fasta"

    # Construct the MAFFT command with conda environment
    command = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "mafft", "--localpair", "--maxiterate", "1000", #L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)
        "--clustalout", # Output format: clustal format
        "--reorder", # Output order: aligned.
        "--treeout", # Guide tree is output to the input.tree file
        "--amino", # Assume the sequences are amino acid
        str(input_fasta), ">", str(output_file)
    ]

    try:
        # Run the command
        subprocess.run(" ".join(command), shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"✅ MAFFT alignment completed. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error("❌ Error running MAFFT:")
        logging.error(e.stderr.decode())
        raise
