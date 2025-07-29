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
from tools import common

def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int):
    """
    Runs MAFFT using the external_tools_env conda environment.
    Args:
    - tool_path: Path to the MAFFT executable (ignored, required by interface).
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    - batch_size: Batch size for processing (not used here, kept for compatibility).
    Raises:
    - RuntimeError: If Conda is not available or if the MAFFT command fails.
    Outputs:
    - <output_dir>/<input_fasta_stem>_aligned.fasta: Aligned sequences in Clustal format.
    - <output_dir>/<input_fasta_stem>.tree: Guide tree in Newick format (if generated).
    """
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed()

    output_file = output_dir / f"{input_fasta.stem}_aligned.fasta"

    command = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "mafft", "--localpair", "--maxiterate", "1000", #L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)
        "--clustalout", # Output format: clustal format
        "--reorder", # Output order: aligned.
        "--treeout", # Guide tree is output to the input.tree file
        "--amino", # Assume the sequences are amino acid
        str(input_fasta)
    ]

    logging.info(f"🔍 Running MAFFT on {input_fasta}...")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Save the aligned output
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        # Log stderr if any
        if result.stderr:
            logging.debug(f"MAFFT stderr:\n{result.stderr}")

        # Move the tree file
        tree_file = input_fasta.with_suffix(".tree")
        if tree_file.exists():
            shutil.move(tree_file, output_dir / tree_file.name)
            tree_output = output_dir / tree_file.name
        else:
            tree_output = "not generated"

        logging.info(f"✅ MAFFT alignment completed. Output saved to {output_file}, tree file saved to {tree_output}")

    except subprocess.CalledProcessError as e:
        logging.error("❌ Error running MAFFT:")
        logging.error(e.stderr)
        raise
