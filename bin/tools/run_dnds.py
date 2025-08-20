import subprocess
import os
from pathlib import Path
import logging
from tools import run_mafft

def run_hyphy(method, alignment, tree, output_dir):
    """
    Run a HyPhy method with the given alignment and tree files.
    
    :param method: The HyPhy method to run (e.g., FEL, FUBAR, SLAC)
    :param alignment: Path to the alignment file
    :param tree: Path to the tree file
    :param output_dir: Directory to save the output
    """
    output_file = os.path.join(output_dir, f"{method}_results.json")
    command = [
        "hyphy", method.lower(),
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--output", output_file
    ]
    try:
        print(f"Running {method}...")
        subprocess.run(command, check=True)
        print(f"{method} completed. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {method}: {e}")


def run(tool_path: Path, input_fasta: Path, output_dir: Path, run_hyphy_analysis: bool = True):
    """
    Runs MAFFT to generate alignments and tree, then optionally runs HyPhy analysis.
    Args:
    - tool_path: Path to the MAFFT executable (ignored, required by interface).
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    - run_hyphy_analysis: Boolean flag to indicate if HyPhy analysis should be run.
    
    Raises:
    - RuntimeError: If MAFFT or HyPhy commands fail.
    """
    # Step 1: Run MAFFT to generate alignment and tree
    logging.info("🔍 Running MAFFT to generate alignment and tree...")
    run_mafft.run(tool_path, input_fasta, output_dir/"mafft", rate4site=False)

    # Paths to the generated alignment and tree files
    alignment_file = output_dir / "mafft" / f"{input_fasta.stem}_aligned.fasta"
    tree_file = output_dir / "mafft" / f"{input_fasta.stem}.tree"

    if not alignment_file.exists() or not tree_file.exists():
        logging.error("❌ MAFFT did not produce the required alignment or tree files.")
        raise RuntimeError("MAFFT failed to generate alignment or tree.")

    logging.info(f"✅ MAFFT completed. Alignment: {alignment_file}, Tree: {tree_file}")

    # Step 2: Run HyPhy analysis if enabled
    if run_hyphy_analysis:
        logging.info("🔍 Running HyPhy analysis...")
        for method in ["FEL", "FUBAR", "SLAC"]:
            run_hyphy(method, alignment_file, tree_file, output_dir)
        logging.info("✅ HyPhy analysis completed.")
