import subprocess
import os
import logging
from pathlib import Path
from tools import common


def run_hyphy(method, alignment, tree, output_dir):
    """
    Run a HyPhy method with the given alignment and tree files.
    Args:
        method: The HyPhy method to run (e.g., FEL, FUBAR, SLAC)
        alignment: Path to the alignment file
        tree: Path to the tree file
        output_dir: Directory to save the output
    """
    input_stem = alignment.stem
    output_file = os.path.join(output_dir, f"{input_stem}_{method}_results.json")
    command = [
        "hyphy", method.lower(),
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--output", output_file
    ]
    try:
        logging.info(f"▶️ Running {method}...")
        subprocess.run(command, check=True)
        logging.info(f"✅ {method} completed. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error running {method}: {e}")


def run(tool_path: Path, input_fasta: Path, output_dir: Path, run_hyphy_analysis: bool = True):
    """
    Runs MACSE to generate codon-aware alignments,
    builds a tree, then optionally runs HyPhy analysis.
    """
    common.create_conda_env_if_needed()
    msa_path = output_dir / "msa"
    msa_path.mkdir(parents=True, exist_ok=True)

    alignment_file = msa_path / f"{input_fasta.stem}_codon_aligned.fasta"

    # Step 1: Run MACSE for codon-aware alignment
    logging.info("🔍 Running MACSE for codon-aware alignment...")
    subprocess.run([
        "macse", "-prog", "alignSequences",
        "-seq", str(input_fasta),
        "-out_aln", str(alignment_file)
    ], check=True)

    if not alignment_file.exists():
        raise RuntimeError("❌ MACSE failed to generate codon alignment.")
    logging.info(f"✅ MACSE completed. Codon alignment: {alignment_file}")

    # Step 2: Build a tree from the codon alignment
    logging.info("🌳 Building phylogenetic tree from codon alignment...")
    tree_file = msa_path / f"{input_fasta.stem}.tree"
    with open(tree_file, "w") as tree_out:
        subprocess.run(["fasttree", "-nt", str(alignment_file)], check=True, stdout=tree_out)
    logging.info(f"✅ Tree built: {tree_file}")

    # Step 3: Run HyPhy analysis if enabled
    if run_hyphy_analysis:
        logging.info("🔍 Running HyPhy analysis...")
        for method in ["FEL", "FUBAR", "SLAC"]:
            run_hyphy(method, alignment_file, tree_file, output_dir)
        logging.info("✅ HyPhy analysis completed.")

    # Step 4: Clean up temporary files
    logging.info("🧹 Cleaning up temporary files...")
    for temp_file in [alignment_file, tree_file]:
        if temp_file.exists():
            temp_file.unlink()
    logging.info("✅ Cleanup completed.")
