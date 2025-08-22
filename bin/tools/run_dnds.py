import subprocess
import os
import logging
import tempfile
from pathlib import Path
from tools import common


def run_hyphy(method, alignment, tree, output_dir):
    """
    Run a HyPhy method with the given alignment and tree files.
    Args:
        Method: The HyPhy method to run (e.g., FEL, FUBAR, SLAC)
        Alignment: Path to the alignment file
        Tree: Path to the tree file
        Output_dir: Directory to save the output
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
    Runs TranslatorX to generate codon-aware alignments (temporary),
    builds a tree, then optionally runs HyPhy analysis.
    """
    common.create_conda_env_if_needed()
    msa_path = output_dir / "msa"
    aln_prefix = msa_path / input_fasta.stem

    # Step 1: Run TranslatorX
    msa_path.mkdir(parents=True, exist_ok=True)  # Ensure msa_path exists
    logging.info("🔍 Running TranslatorX for codon-aware alignment...")
    subprocess.run([
        "translatorx",
        "-i", str(input_fasta),
        "-o", str(aln_prefix),
        "-p", "F"  # use MAFFT for protein alignment
    ], check=True)

    alignment_file = aln_prefix.with_suffix(".nt_ali.fasta")

    if not alignment_file.exists():
        raise RuntimeError("❌ TranslatorX failed to generate codon alignment.")
    logging.info(f"✅ TranslatorX completed. Codon alignment: {alignment_file}")
    
    # Step 2: Build a tree from the codon alignment
    logging.info("🌳 Building phylogenetic tree from codon alignment...")
    tree_file = aln_prefix.with_suffix(".tree")

    with open(tree_file, "w") as tree_out:
        subprocess.run(["fasttree", "-nt", str(alignment_file)], check=True, stdout=tree_out)
    logging.info(f"✅ Tree built: {tree_file}")

    # Step 3: Run HyPhy analysis if enabled
    if run_hyphy_analysis:
        logging.info("🔍 Running HyPhy analysis...")
        for method in ["FEL", "FUBAR", "SLAC"]:
            run_hyphy(method, alignment_file, tree_file, output_dir)
        logging.info("✅ HyPhy analysis completed.")

    # Clean up temporary files
    logging.info("🧹 Cleaning up temporary files...")
    for tmp_file in [alignment_file, tree_file]:
        if tmp_file.exists():
            tmp_file.unlink()
            logging.info(f"✅ Deleted temporary file: {tmp_file}")