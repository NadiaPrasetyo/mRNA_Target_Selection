"""
run_dnds.py
Runner for MACSE and HyPhy analysis pipeline.
Overview:
    - Generates codon-aware alignments using the MACSE tool.
    - Builds a phylogenetic tree from the codon alignment using FastTree.
    - Optionally runs HyPhy analysis (FEL, FUBAR, SLAC) on the alignment and tree.
    - Cleans up temporary files after processing.
Arguments:
    tool_path (Path): Path to the directory containing external tools (e.g., /path/to/tools/).
    input_fasta (Path): Path to the input FASTA file containing nucleotide sequences.
    output_dir (Path): Directory where results will be saved.
    run_hyphy_analysis (bool): Whether to run HyPhy analysis after alignment and tree generation (default: True).
Requirements:
    - Python packages: subprocess, pathlib, logging, os, Bio.
    - A Conda environment named 'ext_tools_env' must be created and configured with the required dependencies.
    - The MACSE, FastTree, and HyPhy tools must be installed and accessible.
Outputs:
    <output_dir>/msa/<input_fasta_stem>_codon_aligned_raw.fasta  # Raw codon alignment from MACSE.
    <output_dir>/msa/<input_fasta_stem>_codon_aligned.fasta      # Cleaned codon alignment for HyPhy.
    <output_dir>/msa/<input_fasta_stem>.tree                    # Phylogenetic tree in Newick format.
    <output_dir>/<input_fasta_stem>_<method>_results.json       # HyPhy analysis results (if enabled).
Notes:
    - This script ensures that the MACSE alignment is cleaned and compatible with HyPhy.
    - Logs are generated to provide detailed information about the execution process.
    - Ensure that the Conda environment is activated and accessible before running the script.
Author: Nadia
"""

import subprocess
import os
import logging
from pathlib import Path
from tools import common
from Bio import SeqIO

def clean_macse_alignment(input_fasta: Path, output_fasta: Path):
    """
    Clean MACSE codon alignment for HyPhy:
    - Replace MACSE special characters with '-'
    - Mask stop codons (TAA, TAG, TGA) with gaps
    - Ensure length is divisible by 3
    """
    valid_codons = {"TAA", "TAG", "TGA"}
    cleaned_records = []

    for record in SeqIO.parse(str(input_fasta), "fasta"):
        seq = str(record.seq).upper()
        # Replace MACSE symbols
        seq = seq.replace("!", "-").replace("~", "-").replace("#", "-")
        # Codon-by-codon masking
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        new_codons = []
        for codon in codons:
            if len(codon) == 3:
                if codon in valid_codons:
                    new_codons.append("---")  # mask stop codon
                else:
                    new_codons.append(codon)
        new_seq = "".join(new_codons)
        # Ensure multiple of 3
        new_seq = new_seq[:len(new_seq) - (len(new_seq) % 3)]
        record.seq = type(record.seq)(new_seq)
        cleaned_records.append(record)

    SeqIO.write(cleaned_records, str(output_fasta), "fasta")


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
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
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
    Args:
        tool_path (unused): Path to the directory containing external tools kept for interface consistency
        input_fasta: Path to the input FASTA file containing nucleotide sequences
        output_dir: Directory where results will be saved
        run_hyphy_analysis: Whether to run HyPhy analysis after alignment and tree generation (default: True)
    """
    msa_path = output_dir / "msa"
    msa_path.mkdir(parents=True, exist_ok=True)
    raw_alignment = msa_path / f"{input_fasta.stem}_codon_aligned_raw.fasta"
    alignment_file = msa_path / f"{input_fasta.stem}_codon_aligned.fasta"

    # Run MACSE
    subprocess.run([
        "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
        "macse", "-prog", "alignSequences",
        "-seq", str(input_fasta),
        "-out_NT", str(raw_alignment)
    ], check=True)

    if not raw_alignment.exists():
        raise RuntimeError("❌ MACSE failed to generate codon alignment.")
    logging.info(f"✅ MACSE completed. Codon alignment: {raw_alignment}")

    # Clean MACSE alignment for HyPhy
    clean_macse_alignment(raw_alignment, alignment_file)
    
        # Step 2: Build a tree from the codon alignment
    logging.info("🌳 Building phylogenetic tree from codon alignment...")
    tree_file = msa_path / f"{input_fasta.stem}.tree"
    with open(tree_file, "w") as tree_out:
        subprocess.run([
            "conda", "run", "-n", common.EXT_TOOLS_ENV_NAME,
            "fasttree", "-nt", str(alignment_file)
        ], check=True, stdout=tree_out)

    # Validate tree content
    tree_str = tree_file.read_text().strip()
    if not tree_str or tree_str in {"()", ";", "();"}:
        logging.warning(f"⚠️ FastTree produced a degenerate tree for {input_fasta}. Skipping HyPhy.")
        return
    if not tree_str.startswith("("):
        raise RuntimeError(f"❌ Invalid Newick tree generated: {tree_file}")

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
