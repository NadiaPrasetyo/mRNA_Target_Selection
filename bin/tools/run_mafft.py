"""
run_mafft.py
Runner for MAFFT multiple sequence alignment tool.
Overview:
    - Ensures the required Conda environment and dependencies are present.
    - Runs MAFFT (L-INS-i algorithm) on the provided input FASTA file using the external_tools_env Conda environment.
    - Outputs the aligned sequences in Clustal format, along with a guide tree.
    - Restores original headers in the alignment and tree files.
    - Optionally runs Rate4Site on the aligned sequences and guide tree if both are successfully generated.
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
    <output_dir>/<input_fasta_stem>_aligned.fasta   # Aligned sequences in Clustal format.
    <output_dir>/<input_fasta_stem>.tree            # Guide tree in Newick format (if generated).
    <output_dir>/rate4site_results/                 # Rate4Site results directory (if Rate4Site is run).
Author: Nadia
"""
import subprocess
from pathlib import Path
import logging
import shutil
from tools import common
import re

############################ HELPER FUNCTIONS ############################

def count_sequences(input_fasta: Path) -> int:
    """
    Count the number of sequences in a FASTA file.
    Args:
    - input_fasta: Path to the input FASTA file.
    Returns:
    - int: Number of sequences in the FASTA file.
    """
    with open(input_fasta, 'r') as f:
        return sum(1 for line in f if line.startswith('>'))
    

def rename_tree_headers_from_files(tree_file, fasta_file):
    """
    Rename tree headers in the Newick tree file to match the FASTA headers
    to ensure consistency between the alignment and the tree.
    Args:
    - tree_file (Path): Path to the Newick tree file.
    - fasta_file (Path): Path to the FASTA file containing sequence headers.
    Returns:
    - None: The tree file is modified in place.
    """
    # Read tree and FASTA file contents
    with open(tree_file, 'r') as tf:
        tree_str = tf.read()
    with open(fasta_file, 'r') as ff:
        fasta_str = ff.read()

    # Step 1: Extract FASTA headers
    fasta_headers = re.findall(r'^>(\S+)', fasta_str, re.MULTILINE)
    header_map = {}

    # Step 2: Build mapping: tree-style ID → FASTA header
    for header in fasta_headers:
        if '|' not in header:
            continue
        uniprot, accver = header.split('|', 1)
        if '.' not in accver:
            continue
        accession, version = accver.split('.', 1)
        tree_id = f"{uniprot}_{accession}_{version}"
        header_map[tree_id] = f"{uniprot}|{accession}.{version}"

    # Step 3: Replace all matching IDs in tree
    def replacer(match):
        """Replace matched tree ID with corresponding FASTA header."""
        original = match.group(0)
        parts = original.split('_', 1)
        if len(parts) != 2:
            return original
        id_body = parts[1]
        return header_map.get(id_body, original)

    # Match tree IDs like 1_A0A0H3K6Z9_CP002114_3
    updated_tree = re.sub(r'\b\d+_[A-Z0-9]+\_[A-Z0-9]+\_\d+\b', replacer, tree_str)

    # Step 4: Overwrite the original tree file
    with open(tree_file, 'w') as tf:
        tf.write(updated_tree)
    
############################ RUN MAFFT FUNCTION ############################


def run(tool_path: Path, input_fasta: Path, output_dir: Path, rate4site: bool = True):
    """
    Runs MAFFT using the external_tools_env conda environment.
    Args:
    - tool_path: Path to the MAFFT executable (ignored, required by interface).
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    - rate4site: Boolean flag to indicate if Rate4Site analysis should be run.
    Raises:
    - RuntimeError: If Conda is not available or if the MAFFT command fails.
    Outputs:
    - <output_dir>/<input_fasta_stem>_aligned.fasta: Aligned sequences in Clustal format.
    - <output_dir>/<input_fasta_stem>.tree: Guide tree in Newick format (if generated).
    """
    # Check if the input FASTA file contains more than one sequence
    sequence_count = count_sequences(input_fasta)
    if sequence_count <= 1:
        logging.info(f"ℹ️ Input FASTA file {input_fasta.name} contains only one sequence. Alignment is not required. Exiting gracefully.")
        return


    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    common.create_conda_env_if_needed()
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / (input_fasta.stem + "_aligned.fasta")

    command = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "mafft", "--localpair", "--maxiterate", "1000", #L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)
        "--reorder", # Output order: aligned.
        "--treeout", # Guide tree is output to the input.tree file
        "--amino", # Assume the sequences are amino acid
        str(input_fasta)  # Input FASTA file
    ]

    logging.info(f"🔍 Running MAFFT on {input_fasta}...")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Save the aligned output
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        logging.info(f"✅ MAFFT alignment completed. Output saved to {output_file}")

        # Check if a tree file was generated
        temp_file = input_fasta.with_suffix(input_fasta.suffix + ".tree")  # Input file with .tree suffix
        if not temp_file.exists():
            logging.warning("⚠️ MAFFT did not generate a tree file.")
        
        #move the tree file to the output directory if it exists
        final_output_tree_file = output_dir / f"{input_fasta.stem}.tree"
        shutil.move(temp_file, final_output_tree_file)
        rename_tree_headers_from_files(final_output_tree_file, output_file)

        logging.info("🔁 Restored original headers in alignment and tree.")
        logging.info("✅ MAFFT alignment and tree restoration completed. Output files: "
                    f"{output_file}, {final_output_tree_file}")
        
        cleanup_files = [temp_file]
        for f in cleanup_files:
            if f.exists():
                f.unlink()

    except subprocess.CalledProcessError as e:
        logging.error("❌ Error running MAFFT:")
        logging.error(e.stderr)
        raise

    # If Rate4Site is enabled, run it on the aligned sequences and tree
    if rate4site:
        from tools import run_rate4site
        # call run_rate4site to run rate4site on the aligned sequences if a the ouput files were created
        if output_file.exists():
            rate4site_output_dir = output_dir / "rate4site_results"
            tree_file = output_dir / f"{input_fasta.stem}.tree"
            if tree_file.exists():
                rate4site_output_dir.mkdir(parents=True, exist_ok=True)
                # Run Rate4Site with the aligned output and the tree file
                run_rate4site.run(input_fasta=output_file, input_tree=tree_file, output_dir=rate4site_output_dir)
            else:
                logging.error("❌ MAFFT did not produce a tree file. Rate4Site will not be run.")
                raise RuntimeError("MAFFT alignment completed but no tree file was generated.")
        else:
            logging.error("❌ MAFFT did not produce the expected output file. Rate4Site will not be run.")
            raise RuntimeError("MAFFT alignment failed, output file not created.")