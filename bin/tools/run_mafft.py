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

############################ HELPER FUNCTIONS ############################
def rename_fasta_headers(original_fasta: Path, renamed_fasta: Path) -> dict:
    """
    Replace FASTA headers with generic seq1, seq2, ... and save the mapping.
    Returns a dictionary mapping new -> original headers.
    """
    from Bio import SeqIO

    records = list(SeqIO.parse(original_fasta, "fasta"))
    mapping = {}
    renamed_records = []

    for i, record in enumerate(records, 1):
        new_id = f"seq{i}"
        mapping[new_id] = record.id
        record.id = new_id
        record.name = ""
        record.description = ""
        renamed_records.append(record)

    SeqIO.write(renamed_records, renamed_fasta, "fasta")
    return mapping


def restore_fasta_headers(renamed_fasta: Path, mapping: dict, restored_fasta: Path):
    """
    Replace generic headers (seq1, ...) back with original headers using the mapping.
    """
    from Bio import SeqIO

    records = []
    for record in SeqIO.parse(renamed_fasta, "clustal"):
        original_id = mapping.get(record.id)
        if not original_id:
            raise ValueError(f"Missing mapping for {record.id}")
        record.id = original_id
        record.name = ""
        record.description = ""
        records.append(record)

    SeqIO.write(records, restored_fasta, "fasta")

def restore_tree_names(tree_file: Path, mapping: dict):
    from Bio import Phylo

    tree = Phylo.read(tree_file, "newick")

    for terminal in tree.get_terminals():
        if terminal.name in mapping:
            terminal.name = mapping[terminal.name]

    Phylo.write(tree, tree_file, "newick")

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
    
############################ RUN MAFFT FUNCTION ############################

def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int):
    from tools import run_rate4site

    sequence_count = count_sequences(input_fasta)
    if sequence_count <= 1:
        logging.info(f"ℹ️ Input FASTA {input_fasta.name} has one sequence. Skipping alignment.")
        return

    if not shutil.which("conda"):
        logging.error("❌ Conda is not available.")
        raise RuntimeError("Conda is required.")

    common.create_conda_env_if_needed()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Rename headers
    renamed_fasta = output_dir / f"{input_fasta.stem}_renamed.fasta"
    mapping = rename_fasta_headers(input_fasta, renamed_fasta)

    mafft_output_file = output_dir / f"{input_fasta.stem}_aligned_simplified.fasta"

    # Run MAFFT
    command = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "mafft", "--localpair", "--maxiterate", "1000",
        "--clustalout", "--reorder", "--treeout", "--amino",
        str(renamed_fasta)
    ]

    logging.info("🔍 Running MAFFT...")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        with open(mafft_output_file, 'w') as f:
            f.write(result.stdout)
        logging.info(f"✅ MAFFT output saved to {mafft_output_file}")
    except subprocess.CalledProcessError as e:
        logging.error("❌ MAFFT failed:")
        logging.error(e.stderr)
        raise

    mafft_tree_file = renamed_fasta.with_suffix(".fasta.tree")

    if not mafft_tree_file.exists():
        logging.error("❌ MAFFT did not generate a tree file.")
        raise RuntimeError("MAFFT failed to produce guide tree.")
    
    # Restore headers
    restored_fasta = output_dir / f"{input_fasta.stem}_aligned.fasta"
    tree_output_file = output_dir / f"{input_fasta.stem}.tree"
    restore_fasta_headers(mafft_output_file, mapping, restored_fasta)
    shutil.move(mafft_tree_file, tree_output_file) # Rename tree file to match output .tree file
    restore_tree_names(tree_output_file, mapping)
    logging.info("🔁 Restored original headers in alignment and tree.")
    logging.info("✅ MAFFT alignment and tree restoration completed. Output files: "
                 f"{restored_fasta}, {tree_output_file}")
    
    cleanup_files = [renamed_fasta, mafft_output_file, mafft_tree_file]
    for file in cleanup_files:
        if file.exists():
            file.unlink()
            logging.info(f"🗑️ Cleaned up temporary file: {file}")

    ##################### RUN RATE4SITE ############################
    if restored_fasta.exists() and tree_output_file.exists():
        logging.info("🔍 Running Rate4Site on the MAFFT output...")
    rate4site_output_dir = output_dir / "rate4site_results"
    rate4site_output_dir.mkdir(parents=True, exist_ok=True)
    run_rate4site.run(input_fasta=restored_fasta, input_tree=tree_output_file, output_dir=rate4site_output_dir)
