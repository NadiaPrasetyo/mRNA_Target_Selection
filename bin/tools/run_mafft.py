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


def restore_clustal_headers(clustal_file: Path, mapping: dict, output_file: Path):
    """
    Replace seq1, seq2, ... headers in a MAFFT CLUSTAL alignment file with original headers.
    Preserves alignment format.
    """
    import re

    # Compile regex to match alignment lines
    header_regex = re.compile(r"^(seq\d+)(\s+)")

    with open(clustal_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            match = header_regex.match(line)
            if match:
                seq_id, spacing = match.groups()
                original_id = mapping.get(seq_id)
                if not original_id:
                    raise ValueError(f"Missing mapping for {seq_id}")

                # Truncate/pad the original ID to match the width of seq_id
                new_id = original_id[:len(seq_id)].ljust(len(seq_id))
                line = line.replace(seq_id, new_id, 1)

            outfile.write(line)

    logging.info(f"🪄 Restored headers in CLUSTAL alignment saved to {output_file}")

def restore_tree_names(tree_file: Path, mapping: dict):
    from Bio import Phylo
    import re

    tree = Phylo.read(tree_file, "newick")

    for terminal in tree.get_terminals():
        # Extract seqN from "1_seq1", "2_seq2", etc.
        match = re.match(r"\d+_(seq\d+)", terminal.name)
        if match:
            simplified_name = match.group(1)
        else:
            simplified_name = terminal.name

        original = mapping.get(simplified_name)
        if original:
            terminal.name = original
        else:
            logging.warning(f"⚠️ No mapping found for terminal: {terminal.name}")

    Phylo.write(tree, tree_file, "newick")
    logging.info(f"🌳 Tree tip names restored in {tree_file}")


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
    # Rename headers
    renamed_temp_fasta = output_dir / f"{input_fasta.stem}_renamed.fasta"
    mapping = rename_fasta_headers(input_fasta, renamed_temp_fasta)

    renamed_temp_output_file = output_dir / f"{input_fasta.stem}_aligned_simplified.fasta"

    command = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "mafft", "--localpair", "--maxiterate", "1000", #L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)
        "--clustalout", # Output format: clustal format
        "--reorder", # Output order: aligned.
        "--treeout", # Guide tree is output to the input.tree file
        "--amino", # Assume the sequences are amino acid
        str(renamed_temp_fasta)
    ]

    logging.info(f"🔍 Running MAFFT on {renamed_temp_fasta}...")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        # Save the aligned output
        with open(renamed_temp_output_file, 'w') as f:
            f.write(result.stdout)

        logging.info(f"✅ MAFFT alignment completed. Output saved to {renamed_temp_output_file}")

        # Move the generated tree file to the output directory
        tree_file = renamed_temp_fasta.with_suffix(input_fasta.suffix + ".tree")  # Input file with .tree suffix
        if tree_file.exists():
            temp_simplified_tree_output_file = output_dir / f"{input_fasta.stem}_simplified.tree"
            shutil.move(tree_file, temp_simplified_tree_output_file)
            logging.info(f"✅ Guide tree saved to {temp_simplified_tree_output_file}")
        else:
            logging.error("❌ MAFFT did not generate a tree file.")
            raise RuntimeError("MAFFT failed to produce guide tree.")
        
        # Restore headers
        final_output_restored_fasta = output_dir / f"{input_fasta.stem}_aligned.fasta"
        final_output_tree_file = output_dir / f"{input_fasta.stem}.tree"

        restore_clustal_headers(renamed_temp_output_file, mapping, final_output_restored_fasta)
        shutil.move(temp_simplified_tree_output_file, final_output_tree_file)
        restore_tree_names(final_output_tree_file, mapping)
        
        logging.info("🔁 Restored original headers in alignment and tree.")
        logging.info("✅ MAFFT alignment and tree restoration completed. Output files: "
                    f"{final_output_restored_fasta}, {final_output_tree_file}")
        
        cleanup_files = [renamed_temp_fasta, renamed_temp_output_file, temp_simplified_tree_output_file]

    except subprocess.CalledProcessError as e:
        logging.error("❌ Error running MAFFT:")
        logging.error(e.stderr)
        raise

    from tools import run_rate4site
    # call run_rate4site to run rate4site on the aligned sequences if a the ouput files were created
    if final_output_restored_fasta.exists():
        rate4site_output_dir = output_dir / "rate4site_results"
        tree_file = output_dir / f"{input_fasta.stem}.tree"
        if tree_file.exists():
            rate4site_output_dir.mkdir(parents=True, exist_ok=True)
            # Run Rate4Site with the aligned output and the tree file
            run_rate4site.run(input_fasta=final_output_restored_fasta, input_tree=tree_file, output_dir=rate4site_output_dir)
        else:
            logging.error("❌ MAFFT did not produce a tree file. Rate4Site will not be run.")
            raise RuntimeError("MAFFT alignment completed but no tree file was generated.")
    else:
        logging.error("❌ MAFFT did not produce the expected output file. Rate4Site will not be run.")
        raise RuntimeError("MAFFT alignment failed, output file not created.")