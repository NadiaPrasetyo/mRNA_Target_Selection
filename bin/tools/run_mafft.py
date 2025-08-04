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
    
def correct_tree_file_headers(tree_file: Path):
    """
    Corrects the headers in a tree file to match the sequence names in a FASTA file
    by removing the num_suffix and replacing the '_' with '.' and '|' in the tree file.
    Args:
    - tree_file: Path to the tree file.
    """
    with open(tree_file, 'r') as f:
        tree_content = f.read()

    for line in tree_content.splitlines():
        if line.startswith('('):
            # Split the line by commas to get individual sequence names
            seq_names = line.strip('()').split(',')
            corrected_names = []
            # 4_A0A0H3K6Z9_BA000018_3 -> A0A0H3K6Z9|BA000018.3
            for name in seq_names:
                # Remove the numeric suffix and replace '_' with '.'
                parts = name.split('_')
                if len(parts) > 1:
                    # Join all parts except the last one with '_'
                    corrected_name = '.'.join(parts[:-1]) + '|' + parts[-1]
                    corrected_names.append(corrected_name)
                    print(f"Corrected name: {corrected_name}")
                else:
                    # If no underscore, just append the name as is
                    corrected_names.append(name)
            # Join the corrected names back into a single string
            corrected_line = '(' + ','.join(corrected_names) + ')'
            tree_content = tree_content.replace(line, corrected_line)
    # Write the corrected content back to the tree file
    with open(tree_file, 'w') as f:
        f.write(tree_content)
        
    
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
        correct_tree_file_headers(final_output_tree_file)

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