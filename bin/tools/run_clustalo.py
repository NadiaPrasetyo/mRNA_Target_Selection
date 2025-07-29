from pathlib import Path
import subprocess


def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int):
    """
    Run Clustal Omega on the given input FASTA file.

    Parameters:
        tool_path (Path): Path to the Clustal Omega executable.
        input_fasta (Path): Path to the input FASTA file.
        output_dir (Path): Directory to store output files.
        batch_size (int): Unused, kept for interface consistency.
    """
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output file paths
    msa_output = output_dir / f"{input_fasta.stem}.aln"
    guide_tree_output = output_dir / f"{input_fasta.stem}.dnd"

    # Build the command
    command = [
        "chmod", "u+x", str(tool_path),  # Ensure the tool is executable
        str(tool_path),
        "-i", str(input_fasta),
        "-o", str(msa_output),
        "--outfmt=clu",
        "--guidetree-out", str(guide_tree_output),
        "--force"
    ]

    # Run Clustal Omega
    try:
        subprocess.run(command, check=True)
        print(f"Clustal Omega alignment completed. Output saved to {msa_output} and guide tree to {guide_tree_output}.")
    except subprocess.CalledProcessError as e:
        print("Error running Clustal Omega:")
        print(e.stderr.decode())
        raise