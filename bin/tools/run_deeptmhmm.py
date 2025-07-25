"""
run_deeptmhmm.py
Run DeepTMHMM to predict membrane topology of proteins.
This script uses the BioLib library to run DeepTMHMM predictions on a given FASTA file.

Overview:
    - Ensures the required BioLib library is available.
    - Runs DeepTMHMM on the provided input FASTA file.
    - Outputs membrane topology predictions as files in the specified output directory.
Arguments:
    tool_path (Path): Directory containing tools (not used in this script).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
Requirements:
    - BioLib library installed in the Python environment.
Outputs:
    <output_dir>/<input_fasta_stem>_deeptmhmm_results/   # Membrane topology prediction results
"""
import logging
from pathlib import Path
import biolib


def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 0):
    """
    Run DeepTMHMM via BioLib to predict membrane topology.

    Args:
        tool_path (Path): Not used for DeepTMHMM via BioLib but required for signature consistency.
        input_fasta (Path): Path to the input FASTA file.
        output_dir (Path): Directory to save the result files.
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"🚀 Running DeepTMHMM on {input_fasta.name}")

        # Load DeepTMHMM tool from BioLib
        deeptmhmm = biolib.load('DTU/DeepTMHMM')

        # Run the prediction job
        deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {input_fasta.as_posix()}')

        # Collect files before saving
        existing_files = set(output_dir.glob("*"))

        # Save result files to the output directory
        deeptmhmm_job.save_files(output_dir.as_posix())

        # rename only the files that was saved in the output files to match the input file stem
        input_stem = input_fasta.stem
        # Collect files after saving
        new_files = set(output_dir.glob("*")) - existing_files
        for file in new_files:
            if file.is_file():
                new_name = output_dir / f"{input_stem}_{file.name}"
                file.rename(new_name)

        # Log the saved files
        saved_files = list(output_dir.glob(f"{input_stem}_*"))
        if saved_files:
            logging.info(f"Files saved: {', '.join(str(f.name) for f in saved_files)}")
        else:
            logging.warning(f"No files were saved in {output_dir} for {input_fasta.name}")

        logging.info(f"✅ DeepTMHMM results saved to {output_dir}")
    except Exception as e:
        logging.error(f"❌ DeepTMHMM failed on {input_fasta.name}: {e}")
        raise
