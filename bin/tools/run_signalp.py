"""
run_signalp.py
Command-line tool to run SignalP predictions on protein FASTA files.

Overview:
    - Executes the SignalP executable on a given input FASTA file.
    - Supports batch processing of sequences for large input files.
    - Collects and organizes SignalP output files, including prediction results and artifacts.
    - Handles temporary directories for intermediate files and moves generated artifacts for easy access.

Arguments:
    input_fasta (Path): Path to the input protein FASTA file.
    output_dir (Path): Directory where SignalP results and artifacts will be stored.
    --batch-size (int, optional): Number of sequences to process per batch (default: 10000).
    --signalp-path (Path, required): Path to the SignalP executable.

Requirements:
    - SignalP executable available at the specified path.
    - Python packages: argparse, pathlib, subprocess, shutil, os, sys.

Usage Example:
    python run_signalp.py proteins.fasta results/ --signalp-path /opt/signalp/signalp --batch-size 5000

Outputs:
    Writes SignalP prediction results to a text file in the output directory.
    Moves generated artifact files (plots, predictions, mature sequences) to a temporary subdirectory.
    Prints informative messages about progress and file locations.

Author: Nadia
"""
import subprocess
from pathlib import Path
import sys
import shutil
from os.path import abspath

def run(signalp_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 10000, organism: str = "gram+"):
    """
    Run SignalP on the given input FASTA file and save results to the specified output directory.
    Args:
        signalp_path (Path): Path to the SignalP executable.
        input_fasta (Path): Path to the input protein FASTA file.
        output_dir (Path): Directory where results will be saved.
        batch_size (int): Number of sequences to process per batch (default: 10000).
    """
    input_fasta_abs = Path(abspath(str(input_fasta)))
    output_dir_abs = Path(abspath(str(output_dir)))
    tmp_dir_abs = output_dir_abs / "tmp"
    tmp_dir_abs.mkdir(exist_ok=True)

    if not signalp_path.exists():
        print(f"❌ SignalP executable not found at: {signalp_path}")
        sys.exit(1)

    if not input_fasta.exists():
        print(f"❌ Input FASTA file does not exist: {input_fasta}")
        sys.exit(1)

    print(f"[INFO] Running SignalP on: {input_fasta}")
    print(f"[INFO] Output directory: {output_dir}")

    basename = input_fasta.stem
    output_file = output_dir_abs / f"{basename}_signalp_phred.txt"

    cmd = [
        "./signalp",
        "-fasta", str(input_fasta_abs),
        "-format", "long",
        "-mature",
        "-batch", str(batch_size),
        "-stdout",
        "-tmp", str(tmp_dir_abs),
        "-org", organism
    ]

    with output_file.open("w") as outfile:
        try:
            subprocess.run(cmd, cwd=signalp_path.parent, stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ SignalP failed: {e.stderr.decode()}")
            sys.exit(1)

    print(f"[INFO] SignalP results written to: {output_file}")

    # Move all *_plot.png and *_pred.txt files generated in the signalp working directory
    # Move generated artifacts from tmp dir, sanitize filenames
    for ext in ("*_plot.png", "*_pred.txt", "*_mature.fasta"):
        for artifact in tmp_dir_abs.glob(ext):
            try:
                safe_name = artifact.name.replace(":", "_")
                dest = tmp_dir_abs / safe_name
                artifact.rename(artifact.with_name(safe_name))  # rename in place
                shutil.move(str(artifact.with_name(safe_name)), dest)
                print(f"[INFO] Moved artifact to tmp: {dest}")
            except FileNotFoundError:
                print(f"[WARNING] Artifact not found: {artifact}")
            except Exception as e:
                print(f"[WARNING] Failed to move artifact {artifact.name}: {e}")


if __name__ == "__main__":
    """Main entry point for the script to enable command-line execution.
    Parses command-line arguments and runs the SignalP prediction.
    Usage:
        python run_signalp.py <input_fasta> <output_dir> --batch-size <batch_size> --signalp-path <signalp_path>
    """
    import argparse
    parser = argparse.ArgumentParser(description="Run SignalP prediction")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size")
    parser.add_argument("--signalp-path", type=Path, required=True, help="Path to signalp executable")
    parser.add_argument("--organism", type=str, default="gram+", help="Organism name, options: arch, gram+, gram-, euk")

    args = parser.parse_args()
    run(args.signalp_path, args.input_fasta, args.output_dir, args.batch_size, args.organism)
