# tools/run_tmhmm.py
"""
run_tmhmm.py
Command-line tool to run TMHMM transmembrane helix prediction on protein FASTA files.

Overview:
    - Ensures TMHMM Perl scripts have correct shebang lines for the current environment.
    - Validates presence of required TMHMM files and scripts.
    - Executes TMHMM prediction on the provided input FASTA file.
    - Processes and annotates TMHMM output, appending '_tmhmm' to sequence identifiers.
    - Saves results to a specified output directory and organizes TMHMM-generated artifacts.

Arguments:
    tmhmm_dir (Path): Path to the TMHMM2.0a installation directory.
    input_fasta (Path): Input protein FASTA file for TMHMM prediction.
    output_dir (Path): Directory to save TMHMM prediction results and artifacts.

Requirements:
    - TMHMM2.0a installed and accessible at the specified path.
    - Perl interpreter available in the system PATH.
    - Python packages: argparse, pathlib, subprocess, shutil, os.

Usage Example:
    python run_tmhmm.py /opt/TMHMM2.0a/bin/tmhmm input.fasta results/tmhmm

Outputs:
    - Writes TMHMM prediction results to <output_dir>/<input_basename>_tmhmm_result.txt.
    - Moves TMHMM_* artifact directories to <output_dir>/tmp/.
    - Prints status messages and error diagnostics to the console.

Author: Nadia
"""
import subprocess
from pathlib import Path
import sys
import shutil

def patch_shebang(file_path: Path, perl_path: str):
    """Ensure the shebang line of file_path points to perl_path.
    If the file does not exist or is empty, it will skip patching.
    Args:
        file_path (Path): Path to the script file to patch.
        perl_path (str): Path to the Perl interpreter."""
    if not file_path.is_file():
        print(f"Warning: {file_path} not found, skipping shebang patch.")
        return
    with file_path.open("r") as f:
        lines = f.readlines()
    if not lines:
        print(f"Warning: {file_path} is empty, skipping shebang patch.")
        return
    current_shebang = lines[0].strip()
    expected_shebang = f"#!{perl_path}"
    if current_shebang != expected_shebang:
        print(f"Patching shebang in {file_path} from {current_shebang} to {expected_shebang}")
        lines[0] = expected_shebang + "\n"
        file_path.write_text("".join(lines))
    else:
        print(f"Shebang in {file_path} is already correct.")


def run(tmhmm_path: Path, input_fasta: Path, output_dir: Path): 
    """Run TMHMM prediction on the provided input FASTA file.
    Args:
        tmhmm_path (Path): Path to the TMHMM2.0a installation directory.
        input_fasta (Path): Input protein FASTA file for TMHMM prediction.
        output_dir (Path): Directory to save TMHMM prediction results and artifacts.
    """
    # Resolve absolute paths
    tmhmm_dir = tmhmm_path.parent.parent.resolve() 
    input_fasta = input_fasta.resolve()
    output_dir = output_dir.resolve()

    tmhmm_script = tmhmm_path.resolve()

    if not tmhmm_script.is_file():
        print(f"❌ Error: tmhmm script not found in {tmhmm_dir / 'bin'}")
        sys.exit(1)

    perl_path = shutil.which("perl")
    if perl_path is None:
        print("❌ Error: perl not found in PATH")
        sys.exit(1)

    # Patch shebangs
    patch_shebang(tmhmm_dir / "bin" / "tmhmm", perl_path)
    patch_shebang(tmhmm_dir / "bin" / "tmhmmformat.pl", perl_path)

    required_files = [
        tmhmm_dir / "lib" / "TMHMM2.0.model",
        tmhmm_dir / "bin" / "tmhmmformat.pl"
    ]
    for f in required_files:
        if not f.is_file():
            print(f"❌ Missing required file: {f}")
            sys.exit(1)

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    basename = input_fasta.stem
    output_file = output_dir / f"{basename}_tmhmm_result.txt"

    cmd = [str(tmhmm_script), str(input_fasta)]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running TMHMM: {e.stderr}")
        sys.exit(1)

    lines = []
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            lines.append(line)
        else:
            fields = line.split("\t")
            if fields:
                fields[0] = fields[0] + "_tmhmm"
                lines.append("\t".join(fields))
            else:
                lines.append(line)

    output_file.write_text("\n".join(lines) + "\n")
    print(f"✔ TMHMM run complete. Output saved to {output_file}")

    # move any artifacts to output directory
    # Move any TMHMM_* dirs into tmp/
    for item in Path.cwd().glob("TMHMM_*"):
        if item.is_dir():
            tmp_dest = output_dir / "tmp"
            tmp_dest.mkdir(exist_ok=True)
            shutil.move(str(item), tmp_dest / item.name)

if __name__ == "__main__":
    """Main entry point for the script to enable command-line execution.
    Parses command-line arguments and runs the TMHMM prediction.
    Usage:
    python run_tmhmm.py /path/to/tmhmm2.0a/bin/tmhmm input.fasta output_dir
    """
    import argparse
    parser = argparse.ArgumentParser(description="Run TMHMM prediction")
    parser.add_argument("tmhmm_dir", type=Path, help="Path to TMHMM2.0a directory")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")

    args = parser.parse_args()
    run(args.tmhmm_dir, args.input_fasta, args.output_dir)
