# tools/run_tmhmm.py
import subprocess
from pathlib import Path
import sys
import shutil
import os

def patch_shebang(file_path: Path, perl_path: str):
    """Ensure the shebang line of file_path points to perl_path."""
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
    # Resolve absolute paths
    tmhmm_dir = tmhmm_path.resolve()
    input_fasta = input_fasta.resolve()
    output_dir = output_dir.resolve()

    tmhmm_script = tmhmm_dir / "bin" / "tmhmm"

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

    cmd = [str(tmhmm_script), "-long", str(input_fasta)]

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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run TMHMM prediction")
    parser.add_argument("tmhmm_dir", type=Path, help="Path to TMHMM2.0a directory")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")

    args = parser.parse_args()
    run(args.tmhmm_dir, args.input_fasta, args.output_dir)
