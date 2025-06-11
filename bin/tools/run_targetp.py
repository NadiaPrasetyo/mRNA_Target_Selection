# tools/run_targetp.py
import subprocess
from pathlib import Path
import sys
import shutil
import os

def run(targetp_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 100):
    if not input_fasta.exists():
        print(f"❌ Input FASTA file does not exist: {input_fasta}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Extract first FASTA header prefix
    prefix = None
    with input_fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                prefix = line[1:].split('|')[0]
                break
    if prefix is None:
        print("❌ No FASTA header found in input file.")
        sys.exit(1)

    if not targetp_path.exists() or not targetp_path.is_file() or not os.access(targetp_path, os.X_OK):
        print(f"❌ Error: TargetP executable not found or not executable: {targetp_path}")
        sys.exit(1)

    print(f"[INFO] Running TargetP on: {input_fasta}")
    print(f"[INFO] Output prefix: {prefix}")
    print(f"[INFO] Output dir: {output_dir}")

    output_file = output_dir / f"{prefix}_targetp.txt"

    cmd = [
        str(targetp_path),
        "-fasta", str(input_fasta),
        "-org", "non-pl",
        "-format", "short",
        "-batch", str(batch_size),
        "-gff3",
        "-mature",
        "-prefix", prefix,
        "-tmp", str(tmp_dir),
        "-stdout"
    ]

    with output_file.open("w") as outfile:
        try:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ TargetP failed: {e.stderr.decode()}")
            sys.exit(1)

    # Move plot file if it exists
    plot_files = list(output_dir.glob("*_plot.png"))
    if plot_files:
        plot_file = output_dir / f"{prefix}_targetp_plot.png"
        shutil.move(str(plot_files[0]), plot_file)
        print(f"[INFO] Plot saved to: {plot_file}")

    print(f"[INFO] TargetP results written to: {output_file}")

    # move any artifacts to output directory
    # Move extra generated files
    for ext in ["_mature.fasta", ".gff3"]:
        f = Path(f"{prefix}{ext}")
        if f.exists():
            shutil.move(str(f), output_dir)

# Optional CLI interface for direct use
if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser(description="Run TargetP prediction")
    parser.add_argument("--targetp_path", type=Path, help="Path to TargetP binary")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size")

    args = parser.parse_args()
    run(args.targetp_path, args.input_fasta, args.output_dir, args.batch_size)
