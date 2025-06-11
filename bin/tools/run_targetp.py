# tools/run_targetp.py
import subprocess
from pathlib import Path
import sys
import shutil

def run_targetp(input_fasta: Path, output_dir: Path, batch_size: int = 100):
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

    # Locate TargetP binary
    script_dir = Path(__file__).parent.resolve()
    targetp_bin = script_dir / "targetp-2.0" / "bin" / "targetp"

    if not targetp_bin.exists() or not targetp_bin.is_file() or not os.access(targetp_bin, os.X_OK):
        # fallback to system targetp in PATH
        from shutil import which
        system_targetp = which("targetp")
        if system_targetp is None:
            print(f"❌ Error: TargetP executable not found at {targetp_bin} or in PATH.")
            sys.exit(1)
        targetp_bin = Path(system_targetp)

    print(f"[INFO] Running TargetP on: {input_fasta}")
    print(f"[INFO] Output prefix: {prefix}")
    print(f"[INFO] Output dir: {output_dir}")

    output_file = output_dir / f"{prefix}_targetp.txt"

    cmd = [
        str(targetp_bin),
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


if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser(description="Run TargetP prediction")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size")

    args = parser.parse_args()
    run_targetp(args.input_fasta, args.output_dir, args.batch_size)
