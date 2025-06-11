# tools/run_signalp.py
import subprocess
from pathlib import Path
import sys
import shutil

def run_signalp(input_fasta: Path, output_dir: Path, batch_size: int = 10000):
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

    # Paths to SignalP binaries
    script_dir = Path(__file__).parent.resolve()
    signalp_dir = script_dir / "signalp-5.0b"
    signalp_bin = signalp_dir / "bin" / "signalp"
    expected_bin_dir = signalp_dir / "bin" / "bin"
    expected_bin = expected_bin_dir / "signalp"

    if not expected_bin.exists():
        expected_bin_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(signalp_bin, expected_bin)

    print(f"[INFO] Running SignalP on: {input_fasta}")
    print(f"[INFO] Output directory: {output_dir}")
    print(f"[INFO] Output prefix: {prefix}")

    output_file = output_dir / f"{prefix}_signalp_phred.txt"

    cmd = [
        str(expected_bin),
        "-fasta", str(input_fasta),
        "-format", "long",
        "-mature",
        "-prefix", prefix,
        "-batch", str(batch_size),
        "-stdout",
        "-tmp", str(tmp_dir)
    ]

    with output_file.open("w") as outfile:
        try:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ SignalP failed: {e.stderr.decode()}")
            sys.exit(1)

    # Move plot file if it exists
    plot_files = list(output_dir.glob("*_plot.png"))
    if plot_files:
        plot_file = output_dir / f"{prefix}_plot.png"
        shutil.move(str(plot_files[0]), plot_file)
        print(f"[INFO] Plot saved to: {plot_file}")

    print(f"[INFO] SignalP results written to: {output_file}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run SignalP prediction")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size")

    args = parser.parse_args()
    run_signalp(args.input_fasta, args.output_dir, args.batch_size)
