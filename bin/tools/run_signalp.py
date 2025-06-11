# tools/run_signalp.py
import subprocess
from pathlib import Path
import sys
import shutil

def patch_path_error(signalp_path: Path):
    # Paths to SignalP binaries
    script_dir = signalp_path.parent.resolve() #script directory is signalp/bin/signalp
    print(f"[CHECK!!!!!!!!!] Script directory: {script_dir}")
    signalp_bin = script_dir / "signalp"
    expected_bin_dir = script_dir / "bin"
    expected_bin = expected_bin_dir / "signalp"

    # Patch: ensure expected_bin exists for tool compatibility
    try:
        if not expected_bin.exists():
            print("[PATCH] signalp executable not found in nested bin/bin. Creating patched structure...")
            expected_bin_dir.mkdir(parents=True, exist_ok=True)

            if not signalp_bin.exists():
                raise FileNotFoundError(f"Original signalp binary not found at {signalp_bin}")

            shutil.copy(signalp_bin, expected_bin)
            print(f"[PATCH] Copied signalp binary to: {expected_bin}")
    except Exception as e:
        print(f"❌ Failed to patch SignalP binary structure: {e}")
        sys.exit(1)


def run(signalp_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 10000):
    if not signalp_path.exists():
        print(f"❌ SignalP executable not found at: {signalp_path}")
        sys.exit(1)

    if not input_fasta.exists():
        print(f"❌ Input FASTA file does not exist: {input_fasta}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Extract FASTA header prefix
    prefix = None
    with input_fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                prefix = line[1:].split('|')[0]
                break
    if prefix is None:
        print("❌ No FASTA header found in input file.")
        sys.exit(1)

    # Ensure SignalP binary is patched correctly
    patch_path_error(signalp_path)

    print(f"[INFO] Running SignalP on: {input_fasta}")
    print(f"[CHECK!!!!!!!!!!!!] SignalP path: {signalp_path}")
    print(f"[INFO] Output directory: {output_dir}")
    print(f"[INFO] Output prefix: {prefix}")

    output_file = output_dir / f"{prefix}_signalp_phred.txt"

    cmd = [
        str(signalp_path),
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
    parser.add_argument("--signalp-path", type=Path, required=True, help="Path to signalp executable")
    # then call:
    args = parser.parse_args()
    run(args.signalp_path, args.input_fasta, args.output_dir, args.batch_size)
