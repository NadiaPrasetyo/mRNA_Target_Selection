# tools/run_signalp.py
import subprocess
from pathlib import Path
import sys
import shutil
from os.path import abspath

def run(signalp_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 10000):
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
        "-tmp", str(tmp_dir_abs)
    ]

    with output_file.open("w") as outfile:
        try:
            subprocess.run(cmd, cwd=signalp_path.parent, stdout=outfile, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(f"❌ SignalP failed: {e.stderr.decode()}")
            sys.exit(1)

    print(f"[INFO] SignalP results written to: {output_file}")

    # Move all *_plot.png and *_pred.txt files generated in the signalp working directory
    for ext in ("*_plot.png", "*_pred.txt"):
        for artifact in signalp_path.parent.glob(ext):
            try:
                dest = tmp_dir_abs / artifact.name
                shutil.move(str(artifact), dest)
                print(f"[INFO] Moved artifact to tmp: {dest}")
            except Exception as e:
                print(f"[WARNING] Failed to move artifact {artifact.name}: {e}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run SignalP prediction")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size")
    parser.add_argument("--signalp-path", type=Path, required=True, help="Path to signalp executable")

    args = parser.parse_args()
    run(args.signalp_path, args.input_fasta, args.output_dir, args.batch_size)
