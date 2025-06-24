import subprocess
from pathlib import Path
import shutil

def run(tool_path: Path, input_fasta: Path, output_dir: Path,
        output_prefix: str = None):
    """
    Run MMseqs2 clustering and export cluster sequences to a FASTA file.

    Args:
        tool_path (Path): Path to MMseqs2 binary directory, or None if in PATH.
        input_fasta (Path): Input FASTA file.
        output_dir (Path): Directory where results will be saved.
        output_prefix (str): Output filename prefix (defaults to input_fasta name).
    """
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = output_prefix or input_fasta.stem
    tmp_dir = output_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    # Paths for intermediate files
    db_path = tmp_dir / f"{prefix}_db"
    clu_path = tmp_dir / f"{prefix}_clu"
    clu_seq_path = tmp_dir / f"{prefix}_clu_seq"

    # Final outputs
    output_tsv = output_dir / f"{prefix}.tsv"
    output_fasta = output_dir / f"{prefix}_clusters.fasta"

    mmseqs = str(tool_path / "mmseqs") if tool_path else "mmseqs"

    cmds = [
        [mmseqs, "createdb", str(input_fasta), str(db_path)],
        [mmseqs, "cluster", str(db_path), str(clu_path), str(tmp_dir)],
        [mmseqs, "createtsv", str(db_path), str(db_path), str(clu_path), str(output_tsv)],
        [mmseqs, "createseqfiledb", str(db_path), str(clu_path), str(clu_seq_path)],
        [mmseqs, "result2flat", str(db_path), str(db_path), str(clu_seq_path), str(output_fasta)],
    ]

    print(f"🧬 Running MMseqs2 clustering: {input_fasta.name} -> {output_tsv.name}, {output_fasta.name}")

    try:
        for cmd in cmds:
            subprocess.run(cmd, check=True)

        if output_tsv.exists() and output_fasta.exists():
            print(f"✅ Clustering completed: {output_tsv.name}, {output_fasta.name}")
        else:
            print(f"⚠️ Output missing: TSV = {output_tsv.exists()}, FASTA = {output_fasta.exists()}")

    except subprocess.CalledProcessError as e:
        print(f"❌ MMseqs2 command failed: {e}")

    # Cleanup tmp
    try:
        shutil.rmtree(tmp_dir)
        print(f"🧹 Temporary files cleaned: {tmp_dir}")
    except Exception as e:
        print(f"⚠️ Failed to clean temporary directory: {e}")