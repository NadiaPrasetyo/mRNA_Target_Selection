import subprocess
from pathlib import Path

def run(tool_path: Path, input_fasta: Path, output_dir: Path,
        output_prefix: str = None, output_format: str = "tsv"):
    """
    Run MMseqs2 clustering on the input FASTA file.

    Args:
        tool_path (Path): Path to the directory containing MMseqs2 (optional, assumed in PATH).
        input_fasta (Path): Input FASTA file.
        output_dir (Path): Directory where clustering results will be saved.
        output_prefix (str): Prefix for output files (defaults to input_fasta name).
        output_format (str): Output format ('tsv' or 'json').

    Raises:
        FileNotFoundError: If input FASTA doesn't exist.
        ValueError: If output format is unsupported.
    """
    if output_format not in ("tsv", "json"):
        raise ValueError("Unsupported output format for MMseqs2. Only 'tsv' and 'json' supported.")

    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    output_dir = output_dir / "cluster"
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = output_prefix or input_fasta.stem
    db_path = output_dir / f"{prefix}_db"
    clu_path = output_dir / f"{prefix}_clu"
    tmp_path = output_dir / f"{prefix}_tmp"
    tsv_path = output_dir / f"{prefix}.tsv"

    cmds = [
        ["mmseqs", "createdb", str(input_fasta), str(db_path)],
        ["mmseqs", "cluster", str(db_path), str(clu_path), str(tmp_path)],
        ["mmseqs", "createtsv", str(db_path), str(db_path), str(clu_path), str(tsv_path)]
    ]

    print(f"🧬 Running MMseqs2 clustering: {input_fasta.name} -> {tsv_path.name}")
    try:
        for cmd in cmds:
            subprocess.run(cmd, check=True)
        if tsv_path.exists():
            print(f"✅ MMseqs2 clustering completed: {tsv_path.name}")
        else:
            print(f"⚠️ MMseqs2 output missing: {tsv_path.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MMseqs2 clustering failed: {e}")
