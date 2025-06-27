"""
run_cluster.py
Run MMseqs2 clustering tool on a FASTA file.
Compatible with tool_runners[tool] pattern and concurrent.futures.
Author: Nadia
"""

import subprocess
from pathlib import Path
import shutil

def validate_fasta(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        print(f"❌ FASTA file missing or empty: {path}")
        return False

    seen_ids = set()
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

    with path.open() as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq_id = line[1:].strip()
                if seq_id in seen_ids:
                    print(f"❌ Duplicate sequence ID at line {i}: {seq_id}")
                    return False
                seen_ids.add(seq_id)
            else:
                invalid_chars = set(line.upper()) - valid_amino_acids
                if invalid_chars:
                    print(f"❌ Invalid characters in sequence at line {i}: {line}")
                    print(f"   → Invalid: {''.join(sorted(invalid_chars))}")
                    return False
    return True

def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size=None):
    """
    Run MMseqs2 clustering on the given FASTA file.
    Only outputs final TSV and cluster FASTA. Cleans intermediate files.
    """
    if not input_fasta.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    if not validate_fasta(input_fasta):
        print(f"🚫 Validation failed for: {input_fasta}")
        raise ValueError(f"Invalid FASTA format or empty sequences in: {input_fasta}")

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = input_fasta.stem

    # Intermediate file names
    db_path = output_dir / f"{prefix}_db"
    clu_path = output_dir / f"{prefix}_clu"
    tmp_path = output_dir / f"{prefix}_tmp"
    clu_seq_path = output_dir / f"{prefix}_clu_seq"

    # Final outputs
    output_tsv = output_dir / f"{prefix}.tsv"
    output_fasta = output_dir / f"{prefix}_clusters.fasta"

    cmds = [
        ["mmseqs", "createdb", str(input_fasta), str(db_path)],
        ["mmseqs", "cluster", str(db_path), str(clu_path), str(tmp_path)],
        ["mmseqs", "createtsv", str(db_path), str(db_path), str(clu_path), str(output_tsv)],
        ["mmseqs", "createseqfiledb", str(db_path), str(clu_path), str(clu_seq_path)],
        ["mmseqs", "result2flat", str(db_path), str(db_path), str(clu_seq_path), str(output_fasta)],
    ]

    print(f"\n🧬 Running MMseqs2 clustering:")
    print(f"   Input:  {input_fasta}")
    print(f"   Output: {output_fasta}")

    try:
        for cmd in cmds:
            print(f"⚙️  Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True)  # no text=True

            # Decode with 'replace' to handle invalid UTF-8 bytes gracefully
            stdout = result.stdout.decode('utf-8', errors='replace') if result.stdout else ''
            stderr = result.stderr.decode('utf-8', errors='replace') if result.stderr else ''

            if result.returncode != 0:
                print(f"❌ Command failed: {' '.join(cmd)}")
                print(f"🔻 STDERR:\n{stderr}")
                print(f"🔻 STDOUT:\n{stdout}")
                raise subprocess.CalledProcessError(result.returncode, cmd)

        if output_tsv.exists() and output_fasta.exists():
            print(f"✅ Clustering completed: {output_tsv.name}, {output_fasta.name}")
        else:
            print(f"⚠️ Output missing: TSV = {output_tsv.exists()}, FASTA = {output_fasta.exists()}")

    except subprocess.CalledProcessError as e:
        print(f"❌ MMseqs2 pipeline failed.")
        raise e  # Optionally keep intermediate files for debugging

    # Clean up intermediate files
    for path in [db_path, clu_path, tmp_path, clu_seq_path]:
        try:
            if path.exists():
                shutil.rmtree(path)
        except Exception as e:
            print(f"⚠️ Failed to delete {path}: {e}")
