"""
run_cluster.py

Command-line tool to cluster input FASTA sequences using MMseqs2 and clean up intermediate files.

Overview:
    - Runs MMseqs2 clustering on a given FASTA file.
    - Generates cluster assignments as TSV and representative sequences as FASTA.
    - Cleans up all intermediate MMseqs2 files after processing.

Arguments:
    _ (unused): Placeholder for compatibility with run_tool interface.
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Output directory for results.
    _batch_size (unused): Ignored for MMseqs clustering.

Requirements:
    - MMseqs2 installed and available in PATH.
    - Input FASTA file with sequences to cluster.

Outputs:
    <output_dir>/<input_basename>_clu.tsv      # Cluster assignments (TSV)
    <output_dir>/<input_basename>_clu.fasta    # Representative sequences (FASTA)

Author: Nadia
"""
import subprocess
from pathlib import Path
import logging
import shutil

def run(_, input_fasta: Path, output_dir: Path, _batch_size: int):
    """
    Run MMseqs2 clustering on the input FASTA file and clean up intermediate files.

    Args:
        _ (unused): Placeholder for compatibility with run_tool interface.
        input_fasta (Path): Path to the input FASTA file.
        output_dir (Path): Output directory for results.
        _batch_size (unused): Ignored for MMseqs clustering.
    """

    db_name = input_fasta.stem
    work_dir = output_dir / "mmseqdb"
    tmp_dir = output_dir / "tmp"
    work_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    db_path = work_dir / db_name
    clu_path = work_dir / f"{db_name}_clu"
    tsv_path = output_dir / f"{db_name}_clu.tsv"
    fasta_output = output_dir / f"{db_name}_clu.fasta"
    seqfiledb_path = clu_path.with_name(f"{db_name}_clu_seq")

    try:
        logging.info(f"📦 Creating MMseqs DB for {input_fasta.name}")
        subprocess.run(["mmseqs", "createdb", str(input_fasta), str(db_path)], check=True)

        logging.info(f"🔗 Running clustering on {input_fasta.name}")
        subprocess.run(["mmseqs", "cluster", str(db_path), str(clu_path), str(tmp_dir)], check=True)

        logging.info(f"📄 Generating TSV output")
        subprocess.run(["mmseqs", "createtsv", str(db_path), str(db_path), str(clu_path), str(tsv_path)], check=True)

        logging.info(f"📁 Creating representative FASTA output")
        subprocess.run(["mmseqs", "createseqfiledb", str(db_path), str(clu_path), str(seqfiledb_path)], check=True)
        subprocess.run(["mmseqs", "result2flat", str(db_path), str(db_path), str(seqfiledb_path), str(fasta_output)], check=True)

        logging.info(f"✅ Clustering completed for {input_fasta.name}")

    finally:
        # Cleanup all intermediate MMseqs files and directories
        logging.info(f"🧹 Cleaning up intermediate files")

        # Delete MMseqs DB and cluster files
        for file in work_dir.glob(f"{db_name}*"):
            try:
                file.unlink()
            except Exception as e:
                logging.warning(f"⚠️ Failed to delete file {file}: {e}")

            finally:
                logging.info(f"🧹 Cleaning up intermediate files")

                # Delete MMseqs DB and cluster files
                for file in work_dir.glob(f"{db_name}*"):
                    if file.exists() and file.is_file():
                        try:
                            file.unlink()
                        except Exception as e:
                            logging.warning(f"⚠️ Failed to delete file {file}: {e}")

                # Delete intermediate seqfile DB files
                for file in work_dir.glob(f"{db_name}_clu_seq*"):
                    if file.exists() and file.is_file():
                        try:
                            file.unlink()
                        except Exception as e:
                            logging.warning(f"⚠️ Failed to delete file {file}: {e}")

                # Remove tmp dir if it exists
                if tmp_dir.exists() and tmp_dir.is_dir():
                    try:
                        shutil.rmtree(tmp_dir)
                    except Exception as e:
                        logging.warning(f"⚠️ Failed to remove tmp directory: {e}")
                else:
                    logging.debug(f"🟡 Skipping tmp directory cleanup; not found: {tmp_dir}")

                # Remove mmseqdb dir if empty
                if work_dir.exists() and work_dir.is_dir():
                    try:
                        if not any(work_dir.iterdir()):
                            work_dir.rmdir()
                    except Exception as e:
                        logging.warning(f"⚠️ Failed to remove mmseqdb directory: {e}")
                else:
                    logging.debug(f"🟡 Skipping mmseqdb directory cleanup; not found: {work_dir}")
                
                logging.info(f"🧹 Cleanup completed")
