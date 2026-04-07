"""
run_deeptmhmm.py
Run DeepTMHMM to predict membrane topology of proteins.
This script uses the BioLib library to run DeepTMHMM predictions on a given FASTA file.

Overview:
    - Ensures the required BioLib library is available.
    - Runs DeepTMHMM on the provided input FASTA file.
    - Outputs membrane topology predictions as files in the specified output directory.
    - Automatically splits large FASTA files and merges outputs.
Arguments:
    tool_path (Path): Directory containing tools (not used in this script).
    input_fasta (Path): Path to the input FASTA file.
    output_dir (Path): Directory where output will be saved.
Requirements:
    - BioLib library installed in the Python environment.
Outputs:
    <output_dir>/<input_fasta_stem>_deeptmhmm_results/   # Membrane topology prediction results
"""
import logging
from pathlib import Path
import tempfile
import shutil
import traceback
import biolib
from Bio import SeqIO

MAX_SEQS_PER_SPLIT = 300
MAX_FILE_SIZE = 2 * 1024 * 1024  # 5 MB


def split_fasta_file(input_fasta: Path, output_dir: Path, max_seqs=MAX_SEQS_PER_SPLIT):
    """Split a FASTA into smaller files (returns list of Paths)."""
    split_dir = Path(output_dir) / f"{input_fasta.stem}_splits"
    split_dir.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    split_files = []
    for i in range(0, len(records), max_seqs):
        split_path = split_dir / f"{input_fasta.stem}_part{i//max_seqs+1}.fasta"
        SeqIO.write(records[i:i+max_seqs], str(split_path), "fasta")
        split_files.append(split_path)
    return split_files, split_dir


def merge_files_by_extension(files, output_dir: Path, base_stem: str):
    """
    Merge multiple DeepTMHMM output files by extension.
    Returns a dict mapping extension -> merged file path.
    """
    ext_groups = {}
    for f in files:
        ext = f.suffix  # includes the dot, e.g., ".md", ".3line", ".gff3"
        ext_groups.setdefault(ext, []).append(f)

    merged_paths = {}
    for ext, group_files in ext_groups.items():
        merged_path = output_dir / f"{base_stem}_deeptmhmm_merged{ext}"
        with open(merged_path, "w") as out_f:
            for i, f in enumerate(group_files):
                with open(f) as in_f:
                    for line in in_f:
                        # skip duplicate headers for .md files
                        if i > 0 and ext == ".md" and line.startswith("#"):
                            continue
                        out_f.write(line)
        logging.info(f"Merged {len(group_files)} {ext} files -> {merged_path.name}")
        merged_paths[ext] = merged_path
    return merged_paths


def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int = 0):
    """Run DeepTMHMM via BioLib, with automatic splitting and type-aware merging."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"🚀 Running DeepTMHMM on {input_fasta.name}")

        # --- DEBUG LOGGING ---
        logging.info(f"Input FASTA path: {input_fasta.resolve()}")
        logging.info(f"Output directory: {output_dir.resolve()}")
        logging.info(f"Tool path (unused): {tool_path}")
        logging.info(f"Batch size: {batch_size}")

        # Check file size and sequence count
        fasta_size = input_fasta.stat().st_size
        seq_count = sum(1 for _ in SeqIO.parse(str(input_fasta), "fasta"))
        logging.info(f"{input_fasta.name} has {seq_count} sequences ({fasta_size} bytes)")

        if fasta_size > MAX_FILE_SIZE or seq_count > MAX_SEQS_PER_SPLIT:
            logging.warning(f"{input_fasta.name} is large — splitting into smaller chunks.")
            split_files, split_dir = split_fasta_file(input_fasta, output_dir)
        else:
            split_files, split_dir = [input_fasta], None

        # Load DeepTMHMM
        logging.info("Loading DeepTMHMM from BioLib...")
        deeptmhmm = biolib.load('DTU/DeepTMHMM')
        logging.info("DeepTMHMM successfully loaded.")

        all_result_files = []
        #  Run each split
        for split_fasta in split_files:
            logging.info(f"Running DeepTMHMM on chunk: {split_fasta.name}")

            # Use BioLib's file input dict instead of passing a local path string.
            with open(split_fasta, 'rb') as fasta_handle:
                deeptmhmm_job = deeptmhmm.cli(
                    args='--fasta sequences.fasta',
                    files={'sequences.fasta': fasta_handle.read()}  # .read() to pass bytes, not a handle
                )

            logging.info(f"DeepTMHMM job object: {deeptmhmm_job}")
            logging.info(f"Job status: {getattr(deeptmhmm_job, 'status', 'unknown')}")
            logging.info(f"DeepTMHMM stdout:\n{deeptmhmm_job.get_stdout()}")
            logging.info(f"DeepTMHMM stderr:\n{deeptmhmm_job.get_stderr()}")

            if b"too large for the app" in (deeptmhmm_job.get_stdout() or b''):
                logging.error(f"{split_fasta.name} too large for DeepTMHMM. Skipping.")
                continue

            existing_files = set(output_dir.glob("*"))
            logging.info(f"Existing files before save: {[f.name for f in existing_files]}")

            logging.info("Saving DeepTMHMM output files...")
            deeptmhmm_job.save_files(str(output_dir))

            new_files = set(output_dir.glob("*")) - existing_files
            logging.info(f"New files detected after save: {[f.name for f in new_files]}")

            # Rename new files to include split stem
            split_stem = split_fasta.stem
            renamed_files = []
            for file in new_files:
                if file.is_file():
                    new_name = output_dir / f"{split_stem}_{file.name}"
                    logging.info(f"Renaming {file.name} -> {new_name.name}")
                    file.rename(new_name)
                    renamed_files.append(new_name)
            all_result_files.extend(renamed_files)

        # Merge results by file type if multiple chunks
        if len(split_files) > 1 and all_result_files:
            merged = merge_files_by_extension(all_result_files, output_dir, input_fasta.stem)
            logging.info(f"✅ Merged results for {input_fasta.name}: {', '.join(f.name for f in merged.values())}")
            # delete all the pre-merged files from splits
            for f in all_result_files:
                f.unlink()
                logging.info(f"Deleted pre-merged file: {f.name}")

        elif all_result_files:
            logging.info(f"✅ Single result (no merge needed): {', '.join(f.name for f in all_result_files)}")
        else:
            logging.warning(f"No results produced for {input_fasta.name}")

        # Cleanup temporary splits
        if split_dir:
            shutil.rmtree(split_dir)

    except Exception as e:
        logging.error(f"❌ DeepTMHMM failed on {input_fasta.name}: {e}")
        logging.info("Full traceback:\n" + traceback.format_exc())
        raise



