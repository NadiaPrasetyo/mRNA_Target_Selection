"""
antigen_analysis.py
Command-line tool to run SignalP and TargetP predictors on input FASTA files for antigen analysis.
Overview:
    - Scans a specified pathogen sequence directory for FASTA files.
    - Runs selected prediction tools (SignalP, TargetP, Cluster, Allergenicity) on each FASTA file.
    - Supports parallel execution of jobs for efficient processing.
    - Organizes results into structured output directories.
Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    sequence_dir (str): Subdirectory under `pathogen_dir` containing FASTA files.
    --tool-root (str, required): Root directory containing tool wrappers and executables.
    --threads (int, optional): Number of parallel threads to use (default: 4).
    --tools (list, optional): List of tools to run (choices: SIGNALP, TARGETP; default: both).
    --batch-size (int, optional): Batch size for SignalP/TargetP (default: 10000).
    --output-dir (str, optional): Output directory for results (default: epitope_outputs).
Requirements:
    - Tool wrappers and executables for SignalP and TargetP available under `tool-root`.
    - Input FASTA files present in the specified sequence directory.
    - Python packages: argparse, pathlib, concurrent.futures.
Usage Example:
    python antigen_analysis.py sars_cov_2 proteins --tool-root /opt/bio_tools --threads 8 --tools SIGNALP
Outputs:
    data/<pathogen_dir>/<output_dir>/<tool>/<input_file>_<tool>.out   # Prediction results for each tool and input
Author: Nadia
"""

import argparse
import sys
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_signalp, run_targetp, run_cluster, run_algpred, run_deeplocpro, common
import shutil

# Define the mapping of tool names to their runner functions
TOOL_RUNNERS = {
    "SIGNALP": run_signalp.run,
    "TARGETP": run_targetp.run,
    "CLUSTER": run_cluster.run, 
    "ALGPRED": run_algpred.run,
    "DEEPLOC": run_deeplocpro.run_deeplocpro
}

# List of valid tools that can be run
VALID_TOOLS = list(TOOL_RUNNERS.keys())


def run_tool(tool_name, runner_func, input_file, output_dir, batch_size, tool_path):
    """
    Run a specific tool on the input file and save the output.
    Args:
        tool_name (str): Name of the tool to run (e.g., SIGNALP, TARGETP, CLUSTER).
        runner_func (callable): Function to run the tool.
        input_file (Path): Input FASTA file to process.
        output_dir (Path): Directory to save the output.
        batch_size (int): Batch size for processing.
        tool_path (Path): Path to the tool executable or script.
    """
    try:
        runner_func(tool_path, input_file, output_dir, batch_size)
        logging.info(f"✅ {tool_name} completed for {input_file.name}")
    except Exception as e:
        logging.error(f"❌ {tool_name} failed for {input_file.name}: {e}")


def run_parallel_jobs(jobs, threads):
    """ Run a list of jobs in parallel using ThreadPoolExecutor.
    Args:
        jobs (list): List of tuples containing (tool_name, runner_func, input_file, output_dir, batch_size, tool_path).
        threads (int): Number of threads to use for parallel execution.
    """

    allergenicity_jobs = [job for job in jobs if job[0] == "ALGPRED"]
    other_jobs = [job for job in jobs if job[0] != "ALGPRED"]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Run non-Allergenicity jobs in parallel
        logging.info(f"Running {len(jobs) - len(allergenicity_jobs)} jobs in parallel with {threads} threads")
        futures = [executor.submit(run_tool, *job) for job in other_jobs]
        for f in futures:
            try:
                f.result()
            except Exception as e:
                logging.error(f"❌ Job failed: {e}")

    # Run Allergenicity jobs serially
    for job in allergenicity_jobs:
        tool_name, _, input_file, _, _, _ = job
        try:
            logging.info(f"Running Allergenicity for {input_file.name}")
            run_tool(*job)
        except Exception as e:
            logging.error(f"❌ Allergenicity failed for {input_file.name}: {e}")


    # Clean up temporary and intermediate files after CLUSTER jobs
    for job in jobs:
        tool_name, _, input_file, output_dir, _, _ = job
        if tool_name != "CLUSTER":
            continue

        db_name = input_file.stem
        work_dir = output_dir / "mmseqdb"
        cluster_input_dir = output_dir.parent / "cluster_inputs"

        logging.info(f"🧹 Cleaning up intermediate files for {db_name}")

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
        
        # Delete alignment results
        aln_result_path = work_dir / f"{db_name}_aln"
        if aln_result_path.exists() and aln_result_path.is_file():
            try:
                aln_result_path.unlink()
            except Exception as e:
                logging.warning(f"⚠️ Failed to delete alignment result {aln_result_path}: {e}")

        # Remove mmseqdb dir if empty
        if work_dir.exists() and work_dir.is_dir():
            try:
                if not any(work_dir.iterdir()):
                    work_dir.rmdir()
            except Exception as e:
                logging.warning(f"⚠️ Failed to remove mmseqdb directory: {e}")
        else:
            logging.debug(f"🟡 mmseqdb directory cleanup; not found: {work_dir}")

        # Remove cluster_inputs dir
        if cluster_input_dir.exists() and cluster_input_dir.is_dir():
            try:
                shutil.rmtree(cluster_input_dir)
            except Exception as e:
                logging.warning(f"⚠️ Failed to clean cluster_inputs directory {cluster_input_dir}: {e}")

        #remove tmp in cluster output dir
        tmp_dir = output_dir / "tmp"
        if tmp_dir.exists() and tmp_dir.is_dir():
            try:
                shutil.rmtree(tmp_dir)
            except Exception as e:
                logging.warning(f"⚠️ Failed to clean tmp directory {tmp_dir}: {e}")

        logging.info(f"✅ Cleanup completed for {db_name}")


def is_job_done(fasta_path: Path, output_dir: Path, mode: str = "strain") -> bool:
    """
    Check if a job has already been processed for a given FASTA file.
    
    Modes:
        - "strain": checks if file.stem (e.g. 'HO_5096_0412_matched_antigens') appears in any output.
        - "accession": checks if accession (e.g. 'ABC123') from 'ABC123_combined.fasta' appears in output names.
    
    Args:
        fasta_path (Path): Input FASTA file path.
        output_dir (Path): Output directory to search for existing results.
        mode (str): Type of identifier to match ('strain' or 'accession').
    
    Returns:
        bool: True if matching output exists, False otherwise.
    """
    stem = fasta_path.stem

    if mode == "strain":
        identifier = stem  # e.g., 'HO_5096_0412_matched_antigens'
    elif mode == "accession":
        identifier = stem.split("_")[0]  # e.g., 'ABC123' from 'ABC123_combined'
    else: 
        raise ValueError(f"Unknown mode '{mode}' passed to is_job_done()")
    
    extensions = [".tsv", ".fasta", ".txt", ".gff3", "_algpred.csv"]

    for ext in extensions:
        for f in output_dir.glob(f"*{ext}"):
            if identifier in f.name:
                return True
    return False


def main():
    """Main function to parse arguments and run the antigen analysis pipeline."""
    parser = argparse.ArgumentParser(description="Run SignalP, TargetP, and Cluster on input FASTA files")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", help="Root directory for tools, required for SignalP and TargetP", default="none")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tools", nargs="+", choices=VALID_TOOLS, default=VALID_TOOLS)
    parser.add_argument("--batch-size", type=int, default=10000)
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"))
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for debugging")
    parser.add_argument("--group", type=str, default="any", help="Group name for DeepLocPro: [any, archaea, positive, negative]; default is 'any'")
    args = parser.parse_args()
    
    data_path = Path("data") / args.pathogen_dir
    
    output_root = data_path / args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    if args.verbose:
        log_file = output_root / "antigen_analysis.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers = [
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode='a')
        ]
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=handlers,
            force=True
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s: %(message)s",
            force=True
        )

    sequence_path = data_path / args.sequence_dir
    if not sequence_path.exists():
        logging.error(f"❌ Invalid input directory: {sequence_path}")
        sys.exit(1)
            
    fasta_files = common.get_fasta_files(data_path, args.sequence_dir)
    if not fasta_files:
        logging.error("❌ No FASTA files found.")
        sys.exit(1)

    tool_paths = {}

    print("DEBUG: fasta_files =", fasta_files)
    print("DEBUG: checking tool paths=", args.tools)
    # Only resolve tool_root if required
    if any(tool in args.tools for tool in ["SIGNALP", "TARGETP"]):
        if args.tool_root == "none":
            logging.error("❌ --tool-root is required for running SignalP and TargetP.")
            sys.exit(1)
    
    tool_root = Path(args.tool_root).resolve()
    print("DEBUG: tool_root =", tool_root)
    
    try:
        tool_paths = common.check_antigen_tools(args.tools, tool_root)
    except (FileNotFoundError, ImportError) as e:
        logging.error(f"❌ {e}")
        sys.exit(1)

    print("DEBUG: tool_paths =", tool_paths)

    jobs = []

    # Handle jobs
    for tool_name in args.tools:
        # Handle Cluster jobs
        if tool_name == "CLUSTER":
            cluster_output = output_root / "cluster"
            cluster_input_dir = output_root / "cluster_inputs"
            grouped_fastas = common.group_cluster_inputs(fasta_files, cluster_input_dir)

            for _, fasta_path in grouped_fastas.items():
                if is_job_done(fasta_path, cluster_output, mode="accession"):
                    logging.info(f"⏭️ Skipping CLUSTER for {fasta_path.name} (already processed)")
                    continue

                jobs.append(("CLUSTER", TOOL_RUNNERS["CLUSTER"], fasta_path, cluster_output, 0, tool_paths["CLUSTER"]))

        # Handle SignalP and TargetP jobs
        elif tool_name in ["SIGNALP", "TARGETP", "ALGPRED"]:
            output_dir = output_root / tool_name.lower()
            output_dir.mkdir(parents=True, exist_ok=True)
            for fasta_file in fasta_files:
                if is_job_done(fasta_file, output_dir, mode="strain"):
                    logging.info(f"⏭️ Skipping {tool_name} for {fasta_file.name} (already processed)")
                    continue

                jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_file, output_dir, args.batch_size, tool_paths[tool_name]))

        # Handle DeepLocPro jobs
        elif tool_name == "DEEPLOC":
            output_dir = output_root / "deeplocpro"
            output_dir.mkdir(parents=True, exist_ok=True)
            for fasta_file in fasta_files:
                if is_job_done(fasta_file, output_dir, mode="strain"):
                    logging.info(f"⏭️ Skipping DEEPLOC for {fasta_file.name} (already processed)")
                    continue

                jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_file, output_dir, args.group, tool_paths["DEEPLOC"]))

    if not jobs:
        logging.info("No jobs to run. All files already processed.")
        # clean up any empty directories created
        for tool_name in args.tools:
            output_dir = output_root / tool_name.lower()
            if output_dir.exists() and not any(output_dir.iterdir()):
                try:
                    output_dir.rmdir()
                except OSError as e:
                    logging.warning(f"⚠️ Failed to remove empty directory {output_dir}: {e}")
                    
        #clean up cluster inputs if they exist
        cluster_input_dir = output_root / "cluster_inputs"
        if cluster_input_dir.exists() and not any(cluster_input_dir.iterdir()):
            try:
                cluster_input_dir.rmdir()
            except OSError as e:
                logging.warning(f"⚠️ Failed to remove empty directory {cluster_input_dir}: {e}")
        return
    
    run_parallel_jobs(jobs, args.threads)


if __name__ == "__main__":
    main()
