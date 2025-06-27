"""
antigen_analysis.py
Command-line tool to run SignalP and TargetP predictors on input FASTA files for antigen analysis.
Overview:
    - Scans a specified pathogen sequence directory for FASTA files.
    - Runs selected prediction tools (SignalP, TargetP) on each FASTA file.
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
from tools import run_signalp, run_targetp, run_cluster, common
import shutil

# Define the mapping of tool names to their runner functions
TOOL_RUNNERS = {
    "SIGNALP": run_signalp.run,
    "TARGETP": run_targetp.run,
    "CLUSTER": run_cluster.run
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
    output_file = output_dir / f"{input_file.stem}_{tool_name.lower()}.out"
    if output_file.exists():
        logging.info(f"⏭️ Skipping {tool_name} for {input_file.name} (output exists)")
        return
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
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(run_tool, *job) for job in jobs]
        for f in futures:
            try:
                f.result()
            except Exception as e:
                logging.error(f"❌ Job failed: {e}")

    # Clean up temporary directories if needed
    for job in jobs:
        tool_name, _, input_file, output_dir, _, _ = job
        if tool_name == "CLUSTER":
            cluster_input_dir = output_dir / "cluster_inputs"
            if cluster_input_dir.exists():
                logging.info(f"🗑️ Cleaning up temporary directory: {cluster_input_dir}")
                for f in cluster_input_dir.glob("*.fasta"):
                    f.unlink()
                cluster_input_dir.rmdir()

            # Optional: cleanup mmseqdb directory
            mmseqdb_dir = output_dir / "mmseqdb"
            if mmseqdb_dir.exists():
                logging.info(f"🗑️ Cleaning up mmseqdb directory: {mmseqdb_dir}")
                shutil.rmtree(mmseqdb_dir)



def main():
    """Main function to parse arguments and run the antigen analysis pipeline."""
    parser = argparse.ArgumentParser(description="Run SignalP, TargetP, and Cluster on input FASTA files")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory for tools")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tools", nargs="+", choices=VALID_TOOLS, default=VALID_TOOLS)
    parser.add_argument("--batch-size", type=int, default=10000)
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"))
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for debugging")
    args = parser.parse_args()

    data_path = Path("data") / args.pathogen_dir
    sequence_path = data_path / args.sequence_dir
    if not sequence_path.exists():
        logging.error(f"❌ Invalid input directory: {sequence_path}")
        sys.exit(1)

    fasta_files = common.get_fasta_files(data_path, args.sequence_dir)
    if not fasta_files:
        logging.error("❌ No FASTA files found.")
        sys.exit(1)

    tool_root = Path(args.tool_root)
    if not tool_root.exists():
        logging.error(f"❌ Tool root does not exist: {tool_root}")
        sys.exit(1)

    try:
        tool_paths = common.check_antigen_tools(tool_root)
    except FileNotFoundError as e:
        logging.error(f"❌ {e}")
        sys.exit(1)

    output_root = data_path / args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    
    if args.verbose:
        # Prepare output directories first to get the correct output_dir
        log_file = output_root / "antigen_analysis.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode='a')]
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s: %(message)s", handlers=handlers)
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


    jobs = []

    # Handle jobs
    for tool_name in args.tools:
        # Handle Cluster jobs
        if tool_name == "CLUSTER":
            cluster_output = output_root / "cluster"
            cluster_input_dir = output_root / "cluster_inputs"
            grouped_fastas = common.group_cluster_inputs(fasta_files, cluster_input_dir)

            for accession, fasta_path in grouped_fastas.items():
                jobs.append(("CLUSTER", TOOL_RUNNERS["CLUSTER"], fasta_path, cluster_output, args.batch_size, tool_paths["CLUSTER"]
                ))

        # Handle SignalP and TargetP jobs
        elif tool_name in ["SIGNALP", "TARGETP"]:
            output_dir = output_root / tool_name.lower()
            output_dir.mkdir(parents=True, exist_ok=True)
            for fasta_file in fasta_files:
                jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_file, output_dir, args.batch_size, tool_paths[tool_name]))


    run_parallel_jobs(jobs, args.threads)


if __name__ == "__main__":
    main()
