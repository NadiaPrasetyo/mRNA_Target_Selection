"""
antigen_analysis.py
Runner for the antigen analysis pipeline.

Overview:
    - Scans a specified pathogen sequence directory for FASTA files.
    - Runs selected prediction tools (SignalP, TargetP, Cluster, AlgPred, DeepLocPro, IfNePitope2, DeepTMHMM, MAFFT, MAFFT_RATE4SITE) on each FASTA file.
    - Supports parallel execution for efficient processing, with serial execution for tools requiring it (e.g., AlgPred, IfNePitope2).
    - Organizes results into structured output directories.
    - Cleans up intermediate and temporary files after processing.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    sequence_dir (str): Subdirectory under `pathogen_dir` containing FASTA files.
    --tool-root (str): Root directory containing tool wrappers and executables (required for SignalP, TargetP, DeepLocPro).
    --threads (int): Number of parallel threads to use (default: 4).
    --tools (list): List of tools to run (choices: SIGNALP, TARGETP, CLUSTER, ALGPRED, DEEPLOC, IFNEPITOPE2, DEEPTMHMM, MAFFT, MAFFT_RATE4SITE; default: all).
    --batch-size (int): Batch size for SignalP/TargetP (default: 10000).
    --output-dir (str): Output directory for results (default: epitope_outputs).
    --verbose: Enable verbose output for debugging.
    --group (str): Group name for DeepLocPro: [any, archaea, positive, negative] (default: any).

Requirements:
    - Tool wrappers and executables for all selected tools available under `tool-root` as needed.
    - Input FASTA files present in the specified sequence directory.
    - Python packages: argparse, pathlib, concurrent.futures, logging, shutil.

Outputs:
    data/<pathogen_dir>/<output_dir>/<tool>/<input_file>_<tool>.out   # Prediction results for each tool and input
    data/<pathogen_dir>/<output_dir>/cluster/<accession>_clu.tsv      # Cluster results
    data/<pathogen_dir>/<output_dir>/algpred/<input_file>_algpred.csv # AlgPred results
    data/<pathogen_dir>/<output_dir>/ifnepitope2/<input_file>_ifnepitope2.csv # IfNePitope2 results
    data/<pathogen_dir>/<output_dir>/deeptmhmm/<input_file>_TMRs.gff3 # DeepTMHMM results
    data/<pathogen_dir>/<output_dir>/mafft/<input_file>.fasta.tree    # MAFFT phylogenetic tree
    data/<pathogen_dir>/<output_dir>/mafft_rate4site/<input_file>.fasta # MAFFT_RATE4SITE alignment
    data/<pathogen_dir>/<output_dir>/mafft_rate4site/<input_file>.tree  # MAFFT_RATE4SITE tree
    data/<pathogen_dir>/<output_dir>/mafft_rate4site/rate4site_results/ # Rate4Site results directory

Author: Nadia
"""

import argparse
import sys
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_signalp, run_targetp, run_cluster, run_algpred, run_deeplocpro, run_ifnepitope2, run_deeptmhmm, run_mafft
from tools import common
import shutil

# Define the mapping of tool names to their runner functions
TOOL_RUNNERS = {
    "SIGNALP": run_signalp.run,
    "TARGETP": run_targetp.run,
    "CLUSTER": run_cluster.run, 
    "ALGPRED": run_algpred.run,
    "DEEPLOC": run_deeplocpro.run,
    "IFNEPITOPE2": run_ifnepitope2.run,
    "DEEPTMHMM": run_deeptmhmm.run,
    "MAFFT": run_mafft.run,
    "MAFFT_RATE4SITE": run_mafft.run
}

# List of valid tools that can be run
VALID_TOOLS = list(TOOL_RUNNERS.keys())


def run_tool(tool_name, runner_func, input_file, output_dir, batch_size, tool_path):
    """
    Run a specific tool on the input file and save the output.
    Args:
        tool_name (str): Name of the tool to run (e.g., SIGNALP, TARGETP, CLUSTER, ALGPRED, DEEPLOC, IFNEPITOPE2).
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

    serial_jobs = [job for job in jobs if job[0] == "ALGPRED" or job[0] == "IFNEPITOPE2"]
    other_jobs = [job for job in jobs if job[0] != "ALGPRED" and job[0] != "IFNEPITOPE2"]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Run non-Allergenicity jobs in parallel
        logging.info(f"Running {len(jobs) - len(serial_jobs)} jobs in parallel with {threads} threads")
        futures = [executor.submit(run_tool, *job) for job in other_jobs]
        for f in futures:
            try:
                f.result()
            except Exception as e:
                logging.error(f"❌ Job failed: {e}")

    # Run Allergenicity and IFNepitope2 jobs serially
    for job in serial_jobs:
        tool_name, _, input_file, _, _, _ = job
        try:
            logging.info(f"Running {tool_name} for {input_file.name}")
            run_tool(*job)
        except Exception as e:
            logging.error(f"❌ {tool_name} failed for {input_file.name}: {e}")


    # Clean up temporary and intermediate files after CLUSTER, MAFFT, and MAFFT_RATE4SITE jobs
    for job in jobs:
        tool_name, _, input_file, output_dir, _, _ = job
        if tool_name not in ["CLUSTER", "MAFFT_RATE4SITE", "MAFFT"]:
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
        
        empty_res = output_dir.parent.parent.parent.parent / "r4s.res"
        if empty_res.exists():
            try:
                empty_res.unlink()
            except Exception as e:
                logging.warning(f"⚠️ Failed to delete empty r4s.res file: {e}")

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

    extensions = [".tsv", ".fasta", ".txt", ".gff3", "_algpred.csv", "_ifnepitope2.csv", "_TMRs.gff3", "_probs.csv", 
                  "_results.md", "_deeptmhmm_results.md", "_predicted_topologies.3line", "_plot.png", "_combined_aligned.tree", "_combined_aligned.fasta",
                  "_combined_aligned.out", "_combined_clu.tsv", "_combined_clu.fasta", "_combined_scores.m8"]

    for ext in extensions:
        for f in output_dir.glob(f"*{ext}"):
            if identifier in f.name:
                return True
    return False


def main():
    """Main function to parse arguments and run the antigen analysis pipeline."""
    parser = argparse.ArgumentParser(description="Run selected antigen prediction tools on pathogen sequences.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", help="Root directory for tools, required for SignalP, TargetP, and DeeplocPro", default="none")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tools", nargs="+", choices=VALID_TOOLS, default=VALID_TOOLS, help="List of tools to run (default: all tools)")
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
    # Only resolve tool_root if required
    if any(tool in args.tools for tool in ["SIGNALP", "TARGETP", "DEEPLOC"]):
        if args.tool_root == "none":
            logging.error("❌ --tool-root is required for running SignalP, TargetP, and DeepLocPro.")
            sys.exit(1)
    
    tool_root = Path(args.tool_root).resolve()
    
    try:
        tool_paths = common.check_antigen_tools(args.tools, tool_root)
    except (FileNotFoundError, ImportError) as e:
        logging.error(f"❌ {e}")
        sys.exit(1)

    jobs = []

    # Handle jobs
    for tool_name in args.tools:
        # Handle Cluster jobs
        if tool_name in ["CLUSTER", "MAFFT_RATE4SITE", "MAFFT"]:
            output_dir = output_root / tool_name.lower()
            cluster_input_dir = output_root / "cluster_inputs"
            grouped_fastas = common.group_cluster_inputs(fasta_files, cluster_input_dir)

            for _, fasta_path in grouped_fastas.items():
                if is_job_done(fasta_path, output_dir, mode="accession"):
                    logging.info(f"⏭️ Skipping CLUSTER for {fasta_path.name} (already processed)")
                    continue

                if tool_name == "CLUSTER":
                    jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_path, output_dir, args.batch_size, tool_paths[tool_name]))
                elif tool_name == "MAFFT_RATE4SITE":
                    jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_path, output_dir, True, tool_paths[tool_name]))
                elif tool_name == "MAFFT":
                    jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_path, output_dir, False, tool_paths[tool_name]))

        # Handle SignalP and TargetP jobs
        elif tool_name in ["SIGNALP", "TARGETP", "ALGPRED", "DEEPTMHMM"]:
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

                jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_file, output_dir, args.group, tool_paths[tool_name]))

        # Handle IfNePitope2 jobs
        elif tool_name == "IFNEPITOPE2":
            output_dir = output_root / "ifnepitope2"
            output_dir.mkdir(parents=True, exist_ok=True)
            for fasta_file in fasta_files:
                if is_job_done(fasta_file, output_dir, mode="strain"):
                    logging.info(f"⏭️ Skipping IFNEPITOPE2 for {fasta_file.name} (already processed)")
                    continue

                jobs.append((tool_name, TOOL_RUNNERS[tool_name], fasta_file, output_dir, 3, tool_paths[tool_name]))

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
