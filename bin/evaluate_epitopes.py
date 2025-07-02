"""
evaluate_epitopes.py
Command-line tool to evaluate predicted epitopes using immunoinformatics tools (e.g., AlgPred, PopCoverage).

Overview:
    - Scans a specified pathogen epitope directory for predicted epitope files (mhci, mhcii, bcell).
    - Runs selected evaluation tools (Allergenicity/AlgPred, PopCoverage) on each epitope file.
    - Validates outputs and skips jobs if results already exist and are valid.
    - Supports parallel execution for efficient processing.
    - Organizes results into structured output directories and cleans up temporary files.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    epitope_dir (Path): Directory under `pathogen_dir` containing epitope predictions (mhci/, mhcii/, bcell/).
    --tool-root (str, required): Root directory containing tool wrappers and executables.
    --threads (int, optional): Number of parallel threads to use (default: 4).
    --tools (list, optional): List of tools to run (choices: Allergenicity, PopCoverage; default: all available).
    --output-dir (Path, optional): Output directory for results (default: evaluation_outputs).
    --verbose (flag, optional): Enable verbose logging.

Requirements:
    - Tool wrappers and executables for AlgPred and PopCoverage available under `tool-root`.
    - Input epitope files present in the specified epitope directory (JSON for mhci/mhcii, CSV for bcell).
    - Python packages: argparse, pathlib, concurrent.futures, logging.

Usage Example:
    python evaluate_epitopes.py sars_cov_2 epitopes --tool-root /opt/bio_tools --threads 8 --tools Allergenicity PopCoverage

Outputs:
    data/<pathogen_dir>/<output_dir>/<tool>/<input_file>_<tool>.*   # Prediction results for each tool and input
Author: Nadia
"""
import argparse
import logging
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import shutil

from tools import run_algpred, run_popcoverage, common, extract_epitopes

tool_runners = {
    "Allergenicity": run_algpred.run,
    "PopCoverage": run_popcoverage.run,
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Evaluate predicted epitopes using immunoinformatics tools.")
    parser.add_argument("pathogen_dir", help="Subdirectory inside data/ for pathogen")
    parser.add_argument("epitope_dir", type=Path, help="Epitope directory containing mhci/, mhcii/, and/or bcell/")
    parser.add_argument("--tool-root", required=True, help="Root directory with tool scripts")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=list(tool_runners.keys()), help="Tools to run")
    parser.add_argument("--output-dir", type=Path, default=Path("evaluation_outputs"), help="Output directory")
    return parser.parse_args()

def check_directories(pathogen_path, tool_root, epitope_path):
    """Check if the required directories exist.
    Args:
        pathogen_path (Path): Path to the pathogen directory.
        tool_root (str): Path to the root directory of tools.
        epitope_path (Path): Path to the epitope directory.
    Raises:
        SystemExit: If any directory does not exist.
    """
    for path in [pathogen_path, tool_root, epitope_path]:
        if not Path(path).exists():
            logging.error(f"❌ Directory does not exist: {path}")
            sys.exit(1)

def discover_epitope_files(epitope_dir):
    """Discover epitope files in the specified directory.
    Args:
        epitope_dir (Path): Path to the epitope directory.
    Returns:
        list: List of discovered epitope files.
    """
    files = []
    for subdir in ["mhci", "mhcii", "bcell"]:
        path = epitope_dir / subdir
        if path.exists():
            files.extend(path.glob("*.json" if subdir != "bcell" else "*.csv"))
    return files

def is_output_valid(tool, input_file, output_dir):
    tool = tool.lower()
    stem = input_file.stem
    out_dir = output_dir / tool

    try:
        if tool == "popcoverage":
            files = list(out_dir.glob(f"{stem}*.txt")) + list(out_dir.glob(f"{stem}*.png"))
            return any(f.exists() and f.stat().st_size > 0 for f in files)

        elif tool == "allergenicity":
            # Allergenicity uses FASTA file name based on the epitope input file's stem
            expected_output = out_dir / f"{stem}_algpred.csv"
            return expected_output.exists() and expected_output.stat().st_size > 0

    except Exception as e:
        logging.warning(f"⚠️ Error validating output for {tool} / {input_file.name}: {e}")
    return False

def prepare_jobs(epitope_files, tools_to_run, output_dir):
    """
    Prepare jobs for parallel execution based on available tools and epitope files.
    Args:
        epitope_files (list): List of epitope files to process.
        tools_to_run (dict): Dictionary of tools to run with their paths.
        output_dir (Path): Directory where outputs will be saved.
    Returns:
        tuple: (jobs_to_run, skipped_info)
    """
    jobs = []
    for tool, tool_path in tools_to_run.items():
        for file in epitope_files:
            if not is_output_valid(tool, file, output_dir):
                jobs.append((tool, tool_path, file))
            else:
                logging.info(f"Skipping {tool} for {file.name}: output already exists and is valid.")
    return jobs

def run_jobs_parallel(jobs, output_dir, epitope_dir, max_threads):
    """
    Run the prepared jobs in parallel using a thread pool.

    Args:
        jobs (list): List of jobs to run, each as a tuple (tool, tool_path, file).
        output_dir (Path): Directory where outputs will be saved.
        epitope_dir (Path): Directory containing epitope data.
        max_threads (int): Maximum number of threads to use for parallel execution.
    """
    output_dir = Path(output_dir)
    epitope_dir = Path(epitope_dir)
    fasta_inputs_dir = output_dir / "fasta_inputs"
    popcov_inputs_dir = output_dir / "popcov_inputs"
    temp_dirs = [fasta_inputs_dir, popcov_inputs_dir]

    allergenicity_jobs = [job for job in jobs if job[0] == "Allergenicity"]
    other_jobs = [job for job in jobs if job[0] == "PopCoverage"]

    epitope_map = {}
    if other_jobs:
        try:
            epitope_map = extract_epitopes.extract_all_epitopes_by_file(epitope_dir)
            extract_epitopes.write_allele_epitopes(epitope_map, popcov_inputs_dir)
        except Exception as e:
            logging.error(f"❌ Failed to extract/write PopCoverage epitopes: {e}")

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for _, tool_path, _ in allergenicity_jobs:
            fasta_inputs_dir.mkdir(parents=True, exist_ok=True)

        for tool, tool_path, file in other_jobs:
            if "bcell" in str(file).lower():
                continue  # Skip B-cell for PopCoverage
            for txt_file in popcov_inputs_dir.glob("*.txt"):
                if txt_file.exists() and txt_file.stat().st_size > 0:
                    out_dir = output_dir / tool.lower()
                    futures.append(executor.submit(tool_runners[tool], tool_path, txt_file, out_dir))

        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"❌ Job failed: {e}")

    # Run Allergenicity jobs serially
    for _, tool_path, file in allergenicity_jobs:
        out_dir = output_dir / "allergenicity"
        try:
            fasta_file = (
                common.parse_csv_to_fasta(file, fasta_inputs_dir, file.stem)
                if "bcell" in str(file).lower()
                else common.parse_json_to_fasta(file, fasta_inputs_dir, file.stem)
            )
            if fasta_file and fasta_file.exists():
                logging.info(f"🚀 Running Allergenicity: {file.name}")
                tool_runners["Allergenicity"](tool_path, fasta_file, out_dir)
            else:
                logging.warning(f"⚠️ FASTA not created for: {file.name}")
        except Exception as e:
            logging.error(f"❌ Allergenicity job failed for {file.name}: {e}")

    # Cleanup
    for temp_dir in temp_dirs:
        try:
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
                logging.info(f"🧹 Cleaned temporary directory: {temp_dir}")
        except Exception as e:
            logging.warning(f"⚠️ Cleanup failed for {temp_dir}: {e}")

def main():
    """Main function to run the evaluation pipeline.
    Parses arguments, checks directories, discovers epitope files, prepares jobs, and runs them in parallel.
    """
    args = parse_arguments()

    if args.verbose:
        # Prepare output directories first to get the correct output_dir
        _, output_dir = common.prepare_output_dirs(Path("data") / args.pathogen_dir, args.output_dir, [])
        log_file = output_dir / "pipeline.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers = [logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode='a')]
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s: %(message)s", handlers=handlers)
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    pathogen_path = Path("data") / args.pathogen_dir
    epitope_path = pathogen_path / args.epitope_dir
    check_directories(pathogen_path, args.tool_root, epitope_path)

    available_tools = common.check_epitope_evaluation_tools(Path(args.tool_root))
    requested_tools = set(args.tools) if args.tools else set(available_tools.keys())
    tools_to_run = {t: available_tools[t] for t in requested_tools if t in available_tools}

    if not tools_to_run:
        logging.error("No valid tools found to run.")
        sys.exit(1)

    epitope_files = discover_epitope_files(epitope_path)
    if not epitope_files:
        logging.error(f"No epitope files found in {epitope_path}")
        sys.exit(1)

    _, output_dir = common.prepare_output_dirs(pathogen_path, args.output_dir, tools_to_run.keys())
    jobs = prepare_jobs(epitope_files, tools_to_run, output_dir)

    if not jobs:
        logging.info("No jobs to run. All tasks are up-to-date.")
        return

    run_jobs_parallel(jobs, output_dir, epitope_path, args.threads)
    common.cleanup_temp(output_dir / "json_inputs")
    common.cleanup_temp(output_dir / "popcov_inputs")
    logging.info("✅ Evaluation pipeline complete.")

if __name__ == "__main__":
    main()
