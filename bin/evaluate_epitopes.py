"""
evaluate_epitopes.py
Command-line tool to evaluate predicted epitopes using immunoinformatics tools (e.g., PopCoverage).

Overview:
    - Scans a specified pathogen epitope directory for predicted epitope files (mhci, mhcii, bcell).
    - Runs selected evaluation tools (e.g., PopCoverage) on each epitope file.
    - Validates outputs and skips jobs if results already exist and are valid.
    - Supports parallel execution for efficient processing.
    - Organizes results into structured output directories and cleans up temporary files.

Usage Example:
    python evaluate_epitopes.py sars_cov_2 epitopes --tool-root /opt/bio_tools --threads 8 --tools PopCoverage

Requirements:
    - Tool wrappers (run_popcoverage) available under `tool-root`.
    - Input epitope files: JSON (mhci, mhcii), CSV (bcell).
    - Python packages: argparse, pathlib, concurrent.futures, logging.

Author: Nadia
"""

import argparse
import logging
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import shutil

from tools import run_popcoverage, common, extract_epitopes

tool_runners = {
    "PopCoverage": run_popcoverage.run,
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Evaluate predicted epitopes using immunoinformatics tools.")
    parser.add_argument("pathogen_dir", help="Subdirectory inside data/ for pathogen")
    parser.add_argument("epitope_dir", type=Path, help="Directory under pathogen_dir with mhci/, mhcii/, bcell/")
    parser.add_argument("--tool-root", required=True, help="Root directory with tool scripts")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=list(tool_runners.keys()), help="Tools to run")
    parser.add_argument("--output-dir", type=Path, default=Path("evaluation_outputs"), help="Directory for outputs, defaults to 'evaluation_outputs'")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
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
            ext = "*.csv" if subdir == "bcell" else "*.json"
            files.extend(path.glob(ext))
    return files

def is_output_valid(tool, input_file, output_dir):
    tool = tool.lower()
    stem = input_file.stem
    out_dir = output_dir / tool
    try:
        if tool == "popcoverage":
            return any((f.exists() and f.stat().st_size > 0)
                       for f in out_dir.glob(f"{stem}*.txt")) or \
                   any((f.exists() and f.stat().st_size > 0)
                       for f in out_dir.glob(f"{stem}*.png"))
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
        list: jobs_to_run
    """
    jobs = []
    for tool, tool_path in tools_to_run.items():
        for file in epitope_files:
            if not is_output_valid(tool, file, output_dir):
                jobs.append((tool, tool_path, file))
            else:
                logging.info(f"⏩ Skipping {tool} for {file.name}: output already exists.")
    return jobs

def run_jobs_parallel(jobs, output_dir, epitope_dir, max_threads):
    """
    Run the prepared jobs. All jobs run in parallel.
    """
    output_dir = Path(output_dir)
    epitope_dir = Path(epitope_dir)
    popcov_inputs_dir = output_dir / "popcov_inputs"

    temp_dirs = [popcov_inputs_dir]

    popcov_jobs = [job for job in jobs if job[0] == "PopCoverage"]
    other_jobs = [job for job in jobs if job[0] not in {"PopCoverage"}]

    # Prepare PopCoverage inputs
    if popcov_jobs:
        try:
            epitope_map = extract_epitopes.extract_all_epitopes_by_file(epitope_dir)
            extract_epitopes.write_allele_epitopes(epitope_map, popcov_inputs_dir)
            temp_dirs.append(popcov_inputs_dir)
        except Exception as e:
            logging.error(f"❌ Failed to prepare PopCoverage inputs: {e}")

    # Parallel for PopCoverage and other tools
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []

        for tool, tool_path, _ in popcov_jobs:
            for txt_file in popcov_inputs_dir.glob("*.txt"):
                if txt_file.exists() and txt_file.stat().st_size > 0:
                    out_dir = output_dir / tool.lower()
                    futures.append(executor.submit(tool_runners[tool], tool_path, txt_file, out_dir))

        for tool, tool_path, file in other_jobs:
            out_dir = output_dir / tool.lower()
            futures.append(executor.submit(tool_runners[tool], tool_path, file, out_dir))

        # Wait for parallel jobs to complete
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"❌ Parallel job failed: {e}")

    # Cleanup temp dirs
    for temp_dir in temp_dirs:
        try:
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
                logging.info(f"🧹 Removed temporary directory: {temp_dir}")
        except Exception as e:
            logging.warning(f"⚠️ Failed to remove {temp_dir}: {e}")

def main():
    """Main function to run the evaluation pipeline.
    Parses arguments, checks directories, discovers epitope files, prepares jobs, and runs them in parallel.
    """
    args = parse_arguments()

    # Logging configuration
    if args.verbose:
        output_dir = common.prepare_output_dirs(Path("data") / args.pathogen_dir, args.output_dir, [])
        log_file = output_dir / "pipeline.log"
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
        logging.error("❌ No valid tools found to run.")
        sys.exit(1)

    epitope_files = discover_epitope_files(epitope_path)
    if not epitope_files:
        logging.error(f"❌ No epitope files found in {epitope_path}")
        sys.exit(1)

    output_dir = common.prepare_output_dirs(pathogen_path, args.output_dir, tools_to_run.keys())
    jobs = prepare_jobs(epitope_files, tools_to_run, output_dir)

    if not jobs:
        logging.info("✅ All tasks are already complete.")
        return

    run_jobs_parallel(jobs, output_dir, epitope_path, args.threads)
    logging.info("✅ Evaluation pipeline complete.")

if __name__ == "__main__":
    main()
