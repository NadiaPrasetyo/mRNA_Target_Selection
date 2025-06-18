import argparse
import logging
import sys
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from tools import run_algpred, run_popcoverage, run_cluster, common, extract_epitopes

# Map tool names to their corresponding runner functions
tool_runners = {
    "Allergenicity": run_algpred.run,
    "PopCoverage": run_popcoverage.run,
    "Cluster": run_cluster.run,
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Evaluate predicted epitopes using immunoinformatics tools.")
    parser.add_argument("pathogen_dir", help="Subdirectory inside data/ for pathogen")
    parser.add_argument("epitope_dir", type=Path, help="Epitope directory containing mhci/, mhcii/, and/or bcell/")
    parser.add_argument("--tool-root", required=True, help="Root directory with tool scripts")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["Allergenicity", "PopCoverage", "Cluster"], help="Tools to run")
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

def is_output_valid(tool: str, input_file: Path, output_dir: Path) -> bool:
    """
    Check if the output for a given tool and input file is valid.
    
    Args:
        tool (str): Tool name ("Cluster", "PopCoverage", "AlgPred").
        input_file (Path): Input file path.
        output_dir (Path): Base output directory.
    
    Returns:
        bool: True if the output is valid and complete, False otherwise.
    """
    tool = tool.lower()
    stem = input_file.stem
    out_dir = output_dir / tool

    try:
        if tool == "cluster":
            # Expecting a JSON file with a "clusters" key
            for file in out_dir.glob(f"{stem}*.json"):
                if file.exists() and file.stat().st_size > 0:
                    with open(file) as f:
                        data = json.load(f)
                    if "clusters" in data:
                        return True
            return False

        elif tool == "popcoverage":
            # Expecting .txt summary + .png plots
            txt_ok = False
            png_ok = False

            for file in out_dir.glob(f"{stem}*.txt"):
                if file.exists() and file.stat().st_size > 0:
                    txt_ok = True
                    break

            for file in out_dir.glob(f"{stem}*.png"):
                if file.exists() and file.stat().st_size > 0:
                    png_ok = True
                    break

            return txt_ok and png_ok

        elif tool == "algpred":
            # Expecting .csv result, either named by input or fallback as outfile.csv
            candidates = list(out_dir.glob(f"{stem}*.csv"))
            if not candidates:
                # Fallback check for default name
                fallback = out_dir / "outfile.csv"
                if fallback.exists() and fallback.stat().st_size > 0:
                    return True
                return False

            for file in candidates:
                if file.stat().st_size > 0:
                    return True
            return False

        else:
            # Unknown tool – consider invalid
            return False

    except Exception as e:
        print(f"⚠️ Error validating output for {tool} / {input_file.name}: {e}")
        return False

def prepare_jobs(epitope_files, tools_to_run, output_dir):
    """
    Prepare jobs for parallel execution based on available tools and epitope files.
    Args:
        epitope_files (list): List of epitope files to process.
        tools_to_run (dict): Dictionary of tools to run with their paths.
        output_dir (Path): Directory where outputs will be saved.
    Returns:
        list: List of jobs to run, each as a tuple (tool, tool_path, file).
    """
    jobs = []
    for tool, tool_path in tools_to_run.items():
        for file in epitope_files:
            if is_output_valid(tool, file, output_dir):
                logging.info(f"Skipping {file.name} for {tool}, already processed.")
                continue
            jobs.append((tool, tool_path, file))
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
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for tool, tool_path, file in jobs:
            out_subdir = output_dir / tool.lower()
            if tool == "Cluster":
                temp_dir = output_dir / "json_inputs"
                epitope_data = extract_epitopes.extract_all_epitopes(file)
                alleles = sorted({allele for _, allele in epitope_data})
                json_inputs = common.reformat_epitope_json_for_cluster(file, temp_dir, alleles, file.stem)
                for inp in json_inputs:
                    futures.append(executor.submit(tool_runners[tool], tool_path, inp, out_subdir))
            elif tool == "PopCoverage":
                if "bcell" in str(file):
                    continue  # skip bcell for PopCoverage
                temp_dir = output_dir / "popcov_inputs"
                epitope_data = extract_epitopes.extract_all_epitopes(epitope_dir)
                grouped = {}
                for peptide, allele in epitope_data:
                    grouped.setdefault(allele, []).append(peptide)
                for allele, peptides in grouped.items():
                    temp_file = temp_dir / f"{allele.replace('*', '').replace(':', '')}.txt"
                    with open(temp_file, "w") as f:
                        for p in peptides:
                            f.write(f"{p} {allele}\n")
                    futures.append(executor.submit(tool_runners[tool], tool_path, temp_file, out_subdir))
            else:
                futures.append(executor.submit(tool_runners[tool], tool_path, file, out_subdir))
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"Job failed: {e}")

def main():
    """Main function to run the evaluation pipeline.
    Parses arguments, checks directories, discovers epitope files, prepares jobs, and runs them in parallel.
    """
    args = parse_arguments()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")

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