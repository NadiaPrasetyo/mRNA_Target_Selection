"""
evaluate_epitopes.py

Command-line tool to evaluate predicted epitope sequences using multiple immunoinformatics tools.

Overview:
    - Runs evaluation tools (Allergenicity, Population Coverage, Cluster) on epitope JSON files.
    - Supports parallel execution for efficient processing of multiple jobs.
    - Handles tool-specific input preparation and output validation.
    - Skips jobs if output already exists and passes validation checks.
    - Cleans up temporary files generated during processing.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    epitope_dir (Path, required): Directory containing epitope predictions with mhci/, mhcii/, and bcell/ subdirs.
    --tool-root (str, required): Root directory containing the evaluation tools.
    --verbose (flag, optional): Enable verbose logging.
    --threads (int, optional): Number of parallel threads (default: 4).
    --tools (list, optional): Specify which tools to run (Allergenicity, PopCoverage, Cluster).
    --output-dir (Path, optional): Directory to save output files (default: evaluation_outputs).

Requirements:
    - Evaluation tools (Allergenicity, Population Coverage, Cluster) installed and available in tool-root.
    - Epitope JSON files present in the specified epitope directory.
    - Python packages: argparse, concurrent.futures, pathlib, logging, json.

Usage Example:
    python evaluate_epitopes.py sars_cov_2 --epitope-dir data/sars_cov_2/epitopes --tool-root tools/ --threads 8 --tools Allergenicity PopCoverage

Outputs:
    <output_dir>/<tool_type>/*.txt or *.json   # Output files from each evaluation tool
    <output_dir>/json_inputs/                  # Temporary JSON files for Cluster tool
    <output_dir>/popcov_inputs/                # Temporary input files for PopCoverage tool

Author: Nadia
"""
import argparse
import sys
import logging
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from tools import run_algpred, run_popcoverage, run_cluster, common, extract_epitopes

# Dictionary mapping tool names to their respective runner functions
tool_runners = {
    "Allergenicity": run_algpred.run,
    "PopCoverage": run_popcoverage.run,
    "Cluster": run_cluster.run,
}

def is_job_completed(tool_type, input_path, base_output_dir):
    """
    Check if a job has already been processed by verifying if a result file exists
    in the tool-specific output directory.
    Looks for any file that contains the input's base name and ends with the tool-specific suffix.
    Args:
        tool_type (str): Type of the tool (e.g., "Allergenicity", "PopCoverage", "Cluster").
        input_path (Path): Path to the input JSON file.
        base_output_dir (Path): Base output directory where results are stored.
    Returns:
        bool: True if the job is completed (output file exists), False otherwise.
    """
    subdir = base_output_dir / tool_type.lower()
    expected_suffix = ".json" if tool_type == "Cluster" else ".txt"

    for file in subdir.glob(f"{input_path.stem}*{expected_suffix}"):
        try:
            if not file.exists() or file.stat().st_size == 0:
                logging.warning(f"{tool_type}: Output file {file} is missing or empty.")
                return False

            if file.suffix == ".json":
                with open(file) as f:
                    data = json.load(f)
                if tool_type == "Cluster" and "clusters" not in data:
                    logging.warning(f"{tool_type}: Missing 'clusters' in {file}.")
                    return False
                if tool_type != "Cluster" and "results" not in data:
                    logging.warning(f"{tool_type}: Missing 'results' in {file}.")
                    return False

            return True
        except Exception as e:
            logging.error(f"Failed to validate output {file}: {e}")
            return False

    return False

def run_predictions_parallel(job_list, output_dir, max_threads, args):
    """
    Run the specified jobs in parallel using a thread pool executor.
    Args:
        job_list (list): List of tuples containing (tool_type, tool_path, input_file).
        output_dir (Path): Directory to save output files.
        max_threads (int): Maximum number of threads to use for parallel execution.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    logging.info(f"Starting parallel execution with {max_threads} thread(s) on {len(job_list)} jobs")
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for tool_type, tool_path, input_file in job_list:
            sub_out = output_dir / tool_type.lower()
            if tool_type == "Cluster":
                temp_dir = output_dir / "json_inputs"
                alleles = extract_epitopes.get_alleles_from_epitope_file(input_file, "mhci") + \
                          extract_epitopes.get_alleles_from_epitope_file(input_file, "mhcii")
                jsons = common.parse_fasta_to_jsons(input_file, temp_dir, alleles, [], "cluster", input_file.stem)
                for jp in jsons:
                    futures.append(executor.submit(tool_runners[tool_type], tool_path, jp, sub_out))
            elif tool_type == "PopCoverage":
                temp_txt = output_dir / "popcov_inputs"
                mhci_ep = list((args.epitope_dir / "mhci").glob("*.json"))
                mhcii_ep = list((args.epitope_dir / "mhcii").glob("*.json"))
                for tool_class, files in [("MHCI", mhci_ep), ("MHCII", mhcii_ep)]:
                    for ep in files:
                        alleles = extract_epitopes.get_alleles_from_epitope_file(ep, tool_class.lower())
                        out = temp_txt / f"{ep.stem}_{tool_class.lower()}.txt"
                        with open(ep) as fin, open(out, "w") as fout:
                            for l in fin:
                                fq = l.strip()
                                if fq:
                                    fout.write(f"{fq} {','.join(alleles)}\n")
                        futures.append(executor.submit(tool_runners[tool_type], tool_path, out, sub_out))
            else:
                futures.append(executor.submit(tool_runners[tool_type], tool_path, input_file, sub_out))

        for f in futures:
            try:
                f.result()
            except Exception as e:
                logging.error(f"Error in job: {e}")

def main():
    """
    Main function to parse command-line arguments and run the epitope evaluation pipeline.
    It checks for required directories, prepares output directories, and runs the specified tools
    on the provided epitope JSON files.
    """
    parser = argparse.ArgumentParser(description="Run evaluation tools: Allergenicity, Population Coverage, Cluster")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("epitope_dir", type=Path,
                        help="Directory containing epitope predictions with mhci/, mhcii/, and bcell/ subdirs")
    parser.add_argument("--tool-root", required=True, help="Root directory containing analysis tools")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["Allergenicity", "PopCoverage", "Cluster"], default=None,
                        help="Specify which tools to run (default: all detected tools)")
    parser.add_argument("--output-dir", type=Path, default=Path("evaluation_outputs"),
                        help="Directory to save output files (default: 'evaluation_outputs')")

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    logging.info("Starting epitope evaluation pipeline")

    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    output_dir = args.output_dir

    for p in [pathogen_path, args.tool_root, args.epitope_dir]:
        if not Path(p).exists():
            print(f"❌ Directory does not exist: {p}")
            sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    tool_map = common.check_epitope_evaluation_tools(Path(args.tool_root))
    selected_tools = set(args.tools) if args.tools else set(tool_map.keys())
    missing_tools = selected_tools - set(tool_map.keys())
    if missing_tools:
        print(f"⚠️ Requested tools not found: {', '.join(missing_tools)}")

    final_tools = {t: tool_map[t] for t in selected_tools if t in tool_map}
    if not final_tools:
        print("❌ No valid tools available. Exiting.")
        sys.exit(1)

    # Gather all JSON files from mhci, mhcii, bcell subdirectories
    epitope_json_files = []
    for sub in ["mhci", "mhcii", "bcell"]:
        subdir = args.epitope_dir / sub
        if subdir.exists():
            epitope_json_files.extend(subdir.glob("*.json"))
    if not epitope_json_files:
        print(f"❌ No epitope JSON files found in {args.epitope_dir}")
        sys.exit(1)

    _, output_dir = common.prepare_output_dirs(pathogen_path, output_dir, final_tools.keys())

    jobs = []
    for t, p in final_tools.items():
        for f in epitope_json_files:
            if is_job_completed(t, f, output_dir):
                logging.info(f"Skipping {f.name} for {t}")
                continue
            jobs.append((t, p, f))
    if not jobs:
        logging.info("No jobs to run; exiting.")
        sys.exit(0)

    logging.info(f"Launching {len(jobs)} job(s)")
    run_predictions_parallel(jobs, output_dir, args.threads, args)

    common.cleanup_temp(output_dir / "json_inputs")
    common.cleanup_temp(output_dir / "popcov_inputs")
    logging.info("All done!")

if __name__ == "__main__":
    main()
