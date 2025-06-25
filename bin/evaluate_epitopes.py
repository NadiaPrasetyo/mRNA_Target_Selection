import argparse
import logging
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import shutil
import re

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
            for file in out_dir.glob(f"{stem}*.tsv"):
                if file.exists() and file.stat().st_size > 0:
                    return True
            return False

        elif tool == "popcoverage":
            basename = input_file.stem
            out_dir = output_dir / "popcoverage"
            files = list(out_dir.glob(f"{basename}*.txt")) + list(out_dir.glob(f"{basename}*.png"))

            return any(f.exists() and f.stat().st_size > 0 for f in files)

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
        return 
    

def group_cluster_inputs(files, fasta_inputs_dir: Path) -> dict:
    """
    Groups input files (JSON or CSV), converts them to FASTA using helper functions,
    and merges all FASTAs in each group into a single group FASTA file.

    Args:
        files (list[Path]): List of input JSON/CSV files.
        fasta_inputs_dir (Path): Directory where output FASTAs will be stored.

    Returns:
        dict: Mapping from group name -> Path to combined FASTA file.
    """
    grouped = defaultdict(list)
    combined_fastas = {}

    BCELL_METHODS = [
        "Chou-Fasman", "Emini", "Karplus-Schulz",
        "Kolaskar-Tongaonkar", "Parker", "Bepipred"
    ]

    # Step 1: Group input files
    for f in files:
        name = f.name.lower()
        lower_path = str(f).lower()

        if "bcell" in lower_path:
            matched = False
            for method in BCELL_METHODS:
                if method.lower().replace(" ", "") in name.replace("_", "").replace("-", "").lower():
                    grouped[f"bcell_{method}"].append(f)
                    matched = True
                    break
            if not matched:
                grouped["bcell_unknown"].append(f)

        elif "mhci" in lower_path or "mhcii" in lower_path:
            directory = "mhcii" if "mhcii" in lower_path else "mhci"
            parts = name.split("_")
            antigen_number = parts[1] if len(parts) > 1 else "unknown"
            antigen_id = parts[2] if len(parts) > 2 else "unknown"
            key = f"{directory}_antigen{antigen_number}_{antigen_id}"
            grouped[key].append(f)            

    # Step 2: Convert files to FASTA and merge per group
    for group_name, file_list in grouped.items():
        group_fasta_dir = fasta_inputs_dir / group_name
        group_fasta_dir.mkdir(parents=True, exist_ok=True)
        all_fasta_lines = []

        for i, file_path in enumerate(file_list):
            basename_prefix = f"{group_name}_{i}"
            try:
                if file_path.suffix == ".json":
                    fasta_path = common.parse_json_to_fasta(file_path, group_fasta_dir, basename_prefix)
                elif file_path.suffix == ".csv":
                    fasta_path = common.parse_csv_to_fasta(file_path, group_fasta_dir, basename_prefix)
                else:
                    logging.warning(f"⚠️ Skipping unsupported file type: {file_path}")
                    continue

                if not fasta_path or not fasta_path.exists():
                    logging.error(f"❌ FASTA not created for: {file_path.name}")
                    continue

                with open(fasta_path, "r") as f:
                    lines = f.read().strip().splitlines()
                if not lines:
                    logging.warning(f"⚠️ FASTA file is empty: {fasta_path.name}")
                    continue

                all_fasta_lines.extend(lines)
                logging.info(f"✅ FASTA added: {fasta_path.name}")

            except Exception as e:
                logging.error(f"❌ Failed to parse and convert {file_path.name}: {e}")

        if all_fasta_lines:
            combined_path = fasta_inputs_dir / f"{group_name}.fasta"
            with open(combined_path, "w") as out_f:
                out_f.write("\n".join(all_fasta_lines))
            logging.info(f"🔗 Combined FASTA written: {combined_path}")
            combined_fastas[group_name] = combined_path
        else:
            logging.warning(f"⚠️ No valid FASTA entries for group {group_name}")

    return combined_fastas

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
    skipped = {tool: [] for tool in tools_to_run.keys()}
    unprocessed = {tool: [] for tool in tools_to_run.keys()}

    for tool, tool_path in tools_to_run.items():
        for file in epitope_files:
            if is_output_valid(tool, file, output_dir):
                logging.info(f"Skipping {file.name} for {tool}, already processed.")
                skipped[tool].append(file)
                continue
            jobs.append((tool, tool_path, file))
            unprocessed[tool].append(file)

    return jobs, unprocessed


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

    popcov_inputs_dir = output_dir / "popcov_inputs"
    fasta_inputs_dir = output_dir / "fasta_inputs"
    temp_dirs = [popcov_inputs_dir, fasta_inputs_dir]

    futures = []

    # Extract epitopes for PopCoverage
    epitope_map = {}
    if any(tool == "PopCoverage" for tool, _, _ in jobs):
        try:
            epitope_map = extract_epitopes.extract_all_epitopes_by_file(epitope_dir)
            extract_epitopes.write_allele_epitopes(epitope_map, popcov_inputs_dir)
        except Exception as e:
            logging.error(f"❌ Failed to extract/write PopCoverage epitopes: {e}")
    else:
        logging.info("ℹ️ No PopCoverage jobs found, skipping epitope extraction.")

    with ThreadPoolExecutor(max_workers=max_threads) as executor:

        # 🧬 Handle Cluster jobs
        cluster_jobs = [file for tool, _, file in jobs if tool == "Cluster"] # filter for Cluster jobs
        if cluster_jobs: # if there are Cluster jobs
            grouped_fastas = group_cluster_inputs(cluster_jobs, fasta_inputs_dir)
            cluster_out_dir = output_dir / "cluster"
            cluster_out_dir.mkdir(parents=True, exist_ok=True)

            for group_name, combined_fasta in grouped_fastas.items():
                if combined_fasta.exists():
                    logging.info(f"🚀 Submitting Cluster job: {group_name} ({combined_fasta})")

                    futures.append(
                        executor.submit(
                            tool_runners["Cluster"], None, combined_fasta, cluster_out_dir, group_name
                        )
                    )
                else:
                    logging.warning(f"⚠️ Missing combined FASTA for group: {group_name}")

        # 🧪 Handle all other jobs
        for tool, tool_path, file in jobs:
            if tool == "Cluster":
                continue  # Already handled

            file = Path(file)
            out_dir = output_dir / tool.lower()

            try:
                if tool == "Allergenicity":
                    if "bcell" in str(file).lower():
                        fasta_file = common.parse_csv_to_fasta(file, fasta_inputs_dir, file.stem)
                    else:
                        fasta_file = common.parse_json_to_fasta(file, fasta_inputs_dir, file.stem)

                    if fasta_file and fasta_file.exists():
                        futures.append(executor.submit(tool_runners[tool], tool_path, fasta_file, out_dir))
                    else:
                        logging.warning(f"⚠️ FASTA not created for: {file.name}")

                elif tool == "PopCoverage":
                    if "bcell" in str(file).lower():
                        continue  # Skip B-cell for PopCoverage

                    for txt_file in popcov_inputs_dir.glob("*.txt"):
                        if txt_file.exists() and txt_file.stat().st_size > 0:
                            futures.append(executor.submit(tool_runners[tool], tool_path, txt_file, out_dir))
                        else:
                            logging.warning(f"⚠️ Skipping empty PopCoverage input: {txt_file.name}")

                else:
                    futures.append(executor.submit(tool_runners[tool], tool_path, file, out_dir))

            except Exception as e:
                logging.error(f"❌ Error preparing job for {file.name}: {e}")

        # ✅ Wait for all submitted jobs
        for future in futures:
            try:
                future.result()
            except Exception as e:
                logging.error(f"❌ Job failed: {e}")

    # # 🧹 Clean up temporary directories
    # for temp_dir in temp_dirs:
    #     try:
    #         if temp_dir.exists():
    #             shutil.rmtree(temp_dir)
    #             logging.info(f"🧹 Cleaned temporary directory: {temp_dir}")
    #     except Exception as e:
    #         logging.warning(f"⚠️ Cleanup failed for {temp_dir}: {e}")
            
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
        handlers = [logging.StreamHandler(sys.stdout)]
        handlers.append(logging.FileHandler(log_file, mode='a'))
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=handlers
        )
    else:
        # Default logging to console only
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s: %(message)s"
        )

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
    jobs, unprocessed = prepare_jobs(epitope_files, tools_to_run, output_dir)

    # # Log unprocessed files clearly for debugging
    # for tool, files in unprocessed.items():
    #     if files:
    #         logging.info(f"🔍 Unprocessed files for {tool}:")
    #         for f in files:
    #             logging.info(f"  - {f}")

    if not jobs:
        logging.info("No jobs to run. All tasks are up-to-date.")
        return

    run_jobs_parallel(jobs, output_dir, epitope_path, args.threads)

    common.cleanup_temp(output_dir / "json_inputs")
    common.cleanup_temp(output_dir / "popcov_inputs")
    logging.info("✅ Evaluation pipeline complete.")

if __name__ == "__main__":
    main()