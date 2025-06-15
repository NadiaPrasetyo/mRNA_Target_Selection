import argparse
from collections import Counter
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_allergenicity, run_population_coverage, run_conservation, common
import sys

# Define which function runs which tool
tool_runners = {
    "Allergenicity": run_allergenicity.run,
    "PopulationCoverage": run_population_coverage.run,
    "Conservation": run_conservation.run
}

def is_job_completed(tool_type, input_path, base_output_dir):
    """
    Check if a job for the tool and input has already been processed.
    """
    subdir = base_output_dir / tool_type.lower()
    subdir.mkdir(parents=True, exist_ok=True)

    base_name = input_path.stem
    expected_suffix = ".json" if tool_type != "Allergenicity" else ".txt"
    
    for file in subdir.glob(f"*{expected_suffix}"):
        if base_name in file.stem:
            return True
    return False

def run_predictions_parallel(job_list, output_dir, max_threads):
    print(f"\n⚙️ Running {len(job_list)} jobs in parallel using {max_threads} thread(s)...")
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [
            executor.submit(tool_runners[tool], tool_path, input_file, output_dir)
            for tool, tool_path, input_file in job_list
        ]
        for f in futures:
            f.result()

def main():
    parser = argparse.ArgumentParser(description="Run allergenicity, population coverage, and conservation analysis")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Path to directory containing tool executables")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["Allergenicity", "PopulationCoverage", "Conservation"],
                        default=["Allergenicity", "PopulationCoverage", "Conservation"],
                        help="Tools to run")
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_post_analysis_outputs"),
                        help="Output directory to save results")

    args = parser.parse_args()

    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    sequence_path = pathogen_path / args.sequence_dir
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validations
    if not pathogen_path.exists():
        print(f"❌ Pathogen directory does not exist: {pathogen_path}")
        sys.exit(1)
    if not sequence_path.exists():
        print(f"❌ Sequence directory does not exist: {sequence_path}")
        sys.exit(1)
    if not Path(args.tool_root).exists():
        print(f"❌ Tool root directory does not exist: {args.tool_root}")
        sys.exit(1)

    # Fetch fasta files
    fasta_files = common.get_fasta_files(pathogen_path, args.sequence_dir)
    if not fasta_files:
        print(f"❌ No FASTA files found in {sequence_path}")
        sys.exit(1)

    # Detect and map available tools
    selected_tools = set(args.tools)
    tool_map = {
        "Allergenicity": Path(args.tool_root) / "algpred2" / "algpred2.py",
        "PopulationCoverage": Path(args.tool_root) / "population_coverage" / "population_coverage.py",
        "Conservation": Path(args.tool_root) / "cluster" / "cluster.py"
    }

    missing_tools = [tool for tool in selected_tools if not tool_map[tool].exists()]
    if missing_tools:
        print(f"❌ Tool(s) not found in tool-root: {', '.join(missing_tools)}")
        sys.exit(1)

    final_tools = {tool: tool_map[tool] for tool in selected_tools}

    # Create output folders and prepare jobs
    _, output_dir = common.prepare_output_dirs(pathogen_path, output_dir, final_tools.keys())
    all_jobs = []

    for tool, path in final_tools.items():
        print(f"\n🧪 Preparing {tool} analysis")
        for fasta_file in fasta_files:
            if is_job_completed(tool, fasta_file, output_dir):
                print(f"⏩ Skipping {fasta_file.name} — {tool} already done.")
                continue
            all_jobs.append((tool, path, fasta_file))

    if not all_jobs:
        print("✅ All jobs are already completed.")
        sys.exit(0)

    # Summary
    print(f"\n📋 Job Summary:")
    job_counter = Counter(job[0] for job in all_jobs)
    for tool, count in job_counter.items():
        print(f"  - {tool}: {count} job(s)")

    run_predictions_parallel(all_jobs, output_dir, args.threads)

    print("\n✅ All post-analysis tasks completed successfully.")

if __name__ == "__main__":
    main()
