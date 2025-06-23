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
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_signalp, run_targetp, common

# Define the mapping of tool names to their runner functions
TOOL_RUNNERS = {
    "SIGNALP": run_signalp.run,
    "TARGETP": run_targetp.run
}

# List of valid tools that can be run
VALID_TOOLS = list(TOOL_RUNNERS.keys())


def run_tool(tool_name: str, runner_func, input_file: Path, output_dir: Path, batch_size: int, tool_path: Path) -> None:
    """
    Run a specific tool on the input file and save the output to the specified directory.
    Args:
        tool_name (str): Name of the tool to run (e.g., SIGNALP, TARGETP).
        runner_func (function): Function to run the tool.
        input_file (Path): Path to the input FASTA file.
        output_dir (Path): Directory to save the output files.
        batch_size (int): Batch size for tools that support batching (e.g., SignalP, TargetP).
        tool_path (Path): Path to the tool executable.
    """
    output_file = output_dir / f"{input_file.stem}_{tool_name.lower()}.out"
    if output_file.exists():
        print(f"⏭️ Skipping {tool_name} for {input_file.name} (output already exists)")
        return

    try:
        runner_func(tool_path, input_file, output_dir, batch_size)
        print(f"✅ {tool_name} completed for {input_file.name}")
    except Exception as e:
        print(f"❌ {tool_name} failed for {input_file.name}: {e}")


def run_parallel_jobs(jobs, threads: int) -> None:
    """
    Run a list of jobs in parallel using a thread pool.
    Args:
        jobs (list): List of tuples containing job parameters (tool_name, runner_func, input_file, output_dir, batch_size, tool_path).
        threads (int): Number of threads to use for parallel execution.
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(run_tool, *job) for job in jobs]
        for f in futures:
            try:
                f.result()
            except TypeError as e:
                print(f"❌ Job failed due to argument mismatch: {e}")
            except Exception as e:
                print(f"❌ Job failed with unexpected error: {e}")


def main():
    """
    Main function to parse arguments and run the antigen analysis pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Run SignalP and TargetP on input FASTA files",
        usage="run_predictors.py <pathogen_dir> <sequence_dir> --tool-root <tool_root> [options]"
    )
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing tool wrappers and executables")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=VALID_TOOLS, default=VALID_TOOLS,
                        help="Specify which tools to run (default: both)")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size for SignalP/TargetP (default: 10000)")
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"),
                        help="Base output directory for results (default: epitope_outputs)")

    args = parser.parse_args()

    data_path = Path("data") / args.pathogen_dir
    sequence_path = data_path / args.sequence_dir
    if not sequence_path.exists():
        print(f"❌ Invalid input directory: {sequence_path}")
        sys.exit(1)

    fasta_files = common.get_fasta_files(data_path, args.sequence_dir)
    if not fasta_files:
        print("❌ No FASTA files found.")
        sys.exit(1)

    tool_root = Path(args.tool_root)
    if not tool_root.exists():
        print(f"❌ Tool root directory does not exist: {tool_root}")
        sys.exit(1)

    try:
        tool_paths = common.check_signalp_targetp(tool_root)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    epitope_root = data_path / args.output_dir
    epitope_root.mkdir(parents=True, exist_ok=True)

    jobs = []
    for tool_name in args.tools:
        tool_out_dir = epitope_root / tool_name.lower()

        if not common.ensure_writable_dir(tool_out_dir):
            print(f"❌ Skipping {tool_name} due to output directory issue.")
            continue

        runner_func = TOOL_RUNNERS[tool_name]

        for fasta in fasta_files:
            output_file = tool_out_dir / f"{fasta.stem}_{tool_name.lower()}.out"
            if output_file.exists():
                print(f"⏭️ Skipping {tool_name} for {fasta.name} (output already exists)")
                continue

            tool_path = tool_paths.get(tool_name)
            jobs.append((tool_name, runner_func, fasta, tool_out_dir, args.batch_size, tool_path))

    if not jobs:
        print("❌ No jobs to run after validation. Exiting.")
        sys.exit(1)

    print(f"\n🚀 Running {len(jobs)} jobs using {args.threads} threads...")
    run_parallel_jobs(jobs, args.threads)
    print("\n✅ All predictions complete.")


if __name__ == "__main__":
    main()
