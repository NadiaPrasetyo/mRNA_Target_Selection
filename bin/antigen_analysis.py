# run signalp, targetp, and tmhmm tools to predict signal peptides, targeting peptides, and transmembrane helices within the antigen sequences.
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
import sys

# Map of tool runners and their expected script paths relative to tool-root
tool_runners = {
    "SIGNALP": "signalp_wrapper.sh",
    "TARGETP": "targetp/targetp",
    "TMHMM": "tmhmm/bin/tmhmm"
}

def run_tool(tool_name, tool_script_path, input_file, output_dir):
    output_file = output_dir / f"{input_file.stem}_{tool_name.lower()}.out"
    try:
        with open(output_file, "w") as outfile:
            subprocess.run([
                tool_script_path,
                str(input_file)
            ], stdout=outfile, stderr=subprocess.PIPE, check=True)
        print(f"✅ {tool_name} completed for {input_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ {tool_name} failed for {input_file.name}: {e.stderr.decode()}")

def run_parallel_jobs(jobs, threads):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(run_tool, *job) for job in jobs]
        for f in futures:
            f.result()

def main():
    parser = argparse.ArgumentParser(
        description="Run SignalP, TargetP, and TMHMM on input FASTA files",
        usage="run_predictors.py <pathogen_dir> <sequence_dir> --tool-root <tool_root> [options]"
    )
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing tool wrappers and executables")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["SIGNALP", "TARGETP", "TMHMM"],
                        default=["SIGNALP", "TARGETP", "TMHMM"],
                        help="Specify which tools to run (default: all)")

    args = parser.parse_args()

    base_path = Path("data") / args.pathogen_dir / args.sequence_dir
    if not base_path.exists():
        print(f"❌ Invalid input directory: {base_path}")
        sys.exit(1)

    fasta_files = list(base_path.glob("*.fasta")) + list(base_path.glob("*.fa"))
    if not fasta_files:
        print(f"❌ No FASTA files found in: {base_path}")
        sys.exit(1)

    tool_root = Path(args.tool_root)
    if not tool_root.exists():
        print(f"❌ Tool root directory does not exist: {tool_root}")
        sys.exit(1)

    output_dir = base_path / "tool_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for tool_name in args.tools:
        script_rel_path = tool_runners[tool_name]
        script_abs_path = tool_root / script_rel_path
        if not script_abs_path.exists():
            print(f"⚠️ Tool script not found: {script_abs_path}")
            continue

        for fasta in fasta_files:
            jobs.append((tool_name, str(script_abs_path), fasta, output_dir))

    if not jobs:
        print("❌ No jobs to run. Exiting.")
        sys.exit(1)

    print(f"\n🚀 Running {len(jobs)} jobs with {args.threads} threads...")
    run_parallel_jobs(jobs, args.threads)
    print("\n✅ All predictions complete.")

if __name__ == "__main__":
    main()
