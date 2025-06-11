# antigen_analysis.py
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
import sys

tool_runners = {
    "SIGNALP": "run_signalp.sh",
    "TARGETP": "run_targetp.sh",
    "TMHMM": "run_tmhmm.sh"
}

# Define required directories and files per tool inside the tool root directory
tool_requirements = {
    "SIGNALP": {
        "dirs": ["signalp/bin"],
        "files": ["signalp/bin/signalp"]
    },
    "TARGETP": {
        "dirs": ["targetp/bin"],
        "files": ["targetp/bin/targetp"]
    },
    "TMHMM": {
        "dirs": ["TMHMM2.0a/bin", "TMHMM2.0a/lib"],
        "files": [
            "TMHMM2.0a/bin/tmhmm",
            "TMHMM2.0a/bin/decodeanhmm",
            "TMHMM2.0a/bin/tmhmmformat.pl",
            "TMHMM2.0a/lib/TMHMM2.0.model",
            "TMHMM2.0a/lib/TMHMM2.0.options"
        ]
    }
}


def check_tool_environment(tool_name, tool_root):
    """Check existence of required dirs and files for a given tool."""
    reqs = tool_requirements.get(tool_name, {})
    dirs = reqs.get("dirs", [])
    files = reqs.get("files", [])

    missing_dirs = [d for d in dirs if not (tool_root / d).is_dir()]
    missing_files = [f for f in files if not (tool_root / f).is_file()]

    if missing_dirs:
        print(f"❌ Missing required directories for {tool_name}: {missing_dirs}")
    if missing_files:
        print(f"❌ Missing required files for {tool_name}: {missing_files}")

    return not (missing_dirs or missing_files)

def run_tool(tool_name, tool_script_path, input_file, output_dir, batch_size):
    output_file = output_dir / f"{input_file.stem}_{tool_name.lower()}.out"
    if output_file.exists():
        print(f"⏭️ Skipping {tool_name} for {input_file.name} (output already exists)")
        return

    try:
        with open(output_file, "w") as outfile:
            if tool_name in ("SIGNALP", "TARGETP"):
                cmd = [tool_script_path, str(input_file), str(output_dir), str(batch_size)]
            elif tool_name == "TMHMM":
                cmd = [tool_script_path, str(input_file), str(output_dir)]
            else:
                print(f"⚠️ Unknown tool: {tool_name}")
                return

            subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, check=True)
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
    parser.add_argument("--batch-size", type=int, default=10000,
                        help="Batch size to use for SignalP (default: 10000)")

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

    # Check tools environment & scripts
    valid_tools = []
    for tool_name in args.tools:
        # Check scripts exist
        script_rel_path = tool_runners.get(tool_name)
        script_abs_path = tool_root / script_rel_path
        if not script_abs_path.exists():
            print(f"❌ Tool script not found for {tool_name}: {script_abs_path}")
            continue

        # Check required dirs/files for tool
        if not check_tool_environment(tool_name, tool_root):
            print(f"❌ Environment check failed for tool {tool_name}, skipping.")
            continue

        valid_tools.append(tool_name)

    if not valid_tools:
        print("❌ No valid tools to run after environment checks. Exiting.")
        sys.exit(1)

    # Prepare output dirs per tool
    epitope_root = base_path / "epitope_outputs"
    epitope_root.mkdir(parents=True, exist_ok=True)

    jobs = []
    for tool_name in valid_tools:
        tool_out_dir = epitope_root / tool_name.lower()
        tool_out_dir.mkdir(parents=True, exist_ok=True)

        if not os.access(tool_out_dir, os.W_OK):
            print(f"❌ Output directory not writable: {tool_out_dir}")
            continue

        script_rel_path = tool_runners[tool_name]
        script_abs_path = tool_root / script_rel_path

        for fasta in fasta_files:
            # Determine expected output filename
            output_file = tool_out_dir / f"{fasta.stem}_{tool_name.lower()}.out"
            if output_file.exists():
                print(f"⏭️ Skipping {tool_name} for {fasta.name} (output already exists)")
                continue

            jobs.append((tool_name, str(script_abs_path), fasta, tool_out_dir, args.batch_size))

    if not jobs:
        print("❌ No jobs to run after checks. Exiting.")
        sys.exit(1)

    print(f"\n🚀 Running {len(jobs)} jobs with {args.threads} threads...")
    run_parallel_jobs(jobs, args.threads)
    print("\n✅ All predictions complete.")

if __name__ == "__main__":
    import os
    main()
