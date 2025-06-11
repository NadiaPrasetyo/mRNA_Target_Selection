# antigen_analysis.py
import argparse
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_signalp, run_targetp, run_tmhmm, common

TOOL_RUNNERS = {
    "SIGNALP": run_signalp.run,
    "TARGETP": run_targetp.run,
    "TMHMM": run_tmhmm.run
}

VALID_TOOLS = list(TOOL_RUNNERS.keys())


def run_tool(tool_name: str, runner_func, input_file: Path, output_dir: Path, batch_size: int, tool_root: Path) -> None:
    output_file = output_dir / f"{input_file.stem}_{tool_name.lower()}.out"
    if output_file.exists():
        print(f"⏭️ Skipping {tool_name} for {input_file.name} (output already exists)")
        return

    try:
        if tool_name in ("SIGNALP", "TARGETP"):
            runner_func(input_file, output_dir, batch_size)
        elif tool_name == "TMHMM":
            runner_func(tool_root / "tmhmm-2.0a", input_file, output_dir)
        else:
            print(f"⚠️ Unknown tool: {tool_name}")
            return
        print(f"✅ {tool_name} completed for {input_file.name}")
    except Exception as e:
        print(f"❌ {tool_name} failed for {input_file.name}: {e}")

def run_parallel_jobs(jobs, threads: int) -> None:
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
    parser = argparse.ArgumentParser(
        description="Run SignalP, TargetP, and TMHMM on input FASTA files",
        usage="run_predictors.py <pathogen_dir> <sequence_dir> --tool-root <tool_root> [options]"
    )
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing tool wrappers and executables")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=VALID_TOOLS, default=VALID_TOOLS,
                        help="Specify which tools to run (default: all)")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size for SignalP/TargetP (default: 10000)")

    args = parser.parse_args()

    base_path = Path("data") / args.pathogen_dir / args.sequence_dir
    if not base_path.exists():
        print(f"❌ Invalid input directory: {base_path}")
        sys.exit(1)

    fasta_files = common.get_fasta_files(Path("data") / args.pathogen_dir, args.sequence_dir)
    if not fasta_files:
        print("❌ No FASTA files found.")
        sys.exit(1)


    tool_root = Path(args.tool_root)
    if not tool_root.exists():
        print(f"❌ Tool root directory does not exist: {tool_root}")
        sys.exit(1)

    found = common.check_signalp_targetp_tmhmm(tool_root)
    if not found:
        print(f"❌ Required components missing under {tool_root}. Exiting.")
        sys.exit(1)



    epitope_root = base_path / "epitope_outputs"
    epitope_root.mkdir(parents=True, exist_ok=True)

    epitope_root = Path("data") / args.pathogen_dir / args.sequence_dir / "epitope_outputs"

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

            jobs.append((tool_name, runner_func, fasta, tool_out_dir, args.batch_size, tool_root))


    if not jobs:
        print("❌ No jobs to run after validation. Exiting.")
        sys.exit(1)

    print(f"\n🚀 Running {len(jobs)} jobs using {args.threads} threads...")
    run_parallel_jobs(jobs, args.threads)
    print("\n✅ All predictions complete.")


if __name__ == "__main__":
    main()
