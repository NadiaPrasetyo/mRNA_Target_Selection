import argparse
import sys
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from tools import run_algpred, run_popcoverage, run_cluster, common

tool_runners = {
    "Allergenicity": run_algpred.run,
    "PopCoverage": run_popcoverage.run,
    "Cluster": run_cluster.run,
}

def is_job_completed(tool_type, input_path, base_output_dir):
    subdir = base_output_dir / tool_type.lower()
    subdir.mkdir(parents=True, exist_ok=True)
    base_name = input_path.stem
    expected_suffix = ".json" if tool_type == "Cluster" else ".txt"
    for file in subdir.glob(f"*{expected_suffix}"):
        if base_name in file.stem:
            return True
    return False

def run_predictions_parallel(job_list, output_dir, max_threads):
    print(f"\n⚙️ Starting parallel execution of {len(job_list)} job(s) using {max_threads} thread(s)...")
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = []
        for tool_type, tool_path, input_file in job_list:
            sub_output_dir = output_dir / tool_type.lower()

            if tool_type == "Cluster" :
                temp_json_dir = output_dir / "json_inputs"
                temp_json_dir.mkdir(parents=True, exist_ok=True)
                json_paths = common.parse_fasta_to_jsons(
                    input_file, temp_json_dir,
                    alleles=[], peptide_lengths=[],
                    tool_type="cluster", strain_name=input_file.stem[:-len("_matched_antigen")]
                )
                for json_path in json_paths:
                    futures.append(
                        executor.submit(tool_runners[tool_type], tool_path, json_path, sub_output_dir)
                    )

            elif tool_type == "PopCoverage":
                temp_txt_dir = output_dir / "popcov_inputs"
                temp_txt_dir.mkdir(parents=True, exist_ok=True)
                txt_paths = common.parse_fasta_to_jsons(
                    input_file, temp_txt_dir,
                    alleles=[], peptide_lengths=[],
                    tool_type="popcoverage", strain_name=input_file.stem[:-len("_matched_antigen")]
                )
                for txt_path in txt_paths:
                    futures.append(
                        executor.submit(tool_runners[tool_type], tool_path, txt_path, sub_output_dir)
                    )

            else:
                futures.append(
                    executor.submit(tool_runners[tool_type], tool_path, input_file, sub_output_dir)
                )

        for f in futures:
            f.result()

def main():
    parser = argparse.ArgumentParser(description="Run evaluation tools: Allergenicity, Population Coverage, Cluster")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing analysis tools")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["Allergenicity", "PopCoverage", "Cluster"], default=None,
                        help="Specify which tools to run (default: all detected tools)")
    parser.add_argument("--output-dir", type=Path, default=Path("evaluation_outputs"),
                        help="Directory to save output files (default: 'evaluation_outputs')")
    args = parser.parse_args()

    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    sequence_path = pathogen_path / args.sequence_dir
    output_dir = args.output_dir

    for p in [pathogen_path, sequence_path, args.tool_root]:
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

    fasta_files = common.get_fasta_files(pathogen_path, args.sequence_dir)
    if not fasta_files:
        print(f"❌ No FASTA files found in {sequence_path}")
        sys.exit(1)

    _, output_dir = common.prepare_output_dirs(pathogen_path, output_dir, final_tools.keys())

    jobs = []
    for tool_type, tool_path in final_tools.items():
        print(f"\n🧪 Preparing {tool_type} analysis...")
        for input_file in fasta_files:
            if is_job_completed(tool_type, input_file, output_dir):
                print(f"⏩ Skipping {input_file.name} — already processed by {tool_type}")
                continue
            jobs.append((tool_type, tool_path, input_file))

    if not jobs:
        print("❌ No jobs to run. Exiting.")
        # clean up any temporary directories created
        common.cleanup_temp_dirs(temp_dir=output_dir / "json_inputs")
        common.cleanup_temp_dirs(temp_dir=output_dir / "popcov_inputs")
        sys.exit(1)

    print(f"\n🚀 Running {len(jobs)} evaluations with {args.threads} thread(s)...")
    job_counter = Counter(j[0] for j in jobs)
    print("\n📋 Job Summary:")
    for k, v in job_counter.items():
        print(f"  - {k}: {v} job(s)")

    run_predictions_parallel(jobs, output_dir, args.threads)
    print("\n✅ Evaluation complete.")

    # clean up any temporary directories created
    common.cleanup_temp_dirs(temp_dir=output_dir / "json_inputs")
    common.cleanup_temp_dirs(temp_dir=output_dir / "popcov_inputs")

if __name__ == "__main__":
    main()
