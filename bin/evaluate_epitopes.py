import argparse
import sys
import logging
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from tools import run_algpred, run_popcoverage, run_cluster, common, extract_epitopes

tool_runners = {
    "Allergenicity": run_algpred.run,
    "PopCoverage": run_popcoverage.run,
    "Cluster": run_cluster.run,
}

def is_job_completed(tool_type, input_path, base_output_dir):
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
                mhci_ep = list((args.epitope_dir / "mhci").glob("*.txt"))
                mhcii_ep = list((args.epitope_dir / "mhcii").glob("*.txt"))
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
    parser = argparse.ArgumentParser(description="Run evaluation tools: Allergenicity, Population Coverage, Cluster")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing analysis tools")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--tools", nargs="+", choices=["Allergenicity", "PopCoverage", "Cluster"], default=None,
                        help="Specify which tools to run (default: all detected tools)")
    parser.add_argument("--output-dir", type=Path, default=Path("evaluation_outputs"),
                        help="Directory to save output files (default: 'evaluation_outputs')")
    parser.add_argument("--epitope-dir", type=Path, required=False,
                        help="Directory containing epitope predictions with mhci/, mhcii/, and bcell/ subdirs")

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    logging.info("Starting epitope evaluation pipeline")

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
    for t, p in final_tools.items():
        for f in fasta_files:
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
