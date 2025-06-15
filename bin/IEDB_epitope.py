import argparse
from collections import Counter
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_mhci, run_mhcii, run_bcell, common
import sys

tool_runners = {
    "MHCI": run_mhci.run,
    "MHCII": run_mhcii.run,
    "BCell": run_bcell.run
}

def is_job_completed(tool_type, input_path, base_output_dir):
    """
    Check if a job has already been processed by verifying if a result
    JSON file exists in the tool-specific output directory.
    Looks for any JSON file that contains the input's base name
    and ends with the tool suffix.
    """
    subdir = base_output_dir / tool_type.lower()
    subdir.mkdir(parents=True, exist_ok=True)
    print(f"🔍 Checking if {tool_type} job for {input_path.name} is completed in {base_output_dir}")

    base_name = input_path.stem
    expected_suffix = f"{tool_type.upper()}.json"

    for file in subdir.glob(f"*{expected_suffix}"):
        if base_name in file.stem:
            return True
    return False

def run_predictions_parallel(job_list, output_dir, max_threads):
    print(f"\n⚙️ Starting parallel execution of {len(job_list)} job(s) using {max_threads} thread(s)...")

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [
            executor.submit(tool_runners[tool_type], json_file, tool_path, output_dir)
            for tool_type, tool_path, json_file in job_list
        ]
        for f in futures:
            f.result()

def main():
    parser = argparse.ArgumentParser(description="Run epitope predictions (MHCI, MHCII, BCell)",
                                     usage="iedb_epitope.py <pathogen_dir> <sequence_dir> --tool-root <tool_root> [options]")
    parser.epilog = (
        f"Example usage:\n"
        f"  iedb_epitope.py influenza sequences --tool-root /path/to/iedb/tools --threads 8 --peptide-lengths 8 11 --tools MHCI MHCII\n\n"
        "This script will run epitope predictions for the specified pathogen and sequence directories using the IEDB tools MHCI and MHCII with peptide length between 8-11 and with 8 parallel threads."
    )
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing IEDB tools")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--mhci-peptide-lengths", "-mhci-pl", nargs=2, type=int, metavar=('MIN', 'MAX'),
                        default=[8, 11], help="Min and max peptide lengths for MHCI (default 8-11)")
    parser.add_argument("--mhcii-peptide-lengths", "-mhcii-pl", nargs=2, type=int, metavar=('MIN', 'MAX'),
                        default=[11, 25], help="Min and max peptide lengths for MHCII (default 11-25)")
    parser.add_argument("--tools", nargs="+", choices=["MHCI", "MHCII", "BCell"], default=None,
                        help="Specify which tools to run (default: all detected tools)")
    parser.add_argument("--mhci-allele-panel", choices=["default", "extended", "custom"], default="default",
                        help="Allele panel for MHCI. Choose 'default', 'extended', or 'custom'")
    parser.add_argument("--mhci-custom-alleles", nargs="+", default=None,
                        help="Custom alleles list for MHCI if panel is 'custom'")
    parser.add_argument("--mhcii-allele-panel", choices=["default", "extended", "custom"], default="default",
                        help="Allele panel for MHCII. Choose 'default', 'extended', or 'custom'")
    parser.add_argument("--mhcii-custom-alleles", nargs="+", default=None,
                        help="Custom alleles list for MHCII if panel is 'custom'")
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"),
                        help="Directory to save output files (default: 'epitope_outputs')")

    args = parser.parse_args()

    data_dir = Path("data")
    
    pathogen_path = data_dir / args.pathogen_dir
    if not pathogen_path.exists() or not pathogen_path.is_dir():
        print(f"❌ Pathogen directory does not exist or is not a directory: {pathogen_path}")
        sys.exit(1)

    # check that output directory exists or create it
    output_dir = output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)    

    sequence_path = pathogen_path / args.sequence_dir
    if not sequence_path.exists() or not sequence_path.is_dir():
        print(f"❌ Sequence directory does not exist or is not a directory: {sequence_path}")
        sys.exit(1)

    tool_root = Path(args.tool_root)
    if not tool_root.exists() or not tool_root.is_dir():
        print(f"❌ Tool root directory does not exist or is not a directory: {tool_root}")
        sys.exit(1)

    tool_map = common.check_iedb_tool(tool_root)
    if not tool_map:
        print(f"❌ No valid IEDB tools found in: {tool_root}")
        sys.exit(1)

    selected_tools = set(args.tools) if args.tools else set(tool_map.keys())
    missing_tools = selected_tools - set(tool_map.keys())
    if missing_tools:
        print(f"⚠️ Warning: Tools requested but not found: {', '.join(missing_tools)}")

    final_tools = {t: tool_map[t] for t in selected_tools if t in tool_map}
    if not final_tools:
        print("❌ No tools available to run after filtering.")
        sys.exit(1)

    fasta_files = common.get_fasta_files(pathogen_path, args.sequence_dir)
    if not fasta_files:
        print(f"❌ No FASTA files found in {sequence_path}")
        sys.exit(1)

    # Convert fasta to txt if BCell is selected
    txt_files = []
    if "BCell" in final_tools:
        temp_txt_dir = pathogen_path / "temp_txt"
        txt_files = common.convert_fasta_to_txt(fasta_files, temp_txt_dir)

    # Prepare temporary JSON directory for MHCI and MHCII and output directories for all tools
    temp_json_dir, output_dir = common.prepare_output_dirs(pathogen_path, output_dir, final_tools.keys())

    all_jobs = []
    for tool_type, tool_path in final_tools.items():
        print(f"\n🧪 Preparing {tool_type} predictions")

        if tool_type == "MHCI":
            alleles = common.get_alleles(tool_type, args.mhci_allele_panel, args.mhci_custom_alleles)
            peptide_lengths = args.mhci_peptide_lengths
        elif tool_type == "MHCII":
            alleles = common.get_alleles(tool_type, args.mhcii_allele_panel, args.mhcii_custom_alleles)
            peptide_lengths = args.mhcii_peptide_lengths
        else:
            alleles = []
            peptide_lengths = None

        input_files = txt_files if tool_type == "BCell" else fasta_files
        for input_file in input_files:
            print(f"🧬 Processing {input_file.name}")

            if is_job_completed(tool_type, input_file, output_dir.parent):
                print(f"⏩ Skipping {input_file.name} — {tool_type} result already exists.")
                continue
            
            if tool_type == "BCell":
                all_jobs.append((tool_type, tool_path, input_file))
            else:
                json_paths = common.parse_fasta_to_jsons(
                    input_file,
                    temp_json_dir,
                    alleles,
                    peptide_lengths,
                    tool_type,
                    args.sequence_dir
                )
                all_jobs.extend([(tool_type, tool_path, jp) for jp in json_paths])

    if not all_jobs:
        print("❌ No jobs to run. Exiting.")
        sys.exit(1)

    print(f"\n🚀 Running predictions with {args.threads} threads...")

    job_counter = Counter([job[0] for job in all_jobs])
    print("\n📋 Job Summary:")
    for tool, count in job_counter.items():
        print(f"  - {tool}: {count} job(s)")

    run_predictions_parallel(all_jobs, output_dir, args.threads)

    common.cleanup_temp(temp_json_dir)
    print("\n✅ Prediction complete.")

if __name__ == "__main__":
    main()
