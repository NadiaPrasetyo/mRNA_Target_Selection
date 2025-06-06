import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_mhci, run_mhcii, run_bcell, common
import sys

tool_runners = {
    "MHCI": run_mhci.run,
    "MHCII": run_mhcii.run,
    "BCell": run_bcell.run
}

def run_predictions_parallel(job_list, output_dir, max_threads):
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
    example_with_options = ("iedb_epitope.py influenza sequences --tool-root /path/to/iedb/tools --threads 8 --peptide-lengths 8 11 --tools MHCI MHCII")
    parser.epilog = (
        f"Example usage:\n"
        f"  {example_with_options}\n\n"
        "This script will run epitope predictions for the specified pathogen and sequence directories using the IEDB tools MHCI and MHCII with peptide length between 8-11 and with 8 parallel threads."
    )
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing IEDB tools")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--peptide-lengths", "-pl", nargs=2, type=int, metavar=('MIN', 'MAX'),
                    default=[9, 11],
                    help="Minimum and maximum peptide lengths to consider")

    parser.add_argument("--tools", nargs="+", choices=["MHCI", "MHCII", "BCell"], default=None,
                        help="Specify which tools to run (default: all detected tools)")

    # Allele panel and custom alleles for MHCI and MHCII
    parser.add_argument("--mhci-allele-panel", choices=["default", "extended", "custom"], default="default",
                        help="Allele panel for MHCI")
    parser.add_argument("--mhci-custom-alleles", nargs="+", default=None,
                        help="Custom alleles list for MHCI if panel is 'custom'")
    parser.add_argument("--mhcii-allele-panel", choices=["default", "extended", "custom"], default="default",
                        help="Allele panel for MHCII")
    parser.add_argument("--mhcii-custom-alleles", nargs="+", default=None,
                        help="Custom alleles list for MHCII if panel is 'custom'")

    args = parser.parse_args()

    # Check base directories exist
    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    if not pathogen_path.exists() or not pathogen_path.is_dir():
        print(f"❌ Pathogen directory does not exist or is not a directory: {pathogen_path}")
        sys.exit(1)

    sequence_path = pathogen_path / args.sequence_dir
    if not sequence_path.exists() or not sequence_path.is_dir():
        print(f"❌ Sequence directory does not exist or is not a directory: {sequence_path}")
        sys.exit(1)

    tool_root = Path(args.tool_root)
    if not tool_root.exists() or not tool_root.is_dir():
        print(f"❌ Tool root directory does not exist or is not a directory: {tool_root}")
        sys.exit(1)

    # Check available tools
    tool_map = common.check_iedb_tool(tool_root)
    if not tool_map:
        print(f"❌ No valid IEDB tools found in: {tool_root}")
        sys.exit(1)

    # If user specified a subset of tools, filter
    selected_tools = set(args.tools) if args.tools else set(tool_map.keys())
    missing_tools = selected_tools - set(tool_map.keys())
    if missing_tools:
        print(f"⚠️ Warning: Tools requested but not found: {', '.join(missing_tools)}")
    # Only keep detected + requested
    final_tools = {t: tool_map[t] for t in selected_tools if t in tool_map}
    if not final_tools:
        print("❌ No tools available to run after filtering.")
        sys.exit(1)

    # Find FASTA files
    fasta_files = common.get_fasta_files(pathogen_path, args.sequence_dir)
    if not fasta_files:
        print(f"❌ No FASTA files found in {sequence_path}")
        sys.exit(1)

    # Prepare output directories
    temp_json_dir, output_dir = common.prepare_output_dirs(pathogen_path)

    peptide_lengths = args.peptide_lengths  # Already [min, max]

    all_jobs = []
    for tool_type, tool_path in final_tools.items():
        print(f"\n🧪 Preparing {tool_type} predictions")

        if tool_type == "MHCI":
            alleles = common.get_alleles(tool_type, args.mhci_allele_panel, args.mhci_custom_alleles)
        elif tool_type == "MHCII":
            alleles = common.get_alleles(tool_type, args.mhcii_allele_panel, args.mhcii_custom_alleles)
        else:  # BCell or others
            alleles = []

        for fasta_file in fasta_files:
            print(f"🧬 Processing {fasta_file.name}")
            json_paths = common.parse_fasta_to_jsons(
                fasta_file,
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
    run_predictions_parallel(all_jobs, output_dir, args.threads)
    common.cleanup_temp(temp_json_dir)
    print("\n✅ Prediction complete.")

if __name__ == "__main__":
    main()
