"""
IEDB_epitope.py
Command-line tool to run IEDB epitope prediction tools (MHCI, MHCII, BCell) on input FASTA files for immunoinformatics analysis.
Overview:
    - Scans a specified pathogen sequence directory for FASTA files.
    - Runs selected IEDB prediction tools (MHCI, MHCII, BCell) on each input file.
    - Supports custom allele panels and peptide length ranges for MHCI and MHCII.
    - Converts FASTA to TXT for BCell predictions as needed.
    - Supports parallel execution of jobs for efficient processing.
    - Organizes results into structured output directories and manages temporary files.
Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    sequence_dir (str): Subdirectory under `pathogen_dir` containing FASTA files.
    --tool-root (str, required): Root directory containing IEDB tool wrappers and executables.
    --threads (int, optional): Number of parallel threads to use (default: 4).
    --mhci-peptide-lengths (int, int, optional): Min and max peptide lengths for MHCI (default: 8 11).
    --mhcii-peptide-lengths (int, int, optional): Min and max peptide lengths for MHCII (default: 11 25).
    --tools (list, optional): List of tools to run (choices: MHCI, MHCII, BCell; default: all detected).
    --mhci-allele-panel (str, optional): Allele panel for MHCI (choices: default, extended, custom; default: default).
    --mhci-custom-alleles (list, optional): Custom alleles for MHCI if panel is 'custom'.
    --mhcii-allele-panel (str, optional): Allele panel for MHCII (choices: default, extended, custom; default: default).
    --mhcii-custom-alleles (list, optional): Custom alleles for MHCII if panel is 'custom'.
    --output-dir (str, optional): Output directory for results (default: epitope_outputs).
Requirements:
    - IEDB tool wrappers and executables for MHCI, MHCII, and BCell available under `tool-root`.
    - Input FASTA files present in the specified sequence directory.
    - Python packages: argparse, pathlib, concurrent.futures, collections.
Usage Example:
    python IEDB_epitope.py influenza sequences --tool-root /opt/iedb_tools --threads 8 --mhci-peptide-lengths 8 11 --tools MHCI MHCII
Outputs:
    data/<pathogen_dir>/<output_dir>/<tool>/<input_file>_<TOOL>.json   # Prediction results for MHCI and MHCII
    data/<pathogen_dir>/<output_dir>/bcell/<input_file>.txt            # Prediction results for BCell
Author: Nadia
"""
import argparse
from collections import Counter
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tools import run_mhci, run_mhcii, run_bcell, run_ellipro, common
import sys
import logging

# Mapping of tool types to their respective runner functions
tool_runners = {
    "MHCI": run_mhci.run,
    "MHCII": run_mhcii.run,
    "BCell": run_bcell.run,
    "Ellipro": run_ellipro.run,  
}

def is_job_completed(tool_type, input_path, base_output_dir):
    """
    Check if a job has already been processed by verifying if a result
    JSON file exists in the tool-specific output directory.
    Looks for any JSON file that contains the input's base name
    and ends with the tool suffix.
    Args:
        tool_type (str): Type of tool (e.g., "MHCI", "MHCII", "BCell").
        input_path (Path): Path to the input file.
        base_output_dir (Path): Base output directory where results are stored.
    Returns:
        bool: True if the job has been completed (result file exists), False otherwise.
    """
    subdir = base_output_dir / tool_type.lower()
    if not subdir.exists():
        logging.warning(f"❌ Output directory for {tool_type} does not exist: {subdir}")
        return False

    base_name = input_path.stem
    expected_suffix = f"{tool_type.upper()}.json" # e.g., "MHCI.json", "MHCII.json", "BCELL.txt", "Ellipro.txt"
    # Check for any file that matches the base name and expected suffix
    
    # Check for any file that matches the base name and expected suffix
    if tool_type in ["BCell", "Ellipro"]:
        expected_suffix = ".txt"

    for file in subdir.glob(f"*{expected_suffix}"):
        if base_name in file.stem:
            return True
    return False

def run_predictions_parallel(job_list, output_dir, max_threads):
    """
    Run the prediction jobs in parallel using a thread pool executor.
    Args:
        job_list (list): List of tuples containing (tool_type, tool_path, input_file).
        output_dir (Path): Base output directory for results.
        max_threads (int): Maximum number of threads to use for parallel execution.
    """
    logging.info(f"\n⚙️ Starting parallel execution of {len(job_list)} job(s) using {max_threads} thread(s)...")

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = [
            executor.submit(tool_runners[tool_type], json_file, tool_path, output_dir)
            for tool_type, tool_path, json_file in job_list
        ]
        for f in futures:
            f.result()


def main():
    parser = argparse.ArgumentParser(description="Run epitope predictions (MHCI, MHCII, BCell, Ellipro)")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence subdirectory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing IEDB tools")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--mhci-peptide-lengths", "-mhci-pl", nargs=2, type=int, default=[8, 11])
    parser.add_argument("--mhcii-peptide-lengths", "-mhcii-pl", nargs=2, type=int, default=[11, 25])
    parser.add_argument("--tools", nargs="+", choices=["MHCI", "MHCII", "BCell", "Ellipro"], default=None)
    parser.add_argument("--mhci-allele-panel", choices=["default", "extended", "custom"], default="default")
    parser.add_argument("--mhci-custom-alleles", nargs="+", default=None)
    parser.add_argument("--mhcii-allele-panel", choices=["default", "extended", "custom"], default="default")
    parser.add_argument("--mhcii-custom-alleles", nargs="+", default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("epitope_outputs"))
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    data_dir = Path("data")
    pathogen_path = data_dir / args.pathogen_dir
    sequence_path = pathogen_path / args.sequence_dir
    output_dir = args.output_dir
    tool_root = Path(args.tool_root).resolve()

    if args.verbose:
        log_file = pathogen_path / output_dir / "antigen_analysis.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers = [
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file, mode='a')
        ]
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=handlers,
            force=True
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s: %(message)s",
            force=True
        )

    if not pathogen_path.exists() or not sequence_path.exists() or not tool_root.exists():
        logging.error("❌ One or more required directories do not exist.")
        sys.exit(1)

    tool_map = common.check_iedb_tool(tool_root)
    if not tool_map:
        logging.warning(f"❌ No valid IEDB tools found in: {tool_root}")
        sys.exit(1)

    # Default tools: all except Ellipro
    default_tools = ["MHCI", "MHCII", "BCell"]
    if args.tools is None:
        selected_tools = set([t for t in default_tools if t in tool_map])
    else:
        selected_tools = set(args.tools)

    # Ellipro can only be run by itself
    if "Ellipro" in selected_tools and len(selected_tools) > 1:
        logging.error("❌ Ellipro can only be run by itself due to input type differences. Please select only Ellipro or other tools.")
        sys.exit(1)

    missing_tools = selected_tools - set(tool_map.keys())
    if missing_tools:
        logging.warning(f"⚠️ Warning: Tools requested but not found: {', '.join(missing_tools)}")

    final_tools = {t: tool_map[t] for t in selected_tools if t in tool_map}
    if not final_tools:
        logging.warning("❌ No tools available to run after filtering.")
        sys.exit(1)

    # ONLY Get PDB files if Ellipro is selected
    pdb_files = []
    if "Ellipro" in final_tools:
        pdb_files = common.get_pdb_files(pathogen_path, args.sequence_dir)
        logging.info(f"📂 Found {len(pdb_files)} structure file(s): {[p.name for p in pdb_files]}")
    else:
        fasta_files = common.get_fasta_files(pathogen_path, args.sequence_dir)

    # Convert fasta to txt if BCell is selected
    txt_files = []
    if "BCell" in final_tools:
        temp_txt_dir = pathogen_path / "temp_txt"
        txt_files = common.convert_fasta_to_txt(fasta_files, temp_txt_dir)

    # Prepare temporary JSON directory for MHCI and MHCII and output directories for all tools
    temp_json_dir, output_dir = common.prepare_output_dirs(pathogen_path, output_dir, final_tools.keys(), temp=True)

    all_jobs = []
    for tool_type, tool_path in final_tools.items():
        logging.info(f"\n🧪 Preparing {tool_type} predictions")

        if tool_type == "MHCI":
            alleles = common.get_alleles(tool_type, args.mhci_allele_panel, args.mhci_custom_alleles)
            peptide_lengths = args.mhci_peptide_lengths
        elif tool_type == "MHCII":
            alleles = common.get_alleles(tool_type, args.mhcii_allele_panel, args.mhcii_custom_alleles)
            peptide_lengths = args.mhcii_peptide_lengths
        else:
            alleles = []
            peptide_lengths = None

        if tool_type == "BCell":
            input_files_tool = txt_files
        elif tool_type == "Ellipro":
            input_files_tool = pdb_files
        else:
            input_files_tool = fasta_files

        for input_file in input_files_tool:
            logging.info(f"🧬 Processing {input_file.name}")

            if is_job_completed(tool_type, input_file, output_dir):
                logging.info(f"⏩ Skipping {input_file.name} — {tool_type} result already exists.")
                continue

            if tool_type in ["MHCI", "MHCII"]:
                json_paths = common.parse_fasta_to_jsons(input_file, temp_json_dir, alleles, peptide_lengths, tool_type, args.sequence_dir)
                for jp in json_paths:
                    if is_job_completed(tool_type, jp, output_dir):
                        logging.info(f"⏩ Skipping {jp.name} — already processed.")
                        continue
                    all_jobs.append((tool_type, tool_path, jp))
            else:
                all_jobs.append((tool_type, tool_path, input_file))

    if not all_jobs:
        logging.warning("❌ No jobs to run.")
        common.cleanup_temp(temp_json_dir)
        if "BCell" in final_tools:
            common.cleanup_temp(temp_txt_dir)
        sys.exit(1)

    logging.info(f"\n🚀 Running predictions with {args.threads} threads...")

    job_counter = Counter([job[0] for job in all_jobs])
    logging.info("\n📋 Job Summary:")
    for tool, count in job_counter.items():
        logging.info(f"  - {tool}: {count} job(s)")

    run_predictions_parallel(all_jobs, output_dir, args.threads)

    common.cleanup_temp(temp_json_dir)
    if "BCell" in final_tools:
        common.cleanup_temp(temp_txt_dir)

    logging.info("\n✅ Prediction complete.")

if __name__ == "__main__":
    main()
