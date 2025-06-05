import argparse
import os
import subprocess
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import tempfile
import getpass

def check_iedb_tool(iedb_dir):
    tool_path = Path(iedb_dir) / "src" / "tcell_mhci.py"
    if not tool_path.exists():
        raise FileNotFoundError(f"IEDB tool not found at {tool_path}")
    return str(tool_path)

def find_fasta_files(pathogen_dir, sequence_dir):
    search_path = Path("data") / pathogen_dir / sequence_dir
    if not search_path.exists():
        raise FileNotFoundError(f"Sequence directory does not exist: {search_path}")
    return list(search_path.rglob("*.fasta"))

def generate_input_json(fasta_file, allele_list, peptide_lengths):
    return {
        "input_sequence_text_file_path": str(fasta_file),
        "peptide_length_range": peptide_lengths,
        "alleles": ",".join(allele_list),
        "predictors": [
            {
                "type": "binding",
                "method": "netmhcpan_ba"
            }
        ]
    }

def run_prediction(fasta_file, tool_path, output_dir, allele_list, peptide_lengths):
    basename = fasta_file.stem
    json_input = generate_input_json(fasta_file, allele_list, peptide_lengths)

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".json") as temp_json:
        json.dump(json_input, temp_json)
        temp_json.flush()
        temp_json_path = temp_json.name

    output_prefix = Path(output_dir) / basename
    cmd = [
        "python3", tool_path,
        "-j", temp_json_path,
        "-o", str(output_prefix),
        "-f", "json"
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Processed: {fasta_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed on {fasta_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run IEDB Class I T Cell predictions on FASTA files.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence directory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    args = parser.parse_args()

    # Ask user for IEDB tool path
    iedb_dir = input("Enter full path to IEDB tool folder (<50 chars): ").strip()
    if len(iedb_dir) > 50:
        raise ValueError("IEDB path must be under 50 characters.")

    tool_path = check_iedb_tool(iedb_dir)
    fasta_files = find_fasta_files(args.pathogen_dir, args.sequence_dir)

    if not fasta_files:
        print("No FASTA files found.")
        return

    output_dir = Path("outputs") / args.pathogen_dir / args.sequence_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Customize alleles and peptide lengths if needed
    allele_list = ["HLA-A*02:01", "HLA-A*01:01"]
    peptide_lengths = [8, 11]

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                run_prediction, fasta, tool_path, output_dir, allele_list, peptide_lengths
            ) for fasta in fasta_files
        ]

        for future in futures:
            future.result()

if __name__ == "__main__":
    main()
