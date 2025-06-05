import argparse
import os
import subprocess
import json
from pathlib import Path

def check_iedb_tool(iedb_dir):
    tool_path = Path(iedb_dir) / "src" / "tcell_mhci.py"
    if not tool_path.exists():
        raise FileNotFoundError(f"IEDB tool not found at {tool_path}")
    return str(tool_path)

def find_first_fasta_file(pathogen_dir, sequence_dir):
    search_path = Path("data") / pathogen_dir / sequence_dir
    if not search_path.exists():
        raise FileNotFoundError(f"Directory not found: {search_path}")
    fasta_files = list(search_path.rglob("*.fasta"))
    return fasta_files[0] if fasta_files else None

def clean_fasta_file(original_fasta, cleaned_fasta):
    with open(original_fasta, "r") as infile, open(cleaned_fasta, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(line)
            else:
                # Remove asterisk characters
                cleaned_seq = line.replace("*", "").strip()
                outfile.write(cleaned_seq + "\n")

def generate_input_json(cleaned_fasta_path, output_json_path, allele_list, peptide_lengths):
    json_data = {
        "input_sequence_text_file_path": str(cleaned_fasta_path),
        "peptide_length_range": peptide_lengths,
        "alleles": ",".join(allele_list),
        "predictors": [
            {
                "type": "binding",
                "method": "netmhcpan_ba"
            }
        ]
    }
    with open(output_json_path, "w") as f:
        json.dump(json_data, f, indent=2)

def run_prediction(tool_path, json_file, output_prefix):
    cmd = [
        "python3", tool_path,
        "-j", json_file,
        "-o", output_prefix,
        "-f", "json"
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Prediction succeeded: {output_prefix}.json")
    except subprocess.CalledProcessError as e:
        print(f"❌ Prediction failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Test IEDB Class I T Cell prediction on a single FASTA file.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence directory inside pathogen_dir/")
    args = parser.parse_args()

    iedb_dir = input("Enter full path to IEDB tool folder (<50 chars): ").strip()
    if len(iedb_dir) > 50:
        raise ValueError("IEDB path must be under 50 characters.")
    tool_path = check_iedb_tool(iedb_dir)

    first_fasta = find_first_fasta_file(args.pathogen_dir, args.sequence_dir)
    if not first_fasta:
        print("No FASTA file found.")
        return

    print(f"🧪 Testing with: {first_fasta.name}")

    output_dir = Path("outputs") / args.pathogen_dir / args.sequence_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cleaned_fasta = output_dir / f"{first_fasta.stem}_cleaned.fasta"
    json_file = output_dir / f"{first_fasta.stem}_input.json"
    output_prefix = output_dir / f"{first_fasta.stem}"

    clean_fasta_file(first_fasta, cleaned_fasta)

    allele_list = ["HLA-A*02:01", "HLA-A*01:01"]
    peptide_lengths = [8, 11]

    generate_input_json(cleaned_fasta, json_file, allele_list, peptide_lengths)
    run_prediction(tool_path, str(json_file), str(output_prefix))

if __name__ == "__main__":
    main()
