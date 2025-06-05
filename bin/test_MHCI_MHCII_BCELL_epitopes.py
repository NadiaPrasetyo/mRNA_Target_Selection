import argparse
import os
import subprocess
import json
import time
import traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def check_iedb_tool(iedb_dir):
    tool_path = Path(iedb_dir) / "src" / "tcell_mhci.py"
    if not tool_path.exists():
        raise FileNotFoundError(f"IEDB tool not found at {tool_path}")
    return str(tool_path)

def parse_fasta_to_jsons(fasta_path, output_dir, alleles, peptide_lengths):
    json_paths = []
    total = 0

    with open(fasta_path, 'r') as infile:
        seq_id = None
        seq_data = []

        for line in infile:
            if line.startswith(">"):
                if seq_id:
                    json_path = write_json(seq_id, seq_data, output_dir, alleles, peptide_lengths)
                    if json_path:
                        json_paths.append(json_path)
                        print(f"📝 Wrote JSON for {seq_id}")
                        total += 1
                seq_id = line.strip()
                seq_data = []
            else:
                seq_data.append(line.strip())

        if seq_id and seq_data:
            json_path = write_json(seq_id, seq_data, output_dir, alleles, peptide_lengths)
            if json_path:
                json_paths.append(json_path)
                print(f"📝 Wrote JSON for {seq_id}")
                total += 1

    print(f"📦 Total sequences parsed: {total}")
    return json_paths

def write_json(seq_id_line, seq_lines, output_dir, alleles, peptide_lengths):
    header = seq_id_line.strip()
    if not header.startswith(">"):
        print(f"⚠️ Invalid FASTA header: {header}")
        return None

    header = header[1:]
    sequence = "".join(seq_lines).replace("*", "").strip()

    if not sequence:
        print(f"⚠️ Empty sequence for {header}")
        return None

    input_sequence_text = f">{header}\n{sequence}"

    json_data = {
        "input_sequence_text": input_sequence_text,
        "peptide_length_range": peptide_lengths,
        "alleles": ",".join(alleles),
        "predictors": [
            {
                "type": "binding",
                "method": "netmhcpan_ba"
            }
        ]
    }

    name = header.replace(" ", "_").replace("/", "_").replace("|", "_")[:40]
    json_path = Path(output_dir) / f"{name}.json"

    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)

    return json_path

def run_prediction(tool_path, json_file, output_dir):
    out_base = Path(json_file).stem
    output_prefix = Path(output_dir) / out_base
    cmd = [
        "python3", tool_path,
        "-j", str(json_file),
        "-o", str(output_prefix),
        "-f", "json"
    ]
    print(f"🔄 Running prediction for: {json_file.name}")

    try:
        start_time = time.time()
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        elapsed = time.time() - start_time
        print(f"✅ Success: {json_file.name} ({elapsed:.2f}s)")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {json_file.name}")
        print("📄 Command:", ' '.join(cmd))
        print("📤 STDOUT:")
        print(e.stdout.strip())
        print("📥 STDERR:")
        print(e.stderr.strip())
    except Exception as ex:
        print(f"❌ Exception during processing {json_file.name}")
        traceback.print_exc()

def process_fasta(tool_path, fasta_file, output_dir, allele_list, peptide_lengths):
    print(f"🧬 Starting FASTA: {fasta_file.name}")
    json_files = parse_fasta_to_jsons(fasta_file, output_dir, allele_list, peptide_lengths)
    for json_file in json_files:
        run_prediction(tool_path, json_file, output_dir)

def main():
    parser = argparse.ArgumentParser(description="Run IEDB MHCI predictions split by FASTA sequence.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence directory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=2, help="Number of parallel FASTA files to process")
    args = parser.parse_args()

    iedb_dir = input("Enter full path to IEDB tool folder (<50 chars): ").strip()
    if len(iedb_dir) > 50:
        raise ValueError("IEDB path must be under 50 characters.")
    tool_path = check_iedb_tool(iedb_dir)

    search_path = Path("data") / args.pathogen_dir / args.sequence_dir
    fasta_files = list(search_path.glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {search_path}")
        return

    output_dir = Path("outputs") / args.pathogen_dir / args.sequence_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    allele_list = ["HLA-A*02:01", "HLA-A*01:01"]
    peptide_lengths = [8, 11]

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(process_fasta, tool_path, fasta_file, output_dir, allele_list, peptide_lengths)
            for fasta_file in fasta_files
        ]
        for f in futures:
            f.result()

if __name__ == "__main__":
    main()
