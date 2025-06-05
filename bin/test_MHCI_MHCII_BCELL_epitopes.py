import argparse
import os
import subprocess
import json
import time
import traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def check_iedb_tool(iedb_dir):
    if "bcell" in iedb_dir.lower():
        tool_type = "BCell"
        tool_path = Path(iedb_dir) / "bcell_standalone.py"
    elif "tc1" in iedb_dir.lower():
        tool_type = "MHCI"
        tool_path = Path(iedb_dir) / "src" / "tcell_mhci.py"
    elif "tc2" in iedb_dir.lower():
        tool_type = "MHCII"
        tool_path = Path(iedb_dir) / "src" / "tcell_mhcii.py"
    else:
        raise ValueError("Unable to infer tool type from path. Include 'bcell', 'tc1', or 'tc2' in the path.")

    if not tool_path.exists():
        raise FileNotFoundError(f"Tool not found at {tool_path}")

    return str(tool_path), tool_type

def parse_fasta_to_jsons(fasta_path, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    json_paths = []
    total = 0

    with open(fasta_path, 'r') as infile:
        seq_id = None
        seq_data = []

        for line in infile:
            if line.startswith(">"):
                if seq_id:
                    json_path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
                    if json_path:
                        json_paths.append(json_path)
                        total += 1
                        print(f"📝 Wrote JSON for {seq_id}")
                seq_id = line.strip()
                seq_data = []
            else:
                seq_data.append(line.strip())

        if seq_id and seq_data:
            json_path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
            if json_path:
                json_paths.append(json_path)
                total += 1
                print(f"📝 Wrote JSON for {seq_id}")

    print(f"📦 Total sequences parsed: {total}")
    return json_paths

def write_json(seq_id_line, seq_lines, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    header = seq_id_line.strip()
    if not header.startswith(">"):
        print(f"⚠️ Invalid FASTA header: {header}")
        return None

    antigen_id = header[1:].split()[0]
    sequence = "".join(seq_lines).replace("*", "").strip()

    if not sequence:
        print(f"⚠️ Empty sequence for {antigen_id}")
        return None

    input_sequence_text = f">{antigen_id}\n{sequence}"

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

    safe_antigen_id = antigen_id.replace(" ", "_").replace("/", "_").replace("|", "_")
    filename = f"{safe_antigen_id}_{strain_name}_{tool_type}.json"
    json_path = Path(temp_dir) / filename

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
    except Exception:
        print(f"❌ Exception during processing {json_file.name}")
        traceback.print_exc()

def main():
    parser = argparse.ArgumentParser(description="Run IEDB predictions in parallel for FASTA sequences.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence directory inside pathogen_dir/")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads for sequence prediction")
    args = parser.parse_args()

    iedb_dir = input("Enter full path to IEDB tool folder: ").strip()
    tool_path, tool_type = check_iedb_tool(iedb_dir)

    strain_name = args.sequence_dir
    base_path = Path("data") / args.pathogen_dir

    search_path = base_path / args.sequence_dir
    fasta_files = list(search_path.glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {search_path}")
        return

    temp_json_dir = base_path / "temp_jsons"
    output_dir = base_path / "epitope_outputs"
    temp_json_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    allele_list = ["HLA-A*02:01", "HLA-A*01:01"]
    peptide_lengths = [8, 11]

    all_json_files = []
    for fasta_file in fasta_files:
        print(f"🧬 Processing FASTA file: {fasta_file.name}")
        json_files = parse_fasta_to_jsons(fasta_file, temp_json_dir, allele_list, peptide_lengths, tool_type, strain_name)
        all_json_files.extend(json_files)

    # 🔁 Run predictions in parallel
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(run_prediction, tool_path, json_file, output_dir)
            for json_file in all_json_files
        ]
        for f in futures:
            f.result()

if __name__ == "__main__":
    main()
