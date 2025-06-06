import argparse
import os
import subprocess
import json
import time
import shutil
import traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# Define MHC allele presets
ALLELE_PRESETS = {
    "mhci-27": [
        "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01",
        "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
        "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01",
        "HLA-A*68:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
        "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01",
        "HLA-B*57:01", "HLA-B*58:01"
    ],
    "mhcii-7": [
        "HLA-DRB1*03:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01", "HLA-DRB3*01:01",
        "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01"
    ],
    "mhcii-27": [
        "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05", "HLA-DRB1*07:01",
        "HLA-DRB1*08:02", "HLA-DRB1*09:01", "HLA-DRB1*11:01", "HLA-DRB1*12:01", "HLA-DRB1*13:02",
        "HLA-DRB1*15:01", "HLA-DRB3*01:01", "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01",
        "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*05:01/DQB1*03:01", "HLA-DQA1*03:01/DQB1*03:02",
        "HLA-DQA1*04:01/DQB1*04:02", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02",
        "HLA-DPA1*02:01/DPB1*01:01", "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*01:03/DPB1*04:01",
        "HLA-DPA1*03:01/DPB1*04:02", "HLA-DPA1*02:01/DPB1*05:01", "HLA-DPA1*02:01/DPB1*14:01"
    ]
}

def check_iedb_tool(tool_root):
    tool_root = Path(tool_root)
    tools = {
        "BCell": tool_root / "bcell_standalone" / "bcell_standalone.py",
        "MHCI": tool_root / "ng_tc1-0.1.2-beta" / "src" / "tcell_mhci.py",
        "MHCII": tool_root / "ng_tc2-0.1.1-beta" / "src" / "tcell_mhcii.py"
    }

    available_tools = [(k, v) for k, v in tools.items() if v.exists()]
    if not available_tools:
        raise FileNotFoundError("No valid IEDB tools found in the provided directory.")

    # Return all found tools
    return dict(available_tools)

def parse_fasta_to_jsons(fasta_path, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    json_paths = []
    total = 0
    with open(fasta_path, 'r') as infile:
        seq_id, seq_data = None, []
        for line in infile:
            if line.startswith(">"):
                if seq_id:
                    json_path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
                    if json_path:
                        json_paths.append(json_path)
                        total += 1
                seq_id, seq_data = line.strip(), []
            else:
                seq_data.append(line.strip())
        if seq_id and seq_data:
            json_path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
            if json_path:
                json_paths.append(json_path)
                total += 1
    return json_paths

def write_json(seq_id_line, seq_lines, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    antigen_id = seq_id_line[1:].split()[0]
    sequence = "".join(seq_lines).replace("*", "").strip()
    if not sequence:
        return None
    input_sequence_text = f">{antigen_id}\n{sequence}"
    json_data = {
        "input_sequence_text": input_sequence_text,
        "peptide_length_range": peptide_lengths,
        "alleles": ",".join(alleles),
        "predictors": [{"type": "binding", "method": "netmhcpan_ba"}]
    }
    safe_id = antigen_id.replace(" ", "_").replace("/", "_").replace("|", "_")
    json_path = Path(temp_dir) / f"{safe_id}_{strain_name}_{tool_type}.json"
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    return json_path

def run_prediction(tool_path, json_file, output_dir):
    output_prefix = Path(output_dir) / Path(json_file).stem
    cmd = ["python3", tool_path, "-j", str(json_file), "-o", str(output_prefix), "-f", "json"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        print(f"✅ Success: {json_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {json_file.name}")
        print(e.stderr)

def main():
    parser = argparse.ArgumentParser(description="Run IEDB predictions for epitope analysis.")
    parser.add_argument("pathogen_dir", help="Pathogen directory inside data/")
    parser.add_argument("sequence_dir", help="Sequence directory inside pathogen_dir/")
    parser.add_argument("--tool-root", required=True, help="Root directory containing IEDB tools")
    parser.add_argument("--allele-panel", choices=["mhci-27", "mhcii-7", "mhcii-27", "custom"], default="mhci-27",
                        help="Allele panel to use for prediction")
    parser.add_argument("--custom-alleles", nargs="+", help="List of custom alleles (only if --allele-panel custom)")
    parser.add_argument("--peptide-lengths", "-pl", nargs="+", type=int, default=[9, 10, 11],
                        help="Peptide lengths to consider for prediction")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    args = parser.parse_args()

    tools = check_iedb_tool(args.tool_root)
    tool_type = "MHCI" if "mhci" in args.allele_panel else "MHCII"
    tool_path = tools.get(tool_type)
    if not tool_path:
        raise RuntimeError(f"Tool {tool_type} not found in {args.tool_root}")

    alleles = args.custom_alleles if args.allele_panel == "custom" else ALLELE_PRESETS[args.allele_panel]

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

    all_json_files = []
    for fasta_file in fasta_files:
        jsons = parse_fasta_to_jsons(fasta_file, temp_json_dir, alleles, args.peptide_lengths, tool_type, args.sequence_dir)
        all_json_files.extend(jsons)

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(run_prediction, tool_path, json_file, output_dir) for json_file in all_json_files]
        for f in futures: f.result()

    shutil.rmtree(temp_json_dir)
    print("🧹 Cleaned up temporary JSONs.")

if __name__ == "__main__":
    main()
