import json
from pathlib import Path
import tempfile

def ensure_writable_dir(path: Path) -> bool:
    """
    Ensure the directory exists and is writable.

    Args:
        path (Path): Directory path to check/create.

    Returns:
        bool: True if the directory exists and is writable, False otherwise.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"❌ Could not create directory {path}: {e}")
        return False

    try:
        with tempfile.TemporaryFile(dir=path):
            pass
    except Exception as e:
        print(f"❌ Directory not writable: {path} ({e})")
        return False

    return True

def get_fasta_files(base_path: Path, sequence_subdir: str):
    seq_dir = base_path / sequence_subdir
    if not seq_dir.exists() or not seq_dir.is_dir():
        print(f"❌ Sequence directory {seq_dir} does not exist or is not a directory.")
        return []
    fasta_files = list(seq_dir.glob("*.fasta"))
    return fasta_files

def prepare_output_dirs(pathogen_path, output_subdir, selected_tools):
    output_dir = pathogen_path / output_subdir

    # Create only the needed tool subdirectories
    for tool in selected_tools:
        (output_dir / tool.lower()).mkdir(parents=True, exist_ok=True)

    temp_json_dir = pathogen_path / "temp_json"
    temp_json_dir.mkdir(parents=True, exist_ok=True)

    return temp_json_dir, output_dir

def cleanup_temp(temp_json_dir: Path):
    import shutil
    if temp_json_dir.exists() and temp_json_dir.is_dir():
        shutil.rmtree(temp_json_dir)
        print("🧹 Cleaned up temporary JSON files.")


MHCI_DEFAULT = [
    "HLA-A*02:01", "HLA-A*01:01"
]

MHCI_EXTENDED = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01",
    "HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"
]

MHCII_DEFAULT = [
    "HLA-DRB1*03:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01", "HLA-DRB3*01:01",
    "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01"
]

MHCII_EXTENDED = [
    "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05", "HLA-DRB1*07:01", "HLA-DRB1*08:02",
    "HLA-DRB1*09:01", "HLA-DRB1*11:01", "HLA-DRB1*12:01", "HLA-DRB1*13:02", "HLA-DRB1*15:01", "HLA-DRB3*01:01",
    "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*05:01/DQB1*03:01",
    "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*04:01/DQB1*04:02", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02",
    "HLA-DPA1*02:01/DPB1*01:01", "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*01:03/DPB1*04:01", "HLA-DPA1*03:01/DPB1*04:02",
    "HLA-DPA1*02:01/DPB1*05:01", "HLA-DPA1*02:01/DPB1*14:01"
]

ALLELE_PRESETS = {
    "MHCI": {
        "default": MHCI_DEFAULT,
        "extended": MHCI_EXTENDED
    },
    "MHCII": {
        "default": MHCII_DEFAULT,
        "extended": MHCII_EXTENDED
    }
}

def get_alleles(tool_type, panel="default", custom_alleles=None):
    panel = panel.lower()
    if panel == "custom":
        if not custom_alleles:
            print(f"⚠️ Custom allele panel selected but no alleles provided for {tool_type}. Using default panel.")
            return ALLELE_PRESETS[tool_type]["default"]
        return [a.strip() for a in custom_alleles]
    elif panel in ALLELE_PRESETS[tool_type]:
        return ALLELE_PRESETS[tool_type][panel]
    else:
        print(f"⚠️ Invalid allele panel '{panel}' for {tool_type}, using default.")
        return ALLELE_PRESETS[tool_type]["default"]
    
def check_signalp_targetp_tmhmm(tool_root: Path) -> dict:
    """
    Validates and returns paths to SignalP, TargetP, and TMHMM executables.
    """
    tools = {}

    # Check SignalP
    signalp_candidates = [
        tool_root / "signalp-5.0b" / "bin" / "signalp",
        tool_root / "signalp-5.0b" / "bin" / "bin" / "signalp"
    ]
    tools["SIGNALP"] = next((p for p in signalp_candidates if p.exists()), None)

    # Check TargetP
    targetp_path = tool_root / "targetp-2.0" / "bin" / "targetp"
    tools["TARGETP"] = targetp_path if targetp_path.exists() else None

    # Check TMHMM
    tmhmm_path = tool_root / "tmhmm-2.0c" / "bin" / "tmhmm"
    tools["TMHMM"] = tmhmm_path if tmhmm_path.exists() else None

    # Log missing tools
    for tool, path in tools.items():
        if not path:
            print(f"❌ {tool} not found in expected location.")

    # Raise error if any are missing
    if not all(tools.values()):
        raise FileNotFoundError("One or more required tools were not found under the given tool root.")

    return tools


def check_iedb_tool(base_path):
    base = Path(base_path)
    paths = {
        "BCell": base / "bcell_standalone" / "predict_antibody_epitope.py",
        "MHCI": base / "ng_tc1-0.1.2-beta" / "src" / "tcell_mhci.py",
        "MHCII": base / "ng_tc2-0.1.1-beta" / "src" / "tcell_mhcii.py"
    }

    tools = {}
    for key, path in paths.items():
        if path.exists():
            tools[key] = str(path)
    return tools

def convert_fasta_to_txt(fasta_files, temp_txt_dir: Path):
    temp_txt_dir.mkdir(parents=True, exist_ok=True)
    txt_files = []
    for fasta_file in fasta_files:
        txt_file = temp_txt_dir / (fasta_file.stem + ".txt")
        txt_file.write_text(fasta_file.read_text())
        txt_files.append(txt_file)
    return txt_files

def write_json(seq_id_line, seq_lines, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    header = seq_id_line.strip()
    antigen_id = header[1:].split()[0]
    sequence = "".join(seq_lines).replace("*", "").strip()
    if not sequence:
        print(f"⚠️ Empty sequence for {antigen_id}")
        return None
    if tool_type == "MHCI":
        method = "netmhcpan_ba"
    elif tool_type == "MHCII":
        method = "netmhciipan_el"  # or _ba depending on your choice
    else:
        method = None  # Or skip

    json_data = {
        "input_sequence_text": f">{antigen_id}\n{sequence}",
        "peptide_length_range": peptide_lengths,
        "alleles": ",".join(a.strip() for a in alleles),
        "predictors": [{"type": "binding", "method": method}]
    }

    safe_antigen_id = antigen_id.replace(" ", "_").replace("/", "_").replace("|", "_")
    filename = f"{safe_antigen_id}_{strain_name}_{tool_type}.json"
    json_path = Path(temp_dir) / filename

    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    return json_path


def parse_fasta_to_jsons(fasta_path, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
    json_paths, seq_id, seq_data = [], None, []
    with open(fasta_path, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                if seq_id:
                    path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
                    if path: json_paths.append(path)
                seq_id, seq_data = line.strip(), []
            else:
                seq_data.append(line.strip())
        if seq_id:
            path = write_json(seq_id, seq_data, temp_dir, alleles, peptide_lengths, tool_type, strain_name)
            if path: json_paths.append(path)
    return json_paths

def check_epitope_evaluation_tools(tool_root: Path) -> dict:
    """
    Detects available evaluation tools (Allergenicity, Population Coverage, Cluster).
    Returns a dictionary mapping tool names to their runner paths.
    """
    tool_map = {
        "Allergenicity": tool_root / "algpred2" / "run_algpred2.py",
        "PopCoverage": tool_root / "population_coverage" /"calculate_population_coverage.py",
        "Cluster": tool_root / "epitope_cluster" / "cluster.py",
    }

    found = {}
    for name, path in tool_map.items():
        if path.exists():
            found[name] = str(path.parent)
        else:
            print(f"❌ {name} tool not found at: {path.parent}")

    return found