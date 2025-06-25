"""
common.py
Common utility functions and constants for mRNA Target Selection pipeline.

This module provides shared utilities for directory management, file handling, tool validation,
allele panel selection, and FASTA/JSON conversion used throughout the mRNA Target Selection workflow.

General Function:
    - Ensures output and temporary directories exist and are writable.
    - Locates and validates external tool executables (SignalP, TargetP, TMHMM, IEDB tools, etc.).
    - Provides allele panel presets and selection logic for MHC-I and MHC-II.
    - Converts FASTA files to text and JSON formats for downstream processing.
    - Cleans up temporary directories and files.
    - Includes helper functions for parsing and preparing input/output files.

Constants:
    - MHCI_DEFAULT, MHCI_EXTENDED: Default and extended allele panels for MHC-I.
    - MHCII_DEFAULT, MHCII_EXTENDED: Default and extended allele panels for MHC-II.
    - ALLELE_PRESETS: Dictionary mapping tool types to their allele panels.

Author: Nadia
"""
import json
from pathlib import Path
import shutil
import tempfile
import csv
import logging
import re

def is_valid_peptide(seq: str) -> bool:
    """
    Checks if a given sequence is a valid peptide sequence.
    A valid peptide consists only of standard amino acid characters (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y).
    Args:
        seq (str): The peptide sequence to validate.
    Returns:
        bool: True if the sequence is valid, False otherwise.
    """
    return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq))

def parse_csv_to_fasta(csv_file: Path, output_dir: Path, basename_prefix: str) -> Path:
    """
    Parses a B-cell CSV file with peptide predictions and writes a FASTA file with contextual headers.
    Infers header info from the first valid 'input:' line and applies it to earlier peptides too.
    
    Args:
        csv_file (Path): Path to input CSV file.
        output_dir (Path): Directory where output FASTA will be saved.
        basename_prefix (str): Base prefix for output file name.

    Returns:
        Path: Path to the generated FASTA file, or None if no peptides found.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_lines = []
    seen = set()
    peptides = []
    current_header = None
    strain_acc = "unknown"

    early_peptides = []
    early_header_written = False

    method = csv_file.name.lower()
    is_bepipred = "bepipred" in method
    is_emini_or_kolaskar = any(x in method for x in ["emini", "kolaskar"])
    is_other = not (is_bepipred or is_emini_or_kolaskar)
    in_bepipred_peptide_block = False

    def write_block(header: str, block_peptides: list):
        if not header or not block_peptides:
            return
        header = header.lstrip(">")
        for i, pep in enumerate(block_peptides):
            if not is_valid_peptide(pep):
                logging.warning(f"⚠️ Invalid peptide skipped: {pep} in header: {header}")
                continue
            key = (header, pep)
            if key not in seen:
                fasta_lines.append(f">{header}|seq{i+1}\n{pep}")
                seen.add(key)

    with open(csv_file) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or all(cell.strip() == "" for cell in row):
                continue

            if row[0].startswith("input:"):
                # Finalize previous peptide block
                if peptides:
                    write_block(current_header, peptides)
                    peptides = []

                # Extract header info
                header_str = row[1] if len(row) > 1 else row[0][len("input:"):].strip()
                match = re.match(r"antigen_(\d+)\|([A-Z0-9]+)\|.*?\|([A-Z0-9.]+)", header_str)
                if match:
                    antigen_num, acc_num, strain_acc = match.groups()
                    current_header = f">antigen{antigen_num}|{acc_num}|{strain_acc}|bcell"

                    # If we had early peptides, write them now using inferred strain
                    if early_peptides and not early_header_written:
                        inferred_header = f">antigenUnknown|unknown|{strain_acc}|bcell"
                        write_block(inferred_header, early_peptides)
                        early_header_written = True
                else:
                    logging.warning(f"⚠️ Could not parse B-cell header: {row}")
                    current_header = f">antigenUnknown|unknown|{strain_acc}|bcell"
                in_bepipred_peptide_block = False
                continue

            # Ensure we have a default header if no 'input:' encountered yet
            if current_header is None:
                current_header = f">antigenUnknown|unknown|{strain_acc}|bcell"

            peptide = None

            # --- Peptide Parsing Logic ---
            if is_bepipred:
                if row[0].startswith("No") and any("Peptipe" in col for col in row):
                    in_bepipred_peptide_block = True
                    continue
                elif row[0].startswith("Position") and "Residue" in row:
                    in_bepipred_peptide_block = False
                    continue
                elif in_bepipred_peptide_block and len(row) >= 4:
                    peptide = row[3].strip()

            elif is_emini_or_kolaskar:
                if (row[0] == "Position" and len(row) >= 5 and row[4] == "Peptide") or \
                   (row[0] == "No" and "Peptipe" in row):
                    continue
                elif len(row) >= 5:
                    peptide = row[4].strip()
                elif len(row) >= 3:
                    peptide = row[2].strip()

            elif is_other:
                if row[0] == "Position" and len(row) >= 5 and row[4] == "Peptide":
                    continue
                elif len(row) >= 5:
                    peptide = row[4].strip()

            if peptide and is_valid_peptide(peptide):
                if current_header.startswith(">antigenUnknown"):
                    early_peptides.append(peptide)
                else:
                    peptides.append(peptide)

    # Final block
    if peptides:
        write_block(current_header, peptides)
    elif early_peptides and not early_header_written:
        # No 'input:' line at all? Still write early peptides
        fallback_header = f">antigenUnknown|unknown|{strain_acc}|bcell"
        write_block(fallback_header, early_peptides)

    if not fasta_lines:
        logging.warning(f"⚠️ No peptides found in B-cell file: {csv_file.name}")
        return None

    fasta_path = output_dir / f"{basename_prefix}.fasta"
    with open(fasta_path, "w") as f_out:
        f_out.write("\n".join(fasta_lines))

    print(f"💾 B-cell FASTA written: {fasta_path}")
    return fasta_path

def parse_json_to_fasta(json_file: Path, output_dir: Path, basename_prefix: str) -> Path:
    """
    Parses a JSON file with peptide predictions and writes a FASTA file of unique peptides
    with contextual headers based on the input file name, appending seq1, seq2, etc.

    Args:
        json_file (Path): Path to input JSON file.
        output_dir (Path): Directory where output FASTA will be saved.
        basename_prefix (str): Base prefix for output file name.

    Returns:
        Path: Path to the generated FASTA file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(json_file) as f:
        data = json.load(f)

    results = data.get("results", [])
    peptides = set()

    for entry in results:
        if entry.get("type") != "peptide_table":
            continue

        columns = entry.get("table_columns", [])
        data_rows = entry.get("table_data", [])

        try:
            peptide_idx = columns.index("peptide")
        except ValueError:
            continue  # No peptide column

        for row in data_rows:
            peptides.add(row[peptide_idx])

    if not peptides:
        print(f"⚠️ No peptides found in {json_file}")
        return None

    filename = json_file.name
    match = re.match(
    r"antigen_(\d+)_([A-Z0-9]+)_(.+?)_([A-Z0-9]+\.\d+).*_(MHCI|MHCII)\.json$",
    filename
    )

    if match:
        antigen_num, acc_num, _, strain_acc, mhc_class = match.groups()
        mhc_class = mhc_class.lower()
        fasta_lines = [
            f">antigen{antigen_num}|{acc_num}|{strain_acc}|{mhc_class}|seq{i+1}\n{pep}"
            for i, pep in enumerate(sorted(peptides))
        ]
    else:
        print(f"⚠️ Filename pattern not recognized for T-cell: {filename}")
        raise SystemExit(1)

    # Write to FASTA
    fasta_path = output_dir / f"{basename_prefix}.fasta"
    with open(fasta_path, "w") as fasta_file:
        fasta_file.write("\n".join(fasta_lines))

    print(f"💾 FASTA written: {fasta_path}")
    return fasta_path

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
    """
    Get all FASTA files in the specified sequence subdirectory.
    Args:
        base_path (Path): Base directory path.
        sequence_subdir (str): Subdirectory name containing FASTA files.
    Returns:
        list: List of Path objects for each FASTA file found.
    """
    seq_dir = base_path / sequence_subdir
    if not seq_dir.exists() or not seq_dir.is_dir():
        print(f"❌ Sequence directory {seq_dir} does not exist or is not a directory.")
        return []
    fasta_files = list(seq_dir.glob("*.fasta"))
    return fasta_files

def prepare_output_dirs(pathogen_path, output_subdir, selected_tools):
    """
    Prepare output directories for the specified pathogen and tools.
    Args:
        pathogen_path (Path): Path to the pathogen directory.
        output_subdir (str): Subdirectory name for output files.
        selected_tools (list): List of selected tools to create subdirectories for.
    Returns:
        tuple: Paths to the temporary JSON directory and the output directory.
    """
    output_dir = pathogen_path / output_subdir

    # Create only the needed tool subdirectories
    for tool in selected_tools:
        (output_dir / tool.lower()).mkdir(parents=True, exist_ok=True)

    # Only create temp_json_dir if needed
    if selected_tools in ["MHCI", "MHCII", "BCell"]:
        temp_json_dir = pathogen_path / "temp_json"
        temp_json_dir.mkdir(parents=True, exist_ok=True)
        return temp_json_dir, output_dir
    
    else:
        # If no temp_json_dir is needed, return None
        return None, output_dir

def cleanup_temp(temp_dir: Path):
    """
    Clean up temporary directory if it exists.
    Args:
        temp_dir (Path): Path to the temporary directory to clean up.
    """
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
        print(f"Cleaned temporary: {temp_dir}")

""" Constants for MHC Allele Panels """
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
    """
    Get the list of alleles for the specified tool type and panel.
    Args:
        tool_type (str): Type of tool (e.g., "MHCI", "MHCII").
        panel (str): Name of the allele panel ("default", "extended", or "custom").
        custom_alleles (list, optional): List of custom alleles if panel is "custom".
    Returns:
        list: List of alleles for the specified tool type and panel.
    """
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
    Args:
        tool_root (Path): Root directory where tools are expected to be located.
    Returns:
        dict: Dictionary mapping tool names to their executable paths.
    Raises:
        FileNotFoundError: If any of the required tools are not found in the expected locations.
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
    """
    Check for the presence of IEDB epitope prediction tools and return their paths.
    Args:
        base_path (str): Base path where IEDB tools are expected to be located.
    Returns:
        dict: Dictionary mapping tool names to their executable paths.
    """
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
    """
    Convert FASTA files to text format by copying their contents to .txt files.
    Args:
        fasta_files (list): List of Path objects for FASTA files.
        temp_txt_dir (Path): Directory where the converted text files will be saved.
    Returns:
        list: List of Path objects for the created text files.
    """
    temp_txt_dir.mkdir(parents=True, exist_ok=True)
    txt_files = []
    for fasta_file in fasta_files:
        txt_file = temp_txt_dir / (fasta_file.stem + ".txt")
        txt_file.write_text(fasta_file.read_text())
        txt_files.append(txt_file)
    return txt_files

def write_json(seq_id_line, seq_lines, temp_dir, alleles, peptide_lengths, tool_type, strain_name, method="netmhcpan_el"):
    """
    Write a JSON file for a given sequence ID and its associated sequence lines.
    Args:
        seq_id_line (str): The sequence ID line from the FASTA file.
        seq_lines (list): List of sequence lines corresponding to the ID.
        temp_dir (Path): Directory where the JSON file will be saved.
        alleles (list): List of alleles to include in the JSON.
        peptide_lengths (tuple): Tuple specifying the peptide length range (min, max).
        tool_type (str): Type of tool for which this JSON is being generated.
        strain_name (str): Name of the strain associated with this sequence.
        method (str): Method to use for prediction, default is "netmhcpan_el".
    Returns:
        json_path (Path): Path to the created JSON file, or None if an error occurred.
    """
    header = seq_id_line.strip()
    antigen_id = header[1:].split()[0]
    sequence = "".join(seq_lines).replace("*", "").strip()
    
    if not sequence:
        print(f"⚠️ Empty sequence for {antigen_id}")
        return None

    if method is None:
        print(f"⚠️ Unknown tool type '{tool_type}' for {antigen_id}, skipping JSON generation.")
        return None

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
    """
    Parse a FASTA file and convert each sequence to a JSON file.
    Args:
        fasta_path (str): Path to the input FASTA file.
        temp_dir (Path): Directory where JSON files will be saved.
        alleles (list): List of alleles to include in the JSON.
        peptide_lengths (tuple): Tuple specifying the peptide length range (min, max).
        tool_type (str): Type of tool for which this JSON is being generated.
        strain_name (str): Name of the strain associated with this sequence.
    Returns:
        list: List of paths to the created JSON files.
    """
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
    Args:
        tool_root (Path): Root directory where tools are expected to be located.
    Returns:
        dict: Dictionary mapping tool names to their executable paths.
    """
    tool_map = {
        "Allergenicity": tool_root / "algpred2" / "algpred2.py",
        "PopCoverage": tool_root / "population_coverage" /"calculate_population_coverage.py",
        "Cluster": tool_root / "miniconda3" / "bin" / "mmseqs"
    }

    found = {}
    for name, path in tool_map.items():
        if path.exists():
            found[name] = str(path.parent)
        else:
            print(f"❌ {name} tool not found at: {path.parent}")

    return found

# Add __main__ with unittest
if __name__ == "__main__":
    """ Main function for unit testing of allele retrieval functionality. """
    import unittest
    class AlleleTests(unittest.TestCase):
        def test_default_mhci(self):
            self.assertTrue(len(get_alleles("MHCI")) >= len(MHCI_DEFAULT))

        def test_custom(self):
            custom = ["HLA-X*01", "HLA-Y*02"]
            aa = get_alleles("MHCI", "custom", custom)
            self.assertEqual(aa, custom)

        def test_invalid_panel(self):
            aa = get_alleles("MHCII", "NOT-FOUND")
            self.assertEqual(aa, MHCII_DEFAULT)

    class JsonToFastaTests(unittest.TestCase):
        def test_parse_json_to_fasta(self):
            # Mock JSON data
            json_data = {
                "results": [
                    {
                    "method": "binding.netmhcpan_el",
                    "type": "peptide_table",
                    "table_columns": [
                        "allele",
                        "peptide",
                        "core",
                        "icore",
                        "score",
                        "percentile"
                    ],
                    "table_data": [
                        [
                        "HLA-A*31:01",
                        "RLNKYTLHR",
                        "RLNKYTLHR",
                        "RLNKYTLHR",
                        0.947356,
                        0.01
                        ],
                        [
                        "HLA-A*23:01",
                        "KYCPRLNKYTL",
                        "KYCPNKYTL",
                        "KYCPRLNKYTL",
                        0.408173,
                        0.22
                        ],
                        [
                        "HLA-A*03:01",
                        "RLNKYTLHR",
                        "RLNKYTLHR",
                        "RLNKYTLHR",
                        0.914693,
                        0.03
                        ]
                    ]
                    }
                ]
            }
            # Create a permanent JSON file with a filename matching the required pattern
            json_filename = "antigen_12_A0A2S1FUJ1_something_BA1.2_more_MHCI.json"
            json_file = Path("/tmp") / json_filename
            json_file.write_text(json.dumps(json_data))

            # Use a permanent output directory and filename for inspection
            output_dir = Path("/tmp/test_json_to_fasta")
            output_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = parse_json_to_fasta(json_file, output_dir, "test_output")

            self.assertIsNotNone(fasta_path)
            self.assertTrue(fasta_path.exists())
            self.assertEqual(fasta_path.suffix, ".fasta")
            fasta_content = fasta_path.read_text().strip()
            # The header should match the new contextual header format
            self.assertTrue(fasta_content.startswith(">"))
            self.assertIn("antigen12|A0A2S1FUJ1|BA1.2|mhci|seq1", fasta_content)
            self.assertIn("RLNKYTLHR", fasta_content)
            self.assertIn("KYCPRLNKYTL", fasta_content)
            # Only two unique peptides, so only seq1 and seq2 should be present
            self.assertIn("seq1", fasta_content)
            self.assertIn("seq2", fasta_content)
            self.assertNotIn("seq3", fasta_content)

            #clean up
            if fasta_path.exists():
                fasta_path.unlink() # remove the created FASTA file
            if json_file.exists():
                json_file.unlink() # remove the created JSON file
            if output_dir.exists():
                shutil.rmtree(output_dir)

    
    class CSVToFastaTests(unittest.TestCase):
        def test_parse_csv_to_fasta(self):
            # Mock CSV data with two antigen blocks and some duplicate peptides
            csv_data = """input: antigen_86|A0A2S1FUJ1|Superantigen-like|HE681097.1|tpos:1234
Position,Residue,Score,Length,Peptide Sequence
1,1,0.95,9,RLNKYTLHR
2,2,0.90,10,KYCPRLNKYTL
3,3,0.85,11,RLNKYTLHR
input: antigen_87|A0A2S1FUJ2|Superantigen-like|HE681098.1|tpos:5678
Position,Residue,Score,Length,Peptide Sequence
1,1,0.92,8,RLNKYTLH
2,2,0.88,9,KYCPRLNKYT
3,3,0.95,7,plot
"""

            # Create a temporary CSV file
            csv_file = Path("/tmp/test_bcell.csv")
            csv_file.write_text(csv_data)

            # Output directory and expected fasta output
            output_dir = Path("/tmp/test_csv_to_fasta")
            output_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = parse_csv_to_fasta(csv_file, output_dir, "test_output")

            # Assertions
            self.assertIsNotNone(fasta_path)
            self.assertTrue(fasta_path.exists())
            self.assertEqual(fasta_path.suffix, ".fasta")

            fasta_content = fasta_path.read_text().strip().splitlines()
            headers = [line for line in fasta_content if line.startswith(">")]
            sequences = [line for line in fasta_content if not line.startswith(">")]

            # Check expected structure
            self.assertTrue(headers[0].startswith(">antigen86|A0A2S1FUJ1|HE681097.1|bcell|seq1"))
            self.assertIn("RLNKYTLHR", sequences)
            self.assertIn("KYCPRLNKYTL", sequences)
            self.assertIn("RLNKYTLH", sequences)
            self.assertIn("KYCPRLNKYT", sequences)
            self.assertNotIn("plot", sequences)
            self.assertNotIn("seq4", sequences)

            # Ensure all sequences are valid
            for seq in sequences:
                self.assertTrue(is_valid_peptide(seq), f"Invalid peptide in FASTA: {seq}")

            # Make sure duplicate peptide is not included twice
            self.assertEqual(sequences.count("RLNKYTLHR"), 1)

            # # Clean up files and directory
            if fasta_path.exists():
                fasta_path.unlink()
            if csv_file.exists():
                csv_file.unlink()
            if output_dir.exists():
                shutil.rmtree(output_dir)

        def test_bepipred_csv_to_fasta(self):
            csv_data = """Position,Residue,Score,Assignment
1,M,-0.481,.
2,A,-0.099,.
3,I,0.070,.
4,S,0.241,.
input:,antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462
Predicted,peptides
No,Start,End,Peptipe,Length
1,10,13,TFNK,4
2,39,44,YSGAGK,6
3,56,59,AASS,4
4,78,81,been,4
Position,Residue,Score,Assignment
        """

            csv_file = Path("/tmp/bepipred.csv")
            csv_file.write_text(csv_data)
            output_dir = Path("/tmp/test_bepipred_output")
            output_dir.mkdir(parents=True, exist_ok=True)

            fasta_path = parse_csv_to_fasta(csv_file, output_dir, "bepipred_test")
            self.assertIsNotNone(fasta_path)
            fasta_lines = fasta_path.read_text().strip().splitlines()

            headers = [l for l in fasta_lines if l.startswith(">")]
            sequences = [l for l in fasta_lines if not l.startswith(">")]

            self.assertEqual(len(headers), 3)
            self.assertIn("TFNK", sequences)
            self.assertIn("YSGAGK", sequences)
            self.assertIn("AASS", sequences)
            self.assertNotIn("been", sequences)

            # Cleanup
            fasta_path.unlink()
            csv_file.unlink()
            shutil.rmtree(output_dir)

        def test_emini_csv_to_fasta(self):
            csv_data = """Position,Residue,Start,End,Peptide,Score
3,I,1,6,MAISQE,0.396
4,S,2,7,AISQER,0.783
5,Q,3,8,ISQERK,1.551
6,E,4,9,SQERKN,3.557
input:,antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462
Predicted,peptides
No,Start,End,Peptipe,Length
10,17,TFNKKKQK,4.009875
69,81,ITNYSEKGMREIK,1.9206153846153846
111,121,KKSKTEIKQRV,2.6391818181818176
131,138,SDKKDQFP,2.94825
Position,Residue,Start,End,Peptide,Score
3,E,1,6,MIEFRQ,0.775
4,F,2,7,IEFRQV,0.581
5,R,3,8,error,0.514
        """

            csv_file = Path("/tmp/emini.csv")
            csv_file.write_text(csv_data)
            output_dir = Path("/tmp/test_emini_output")
            output_dir.mkdir(parents=True, exist_ok=True)

            fasta_path = parse_csv_to_fasta(csv_file, output_dir, "emini_test")
            self.assertIsNotNone(fasta_path)
            fasta_lines = fasta_path.read_text().strip().splitlines()

            headers = [l for l in fasta_lines if l.startswith(">")]
            sequences = [l for l in fasta_lines if not l.startswith(">")]

            self.assertEqual(len(headers), 10)
            self.assertIn("MAISQE", sequences)
            self.assertIn("AISQER", sequences)
            self.assertIn("ISQERK", sequences)
            self.assertIn("SQERKN", sequences)
            self.assertIn("TFNKKKQK", sequences)
            self.assertIn("ITNYSEKGMREIK", sequences)
            self.assertIn("KKSKTEIKQRV", sequences)
            self.assertIn("SDKKDQFP", sequences)
            self.assertIn("MIEFRQ", sequences)
            self.assertIn("IEFRQV", sequences)
            self.assertNotIn("error", sequences)

            # Cleanup
            fasta_path.unlink()
            csv_file.unlink()
            shutil.rmtree(output_dir)

        def test_missing_initial_input_header(self):
            csv_data = """Position,Residue,Start,End,Peptide,Score
4,S,1,7,MAISQER,2.057
5,Q,2,8,AISQERK,3.471
6,E,3,9,ISQERKN,4.171
input:,antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462
Position,Residue,Start,End,Peptide,Score
4,F,1,7,MIEFRQV,-1.014
5,R,2,8,IEFRQVS,0.514
6,V,3,9,test,1.214
        """

            csv_file = Path("/tmp/missing_header.csv")
            csv_file.write_text(csv_data)
            output_dir = Path("/tmp/test_missing_header_output")
            output_dir.mkdir(parents=True, exist_ok=True)

            fasta_path = parse_csv_to_fasta(csv_file, output_dir, "missing_header_test")
            self.assertIsNotNone(fasta_path)
            fasta_lines = fasta_path.read_text().strip().splitlines()

            headers = [l for l in fasta_lines if l.startswith(">")]
            sequences = [l for l in fasta_lines if not l.startswith(">")]

            # Should contain 5 entries in total
            self.assertEqual(len(headers), 5)
            self.assertIn("MAISQER", sequences)
            self.assertIn("AISQERK", sequences)
            self.assertIn("ISQERKN", sequences)
            self.assertIn("MIEFRQV", sequences)
            self.assertIn("IEFRQVS", sequences)
            self.assertNotIn("test", sequences)
            
            # First header should use "antigenUnknown"
            self.assertTrue(headers[0].startswith(">antigenUnknown|unknown|unknown|bcell"))

            # Second block should use parsed antigen_149 header
            self.assertTrue(headers[3].startswith(">antigen149|Q99QV7|HE681097.1|bcell"))

            # Cleanup
            fasta_path.unlink()
            csv_file.unlink()
            shutil.rmtree(output_dir)



    unittest.main()