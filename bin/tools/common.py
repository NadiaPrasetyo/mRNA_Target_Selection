"""
common.py
Common utility functions and constants for mRNA Target Selection pipeline.

This module provides shared utilities for directory management, file handling, tool validation,
allele panel selection, and FASTA/JSON conversion used throughout the mRNA Target Selection workflow.

General Functionality:
    - Ensures output and temporary directories exist and are writable.
    - Locates and validates external tool executables (SignalP, TargetP, TMHMM, IEDB tools, etc.).
    - Provides allele panel presets and selection logic for MHC-I and MHC-II.
    - Converts FASTA files to text and JSON formats for downstream processing.
    - Cleans up temporary directories and files.
    - Includes helper functions for parsing and preparing input/output files.
    - Groups sequences from multiple FASTA files by accession code and merges them.
    - Splits protein FASTA files into peptide FASTA files using a sliding window.
    - Contains unit tests for key parsing and conversion functions.

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
from collections import defaultdict
from typing import List
import subprocess
from Bio import SeqIO

def is_valid_peptide(seq: str, min_length: int = 8) -> bool:
    """
    Checks if a given sequence is a valid peptide sequence.
    A valid peptide:
    - Contains only standard amino acids (ACDEFGHIKLMNPQRSTVWY)
    - Is at least `min_length` amino acids long
    """
    return isinstance(seq, str) and len(seq) >= min_length and bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq))

def parse_csv_to_fasta(csv_file: Path, output_dir: Path, basename_prefix: str, min_length=8) -> Path:
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
    in_bepipred_peptide_block = False

    def write_block(header: str, block_peptides: list):
        """Write a block of peptides to the FASTA lines.
        Args:
            header (str): Header for the peptide block.
            block_peptides (list): List of peptide sequences to write.
        """
        if not header or not block_peptides:
            return
        header = header.lstrip(">")
        for i, pep in enumerate(block_peptides):
            if not is_valid_peptide(pep, min_length):
                logging.warning(f"⚠️ Peptide skipped (invalid or too short): {pep} in header: {header}")
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
                    if not is_valid_peptide(peptide, min_length):
                        peptide = None
            else:
                if row[0] == "Position" and len(row) >= 5 and row[4] == "Peptide":
                    continue
                elif len(row) >= 5:
                    peptide = row[4].strip()
            if peptide and is_valid_peptide(peptide, min_length):
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

def group_cluster_inputs(fasta_files: List[Path], fasta_inputs_dir: Path) -> dict:
    """
    Groups sequences from multiple FASTA files by accession code and merges them.

    Args:
        fasta_files (List[Path]): List of input FASTA files.
        fasta_inputs_dir (Path): Output directory for grouped antigen FASTAs.

    Returns:
        dict: Mapping of accession code to combined FASTA Path.
    """
    grouped_by_accession = defaultdict(list)
    header_pattern = re.compile(
        r"^>[^|]*\|(?P<accession>[A-Z0-9_.-]+)\|[^|]*\|(?P<strain_acc>[A-Z0-9_.-]+)\|.*$"
    )
    # example header formats:
    # >antigen_77|Q5HDD7|Immunoglobulin-binding|HE681097.1|tpos:773380-773815
    # >Q5HDD7|HE681097.1

    for fasta_file in fasta_files:
        with open(fasta_file, "r") as fh:
            lines = fh.readlines()

        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                header = lines[i].strip()
                sequence = []
                i += 1
                while i < len(lines) and not lines[i].startswith(">"):
                    sequence.append(lines[i].strip())
                    i += 1

                match = header_pattern.match(header)
                if match:
                    accession = match.group("accession")
                    strain_acc = match.group("strain_acc") or ""  # Default to empty if not present
                    header = f">{accession}|{strain_acc}"  # Reformat header to only have the accession and strain_acc
                    grouped_by_accession[accession].append((header, ''.join(sequence)))
                else:
                    raise ValueError(f"Could not parse accession from header: {header}")
            else:
                i += 1

    output_paths = {}
    for accession, records in grouped_by_accession.items():
        output_path = fasta_inputs_dir / f"{accession}_combined.fasta"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as out_f:
            for header, seq in records:
                out_f.write(f"{header}\n")
                out_f.write(f"{seq}\n")
        output_paths[accession] = output_path

    return output_paths

def parse_json_to_fasta(json_file: Path, output_dir: Path, basename_prefix: str, min_length: int = 8) -> Path:
    """
    Parses a JSON file with peptide predictions and writes a FASTA file of unique, valid peptides
    with contextual headers based on the input file name.

    Args:
        json_file (Path): Path to input JSON file.
        output_dir (Path): Directory where output FASTA will be saved.
        basename_prefix (str): Base prefix for output file name.
        min_length (int): Minimum valid peptide length.

    Returns:
        Path: Path to the generated FASTA file, or None if no valid peptides found.
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
            pep = row[peptide_idx] if peptide_idx < len(row) else None
            if is_valid_peptide(pep, min_length):
                peptides.add(pep)
            else:
                logging.warning(f"⚠️ Skipped invalid peptide: '{pep}' in {json_file.name}")

    if not peptides:
        logging.warning(f"⚠️ No valid peptides found in {json_file}")
        return None

    # Extract metadata from filename
    filename = json_file.name
    match = re.match(
        r"antigen_(\d+)_([A-Z0-9]+)_(.+?)_([A-Z0-9]+\.\d+).*_(MHCI|MHCII)\.json$",
        filename
    )

    if not match:
        logging.error(f"❌ Filename pattern not recognized for T-cell JSON: {filename}")
        return None

    antigen_num, acc_num, _, strain_acc, mhc_class = match.groups()
    mhc_class = mhc_class.lower()

    # Write FASTA lines
    fasta_lines = [
        f">antigen{antigen_num}|{acc_num}|{strain_acc}|{mhc_class}|seq{i+1}\n{pep}"
        for i, pep in enumerate(sorted(peptides))
    ]

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

CONDA_ENV_NAME = "external_tools_env"
CONDA_ENV_YML = Path("ext_tools_dependencies.yml")

def create_conda_env_if_needed():
    """Create Conda environment if it doesn't exist."""
    logging.info(f"🔍 Checking for Conda environment '{CONDA_ENV_NAME}'...")
    result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    if CONDA_ENV_NAME not in result.stdout:
        logging.info("📦 Conda environment not found. Creating from YAML...")
        subprocess.run(["conda", "env", "create", "-f", str(CONDA_ENV_YML)], check=True)
    else:
        logging.info("✅ Conda environment already exists.")

def get_pdb_files(pathogen_path, sequence_dir):
    """
    Get a list of structure files (.pdb, .cif, .cif.gz) from the given sequence directory inside the pathogen path.
    Args:
        pathogen_path (Path): Base path to the pathogen.
        sequence_dir (str): Subdirectory under pathogen_path where structure files are located.
    Returns:
        list[Path]: List of structure file paths.
    """
    search_path = pathogen_path / sequence_dir
    extensions = ["*.pdb", "*.cif", "*.cif.gz"]

    structure_files = []
    for ext in extensions:
        structure_files.extend(search_path.glob(ext))

    return structure_files


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
        temp (bool): Whether to create a temporary JSON directory (default: False).
    
    Returns:
        tuple: (output_dir)
    """
    output_dir = pathogen_path / output_subdir

    # Create tool-specific subdirectories under output_dir
    for tool in selected_tools:
        (output_dir / tool.lower()).mkdir(parents=True, exist_ok=True)

    return output_dir

def cleanup_temp(temp_dirs):
    """
    Clean up temporary directories created during processing.
    Args:
        temp_dirs (list): List of temporary Path objects to clean up.
    """
    for temp_dir in temp_dirs:
        if temp_dir.exists():
            try:
                logging.info(f"🗑️ Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
            except Exception as e:
                logging.error(f"❌ Failed to delete {temp_dir}: {e}")
        else:
            logging.debug(f"⚠️ Temporary directory does not exist: {temp_dir}")


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

def check_antigen_tools(tools: list[str], tool_root: Path) -> dict:
    """
    Detects and validates selected antigen tools.
    
    Args:
        tools (list): List of requested tools (e.g., ["SIGNALP", "TARGETP", "CLUSTER", "DEEPLOC", "ALGPRED"])
        tool_root (Path): Root directory where SignalP/TargetP are expected

    Returns:
        dict: Mapping of tool names to either their executable path or placeholder string
    """
    tool_paths = {}

    if "SIGNALP" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")  
        signalp_path = tool_root / "signalp-5.0b" / "bin" / "signalp"
        if not signalp_path.exists():
            raise FileNotFoundError(
                f"SignalP not found at {signalp_path}. "
                "Install via https://services.healthtech.dtu.dk/services/SignalP-5.0/"
            )
        tool_paths["SIGNALP"] = signalp_path

    if "TARGETP" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")
        targetp_path = tool_root / "targetp-2.0" / "bin" / "targetp"
        if not targetp_path.exists():
            raise FileNotFoundError(
                f"TargetP not found at {targetp_path}. "
                "Install via https://services.healthtech.dtu.dk/services/TargetP-2.0/"
            )
        tool_paths["TARGETP"] = targetp_path

    if "CLUSTER" in tools:
        try:
            subprocess.run(["mmseqs", "version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            tool_paths["CLUSTER"] = "mmseqs"  # placeholder
        except Exception:
            raise FileNotFoundError("MMseqs2 not found in PATH. Install via: conda install -c bioconda mmseqs2")

    if "ALGPRED" in tools:
        algpred_path = Path("ext_tools_dependencies.yml")
        if not algpred_path.exists():
            logging.warning(
                f"⚠️ Algpred2 dependencies file not found at {algpred_path}. "
                "Please get the dependency file from GitHub."
            )
        tool_paths["ALGPRED"] = "algpred2"  # placeholder

    if "DEEPLOC" in tools:
        if not tool_root.exists():
            raise FileNotFoundError(f"Tool root directory {tool_root} does not exist.")
        if not tool_root.is_dir():
            raise NotADirectoryError(f"Tool root {tool_root} is not a directory.")  
        deeploc_path = tool_root / "deeplocpro"
        if not deeploc_path.exists():
            raise FileNotFoundError(
                f"Deeplocpro not found at {deeploc_path}. "
                "Install via https://github.com/Jaimomar99/deeplocpro"
            )
        tool_paths["DEEPLOC"] = deeploc_path

    if "IFNEPITOPE2":
        algpred_path = Path("ext_tools_dependencies.yml")
        if not algpred_path.exists():
            logging.warning(
                f"⚠️ Dependencies file not found at {algpred_path}. "
                "Please get the dependency file from GitHub."
            )
        tool_paths["IFNEPITOPE2"] = "ifnepitope2"  # placeholder

    if "DEEPTMHMM":
        # check that pybiolib is installed
        try:
            import biolib
        except ImportError:
            raise ImportError(
                "PyBioLib is not installed. Please install it via: pip install pybiolib"
            )
        tool_paths["DEEPTMHMM"] = "deetmhmm"  # placeholder

    if "MAFFT_RATE4SITE" or "MAFFT" in tools:
        try:
            algpred_path = Path("ext_tools_dependencies.yml")
            if not algpred_path.exists():
                logging.warning(
                    f"⚠️ Dependencies file not found at {algpred_path}. "
                    "Please get the dependency file from GitHub."
                )
            tool_paths["MAFFT_RATE4SITE"] = "mafft"  # placeholder
            tool_paths["MAFFT"] = "mafft"  # placeholder
        except Exception as e:
            raise FileNotFoundError(f"MAFFT_RATE4SITE not found: {e}")
        
    return tool_paths

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
        "MHCII": base / "ng_tc2-0.2.2-beta" / "src" / "tcell_mhcii.py",
        "Ellipro": base / "ElliPro.jar",
        "MixMHC2pred": base / "MixMHC2pred-2.0" / "MixMHC2pred_unix"
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

def write_json(seq_id_line, seq_lines, temp_dir, alleles, peptide_lengths, tool_type, strain_name):
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
    method = "netmhcpan_el" if tool_type.lower() == "mhci" else "netmhciipan_el"
    
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
        "PopCoverage": tool_root / "population_coverage" /"calculate_population_coverage.py",
        "IFNepitope2": tool_root #doesn't need it
    }

    found = {}
    for name, path in tool_map.items():
        if path.exists():
            found[name] = str(path.parent)
        else:
            print(f"❌ {name} tool not found at: {path.parent}")

    return found

def split_protein_fasta_to_peptides(input_fasta, output_dir, peptide_length=15):
    """
    Splits one or more protein FASTA files into peptide FASTA files using a sliding window.

    Args:
        input_fasta (str, Path, or list): Path(s) to protein FASTA file(s).
        output_dir (str or Path): Directory to save peptide FASTA files.
        peptide_length (int): Length of peptides (default=15).

    Returns:
        list[Path]: List of generated peptide FASTA file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Allow input_fasta to be a single path or a list of paths
    if isinstance(input_fasta, (str, Path)):
        input_fasta = [Path(input_fasta)]
    else:
        input_fasta = [Path(f) for f in input_fasta]

    output_files = []

    for fasta_file in input_fasta:
        output_fasta = output_dir / (fasta_file.stem + "_peptides.fasta")
        with open(output_fasta, "w") as out_f:
            for record in SeqIO.parse(fasta_file, "fasta"):
                protein_seq = str(record.seq)
                header = record.id

                for i in range(len(protein_seq) - peptide_length + 1):
                    peptide = protein_seq[i:i + peptide_length]
                    peptide_id = f"{header}_seq{i+1}"
                    out_f.write(f">{peptide_id}\n{peptide}\n")

        print(f"✅ Peptide FASTA written to: {output_fasta}")
        output_files.append(output_fasta)

    return output_files


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

            fasta_path = Path("/tmp/test_json_to_fasta")
            fasta_path.mkdir(parents=True, exist_ok=True)
            fasta_file = parse_json_to_fasta(json_file, fasta_path, "test_output")
            logging.info(f"FASTA file created at: {fasta_file}")

            self.assertIsNotNone(fasta_file)
            self.assertTrue(fasta_path.exists())
            fasta_content = fasta_file.read_text().strip()
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
            if fasta_file.exists():
                fasta_file.unlink() # remove the created FASTA file
            if fasta_path.exists():
                fasta_path.rmdir()  # remove the FASTA directory
            if json_file.exists():
                json_file.unlink() # remove the created JSON file
            if fasta_file.exists():
                fasta_file.unlink()

    
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
1,10,13,TFNKAATTTT,4
2,39,44,YSGAGKSSS,6
3,56,59,ASTAGA,4
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

            self.assertEqual(len(headers), 2)
            self.assertIn("TFNKAATTTT", sequences)
            self.assertIn("YSGAGKSSS", sequences)
            self.assertNotIn("ASTAGA", sequences)
            self.assertNotIn("been", sequences)

            # Cleanup
            fasta_path.unlink()
            csv_file.unlink()
            shutil.rmtree(output_dir)

    
    class GroupClusterInputsTests(unittest.TestCase):
        def setUp(self):
            self.test_dir = tempfile.TemporaryDirectory()
            self.input_dir = Path(self.test_dir.name) / "input"
            self.output_dir = Path(self.test_dir.name) / "output"
            self.input_dir.mkdir(parents=True)
            self.output_dir.mkdir(parents=True)

            self.fasta1 = self.input_dir / "strain1.fasta"
            self.fasta2 = self.input_dir / "strain2.fasta"

            self.fasta1.write_text(
                ">antigen_153|Q6GHG2|Small|HE681097.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTHKKDHHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
            )

            self.fasta2.write_text(
                ">antigen_153|Q6GHG2|Small|HE681098.1|tpos:380962-381050\n"
                "MAISQERKNEIIKEYRVHETDTGSPEVQIAVLTAEINAVNEHLRTTESTHSRRGLLKMVGRRRHLLNYLRSKDIQRYRELIKSLGIRR\n"
                ">antigen_149|Q99QV7|Putative|HE681097.1|tpos:139239-139462\n"
                "MIEFRQVSKTFNKKKQKIHALKDVSFKVNRNDIFGVIGYSGAGKSTLVRLVNHLEAASSGQVLVDGHDITNY\n"
            )

        def tearDown(self):
            self.test_dir.cleanup()

        def test_group_by_accession_combined_output(self):
            result = group_cluster_inputs(
                fasta_files=[self.fasta1, self.fasta2],
                fasta_inputs_dir=self.output_dir
            )

            expected_files = {
                "Q6GHG2": self.output_dir / "Q6GHG2_combined.fasta",
                "Q99QV7": self.output_dir / "Q99QV7_combined.fasta"
            }

            for accession, path in expected_files.items():
                self.assertIn(accession, result)
                self.assertTrue(path.exists())

            # Check Q6GHG2 file for both strain accessions
            with open(expected_files["Q6GHG2"], "r") as f:
                content = f.read()
                self.assertEqual(content.count(">"), 2)
                self.assertIn(">Q6GHG2|HE681097.1", content)
                self.assertIn(">Q6GHG2|HE681098.1", content)

            # Check Q99QV7 file for correct strain accession
            with open(expected_files["Q99QV7"], "r") as f:
                content = f.read()
                self.assertEqual(content.count(">"), 1)
                self.assertIn(">Q99QV7|HE681097.1", content)


    unittest.main()
