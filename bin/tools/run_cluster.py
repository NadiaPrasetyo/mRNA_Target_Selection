"""
run_cluster.py
Utility to run a clustering tool on preprocessed JSON input files.

Overview:
    - Executes the clustering script (`run_cluster.py`) with specified input and output parameters.
    - Handles output directory creation and command construction.
    - Reports success or failure of the clustering process.

Arguments:
    tool_path (Path): Path to the directory containing `run_cluster.py`.
    input_json (Path): Preprocessed JSON input file.
    output_dir (Path): Directory where clustering results will be saved.
    output_prefix (str, optional): Prefix for output files (default: derived from input_json name).
    output_format (str, optional): Output format ('tsv' or 'json', default: 'tsv').

Requirements:
    - Python packages: subprocess, pathlib.
    - The clustering script (`run_cluster.py`) must exist at the specified tool_path.

Usage Example:
    python run_cluster.py /path/to/tools results/input.json results/output --output_prefix clustered --output_format json

Outputs:
    Prints status messages indicating the progress and result of the clustering operation.
    Writes clustering results to the specified output directory.

Author: Nadia
"""
import subprocess
from pathlib import Path

def run(tool_path: Path, input_json: Path, output_dir: Path,
        output_prefix: str = None,
        output_format: str = "tsv"):
    """
    Runs the cluster tool using a preprocessed JSON input.

    Args:
        tool_path (Path): Path to `run_cluster.py` (the clustering script)
        input_json (Path): Pre-converted JSON input (via `common.parse_fasta_to_jsons`)
        output_dir (Path): Directory where results will be saved
        output_prefix (str): Optional output prefix (default: based on input_json name)
        output_format (str): Output format ('tsv' or 'json')
    """
    tool_path = Path(tool_path)  # convert string to Path object
    script_path = tool_path/"run_cluster.py"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix or input_json.stem.replace("_input", "")
    output_path = output_dir / output_prefix

    cmd = [
        "python3", str(script_path),
        "-j", str(input_json),
        "-o", str(output_path),
        "-f", output_format
    ]

    print(f"🚀 Running cluster: {input_json.name}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Clustering failed: {e}")
    else:
        print(f"✅ Output written: {output_path}.{output_format}")
