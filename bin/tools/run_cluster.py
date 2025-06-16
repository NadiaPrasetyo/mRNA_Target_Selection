import subprocess
from pathlib import Path

def run(tool_script: Path, input_json: Path, output_dir: Path,
        output_prefix: str = None,
        output_format: str = "tsv"):
    """
    Runs the cluster tool using a preprocessed JSON input.

    Args:
        tool_script (Path): Path to `run_cluster.py` (the clustering script)
        input_json (Path): Pre-converted JSON input (via `common.parse_fasta_to_jsons`)
        output_dir (Path): Directory where results will be saved
        output_prefix (str): Optional output prefix (default: based on input_json name)
        output_format (str): Output format ('tsv' or 'json')
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_prefix or input_json.stem.replace("_input", "")
    output_path = output_dir / output_prefix

    cmd = [
        "python3", str(tool_script),
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
