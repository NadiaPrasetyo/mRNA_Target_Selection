import subprocess
from pathlib import Path

def run(json_file, tool_path, output_dir):
    output_dir = Path(output_dir) / "mhcii"
    output_dir.mkdir(parents=True, exist_ok=True)

    out_base = Path(json_file).stem
    output_prefix = output_dir / out_base

    cmd = ["python3", tool_path, "-j", str(json_file), "-o", str(output_prefix), "-f", "json"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        print(f"✅ MHCII done: {json_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MHCII error: {json_file.name}")
        print(e.stderr)
