import subprocess
from pathlib import Path

def run(json_file, tool_path, output_dir):
    out_base = Path(json_file).stem
    output_prefix = Path(output_dir) / out_base
    cmd = ["python3", tool_path, "-j", str(json_file), "-o", str(output_prefix), "-f", "json"]
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        print(f"✅ MHCI done: {json_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MHCI error: {json_file.name}")
        print(e.stderr)
