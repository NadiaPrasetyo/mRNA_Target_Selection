# tools/run_mixmhc2pred.py
from pathlib import Path
import subprocess

MHCII_DEFAULT = [
    "DRB1_03_01", "DRB1_07_01", "DRB1_15_01", "DRB3_01_01",
    "DRB3_02_02", "DRB4_01_01", "DRB5_01_01"
]

MHCII_EXTENDED = [
    "DRB1_01_01", "DRB1_03_01", "DRB1_04_01", "DRB1_04_05", "DRB1_07_01", "DRB1_08_02",
    "DRB1_09_01", "DRB1_11_01", "DRB1_12_01", "DRB1_13_02", "DRB1_15_01", "DRB3_01_01",
    "DRB3_02_02", "DRB4_01_01", "DRB5_01_01", "DQA1_05_01__DQB1_02_01", "DQA1_05_01__DQB1_03_01",
    "DQA1_03_01__DQB1_03_02", "DQA1_04_01__DQB1_04_02", "DQA1_01_01__DQB1_05_01", "DQA1_01_02__DQB1_06_02",
    "DPA1_02_01__DPB1_01_01", "DPA1_01_03__DPB1_02_01", "DPA1_01_03__DPB1_04_01", "DPA1_03_01__DPB1_04_02",
    "DPA1_02_01__DPB1_05_01", "DPA1_02_01__DPB1_14_01"
]

def run(input_file, tool_path, output_dir, alleles=MHCII_DEFAULT, no_context=True):
    input_file = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{input_file.stem}_mixmhc2pred_output.txt"

    cmd = [
        str(tool_path),
        "-i", str(input_file),
        "-o", str(output_file),
        "-a", *alleles
    ]

    if no_context:
        cmd.append("--no_context")

    print(f"Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Output: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MixMHC2pred failed: {e}")
