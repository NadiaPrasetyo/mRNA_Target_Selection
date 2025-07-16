"""
run_ellipro.py
Utility to execute ElliPro conformational B-cell epitope prediction tool on a given structure file (PDB, CIF, or CIF.GZ).

Overview:
    - Supports input: .pdb, .cif, .cif.gz
    - Automatically unzips .gz and converts .cif to .pdb
    - Runs ElliPro on resulting .pdb file
    - Stores output in specified directory

Arguments:
    input_file (str or Path): Structure file to analyze. Can be PDB, CIF, or gzipped CIF.
    tool_root (str or Path): Path containing 'ElliPro.jar' and 'Convert_CIF_to_PDB/cif_2_pdb.py'.
    output_dir (str or Path): Output directory for prediction results.

Requirements:
    - Java must be available on PATH
    - tool_root must include:
        - ElliPro.jar
        - Convert_CIF_to_PDB/cif_2_pdb.py

Author: Nadia
"""

import subprocess
import gzip
import shutil
import logging
from pathlib import Path


def run(input_file, tool_root, output_dir):
    input_file = Path(input_file)
    tool_root = Path(tool_root)
    output_dir = Path(output_dir)/"ellipro"
    output_dir.mkdir(parents=True, exist_ok=True)

    ellipro_jar = tool_root
    cif_converter = tool_root.parent / "Convert_CIF_to_PDB" / "cif_2_PDB.py"

    if not ellipro_jar.exists():
        logging.error(f"❌ ElliPro.jar not found at {ellipro_jar}")
        return
    if not cif_converter.exists():
        logging.error(f"❌ CIF to PDB converter not found at {cif_converter}")
        return

    stem = input_file.stem.replace(".cif", "").replace(".pdb", "")
    working_file = input_file

    # Step 1: Unzip if input is .gz
    if input_file.suffix == ".gz":
        logging.info(f"🔍 Detected gzipped file: {input_file.name}, unzipping")
        unzipped = input_file.with_suffix("")
        try:
            with gzip.open(input_file, 'rb') as f_in, open(unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            working_file = unzipped
            logging.info(f"✅ Unzipped {input_file.name} → {unzipped.name}")
        except Exception as e:
            logging.error(f"❌ Failed to unzip {input_file.name}: {e}")
            return

    # Step 2: Convert CIF to PDB
    if working_file.suffix == ".cif":
        logging.info(f"🔍 Detected CIF file: {working_file.name}, converting to PDB")
        try:
            cmd = [
                "python3", str(cif_converter),
                "--folder", str(working_file.parent),
                "--output", str(output_dir)
            ]
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            logging.info(f"✅ Converted CIF to PDB: {stem}.pdb")
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ CIF to PDB conversion failed: {working_file.name}")
            logging.error(e.stderr)
            return
        working_file = output_dir / f"{stem}.pdb"

    # Step 3: Run ElliPro (no --chains)
    logging.info(f"🔍 Running ElliPro on {working_file.name}")
    output_file = output_dir / f"{stem}.txt"
    try:
        cmd = [
            "java", "-jar", str(ellipro_jar),
            "--input-file", str(working_file),
            "--output", str(output_file)
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"✅ ElliPro completed: {output_file.name}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ ElliPro failed: {working_file.name}")
        logging.error(e.stderr)
