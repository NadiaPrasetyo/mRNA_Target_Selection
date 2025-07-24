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
import gemmi

def convert_cif_to_pdb_with_gemmi(input_cif, output_pdb):
    """Convert CIF file to PDB format using Gemmi.
    Args:
        input_cif (Path): Path to the input CIF file.
        output_pdb (Path): Path to save the output PDB file.
    Returns:
        output_pdb (Path): Path to the generated PDB file.
    """
    try:
        doc = gemmi.cif.read_file(str(input_cif))
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)

        # Truncate 2-letter chain IDs to 1-letter
        for model in structure:
            for chain in model:
                if len(chain.name) > 1:
                    chain.name = chain.name[0]

        structure.write_pdb(str(output_pdb))
        logging.info(f"✅ Gemmi converted: {input_cif.name} → {output_pdb.name}")
        return output_pdb
    except Exception as e:
        logging.error(f"❌ Gemmi CIF→PDB failed: {e}")
        return None

def run(input_file, tool_root, output_dir):
    """
    Run ElliPro on the specified input file after necessary conversions.
    Args:
        input_file (str or Path): Path to the input structure file (PDB, CIF, or CIF.GZ).
        tool_root (str or Path): Path to the ElliPro tool directory containing ElliPro.jar.
        output_dir (str or Path): Directory to save the ElliPro output.
    """
    input_file = Path(input_file)
    tool_root = Path(tool_root)
    output_dir = Path(output_dir)/"ellipro"
    output_dir.mkdir(parents=True, exist_ok=True)

    ellipro_jar = tool_root

    if not ellipro_jar.exists():
        logging.error(f"❌ ElliPro.jar not found at {ellipro_jar}")
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
    
    # Step 2: Convert CIF to PDB (using Gemmi)
    if working_file.suffix == ".cif":
        logging.info(f"🔍 Detected CIF file: {working_file.name}, converting to PDB with Gemmi")
        pdb_output = output_dir / f"{stem}.pdb"
        working_file = convert_cif_to_pdb_with_gemmi(working_file, pdb_output)

        if not working_file or not working_file.exists():
            logging.error(f"❌ Failed to generate PDB from CIF: {pdb_output}")
            return


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
        
