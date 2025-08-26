"""
run_dssp.py
Convert PDB -> mmCIF (preferably with gemmi), then run mkdssp with --mmcif-dictionary.

Author: Nadia
"""
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import gemmi
from tools import common


def _find_mmcif_dictionary(conda_prefix: Path) -> Path | None:
    """Search common locations for mmcif_pdbx.dic inside a conda env."""
    candidates = [
        conda_prefix / "share/libcifpp/mmcif_pdbx.dic",
        conda_prefix / "share/mmcif_pdbx.dic",
        conda_prefix / "share/dssp/mmcif_pdbx.dic",
        conda_prefix / "lib/libcifpp/mmcif_pdbx.dic",
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def _pdb_to_mmcif_with_gemmi(pdb_path: Path) -> str:
    """Use gemmi to read the PDB and write a standards-compliant mmCIF file."""
    # gemmi.read_structure handles many PDB oddities robustly
    structure = gemmi.read_structure(str(pdb_path))
    tmp_cif = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
    tmp_cif.close()
    # make_mmcif_document() returns a gemmi.cif.Document; write_file writes mmCIF
    doc = structure.make_mmcif_document()
    doc.write_file(tmp_cif.name)
    return tmp_cif.name


def _pdb_to_mmcif_tmp(pdb_path: Path) -> str:
    """
    Convert a PDB to a temporary mmCIF file using gemmi (better and faster)
    Caller is responsible for deleting the returned path.
    """
    try:
        return _pdb_to_mmcif_with_gemmi(pdb_path)
    except Exception as e:
        logging.warning("gemmi conversion failed: %s", e)


def run(input_file, tool_root, output_dir):
    """
    Convert input PDB to mmCIF and run mkdssp on the mmCIF using --mmcif-dictionary.

    Args:
        input_file (str or Path): Path to the input PDB structure file.
        tool_root (str or Path): Unused, kept for API consistency.
        output_dir (str or Path): Directory to save the DSSP output.
    """
    common.create_conda_env_if_needed()

    input_file = Path(input_file)
    output_dir = Path(output_dir) / "dssp"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{input_file.stem}.dssp"

    # find mkdssp executable
    mkdssp_path = shutil.which("mkdssp")
    if mkdssp_path is None:
        logging.error("❌ mkdssp not found on PATH. Please ensure mkdssp is installed and on PATH.")
        return False

    # locate the mmcif dictionary from CONDA_PREFIX
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        logging.error("❌ CONDA_PREFIX not set. Cannot locate DSSP/mmCIF dictionary.")
        return False

    dic_path = _find_mmcif_dictionary(Path(conda_prefix))
    if not dic_path:
        logging.error("❌ Could not find 'mmcif_pdbx.dic' in the conda environment. Checked common locations under CONDA_PREFIX.")
        return False

    tmp_cif_path = None
    try:
        # 1) Convert PDB -> mmCIF (gemmi preferred)
        tmp_cif_path = _pdb_to_mmcif_tmp(input_file)

        # 2) Run mkdssp on mmCIF with explicit dictionary
        cmd = [
            mkdssp_path,
            "--mmcif-dictionary", str(dic_path),
            tmp_cif_path,
            str(output_file),
            "--verbose",
        ]
        logging.debug("Running mkdssp: %s", " ".join(cmd))
        # Capture output for diagnostics; mkdssp often prints diagnostic messages to stderr.
        proc = subprocess.run(cmd, check=True, text=True, capture_output=True)

        if proc.stdout:
            logging.debug("mkdssp stdout: %s", proc.stdout.strip())
        if proc.stderr:
            logging.debug("mkdssp stderr: %s", proc.stderr.strip())

        # optionally verify output file exists and non-empty
        if not Path(output_file).is_file() or Path(output_file).stat().st_size == 0:
            logging.error("❌ mkdssp did not produce a DSSP file (empty or missing): %s", output_file)
            return False

        logging.info(f"✅ DSSP completed: {output_file.name}")
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ mkdssp failed: {input_file.name}")
        if e.stdout:
            logging.error("mkdssp stdout:\n%s", e.stdout.strip())
        if e.stderr:
            logging.error("mkdssp stderr:\n%s", e.stderr.strip())
        logging.exception(e)
        return False

    except Exception as e:
        logging.error(f"❌ DSSP processing failed: {input_file.name}")
        logging.exception(e)
        return False

    finally:
        # Clean up temp mmCIF
        if tmp_cif_path:
            try:
                Path(tmp_cif_path).unlink()
            except Exception:
                pass
