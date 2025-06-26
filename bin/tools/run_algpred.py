"""
run_algpred.py
Utility to patch and run the AlgPred2.0 allergenicity prediction tool on protein FASTA files.

Overview:
    - Automatically patches AlgPred2.0's output CSV separator if set incorrectly (e.g., sep='\n').
    - Fixes outdated sklearn.externals.joblib import by replacing it with direct joblib import.
    - Only applies patches when needed (avoids unnecessary modifications).
    - Executes AlgPred2.0 on a given input FASTA file using a specified model and output configuration.
    - Organizes results into a dedicated output directory and provides backup of patched scripts.

Functions:
    patch_to_csv_sep(file_path: Path):
        Scans and patches the specified Python script to ensure pandas.DataFrame.to_csv() uses a comma separator.
        Creates a backup of the original script before modification.

    run(tool_path: Path, input_fasta: Path, output_dir: Path):
        Runs AlgPred2.0 on the provided FASTA file.
            tool_path (Path): Path to the AlgPred2.0 directory containing algpred2.py.
            input_fasta (Path): Path to the input protein FASTA file.
            output_dir (Path): Directory where results will be stored.
        - Patches algpred2.py if necessary.
        - Executes AlgPred2.0 with model 1 (AAC + RF) and displays allergen peptides only.
        - Saves results in an 'algpred' subdirectory under the specified output directory.

Requirements:
    - AlgPred2.0 source code (algpred2.py) available in the specified tool_path.
    - Python packages: subprocess, pathlib, ast, astor.

Usage Example:
    python run_algpred.py /path/to/algpred2 /path/to/input.fasta /path/to/output_dir

Outputs:
    - Patched algpred2.py script (with backup).
    - AlgPred2.0 prediction results in the output directory.
    - Console messages indicating patching and execution status.

Author: Nadia
"""
import subprocess
from pathlib import Path
import ast
import astor
import logging

class PatchDetector(ast.NodeVisitor):
    """
    Detects if a patch is needed:
        - Deprecated sklearn.externals.joblib
        - Invalid .to_csv(sep='\\n')
    """
    needs_joblib_patch = False
    needs_sep_patch = False

    def visit_ImportFrom(self, node):
        if node.module == "sklearn.externals":
            for alias in node.names:
                if alias.name == "joblib":
                    self.needs_joblib_patch = True

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'to_csv':
            for kw in node.keywords:
                if kw.arg == 'sep' and isinstance(kw.value, ast.Constant) and kw.value.value == '\n':
                    self.needs_sep_patch = True

class PatchFixer(ast.NodeTransformer):
    """
    Applies both:
        - sklearn.externals.joblib -> import joblib
        - .to_csv(sep='\\n') -> .to_csv(sep=',')
    """
    def visit_ImportFrom(self, node):
        if node.module == "sklearn.externals":
            for alias in node.names:
                if alias.name == "joblib":
                    logging(f"⚠️ Replacing 'from sklearn.externals import joblib' with 'import joblib' (line {node.lineno})")
                    return ast.Import(names=[ast.alias(name="joblib", asname=None)])
        return node

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'to_csv':
            for kw in node.keywords:
                if kw.arg == 'sep' and isinstance(kw.value, ast.Constant) and kw.value.value == '\n':
                    logging(f"⚠️ Replacing sep='\\n' with sep=',' in .to_csv() (line {node.lineno})")
                    kw.value = ast.Constant(value=',')
        return self.generic_visit(node)

def patch_algpred2(file_path: Path):
    """
    Checks for known issues and applies patches if needed.
    """
    source = file_path.read_text()
    tree = ast.parse(source)

    detector = PatchDetector()
    detector.visit(tree)

    if detector.needs_joblib_patch or detector.needs_sep_patch:
        logging("🔧 Patching algpred2.py...")
        patched_tree = PatchFixer().visit(tree)
        patched_code = astor.to_source(patched_tree)

        backup_path = file_path.with_suffix(file_path.suffix + ".bak")
        backup_path.write_text(source)
        file_path.write_text(patched_code)

        logging(f"🛠️ Patched algpred2.py, backup saved as {backup_path.name}")
    else:
        logging("✅ algpred2.py is already patched — no changes made.")

def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    """
    Runs AlgPred2.0 on the provided FASTA file.
    """
    tool_path = Path(tool_path)
    script_path = tool_path / "algpred2.py"
    if not script_path.exists():
        raise FileNotFoundError(f"AlgPred2 script not found: {script_path}")

    try:
        patch_algpred2(script_path)
    except Exception as e:
        logging(f"⚠️ Patch failed: {e}")

    output_subdir = output_dir / "algpred"
    output_subdir.mkdir(parents=True, exist_ok=True)

    output_csv = output_subdir / f"{input_fasta.stem}_result.csv"

    cmd = [
        "python3", str(script_path),
        "-i", str(input_fasta),
        "-o", str(output_csv),
        "-m", "1",
        "-d", "1"
    ]

    logging(f"🚀 Running AlgPred2 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging(f"❌ AlgPred2 failed: {e}")
    else:
        if output_csv.exists():
            logging(f"✅ AlgPred2 output saved: {output_csv.name}")
        else:
            fallback = output_subdir / "outfile.csv"
            if fallback.exists():
                logging(f"⚠️ Output fallback to: {fallback.name}")
            else:
                logging(f"⚠️ No output CSV found.")