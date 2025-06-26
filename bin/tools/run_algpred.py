"""
run_algpred.py
Utility to patch and run the AlgPred2.0 allergenicity prediction tool on protein FASTA files.

Patches:
    - Fixes deprecated pandas `.str.replace('>', '')` to `.str.replace('>', '', regex=False)`
    - Ensures `.to_csv(sep=',')` (fixes accidental newline separators)
    - Replaces deprecated sklearn.externals.joblib import with direct joblib import

Author: Nadia
"""

import subprocess
from pathlib import Path
import ast
import astor
import logging

# Configure logging (if not configured externally)
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

### Patch 1: Fix to_csv(sep='\n') ###
class CSVSeparatorFixer(ast.NodeTransformer):
    """
    AST transformer to patch pandas DataFrame.to_csv() calls with sep='\\n' to use sep=','.
    """
    def __init__(self):
        self.patched = False

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'to_csv':
            for kw in node.keywords:
                if kw.arg == 'sep' and isinstance(kw.value, ast.Constant) and kw.value.value == '\n':
                    logging.warning(f"Patching sep='\\n' on line {node.lineno}")
                    kw.value = ast.Constant(value=',')
                    self.patched = True
        return self.generic_visit(node)

def patch_to_csv_sep(file_path: Path):
    """
    Patches the specified Python script to ensure pandas DataFrame.to_csv() uses a comma separator.
    """
    source = file_path.read_text()
    tree = ast.parse(source)
    fixer = CSVSeparatorFixer()
    fixer.visit(tree)

    if fixer.patched:
        patched_code = astor.to_source(tree)
        backup_path = file_path.with_suffix(file_path.suffix + ".bak")
        backup_path.write_text(source)
        file_path.write_text(patched_code)
        logging.info(f"Patched .to_csv() in {file_path.name}, backup saved at {backup_path.name}")
    else:
        logging.info(".to_csv() sep is fine — no patch needed.")

### Patch 2: Fix sklearn.externals.joblib ###
class JoblibFixDetector(ast.NodeVisitor):
    """
    Detects if 'from sklearn.externals import joblib' is present.
    """
    def __init__(self):
        self.needs_patch = False

    def visit_ImportFrom(self, node):
        if node.module == "sklearn.externals":
            for alias in node.names:
                if alias.name == "joblib":
                    self.needs_patch = True
        self.generic_visit(node)

class JoblibFixer(ast.NodeTransformer):
    """
    Replaces 'from sklearn.externals import joblib' with 'import joblib'.
    """
    def visit_ImportFrom(self, node):
        if node.module == "sklearn.externals":
            for alias in node.names:
                if alias.name == "joblib":
                    logging.warning(f"Replacing sklearn.externals.joblib with import joblib (line {node.lineno})")
                    return ast.Import(names=[ast.alias(name="joblib", asname=None)])
        return node

def patch_joblib_import(file_path: Path):
    """
    Patches algpred2.py if sklearn.externals.joblib is used.
    """
    source = file_path.read_text()
    tree = ast.parse(source)

    detector = JoblibFixDetector()
    detector.visit(tree)

    if detector.needs_patch:
        tree = JoblibFixer().visit(tree)
        patched_code = astor.to_source(tree)
        backup_path = file_path.with_suffix(file_path.suffix + ".bak")
        backup_path.write_text(source)
        file_path.write_text(patched_code)
        logging.info(f"Patched joblib import in {file_path.name}, backup saved at {backup_path.name}")
    else:
        logging.info("joblib import is already correct — no patch needed.")

### Patch 3: Fix .str.replace('>', '') -> .str.replace('>', '', regex=False) ###
def patch_str_replace_gt(file_path: Path):
    original = file_path.read_text()
    if ".str.replace('>', '')" in original:
        patched = original.replace(".str.replace('>', '')", ".str.replace('>', '', regex=False)")
        backup_path = file_path.with_suffix(file_path.suffix + ".bak")
        backup_path.write_text(original)
        file_path.write_text(patched)
        logging.info("Patched .str.replace('>', '') to .str.replace('>', '', regex=False)")
    else:
        logging.info(".str.replace('>', '') not found — no patch needed.")

### Wrapper function ###
def patch_algpred_script(script_path: Path):
    patch_to_csv_sep(script_path)
    patch_joblib_import(script_path)
    patch_str_replace_gt(script_path)

### Runner ###
def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    """
    Runs AlgPred2.0 on the provided FASTA file.
    """
    tool_path = Path(tool_path)
    script_path = tool_path / "algpred2.py"
    if not script_path.exists():
        raise FileNotFoundError(f"{script_path} not found!")

    patch_algpred_script(script_path)

    output_dir = output_dir / "algpred"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "result.csv"
    command = [
        "python3", str(script_path),
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "1",  # Model 1 (AAC + RF)
        "-d", "1"   # Only show Allergen predictions
    ]

    logging.info("Running AlgPred2.0...")
    subprocess.run(command, check=True)
    logging.info(f"AlgPred2.0 finished. Output saved to: {output_file}")

