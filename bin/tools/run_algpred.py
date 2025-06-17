"""
run_algpred.py
Utility to patch and run the AlgPred2.0 allergenicity prediction tool on protein FASTA files.

Overview:
    - Automatically patches AlgPred2.0's output CSV separator if set incorrectly (e.g., sep='\n').
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

class CSVSeparatorFixer(ast.NodeTransformer):
    """
    AST transformer to patch pandas DataFrame.to_csv() calls with sep='\n' to use sep=','.
    This is necessary to ensure proper CSV formatting.
    """
    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'to_csv':
            for kw in node.keywords:
                if kw.arg == 'sep' and isinstance(kw.value, ast.Constant) and kw.value.value == '\n':
                    print(f"⚠️ Patching sep='\\n' on line {node.lineno}")
                    kw.value = ast.Constant(value=',')  # Patch to comma
        return self.generic_visit(node)

def patch_to_csv_sep(file_path: Path):
    """
    Patches the specified Python script to ensure pandas DataFrame.to_csv() uses a comma separator.
    Args:
        file_path (Path): Path to the algpred2.py script to be patched.
    """
    source = file_path.read_text()
    tree = ast.parse(source)
    tree = CSVSeparatorFixer().visit(tree)
    patched_code = astor.to_source(tree)

    backup_path = file_path.with_suffix(file_path.suffix + ".bak")
    backup_path.write_text(source)
    file_path.write_text(patched_code)

    print(f"🛠️ Patched .to_csv() in {file_path.name}, backup saved at {backup_path.name}")

def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    """
    Runs AlgPred2.0 tool on the given FASTA file.
    Args:
        tool_path (Path): Path to algpred2.py
        input_fasta (Path): Path to the input FASTA file
        output_dir (Path): Base output directory
    Raises:
        FileNotFoundError: If algpred2.py script is not found in the specified tool_path.
        subprocess.CalledProcessError: If AlgPred2.0 execution fails.
    """
    tool_path = Path(tool_path)  # convert string to Path object
    script_path = tool_path / "algpred2.py"
    if not script_path.exists():
        raise FileNotFoundError(f"AlgPred2 script not found: {script_path}")

    # Pre-patch script if necessary
    try:
        patch_to_csv_sep(script_path)
    except Exception as e:
        print(f"⚠️ Could not patch CSV sep: {e}")

    output_subdir = output_dir / "algpred"
    output_subdir.mkdir(parents=True, exist_ok=True)
    output_file = output_subdir / f"{input_fasta.stem}_ALGPRED.txt"

    cmd = [
        "python3", str(script_path),
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "1",          # Model 1 (AAC + RF)
        "-d", "1"           # Display allergen peptides only
    ]

    print(f"🚀 Running AlgPred2 on {input_fasta.name}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ AlgPred2 failed: {e}")
    else:
        print(f"✅ AlgPred2 completed: {output_file.name}")
