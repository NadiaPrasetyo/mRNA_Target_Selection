import subprocess
from pathlib import Path
import ast
import astor

class CSVSeparatorFixer(ast.NodeTransformer):
    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'to_csv':
            for kw in node.keywords:
                if kw.arg == 'sep' and isinstance(kw.value, ast.Constant) and kw.value.value == '\n':
                    print(f"⚠️ Patching sep='\\n' on line {node.lineno}")
                    kw.value = ast.Constant(value=',')  # Patch to comma
        return self.generic_visit(node)

def patch_to_csv_sep(file_path: Path):
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
    Arguments:
        tool_path (Path): Path to algpred2.py
        input_fasta (Path): Path to the input FASTA file
        output_dir (Path): Base output directory
    """
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
