"""
run_algpred.py
Utility to patch and run the AlgPred2.0 allergenicity prediction tool on protein FASTA files.

Patches:
    - Fixes deprecated pandas `.str.replace('>', '')` to `.str.replace('>', '', regex=False)`
    - Fixes `.to_csv(sep='\n')` to `.to_csv(sep=',')`
    - Replaces deprecated sklearn.externals.joblib import with direct joblib import
    - Injects compatibility for old pickled RandomForestClassifier paths
    - Ensures `rf_model` is loaded via absolute path from the script directory

Author: Nadia (refined for structural clarity)
"""


import subprocess
import logging
from pathlib import Path
import ast
import astor

logging.basicConfig(level=logging.INFO)

### Patch 1: Fix .to_csv(sep='\n') ###
class CSVSeparatorFixer(ast.NodeTransformer):
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

def patch_to_csv_sep(script: str) -> tuple[str, bool]:
    tree = ast.parse(script)
    fixer = CSVSeparatorFixer()
    fixer.visit(tree)
    return (astor.to_source(tree), fixer.patched) if fixer.patched else (script, False)

### Patch 2: Fix joblib import ###
def patch_joblib_import(script: str) -> tuple[str, bool]:
    if "sklearn.externals" in script:
        patched = (
            script.replace("from sklearn.externals import joblib", "import joblib")
                  .replace("import sklearn.externals.joblib", "import joblib")
        )
        return patched, True
    return script, False

### Patch 3: Fix .str.replace('>', '') ###
def patch_str_replace(script: str) -> tuple[str, bool]:
    return (script.replace(".str.replace('>', '')", ".str.replace('>', '', regex=False)"), True) \
        if ".str.replace('>', '')" in script else (script, False)

def patch_rf_pickle(script: str) -> tuple[str, bool]:
    injection = (
        "# Patch for backward compatibility with old sklearn model paths\n"
        "import sys\n"
        "import types\n"
        "import sklearn.ensemble._forest\n"
        "import sklearn.tree._tree\n"
        "import sklearn.tree._classes\n"
        "import sklearn.ensemble._gb\n"
        "import sklearn.ensemble._base\n"
        "sys.modules['sklearn.ensemble.forest'] = sklearn.ensemble._forest\n"
        "sys.modules['sklearn.tree.tree'] = sklearn.tree._tree\n"
        "sys.modules['sklearn.tree._tree'] = sklearn.tree._tree\n"
        "sys.modules['sklearn.tree._classes'] = sklearn.tree._classes\n"
        "sys.modules['sklearn.ensemble._gradient_boosting'] = sklearn.ensemble._gb\n"
        "sys.modules['sklearn.ensemble.base'] = sklearn.ensemble._base\n"
    )

    if "sys.modules['sklearn.tree.tree']" in script:
        return script, False  # Already patched

    lines = script.splitlines()
    insert_idx = 0

    for i, line in enumerate(lines):
        if not line.strip().startswith(("import", "from")):
            insert_idx = i
            break

    lines.insert(insert_idx, injection)
    return "\n".join(lines), True


### Patch 5: Ensure rf_model is loaded via full path ###
class RFModelPathFixer(ast.NodeTransformer):
    def __init__(self):
        self.patched = False

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'load':
            if any(isinstance(arg, ast.Name) and arg.id == 'file_name2' for arg in node.args):
                node.args[0] = ast.Call(
                    func=ast.Name(id='str', ctx=ast.Load()),
                    args=[ast.BinOp(
                        left=ast.Attribute(
                            value=ast.Call(
                                func=ast.Name(id='Path', ctx=ast.Load()),
                                args=[ast.Name(id='__file__', ctx=ast.Load())],
                                keywords=[]
                            ),
                            attr='parent',
                            ctx=ast.Load()
                        ),
                        op=ast.Div(),
                        right=ast.Constant(value='rf_model')
                    )],
                    keywords=[]
                )
                self.patched = True
        return self.generic_visit(node)

def patch_rf_model_path(script: str) -> tuple[str, bool]:
    if "'rf_model'" not in script and 'file_name2' not in script:
        return script, False

    tree = ast.parse(script)
    fixer = RFModelPathFixer()
    tree.body.insert(0, ast.ImportFrom(module='pathlib', names=[ast.alias(name='Path')], level=0))
    fixer.visit(tree)
    return (astor.to_source(tree), fixer.patched) if fixer.patched else (script, False)


### Backup, patch, and write the script ###
def patch_script(script_path: Path):
    original = script_path.read_text()
    modified = original

    modified, p1 = patch_to_csv_sep(modified)
    modified, p2 = patch_joblib_import(modified)
    modified, p3 = patch_str_replace(modified)
    modified, p4 = patch_rf_pickle(modified)
    modified, p5 = patch_rf_model_path(modified)  # 👈 new patch

    if modified != original:
        backup_path = script_path.with_suffix(".bak")
        if not backup_path.exists():
            backup_path.write_text(original)
            logging.info(f"🛡️  Backup created: {backup_path}")
        else:
            logging.info("ℹ️  Backup already exists.")
        script_path.write_text(modified)
        logging.info("🛠️  Script patched successfully.")
    else:
        logging.info("✅ Script is already up-to-date. No patching needed.")


### Check if rf_model is a Git LFS pointer ###
def ensure_real_model(model_path: Path):
    if not model_path.exists():
        raise FileNotFoundError(f"Model file {model_path} not found.")

    # Pointer files are typically <1MB and contain LFS metadata
    if model_path.stat().st_size < 1_000_000:
        content = model_path.read_text(errors='ignore')
        if "git-lfs.github.com" in content:
            logging.warning("⚠️  Model appears to be a Git LFS pointer. Attempting `git lfs pull`...")
            result = subprocess.run(["git", "lfs", "pull"], cwd=model_path.parent, capture_output=True, text=True)
            if result.returncode != 0:
                logging.error(result.stderr)
                raise RuntimeError("❌ Git LFS pull failed.")
            logging.info("✅ Model successfully pulled from Git LFS.")
        else:
            logging.info("✅ Model file is small but valid.")
    else:
        logging.info("✅ Model file is present and looks complete.")

### Main runner ###
def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    tool_path = Path(tool_path)
    script_path = tool_path / "algpred2.py"
    model_path = tool_path /  "rf_model"

    # Step 1: Check and backup original script
    if not script_path.with_suffix(".bak").exists():
        logging.info("🔍 Detected original script — creating backup and preparing to patch.")
    else:
        logging.info("🔍 Script has already been patched previously.")

    # Step 2: If original, ensure model is valid
    ensure_real_model(model_path)

    # Step 3: Patch script
    patch_script(script_path)

    # Step 4: Run AlgPred
    output_dir = output_dir / "algpred"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "result.csv"

    cmd = [
        "python3", str(script_path),
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "1",
        "-d", "1"
    ]

    logging.info("🚀 Running AlgPred2.0...")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"✅ Finished. Output saved to: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ AlgPred2.0 failed with return code {e.returncode}")
        raise
