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


### Patch 1: Fix to_csv(sep='\n') ###
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

def patch_to_csv_sep(file_path: Path, original_text: str) -> str:
    tree = ast.parse(original_text)
    fixer = CSVSeparatorFixer()
    fixer.visit(tree)

    if fixer.patched:
        logging.info(f"✅ Will patch .to_csv(sep='\\n') to sep=',' in {file_path.name}")
        return astor.to_source(tree)
    else:
        logging.info("ℹ️  No .to_csv(sep='\\n') found — skipping.")
        return original_text


### Patch 2: Fix joblib import ###
def patch_joblib_import(file_path: Path, text: str) -> str:
    if ("from sklearn.externals import joblib" in text) or ("import sklearn.externals.joblib" in text):
        logging.info(f"✅ Will patch joblib import in {file_path.name}")
        return (
            text.replace("from sklearn.externals import joblib", "import joblib")
                .replace("import sklearn.externals.joblib", "import joblib")
        )
    logging.info("ℹ️  No deprecated joblib import found — skipping.")
    return text


### Patch 3: Fix .str.replace('>', '') ###
def patch_str_replace_gt(file_path: Path, text: str) -> str:
    if ".str.replace('>', '')" in text:
        logging.info(f"✅ Will patch .str.replace('>', '') in {file_path.name}")
        return text.replace(".str.replace('>', '')", ".str.replace('>', '', regex=False)")
    logging.info("ℹ️  .str.replace('>', '') not found — skipping.")
    return text


### Patch 4: Fix old RandomForestClassifier pickle path ###
def patch_random_forest_import(file_path: Path, text: str) -> str:
    inject_code = (
        "import sys\n"
        "import types\n"
        "import sklearn.ensemble._forest\n"
        "# Compatibility for old pickled RandomForestClassifier\n"
        "sys.modules['sklearn.ensemble.forest'] = types.ModuleType('sklearn.ensemble.forest')\n"
        "sys.modules['sklearn.ensemble.forest'].RandomForestClassifier = sklearn.ensemble._forest.RandomForestClassifier\n"
    )

    if "RandomForestClassifier" in text and "sklearn.ensemble.forest" not in text:
        # Insert before first joblib.load call
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if "joblib.load" in line:
                lines.insert(i, inject_code)
                logging.info(f"✅ Injected backward compatibility for RandomForestClassifier path.")
                return "\n".join(lines)
        # fallback
        lines.insert(0, inject_code)
        return "\n".join(lines)

    logging.info("ℹ️  No RandomForestClassifier compatibility patch needed.")
    return text


### Ensure Git LFS model is real ###
def ensure_real_model(model_path: Path):
    """
    Ensures that the rf_model file is not a Git LFS pointer. If it is, runs 'git lfs pull'.
    """
    if not model_path.exists():
        raise FileNotFoundError(f"Model file {model_path} not found!")

    content = model_path.read_text(errors='ignore')
    if "git-lfs.github.com" in content and content.strip().startswith("version https://git-lfs.github.com/spec/v1"):
        logging.warning(f"Model file {model_path} appears to be a Git LFS pointer.")
        logging.info("Attempting to fetch actual model via 'git lfs pull'...")

        result = subprocess.run(["git", "lfs", "pull"], cwd=model_path.parent, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error("❌ Failed to fetch model using 'git lfs pull'.")
            logging.error(result.stderr)
            raise RuntimeError("Git LFS pull failed.")
        logging.info("✅ Git LFS pull successful.")
    else:
        logging.info(f"✅ Model file {model_path.name} looks valid — not a Git LFS pointer.")


### Wrapper to coordinate all patches safely ###
def patch_algpred_script(script_path: Path):
    original_text = script_path.read_text()
    patched_text = original_text

    patched_text = patch_to_csv_sep(script_path, patched_text)
    patched_text = patch_joblib_import(script_path, patched_text)
    patched_text = patch_str_replace_gt(script_path, patched_text)
    patched_text = patch_random_forest_import(script_path, patched_text)

    if patched_text != original_text:
        backup_path = script_path.with_suffix(script_path.suffix + ".bak")
        if not backup_path.exists():
            backup_path.write_text(original_text)
            logging.info(f"🛡️  Original script backed up to {backup_path.name}")
        else:
            logging.info("🛡️  Backup already exists — not overwriting.")

        script_path.write_text(patched_text)
        logging.info(f"🛠️  Patched script saved to {script_path.name}")
    else:
        logging.info("✅ Script already patched — nothing changed.")


### Runner ###
def run(tool_path: Path, input_fasta: Path, output_dir: Path):
    """
    Runs AlgPred2.0 on the provided FASTA file.
    """
    tool_path = Path(tool_path)
    script_path = tool_path / "algpred2.py"
    model_path = tool_path / "rf_model"

    if not script_path.exists():
        raise FileNotFoundError(f"{script_path} not found!")
    if not model_path.exists():
        raise FileNotFoundError(f"Model file {model_path} not found!")

    # patch the script to ensure compatibility
    patch_algpred_script(script_path)
    # Ensure real model is present (not LFS pointer)
    ensure_real_model(model_path)

    output_dir = output_dir / "algpred"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "result.csv"

    command = [
        "python3", str(script_path),
        "-i", str(input_fasta),
        "-o", str(output_file),
        "-m", "1",
        "-d", "1"
    ]

    logging.info("🚀 Running AlgPred2.0...")
    try:
        subprocess.run(command, check=True)
        logging.info(f"✅ AlgPred2.0 completed successfully. Output saved to: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ AlgPred2.0 failed with return code {e.returncode}")
        raise