from pathlib import Path
import argparse
import logging
import subprocess

def patch_error_msg(tool_path: Path):
    """
    Patch the errorMsg.cpp file in the rate4site source to fix a bug.
    This patch modifies the condition to check if _errorOut is not equal to &std::cerr.
    
    Args:
    - tool_path: Path to the directory containing the rate4site source code.
    
    Raises:
    - FileNotFoundError: If the errorMsg.cpp file does not exist at the expected location.
    """
    source_dir = tool_path / "sourceMar09"
    file_path = source_dir / "errorMsg.cpp"
    print(f"🔍 Looking for errorMsg.cpp at: {file_path}")

    if not file_path.exists():
        raise FileNotFoundError(f"Could not find errorMsg.cpp at expected location: {file_path}")

    # Read and patch the source file
    patched_lines = []
    with open(file_path, "r") as f:
        for line in f:
            if "*_errorOut != cerr" in line:
                line = line.replace("*_errorOut != cerr", "_errorOut != &std::cerr")
            patched_lines.append(line)

    # Write back the patched file
    with open(file_path, "w") as f:
        f.writelines(patched_lines)
    logging.info("🔧 Patching errorMsg.cpp to fix bug in Rate4Site...")
    logging.info("✅ Patch applied successfully.")

def patch_some_util(tool_path: Path):
    """
    Patch someUtil.cpp to replace invalid 'ifstream == NULL' checks.
    """
    logging.info("🔧 Patching someUtil.cpp to fix file open check...")
    file_path = tool_path / "sourceMar09" / "someUtil.cpp"
    
    if not file_path.exists():
        raise FileNotFoundError(f"someUtil.cpp not found at: {file_path}")

    patched_lines = []
    modified = False
    with file_path.open("r") as f:
        for line in f:
            if "== NULL" in line and "ifstream" in line:
                new_line = line.replace("== NULL", ".is_open() == false")
                patched_lines.append(new_line)
                modified = True
            elif "f == NULL" in line:
                new_line = line.replace("f == NULL", "!f")
                patched_lines.append(new_line)
                modified = True
            else:
                patched_lines.append(line)

    if modified:
        with file_path.open("w") as f:
            f.writelines(patched_lines)
        logging.info("✅ Patch applied successfully.")
    else:
        logging.info("ℹ️ No changes made to someUtil.cpp (already patched?).")



def run(tool_path: Path, input_fasta: Path, output_dir: Path, batch_size: int):
    """
    Runs Rate4Site using the patched source code.
    
    Args:
    - tool_path: Path to the directory containing the rate4site source code.
    - input_fasta: Path to the input FASTA file.
    - output_dir: Path to the output directory.
    - batch_size: Unused, present for interface compatibility.
    
    Raises:
    - RuntimeError: If the patching fails or if the command execution fails.
    """
    try:
        patch_error_msg(tool_path)
        patch_some_util(tool_path)  # Assume this is your patch for someUtil.cpp

        logging.info("🔨 Recompiling Rate4Site source code with patched files...")
        command = ["make", "-C", str(tool_path / "sourceMar09")]
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error("❌ Failed to recompile Rate4Site source code:")
            logging.error(result.stderr)
            raise RuntimeError("Compilation failed.")

        logging.info("✅ Compilation successful. Ready to run Rate4Site.")
        # Proceed to run Rate4Site...

    except Exception as e:
        raise RuntimeError(f"Failed to run Rate4Site: {e}")

    
def main():
    """
    Main function to execute the script.
    Take input from command line arguments or predefined paths.
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler()])
    parser = argparse.ArgumentParser(description="Run Rate4Site with patched source code.")
    parser.add_argument("--tool-path", required=True, type=Path, help="Path to the directory containing the rate4site source code. e.g. /home/usr/rate4site.3.2.source/")
    parser.add_argument("--input-fasta", required=False, type=Path, help="Path to the input FASTA file.")
    parser.add_argument("--output-dir", required=False, type=Path, help="Path to the output directory.")

    args = parser.parse_args()
    run(args.tool_path, args.input_fasta, args.output_dir, batch_size=0)

if __name__ == "__main__":
    main()
