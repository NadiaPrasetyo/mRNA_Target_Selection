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


def patch_source_code(tool_path: Path):
    source_file = tool_path / "sourceMar09/someUtil.cpp"
    if not source_file.exists():
        raise FileNotFoundError(f"Could not find someUtil.cpp at {source_file}")
    
    with source_file.open("r") as f:
        lines = f.readlines()

    # Patch the line containing 'if (f == NULL)'
    new_lines = []
    for line in lines:
        if "if (f == NULL)" in line:
            new_line = line.replace("if (f == NULL)", "if (!f.is_open())")
            print(f"Patching line:\n  {line.strip()}\n→ {new_line.strip()}")
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    with source_file.open("w") as f:
        f.writelines(new_lines)
    
    logging.info("🔧 Patching someUtil.cpp to fix file open check...")
    logging.info("✅ Patch applied successfully.")


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
        patch_source_code(tool_path)
        
        # Here you would typically run the command to execute Rate4Site
        # For example:
        # subprocess.run(["rate4site", str(input_fasta), "-o", str(output_dir)], check=True)
    except Exception as e:
        raise RuntimeError(f"Failed to run Rate4Site: {e}")
    
     # recompile the patched source code
    try:
        result = subprocess.run(["make", tool_path/"sourceMar09"], check=True, capture_output=True, text=True)
        # print all output from the make command
        logging.info(result.stdout)
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Failed to recompile Rate4Site source code: {e}")
        raise
    
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
