import biolib
import logging
from pathlib import Path

def run_deeplocpro(tool_path, input_file, output_dir, group):
    """
    Runner for DeepLocPro using BioLib.

    Args:
        tool_path (Path): Not used.
        input_file (Path): Input FASTA file.
        output_dir (Path): Output directory.
        group (str): Group argument for DeepLocPro (any, archaea, positive, negative).
    """
    output_dir = Path(output_dir).resolve()
    input_file = Path(input_file).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        logging.info(f"🔬 Running DeepLocPro on {input_file.name} with group: {group}")
        deeplocpro = biolib.load("KU/DeepLocPro")

        result = deeplocpro.cli(
            args=f"-f {input_file.name} -o output -p -d cpu -g {group}",
            input_files=[str(input_file)],
            output_files=["output/"],  # This makes a subdir inside sandbox for output
        )

        logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
        logging.debug(f"STDERR:\n{result.get_stderr().decode()}")
        logging.info(f"✅ DeepLocPro completed: {input_file.name}")

        # Copy files from BioLib's output dir to actual `output_dir`
        result.save_files(output_dir)

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise

