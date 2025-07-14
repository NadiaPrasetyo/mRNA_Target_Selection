import biolib
import logging
import shutil
import tempfile
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
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    plots_dir = output_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        logging.info(f"🔬 Running DeepLocPro on {input_file.name} with group: {group}")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            tmp_input = tmpdir / input_file.name
            shutil.copy(input_file, tmp_input)

            deeplocpro = biolib.load("KU/DeepLocPro")

            # Run DeepLocPro inside sandbox; all output goes into 'output/' dir inside sandbox
            result = deeplocpro.cli(args=f"-f {input_file.name} -o output -p -d cpu -g {group}")

            # Retrieve output files from sandbox
            for file_path in result.output_files:
                file_path = Path(file_path)
                if file_path.suffix.lower() in [".png", ".jpg", ".svg"]:
                    shutil.copy(file_path, plots_dir / file_path.name)
                else:
                    shutil.copy(file_path, output_dir / file_path.name)

            logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
            logging.debug(f"STDERR:\n{result.get_stderr().decode()}")
            logging.info(f"✅ DeepLocPro completed: {input_file.name}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
