import biolib
import logging
import shutil
from pathlib import Path
import tempfile

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

        # Create a temp working dir
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tmp_input = tmp / input_file.name
            shutil.copy(input_file, tmp_input)

            deeplocpro = biolib.load("KU/DeepLocPro")

            # Tell DeepLocPro to write output into 'output/' subdir inside the sandbox
            result = deeplocpro.cli(
                args=f"-f {input_file.name} -o output -p -d cpu -g {group}",
                working_dir=tmp
            )

            # Copy back result files
            for f in result.output_files:
                fpath = Path(f)
                if fpath.suffix.lower() in [".png", ".jpg", ".svg"]:
                    shutil.copy(f, plots_dir / fpath.name)
                else:
                    shutil.copy(f, output_dir / fpath.name)

            logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
            logging.debug(f"STDERR:\n{result.get_stderr().decode()}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
