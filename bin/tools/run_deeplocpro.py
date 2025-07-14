import biolib
import logging
import shutil
import tempfile
from pathlib import Path

def run_deeplocpro(tool_path, input_file, output_dir, group):
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
            result = deeplocpro.cli(args=f"-f {input_file.name} -o output -p -d cpu -g {group}", working_dir=str(tmpdir))

            # Save all sandbox output into a temporary directory
            sandbox_output = tmpdir / "output_files"
            sandbox_output.mkdir()
            result.save_files(sandbox_output)

            # Move files: plots to /plots, others to /output_dir
            for file in sandbox_output.iterdir():
                if file.suffix.lower() in [".png", ".jpg", ".svg"]:
                    shutil.move(str(file), plots_dir / file.name)
                else:
                    shutil.move(str(file), output_dir / file.name)

            logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
            logging.debug(f"STDERR:\n{result.get_stderr().decode()}")
            logging.info(f"✅ DeepLocPro completed: {input_file.name}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
