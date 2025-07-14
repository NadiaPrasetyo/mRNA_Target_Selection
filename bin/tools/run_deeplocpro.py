import biolib
import logging
import shutil
from pathlib import Path

def run_deeplocpro(tool_path, input_file, output_dir, group):
    """Run DeepLocPro on a given input file and save the results to the specified output directory.
    Args:
        tool_path (str): unused, path to the DeepLocPro tool (not required for biolib).
        input_file (str): Path to the input FASTA file.
        output_dir (str): Directory where the results will be saved.
        group (str): Group name for DeepLocPro, e.g., 'any', 'archaea', 'positive', 'negative'.
    """

    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    plots_dir = output_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        logging.info(f"🔬 Running DeepLocPro on {input_file.name} with group: {group}")
        deeplocpro = biolib.load("KU/DeepLocPro")

        result = deeplocpro.cli(
            args=f"-f {input_file.name} -o output -p -d cpu -g {group}",
            input_files=[str(input_file)]
        )

        result.save_files(output_dir)

        # Move plot images into the plots/ folder
        for f in output_dir.iterdir():
            if f.suffix.lower() in {".png", ".jpg", ".svg"}:
                shutil.move(str(f), plots_dir / f.name)
            else:
                shutil.move(str(f), output_dir / f.name)

        logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
        logging.debug(f"STDERR:\n{result.get_stderr().decode()}")
        logging.info(f"✅ DeepLocPro completed: {input_file.name}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
