import biolib
import logging
from pathlib import Path

# def run_deeplocpro(tool_path, input_file, output_dir, group):
#     """
#     Runner for DeepLocPro using BioLib.

#     Args:
#         tool_path (Path): Not used.
#         input_file (Path): Input FASTA file.
#         output_dir (Path): Output directory.
#         group (str): Group argument for DeepLocPro (any, archaea, positive, negative).
#     """
#     output_dir = Path(output_dir).resolve()
#     input_file = Path(input_file).resolve()
#     plots_dir = output_dir / "plots"
#     output_dir.mkdir(parents=True, exist_ok=True)
#     plots_dir.mkdir(parents=True, exist_ok=True)
    

#     try:
#         logging.info(f"🔬 Running DeepLocPro on {input_file.name} with group: {group}")
#         deeplocpro = biolib.load("KU/DeepLocPro")
#         result = deeplocpro.cli(
#             args=f"-f {input_file} -o {output_dir} -p -d cpu -g {group}"
#         )
#         logging.debug(f"DeepLocPro result: {result}")
#         logging.info(f"✅ DeepLocPro completed: {input_file.name}")
        
#     except Exception as e:
#         logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
#         raise

def run_deeplocpro(tool_path, input_file, output_dir, group):
    deeplocpro = biolib.load('KU/DeepLocPro')
    print(deeplocpro.cli(args='--help')) 
    