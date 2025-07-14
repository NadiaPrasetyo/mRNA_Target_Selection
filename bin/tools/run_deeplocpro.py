import biolib
import logging
import shutil
from pathlib import Path
import os

def force_model_to_cpu(monkeypatch_result):
    import torch

    # Locate and patch the embed_batch method at runtime
    try:
        model_mod = __import__("DeepLocPro.model", fromlist=["Model"])
        original_embed_batch = model_mod.Model.embed_batch

        def patched_embed_batch(self, sequences):
            sequences = [s.to("cpu") for s in sequences]
            toks = self.batch_converter(sequences)[2].to("cpu")
            out = self.esm_model(toks, repr_layers=[33], return_contacts=False)
            return out["representations"][33].to("cpu"), toks.to("cpu")

        model_mod.Model.embed_batch = patched_embed_batch
        logging.info("✅ Patched embed_batch to force CPU usage.")
    except Exception as e:
        logging.warning(f"⚠️ Failed to patch embed_batch: {e}")


def run_deeplocpro(tool_path, input_file, output_dir, group):
    """
    Run DeepLocPro on a given input file and save the results to the specified output directory.
    """
    input_file = Path(input_file).resolve()
    output_dir = Path(output_dir).resolve()
    plots_dir = output_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""  # Hide GPUs
        force_model_to_cpu()      # Patch only if needed

        deeplocpro = biolib.load("KU/DeepLocPro")

        result = deeplocpro.cli(
            args=f"-f {input_file} -o output -p -d cpu -g {group}"
        )

        result.save_files(output_dir)

        for f in output_dir.iterdir():
            if f.suffix.lower() in {".png", ".jpg", ".svg"}:
                shutil.move(str(f), plots_dir / f.name)
            else:
                shutil.move(str(f), output_dir / f.name)

        logging.debug(f"STDOUT:\n{result.get_stdout().decode()}")
        logging.debug(f"STDERR:\n{result.get_stderr().decode()}")

    except Exception as e:
        logging.error(f"❌ DeepLocPro failed on {input_file.name}: {e}")
        raise
