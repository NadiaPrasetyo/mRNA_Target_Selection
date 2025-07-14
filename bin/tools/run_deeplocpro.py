import biolib
import logging
import shutil
from pathlib import Path
import os
import importlib.util
import inspect


def patch_deeplocpro_to_cpu_if_needed():
    """
    Monkey-patch DeepLocPro's embed_batch method to force CPU usage,
    but only if it's not already patched.
    """
    import torch

    try:
        # Locate model.py
        deeplocpro_path = None
        for path in Path("/opt/conda/lib/python3.10/site-packages/DeepLocPro").rglob("model.py"):
            deeplocpro_path = path
            break

        if not deeplocpro_path:
            logging.warning("Could not locate DeepLocPro model.py for patching.")
            return

        # Dynamically import the model module
        spec = importlib.util.spec_from_file_location("deeploc_model", deeplocpro_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        model_class = module.DeepLocModel
        original_fn = model_class.embed_batch

        # Check if it's already patched
        source = inspect.getsource(original_fn)
        if ".to(" in source and "cpu" in source:
            logging.debug("DeepLocPro already patched for CPU. Skipping patch.")
            return

        # Define patched method
        def patched_embed_batch(self, sequences):
            toks = self.batch_converter(sequences)[2]
            device = torch.device("cpu")
            self.esm_model = self.esm_model.to(device)
            toks = toks.to(device)

            with torch.no_grad():
                out = self.esm_model(toks, repr_layers=[33], return_contacts=False)["representations"][33]

            return out, toks.ne(self.batch_converter.alphabet.padding_idx)

        model_class.embed_batch = patched_embed_batch
        logging.info("✅ DeepLocPro patched to use CPU.")

    except Exception as e:
        logging.warning(f"⚠️ DeepLocPro patch check/patch failed: {e}")


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
        patch_deeplocpro_to_cpu_if_needed()      # Patch only if needed

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
