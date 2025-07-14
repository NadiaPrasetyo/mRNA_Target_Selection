import biolib
import logging
import shutil
from pathlib import Path
import os

def force_embed_batch_cpu(deeplocpro):
    try:
        model = deeplocpro.model  # assuming deeplocpro has a 'model' attribute
        original_embed_batch = model.embed_batch

        def patched_embed_batch(self, sequences):
            sequences = [s.to("cpu") for s in sequences]
            toks = self.batch_converter(sequences)[2].to("cpu")
            out = self.esm_model(toks, repr_layers=[33], return_contacts=False)
            return out["representations"][33].to("cpu"), toks.to("cpu")

        # Replace method on the model instance
        import types
        model.embed_batch = types.MethodType(patched_embed_batch, model)
        print("✅ Patched embed_batch to force CPU usage.")
    except Exception as e:
        print(f"⚠️ Failed to patch embed_batch: {e}")



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

        deeplocpro = biolib.load("KU/DeepLocPro")
        print(dir(deeplocpro))

        force_embed_batch_cpu(deeplocpro)


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
