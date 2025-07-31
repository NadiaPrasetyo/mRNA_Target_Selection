import subprocess
import logging
import shutil
from pathlib import Path
from tools import common

def run(input_fasta: Path, input_tree: Path, output_dir: Path):
    if not shutil.which("conda"):
        logging.error("❌ Conda is not available in PATH.")
        raise RuntimeError("Conda is required but not found.")

    # Ensure conda env is ready (assumes this function is defined elsewhere)
    common.create_conda_env_if_needed()

    input_stem = input_fasta.stem

    # Prepare output paths
    output_file = output_dir / f"{input_stem}.out"
    tree_out_file = output_dir / f"{input_stem}.tree"
    unnormalized_rates_file = output_dir / f"{input_stem}.unnormalized"

    # Make sure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Command construction
    cmd = [
        "conda", "run", "-n", common.CONDA_ENV_NAME,
        "rate4site",      # Assuming rate4site is in the conda environment
        "-s", str(input_fasta),
        "-t", str(input_tree),
        "-o", str(output_file),
        "-x", str(tree_out_file),
        "-y", str(unnormalized_rates_file),
        "-ib",                  # Rate inference method flag:
                                # -im = rates are inferred using the maximum likelihood method
                                # -ib = rates are inferred using the empirical Bayes method 
        "-Mw"                   # The following amino-acids models are supported:
                                # DAY (-md), JTT (-mj), REV (-mr), aaJC (-ma), LG (-Ml), WAG (-Mw)
    ]

    logging.info(f"🚀 Running Rate4Site for {input_fasta.name} with tree {input_tree.name}")
    try:
        subprocess.run(cmd, check=True)
        logging.info("✅ Rate4Site completed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Rate4Site failed with error: {e}")
        raise
