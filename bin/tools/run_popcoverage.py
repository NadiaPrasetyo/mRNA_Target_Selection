"""
run_popcoverage.py
Command-line utility to run the population coverage analysis tool on epitope-allele input files.

Overview:
    - Executes the calculate_population_coverage.py script for a given set of epitopes and HLA alleles.
    - Supports analysis for specific populations and MHC classes (I, II, or combined).
    - Optionally generates population coverage plots.
    - Saves the population coverage results to a specified output directory.

Arguments:
    tool_path (Path): Path to the directory containing calculate_population_coverage.py.
    input_file (Path): Path to the input .txt file with epitope and allele information.
    output_dir (Path): Directory to save the population coverage results.
    population (str, optional): Population name to analyze (default: 'Global').
    mhc_class (str, optional): MHC class to analyze ('I', 'II', or 'combined'; default: 'combined').
    plot (bool, optional): Whether to generate a population coverage plot (default: False).

Requirements:
    - calculate_population_coverage.py script available in the specified tool_path.
    - Python packages: subprocess, pathlib.

Usage Example:
    python run_popcoverage.py /path/to/tools input_epitopes.txt results/ --population EastAsia --mhc_class I --plot

Outputs:
    - Writes the population coverage analysis results to a text file in the output directory.
    - Optionally generates and saves a plot of the population coverage.

Author: Nadia
"""
import subprocess
from pathlib import Path
import shutil
import logging


def run(tool_path: Path, input_file: Path, output_dir: Path,
        population="World", mhc_class="combined", plot=False):
    """
    Run population coverage tool on a given input file.

    Args:
        tool_path (Path): Path to the directory containing calculate_population_coverage.py
        input_file (Path): Path to the input .txt file (epitope, allele format)
        output_dir (Path): Directory to save the results
        population (str): Population name (default: 'World')
        mhc_class (str): MHC class ('I', 'II', or 'combined')
        plot (bool): Whether to generate a plot (default: False)
    """
    tool_path = Path(tool_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    script_path = tool_path / "calculate_population_coverage.py"
    output_prefix = input_file.stem

    # Output file
    output_txt = output_dir / f"{output_prefix}.txt"

    # Build base command
    cmd = [
        "python", str(script_path),
        "-p", population,
        "-c", mhc_class,
        "-f", str(input_file)
    ]

    # Plot output handling
    if plot:
        plot_path = output_dir / "plot" / output_prefix
        plot_path.parent.mkdir(parents=True, exist_ok=True)  # ensure plot/ exists
        cmd += ["--plot", str(plot_path)]  # Script will add .png internally

    try:
        print(f"📊 Running PopCoverage: {input_file.name} for {population} (MHC-{mhc_class})")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(output_txt, "w") as f:
            f.write(result.stdout)

        print(f"✅ PopCoverage results saved: {output_txt.name}")

        # If plots were generated, flatten them into plot/ directory
        if plot:
            plot_dir = output_dir / "plot"
            for subdir in list(plot_dir.iterdir()):
                if subdir.is_dir():
                    for png_file in subdir.glob("*.png"):
                        newname = f"{subdir.name}_{png_file.name}"
                        dest = plot_dir / newname
                        shutil.move(str(png_file), str(dest))
                    try:
                        subdir.rmdir()
                    except Exception as e:
                        logging.warning(f"⚠️ Failed to remove empty directory {subdir}: {e}")

    except subprocess.CalledProcessError as e:
        print(f"❌ PopCoverage failed for {input_file.name}")
        if e.stderr:
            print(e.stderr)