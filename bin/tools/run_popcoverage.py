import subprocess
from pathlib import Path

def run(tool_path, input_file, output_dir, population="Global", mhc_class="combined", plot=False):
    """
    Run population coverage tool on a given input file.

    Parameters:
    - tool_path (Path): Path to the directory containing calculate_population_coverage.py
    - input_file (Path): Path to the input .txt file (epitope, allele format)
    - output_dir (Path): Directory to save the results
    - population (str): Population name (default: 'Global')
    - mhc_class (str): MHC class ('I', 'II', or 'combined')
    - plot (bool): Whether to generate a plot (default: False)
    """
    script_path = tool_path / "calculate_population_coverage.py"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = input_file.stem
    output_file = output_dir / f"{output_prefix}.txt"

    # Build command
    cmd = [
        "python", str(script_path),
        "-p", population,
        "-c", mhc_class,
        "-f", str(input_file)
    ]

    if plot:
        plot_path = output_dir / f"{output_prefix}_plot"
        cmd += ["--plot", str(plot_path)]

    try:
        print(f"📊 Running PopCoverage: {input_file.name} for {population} (MHC-{mhc_class})")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Save stdout to output file
        with open(output_file, "w") as f:
            f.write(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"❌ PopCoverage failed for {input_file.name}")
        print(e.stderr)
