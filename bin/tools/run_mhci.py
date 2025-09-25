"""
run_mhci.py
Utility to execute NetMHCpan-4.2 epitope prediction tool on a given FASTA input file.

Overview:
    - Runs the external NetMHCpan-4.2 tool using subprocess.
    - Uses a predefined extended set of MHCI alleles.
    - Organizes output files into a dedicated 'mhci' subdirectory within the specified output directory.
    - Handles process execution and error reporting.

Arguments:
    fasta_file (str or Path): Path to the input FASTA file containing peptide/protein sequences.
    tool_path (str or Path): Path to the NetMHCpan-4.2 executable.
    output_dir (str or Path): Directory where output files will be stored. Results are placed in a 'mhci' subfolder.

Requirements:
    - NetMHCpan-4.2 must be installed and executable at `tool_path`.
    - Python packages: subprocess, pathlib.

Usage Example:
    run("input/protein.fasta",
        "/projects/.../netMHCpan-4.2/Linux_x86_64/bin/netMHCpan-4.2",
        "results/epitopes")

Outputs:
    - Writes tool output to `output_dir/mhci/<basename>.out`
    - Writes Excel output to `output_dir/mhci/<basename>.xls`
    - Prints success or error messages accordingly.

Author: Nadia
"""
import subprocess
from pathlib import Path
import logging

MHCI_EXTENDED = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
    "HLA-A*26:01", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01", "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01",
    "HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"
]

def run(fasta_file, tool_path, output_dir):
    # Default parameters
    rank_threshold = 0.5
    weak_threshold = 2.0
    allele_string = ",".join(MHCI_EXTENDED)

    fasta_file = Path(fasta_file)
    tool_path = Path(tool_path)
    output_dir = Path(output_dir) / "mhci"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{fasta_file.stem}_mhci_out"

    cmd = [
        str(tool_path),
        "-a", allele_string,
        "-f", str(fasta_file),
        "-rankS", str(rank_threshold),
        "-rankW", str(weak_threshold),
        "-xls",
        "-xlsfile", str(output_file.with_suffix(".xls")),
    ]

    logging.info(f"Running NetMHCpan-4.2 on {fasta_file} with {len(MHCI_EXTENDED)} alleles...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        # Save stdout into .out file
        output_file.write_text(result.stdout)

        logging.info(f"✅ NetMHCpan-4.2 completed successfully.")
        logging.info(f"Results saved in: {output_file} and {output_file.with_suffix('.xls')}")
    except subprocess.CalledProcessError as e:
        logging.info("❌ NetMHCpan-4.2 failed.")
        logging.info("Command:", " ".join(cmd))
        logging.info("Stdout:", e.stdout)
        logging.info("Stderr:", e.stderr)
        raise
