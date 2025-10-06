#!/usr/bin/env python3
"""
run_whole_pipeline.py

Master pipeline script to run the complete mRNA target selection workflow.

This script orchestrates the entire pipeline from IEDB data fetching through final analysis,
running all tools in the correct order with appropriate dependencies, including random sequence
generation and processing for statistical comparison.

Usage:
    python run_whole_pipeline.py <pathogen_dir> <pathogen_name> [options]

Example:
    python run_whole_pipeline.py sars_cov_2 "SARS-CoV-2" --tool-root /opt/bio_tools --threads 8
"""

import argparse
import os
import sys
import subprocess
import logging
from pathlib import Path
import shutil

def setup_logging(verbose=False, log_file=None):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"
    
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(level)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(format_str))
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is provided
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(format_str))
        logger.addHandler(file_handler)

def run_command(cmd, description, failed_steps, success_exit_codes=(0,)):
    """
    Run a shell command, log outcome, and track failures.

    Args:
        cmd (list[str]): Command and arguments.
        description (str): Human-readable step description.
        failed_steps (list[str]): Collector for failed step descriptions.
        success_exit_codes (tuple[int]): Acceptable exit codes.
    """
    logging.info(f"🚀 {description}")
    logging.debug(f"Command: {' '.join(shlex.quote(c) for c in cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=False, # handle errors manually
            text=True,
            capture_output=True
        )

        # Log stdout/stderr if present
        if result.stdout.strip():
            logging.debug("STDOUT:\n" + result.stdout.strip())
        if result.stderr.strip():
            logging.debug("STDERR:\n" + result.stderr.strip())

        # Handle non-zero exit codes
        if result.returncode not in success_exit_codes:
            logging.error(f"❌ {description} failed (exit {result.returncode})")

            # Show a small preview of the stderr for quick diagnosis
            if result.stderr:
                preview = "\n".join(result.stderr.strip().splitlines()[-10:])
                logging.error("Last few stderr lines:\n" + preview)

            failed_steps.append(description)
            return False

        logging.info(f"✅ {description} completed")
        return True

    except FileNotFoundError:
        logging.error(f"❌ Command not found: {cmd[0]}")
        failed_steps.append(description)
        return False
    except Exception as e:
        logging.exception(f"❌ Unexpected error running {description}: {e}")
        failed_steps.append(description)
        return False
    
def check_dependencies():
    """Check if required external tools are available."""
    required_tools = ["mmseqs", "wget"]
    missing = []
    
    for tool in required_tools:
        if not shutil.which(tool):
            missing.append(tool)
    
    if missing:
        logging.error(f"Missing required tools: {', '.join(missing)}")
        return False
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Run the complete mRNA target selection pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("-pd","--pathogen_dir", help="Pathogen directory name (e.g., sars_cov_2)")
    parser.add_argument("-n","--pathogen_name", help="Full pathogen name (e.g., 'SARS-CoV-2')")
    parser.add_argument("-tr","--tool-root", required=True, help="Root directory containing tool executables")
    parser.add_argument("-t","--threads", type=int, default=4, help="Number of threads to use (default: 4)")
    parser.add_argument("-r","--random-genomes", type=int, default=5, help="Number of random genomes to fetch (default: 5)")
    parser.add_argument("--pfam-hmm", help="Path to Pfam-A.hmm file (required for Pfam analysis)")
    parser.add_argument("--skip", nargs="+", choices=[
        "iedb", "compile", "uniprot", "random", "genomes", "align", "pdb", "pfam", "analysis", "epitopes", "features"
    ], help="Skip specific steps (default: none)")
    parser.add_argument("--steps", nargs="+", choices=[
        "iedb", "compile", "uniprot", "random", "genomes", "align", "pdb", "pfam", "analysis", "epitopes", "features"
    ], help="Run only specific steps (default: all)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--dry-run", action="store_true", help="Show commands without executing")
    parser.add_argument("--human-negative", action ="store_true", help="Use human instead of random as a negative set")

    args = parser.parse_args()
    
    # Setup directories and logging
    base_dir = Path("data") / args.pathogen_dir
    base_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = base_dir / "pipeline.log" if args.verbose else None
    setup_logging(args.verbose, log_file)
    
    logging.info(f"Starting pipeline for {args.pathogen_name}")
    logging.info(f"Output directory: {base_dir}")
    
    # Check dependencies
    if not args.dry_run and not check_dependencies():
        sys.exit(1)
    
    # Define all steps
    all_steps = ["iedb", "compile", "random", "genomes", "uniprot", "align", "analysis", "epitopes", "pdb", "pfam", "evaluate", "features"]
    steps_to_run = [step for step in (args.steps if args.steps else all_steps) if not args.skip or step not in args.skip]
    
    success_count = 0           # successful sub-commands
    attempted_count = 0         # attempted sub-commands
    failed_steps = []      # <<< NEW: Track failed steps


    prefix = "human" if args.human_negative else "random"

    try:
        # Step 1: IEDB Fetch
        if "iedb" in steps_to_run:
            if args.skip and "iedb" in args.skip or (base_dir / f"{args.pathogen_dir}_IEDB_antigens.csv").exists():
                logging.info("⏭️  Skipping IEDB fetch - data already exists")
            else:
                cmd = ["python", "bin/IEDB_fetch.py", args.pathogen_dir, f"{args.pathogen_name}"]
                attempted_count += 1
                if args.dry_run:
                    logging.info(f"Would run: {' '.join(cmd)}")
                elif run_command(cmd, "IEDB data fetch", failed_steps):
                    success_count += 1
        
        # Step 2: Compile Antigens
        if "compile" in steps_to_run:
            cmd = ["python", "bin/compile_antigens.py", args.pathogen_dir, f"{args.pathogen_name}"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Compile antigens", failed_steps):
                success_count += 1
        
        # Step 3: Fetch UniProt Sequences
        if "uniprot" in steps_to_run:
            cmd = ["python", "bin/fetch_sequences_Uniprot.py", args.pathogen_dir, f"{args.pathogen_name}"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Fetch UniProt sequences", failed_steps):
                success_count += 1

        # Step 4: Generate Random Sequences
        if "random" in steps_to_run:
            if args.human_negative:
                cmd = ["python", "bin/generate_random_sequences.py", args.pathogen_dir, f"{args.pathogen_name}", "--human"]
                attempted_count += 1
            else:
                cmd = ["python", "bin/generate_random_sequences.py", args.pathogen_dir, f"{args.pathogen_name}"]
                attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"Generate {prefix} sequences", failed_steps):
                success_count += 1
           
        # Step 5: Fetch NCBI Strain Genomes
        if "genomes" in steps_to_run:
            if args.skip and "genomes" in args.skip and (base_dir / "strain_genomes").exists():
                logging.info("⏭️  Skipping genome fetch - data already exists")
            else:
                # chmod the script first
                cmd = ["chmod", "+x", "bin/fetch_NCBI_strain_genome.sh"]
                attempted_count += 1
                if args.dry_run:
                    logging.info(f"Would run: {' '.join(cmd)}")
                elif run_command(cmd, "Make bin/fetch_NCBI_strain_genome.sh executable", failed_steps):
                    success_count += 1

                cmd = ["bash", "bin/fetch_NCBI_strain_genome.sh", "--random", f"{args.pathogen_name}", 
                       "--random-num", str(args.random_genomes), "--threads", str(args.threads), args.pathogen_dir]
                attempted_count += 1
                if args.dry_run:
                    logging.info(f"Would run: {' '.join(cmd)}")
                elif run_command(cmd, "Fetch strain genomes", failed_steps):
                    success_count += 1
        
        # Step 6: Align Antigens with MMseqs2
        if "align" in steps_to_run:
            # Protein alignment (antigens)
            cmd = ["python", "bin/align_antigens_mmseqs.py", args.pathogen_dir, f"{args.pathogen_name}", 
                   "--threads", str(args.threads), "--output-dir", "mmseqs_protein"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Protein alignment with MMseqs2", failed_steps):
                success_count += 1
            
            # Nucleotide alignment (antigens)
            cmd = ["python", "bin/align_antigens_mmseqs.py", args.pathogen_dir, f"{args.pathogen_name}", 
                   "--threads", str(args.threads), "--output-dir", "mmseqs_nucleotide", "--mode", "nucleotide"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Nucleotide alignment with MMseqs2", failed_steps):
                success_count += 1
            
            # Random protein alignment
            cmd = ["python", "bin/align_antigens_mmseqs.py", args.pathogen_dir, f"{prefix}", 
                   "--threads", str(args.threads), "--output-dir", f"{prefix}_mmseqs_protein"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"{prefix} protein alignment with MMseqs2", failed_steps):
                success_count += 1
            
            # Random nucleotide alignment
            cmd = ["python", "bin/align_antigens_mmseqs.py", args.pathogen_dir, f"{prefix}", 
                   "--threads", str(args.threads), "--output-dir", f"{prefix}_mmseqs_nucleotide", "--mode", "nucleotide"]
            attempted_count += 1
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"{prefix} nucleotide alignment with MMseqs2", failed_steps):
                success_count += 1
        
        # Step 7: Fetch PDB Structures
        if "pdb" in steps_to_run:
            # PDB structures for antigens
            cmd = ["python", "bin/fetch_PDB_structure.py", args.pathogen_dir, "mmseqs_protein",
                   "--threads", str(args.threads), "--output-dir", "pdb_sequences"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Fetch PDB structures", failed_steps):
                success_count += 1
            
            # PDB structures for random sequences
            cmd = ["python", "bin/fetch_PDB_structure.py", args.pathogen_dir, f"{prefix}_mmseqs_protein",
                   "--threads", str(args.threads), "--output-dir", f"{prefix}_pdb_sequences"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"Fetch PDB structures for {prefix} sequences", failed_steps):
                success_count += 1
            
            # IEDB epitope prediction on PDB structures (non-default tools) - antigens
            nondefault_tools = ["DSSP", "ProtLearn", "Ellipro"]
            cmd = ["python", "bin/IEDB_epitope.py", args.pathogen_dir, "pdb_sequences",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--tools"] + nondefault_tools + ["--output-dir", "epitope_outputs"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "IEDB epitope prediction (structure-based tools)", failed_steps):
                success_count += 1
            
            # IEDB epitope prediction on random PDB structures (non-default tools)
            cmd = ["python", "bin/IEDB_epitope.py", args.pathogen_dir, f"{prefix}_pdb_sequences",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--tools"] + nondefault_tools + ["--output-dir", f"{prefix}_analysis"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"IEDB epitope prediction on {prefix} structures (structure-based tools)", failed_steps):
                success_count += 1
        
        # Step 8: Fetch Pfam HMMs
        if "pfam" in steps_to_run:
            if not args.pfam_hmm:
                logging.warning("⚠️  Skipping Pfam analysis - no Pfam-A.hmm file provided")
            else:
                pathogen_name_clean = args.pathogen_name.replace(" ", "_").lower()
                
                # Pfam analysis for antigens
                cmd = ["python", "bin/fetch_pfam_hmmer.py", args.pathogen_dir,
                       "--pathogen_name", pathogen_name_clean, "--pfam_hmm", args.pfam_hmm]
                attempted_count += 1
                if args.dry_run:
                    logging.info(f"Would run: {' '.join(cmd)}")
                elif run_command(cmd, "Fetch Pfam HMMs", failed_steps):
                    success_count += 1
                
                # Pfam analysis for random sequences
                cmd = ["python", "bin/fetch_pfam_hmmer.py", args.pathogen_dir,
                       "--pathogen_name", f"{prefix}", "--pfam_hmm", args.pfam_hmm, "--output-dir", f"{prefix}_pfam"]
                attempted_count += 1
                if args.dry_run:
                    logging.info(f"Would run: {' '.join(cmd)}")
                elif run_command(cmd, f"Fetch Pfam HMMs for {prefix} sequences", failed_steps):
                    success_count += 1

        
        # Step 9: Antigen Analysis (protein-based tools)
        if "analysis" in steps_to_run:
            
            # Antigen analysis
            cmd = ["python", "bin/antigen_analysis.py", args.pathogen_dir, "mmseqs_protein",
                   "--tool-root", args.tool_root, "--threads", str(args.threads), 
                   "--output-dir", "epitope_outputs"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Antigen analysis (protein-based)", failed_steps):
                success_count += 1
            
            # Random analysis
            cmd = ["python", "bin/antigen_analysis.py", args.pathogen_dir, f"{prefix}_mmseqs_protein",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--output-dir", f"{prefix}_analysis"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"{prefix} sequence analysis (protein-based)", failed_steps):
                success_count += 1
            
            # dN/dS analysis on nucleotide alignments (antigens)
            cmd = ["python", "bin/antigen_analysis.py", args.pathogen_dir, "mmseqs_nucleotide",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--tools", "DNDS", "--output-dir", "epitope_outputs"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "dN/dS analysis", failed_steps):
                success_count += 1
            
            # dN/dS analysis on random nucleotide alignments
            cmd = ["python", "bin/antigen_analysis.py", args.pathogen_dir, f"{prefix}_mmseqs_nucleotide",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--tools", "DNDS", "--output-dir", f"{prefix}_analysis"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"{prefix} dN/dS analysis", failed_steps):
                success_count += 1
        
        # Step 10: IEDB Epitope Prediction (default tools)
        if "epitopes" in steps_to_run:
            
            # Epitope prediction on antigens
            cmd = ["python", "bin/IEDB_epitope.py", args.pathogen_dir, "mmseqs_protein",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--output-dir", "epitope_outputs"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "IEDB epitope prediction (default tools)", failed_steps):
                success_count += 1
            
            # Epitope prediction on random sequences
            cmd = ["python", "bin/IEDB_epitope.py", args.pathogen_dir, f"{prefix}_mmseqs_protein",
                   "--tool-root", args.tool_root, "--threads", str(args.threads),
                   "--output-dir", f"{prefix}_analysis"]
            attempted_count += 1
            if args.verbose:
                cmd.append("--verbose")
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, f"IEDB epitope prediction on {prefix} sequences (default tools)", failed_steps):
                success_count += 1
        
        # Step 11: Feature Analysis and KS Tests (Final step after all others complete)
        if "features" in steps_to_run:
            if args.human_negative:
                cmd = ["python", "bin/calculate_features_kstest.py", args.pathogen_dir,
                       "--threads", str(args.threads), "--human"]
                attempted_count += 1
            else:
                cmd = ["python", "bin/calculate_features_kstest.py", args.pathogen_dir,
                   "--threads", str(args.threads)]
                attempted_count += 1
            if args.verbose:
                cmd.extend(["--verbose", "--write-raw"])
            if args.dry_run:
                logging.info(f"Would run: {' '.join(cmd)}")
            elif run_command(cmd, "Feature analysis and statistical tests", failed_steps):
                success_count += 1
    
    except KeyboardInterrupt:
        logging.info("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        sys.exit(1)
    
    # Final summary
    logging.info("=" * 60)
    if args.dry_run:
        logging.info("Dry run completed - no commands were executed")
    else:
        logging.info(
            f"Pipeline completed: {success_count}/{attempted_count} sub-commands succeeded"
        )
        if failed_steps:
            logging.warning("⚠️ Some steps failed:")
            for step in failed_steps:
                logging.warning(f"   - {step}")
            logging.warning("Check pipeline.log for full details.")
            sys.exit(1)
        else:
            logging.info("🎉 All steps completed successfully!")
            logging.info(f"Results available in: {base_dir.absolute()}")


if __name__ == "__main__":
    main()