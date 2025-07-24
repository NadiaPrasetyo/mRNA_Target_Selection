"""
calculate_features_kstest.py

Command-line tool to extract immunological and sequence features from epitope and random protein sets,
and compare their distributions using the Kolmogorov-Smirnov (KS) test.

Overview:
    - Parses feature outputs (B-cell, MHC I/II, SignalP, TargetP, allergenicity, cluster conservation, population coverage)
      from specified directories for both epitope (positive) and random sets.
    - Aggregates and structures feature data for statistical comparison.
    - Performs KS tests for each feature/subfeature between positive and random sets.
    - Optionally writes raw feature data to disk for further analysis.
    - Outputs a CSV summary of KS test statistics and p-values.

Arguments:
    pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
    --threads (int, optional): Number of parallel workers for feature extraction (default: 1).
    --verbose (flag, optional): If set, enables verbose logging to file.
    --write-raw (flag, optional): If set, writes raw feature data to disk (can be large).

Requirements:
    - Feature output files in expected formats under:
        data/<pathogen_dir>/epitope_outputs/
        data/<pathogen_dir>/random_analysis/
        data/<pathogen_dir>/evaluation_outputs/
        data/<pathogen_dir>/random_evaluation/
    - Python packages: pandas, scipy

Usage Example:
    python calculate_features_kstest.py sars_cov_2 --threads 4 --verbose --write-raw

Outputs:
    data/<pathogen_dir>/ks_test_results.csv         # KS statistics and p-values for each feature
    data/<pathogen_dir>/raw_positive_features/      # (optional) Raw feature CSVs for positive set
    data/<pathogen_dir>/raw_random_features/        # (optional) Raw feature CSVs for random set
    data/<pathogen_dir>/ks_test.log                 # (optional) Verbose log file

Author: Nadia
"""
import os
import json
import csv
import sys
import pandas as pd
from scipy.stats import ks_2samp
import collections
from collections import defaultdict
import logging
import argparse
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------- Utility Functions -----------------------------

def init_logging(verbose=False, pathogen="unknown"):
    """
    Initialize logging configuration.
    If verbose is True, logs will be written to a file in the specified pathogen directory.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if verbose:
        os.makedirs(f"data/{pathogen}", exist_ok=True)
        fh = logging.FileHandler(f"data/{pathogen}/ks_test.log", mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

def sizeof_fmt(num, suffix="B"):
    """
    Convert a number of bytes into a human-readable format with appropriate suffix.
    Examples:
        1024 -> "1.0KB"
        1048576 -> "1.0MB"
        1073741824 -> "1.0GB"
    """
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return f"{num:.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}P{suffix}"

# ----------------------------- Feature Parsers -----------------------------

def parse_bcell_dir(directory):
    """
    Parse B-cell epitope prediction files in the specified directory.
    Expected files are CSVs with various B-cell prediction methods.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing B-cell prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a B-cell feature.
    Each dictionary contains:
        - "feature": "bcell"
        - "subfeature": method name (e.g., "bepipred", "chou-fasman", etc.)
        - "value": numerical score for the epitope
    """
    logging.info(f"Parsing B-cell features in {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".csv")]
        logging.debug(f"Found {len(files)} CSV files in B-cell dir")
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        method = os.path.basename(file).split("_")[-1].replace(".csv", "").lower()
        path = os.path.join(directory, file)
        logging.debug(f"Parsing B-cell file: {file} with method {method}")
        try:
            with open(path) as f:
                reader = csv.reader(f)
                headers = next(reader)
                if method == "bepipred":
                    idx = headers.index("Score")
                    for i, row in enumerate(reader):
                        try:
                            val = float(row[idx])
                            results.append({"feature": "bcell", "subfeature": method, "value": val})
                        except Exception as e:
                            logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")
                else:
                    if method in ["chou-fasman", "emini", "karplus-schulz", "kolaskar-tongaonkar", "parker"]:
                        idx = headers.index("Score")
                        for i, row in enumerate(reader):
                            try:
                                val = float(row[idx])
                                results.append({"feature": "bcell", "subfeature": method, "value": val})
                            except Exception as e:
                                logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing B-cell file {file}: {e}")
    logging.info(f"Completed B-cell parsing with {len(results)} results")
    return results

def parse_mhc_dir(directory):
    """
    Parse MHC epitope prediction files in the specified directory.
    Expected files are JSONs with MHC I/II prediction results.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing MHC prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an MHC feature.
    Each dictionary contains:
        - "feature": "mhci" or "mhcii" based on the directory name
        - "subfeature": "score", "percentile", or "peptide_length"
        - "value": numerical score or length value
    """
    prefix = "mhcii" if "mhcii" in directory else "mhci"
    logging.info(f"Parsing MHC dir {directory} with prefix {prefix}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".json")]
        logging.debug(f"Found {len(files)} JSON files in MHC dir")
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing MHC file: {file}")
        try:
            with open(path) as f:
                data = json.load(f)
                for result in data.get("results", []):
                    if result.get("type") == "peptide_table":
                        cols = result.get("table_columns", [])
                        table = result.get("table_data", [])
                        try:
                            idx_score = cols.index("score")
                            idx_percentile = cols.index("percentile")
                            idx_peptide = cols.index("peptide")
                        except ValueError as e:
                            logging.warning(f"Missing expected columns in {file}: {e}")
                            continue

                        for i, row in enumerate(table):
                            try:
                                results.append({"feature": prefix, "subfeature": "score", "value": float(row[idx_score])})
                                results.append({"feature": prefix, "subfeature": "percentile", "value": float(row[idx_percentile])})
                                results.append({"feature": prefix, "subfeature": "peptide_length", "value": len(row[idx_peptide])})
                            except Exception as e:
                                logging.debug(f"Skipping row {i} in {file} due to error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing MHC file {file}: {e}")
    logging.info(f"Completed MHC parsing with {len(results)} results")
    return results

def parse_signalp_dir(directory):
    """
    Parse SignalP prediction files in the specified directory.
    Expected files are text files with SignalP results.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing SignalP prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a SignalP feature.
    Each dictionary contains:
        - "feature": "signalp"
        - "subfeature": "prob_signalp" or "prob_other"
        - "value": numerical probability value
    """
    logging.info(f"Parsing SignalP dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".txt")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing SignalP file: {file}")
        try:
            with open(path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        try:
                            results.append({"feature": "signalp", "subfeature": "prob_signalp", "value": float(parts[2])})
                            results.append({"feature": "signalp", "subfeature": "prob_other", "value": float(parts[3])})
                        except Exception as e:
                            logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing SignalP file {file}: {e}")
    logging.info(f"Completed SignalP parsing with {len(results)} results")
    return results

def parse_targetp_dir(directory):
    """
    Parse TargetP prediction files in the specified directory.
    Expected files are text files with TargetP results.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing TargetP prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a TargetP feature.
    Each dictionary contains:
        - "feature": "targetp"
        - "subfeature": "prob_noTP", "prob_SP", or "prob_mTP"
        - "value": numerical probability value
    """
    logging.info(f"Parsing TargetP dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".txt")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing TargetP file: {file}")
        try:
            with open(path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        try:
                            results.append({"feature": "targetp", "subfeature": "prob_noTP", "value": float(parts[2])})
                            results.append({"feature": "targetp", "subfeature": "prob_SP", "value": float(parts[3])})
                            results.append({"feature": "targetp", "subfeature": "prob_mTP", "value": float(parts[4])})
                        except Exception as e:
                            logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing TargetP file {file}: {e}")
    logging.info(f"Completed TargetP parsing with {len(results)} results")
    return results

def parse_allergenicity_dir(directory):
    """
    Parse allergenicity prediction files in the specified directory.
    Expected files are CSVs with MERCI, BLAST, and Hybrid scores.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing allergenicity prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an allergenicity feature.
    Each dictionary contains:
        - "feature": "allergenicity"
        - "subfeature": "merci_score", "blast_score", or "hybrid_score"
        - "value": numerical score value
    """
    logging.info(f"Parsing allergenicity dir {directory}")
    results = []

    # Step 1: List all allergenicity CSV files
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".csv")]
        logging.debug(f"Found {len(files)} allergenicity files")
    except Exception as e:
        logging.error(f"Failed listing allergenicity directory {directory}: {e}")
        return []

    # Step 2: Parse each CSV file
    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing allergenicity file {file}")
        try:
            with open(path, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for i, row in enumerate(reader):
                    try:
                        merci_score = float(row["MERCI Score"])
                        blast_score = float(row["BLAST Score"])
                        hybrid_score = float(row["Hybrid Score"])

                        results.extend([
                            {
                                "feature": "allergenicity",
                                "subfeature": "merci_score",
                                "value": merci_score
                            },
                            {
                                "feature": "allergenicity",
                                "subfeature": "blast_score",
                                "value": blast_score
                            },
                            {
                                "feature": "allergenicity",
                                "subfeature": "hybrid_score",
                                "value": hybrid_score
                            }
                        ])
                    except (ValueError, KeyError) as e:
                        logging.debug(f"Skipping malformed row {i} in {file}: {e}")
        except Exception as e:
            logging.error(f"Failed parsing allergenicity file {file}: {e}")

    logging.info(f"Completed allergenicity parsing with {len(results)} results")
    return results

def parse_cluster_dir(directory):
    """
    Parse cluster conservation files in the specified directory.
    Expected files are m8 format files with cluster conservation scores.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing cluster conservation files.
    Returns:
        List of dictionaries, where each dictionary represents a cluster conservation feature.
    Each dictionary contains:
        - "feature": "cluster_conservation"
        - "subfeature": "conservation_score"
        - "value": numerical conservation score
    """
    logging.info(f"Parsing cluster conservation dir {directory}")
    clusters = defaultdict(list)
    unique_strains = set()

    # Step 1: List all cluster files
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".m8")]
        logging.debug(f"Found {len(files)} cluster files")
    except Exception as e:
        logging.error(f"Failed listing cluster directory {directory}: {e}")
        return []

    # Step 2: Parse each cluster file
    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing cluster file {file}")
        try:
            with open(path) as f:
                for i, line in enumerate(f):
                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue
                    query_id = parts[0]
                    member_id = parts[1]
                    unique_strains.add(member_id) #collect unique strains from ALL the clusters
                    try:
                        percent_identity = float(parts[2])
                        clusters[query_id].append(percent_identity)
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing cluster file {file}: {e}")

    results = []
    cluster_scores = []

    # Step 3: Compute conservation scores and store them as subfeatures
    for identities in clusters.values():
        if identities:
            conservation_score = sum(identities) / len(unique_strains)  # Average over unique strains
            cluster_scores.append(conservation_score)
            results.append({
                "feature": "cluster_conservation",
                "subfeature": "conservation_score",  
                "value": conservation_score
            })

    logging.info(f"Completed cluster parsing with {len(results)} results")
    return results

def parse_popcov_dir(directory):
    """
    Parse population coverage files in the specified directory.
    Expected files are text files with population coverage data.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing population coverage files.
    Returns:
        List of dictionaries, where each dictionary represents a population coverage feature.
    Each dictionary contains:
        - "feature": "popcov"
        - "subfeature": "coverage_average"
        - "value": numerical average coverage value
    """
    logging.info(f"Parsing population coverage dir {directory}")
    results = []

    # Step 1: List all popcov files
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".txt")]
        logging.debug(f"Found {len(files)} popcov files")
    except Exception as e:
        logging.error(f"Failed listing popcov directory {directory}: {e}")
        return []

    # Step 2: Parse each popcov file
    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing popcov file {file}")
        try:
            with open(path) as f:
                lines = f.readlines()
        except Exception as e:
            logging.error(f"Failed reading popcov file {file}: {e}")
            continue

        # Step 3: Extract "average" coverage from the first table
        data_started = False
        for i, line in enumerate(lines):
            stripped = line.strip()

            # Skip empty or header lines
            if not stripped or stripped.lower().startswith("class combined"):
                continue

            if not data_started:
                if stripped.startswith("population/area") and "coverage" in stripped:
                    data_started = True
                continue

            # We're now in the first table
            if data_started:
                parts = stripped.split('\t')
                if len(parts) < 4:
                    continue

                region = parts[0].strip()
                coverage_str = parts[1].strip()

                if region.lower() == "average":
                    try:
                        coverage_val = float(coverage_str.strip('%'))
                        results.append({
                            "feature": "popcov",
                            "subfeature": "coverage_average",
                            "value": coverage_val
                        })
                    except ValueError as e:
                        logging.debug(f"Failed to parse coverage value in {file}, line {i}: {e}")
                    break  # only need the "average" row from the first table

    logging.info(f"Completed popcov parsing with {len(results)} results")
    return results

def parse_deeplocpro_dir(directory):
    """
    Parse DeeplocPro prediction files in the specified directory.
    Expected files are text files with DeeplocPro results.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing DeeplocPro prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a DeeplocPro feature.
    Each dictionary contains:
        - "feature": "deeplocpro"
        - "subfeature": "prob_location1", "prob_location2", etc.
        - "value": numerical probability value
    """
    logging.info(f"Parsing DeeplocPro dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".csv")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing DeeplocPro file: {file}")
        try:
            
                df = pd.read_csv(path)
                # Assume columns: ...,"Cell wall/surface","Extracellular","Cytoplasmic","Cytoplasmic membrane","Outer membrane","Periplasmic"
                prob_cols = [
                    ("cell_wall_surface", "Cell wall/surface"),
                    ("extracellular", "Extracellular"),
                    ("cytoplasmic", "Cytoplasmic"),
                    ("cytoplasmic_membrane", "Cytoplasmic membrane"),
                    ("outer_membrane", "Outer membrane"),
                    ("periplasmic", "Periplasmic"),
                ]
                for i, row in df.iterrows():
                    probs = []
                    for loc, col in prob_cols:
                        try:
                            prob = float(row[col])
                        except Exception:
                            prob = None
                        if prob is not None:
                            results.append({
                                "feature": "deeplocpro",
                                "subfeature": f"prob_{loc}",
                                "value": prob
                            })
                            probs.append(prob)
                    # Add max probability as "max"
                    if probs:
                        results.append({
                            "feature": "deeplocpro",
                            "subfeature": "prob_max",
                            "value": max(probs)
                        })
        except Exception as e:
            logging.error(f"Failed parsing DeeplocPro file {file}: {e}")
    logging.info(f"Completed DeeplocPro parsing with {len(results)} results")
    return results


def parse_ellipro_dir(directory):
    """
    Parse Ellipro prediction files in the specified directory.
    Expected files are text files with linear and discontinuous epitope scores.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing Ellipro prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an Ellipro feature.
    Each dictionary contains:
        - "feature": "ellipro"
        - "subfeature": "linear_epitope_score" or "discontinuous_epitope_score"
        - "value": numerical score value
    """
    logging.info(f"Parsing Ellipro dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".txt")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results
    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing Ellipro file: {file}")
        try:
            with open(path) as f:
                lines = f.readlines()
                # Parse linear epitopes
                for i, line in enumerate(lines):
                    if line.startswith("No.,Structure,Chain,Start Position,End Position,Peptide,Number of Residues,Score,Type"):
                        continue  # Skip header
                    parts = line.strip().split(',')
                    if len(parts) < 8:
                        continue  # Skip malformed lines
                    try:
                        score = float(parts[7])
                        results.append({
                            "feature": "ellipro",
                            "subfeature": "linear_epitope_score",
                            "value": score
                        })
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

                # Parse discontinuous epitopes
                for i, line in enumerate(lines):
                    if line.startswith("No.,Structure,Residues,Number of Residues,Score,Type"):
                        continue  # Skip header
                    parts = line.strip().split(',')
                    if len(parts) < 6:
                        continue  # Skip malformed lines
                    try:
                        score = float(parts[4])
                        results.append({
                            "feature": "ellipro",
                            "subfeature": "discontinuous_epitope_score",
                            "value": score
                        })
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing Ellipro file {file}: {e}")
    logging.info(f"Completed Ellipro parsing with {len(results)} results")
    return results


def parse_ifnepitope2_dir(directory):
    """
    Parse IFNepitope2 prediction files in the specified directory.
    Expected files are CSVs with ML, BLAST, and total scores.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing IFNepitope2 prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an IFNepitope2 feature.
    Each dictionary contains:
        - "feature": "ifnepitope2"
        - "subfeature": "ml_score", "blast_score", or "total_score"
        - "value": numerical score value
    """
    logging.info(f"Parsing IFNepitope2 dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".csv")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing IFNepitope2 file: {file}")
        try:
            with open(path) as f:
                reader = csv.DictReader(f)
                for i, row in enumerate(reader):
                    try:
                        ml_score = float(row["ML_Score"])
                        blast_score = float(row["BLAST_Score"])
                        total_score = float(row["Total_Score"])
                        results.append({
                            "feature": "ifnepitope2",
                            "subfeature": "ml_score",
                            "value": ml_score
                        })
                        results.append({
                            "feature": "ifnepitope2",
                            "subfeature": "blast_score",
                            "value": blast_score
                        })
                        results.append({
                            "feature": "ifnepitope2",
                            "subfeature": "total_score",
                            "value": total_score
                        })
                    except ValueError as e:
                        logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing IFNepitope2 file {file}: {e}")
    logging.info(f"Completed IFNepitope2 parsing with {len(results)} results")
    return results


def parse_mixmhc2pred_dir(directory):
    """
    Parse MixMHC2Pred prediction files in the specified directory.
    Expected files are text files with MixMHC2Pred results.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing MixMHC2Pred prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a MixMHC2Pred feature.
    Each dictionary contains:
        - "feature": "mixmhc2pred"
        - "subfeature": "rank_best" or "best_allele"
        - "value": numerical rank value or allele name
    """
    logging.info(f"Parsing MixMHC2Pred dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".txt")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing MixMHC2Pred file: {file}")
        try:
            with open(path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#") or not line.strip() or line.startswith("Peptide"):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 8:
                        continue
                    try:
                        rank_best = parts[7]  # Assuming "%Rank_best" is the 8th column (index 7)
                        results.append({
                            "feature": "mixmhc2pred",
                            "subfeature": "rank_best",
                            "value": float(rank_best)
                        })
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing MixMHC2Pred file {file}: {e}")
    logging.info(f"Completed MixMHC2Pred parsing with {len(results)} results")
    return results


# ----------------------------- Orchestration Functions -----------------------------

def extract_all_features(base_dir, eval_dir, threads=1):
    """
    Extract all features from the specified base and evaluation directories using multiple threads.
    Args:
        base_dir (str): Path to the base directory containing epitope outputs.
        eval_dir (str): Path to the evaluation directory containing random analysis outputs.
        threads (int): Number of threads to use for parallel parsing.
    Returns:
        List of dictionaries, where each dictionary represents a feature.
    Each dictionary contains:
        - "feature": feature type (e.g., "bcell", "mhci", "mhcii", "signalp", "targetp", "allergenicity", "cluster", "popcov")
        - "subfeature": specific subfeature name (e.g., "bepipred", "score", "prob_signalp", etc.)
        - "value": numerical value for the feature
    """
    logging.info(f"Extracting features from base_dir: {base_dir} and eval_dir: {eval_dir} using {threads} threads")
    parsers = {
        "bcell": lambda: parse_bcell_dir(os.path.join(base_dir, "bcell")),
        "mhci": lambda: parse_mhc_dir(os.path.join(base_dir, "mhci")),
        "mhcii": lambda: parse_mhc_dir(os.path.join(base_dir, "mhcii")),
        "signalp": lambda: parse_signalp_dir(os.path.join(base_dir, "signalp")),
        "targetp": lambda: parse_targetp_dir(os.path.join(base_dir, "targetp")),
        "allergenicity": lambda: parse_allergenicity_dir(os.path.join(base_dir, "algpred")),
        "cluster": lambda: parse_cluster_dir(os.path.join(eval_dir, "cluster")),
        "popcov": lambda: parse_popcov_dir(os.path.join(eval_dir, "popcoverage")),
        "deeplocpro": lambda: parse_deeplocpro_dir(os.path.join(base_dir, "deeplocpro")),
        "ellipro": lambda: parse_ellipro_dir(os.path.join(base_dir, "ellipro")),
        "ifnepitope2": lambda: parse_ifnepitope2_dir(os.path.join(base_dir, "ifnepitope2")),
        "mixmhc2pred": lambda: parse_mixmhc2pred_dir(os.path.join(base_dir, "mixmhc2pred")),
    }

    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_name = {executor.submit(parser): name for name, parser in parsers.items()}
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                result = future.result()
                if result is None:
                    logging.warning(f"{name} parser returned None; skipping")
                    continue
                logging.debug(f"{name} returned {len(result)} items")
                results.extend(result)
                logging.info(f"{name} feature extraction complete")
            except Exception as e:
                logging.error(f"{name} feature extraction failed: {e}")

    return results


def calculate_auroc(pos_vals, rand_vals):
    """
    Compute AUROC from positive and random value lists.
    Returns AUROC or None if computation fails.
    """
    try:
        y_true = [1] * len(pos_vals) + [0] * len(rand_vals)
        y_scores = pos_vals + rand_vals

        # Check for numeric values
        if all(isinstance(v, (int, float)) for v in y_scores):
            return roc_auc_score(y_true, y_scores)
        else:
            return None
    except Exception as e:
        logging.debug(f"AUROC calculation failed: {e}")
        return None
    
def plot_roc_curve(pos_vals, rand_vals, feature, subfeature, output_dir):
    """
    Plot ROC curve and save as PNG.
    """
    output_dir = os.path.join(output_dir, "roc_plots")
    logging.info(f"Plotting ROC curve for {feature}/{subfeature} in {output_dir}")
    try:
        y_true = [1] * len(pos_vals) + [0] * len(rand_vals)
        y_scores = pos_vals + rand_vals

        # Only plot if numeric
        if not all(isinstance(v, (int, float)) for v in y_scores):
            return

        fpr, tpr, _ = roc_curve(y_true, y_scores)
        auc = roc_auc_score(y_true, y_scores)

        # Make output directory if needed
        os.makedirs(output_dir, exist_ok=True)
        filename = f"{feature}_{subfeature}_roc.png".replace("/", "_")
        filepath = os.path.join(output_dir, filename)

        # Plot
        plt.figure()
        plt.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
        plt.plot([0, 1], [0, 1], "k--", label="Random Classifier")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"ROC Curve: {feature} / {subfeature}")
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()

    except Exception as e:
        logging.debug(f"Failed to plot ROC curve for {feature}/{subfeature}: {e}")

def plot_auroc_summary(results_df, output_dir):
    """
    Create a bar plot of AUROC values for all features, corrected and sorted.
    AUROCs < 0.5 are adjusted as 1 - AUROC.
    """
    output_path = os.path.join(output_dir, "auroc_summary.png")
    os.makedirs(output_dir, exist_ok=True)
    try:
        # Drop missing AUROCs
        df = results_df.dropna(subset=["auroc"]).copy()

        # Adjust AUROCs < 0.5
        df["adjusted_auroc"] = df["auroc"].apply(lambda x: x if x >= 0.5 else 1 - x)

        # Combine feature + subfeature for labeling
        df["label"] = df["feature"] + " / " + df["subfeature"]

        # Sort
        df = df.sort_values("adjusted_auroc", ascending=False)

        # Plot
        plt.figure(figsize=(10, max(4, 0.3 * len(df))))
        bars = plt.barh(df["label"], df["adjusted_auroc"], color="skyblue")
        plt.xlabel("AUROC (adjusted, min = 0.5)")
        plt.title("AUROC Summary (Sorted High to Low)")
        plt.xlim(0.5, 1.0)
        plt.gca().invert_yaxis()  # Highest on top

        # Optional: Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.01, bar.get_y() + bar.get_height()/2,
                     f"{width:.3f}", va="center")

        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()

        logging.info(f"AUROC summary plot saved to {output_path}")

    except Exception as e:
        logging.error(f"Failed to generate AUROC summary plot: {e}")


def compare_ks(pos_features, rand_features, output_dir):
    """
    Compare distributions using KS test and AUROC.
    Returns a DataFrame with all results.
    """
    logging.info("Starting KS and AUROC comparison")
    results = []

    def group_values(data):
        grouped = defaultdict(list)
        for row in data:
            if "feature" in row and "subfeature" in row and "value" in row:
                key = (row["feature"], row["subfeature"])
                grouped[key].append(row["value"])
        return grouped

    pos_grouped = group_values(pos_features)
    rand_grouped = group_values(rand_features)

    all_keys = set(pos_grouped.keys()).union(rand_grouped.keys())

    for (feature, subfeature) in all_keys:
        pos_vals = pos_grouped.get((feature, subfeature), [])
        rand_vals = rand_grouped.get((feature, subfeature), [])

        ks_stat, pval, auroc = None, None, None

        if pos_vals and rand_vals:
            try:
                ks_stat, pval = ks_2samp(pos_vals, rand_vals)
                auroc = calculate_auroc(pos_vals, rand_vals)

                # Optional: Generate ROC curve plot if AUROC is computable
                if auroc is not None:
                    plot_roc_curve(pos_vals, rand_vals, feature, subfeature, output_dir)

            except Exception as e:
                logging.debug(f"KS/AUROC failed for {feature}/{subfeature}: {e}")

        results.append({
            "feature": feature,
            "subfeature": subfeature,
            "ks_statistic": ks_stat,
            "p_value": pval,
            "auroc": auroc,
            "positive_n": len(pos_vals),
            "random_n": len(rand_vals)
        })

    logging.info("KS and AUROC comparison complete")
    return pd.DataFrame(results)


def write_features_by_feature(features, label, output_dir):
    """
    Write features to disk, grouped by "feature" field.
    Each feature will be written to a separate CSV file named <feature>_<label>_raw_data.csv
    Args:
        features (list): List of dictionaries representing features.
        label (str): Label for the features (e.g., "positive", "random").
        output_dir (str): Directory to write the feature files to.
    """
    logging.info(f"Writing features to disk for label {label} in {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    # Group features by "feature" field
    grouped = collections.defaultdict(list)
    for row in features:
        grouped[row["feature"]].append(row)

    for feature, rows in grouped.items():
        filepath = os.path.join(output_dir, f"{feature}_{label}_raw_data.csv")

        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["label", "feature", "subfeature", "value"])
            writer.writeheader()

            for row in rows:
                writer.writerow({
                    "label": label,
                    "feature": row.get("feature", ""),
                    "subfeature": row.get("subfeature", ""),
                    "value": row.get("value", "")
                })
        logging.debug(f"Wrote feature file: {filepath}")

# ----------------------------- Entry Point -----------------------------

def main(pathogen_dir, threads, verbose=False, write_raw=False):
    """
    Main entry point for the script.
    Initializes logging, extracts features, performs KS tests, and writes results.
    Args:
        pathogen_dir (str): Pathogen directory name under data/
        threads (int): Optional -number of threads to use for feature extraction.
        verbose (bool): Optional -if True, enables verbose logging to file.
        write_raw (bool): Optional -if True, writes raw feature data to disk (can be large).
    """
    init_logging(verbose, pathogen_dir)
    logger = logging.getLogger()

    logger.info(f"Starting processing for pathogen: {pathogen_dir}")

    pos_dir = os.path.join("data", pathogen_dir, "epitope_outputs")
    rand_dir = os.path.join("data", pathogen_dir, "random_analysis")
    pos_eval_dir = os.path.join("data", pathogen_dir, "evaluation_outputs")
    rand_eval_dir = os.path.join("data", pathogen_dir, "random_evaluation")

    logger.info(f"Extracting positive features from {pos_dir}")
    pos_features = extract_all_features(pos_dir, pos_eval_dir, threads)

    logger.info(f"Extracting random features from {rand_dir}")
    rand_features = extract_all_features(rand_dir, rand_eval_dir, threads)

    logger.info("Estimating memory usage of extracted features")
    logger.info(f"Positive features: {sizeof_fmt(sys.getsizeof(pos_features))}")
    logger.info(f"Random features: {sizeof_fmt(sys.getsizeof(rand_features))}")
    output_dir = os.path.join("results", pathogen_dir)

    raw_out_dir = os.path.join("results", pathogen_dir, "raw_data")

    if write_raw:
        logger.info(f"Writing positive features to {raw_out_dir}")
        write_features_by_feature(pos_features, "positive", raw_out_dir)

        logger.info(f"Writing random features to {raw_out_dir}")
        write_features_by_feature(rand_features, "random", raw_out_dir)


    logger.info("Running KS test on features")
    result_df = compare_ks(pos_features, rand_features, output_dir)
    plot_auroc_summary(result_df, output_dir)

    # Sort the DataFrame alphabetically by the first column
    result_df = result_df.sort_values(by=result_df.columns[0])
    logger.info("\n" + result_df.to_string(index=False))

    ks_out_path = os.path.join("results", pathogen_dir, "ks_test_results.csv")
    logger.info(f"Writing KS test results to {ks_out_path}")
    result_df.to_csv(ks_out_path, index=False)

    logger.info("Processing complete.")

if __name__ == "__main__":
    """ Entry point for command-line execution.
    Parses command-line arguments and calls the main function.
    Usage:
        python calculate_features_kstest.py <pathogen_dir> [--threads <num_threads>]
    Example:
        python calculate_features_kstest.py sars_cov_2 --threads 4 --verbose --write-raw
    Arguments:
        pathogen_dir (str): Subdirectory under `data/` containing pathogen data.
        --threads (int, optional): Number of threads to use for feature extraction (default: 1).
        --verbose (flag, optional): If set, enables verbose logging to file.
        --write-raw (flag, optional): If set, writes raw feature data to disk (can be large).
    Outputs:
        - data/<pathogen_dir>/ks_test_results.csv         # KS statistics and p-values for each feature
        - data/<pathogen_dir>/raw_positive_features/      # (optional) Raw feature CSVs for positive set
        - data/<pathogen_dir>/raw_random_features/        # (optional) Raw feature CSVs for random set
        - data/<pathogen_dir>/ks_test.log                 # (optional) Verbose log file
    """
    parser = argparse.ArgumentParser(description="KS-test comparison of epitope vs. random features.")
    parser.add_argument("pathogen_dir", help="Pathogen directory name under data/")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for parsing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to file")
    parser.add_argument("--write-raw", action="store_true", help="Write raw feature data to disk (optional, large output)")
    args = parser.parse_args()

    main(args.pathogen_dir, args.threads, args.verbose, args.write_raw)
