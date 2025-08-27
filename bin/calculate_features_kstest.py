"""
calculate_features_kstest.py
Command-line tool to extract immunological and sequence features from epitope and random protein sets, and compare their distributions using KS-test and AUROC.

Overview:
    - Extracts features (B-cell, MHC I/II, SignalP, TargetP, allergenicity, cluster conservation, DeeplocPro, Ellipro, IFNepitope2, MixMHC2Pred, DeepTMHMM, Rate4Site, Rate4Site_Mafft_DeepTMHMM) from epitope (positive) and random sets.
    - Aggregates and structures feature data for statistical comparison.
    - Performs Kolmogorov-Smirnov (KS) tests and computes AUROC for each feature/subfeature between positive and random sets.
    - Optionally writes raw feature data to disk for further analysis.
    - Generates ROC curve plots and an AUROC summary bar plot.

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
    - Python packages: pandas, scipy, scikit-learn, matplotlib, os, csv, json, argparse.

Usage Example:
    python calculate_features_kstest.py sars_cov_2 --threads 4 --verbose --write-raw

Outputs:
    results/<pathogen_dir>/ks_test_results.csv      # KS statistics, p-values, and AUROC for each feature
    results/<pathogen_dir>/raw_data/                # (optional) Raw feature CSVs for positive and random sets
    results/<pathogen_dir>/roc_plots/               # ROC curve plots for each feature/subfeature
    results/<pathogen_dir>/auroc_summary.png        # AUROC summary bar plot
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
from math import log
import argparse
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from concurrent.futures import ThreadPoolExecutor, as_completed
import statistics
from Bio import AlignIO
import glob
import re
import traceback
import tempfile
from pathlib import Path

# ----------------------------- Utility Functions -----------------------------

def init_logging(verbose=False, pathogen="unknown"):
    """
    Initialize logging configuration.
    If verbose is True, logs will be written to a file in the specified pathogen directory.
    Args:
        verbose (bool): If True, enable verbose logging to file.
        pathogen (str): Name of the pathogen for log file naming.
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
        os.makedirs(f"results/{pathogen}", exist_ok=True)
        fh = logging.FileHandler(f"results/{pathogen}/ks_test.log", mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

def sizeof_fmt(num, suffix="B"):
    """
    Convert a number of bytes into a human-readable format with appropriate suffix.
    Args:
        num (int): Number of bytes.
        suffix (str): Suffix to append (default: "B").
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

    - add the percent identities and divide by number of strains (simple)
    - sum of bit scores (length confounding) / sum of the length
    - sum the log₁₀(e-value)

    0 	Query sequence ID
    1 	Subject (database) sequence ID
    2 	Percent Identity
    3 	Alignment Length
    4 	Number of gaps
    5 	Number of mismatches
    6 	Start on the query sequence
    7 	End on the query sequence
    8 	Start on the database sequence
    9 	End on the database sequence
    10 	E value - the expectation that this alignment is random given the length of the sequence and length of the database
    11 	bit score - the score of the alignment itself
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
                        length = int(parts[3])
                        e_value = float(parts[10])
                        bit_score = float(parts[11])
                        clusters[query_id].append(percent_identity)
                        clusters[query_id].append(bit_score / length)  # Normalize bit score by length
                        if e_value > 0:
                            clusters[query_id].append(log(e_value))  # Normalize e-value by log base 10
                        else:
                            logging.debug(f"Skipping line {i} in {file} due to non-positive e_value: {e_value}")
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing cluster file {file}: {e}")

    results = []

    # Step 3: Compute conservation scores and store them as subfeatures
    for cluster_id, scores in clusters.items():
        if not scores:
            continue
        try:
            percent_identity_num_strain = sum(scores[0::3]) / len(unique_strains)  #sum of percent identities divided by number of strains
            avg_bit_score = sum(scores[1::3]) / len(scores[1::3]) # sum of bit scores divided by number of scores
            avg_log_e_value = sum(scores[2::3]) / len(scores[2::3]) # sum of e-values divided by number of e-values
            results.append({
                "feature": "cluster_conservation",
                "subfeature": "percent_identity/num_strain",
                "value": percent_identity_num_strain
            })
            results.append({
                "feature": "cluster_conservation",
                "subfeature": "bit_score_normalized",
                "value": avg_bit_score
            })
            results.append({
                "feature": "cluster_conservation",
                "subfeature": "e_value_average",
                "value": avg_log_e_value
            })
        except ZeroDivisionError as e:
            logging.debug(f"Skipping cluster {cluster_id} due to division by zero: {e}")

    logging.info(f"Completed cluster parsing with {len(results)} results")
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
                    parts = line.strip().rsplit(',', 2)  # Split into three parts from the right
                    if len(parts) < 3:
                        continue  # Skip malformed lines
                    try:
                        score = float(parts[-2])  # Take the second-to-last field as the score
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
                    parts = line.strip().rsplit(',', 2)  # Split into three parts from the right
                    if len(parts) < 3:
                        continue  # Skip malformed lines
                    try:
                        score = float(parts[-2])  # Take the second-to-last field as the score
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

def parse_deeptmhmm_dir(directory):
    """
    Parse Deeptmhmm prediction files in the specified directory.
    Expected files are 3-line format files with sequence and topology.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing Deeptmhmm prediction files.
    Returns:
        List of dictionaries, where each dictionary represents a Deeptmhmm feature.
    Each dictionary contains:
        - "feature": "deeptmhmm"
        - "subfeature": "proportion_outside"
        - "value": numerical proportion of outside residues in the topology
    """
    logging.info(f"Parsing Deeptmhmm dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".3line")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results
    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing Deeptmhmm file: {file}")
        try:
            with open(path) as f:
                lines = f.readlines()
                i = 0
                while i < len(lines):
                    # Look for header line
                    if lines[i].startswith(">"):
                        if i + 2 >= len(lines):
                            logging.warning(f"File {file} does not have enough lines after header at line {i}")
                            break
                        sequence = lines[i + 1].strip()
                        topology = lines[i + 2].strip()
                        if len(sequence) != len(topology):
                            logging.warning(f"Sequence and topology lengths do not match in {file} at header {lines[i].strip()}")
                        # Calculate proportion of outside residues
                        outside_count = sum(1 for char in topology if char == 'O')
                        total_count = len(sequence)
                        if total_count > 0:
                            proportion_outside = outside_count / total_count
                        else:
                            proportion_outside = 0.0
                        results.append({
                            "feature": "deeptmhmm",
                            "subfeature": "proportion_outside",
                            "value": proportion_outside
                        })
                        i += 3  # Move to next header or end
                    else:
                        i += 1  # Skip lines until next header
        except Exception as e:
            logging.error(f"Failed parsing Deeptmhmm file {file}: {e}")
    logging.info(f"Completed Deeptmhmm parsing with {len(results)} results")
    return results

def parse_rate4site_dir(directory):
    """
    Parse Rate4Site output files in the specified directory.
    Extracts per-site conservation scores and computes rolling window statistics.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing Rate4Site output files.
    Returns:
        List of dictionaries, where each dictionary represents a Rate4Site feature.
    Each dictionary contains:
        - "feature": "rate4site"
        - "subfeature": specific subfeature name (e.g., "per_site_score", "rolling_mean", etc.)
        - "value": numerical value for the feature
    """
    logging.info(f"Parsing Rate4Site dir {directory}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".out")]
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing Rate4Site file: {file}")
        try:
            scores = []
            with open(path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split()  # Split by whitespace
                    if len(parts) < 4:
                        continue
                    try:
                        score = float(parts[2])
                        scores.append(score)
                        results.append({
                            "feature": "rate4site",
                            "subfeature": "per_site_score",
                            "value": score
                        })
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

            # Compute rolling window statistics (15-peptide window)
            window_size = 15
            for j in range(len(scores)):
                window = scores[max(0, j - window_size + 1):j + 1]
                if len(window) > 0:
                    results.append({
                        "feature": "rate4site",
                        "subfeature": "rolling_mean",
                        "value": sum(window) / len(window)
                    })
                    results.append({
                        "feature": "rate4site",
                        "subfeature": "rolling_median",
                        "value": statistics.median(window)
                    })
                    results.append({
                        "feature": "rate4site",
                        "subfeature": "rolling_max",
                        "value": max(window)
                    })
                    results.append({
                        "feature": "rate4site",
                        "subfeature": "rolling_min",
                        "value": min(window)
                    })
        except Exception as e:
            logging.error(f"Failed parsing Rate4Site file {file}: {e}")
    logging.info(f"Completed Rate4Site parsing with {len(results)} results")
    return results


def parse_rate4site_mafft_deeptmhmm_dir(rate4site_dir, mafft_dir, deeptmhmm_dir):
    """
    Parse the output directories for Rate4Site, MAFFT, and DeepTMHMM results by:
    - Reading Rate4Site output files to get per-site conservation scores.
    - Reading MAFFT alignment files to get the reference sequence and alignment.
    - Reading DeepTMHMM topology files to get strain-antigen topologies.
    Returns a list of dictionaries with feature values.
    Args:
        rate4site_dir (str): Path to the directory containing Rate4Site output files.
        mafft_dir (str): Path to the directory containing MAFFT alignment files.
        deeptmhmm_dir (str): Path to the directory containing DeepTMHMM topology files.
    Returns:
        List of dictionaries, where each dictionary represents a Rate4Site feature.
    Each dictionary contains:
        - "feature": "rate4site_deeptmhmm"
        - "subfeature": "outside_score_per_site" or "rolling_avg" or "rolling_max" or "rolling_min"
        - "strain": strain name
        - "antigen": antigen accession
        - "value": numerical value for the feature
    """
    def parse_rate4site_out(filepath):
        """        
        Parse a Rate4Site output file to extract per-site conservation scores.
        Args:
            filepath (str): Path to the Rate4Site output file.
        Returns:
            List of tuples (position, amino acid, score) for each site.
        """
        scores = []
        with open(filepath) as f:
            for line in f:
                if line.strip().startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    pos, aa, score, msa_data = parts
                    scores.append((int(pos), aa, float(score)))
        logging.debug(f"Parsed {len(scores)} scores from {filepath}")
        return scores

    def load_all_deeptmhmm_topologies(deeptmhmm_dir):
        """
        Load all DeepTMHMM topologies from the specified directory.
        Args:
            deeptmhmm_dir (str): Path to the directory containing DeepTMHMM topology files.
        Returns:
            dict: A mapping of strain names to their antigen topologies.
        """
        topologies = {}  # strain -> antigen -> topology
        for path in glob.glob(os.path.join(deeptmhmm_dir, "*.3line")):
            try:
                with open(path) as f:
                    lines = f.read().splitlines()
                    for i in range(0, len(lines), 3):
                        try:
                            header = lines[i]
                            seq = lines[i + 1]
                            topology = lines[i + 2]
                            fields = header.lstrip(">").split("|")
                            if len(fields) < 4:
                                logging.warning(f"Unexpected header format: {header}")
                                continue
                            antigen_accession = fields[1]
                            strain = fields[3]
                            topologies.setdefault(strain, {})[antigen_accession] = topology
                        except IndexError:
                            logging.warning(f"Malformed 3-line record in {path} at lines {i}-{i+2}")
            except Exception:
                logging.exception(f"Error reading {path}")
        logging.info(f"Loaded topologies for {len(topologies)} strains from {deeptmhmm_dir}")
        return topologies

    def get_reference_sequence_index_map(alignment, accession):
        """
        Create a mapping of reference sequence positions to alignment indices.
        """
        ref_seq = None
        for record in alignment:
            if accession in record.id:
                ref_seq = record.seq
                break
        if ref_seq is None:
            raise ValueError(f"Reference sequence {accession} not found in alignment.")

        mapping = {}
        ref_pos = 0
        for i, res in enumerate(ref_seq):
            if res != '-':
                ref_pos += 1
                mapping[ref_pos] = i
        logging.debug(f"Reference map created for {accession}: {len(mapping)} positions")
        return mapping

    def get_strain_sequence_from_alignment(alignment, strain):
        """
        Get the sequence for a specific strain from the alignment.
        """
        for record in alignment:
            if strain in record.id:
                return str(record.seq)
        return None

    def map_topology_to_alignment(strain_seq, topology):
        """
        Map the topology to the strain sequence alignment.
        """
        mapped = []
        topo_index = 0
        for res in strain_seq:
            if res == "-":
                mapped.append(None)
            elif topo_index < len(topology):
                mapped.append(topology[topo_index])
                topo_index += 1
            else:
                mapped.append(None)
        return mapped

    results = []
    strain_antigen_topos = load_all_deeptmhmm_topologies(deeptmhmm_dir)

    for rate4site_file in glob.glob(os.path.join(rate4site_dir, "*.out")):
        try:
            accession = os.path.basename(rate4site_file).split("_combined")[0]
            mafft_file = os.path.join(mafft_dir, f"{accession}_combined_aligned.fasta")

            if not os.path.exists(mafft_file):
                logging.warning(f"Missing MAFFT file: {mafft_file}")
                continue

            logging.info(f"Processing accession: {accession}")
            scores = parse_rate4site_out(rate4site_file)
            if not scores:
                logging.warning(f"No scores found in {rate4site_file}")
                continue

            alignment = AlignIO.read(mafft_file, "fasta")
            ref_map = get_reference_sequence_index_map(alignment, accession)

            for strain, antigen_map in strain_antigen_topos.items():
                if accession not in antigen_map:
                    continue
                topology = antigen_map[accession]
                strain_seq = get_strain_sequence_from_alignment(alignment, strain)
                if strain_seq is None:
                    logging.debug(f"Strain {strain} not in alignment for {accession}")
                    continue

                topo_by_pos = map_topology_to_alignment(strain_seq, topology)
                strain_topology_positions = []

                for pos, aa, score in scores:
                    aln_index = ref_map.get(pos)
                    if aln_index is None or aln_index >= len(topo_by_pos):
                        continue
                    topo = topo_by_pos[aln_index]
                    if topo is None:
                        continue

                    results.append({
                        "feature": "rate4site_deeptmhmm",
                        "subfeature": "outside_score_per_site",
                        "strain": strain,
                        "antigen": accession,
                        "value": score
                    })

                    strain_topology_positions.append((pos, aa, score, topo))

                # Rolling window stats for this strain-antigen pair
                outside_scores = [s for _, _, s, t in strain_topology_positions if t == "O"]
                window_size = 15
                for i in range(len(outside_scores) - window_size + 1):
                    window = outside_scores[i:i + window_size]
                    results.extend([
                        {
                            "strain": strain,
                            "accession": accession,
                            "feature": "rate4site_deeptmhmm",
                            "subfeature": "rolling_avg",
                            "value": sum(window) / window_size
                        },
                        {
                            "strain": strain,
                            "accession": accession,
                            "feature": "rate4site_deeptmhmm",
                            "subfeature": "rolling_median",
                            "value": statistics.median(window)
                        },
                        {
                            "strain": strain,
                            "accession": accession,
                            "feature": "rate4site_deeptmhmm",
                            "subfeature": "rolling_min",
                            "value": min(window)
                        },
                        {
                            "strain": strain,
                            "accession": accession,
                            "feature": "rate4site_deeptmhmm",
                            "subfeature": "rolling_max",
                            "value": max(window)
                        }
                    ])

        except Exception as e:
            logging.exception(f"Error processing {rate4site_file}: {str(e)}")

    logging.info(f"Total features returned: {len(results)}")
    return results

def parse_dssp_dir(dssp_dir):
    """
    Parse the DSSP directory for structural features.
    data/S.pyogenes/epitope_outputs/dssp/1I5K_P49054.dssp
    --E------------EE----EEE-----------------HHH-------E--------EEEE------EEEE---E----E------------EE----EE--------------HHH-HHH-------E--------EEEE------EEEE---E-----HHHHHHHHHHHHHHHHHHHHH-----HHHHHHHHHHHHHHHHHHHHH-  1I5K_P49054.pdb

    """
    results = []
    for dssp_file in Path(dssp_dir).glob("*.dssp"):
        contents = dssp_file.read_text()
        percent_helix = contents.count('H') / len(contents) if len(contents) > 0 else 0
        percent_sheet = contents.count('E') / len(contents) if len(contents) > 0 else 0
        percent_loop = contents.count('-') / len(contents) if len(contents) > 0 else 0

        results.append({
            "feature": "dssp",
            "subfeature": "percent_helix",
            "value": percent_helix
        })
        results.append({
            "feature": "dssp",
            "subfeature": "percent_sheet",
            "value": percent_sheet
        })
        results.append({
            "feature": "dssp",
            "subfeature": "percent_loop",
            "value": percent_loop
        })

    return results

def parse_dnds_dir(directory):
    """
    Parse the DNDS directory for structural features.

    """
    results = []
    for dnds_file in Path(directory).glob("*.json"):
        # e.g.: A0A0H2UTN5_combined_codon_aligned_FEL_results.json
        type = dnds_file.stem.split("_")[-2] # Extract the type from the filename
        content = dnds_file.read_json().get("MLE").get("content").get("0")
        match type:
            case "FEL":
                for each in content:
                    sum_n += each[1]
                    sum_s += each[0]

                results.append({
                    "feature": "FEL dN/dS",
                    "subfeature": "sumN/sumS",
                    "value": sum_n / sum_s if sum_s > 0 else 0
                })
            case "SLAC":
                content = content.get("by-site").get("AVERAGED")
                for each in content:
                    sum_n += each[3]
                    sum_s += each[2]

                    dn = each[6]
                    ds = each[5]
                    results.append({
                        "feature": "SLAC dN/dS",
                        "subfeature": "dN/dS",
                        "value": dn / ds if ds > 0 else 0
                    })

                results.append({
                    "feature": "SLAC dN/dS",
                    "subfeature": "sumN/sumS",
                    "value": sum_n / sum_s if sum_s > 0 else 0
                })
            case "FUBAR":
                for each in content:
                    dn = each[1]
                    ds = each[0]

                    prob_neg_selection = each[3]
                    prob_pos_selection = each[4]

                    results.append({
                        "feature": "FUBAR dN/dS",
                        "subfeature": "dN/dS",
                        "value": dn / ds if ds > 0 else 0
                    })

                    results.append({
                        "feature": "FUBAR dN/dS",
                        "subfeature": "Probability negative selection",
                        "value": prob_neg_selection
                    })
                    results.append({
                        "feature": "FUBAR dN/dS",
                        "subfeature": "Probability positive selection",
                        "value": prob_pos_selection
                    })
            
            case _:
                logging.warning(f"Unknown dN/dS type: {type}")
                continue

    return results

def parse_protlearn_dir(directory):
    """
    Parse the ProtLearn directory for structural features.
    feature,value
length,[211.]
aac_A,0.04739336492890995
aac_C,0.05687203791469194
aac_D,0.06635071090047394
aac_E,0.10900473933649289
aac_F,0.018957345971563982

    """
    results = []
    for protlearn_file in Path(directory).glob("*.csv"):
        with open(protlearn_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                feature = row["feature"]
                value = row["value"]
                results.append({
                    "feature": "ProtLearn",
                    "subfeature": feature,
                    "value": int(value)
                })
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
        - "feature": feature type (e.g., "bcell", "mhci", "mhcii", "signalp", "targetp", "allergenicity", "cluster", "popcov", "deeplocpro", "ellipro",
        "ifnepitope2", "mixmhc2pred",  "deeptmhmm", "rate4site", "rate4site_mafft_deeptmhmm")
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
        "cluster": lambda: parse_cluster_dir(os.path.join(base_dir, "cluster")),
        "deeplocpro": lambda: parse_deeplocpro_dir(os.path.join(base_dir, "deeplocpro")),
        "ellipro": lambda: parse_ellipro_dir(os.path.join(base_dir, "ellipro")),
        "ifnepitope2": lambda: parse_ifnepitope2_dir(os.path.join(base_dir, "ifnepitope2")),
        "mixmhc2pred": lambda: parse_mixmhc2pred_dir(os.path.join(base_dir, "mixmhc2pred")),
        "deeptmhmm": lambda: parse_deeptmhmm_dir(os.path.join(base_dir, "deeptmhmm")),
        "rate4site": lambda: parse_rate4site_dir(os.path.join(base_dir, "mafft_rate4site/rate4site_results")),
        "rate4site_mafft_deeptmhmm": lambda: parse_rate4site_mafft_deeptmhmm_dir(
            os.path.join(base_dir, "mafft_rate4site/rate4site_results"),
            os.path.join(base_dir, "mafft_rate4site"),
            os.path.join(base_dir, "deeptmhmm")
        ),
        "dnds": lambda: parse_dnds_dir(os.path.join(base_dir, "dnds")),
        "protlearn": lambda: parse_protlearn_dir(os.path.join(base_dir, "protlearn")),
        "dssp": lambda: parse_dssp_dir(os.path.join(base_dir, "dssp"))
    }

    import traceback

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
                logging.debug(traceback.format_exc())

    return results


def calculate_auroc(pos_vals, rand_vals):
    """
    Compute AUROC from positive and random value lists.
    Args:
        pos_vals (list): List of positive values (e.g., epitope scores).
        rand_vals (list): List of random values (e.g., background scores).
    Returns:
        float: AUROC value or None if computation fails.
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
    Args:
        pos_vals (list): List of positive values (e.g., epitope scores).
        rand_vals (list): List of random values (e.g., background scores).
        feature (str): Feature name (e.g., "bcell", "mhci").
        subfeature (str): Subfeature name (e.g., "bepipred", "score").
        output_dir (str): Directory to save the ROC plot.
    This function will create a directory "roc_plots" inside output_dir if it doesn't exist
    and save the ROC plot as a PNG file named "{feature}_{subfeature}_roc.png".
    If the values are not numeric, it will skip plotting.
    If any error occurs during plotting, it will log the error but not raise an exception.
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

def categorize_feature(feature, subfeature):
    """
    Categorize features for coloring in AUROC summary plot.
    Returns a string category.
    """
    # Subcellular localisation
    if feature in ["signalp", "targetp", "deeplocpro", "deeptmhmm"]:
        return "Subcellular localisation"
    # Allergenicity
    if feature == "allergenicity":
        return "Allergenicity"
    # Immunogenicity
    if feature == "ifnepitope2":
        return "Immunogenicity"
    # Conservation Analysis
    if feature in ["cluster_conservation", "rate4site", "rate4site_deeptmhmm", "dnds"] or subfeature in [
        "Percent identity / number of strains",
        "Average Log₁₀ e-value",
        "Average bit-score / length",
        "percent_identity/num_strain",
        "e_value_average",
        "bit_score_normalized",
        "rolling_mean",
        "rolling_median",
        "rolling_max",
        "rolling_min",
        "outside_score_per_site",
        "per_site_score"
    ]:
        return "Conservation Analysis Across Strains"
    # Epitope Prediction
    if feature in ["bcell", "ellipro", "mhci", "mhcii", "mixmhc2pred"]:
        return "Epitope Prediction"
    # Epitope evaluation
    if feature in ["dssp", "protlearn"]:
        return "Structure Analysis"
    return "Other"

def plot_auroc_summary(results_df, output_dir):
    """
    Create a bar plot of AUROC values for all features, categorized and colored.
    AUROCs < 0.5 are adjusted as 1 - AUROC.
    Args:
        results_df (pd.DataFrame): DataFrame containing feature results with columns:
            - "feature": feature name (e.g., "bcell", "mhci")
            - "subfeature": subfeature name (e.g., "bepipred", "score")
            - "auroc": AUROC value (float)
        output_dir (str): Directory to save the AUROC summary plot.
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

        # Categorize features
        df["category"] = df.apply(lambda row: categorize_feature(row["feature"], row["subfeature"]), axis=1)

        # Sort
        df = df.sort_values("adjusted_auroc", ascending=False)

        # Define color palette for categories
        category_palette = {
            "Subcellular localisation": "#1b9e77",
            "Allergenicity": "#d95f02",
            "Immunogenicity": "#7570b3",
            "Conservation Analysis Across Strains": "#e7298a",
            "Epitope Prediction": "#66a61e",
            "Epitope evaluation": "#e6ab02",
            "Other": "#a6761d"
        }
        # Map colors
        colors = df["category"].map(category_palette).fillna("#a6761d")

        # Plot
        plt.figure(figsize=(10, max(4, 0.3 * len(df))))
        bars = plt.barh(df["label"], df["adjusted_auroc"], color=colors)
        plt.xlabel("AUROC (adjusted, min = 0.5)")
        plt.title("AUROC Summary (Sorted High to Low)")
        plt.xlim(0.5, 1.0)
        plt.gca().invert_yaxis()  # Highest on top

        # Add value labels
        for bar in bars:
            width = bar.get_width()
            plt.text(width + 0.01, bar.get_y() + bar.get_height()/2,
                     f"{width:.3f}", va="center")

        # Add legend for categories
        handles = []
        for cat, color in category_palette.items():
            handles.append(plt.Rectangle((0,0),1,1, color=color, label=cat))
        plt.legend(handles=handles, title="Category", loc="lower right", fontsize=10)

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
    Args:
        pos_features (list): List of dictionaries representing positive features.
        rand_features (list): List of dictionaries representing random features.
        output_dir (str): Directory to save the results.
    Returns:
        pd.DataFrame: DataFrame containing KS statistics, p-values, AUROC values, and counts.
    Each row contains:
        - "feature": feature name (e.g., "bcell", "mhci")
        - "subfeature": subfeature name (e.g., "bepipred", "score")
        - "ks_statistic": KS statistic value (float)
        - "p_value": p-value from KS test (float)
        - "auroc": AUROC value (float)
        - "positive_n": number of positive samples (int)
        - "random_n": number of random samples (int)
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
