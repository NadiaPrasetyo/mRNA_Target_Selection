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
from scipy.stats import ks_2samp, ttest_ind
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
import traceback
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
    Parse B-cell epitope prediction CSV files in a directory.

    Expected format:
    - A metadata line starting with "input:".
    - A peptide table between:
        "No,Start,End,Peptipe,Length"
        ... and ...
        "Position,Residue,Score,Assignment"

    Extracted features:
        - Number of peptides
        - Average peptide length

    Args:
        directory (str): Path to directory containing B-cell CSV files.

    Returns:
        List[dict]: Each dictionary contains:
            - "feature": always "bcell"
            - "subfeature": method-specific metric (e.g., "bepipred_num_peptides")
            - "value": numerical value
    """
    logging.info(f"Parsing B-cell features in {directory}")
    results = []

    try:
        subdirs = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
        logging.debug(f"Found {len(subdirs)} subdirectories in B-cell dir")
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for subdir in subdirs:
        subdir_path = os.path.join(directory, subdir)
        files = [f for f in os.listdir(subdir_path)
                if "Bcell_epitope_preds" in f and f.endswith(".fasta")]

        for file in files:
            path = os.path.join(subdir_path, file)  
            logging.debug(f"Parsing B-cell file: {path}")
            # parse_file(path)
            try:
                with open(path) as f:
                    peptide_lengths = defaultdict(list)
                    num_peptide = defaultdict(int)
                    lines = f.readlines()
                    i = 0
                    while i < len(lines):
                        # Look for header line
                        if lines[i].startswith(">"):
                            if i + 1 >= len(lines):
                                logging.warning(f"File {file} does not have enough lines after header at line {i}")
                                break
                            header = lines[i].strip()
                            accession = header.replace(">", "")  # Extract accession from header
                            sequence = lines[i + 1].strip()
                            peptide_length = 0
                            # get peptides from sequence, peptides have capitalized letters but we only want peptides with a minimum of 5 aa
                            for char in sequence:
                                if char.isupper():
                                    peptide_length += 1
                                else:
                                    if peptide_length >= 5:
                                        peptide_lengths[accession].append(peptide_length)
                                        num_peptide[accession]+=1

                                    peptide_length = 0

                            i += 2  # Move to next header or end
                        else:
                            i += 1  # Skip lines until next header

                    for accession in peptide_lengths:
                        avg_length = float(statistics.mean(peptide_lengths[accession])) if peptide_lengths[accession] else 0
                        total_peptides = num_peptide[accession] if accession in num_peptide else 0
                        results.append({"accession": accession, "feature": "bcell", "subfeature": f"num_peptides", "value": total_peptides})
                        results.append({"accession": accession, "feature": "bcell", "subfeature": f"avg_peptide_length", "value": avg_length})
            
            except Exception as e:
                logging.error(f"Failed parsing B-cell file {file}: {e}")

    logging.info(f"Completed B-cell parsing with {len(results)} results")
    return results

def parse_mhc_dir(directory):
    """
    Parse MHC epitope prediction files in the specified directory.
    Expected files are JSONs with MHC I/II prediction results.
    Filters peptides based on %Rank thresholds for SBs and WBs.
    Returns a list of dictionaries with feature values.
    Args:
        directory (str): Path to the directory containing MHC prediction files.
    Returns:
        List of dictionaries, where each dictionary represents an MHC feature.
    Each dictionary contains:
        - "feature": "mhci" or "mhcii" based on the directory name
        - "subfeature": "score", "percentile", "peptide_length", or "num_peptides"
        - "value": numerical score, length, or count value
    """
    prefix = "mhcii" if "mhcii" in directory else "mhci"
    logging.info(f"Parsing MHC dir {directory} with prefix {prefix}")
    results = []
    try:
        files = [f for f in os.listdir(directory) if "matched_antigens" in f.lower()]
        logging.debug(f"Found {len(files)} out files in MHC dir")
    except Exception as e:
        logging.error(f"Failed listing directory {directory}: {e}")
        return results

    for file in files:
        path = os.path.join(directory, file)
        logging.debug(f"Parsing MHC file: {file}")
        try:
            strain = Path(file).stem.split("_")[0]
            with open(path) as f:
                data = f.readlines()
                num_peptides = defaultdict(int)
                scores = defaultdict(list)
                percentiles = defaultdict(list)
                num_sb = defaultdict(int)
                num_wb = defaultdict(int)

                for i, line in enumerate(data):
                    if line.startswith("#") or line.startswith("-") or not line.strip() or "Pos" in line:
                        continue
                    # skip all lines without <=
                    if "<=" not in line:
                        continue

                    parts = line.split()
                    if len(parts) < 10:
                        logging.debug(f"Skipping malformed line {i} in {file}: {line.strip()}")
                        continue
                    try:
                        id = parts[10 if prefix == "mhci" else 7]
                        accession = f'{id.split("_")[0]}_{strain}'
                        score = float(parts[11 if prefix == "mhci" else 8])
                        percentile = float(parts[12 if prefix == "mhci" else 9])
                        binding_strength = parts[14 if prefix == "mhci" else 12] if len(parts) > (13 if prefix == "mhci" else 11) else "NA"

                        logging.debug(f"Parsed line {i} in {file}: accession={accession}, score={score}, percentile={percentile}, binding_strength={binding_strength}")
                        num_peptides[accession] += 1
                        # Filter for Strong Binders only
                        if "SB" in binding_strength:
                            num_sb[accession] += 1
                        elif "WB" in binding_strength:
                            num_wb[accession] += 1

                        scores[accession].append(score)
                        percentiles[accession].append(percentile)
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

                for accession in num_peptides:
                    avg_score = float(statistics.mean(scores[accession])) if scores[accession] else 0
                    avg_percentile = float(statistics.mean(percentiles[accession])) if percentiles[accession] else 0
                    results.append({"accession": accession, "feature": prefix, "subfeature": "score", "value": avg_score})
                    results.append({"accession": accession, "feature": prefix, "subfeature": "percentile", "value": avg_percentile})
                    results.append({"accession": accession, "feature": prefix, "subfeature": "num_peptides", "value": num_peptides[accession]})
                    results.append({"accession": accession, "feature": prefix, "subfeature": "num_strong_binders", "value": num_sb[accession]})
                    results.append({"accession": accession, "feature": prefix, "subfeature": "num_weak_binders", "value": num_wb[accession]})

                
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
                    accession = parts[0].split('|')[1] + "_" + parts[0].split('|')[3]	
                    if len(parts) >= 4:
                        try:
                            results.append({"accession": accession, "feature": "signalp", "subfeature": "prob_signalp", "value": float(parts[2])})
                            results.append({"accession": accession, "feature": "signalp", "subfeature": "prob_other", "value": float(parts[3])})
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
                        accession = parts[0].split('|')[1] + "_" + parts[0].split('|')[3]
                        try:
                            results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_noTP", "value": float(parts[2])})
                            results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_SP", "value": float(parts[3])})
                            results.append({"accession": accession, "feature": "targetp", "subfeature": "prob_mTP", "value": float(parts[4])})
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
                        accession = row["Subject"].split('|')[1] + "_" + row["Subject"].split('|')[3]
                        merci_score = float(row["MERCI Score"])
                        blast_score = float(row["BLAST Score"])
                        hybrid_score = float(row["Hybrid Score"])

                        results.extend([
                            {
                                "accession": accession,
                                "feature": "allergenicity",
                                "subfeature": "merci_score",
                                "value": merci_score
                            },
                            {
                                "accession": accession,
                                "feature": "allergenicity",
                                "subfeature": "blast_score",
                                "value": blast_score
                            },
                            {
                                "accession": accession,
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
            accession = cluster_id
            percent_identity_num_strain = sum(scores[0::3]) / len(unique_strains)  #sum of percent identities divided by number of strains
            avg_bit_score = sum(scores[1::3]) / len(scores[1::3]) # sum of bit scores divided by number of scores
            avg_log_e_value = sum(scores[2::3]) / len(scores[2::3]) # sum of e-values divided by number of e-values
            results.append({
                "accession": accession,
                "feature": "cluster_conservation",
                "subfeature": "percent_identity/num_strain",
                "value": percent_identity_num_strain
            })
            results.append({
                "accession": accession,
                "feature": "cluster_conservation",
                "subfeature": "bit_score_normalized",
                "value": avg_bit_score
            })
            results.append({
                "accession": accession,
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
                    try:
                        accession = f"{row['ACC'].split('|')[1]}_{row['ACC'].split('|')[3]}"
                    except IndexError:
                        logging.warning(f"Unexpected ACC format: {row['ACC']}")
                        accession = "unknown"

                    for loc, col in prob_cols:
                        try:
                            prob = float(row[col])
                            results.append({
                                "accession": accession,
                                "feature": "deeplocpro",
                                "subfeature": f"prob_{loc}",
                                "value": prob
                            })
                        except ValueError:
                            logging.debug(f"Invalid probability value in column {col} for row {i}")
                        except KeyError:
                            logging.debug(f"Missing column {col} in row {i}")
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
        linear_score = []
        discontinuous_score = []
        # Determine accession from filename
        if "_" in os.path.basename(file):
            accession = os.path.basename(file).split("_")[0].replace(".txt", "") if os.path.basename(file).split("_")[1] != "AF.txt" else os.path.basename(file).split("_")[1].replace(".txt", "")
        else:
            accession = os.path.basename(file).replace(".txt", "")
        try:
            with open(path) as f:
                lines = f.readlines()
                # Parse linear epitopes
                for i, line in enumerate(lines):
                    if i == 0 or "No." in line:
                        continue  # Skip header
                    parts = line.strip().rsplit(',', 2)  # Split into three parts from the right
                    if len(parts) < 3:
                        continue  # Skip malformed lines
                    try:
                        # Take the second-to-last field as the score
                        linear_score.append(float(parts[-2]))
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

                # Parse discontinuous epitopes
                for i, line in enumerate(lines):
                    if i == 0 or "No." in line:
                        continue  # Skip header
                    parts = line.strip().rsplit(',', 2)  # Split into three parts from the right
                    if len(parts) < 3:
                        continue  # Skip malformed lines
                    try:
                        discontinuous_score.append(float(parts[-2]))  # Take the second-to-last field as the score
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
            # Store results
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "mean_linear_score",
                "value": statistics.mean(linear_score) if linear_score else 0
            })
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "mean_discontinuous_score",
                "value": statistics.mean(discontinuous_score) if discontinuous_score else 0
            })
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "med_linear_score",
                "value": statistics.median(linear_score) if linear_score else 0
            })
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "med_discontinuous_score",
                "value": statistics.median(discontinuous_score) if discontinuous_score else 0
            })
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "num_linear_epitopes",
                "value": len(linear_score)
            })
            results.append({
                "accession": accession,
                "feature": "ellipro",
                "subfeature": "num_discontinuous_epitopes",
                "value": len(discontinuous_score)
            })
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
        ml_score = defaultdict(list)
        blast_score = defaultdict(list)
        total_score = defaultdict(list)
        try:
            with open(path) as f:
                reader = csv.DictReader(f)
                for i, row in enumerate(reader):
                    try:
                        accession = row["Seq_ID"].split('|')[1] + "_" + row["Seq_ID"].split('|')[3]
                        ml_score[accession].append(float(row["ML_Score"]))
                        blast_score[accession].append(float(row["BLAST_Score"]))
                        total_score[accession].append(float(row["Total_Score"]))
                        
                    except ValueError as e:
                        logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")

            for accession in ml_score:
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "mean_ml_score",
                    "value": statistics.mean(ml_score[accession]) if ml_score[accession] else 0
                })
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "mean_blast_score",
                    "value": statistics.mean(blast_score[accession]) if blast_score[accession] else 0
                })
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "mean_total_score",
                    "value": statistics.mean(total_score[accession]) if total_score[accession] else 0
                })
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "med_ml_score",
                    "value": statistics.median(ml_score[accession]) if ml_score[accession] else 0
                })
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "med_blast_score",
                    "value": statistics.median(blast_score[accession]) if blast_score[accession] else 0
                })
                results.append({
                    "accession": accession,
                    "feature": "ifnepitope2",
                    "subfeature": "med_total_score",
                    "value": statistics.median(total_score[accession]) if total_score[accession] else 0
                })
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
                        header = lines[i].strip()
                        accession = header.split('|')[1] + "_" + header.split('|')[3]  # Extract accession from header
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
                            "accession": accession,
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
        accession = os.path.basename(file).split("_")[0] if "_" in os.path.basename(file) else os.path.basename(file).replace(".out", "")
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
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")

            results.append({
                "accession": accession,
                "feature": "rate4site",
                "subfeature": "mean_per_site_score",
                "value": statistics.mean(scores) if scores else 0
            })
            results.append({
                "accession": accession,
                "feature": "rate4site",
                "subfeature": "median_per_site_score",
                "value": statistics.median(scores) if scores else 0
            })
            results.append({
                "accession": accession,
                "feature": "rate4site",
                "subfeature": "min_per_site_score",
                "value": min(scores) if scores else 0
            })
            results.append({
                "accession": accession,
                "feature": "rate4site",
                "subfeature": "max_per_site_score",
                "value": max(scores) if scores else 0
            })

            # # Compute rolling window statistics (15-peptide window)
            # window_size = 15
            # for j in range(len(scores)):
            #     window = scores[max(0, j - window_size + 1):j + 1]
            #     if len(window) > 0:
            #         results.append({
            #             "feature": "rate4site",
            #             "subfeature": "rolling_mean",
            #             "value": sum(window) / len(window)
            #         })
            #         results.append({
            #             "feature": "rate4site",
            #             "subfeature": "rolling_median",
            #             "value": statistics.median(window)
            #         })
            #         results.append({
            #             "feature": "rate4site",
            #             "subfeature": "rolling_max",
            #             "value": max(window)
            #         })
            #         results.append({
            #             "feature": "rate4site",
            #             "subfeature": "rolling_min",
            #             "value": min(window)
            #         })
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
            outside_per_site_scores = []
            accession = os.path.basename(rate4site_file).split("_")[0]
            mafft_file = os.path.join(mafft_dir, f"{accession}_combined_aligned.fasta")

            if not os.path.exists(mafft_file):
                logging.warning(f"Missing MAFFT file: {mafft_file}")
                continue

            logging.info(f"Processing rate4site accession: {accession}")
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

                    outside_per_site_scores.append(score)

                    strain_topology_positions.append((pos, aa, score, topo))

                results.append({
                    "accession": accession,
                    "feature": "rate4site_deeptmhmm",
                    "subfeature": "mean_outside_scores",
                    "value": statistics.mean(outside_per_site_scores) if outside_per_site_scores else None
                })
                results.append({
                    "accession": accession,
                    "feature": "rate4site_deeptmhmm",
                    "subfature": "median_outside_scores",
                    "value": statistics.median(outside_per_site_scores) if outside_per_site_scores else None
                })
                results.append({
                    "accession": accession,
                    "feature": "rate4site_deeptmhmm",
                    "subfeature": "min_outside_scores",
                    "value": min(outside_per_site_scores) if outside_per_site_scores else None
                })
                results.append({
                    "accession": accession,
                    "feature": "rate4site_deeptmhmm",
                    "subfeature": "max_outside_scores",
                    "value": max(outside_per_site_scores) if outside_per_site_scores else None
                })

                # # Rolling window stats for this strain-antigen pair
                # outside_scores = [s for _, _, s, t in strain_topology_positions if t == "O"]
                # window_size = 15
                # for i in range(len(outside_scores) - window_size + 1):
                #     window = outside_scores[i:i + window_size]
                #     results.extend([
                #         {
                #             "strain": strain,
                #             "accession": accession,
                #             "feature": "rate4site_deeptmhmm",
                #             "subfeature": "rolling_avg",
                #             "value": sum(window) / window_size
                #         },
                #         {
                #             "strain": strain,
                #             "accession": accession,
                #             "feature": "rate4site_deeptmhmm",
                #             "subfeature": "rolling_median",
                #             "value": statistics.median(window)
                #         },
                #         {
                #             "strain": strain,
                #             "accession": accession,
                #             "feature": "rate4site_deeptmhmm",
                #             "subfeature": "rolling_min",
                #             "value": min(window)
                #         },
                #         {
                #             "strain": strain,
                #             "accession": accession,
                #             "feature": "rate4site_deeptmhmm",
                #             "subfeature": "rolling_max",
                #             "value": max(window)
                #         }
                #     ])

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
        if "_" in dssp_file.stem:
            accession = dssp_file.stem.split("_")[0] if dssp_file.stem.split("_")[1] != "AF" else dssp_file.stem.split("_")[1]
        else:
            accession = dssp_file.stem
        contents = dssp_file.read_text()
        percent_helix = contents.count('H') / len(contents) if len(contents) > 0 else 0
        percent_sheet = contents.count('E') / len(contents) if len(contents) > 0 else 0
        percent_loop = contents.count('-') / len(contents) if len(contents) > 0 else 0

        results.append({
            "accession": accession,
            "feature": "dssp",
            "subfeature": "percent_helix",
            "value": percent_helix
        })
        results.append({
            "accession": accession,
            "feature": "dssp",
            "subfeature": "percent_sheet",
            "value": percent_sheet
        })
        results.append({
            "accession": accession,
            "feature": "dssp",
            "subfeature": "percent_loop",
            "value": percent_loop
        })

    return results


def parse_dnds_dir(directory):
    """
    Parse the DNDS directory for structural features.
    Handles FEL, SLAC, and FUBAR result files.
    """

    def safe_div(n, d):
        """Safely divide n / d, return 0 if invalid or denominator <= 0."""
        try:
            if d is None or n is None or d <= 0:
                return 0
            return n / d
        except Exception:
            return 0

    results = []
    for dnds_file in Path(directory).glob("*.json"):
        try:
            # Skip empty files quickly
            if dnds_file.stat().st_size == 0:
                logging.warning(f"Skipping empty JSON file: {dnds_file}")
                continue
            accession = dnds_file.stem.split("_")[0] if "_" in dnds_file.stem else dnds_file.stem
            type = dnds_file.stem.split("_")[-2]  # Extract method type (FEL, SLAC, FUBAR)
            with open(dnds_file, "r") as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    logging.warning(f"Skipping invalid JSON file: {dnds_file}")
                    continue

            content = data.get("MLE", {}).get("content", {}).get("0")
            if not content:
                logging.warning(f"Missing or invalid content in {dnds_file}")
                continue

            match type:
                case "FEL":
                    sum_alpha, sum_beta, sum_lrt, count = 0, 0, 0, 0
                    for each in content:
                        alpha, beta, lrt = each[0], each[1], each[3]
                        if alpha is not None and alpha != 0 and beta is not None and beta != 0:
                            sum_alpha += alpha
                            sum_beta += beta
                            sum_lrt += lrt
                            count += 1

                    results.append({
                        "accession": accession,
                        "feature": "FEL",
                        "subfeature": "mean_alpha",
                        "value": safe_div(sum_alpha, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "FEL",
                        "subfeature": "mean_beta",
                        "value": safe_div(sum_beta, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "FEL",
                        "subfeature": "mean_lrt",
                        "value": safe_div(sum_lrt, count)
                    })

                case "SLAC":
                    content = content.get("by-site", {}).get("AVERAGED", [])
                    if not content:
                        logging.warning(f"Missing or invalid SLAC content in {dnds_file}")
                        continue

                    sum_n, sum_s, sum_dn, sum_ds, count = 0, 0, 0, 0, 0
                    for each in content:
                        s, n, ds, dn = each[2], each[3], each[5], each[6]
                        if n is not None and s is not None and ds is not None and dn is not None:
                            if n != s and s != 0:
                                sum_n += n
                                sum_s += s
                                sum_dn += dn
                                sum_ds += ds
                                count += 1
                        
                    results.append({
                        "accession": accession,
                        "feature": "SLAC",
                        "subfeature": "mean_dN/dS",
                        "value": safe_div(sum_dn, sum_ds)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "SLAC",
                        "subfeature": "mean_n",
                        "value": safe_div(sum_n, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "SLAC",
                        "subfeature": "mean_s",
                        "value": safe_div(sum_s, count)
                    })

                case "FUBAR":
                    sum_alpha, sum_beta, sum_prob_pos, sum_prob_neg, count = 0, 0, 0, 0, 0
                    for each in content:
                        alpha, beta, difference, prob_pos, prob_neg = each[0], each[1], each[2], each[3], each[4]
                        # Only counting the positive difference (non synonymous > synonymous changes) and if they are not none
                        if difference>0 and alpha is not None and beta is not None and prob_pos is not None and prob_neg is not None:
                            sum_alpha += alpha
                            sum_beta += beta
                            sum_prob_pos += prob_pos
                            sum_prob_neg += prob_neg
                            count += 1

                    results.append({
                        "accession": accession,
                        "feature": "FUBAR",
                        "subfeature": "mean_alpha",
                        "value": safe_div(sum_alpha, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "FUBAR",
                        "subfeature": "mean_beta",
                        "value": safe_div(sum_beta, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "FUBAR",
                        "subfeature": "mean_prob_pos",
                        "value": safe_div(sum_prob_pos, count)
                    })

                    results.append({
                        "accession": accession,
                        "feature": "FUBAR",
                        "subfeature": "mean_prob_neg",
                        "value": safe_div(sum_prob_neg, count)
                    })

                case _:
                    logging.warning(f"Unknown dN/dS type: {type}")
                    continue

        except Exception as e:
            logging.error(f"Failed parsing {dnds_file}: {e}", exc_info=True)
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
        if "_" in protlearn_file.stem:
            accession = protlearn_file.stem.split("_")[0] if protlearn_file.stem.split("_")[1] != "AF" else protlearn_file.stem.split("_")[1]
        with open(protlearn_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                feature = row["feature"]
                value = row["value"].strip("[]") if row["value"] is not None else None  # Remove brackets if present
                try:
                    if value is not None:  # Ensure value is not None before conversion
                        value = int(float(value))  # Convert to float first, then to int
                        results.append({
                        "accession": accession,
                        "feature": "ProtLearn",
                        "subfeature": feature,
                        "value": value
                    })
                except ValueError as e:
                    logging.debug(f"Skipping row due to conversion error: {e}")
    return results

# ----------------------------- Orchestration Functions -----------------------------

def extract_all_features(base_dir, threads=1):
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
    logging.info(f"Extracting features from base_dir: {base_dir} using {threads} threads")
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
    If the AUROC is exactly 0.5, it will skip plotting.
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

        auc = roc_auc_score(y_true, y_scores)

        # Skip plotting if AUROC is exactly 0.5
        if auc == 0.5 or auc == 0.500:
            logging.info(f"Skipping ROC plot for {feature}/{subfeature} as AUROC is 0.5")
            return

        fpr, tpr, _ = roc_curve(y_true, y_scores)

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
    if feature in ["cluster_conservation", "rate4site", "rate4site_deeptmhmm", "FEL", "SLAC", "FUBAR"] or subfeature in [
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
    if feature in ["dssp", "ProtLearn"]:
        return "Structure Analysis"
    return "Other"

def plot_auroc_summary(results_df, output_dir, prefix):
    """
    Create a bar plot of AUROC values for all features, categorized and colored.
    AUROCs < 0.5 are adjusted as 1 - AUROC. Excludes data where AUROC is exactly 0.5.
    Bars are solid for positive t-test results and striped for negative t-test results.
    Args:
        results_df (pd.DataFrame): DataFrame containing feature results with columns:
            - "feature": feature name (e.g., "bcell", "mhci")
            - "subfeature": subfeature name (e.g., "bepipred", "score")
            - "auroc": AUROC value (float)
            - "t_statistic": t-test statistic (float)
        output_dir (str): Directory to save the AUROC summary plot.
    """
    output_path = os.path.join(output_dir, f"auroc_summary_{prefix}.png")
    os.makedirs(output_dir, exist_ok=True)
    try:
        # Drop missing AUROCs and exclude AUROCs exactly 0.5
        df = results_df.dropna(subset=["auroc", "t_statistic"]).copy()
        df = df[df["auroc"] != 0.5]
        # filter for p value <0.05
        df = df[df["p_value"] < 0.05]

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
            "Structure Analysis": "#d010e1",
            "Other": "#a6761d"
        }
        # Map colors
        colors = df["category"].map(category_palette).fillna("#a6761d")

        # Determine bar patterns: solid if t_statistic >=0, striped if <0
        hatches = ['' if t >= 0 else '////' for t in df["t_statistic"]]

        # Plot
        plt.figure(figsize=(10, max(4, 0.3 * len(df))))
        bars = plt.barh(df["label"], df["adjusted_auroc"], color=colors, hatch=None)
        # Apply hatches
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)

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


def compare_ks(pos_features, rand_features, output_dir, prefix):
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
        - f"{prefix}_n": number of random samples (int)
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
            f"{prefix}_n": len(rand_vals)
        })

    logging.info("KS and AUROC comparison complete")
    return pd.DataFrame(results)

    
def write_features_by_feature(features, label, output_dir):
    """
    Write features to disk, grouped by "feature" field.
    Each feature will be written to a separate CSV file named <feature>_<label>_raw_data.csv
    Args:
        features (list): List of dictionaries representing features.
        label (str): Label for the features (e.g., "positive", "random or human").
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
            writer = csv.DictWriter(f, fieldnames=["accession", "label", "feature", "subfeature", "value"])
            writer.writeheader()

            for row in rows:
                writer.writerow({
                    "accession": row.get("accession", ""),
                    "label": label,
                    "feature": row.get("feature", ""),
                    "subfeature": row.get("subfeature", ""),
                    "value": row.get("value", "")                    
                })
        logging.debug(f"Wrote feature file: {filepath}")


def add_ttest_results(result_df, pos_features, rand_features, logger):
    """
    Runs t-tests on features to determine directionality and updates result_df in place.
    
    Args:
        result_df (pd.DataFrame): DataFrame containing KS test results
        pos_features (list[dict]): Extracted positive feature data
        rand_features (list[dict]): Extracted random feature data
        logger (logging.Logger): Logger instance for debug/info messages
    
    Returns:
        pd.DataFrame: Updated DataFrame with t-test results
    """

    def group_values(data):
        grouped = defaultdict(list)
        for row in data:
            if "feature" in row and "subfeature" in row and "value" in row:
                key = (row["feature"], row["subfeature"])
                grouped[key].append(row["value"])
        return grouped

    pos_grouped = group_values(pos_features)
    rand_grouped = group_values(rand_features)

    # add empty columns if missing
    if "t_statistic" not in result_df.columns:
        result_df["t_statistic"] = None
    if "t_p_value" not in result_df.columns:
        result_df["t_p_value"] = None

    for index, row in result_df.iterrows():
        feature = row["feature"]
        subfeature = row["subfeature"]
        pos_vals = pos_grouped.get((feature, subfeature), [])
        rand_vals = rand_grouped.get((feature, subfeature), [])

        if pos_vals and rand_vals:
            try:
                t_stat, t_pval = ttest_ind(pos_vals, rand_vals, equal_var=False)
                result_df.at[index, "t_statistic"] = t_stat
                result_df.at[index, "t_p_value"] = t_pval
            except Exception as e:
                logger.debug(f"T-test failed for {feature}/{subfeature}: {e}")

    return result_df


# ----------------------------- Entry Point -----------------------------

def main(pathogen_dir, threads, verbose=False, write_raw=False, human_negative=False, raw_out_dir=None):
    """
    Main entry point for the script.
    Initializes logging, extracts features, performs KS tests, t-tests, and writes results.
    Args:
        pathogen_dir (str): Pathogen directory name under data/
        threads (int): Optional -number of threads to use for feature extraction.
        verbose (bool): Optional -if True, enables verbose logging to file.
        write_raw (bool): Optional -if True, writes raw feature data to disk (can be large).
    """
    init_logging(verbose, pathogen_dir)
    logger = logging.getLogger()

    prefix = "human" if human_negative else "random"

    logger.info(f"Starting processing for pathogen: {pathogen_dir}")

    pos_dir = os.path.join("data", pathogen_dir, "epitope_outputs")
    rand_dir = os.path.join("data", pathogen_dir, f"{prefix}_analysis")

    logger.info(f"Extracting positive features from {pos_dir}")
    pos_features = extract_all_features(pos_dir, threads)

    logger.info(f"Extracting random features from {rand_dir}")
    rand_features = extract_all_features(rand_dir, threads)

    logger.info("Estimating memory usage of extracted features")
    logger.info(f"Positive features: {sizeof_fmt(sys.getsizeof(pos_features))}")
    logger.info(f"{prefix} features: {sizeof_fmt(sys.getsizeof(rand_features))}")
    output_dir = os.path.join("results", pathogen_dir)



    if raw_out_dir is None:
        raw_out_dir = os.path.join("results", pathogen_dir, "raw_data")

    if write_raw:
        logger.info(f"Writing positive features to {raw_out_dir}")
        write_features_by_feature(pos_features, "positive", raw_out_dir)

        logger.info(f"Writing {prefix} features to {raw_out_dir}")
        write_features_by_feature(rand_features, prefix, raw_out_dir)

    logger.info("Running KS test on features")
    result_df = compare_ks(pos_features, rand_features, output_dir, prefix)

    logger.info("Running t-test on features to determine directionality")
    result_df = add_ttest_results(result_df, pos_features, rand_features, logger)

    plot_auroc_summary(result_df, output_dir, prefix)

    # Sort the DataFrame alphabetically by the first column
    result_df = result_df.sort_values(by=result_df.columns[0])
    logger.info("\n" + result_df.to_string(index=False))

    ks_out_path = os.path.join("results", pathogen_dir, f"ks_test_results_{prefix}.csv")
    logger.info(f"Writing KS test and t-test results to {ks_out_path}")
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
    parser.add_argument("--raw-output-dir", help="Write raw features to this output folder. Default: results/pathogen_dir/raw_data")
    parser.add_argument("--human", action="store_true", help="Use human negative set instead of random")
    args = parser.parse_args()

    main(args.pathogen_dir, args.threads, args.verbose, args.write_raw, args.human, args.raw_output_dir)
