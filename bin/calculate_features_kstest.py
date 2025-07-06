import os
import json
import csv
import sys
import pandas as pd
from scipy.stats import ks_2samp
from statistics import mean, median
from collections import defaultdict
from Bio import SeqIO
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------- Utility Functions -----------------------------

def safe_mean(lst):
    """Compute mean safely (returns 0.0 if list is empty)."""
    return mean(lst) if lst else 0.0

def safe_median(lst):
    """Compute median safely (returns 0.0 if list is empty)."""
    return median(lst) if lst else 0.0

def init_logging(verbose=False, pathogen="unknown"):
    """
    Initializes logging to stdout and optionally to a file.

    Args:
        verbose (bool): If True, logs to both console and 'log.txt'.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if verbose:
        fh = logging.FileHandler(f"data/{pathogen}/log.txt", mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

def sizeof_fmt(num, suffix="B"):
    """Return human-readable file size string."""
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return f"{num:.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}P{suffix}"

# ----------------------------- Feature Parsers -----------------------------

def parse_bcell_dir(directory):
    """
    Parses B-cell epitope prediction scores.

    Returns:
        dict: Keys are method names (e.g., 'bcell_bepipred'), values are dicts with 'mean' and 'median'.
    """
    scores = defaultdict(list)

    for file in os.listdir(directory):
        if not file.endswith(".csv"):
            continue
        method = os.path.basename(file).split("_")[-1].replace(".csv", "")
        path = os.path.join(directory, file)
        with open(path) as f:
            reader = csv.reader(f)
            headers = next(reader)
            idx = headers.index("Score")
            for row in reader:
                try:
                    score = float(row[idx])
                    method_key = f"bcell_{method.lower()}"
                    scores[method_key].append(score)
                except Exception:
                    continue

    return {key: {"mean": safe_mean(vals), "median": safe_median(vals)} for key, vals in scores.items()}

def parse_mhc_dir(directory):
    """
    Parses MHC class I/II binding predictions.

    Returns:
        dict: Contains scores for 'mhc_score', 'mhc_percentile', 'mhc_peptide_length'.
    """
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".json"):
            continue
        path = os.path.join(directory, file)
        try:
            with open(path) as f:
                data = json.load(f)
                for result in data.get("results", []):
                    if result.get("type") != "peptide_table":
                        continue
                    cols = result["table_columns"]
                    idx_score = cols.index("score")
                    idx_percentile = cols.index("percentile")
                    idx_peptide = cols.index("peptide")
                    for row in result["table_data"]:
                        try:
                            scores["score"].append(float(row[idx_score]))
                            scores["percentile"].append(float(row[idx_percentile]))
                            scores["peptide_length"].append(len(row[idx_peptide]))
                        except Exception:
                            continue
        except Exception:
            continue

    return {f"mhc_{k}": {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_signalp_dir(directory):
    """
    Parses SignalP results.

    Returns:
        dict: Includes 'signalp_prob_signalp', 'signalp_prob_other'.
    """
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    try:
                        scores["signalp_prob_signalp"].append(float(parts[2]))
                        scores["signalp_prob_other"].append(float(parts[3]))
                    except Exception:
                        continue

    return {k: {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_targetp_dir(directory):
    """
    Parses TargetP results.

    Returns:
        dict: Includes 'targetp_prob_noTP', 'targetp_prob_SP', 'targetp_prob_mTP'.
    """
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    try:
                        scores["targetp_prob_noTP"].append(float(parts[2]))
                        scores["targetp_prob_SP"].append(float(parts[3]))
                        scores["targetp_prob_mTP"].append(float(parts[4]))
                    except Exception:
                        continue

    return {k: {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_allergenicity_dir(directory):
    """
    Parses allergenicity results from AlgPred.

    Returns:
        dict: Includes hybrid, ML, MERCI, and BLAST scores under 'allergenicity_*'.
    """
    hybrid, ml, merci, blast = [], [], [], []

    for file in os.listdir(directory):
        if not file.endswith("_algpred.csv"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    hybrid.append(float(row.get("Hybrid Score", 0)))
                    ml.append(float(row.get("ML Score", 0)))
                    merci.append(float(row.get("MERCI Score", 0)))
                    blast.append(float(row.get("BLAST Score", 0)))
                except Exception:
                    continue

    return {
        "allergenicity_hybrid": {"mean": safe_mean(hybrid), "median": safe_median(hybrid)},
        "allergenicity_ml": {"mean": safe_mean(ml), "median": safe_median(ml)},
        "allergenicity_merci": {"mean": safe_mean(merci), "median": safe_median(merci)},
        "allergenicity_blast": {"mean": safe_mean(blast), "median": safe_median(blast)}
    }

def parse_cluster_dir(directory):
    """
    Parses sequence similarity clustering files (.m8).

    Returns:
        dict: Contains 'cluster_conservation_score' with mean/median percent identity.
    """
    clusters = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".m8"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                try:
                    percent_identity = float(parts[2])
                    clusters[parts[0]].append(percent_identity)
                except Exception:
                    continue

    conservation_scores = [safe_mean(vals) for vals in clusters.values() if vals]

    return {
        "cluster_conservation_score": {
            "mean": safe_mean(conservation_scores),
            "median": safe_median(conservation_scores)
        }
    }

def parse_popcov_dir(directory):
    """
    Parses population coverage results.

    Returns:
        dict: Includes 'popcov_percent_individuals' and 'popcov_cumulative_coverage'.
    """
    individuals, coverage = [], []

    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            continue
        path = os.path.join(directory, file)
        try:
            with open(path) as f:
                in_table = False
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith("population/area") and "cumulative_coverage" in line:
                        in_table = True
                        continue
                    if in_table:
                        parts = line.split('\t')
                        if len(parts) == 4:
                            try:
                                individuals.append(float(parts[2]))
                                coverage.append(float(parts[3]))
                            except Exception:
                                continue
        except Exception:
            continue

    return {
        "popcov_percent_individuals": {"mean": safe_mean(individuals), "median": safe_median(individuals)},
        "popcov_cumulative_coverage": {"mean": safe_mean(coverage), "median": safe_median(coverage)}
    }

# ----------------------------- Orchestration Functions -----------------------------

def extract_all_features(base_dir, eval_dir, threads=1):
    """
    Extracts features from all prediction directories.

    Args:
        base_dir (str): Base path to look for subdirectories.
        threads (int): Number of threads to use.

    Returns:
        dict: Key = feature group name, value = parsed result.
    """
    parsers = {
        "bcell": lambda: parse_bcell_dir(os.path.join(base_dir, "bcell")),
        "mhci": lambda: parse_mhc_dir(os.path.join(base_dir, "mhci")),
        "mhcii": lambda: parse_mhc_dir(os.path.join(base_dir, "mhcii")),
        "signalp": lambda: parse_signalp_dir(os.path.join(base_dir, "signalp")),
        "targetp": lambda: parse_targetp_dir(os.path.join(base_dir, "targetp")),
        "allergenicity": lambda: parse_allergenicity_dir(os.path.join(eval_dir, "allergenicity")),
        "cluster": lambda: parse_cluster_dir(os.path.join(eval_dir, "cluster")),
        "popcoverage": lambda: parse_popcov_dir(os.path.join(eval_dir, "popcoverage")),
    }

    results = {}
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fn): name for name, fn in parsers.items()}
        for future in as_completed(futures):
            name = futures[future]
            try:
                results[name] = future.result()
            except Exception as e:
                print(f"Error parsing {name}: {e}")

    return results

def compare_ks(pos_features, rand_features):
    """
    Performs Kolmogorov–Smirnov test between positive and random features.

    Returns:
        pd.DataFrame: Results with KS statistic and p-value per (feature, subfeature).
    """
    results = []

    for feature, pos_data in pos_features.items():
        rand_data = rand_features.get(feature, {})

        if isinstance(pos_data, dict):
            for subfeature, val in pos_data.items():
                pos_vals = [val["mean"]]
                rand_vals = [rand_data.get(subfeature, {}).get("mean", 0)]

                try:
                    stat, pval = ks_2samp(pos_vals, rand_vals)
                except Exception:
                    stat, pval = None, None

                results.append({
                    "feature": feature,
                    "subfeature": subfeature,
                    "ks_statistic": stat,
                    "p_value": pval,
                    "positive_n": len(pos_vals),
                    "random_n": len(rand_vals)
                })

    return pd.DataFrame(results)

def write_features_by_feature(features, label, output_dir):
    """
    Writes raw feature data to CSV by feature.

    Args:
        features (dict): Dictionary of features and subfeatures.
        label (str): 'positive' or 'random'.
        output_dir (str): Directory to write CSV files.
    """
    os.makedirs(output_dir, exist_ok=True)

    for feature, subdata in features.items():
        filepath = os.path.join(output_dir, f"{feature}_{label}_raw_data.csv")

        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["label", "feature", "subfeature", "value"])
            writer.writeheader()

            if isinstance(subdata, dict):
                for subfeature, stats in subdata.items():
                    for k in ("mean", "median"):
                        writer.writerow({
                            "label": label,
                            "feature": feature,
                            "subfeature": f"{subfeature}_{k}",
                            "value": stats[k]
                        })

# ----------------------------- Entry Point -----------------------------

def main(pathogen_dir, threads, verbose=False):
    init_logging(verbose, pathogen_dir)
    logger = logging.getLogger()

    pos_dir = os.path.join("data", pathogen_dir, "epitope_outputs")
    rand_dir = os.path.join("data", pathogen_dir, "random_analysis")
    pos_eval_dir = os.path.join("data", pathogen_dir, "evaluation_outputs")
    rand_eval_dir = os.path.join("data", pathogen_dir, "random_evaluation")

    logger.info(f"Extracting features for {pathogen_dir} using {threads} thread(s)")
    pos_features = extract_all_features(pos_dir, pos_eval_dir, threads)
    rand_features = extract_all_features(rand_dir, rand_eval_dir, threads)

    logger.info("Estimating memory usage...")
    logger.info(f"Positive features: {sizeof_fmt(sys.getsizeof(pos_features))}")
    logger.info(f"Random features: {sizeof_fmt(sys.getsizeof(rand_features))}")

    raw_out_dir_pos = os.path.join("data", pathogen_dir, "raw_positive_features")
    raw_out_dir_rand = os.path.join("data", pathogen_dir, "raw_random_features")

    write_features_by_feature(pos_features, "positive", raw_out_dir_pos)
    write_features_by_feature(rand_features, "random", raw_out_dir_rand)

    logger.info("Running KS test...")
    result_df = compare_ks(pos_features, rand_features)
    logger.info("\n" + result_df.to_string(index=False))

    result_df.to_csv(os.path.join("data", pathogen_dir, "ks_test_results.csv"), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KS-test comparison of epitope vs. random features.")
    parser.add_argument("pathogen_dir", help="Pathogen directory name under data/")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for parsing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to file")
    args = parser.parse_args()

    main(args.pathogen_dir, args.threads, args.verbose)
