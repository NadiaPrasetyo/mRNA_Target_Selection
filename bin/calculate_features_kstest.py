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

    # Clear existing handlers to avoid duplicate logs if re-run in notebook or similar
    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if verbose:
        os.makedirs(f"data/{pathogen}", exist_ok=True)
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
    logging.debug(f"Starting parse_bcell_dir on directory: {directory}")
    scores = defaultdict(list)

    for file in os.listdir(directory):
        if not file.endswith(".csv"):
            logging.debug(f"Skipping non-csv file in bcell dir: {file}")
            continue
        method = os.path.basename(file).split("_")[-1].replace(".csv", "")
        path = os.path.join(directory, file)
        logging.debug(f"Parsing B-cell file: {file} with method: {method}")
        with open(path) as f:
            reader = csv.reader(f)
            headers = next(reader)
            idx = headers.index("Score")
            for i, row in enumerate(reader, start=1):
                try:
                    score = float(row[idx])
                    method_key = f"bcell_{method.lower()}"
                    scores[method_key].append(score)
                except Exception as e:
                    logging.debug(f"Skipping bad score at line {i} in {file}: {e}")
                    continue
    logging.debug("Finished parse_bcell_dir")
    return {key: {"mean": safe_mean(vals), "median": safe_median(vals)} for key, vals in scores.items()}

def parse_mhc_dir(directory):
    logging.debug(f"Starting parse_mhc_dir on directory: {directory}")
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".json"):
            logging.debug(f"Skipping non-json file in mhc dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing MHC file: {file}")
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
                    for i, row in enumerate(result["table_data"], start=1):
                        try:
                            scores["score"].append(float(row[idx_score]))
                            scores["percentile"].append(float(row[idx_percentile]))
                            scores["peptide_length"].append(len(row[idx_peptide]))
                        except Exception as e:
                            logging.debug(f"Skipping bad row {i} in {file}: {e}")
                            continue
        except Exception as e:
            logging.warning(f"Failed to parse MHC file {file}: {e}")
            continue
    logging.debug("Finished parse_mhc_dir")
    return {f"mhc_{k}": {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_signalp_dir(directory):
    logging.debug(f"Starting parse_signalp_dir on directory: {directory}")
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            logging.debug(f"Skipping non-txt file in signalp dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing SignalP file: {file}")
        with open(path) as f:
            for i, line in enumerate(f, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    try:
                        scores["signalp_prob_signalp"].append(float(parts[2]))
                        scores["signalp_prob_other"].append(float(parts[3]))
                    except Exception as e:
                        logging.debug(f"Skipping bad line {i} in {file}: {e}")
                        continue
    logging.debug("Finished parse_signalp_dir")
    return {k: {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_targetp_dir(directory):
    logging.debug(f"Starting parse_targetp_dir on directory: {directory}")
    scores = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            logging.debug(f"Skipping non-txt file in targetp dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing TargetP file: {file}")
        with open(path) as f:
            for i, line in enumerate(f, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    try:
                        scores["targetp_prob_noTP"].append(float(parts[2]))
                        scores["targetp_prob_SP"].append(float(parts[3]))
                        scores["targetp_prob_mTP"].append(float(parts[4]))
                    except Exception as e:
                        logging.debug(f"Skipping bad line {i} in {file}: {e}")
                        continue
    logging.debug("Finished parse_targetp_dir")
    return {k: {"mean": safe_mean(v), "median": safe_median(v)} for k, v in scores.items()}

def parse_allergenicity_dir(directory):
    logging.debug(f"Starting parse_allergenicity_dir on directory: {directory}")
    hybrid, ml, merci, blast = [], [], [], []

    for file in os.listdir(directory):
        if not file.endswith("_algpred.csv"):
            logging.debug(f"Skipping non-algpred csv file in allergenicity dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing allergenicity file: {file}")
        with open(path) as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader, start=1):
                try:
                    hybrid.append(float(row.get("Hybrid Score", 0)))
                    ml.append(float(row.get("ML Score", 0)))
                    merci.append(float(row.get("MERCI Score", 0)))
                    blast.append(float(row.get("BLAST Score", 0)))
                except Exception as e:
                    logging.debug(f"Skipping bad row {i} in {file}: {e}")
                    continue

    logging.debug("Finished parse_allergenicity_dir")
    return {
        "allergenicity_hybrid": {"mean": safe_mean(hybrid), "median": safe_median(hybrid)},
        "allergenicity_ml": {"mean": safe_mean(ml), "median": safe_median(ml)},
        "allergenicity_merci": {"mean": safe_mean(merci), "median": safe_median(merci)},
        "allergenicity_blast": {"mean": safe_mean(blast), "median": safe_median(blast)}
    }

def parse_cluster_dir(directory):
    logging.debug(f"Starting parse_cluster_dir on directory: {directory}")
    clusters = defaultdict(list)
    for file in os.listdir(directory):
        if not file.endswith(".m8"):
            logging.debug(f"Skipping non-.m8 file in cluster dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing cluster file: {file}")
        with open(path) as f:
            for i, line in enumerate(f, start=1):
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    logging.debug(f"Skipping incomplete line {i} in {file}")
                    continue
                try:
                    percent_identity = float(parts[2])
                    clusters[parts[0]].append(percent_identity)
                except Exception as e:
                    logging.debug(f"Skipping bad percent identity at line {i} in {file}: {e}")
                    continue

    conservation_scores = [safe_mean(vals) for vals in clusters.values() if vals]

    logging.debug("Finished parse_cluster_dir")
    return {
        "cluster_conservation_score": {
            "mean": safe_mean(conservation_scores),
            "median": safe_median(conservation_scores)
        }
    }

def parse_popcov_dir(directory):
    logging.debug(f"Starting parse_popcov_dir on directory: {directory}")
    individuals, coverage = [], []

    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            logging.debug(f"Skipping non-txt file in popcoverage dir: {file}")
            continue
        path = os.path.join(directory, file)
        logging.debug(f"Parsing population coverage file: {file}")
        try:
            with open(path) as f:
                in_table = False
                for i, line in enumerate(f, start=1):
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
                            except Exception as e:
                                logging.debug(f"Skipping bad data at line {i} in {file}: {e}")
                                continue
        except Exception as e:
            logging.warning(f"Failed to parse population coverage file {file}: {e}")
            continue

    logging.debug("Finished parse_popcov_dir")
    return {
        "popcov_percent_individuals": {"mean": safe_mean(individuals), "median": safe_median(individuals)},
        "popcov_cumulative_coverage": {"mean": safe_mean(coverage), "median": safe_median(coverage)}
    }

# ----------------------------- Orchestration Functions -----------------------------

def extract_all_features(base_dir, eval_dir, threads=1):
    logging.info(f"Extracting features from base_dir: {base_dir} and eval_dir: {eval_dir} using {threads} threads")
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
                logging.info(f"Starting parser for: {name}")
                results[name] = future.result()
                logging.info(f"Completed parser for: {name}")
            except Exception as e:
                logging.error(f"Error parsing {name}: {e}")

    logging.info("Completed extracting all features")
    return results

def compare_ks(pos_features, rand_features):
    logging.info("Starting KS test comparison")
    results = []

    for feature, pos_data in pos_features.items():
        rand_data = rand_features.get(feature, {})

        if isinstance(pos_data, dict):
            for subfeature, val in pos_data.items():
                pos_vals = [val["mean"]]
                rand_vals = [rand_data.get(subfeature, {}).get("mean", 0)]

                try:
                    stat, pval = ks_2samp(pos_vals, rand_vals)
                except Exception as e:
                    logging.debug(f"KS test failed for feature {feature} subfeature {subfeature}: {e}")
                    stat, pval = None, None

                results.append({
                    "feature": feature,
                    "subfeature": subfeature,
                    "ks_statistic": stat,
                    "p_value": pval,
                    "positive_n": len(pos_vals),
                    "random_n": len(rand_vals)
                })

    logging.info("KS test comparison complete")
    return pd.DataFrame(results)

def write_features_by_feature(features, label, output_dir):
    logging.info(f"Writing features to disk for label {label} in {output_dir}")
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
        logging.debug(f"Wrote feature file: {filepath}")

# ----------------------------- Entry Point -----------------------------

def main(pathogen_dir, threads, verbose=False):
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

    raw_out_dir_pos = os.path.join("data", pathogen_dir, "raw_positive_features")
    raw_out_dir_rand = os.path.join("data", pathogen_dir, "raw_random_features")

    logger.info(f"Writing positive features to {raw_out_dir_pos}")
    write_features_by_feature(pos_features, "positive", raw_out_dir_pos)

    logger.info(f"Writing random features to {raw_out_dir_rand}")
    write_features_by_feature(rand_features, "random", raw_out_dir_rand)

    logger.info("Running KS test on features")
    result_df = compare_ks(pos_features, rand_features)
    logger.info("\n" + result_df.to_string(index=False))

    ks_out_path = os.path.join("data", pathogen_dir, "ks_test_results.csv")
    logger.info(f"Writing KS test results to {ks_out_path}")
    result_df.to_csv(ks_out_path, index=False)

    logger.info("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KS-test comparison of epitope vs. random features.")
    parser.add_argument("pathogen_dir", help="Pathogen directory name under data/")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use for parsing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging to file")
    args = parser.parse_args()

    main(args.pathogen_dir, args.threads, args.verbose)
