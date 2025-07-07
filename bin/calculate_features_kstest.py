import os
import json
import csv
import sys
import pandas as pd
from scipy.stats import ks_2samp
from statistics import mean, median
import collections
from collections import defaultdict
import logging
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------- Utility Functions -----------------------------

def safe_mean(lst):
    return mean(lst) if lst else 0.0

def safe_median(lst):
    return median(lst) if lst else 0.0

def init_logging(verbose=False, pathogen="unknown"):
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
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return f"{num:.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}P{suffix}"

# ----------------------------- Feature Parsers -----------------------------

def parse_bcell_dir(directory):
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
                            results.append({"feature": "bcell", "sub_feature": method, "value": val})
                        except Exception as e:
                            logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")
                else:
                    if method in ["chou-fasman", "emini", "karplus-schulz", "kolaskar-tongaonkar", "parker"]:
                        idx = headers.index("Score")
                        for i, row in enumerate(reader):
                            try:
                                val = float(row[idx])
                                results.append({"feature": "bcell", "sub_feature": method, "value": val})
                            except Exception as e:
                                logging.debug(f"Skipping row {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing B-cell file {file}: {e}")
    logging.info(f"Completed B-cell parsing with {len(results)} results")
    return results

def parse_mhc_dir(directory):
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
                                results.append({"feature": prefix, "sub_feature": "score", "value": float(row[idx_score])})
                                results.append({"feature": prefix, "sub_feature": "percentile", "value": float(row[idx_percentile])})
                                results.append({"feature": prefix, "sub_feature": "peptide_length", "value": len(row[idx_peptide])})
                            except Exception as e:
                                logging.debug(f"Skipping row {i} in {file} due to error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing MHC file {file}: {e}")
    logging.info(f"Completed MHC parsing with {len(results)} results")
    return results

def parse_signalp_dir(directory):
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
                            results.append({"feature": "signalp", "sub_feature": "prob_signalp", "value": float(parts[2])})
                            results.append({"feature": "signalp", "sub_feature": "prob_other", "value": float(parts[3])})
                        except Exception as e:
                            logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing SignalP file {file}: {e}")
    logging.info(f"Completed SignalP parsing with {len(results)} results")
    return results

def parse_targetp_dir(directory):
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
                            results.append({"feature": "targetp", "sub_feature": "prob_noTP", "value": float(parts[2])})
                            results.append({"feature": "targetp", "sub_feature": "prob_SP", "value": float(parts[3])})
                            results.append({"feature": "targetp", "sub_feature": "prob_mTP", "value": float(parts[4])})
                        except Exception as e:
                            logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing TargetP file {file}: {e}")
    logging.info(f"Completed TargetP parsing with {len(results)} results")
    return results

def parse_allergenicity_dir(directory):
    result = {}

    if not os.path.isdir(directory):
        logging.warning(f"Allergenicity directory not found: {directory}")
        return result

    for filename in os.listdir(directory):
        if not filename.endswith(".json"):
            continue

        filepath = os.path.join(directory, filename)

        try:
            with open(filepath, "r") as f:
                data = json.load(f)
        except Exception as e:
            logging.warning(f"Failed to read or parse Allergenicity JSON: {filepath} - {e}")
            continue

        allergen_score = data.get("AllergenicityScore")
        if allergen_score is not None:
            result.setdefault("allergenicity", {}).setdefault("score", []).append(allergen_score)

    return result

def parse_cluster_dir(directory):
    logging.info(f"Parsing cluster conservation dir {directory}")
    clusters = defaultdict(list)
    try:
        files = [f for f in os.listdir(directory) if f.endswith(".m8")]
        logging.debug(f"Found {len(files)} cluster files")
    except Exception as e:
        logging.error(f"Failed listing cluster directory {directory}: {e}")
        return []

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
                    try:
                        percent_identity = float(parts[2])
                        clusters[query_id].append(percent_identity)
                    except ValueError as e:
                        logging.debug(f"Skipping line {i} in {file} due to conversion error: {e}")
        except Exception as e:
            logging.error(f"Failed parsing cluster file {file}: {e}")

    results = []
    cluster_scores = []
    for cluster_id, identities in clusters.items():
        if identities:
            conservation_score = sum(identities) / len(identities)
            cluster_scores.append(conservation_score)
            results.append({"feature": "cluster_conservation", "sub_feature": cluster_id, "value": conservation_score})

    if cluster_scores:
        results.append({"feature": "cluster_conservation", "sub_feature": "mean", "value": mean(cluster_scores)})
        results.append({"feature": "cluster_conservation", "sub_feature": "median", "value": median(cluster_scores)})

    logging.info(f"Completed cluster parsing with {len(results)} results")
    return results

def parse_popcov_dir(directory):
    result = {}

    if not os.path.isdir(directory):
        logging.warning(f"PopCov directory not found: {directory}")
        return result

    for filename in os.listdir(directory):
        if not filename.endswith(".json"):
            continue

        filepath = os.path.join(directory, filename)

        try:
            with open(filepath, "r") as f:
                data = json.load(f)
        except Exception as e:
            logging.warning(f"Failed to read or parse PopCov JSON: {filepath} - {e}")
            continue

        popcov_scores = data.get("population_coverage", {})
        for region, score in popcov_scores.items():
            result.setdefault("popcov", {}).setdefault(region, []).append(score)

    return result

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
        "popcov": lambda: parse_popcov_dir(os.path.join(eval_dir, "popcoverage"))
    }

    results = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_name = {executor.submit(parser): name for name, parser in parsers.items()}
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                results.extend(future.result())
                logging.info(f"{name} feature extraction complete")
            except Exception as e:
                logging.error(f"{name} feature extraction failed: {e}")
    return results

def compare_ks(pos_features, rand_features):
    logging.info("Starting KS test comparison")
    results = []

    for feature, pos_subfeatures in pos_features.items():
        rand_subfeatures = rand_features.get(feature, {})

        if isinstance(pos_subfeatures, dict):
            for subfeature, pos_vals in pos_subfeatures.items():
                rand_vals = rand_subfeatures.get(subfeature, [])

                if not pos_vals or not rand_vals:
                    continue  # Skip if either side has no data

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

    if write_raw:
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
    parser.add_argument("--write-raw", action="store_true", help="Write raw feature data to disk (optional, large output)")
    args = parser.parse_args()

    main(args.pathogen_dir, args.threads, args.verbose, args.write_raw)
