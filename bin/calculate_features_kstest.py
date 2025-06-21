import os
import json
import csv
import re
import pandas as pd
from scipy.stats import ks_2samp

def parse_bcell_dir(directory):
    scores = {name: [] for name in [
        "Chou-Fasman", "Emini", "Karplus-Schulz",
        "Kolaskar-Tongaonkar", "Parker", "Bepipred"
    ]}
    for file in os.listdir(directory):
        if not file.endswith(".csv"):
            continue
        method = os.path.basename(file).split("_")[-1].replace(".csv", "")
        path = os.path.join(directory, file)
        with open(path) as f:
            reader = csv.reader(f)
            headers = next(reader)
            if method.lower() == "bepipred":
                idx = headers.index("Score")
                for row in reader:
                    try:
                        scores["Bepipred"].append(float(row[idx]))
                    except:
                        continue
            elif method in scores:
                idx = headers.index("Score")
                for row in reader:
                    try:
                        scores[method].append(float(row[idx]))
                    except:
                        continue
    return scores

def parse_mhc_dir(directory):
    scores = {"score": [], "percentile": [], "peptide_length": []}
    for file in os.listdir(directory):
        if not file.endswith(".json"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            try:
                data = json.load(f)
                for result in data.get("results", []):
                    if result.get("type") == "peptide_table":
                        cols = result["table_columns"]
                        table = result["table_data"]
                        idx_score = cols.index("score")
                        idx_percentile = cols.index("percentile")
                        idx_peptide = cols.index("peptide")
                        for row in table:
                            try:
                                scores["score"].append(float(row[idx_score]))
                                scores["percentile"].append(float(row[idx_percentile]))
                                scores["peptide_length"].append(len(row[idx_peptide]))
                            except:
                                continue
            except:
                continue
    return scores

def parse_signalp_dir(directory):
    features = {"predicted_feature": [], "prob_signalp": [], "prob_other": []}
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
                        features["predicted_feature"].append(parts[1])
                        features["prob_signalp"].append(float(parts[2]))
                        features["prob_other"].append(float(parts[3]))
                    except:
                        continue
    return features

def parse_targetp_dir(directory):
    features = {"predicted_feature": [], "prob_noTP": [], "prob_SP": [], "prob_mTP": []}
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
                        features["predicted_feature"].append(parts[1])
                        features["prob_noTP"].append(float(parts[2]))
                        features["prob_SP"].append(float(parts[3]))
                        features["prob_mTP"].append(float(parts[4]))
                    except:
                        continue
    return features

def parse_tmhmm_dir(directory):
    features = {
        "predicted_thms": [],
        "expected_aa_num": [],
        "expected_aa_num_60": [],
        "prob_n_in": [],
        "has_nterm_signal_seq": []
    }
    for file in os.listdir(directory):
        if not file.endswith(".txt"):
            continue
        path = os.path.join(directory, file)
        with open(path) as f:
            content = f.read()
            matches = {
                "predicted_thms": re.search(r"Number of predicted TMHs:\s+(\d+)", content),
                "expected_aa_num": re.search(r"Exp number of AAs in TMHs:\s+([\d.]+)", content),
                "expected_aa_num_60": re.search(r"first 60 AAs:\s+([\d.]+)", content),
                "prob_n_in": re.search(r"Total prob of N-in:\s+([\d.]+)", content),
                "has_nterm_signal_seq": "POSSIBLE N-term signal sequence" in content
            }
            try:
                features["predicted_thms"].append(int(matches["predicted_thms"].group(1)))
                features["expected_aa_num"].append(float(matches["expected_aa_num"].group(1)))
                features["expected_aa_num_60"].append(float(matches["expected_aa_num_60"].group(1)))
                features["prob_n_in"].append(float(matches["prob_n_in"].group(1)))
                features["has_nterm_signal_seq"].append(int(matches["has_nterm_signal_seq"]))
            except:
                continue
    return features

def extract_all_features(base_dir):
    return {
        "bcell": parse_bcell_dir(os.path.join(base_dir, "bcell")),
        "mhci": parse_mhc_dir(os.path.join(base_dir, "mhci")),
        "mhcii": parse_mhc_dir(os.path.join(base_dir, "mhcii")),
        "signalp": parse_signalp_dir(os.path.join(base_dir, "signalp")),
        "targetp": parse_targetp_dir(os.path.join(base_dir, "targetp")),
        "tmhmm": parse_tmhmm_dir(os.path.join(base_dir, "tmhmm")),
    }

def compare_ks(pos_features, rand_features):
    results = []

    for feature in pos_features:
        pos_data = pos_features[feature]
        rand_data = rand_features.get(feature, {})

        if isinstance(pos_data, dict):
            for subfeature in pos_data:
                pos_vals = pos_data.get(subfeature, [])
                rand_vals = rand_data.get(subfeature, [])

                if pos_vals and rand_vals:
                    stat, pval = ks_2samp(pos_vals, rand_vals)
                else:
                    stat, pval = None, None

                results.append({
                    "feature": feature,
                    "subfeature": subfeature,
                    "ks_statistic": stat,
                    "p_value": pval,
                    "positive_n": len(pos_vals),
                    "random_n": len(rand_vals)
                })

        elif isinstance(pos_data, list):
            if pos_data and rand_data:
                stat, pval = ks_2samp(pos_data, rand_data)
            else:
                stat, pval = None, None

            results.append({
                "feature": feature,
                "subfeature": None,
                "ks_statistic": stat,
                "p_value": pval,
                "positive_n": len(pos_data),
                "random_n": len(rand_data)
            })

    return pd.DataFrame(results)

def write_features_to_csv(features, label, filepath):
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "feature", "subfeature", "value"])
        writer.writeheader()
        for feature, subdata in features.items():
            if isinstance(subdata, dict):
                for subfeature, values in subdata.items():
                    for val in values:
                        writer.writerow({
                            "label": label,
                            "feature": feature,
                            "subfeature": subfeature,
                            "value": val
                        })
            elif isinstance(subdata, list):
                for val in subdata:
                    writer.writerow({
                        "label": label,
                        "feature": feature,
                        "subfeature": None,
                        "value": val
                    })

import sys

def sizeof_fmt(num, suffix="B"):
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return f"{num:.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}P{suffix}"

def main(pathogen_dir):
    pos_dir = os.path.join("data", pathogen_dir, "epitope_outputs")
    rand_dir = os.path.join("data", pathogen_dir, "random_analysis")

    print(f"Extracting features for: {pathogen_dir}")
    pos_features = extract_all_features(pos_dir)
    rand_features = extract_all_features(rand_dir)

    print("Estimating memory usage...")
    pos_size = sys.getsizeof(pos_features)
    rand_size = sys.getsizeof(rand_features)
    print(f"Positive features: {sizeof_fmt(pos_size)}")
    print(f"Random features: {sizeof_fmt(rand_size)}")


    raw_out_path_pos = os.path.join("data", pathogen_dir, "raw_positive_features.csv")
    raw_out_path_rand = os.path.join("data", pathogen_dir, "raw_random_features.csv")

    print("Writing raw feature data to CSV...")
    write_features_to_csv(pos_features, "positive", raw_out_path_pos)
    write_features_to_csv(rand_features, "random", raw_out_path_rand)
    print(f"Positive features -> {raw_out_path_pos}")
    print(f"Random features  -> {raw_out_path_rand}")

    combined_features = {**pos_features, **rand_features}
    combined_out_path = os.path.join("data", pathogen_dir, "combined_features.csv")
    print("Writing combined feature data to CSV...")
    write_features_to_csv(combined_features, "combined", combined_out_path)
    print(f"Combined features -> {combined_out_path}")


    # KS test
    print("Performing KS test...")
    result_df = compare_ks(pos_features, rand_features)
    print(result_df.to_string(index=False))

    ks_out_path = os.path.join("data", pathogen_dir, "ks_test_results.csv")
    result_df.to_csv(ks_out_path, index=False)
    print(f"\nKS test results saved to: {ks_out_path}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="KS-test comparison of epitope vs. random features.")
    parser.add_argument("pathogen_dir", help="Pathogen directory name under data/")
    args = parser.parse_args()
    main(args.pathogen_dir)
