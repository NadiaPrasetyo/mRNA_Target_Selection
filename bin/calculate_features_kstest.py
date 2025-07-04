import os
import json
import csv
import re
import pandas as pd
from scipy.stats import ks_2samp
from statistics import mean, median
from collections import defaultdict
from Bio import SeqIO

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

def extract_all_features(base_dir):
    return {
        "bcell": parse_bcell_dir(os.path.join(base_dir, "bcell")),
        "mhci": parse_mhc_dir(os.path.join(base_dir, "mhci")),
        "mhcii": parse_mhc_dir(os.path.join(base_dir, "mhcii")),
        "signalp": parse_signalp_dir(os.path.join(base_dir, "signalp")),
        "targetp": parse_targetp_dir(os.path.join(base_dir, "targetp")),
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

def extract_evaluation_features(base_dir):
    return {
        "allergenicity": parse_allergenicity_dir(os.path.join(base_dir, "allergenicity")),
        "cluster": parse_cluster_dir(os.path.join(base_dir, "cluster")),
        "popcoverage": parse_popcov_dir(os.path.join(base_dir, "popcoverage")),
    }

def parse_allergenicity_dir(directory):
# input files:{file_stem}_algpred.csv AND associated {file_stem}.fasta in fasta_inputs directory that is in the same directory as the directory i.e. parent/directory, parent/fasta_inputs
# fields of csv: Subject,ML Score,MERCI Score,BLAST Score,Hybrid Score,Prediction
# to do: parse the fasta file to get the sequence of predicted epitopes AND not predicted epitopes
# if not predicted (i.e. the sequence subject is found in the fasta but is not found in the csv), the sequence is Non-Allergen (prediction), leave the hybrid score as 0
# to be extracted: Subject,Hybrid Score, prediction, sequence
    fasta_dir = os.path.join(os.path.dirname(directory), "fasta_inputs")
    results = {
        "hybrid_score": [],
        "is_allergen": [],
        "percent_allergenicity": []
    }

    all_sequences = {}
    for file in os.listdir(fasta_dir):
        if not file.endswith(".fasta") and not file.endswith(".fa"):
            continue
        for record in SeqIO.parse(os.path.join(fasta_dir, file), "fasta"):
            all_sequences[record.id] = str(record.seq)

    for file in os.listdir(directory):
        if not file.endswith("_algpred.csv"):
            continue
        csv_path = os.path.join(directory, file)
        detected = set()
        with open(csv_path, newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                subj = row["Subject"]
                results["hybrid_score"].append(float(row.get("Hybrid Score", 0)))
                results["is_allergen"].append(1)
                detected.add(subj)
        
        # Get corresponding fasta file
        stem = file.replace("_algpred.csv", "")
        fasta_path = os.path.join(fasta_dir, stem + ".fasta")
        if os.path.exists(fasta_path):
            for record in SeqIO.parse(fasta_path, "fasta"):
                if record.id not in detected:
                    results["hybrid_score"].append(0.0)
                    results["is_allergen"].append(0)

    # Allergenicity percentile
    if len(results["is_allergen"]) > 0:
        percent_allergenicity = sum(results["is_allergen"]) / len(results["is_allergen"])
        results["percent_allergenicity"] = [percent_allergenicity] * len(results["is_allergen"])
    
    return {"hybrid_score": results["hybrid_score"], "percent_allergenicity": results["percent_allergenicity"]}


def parse_cluster_dir(directory):
# input files: {antigen accession}_combined_scores.m8
# format:
# Column	Content
# 0	Query sequence ID
# 1	Subject (database) sequence ID
# 2	Percent Identity
# 3	Alignment Length
# 4	Number of gaps
# 5	Number of mismatches
# 6	Start on the query sequence
# 7	End on the query sequence
# 8	Start on the database sequence
# 9	End on the database sequence
# 10	E value - the expectation that this alignment is random given the length of the sequence and length of the database
# 11	bit score - the score of the alignment itself
# to be extracted: query sequence ID, Subject sequence ID, Percent Identity, Alignment Length, E value, bit score
    """
    Parse .m8 files to compute per-cluster conservation scores (mean and median percent identity).
    Each unique query sequence is considered a separate cluster.
    Returns:
        {
            "cluster_conservation_mean": [mean1, mean2, ...],
            "cluster_conservation_median": [median1, median2, ...]
        }
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
                query_id = parts[0]
                try:
                    percent_identity = float(parts[2])
                    clusters[query_id].append(percent_identity)
                except ValueError:
                    continue

    # Aggregate percent identity scores per cluster
    conservation_means = []
    conservation_medians = []
    for cluster_id, identities in clusters.items():
        if identities:
            conservation_means.append(mean(identities))
            conservation_medians.append(median(identities))

    return {
        "cluster_conservation_mean": conservation_means,
        "cluster_conservation_median": conservation_medians
    }

def parse_popcov_dir(directory):
# files: {file_stem}.csv and associated {file_stem}.fasta in popcov_inputs directory that is in the same directory as the directory i.e. parent/directory, parent/popcov_inputs
# format: population/area	epitope_hits	percent_individuals	cumulative_coverage
# skip the first few lines until the header: population/area	epitope_hits	percent_individuals	cumulative_coverage
# skip empty lines and lines with more than 4 columns
# to do: parse the fasta file to get the sequence of predicted epitopes 0 for the first, 1, 2, 3, etc. for the rest
# to be extracted: sequence, percent_individuals, cumulative_coverage
    fasta_dir = os.path.join(os.path.dirname(directory), "popcov_inputs")
    cumulative_coverage = []
    percent_individuals = []

    for file in os.listdir(directory):
        if not file.endswith(".csv"):
            continue
        path = os.path.join(directory, file)

        with open(path) as f:
            for line in f:
                if line.strip().startswith("population/area"):
                    break  # Skip lines until header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) != 4 or not parts[0] or not parts[1]:
                    continue
                try:
                    percent_individuals.append(float(parts[2]))
                    cumulative_coverage.append(float(parts[3]))
                except:
                    continue

    return {
        "popcov_cumulative_mean": [mean(cumulative_coverage)] * len(cumulative_coverage),
        "popcov_cumulative_median": [median(cumulative_coverage)] * len(cumulative_coverage)
    }

def main(pathogen_dir):
    pos_dir = os.path.join("data", pathogen_dir, "epitope_outputs")
    pos_eval_dir = os.path.join("data", pathogen_dir, "evaluation_outputs")
    rand_dir = os.path.join("data", pathogen_dir, "random_analysis")
    rand_eval_dir = os.path.join("data", pathogen_dir, "random_evaluation")

    print(f"Extracting features for: {pathogen_dir}")
    pos_features = extract_all_features(pos_dir)
    rand_features = extract_all_features(rand_dir)
    pos_eval_features = extract_evaluation_features(pos_eval_dir)
    rand_eval_features = extract_evaluation_features(rand_eval_dir)

    # Merge evaluation features
    pos_features.update(pos_eval_features)
    rand_features.update(rand_eval_features)

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
