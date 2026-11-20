#!/usr/bin/env python3
"""
Random Forest antigen-prediction pipeline -- NEGATIVE-SUBSAMPLING SENSITIVITY CHECK.

This is a copy of the original pipeline with ONE change: a random fraction (default 75%) of
the available negative (non-antigen) training samples is kept; the rest are discarded.
Everything else (features, model, MCC-optimised OOB threshold, validation, bootstrap CIs,
feature importances, target prediction) is identical to the original, so the outputs can be
compared directly with the original run to show how robust the model is to a modest
perturbation of the negative set.

Where the subsampling happens (--subsample-stage):
  - pool  (default): negatives are subsampled from the combined training-bacteria pool BEFORE
                     the train/validation split. Both the training and the validation sets
                     therefore contain fewer negatives (validation prevalence goes up, so
                     precision / PR-AUC are not directly comparable with the original; compare
                     sensitivity, specificity, MCC, ROC-AUC, or use 'train' below).
  - train          : the train/validation split is done exactly as in the original (same
                     random_state, same validation samples), and only the negatives in the
                     TRAINING split are subsampled. The validation set is identical to the
                     original run, so ALL metrics are directly comparable.

The sampling is controlled by --seed (independent from the fixed split/model seeds), so you can
repeat the check with several seeds, e.g.:

    for s in 1 2 3 4 5; do
        python rf_antigen_pipeline_neg75.py --input-dir E.coli S.aureus ... \
            --seed $s --output-dir ./results_neg75/seed_$s
    done

Reads raw CSVs from <base_dir>/<bacterium>/raw_data/*_raw_data.csv
Each raw CSV is expected to contain at least these columns:
 - accession
 - feature
 - subfeature
 - value
 - label (required for training bacteria; not needed for the target bacterium)

Pipeline:
 1) Load and concatenate data per bacterium
 2) Pivot to accession x feature_subfeature matrix (mean of replicates)
 3) Align feature columns across all bacteria (union of features, fill missing with 0)
 4) Train RandomForest on all bacteria EXCEPT the target bacterium (default: S.aureus),
    after randomly keeping only --neg-fraction of the negative samples
 5) Choose the classification threshold that maximises MCC on the out-of-bag (OOB)
    predictions of the training split (or use a fixed --threshold if given), then validate
    on a hold-out set drawn from the training bacteria at that threshold
    (full metrics + bootstrap confidence intervals + confusion matrix)
 6) Compute feature importances (feature usefulness)
 7) Predict probabilities on the target bacterium (prediction only, no metrics because
    the target bacterium is unlabelled) and save CSV sorted by probability (desc)

Outputs (saved in output_dir):
 - <target_bacterium>_predictions.csv  (accession, prob_antigen, pred_label, protein/gene names)
 - <target_bacterium>_features_with_probs.csv (full feature matrix with probabilities)
 - feature_importances.csv (feature, importance, rank)
 - validation_metrics.csv  (metric, estimate, CI lower, CI upper)
 - threshold_mcc_curve.csv (OOB MCC at each candidate threshold; only when threshold is MCC-optimised)
 - negatives_retained.csv  (accession and whether each negative was kept by the subsampling)
 - model_report.txt        (train/validation sizes, subsampling details, confusion matrix,
                            all validation metrics with CIs, target prediction summary)

UniProt name cache: protein/gene names are cached in --name-cache (default
./uniprot_name_cache.json, NOT inside output_dir, so all runs of a loop share it). Only the
first run queries UniProt; later runs read the names from the cache.
"""

import os
import glob
import argparse
import json
import tempfile
import logging
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import requests
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_auc_score, average_precision_score, confusion_matrix
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split

# ----------------------
# Logging setup
# ----------------------
def setup_logging(verbose: bool) -> None:
    """
    Setup logging configuration.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

# ----------------------
# Data loading utils
# ----------------------
def load_bacterium_data(base_dir: str, bacterium: str) -> pd.DataFrame:
    """Load and concatenate all raw data CSVs for a bacterium.
    Args:
        base_dir (str): Base directory containing bacterium subdirectories.
        bacterium (str): Name of the bacterium (subdirectory name).
    Returns:
        pd.DataFrame: Concatenated DataFrame of all raw data for the bacterium.
    """
    folder = os.path.join(base_dir, bacterium, "raw_data")
    logging.info(f"Loading data for bacterium '{bacterium}' from: {folder}")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    if not files:
        logging.warning(f"No files found for {bacterium} at {folder}")
        return pd.DataFrame()
    dfs = []
    for f in files:
        logging.debug(f"Reading {f}")
        try:
            df = pd.read_csv(f)
        except Exception as e:
            logging.error(f"Failed to read {f}: {e}")
            continue
        df["bacterium"] = bacterium
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)
    # Strip whitespace from column names
    out.columns = [c.strip() for c in out.columns]
    return out

# ----------------------
# Preprocess & pivot
# ----------------------
def preprocess_and_pivot(df: pd.DataFrame):
    """
    Build feature matrix:
      - feature_subfeature = feature_subfeature
      - pivot so each accession is a row, each column is a feature_subfeature (mean)

    Returns (X_df, labels_series) where:
      - X_df.index = accession
      - labels_series.index = accession (if label present, else None), dtype object

    Filters out accessions with <=4 characters.
    """
    required_cols = {"accession", "feature", "subfeature", "value"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Input dataframe missing required cols: {missing_cols}")

    df = df.copy()

    # Filter accessions with more than 4 characters
    df = df[df["accession"].astype(str).str.len() > 4]
    if df.empty:
        logging.warning("No accessions with more than 4 characters found.")
        return pd.DataFrame(), pd.Series(dtype=object)

    # Compose feature ID
    df["feature_subfeature"] = df["feature"].astype(str) + "_" + df["subfeature"].astype(str)

    # Ensure numeric values
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["value"].notna()]

    # Aggregate mean across replicates
    agg = (
        df.groupby(["accession", "feature_subfeature"], observed=True)["value"]
          .mean()
          .reset_index()
    )
    # Pivot to wide format
    X = agg.pivot(index="accession", columns="feature_subfeature", values="value")
    X = X.fillna(0.0)  # missing features -> 0
    # Labels (if present)
    labels = None
    if "label" in df.columns:
        labels = (
            df.groupby("accession", observed=True)["label"]
              .agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0])
        )
        # Align labels to pivoted X
        labels = labels.loc[X.index]

    return X, labels

# ----------------------
# Label handling heuristic -> binary
# ----------------------
def make_binary_labels(labels: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    Convert labels to binary 0/1.
    Args:
      labels (pd.Series): Input labels.
    Heuristic:
      - If labels are numeric and only {0,1}, use directly
      - If any label contains 'antigen'/'positive'/'pos' (case-ins), map those to 1
      - Else if two classes, map the second label to 1 (LabelEncoder) and warn
      - Else raise
    Returns: (y_array, label_encoder) where label_encoder.inverse_transform maps 0/1 back.
    """
    if labels is None:
        raise ValueError("No labels provided to make_binary_labels()")

    labs = labels.astype(str).copy()
    unique = sorted(set(labs))
    logging.info(f"Unique label values found: {unique}")

    # Try numeric
    try:
        numeric = pd.to_numeric(labels.dropna().unique())
        unique_numeric = set(numeric.tolist())
        if unique_numeric.issubset({0, 1}):
            y = labels.astype(int).to_numpy()
            le = LabelEncoder()
            le.fit([0, 1])
            return y, le
    except Exception:
        pass

    # check for positive tokens
    pos_tokens = ("antigen", "positive", "pos", "1", "true", "yes")
    positive_vals = [u for u in unique if any(tok in u.lower() for tok in pos_tokens)]
    if positive_vals:
        mapping = {}
        for u in unique:
            mapping[u] = 1 if u in positive_vals else 0
        y = labels.map(mapping).astype(int).to_numpy()
        le = LabelEncoder()
        le.fit(["0", "1"])
        return y, le

    # fallback: two classes -> map second encoded class to 1
    if len(unique) == 2:
        le = LabelEncoder()
        le.fit(unique)
        y_encoded = le.transform(labs)
        # we will map encoded==1 to positive (so second class -> 1)
        return y_encoded, le

    raise ValueError(f"Unable to convert labels to binary. Found >2 classes: {unique}")

# ----------------------
# Negative subsampling (sensitivity check)
# ----------------------
def subsample_negatives(X: pd.DataFrame, y: np.ndarray, fraction: float,
                        seed: int) -> Tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """Randomly keep `fraction` of the negative (y == 0) samples; all positives are kept.

    Sampling is without replacement and uses positional indices, so duplicated index values
    are not a problem. Row order of the retained samples is preserved.

    Returns (X_sub, y_sub, retained_df) where retained_df lists every original negative
    (accession, kept True/False) for traceability.
    """
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"neg fraction must be in (0, 1], got {fraction}")

    rng = np.random.default_rng(seed)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    n_keep = max(1, int(round(fraction * len(neg_idx))))
    keep_neg = rng.choice(neg_idx, size=n_keep, replace=False)

    keep = np.sort(np.concatenate([pos_idx, keep_neg]))
    kept_set = set(keep_neg.tolist())
    retained_df = pd.DataFrame({
        "accession": X.index[neg_idx],
        "kept": [i in kept_set for i in neg_idx],
    })
    return X.iloc[keep], y[keep], retained_df

# ----------------------
# Train Random Forest
# ----------------------
def train_rf(X_train: pd.DataFrame, y_train: np.ndarray, n_estimators: int = 500, random_state: int = 42) -> RandomForestClassifier:
    logging.info(f"Training RandomForest on {X_train.shape[0]} samples and {X_train.shape[1]} features...")
    rf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=-1, random_state=random_state, class_weight='balanced', oob_score=True)
    rf.fit(X_train, y_train)
    return rf

# ----------------------
# Evaluation metrics (validation set only)
# ----------------------
METRIC_ORDER = [
    "sensitivity_recall", "specificity", "precision", "f1", "mcc",
    "accuracy", "roc_auc", "pr_auc",
]
METRIC_LABELS = {
    "sensitivity_recall": "Sensitivity / Recall (TPR)",
    "specificity": "Specificity (TNR)",
    "precision": "Precision (PPV)",
    "f1": "F1 score",
    "mcc": "Matthews correlation coefficient",
    "accuracy": "Accuracy",
    "roc_auc": "ROC-AUC",
    "pr_auc": "PR-AUC (average precision)",
}


def _safe_div(num: float, den: float) -> float:
    return num / den if den > 0 else np.nan


def point_metrics(y_true: np.ndarray, probs: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute all metrics at a fixed classification threshold (AUCs are threshold-free).

    Undefined values (e.g. precision when nothing is predicted positive) are returned as NaN.
    """
    preds = (probs >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, preds, labels=[0, 1]).ravel()
    tn, fp, fn, tp = float(tn), float(fp), float(fn), float(tp)

    mcc_den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    metrics = {
        "tn": tn, "fp": fp, "fn": fn, "tp": tp,
        "sensitivity_recall": _safe_div(tp, tp + fn),
        "specificity": _safe_div(tn, tn + fp),
        "precision": _safe_div(tp, tp + fp),
        "f1": _safe_div(2 * tp, 2 * tp + fp + fn),
        "mcc": _safe_div(tp * tn - fp * fn, mcc_den),
        "accuracy": _safe_div(tp + tn, tp + tn + fp + fn),
        "roc_auc": np.nan,
        "pr_auc": np.nan,
    }
    if len(np.unique(y_true)) == 2:
        metrics["roc_auc"] = roc_auc_score(y_true, probs)
        metrics["pr_auc"] = average_precision_score(y_true, probs)
    return metrics


def select_threshold_mcc(y_true: np.ndarray, probs: np.ndarray,
                         grid: Optional[np.ndarray] = None) -> Tuple[float, float, pd.DataFrame]:
    """Pick the probability threshold that maximises MCC.

    Should be run on predictions that were NOT used to fit the model (e.g. OOB predictions),
    never on the validation set that is later used for reporting.
    If several thresholds tie for the best MCC, the median of the tied thresholds is used.

    Returns (best_threshold, best_mcc, curve_df[threshold, mcc]).
    """
    if grid is None:
        grid = np.round(np.arange(0.01, 1.0, 0.01), 2)
    mccs = []
    for t in grid:
        preds = probs >= t
        tp = float(np.sum(preds & (y_true == 1)))
        fp = float(np.sum(preds & (y_true == 0)))
        fn = float(np.sum(~preds & (y_true == 1)))
        tn = float(np.sum(~preds & (y_true == 0)))
        den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mccs.append((tp * tn - fp * fn) / den if den > 0 else 0.0)
    mccs = np.asarray(mccs)
    best = mccs.max()
    best_thr = float(np.median(grid[np.isclose(mccs, best)]))
    return best_thr, float(best), pd.DataFrame({"threshold": grid, "mcc": mccs})


def bootstrap_ci(y_true: np.ndarray, probs: np.ndarray, threshold: float,
                 n_boot: int = 1000, ci_level: float = 0.95,
                 random_state: int = 42) -> Dict[str, Tuple[float, float]]:
    """Percentile bootstrap confidence intervals (resampling validation samples with replacement).

    Resamples containing only one class are skipped (AUCs undefined).
    """
    rng = np.random.default_rng(random_state)
    n = len(y_true)
    draws = {m: [] for m in METRIC_ORDER}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y_true[idx]
        if len(np.unique(yb)) < 2:
            continue
        m = point_metrics(yb, probs[idx], threshold)
        for k in METRIC_ORDER:
            draws[k].append(m[k])

    alpha = (1.0 - ci_level) / 2.0
    cis = {}
    for k, vals in draws.items():
        arr = np.asarray(vals, dtype=float)
        if arr.size == 0 or np.all(np.isnan(arr)):
            cis[k] = (np.nan, np.nan)
        else:
            cis[k] = (float(np.nanpercentile(arr, 100 * alpha)),
                      float(np.nanpercentile(arr, 100 * (1 - alpha))))
    return cis


def _fmt(x: float) -> str:
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.4f}"


def format_metrics_report(metrics: Dict[str, float], cis: Dict[str, Tuple[float, float]],
                          threshold: float, ci_level: float, n_boot: int,
                          n_pos: int, n_neg: int) -> str:
    """Build the text block for the validation section of the report."""
    tn, fp, fn, tp = (int(metrics[k]) for k in ("tn", "fp", "fn", "tp"))
    prevalence = n_pos / (n_pos + n_neg) if (n_pos + n_neg) else np.nan
    lines = [
        "=" * 70,
        "VALIDATION METRICS (hold-out set from training bacteria)",
        "=" * 70,
        f"Classification threshold: P(antigen) >= {threshold}",
        f"Validation samples: {n_pos + n_neg} (pos={n_pos}, neg={n_neg}, prevalence={_fmt(prevalence)})",
        f"Confidence intervals: {int(ci_level * 100)}% percentile bootstrap, {n_boot} resamples",
        "",
        f"Confusion matrix at threshold {threshold} (rows = true, columns = predicted):",
        f"{'':>18}{'Pred neg':>12}{'Pred pos':>12}",
        f"{'True neg':>18}{tn:>12d}{fp:>12d}   (TN, FP)",
        f"{'True pos':>18}{fn:>12d}{tp:>12d}   (FN, TP)",
        "",
        f"{'Metric':<36}{'Estimate':>10}   {int(ci_level * 100)}% CI",
        "-" * 70,
    ]
    for k in METRIC_ORDER:
        lo, hi = cis[k]
        lines.append(f"{METRIC_LABELS[k]:<36}{_fmt(metrics[k]):>10}   [{_fmt(lo)}, {_fmt(hi)}]")
    lines += [
        "",
        f"Note: ROC-AUC measures ranking quality across all thresholds; it is not precision or recall.",
        f"PR-AUC should be compared with the no-skill baseline (= prevalence = {_fmt(prevalence)}).",
        "",
    ]
    return "\n".join(lines)

# ----------------------
# Accession to Name mapping (optional)
# ----------------------

def _load_name_cache(path: Optional[str]) -> Dict[str, dict]:
    """Load the on-disk UniProt name cache ({accession: {protein_names, gene_names, found}})."""
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as e:
        logging.warning(f"Could not read name cache {path} ({e}); starting with an empty cache.")
        return {}


def _save_name_cache(path: Optional[str], new_entries: Dict[str, dict]) -> None:
    """Merge new entries into the on-disk cache and write it atomically.

    The file is re-read just before writing, so runs executing in parallel do not discard
    each other's entries (worst case: the same accession is fetched twice).
    """
    if not path or not new_entries:
        return
    merged = _load_name_cache(path)
    merged.update(new_entries)
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=folder, suffix=".tmp")
    with os.fdopen(fd, "w") as fh:
        json.dump(merged, fh)
    os.replace(tmp, path)  # atomic: readers never see a half-written file


def load_accession_name_and_gene_fetch(df: pd.DataFrame, cache_path: Optional[str] = None) -> pd.DataFrame:
    """
    Fetch protein names and gene names from UniProt for given accessions.

    If 'protein_names' (or legacy 'name') or 'gene_names' columns exist, values are preserved
    and only missing entries are fetched. Fetched results (including "accession not found")
    are stored in a persistent JSON cache, so repeated runs only query UniProt for accessions
    never seen before. Failed requests (network errors, timeouts) are NOT cached and will be
    retried next time.

    Args:
        df (pd.DataFrame): DataFrame containing at least an 'accession' column,
            and optionally 'protein_names'/'name' and/or 'gene_names' columns.
        cache_path (Optional[str]): Path of the JSON cache file. None/empty disables caching.

    Returns:
        pd.DataFrame: DataFrame indexed by 'accession' with columns ['protein_names', 'gene_names'].
    """
    if "accession" not in df.columns:
        raise ValueError("Input DataFrame must contain an 'accession' column.")

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    df = df.copy()

    if "name" in df.columns and "protein_names" not in df.columns:
        df = df.rename(columns={"name": "protein_names"})

    # Ensure required columns exist
    for col in ["protein_names", "gene_names"]:
        if col not in df.columns:
            df[col] = None

    headers = {"accept": "application/json"}

    cache = _load_name_cache(cache_path)
    new_entries: Dict[str, dict] = {}
    n_hits = 0
    logging.info(f"UniProt name cache: {len(cache)} entries loaded"
                 f"{f' from {cache_path}' if cache_path else ' (caching disabled)'}")

    for accession in df["accession"].unique():
        row = df.loc[df["accession"] == accession].iloc[0]

        if pd.notna(row["protein_names"]) and pd.notna(row["gene_names"]):
            continue  # Skip already known entries

        if accession in cache:
            hit = cache[accession]
            df.loc[df["accession"] == accession, "protein_names"] = hit.get("protein_names")
            df.loc[df["accession"] == accession, "gene_names"] = hit.get("gene_names")
            n_hits += 1
            continue

        params = {
            "query": f"accession:{accession}",
            "fields": "accession,protein_name,gene_names",
            "format": "json",
            "size": 50  # get multiple entries if available
        }

        try:
            response = requests.get(base_url, headers=headers, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()

            if not data.get("results"):
                logging.warning(f"No results found for accession {accession}")
                new_entries[accession] = {"protein_names": None, "gene_names": None, "found": False}
                continue

            # Gather all protein and gene names
            protein_names = set()
            gene_names = set()

            for entry in data["results"]:
                desc = entry.get("proteinDescription", {})

                # --- Protein names ---
                if "recommendedName" in desc:
                    full_name = desc["recommendedName"].get("fullName", {}).get("value", "")
                    if full_name:
                        protein_names.add(full_name)

                if "submissionNames" in desc:
                    for sub in desc["submissionNames"]:
                        full_name = sub.get("fullName", {}).get("value", "")
                        if full_name:
                            protein_names.add(full_name)

                if "alternativeNames" in desc:
                    for alt in desc["alternativeNames"]:
                        full_name = alt.get("fullName", {}).get("value", "")
                        if full_name:
                            protein_names.add(full_name)

                # --- Gene names ---
                genes = entry.get("genes", [])
                for g in genes:
                    # main gene name
                    main = g.get("geneName", {}).get("value")
                    if main:
                        gene_names.add(main)
                    # synonyms
                    for syn in g.get("synonyms", []):
                        val = syn.get("value")
                        if val:
                            gene_names.add(val)
                    # ORF names
                    for orf in g.get("orfNames", []):
                        val = orf.get("value")
                        if val:
                            gene_names.add(val)

            # Prepare string representations
            name_str = ", ".join(sorted(protein_names)) if protein_names else None
            gene_str = ", ".join(sorted(gene_names)) if gene_names else None

            # Update dataframe
            df.loc[df["accession"] == accession, "protein_names"] = name_str
            df.loc[df["accession"] == accession, "gene_names"] = gene_str

            new_entries[accession] = {"protein_names": name_str, "gene_names": gene_str, "found": True}
            logging.debug(f"Fetched {accession} -> names: {name_str}; genes: {gene_str}")

            # Periodic save so an interrupted run keeps its progress
            if cache_path and len(new_entries) % 50 == 0:
                _save_name_cache(cache_path, new_entries)

        except Exception as e:
            logging.warning(f"Failed to fetch data for accession {accession}: {e}")

    _save_name_cache(cache_path, new_entries)
    logging.info(f"UniProt names: {n_hits} from cache, {len(new_entries)} newly fetched")

    return df[["accession", "protein_names", "gene_names"]].drop_duplicates().set_index("accession")


# ----------------------
# Main pipeline
# ----------------------
def main(base_dir: str, output_dir: str, input_dirs: List[str], target_bacterium: str,
         verbose: bool, threshold: Optional[float] = None, n_bootstrap: int = 1000,
         ci_level: float = 0.95, neg_fraction: float = 0.75, seed: int = 42,
         subsample_stage: str = "pool", name_cache: Optional[str] = "uniprot_name_cache.json") -> None:
    """
    Main pipeline function.
    Args:
        base_dir (str): Base directory containing bacterium subdirectories.
        output_dir (str): Directory to save outputs.
        input_dirs (List[str]): List of bacterium names to include.
        target_bacterium (str): Bacterium to predict on (excluded from training; unlabelled).
        verbose (bool): Whether to enable verbose logging.
        threshold (Optional[float]): Fixed probability threshold. If None (default), the
            threshold is chosen by maximising MCC on the training split's OOB predictions.
        n_bootstrap (int): Number of bootstrap resamples for confidence intervals.
        ci_level (float): Confidence level for intervals (e.g. 0.95).
        neg_fraction (float): Fraction of negative samples to keep (default 0.75).
        seed (int): Random seed for the negative subsampling only.
        subsample_stage (str): 'pool' (subsample before the train/val split) or
            'train' (subsample only the training split; validation set identical to original).
        name_cache (Optional[str]): JSON file caching UniProt protein/gene names across runs
            (shared by all runs that point to the same path). None/empty disables the cache.
    """
    setup_logging(verbose)
    os.makedirs(output_dir, exist_ok=True)

    # 1) Load all requested bacteria
    per_bacteria = {}
    for bacteria in input_dirs:
        df = load_bacterium_data(base_dir, bacteria)
        if df.empty:
            logging.warning(f"No data for {bacteria}; skipping.")
            continue
        per_bacteria[bacteria] = df
    if not per_bacteria:
        logging.error("No data loaded. Exiting.")
        return

    # 2) Build full matrix per bacterium
    matrices = {}
    labels_per_bacteria = {}
    for bacteria, df in per_bacteria.items():
        try:
            feature_matrix, labels = preprocess_and_pivot(df)
        except Exception as e:
            logging.error(f"Failed to pivot data for {bacteria}: {e}")
            continue
        if feature_matrix.empty:
            logging.warning(f"Pivot produced empty matrix for {bacteria}")
            continue
        matrices[bacteria] = feature_matrix
        labels_per_bacteria[bacteria] = labels

    if target_bacterium not in matrices:
        logging.error(f"Target bacterium '{target_bacterium}' not found among loaded matrices. Available: {list(matrices.keys())}")
        return

    # 3) Create training set = all other bacteria combined
    train_bacts = [bacteria for bacteria in matrices.keys() if bacteria != target_bacterium]
    if not train_bacts:
        logging.error("No bacteria left for training after excluding target bacterium.")
        return

    logging.info(f"Training on: {train_bacts}. Predicting on target: {target_bacterium}")

    # Align feature columns across all bacteria: union of columns, fill missing with 0
    all_features = sorted(set().union(*(matrices[bacteria].columns.tolist() for bacteria in matrices)))
    logging.info(f"Union feature count: {len(all_features)}")

    def align_df(df: pd.DataFrame) -> pd.DataFrame:
        # ensure all_features columns present
        missing = [c for c in all_features if c not in df.columns]
        if missing:
            for m in missing:
                df[m] = 0.0
        # Reorder columns
        return df[all_features]

    X_train_list = []
    y_train_list = []
    for bacteria in train_bacts:
        feature_matrix = align_df(matrices[bacteria])
        labels = labels_per_bacteria.get(bacteria)
        if labels is None:
            logging.warning(f"No labels for bacterium {bacteria}; skipping those samples from training.")
            continue
        # Align labels index/order to feature_matrix index
        labels = labels.reindex(feature_matrix.index).ffill()  # best-effort
        try:
            yb, le_b = make_binary_labels(labels)
        except Exception as e:
            logging.error(f"Failed converting labels for training bacterium {bacteria}: {e}. Skipping.")
            continue
        X_train_list.append(feature_matrix)
        y_train_list.append(pd.Series(yb, index=feature_matrix.index))

    if not X_train_list:
        logging.error("No labeled training data available after processing. Exiting.")
        return

    X_train = pd.concat(X_train_list, axis=0)
    y_train = pd.concat(y_train_list, axis=0).to_numpy()
    logging.info(f"Final training matrix: {X_train.shape}, labels: {np.bincount(y_train)}")
    n_neg_available = int((y_train == 0).sum())
    n_pos_available = int((y_train == 1).sum())

    # ------------------------------------------------------------------
    # SENSITIVITY CHECK: random subsample of the negatives
    # ------------------------------------------------------------------
    if subsample_stage not in ("pool", "train"):
        raise ValueError("subsample_stage must be 'pool' or 'train'")

    retained_df = None
    if subsample_stage == "pool":
        X_train, y_train, retained_df = subsample_negatives(X_train, y_train, neg_fraction, seed)
        logging.info(f"[pool] Kept {int(round(neg_fraction * 100))}% of negatives "
                     f"(seed={seed}): {n_neg_available} -> {int((y_train == 0).sum())} negatives "
                     f"({n_pos_available} positives untouched)")

    # Split a validation set from the combined training data (stratify if possible)
    try:
        X_tr, X_val, y_tr, y_val = train_test_split(X_train, y_train, test_size=0.15, stratify=y_train, random_state=42)
    except Exception:
        X_tr, X_val, y_tr, y_val = train_test_split(X_train, y_train, test_size=0.15, random_state=42)

    if subsample_stage == "train":
        n_neg_tr_before = int((y_tr == 0).sum())
        X_tr, y_tr, retained_df = subsample_negatives(X_tr, y_tr, neg_fraction, seed)
        logging.info(f"[train] Kept {int(round(neg_fraction * 100))}% of training-split negatives "
                     f"(seed={seed}): {n_neg_tr_before} -> {int((y_tr == 0).sum())}; "
                     f"validation set unchanged ({len(y_val)} samples)")

    retained_df.to_csv(os.path.join(output_dir, "negatives_retained.csv"), index=False)

    # Impute if necessary (shouldn't be needed if we filled with zeros)
    imputer = SimpleImputer(strategy="median")
    X_tr_imputed = pd.DataFrame(imputer.fit_transform(X_tr), index=X_tr.index, columns=X_tr.columns)
    X_val_imputed = pd.DataFrame(imputer.transform(X_val), index=X_val.index, columns=X_val.columns)

    # 4) Train Random Forest
    rf = train_rf(X_tr_imputed, y_tr)

    # Choose threshold (MCC-optimised on OOB predictions, i.e. without touching the validation set)
    if threshold is None:
        oob_probs = rf.oob_decision_function_[:, 1]
        ok = ~np.isnan(oob_probs)
        threshold, oob_mcc, curve_df = select_threshold_mcc(y_tr[ok], oob_probs[ok])
        curve_df.to_csv(os.path.join(output_dir, "threshold_mcc_curve.csv"), index=False)
        threshold_source = (f"MCC-optimised on out-of-bag training predictions "
                            f"(OOB MCC = {oob_mcc:.4f}, {int(ok.sum())} samples)")
        logging.info(f"Selected threshold {threshold:.2f} by max OOB MCC ({oob_mcc:.4f})")
    else:
        threshold_source = "fixed value supplied via --threshold"
        logging.info(f"Using fixed threshold {threshold}")

    # Validate
    val_probs = rf.predict_proba(X_val_imputed)[:, 1]
    val_metrics = point_metrics(y_val, val_probs, threshold)
    val_cis = bootstrap_ci(y_val, val_probs, threshold, n_boot=n_bootstrap, ci_level=ci_level)

    logging.info(
        f"Validation @thr={threshold}: "
        f"Sens={_fmt(val_metrics['sensitivity_recall'])} Spec={_fmt(val_metrics['specificity'])} "
        f"Prec={_fmt(val_metrics['precision'])} F1={_fmt(val_metrics['f1'])} "
        f"MCC={_fmt(val_metrics['mcc'])} ROC-AUC={_fmt(val_metrics['roc_auc'])} "
        f"PR-AUC={_fmt(val_metrics['pr_auc'])}"
    )

    # Compute class splits
    def class_split(y: np.ndarray) -> Tuple[int, int]:
        pos = int((y == 1).sum())
        neg = int((y == 0).sum())
        return pos, neg

    train_pos, train_neg = class_split(y_tr)
    val_pos, val_neg = class_split(y_val)

    # Save model report
    report_path = os.path.join(output_dir, "model_report.txt")
    with open(report_path, "w") as fh:
        fh.write("*** SENSITIVITY CHECK: random negative subsampling ***\n")
        fh.write(f"Negatives available in training bacteria: {n_neg_available} "
                 f"(positives: {n_pos_available})\n")
        fh.write(f"Fraction of negatives kept: {neg_fraction:.2f} | sampling seed: {seed} | "
                 f"stage: {subsample_stage}\n")
        if subsample_stage == "train":
            fh.write("Validation set is identical to the original (un-subsampled) run.\n")
        else:
            fh.write("Negatives were subsampled before the train/validation split, so the "
                     "validation set has a different prevalence than in the original run "
                     "(precision / PR-AUC are not directly comparable).\n")
        fh.write("\n")
        fh.write(f"Training bacteria: {train_bacts}\n")
        fh.write(f"Target bacterium (prediction only, unlabelled): {target_bacterium}\n")
        fh.write(f"Training samples: {X_tr.shape[0]} (pos={train_pos}, neg={train_neg})\n")
        fh.write(f"Validation samples: {X_val.shape[0]} (pos={val_pos}, neg={val_neg})\n")
        fh.write(f"Threshold: {threshold:.2f} -- {threshold_source}\n")
        fh.write("Note: the threshold was chosen without using the validation set; the CIs below "
                 "do not include uncertainty in the threshold itself.\n\n")
        fh.write(format_metrics_report(val_metrics, val_cis, threshold, ci_level,
                                       n_bootstrap, val_pos, val_neg))

    # Save metrics as CSV too
    rows = []
    for k in METRIC_ORDER:
        rows.append({"metric": METRIC_LABELS[k], "estimate": val_metrics[k],
                     "ci_lower": val_cis[k][0], "ci_upper": val_cis[k][1],
                     "threshold": threshold if k not in ("roc_auc", "pr_auc") else np.nan})
    for k in ("tn", "fp", "fn", "tp"):
        rows.append({"metric": k.upper(), "estimate": val_metrics[k],
                     "ci_lower": np.nan, "ci_upper": np.nan, "threshold": threshold})
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "validation_metrics.csv"), index=False)
    logging.info("Saved validation metrics.")

    # 5) Feature importances
    importances = rf.feature_importances_
    feat_imp_df = pd.DataFrame({"feature": X_tr_imputed.columns, "importance": importances})
    feat_imp_df = feat_imp_df.sort_values("importance", ascending=False).reset_index(drop=True)
    feat_imp_df["rank"] = feat_imp_df.index + 1
    feat_imp_df.to_csv(os.path.join(output_dir, "feature_importances.csv"), index=False)
    logging.info(f"Saved feature importances to {os.path.join(output_dir, 'feature_importances.csv')}")

    # 6) Predict on the target bacterium (no labels / no AUC: prediction only)
    X_target_raw = align_df(matrices[target_bacterium])
    X_target_imputed = pd.DataFrame(imputer.transform(X_target_raw),
                                    index=X_target_raw.index, columns=X_target_raw.columns)

    target_probs = rf.predict_proba(X_target_imputed)[:, 1]
    target_preds = (target_probs >= threshold).astype(int)
    out_df = pd.DataFrame({
        "accession": X_target_imputed.index,
        "prob_antigen": target_probs,
        "pred_label": target_preds
    }).set_index("accession")

    # --- Add protein names using the existing function ---
    try:
        logging.info("Fetching protein names for target accessions...")
        name_map = load_accession_name_and_gene_fetch(out_df.reset_index()[["accession"]],
                                                     cache_path=name_cache)
        out_df = out_df.merge(name_map, left_index=True, right_index=True, how="left")
        logging.info("Added protein names to predictions.")
    except Exception as e:
        logging.warning(f"Failed to add protein names: {e}")

    # Order from most likely antigenic -> least likely
    out_df = out_df.sort_values("prob_antigen", ascending=False)
    out_csv = os.path.join(output_dir, f"{target_bacterium.replace('.', '_')}_predictions.csv")
    out_df.to_csv(out_csv)
    logging.info(f"Saved target predictions to {out_csv} (rows: {len(out_df)})")

    # Save full target feature matrix with probabilities
    X_target_export = X_target_imputed.copy()
    X_target_export["prob_antigen"] = target_probs
    X_target_export = X_target_export.sort_values("prob_antigen", ascending=False)
    X_target_export.to_csv(os.path.join(output_dir, f"{target_bacterium.replace('.', '_')}_features_with_probs.csv"))
    logging.info("Saved target features with probabilities.")

    # Target section of the report: prediction summary only (target is unlabelled)
    n_pred_pos = int(target_preds.sum())
    with open(report_path, "a") as fh:
        fh.write("=" * 70 + "\n")
        fh.write(f"TARGET BACTERIUM PREDICTIONS: {target_bacterium}\n")
        fh.write("=" * 70 + "\n")
        fh.write("Target data are unlabelled, so no accuracy/AUC is reported for this bacterium.\n")
        fh.write(f"Target proteins scored: {len(out_df)}\n")
        fh.write(f"Predicted antigens at threshold {threshold}: {n_pred_pos} "
                 f"({n_pred_pos / len(out_df):.1%})\n")
        fh.write(f"Mean predicted probability: {float(np.mean(target_probs)):.4f}\n")

    logging.info("Pipeline finished successfully.")


# ----------------------
# CLI
# ----------------------
if __name__ == "__main__":
    """ Main entry point """
    parser = argparse.ArgumentParser(
        description="Random Forest antigenicity pipeline -- negative-subsampling sensitivity check")
    parser.add_argument("--base-dir", type=str, default="./results", help="Base directory containing <bacterium>/raw_data/*_raw_data.csv")
    parser.add_argument("--output-dir", type=str, default="./results_neg75", help="Directory to save outputs (default differs from the original so results are not overwritten)")
    parser.add_argument("--input-dir", nargs="+", required=True, help="List of bacteria to load (folders under base_dir). e.g. E.coli S.aureus")
    parser.add_argument("--target-bacterium", type=str, default="S.aureus", help="Bacterium to hold out of training and predict on (default: S.aureus)")
    parser.add_argument("--threshold", type=float, default=None, help="Fixed probability threshold. Default: choose the MCC-optimal threshold from out-of-bag training predictions")
    parser.add_argument("--n-bootstrap", type=int, default=1000, help="Bootstrap resamples for confidence intervals (default: 1000)")
    parser.add_argument("--ci-level", type=float, default=0.95, help="Confidence level for intervals (default: 0.95)")
    parser.add_argument("--neg-fraction", type=float, default=0.75, help="Fraction of negative samples to keep, sampled at random without replacement (default: 0.75)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for the negative subsampling (change it to repeat the check with a different random subset)")
    parser.add_argument("--subsample-stage", choices=["pool", "train"], default="pool",
                        help="'pool': subsample negatives before the train/val split (default). "
                             "'train': keep the original validation set and subsample only the training split")
    parser.add_argument("--name-cache", type=str, default="uniprot_name_cache.json",
                        help="JSON file caching UniProt protein/gene names across runs (default: ./uniprot_name_cache.json). "
                             "Use the same path for every run in a loop. Pass an empty string ('') to disable; delete the file to reset")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    main(base_dir=args.base_dir, output_dir=args.output_dir, input_dirs=args.input_dir,
         target_bacterium=args.target_bacterium, verbose=args.verbose,
         threshold=args.threshold, n_bootstrap=args.n_bootstrap, ci_level=args.ci_level,
         neg_fraction=args.neg_fraction, seed=args.seed, subsample_stage=args.subsample_stage,
         name_cache=args.name_cache or None)