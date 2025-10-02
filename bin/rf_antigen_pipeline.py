#!/usr/bin/env python3
"""
Random Forest antigen-prediction pipeline.

Reads raw CSVs from <base_dir>/<bacterium>/raw_data/*_raw_data.csv
Each raw CSV is expected to contain at least these columns:
 - accession
 - feature
 - subfeature
 - value
 - label

Pipeline:
 1) Load and concatenate data per bacterium
 2) Pivot to accession x feature_subfeature matrix (mean of replicates)
 3) Train RandomForest on all bacteria EXCEPT the test bacterium (default: S.aureus)
 4) Compute feature importances (feature usefulness)
 5) Predict probabilities on test bacterium and save CSV sorted by probability (desc)

Outputs (saved in output_dir):
 - saureus_predictions.csv  (accession, prob_antigen, pred_label, true_label if available)
 - feature_importances.csv (feature, importance, rank)
 - model_report.txt        (train/test sizes, AUC if labels exist on test)
"""

import os
import glob
import argparse
import logging
from typing import List, Tuple
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, classification_report
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split

# ----------------------
# Logging setup
# ----------------------
def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

# ----------------------
# Data loading utils
# ----------------------
def load_bacterium_data(base_dir: str, bacterium: str) -> pd.DataFrame:
    """Load and concatenate all raw data CSVs for a bacterium."""
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
    # Standardize column names lower-case
    out.columns = [c.strip() for c in out.columns]
    return out

# ----------------------
# Preprocess & pivot
# ----------------------
def preprocess_and_pivot(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Build feature matrix:
      - feature_subfeature = feature_subfeature
      - pivot so each accession is a row, each column is a feature_subfeature (mean)
    Returns (X_df, labels_series) where:
      - X_df.index = accession
      - labels_series.index = accession (if label present), dtype object
    """
    required_cols = {"accession", "feature", "subfeature", "value"}
    if not required_cols.issubset(set(df.columns)):
        raise ValueError(f"Input dataframe missing required cols: {required_cols - set(df.columns)}")
    df = df.copy()
    # Compose feature id
    df["feature_subfeature"] = df["feature"].astype(str) + "_" + df["subfeature"].astype(str)
    # Keep only positive numeric values if possible
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["value"].notna()]
    # Aggregate mean across replicates
    agg = (
        df.groupby(["accession", "feature_subfeature"], observed=True)["value"]
          .mean()
          .reset_index()
    )
    # Pivot
    X = agg.pivot(index="accession", columns="feature_subfeature", values="value")
    # Fill missing with 0 (absent signal). You can change to median/imputer if desired.
    X = X.fillna(0.0)
    # Labels (if present)
    labels = None
    if "label" in df.columns:
        # choose majority label per accession
        labels = (df.groupby("accession", observed=True)["label"]
                    .agg(lambda x: x.mode().iat[0] if not x.mode().empty else x.iloc[0]))
        labels = labels.loc[X.index]  # align to X rows
    return X, labels

# ----------------------
# Label handling heuristic -> binary
# ----------------------
def make_binary_labels(labels: pd.Series) -> Tuple[np.ndarray, LabelEncoder]:
    """
    Convert labels to binary 0/1.
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
        if unique_numeric.issubset({0,1}):
            y = labels.astype(int).to_numpy()
            le = LabelEncoder()
            le.fit([0,1])
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
        le.fit(["0","1"])
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
# Train Random Forest
# ----------------------
def train_rf(X_train: pd.DataFrame, y_train: np.ndarray, n_estimators:int=500, random_state:int=42) -> RandomForestClassifier:
    logging.info(f"Training RandomForest on {X_train.shape[0]} samples and {X_train.shape[1]} features...")
    rf = RandomForestClassifier(n_estimators=n_estimators, n_jobs=-1, random_state=random_state, class_weight='balanced')
    rf.fit(X_train, y_train)
    return rf

# ----------------------
# Main pipeline
# ----------------------
def main(base_dir: str, output_dir: str, input_dirs: List[str], test_bacterium: str, verbose: bool) -> None:
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
            Xb, lb = preprocess_and_pivot(df)
        except Exception as e:
            logging.error(f"Failed to pivot data for {bacteria}: {e}")
            continue
        if Xb.empty:
            logging.warning(f"Pivot produced empty matrix for {bacteria}")
            continue
        matrices[bacteria] = Xb
        labels_per_bacteria[bacteria] = lb

    if test_bacterium not in matrices:
        logging.error(f"Test bacterium '{test_bacterium}' not found among loaded matrices. Available: {list(matrices.keys())}")
        return

    # 3) Create training set = all other bacteria combined
    train_bacts = [bacteria for bacteria in matrices.keys() if bacteria != test_bacterium]
    if not train_bacts:
        logging.error("No bacteria left for training after excluding test bacterium.")
        return

    logging.info(f"Training on: {train_bacts}. Testing on: {test_bacterium}")

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
        Xb = align_df(matrices[bacteria])
        lb = labels_per_bacteria.get(bacteria)
        if lb is None:
            logging.warning(f"No labels for bacterium {bacteria}; skipping those samples from training.")
            continue
        # Align labels index/order to Xb index
        lb = lb.reindex(Xb.index).fillna(method="ffill")  # best-effort
        try:
            yb, le_b = make_binary_labels(lb)
        except Exception as e:
            logging.error(f"Failed converting labels for training bacterium {bacteria}: {e}. Skipping.")
            continue
        X_train_list.append(Xb)
        y_train_list.append(pd.Series(yb, index=Xb.index))

    if not X_train_list:
        logging.error("No labeled training data available after processing. Exiting.")
        return

    X_train = pd.concat(X_train_list, axis=0)
    y_train = pd.concat(y_train_list, axis=0).to_numpy()
    logging.info(f"Final training matrix: {X_train.shape}, labels: {np.bincount(y_train)}")

    # Optional: split a small validation set from the combined training for quick check (stratify if possible)
    try:
        X_tr, X_val, y_tr, y_val = train_test_split(X_train, y_train, test_size=0.15, stratify=y_train, random_state=42)
    except Exception:
        X_tr, X_val, y_tr, y_val = train_test_split(X_train, y_train, test_size=0.15, random_state=42)

    # Impute if necessary (shouldn't be needed if we filled with zeros)
    imputer = SimpleImputer(strategy="median")
    X_tr_imputed = pd.DataFrame(imputer.fit_transform(X_tr), index=X_tr.index, columns=X_tr.columns)
    X_val_imputed = pd.DataFrame(imputer.transform(X_val), index=X_val.index, columns=X_val.columns)

    # 4) Train Random Forest
    rf = train_rf(X_tr_imputed, y_tr)

    # Validate
    val_probs = rf.predict_proba(X_val_imputed)[:, 1]
    val_preds = (val_probs >= 0.5).astype(int)
    val_auc = None
    try:
        val_auc = roc_auc_score(y_val, val_probs)
    except Exception:
        val_auc = None
    val_acc = accuracy_score(y_val, val_preds)
    logging.info(f"Validation ACC: {val_acc:.4f} AUC: {val_auc if val_auc is None else val_auc:.4f}")

    # Save model report
    report_path = os.path.join(output_dir, "model_report.txt")
    with open(report_path, "w") as fh:
        fh.write(f"Training bacteria: {train_bacts}\n")
        fh.write(f"Training samples: {X_tr.shape[0]}\n")
        fh.write(f"Validation samples: {X_val.shape[0]}\n")
        fh.write(f"Validation accuracy: {val_acc:.4f}\n")
        fh.write(f"Validation AUC: {val_auc}\n")

    # 5) Feature importances
    importances = rf.feature_importances_
    feat_imp_df = pd.DataFrame({"feature": X_tr_imputed.columns, "importance": importances})
    feat_imp_df = feat_imp_df.sort_values("importance", ascending=False).reset_index(drop=True)
    feat_imp_df["rank"] = feat_imp_df.index + 1
    feat_imp_df.to_csv(os.path.join(output_dir, "feature_importances.csv"), index=False)
    logging.info(f"Saved feature importances to {os.path.join(output_dir, 'feature_importances.csv')}")

    # 6) Test on S.aureus (test_bacterium)
    X_test_raw = align_df(matrices[test_bacterium])
    # if labels exist for test, attempt to make binary labels using same heuristic
    y_test_series = labels_per_bacteria.get(test_bacterium)
    y_test = None
    label_encoder_for_test = None
    if y_test_series is not None:
        try:
            y_test, label_encoder_for_test = make_binary_labels(y_test_series.reindex(X_test_raw.index).fillna(method="ffill"))
        except Exception as e:
            logging.warning(f"Couldn't create binary labels for test bacterium: {e}. Continuing without true labels.")
            y_test = None

    X_test_imputed = pd.DataFrame(imputer.transform(X_test_raw), index=X_test_raw.index, columns=X_test_raw.columns)

    # Predict probabilities
    test_probs = rf.predict_proba(X_test_imputed)[:, 1]
    test_preds = (test_probs >= 0.5).astype(int)
    out_df = pd.DataFrame({
        "accession": X_test_imputed.index,
        "prob_antigen": test_probs,
        "pred_label": test_preds
    }).set_index("accession")

    if y_test is not None:
        out_df["true_label"] = y_test

    # Order from most likely antigenic -> least likely
    out_df = out_df.sort_values("prob_antigen", ascending=False)
    out_csv = os.path.join(output_dir, f"{test_bacterium.replace('.','_')}_predictions.csv")
    out_df.to_csv(out_csv)
    logging.info(f"Saved test predictions to {out_csv} (rows: {len(out_df)})")

    # Also save the feature values for the top N proteins (optional small sample)
    # Save full X_test with probabilities
    X_test_export = X_test_imputed.copy()
    X_test_export["prob_antigen"] = test_probs
    X_test_export = X_test_export.sort_values("prob_antigen", ascending=False)
    X_test_export.to_csv(os.path.join(output_dir, f"{test_bacterium.replace('.','_')}_features_with_probs.csv"))
    logging.info("Saved test features with probabilities.")

    # If test labels exist compute metrics
    if y_test is not None:
        try:
            test_auc = roc_auc_score(y_test, test_probs)
        except Exception:
            test_auc = None
        test_acc = accuracy_score(y_test, test_preds)
        logging.info(f"Test ACC: {test_acc:.4f} AUC: {test_auc}")
        with open(report_path, "a") as fh:
            fh.write(f"Test bacterium: {test_bacterium}\n")
            fh.write(f"Test samples: {X_test_imputed.shape[0]}\n")
            fh.write(f"Test accuracy: {test_acc:.4f}\n")
            fh.write(f"Test AUC: {test_auc}\n")

    logging.info("Pipeline finished successfully.")


# ----------------------
# CLI
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Random Forest antigenicity pipeline")
    parser.add_argument("--base-dir", type=str, default="./results", help="Base directory containing <bacterium>/raw_data/*_raw_data.csv")
    parser.add_argument("--output-dir", type=str, default="./results", help="Directory to save outputs")
    parser.add_argument("--input-dir", nargs="+", required=True, help="List of bacteria to load (folders under base_dir). e.g. E.coli S.aureus")
    parser.add_argument("--test-bacterium", type=str, default="S.aureus", help="Bacterium to hold-out and test (default: S.aureus)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    main(base_dir=args.base_dir, output_dir=args.output_dir, input_dirs=args.input_dir, test_bacterium=args.test_bacterium, verbose=args.verbose)
