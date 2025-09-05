import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from scipy.sparse import coo_matrix
import argparse
import logging
from adjustText import adjust_text


# ----------------------
# Configuration
# ----------------------
def setup_logging(verbose: bool) -> None:
    """Configure logging to console and optionally file."""
    log_file = "pca_analysis.log" if verbose else None
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"

    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(level)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(format_str))
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(format_str))
        logger.addHandler(file_handler)


# ----------------------
# Data loading
# ----------------------
def load_bacterium_data(base_dir: str, bacterium: str) -> pd.DataFrame:
    """Load and concatenate all raw data CSVs for a bacterium."""
    folder = os.path.join(base_dir, bacterium, "raw_data")
    logging.info(f"Loading data for bacterium: {bacterium} from {folder}")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    dfs = []
    for f in files:
        logging.debug(f"Reading file: {f}")
        df = pd.read_csv(f)
        df["bacterium"] = bacterium
        dfs.append(df)
    if not dfs:
        logging.warning(f"No files found for {bacterium}")
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)

    # Optimize dtypes
    for col in ["label", "feature", "subfeature", "bacterium"]:
        out[col] = out[col].astype("category")
    out["value"] = pd.to_numeric(out["value"], downcast="float")
    return out


def preprocess_data(df: pd.DataFrame) -> pd.DataFrame:
    """Drop missing values and non-positive values."""
    return df.dropna()[df["value"] > 0]


# ----------------------
# Plotting utilities
# ----------------------
def plot_scree(ipca, output_dir: str):
    """Plot scree plot with variance explained and cumulative variance."""
    explained_var = ipca.explained_variance_ratio_ * 100
    cum_var = np.cumsum(explained_var)

    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(explained_var) + 1), explained_var, 'o-', label="Individual")
    plt.plot(range(1, len(cum_var) + 1), cum_var, 's--', label="Cumulative")
    plt.axhline(1, color="gray", linestyle="--", alpha=0.5)  # Kaiser criterion
    plt.axhline(90, color="red", linestyle="--", alpha=0.6, label="90% threshold")

    plt.xlabel("Principal Component")
    plt.ylabel("Variance Explained (%)")
    plt.title("Scree Plot")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(os.path.join(output_dir, "scree_plot.png"), dpi=300)
    plt.close()

    # Save explained variance table
    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(explained_var))],
        "Variance (%)": explained_var,
        "Cumulative (%)": cum_var
    }).to_csv(os.path.join(output_dir, "explained_variance.csv"), index=False)


def plot_pca_biplot(pca_df, ipca, feature_enc, output_dir: str, top_n=10, scale=2.5):
    """PCA biplot with samples and top feature loadings, improved label placement."""
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="label", s=60, alpha=0.7)

    # Loadings
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    # Top features by vector length
    norms = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(norms)[-top_n:]

    texts = []
    for i in top_idx:
        x, y = loadings[i, 0] * scale, loadings[i, 1] * scale
        plt.arrow(0, 0, x, y, color='red', alpha=0.5, head_width=0.02*scale)
        texts.append(plt.text(x, y, feature_names[i], fontsize=9, color="darkred"))

    # Adjust labels with repelling forces
    adjust_text(
        texts,
        only_move={'points':'y', 'text':'y'},  # move vertically for cleaner alignment
        arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
        expand_points=(1.2, 1.4),
        expand_text=(1.2, 1.4),
        force_points=0.2,
        force_text=0.2
    )

    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.title("PCA Biplot (Top Features)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_biplot.png"), dpi=300)
    plt.close()
    

def plot_loading_scatter(ipca, feature_enc, output_dir: str, top_n=10):
    """Scatter plot of feature loadings on PC1 vs PC2."""
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    norms = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(norms)[-top_n:]

    plt.figure(figsize=(10, 8))
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.scatter(loadings[:, 0], loadings[:, 1], alpha=0.3, color="lightgray", s=20, label="Other features")

    for i in top_idx:
        x, y = loadings[i, 0], loadings[i, 1]
        plt.scatter(x, y, color="red", s=60)
        plt.text(x * 1.1, y * 1.1, feature_names[i], fontsize=9)

    plt.xlabel("Loading on PC1")
    plt.ylabel("Loading on PC2")
    plt.title(f"PCA Feature Loadings (Top {top_n})")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(os.path.join(output_dir, "pca_loading_scatter.png"), dpi=300)
    plt.close()


def plot_correlation_matrix(X_scaled, feature_enc, output_dir: str, max_features=2000):
    """Plot correlation matrix heatmap with clustering."""
    n_features = X_scaled.shape[1]

    if n_features > max_features:
        logging.warning(f"Too many features ({n_features}), downsampling to top {max_features}.")
        col_var = np.array(X_scaled.power(2).mean(axis=0) - np.power(X_scaled.mean(axis=0), 2)).ravel()
        top_idx = np.argpartition(col_var, -max_features)[-max_features:]
        X_scaled = X_scaled[:, top_idx]
        feature_names = feature_enc.inverse_transform(top_idx)
    else:
        feature_names = feature_enc.inverse_transform(np.arange(n_features))

    X_dense = X_scaled.toarray().astype(np.float32)
    corr_matrix = np.corrcoef(X_dense, rowvar=False)

    corr_df = pd.DataFrame(corr_matrix, index=feature_names, columns=feature_names)
    sns.clustermap(corr_df, cmap="coolwarm", center=0, figsize=(14, 12))
    plt.savefig(os.path.join(output_dir, "correlation_matrix.png"), dpi=300)
    plt.close()

    # Save values
    corr_df.to_csv(os.path.join(output_dir, "correlation_matrix.csv"))


def plot_covariance_matrix(X_scaled, feature_enc, output_dir: str, max_features=2000):
    """Plot covariance matrix heatmap (optionally subset to top variable features)."""
    n_features = X_scaled.shape[1]

    if n_features > max_features:
        logging.warning(f"Too many features ({n_features}), downsampling to top {max_features}.")
        col_var = np.array(X_scaled.power(2).mean(axis=0) - np.power(X_scaled.mean(axis=0), 2)).ravel()
        top_idx = np.argpartition(col_var, -max_features)[-max_features:]
        X_scaled = X_scaled[:, top_idx]
        feature_names = feature_enc.inverse_transform(top_idx)
    else:
        feature_names = feature_enc.inverse_transform(np.arange(n_features))

    X_dense = X_scaled.toarray().astype(np.float32)
    cov_matrix = np.cov(X_dense, rowvar=False)

    cov_df = pd.DataFrame(cov_matrix, index=feature_names, columns=feature_names)
    plt.figure(figsize=(12, 10))
    sns.heatmap(cov_df, cmap="coolwarm", center=0)
    plt.title("Feature Covariance Matrix")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "covariance_matrix.png"), dpi=300)
    plt.close()

    # Save values
    cov_df.to_csv(os.path.join(output_dir, "covariance_matrix.csv"))



# ----------------------
# Main analysis
# ----------------------
def main(base_dir: str, output_dir: str, input_dirs: list[str]) -> None:
    logging.info("Starting PCA analysis...")
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    dfs = [load_bacterium_data(base_dir, b) for b in input_dirs]
    dfs = [df for df in dfs if not df.empty]
    if not dfs:
        logging.error("No data loaded. Exiting.")
        return
    all_data = pd.concat(dfs, ignore_index=True)
    logging.info(f"Loaded data for {len(input_dirs)} bacteria. Rows: {len(all_data):,}")

    # ----------------------
    # Assign sample IDs
    # ----------------------
    logging.info("Assigning sample IDs...")
    label_codes, label_indices = pd.factorize(all_data["label"])
    counters = np.zeros(len(np.unique(label_codes)), dtype=int)
    sample_ids = []

    for code in label_codes:
        sample_ids.append(f"{counters[code]}_{all_data['label'].iloc[len(sample_ids)]}")
        counters[code] += 1

    all_data["sample_id"] = sample_ids
    all_data["feature_subfeature"] = all_data["feature"].astype(str) + "_" + all_data["subfeature"].astype(str)


    all_data = preprocess_data(all_data)

    # ----------------------
    # Sparse matrix
    # ----------------------
    sample_enc = LabelEncoder()
    feature_enc = LabelEncoder()
    row_idx = sample_enc.fit_transform(all_data["sample_id"])
    col_idx = feature_enc.fit_transform(all_data["feature_subfeature"])

    X_sparse = coo_matrix(
        (all_data["value"].astype(np.float32), (row_idx, col_idx)),
        shape=(len(sample_enc.classes_), len(feature_enc.classes_))
    ).tocsr()

    logging.info(f"Sparse matrix shape: {X_sparse.shape}, nnz={X_sparse.nnz:,}")


    meta = pd.DataFrame({"sample_id": sample_enc.classes_, "label": [s.split("_")[-1] for s in sample_enc.classes_]})
    meta = meta.set_index("sample_id")

    # PCA
    scaler = StandardScaler(with_mean=False)
    X_scaled = scaler.fit_transform(X_sparse)

    ipca = IncrementalPCA(n_components=50, batch_size=10000)
    pcs = ipca.fit_transform(X_scaled)
    pca_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"], index=sample_enc.classes_).join(meta)

    # Plots
    plot_scree(ipca, output_dir)
    plot_pca_biplot(pca_df, ipca, feature_enc, output_dir)
    plot_loading_scatter(ipca, feature_enc, output_dir)
    plot_covariance_matrix(X_scaled, feature_enc, output_dir)
    plot_correlation_matrix(X_scaled, feature_enc, output_dir)

    logging.info(f"✅ Analysis complete. Plots and CSVs saved in {output_dir}")


# ----------------------
# Entry point
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run PCA analysis on bacterial raw data.")
    parser.add_argument("--base-dir", type=str, default="./results")
    parser.add_argument("--output-dir", type=str, default="./results")
    parser.add_argument("--input-dir", nargs="+", required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    main(args.base_dir, args.output_dir, args.input_dir)
