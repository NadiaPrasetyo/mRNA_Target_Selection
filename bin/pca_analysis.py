import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
import argparse
import logging

# ----------------------
# Configuration
# ----------------------
def setup_logging(verbose: bool) -> None:
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
# Load data
# ----------------------
def load_bacterium_data(base_dir: str, bacterium: str) -> pd.DataFrame:
    folder = os.path.join(base_dir, bacterium, "raw_data")
    logging.info(f"Loading data for bacterium: {bacterium} from {folder}")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    logging.debug(f"Found {len(files)} files for {bacterium}")
    dfs = []
    for f in files:
        logging.debug(f"Reading file: {f}")
        df = pd.read_csv(f)
        df["bacterium"] = bacterium
        dfs.append(df)
    if not dfs:
        logging.warning(f"No files found for {bacterium} in {folder}")
        return pd.DataFrame()
    out = pd.concat(dfs, ignore_index=True)

    # Optimize dtypes
    for col in ["label", "feature", "subfeature", "bacterium"]:
        out[col] = out[col].astype("category")
    out["value"] = pd.to_numeric(out["value"], downcast="float")
    return out


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
    logging.info(f"Raw dataframe memory usage: {all_data.memory_usage(deep=True).sum()/1e9:.2f} GB")

    # ----------------------
    # Pivot → Sparse matrix
    # ----------------------
    logging.info("Pivoting data to sparse wide format...")
    all_data["sample_id"] = (
        all_data.groupby(["label"]).cumcount().astype(str)
        + "_" + all_data["label"].astype(str)
    )
    all_data["feature_subfeature"] = (
        all_data["feature"].astype(str) + "_" + all_data["subfeature"].astype(str)
    )

    wide_df = all_data.pivot_table(
        index="sample_id",
        columns="feature_subfeature",
        values="value",
        aggfunc="mean"
    ).astype(pd.SparseDtype("float32", fill_value=0))

    logging.info(f"Pivoted to shape {wide_df.shape}, sparse density ≈ {wide_df.sparse.density:.6f}")
    logging.info(f"Approx sparse storage size: {wide_df.memory_usage(deep=True).sum()/1e9:.2f} GB")

    meta = all_data.drop_duplicates(subset=["sample_id"])[["sample_id", "label"]].set_index("sample_id")

    # ----------------------
    # Incremental PCA
    # ----------------------
    logging.info("Performing Incremental PCA...")
    X_sparse = wide_df.sparse.to_coo()   # convert to scipy.sparse
    scaler = StandardScaler(with_mean=False)  # works with sparse
    X_scaled = scaler.fit_transform(X_sparse)

    ipca = IncrementalPCA(n_components=50, batch_size=10000)
    pcs = ipca.fit_transform(X_scaled)
    logging.info("Incremental PCA completed")

    pca_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"], index=wide_df.index).join(meta)

    # ----------------------
    # Scree plot
    # ----------------------
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(ipca.explained_variance_ratio_) + 1),
             ipca.explained_variance_ratio_, 'o-', markersize=6)
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    plt.title('Scree Plot')
    plt.grid(True)
    scree_plot_path = os.path.join(output_dir, "scree_plot.png")
    plt.savefig(scree_plot_path, dpi=300)
    plt.close()
    logging.info(f"Scree plot saved to {scree_plot_path}")

    # ----------------------
    # PCA biplot (skip arrows if too many features)
    # ----------------------
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="label", s=60, alpha=0.8)
    n_features = wide_df.shape[1]
    if n_features <= 200:  # safety cutoff
        loadings = ipca.components_[:2].T * np.sqrt(ipca.explained_variance_[:2])
        for i, col in enumerate(wide_df.columns):
            plt.arrow(0, 0, loadings[i, 0], loadings[i, 1], color='gray', alpha=0.5, head_width=0.02)
            plt.text(loadings[i, 0]*1.1, loadings[i, 1]*1.1, col, fontsize=8, alpha=0.7)
    plt.title("PCA Biplot")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    biplot_path = os.path.join(output_dir, "pca_biplot.png")
    plt.savefig(biplot_path, dpi=300)
    plt.close()
    logging.info(f"PCA biplot saved to {biplot_path}")

    # ----------------------
    # Skip huge covariance
    # ----------------------
    if wide_df.shape[1] <= 2000:
        logging.info("Generating covariance matrix heatmap...")
        cov_matrix = np.cov(X_scaled.T.toarray())
        cov_df = pd.DataFrame(cov_matrix, index=wide_df.columns, columns=wide_df.columns)
        plt.figure(figsize=(12, 10))
        sns.heatmap(cov_df, cmap="coolwarm", center=0)
        covariance_path = os.path.join(output_dir, "covariance_matrix.png")
        plt.savefig(covariance_path, dpi=300)
        plt.close()
        logging.info(f"Covariance matrix heatmap saved to {covariance_path}")
    else:
        logging.warning("Covariance matrix skipped (too many features)")

    logging.info(f"✅ Analysis complete. Plots saved in {output_dir}")


# ----------------------
# Entry point
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run PCA analysis on bacterial raw data.")
    parser.add_argument("--base-dir", type=str, default="./results", help="Base directory containing bacteria data , default: ./results")
    parser.add_argument("--output-dir", type=str, default="./results", help="Directory to save output plots, default: ./results")
    parser.add_argument("--input-dir", nargs="+", required=True, help="List of bacteria names (folders under base-dir)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    setup_logging(args.verbose)
    main(args.base_dir, args.output_dir, args.input_dir)
