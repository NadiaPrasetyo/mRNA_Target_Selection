import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse
import logging

# ----------------------
# Configuration
# ----------------------
def setup_logging(verbose: bool) -> None:
    """
    Set up logging configuration.
    If verbose is True, logs will be printed to console and written to a file.
    """
    log_file = "pca_analysis.log" if verbose else None
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"
    
    # Clear previous handlers
    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(level)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(format_str))
    logger.addHandler(console_handler)
    
    # File handler if verbose
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
    logging.debug(f"Found {len(files)} files for bacterium: {bacterium}")
    dfs = []
    for f in files:
        logging.debug(f"Reading file: {f}")
        df = pd.read_csv(f)
        df["bacterium"] = bacterium
        dfs.append(df)
    if not dfs:
        logging.warning(f"No files found for {bacterium} in {folder}")
        return pd.DataFrame()
    logging.info(f"Loaded {len(dfs)} dataframes for bacterium: {bacterium}")
    return pd.concat(dfs, ignore_index=True)


# ----------------------
# Main analysis
# ----------------------
def main(base_dir: str, output_dir: str, input_dirs: list[str]) -> None:
    logging.info("Starting PCA analysis...")
    os.makedirs(output_dir, exist_ok=True)

    # Load all bacteria data into one combined dataset
    logging.info("Loading data for all bacteria...")
    dfs = [load_bacterium_data(base_dir, b) for b in input_dirs]
    dfs = [df for df in dfs if not df.empty]
    if not dfs:
        logging.error("No data loaded. Exiting.")
        return
    all_data = pd.concat(dfs, ignore_index=True)
    logging.info(f"Loaded data for {len(input_dirs)} bacteria. Total rows: {len(all_data)}")

    # ----------------------
    # Pivot to wide format
    # ----------------------
    logging.info("Pivoting data to wide format...")
    all_data["sample_id"] = (
        all_data.groupby(["label"]).cumcount().astype(str)
        + "_" + all_data["label"]
    )
    all_data["feature_subfeature"] = all_data["feature"] + "_" + all_data["subfeature"]
    wide_df = all_data.pivot_table(
        index="sample_id",
        columns="feature_subfeature",
        values="value",
        aggfunc="mean"
    )
    logging.info(f"Pivoted data to wide format. Shape: {wide_df.shape}")

    # Metadata for labels
    meta = all_data.drop_duplicates(subset=["sample_id"])[["sample_id", "label"]].set_index("sample_id")
    logging.debug(f"Metadata shape: {meta.shape}")

    # ----------------------
    # PCA
    # ----------------------
    logging.info("Performing PCA...")
    X = wide_df.fillna(0)
    logging.debug("Missing values filled with 0")
    X_scaled = StandardScaler().fit_transform(X)
    logging.debug("Data standardized")

    pca = PCA()
    pcs = pca.fit_transform(X_scaled)
    logging.info("PCA completed")

    pca_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"], index=wide_df.index)
    pca_df = pca_df.join(meta)
    logging.debug(f"PCA DataFrame shape: {pca_df.shape}")

    # ----------------------
    # Scree plot
    # ----------------------
    logging.info("Generating scree plot...")
    plt.figure(figsize=(8, 6))
    plt.plot(
        range(1, len(pca.explained_variance_ratio_) + 1),
        pca.explained_variance_ratio_, 'o-', markersize=6
    )
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    plt.title('Scree Plot')
    plt.grid(True)
    scree_plot_path = os.path.join(output_dir, "scree_plot.png")
    plt.savefig(scree_plot_path, dpi=300)
    plt.close()
    logging.info(f"Scree plot saved to {scree_plot_path}")

    # ----------------------
    # PCA biplot
    # ----------------------
    logging.info("Generating PCA biplot...")
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=pca_df,
        x="PC1", y="PC2",
        hue="label",
        s=60, alpha=0.8
    )
    loadings = pca.components_[:2].T * np.sqrt(pca.explained_variance_[:2])
    for i, col in enumerate(wide_df.columns):
        plt.arrow(0, 0, loadings[i, 0], loadings[i, 1],
                  color='gray', alpha=0.5, head_width=0.02)
        plt.text(loadings[i, 0]*1.1, loadings[i, 1]*1.1,
                 col, fontsize=8, alpha=0.7)
    plt.title("PCA Biplot (All Bacteria Combined)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    biplot_path = os.path.join(output_dir, "pca_biplot.png")
    plt.savefig(biplot_path, dpi=300)
    plt.close()
    logging.info(f"PCA biplot saved to {biplot_path}")

    # ----------------------
    # Loading plot (Top 2 PCs)
    # ----------------------
    logging.info("Generating PCA loadings heatmap...")
    loadings_df = pd.DataFrame(pca.components_[:2].T, index=wide_df.columns, columns=['PC1', 'PC2'])
    plt.figure(figsize=(12, 8))
    sns.heatmap(loadings_df, annot=True, cmap="coolwarm", center=0)
    plt.title("PCA Loadings (PC1 & PC2)")
    plt.ylabel("Features")
    plt.xlabel("Principal Components")
    loadings_path = os.path.join(output_dir, "pca_loadings.png")
    plt.savefig(loadings_path, dpi=300)
    plt.close()
    logging.info(f"PCA loadings heatmap saved to {loadings_path}")

    # ----------------------
    # Covariance matrix heatmap
    # ----------------------
    logging.info("Generating covariance matrix heatmap...")
    cov_matrix = np.cov(X_scaled.T)
    cov_df = pd.DataFrame(cov_matrix, index=wide_df.columns, columns=wide_df.columns)
    plt.figure(figsize=(12, 10))
    sns.heatmap(cov_df, cmap="coolwarm", center=0)
    plt.title("Feature Covariance Matrix (All Bacteria Combined)")
    covariance_path = os.path.join(output_dir, "covariance_matrix.png")
    plt.savefig(covariance_path, dpi=300)
    plt.close()
    logging.info(f"Covariance matrix heatmap saved to {covariance_path}")

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
