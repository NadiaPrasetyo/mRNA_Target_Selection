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
def setup_logging(verbose):
    """
    Set up logging configuration.
    if Verbose is true, it will print to both the CLI and a file
    """
    log_file = "pca_analysis.log" if verbose else None
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s - %(levelname)s - %(message)s"
    
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(level)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(format_str))
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is provided
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(format_str))
        logger.addHandler(file_handler)


# ----------------------
# Load data
# ----------------------
def load_bacterium_data(base_dir, bacterium):
    folder = os.path.join(base_dir, bacterium, "raw_data")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    dfs = []
    for f in files:
        df = pd.read_csv(f)
        df["bacterium"] = bacterium
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

def main(base_dir, output_dir, input_dirs):
    input_dirs = input_dirs if isinstance(input_dirs, list) else [input_dirs]

    # Load all bacteria data into one combined dataset
    all_data = pd.concat([load_bacterium_data(base_dir, b) for b in input_dirs], ignore_index=True)

    # ----------------------
    # Pivot to wide format
    # ----------------------
    # Each sample = combination of label + (row index across all bacteria)
    all_data["sample_id"] = (
        all_data.groupby(["label"]).cumcount().astype(str)
        + "_" + all_data["label"]
    )

    # Join feature and subfeature for column names
    all_data["feature_subfeature"] = all_data["feature"] + "_" + all_data["subfeature"]

    # Wide format: rows = samples, columns = feature_subfeature
    wide_df = all_data.pivot_table(
        index="sample_id",
        columns="feature_subfeature",
        values="value",
        aggfunc="mean"
    )

    # Metadata for labels
    meta = all_data.drop_duplicates(subset=["sample_id"])[["sample_id", "label"]].set_index("sample_id")

    # ----------------------
    # PCA
    # ----------------------
    X = wide_df.fillna(0)  # replace missing values
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA()
    pcs = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"], index=wide_df.index)
    pca_df = pca_df.join(meta)

    # ----------------------
    # Scree plot
    # ----------------------
    plt.figure(figsize=(8, 6))
    plt.plot(range(1, len(pca.explained_variance_ratio_) + 1), pca.explained_variance_ratio_, 'o-', markersize=6)
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    plt.title('Scree Plot')
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "scree_plot.png"), dpi=300)
    plt.close()

    # ----------------------
    # PCA biplot
    # ----------------------
    plt.figure(figsize=(10, 8))

    # Scatter plot of samples
    sns.scatterplot(
        data=pca_df,
        x="PC1", y="PC2",
        hue="label",
        s=60, alpha=0.8
    )

    # Add feature loadings for first 2 PCs
    loadings = pca.components_[:2].T * np.sqrt(pca.explained_variance_[:2])
    for i, col in enumerate(wide_df.columns):
        plt.arrow(0, 0, loadings[i, 0], loadings[i, 1],
                  color='gray', alpha=0.5, head_width=0.02)
        plt.text(loadings[i, 0]*1.1, loadings[i, 1]*1.1,
                 col, fontsize=8, alpha=0.7)

    plt.title("PCA Biplot (All Bacteria Combined)")
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.savefig(os.path.join(output_dir, "pca_biplot.png"), dpi=300)
    plt.close()

    # ----------------------
    # Loading plot (Top 2 PCs)
    # ----------------------
    loadings_df = pd.DataFrame(pca.components_[:2].T, index=wide_df.columns, columns=['PC1', 'PC2'])

    plt.figure(figsize=(12, 8))
    sns.heatmap(loadings_df, annot=True, cmap="coolwarm", center=0)
    plt.title("PCA Loadings (PC1 & PC2)")
    plt.ylabel("Features")
    plt.xlabel("Principal Components")
    plt.savefig(os.path.join(output_dir, "pca_loadings.png"), dpi=300)
    plt.close()

    # ----------------------
    # Covariance matrix heatmap
    # ----------------------
    cov_matrix = np.cov(X_scaled.T)
    cov_df = pd.DataFrame(cov_matrix, index=wide_df.columns, columns=wide_df.columns)

    plt.figure(figsize=(12, 10))
    sns.heatmap(cov_df, cmap="coolwarm", center=0)
    plt.title("Feature Covariance Matrix (All Bacteria Combined)")
    plt.savefig(os.path.join(output_dir, "covariance_matrix.png"), dpi=300)
    plt.close()

    print(f"✅ Analysis complete. Plots saved in {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCA Analysis of Bacterial Data")
    parser.add_argument("--output-dir", type=str, default="./results/pca_outputs", help="Directory to save output plots, default: ./results/pca_outputs")
    parser.add_argument("--base-dir", type=str, default="./results", help="Base directory for input data, default: ./results")
    parser.add_argument("--input-dir", required=True, nargs="+", type=str, help="List of bacteria directories (e.g. S.aureus S.pneumoniae S.pyogenes)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()
    setup_logging(args.verbose)
    os.makedirs(args.output_dir, exist_ok=True)
    main(args.base_dir, args.output_dir, args.input_dir)
