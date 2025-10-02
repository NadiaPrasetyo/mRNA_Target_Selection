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
from scipy.stats import ks_2samp, ttest_ind
from sklearn.metrics import roc_auc_score
import matplotlib.patches as mpatches

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
    """Drop missing and non-positive values."""
    return df.dropna()[df["value"] > 0]


def aggregate_by_accession(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot so each accession is a single row, each column is a feature_subfeature.
    The value is typically the mean across replicates.
    """
    logging.info("Aggregating data by accession...")

    # Make a combined feature identifier
    df = df.copy()
    df["feature_subfeature"] = df["feature"].astype(str) + "_" + df["subfeature"].astype(str)

    # Aggregate values by accession + feature
    agg_df = (
        df.groupby(["accession", "feature_subfeature"], observed=True)["value"]
          .mean()                       # <-- choose mean or sum as needed
          .reset_index()
    )

    # Keep a label for each accession (e.g. majority vote)
    label_map = (
        df.groupby("accession")["label"]
          .agg(lambda x: x.mode().iat[0])  # pick most common label per accession
    )
    return agg_df, label_map

# ----------------------
# Feature categorization
# ----------------------
def categorize_feature(feature, subfeature):
    """Categorize features for AUROC/KS summary plots."""
    if feature in ["signalp", "targetp", "deeplocpro", "deeptmhmm"]:
        return "Subcellular localisation"
    if feature == "allergenicity":
        return "Allergenicity"
    if feature == "ifnepitope2":
        return "Immunogenicity"
    if feature in ["cluster_conservation", "rate4site", "rate4site_deeptmhmm", "FEL", "SLAC", "FUBAR"] or subfeature in [
        "Percent identity / number of strains",
        "Average Log₁₀ e-value",
        "Average bit-score / length",
        "percent_identity/num_strain",
        "e_value_average",
        "bit_score_normalized",
        "rolling_mean",
        "rolling_median",
        "rolling_max",
        "rolling_min",
        "outside_score_per_site",
        "per_site_score"
    ]:
        return "Conservation Analysis Across Strains"
    if feature in ["bcell", "ellipro", "mhci", "mhcii", "mixmhc2pred"]:
        return "Epitope Prediction"
    if feature in ["dssp", "ProtLearn"]:
        return "Structure Analysis"
    return "Other"

# ----------------------
# KS Test, t-test, and AUROC
# ----------------------
def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute KS, t-test, and AUROC for each feature/subfeature."""
    results = []
    grouped = df.groupby(["feature", "subfeature"])
    for (feature, subfeature), group in grouped:
        pos_vals = group.loc[group["label"] != "random", "value"].values
        rand_vals = group.loc[group["label"] == "random", "value"].values
        if len(pos_vals) == 0 or len(rand_vals) == 0:
            continue
        ks_stat, ks_p = ks_2samp(pos_vals, rand_vals)
        t_stat, t_p = ttest_ind(pos_vals, rand_vals, equal_var=False)
        y_true = np.concatenate([np.ones(len(pos_vals)), np.zeros(len(rand_vals))])
        y_scores = np.concatenate([pos_vals, rand_vals])
        try:
            auroc = roc_auc_score(y_true, y_scores)
        except ValueError:
            auroc = np.nan
        results.append({
            "feature": feature,
            "subfeature": subfeature,
            "ks_statistic": ks_stat,
            "ks_pvalue": ks_p,
            "t_statistic": t_stat,
            "t_pvalue": t_p,
            "auroc": auroc
        })
    return pd.DataFrame(results)

# ----------------------
# AUROC and KS plots
# ----------------------
def plot_auroc_summary(results_df, output_dir, prefix="all"):
    output_path = os.path.join(output_dir, f"auroc_summary_{prefix}.png")
    os.makedirs(output_dir, exist_ok=True)

    df = results_df.dropna(subset=["auroc", "t_pvalue", "t_statistic"]).copy()
    df = df[(df["auroc"] != 0.5) & (df["t_pvalue"] < 0.05)]

    if df.empty:
        logging.warning("No significant AUROC values to plot.")
        return

    # Adjust AUROCs < 0.5
    df["adjusted_auroc"] = df["auroc"].apply(lambda x: x if x >= 0.5 else 1 - x)
    df["label"] = df["feature"] + " / " + df["subfeature"]
    df["category"] = df.apply(lambda row: categorize_feature(row["feature"], row["subfeature"]), axis=1)
    df = df.sort_values("adjusted_auroc", ascending=False)

    category_palette = {
        "Subcellular localisation": "#1b9e77",
        "Allergenicity": "#d95f02",
        "Immunogenicity": "#7570b3",
        "Conservation Analysis Across Strains": "#e7298a",
        "Epitope Prediction": "#66a61e",
        "Epitope evaluation": "#e6ab02",
        "Structure Analysis": "#d010e1",
        "Other": "#a6761d"
    }

    colors = df["category"].map(category_palette).fillna("#a6761d")
    hatches = ['' if t >= 0 else '////' for t in df["t_statistic"]]

    plt.figure(figsize=(10, max(4, 0.3 * len(df))))
    bars = plt.barh(df["label"], df["adjusted_auroc"], color=colors)

    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.01, bar.get_y() + bar.get_height()/2, f"{width:.3f}", va="center", fontsize=9)

    handles = [mpatches.Patch(color=color, label=cat) for cat, color in category_palette.items()]
    handles += [
        mpatches.Patch(facecolor='white', edgecolor='black', hatch='////', label='Enriched in Random'),
        mpatches.Patch(facecolor='white', edgecolor='black', label='Enriched in Positive')
    ]

    plt.legend(handles=handles, title="Category / Directionality", loc="lower right", fontsize=9)
    plt.xlabel("AUROC (adjusted, min=0.5)")
    plt.title("AUROC Summary (Significant Features)")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    logging.info(f"AUROC summary plot saved to {output_path}")


def plot_ks_summary(results_df, output_dir, prefix="all"):
    output_path = os.path.join(output_dir, f"ks_summary_{prefix}.png")
    os.makedirs(output_dir, exist_ok=True)

    df = results_df.dropna(subset=["ks_statistic", "ks_pvalue", "t_statistic"]).copy()
    df = df[df["ks_pvalue"] < 0.05]

    if df.empty:
        logging.warning("No significant KS statistics to plot.")
        return

    df["label"] = df["feature"] + " / " + df["subfeature"]
    df["category"] = df.apply(lambda row: categorize_feature(row["feature"], row["subfeature"]), axis=1)
    df = df.sort_values("ks_statistic", ascending=False)

    category_palette = {
        "Subcellular localisation": "#1b9e77",
        "Allergenicity": "#d95f02",
        "Immunogenicity": "#7570b3",
        "Conservation Analysis Across Strains": "#e7298a",
        "Epitope Prediction": "#66a61e",
        "Epitope evaluation": "#e6ab02",
        "Structure Analysis": "#d010e1",
        "Other": "#a6761d"
    }

    colors = df["category"].map(category_palette).fillna("#a6761d")
    hatches = ['' if t >= 0 else '////' for t in df["t_statistic"]]

    plt.figure(figsize=(10, max(4, 0.3 * len(df))))
    bars = plt.barh(df["label"], df["ks_statistic"], color=colors)

    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.01, bar.get_y() + bar.get_height()/2, f"{width:.3f}", va="center", fontsize=9)

    handles = [mpatches.Patch(color=color, label=cat) for cat, color in category_palette.items()]
    handles += [
        mpatches.Patch(facecolor='white', edgecolor='black', hatch='////', label='Enriched in Random'),
        mpatches.Patch(facecolor='white', edgecolor='black', label='Enriched in Positive')
    ]

    plt.legend(handles=handles, title="Category / Directionality", loc="lower right", fontsize=9)
    plt.xlabel("KS Statistic")
    plt.title("KS Statistics Summary (Significant Features)")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    logging.info(f"KS summary plot saved to {output_path}")


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


def plot_pca_biplot(pca_df, ipca, feature_enc, output_dir: str, top_n=20, scale=2.5):
    """PCA biplot with samples and top feature loadings using manual label placement."""
    plt.figure(figsize=(12, 10))
    
    # Create a custom palette where 'random' is red and others use viridis
    unique_labels = pca_df['label'].unique()
    if 'random' in unique_labels:
        # Create custom palette: random=red, others=viridis
        n_other_labels = len(unique_labels) - 1
        other_colors = sns.color_palette("viridis", n_other_labels)
        
        custom_palette = {}
        other_idx = 0
        for label in unique_labels:
            if label == 'random':
                custom_palette[label] = 'red'
            else:
                custom_palette[label] = other_colors[other_idx]
                other_idx += 1
    else:
        custom_palette = "viridis"
    
    # Create scatter plot with samples using custom palette
    scatter = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="label", 
                             s=60, alpha=0.7, palette=custom_palette)
    
    # Store legend handles and labels
    handles, labels = scatter.get_legend_handles_labels()
    
    # Loadings
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    # Top features by vector length
    norms = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(norms)[-top_n:]
    
    # Calculate dynamic scaling based on data range
    x_range = pca_df["PC1"].max() - pca_df["PC1"].min()
    y_range = pca_df["PC2"].max() - pca_df["PC2"].min()
    avg_range = (x_range + y_range) / 2
    dynamic_scale = scale * (avg_range / 5)  # Adjust scale based on data range

    texts = []
    
    # Plot feature vectors and labels
    for i in top_idx:
        x, y = loadings[i, 0] * dynamic_scale, loadings[i, 1] * dynamic_scale
        
        # Draw arrow using plt.arrow instead of FancyArrowPatch
        plt.arrow(0, 0, x, y, color='red', alpha=0.7, 
                  head_width=0.03*dynamic_scale, 
                  length_includes_head=True, 
                  linewidth=1.5,
                  overhang=0.3)  # Added overhang to improve arrow appearance
        
        # Add text label with initial positioning
        text = plt.text(x * 1.15, y * 1.15, feature_names[i], fontsize=10, 
                        color="darkred", weight='bold',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", 
                                 alpha=0.9, edgecolor="red", linewidth=0.5))
        texts.append(text)

    # Adjust labels if adjust_text is available
    try:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color='gray', lw=0.8, alpha=0.7),
            expand_points=(1.5, 1.8),
            expand_text=(1.3, 1.6),
            force_points=(0.5, 0.8),
            force_text=(0.8, 1.2),
            va='center',
            ha='center',
            only_move={'points':'xy', 'text':'xy', 'objects':'xy'},
            avoid_points=True,
            avoid_text=True,
            lim=100  # Increase iteration limit for better convergence
        )
    except ImportError:
        # Fallback: manual label placement if adjust_text is not available
        print("adjustText package not available, using manual label placement")
        for text in texts:
            # Simple manual adjustment to reduce overlaps
            pos = text.get_position()
            text.set_position((pos[0] + 0.02, pos[1] + 0.02))

    # Add center lines
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.7)
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.7)
    
    # Add labels and title
    plt.xlabel(f"PC1 ({ipca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=12)
    plt.ylabel(f"PC2 ({ipca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=12)
    plt.title("PCA Biplot (Top Features)", fontsize=14, pad=20)
    
    # Add grid for better readability
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Add legend using the stored handles and labels
    plt.legend(handles=handles, labels=labels, 
               title="Label",
               loc='upper left',
               bbox_to_anchor=(1.05, 1),
               borderaxespad=0.)
    
    plt.tight_layout()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "pca_biplot.png"), dpi=300, bbox_inches='tight')
    plt.close()


def plot_loading_scatter(ipca, feature_enc, output_dir: str, top_n=20):
    """Scatter plot of feature loadings on PC1 vs PC2 with non-overlapping labels."""
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    norms = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(norms)[-top_n:]
    
    # Get the top feature names and their loadings
    top_features = [feature_names[i] for i in top_idx]
    top_loadings = loadings[top_idx]

    plt.figure(figsize=(12, 10))
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    
    # Plot all points in light gray
    plt.scatter(loadings[:, 0], loadings[:, 1], alpha=0.3, color="lightgray", s=20, label="Other features")
    
    # Plot top points in red
    plt.scatter(top_loadings[:, 0], top_loadings[:, 1], color="red", s=80, label=f"Top {top_n} features")
    
    # Add labels with adjustments to prevent overlap
    text_objects = []
    for i, (x, y) in enumerate(top_loadings):
        # Determine the optimal text position
        offset_x = 0.02 * (1 if x >= 0 else -1)
        offset_y = 0.02 * (1 if y >= 0 else -1)
        
        # Create annotation with arrow
        ann = plt.annotate(
            top_features[i], 
            xy=(x, y), 
            xytext=(x + offset_x, y + offset_y),
            fontsize=10,
            ha='left' if x >= 0 else 'right',
            va='bottom' if y >= 0 else 'top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0", color="black", alpha=0.6, lw=0.8)
        )
        text_objects.append(ann)
    
    # Adjust text positions to minimize overlaps
    adjust_text(text_objects, 
                arrowprops=dict(arrowstyle="->", color='black', lw=0.5),
                expand_points=(1.5, 1.5),
                expand_text=(1.2, 1.2),
                force_text=(0.5, 0.8),
                only_move={'points':'y', 'text':'xy', 'objects':'xy'})
    
    plt.xlabel("Loading on PC1", fontsize=12)
    plt.ylabel("Loading on PC2", fontsize=12)
    plt.title(f"PCA Feature Loadings (Top {top_n})", fontsize=14, pad=20)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    
    # Save with high DPI for better quality
    plt.savefig(os.path.join(output_dir, "pca_loading_scatter.png"), dpi=300, bbox_inches='tight')
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
    logging.info("Starting PCA + KS/AUROC analysis...")
    os.makedirs(output_dir, exist_ok=True)

    # ----------------------
    # Load data
    # ----------------------
    dfs = [load_bacterium_data(base_dir, b) for b in input_dirs]
    dfs = [df for df in dfs if not df.empty]
    if not dfs:
        logging.error("No data loaded. Exiting.")
        return
    all_data = pd.concat(dfs, ignore_index=True)
    logging.info(f"Loaded data for {len(input_dirs)} bacteria. Rows: {len(all_data):,}")

    # ----------------------
    # Preprocess
    # ----------------------
    all_data = all_data[all_data["value"].notna() & (all_data["value"] > 0)]
    all_data["feature_subfeature"] = all_data["feature"].astype(str) + "_" + all_data["subfeature"].astype(str)

    # ----------------------
    # Aggregate by accession
    # ----------------------
    agg_df, accession_labels = aggregate_by_accession(all_data)

    # ----------------------
    # Encode rows (accessions) and columns (features)
    # ----------------------
    accession_enc = LabelEncoder()
    feature_enc = LabelEncoder()

    row_idx = accession_enc.fit_transform(agg_df["accession"])
    col_idx = feature_enc.fit_transform(agg_df["feature_subfeature"])

    X_sparse = coo_matrix(
        (agg_df["value"].astype(np.float32), (row_idx, col_idx)),
        shape=(len(accession_enc.classes_), len(feature_enc.classes_))
    ).tocsr()

    logging.info(f"Sparse matrix shape: {X_sparse.shape}, nnz={X_sparse.nnz:,}")

    # ----------------------
    # Metadata for plotting
    # ----------------------
    meta = pd.DataFrame({
        "accession": accession_enc.classes_,
        "label": [accession_labels[a] for a in accession_enc.classes_]
    }).set_index("accession")

    # ----------------------
    # PCA
    # ----------------------
    scaler = StandardScaler(with_mean=False)
    X_scaled = scaler.fit_transform(X_sparse)

    ipca = IncrementalPCA(n_components=50, batch_size=10000)
    pcs = ipca.fit_transform(X_scaled)
    pca_df = pd.DataFrame(pcs[:, :2], columns=["PC1", "PC2"], index=accession_enc.classes_).join(meta)

    # ----------------------
    # Plots
    # ----------------------
    plot_scree(ipca, output_dir)
    plot_pca_biplot(pca_df, ipca, feature_enc, output_dir)
    plot_loading_scatter(ipca, feature_enc, output_dir)
    plot_covariance_matrix(X_scaled, feature_enc, output_dir)
    plot_correlation_matrix(X_scaled, feature_enc, output_dir)

    # KS/t-test/AUROC computations
    results_df = compute_stats(all_data)
    results_df.to_csv(os.path.join(output_dir, "ks_auroc_results.csv"), index=False)

    # KS and AUROC plots
    plot_ks_summary(results_df, output_dir)
    plot_auroc_summary(results_df, output_dir)

    logging.info(f"✅ Analysis complete. All plots and CSVs saved in {output_dir}")

# ----------------------
# Entry point
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run PCA and KS/AUROC analysis on bacterial data.")
    parser.add_argument("--base-dir", type=str, default="./results")
    parser.add_argument("--output-dir", type=str, default="./results")
    parser.add_argument("--input-dir", nargs="+", required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(args.verbose)
    main(args.base_dir, args.output_dir, args.input_dir)
