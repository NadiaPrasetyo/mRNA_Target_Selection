"""
pca_analysis.py
Runner for performing PCA and statistical analysis on bacterial data.
Overview:
    - Loads raw data for specified bacteria from CSV files.
    - Preprocesses and aggregates data by accession.
    - Encodes features and accessions into sparse matrices for efficient computation.
    - Performs Z-score normalization and PCA on the data.
    - Computes statistical tests (KS test, t-test, AUROC) for feature significance.
    - Generates various plots for PCA results, feature loadings, and statistical summaries.
    - Saves all results and plots in a structured directory format for downstream analysis.
Arguments:
    --base-dir (str): Base directory containing the data for each bacterium (default: "./results").
    --output-dir (str): Directory where results and plots will be saved (default: "./results").
    --input-dir (list[str]): List of subdirectories under `base-dir` for each bacterium (required).
    --verbose: If set, enables verbose logging to console and log file.
Requirements:
    - Python packages: argparse, pandas, numpy, matplotlib, seaborn, sklearn, scipy, adjustText.
    - Input data must be in CSV format with specific columns (e.g., accession, feature, subfeature, value, label).
    - Ensure that the `adjustText` package is installed for better label placement in plots.
Outputs:
    <output_dir>/                                      # Directory containing all results and plots.
    <output_dir>/scree_plot.png                       # Scree plot showing variance explained by PCA components.
    <output_dir>/pca_biplot.png                       # PCA biplot with samples and representative feature vectors.
    <output_dir>/feature_loadings.csv                 # CSV file with PCA feature loadings.
    <output_dir>/representative_features.csv          # CSV file with representative features for PCA biplot.
    <output_dir>/ks_auroc_results.csv                 # CSV file with KS test, t-test, and AUROC results.
    <output_dir>/auroc_summary_all.png                # AUROC summary plot for significant features.
    <output_dir>/ks_summary_all.png                   # KS statistics summary plot for significant features.
    <output_dir>/correlation_matrix.png               # Heatmap of feature correlation matrix.
    <output_dir>/covariance_matrix.png                # Heatmap of feature covariance matrix.
Notes:
    - The script assumes that each bacterium's data is stored in a subdirectory under `base-dir`.
    - Aggregated data and labels are saved as debug CSVs for inspection.
    - Statistical tests are performed to compare "positive" labels against "random" labels.
    - Ensure that the input data is properly formatted and contains the required columns.
Author: Nadia
"""
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
from scipy.stats import ks_2samp, ttest_ind, spearmanr
from sklearn.metrics import roc_auc_score
import matplotlib.patches as mpatches
from scipy.sparse import issparse

# ----------------------
# Configuration
# ----------------------
def setup_logging(verbose: bool) -> None:
    """Configure logging to console and optionally file."""
    log_file = "pca_analysis.log" if verbose else None
    level = logging.INFO
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
    """Load and concatenate all raw data CSVs for a bacterium.
    Args:
        base_dir (str): Base directory containing bacterium subdirectories.
        bacterium (str): Name of the bacterium (subdirectory name).
    Returns:
        pd.DataFrame: Concatenated DataFrame of all raw data for the bacterium.
    """
    folder = os.path.join(base_dir, bacterium, "raw_data")
    logging.info(f"Loading data for bacterium: {bacterium} from {folder}")
    files = glob.glob(os.path.join(folder, "*_raw_data.csv"))
    dfs = []
    for f in files:
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


def aggregate_by_accession(df: pd.DataFrame, output_prefix: str = "debug") -> tuple[pd.DataFrame, pd.Series]:
    """
    Pivot so each accession is a single row, each column is a feature_subfeature.
    The value is typically the mean across replicates.

    Only includes accessions with more than 4 characters.
    Also writes debug CSVs for the aggregated data and labels.
    """
    logging.info("Aggregating data by accession...")

    df = df.copy()

    # Filter accessions with more than 4 characters
    df = df[df["accession"].astype(str).str.len() > 4]
    if df.empty:
        logging.warning("No accessions with more than 4 characters found.")
        return pd.DataFrame(), pd.Series(dtype=object)

    # Make a combined feature identifier
    df["feature_subfeature"] = df["feature"].astype(str) + "_" + df["subfeature"].astype(str)

    # Aggregate values by accession + feature
    agg_df = (
        df.groupby(["accession", "feature_subfeature"], observed=True)["value"]
          .mean()
          .reset_index()
    )

    # Keep a label for each accession (e.g. majority vote)
    label_map = (
        df.groupby("accession")["label"]
          .agg(lambda x: x.mode().iat[0])
    )

    # Write debug CSVs
    agg_df.to_csv(f"results/{output_prefix}_aggregated.csv", index=False)
    label_map.to_csv(f"results/{output_prefix}_labels.csv", header=True)

    logging.info(f"Aggregated data written to results/{output_prefix}_aggregated.csv")
    logging.info(f"Label map written to results/{output_prefix}_labels.csv")

    return agg_df, label_map

# ----------------------
# Feature categorization
# ----------------------
def categorize_feature(feature, subfeature):
    """Categorize features for AUROC/KS summary plots.
    Categories:
        - Subcellular localisation
        - Allergenicity
        - Immunogenicity
        - Conservation Analysis Across Strains
        - Epitope Prediction
        - Structure Analysis
        - Other
    """
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
        return "Conservation"
    if feature in ["bcell", "ellipro", "mhci", "mhcii", "mixmhc2pred", "discotope"]:
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
import re

# Extend this with your actual feature keys for full coverage.
# Falls back to a title-cased version of the raw feature name if not found.
FEATURE_DISPLAY_NAMES = {
    "deeplocpro": "DeepLocPro -",
    "deeptmhmm": "DeepTMHMM -",
    "signalp": "SignalP -",
    "bcell": "BepiPred -",
    "targetp": "TargetP -",
    "allergenicity": "AlgPred -",
    "discotope": "DiscoTope -",
    "cluster": "Cluster -",
    "mhci": "NetMHCPan -",
    "mhcii": "NetMHCIIPan -"
    
}


def _display_feature_name(feature):
    return FEATURE_DISPLAY_NAMES.get(feature.lower(), feature.replace("_", " ").title())


def _shorten_subfeature(subfeature):
    """Clean up a subfeature string: strip common prefixes, title-case,
    and convert 'proportion_x' patterns into 'X %'."""
    s = subfeature.strip()
    is_proportion = bool(re.match(r'^(proportion|prop|)_', s, flags=re.IGNORECASE))
    s = re.sub(r'^(prob|proportion|prop|num)_', '', s, flags=re.IGNORECASE)
    s = s.replace('_', ' ').strip()
    s = s.replace('percent identity/num strain', 'PID Across Strains')
    if s.lower() not in ["sp", "notp", "pid across strains"]:
        s = s.title()
    if is_proportion:
        s = f"{s} %"
    return s


def _shorten_label(feature, subfeature):
    return f"{_display_feature_name(feature)} {_shorten_subfeature(subfeature)}"

def plot_auroc_summary(results_df, output_dir, prefix="all"):
    """
    Plot AUROC summary and Top 20 plot for significant features.
    """
    output_path_all = os.path.join(output_dir, f"auroc_summary_{prefix}.png")
    output_path_top20 = os.path.join(output_dir, f"auroc_summary_top20_{prefix}.png")
    os.makedirs(output_dir, exist_ok=True)

    # Set Times New Roman globally for this plot
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]

    df = results_df.dropna(subset=["auroc", "t_pvalue", "t_statistic"]).copy()
    df = df[(df["auroc"] != 0.5) & (df["t_pvalue"] < 0.05)]

    if df.empty:
        logging.warning("No significant AUROC values to plot.")
        return

    df["adjusted_auroc"] = df["auroc"].apply(lambda x: x if x >= 0.5 else 1 - x)
    df["label"] = df["feature"] + " / " + df["subfeature"]
    df["category"] = df.apply(lambda row: categorize_feature(row["feature"], row["subfeature"]), axis=1)
    df = df.sort_values("adjusted_auroc", ascending=False)

    category_palette = {
        "Subcellular localisation": "#0072B2",
        "Allergenicity": "#D55E00",
        "Immunogenicity": "#56B6E9" if False else "#56B4E9",
        "Conservation": "#CC79A7",
        "Epitope Prediction": "#009E73",
        "Structure Analysis": "#E69F00",
        "Other": "#999999"
    }

    def _plot(df_subset, save_path, title_suffix="", short_labels=False, bar_height=0.6,
              title_fs=24, label_fs=16, tick_fs=13, legend_fs=12, legend_title_fs=13, value_fs=13):
        colors = df_subset["category"].map(category_palette).fillna("#a6761d")
        hatches = ['' if t >= 0 else '////' for t in df_subset["t_statistic"]]

        labels = (
            df_subset.apply(lambda row: _shorten_label(row["feature"], row["subfeature"]), axis=1)
            if short_labels
            else df_subset["label"]
        )

        fig_height = max(4, bar_height * len(df_subset))
        plt.figure(figsize=(12, fig_height))
        bars = plt.barh(labels, df_subset["adjusted_auroc"], color=colors)

        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)

        for bar in bars:
            width = bar.get_width()
            plt.text(
                width + 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{width:.3f}",
                va="center",
                ha="left",
                fontsize=value_fs,
                color="black",
                backgroundcolor="white",
                bbox=dict(
                    facecolor="white", edgecolor="none", boxstyle="square,pad=0.1"
                ),
                clip_on=False
            )

        handles = [mpatches.Patch(color=color, label=cat)
                   for cat, color in category_palette.items()]
        handles += [
            mpatches.Patch(facecolor='white', edgecolor='black', hatch='////', label='Enriched in Random'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='Enriched in Positive')
        ]

        plt.legend(handles=handles, title="Category / Directionality", loc="lower right",
                   fontsize=legend_fs, title_fontsize=legend_title_fs)
        plt.xlabel("AUROC", fontsize=label_fs)
        plt.ylabel("")
        plt.title(f"AUROC Summary {title_suffix}".strip(), fontsize=title_fs, fontweight='bold')
        plt.xticks(fontsize=tick_fs)
        plt.yticks(fontsize=tick_fs)

        x_max = max(df_subset["adjusted_auroc"]) + 0.1
        plt.xlim(0.5, x_max)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

    # "All features": full labels, original spacing/fonts
    _plot(df, output_path_all, "(Significant Features)", short_labels=False, bar_height=0.3)
    logging.info(f"AUROC summary plot saved to {output_path_all}")

    # Top 20: shortened labels, taller bars, larger fonts
    df_top20 = df.head(20)
    if not df_top20.empty:
        _plot(
            df_top20, output_path_top20, "(Top 20 Features)",
            short_labels=True, bar_height=0.9,
            title_fs=24, label_fs=20, tick_fs=17, legend_fs=15, legend_title_fs=16, value_fs=16
        )
        logging.info(f"AUROC top 20 plot saved to {output_path_top20}")


def plot_ks_summary(results_df, output_dir, prefix="all"):
    """
    Plot KS summary and Top 20 plot for significant features.
    """
    output_path_all = os.path.join(output_dir, f"ks_summary_{prefix}.png")
    output_path_top20 = os.path.join(output_dir, f"ks_summary_top20_{prefix}.png")
    os.makedirs(output_dir, exist_ok=True)

    # Set Times New Roman globally for this plot
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]

    df = results_df.dropna(subset=["ks_statistic", "ks_pvalue", "t_statistic"]).copy()
    df = df[df["ks_pvalue"] < 0.05]

    if df.empty:
        logging.warning("No significant KS statistics to plot.")
        return

    df["label"] = df["feature"] + " / " + df["subfeature"]
    df["category"] = df.apply(lambda row: categorize_feature(row["feature"], row["subfeature"]), axis=1)
    df = df.sort_values("ks_statistic", ascending=False)

    category_palette = {
        "Subcellular localisation": "#0072B2",
        "Allergenicity": "#D55E00",
        "Immunogenicity": "#56B4E9",
        "Conservation": "#CC79A7",
        "Epitope Prediction": "#009E73",
        "Structure Analysis": "#E69F00",
        "Other": "#999999"
    }

    def _plot(df_subset, save_path, title_suffix="", short_labels=False, bar_height=0.6,
              title_fs=24, label_fs=16, tick_fs=13, legend_fs=12, legend_title_fs=13, value_fs=13):
        colors = df_subset["category"].map(category_palette).fillna("#a6761d")
        hatches = ['' if t >= 0 else '////' for t in df_subset["t_statistic"]]

        labels = (
            df_subset.apply(lambda row: _shorten_label(row["feature"], row["subfeature"]), axis=1)
            if short_labels
            else df_subset["label"]
        )

        fig_height = max(4, bar_height * len(df_subset))
        plt.figure(figsize=(12, fig_height))

        bars = plt.barh(labels, df_subset["ks_statistic"], color=colors)

        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)

        for bar in bars:
            width = bar.get_width()
            plt.text(
                width + 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{width:.3f}",
                va="center",
                ha="left",
                fontsize=value_fs,
                color="black",
                backgroundcolor="white",
                bbox=dict(
                    facecolor="white", edgecolor="none", boxstyle="square,pad=0.1"
                ),
                clip_on=False
            )

        handles = [mpatches.Patch(color=color, label=cat)
                   for cat, color in category_palette.items()]
        handles += [
            mpatches.Patch(facecolor='white', edgecolor='black', hatch='////', label='Enriched in Random'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='Enriched in Positive')
        ]

        plt.legend(handles=handles, title="Category / Directionality", loc="lower right",
                   fontsize=legend_fs, title_fontsize=legend_title_fs)
        plt.xlabel("KS Statistic", fontsize=label_fs)
        plt.ylabel("")
        plt.title(f"KS Statistics Summary {title_suffix}".strip(), fontsize=title_fs, fontweight='bold')
        plt.xticks(fontsize=tick_fs)
        plt.yticks(fontsize=tick_fs)

        x_max = max(df_subset["ks_statistic"]) + 0.1
        plt.xlim(0, x_max)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

    # "All features": full labels, original spacing/fonts
    _plot(df, output_path_all, "(Significant Features)", short_labels=False, bar_height=0.3)
    logging.info(f"KS summary plot saved to {output_path_all}")

    # Top 20: shortened labels, taller bars, larger fonts
    df_top20 = df.head(20)
    if not df_top20.empty:
        _plot(
            df_top20, output_path_top20, "(Top 20 Features)",
            short_labels=True, bar_height=0.9,
            title_fs=24, label_fs=20, tick_fs=17, legend_fs=15, legend_title_fs=16, value_fs=16
        )
        logging.info(f"KS top 20 plot saved to {output_path_top20}")

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



def plot_pca_biplot(pca_df, ipca, feature_enc, output_dir: str, scale=2.5, n_directions=12):
    """
    PCA biplot with samples and representative feature vectors.
    Instead of labeling all features, we choose one representative
    per direction (with the largest vector length).
    
    Also exports a CSV of all feature loadings.
    
    Parameters
    ----------
    pca_df : pd.DataFrame
        DataFrame containing PC1, PC2, and sample labels.
    ipca : fitted PCA object (e.g., sklearn.decomposition.PCA)
    feature_enc : encoder with inverse_transform to retrieve feature names
    output_dir : str
        Directory to save plot and CSV
    scale : float
        Base scaling factor for arrows
    n_directions : int
        Number of angular sectors to partition the unit circle.
        One representative feature is selected from each sector.
    """
    plt.figure(figsize=(12, 10))

    # =========================
    # Color palette for samples
    # =========================
    unique_labels = pca_df['label'].unique()
    if 'random' in unique_labels:
        n_other_labels = len(unique_labels) - 1
        
        # Use a perceptually uniform blue palette for non-random samples
        other_colors = sns.light_palette("#0072B2", n_other_labels, reverse=True)  # Okabe–Ito blue
        
        custom_palette = {}
        other_idx = 0
        for label in unique_labels:
            if label == 'random':
                # Use a color-blind–safe red/orange (high contrast, distinct in grayscale)
                custom_palette[label] = "#D55E00"   # Okabe–Ito vermilion
            else:
                custom_palette[label] = other_colors[other_idx]
                other_idx += 1
    else:
        custom_palette = sns.light_palette("#0072B2", as_cmap=True)

    scatter = sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="label",
        s=60, alpha=0.7, palette=custom_palette
    )
    
    handles, labels = scatter.get_legend_handles_labels()

    # =========================
    # Feature loadings
    # =========================
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    # Compute direction (angle) and distance (norm) for each feature
    norms = np.linalg.norm(loadings, axis=1)
    angles = np.arctan2(loadings[:, 1], loadings[:, 0])  # radians [-pi, pi]

    # Normalize angles to [0, 2π)
    angles = (angles + 2 * np.pi) % (2 * np.pi)

    # DataFrame for all features
    loadings_df = pd.DataFrame({
        "feature": feature_names,
        "PC1_loading": loadings[:, 0],
        "PC2_loading": loadings[:, 1],
        "angle_rad": angles,
        "distance_from_origin": norms
    })

    # =========================
    # Select representatives by direction
    # =========================
    sector_edges = np.linspace(0, 2 * np.pi, n_directions + 1)
    representatives = []

    for i in range(n_directions):
        start, end = sector_edges[i], sector_edges[i + 1]
        in_sector = loadings_df[(loadings_df['angle_rad'] >= start) &
                                (loadings_df['angle_rad'] < end)]
        if not in_sector.empty:
            # pick feature with max distance in this direction
            rep = in_sector.loc[in_sector['distance_from_origin'].idxmax()]
            representatives.append(rep)

    rep_df = pd.DataFrame(representatives)

    # =========================
    # Dynamic scaling
    # =========================
    x_range = pca_df["PC1"].max() - pca_df["PC1"].min()
    y_range = pca_df["PC2"].max() - pca_df["PC2"].min()
    avg_range = (x_range + y_range) / 2
    dynamic_scale = scale * (avg_range / 5)

    texts = []

    # =========================
    # Axes and formatting
    # =========================
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.7)
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--', alpha=0.7)

    plt.xlabel(f"PC1 ({ipca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=12)
    plt.ylabel(f"PC2 ({ipca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=12)
    plt.title(f"PCA Biplot (Representative Features by Direction, {n_directions} sectors)", fontsize=14, pad=20)

    plt.grid(True, linestyle='--', alpha=0.3)

    plt.legend(
        handles=handles, labels=labels,
        title="Label",
        loc='upper left',
        bbox_to_anchor=(1.05, 1),
        borderaxespad=0.
    )

    plt.tight_layout()

    # =========================
    # Save outputs
    # =========================
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "pca_biplot.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # Save all loadings and representatives
    loadings_df.to_csv(os.path.join(output_dir, "feature_loadings.csv"), index=False)
    rep_df.to_csv(os.path.join(output_dir, "representative_features.csv"), index=False)

    print(f"Saved PCA biplot to {os.path.join(output_dir, 'pca_biplot.png')}")
    print(f"Saved full loadings to {os.path.join(output_dir, 'feature_loadings.csv')}")
    print(f"Saved representative features to {os.path.join(output_dir, 'representative_features.csv')}")



def plot_loading_scatter(ipca, feature_enc, output_dir: str, top_n=50):
    """Scatter plot of feature loadings on PC1 vs PC2 with non-overlapping labels."""
    loadings = ipca.components_[:2].T
    feature_names = feature_enc.inverse_transform(np.arange(loadings.shape[0]))

    norms = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(norms)[-top_n:]
    
    top_features = [feature_names[i] for i in top_idx]
    top_loadings = loadings[top_idx]

    plt.figure(figsize=(12, 10))
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    
    plt.scatter(loadings[:, 0], loadings[:, 1], alpha=0.3, color="lightgray", s=20, label="Other features")
    plt.scatter(top_loadings[:, 0], top_loadings[:, 1], color="red", s=80, label=f"Top {top_n} features")
    
    text_objects = []
    for i, (x, y) in enumerate(top_loadings):
        offset_x = 0.02 * (1 if x >= 0 else -1)
        offset_y = 0.02 * (1 if y >= 0 else -1)
        
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
    
    plt.savefig(os.path.join(output_dir, "pca_loading_scatter.png"), dpi=300, bbox_inches='tight')
    plt.close()

# ----------------------
# Correlation matrix plotting
# ----------------------
def plot_correlation_matrix(
    X_scaled, 
    feature_enc, 
    output_dir: str, 
    max_features=2000,
    subset_prefixes=("discotope_", "ellipro_", "dssp_", "ProtLearn_", "FEL_", "FUBAR_", "SLAC_", "cluster_percent_")
):
    """
    Plot Spearman correlation heatmap (clustered) for all features,
    and a separate heatmap for selected subsets of features.
    """
    # Set Times New Roman globally for these plots
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]

    os.makedirs(output_dir, exist_ok=True)
    n_features = X_scaled.shape[1]

    if n_features > max_features:
        logging.warning(f"Too many features ({n_features}), downsampling to top {max_features} by variance.")
        if issparse(X_scaled):
            col_var = np.array(
                X_scaled.power(2).mean(axis=0) - np.power(X_scaled.mean(axis=0), 2)
            ).ravel()
        else:
            col_var = X_scaled.var(axis=0)

        # Always keep features in subset_prefixes (the ones the subset
        # heatmap is built from) so variance-based downsampling can't
        # crowd them out with high-variance, unrelated features
        # (e.g. thousands of per-site conservation columns).
        all_feature_names = feature_enc.inverse_transform(np.arange(n_features))
        protected_mask = np.array([
            any(f.startswith(prefix) for prefix in subset_prefixes)
            for f in all_feature_names
        ])
        protected_idx = np.where(protected_mask)[0]
        n_protected = len(protected_idx)

        if n_protected >= max_features:
            logging.warning(
                f"{n_protected} protected subset features already meet/exceed "
                f"max_features={max_features}; keeping protected features only."
            )
            top_idx = protected_idx
        else:
            remaining_budget = max_features - n_protected
            candidate_idx = np.where(~protected_mask)[0]
            candidate_var = col_var[candidate_idx]
            if len(candidate_idx) > remaining_budget:
                top_candidates = candidate_idx[
                    np.argpartition(candidate_var, -remaining_budget)[-remaining_budget:]
                ]
            else:
                top_candidates = candidate_idx
            top_idx = np.concatenate([protected_idx, top_candidates])

        top_idx = np.sort(top_idx)
        X_scaled = X_scaled[:, top_idx]
        feature_names = feature_enc.inverse_transform(top_idx)
    else:
        feature_names = feature_enc.inverse_transform(np.arange(n_features))

    if issparse(X_scaled):
        X_dense = X_scaled.toarray().astype(np.float32)
    else:
        X_dense = X_scaled.astype(np.float32)

    logging.info("Computing Spearman correlation matrix...")
    corr_matrix, _ = spearmanr(X_dense, axis=0)
    corr_df = pd.DataFrame(corr_matrix, index=feature_names, columns=feature_names)

    nan_cols = corr_df.columns[corr_df.isna().all()]
    if len(nan_cols) > 0:
        logging.warning(f"Constant features with NaN correlation dropped: {list(nan_cols)}")
        corr_df = corr_df.drop(index=nan_cols, columns=nan_cols)

    logging.info("Plotting clustered correlation heatmap...")
    sns.set(style="white")
    # sns.set() resets rcParams, so re-apply the font after it
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
    g = sns.clustermap(
        corr_df,
        cmap="coolwarm",
        center=0,
        figsize=(21, 18),
        xticklabels=False,
        yticklabels=False,
        row_cluster=True,
        col_cluster=True,
        dendrogram_ratio=(0.1, 0),
        cbar_pos=None
    )

    clustered_features = corr_df.index[g.dendrogram_row.reordered_ind]

    max_labels = 80
    step = max(1, len(clustered_features) // max_labels)

    g.ax_heatmap.set_yticks(np.arange(0, len(clustered_features), step))
    g.ax_heatmap.set_yticklabels(
        clustered_features[::step],
        fontsize=8,
        rotation=0
    )

    g.ax_heatmap.yaxis.set_label_position("right")
    g.ax_heatmap.yaxis.tick_right()

    g.ax_heatmap.tick_params(
        left=False,
        labelleft=False,
        right=True,
        labelright=True,
        bottom=False,
        labelbottom=False,
        labelsize=8,
        pad=2
    )

    g.figure.suptitle(
        "Spearman Correlation Matrix (Clustered)",
        fontsize=26, x=0.5, ha="center", y=1.02
    )
    plt.savefig(os.path.join(output_dir, "correlation_matrix_spearman.png"), dpi=600, bbox_inches="tight")
    plt.close()

    corr_df.to_csv(os.path.join(output_dir, "correlation_matrix_spearman.csv"))

    logging.info("Saved full correlation heatmap.")

    logging.info("Creating subset correlation heatmap...")

    selected_features = [
        f for f in corr_df.columns
        if any(f.startswith(prefix) for prefix in subset_prefixes)
    ]

    if len(selected_features) < 2:
        logging.warning("Not enough subset features found to compute correlation heatmap.")
        return

    corr_subset = corr_df.loc[selected_features, selected_features]

    nan_cols_subset = corr_subset.columns[corr_subset.isna().all()]
    if len(nan_cols_subset) > 0:
        logging.warning(f"Subset constant features dropped: {list(nan_cols_subset)}")
        corr_subset = corr_subset.drop(index=nan_cols_subset, columns=nan_cols_subset)

    if corr_subset.shape[0] < 2:
        logging.warning("Subset correlation matrix has < 2 features after cleaning.")
        return

    sns.set(style="white")
    # sns.set() resets rcParams, so re-apply the font after it
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
    g_sub = sns.clustermap(
        corr_subset,
        cmap="coolwarm",
        center=0,
        figsize=(18, 14),
        xticklabels=False,
        yticklabels=False,
        row_cluster=True,
        col_cluster=True,
        dendrogram_ratio=(0.1, 0),
        cbar_pos=(-0.02, 0.8, 0.05, 0.18),
    )

    clustered_features = corr_subset.index[g_sub.dendrogram_row.reordered_ind]

    max_labels = 80
    step = max(1, len(clustered_features) // max_labels)

    g_sub.ax_heatmap.set_yticks(np.arange(0, len(clustered_features), step))
    g_sub.ax_heatmap.set_yticklabels(
        clustered_features[::step],
        fontsize=20,
        rotation=0
    )

    g_sub.ax_heatmap.yaxis.set_label_position("right")
    g_sub.ax_heatmap.yaxis.tick_right()

    g_sub.ax_heatmap.tick_params(
        left=False,
        labelleft=False,
        right=True,
        labelright=True,
        bottom=False,
        labelbottom=False,
        labelsize=20,
        pad=2
    )
    
    if g_sub.cax is not None:
        g_sub.cax.tick_params(labelsize=20)


    g_sub.figure.suptitle(
        "Spearman Correlation (Selected Feature Subset)",
        fontsize=32, x=0.5, ha="center", y=1.02
    )

    plt.savefig(
        os.path.join(output_dir, "correlation_matrix_spearman_subset.png"),
        dpi=600, bbox_inches="tight"
    )
    plt.close()

    corr_subset.to_csv(
        os.path.join(output_dir, "correlation_matrix_spearman_subset.csv")
    )

    logging.info("Saved subset correlation heatmap and matrix.")

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
    logging.info(f"{all_data['feature'].nunique()} features detected: \n{all_data['feature'].unique()}")

    # ----------------------
    # Preprocess
    # ----------------------
    all_data = all_data[all_data["value"].notna()]
    all_data["feature_subfeature"] = (
        all_data["feature"].astype(str) + "_" + all_data["subfeature"].astype(str)
    )

    # ----------------------
    # Compute KS/AUROC BEFORE PCA (feature selection happens here)
    # ----------------------
    logging.info("Computing KS/AUROC statistics (before PCA)...")
    results_df = compute_stats(all_data)
    results_df.to_csv(os.path.join(output_dir, "ks_auroc_results.csv"), index=False)

    # Threshold for PCA feature selection
    ks_threshold = 0.15
    significant = results_df.loc[results_df["ks_statistic"] > ks_threshold].copy()

    significant["feature_subfeature"] = (
        significant["feature"].astype(str) + "_" + significant["subfeature"].astype(str)
    )
    selected_features = set(significant["feature_subfeature"])

    logging.info(f"{len(selected_features):,} features passed KS > {ks_threshold}")

    # ----------------------
    # Aggregate by accession
    # (done unconditionally now so the full, unfiltered feature set is
    #  always available for the correlation matrix, matching the original
    #  behaviour, even when the KS filter yields no PCA-eligible features.)
    # ----------------------
    agg_df, accession_labels = aggregate_by_accession(all_data)

    logging.info(f"Aggregated data shape: {agg_df.shape}")
    logging.info(f"Unique features: {agg_df['feature_subfeature'].nunique():,}")

    # ----------------------
    # Encode rows (accessions) and columns (features) on the FULL,
    # unfiltered feature set -- this is what the correlation matrix
    # is built from, reverted back to original (VER1) behaviour.
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
    logging.info(f"X_sparse unique features: {len(feature_enc.classes_):,}")

    # ----------------------
    # Metadata for plotting
    # ----------------------
    meta = pd.DataFrame({
        "accession": accession_enc.classes_,
        "label": [accession_labels[a] for a in accession_enc.classes_]
    }).set_index("accession")

    # ----------------------
    # Z-score Normalization (per feature) on the FULL feature set
    # ----------------------
    logging.info("Applying Z-score normalization per feature...")
    scaler = StandardScaler(with_mean=False)  # works with sparse
    X_scaled = scaler.fit_transform(X_sparse)

    logging.info(f"X_scaled shape after scaling: {X_scaled.shape}, nnz={X_scaled.nnz:,}")

    # ----------------------
    # Correlation matrix -- reverted to use the full, unfiltered
    # X_scaled/feature_enc, exactly as in the original script.
    # ----------------------
    plot_correlation_matrix(X_scaled, feature_enc, output_dir)

    if len(selected_features) == 0:
        logging.warning("No features exceed KS threshold. Skipping PCA.")
        plot_ks_summary(results_df, output_dir)
        plot_auroc_summary(results_df, output_dir)
        return

    # ----------------------
    # Build a KS-filtered matrix for PCA only (does not affect the
    # correlation matrix, which already ran above on the full data).
    # ----------------------
    agg_df_filtered = agg_df[agg_df["feature_subfeature"].isin(selected_features)]
    logging.info(f"After KS filtering: {len(agg_df_filtered['feature_subfeature'].unique()):,} features remain for PCA")

    accession_enc_pca = LabelEncoder()
    feature_enc_pca = LabelEncoder()

    row_idx_pca = accession_enc_pca.fit_transform(agg_df_filtered["accession"])
    col_idx_pca = feature_enc_pca.fit_transform(agg_df_filtered["feature_subfeature"])

    X_sparse_pca = coo_matrix(
        (agg_df_filtered["value"].astype(np.float32), (row_idx_pca, col_idx_pca)),
        shape=(len(accession_enc_pca.classes_), len(feature_enc_pca.classes_))
    ).tocsr()

    meta_pca = pd.DataFrame({
        "accession": accession_enc_pca.classes_,
        "label": [accession_labels[a] for a in accession_enc_pca.classes_]
    }).set_index("accession")

    scaler_pca = StandardScaler(with_mean=False)
    X_scaled_pca = scaler_pca.fit_transform(X_sparse_pca)

    # ----------------------
    # PCA
    # ----------------------
    logging.info("Running PCA on KS-filtered data...")
    num_features = X_scaled_pca.shape[1]
    n_components = min(50, num_features)

    if n_components < 2:
        logging.error(f"Not enough features ({num_features}) for PCA. Need at least 2.")
        plot_ks_summary(results_df, output_dir)
        plot_auroc_summary(results_df, output_dir)
        return

    logging.info(f"Using n_components={n_components} for PCA (num_features={num_features})")
    ipca = IncrementalPCA(n_components=n_components, batch_size=10000)
    pcs = ipca.fit_transform(X_scaled_pca)
    pca_df = pd.DataFrame(
        pcs[:, :2], columns=["PC1", "PC2"], index=accession_enc_pca.classes_
    ).join(meta_pca)
    logging.info(f"PCA output shape: {pcs.shape}")

    # ----------------------
    # Plots
    # ----------------------
    plot_scree(ipca, output_dir)
    plot_pca_biplot(pca_df, ipca, feature_enc_pca, output_dir)
    plot_loading_scatter(ipca, feature_enc_pca, output_dir)

    # KS and AUROC summary plots (computed earlier)
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

    #B.melitensis C.burnetii C.trachomatis H.influenzae H.pylori P.aeruginosa N.gonorrhoeae S.enteritidis S.pneumoniae S.pyogenes T.pallidum V.cholerae