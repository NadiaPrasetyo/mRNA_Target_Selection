# visualize_features.R
#
# Command-line R script to visualize distributions of significant features identified by KS test
# in mRNA target selection pipelines.
#
# Overview:
#   - Loads Kolmogorov-Smirnov (KS) test results to identify significant features/subfeatures (p < 0.05).
#   - Aggregates raw feature data for both positive and random sets for each significant feature/subfeature.
#   - Generates violin plots comparing the distributions of each feature/subfeature between positive and random sets.
#   - Saves plots as PNG files for downstream interpretation.
#
# Arguments/Paths:
#   ks_file (str): Path to CSV file containing KS test results.
#   base_path (str): Directory containing raw feature data CSVs for positive and random sets.
#   output_dir (str): Directory to save output distribution plots.
#
# Requirements:
#   - Raw feature data files named as <feature>_positive_raw_data.csv and <feature>_random_raw_data.csv
#     under the specified base_path.
#   - KS test results CSV with columns: feature, subfeature, p_value.
#   - R packages: tidyverse (ggplot2, dplyr, readr, etc.)
#
# Usage Example:
#   Rscript visualize_features.R
#
# Outputs:
#   <output_dir>/<feature>_<subfeature>.png   # Violin plots for each significant feature/subfeature
#
# Author: Nadia
# This file visualizes the features collected by the pipeline
# Load libraries
library(tidyverse)

# Define paths
ks_file <- "../../results/C.trachomatis/ks_test_results.csv"
base_path <- "../../results/C.trachomatis/raw_data"
output_dir <- "../../results/C.trachomatis/feature_distributions"

# Read KS test results and filter for significant features with p-value < 0.05
ks_results <- read_csv(ks_file, show_col_types = FALSE)
feature_subfeature <- ks_results %>%
  filter(p_value < 0.05) %>%
  select(feature, subfeature) %>%
  distinct()

# Function to load all data
load_all_data <- function(feature_subfeature) {
  data_list <- list()
  
  for (i in 1:nrow(feature_subfeature)) {
    feature <- feature_subfeature$feature[i]
    subfeature <- feature_subfeature$subfeature[i]
    
    for (set_label in c("positive", "random")) {
      file_path <- file.path(base_path, paste0(feature, "_", set_label, "_raw_data.csv"))
      
      if (file.exists(file_path)) {
        df <- read_csv(file_path, show_col_types = FALSE) %>%
          filter(subfeature == !!subfeature) %>%
          mutate(set = set_label)
        
        data_list[[length(data_list) + 1]] <- df
      } else {
        warning(paste("Missing file:", file_path))
      }
    }
  }
  
  bind_rows(data_list)
}

# Function to plot distributions
plot_distributions <- function(data, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  unique_combos <- data %>% distinct(feature, subfeature)
  
  for (i in 1:nrow(unique_combos)) {
    feat <- unique_combos$feature[i]
    subf <- unique_combos$subfeature[i]
    
    df_plot <- data %>% filter(feature == feat, subfeature == subf)
    
    if (nrow(df_plot) == 0) next
    
    # Sanitize filename: replace slashes with hyphens
    safe_feat <- str_replace_all(feat, "/", "-")
    safe_subf <- str_replace_all(subf, "/", "-")
    
    # Set y-axis limits if values should be bounded
    y_max <- max(df_plot$value, na.rm = TRUE)
    y_limits <- if (y_max <= 1.1) c(0, 1) else NULL
    
    p <- ggplot(df_plot, aes(x = set, y = value, fill = set)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      labs(title = paste(feat, "-", subf), x = "", y = "Value") +
      theme_minimal(base_size = 14) +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "none")
    
    if (!is.null(y_limits)) {
      p <- p + coord_cartesian(ylim = y_limits)
    }
    
    ggsave(
      filename = file.path(output_dir, paste0(safe_feat, "_", safe_subf, ".png")),
      plot = p,
      width = 8,
      height = 5
    )
  }
}

# Function to categorize features for coloring
categorize_feature <- function(feature, subfeature) {
  # Subcellular localisation
  if (feature %in% c("signalp", "targetp", "deeplocpro", "deeptmhmm")) {
    return("Subcellular localisation")
  }
  # Allergenicity
  if (feature == "allergenicity") {
    return("Allergenicity")
  }
  # Immunogenicity
  if (feature == "ifnepitope2") {
    return("Immunogenicity")
  }
  # Conservation Analysis
  if (feature %in% c("cluster_conservation", "rate4site", "rate4site_deeptmhmm", "dnds", "FEL", "FUBAR", "SLAC") ||
      subfeature %in% c("Percent identity / number of strains", "Average Log₁₀ e-value", "Average bit-score / length")) {
    return("Conservation Analysis Across Strains")
  }
  # Epitope Prediction
  if (feature %in% c("bcell", "ellipro", "mchi", "mhcii", "mixmhc2pred")) {
    return("Epitope Prediction")
  }
  # Epitope evaluation
  if (feature == "popcov") {
    return("Epitope evaluation")
  }

  if (feature %in% c("ProtLearn", "dssp")) {
    return("Structure Analysis")
  }
  return("Other")
}

# Function to visualize KS statistics (with value labels and categorized colors)
plot_ks_statistics <- function(ks_df, output_dir) {
  if (!"ks_statistic" %in% colnames(ks_df)) {
    stop("Column 'ks_statistic' not found in KS results.")
  }
  
  ks_filtered <- ks_df %>%
    filter(p_value < 0.05) %>%
    mutate(
      label = paste(feature, subfeature, sep = " / "),
      label_safe = str_replace_all(label, "/", "-"),
      category = mapply(categorize_feature, feature, subfeature)
    ) %>%
    arrange(desc(ks_statistic))
  
  # Define color palette for categories
  category_palette <- c(
    "Subcellular localisation" = "#1b9e77",
    "Allergenicity" = "#d95f02",
    "Immunogenicity" = "#7570b3",
    "Conservation Analysis Across Strains" = "#e7298a",
    "Epitope Prediction" = "#66a61e",
    "Epitope evaluation" = "#e6ab02",
    "Structure Analysis" = "#d010e1",
    "Other" = "#a6761d"
  )
  
  p <- ggplot(ks_filtered, aes(x = reorder(label_safe, ks_statistic), y = ks_statistic, fill = category)) +
    geom_col() +
    geom_text(aes(label = round(ks_statistic, 3)), 
              hjust = -0.1, size = 3.5, color = "black") +
    coord_flip(ylim = c(0, 1)) +
    labs(
      title = "KS Statistic Summary (p < 0.05)",
      x = "Feature / Subfeature",
      y = "KS Statistic",
      fill = "Category"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold")) +
    scale_fill_manual(values = category_palette)
  
  ggsave(
    filename = file.path(output_dir, "ks_statistics_summary.png"),
    plot = p,
    width = 10,
    height = max(6, nrow(ks_filtered) * 0.25)  # Adjust height based on number of bars
  )
}

# Run pipeline
all_data <- load_all_data(feature_subfeature)
plot_distributions(all_data, output_dir)
plot_ks_statistics(ks_results, output_dir)

message("✅ Plots saved to: ", output_dir)
