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
ks_file <- "../../results/S.aureus/ks_test_results.csv"
base_path <- "../../results/S.aureus/raw_data"
output_dir <- "../../results/S.aureus/feature_distributions"

# Read KS test results and filter for significant features with p-value < 0.05
feature_subfeature <- read_csv(ks_file, show_col_types = FALSE) %>%
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
    
    p <- ggplot(df_plot, aes(x = set, y = value, fill = set)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      labs(title = paste(feat, "-", subf), x = "", y = "Value") +
      theme_minimal(base_size = 14) +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "none")
    
    ggsave(filename = file.path(output_dir, paste0(feat, "_", subf, ".png")),
           plot = p, width = 8, height = 5)
  }
}

# Run pipeline
all_data <- load_all_data(feature_subfeature)
plot_distributions(all_data, output_dir)

message("✅ Plots saved to: ", output_dir)
