# visualize_features.R
#
# Command-line R script to visualize distributions of significant features identified by KS test
# in mRNA target selection pipelines.
#
# Overview:
#   - Loads Kolmogorov-Smirnov (KS) test results to identify significant features/subfeatures (p < 0.05).
#   - Aggregates raw feature data for both positive and random sets for each significant feature/subfeature.
#   - Generates violin plots comparing distributions between positive and random sets.
#   - Saves plots as PNG files for downstream interpretation.
#
# Arguments/Paths:
#   ks_file (str): Path to CSV file containing KS test results.
#   base_path (str): Directory containing raw feature data CSVs.
#   output_dir (str): Directory to save output distribution plots.
#
# Requirements:
#   - Raw feature data files named as <feature>_positive_raw_data.csv and <feature>_random_raw_data.csv
#     under the specified base_path.
#   - KS test results CSV with columns: feature, subfeature, p_value.
#   - R packages: tidyverse, glue
#
# Usage Example:
#   Rscript visualize_features.R
#
# Outputs:
#   <output_dir>/<feature>_<subfeature>.png   # Violin plots
#   <summary_out>/ks_statistics_summary.png   # KS statistics summary plot
#
# Author: Nadia

library(tidyverse)
library(glue)

# -------------------------------
# Define paths and parameters
# -------------------------------
# List of pathogens to process
#pathogens <- c("S.aureus", "S.pneumoniae", "S.pyogenes", "C.trachomatis", 
#               "P.aeruginosa", "H.pylori", "N.gonorrhoeae", "C.burnetii", 
#               "B.melitensis")
pathogens <- c("B.melitensis")

# -------------------------------
# Run pipeline function
# -------------------------------

run_visualization <- function(pathogen_dir) {
  ks_file <- glue("../../results/{pathogen_dir}/ks_test_results_random.csv")
  base_path <- glue("../../results/{pathogen_dir}/raw_data")
  output_dir <- glue("../../results/{pathogen_dir}/feature_distributions")
  summary_out <- glue("../../results/{pathogen_dir}")
  
  # Read KS test results
  ks_results <- read_csv(ks_file, show_col_types = FALSE)
  feature_subfeature <- ks_results %>%
    filter(p_value < 0.05) %>%
    select(feature, subfeature) %>%
    distinct()
  
  # Load and plot distributions
  all_data <- load_all_data(feature_subfeature, base_path)
  plot_distributions(all_data, output_dir)
  plot_ks_statistics(ks_results, summary_out)
  
  message(glue("✅ Plots saved for {pathogen_dir} to {output_dir}"))
}

# -------------------------------
# Function to load raw data
# -------------------------------
load_all_data <- function(feature_subfeature, base_path) {
  data_list <- list()
  
  for (i in seq_len(nrow(feature_subfeature))) {
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
        warning(glue("Missing file: {file_path}"))
      }
    }
  }
  
  bind_rows(data_list)
}

# -------------------------------
# Function to plot distributions
# -------------------------------
plot_distributions <- function(data, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  unique_combos <- data %>% distinct(feature, subfeature)
  
  for (i in seq_len(nrow(unique_combos))) {
    feat <- unique_combos$feature[i]
    subf <- unique_combos$subfeature[i]
    
    df_plot <- data %>% filter(feature == feat, subfeature == subf)
    if (nrow(df_plot) == 0) next
    
    # Sanitize filename
    safe_feat <- str_replace_all(feat, "/", "-")
    safe_subf <- str_replace_all(subf, "/", "-")
    
    # Optional y-axis limits
    y_max <- max(df_plot$value, na.rm = TRUE)
    y_limits <- if (y_max <= 1.1) c(0, 1) else NULL
    
    p <- ggplot(df_plot, aes(x = set, y = value, fill = set)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      labs(title = paste(feat, "-", subf), x = NULL, y = "Value") +
      theme_minimal(base_size = 14) +
      scale_fill_brewer(palette = "Set2") +
      theme(legend.position = "none")
    
    if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
    
    ggsave(
      filename = file.path(output_dir, paste0(safe_feat, "_", safe_subf, ".png")),
      plot = p,
      width = 8,
      height = 5
    )
  }
}

# -------------------------------
# Categorize features for KS plot
# -------------------------------
categorize_feature <- function(feature, subfeature) {
  case_when(
    feature %in% c("signalp", "targetp", "deeplocpro", "deeptmhmm") ~ "Subcellular localisation",
    feature == "allergenicity" ~ "Allergenicity",
    feature == "ifnepitope2" ~ "Immunogenicity",
    feature %in% c("cluster_conservation", "rate4site", "rate4site_deeptmhmm", "dnds", "FEL", "FUBAR", "SLAC") |
      subfeature %in% c("Percent identity / number of strains", "Average Log₁₀ e-value", "Average bit-score / length") ~ "Conservation Analysis Across Strains",
    feature %in% c("bcell", "ellipro", "mchi", "mhcii", "mixmhc2pred") ~ "Epitope Prediction",
    feature == "popcov" ~ "Epitope evaluation",
    feature %in% c("ProtLearn", "dssp") ~ "Structure Analysis",
    TRUE ~ "Other"
  )
}

# -------------------------------
# Function to plot KS statistics with t-test directionality
# -------------------------------
plot_ks_statistics <- function(ks_df, output_dir) {
  if (!"ks_statistic" %in% colnames(ks_df)) stop("Column 'ks_statistic' not found in KS results.")
  
  # Filter significant features and set bar patterns
  ks_filtered <- ks_df %>%
    filter(p_value < 0.05) %>%
    mutate(
      label = paste(feature, subfeature, sep = " / "),
      label_safe = str_replace_all(label, "/", "-"),
      category = mapply(categorize_feature, feature, subfeature),
      t_direction = case_when(
        is.na(t_statistic) ~ "Unknown",
        t_statistic >= 0 ~ "Positive t",
        t_statistic < 0 ~ "Negative t"
      ),
      pattern = case_when(
        t_direction == "Positive t" ~ "solid",
        t_direction == "Negative t" ~ "striped",
        TRUE ~ "none"
      )
    ) %>%
    arrange(desc(ks_statistic))
  
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
  
  if (!requireNamespace("ggpattern", quietly = TRUE)) {
    install.packages("ggpattern")
  }
  library(ggpattern)
  
  p <- ggplot(ks_filtered, aes(
    x = reorder(label_safe, ks_statistic),
    y = ks_statistic,
    fill = category,
    pattern = t_direction
  )) +
    geom_col_pattern(
      color = "black",
      pattern_fill = "black",
      pattern_angle = 45,
      pattern_density = 0.05,
      pattern_spacing = 0.05
    ) +
    geom_text(aes(label = round(ks_statistic, 3)), hjust = -0.1, size = 3.5, color = "black") +
    coord_flip(ylim = c(0, 1)) +
    labs(
      title = "KS Statistic Summary (p < 0.05)",
      x = "Feature / Subfeature",
      y = "KS Statistic",
      fill = "Category",
      pattern = "T-test Direction"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold")) +
    scale_fill_manual(values = category_palette) +
    scale_pattern_manual(values = c("Positive t" = "none", "Negative t" = "stripe"))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(
    filename = file.path(output_dir, "ks_statistics_summary.png"),
    plot = p,
    width = 10,
    height = max(6, nrow(ks_filtered) * 0.25)
  )
}



# -------------------------------
# Run pipeline
# -------------------------------

for (pathogen in pathogens) {
  run_visualization(pathogen)
}
