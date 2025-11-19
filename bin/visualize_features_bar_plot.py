import pandas as pd
import matplotlib.pyplot as plt
import argparse

# -----------------------------
# Load your CSV file
# -----------------------------
parser = argparse.ArgumentParser(description="Visualize feature averages from a CSV file.")
parser.add_argument("--input-csv", help="Path to the input CSV file")
parser.add_argument("--output-path", help="Path to save the output bar plot image, default to results/feature_means_bar_plot.png", default=None)
args = parser.parse_args()

df = pd.read_csv(args.input_csv)
if output_path := args.output_path is None:
    output_path = "results/feature_means_bar_plot.png"

# -----------------------------
# Select numeric columns
# -----------------------------
numeric_df = df.select_dtypes(include='number')

# -----------------------------
# Compute averages of features
# -----------------------------
feature_means = numeric_df.mean().sort_values(ascending=False)

# -----------------------------
# Print sorted feature means
# -----------------------------
print("\n=== Sorted Feature Means ===")
print(feature_means)

# -----------------------------
# Visualization
# -----------------------------
plt.figure(figsize=(10, len(feature_means) * 0.2))

plt.barh(feature_means.index, feature_means.values)
plt.xlabel("Mean Value")
plt.title("Average of Features (Sorted)")
plt.gca().invert_yaxis()  # Highest value on top

plt.tight_layout()
plt.savefig(output_path)
