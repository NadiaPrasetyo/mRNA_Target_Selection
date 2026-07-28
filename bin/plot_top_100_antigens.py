"""
Plot: Top 100 Predicted Staphylococcus aureus Antigen Summary
Boxplot of Antigenicity by functional Category, sorted left-to-right by
descending median antigenicity, colored with a colorblind-safe palette.

Reads data from a CSV file with columns: Category, Antigenicity, n
  - One row per data point.
  - The "n" column only needs to be filled in on the FIRST row of each
    category (the reported sample size for that group); it can be left
    blank on subsequent rows of the same category.

Usage:
    python plot_antigenicity.py [path_to_csv]

If no path is given, defaults to "antigenicity_data.csv" in the same
folder as this script.

Requirements: matplotlib, numpy, pandas
    pip install matplotlib numpy pandas
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import argparse

# ---------------------------------------------------------------
# 1. Load data from CSV
# ---------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("csv_path", nargs="?", default=Path(__file__).parent / "antigenicity_data.csv")
args = parser.parse_args()
csv_path = args.csv_path

df = pd.read_csv(csv_path)
df.columns = [c.strip() for c in df.columns]  # tolerate stray whitespace in headers

required_cols = {"Category", "Antigenicity", "n"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"CSV is missing required column(s): {missing}")

# Forward-fill n within each category's block of rows (n is only given
# on the first row of each category in the CSV)
df["n"] = df.groupby("Category", sort=False)["n"].ffill()

# Build: Category -> list of antigenicity values, and Category -> n
data = {
    cat: group["Antigenicity"].tolist()
    for cat, group in df.groupby("Category", sort=False)
}
n_reported = {
    cat: int(group["n"].iloc[0])
    for cat, group in df.groupby("Category", sort=False)
}

# ---------------------------------------------------------------
# 2. Sort categories by descending median antigenicity
# ---------------------------------------------------------------
categories_sorted = sorted(
    data.keys(), key=lambda cat: np.median(data[cat]), reverse=True
)
values_sorted = [data[cat] for cat in categories_sorted]
n_labels = [f"n={n_reported[cat]}" for cat in categories_sorted]

# ---------------------------------------------------------------
# 2b. Color palette (kept consistent with other plots)
# ---------------------------------------------------------------
# NOTE: the palette's keys don't correspond to these antigen categories,
# so we use its color VALUES, one per category, in the sorted
# (left-to-right) order -- keeping the same look/feel as your other
# figures. Since there can be more categories than colors in the
# original palette, it's extended with additional colorblind-safe
# colors from the same family (Okabe-Ito / Paul Tol qualitative
# palettes) so every category gets its own distinct color.
extended_palette_colors = [
    "#0072B2",
    "#D55E00",
    "#56B4E9",
    "#CC79A7",
    "#009E73",
    "#E69F00",
    "#999999",
    "#0072B2",
    "#D55E00",
    "#56B4E9",
    "#CC79A7",
    "#009E73",
    "#E69F00",
    "#999999",
]
palette_colors = extended_palette_colors
if len(categories_sorted) > len(palette_colors):
    raise ValueError(
        f"Not enough colors in palette ({len(palette_colors)}) for "
        f"{len(categories_sorted)} categories. Add more colors to "
        f"extended_palette_colors."
    )
box_colors = [palette_colors[i] for i in range(len(categories_sorted))]

# ---------------------------------------------------------------
# 3. Font setup (Times New Roman, falling back to Liberation Serif
#    which is metrically identical, in case TNR isn't installed)
# ---------------------------------------------------------------
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = [
    "Times New Roman",
    "Liberation Serif",
    "Times",
    "DejaVu Serif",
]
matplotlib.rcParams["mathtext.fontset"] = "stix"  # Times-like serif math/italic font

BASE_FONT_SIZE = 20

# ---------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(15, 10))

box = ax.boxplot(
    values_sorted,
    showmeans=True,
    meanline=False,
    meanprops=dict(marker="x", markeredgecolor="black", markersize=9, markeredgewidth=1.5),
    medianprops=dict(color="black", linewidth=1.2),
    boxprops=dict(color="black", linewidth=1.2),
    whiskerprops=dict(color="black", linewidth=1.2),
    capprops=dict(color="black", linewidth=1.2),
    flierprops=dict(
        marker="o",
        markerfacecolor="none",
        markeredgecolor="black",
        markersize=7,
    ),
    widths=0.5,
    patch_artist=True,
)

# Color the box faces; keep edges/whiskers/etc. black for readability
for patch, color in zip(box["boxes"], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
    patch.set_edgecolor("black")
    patch.set_linewidth(1.2)

# Title (with italic species name, like the original)
ax.set_title(
    r"Top 100 Predicted $\it{Staphylococcus\ aureus}$ Antigen Summary",
    fontsize=BASE_FONT_SIZE + 12,
    pad=20,
)

ax.set_ylabel("Antigenicity", fontsize=BASE_FONT_SIZE + 2)

# X tick labels: category names, rotated, colored to match their box
ax.set_xticks(range(1, len(categories_sorted) + 1))
ax.set_xticklabels(categories_sorted, rotation=45, ha="right", fontsize=BASE_FONT_SIZE)


# n= labels just under the axis, above the rotated category names
for i, label in enumerate(n_labels, start=1):
    ax.text(
        i,
        -0.028,
        label,
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=BASE_FONT_SIZE,
    )

# Push category-name tick labels down to make room for the n= row
ax.tick_params(axis="x", pad=22)

y_min = min(min(v) for v in values_sorted)
y_max = max(max(v) for v in values_sorted)
ax.set_ylim(np.floor(y_min * 100) / 100 - 0.005, np.ceil(y_max * 100) / 100 + 0.005)
ax.tick_params(axis="y", labelsize=BASE_FONT_SIZE)

ax.yaxis.grid(True, color="lightgray", linewidth=0.8)
ax.set_axisbelow(True)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()

# ---------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------
plt.savefig("results/antigenicity_boxplot.png", dpi=300, bbox_inches="tight")
print("Saved plot to results/antigenicity_boxplot.png")
# plt.show()