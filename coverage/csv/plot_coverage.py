#!/usr/bin/env python3
"""
Coverage visualization script for CRISPR-Cas9 experiment
Generates plots from coverage_summary_all_samples.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_CSV  = "coverage_summary_all_samples.csv"
OUTPUT_DIR = "coverage_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Color palette per group ───────────────────────────────────────────────────
GROUP_COLORS = {
    "Control":    "#2196F3",   # blue
    "Only_MNP":   "#4CAF50",   # green
    "Plasmid_Ko": "#FF9800",   # orange
    "RNP_Cas":    "#F44336",   # red
}

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)

# Shorten sample names for plotting
df["Sample_short"] = df["Sample"].str.replace("_L002", "").str.replace("_S\d+", "", regex=True)

# Assign colors
df["Color"] = df["Group"].map(GROUP_COLORS)

print(f"Loaded {len(df)} samples")
print(f"Groups: {df['Group'].unique()}")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Mean Depth per Sample (Bar plot)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 6))

bars = ax.bar(
    range(len(df)),
    df["Mean_Depth"],
    color=df["Color"],
    edgecolor="white",
    linewidth=0.5,
    alpha=0.85
)

# Add value labels on bars
for bar, val in zip(bars, df["Mean_Depth"]):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.5,
        f"{val:.1f}×",
        ha="center", va="bottom",
        fontsize=7, fontweight="bold", color="#333333"
    )

# Add reference line at mean
mean_depth = df["Mean_Depth"].mean()
ax.axhline(mean_depth, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.text(len(df) - 0.5, mean_depth + 0.5, f"Mean: {mean_depth:.1f}×",
        ha="right", fontsize=9, color="black", alpha=0.7)

# Formatting
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["Sample_short"], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Mean Depth (×)", fontsize=11)
ax.set_title("Mean Sequencing Depth at CRISPR-Cas9 Target Site (bdnf)\n"
             "NC_024333.1:15,921,539–15,922,558 (±500bp window)",
             fontsize=12, fontweight="bold", pad=15)
ax.set_ylim(0, df["Mean_Depth"].max() * 1.15)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")

# Legend
legend_patches = [
    mpatches.Patch(color=c, label=g)
    for g, c in GROUP_COLORS.items()
]
ax.legend(handles=legend_patches, loc="upper right",
          framealpha=0.9, fontsize=9)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot1_mean_depth_per_sample.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 1 saved: mean depth per sample")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Mean Depth by Group (Boxplot)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 6))

groups     = list(GROUP_COLORS.keys())
group_data = [df[df["Group"] == g]["Mean_Depth"].values for g in groups]
colors     = [GROUP_COLORS[g] for g in groups]

bp = ax.boxplot(
    group_data,
    patch_artist=True,
    widths=0.5,
    medianprops=dict(color="white", linewidth=2)
)

for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

for whisker in bp["whiskers"]:
    whisker.set(color="#555555", linewidth=1.2)
for cap in bp["caps"]:
    cap.set(color="#555555", linewidth=1.2)
for flier in bp["fliers"]:
    flier.set(marker="o", color="#555555", alpha=0.5)

# Overlay individual points
for i, (gdata, color) in enumerate(zip(group_data, colors), start=1):
    x = np.random.normal(i, 0.07, size=len(gdata))
    ax.scatter(x, gdata, color=color, edgecolor="white",
               linewidth=0.8, s=50, zorder=3, alpha=0.9)

ax.set_xticks(range(1, len(groups) + 1))
ax.set_xticklabels(groups, fontsize=10)
ax.set_ylabel("Mean Depth (×)", fontsize=11)
ax.set_title("Mean Depth Distribution by Experimental Group\nat CRISPR-Cas9 Target Site (bdnf)",
             fontsize=12, fontweight="bold", pad=15)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot2_depth_by_group_boxplot.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 2 saved: depth by group boxplot")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Coverage % and Total Reads per Sample
# ═════════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Top: Coverage %
ax1.bar(range(len(df)), df["Coverage_Pct"],
        color=df["Color"], edgecolor="white", linewidth=0.5, alpha=0.85)
ax1.axhline(100, color="green", linestyle="--", linewidth=1, alpha=0.6,
            label="100% coverage")
ax1.set_ylabel("Coverage (%)", fontsize=11)
ax1.set_title("Coverage % and Total Reads at CRISPR-Cas9 Target Site (bdnf)",
              fontsize=12, fontweight="bold", pad=15)
ax1.set_ylim(0, 115)
ax1.spines[["top", "right"]].set_visible(False)
ax1.set_facecolor("#f9f9f9")
ax1.legend(fontsize=9)

for i, val in enumerate(df["Coverage_Pct"]):
    ax1.text(i, val + 1, f"{val:.0f}%",
             ha="center", va="bottom", fontsize=7, fontweight="bold")

# Bottom: Total reads
ax2.bar(range(len(df)), df["Total_Reads"],
        color=df["Color"], edgecolor="white", linewidth=0.5, alpha=0.85)
ax2.set_ylabel("Total Reads", fontsize=11)
ax2.set_xticks(range(len(df)))
ax2.set_xticklabels(df["Sample_short"], rotation=45, ha="right", fontsize=8)
ax2.spines[["top", "right"]].set_visible(False)
ax2.set_facecolor("#f9f9f9")

for i, val in enumerate(df["Total_Reads"]):
    ax2.text(i, val + 2, str(val),
             ha="center", va="bottom", fontsize=7, fontweight="bold")

# Shared legend
legend_patches = [
    mpatches.Patch(color=c, label=g)
    for g, c in GROUP_COLORS.items()
]
ax2.legend(handles=legend_patches, loc="upper right",
           framealpha=0.9, fontsize=9)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot3_coverage_pct_and_reads.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 3 saved: coverage % and total reads")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Heatmap: Mean Depth per Sample and Group
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(5, 8))

# Pivot para heatmap
pivot = df.set_index("Sample_short")[["Mean_Depth", "Total_Reads",
                                      "Coverage_Pct", "Mean_BaseQ",
                                      "Mean_MapQ"]].copy()

# Normalizar por columna para el heatmap
pivot_norm = (pivot - pivot.min()) / (pivot.max() - pivot.min())

im = ax.imshow(pivot_norm.values, cmap="YlOrRd", aspect="auto",
               vmin=0, vmax=1)

# Etiquetas
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(["Mean\nDepth", "Total\nReads", "Coverage\n%",
                     "Mean\nBaseQ", "Mean\nMapQ"],
                   fontsize=9)
ax.set_yticks(range(len(pivot)))
ax.set_yticklabels(pivot.index, fontsize=8)

# Valores en celdas
for i in range(len(pivot)):
    for j, col in enumerate(pivot.columns):
        val = pivot.iloc[i, j]
        text = f"{val:.1f}" if col != "Total_Reads" else str(int(val))
        ax.text(j, i, text, ha="center", va="center",
                fontsize=7, color="black" if pivot_norm.iloc[i, j] < 0.7
                else "white")

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Normalized value", fontsize=9)

# Color de fila por grupo
for i, (_, row) in enumerate(df.iterrows()):
    ax.add_patch(mpatches.Rectangle(
        (-0.5 + len(pivot.columns), i - 0.5), 0.3, 1,
        color=GROUP_COLORS[row["Group"]], clip_on=False
    ))

ax.set_title("Coverage Metrics Heatmap\nCRISPR-Cas9 Target Site (bdnf)",
             fontsize=11, fontweight="bold", pad=15)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot4_heatmap_coverage_metrics.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 4 saved: heatmap")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Mean Depth vs Mean MapQ (scatter)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 6))

for group, color in GROUP_COLORS.items():
    subset = df[df["Group"] == group]
    ax.scatter(
        subset["Mean_MapQ"],
        subset["Mean_Depth"],
        c=color, label=group,
        s=100, edgecolor="white",
        linewidth=0.8, alpha=0.9, zorder=3
    )
    # Etiquetas de muestra
    for _, row in subset.iterrows():
        ax.annotate(
            row["Sample_short"],
            (row["Mean_MapQ"], row["Mean_Depth"]),
            textcoords="offset points",
            xytext=(5, 3),
            fontsize=6, color="#555555"
        )

ax.set_xlabel("Mean Mapping Quality (MapQ)", fontsize=11)
ax.set_ylabel("Mean Depth (×)", fontsize=11)
ax.set_title("Mean Depth vs Mapping Quality\nat CRISPR-Cas9 Target Site (bdnf)",
             fontsize=12, fontweight="bold", pad=15)
ax.legend(fontsize=9, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")
ax.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot5_depth_vs_mapq.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 5 saved: depth vs mapq scatter")

# ═════════════════════════════════════════════════════════════════════════════
# Summary stats
# ═════════════════════════════════════════════════════════════════════════════
print("\n=== Summary Statistics by Group ===")
summary = df.groupby("Group").agg(
    N=("Sample", "count"),
    Mean_Depth_mean=("Mean_Depth", "mean"),
    Mean_Depth_std=("Mean_Depth", "std"),
    Total_Reads_mean=("Total_Reads", "mean"),
    Coverage_Pct_mean=("Coverage_Pct", "mean"),
    Mean_MapQ_mean=("Mean_MapQ", "mean")
).round(2)
print(summary.to_string())

summary.to_csv(f"{OUTPUT_DIR}/coverage_summary_by_group.csv")
print(f"\n✅ Summary stats saved: {OUTPUT_DIR}/coverage_summary_by_group.csv")
print(f"\n✅ All plots saved to: {OUTPUT_DIR}/")
print("   plot1_mean_depth_per_sample.png")
print("   plot2_depth_by_group_boxplot.png")
print("   plot3_coverage_pct_and_reads.png")
print("   plot4_heatmap_coverage_metrics.png")
print("   plot5_depth_vs_mapq.png")