#!/usr/bin/env python3
"""
Depth per position visualization script for CRISPR-Cas9 experiment
Generates plots from depth_per_position_all_samples.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_CSV  = "depth_per_position_all_samples.csv"
OUTPUT_DIR = "depth_position_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── CRISPR site coordinates ───────────────────────────────────────────────────
SGRNA_START = 15922039
SGRNA_END   = 15922058
CUT_SITE    = (SGRNA_START + SGRNA_END) // 2   # ~15922048

# ── Sample → Group mapping ────────────────────────────────────────────────────
def get_group(sample):
    if "Control" in sample:   return "Control"
    if "Only_MNP" in sample:  return "Only_MNP"
    if "Plasmid" in sample:   return "Plasmid_Ko"
    if "RNP" in sample:       return "RNP_Cas"
    return "Unknown"

GROUP_COLORS = {
    "Control":    "#2196F3",
    "Only_MNP":   "#4CAF50",
    "Plasmid_Ko": "#FF9800",
    "RNP_Cas":    "#F44336",
}

GROUP_ALPHA = {
    "Control":    0.8,
    "Only_MNP":   0.8,
    "Plasmid_Ko": 0.8,
    "RNP_Cas":    0.8,
}

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"Loading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV)
df = df.set_index("Position")

# Assign groups
sample_groups = {col: get_group(col) for col in df.columns}
print(f"Loaded {len(df)} positions × {len(df.columns)} samples")
print(f"Position range: {df.index.min()} → {df.index.max()}")

# Smooth data for cleaner plots (rolling window)
df_smooth = df.rolling(window=10, center=True, min_periods=1).mean()

# ── Helper: add CRISPR site annotations ──────────────────────────────────────
def annotate_crispr_site(ax, ymin, ymax):
    """Add sgRNA site and cut site annotations to a plot"""
    # sgRNA region
    ax.axvspan(SGRNA_START, SGRNA_END,
               alpha=0.15, color="red", label="sgRNA site")
    # Cut site
    ax.axvline(CUT_SITE, color="red", linestyle="--",
               linewidth=1.2, alpha=0.7, label="Cut site")
    ax.text(CUT_SITE, ymax * 0.97, "✂ Cut",
            ha="center", va="top", fontsize=8,
            color="red", fontweight="bold")
    # Zone boundaries
    for pos, label in [
        (SGRNA_START - 100, "←100bp"),
        (SGRNA_END  + 100, "100bp→"),
    ]:
        ax.axvline(pos, color="gray", linestyle=":",
                   linewidth=0.8, alpha=0.5)

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Coverage profile per sample (all samples overlaid)
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)

groups_ordered = ["Control", "Only_MNP", "Plasmid_Ko", "RNP_Cas"]

for ax, group in zip(axes, groups_ordered):
    group_samples = [s for s, g in sample_groups.items() if g == group]
    color = GROUP_COLORS[group]

    for sample in group_samples:
        short = (sample.replace("_L002", "")
                       .replace(f"_{sample.split('_S')[1].split('_')[0]}", "")
                 if "_S" in sample else sample)
        ax.plot(
            df_smooth.index,
            df_smooth[sample],
            color=color, linewidth=1,
            alpha=0.6, label=short
        )

    # Group mean
    group_mean = df_smooth[group_samples].mean(axis=1)
    ax.plot(df_smooth.index, group_mean,
            color=color, linewidth=2.5,
            alpha=1.0, label=f"{group} mean", zorder=5)

    ymax = df_smooth[group_samples].max().max()
    annotate_crispr_site(ax, 0, ymax)

    ax.set_ylabel("Depth (×)", fontsize=10)
    ax.set_title(f"{group}", fontsize=11, fontweight="bold",
                 color=color)
    ax.legend(fontsize=7, loc="upper left",
              framealpha=0.7, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")

axes[-1].set_xlabel("Genomic Position (NC_024333.1)", fontsize=10)
fig.suptitle(
    "Sequencing Depth Profile Along CRISPR-Cas9 Target Region (bdnf)\n"
    f"NC_024333.1 | sgRNA: TGAGAGACGCCCCGGGCATG | Cut site: ~{CUT_SITE:,}",
    fontsize=12, fontweight="bold", y=1.01
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot1_depth_profile_per_group.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 1 saved: depth profile per group")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Mean ± SEM per group (key comparison plot)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 6))

for group in groups_ordered:
    group_samples = [s for s, g in sample_groups.items() if g == group]
    color = GROUP_COLORS[group]

    group_data = df_smooth[group_samples]
    mean = group_data.mean(axis=1)
    sem  = group_data.sem(axis=1)

    ax.plot(df_smooth.index, mean,
            color=color, linewidth=2,
            label=group, zorder=3)
    ax.fill_between(
        df_smooth.index,
        mean - sem, mean + sem,
        color=color, alpha=0.2
    )

ymax = df_smooth.max().max()
annotate_crispr_site(ax, 0, ymax)

ax.set_xlabel("Genomic Position (NC_024333.1)", fontsize=11)
ax.set_ylabel("Mean Depth (×) ± SEM", fontsize=11)
ax.set_title(
    "Mean Sequencing Depth by Experimental Group\n"
    "Around CRISPR-Cas9 Cut Site (bdnf) — 10bp Rolling Average",
    fontsize=12, fontweight="bold", pad=15
)
ax.legend(fontsize=10, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")
ax.grid(True, axis="y", linestyle="--", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot2_mean_depth_by_group.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 2 saved: mean depth by group")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Heatmap: position × sample
# ═════════════════════════════════════════════════════════════════════════════
# Downsample positions for readability (every 5 positions)
df_ds = df_smooth.iloc[::5, :]

# Sort samples by group
sorted_samples = []
for group in groups_ordered:
    sorted_samples += [s for s, g in sample_groups.items() if g == group]

df_heatmap = df_ds[sorted_samples].T

fig, ax = plt.subplots(figsize=(16, 8))

im = ax.imshow(
    df_heatmap.values,
    cmap="YlOrRd",
    aspect="auto",
    interpolation="nearest"
)

# X axis — positions (sampled)
n_xticks = 10
xtick_idx = np.linspace(0, len(df_ds) - 1, n_xticks, dtype=int)
ax.set_xticks(xtick_idx)
ax.set_xticklabels(
    [f"{df_ds.index[i]:,}" for i in xtick_idx],
    rotation=45, ha="right", fontsize=8
)

# Y axis — samples
short_names = [
    (s.replace("_L002", "")
      .replace(f"_S{s.split('_S')[1].split('_')[0]}", "")
     if "_S" in s else s)
    for s in sorted_samples
]
ax.set_yticks(range(len(sorted_samples)))
ax.set_yticklabels(short_names, fontsize=8)

# Mark sgRNA site on x axis
pos_array = df_ds.index.values
sgrna_start_idx = np.searchsorted(pos_array, SGRNA_START)
sgrna_end_idx   = np.searchsorted(pos_array, SGRNA_END)
cut_idx         = np.searchsorted(pos_array, CUT_SITE)

ax.axvline(cut_idx, color="white", linewidth=2,
           linestyle="--", alpha=0.8, label="Cut site")
ax.axvspan(sgrna_start_idx, sgrna_end_idx,
           alpha=0.3, color="white")

ax.text(cut_idx, -1, "✂ Cut", ha="center", va="bottom",
        fontsize=8, color="red", fontweight="bold",
        transform=ax.get_xaxis_transform())

# Group separators
group_boundaries = []
count = 0
for group in groups_ordered:
    n = sum(1 for g in sample_groups.values() if g == group)
    if count > 0:
        group_boundaries.append(count - 0.5)
    count += n

for boundary in group_boundaries:
    ax.axhline(boundary, color="white", linewidth=2)

# Group labels on right
count = 0
for group in groups_ordered:
    n = sum(1 for g in sample_groups.values() if g == group)
    ax.text(
        len(df_ds) + 1, count + n / 2 - 0.5,
        group, va="center", fontsize=8,
        color=GROUP_COLORS[group], fontweight="bold",
        transform=ax.transData
    )
    count += n

cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.12)
cbar.set_label("Depth (×)", fontsize=9)

ax.set_xlabel("Genomic Position (NC_024333.1)", fontsize=11)
ax.set_ylabel("Sample", fontsize=11)
ax.set_title(
    "Depth Heatmap Across CRISPR-Cas9 Target Region (bdnf)\n"
    "Samples sorted by experimental group",
    fontsize=12, fontweight="bold", pad=15
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot3_heatmap_depth_per_position.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 3 saved: heatmap depth per position")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Normalized depth (Z-score per sample)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(14, 6))

# Normalize each sample by its own mean (ratio to mean)
df_norm = df_smooth.div(df_smooth.mean(axis=0), axis=1)

for group in groups_ordered:
    group_samples = [s for s, g in sample_groups.items() if g == group]
    color = GROUP_COLORS[group]

    norm_mean = df_norm[group_samples].mean(axis=1)
    norm_sem  = df_norm[group_samples].sem(axis=1)

    ax.plot(df_norm.index, norm_mean,
            color=color, linewidth=2,
            label=group, zorder=3)
    ax.fill_between(
        df_norm.index,
        norm_mean - norm_sem,
        norm_mean + norm_sem,
        color=color, alpha=0.2
    )

ax.axhline(1.0, color="black", linestyle="--",
           linewidth=1, alpha=0.5, label="Mean baseline")

annotate_crispr_site(ax, 0, df_norm.max().max())

ax.set_xlabel("Genomic Position (NC_024333.1)", fontsize=11)
ax.set_ylabel("Normalized Depth\n(ratio to sample mean)", fontsize=11)
ax.set_title(
    "Normalized Depth Profile by Experimental Group\n"
    "Values < 1.0 indicate reduced coverage relative to sample mean",
    fontsize=12, fontweight="bold", pad=15
)
ax.legend(fontsize=10, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")
ax.grid(True, axis="y", linestyle="--", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot4_normalized_depth_profile.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 4 saved: normalized depth profile")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Depth ratio: sgRNA window vs flanking regions per sample
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(12, 5))

results = []
for sample in sorted_samples:
    group = sample_groups[sample]

    sgrna_mask   = (df.index >= SGRNA_START) & (df.index <= SGRNA_END)
    flank_mask   = ((df.index >= SGRNA_START - 200) &
                    (df.index < SGRNA_START)) | \
                   ((df.index > SGRNA_END) &
                    (df.index <= SGRNA_END + 200))

    sgrna_mean = df.loc[sgrna_mask, sample].mean()
    flank_mean = df.loc[flank_mask, sample].mean()
    ratio      = sgrna_mean / flank_mean if flank_mean > 0 else np.nan

    short = (sample.replace("_L002", "")
                   .replace(f"_S{sample.split('_S')[1].split('_')[0]}", "")
             if "_S" in sample else sample)

    results.append({
        "sample":     short,
        "group":      group,
        "sgrna_depth": sgrna_mean,
        "flank_depth": flank_mean,
        "ratio":       ratio
    })

results_df = pd.DataFrame(results)

colors = [GROUP_COLORS[g] for g in results_df["group"]]
bars = ax.bar(
    range(len(results_df)),
    results_df["ratio"],
    color=colors, edgecolor="white",
    linewidth=0.5, alpha=0.85
)

# Reference line at 1.0
ax.axhline(1.0, color="black", linestyle="--",
           linewidth=1.5, alpha=0.7, label="No change (ratio = 1.0)")
ax.axhline(0.8, color="orange", linestyle=":",
           linewidth=1, alpha=0.7, label="20% reduction threshold")

# Value labels
for bar, val in zip(bars, results_df["ratio"]):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.005,
        f"{val:.3f}",
        ha="center", va="bottom",
        fontsize=7, fontweight="bold"
    )

ax.set_xticks(range(len(results_df)))
ax.set_xticklabels(results_df["sample"],
                   rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Depth Ratio\n(sgRNA site / ±200bp flanking)", fontsize=11)
ax.set_title(
    "Depth Ratio at sgRNA Site vs Flanking Regions\n"
    "Ratio < 1.0 may indicate reduced coverage at CRISPR cut site",
    fontsize=12, fontweight="bold", pad=15
)
ax.legend(fontsize=9, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")

group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
ax.legend(handles=group_patches + [
    mpatches.Patch(color="none", label=""),
    plt.Line2D([0], [0], color="black", linestyle="--",
               label="No change (1.0)"),
    plt.Line2D([0], [0], color="orange", linestyle=":",
               label="20% reduction")
], fontsize=8, framealpha=0.9)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot5_depth_ratio_sgrna_vs_flanking.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 5 saved: depth ratio sgRNA vs flanking")

# Save ratio table
results_df.to_csv(f"{OUTPUT_DIR}/depth_ratio_per_sample.csv", index=False)
print(f"\n=== Depth Ratio Summary ===")
print(results_df[["sample", "group", "sgrna_depth",
                   "flank_depth", "ratio"]].to_string(index=False))
print(f"\n✅ All outputs saved to: {OUTPUT_DIR}/")
print("   plot1_depth_profile_per_group.png")
print("   plot2_mean_depth_by_group.png")
print("   plot3_heatmap_depth_per_position.png")
print("   plot4_normalized_depth_profile.png")
print("   plot5_depth_ratio_sgrna_vs_flanking.png")
print("   depth_ratio_per_sample.csv")