#!/usr/bin/env python3
"""
Coverage by zone visualization script for CRISPR-Cas9 experiment
Generates plots from coverage_by_zone_all_samples.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import os

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_CSV  = "coverage_by_zone_all_samples.csv"
OUTPUT_DIR = "coverage_zone_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Color palette ─────────────────────────────────────────────────────────────
GROUP_COLORS = {
    "Control":    "#2196F3",
    "Only_MNP":   "#4CAF50",
    "Plasmid_Ko": "#FF9800",
    "RNP_Cas":    "#F44336",
}

ZONE_ORDER = [
    "upstream_500bp",
    "upstream_100bp",
    "sgRNA_site",
    "downstream_100bp",
    "downstream_500bp"
]

ZONE_LABELS = {
    "upstream_500bp":    "Upstream\n500bp",
    "upstream_100bp":    "Upstream\n100bp",
    "sgRNA_site":        "sgRNA\nSite",
    "downstream_100bp":  "Downstream\n100bp",
    "downstream_500bp":  "Downstream\n500bp"
}

ZONE_COLORS = {
    "upstream_500bp":    "#90CAF9",
    "upstream_100bp":    "#1565C0",
    "sgRNA_site":        "#F44336",
    "downstream_100bp":  "#E65100",
    "downstream_500bp":  "#A5D6A7"
}

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)

# Exclude full_region rows for zone-specific plots
df_zones = df[df["Zone"] != "full_region"].copy()
df_full  = df[df["Zone"] == "full_region"].copy()

# Ordered zone category
df_zones["Zone"] = pd.Categorical(
    df_zones["Zone"], categories=ZONE_ORDER, ordered=True
)
df_zones = df_zones.sort_values(["Sample", "Zone"])

# Shorten sample names
df_zones["Sample_short"] = (df_zones["Sample"]
    .str.replace("_L002", "")
    .str.replace("_S\d+", "", regex=True))
df_full["Sample_short"] = (df_full["Sample"]
    .str.replace("_L002", "")
    .str.replace("_S\d+", "", regex=True))

print(f"Loaded {len(df_zones)} zone rows across "
      f"{df_zones['Sample'].nunique()} samples")

groups_ordered = ["Control", "Only_MNP", "Plasmid_Ko", "RNP_Cas"]

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Mean Depth per Zone and Group (line plot — main comparison)
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

metrics = [
    ("Mean_Depth",    "Mean Depth (×)",          "Mean Depth by Zone and Group"),
    ("Coverage_Pct",  "Coverage (%)",             "Coverage % by Zone and Group"),
]

for ax, (metric, ylabel, title) in zip(axes, metrics):
    for group in groups_ordered:
        gdf   = df_zones[df_zones["Group"] == group]
        color = GROUP_COLORS[group]

        zone_stats = gdf.groupby("Zone")[metric].agg(["mean", "sem"])
        means = [zone_stats.loc[z, "mean"] if z in zone_stats.index else np.nan
                 for z in ZONE_ORDER]
        sems  = [zone_stats.loc[z, "sem"]  if z in zone_stats.index else np.nan
                 for z in ZONE_ORDER]

        x = np.arange(len(ZONE_ORDER))
        ax.plot(x, means, color=color, marker="o",
                linewidth=2.5, markersize=8,
                label=group, zorder=3)
        ax.fill_between(
            x,
            [m - s if not np.isnan(m) else np.nan
             for m, s in zip(means, sems)],
            [m + s if not np.isnan(m) else np.nan
             for m, s in zip(means, sems)],
            color=color, alpha=0.15
        )
        ax.errorbar(x, means, yerr=sems,
                    fmt="none", color=color,
                    capsize=4, linewidth=1.5, zorder=4)

    # Highlight sgRNA site
    ax.axvspan(1.5, 2.5, alpha=0.08, color="red")
    ax.axvline(2, color="red", linestyle="--",
               linewidth=1.2, alpha=0.6)
    ax.text(2, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 60,
            "✂ Cut site", ha="center", va="top",
            fontsize=8, color="red", style="italic")

    ax.set_xticks(range(len(ZONE_ORDER)))
    ax.set_xticklabels([ZONE_LABELS[z] for z in ZONE_ORDER], fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)

fig.suptitle(
    "Coverage Metrics by Zone Around CRISPR-Cas9 Cut Site (bdnf)\n"
    "sgRNA: TGAGAGACGCCCCGGGCATG | NC_024333.1:15,922,039–15,922,058",
    fontsize=12, fontweight="bold", y=1.02
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot1_coverage_by_zone_lineplot.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 1 saved: coverage by zone line plot")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Boxplot per zone for Mean Depth and Coverage %
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, len(ZONE_ORDER),
                         figsize=(16, 10), sharey="row")

for col, zone in enumerate(ZONE_ORDER):
    zone_df = df_zones[df_zones["Zone"] == zone]

    for row, (metric, ylabel) in enumerate([
        ("Mean_Depth",   "Mean Depth (×)"),
        ("Coverage_Pct", "Coverage (%)")
    ]):
        ax = axes[row][col]
        group_data = [zone_df[zone_df["Group"] == g][metric].values
                      for g in groups_ordered]
        colors     = [GROUP_COLORS[g] for g in groups_ordered]

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
            whisker.set(color="#555555", linewidth=1)
        for cap in bp["caps"]:
            cap.set(color="#555555", linewidth=1)

        # Overlay individual points
        for i, (gdata, color) in enumerate(zip(group_data, colors), start=1):
            jitter = np.random.normal(i, 0.06, size=len(gdata))
            ax.scatter(jitter, gdata, color=color,
                       edgecolor="white", linewidth=0.5,
                       s=35, zorder=3, alpha=0.9)

        ax.set_xticks(range(1, len(groups_ordered) + 1))
        ax.set_xticklabels(
            [g.replace("_", "\n") for g in groups_ordered],
            fontsize=6
        )

        # Title only on top row
        if row == 0:
            title_color = "red" if zone == "sgRNA_site" else "black"
            ax.set_title(ZONE_LABELS[zone], fontsize=10,
                         fontweight="bold", color=title_color)

        # Y label only on leftmost column
        if col == 0:
            ax.set_ylabel(ylabel, fontsize=10)

        # Highlight sgRNA site
        if zone == "sgRNA_site":
            ax.set_facecolor("#FFF3F3")
            for spine in ax.spines.values():
                spine.set_edgecolor("#F44336")
                spine.set_linewidth(1.5)
        else:
            ax.set_facecolor("#f9f9f9")

        ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    "Coverage Distribution per Zone and Experimental Group\n"
    "Around CRISPR-Cas9 Cut Site (bdnf)",
    fontsize=12, fontweight="bold"
)

# Group legend
group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
fig.legend(handles=group_patches, loc="lower center",
           ncol=4, fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot2_boxplot_per_zone.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 2 saved: boxplot per zone")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Heatmaps: Mean Depth and Coverage % (sample × zone)
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 9))

# Sort samples by group
sample_order = []
for group in groups_ordered:
    samples = df_zones[df_zones["Group"] == group]["Sample_short"].unique()
    sample_order += list(samples)

for ax, (metric, title, fmt) in zip(axes, [
    ("Mean_Depth",   "Mean Depth (×)",  "{:.1f}×"),
    ("Coverage_Pct", "Coverage (%)",    "{:.1f}%"),
]):
    pivot = df_zones.pivot_table(
        index="Sample_short",
        columns="Zone",
        values=metric
    )[ZONE_ORDER]

    # Reorder rows by group
    pivot = pivot.reindex(
        [s for s in sample_order if s in pivot.index]
    )

    im = ax.imshow(pivot.values, cmap="YlOrRd",
                   aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(ZONE_ORDER)))
    ax.set_xticklabels([ZONE_LABELS[z] for z in ZONE_ORDER],
                       fontsize=9)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=8)

    # Cell values
    vmax = np.nanmax(pivot.values)
    for i in range(len(pivot)):
        for j in range(len(ZONE_ORDER)):
            val = pivot.iloc[i, j]
            if not np.isnan(val):
                text_color = "white" if val > vmax * 0.65 else "black"
                ax.text(j, i, fmt.format(val),
                        ha="center", va="center",
                        fontsize=7, color=text_color,
                        fontweight="bold")

    # Group color indicators
    sample_group_map = df_zones[["Sample_short", "Group"]].drop_duplicates()
    sample_group_map = sample_group_map.set_index("Sample_short")

    for i, sample in enumerate(pivot.index):
        group = sample_group_map.loc[sample, "Group"]
        ax.add_patch(mpatches.Rectangle(
            (len(ZONE_ORDER) - 0.5, i - 0.5), 0.3, 1,
            color=GROUP_COLORS[group],
            clip_on=False, zorder=5
        ))

    # Group separators
    prev_group = None
    for i, sample in enumerate(pivot.index):
        group = sample_group_map.loc[sample, "Group"]
        if prev_group and group != prev_group:
            ax.axhline(i - 0.5, color="white", linewidth=2)
        prev_group = group

    # Highlight sgRNA site column
    sgrna_idx = ZONE_ORDER.index("sgRNA_site")
    ax.add_patch(mpatches.Rectangle(
        (sgrna_idx - 0.5, -0.5), 1, len(pivot),
        fill=False, edgecolor="red",
        linewidth=2, clip_on=False
    ))

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.15)
    cbar.set_label(title, fontsize=9)
    ax.set_title(f"Heatmap: {title}", fontsize=11,
                 fontweight="bold", pad=10)

fig.suptitle(
    "Coverage Heatmaps by Zone\nCRISPR-Cas9 Target Site (bdnf)",
    fontsize=12, fontweight="bold"
)

group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
fig.legend(handles=group_patches, loc="lower center",
           ncol=4, fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot3_heatmap_coverage_by_zone.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 3 saved: heatmap coverage by zone")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Ratio: sgRNA site vs flanking zones per sample
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

for ax, (metric, ylabel, title) in zip(axes, [
    ("Mean_Depth",   "Depth Ratio\n(sgRNA / flanking)",   "Mean Depth Ratio"),
    ("Coverage_Pct", "Coverage Ratio\n(sgRNA / flanking)", "Coverage % Ratio"),
]):
    ratios = []
    for sample in df_zones["Sample_short"].unique():
        sdf   = df_zones[df_zones["Sample_short"] == sample]
        group = sdf["Group"].iloc[0]

        sgrna_val = sdf[sdf["Zone"] == "sgRNA_site"][metric]
        flank_val = sdf[sdf["Zone"].isin([
            "upstream_100bp", "downstream_100bp"
        ])][metric]

        if len(sgrna_val) > 0 and len(flank_val) > 0:
            ratio = sgrna_val.values[0] / flank_val.mean()
            ratios.append({
                "sample": sample,
                "group":  group,
                "ratio":  ratio,
                "sgrna":  sgrna_val.values[0],
                "flank":  flank_val.mean()
            })

    ratios_df = pd.DataFrame(ratios)
    ratios_df = ratios_df.sort_values("group")
    colors = [GROUP_COLORS[g] for g in ratios_df["group"]]

    bars = ax.bar(
        range(len(ratios_df)),
        ratios_df["ratio"],
        color=colors, edgecolor="white",
        linewidth=0.5, alpha=0.85
    )

    # Reference lines
    ax.axhline(1.0, color="black", linestyle="--",
               linewidth=1.5, alpha=0.6, label="No change (1.0)")
    ax.axhline(0.8, color="orange", linestyle=":",
               linewidth=1.2, alpha=0.7, label="20% reduction")
    ax.axhline(1.2, color="purple", linestyle=":",
               linewidth=1.2, alpha=0.7, label="20% increase")

    # Value labels
    for bar, val in zip(bars, ratios_df["ratio"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.005,
            f"{val:.3f}",
            ha="center", va="bottom",
            fontsize=7, fontweight="bold"
        )

    ax.set_xticks(range(len(ratios_df)))
    ax.set_xticklabels(ratios_df["sample"],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")

    # Group color indicators
    for i, (_, row) in enumerate(ratios_df.iterrows()):
        ax.add_patch(mpatches.Rectangle(
            (i - 0.45, -0.06), 0.9, 0.025,
            color=GROUP_COLORS[row["group"]],
            clip_on=False, zorder=5
        ))

    # Legends
    group_patches = [mpatches.Patch(color=c, label=g)
                     for g, c in GROUP_COLORS.items()]
    line_handles = [
        plt.Line2D([0], [0], color="black", linestyle="--",
                   label="No change (1.0)"),
        plt.Line2D([0], [0], color="orange", linestyle=":",
                   label="20% reduction"),
        plt.Line2D([0], [0], color="purple", linestyle=":",
                   label="20% increase"),
    ]
    legend_zones = ax.legend(
        handles=line_handles,
        loc="upper right", fontsize=8, framealpha=0.9
    )
    ax.add_artist(legend_zones)
    ax.legend(handles=group_patches, loc="upper left",
              fontsize=8, framealpha=0.9)

fig.suptitle(
    "Coverage Ratio: sgRNA Site vs ±100bp Flanking Regions\n"
    "Ratio < 1.0 = reduced coverage at cut site | "
    "Ratio > 1.0 = increased coverage",
    fontsize=12, fontweight="bold"
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot4_coverage_ratio_sgrna_vs_flanking.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 4 saved: coverage ratio sgRNA vs flanking")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 5 — Summary stats table as figure
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(16, 5))

for ax, (metric, title) in zip(axes, [
    ("Mean_Depth",   "Mean Depth (×) by Group and Zone"),
    ("Coverage_Pct", "Coverage (%) by Group and Zone"),
]):
    summary = df_zones.groupby(["Group", "Zone"])[metric].agg(
        ["mean", "std"]
    ).round(2).reset_index()
    summary["Zone"] = pd.Categorical(
        summary["Zone"], categories=ZONE_ORDER, ordered=True
    )
    summary = summary.sort_values(["Group", "Zone"])

    x     = np.arange(len(ZONE_ORDER))
    n_grp = len(groups_ordered)
    width = 0.18

    for i, group in enumerate(groups_ordered):
        gdata  = summary[summary["Group"] == group].set_index("Zone")
        means  = [gdata.loc[z, "mean"] if z in gdata.index else 0
                  for z in ZONE_ORDER]
        stds   = [gdata.loc[z, "std"]  if z in gdata.index else 0
                  for z in ZONE_ORDER]
        offset = (i - n_grp / 2 + 0.5) * width

        ax.bar(x + offset, means, width=width,
               color=GROUP_COLORS[group],
               label=group, alpha=0.85,
               edgecolor="white", linewidth=0.4)
        ax.errorbar(x + offset, means, yerr=stds,
                    fmt="none", color="black",
                    capsize=3, linewidth=1)

    # Highlight sgRNA zone
    ax.axvspan(1.5, 2.5, alpha=0.07, color="red")
    ax.axvline(2, color="red", linestyle="--",
               linewidth=1, alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([ZONE_LABELS[z] for z in ZONE_ORDER], fontsize=10)
    ax.set_ylabel(f"Mean ± SD {title.split('(')[1].split(')')[0]}",
                  fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")
    ax.grid(True, axis="y", linestyle="--", alpha=0.3)

fig.suptitle(
    "Group Summary: Coverage Metrics by Zone\n"
    "Around CRISPR-Cas9 Cut Site (bdnf) — Mean ± SD",
    fontsize=12, fontweight="bold"
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot5_group_summary_by_zone.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 5 saved: group summary by zone")

# ── Summary CSV ───────────────────────────────────────────────────────────────
summary_out = df_zones.groupby(["Group", "Zone"]).agg(
    N=("Sample", "count"),
    Mean_Depth_mean=("Mean_Depth", "mean"),
    Mean_Depth_sd=("Mean_Depth", "std"),
    Coverage_Pct_mean=("Coverage_Pct", "mean"),
    Coverage_Pct_sd=("Coverage_Pct", "std"),
    NumReads_mean=("NumReads", "mean"),
    Mean_MapQ_mean=("Mean_MapQ", "mean")
).round(2)

summary_out.to_csv(f"{OUTPUT_DIR}/coverage_by_zone_summary.csv")

print(f"\n=== Summary by Group and Zone ===")
print(summary_out.to_string())
print(f"\n✅ All outputs saved to: {OUTPUT_DIR}/")
print("   plot1_coverage_by_zone_lineplot.png")
print("   plot2_boxplot_per_zone.png")
print("   plot3_heatmap_coverage_by_zone.png")
print("   plot4_coverage_ratio_sgrna_vs_flanking.png")
print("   plot5_group_summary_by_zone.png")
print("   coverage_by_zone_summary.csv")