#!/usr/bin/env python3
"""
Depth by zone visualization script for CRISPR-Cas9 experiment
Generates plots from depth_by_zone.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import os

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_CSV  = "depth_by_zone.csv"
OUTPUT_DIR = "depth_zone_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Color palette per group ───────────────────────────────────────────────────
GROUP_COLORS = {
    "Control":    "#2196F3",   # blue
    "Only_MNP":   "#4CAF50",   # green
    "Plasmid_Ko": "#FF9800",   # orange
    "RNP_Cas":    "#F44336",   # red
}

# Zone order and labels
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

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_CSV)
df["Zone"] = pd.Categorical(df["Zone"], categories=ZONE_ORDER, ordered=True)
df = df.sort_values(["Sample", "Zone"])

# Shorten sample names
df["Sample_short"] = (df["Sample"]
                      .str.replace("_L002", "")
                      .str.replace("_S\d+", "", regex=True))

print(f"Loaded {len(df)} rows")
print(f"Samples: {df['Sample'].nunique()}")
print(f"Groups:  {df['Group'].unique()}")
print(f"Zones:   {df['Zone'].unique()}")

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 1 — Mean Depth per Zone per Sample (grouped bar chart) — CORREGIDO
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(16, 7))

samples   = df["Sample_short"].unique()
n_samples = len(samples)
n_zones   = len(ZONE_ORDER)
bar_width = 0.15
x         = np.arange(n_samples)

zone_colors = ["#90CAF9", "#1565C0", "#F44336", "#E65100", "#A5D6A7"]

for i, zone in enumerate(ZONE_ORDER):
    zone_data = []
    for sample in samples:
        val = df[(df["Sample_short"] == sample) &
                 (df["Zone"] == zone)]["Mean_Depth"]
        zone_data.append(val.values[0] if len(val) > 0 else 0)

    offset = (i - n_zones / 2 + 0.5) * bar_width
    ax.bar(
        x + offset, zone_data,
        width=bar_width,
        color=zone_colors[i],
        label=ZONE_LABELS[zone].replace("\n", " "),
        edgecolor="white", linewidth=0.4, alpha=0.9
    )

# Indicadores de grupo debajo del eje x
for i, sample in enumerate(samples):
    group = df[df["Sample_short"] == sample]["Group"].iloc[0]
    ax.add_patch(mpatches.Rectangle(
        (i - 0.45, -8), 0.9, 4,
        color=GROUP_COLORS[group],
        clip_on=False, zorder=5
    ))

ax.set_xticks(x)
ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Mean Depth (×)", fontsize=11)
ax.set_title(
    "Mean Sequencing Depth by Zone Around CRISPR-Cas9 Cut Site (bdnf)\n"
    "sgRNA: TGAGAGACGCCCCGGGCATG | NC_024333.1:15,922,039–15,922,058",
    fontsize=12, fontweight="bold", pad=15
)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")
ax.set_ylim(-10, df["Mean_Depth"].max() * 1.15)

# ✅ Leyenda de zonas (colores de barras) — esquina superior derecha
zone_patches = [mpatches.Patch(color=zone_colors[i],
                               label=ZONE_LABELS[z].replace("\n", " "))
                for i, z in enumerate(ZONE_ORDER)]
legend_zones = ax.legend(
    handles=zone_patches,
    loc="upper right",
    title="Zone",
    fontsize=8,
    title_fontsize=9,
    framealpha=0.9
)

# ✅ Leyenda de grupos (colores debajo del eje) — esquina superior izquierda
group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
ax.add_artist(legend_zones)   # ← preserva la leyenda de zonas
ax.legend(
    handles=group_patches,
    loc="upper left",
    title="Experimental Group",
    fontsize=8,
    title_fontsize=9,
    framealpha=0.9
)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot1_depth_by_zone_per_sample.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 1 saved: depth by zone per sample")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Mean Depth by Zone and Group (line plot — key comparison)
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))

group_zone = df.groupby(["Group", "Zone"])["Mean_Depth"].agg(
    ["mean", "std", "sem"]
).reset_index()

x_pos = np.arange(len(ZONE_ORDER))

for group, color in GROUP_COLORS.items():
    gdata = group_zone[group_zone["Group"] == group].set_index("Zone")
    means = [gdata.loc[z, "mean"] if z in gdata.index else np.nan
             for z in ZONE_ORDER]
    sems  = [gdata.loc[z, "sem"]  if z in gdata.index else np.nan
             for z in ZONE_ORDER]

    ax.plot(x_pos, means, color=color, marker="o",
            linewidth=2.5, markersize=8, label=group, zorder=3)
    ax.fill_between(
        x_pos,
        [m - s for m, s in zip(means, sems)],
        [m + s for m, s in zip(means, sems)],
        color=color, alpha=0.15
    )
    ax.errorbar(x_pos, means, yerr=sems,
                fmt="none", color=color, capsize=4,
                linewidth=1.5, zorder=4)

# Highlight sgRNA site
ax.axvspan(1.5, 2.5, alpha=0.1, color="red", label="sgRNA site")
ax.axvline(2, color="red", linestyle="--", linewidth=1, alpha=0.5)
ax.text(2, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 70,
        "Cut site", ha="center", va="bottom",
        fontsize=8, color="red", alpha=0.7)

ax.set_xticks(x_pos)
ax.set_xticklabels([ZONE_LABELS[z] for z in ZONE_ORDER], fontsize=10)
ax.set_ylabel("Mean Depth (×) ± SEM", fontsize=11)
ax.set_title("Mean Depth by Zone and Experimental Group\n"
             "Around CRISPR-Cas9 Cut Site (bdnf)",
             fontsize=12, fontweight="bold", pad=15)
ax.legend(fontsize=9, framealpha=0.9)
ax.spines[["top", "right"]].set_visible(False)
ax.set_facecolor("#f9f9f9")
ax.grid(True, axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot2_depth_by_zone_and_group_lineplot.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 2 saved: depth by zone and group line plot")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Boxplot per Zone colored by Group
# ═════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, len(ZONE_ORDER), figsize=(16, 6), sharey=True)

for ax, zone in zip(axes, ZONE_ORDER):
    zone_df = df[df["Zone"] == zone]
    groups  = list(GROUP_COLORS.keys())

    group_data = [zone_df[zone_df["Group"] == g]["Mean_Depth"].values
                  for g in groups]
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
        whisker.set(color="#555555", linewidth=1)
    for cap in bp["caps"]:
        cap.set(color="#555555", linewidth=1)

    # Overlay points
    for i, (gdata, color) in enumerate(zip(group_data, colors), start=1):
        x_jitter = np.random.normal(i, 0.06, size=len(gdata))
        ax.scatter(x_jitter, gdata, color=color,
                   edgecolor="white", linewidth=0.5,
                   s=40, zorder=3, alpha=0.9)

    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([g.replace("_", "\n") for g in groups],
                       fontsize=7)
    ax.set_title(ZONE_LABELS[zone], fontsize=10, fontweight="bold",
                 color="red" if zone == "sgRNA_site" else "black")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")

    # Highlight sgRNA site
    if zone == "sgRNA_site":
        ax.set_facecolor("#FFF3F3")
        for spine in ax.spines.values():
            spine.set_edgecolor("#F44336")
            spine.set_linewidth(1.5)

axes[0].set_ylabel("Mean Depth (×)", fontsize=11)

fig.suptitle("Depth Distribution per Zone and Experimental Group\n"
             "Around CRISPR-Cas9 Cut Site (bdnf)",
             fontsize=12, fontweight="bold", y=1.02)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot3_boxplot_per_zone.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 3 saved: boxplot per zone")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Heatmap: Mean Depth por muestra y zona
# ═════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 10))

# Pivot: filas = muestras, columnas = zonas
pivot = df.pivot_table(
    index="Sample_short",
    columns="Zone",
    values="Mean_Depth"
)[ZONE_ORDER]

# Ordenar por grupo
sample_groups = df[["Sample_short", "Group"]].drop_duplicates()
sample_groups = sample_groups.sort_values("Group")
pivot = pivot.loc[sample_groups["Sample_short"]]

im = ax.imshow(pivot.values, cmap="YlOrRd", aspect="auto")

ax.set_xticks(range(len(ZONE_ORDER)))
ax.set_xticklabels([ZONE_LABELS[z] for z in ZONE_ORDER], fontsize=9)
ax.set_yticks(range(len(pivot)))
ax.set_yticklabels(pivot.index, fontsize=8)

# Valores en celdas
vmax = pivot.values.max()
for i in range(len(pivot)):
    for j in range(len(ZONE_ORDER)):
        val = pivot.iloc[i, j]
        text_color = "white" if val > vmax * 0.6 else "black"
        ax.text(j, i, f"{val:.1f}×",
                ha="center", va="center",
                fontsize=7, color=text_color,
                fontweight="bold")

# Colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.12)
cbar.set_label("Mean Depth (×)", fontsize=9)

# Indicadores de grupo a la derecha
for i, (_, row) in enumerate(sample_groups.iterrows()):
    ax.add_patch(mpatches.Rectangle(
        (len(ZONE_ORDER) - 0.5, i - 0.5), 0.4, 1,
        color=GROUP_COLORS[row["Group"]],
        clip_on=False, zorder=5
    ))

# Línea separadora entre grupos
prev_group = None
for i, (_, row) in enumerate(sample_groups.iterrows()):
    if prev_group and row["Group"] != prev_group:
        ax.axhline(i - 0.5, color="white", linewidth=2)
    prev_group = row["Group"]

# Highlight sgRNA site column
ax.add_patch(mpatches.Rectangle(
    (1.5, -0.5), 1, len(pivot),
    fill=False, edgecolor="red",
    linewidth=2, clip_on=False
))

ax.set_title("Mean Depth Heatmap by Zone\nCRISPR-Cas9 Target Site (bdnf)",
             fontsize=11, fontweight="bold", pad=15)

# Group legend
group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
ax.legend(handles=group_patches, bbox_to_anchor=(1.25, 1),
          loc="upper left", fontsize=8, title="Group",
          title_fontsize=9)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot4_heatmap_depth_by_zone.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 4 saved: heatmap depth by zone")

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 5 — sgRNA site depth vs upstream/downstream (scatter)
# ═════════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

sgrna_depth = df[df["Zone"] == "sgRNA_site"].set_index("Sample_short")["Mean_Depth"]
up500_depth  = df[df["Zone"] == "upstream_500bp"].set_index("Sample_short")["Mean_Depth"]
dn500_depth  = df[df["Zone"] == "downstream_500bp"].set_index("Sample_short")["Mean_Depth"]

for sample in df["Sample_short"].unique():
    group = df[df["Sample_short"] == sample]["Group"].iloc[0]
    color = GROUP_COLORS[group]

    if sample in sgrna_depth.index and sample in up500_depth.index:
        ax1.scatter(up500_depth[sample], sgrna_depth[sample],
                    color=color, s=80, edgecolor="white",
                    linewidth=0.8, alpha=0.9, zorder=3)
        ax1.annotate(sample, (up500_depth[sample], sgrna_depth[sample]),
                     textcoords="offset points", xytext=(4, 3),
                     fontsize=6, color="#555555")

    if sample in sgrna_depth.index and sample in dn500_depth.index:
        ax2.scatter(dn500_depth[sample], sgrna_depth[sample],
                    color=color, s=80, edgecolor="white",
                    linewidth=0.8, alpha=0.9, zorder=3)
        ax2.annotate(sample, (dn500_depth[sample], sgrna_depth[sample]),
                     textcoords="offset points", xytext=(4, 3),
                     fontsize=6, color="#555555")

# Diagonal reference line
for ax in [ax1, ax2]:
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "k--", linewidth=1, alpha=0.4, label="x=y")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_facecolor("#f9f9f9")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_ylabel("sgRNA Site Mean Depth (×)", fontsize=10)

ax1.set_xlabel("Upstream 500bp Mean Depth (×)", fontsize=10)
ax1.set_title("sgRNA Site vs Upstream 500bp", fontsize=11, fontweight="bold")

ax2.set_xlabel("Downstream 500bp Mean Depth (×)", fontsize=10)
ax2.set_title("sgRNA Site vs Downstream 500bp", fontsize=11, fontweight="bold")

# Shared legend
group_patches = [mpatches.Patch(color=c, label=g)
                 for g, c in GROUP_COLORS.items()]
fig.legend(handles=group_patches, loc="lower center",
           ncol=4, fontsize=9, framealpha=0.9,
           bbox_to_anchor=(0.5, -0.05))

fig.suptitle("sgRNA Site Depth vs Flanking Regions\n"
             "Points below diagonal = reduced coverage at cut site",
             fontsize=11, fontweight="bold")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plot5_sgrna_vs_flanking_scatter.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("✅ Plot 5 saved: sgRNA vs flanking scatter")

# ═════════════════════════════════════════════════════════════════════════════
# Summary statistics
# ═════════════════════════════════════════════════════════════════════════════
print("\n=== Mean Depth by Group and Zone ===")
summary = df.groupby(["Group", "Zone"])["Mean_Depth"].agg(
    ["mean", "std", "min", "max"]
).round(2)
print(summary.to_string())
summary.to_csv(f"{OUTPUT_DIR}/depth_by_zone_summary.csv")

print(f"\n✅ All plots saved to: {OUTPUT_DIR}/")
print("   plot1_depth_by_zone_per_sample.png")
print("   plot2_depth_by_zone_and_group_lineplot.png")
print("   plot3_boxplot_per_zone.png")
print("   plot4_heatmap_depth_by_zone.png")
print("   plot5_sgrna_vs_flanking_scatter.png")
print("   depth_by_zone_summary.csv")