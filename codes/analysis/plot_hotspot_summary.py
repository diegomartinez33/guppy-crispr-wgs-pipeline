#!/usr/bin/env python3
"""
Hotspot analysis summary plots.

Reads window_counts_annotated.csv and hotspots.bed produced by
hotspot_analysis.py and generates five publication-ready figures:

  hotspot_summary_01_variants_per_chrom.png  — total SNP+INDEL counts per LG
  hotspot_summary_02_snp_indel_breakdown.png — stacked SNP vs INDEL per LG
  hotspot_summary_03_hotspot_regions.png     — merged hotspot regions ranked by Z-score
  hotspot_summary_04_densest_zoom.png        — zoomed density plot for the densest LG
  hotspot_summary_05_density_grid.png        — grid of density plots for all active LGs

Usage:
    python codes/analysis/plot_hotspot_summary.py
    REF_VERSION=v2 python codes/analysis/plot_hotspot_summary.py   # new reference
"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
REF_VERSION = os.environ.get("REF_VERSION", "v1")
OUT_SUFFIX = "" if REF_VERSION == "v1" else f"_{REF_VERSION}"
HOT_DIR     = os.path.join(PROJECT_DIR, f"gatk/trimmomatic{OUT_SUFFIX}/hotspots")
OUT_DIR     = os.path.join(HOT_DIR, "summary_plots")
os.makedirs(OUT_DIR, exist_ok=True)

ANNOTATED_CSV  = os.path.join(HOT_DIR, "window_counts_annotated.csv")
HOTSPOTS_BED   = os.path.join(HOT_DIR, "hotspots.bed")

# ── Colour palette ─────────────────────────────────────────────────────────────
PALETTE = ["#4C8EBF", "#E8A838"]
HOTSPOT_RED  = "#D94040"
SNP_COLOR    = "#3A7AB5"
INDEL_COLOR  = "#E07B30"

# Genome-specific accession ranges (not a formula that generalizes across
# assemblies) - see hotspot_analysis.py for the full rationale.
CHROM_TO_LG = {
    "v1": {f"NC_0243{31 + i:02d}.1": f"LG{i + 1}" for i in range(23)},
    "v2": {f"NC_0888{30 + i:02d}.1": f"LG{i + 1}" for i in range(23)},
}[REF_VERSION]
MAIN_CHROMS = list(CHROM_TO_LG) + (["NW_007615013.1"] if REF_VERSION == "v1" else [])

CHROM_COLORS = {}
for _i, _c in enumerate(MAIN_CHROMS):
    CHROM_COLORS[_c] = PALETTE[_i % 2]


def chrom_label(chrom):
    return CHROM_TO_LG.get(chrom, "Un")


# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
df = pd.read_csv(ANNOTATED_CSV)
df_main = df[df["chrom"].isin(MAIN_CHROMS)].copy()

by_chrom = (
    df_main.groupby("chrom")
    .agg(
        total_variants=("total_count", "sum"),
        snp_count=("snp_count", "sum"),
        indel_count=("indel_count", "sum"),
        hotspot_windows=("is_hotspot", "sum"),
        max_z=("z_score", "max"),
    )
    .reindex(MAIN_CHROMS)
    .reset_index()
    .rename(columns={"index": "chrom"})
    .fillna(0)
)
by_chrom["label"] = by_chrom["chrom"].apply(chrom_label)

hotspots = pd.read_csv(HOTSPOTS_BED, sep="\t")
hotspots["label"] = hotspots["chrom"].apply(chrom_label)
hotspots["width_kb"] = (hotspots["end"] - hotspots["start"]) / 1000
hotspots["mid_mb"] = (hotspots["start"] + hotspots["end"]) / 2 / 1e6

# Chromosomes that actually have variant data
active_chroms = [
    c for c in MAIN_CHROMS
    if by_chrom.loc[by_chrom["chrom"] == c, "total_variants"].values[0] > 0
]
n_active = len(active_chroms)

# Chromosome with highest single-window Z-score (for the zoom plot)
densest_chrom = hotspots.loc[hotspots["max_z"].idxmax(), "chrom"]
densest_label = chrom_label(densest_chrom)

print(f"  {len(df_main):,} windows  |  {len(hotspots)} hotspot regions  |  {n_active} active LGs")
print(f"  Densest LG: {densest_label} ({densest_chrom}), max Z = {hotspots['max_z'].max():.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — Total variants per chromosome (bar chart)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 1: variants per chromosome...")

fig, ax = plt.subplots(figsize=(14, 5))
x = np.arange(len(by_chrom))
colors = [CHROM_COLORS.get(c, "#888888") for c in by_chrom["chrom"]]
bars = ax.bar(x, by_chrom["total_variants"] / 1e6, color=colors, edgecolor="white", linewidth=0.4)

for bar, val in zip(bars, by_chrom["total_variants"]):
    if val > 0:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.04,
                f"{val/1e6:.1f}", ha="center", va="bottom", fontsize=7)

ax.set_xticks(x)
ax.set_xticklabels(by_chrom["label"], fontsize=8, rotation=45, ha="right")
ax.set_ylabel("Total PASS variants (millions)", fontsize=10)
ax.set_xlabel("Linkage group", fontsize=10)
ax.set_title(
    "Total PASS variant count per linkage group\n"
    "Colombian guppy WGS — SNPs + INDELs combined",
    fontsize=11, fontweight="bold"
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_01_variants_per_chrom.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 1 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — Stacked SNP vs INDEL per chromosome
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 2: SNP vs INDEL breakdown...")

fig, ax = plt.subplots(figsize=(14, 5))
x = np.arange(len(by_chrom))
ax.bar(x, by_chrom["snp_count"] / 1e6, color=SNP_COLOR, label="SNPs", edgecolor="white", linewidth=0.4)
ax.bar(x, by_chrom["indel_count"] / 1e6, bottom=by_chrom["snp_count"] / 1e6,
       color=INDEL_COLOR, label="INDELs", edgecolor="white", linewidth=0.4)

ax.set_xticks(x)
ax.set_xticklabels(by_chrom["label"], fontsize=8, rotation=45, ha="right")
ax.set_ylabel("PASS variants (millions)", fontsize=10)
ax.set_xlabel("Linkage group", fontsize=10)
ax.set_title(
    "SNP vs INDEL breakdown per linkage group\n"
    "Colombian guppy WGS — PASS variants only",
    fontsize=11, fontweight="bold"
)
ax.legend(fontsize=9, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_02_snp_indel_breakdown.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 2 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 — Hotspot regions ranked by Z-score (lollipop chart)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 3: hotspot regions ranked by Z-score...")

top = hotspots.sort_values("max_z", ascending=False).head(40).reset_index(drop=True)
top["region_label"] = (
    top["label"] + ":" +
    (top["start"] / 1e6).round(2).astype(str) + "–" +
    (top["end"] / 1e6).round(2).astype(str) + " Mb"
)
colors_top = [CHROM_COLORS.get(c, "#888888") for c in top["chrom"]]

# Build legend from unique chroms in top 40
top_chroms_ordered = list(dict.fromkeys(top["chrom"].tolist()))
legend_patches = [
    mpatches.Patch(color=CHROM_COLORS.get(c, "#888888"), label=chrom_label(c))
    for c in top_chroms_ordered
]

fig, ax = plt.subplots(figsize=(10, 12))
y = np.arange(len(top))[::-1]
ax.hlines(y, 0, top["max_z"], colors=colors_top, linewidth=1.5, alpha=0.7)
ax.scatter(top["max_z"], y, color=colors_top, s=60, zorder=3)
ax.axvline(4.0, color="red", linewidth=1, linestyle="--", alpha=0.6, label="Z = 4.0")

ax.set_yticks(y)
ax.set_yticklabels(top["region_label"], fontsize=7.5)
ax.set_xlabel("Max Z-score (variant density)", fontsize=10)
ax.set_title(
    "Top 40 hotspot regions ranked by Z-score\n"
    "Merged windows (bedtools merge -d 2000 bp)",
    fontsize=11, fontweight="bold"
)
ax.legend(handles=legend_patches, fontsize=8, frameon=False, loc="lower right",
          ncol=2 if len(legend_patches) > 8 else 1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_03_hotspot_regions.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 3 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — Densest LG zoom (full chromosome + zoomed hotspot cluster)
# ══════════════════════════════════════════════════════════════════════════════
print(f"Plot 4: {densest_label} zoom ({densest_chrom})...")

lg_data = df_main[df_main["chrom"] == densest_chrom].copy().sort_values("start")
lg_hot  = hotspots[hotspots["chrom"] == densest_chrom]
mu      = lg_data["density"].mean()
sd      = lg_data["density"].std(ddof=1)
mids    = (lg_data["start"] + lg_data["end"]) / 2 / 1e6
smoothed = lg_data["density"].rolling(5, center=True, min_periods=1).mean()
max_z_val = lg_hot["max_z"].max() if len(lg_hot) > 0 else 0

# Find the cluster: region around the hotspot with the highest Z-score
if len(lg_hot) > 0:
    top_hot = lg_hot.loc[lg_hot["max_z"].idxmax()]
    cluster_center = (top_hot["start"] + top_hot["end"]) / 2
    cluster_half   = max(5_000_000, (top_hot["end"] - top_hot["start"]) * 5)
    zoom_start = max(0, cluster_center - cluster_half)
    zoom_end   = cluster_center + cluster_half
else:
    zoom_start, zoom_end = 0, 5_000_000

fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=False)
chrom_color = CHROM_COLORS.get(densest_chrom, "#4C8EBF")

# Full chromosome
ax = axes[0]
ax.plot(mids, lg_data["density"].values, color="#aaaaaa", linewidth=0.4, alpha=0.5)
ax.plot(mids, smoothed.values, color=chrom_color, linewidth=1.2)
ax.axhline(mu, color="gray", linewidth=0.8, linestyle="--")
ax.fill_between(mids, mu - 2 * sd, mu + 2 * sd, color="gray", alpha=0.12)
for _, row in lg_hot.iterrows():
    ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=HOTSPOT_RED, alpha=0.25)
ax.set_ylabel("Variants per kb", fontsize=9)
ax.set_xlabel(f"Position on {densest_label} (Mb)", fontsize=9)
ax.set_title(f"{densest_label} ({densest_chrom}) — full chromosome", fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Zoomed: around the densest hotspot cluster
ax = axes[1]
zoom = lg_data[(lg_data["start"] >= zoom_start) & (lg_data["end"] <= zoom_end)].copy()
mid_z = (zoom["start"] + zoom["end"]) / 2 / 1e6
sm_z  = zoom["density"].rolling(5, center=True, min_periods=1).mean()
ax.plot(mid_z, zoom["density"].values, color="#aaaaaa", linewidth=0.4, alpha=0.5)
ax.plot(mid_z, sm_z.values, color=chrom_color, linewidth=1.5)
ax.axhline(mu, color="gray", linewidth=0.8, linestyle="--", label=f"chrom mean ({mu:.1f})")
ax.fill_between(mid_z, mu - 2 * sd, mu + 2 * sd, color="gray", alpha=0.12)
for _, row in lg_hot[(lg_hot["start"] >= zoom_start) & (lg_hot["end"] <= zoom_end)].iterrows():
    ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=HOTSPOT_RED, alpha=0.30,
               label=f"Hotspot (Z={row['max_z']:.1f})")
ax.set_ylabel("Variants per kb", fontsize=9)
ax.set_xlabel(f"Position on {densest_label} (Mb)", fontsize=9)
zoom_label = f"{zoom_start/1e6:.0f}–{zoom_end/1e6:.0f} Mb"
ax.set_title(f"{densest_label} — zoomed {zoom_label} (hotspot cluster, max Z = {max_z_val:.1f})", fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.suptitle(f"{densest_label} variation hotspot cluster", fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_04_densest_zoom.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 4 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 5 — Grid: variant density + hotspot shading for ALL active LGs
# ══════════════════════════════════════════════════════════════════════════════
print(f"Plot 5: density grid ({n_active} active LGs)...")

ncols  = 4
nrows  = math.ceil(n_active / ncols)
fw     = ncols * 5.5
fh     = nrows * 4.0

fig = plt.figure(figsize=(fw, fh))
gs  = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.55, wspace=0.35)

for idx, chrom in enumerate(active_chroms):
    row_i, col_i = divmod(idx, ncols)
    ax = fig.add_subplot(gs[row_i, col_i])

    sub  = df_main[df_main["chrom"] == chrom].copy().sort_values("start")
    hot  = hotspots[hotspots["chrom"] == chrom]
    mu   = sub["density"].mean()
    sd   = sub["density"].std(ddof=1)
    mids = (sub["start"] + sub["end"]) / 2 / 1e6
    sm   = sub["density"].rolling(5, center=True, min_periods=1).mean()

    ax.plot(mids, sub["density"].values, color="#cccccc", linewidth=0.4)
    ax.plot(mids, sm.values, color=CHROM_COLORS.get(chrom, "#4C8EBF"), linewidth=1.2)
    ax.axhline(mu, color="gray", linewidth=0.7, linestyle="--")
    if sd > 0:
        ax.fill_between(mids, mu - 2 * sd, mu + 2 * sd, color="gray", alpha=0.10)
    for _, row_h in hot.iterrows():
        ax.axvspan(row_h["start"] / 1e6, row_h["end"] / 1e6, color=HOTSPOT_RED, alpha=0.25)

    n_hot = int(hot.shape[0])
    max_z = hot["max_z"].max() if n_hot > 0 else 0
    ax.set_title(
        f"{chrom_label(chrom)} — {n_hot} hotspot regions  (max Z={max_z:.1f})",
        fontsize=8.5, fontweight="bold"
    )
    ax.set_xlabel("Position (Mb)", fontsize=7.5)
    ax.set_ylabel("Variants / kb", fontsize=7.5)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Hide unused subplot cells
n_cells = nrows * ncols
for empty_idx in range(n_active, n_cells):
    row_i, col_i = divmod(empty_idx, ncols)
    fig.add_subplot(gs[row_i, col_i]).set_visible(False)

fig.suptitle(
    f"Variant density across all {n_active} linkage groups with PASS variants\n"
    "Red shading = merged hotspot regions  |  Grey band = mean ± 2 SD",
    fontsize=12, fontweight="bold"
)
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_05_density_grid.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 5 saved")

print(f"\nAll plots saved to {OUT_DIR}/")
