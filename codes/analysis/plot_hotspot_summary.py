#!/usr/bin/env python3
"""
Hotspot analysis summary plots.

Reads window_counts_annotated.csv and hotspots.bed produced by
hotspot_analysis.py and generates five publication-ready figures:

  hotspot_summary_01_variants_per_chrom.png  — total SNP+INDEL counts per LG
  hotspot_summary_02_snp_indel_breakdown.png — stacked SNP vs INDEL per LG
  hotspot_summary_03_hotspot_regions.png     — merged hotspot regions ranked by Z-score
  hotspot_summary_04_lg5_zoom.png            — zoomed density plot for LG5 (densest region)
  hotspot_summary_05_density_grid.png        — 2×3 grid of density plots for LG1–LG5

Usage:
    python codes/analysis/plot_hotspot_summary.py
"""

import re
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
HOT_DIR     = os.path.join(PROJECT_DIR, "gatk/trimmomatic/hotspots")
OUT_DIR     = os.path.join(HOT_DIR, "summary_plots")
os.makedirs(OUT_DIR, exist_ok=True)

ANNOTATED_CSV  = os.path.join(HOT_DIR, "window_counts_annotated.csv")
HOTSPOTS_BED   = os.path.join(HOT_DIR, "hotspots.bed")

# ── Colour palette ─────────────────────────────────────────────────────────────
PALETTE = ["#4C8EBF", "#E8A838"]
HOTSPOT_RED  = "#D94040"
SNP_COLOR    = "#3A7AB5"
INDEL_COLOR  = "#E07B30"
CHROM_COLORS = {}
for _i in range(23):
    CHROM_COLORS[f"NC_0243{31 + _i:02d}.1"] = PALETTE[_i % 2]
CHROM_COLORS["NW_007615013.1"] = PALETTE[23 % 2]

MAIN_CHROMS = [f"NC_024{331 + i:03d}.1" for i in range(23)] + ["NW_007615013.1"]


def chrom_label(chrom):
    m = re.match(r"NC_024(3\d{2})\.1", chrom)
    return f"LG{int(m.group(1)) - 330}" if m else "Un"


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
)
by_chrom["label"] = by_chrom["chrom"].apply(chrom_label)

hotspots = pd.read_csv(HOTSPOTS_BED, sep="\t")
hotspots["label"] = hotspots["chrom"].apply(chrom_label)
hotspots["width_kb"] = (hotspots["end"] - hotspots["start"]) / 1000
hotspots["mid_mb"] = (hotspots["start"] + hotspots["end"]) / 2 / 1e6

print(f"  {len(df_main):,} windows loaded  |  {len(hotspots)} hotspot regions")


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

legend_patches = [
    mpatches.Patch(color=CHROM_COLORS["NC_024331.1"], label="LG1"),
    mpatches.Patch(color=CHROM_COLORS["NC_024332.1"], label="LG2"),
    mpatches.Patch(color=CHROM_COLORS["NC_024333.1"], label="LG3"),
    mpatches.Patch(color=CHROM_COLORS["NC_024334.1"], label="LG4"),
    mpatches.Patch(color=CHROM_COLORS["NC_024335.1"], label="LG5"),
]
ax.legend(handles=legend_patches, fontsize=8, frameon=False, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_03_hotspot_regions.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 3 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — LG5 zoom (highest Z-score chromosome; cluster at 0–762 kb)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 4: LG5 zoom...")

lg5 = df_main[df_main["chrom"] == "NC_024335.1"].copy().sort_values("start")
lg5_hot = hotspots[hotspots["chrom"] == "NC_024335.1"]
mu5 = lg5["density"].mean()
sd5 = lg5["density"].std(ddof=1)
mid5 = (lg5["start"] + lg5["end"]) / 2 / 1e6
smoothed5 = lg5["density"].rolling(5, center=True, min_periods=1).mean()

fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=False)

# Full chromosome
ax = axes[0]
ax.plot(mid5, lg5["density"].values, color="#aaaaaa", linewidth=0.4, alpha=0.5)
ax.plot(mid5, smoothed5.values, color=CHROM_COLORS["NC_024335.1"], linewidth=1.2)
ax.axhline(mu5, color="gray", linewidth=0.8, linestyle="--")
ax.fill_between(mid5, mu5 - 2 * sd5, mu5 + 2 * sd5, color="gray", alpha=0.12)
for _, row in lg5_hot.iterrows():
    ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=HOTSPOT_RED, alpha=0.25)
ax.set_ylabel("Variants per kb", fontsize=9)
ax.set_xlabel("Position on LG5 (Mb)", fontsize=9)
ax.set_title("LG5 (NC_024335.1) — full chromosome", fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Zoomed: 0–5 Mb where all hotspots sit
ax = axes[1]
zoom = lg5[lg5["start"] <= 5_000_000].copy()
mid_z = (zoom["start"] + zoom["end"]) / 2 / 1e6
sm_z  = zoom["density"].rolling(5, center=True, min_periods=1).mean()
ax.plot(mid_z, zoom["density"].values, color="#aaaaaa", linewidth=0.4, alpha=0.5)
ax.plot(mid_z, sm_z.values, color=CHROM_COLORS["NC_024335.1"], linewidth=1.5)
ax.axhline(mu5, color="gray", linewidth=0.8, linestyle="--", label=f"chrom mean ({mu5:.1f})")
ax.fill_between(mid_z, mu5 - 2 * sd5, mu5 + 2 * sd5, color="gray", alpha=0.12)
for _, row in lg5_hot[lg5_hot["start"] <= 5_000_000].iterrows():
    ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=HOTSPOT_RED, alpha=0.30,
               label=f"Hotspot (Z={row['max_z']:.1f})")
ax.set_ylabel("Variants per kb", fontsize=9)
ax.set_xlabel("Position on LG5 (Mb)", fontsize=9)
ax.set_title("LG5 — zoomed 0–5 Mb (hotspot cluster, max Z = 12.0)", fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.suptitle("LG5 variation hotspot cluster", fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_04_lg5_zoom.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 4 saved")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 5 — 2×3 grid: density + hotspot shading for LG1–LG5
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 5: density grid LG1–LG5...")

active_chroms = [c for c in MAIN_CHROMS if by_chrom.loc[by_chrom["chrom"] == c, "total_variants"].values[0] > 0]
n = len(active_chroms)   # 5
ncols = 3
nrows = 2

fig = plt.figure(figsize=(18, 9))
gs  = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.45, wspace=0.35)

for idx, chrom in enumerate(active_chroms):
    row, col = divmod(idx, ncols)
    ax = fig.add_subplot(gs[row, col])

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

    n_hot  = int(hot.shape[0])
    max_z  = hot["max_z"].max() if n_hot > 0 else 0
    ax.set_title(
        f"{chrom_label(chrom)} — {n_hot} hotspot regions  (max Z={max_z:.1f})",
        fontsize=9, fontweight="bold"
    )
    ax.set_xlabel("Position (Mb)", fontsize=8)
    ax.set_ylabel("Variants / kb", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Hide unused subplot (slot [1,2] is empty since we have 5 chroms in a 2×3 grid)
if n < nrows * ncols:
    fig.add_subplot(gs[1, 2]).set_visible(False)

fig.suptitle(
    "Variant density across linkage groups with PASS variants\n"
    "Red shading = merged hotspot regions  |  Grey band = mean ± 2 SD",
    fontsize=12, fontweight="bold"
)
plt.savefig(os.path.join(OUT_DIR, "hotspot_summary_05_density_grid.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  ✅ Plot 5 saved")

print(f"\nAll plots saved to {OUT_DIR}/")
