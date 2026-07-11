#!/usr/bin/env python3
"""
Variation hotspot detection.

Reads window_counts.tsv produced by hotspot_windows.sh, computes variant
density and per-chromosome Z-scores, calls hotspots (BH-corrected FDR < 0.05
or Z > 4.0 as fallback), merges adjacent windows, and produces:

  window_counts_annotated.csv       — all windows with density, z_score, is_hotspot
  hotspots.bed                      — merged hotspot regions (max_z, max_density)
  hotspot_manhattan_genome.png      — genome-wide Manhattan plot
  hotspot_plots/hotspot_density_*.png — one density plot per chromosome

Usage (after hotspot_windows.sh completes; requires module load bcftools bedtools):
    python codes/analysis/hotspot_analysis.py
"""

import re
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
VCF_DIR     = os.path.join(PROJECT_DIR, "gatk/trimmomatic/vcf_filtered")
HOT_DIR     = os.path.join(PROJECT_DIR, "gatk/trimmomatic/hotspots")
PLOT_DIR    = os.path.join(HOT_DIR, "hotspot_plots")
os.makedirs(PLOT_DIR, exist_ok=True)

WINDOW_TSV = os.path.join(HOT_DIR, "window_counts.tsv")
Z_THRESH   = 4.0

# ── Chromosome colors: alternating blue/amber across all 24 contigs ────────────
_PALETTE = ["#4C8EBF", "#E8A838"]
CHROM_COLORS = {}
for _i in range(23):
    CHROM_COLORS[f"NC_0243{31 + _i:02d}.1"] = _PALETTE[_i % 2]
CHROM_COLORS["NW_007615013.1"] = _PALETTE[23 % 2]


def chrom_label(chrom):
    """NC_024331.1 → LG1;  NW_007615013.1 → Un."""
    m = re.match(r"NC_024(3\d{2})\.1", chrom)
    return f"LG{int(m.group(1)) - 330}" if m else "Un"


def get_contig_lengths(vcf_path):
    """Return (ordered_list, {chrom: length}) from VCF ##contig header lines."""
    result = subprocess.run(
        ["bcftools", "view", "-h", vcf_path],
        capture_output=True, text=True, check=True,
    )
    order, lengths = [], {}
    for line in result.stdout.splitlines():
        if not line.startswith("##contig"):
            continue
        m = re.search(r"ID=([^,>]+),length=(\d+)", line)
        if m:
            cid = m.group(1)
            if cid not in lengths:
                order.append(cid)
            lengths[cid] = int(m.group(2))
    return order, lengths


# ══════════════════════════════════════════════════════════════════════════════
# 1. Contig lengths (for cumulative X positions in Manhattan plot)
# ══════════════════════════════════════════════════════════════════════════════
print("Parsing contig lengths from VCF header...")
chroms_ordered, contig_lengths = get_contig_lengths(
    os.path.join(VCF_DIR, "snps_filtered.vcf.gz")
)
print(f"  {len(contig_lengths)} contigs found")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Load window counts and compute density / Z-scores
# ══════════════════════════════════════════════════════════════════════════════
print(f"Loading {WINDOW_TSV}...")
df = pd.read_csv(
    WINDOW_TSV, sep="\t",
    names=["chrom", "start", "end", "snp_count", "indel_count", "total_count"],
)

df["win_size"] = df["end"] - df["start"]
df["density"]  = df["total_count"] / df["win_size"] * 1000   # variants per kb

# Per-chromosome Z-score (population-specific baseline)
df["z_score"] = np.nan
for chrom, grp in df.groupby("chrom"):
    mu = grp["density"].mean()
    sd = grp["density"].std(ddof=1)
    if sd > 0:
        df.loc[grp.index, "z_score"] = (grp["density"] - mu) / sd
    else:
        df.loc[grp.index, "z_score"] = 0.0

# Hotspot calling
try:
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import norm
    pvals = norm.sf(df["z_score"].fillna(0).values)
    reject, _, _, _ = multipletests(pvals, method="fdr_bh")
    df["is_hotspot"] = reject
    n_method = "BH-corrected FDR < 0.05"
except ImportError:
    df["is_hotspot"] = df["z_score"] > Z_THRESH
    n_method = f"Z > {Z_THRESH}"

n_hot = int(df["is_hotspot"].sum())
print(f"  Hotspot windows: {n_hot:,}  ({n_method})")

csv_path = os.path.join(HOT_DIR, "window_counts_annotated.csv")
df.to_csv(csv_path, index=False)
print(f"✅ Annotated CSV → {csv_path}")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Merge adjacent hotspot windows → hotspots.bed
# ══════════════════════════════════════════════════════════════════════════════
_raw_bed = os.path.join(HOT_DIR, "_hotspot_raw.bed")
hotspot_bed_out = os.path.join(HOT_DIR, "hotspots.bed")

hot_df = df[df["is_hotspot"]].copy().sort_values(["chrom", "start"])
hot_df[["chrom", "start", "end", "z_score", "density"]].to_csv(
    _raw_bed, sep="\t", index=False, header=False,
)

merge_result = subprocess.run(
    ["bedtools", "merge", "-d", "2000", "-c", "4,5", "-o", "max,max",
     "-i", _raw_bed],
    capture_output=True, text=True, check=True,
)

with open(hotspot_bed_out, "w") as fh:
    fh.write("chrom\tstart\tend\tmax_z\tmax_density\n")
    fh.write(merge_result.stdout)

n_regions = merge_result.stdout.strip().count("\n") + (1 if merge_result.stdout.strip() else 0)
print(f"✅ Hotspots BED ({n_regions} merged regions) → {hotspot_bed_out}")
os.remove(_raw_bed)


# ══════════════════════════════════════════════════════════════════════════════
# 4. Cumulative genomic positions for Manhattan X-axis
# ══════════════════════════════════════════════════════════════════════════════
cum_offset = {}
offset = 0
for chrom in chroms_ordered:
    cum_offset[chrom] = offset
    offset += contig_lengths.get(chrom, 0)

df["cum_mid"] = df.apply(
    lambda r: cum_offset.get(r["chrom"], 0) + (r["start"] + r["end"]) / 2,
    axis=1,
)


# ══════════════════════════════════════════════════════════════════════════════
# 5. Manhattan plot (genome-wide)
# ══════════════════════════════════════════════════════════════════════════════
print("Generating Manhattan plot...")

fig, ax = plt.subplots(figsize=(18, 5))

for chrom in chroms_ordered:
    sub = df[df["chrom"] == chrom]
    if sub.empty:
        continue
    base_color = CHROM_COLORS.get(chrom, "#888888")
    hot_mask   = sub["is_hotspot"]
    ax.scatter(
        sub.loc[~hot_mask, "cum_mid"], sub.loc[~hot_mask, "z_score"],
        color=base_color, s=1, alpha=0.4, rasterized=True,
    )
    ax.scatter(
        sub.loc[hot_mask, "cum_mid"], sub.loc[hot_mask, "z_score"],
        color="#D94040", s=4, alpha=0.8, rasterized=True,
    )

ax.axhline(Z_THRESH, color="red", linewidth=1, linestyle="--", label=f"Z = {Z_THRESH}")

xtick_pos, xtick_lab = [], []
for chrom in chroms_ordered:
    if chrom not in contig_lengths:
        continue
    mid = cum_offset[chrom] + contig_lengths[chrom] / 2
    xtick_pos.append(mid)
    xtick_lab.append(chrom_label(chrom))

ax.set_xticks(xtick_pos)
ax.set_xticklabels(xtick_lab, fontsize=7, rotation=45, ha="right")
ax.set_ylabel("Z-score (variant density)", fontsize=10)
ax.set_xlabel("Chromosome (linkage group)", fontsize=10)
ax.set_title(
    "Variation hotspot detection — Guppy WGS Colombian population\n"
    f"10 kb windows, 2 kb step  |  Red points: hotspot windows ({n_hot:,})",
    fontsize=10, fontweight="bold",
)
ax.legend(fontsize=9, frameon=False, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
out_path = os.path.join(HOT_DIR, "hotspot_manhattan_genome.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Manhattan plot → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Per-chromosome density plots
# ══════════════════════════════════════════════════════════════════════════════
print("Generating per-chromosome density plots...")

for chrom in chroms_ordered:
    sub = df[df["chrom"] == chrom].copy().sort_values("start")
    if sub.empty or len(sub) < 5:
        continue

    mu = sub["density"].mean()
    sd = sub["density"].std(ddof=1)

    fig, ax = plt.subplots(figsize=(14, 4))
    midpoints = (sub["start"] + sub["end"]) / 2 / 1e6   # Mb

    ax.plot(midpoints, sub["density"], color="#aaaaaa", linewidth=0.5, alpha=0.6)

    smoothed = sub["density"].rolling(5, center=True, min_periods=1).mean()
    color    = CHROM_COLORS.get(chrom, "#4C8EBF")
    ax.plot(midpoints, smoothed, color=color, linewidth=1.2)

    ax.axhline(mu, color="gray", linewidth=0.8, linestyle="--")
    if sd > 0:
        ax.fill_between(midpoints, mu - 2 * sd, mu + 2 * sd,
                        color="gray", alpha=0.12)

    for _, row in sub[sub["is_hotspot"]].iterrows():
        ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color="#D94040", alpha=0.25)

    ax.set_xlabel("Position (Mb)", fontsize=9)
    ax.set_ylabel("Variants per kb", fontsize=9)
    ax.set_title(f"Variant density — {chrom_label(chrom)} ({chrom})", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    safe_chrom = chrom.replace(".", "_")
    out_path   = os.path.join(PLOT_DIR, f"hotspot_density_{safe_chrom}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

print(f"✅ Per-chromosome plots → {PLOT_DIR}/")
print("\nHotspot analysis complete.")
