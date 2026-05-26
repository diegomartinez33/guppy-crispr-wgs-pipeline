#!/usr/bin/env python3
"""
GATK variant summary statistics.

Runs bcftools stats on the hard-filtered SNP and INDEL VCFs,
parses per-sample counts, and produces:
  - gatk_variant_summary.csv   — per-sample variant count table
  - gatk_summary_barplot.png   — SNP + INDEL counts per sample (grouped)
  - gatk_titv_boxplot.png      — Ti/Tv ratio per group

Usage (requires module load bcftools/1.15.1):
    python codes/analysis/gatk_variant_summary.py
"""

import subprocess
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
VCF_DIR     = os.path.join(PROJECT_DIR, "gatk/trimmomatic/vcf_filtered")
OUT_DIR     = os.path.join(PROJECT_DIR, "codes/analysis/gatk_summary")
os.makedirs(OUT_DIR, exist_ok=True)

SNP_VCF   = os.path.join(VCF_DIR, "snps_filtered.vcf.gz")
INDEL_VCF = os.path.join(VCF_DIR, "indels_filtered.vcf.gz")

# ── Sample / group config ──────────────────────────────────────────────────────
GROUPS = {
    "Control":    ["Control_MNP_I_S54_L002",  "Control_MNP_II_S55_L002",  "Control_MNP_III_S56_L002"],
    "Only_MNP":   ["Only_MNP_C1_S57_L002",   "Only_MNP_C2_S58_L002",    "Only_MNP_C3_S59_L002",   "Only_MNP_C4_S60_L002"],
    "Plasmid_Ko": ["Plasmid_Ko_P1_S61_L002",  "Plasmid_Ko_P2_S62_L002",  "Plasmid_Ko_P3_S63_L002", "Plasmid_Ko_P4_S64_L002"],
    "RNP_Cas":    ["RNP_Cas1_S65_L002",       "RNP_Cas2_S66_L002",       "RNP_Cas3_S67_L002",      "RNP_Cas4_S68_L002"],
}
SAMPLES      = [s for sl in GROUPS.values() for s in sl]
SAMPLE_LABEL = {s: re.sub(r"_S\d+_L002$", "", s) for s in SAMPLES}
GROUP_OF     = {s: g for g, sl in GROUPS.items() for s in sl}
GROUP_COLORS = {
    "Control": "#4C8EBF", "Only_MNP": "#E8A838",
    "Plasmid_Ko": "#60B260", "RNP_Cas": "#D94040",
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Run bcftools stats and parse PSC (per-sample counts)
# ══════════════════════════════════════════════════════════════════════════════
# PSC column layout (tab-separated fields after "PSC"):
#   [1]=id  [2]=sample  [3]=nRefHom  [4]=nHet  [5]=nHomAlt
#   [6]=nTs  [7]=nTv  [8]=nIndel  [9]=avgDepth  [10]=nSingletons  [13]=nMissing
# Note: for an INDEL-only VCF, nHet/nHomAlt are 0 (those are SNP-specific);
#       indel counts are in [8]=nIndel.

def run_stats(vcf_path):
    """Return bcftools stats stdout as a string."""
    result = subprocess.run(
        ["bcftools", "stats", "-s", "-", vcf_path],
        capture_output=True, text=True, check=True
    )
    return result.stdout

def parse_psc(stats_text):
    """Parse PSC lines → dict keyed by sample name."""
    psc = {}
    for line in stats_text.splitlines():
        if not line.startswith("PSC"):
            continue
        f = line.split("\t")
        sample    = f[2]
        n_het     = int(f[4])
        n_hom_alt = int(f[5])
        n_ts      = int(f[6])
        n_tv      = int(f[7])
        n_indel   = int(f[8])
        avg_depth = float(f[9])
        n_missing = int(f[13])
        psc[sample] = {
            "n_het":     n_het,
            "n_hom_alt": n_hom_alt,
            "n_snps":    n_het + n_hom_alt,
            "n_ts":      n_ts,
            "n_tv":      n_tv,
            "titv":      n_ts / n_tv if n_tv > 0 else np.nan,
            "n_indels":  n_indel,
            "avg_depth": avg_depth,
            "n_missing": n_missing,
        }
    return psc

def parse_tstv(stats_text):
    """Return overall Ti/Tv ratio (PASS variants)."""
    for line in stats_text.splitlines():
        if line.startswith("TSTV"):
            f = line.split("\t")
            return float(f[7])   # ratio for PASS variants
    return np.nan


print("Running bcftools stats on SNP VCF...")
snp_stats  = run_stats(SNP_VCF)
snp_psc    = parse_psc(snp_stats)
titv_all   = parse_tstv(snp_stats)
print(f"  Genome-wide Ti/Tv (PASS SNPs): {titv_all:.3f}")

print("Running bcftools stats on INDEL VCF...")
indel_stats = run_stats(INDEL_VCF)
indel_psc   = parse_psc(indel_stats)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Build summary DataFrame
# ══════════════════════════════════════════════════════════════════════════════
rows = []
for s in SAMPLES:
    sp = snp_psc.get(s, {})
    ip = indel_psc.get(s, {})
    rows.append({
        "sample":      SAMPLE_LABEL[s],
        "group":       GROUP_OF[s],
        "n_snps":      sp.get("n_snps",    np.nan),
        "n_snps_het":  sp.get("n_het",     np.nan),
        "n_snps_hom":  sp.get("n_hom_alt", np.nan),
        "titv":        sp.get("titv",      np.nan),
        "n_indels":    ip.get("n_indels",  np.nan),
        "avg_depth":   sp.get("avg_depth", np.nan),
        "n_missing":   sp.get("n_missing", np.nan),
    })

summary = pd.DataFrame(rows)
csv_path = os.path.join(OUT_DIR, "gatk_variant_summary.csv")
summary.to_csv(csv_path, index=False)
print(f"\n✅ Summary CSV → {csv_path}")

print("\n=== Per-group means ===")
for grp in GROUPS:
    sub = summary[summary["group"] == grp]
    print(f"  {grp:<12}  SNPs: {sub['n_snps'].mean():>10,.0f}  "
          f"INDELs: {sub['n_indels'].mean():>8,.0f}  "
          f"Ti/Tv: {sub['titv'].mean():.3f}  "
          f"Depth: {sub['avg_depth'].mean():.1f}×")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Bar plot: SNP and INDEL counts per sample, grouped
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

x        = 0
xtick_pos, xtick_labels = [], []
grp_label_pos = []
dividers  = []

# Collect x positions per sample
xpos_map = {}
for grp, samples in GROUPS.items():
    grp_start = x
    for s in samples:
        xpos_map[s] = x
        xtick_pos.append(x)
        xtick_labels.append(SAMPLE_LABEL[s])
        x += 1
    grp_label_pos.append(((grp_start + x - 1) / 2, grp))
    dividers.append(x - 0.5)
    x += 0.6

# ── Top: SNPs (stacked het / hom alt) ─────────────────────────────────────────
ax = axes[0]
for s in SAMPLES:
    row   = summary[summary["sample"] == SAMPLE_LABEL[s]].iloc[0]
    col   = GROUP_COLORS[GROUP_OF[s]]
    xp    = xpos_map[s]
    ax.bar(xp, row["n_snps_het"],  width=0.7, color=col, alpha=0.85, label="_")
    ax.bar(xp, row["n_snps_hom"],  width=0.7, color=col, alpha=0.45,
           bottom=row["n_snps_het"], label="_")

for d in dividers[:-1]:
    ax.axvline(d, color="#cccccc", linewidth=1, linestyle="--")
ax.set_ylabel("Number of SNPs", fontsize=10)
ax.set_title("Hard-filtered SNPs per sample  (dark = heterozygous, light = homozygous alt)",
             fontsize=10)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v/1e6:.1f}M"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Bottom: INDELs ────────────────────────────────────────────────────────────
ax = axes[1]
for s in SAMPLES:
    row = summary[summary["sample"] == SAMPLE_LABEL[s]].iloc[0]
    col = GROUP_COLORS[GROUP_OF[s]]
    ax.bar(xpos_map[s], row["n_indels"], width=0.7, color=col, alpha=0.8)

for d in dividers[:-1]:
    ax.axvline(d, color="#cccccc", linewidth=1, linestyle="--")
ax.set_ylabel("Number of INDELs", fontsize=10)
ax.set_title("Hard-filtered INDELs per sample", fontsize=10)
ax.set_xticks(xtick_pos)
ax.set_xticklabels(xtick_labels, rotation=40, ha="right", fontsize=8)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v/1e3:.0f}k"))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Group labels on top panel
ymax_snp = axes[0].get_ylim()[1]
for xpos, grp in grp_label_pos:
    axes[0].text(xpos, ymax_snp * 0.97, grp, ha="center", fontsize=9,
                 color=GROUP_COLORS[grp], fontweight="bold")

patches = [mpatches.Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
fig.legend(handles=patches, loc="lower center", ncol=4, fontsize=9,
           bbox_to_anchor=(0.5, -0.01), frameon=False)

fig.suptitle(
    "GATK variant counts — Guppy WGS (Poecilia reticulata, Colombian population)\n"
    f"Reference: GCF_000633615.1 (Guanapo)  ·  Genome-wide Ti/Tv = {titv_all:.2f}",
    fontsize=10, fontweight="bold",
)
plt.savefig(os.path.join(OUT_DIR, "gatk_summary_barplot.png"), dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Barplot → {os.path.join(OUT_DIR, 'gatk_summary_barplot.png')}")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Ti/Tv boxplot per group
# ══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))

group_names = list(GROUPS.keys())
data = [summary[summary["group"] == g]["titv"].dropna().values for g in group_names]
colors = [GROUP_COLORS[g] for g in group_names]

bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                medianprops=dict(color="black", linewidth=2))
for patch, col in zip(bp["boxes"], colors):
    patch.set_facecolor(col)
    patch.set_alpha(0.75)

# overlay individual points
for i, (vals, col) in enumerate(zip(data, colors), start=1):
    ax.scatter([i] * len(vals), vals, color=col, s=60, zorder=3, edgecolors="white")

ax.axhline(titv_all, color="gray", linewidth=1, linestyle="--",
           label=f"Genome-wide mean ({titv_all:.2f})")
ax.set_xticks(range(1, len(group_names) + 1))
ax.set_xticklabels(group_names, fontsize=10)
ax.set_ylabel("Ti/Tv ratio", fontsize=10)
ax.set_title("Transition / Transversion ratio per group\n(hard-filtered SNPs)", fontsize=10)
ax.legend(fontsize=9, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.savefig(os.path.join(OUT_DIR, "gatk_titv_boxplot.png"), dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Ti/Tv boxplot → {os.path.join(OUT_DIR, 'gatk_titv_boxplot.png')}")
