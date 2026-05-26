#!/usr/bin/env python3
"""
On-target vs off-target CRISPR editing efficiency comparison.

Reads CRISPResso2 on-target and CRISPRessoWGS off-target quantification
outputs and produces three figures + one summary CSV.

Usage (from any directory):
    python codes/analysis/plot_editing_comparison.py

Inputs:
    crispresso/ontarget/trimmomatic/{SAMPLE}/CRISPResso_on_{SAMPLE}/
        CRISPResso_quantification_of_editing_frequency.txt
    crispresso/wgs/trimmomatic/{SAMPLE}/CRISPRessoWGS_on_{SAMPLE}/
        SAMPLES_QUANTIFICATION_SUMMARY.txt

Outputs (codes/analysis/editing_comparison/):
    editing_summary.csv       — combined Modified% table (all samples × sites)
    editing_heatmap.png       — heatmap: on-target (left) + off-target (right)
    ontarget_barplot.png      — on-target editing efficiency per group
    offtarget_dotplot.png     — off-target Modified% per site (once WGS data available)
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR  = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
ONTARGET_DIR = os.path.join(PROJECT_DIR, "crispresso/ontarget/trimmomatic")
WGS_DIR      = os.path.join(PROJECT_DIR, "crispresso/wgs/trimmomatic")
OUT_DIR      = os.path.join(PROJECT_DIR, "codes/analysis/editing_comparison")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Samples and groups ─────────────────────────────────────────────────────────
GROUPS = {
    "Control":    ["Control_MNP_I_S54_L002",  "Control_MNP_II_S55_L002",  "Control_MNP_III_S56_L002"],
    "Only_MNP":   ["Only_MNP_C1_S57_L002",   "Only_MNP_C2_S58_L002",    "Only_MNP_C3_S59_L002",   "Only_MNP_C4_S60_L002"],
    "Plasmid_Ko": ["Plasmid_Ko_P1_S61_L002",  "Plasmid_Ko_P2_S62_L002",  "Plasmid_Ko_P3_S63_L002", "Plasmid_Ko_P4_S64_L002"],
    "RNP_Cas":    ["RNP_Cas1_S65_L002",       "RNP_Cas2_S66_L002",       "RNP_Cas3_S67_L002",      "RNP_Cas4_S68_L002"],
}
SAMPLES       = [s for slist in GROUPS.values() for s in slist]
SAMPLE_LABEL  = {s: re.sub(r"_S\d+_L002$", "", s) for s in SAMPLES}
GROUP_OF      = {s: g for g, slist in GROUPS.items() for s in slist}
GROUP_COLORS  = {
    "Control":    "#4C8EBF",
    "Only_MNP":   "#E8A838",
    "Plasmid_Ko": "#60B260",
    "RNP_Cas":    "#D94040",
}

# ── Off-target sites (name as in SAMPLES_QUANTIFICATION_SUMMARY.txt) ──────────
# Key  : site name used by CRISPRessoWGS
# Value: short label for plots
OT_SITES = {
    "NC_024331.1_5708724_3mm_CRISPOR":  "OT1\nNC_024331.1:5708724\n3mm · exon",
    "NC_024331.1_13951199_4mm_CRISPOR": "OT2\nNC_024331.1:13951199\n4mm · exon",
    "NC_024331.1_26228796_4mm_CRISPOR": "OT3\nNC_024331.1:26228796\n4mm · intergenic",
    "NC_024332.1_5810651_4mm_CRISPOR":  "OT4\nNC_024332.1:5810651\n4mm · intron",
    "NC_024338.1_20512932_4mm_CRISPOR": "OT5\nNC_024338.1:20512932\n4mm · intergenic",
    "NC_024339.1_7034820_4mm_CRISPOR":  "OT6\nNC_024339.1:7034820\n4mm · intron",
    "NC_024340.1_12200655_4mm_CRISPOR": "OT7\nNC_024340.1:12200655\n4mm · intergenic",
    "NC_024349.1_24882037_4mm_CRISPOR": "OT8\nNC_024349.1:24882037\n4mm · intron",
}
OT_NAMES = list(OT_SITES.keys())


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Read on-target data
# ═══════════════════════════════════════════════════════════════════════════════
print("Reading on-target CRISPResso results...")

def read_ontarget(sample):
    """Return (weighted_modified_pct, total_reads_aligned) for one sample.
    The on-target run uses two amplicons (bdnf_fwd and bdnf_rc); this function
    computes a weighted average so both strands contribute equally.
    """
    path = os.path.join(
        ONTARGET_DIR, sample,
        f"CRISPResso_on_{sample}",
        "CRISPResso_quantification_of_editing_frequency.txt",
    )
    if not os.path.isfile(path):
        return np.nan, 0
    df = pd.read_csv(path, sep="\t")
    total = df["Reads_aligned"].sum()
    if total == 0:
        return np.nan, 0
    weighted = (df["Modified%"] * df["Reads_aligned"]).sum() / total
    return weighted, int(total)

ontarget = {}
for s in SAMPLES:
    mod_pct, reads = read_ontarget(s)
    ontarget[s] = {"modified_pct": mod_pct, "reads": reads}
    status = f"{mod_pct:.2f}%" if not np.isnan(mod_pct) else "MISSING"
    print(f"  {SAMPLE_LABEL[s]:<22} {status}  (n={reads})")

ontarget_df = pd.DataFrame(ontarget).T  # index = sample


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Read off-target WGS data
# ═══════════════════════════════════════════════════════════════════════════════
print("\nReading CRISPRessoWGS off-target results...")

def read_wgs(sample):
    """Return DataFrame indexed by site name with columns Modified% and Reads_total."""
    path = os.path.join(
        WGS_DIR, sample,
        f"CRISPRessoWGS_on_{sample}",
        "SAMPLES_QUANTIFICATION_SUMMARY.txt",
    )
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path, sep="\t")
    df["Modified%"]   = pd.to_numeric(df["Modified%"],   errors="coerce")
    df["Reads_total"] = pd.to_numeric(df["Reads_total"], errors="coerce")
    return df.set_index("Name")

wgs_mod   = pd.DataFrame(index=SAMPLES, columns=OT_NAMES, dtype=float)
wgs_reads = pd.DataFrame(index=SAMPLES, columns=OT_NAMES, dtype=float)

for s in SAMPLES:
    df = read_wgs(s)
    for site in OT_NAMES:
        if df is not None and site in df.index:
            wgs_mod.loc[s, site]   = df.loc[site, "Modified%"]
            wgs_reads.loc[s, site] = df.loc[site, "Reads_total"]

has_wgs = wgs_mod.notna().any().any()
if has_wgs:
    print(f"  WGS data available for {wgs_mod.notna().any(axis=1).sum()}/{len(SAMPLES)} samples")
else:
    print("  ⚠  All WGS values are NA — CRISPRessoWGS re-run may still be in progress")
    print("     (Re-run with fixed 103-bp BED; re-run this script after job completes)")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Summary CSV
# ═══════════════════════════════════════════════════════════════════════════════
rows = []
for s in SAMPLES:
    row = {
        "sample": SAMPLE_LABEL[s],
        "group":  GROUP_OF[s],
        "on_target_modified_pct": ontarget_df.loc[s, "modified_pct"],
        "on_target_reads":        ontarget_df.loc[s, "reads"],
    }
    for site in OT_NAMES:
        col = site.split("_CRISPOR")[0]   # shorter column name
        row[f"{col}_modified_pct"] = wgs_mod.loc[s, site]
        row[f"{col}_reads"]        = wgs_reads.loc[s, site]
    rows.append(row)

summary_df = pd.DataFrame(rows)
csv_path = os.path.join(OUT_DIR, "editing_summary.csv")
summary_df.to_csv(csv_path, index=False)
print(f"\n✅ Summary CSV → {csv_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4a. Heatmap: on-target (left panel) + off-target (right panel)
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax_on, ax_ot) = plt.subplots(
    1, 2, figsize=(20, 7),
    gridspec_kw={"width_ratios": [1, len(OT_NAMES)]},
)

ylabels = [SAMPLE_LABEL[s] for s in SAMPLES]
group_sizes = [len(v) for v in GROUPS.values()]
dividers = [sum(group_sizes[:i]) - 0.5 for i in range(1, len(group_sizes))]

# ── Left: on-target ────────────────────────────────────────────────────────────
on_vals = ontarget_df.loc[SAMPLES, "modified_pct"].astype(float).values.reshape(-1, 1)
im_on   = ax_on.imshow(on_vals, aspect="auto", cmap="Reds", vmin=0, vmax=100)

ax_on.set_xticks([0])
ax_on.set_xticklabels(["On-target\nbdnf (NC_024333.1)\n:15922039–15922058"], fontsize=8)
ax_on.set_yticks(range(len(SAMPLES)))
ax_on.set_yticklabels(ylabels, fontsize=8)

for tick, s in zip(ax_on.get_yticklabels(), SAMPLES):
    tick.set_color(GROUP_COLORS[GROUP_OF[s]])

for i, s in enumerate(SAMPLES):
    v = ontarget_df.loc[s, "modified_pct"]
    r = int(ontarget_df.loc[s, "reads"])
    if pd.notna(v):
        txt_color = "white" if v > 55 else "black"
        ax_on.text(0, i, f"{v:.1f}%\n(n={r})", ha="center", va="center",
                   fontsize=7, color=txt_color)
    else:
        ax_on.text(0, i, "N/A", ha="center", va="center", fontsize=7, color="#999999")

for d in dividers:
    ax_on.axhline(d, color="white", linewidth=2)

plt.colorbar(im_on, ax=ax_on, label="Modified reads (%)", shrink=0.55, pad=0.02)
ax_on.set_title("On-target\nediting (%)", fontsize=9, fontweight="bold")

# ── Right: off-target ──────────────────────────────────────────────────────────
ot_vals  = wgs_mod.loc[SAMPLES, OT_NAMES].astype(float).values
vmax_ot  = max(5.0, np.nanmax(ot_vals)) if has_wgs else 5.0
im_ot    = ax_ot.imshow(ot_vals, aspect="auto", cmap="Oranges", vmin=0, vmax=vmax_ot)

ax_ot.set_xticks(range(len(OT_NAMES)))
ax_ot.set_xticklabels([OT_SITES[s] for s in OT_NAMES], fontsize=7)
ax_ot.set_yticks(range(len(SAMPLES)))
ax_ot.set_yticklabels([])

for i, s in enumerate(SAMPLES):
    for j, site in enumerate(OT_NAMES):
        v = wgs_mod.loc[s, site]
        r = wgs_reads.loc[s, site]
        if pd.notna(v):
            txt_color = "white" if v > vmax_ot * 0.6 else "black"
            ax_ot.text(j, i, f"{v:.2f}%", ha="center", va="center",
                       fontsize=6, color=txt_color)
        else:
            reads_str = f"n={int(r)}" if pd.notna(r) else ""
            ax_ot.text(j, i, f"NA\n{reads_str}", ha="center", va="center",
                       fontsize=5.5, color="#999999")

for d in dividers:
    ax_ot.axhline(d, color="white", linewidth=2)

plt.colorbar(im_ot, ax=ax_ot, label=f"Modified reads (%, scale 0–{vmax_ot:.1f})",
             shrink=0.55, pad=0.02)
ax_ot.set_title(
    "Off-target editing (%)\n[CRISPRessoWGS · ±50 bp window · CRISPOR predicted sites]",
    fontsize=9, fontweight="bold",
)

if not has_wgs:
    ax_ot.text(
        len(OT_NAMES) / 2 - 0.5, len(SAMPLES) / 2 - 0.5,
        "⚠  WGS re-run in progress\n(re-run this script after job finishes)",
        ha="center", va="center", fontsize=11, color="#555555",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.9),
    )

# Group legend
patches = [mpatches.Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
fig.legend(handles=patches, loc="lower center", ncol=4, fontsize=9,
           bbox_to_anchor=(0.5, -0.02), frameon=False)

fig.suptitle(
    "CRISPR-Cas9 editing: on-target (bdnf) vs predicted off-target sites\n"
    "Guppy (Poecilia reticulata) · SpCas9 · sgRNA: TGAGAGACGCCCCGGGCATG (− strand)",
    fontsize=10, fontweight="bold", y=1.02,
)
plt.tight_layout()
heatmap_path = os.path.join(OUT_DIR, "editing_heatmap.png")
plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Heatmap         → {heatmap_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4b. Bar plot: on-target editing efficiency per sample, grouped
# ═══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 5))

x = 0
xtick_pos, xtick_labels = [], []
group_label_pos = []

for grp, samples in GROUPS.items():
    grp_x_start = x
    for s in samples:
        v = ontarget_df.loc[s, "modified_pct"]
        height = v if pd.notna(v) else 0
        ax.bar(x, height, color=GROUP_COLORS[grp], alpha=0.80, width=0.7, edgecolor="white")
        if pd.notna(v):
            ax.text(x, height + 0.2, f"{v:.1f}%", ha="center", va="bottom", fontsize=7)
        xtick_pos.append(x)
        xtick_labels.append(SAMPLE_LABEL[s])
        x += 1
    group_label_pos.append(((grp_x_start + x - 1) / 2, grp))
    ax.axvline(x - 0.5, color="#cccccc", linewidth=1, linestyle="--")
    x += 0.5

ax.set_xticks(xtick_pos)
ax.set_xticklabels(xtick_labels, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("Modified reads (%)", fontsize=10)
ax.set_title(
    "On-target editing efficiency per sample\n"
    "bdnf sgRNA site · NC_024333.1:15922039–15922058",
    fontsize=10,
)
ax.set_xlim(-0.5, x - 0.5)
ymax = max(10, ontarget_df["modified_pct"].dropna().max() * 1.2 if ontarget_df["modified_pct"].notna().any() else 10)
ax.set_ylim(0, ymax)

for xpos, grp in group_label_pos:
    ax.text(xpos, ymax * 0.95, grp, ha="center", fontsize=9,
            color=GROUP_COLORS[grp], fontweight="bold")

ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
bar_path = os.path.join(OUT_DIR, "ontarget_barplot.png")
plt.savefig(bar_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ On-target barplot → {bar_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# 4c. Off-target dot plot (per site, coloured by group) — only if WGS data exist
# ═══════════════════════════════════════════════════════════════════════════════
dot_path = os.path.join(OUT_DIR, "offtarget_dotplot.png")

print(f"  has_wgs={has_wgs}  (non-null WGS cells: {int(wgs_mod.notna().sum().sum())})")

if has_wgs:
    fig, ax = plt.subplots(figsize=(14, 6))

    for j, site in enumerate(OT_NAMES):
        for s in SAMPLES:
            v = wgs_mod.loc[s, site]
            if pd.notna(v):
                ax.scatter(j, float(v), color=GROUP_COLORS[GROUP_OF[s]], s=60,
                           alpha=0.8, zorder=3)

    ax.set_xticks(range(len(OT_NAMES)))
    ax.set_xticklabels([OT_SITES[n] for n in OT_NAMES], fontsize=7)
    ax.set_ylabel("Modified reads (%)", fontsize=10)
    ax.set_title(
        "Off-target editing rates per predicted site (all samples)\n"
        "CRISPRessoWGS · SpCas9 sgRNA TGAGAGACGCCCCGGGCATG",
        fontsize=10,
    )
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_xlim(-0.5, len(OT_NAMES) - 0.5)
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    patches = [mpatches.Patch(color=c, label=g) for g, c in GROUP_COLORS.items()]
    ax.legend(handles=patches, fontsize=8, frameon=False)

    plt.savefig(dot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"✅ Off-target dotplot → {dot_path}")
else:
    print(f"⏭  Off-target dotplot skipped (no WGS data yet) — will be generated on re-run")

# ── Console summary ────────────────────────────────────────────────────────────
print("\n=== On-target editing rates by group ===")
for grp, samples in GROUPS.items():
    vals = [ontarget_df.loc[s, "modified_pct"] for s in samples]
    valid = [v for v in vals if pd.notna(v)]
    pct_str = "  ".join(f"{v:.2f}%" if pd.notna(v) else "NA" for v in vals)
    mean_str = f"mean={np.mean(valid):.2f}%" if valid else "mean=NA"
    print(f"  {grp:<12} {pct_str}   ({mean_str})")

if has_wgs:
    print("\n=== Off-target Modified% summary (mean across samples) ===")
    for site in OT_NAMES:
        vals = wgs_mod[site].dropna()
        label = OT_SITES[site].split("\n")[0]
        print(f"  {label:<5}  mean={vals.mean():.3f}%  max={vals.max():.3f}%  n_valid={len(vals)}")
