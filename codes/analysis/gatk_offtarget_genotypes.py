#!/usr/bin/env python3
"""
GATK genotype table at CRISPR predicted off-target sites.

Reads offtarget_variants.vcf.gz (produced by gatk SelectVariants restricted
to the 23 bp protospacer windows of the 8 predicted off-target sites) and
produces:
  - offtarget_genotypes.csv          — per-variant, per-sample genotype table
  - offtarget_genotype_heatmap.png   — genotype heatmap (0/0, 0/1, 1/1)

Key result (v1): GATK found only 2 SNPs at OT4 (NC_024332.1:5810651-5810674),
both present in the Control group → pre-existing Colombian-population polymorphisms,
not CRISPR-induced variants. The other 7 off-target protospacer windows have
no germline variants detectable against the Guanapo reference.

Note: the VCF was extracted from the raw joint-genotype VCF (all_samples.vcf.gz)
before VariantFiltration, so FILTER shows "." (unannotated). Both variants pass
all hard-filter thresholds (QD>2, FS<60, MQ>40) and would be PASS if re-filtered.

Dual-genome (v1/v2, see codes/genome_versions.sh): set REF_VERSION=v2 to run
against the GCF_904066995.2 off-target genotype VCF instead. Output goes to a
version-suffixed OUT_DIR so v1 and v2 results never overwrite each other.

Usage:
    python codes/analysis/gatk_offtarget_genotypes.py
    REF_VERSION=v2 python codes/analysis/gatk_offtarget_genotypes.py
"""

import gzip
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
REF_VERSION = os.environ.get("REF_VERSION", "v1")
OUT_SUFFIX  = "" if REF_VERSION == "v1" else f"_{REF_VERSION}"
VCF_PATH    = os.path.join(PROJECT_DIR,
              f"gatk/trimmomatic{OUT_SUFFIX}/vcf_offtargets/offtarget_variants.vcf.gz")
OUT_DIR     = os.path.join(PROJECT_DIR, f"codes/analysis/gatk_summary{OUT_SUFFIX}")
os.makedirs(OUT_DIR, exist_ok=True)
print(f"REF_VERSION={REF_VERSION}  VCF_PATH={VCF_PATH}  OUT_DIR={OUT_DIR}")

# ── Off-target site protospacer windows (23 bp, from offtargets_intervals.list)
# Same 8 predicted sites (identical guide/mismatch design), remapped to each
# reference assembly's own coordinates (crispresso{,_v2}/offtargets/combined/
# combined_offtargets.csv). v1's "locus" is an exon/intron/intergenic
# classification against the Guanapo annotation; v2's combine_offtargets run
# only recorded a distance label (no equivalent classification computed yet),
# so its "locus" is left as that same descriptive string instead of guessing.
OT_SITES_BY_VERSION = {
    "v1": {
        "OT1": {"chrom": "NC_024331.1", "start": 5708724,  "end": 5708747,  "mm": 3, "locus": "exon"},
        "OT2": {"chrom": "NC_024331.1", "start": 13951199, "end": 13951222, "mm": 4, "locus": "exon"},
        "OT3": {"chrom": "NC_024331.1", "start": 26228796, "end": 26228819, "mm": 4, "locus": "intergenic"},
        "OT4": {"chrom": "NC_024332.1", "start": 5810651,  "end": 5810674,  "mm": 4, "locus": "intron"},
        "OT5": {"chrom": "NC_024338.1", "start": 20512932, "end": 20512955, "mm": 4, "locus": "intergenic"},
        "OT6": {"chrom": "NC_024339.1", "start": 7034820,  "end": 7034843,  "mm": 4, "locus": "intron"},
        "OT7": {"chrom": "NC_024340.1", "start": 12200655, "end": 12200678, "mm": 4, "locus": "intergenic"},
        "OT8": {"chrom": "NC_024349.1", "start": 24882037, "end": 24882060, "mm": 4, "locus": "intron"},
    },
    "v2": {
        "OT1": {"chrom": "NC_088830.1", "start": 7340936,  "end": 7340959,  "mm": 3, "locus": "NC_088830.1 7.34 Mbp"},
        "OT2": {"chrom": "NC_088830.1", "start": 14098839, "end": 14098862, "mm": 4, "locus": "NC_088830.1 14.10 Mbp"},
        "OT3": {"chrom": "NC_088830.1", "start": 26619950, "end": 26619973, "mm": 4, "locus": "NC_088830.1 26.62 Mbp"},
        "OT4": {"chrom": "NC_088831.1", "start": 6154659,  "end": 6154682,  "mm": 4, "locus": "NC_088831.1 6.15 Mbp"},
        "OT5": {"chrom": "NC_088837.1", "start": 25867890, "end": 25867913, "mm": 4, "locus": "NC_088837.1 25.87 Mbp"},
        "OT6": {"chrom": "NC_088838.1", "start": 8510394,  "end": 8510417,  "mm": 4, "locus": "NC_088838.1 8.51 Mbp"},
        "OT7": {"chrom": "NC_088839.1", "start": 3391653,  "end": 3391676,  "mm": 4, "locus": "NC_088839.1 3.39 Mbp"},
        "OT8": {"chrom": "NC_088848.1", "start": 20545705, "end": 20545728, "mm": 4, "locus": "NC_088848.1 20.55 Mbp"},
    },
}
OT_SITES = OT_SITES_BY_VERSION[REF_VERSION]

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
# Helpers
# ══════════════════════════════════════════════════════════════════════════════
def assign_ot(chrom, pos):
    for sid, info in OT_SITES.items():
        if chrom == info["chrom"] and info["start"] <= pos <= info["end"]:
            return sid
    return "unknown"

def gt_code(gt_str):
    """0=hom_ref  1=het  2=hom_alt  -1=missing"""
    if gt_str in ("./.", ".|.", "."):
        return -1
    alleles = re.split(r"[/|]", gt_str.split(":")[0])
    try:
        a = [int(x) for x in alleles]
    except ValueError:
        return -1
    if all(x == 0 for x in a):
        return 0
    if all(x > 0 for x in a):
        return 2
    return 1

def parse_field(gt_str, fmt_keys, key):
    """Extract a FORMAT field value by key."""
    try:
        idx = fmt_keys.index(key)
        return gt_str.split(":")[idx]
    except (ValueError, IndexError):
        return "."

def ad_to_af(ad_str):
    """Return alt allele fraction from AD field (ref,alt)."""
    try:
        parts = [int(x) for x in ad_str.split(",")]
        total = sum(parts)
        return parts[1] / total if total > 0 else np.nan
    except (ValueError, IndexError):
        return np.nan


# ══════════════════════════════════════════════════════════════════════════════
# 1. Parse VCF
# ══════════════════════════════════════════════════════════════════════════════
print("Parsing off-target VCF...")

vcf_samples = None
records     = []

with gzip.open(VCF_PATH, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            vcf_samples = fields[9:]
            continue

        chrom  = fields[0]
        pos    = int(fields[1])
        ref    = fields[3]
        alt    = fields[4]
        qual   = fields[5]
        filt   = fields[6]
        fmt    = fields[8].split(":")

        vtype = "SNP" if len(ref) == len(alt) and "," not in alt else "INDEL"
        site  = assign_ot(chrom, pos)
        dist  = pos - OT_SITES[site]["start"] if site != "unknown" else np.nan

        rec = {
            "site":     site,
            "chrom":    chrom,
            "pos":      pos,
            "ref":      ref,
            "alt":      alt,
            "qual":     float(qual) if qual != "." else np.nan,
            "filter":   filt,
            "type":     vtype,
            "dist_from_site_start": int(dist) if not np.isnan(dist) else ".",
        }

        # Per-sample genotype columns
        for sample, gt_str in zip(vcf_samples, fields[9:]):
            code = gt_code(gt_str)
            ad   = parse_field(gt_str, fmt, "AD")
            dp   = parse_field(gt_str, fmt, "DP")
            gq   = parse_field(gt_str, fmt, "GQ")
            af   = ad_to_af(ad)
            rec[f"{SAMPLE_LABEL[sample]}_gt"]  = code
            rec[f"{SAMPLE_LABEL[sample]}_dp"]  = dp
            rec[f"{SAMPLE_LABEL[sample]}_af"]  = f"{af:.2f}" if not np.isnan(af) else "."

        records.append(rec)

print(f"  {len(records)} variant record(s) found")

# Reorder VCF samples to match GROUPS order
ordered_samples = SAMPLES  # already in group order

df = pd.DataFrame(records)

# ── Per-group summary columns ─────────────────────────────────────────────────
for grp, grp_samples in GROUPS.items():
    gt_cols = [f"{SAMPLE_LABEL[s]}_gt" for s in grp_samples]
    codes   = df[gt_cols].values   # shape (n_variants, n_samples_in_group)
    n_alleles = len(grp_samples) * 2
    # count alt alleles: het=1, hom_alt=2; missing=-1 excluded
    alt_count = np.where(codes == 1, 1, np.where(codes == 2, 2, 0)).sum(axis=1)
    missing   = (codes == -1).sum(axis=1)
    valid_alleles = n_alleles - missing * 2
    af_grp    = np.where(valid_alleles > 0, alt_count / valid_alleles, np.nan)
    df[f"{grp}_AF"]      = af_grp.round(3)
    df[f"{grp}_n_het"]   = (codes == 1).sum(axis=1)
    df[f"{grp}_n_homalt"]= (codes == 2).sum(axis=1)

# ── Background classification ─────────────────────────────────────────────────
# "background" = alt allele present in Control; "candidate_edit" = absent in Control
df["classification"] = np.where(df["Control_AF"] > 0, "background_variant", "candidate_edit")

csv_path = os.path.join(OUT_DIR, "offtarget_genotypes.csv")
df.to_csv(csv_path, index=False)
print(f"✅ Genotype CSV → {csv_path}")

# Console summary
print("\n=== Off-target variant summary ===")
for _, row in df.iterrows():
    info = OT_SITES.get(row["site"], {})
    print(f"\n  {row['site']}  {row['chrom']}:{row['pos']}  "
          f"REF={row['ref']} ALT={row['alt']}  "
          f"QUAL={row['qual']:.0f}  {row['type']}  locus={info.get('locus','')}  "
          f"dist_from_start={row['dist_from_site_start']} bp")
    for grp in GROUPS:
        print(f"    {grp:<12} AF={row[f'{grp}_AF']:.3f}  "
              f"het={row[f'{grp}_n_het']}  hom_alt={row[f'{grp}_n_homalt']}")
    print(f"    → {row['classification']}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Genotype heatmap — all 8 sites shown; one row per variant, one placeholder
#    row for sites with no variants detected.
# ══════════════════════════════════════════════════════════════════════════════
# Map: -1=missing, 0=0/0, 1=0/1, 2=1/1
# Colors: gray, white, cornflowerblue, steelblue
CMAP   = ListedColormap(["#CCCCCC", "#F7F7F7", "#87CEEB", "#1E6FA5"])
BOUNDS = [-1.5, -0.5, 0.5, 1.5, 2.5]
NORM   = BoundaryNorm(BOUNDS, CMAP.N)

gt_labels      = []   # y-axis labels
gt_matrix      = []   # numeric codes for cell color
gt_annotations = []   # list of (code, af_str) per cell
row_types      = []   # "variant" | "empty"
row_class      = []   # classification label for right axis
site_dividers  = []   # y-positions for horizontal dividers between sites

y = 0
for site_id, info in OT_SITES.items():
    site_vars = df[df["site"] == site_id]

    if len(site_vars) == 0:
        # All samples are hom-ref at this site — no VCF record exists
        label = (f"{site_id}  {info['chrom']}:{info['start']}–{info['end']}\n"
                 f"no variants detected  ({info['mm']}mm · {info['locus']})")
        gt_labels.append(label)
        gt_matrix.append([0] * len(ordered_samples))
        gt_annotations.append([(0, ".")] * len(ordered_samples))
        row_types.append("empty")
        row_class.append("no variants detected")
        y += 1
    else:
        for _, row in site_vars.iterrows():
            label = (f"{site_id}  {info['chrom']}:{row['pos']}\n"
                     f"REF={row['ref']} → ALT={row['alt']}  "
                     f"({info['mm']}mm · {info['locus']})")
            gt_labels.append(label)
            codes = [int(row[f"{SAMPLE_LABEL[s]}_gt"]) for s in ordered_samples]
            gt_matrix.append(codes)
            gt_annotations.append([
                (codes[si], row[f"{SAMPLE_LABEL[s]}_af"])
                for si, s in enumerate(ordered_samples)
            ])
            row_types.append("variant")
            row_class.append(row["classification"])
            y += 1

    site_dividers.append(y - 0.5)   # horizontal line after last row of this site

gt_matrix = np.array(gt_matrix, dtype=float)
n_rows  = len(gt_labels)
n_samps = len(ordered_samples)

fig_h   = max(6, n_rows * 1.1 + 3.0)
fig, ax = plt.subplots(figsize=(17, fig_h))

im = ax.imshow(gt_matrix, aspect="auto", cmap=CMAP, norm=NORM,
               interpolation="nearest")

# ── x-axis: sample labels colored by group ────────────────────────────────────
ax.set_xticks(range(n_samps))
ax.set_xticklabels([SAMPLE_LABEL[s] for s in ordered_samples],
                   rotation=40, ha="right", fontsize=10)
ax.tick_params(axis='x', pad=6)
for tick, s in zip(ax.get_xticklabels(), ordered_samples):
    tick.set_color(GROUP_COLORS[GROUP_OF[s]])

# ── y-axis: row labels; dim empty-site labels ─────────────────────────────────
ax.set_yticks(range(n_rows))
ax.set_yticklabels(gt_labels, fontsize=10)
for i, (tick, rtype) in enumerate(zip(ax.get_yticklabels(), row_types)):
    tick.set_color("#999999" if rtype == "empty" else "black")

# ── Cell annotations ──────────────────────────────────────────────────────────
gt_str_map = {-1: "./.", 0: "0/0", 1: "0/1", 2: "1/1"}
for ri in range(n_rows):
    for si in range(n_samps):
        code, af_val = gt_annotations[ri][si]
        if row_types[ri] == "empty":
            # show "0/0" dimly for clarity without clutter
            ax.text(si, ri, "0/0", ha="center", va="center",
                    fontsize=7, color="#BBBBBB")
        else:
            txt = gt_str_map[code]
            if code > 0 and af_val != ".":
                txt += f"\n{af_val}"
            color = "white" if code == 2 else ("#888888" if code == -1 else "black")
            ax.text(si, ri, txt, ha="center", va="center", fontsize=8.5, color=color)

# ── Horizontal dividers between OT sites ──────────────────────────────────────
for yd in site_dividers[:-1]:   # skip last (after OT8, at bottom edge)
    ax.axhline(yd, color="#666666", linewidth=0.8, linestyle="--")

# ── Vertical dividers between groups ─────────────────────────────────────────
xpos = -0.5
for grp, grp_samples in GROUPS.items():
    xpos += len(grp_samples)
    if xpos < n_samps - 0.5:
        ax.axvline(xpos, color="white", linewidth=2.5)

# ── Group labels below x-tick labels ─────────────────────────────────────────
xpos = -0.5
for grp, grp_samples in GROUPS.items():
    mid = xpos + len(grp_samples) / 2
    ax.text(mid, -0.22, grp, ha="center", va="top", fontsize=12,
            color=GROUP_COLORS[grp], fontweight="bold",
            transform=ax.get_xaxis_transform())
    xpos += len(grp_samples)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color="#CCCCCC", label="./. (missing)"),
    mpatches.Patch(facecolor="#F7F7F7", edgecolor="#AAAAAA",
                   label="0/0 (hom ref)"),
    mpatches.Patch(color="#87CEEB", label="0/1 (heterozygous)"),
    mpatches.Patch(color="#1E6FA5", label="1/1 (hom alt)"),
]
ax.legend(handles=legend_patches, loc="upper right", fontsize=10,
          bbox_to_anchor=(1.22, 1.08), frameon=True)

# ── Right axis: classification ────────────────────────────────────────────────
CLASS_COLORS = {
    "no variants detected": "#2E8B57",
    "background_variant":   "#555555",
    "candidate_edit":       "#B22222",
}
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(range(n_rows))
ax2.set_yticklabels(
    [c.replace("_", " ") for c in row_class],
    fontsize=10,
)
for tick, cls in zip(ax2.get_yticklabels(), row_class):
    tick.set_color(CLASS_COLORS.get(cls, "#555555"))

ax.set_title(
    "GATK genotypes at predicted CRISPR off-target sites\n"
    "Guppy WGS · 8 protospacer windows (23 bp each) · "
    f"{len(records)} variant(s) detected across all sites",
    fontsize=13, fontweight="bold", pad=14,
)

plt.savefig(os.path.join(OUT_DIR, "offtarget_genotype_heatmap.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print(f"✅ Genotype heatmap → {os.path.join(OUT_DIR, 'offtarget_genotype_heatmap.png')}")

print("\n=== Interpretation ===")
print(f"  REF_VERSION: {REF_VERSION}")
print(f"  Total variants at 8 off-target protospacer windows: {len(records)}")
bg = (df["classification"] == "background_variant").sum()
ce = (df["classification"] == "candidate_edit").sum()
print(f"  Background (present in Control): {bg}")
print(f"  Candidate edits (absent in Control): {ce}")

sites_with_variants = sorted(set(df["site"]) - {"unknown"})
sites_without = [s for s in OT_SITES if s not in sites_with_variants]
if sites_without:
    print(f"  {', '.join(sites_without)}: 0 variants → no confounding germline SNPs")
for s in sites_with_variants:
    sub = df[df["site"] == s]
    control_gt_cols = [f"{SAMPLE_LABEL[x]}_gt" for x in GROUPS["Control"]]
    carriers = [c for c in ["Control_AF", "Only_MNP_AF", "Plasmid_Ko_AF", "RNP_Cas_AF"]
                if (sub[c] > 0).any()]
    cls = sorted(sub["classification"].unique())
    print(f"  {s} ({OT_SITES[s]['chrom']}): {len(sub)} variant(s), "
          f"classification={cls}, AF>0 in groups={[c.replace('_AF','') for c in carriers]}")
