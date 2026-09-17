#!/usr/bin/env python3
"""
Per-hotspot-region gene summary + zebrafish ortholog ENSDARG extraction.

Two independent outputs, reconstructed from v1's existing (previously
unscripted, ad hoc) hotspot_gene_summary.tsv/hotspot_zebrafish_ENSDARG.txt
by reverse-engineering the exact derivation from their content - verified
byte-identical (up to formatting) against v1 before being trusted for v2.

1. hotspot_gene_summary.tsv - one row per hotspot region that has >=1
   overlapping gene (from hotspot_gene_overlap.sh's hotspot_gene_overlaps.tsv),
   with the count of overlapping genes, how many have a real (non-LOC-
   placeholder) guppy gene symbol, and what those symbols are. Purely local
   - no external tool involved, despite living next to the gProfiler files.
   Sorted by max_z descending (matches v1).

2. hotspot_zebrafish_ENSDARG.txt - unique zebrafish Ensembl gene IDs
   (ortholog_ensg column, excluding "N/A") from the g:Orth CSV downloaded
   from gProfiler (https://biit.cs.ut.ee/gprofiler/orth) - see CLAUDE.md,
   "Hotspots under v2", for how that query was run. Only generated if a
   gProfiler_*.csv is present in HOT_DIR.

NOT reproduced here: v1's hotspot_genes_zebrafish.txt (172 lowercase gene
symbols) has no confirmed derivation - a plain lowercase of the CSV's
ortholog_name column does NOT reproduce it (spot-checked, mismatches on
~140/172 entries), so it was likely hand-edited or filtered through a step
that isn't recorded anywhere. Not worth guessing at silently.

Usage:
    python codes/analysis/hotspot_gene_summary.py
    REF_VERSION=v2 python codes/analysis/hotspot_gene_summary.py
"""

import glob
import os
import re

PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
REF_VERSION = os.environ.get("REF_VERSION", "v1")
OUT_SUFFIX = "" if REF_VERSION == "v1" else f"_{REF_VERSION}"
HOT_DIR = os.path.join(PROJECT_DIR, f"gatk/trimmomatic{OUT_SUFFIX}/hotspots")

CHROM_TO_LG = {
    "v1": {f"NC_0243{31 + i:02d}.1": f"LG{i + 1}" for i in range(23)},
    "v2": {f"NC_0888{30 + i:02d}.1": f"LG{i + 1}" for i in range(23)},
}[REF_VERSION]


def chrom_label(chrom):
    return CHROM_TO_LG.get(chrom, "Un")


# ══════════════════════════════════════════════════════════════════════════════
# 1. hotspot_gene_summary.tsv
# ══════════════════════════════════════════════════════════════════════════════
overlaps_path = os.path.join(HOT_DIR, "hotspot_gene_overlaps.tsv")
print(f"Reading {overlaps_path}...")

regions = {}  # (chrom, start, end) -> {"max_z":..., "max_density":..., "n_genes":0, "names":[...]}
with open(overlaps_path) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom, start, end, max_z, max_density = f[0], f[1], f[2], f[3], f[4]
        attrs = f[13]
        m = re.search(r"Name=([^;]+)", attrs)
        name = m.group(1) if m else None

        key = (chrom, start, end)
        if key not in regions:
            regions[key] = {"max_z": float(max_z), "max_density": max_density, "n_genes": 0, "names": []}
        regions[key]["n_genes"] += 1
        if name and not name.startswith("LOC"):
            regions[key]["names"].append(name)

rows = []
for (chrom, start, end), info in regions.items():
    named_sorted = sorted(info["names"])
    rows.append({
        "h_chrom": chrom, "h_start": start, "h_end": end,
        "max_z": info["max_z"], "max_density": float(info["max_density"]),
        "LG": chrom_label(chrom),
        "n_genes": info["n_genes"],
        "named_genes": len(named_sorted),
        "gene_names": ", ".join(named_sorted),
    })

rows.sort(key=lambda r: (-r["max_z"], int(r["h_start"])))

summary_path = os.path.join(HOT_DIR, "hotspot_gene_summary.tsv")
with open(summary_path, "w") as fh:
    fh.write("h_chrom\th_start\th_end\tmax_z\tmax_density\tLG\tn_genes\tnamed_genes\tgene_names\n")
    for r in rows:
        fh.write(f"{r['h_chrom']}\t{r['h_start']}\t{r['h_end']}\t{r['max_z']}\t{r['max_density']}\t"
                  f"{r['LG']}\t{r['n_genes']}\t{r['named_genes']}\t{r['gene_names']}\n")

print(f"  {len(rows)} regions with >=1 overlapping gene")
print(f"✅ {summary_path}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. hotspot_zebrafish_ENSDARG.txt (only if a gProfiler CSV is present)
# ══════════════════════════════════════════════════════════════════════════════
csvs = sorted(glob.glob(os.path.join(HOT_DIR, "gProfiler_*.csv")))
if not csvs:
    print("\nNo gProfiler_*.csv found in HOT_DIR - skipping ENSDARG extraction.")
else:
    csv_path = csvs[-1]  # most recent by name (timestamped)
    print(f"\nReading {csv_path}...")
    ensdarg = set()
    with open(csv_path) as fh:
        next(fh)  # header
        for line in fh:
            fields = [c.strip('"') for c in line.rstrip("\n").split(",")]
            ortholog_ensg = fields[5]
            if ortholog_ensg and ortholog_ensg != "N/A":
                ensdarg.add(ortholog_ensg)

    ensdarg_path = os.path.join(HOT_DIR, "hotspot_zebrafish_ENSDARG.txt")
    with open(ensdarg_path, "w") as fh:
        for e in sorted(ensdarg):
            fh.write(e + "\n")
    print(f"  {len(ensdarg)} unique zebrafish orthologs")
    print(f"✅ {ensdarg_path}")
