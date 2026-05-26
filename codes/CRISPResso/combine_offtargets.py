#!/usr/bin/env python3
"""
Combina off-targets de Cas-OFFinder y CRISPOR
para análisis con CRISPRessoWGS
"""

import pandas as pd
import os

# ── Configuración ─────────────────────────────────────────────────────────────
PROJECT_DIR  = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
CASOFF_FILE  = os.path.join(PROJECT_DIR, "crispresso/offtargets/offtargets.txt")
CRISPOR_FILE = os.path.join(PROJECT_DIR, "data/offtargets/crispor/crispor_326forw_casoffinder_format.txt")
OUTPUT_DIR   = os.path.join(PROJECT_DIR, "crispresso/offtargets/combined")
SGRNA        = "TGAGAGACGCCCCGGGCATG"
TOLERANCE    = 2

# Window padding for the CRISPRessoWGS BED file.
#
# CRISPRessoWGSCORE.write_trimmed_fastq trims each read to [bpstart, bpend]
# by checking that BOTH coordinates appear in the read's reference positions.
# That means a single read must span the entire window end-to-end:
#
#   max usable window = READ_LENGTH - 1  (≈ 149 bp for 150 bp reads)
#
# The BED window = protospacer (23 bp) + 2 × WGS_PADDING.
# With WGS_PADDING = 40 → window = 103 bp  ✓  (< 149 bp)
#
# The original value was PADDING = 100 → window = 223 bp, wider than one read
# → CRISPRessoWGS reported 0 reads for every site despite the BAM having reads.
#
# To change the window, adjust WGS_PADDING here and re-run this script to
# regenerate offtargets_crispresso_wgs.bed before re-submitting crispresso_wgs.sh.
READ_LENGTH  = 150   # Illumina NovaSeq X 2×150 bp paired-end
WGS_PADDING  = 40    # bp of flanking context on each side of the protospacer

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Leer Cas-OFFinder ─────────────────────────────────────────────────────────
print("Reading Cas-OFFinder results...")
casoff_rows = []

with open(CASOFF_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts      = line.strip().split("\t")
        mismatches = int(parts[-1])
        strand     = parts[-2]
        offtarget  = parts[-3]
        position   = int(parts[-4])
        chrom      = parts[1].split()[0]

        casoff_rows.append({
            "chromosome":    chrom,
            "start":         position,
            "end":           position + 23,
            "strand":        strand,
            "offtarget_seq": offtarget,
            "mismatches":    mismatches,
            "mit_score":     None,
            "cfd_score":     None,
            "locus":         None,
            "source":        "CasOFFinder"
        })

casoff_df = pd.DataFrame(casoff_rows)
print(f"Cas-OFFinder off-targets: {len(casoff_df)}")

# ── Leer CRISPOR ──────────────────────────────────────────────────────────────
print("\nReading CRISPOR results...")
crispor_df = pd.read_csv(CRISPOR_FILE, sep="\t", comment="#")
crispor_df = crispor_df.rename(columns={
    "Chromosome": "chromosome",
    "Position":   "start",
    "Strand":     "strand",
    "Sequence":   "offtarget_seq",
    "Mismatches": "mismatches",
    "MIT_Score":  "mit_score",
    "CFD_Score":  "cfd_score",
    "Locus":      "locus"
})
crispor_df["end"]    = crispor_df["start"] + 23
crispor_df["source"] = "CRISPOR"
print(f"CRISPOR off-targets: {len(crispor_df)}")

# ── Combinar ──────────────────────────────────────────────────────────────────
print("\nCombining results...")

# Suprimir FutureWarning de pandas
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

combined = pd.concat([casoff_df, crispor_df], ignore_index=True)
combined = combined.sort_values(
    ["chromosome", "start", "strand"]
).reset_index(drop=True)

# ── Detectar duplicados con tolerancia ±2bp ───────────────────────────────────
print(f"Detecting near-duplicates (tolerance: ±{TOLERANCE}bp)...")

combined["duplicate_group"] = -1
group_id = 0

for i, row_i in combined.iterrows():
    if combined.loc[i, "duplicate_group"] != -1:
        continue
    combined.loc[i, "duplicate_group"] = group_id
    for j, row_j in combined.iterrows():
        if i == j or combined.loc[j, "duplicate_group"] != -1:
            continue
        if (row_i["chromosome"] == row_j["chromosome"] and
            row_i["strand"]     == row_j["strand"] and
            abs(row_i["start"]  - row_j["start"]) <= TOLERANCE):
            combined.loc[j, "duplicate_group"] = group_id
    group_id += 1

# Marcar sitios encontrados por ambas herramientas
group_counts = combined["duplicate_group"].value_counts()
multi_groups = group_counts[group_counts > 1].index
combined["found_by_both"] = combined["duplicate_group"].isin(multi_groups)

# ── Prioridad explícita: CRISPOR primero ──────────────────────────────────────
# CRISPOR tiene MIT score, CFD score y anotación de locus
combined["source_priority"] = combined["source"].map({
    "CRISPOR":     0,   # ← prioridad alta
    "CasOFFinder": 1    # ← prioridad baja
})

combined = combined.sort_values(
    ["duplicate_group", "source_priority"],
    ascending=[True, True]
)

combined_dedup = combined.drop_duplicates(
    subset="duplicate_group", keep="first"
).sort_values(["mismatches", "chromosome", "start"]).reset_index(drop=True)

# Limpiar columna auxiliar
combined_dedup = combined_dedup.drop(columns=["source_priority", "duplicate_group"])

only_casoff  = (~combined_dedup["found_by_both"] &
                (combined_dedup["source"] == "CasOFFinder")).sum()
only_crispor = (~combined_dedup["found_by_both"] &
                (combined_dedup["source"] == "CRISPOR")).sum()
both         = combined_dedup["found_by_both"].sum()

print(f"\nTotal unique off-targets:  {len(combined_dedup)}")
print(f"Only in CasOFFinder:       {only_casoff}")
print(f"Only in CRISPOR:           {only_crispor}")
print(f"Found by both:             {both}")

# ── Excluir on-target ─────────────────────────────────────────────────────────
ontarget_mask   = ((combined_dedup["chromosome"] == "NC_024333.1") &
                   (combined_dedup["mismatches"]  == 0))
offtargets_only = combined_dedup[~ontarget_mask].copy()
print(f"\nOn-target excluded → Off-targets for WGS: {len(offtargets_only)}")

# ── Output 1: CSV completo ────────────────────────────────────────────────────
CSV_OUT = os.path.join(OUTPUT_DIR, "combined_offtargets.csv")
combined_dedup.to_csv(CSV_OUT, index=False)
print(f"\n✅ CSV: {CSV_OUT}")

# ── Output 2: BED para CRISPRessoWGS ─────────────────────────────────────────
WGS_BED_OUT  = os.path.join(OUTPUT_DIR, "offtargets_crispresso_wgs.bed")
PROTOSPACER  = 23                       # sgRNA (20 bp) + PAM (3 bp)
WINDOW_SIZE  = PROTOSPACER + 2 * WGS_PADDING

assert WINDOW_SIZE < READ_LENGTH, (
    f"WGS_PADDING={WGS_PADDING}: window ({WINDOW_SIZE} bp) >= read length "
    f"({READ_LENGTH} bp). CRISPRessoWGS needs window < read length so that "
    f"individual reads can span bpstart–bpend. Reduce WGS_PADDING."
)

with open(WGS_BED_OUT, "w") as f:
    for _, row in offtargets_only.iterrows():
        name         = (f"{row['chromosome']}_"
                        f"{row['start']}_"
                        f"{row['mismatches']}mm_"
                        f"{row['source']}")
        start_padded = max(0, row["start"] - WGS_PADDING)
        end_padded   = row["end"] + WGS_PADDING
        # Use the actual off-target guide (20 bp, no PAM) not the original sgRNA.
        # CRISPRessoWGS passes column-5 as -g to inner CRISPResso; the original
        # sgRNA has 3-4 mismatches vs the off-target amplicon and causes exit code 2.
        guide = row["offtarget_seq"][:20]
        f.write(f"{row['chromosome']}\t"
                f"{start_padded}\t"
                f"{end_padded}\t"
                f"{name}\t"
                f"{guide}\n")

print(f"✅ CRISPRessoWGS BED: {WGS_BED_OUT}  (window: {WINDOW_SIZE} bp)")

# ── Output 3: BED para IGV ────────────────────────────────────────────────────
IGV_BED_OUT = os.path.join(OUTPUT_DIR, "combined_offtargets_igv.bed")

with open(IGV_BED_OUT, "w") as f:
    f.write("# Combined off-targets: Cas-OFFinder + CRISPOR\n")
    f.write("# chr\tstart\tend\tname\tscore\tstrand\n")
    for _, row in offtargets_only.iterrows():
        name  = f"{row['chromosome']}_{row['start']}_{row['mismatches']}mm"
        score = int(row["mit_score"] * 1000) if pd.notna(row["mit_score"]) else 0
        f.write(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t"
                f"{name}\t{score}\t{row['strand']}\n")

print(f"✅ IGV BED: {IGV_BED_OUT}")

# ── Output 4: Intervalos GATK ────────────────────────────────────────────────
INTERVALS_OUT = os.path.join(OUTPUT_DIR, "offtargets_intervals.list")

with open(INTERVALS_OUT, "w") as f:
    for _, row in offtargets_only.iterrows():
        f.write(f"{row['chromosome']}:{row['start']}-{row['end']}\n")

print(f"✅ GATK intervals: {INTERVALS_OUT}")

# ── Resumen final ─────────────────────────────────────────────────────────────
print("\n=== Off-targets para análisis WGS ===")
print(offtargets_only[[
    "chromosome", "start", "strand", "mismatches",
    "offtarget_seq", "mit_score", "cfd_score",
    "source", "found_by_both", "locus"
]].to_string(index=False))

print(f"\n=== Por mismatches ===")
print(offtargets_only["mismatches"].value_counts().sort_index())

print(f"\n=== Por fuente ===")
print(offtargets_only["source"].value_counts())

print(f"\n=== Por cromosoma ===")
print(offtargets_only["chromosome"].value_counts())

print(f"\n=== Archivos generados ===")
print(f"  combined_offtargets.csv          → tabla completa")
print(f"  offtargets_crispresso_wgs.bed    → input CRISPRessoWGS")
print(f"  combined_offtargets_igv.bed      → visualización IGV")
print(f"  offtargets_intervals.list        → intervalos GATK")