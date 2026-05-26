#!/usr/bin/env python3
"""
Convierte output CRISPOR (.xls) a BED y formato compatible
con Cas-OFFinder y CRISPRessoWGS
Filtra por guideId especificado
"""

import pandas as pd
import sys
import os

# ── Configuración ─────────────────────────────────────────────────────────────
PROJECT_DIR="/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
CRISPOR_FILE = os.path.join(PROJECT_DIR, "data/crispor_offtargets.xls")
OUTPUT_DIR   = os.path.join(PROJECT_DIR, "data/offtargets/crispor")
GUIDE_ID     = "326forw"        # guideId a filtrar
MAX_MISMATCHES = 4              # máximo mismatches a incluir

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Leer archivo CRISPOR ──────────────────────────────────────────────────────
# Saltar las filas de metadata (# Name, # Sequence, etc.)
# La fila de headers está en la fila 9 (índice 8)
print(f"Leyendo archivo: {CRISPOR_FILE}")

df = pd.read_excel(
    CRISPOR_FILE,
    header=8,           # fila 9 como header (índice 8)
    engine="xlrd"       # para archivos .xls antiguos
)

print(f"Total de off-targets en archivo: {len(df)}")
print(f"Columnas encontradas: {list(df.columns)}")

# ── Filtrar por guideId ───────────────────────────────────────────────────────
df_filtered = df[df["guideId"] == GUIDE_ID].copy()
print(f"Off-targets para '{GUIDE_ID}': {len(df_filtered)}")

# Filtrar por número de mismatches
df_filtered = df_filtered[df_filtered["mismatchCount"] <= MAX_MISMATCHES]
print(f"Off-targets con <= {MAX_MISMATCHES} mismatches: {len(df_filtered)}")

# ── Limpiar y estandarizar columnas ──────────────────────────────────────────
# Renombrar columnas para mayor claridad
df_filtered = df_filtered.rename(columns={
    "guideId":          "guide_id",
    "guideSeq":         "guide_seq",
    "offtargetSeq":     "offtarget_seq",
    "mismatchCount":    "mismatches",
    "mitOfftargetScore":"mit_score",
    "cfdOfftargetScore":"cfd_score",
    "chrom":            "chromosome",
    "start":            "start",
    "end":              "end",
    "strand":           "strand",
    "locusDesc":        "locus_description"
})

# Asegurar tipos correctos
df_filtered["start"] = pd.to_numeric(df_filtered["start"], errors="coerce").astype(int)
df_filtered["end"]   = pd.to_numeric(df_filtered["end"],   errors="coerce").astype(int)

# Ordenar por score MIT descendente
df_filtered = df_filtered.sort_values("mit_score", ascending=False)

# ── Output 1: CSV completo filtrado ───────────────────────────────────────────
CSV_OUT = os.path.join(OUTPUT_DIR, f"crispor_{GUIDE_ID}_offtargets.csv")
df_filtered.to_csv(CSV_OUT, index=False)
print(f"\n✅ CSV completo: {CSV_OUT}")

# ── Output 2: Archivo BED para CRISPRessoWGS ─────────────────────────────────
# Formato: chr  start  end  name  score  strand
BED_OUT = os.path.join(OUTPUT_DIR, f"crispor_{GUIDE_ID}_offtargets.bed")

with open(BED_OUT, "w") as f:
    f.write("# BED file generado desde CRISPOR para CRISPRessoWGS\n")
    f.write(f"# guideId: {GUIDE_ID}\n")
    f.write(f"# guideSeq: {df_filtered['guide_seq'].iloc[0]}\n")
    f.write(f"# Max mismatches: {MAX_MISMATCHES}\n")

    for _, row in df_filtered.iterrows():
        name = f"{row['chromosome']}_{row['start']}_{row['mismatches']}mm"
        score = round(row["mit_score"] * 1000)  # escalar para BED
        f.write(f"{row['chromosome']}\t"
                f"{row['start']}\t"
                f"{row['end']}\t"
                f"{name}\t"
                f"{score}\t"
                f"{row['strand']}\n")

print(f"✅ BED file: {BED_OUT}")

# ── Output 3: Formato Cas-OFFinder compatible ─────────────────────────────────
# Para combinar con resultados de Cas-OFFinder
CASOFF_OUT = os.path.join(OUTPUT_DIR, f"crispor_{GUIDE_ID}_casoffinder_format.txt")

with open(CASOFF_OUT, "w") as f:
    f.write("# Converted from CRISPOR to Cas-OFFinder format\n")
    f.write("Chromosome\tPosition\tStrand\tSequence\tMismatches\tMIT_Score\tCFD_Score\tLocus\n")

    for _, row in df_filtered.iterrows():
        f.write(f"{row['chromosome']}\t"
                f"{row['start']}\t"
                f"{row['strand']}\t"
                f"{row['offtarget_seq']}\t"
                f"{row['mismatches']}\t"
                f"{row['mit_score']:.6f}\t"
                f"{row['cfd_score']:.6f}\t"
                f"{row['locus_description']}\n")

print(f"✅ Cas-OFFinder format: {CASOFF_OUT}")

# ── Output 4: Intervalos para GATK SelectVariants ────────────────────────────
INTERVALS_OUT = os.path.join(OUTPUT_DIR, f"crispor_{GUIDE_ID}_intervals.list")

with open(INTERVALS_OUT, "w") as f:
    for _, row in df_filtered.iterrows():
        f.write(f"{row['chromosome']}:{row['start']}-{row['end']}\n")

print(f"✅ GATK intervals: {INTERVALS_OUT}")

# ── Resumen estadístico ───────────────────────────────────────────────────────
print(f"\n=== Resumen para {GUIDE_ID} ===")
print(f"Total off-targets:          {len(df_filtered)}")
print(f"\nDistribución por mismatches:")
print(df_filtered["mismatches"].value_counts().sort_index().to_string())
print(f"\nDistribución por cromosoma:")
print(df_filtered["chromosome"].value_counts().to_string())
print(f"\nDistribución por tipo de locus:")
locus_type = df_filtered["locus_description"].str.split(":").str[0]
print(locus_type.value_counts().to_string())
print(f"\nTop 10 off-targets por MIT score:")
print(df_filtered[["chromosome", "start", "end", "strand",
                    "mismatches", "mit_score", "cfd_score",
                    "locus_description"]].head(10).to_string(index=False))