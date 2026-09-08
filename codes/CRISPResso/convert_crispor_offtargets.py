#!/usr/bin/env python3
"""
Convierte el TSV de off-targets producido por crispor.py (contenedor
Singularity, ver crispor_offtarget_scan.sh) a BED y formato compatible con
Cas-OFFinder y CRISPRessoWGS. Filtra por la secuencia de la guía (no por
guideId - ver nota abajo).

Reemplaza el flujo anterior, que leía un .xls descargado manualmente de la
web de crispor.org (data/crispor_offtargets.xls) para el guideId "326forw" -
ese archivo era un artefacto de un solo uso para el genoma viejo, sin forma
de regenerarlo para un genoma nuevo sin repetir el paso manual. Ahora se usa
crispor.py corriendo dentro del contenedor Singularity ya configurado
(mismo genoma registrado que usa ko_guide_scan.py), lo que lo hace
reproducible por script para cualquier REF_VERSION - ver CLAUDE.md,
"Migration to GCF_904066995.2 (v2)".

Filtrado por guideSeq en vez de guideId: el guideId que asigna crispor.py
(ej. "8rev", "447rev") depende de la posición dentro del archivo FASTA de
entrada y de detalles internos de numeración - NO es estable entre
corridas distintas (input diferente, versión del contenedor, etc.). La
columna "guideSeq" del TSV de off-targets guarda la secuencia completa de
23pb (espaciador+PAM) del sitio on-target real, que sí es estable (misma
guía biológica, confirmada idéntica entre v1/v2 vía Cas-OFFinder).

Las columnas del TSV de salida de crispor.py (seqId, guideId, guideSeq,
offtargetSeq, mismatchPos, mismatchCount, mitOfftargetScore,
cfdOfftargetScore, chrom, start, end, strand, locusDesc) coinciden
exactamente con las que ya esperaba este script del .xls - por eso
combine_offtargets.py no necesita cambios de lógica, solo el
parametrizado de rutas.
"""

import os
import sys
import pandas as pd

PROJECT_DIR = "/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data"
sys.path.insert(0, os.path.join(PROJECT_DIR, "codes", "analysis"))

REF_VERSION = os.environ.get("REF_VERSION", "v1")
OUT_SUFFIX = "" if REF_VERSION == "v1" else f"_{REF_VERSION}"

CRISPOR_FILE = os.path.join(
    PROJECT_DIR, f"crispresso{OUT_SUFFIX}/offtargets/crispor_container/bdnf_offtargets.tsv"
)
OUTPUT_DIR = os.path.join(PROJECT_DIR, f"crispresso{OUT_SUFFIX}/offtargets/crispor")
GUIDE_TARGET_SEQ = "TGAGAGACGCCCCGGGCATGCGG"  # bdnf sgRNA spacer+PAM (23bp) at its actual
# genomic on-target site (confirmed via Cas-OFFinder, 0 mismatches, both v1/v2 - see
# CLAUDE.md). The offtargets TSV's "guideSeq" column stores this full 23bp on-target
# sequence (not just the 20bp spacer) for every off-target row belonging to that guide.
MAX_MISMATCHES = 4  # máximo mismatches a incluir

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Leer TSV de crispor.py ────────────────────────────────────────────────────
print(f"Leyendo archivo: {CRISPOR_FILE}")

df = pd.read_csv(CRISPOR_FILE, sep="\t")

print(f"Total de off-targets en archivo: {len(df)}")
print(f"Columnas encontradas: {list(df.columns)}")

# ── Filtrar por secuencia de la guía ──────────────────────────────────────────
df_filtered = df[df["guideSeq"] == GUIDE_TARGET_SEQ].copy()
print(f"Off-targets para guideSeq='{GUIDE_TARGET_SEQ}': {len(df_filtered)}")

# Filtrar por número de mismatches
df_filtered = df_filtered[df_filtered["mismatchCount"] <= MAX_MISMATCHES]
print(f"Off-targets con <= {MAX_MISMATCHES} mismatches: {len(df_filtered)}")

# ── Limpiar y estandarizar columnas ──────────────────────────────────────────
# Renombrar columnas para mayor claridad
df_filtered = df_filtered.rename(columns={
    "guideId":           "guide_id",
    "guideSeq":          "guide_seq",
    "offtargetSeq":      "offtarget_seq",
    "mismatchCount":     "mismatches",
    "mitOfftargetScore": "mit_score",
    "cfdOfftargetScore": "cfd_score",
    "chrom":             "chromosome",
    "start":             "start",
    "end":               "end",
    "strand":            "strand",
    "locusDesc":         "locus_description",
})

# Asegurar tipos correctos
df_filtered["start"] = pd.to_numeric(df_filtered["start"], errors="coerce").astype(int)
df_filtered["end"]   = pd.to_numeric(df_filtered["end"],   errors="coerce").astype(int)

# Ordenar por score MIT descendente
df_filtered = df_filtered.sort_values("mit_score", ascending=False)

# ── Output 1: CSV completo filtrado ───────────────────────────────────────────
CSV_OUT = os.path.join(OUTPUT_DIR, "crispor_bdnf_offtargets.csv")
df_filtered.to_csv(CSV_OUT, index=False)
print(f"\n✅ CSV completo: {CSV_OUT}")

# ── Output 2: Archivo BED para CRISPRessoWGS ─────────────────────────────────
# Formato: chr  start  end  name  score  strand
BED_OUT = os.path.join(OUTPUT_DIR, "crispor_bdnf_offtargets.bed")

with open(BED_OUT, "w") as f:
    f.write("# BED file generado desde CRISPOR (contenedor) para CRISPRessoWGS\n")
    f.write(f"# guideSeq: {GUIDE_TARGET_SEQ}\n")
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
CASOFF_OUT = os.path.join(OUTPUT_DIR, "crispor_bdnf_casoffinder_format.txt")

with open(CASOFF_OUT, "w") as f:
    f.write("# Converted from CRISPOR (container) to Cas-OFFinder format\n")
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
INTERVALS_OUT = os.path.join(OUTPUT_DIR, "crispor_bdnf_intervals.list")

with open(INTERVALS_OUT, "w") as f:
    for _, row in df_filtered.iterrows():
        f.write(f"{row['chromosome']}:{row['start']}-{row['end']}\n")

print(f"✅ GATK intervals: {INTERVALS_OUT}")

# ── Resumen estadístico ───────────────────────────────────────────────────────
print(f"\n=== Resumen para guideSeq={GUIDE_TARGET_SEQ} (REF_VERSION={REF_VERSION}) ===")
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