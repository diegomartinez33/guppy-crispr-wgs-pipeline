#!/bin/bash
# Compare CRISPR knockout guide candidates (SpCas9 NGG) between the NCBI
# Trinidad reference and the Colombian population genome, per gene.
#
# ── MODIFICAR AQUÍ para analizar otros genes ────────────────────────────
GENES=(bdnf agap3 grin1a grin1b gria1a gria1b gria2b nlgn1)
POPULATION=pseudogenome   # "pseudogenome" (recomendado, ver abajo), "scaffolded",
                          # o "pseudogenome_v2" (contra el nuevo genoma de
                          # referencia GCF_904066995.2, una vez construido -
                          # ver CLAUDE.md, "Migration to GCF_904066995.2 (v2)")
# ─────────────────────────────────────────────────────────────────────────
#
# Por qué "pseudogenome" por defecto: preserva la estructura exón/intrón
# casi exactamente (mismas coordenadas, solo sustituciones SNP/indel), así
# que la comparación gen-por-gen es más confiable. "scaffolded" (ensamblaje
# de novo) puede tener reordenamientos o secuencia novedosa que el
# pseudogenoma no representa, pero también carga el riesgo de artefactos de
# ensamblaje (ver README de colombian_scaffolded_genome, sección
# Limitations, punto 1: ~5-7% de genes BUSCO "completos" tienen codones de
# parada internos, probablemente artefactos, no biología real) que podrían
# confundirse con variantes poblacionales genuinas.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCRIPT_DIR=${PROJECT_DIR}/codes/analysis
OUT_DIR=${PROJECT_DIR}/analysis/ko_guide_scan

mkdir -p "$OUT_DIR" logs/

# crispresso2_env's samtools is broken (missing libcrypto.so.1.0.0) and
# this script only needs Python stdlib (no pandas/pysam), so use the
# cluster's own samtools/minimap2 modules instead of any conda env.
# ORDER MATTERS: minimap2's module reloads an older anaconda base that
# shadows samtools/1.16.1's own libs if loaded second - load minimap2
# FIRST, then samtools/1.16.1, so the correct samtools wins the PATH.
module load minimap2
module load samtools/1.16.1
# For the optional CRISPOR (Singularity) scoring step - the pipeline
# degrades gracefully (prints a note, keeps the manual PAM scan results)
# if this module or the container/genomes aren't available.
module load singularity/3.7.1

echo "Start time: $(date)"
echo "Population genome: $POPULATION"
echo "Genes: ${GENES[*]}"
echo ""

for GENE in "${GENES[@]}"; do
    echo "=================================================="
    python3 "${SCRIPT_DIR}/ko_guide_scan.py" --gene "$GENE" --population "$POPULATION"
    echo ""
done

echo "=================================================="
echo "All results in: $OUT_DIR"
echo "End time: $(date)"
