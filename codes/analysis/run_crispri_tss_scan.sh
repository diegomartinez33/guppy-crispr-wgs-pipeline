#!/bin/bash
# Compare CRISPRi TSS-window guide candidates between the NCBI Trinidad
# reference and the Colombian population genome, per gene. Mirrors
# run_ko_guide_scan.sh's structure/conventions exactly (see that script's
# header for the "why pseudogenome by default" rationale) but calls
# crispri_tss_scan.py instead of ko_guide_scan.py - added 2026-09-20 since
# crispri_tss_scan.py previously had no SLURM/bash wrapper of its own.
#
# ── MODIFICAR AQUÍ para analizar otros genes ────────────────────────────
GENES=(bdnf agap3 grin1a grin1b gria1a gria1b gria2b nlgn1)
POPULATION=${POPULATION:-pseudogenome}   # "pseudogenome" (recomendado), "scaffolded",
                          # o "pseudogenome_v2" (ver run_ko_guide_scan.sh)
# ─────────────────────────────────────────────────────────────────────────

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCRIPT_DIR=${PROJECT_DIR}/codes/analysis
OUT_DIR=${PROJECT_DIR}/analysis/ko_guide_scan

mkdir -p "$OUT_DIR" logs/

module load minimap2
module load samtools/1.16.1
module load singularity/3.7.1

echo "Start time: $(date)"
echo "Population genome: $POPULATION"
echo "Genes: ${GENES[*]}"
echo ""

for GENE in "${GENES[@]}"; do
    echo "=================================================="
    python3 "${SCRIPT_DIR}/crispri_tss_scan.py" --gene "$GENE" --population "$POPULATION"
    echo ""
done

echo "=================================================="
echo "All results in: $OUT_DIR"
echo "End time: $(date)"
