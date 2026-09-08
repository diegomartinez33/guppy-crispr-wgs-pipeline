#!/bin/bash
#SBATCH --job-name=busco_qc_polished
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/busco_qc_polished_%j.out
#SBATCH --error=logs/busco_qc_polished_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 5 QC: Assembly completeness assessment (BUSCO) on the
# NextPolish-corrected genome, for direct comparison against the pre-polish
# baseline (assembly/qc_results/busco/ - C:87.1% [S:86.2%,D:0.9%], F:6.2%,
# M:6.7%). This is the more important comparison, since polishing targeted
# the 220 "complete" BUSCO genes with internal stop codons (indel/frameshift
# artifacts) - watch for E% (erroneous) and the split between S/D changing.
#
# Same odb10 fix as busco_qc.sh - see that script + CLAUDE.md Known Issues
# ("BUSCO — Hardcoded odb10 Version Check") for why odb12.2 cannot be used
# with this BUSCO 5.7.1 install.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCAFFOLD=${PROJECT_DIR}/assembly/nextpolish_output/genome.nextpolish.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
BUSCO_DL=${QC_DIR}/busco_downloads

mkdir -p "${QC_DIR}/busco_polished" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: genome.nextpolish.fasta not found at $SCAFFOLD - did Phase 5 (nextpolish_genome.sh) finish?"
    exit 1
fi

module load busco/5.7.1
cd "${QC_DIR}/busco_polished"

echo "Start time: $(date)"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb10 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o busco_colombian_scaffold_polished
echo "End time: $(date)"

echo "=== Summary ==="
find "${QC_DIR}/busco_polished/busco_colombian_scaffold_polished" -name "short_summary*.txt" -exec cat {} \;
