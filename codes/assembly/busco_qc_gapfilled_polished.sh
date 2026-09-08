#!/bin/bash
#SBATCH --job-name=busco_qc_gapfilled_polished
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/busco_qc_gapfilled_polished_%j.out
#SBATCH --error=logs/busco_qc_gapfilled_polished_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Targeted post-gap-fill polishing QC (job 716452, completed 2026-09-06) -
# see quast_qc_gapfilled_polished.sh for full rationale. Compare against
# assembly/qc_results/busco_gapfilled/ (pre-polish gap-filled baseline:
# C:? [S:?,D:?], F:?, M:?, watch specifically whether the 220 genes with
# internal stop codons documented for the pre-gapfill assembly improve).
#
# BUSCO doesn't align against a reference, so Chr0 doesn't need excluding
# here (unlike QUAST) - full genome.nextpolish.fasta used, same as the
# other busco_qc*.sh stages. Same odb10 fix as prior BUSCO stages.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCAFFOLD=${PROJECT_DIR}/assembly/nextpolish_output_gapfilled/genome.nextpolish.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
BUSCO_DL=${QC_DIR}/busco_downloads

mkdir -p "${QC_DIR}/busco_gapfilled_polished" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: genome.nextpolish.fasta not found at $SCAFFOLD - did nextpolish_gapfilled_genome.sh finish?"
    exit 1
fi

module load busco/5.7.1
cd "${QC_DIR}/busco_gapfilled_polished"

echo "Start time: $(date)"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb10 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o busco_colombian_scaffold_gapfilled_polished
echo "End time: $(date)"

echo "=== Summary ==="
find "${QC_DIR}/busco_gapfilled_polished/busco_colombian_scaffold_gapfilled_polished" -name "short_summary*.txt" -exec cat {} \;
