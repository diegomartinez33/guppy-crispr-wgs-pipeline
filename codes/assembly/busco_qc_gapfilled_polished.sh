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
source "${PROJECT_DIR}/codes/genome_versions.sh"
SCAFFOLD=${PROJECT_DIR}/assembly/nextpolish_output_gapfilled${OUT_SUFFIX}/genome.nextpolish.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
# BUSCO_DL (actinopterygii_odb10 lineage data) is species/version-independent
# - shared across both REF_VERSIONs, no need to re-download.
BUSCO_DL=${QC_DIR}/busco_downloads
RUN_DIR=${QC_DIR}/busco_gapfilled_polished${OUT_SUFFIX}

mkdir -p "$RUN_DIR" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: genome.nextpolish.fasta not found at $SCAFFOLD - did nextpolish_gapfilled_genome.sh finish?"
    exit 1
fi

module load busco/5.7.1
cd "$RUN_DIR"

echo "Start time: $(date)"
echo "REF_VERSION=${REF_VERSION}  SCAFFOLD=${SCAFFOLD}"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb10 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o "busco_colombian_scaffold_gapfilled_polished${OUT_SUFFIX}"
echo "End time: $(date)"

echo "=== Summary ==="
find "${RUN_DIR}/busco_colombian_scaffold_gapfilled_polished${OUT_SUFFIX}" -name "short_summary*.txt" -exec cat {} \;
