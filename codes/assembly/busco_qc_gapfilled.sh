#!/bin/bash
#SBATCH --job-name=busco_qc_gapfilled
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/busco_qc_gapfilled_%j.out
#SBATCH --error=logs/busco_qc_gapfilled_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Roadmap option 2 QC: Assembly completeness assessment (BUSCO) on the
# TGS-GapCloser gap-filled genome (job 710348 - 56.3% of gap regions
# filled), for direct comparison against the ORIGINAL baseline
# (assembly/qc_results/busco/ - C:87.1% [S:86.2%,D:0.9%], F:6.2%, M:6.7%,
# 220 genes w/ internal stop codons). Unlike NextPolish (reprocessed the
# same short reads, no new information, no improvement - see "NextPolish —
# No Improvement" Known Issue), gap-filling adds genuinely new sequence
# from long reads that can span what short reads/RagTag gaps couldn't -
# watch specifically for Missing/Fragmented BUSCOs moving to Complete.
#
# Same odb10 fix as busco_qc.sh/busco_qc_polished.sh - see CLAUDE.md Known
# Issues ("BUSCO — Hardcoded odb10 Version Check").

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCAFFOLD=${PROJECT_DIR}/assembly/tgsgapcloser_output/colombian_gapfilled.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
BUSCO_DL=${QC_DIR}/busco_downloads

mkdir -p "${QC_DIR}/busco_gapfilled" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: colombian_gapfilled.fasta not found at $SCAFFOLD - did tgsgapcloser_genome.sh finish?"
    exit 1
fi

module load busco/5.7.1
cd "${QC_DIR}/busco_gapfilled"

echo "Start time: $(date)"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb10 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o busco_colombian_scaffold_gapfilled
echo "End time: $(date)"

echo "=== Summary ==="
find "${QC_DIR}/busco_gapfilled/busco_colombian_scaffold_gapfilled" -name "short_summary*.txt" -exec cat {} \;
