#!/bin/bash
#SBATCH --job-name=quast_qc_gapfilled
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=short
#SBATCH --output=logs/quast_qc_gapfilled_%j.out
#SBATCH --error=logs/quast_qc_gapfilled_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Roadmap option 2 QC: Assembly quality assessment (QUAST) on the
# TGS-GapCloser gap-filled genome (job 710348 - 56.3% of gap regions
# filled, 179,359/318,572), for direct comparison against the ORIGINAL
# pre-polish baseline (assembly/qc_results/quast/ - genome fraction
# 82.807%, N50 28.3Mb, duplication ratio 1.078, 7589 misassemblies) - NOT
# the NextPolish-polished version, which was rejected (see "NextPolish —
# No Improvement" Known Issue).
#
# Same fixes as quast_qc.sh/quast_qc_polished.sh: --fragmented (both REF
# and our scaffold are fragmented), Chr0 excluded
# (colombian_gapfilled.noChr0.fasta - TGS-GapCloser preserves original
# sequence names, unlike NextPolish, so it's still "Chr0_RagTag" here).

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SCAFFOLD=${PROJECT_DIR}/assembly/tgsgapcloser_output/colombian_gapfilled.noChr0.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results

mkdir -p "${QC_DIR}/quast_gapfilled" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: colombian_gapfilled.noChr0.fasta not found at $SCAFFOLD - did tgsgapcloser_genome.sh finish and get Chr0-filtered?"
    exit 1
fi

module load quast/5.0.2
echo "Start time: $(date)"
quast.py "$SCAFFOLD" -r "$REF" -o "${QC_DIR}/quast_gapfilled" -t "${SLURM_CPUS_PER_TASK}" --fragmented
echo "End time: $(date)"

echo "=== Summary ==="
cat "${QC_DIR}/quast_gapfilled/report.txt" 2>/dev/null | head -30
