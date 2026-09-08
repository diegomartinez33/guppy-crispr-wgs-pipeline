#!/bin/bash
#SBATCH --job-name=quast_qc_polished
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=short
#SBATCH --output=logs/quast_qc_polished_%j.out
#SBATCH --error=logs/quast_qc_polished_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 5 QC: Assembly quality assessment (QUAST) on the NextPolish-corrected
# genome, for direct comparison against the pre-polish baseline
# (assembly/qc_results/quast/ - genome fraction 82.807%, N50 28.3Mb,
# duplication ratio 1.078, 7589 misassemblies).
#
# Same fixes as quast_qc.sh (see that script + CLAUDE.md Known Issues):
# --fragmented (both REF and our scaffold are fragmented), and Chr0
# excluded (genome.nextpolish.noChr0.fasta - NextPolish renamed it to
# Chr0_RagTag_np1212, appending the task sequence "1212" to every sequence
# name during its 2 polishing rounds).

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SCAFFOLD=${PROJECT_DIR}/assembly/nextpolish_output/genome.nextpolish.noChr0.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results

mkdir -p "${QC_DIR}/quast_polished" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: genome.nextpolish.noChr0.fasta not found at $SCAFFOLD - did Phase 5 (nextpolish_genome.sh) finish and get Chr0-filtered?"
    exit 1
fi

module load quast/5.0.2
echo "Start time: $(date)"
quast.py "$SCAFFOLD" -r "$REF" -o "${QC_DIR}/quast_polished" -t "${SLURM_CPUS_PER_TASK}" --fragmented
echo "End time: $(date)"

echo "=== Summary ==="
cat "${QC_DIR}/quast_polished/report.txt" 2>/dev/null | head -30
