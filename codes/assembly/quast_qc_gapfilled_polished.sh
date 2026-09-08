#!/bin/bash
#SBATCH --job-name=quast_qc_gapfilled_polished
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=short
#SBATCH --output=logs/quast_qc_gapfilled_polished_%j.out
#SBATCH --error=logs/quast_qc_gapfilled_polished_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Targeted post-gap-fill polishing QC (job 716452, completed 2026-09-06):
# NextPolish re-run on the TGS-GapCloser gap-filled genome specifically, to
# reconcile the newly-filled Nanopore-derived sequence against the Illumina
# short reads (unlike the original whole-genome NextPolish attempt, which
# found nothing new to correct - see "NextPolish - No Improvement" Known
# Issue - this one has new sequence that was never Illumina-corrected).
#
# Compare against BOTH prior stages:
#   - assembly/qc_results/quast_gapfilled/ (pre-polish gap-filled baseline:
#     genome fraction 92.088%, 23532 misassemblies - see "TGS-GapCloser -
#     Result" Known Issue)
#   - assembly/qc_results/quast/ (original pre-gapfill baseline: genome
#     fraction 82.807%, 7589 misassemblies)
# to see whether this targeted polish recovers any of the structural
# precision lost during gap-filling, on top of the completeness gain.
#
# Same fixes as prior QUAST stages: --fragmented, Chr0 excluded
# (genome.nextpolish.noChr0.fasta - NextPolish renamed it to
# Chr0_RagTag_np1212, appending the task sequence "1212" during polishing).

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SCAFFOLD=${PROJECT_DIR}/assembly/nextpolish_output_gapfilled/genome.nextpolish.noChr0.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results

mkdir -p "${QC_DIR}/quast_gapfilled_polished" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: genome.nextpolish.noChr0.fasta not found at $SCAFFOLD - did nextpolish_gapfilled_genome.sh finish and get Chr0-filtered?"
    exit 1
fi

module load quast/5.0.2
echo "Start time: $(date)"
quast.py "$SCAFFOLD" -r "$REF" -o "${QC_DIR}/quast_gapfilled_polished" -t "${SLURM_CPUS_PER_TASK}" --fragmented
echo "End time: $(date)"

echo "=== Summary ==="
cat "${QC_DIR}/quast_gapfilled_polished/report.txt" 2>/dev/null | head -30
