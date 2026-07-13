#!/bin/bash
#SBATCH --job-name=quast_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/quast_qc_%j.out
#SBATCH --error=logs/quast_qc_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 3a: Assembly quality assessment (QUAST).
#
# Split out from the original combined qc_quast_busco.sh: job 683049 ran
# QUAST + BUSCO sequentially in one 8h job and timed out mid-QUAST (BUSCO
# never even started). Root cause was not a hang like the RagTag issue -
# minimap2 alignment itself finished in ~7min, but QUAST's single-threaded
# Python misassembly classification had to process a very large number of
# alignment segments, because BOTH the reference (2768 fragments) and our
# scaffold (RagTag join gaps, ~7220 N's/100kbp) are fragmented. QUAST's own
# log explicitly recommended --fragmented for this. Now runs independently
# of BUSCO (see busco_qc.sh) so neither blocks the other's time budget, with
# a generous 24h limit based on the observed slow-but-real progress rate.
# See "QUAST/BUSCO — slow contig analyzer on fragmented ref+scaffold" in
# CLAUDE.md Known Issues.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SCAFFOLD=${PROJECT_DIR}/assembly/ragtag_output/ragtag.scaffold.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results

mkdir -p "${QC_DIR}/quast" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: ragtag.scaffold.fasta not found at $SCAFFOLD - did Phase 2 (ragtag_scaffold.sh) finish?"
    exit 1
fi

module load quast/5.0.2
echo "Start time: $(date)"
quast.py "$SCAFFOLD" -r "$REF" -o "${QC_DIR}/quast" -t "${SLURM_CPUS_PER_TASK}" --fragmented
echo "End time: $(date)"

echo "=== Summary ==="
cat "${QC_DIR}/quast/report.txt" 2>/dev/null | head -30
