#!/bin/bash
#SBATCH --job-name=spades_coassembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=470G
#SBATCH --time=10-00:00:00
#SBATCH --partition=bigmem
#SBATCH --output=logs/spades_coassembly_%j.out
#SBATCH --error=logs/spades_coassembly_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 1: De novo co-assembly of the 3 Control replicates (SPAdes).
# Input FASTQs total ~47GB compressed - sized for bigmem given the project's
# prior experience underestimating GATK HaplotypeCaller resource needs.
# If this job times out, resume with:
#   spades.py --continue -o ${OUTPUT_DIR}
# instead of restarting from scratch (--checkpoints last enables this).
#
# --isolate is required here: combined depth across the 3 replicates is
# ~160-270x (47GB compressed reads vs a ~700Mb genome). Without --isolate,
# job 646213 crashed after 46min with an mmap ENOMEM at the k-mer counting
# step (BayesHammer tried to build an oversized index for the default
# single-cell-oriented code path). SPAdes's own params.txt output explicitly
# recommends --isolate for "high-coverage isolate and multi-cell data" -
# exactly this dataset - and this was missed in the initial draft.
#
# -m 450 (not 250): with --isolate, job 652298 ran ~24h and got as far as
# the K77 Distance Estimation stage before a genuine OOM (mimalloc ENOMEM,
# peak RSS ~249GB against a 250G limit with zero headroom). 450G leaves
# comfortable margin on the ~565GB bigmem nodes. SBATCH --mem is set higher
# still (470G) than spades.py's own -m 450 so SPAdes's internal self-limit
# triggers a graceful stop before the SLURM cgroup would hard-kill it.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
INPUT_DIR=${PROJECT_DIR}/trimmed_trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/assembly/spades_control_coassembly

mkdir -p "$OUTPUT_DIR" logs/

module load spades/4.0.0

echo "Start time: $(date)"
echo "Available space: $(df -h "$PROJECT_DIR" | awk 'NR==2 {print $4}')"

SAMPLES="Control_MNP_I_S54_L002 Control_MNP_II_S55_L002 Control_MNP_III_S56_L002"
for S in $SAMPLES; do
    for R in R1 R2; do
        F="${INPUT_DIR}/${S}_${R}_paired.fastq.gz"
        if [ ! -f "$F" ]; then
            echo "ERROR: missing input file: $F"
            exit 1
        fi
    done
done

spades.py \
    --pe1-1 "${INPUT_DIR}/Control_MNP_I_S54_L002_R1_paired.fastq.gz" \
    --pe1-2 "${INPUT_DIR}/Control_MNP_I_S54_L002_R2_paired.fastq.gz" \
    --pe2-1 "${INPUT_DIR}/Control_MNP_II_S55_L002_R1_paired.fastq.gz" \
    --pe2-2 "${INPUT_DIR}/Control_MNP_II_S55_L002_R2_paired.fastq.gz" \
    --pe3-1 "${INPUT_DIR}/Control_MNP_III_S56_L002_R1_paired.fastq.gz" \
    --pe3-2 "${INPUT_DIR}/Control_MNP_III_S56_L002_R2_paired.fastq.gz" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -m 450 \
    --isolate \
    --checkpoints last \
    -o "$OUTPUT_DIR"

echo "End time: $(date)"

CONTIGS="${OUTPUT_DIR}/contigs.fasta"
if [ -f "$CONTIGS" ]; then
    echo "Contigs: $(grep -c '^>' "$CONTIGS") sequences"
    ls -lh "$CONTIGS"
else
    echo "ERROR: contigs.fasta not produced"
    exit 1
fi
