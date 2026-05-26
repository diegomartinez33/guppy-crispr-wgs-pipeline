#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegomartinez3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INPUT_DIR=${PROJECT_DIR}/raw_fastq
OUTPUT_DIR=${PROJECT_DIR}/trimmed_fastp
LOG_DIR=${PROJECT_DIR}/logs

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ── Activate conda ───────────────────────────────────────────────────────────
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastp_env

# ── Sample for this task ─────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "Processing sample: $SAMPLE"
echo "Start time: $(date)"

# ── Run fastp ────────────────────────────────────────────────────────────────
fastp \
  -i  "${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz" \
  -I  "${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz" \
  -o  "${OUTPUT_DIR}/${SAMPLE}_R1_filtered.fastq.gz" \
  -O  "${OUTPUT_DIR}/${SAMPLE}_R2_filtered.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 12 \
  --length_required 50 \
  --thread "${SLURM_CPUS_PER_TASK}" \
  --html "${LOG_DIR}/${SAMPLE}_fastp.html" \
  --json "${LOG_DIR}/${SAMPLE}_fastp.json"

echo "End time: $(date)"
echo "Done: $SAMPLE"