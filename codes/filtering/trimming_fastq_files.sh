#!/bin/bash
#SBATCH --job-name=trimmomatic_trim
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=logs/trimmomatic_%A_%a.out
#SBATCH --error=logs/trimmomatic_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegomartinez3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INPUT_DIR=${PROJECT_DIR}/raw_fastq
OUTPUT_DIR=${PROJECT_DIR}/trimmed_trimmomatic
LOG_DIR=${PROJECT_DIR}/logs
ADAPTERS=/hpcfs/apps/conda4.12.0/envs/trimmomatic-0.39/share/trimmomatic/adapters/NexteraPE-PE.fa

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ── Cargar módulo ─────────────────────────────────────────────────────────────
module load trimmomatic    # ← ajusta al nombre exacto del módulo

# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "Processing sample: $SAMPLE"
echo "Trimmomatic version: $(trimmomatic -version)"
echo "Start time: $(date)"

# ── Run Trimmomatic ───────────────────────────────────────────────────────────
trimmomatic PE \
  -threads "${SLURM_CPUS_PER_TASK}" \
  -phred33 \
  "${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz" \
  "${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz" \
  "${OUTPUT_DIR}/${SAMPLE}_R1_paired.fastq.gz" \
  "${OUTPUT_DIR}/${SAMPLE}_R1_unpaired.fastq.gz" \
  "${OUTPUT_DIR}/${SAMPLE}_R2_paired.fastq.gz" \
  "${OUTPUT_DIR}/${SAMPLE}_R2_unpaired.fastq.gz" \
  ILLUMINACLIP:"${ADAPTERS}":2:30:10:8:keepBothReads \
  LEADING:2 \
  TRAILING:2 \
  SLIDINGWINDOW:4:12 \
  MINLEN:50 \
  2> "${LOG_DIR}/${SAMPLE}_trimmomatic.log"

echo "End time: $(date)"
echo "Done: $SAMPLE"