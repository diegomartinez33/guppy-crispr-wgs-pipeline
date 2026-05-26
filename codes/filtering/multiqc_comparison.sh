#!/bin/bash
#SBATCH --job-name=multiqc_comparison
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=20:00:00
#SBATCH --output=logs/multiqc.out
#SBATCH --error=logs/multiqc.err
#SBATCH --partition=short
#SBATCH --mail-user=diegomartinez3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
RAW_FASTQC=${PROJECT_DIR}/raw_fastq/fastqc_results
TRIMMOMATIC_FASTQC=${PROJECT_DIR}/trimmed_trimmomatic/fastqc_results
FASTP_FASTQC=${PROJECT_DIR}/trimmed_fastp/fastqc_results
TRIMMOMATIC_LOGS=${PROJECT_DIR}/trimmed_trimmomatic/logs
FASTP_LOGS=${PROJECT_DIR}/trimmed_fastp/logs
OUTPUT_DIR=${PROJECT_DIR}/multiqc_report

mkdir -p "$OUTPUT_DIR"

# ── Cargar módulo o activar conda ─────────────────────────────────────────────
module load multiqc 2>/dev/null || {
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate fastp_env
    #pip install multiqc -q
}

echo "Start time: $(date)"

# ── Correr MultiQC ────────────────────────────────────────────────────────────
multiqc \
  "$RAW_FASTQC" \
  "$TRIMMOMATIC_FASTQC" \
  "$FASTP_FASTQC" \
  "$TRIMMOMATIC_LOGS" \
  "$FASTP_LOGS" \
  --outdir "$OUTPUT_DIR" \
  --filename multiqc_trimming_comparison \
  --dirs \
  --dirs-depth 2 \
  --title "Trimming Comparison: Raw vs Trimmomatic vs fastp"

echo "End time: $(date)"
echo "MultiQC report: ${OUTPUT_DIR}/multiqc_trimming_comparison.html"