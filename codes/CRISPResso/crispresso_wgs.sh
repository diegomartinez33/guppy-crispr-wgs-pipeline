#!/bin/bash
#SBATCH --job-name=crispresso_wgs
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/crispresso_wgs_%A_%a.out
#SBATCH --error=logs/crispresso_wgs_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
BAM_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/markdup
OUTPUT_DIR=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/wgs/trimmomatic
REGION_FILE=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/offtargets/combined/offtargets_crispresso_wgs.bed

mkdir -p "$OUTPUT_DIR" logs/

# ── Activar ambiente ──────────────────────────────────────────────────────────
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"

which CRISPRessoWGS || { echo "❌ CRISPRessoWGS not found"; exit 1; }

# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BAM="${BAM_DIR}/${SAMPLE}.markdup.bam"

echo "Sample:     $SAMPLE"
echo "BAM:        $BAM"
echo "Start time: $(date)"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

# ── Verificar índice del BAM ──────────────────────────────────────────────────
if [ ! -f "${BAM}.bai" ]; then
    echo "Indexing BAM..."
    samtools index "$BAM"
fi

# ── CRISPRessoWGS ─────────────────────────────────────────────────────────────
CRISPRessoWGS \
  --bam_file "$BAM" \
  --region_file "$REGION_FILE" \
  --reference_file "$REF" \
  --output_folder "${OUTPUT_DIR}/${SAMPLE}" \
  --name "${SAMPLE}" \
  --n_processes "${SLURM_CPUS_PER_TASK}" \
  --min_reads_to_use_region 10 \
  --min_frequency_alleles_around_cut_to_plot 0.05 \
  --expand_ambiguous_alignments \
  --place_report_in_output_folder \
  --skip_failed

echo "End time: $(date)"
echo "Done: $SAMPLE"