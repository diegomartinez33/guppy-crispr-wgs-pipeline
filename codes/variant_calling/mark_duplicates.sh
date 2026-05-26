#!/bin/bash
#SBATCH --job-name=markdup_fastp
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=logs/markdup_fastp_%A_%a.out
#SBATCH --error=logs/markdup_fastp_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INPUT_DIR=${PROJECT_DIR}/mapping/trimmomatic        # cambiar a /mapping/trimmomatic para ese script
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/markdup  # cambiar a /gatk/trimmomatic/markdup
METRICS_DIR=${PROJECT_DIR}/gatk/trimmomatic/metrics

mkdir -p "$OUTPUT_DIR" "$METRICS_DIR" logs/

# ── Configurar TMPDIR ─────────────────────────────────────────────────────────
TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT
echo "Using TMPDIR: $TMPDIR"
echo "Available space: $(df -h /tmp | awk 'NR==2 {print $4}')"

# ── Cargar módulos ────────────────────────────────────────────────────────────
module load gatk4/4.4.0.0
module load samtools/1.16.1

# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "Sample:     $SAMPLE"
echo "Input BAM:  ${INPUT_DIR}/${SAMPLE}.sorted.bam"
echo "Output BAM: ${OUTPUT_DIR}/${SAMPLE}.markdup.bam"
echo "Start time: $(date)"

# ── Verificar que el input existe ────────────────────────────────────────────
if [ ! -f "${INPUT_DIR}/${SAMPLE}.sorted.bam" ]; then
    echo "ERROR: Input BAM not found: ${INPUT_DIR}/${SAMPLE}.sorted.bam"
    exit 1
fi

# ── MarkDuplicates ────────────────────────────────────────────────────────────
gatk MarkDuplicates \
  -I "${INPUT_DIR}/${SAMPLE}.sorted.bam" \
  -O "${OUTPUT_DIR}/${SAMPLE}.markdup.bam" \
  -M "${METRICS_DIR}/${SAMPLE}_dup_metrics.txt" \
  --REMOVE_DUPLICATES false \
  --VALIDATION_STRINGENCY SILENT \
  --TMP_DIR "$TMPDIR"

# ── Indexar BAM resultante ────────────────────────────────────────────────────
gatk BuildBamIndex \
  -I "${OUTPUT_DIR}/${SAMPLE}.markdup.bam" \
  --TMP_DIR "$TMPDIR"

  # ── Indexar BAM con samtools (más explícito que BuildBamIndex) ────────────────
samtools index \
  "${OUTPUT_DIR}/${SAMPLE}.markdup.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.markdup.bam.bai"

# ── Verificar output ──────────────────────────────────────────────────────────
if [ -f "${OUTPUT_DIR}/${SAMPLE}.markdup.bam" ]; then
    echo "✅ BAM created: $(du -sh ${OUTPUT_DIR}/${SAMPLE}.markdup.bam)"
else
    echo "❌ ERROR: Output BAM not created"
    exit 1
fi

if [ -f "${OUTPUT_DIR}/${SAMPLE}.markdup.bam.bai" ]; then
    echo "✅ BAI created: ${OUTPUT_DIR}/${SAMPLE}.markdup.bam.bai"
else
    echo "❌ ERROR: BAI index not created"
    exit 1
fi

echo "End time: $(date)"
echo "Done: $SAMPLE"