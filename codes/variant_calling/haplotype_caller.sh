#!/bin/bash
#SBATCH --job-name=haplotype_trimmomatic
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00          # paso más lento del pipeline
#SBATCH --output=logs/haplotype_trimmomatic_%A_%a.out
#SBATCH --error=logs/haplotype_trimmomatic_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/markdup  # cambiar a /gatk/trimmomatic/markdup
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/gvcf
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna

mkdir -p "$OUTPUT_DIR" logs/

# ── Configurar TMPDIR ─────────────────────────────────────────────────────────
TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT
echo "Using TMPDIR: $TMPDIR"
echo "Available space: $(df -h /tmp | awk 'NR==2 {print $4}')"

# ── Cargar módulos ────────────────────────────────────────────────────────────
module load gatk4/4.4.0.0

# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "Sample:     $SAMPLE"
echo "Start time: $(date)"

# ── HaplotypeCaller ───────────────────────────────────────────────────────────
gatk HaplotypeCaller \
  -R "$REF" \
  -I "${INPUT_DIR}/${SAMPLE}.markdup.bam" \
  -O "${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz" \
  -ERC GVCF \
  --sample-name "${SAMPLE}" \
  --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
  -ploidy 2 \
  --tmp-dir "$TMPDIR"

echo "End time: $(date)"
echo "Done: $SAMPLE"