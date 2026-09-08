#!/bin/bash
#SBATCH --job-name=haplotype_trimmomatic
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/haplotype_trimmomatic_%A_%a.out
#SBATCH --error=logs/haplotype_trimmomatic_%A_%a.err
#SBATCH --partition=medium
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INTERVALS_FILE=${INTERVALS}
INPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/markdup
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/gvcf

mkdir -p "$OUTPUT_DIR" logs/

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
echo "Using TMPDIR: $TMPDIR"
echo "Available space: $(df -h /tmp | awk 'NR==2 {print $4}')"

module load gatk4/4.4.0.0

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "Sample:     $SAMPLE"
echo "Start time: $(date)"

BAM="${INPUT_DIR}/${SAMPLE}.markdup.bam"
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz" \
    -L "$INTERVALS_FILE" \
    -ERC GVCF \
    --sample-name "$SAMPLE" \
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    -ploidy 2 \
    --tmp-dir "$TMPDIR"

echo "End time: $(date)"
echo "Done: $SAMPLE"
