#!/bin/bash
#SBATCH --job-name=merge_gvcfs
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=logs/merge_gvcfs_%A_%a.out
#SBATCH --error=logs/merge_gvcfs_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=FAIL,ARRAY_TASKS

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INTERVALS_FILE=${INTERVALS}
SCATTER_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/gvcf_scatter
GVCF_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/gvcf

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
mkdir -p "$GVCF_DIR" logs/

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT

module load gatk4/4.4.0.0

echo "Merging GVCFs for: $SAMPLE"
echo "Start: $(date)"

# Verify all per-chromosome GVCFs (per $INTERVALS_FILE) exist
MISSING=0
while IFS= read -r CHROM; do
    F="${SCATTER_DIR}/${SAMPLE}/${CHROM}.g.vcf.gz"
    if [ ! -f "$F" ]; then
        echo "ERROR: Missing scatter GVCF: $F"
        MISSING=1
    fi
done < "$INTERVALS_FILE"

if [ "$MISSING" -eq 1 ]; then
    echo "Aborting — not all scatter GVCFs are present."
    exit 1
fi

# Build -I arguments in chromosome order (matches intervals.list order)
I_ARGS=""
while IFS= read -r CHROM; do
    I_ARGS="${I_ARGS} -I ${SCATTER_DIR}/${SAMPLE}/${CHROM}.g.vcf.gz"
done < "$INTERVALS_FILE"

OUT_GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"

gatk MergeVcfs \
    $I_ARGS \
    -O "$OUT_GVCF" \
    --TMP_DIR "$TMPDIR"

gatk IndexFeatureFile -I "$OUT_GVCF"

echo "End: $(date)"
echo "Done: ${SAMPLE} → ${OUT_GVCF}"
