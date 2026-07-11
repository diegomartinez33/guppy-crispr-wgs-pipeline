#!/bin/bash
#SBATCH --job-name=haplo_scatter
#SBATCH --array=1-360%40
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --output=logs/haplo_scatter_%A_%a.out
#SBATCH --error=logs/haplo_scatter_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=FAIL,ARRAY_TASKS

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INTERVALS_FILE=${PROJECT_DIR}/reference/intervals.list
INPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/markdup
OUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/gvcf_scatter
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna

N_CHROMS=24

SAMPLE_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / N_CHROMS + 1 ))
CHROM_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) % N_CHROMS + 1 ))

SAMPLE=$(sed -n "${SAMPLE_IDX}p" "$SAMPLE_LIST")
CHROM=$(sed -n "${CHROM_IDX}p" "$INTERVALS_FILE")

mkdir -p "${OUT_DIR}/${SAMPLE}" logs/

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT

module load gatk4/4.4.0.0

echo "Task ${SLURM_ARRAY_TASK_ID}: sample=${SAMPLE}  chrom=${CHROM}"
echo "Start: $(date)"

BAM="${INPUT_DIR}/${SAMPLE}.markdup.bam"
OUT_GVCF="${OUT_DIR}/${SAMPLE}/${CHROM}.g.vcf.gz"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM" \
    -O "$OUT_GVCF" \
    -L "$CHROM" \
    -ERC GVCF \
    --sample-name "$SAMPLE" \
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    -ploidy 2 \
    --tmp-dir "$TMPDIR"

echo "End: $(date)"
echo "Done: ${SAMPLE} ${CHROM}"
