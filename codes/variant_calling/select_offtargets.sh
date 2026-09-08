#!/bin/bash
#SBATCH --job-name=select_offtargets
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --output=logs/select_offtargets.out
#SBATCH --error=logs/select_offtargets.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
INPUT_VCF=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/vcf/all_samples.vcf.gz
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/vcf_offtargets
OFFTARGET_INTERVALS=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/offtargets/combined/offtargets_intervals.list

mkdir -p "$OUTPUT_DIR" logs/

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT

module load gatk4/4.4.0.0

echo "Start time: $(date)"

# Extract all variants at off-target sites
gatk SelectVariants \
  -R "$REF" \
  -V "$INPUT_VCF" \
  -L "$OFFTARGET_INTERVALS" \
  -O "${OUTPUT_DIR}/offtarget_variants.vcf.gz" \
  --tmp-dir "$TMPDIR"
echo "✅ Off-target variants selected"

# Extract only INDELs at off-target sites (most relevant for CRISPR)
gatk SelectVariants \
  -R "$REF" \
  -V "$INPUT_VCF" \
  -L "$OFFTARGET_INTERVALS" \
  --select-type-to-include INDEL \
  -O "${OUTPUT_DIR}/offtarget_indels.vcf.gz" \
  --tmp-dir "$TMPDIR"
echo "✅ Off-target INDELs selected"

echo "End time: $(date)"