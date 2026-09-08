#!/bin/bash
#SBATCH --job-name=varfilt_trimmomatic
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --output=logs/varfilt_trimmomatic.out
#SBATCH --error=logs/varfilt_trimmomatic.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
INPUT_VCF=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/vcf/all_samples.vcf.gz
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/vcf_filtered

mkdir -p "$OUTPUT_DIR" logs/

# ── Configurar TMPDIR ─────────────────────────────────────────────────────────
TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT
echo "Using TMPDIR: $TMPDIR"
echo "Available space: $(df -h /tmp | awk 'NR==2 {print $4}')"

module load gatk4/4.4.0.0

echo "Start time: $(date)"

# ── Separar SNPs e INDELs ─────────────────────────────────────────────────────
# SNPs
gatk SelectVariants \
  -R "$REF" \
  -V "$INPUT_VCF" \
  --select-type-to-include SNP \
  -O "${OUTPUT_DIR}/snps_raw.vcf.gz" \
  --tmp-dir "$TMPDIR"

# INDELs (más importantes para CRISPR)
gatk SelectVariants \
  -R "$REF" \
  -V "$INPUT_VCF" \
  --select-type-to-include INDEL \
  -O "${OUTPUT_DIR}/indels_raw.vcf.gz" \
  --tmp-dir "$TMPDIR"

# ── Filtrar SNPs ──────────────────────────────────────────────────────────────
gatk VariantFiltration \
  -R "$REF" \
  -V "${OUTPUT_DIR}/snps_raw.vcf.gz" \
  --filter-expression "QD < 2.0"     --filter-name "QD2" \
  --filter-expression "FS > 60.0"    --filter-name "FS60" \
  --filter-expression "MQ < 40.0"    --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "${OUTPUT_DIR}/snps_filtered.vcf.gz" \
  --tmp-dir "$TMPDIR"

# ── Filtrar INDELs ────────────────────────────────────────────────────────────
gatk VariantFiltration \
  -R "$REF" \
  -V "${OUTPUT_DIR}/indels_raw.vcf.gz" \
  --filter-expression "QD < 2.0"     --filter-name "QD2" \
  --filter-expression "FS > 200.0"   --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O "${OUTPUT_DIR}/indels_filtered.vcf.gz" \
  --tmp-dir "$TMPDIR"

echo "End time: $(date)"
echo "Filtered VCFs ready for CRISPResso2 analysis"