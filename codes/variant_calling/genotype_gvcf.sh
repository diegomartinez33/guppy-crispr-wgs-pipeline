#!/bin/bash
#SBATCH --job-name=genotype_trimmomatic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=logs/genotype_trimmomatic.out
#SBATCH --error=logs/genotype_trimmomatic.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
INPUT_DB=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/genomicsdb
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/vcf

mkdir -p "$OUTPUT_DIR" logs/

# ── Configurar TMPDIR ─────────────────────────────────────────────────────────
TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
echo "Using TMPDIR: $TMPDIR"
echo "Available space: $(df -h /tmp | awk 'NR==2 {print $4}')"

# ── Cargar módulos ────────────────────────────────────────────────────────────
module load gatk4/4.4.0.0

echo "Start time: $(date)"

# ── GenotypeGVCFs ─────────────────────────────────────────────────────────────
gatk GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://${INPUT_DB}" \
  -O "${OUTPUT_DIR}/all_samples.vcf.gz" \
  --tmp-dir "$TMPDIR"

# ── create indexes ───────────────────────────────────────────────────────────

GVCF_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/gvcf

find "${GVCF_DIR}" -name "*.g.vcf.gz" | wc -l
find "${GVCF_DIR}" -name "*.g.vcf.gz.tbi" | wc -l

module load gatk4/4.4.0.0
for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
    if [ ! -f "${gvcf}.tbi" ]; then
        echo "Indexing: $gvcf"
        gatk IndexFeatureFile -I "$gvcf"
    fi
done

echo "End time: $(date)"