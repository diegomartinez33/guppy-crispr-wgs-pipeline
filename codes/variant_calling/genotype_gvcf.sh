#!/bin/bash
#SBATCH --job-name=genotype_trimmomatic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=logs/genotype_trimmomatic.out
#SBATCH --error=logs/genotype_trimmomatic.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
INPUT_DB=${PROJECT_DIR}/gatk/trimmomatic/genomicsdb
OUTPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/vcf
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

echo "Start time: $(date)"

# ── GenotypeGVCFs ─────────────────────────────────────────────────────────────
gatk GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://${INPUT_DB}" \
  -O "${OUTPUT_DIR}/all_samples.vcf.gz" \
  --tmp-dir "$TMPDIR"

# ── create indexes ───────────────────────────────────────────────────────────

GVCF_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data/gatk/trimmomatic/gvcf

ls ${GVCF_DIR}/*.g.vcf.gz | wc -l      # debe dar 15
ls ${GVCF_DIR}/*.g.vcf.gz.tbi | wc -l  # debe dar 15

# Si faltan índices
module load gatk4/4.4.0.0
for gvcf in ${GVCF_DIR}/*.g.vcf.gz; do
    if [ ! -f "${gvcf}.tbi" ]; then
        echo "Indexing: $gvcf"
        gatk IndexFeatureFile -I "$gvcf"
    fi
done

echo "End time: $(date)"