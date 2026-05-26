#!/bin/bash
#SBATCH --job-name=reindex_gvcf
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/reindex_gvcf_%A_%a.out
#SBATCH --error=logs/reindex_gvcf_%A_%a.err
#SBATCH --partition=short

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
GVCF_DIR=${PROJECT_DIR}/gatk/trimmomatic/gvcf
SAMPLE_LIST=${PROJECT_DIR}/samples.txt

mkdir -p logs/

module load samtools/1.16.1  # incluye bgzip y tabix
module load gatk4/4.4.0.0

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"

echo "Sample:     $SAMPLE"
echo "Start time: $(date)"

# ── Verificar que el GVCF existe ──────────────────────────────────────────────
if [ ! -f "$GVCF" ]; then
    echo "ERROR: GVCF not found: $GVCF"
    exit 1
fi

# ── Re-comprimir con bgzip ────────────────────────────────────────────────────
echo "Re-compressing with bgzip..."
zcat "$GVCF" | bgzip -c > "${GVCF}.bgz.tmp"

# Verificar que el nuevo archivo es válido
if [ ! -s "${GVCF}.bgz.tmp" ]; then
    echo "ERROR: bgzip compression failed"
    exit 1
fi

# Reemplazar el archivo original
mv "${GVCF}.bgz.tmp" "$GVCF"

echo "Verifying bgzip format..."
file "$GVCF"

# ── Indexar con tabix ─────────────────────────────────────────────────────────
echo "Indexing with tabix..."
tabix -p vcf "$GVCF"

if [ -f "${GVCF}.tbi" ]; then
    echo "✅ Index created: ${GVCF}.tbi"
else
    echo "❌ Index failed — trying with GATK..."
    gatk IndexFeatureFile -I "$GVCF"
fi

echo "End time: $(date)"
echo "Done: $SAMPLE"