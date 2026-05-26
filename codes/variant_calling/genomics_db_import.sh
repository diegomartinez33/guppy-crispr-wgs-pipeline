#!/bin/bash
#SBATCH --job-name=genomicsdb_trimmomatic
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/genomicsdb_trimmomatic.out
#SBATCH --error=logs/genomicsdb_trimmomatic.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
GVCF_DIR=${PROJECT_DIR}/gatk/trimmomatic/gvcf
OUTPUT_DB=${PROJECT_DIR}/gatk/trimmomatic/genomicsdb
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INTERVALS_FILE=${PROJECT_DIR}/reference/intervals.list

mkdir -p "$OUTPUT_DB" logs/

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

# ── Verificar índices .tbi ────────────────────────────────────────────────────
echo "=== Verifying GVCF indices ==="
MISSING_TBI=0
while IFS= read -r SAMPLE; do
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    TBI="${GVCF}.tbi"
    if [ ! -f "$TBI" ]; then
        echo "❌ Missing index: $TBI"
        MISSING_TBI=1
    else
        echo "✅ $(basename $SAMPLE)"
    fi
done < "$SAMPLE_LIST"

if [ "$MISSING_TBI" -eq 1 ]; then
    echo "ERROR: Some .tbi indices are missing. Run reindex_gvcf.sh first."
    exit 1
fi

# ── Crear archivo de intervalos si no existe ──────────────────────────────────
if [ ! -f "$INTERVALS_FILE" ]; then
    echo "Creating intervals file..."
    grep "^>" "$REF" \
      | sed 's/^>//' \
      | awk '{print $1}' \
      | head -24 \
      > "$INTERVALS_FILE"
fi

echo "Intervals:"
cat "$INTERVALS_FILE"

# ── Construir argumentos -V ───────────────────────────────────────────────────
V_ARGS=""
while IFS= read -r SAMPLE; do
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    if [ ! -f "$GVCF" ]; then
        echo "WARNING: GVCF not found: $GVCF"
        continue
    fi
    V_ARGS="${V_ARGS} -V ${GVCF}"
done < "$SAMPLE_LIST"

echo "Number of GVCFs: $(echo $V_ARGS | tr ' ' '\n' | grep -c '.g.vcf.gz')"

# ── Eliminar genomicsdb previo si existe ──────────────────────────────────────
if [ -d "$OUTPUT_DB" ] && [ "$(ls -A $OUTPUT_DB)" ]; then
    echo "WARNING: Removing existing genomicsdb..."
    rm -rf "$OUTPUT_DB"
    mkdir -p "$OUTPUT_DB"
fi

# ── Eliminar genomicsdb previo si existe ──────────────────────────────────────
if [ -d "$OUTPUT_DB" ]; then
    echo "Removing existing genomicsdb workspace..."
    rm -rf "$OUTPUT_DB"
    echo "✅ Removed: $OUTPUT_DB"
fi
# NO crear el directorio — GATK lo crea internamente

# ── GenomicsDBImport ──────────────────────────────────────────────────────────
echo "Running GenomicsDBImport..."
gatk GenomicsDBImport \
  $V_ARGS \
  --genomicsdb-workspace-path "$OUTPUT_DB" \
  -L "$INTERVALS_FILE" \
  --reader-threads "${SLURM_CPUS_PER_TASK}" \
  --tmp-dir "$TMPDIR"

echo "End time: $(date)"