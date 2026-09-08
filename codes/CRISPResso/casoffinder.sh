#!/bin/bash
#SBATCH --job-name=casoffinder
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/casoffinder.out
#SBATCH --error=logs/casoffinder.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Activar ambiente ──────────────────────────────────────────────────────────
# CONDA_BASE was referenced here before being defined (real bug, fixed
# 2026-09 while porting this script to REF_VERSION - see crispresso_wgs.sh
# for the correct pattern this now matches).
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
OUTPUT_DIR=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/offtargets

mkdir -p "$OUTPUT_DIR" logs/

echo "Start time: $(date)"
echo "Reference: $REF"
echo "Output dir: $OUTPUT_DIR"

# ── Crear archivo de input ────────────────────────────────────────────────────
cat > "${OUTPUT_DIR}/casoffinder_input.txt" << EOF
${REF}
NNNNNNNNNNNNNNNNNNNNNGG
TGAGAGACGCCCCGGGCATGNGG 4
EOF

# Verificar que el input se creó correctamente
echo "=== Input file ==="
cat "${OUTPUT_DIR}/casoffinder_input.txt"

# ── Correr Cas-OFFinder ───────────────────────────────────────────────────────
# ✅ Sin comentarios después de \ — todo limpio
cas-offinder "${OUTPUT_DIR}/casoffinder_input.txt" C "${OUTPUT_DIR}/offtargets.txt"

# ── Verificar resultados ──────────────────────────────────────────────────────
if [ -f "${OUTPUT_DIR}/offtargets.txt" ]; then
    echo "✅ Off-targets encontrados: $(wc -l < ${OUTPUT_DIR}/offtargets.txt)"
    echo ""
    echo "=== Distribución por mismatches ==="
    awk '{print $5}' "${OUTPUT_DIR}/offtargets.txt" | sort | uniq -c | sort -k2n
    echo ""
    echo "=== Preview primeros 10 resultados ==="
    head -10 "${OUTPUT_DIR}/offtargets.txt"
else
    echo "❌ ERROR: offtargets.txt no fue creado"
    exit 1
fi

echo "End time: $(date)"