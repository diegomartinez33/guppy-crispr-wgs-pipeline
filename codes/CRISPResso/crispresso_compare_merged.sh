#!/bin/bash
#SBATCH --job-name=crispresso_compare
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crispresso_compare.out
#SBATCH --error=logs/crispresso_compare.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
ONTARGET_DIR=${PROJECT_DIR}/crispresso/ontarget/trimmomatic/merged
OUTPUT_DIR=${PROJECT_DIR}/crispresso/compare/trimmomatic/merged

mkdir -p "$OUTPUT_DIR" logs/

# ── Activar ambiente ──────────────────────────────────────────────────────────
source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
conda activate crispresso2_env

echo "Start time: $(date)"

# ── Función para correr CRISPRessoCompare ─────────────────────────────────────
run_compare() {
    local SAMPLE1_NAME=$1
    local SAMPLE2_NAME=$2
    local FOLDER1="${ONTARGET_DIR}/${SAMPLE1_NAME}/CRISPResso_on_${SAMPLE1_NAME}"
    local FOLDER2="${ONTARGET_DIR}/${SAMPLE2_NAME}/CRISPResso_on_${SAMPLE2_NAME}"
    local COMP_NAME="${SAMPLE2_NAME}_vs_${SAMPLE1_NAME}"

    echo "--- Comparing: $SAMPLE2_NAME vs $SAMPLE1_NAME ---"

    # Verificar que los folders existen
    if [ ! -d "$FOLDER1" ]; then
        echo "❌ ERROR: Folder not found: $FOLDER1"
        return 1
    fi
    if [ ! -d "$FOLDER2" ]; then
        echo "❌ ERROR: Folder not found: $FOLDER2"
        return 1
    fi

    CRISPRessoCompare \
      "$FOLDER1" \
      "$FOLDER2" \
      --sample_1_name "$SAMPLE1_NAME" \
      --sample_2_name "$SAMPLE2_NAME" \
      --output_folder "${OUTPUT_DIR}/${COMP_NAME}" \
      --name "$COMP_NAME"

    echo "✅ Done: $COMP_NAME"
    echo ""
}

# ── Comparaciones vs Control (referencia) ─────────────────────────────────────
echo "=== Control vs tratamientos ==="
run_compare "Control" "RNP_Cas"
run_compare "Control" "Plasmid_Ko"
run_compare "Control" "Only_MNP"

# ── Comparaciones entre tratamientos ─────────────────────────────────────────
echo "=== Entre tratamientos ==="
run_compare "RNP_Cas"    "Plasmid_Ko"
run_compare "RNP_Cas"    "Only_MNP"
run_compare "Plasmid_Ko" "Only_MNP"

echo "End time: $(date)"
echo ""
echo "=== Comparaciones completadas ==="
echo "Control vs tratamientos:"
echo "  ✅ RNP_Cas_vs_Control"
echo "  ✅ Plasmid_Ko_vs_Control"
echo "  ✅ Only_MNP_vs_Control"
echo "Entre tratamientos:"
echo "  ✅ Plasmid_Ko_vs_RNP_Cas"
echo "  ✅ Only_MNP_vs_RNP_Cas"
echo "  ✅ Only_MNP_vs_Plasmid_Ko"