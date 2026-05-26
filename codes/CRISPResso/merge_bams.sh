#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/merge_bams.out
#SBATCH --error=logs/merge_bams.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
BAM_DIR=${PROJECT_DIR}/gatk/trimmomatic/markdup
MERGE_DIR=${PROJECT_DIR}/mapping/trimmomatic/merged

mkdir -p "$MERGE_DIR" logs/

module load samtools/1.16.1

echo "Start time: $(date)"

# ── Función: merge + sort + index ─────────────────────────────────────────────
merge_group() {
    local GROUP_NAME=$1
    shift
    local BAMS=("$@")

    echo "=== Merging group: $GROUP_NAME ==="
    echo "Input BAMs:"
    for b in "${BAMS[@]}"; do echo "  $b"; done

    SORTED_BAM="${MERGE_DIR}/${GROUP_NAME}_merged.sorted.bam"
    TMP_BAM="${MERGE_DIR}/${GROUP_NAME}_merged_tmp.bam"

    # Verificar que todos los BAMs existen
    for b in "${BAMS[@]}"; do
        if [ ! -f "$b" ]; then
            echo "ERROR: BAM not found: $b"
            return 1
        fi
    done

    # Merge
    samtools merge -f -@ "${SLURM_CPUS_PER_TASK}" "$TMP_BAM" "${BAMS[@]}"

    # Sort
    samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$SORTED_BAM" "$TMP_BAM"

    # Index
    samtools index "$SORTED_BAM"

    # Limpiar temporal
    rm -f "$TMP_BAM"

    # Estadísticas
    echo "✅ $GROUP_NAME merged successfully"
    echo "   Output: $SORTED_BAM"
    echo "   Total reads: $(samtools view -c $SORTED_BAM)"
    echo "   File size:   $(du -sh $SORTED_BAM | cut -f1)"
    echo ""
}

# ── Merge por grupo ───────────────────────────────────────────────────────────
merge_group "Control" \
    "${BAM_DIR}/Control_MNP_I_S54_L002.markdup.bam" \
    "${BAM_DIR}/Control_MNP_II_S55_L002.markdup.bam" \
    "${BAM_DIR}/Control_MNP_III_S56_L002.markdup.bam"

merge_group "RNP_Cas" \
    "${BAM_DIR}/RNP_Cas1_S65_L002.markdup.bam" \
    "${BAM_DIR}/RNP_Cas2_S66_L002.markdup.bam" \
    "${BAM_DIR}/RNP_Cas3_S67_L002.markdup.bam" \
    "${BAM_DIR}/RNP_Cas4_S68_L002.markdup.bam"

merge_group "Plasmid_Ko" \
    "${BAM_DIR}/Plasmid_Ko_P1_S61_L002.markdup.bam" \
    "${BAM_DIR}/Plasmid_Ko_P2_S62_L002.markdup.bam" \
    "${BAM_DIR}/Plasmid_Ko_P3_S63_L002.markdup.bam" \
    "${BAM_DIR}/Plasmid_Ko_P4_S64_L002.markdup.bam"

merge_group "Only_MNP" \
    "${BAM_DIR}/Only_MNP_C1_S57_L002.markdup.bam" \
    "${BAM_DIR}/Only_MNP_C2_S58_L002.markdup.bam" \
    "${BAM_DIR}/Only_MNP_C3_S59_L002.markdup.bam" \
    "${BAM_DIR}/Only_MNP_C4_S60_L002.markdup.bam"

# ── Verificación final ────────────────────────────────────────────────────────
echo "=== Final verification ==="
for GROUP in Control RNP_Cas Plasmid_Ko Only_MNP; do
    BAM="${MERGE_DIR}/${GROUP}_merged.sorted.bam"
    BAI="${MERGE_DIR}/${GROUP}_merged.sorted.bam.bai"
    if [ -f "$BAM" ] && [ -f "$BAI" ]; then
        echo "✅ $GROUP: $(samtools view -c $BAM) reads"
    else
        echo "❌ $GROUP: BAM or BAI missing"
    fi
done

echo "End time: $(date)"