#!/bin/bash
#SBATCH --job-name=crispresso_merged
#SBATCH --array=1-4          # ← cambiar de 0-3 a 1-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/crispresso_merged_%A_%a.out
#SBATCH --error=logs/crispresso_merged_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
MERGE_DIR=${PROJECT_DIR}/mapping/trimmomatic/merged
OUTPUT_DIR=${PROJECT_DIR}/crispresso/ontarget/trimmomatic/merged
REGION="NC_024333.1:15922000-15922100"

mkdir -p "$OUTPUT_DIR" logs/

# ── Activar ambiente ──────────────────────────────────────────────────────────
source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
conda activate crispresso2_env

SAMTOOLS=$(which samtools)
AMPLICON_FWD=$(grep -v "^>" "${PROJECT_DIR}/reference/amplicon_bdnf_100bp.fa" | tr -d '\n')
AMPLICON_RC=$(grep -v "^>" "${PROJECT_DIR}/reference/amplicon_bdnf_100bp_rc.fa" | tr -d '\n')
SGRNA="TGAGAGACGCCCCGGGCATG"

echo "Amplicon FWD length: ${#AMPLICON_FWD}bp"
echo "Amplicon RC  length: ${#AMPLICON_RC}bp"
echo "Region: $REGION"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# ── Seleccionar grupo según array task (1-based) ──────────────────────────────
case $SLURM_ARRAY_TASK_ID in
    1) GROUP="Control" ;;
    2) GROUP="RNP_Cas" ;;
    3) GROUP="Plasmid_Ko" ;;
    4) GROUP="Only_MNP" ;;
    *) echo "ERROR: Invalid array task ID: $SLURM_ARRAY_TASK_ID"; exit 1 ;;
esac

BAM="${MERGE_DIR}/${GROUP}_merged.sorted.bam"

echo "Group:      $GROUP"
echo "Input BAM:  $BAM"
echo "Start time: $(date)"

# ── Verificar BAM ─────────────────────────────────────────────────────────────
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    echo "Available BAMs in $MERGE_DIR:"
    ls "$MERGE_DIR"/*.bam 2>/dev/null || echo "No BAMs found"
    exit 1
fi

# ── Filtrar región y duplicados ───────────────────────────────────────────────
FILTERED_BAM="${OUTPUT_DIR}/${GROUP}_filtered.bam"
SORTED_BAM="${OUTPUT_DIR}/${GROUP}_filtered.sorted.bam"

$SAMTOOLS view -b -F 0x400 -h "$BAM" "$REGION" -o "$FILTERED_BAM"

if [ ! -s "$FILTERED_BAM" ]; then
    echo "ERROR: Filtered BAM is empty"
    exit 1
fi

$SAMTOOLS sort -o "$SORTED_BAM" "$FILTERED_BAM"
$SAMTOOLS index "$SORTED_BAM"

echo "Reads in region: $($SAMTOOLS view -c $SORTED_BAM)"

# ── CRISPResso2 ───────────────────────────────────────────────────────────────
CRISPResso \
  --bam_input "$SORTED_BAM" \
  --amplicon_seq "$AMPLICON_FWD,$AMPLICON_RC" \
  --amplicon_name "bdnf_fwd,bdnf_rc" \
  --guide_seq "$SGRNA" \
  --bam_chr_loc "$REGION" \
  --output_folder "${OUTPUT_DIR}/${GROUP}" \
  --name "${GROUP}" \
  --n_processes "${SLURM_CPUS_PER_TASK}" \
  --min_frequency_alleles_around_cut_to_plot 0.05 \
  --expand_ambiguous_alignments \
  --place_report_in_output_folder \
  --write_cleaned_report

# ── Limpiar temporales ────────────────────────────────────────────────────────
rm -f "$FILTERED_BAM" "$SORTED_BAM" "${SORTED_BAM}.bai"

echo "End time: $(date)"
echo "Done: $GROUP"