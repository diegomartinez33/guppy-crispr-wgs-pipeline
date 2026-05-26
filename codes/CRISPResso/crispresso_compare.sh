#!/bin/bash
#SBATCH --job-name=crispresso_compare
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crispresso_compare.out
#SBATCH --error=logs/crispresso_compare.err
#SBATCH --partition=your_partition

source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"
mamba activate crispresso2_env

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
ONTARGET_DIR=${PROJECT_DIR}/crispresso/ontarget/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/crispresso/compare

mkdir -p "$OUTPUT_DIR"

echo "Start time: $(date)"

# ── Comparar cada grupo vs Control_MNP ───────────────────────────────────────

# RNP_Cas1 vs Control
CRISPRessoCompare \
  --sample_1_name "Control_MNP_I" \
  --sample_2_name "RNP_Cas1" \
  "${ONTARGET_DIR}/Control_MNP_I_S54_L002" \
  "${ONTARGET_DIR}/RNP_Cas1_S65_L002" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas1_vs_Control" \
  --name "RNP_Cas1_vs_Control"

# RNP_Cas2 vs Control
CRISPRessoCompare \
  --sample_1_name "Control_MNP_II" \
  --sample_2_name "RNP_Cas2" \
  "${ONTARGET_DIR}/Control_MNP_II_S55_L002" \
  "${ONTARGET_DIR}/RNP_Cas2_S66_L002" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas2_vs_Control" \
  --name "RNP_Cas2_vs_Control"

# Plasmid_Ko vs Control
CRISPRessoCompare \
  --sample_1_name "Control_MNP_III" \
  --sample_2_name "Plasmid_Ko_P1" \
  "${ONTARGET_DIR}/Control_MNP_III_S56_L002" \
  "${ONTARGET_DIR}/Plasmid_Ko_P1_S61_L002" \
  --output_folder "${OUTPUT_DIR}/Plasmid_vs_Control" \
  --name "Plasmid_vs_Control"

echo "End time: $(date)"